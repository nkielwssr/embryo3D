"""Repères d'extrémités pointés par l'utilisateur sur un montage de captures vidéo 360° : vésicule optique = œil (cercle vert) et
pointe de la queue (point vert), une case par stade (étiquette de stade en haut à gauche). Chaque case est recalée sur sa
vidéo comme les étages manuels (marques_video), les deux marques sont ramenées dans le plan médian du volume et écrites dans
<dossier>/out/topographie/video360/reperes_utilisateur.json ; topographie_axe.py les lit en priorité : la pointe de queue ferme Co-4 ; l'œil est un repère crânien
(abscisse le long de l'axe), le bord otique des occipitaux restant calculé sur le label vesicules_otiques.

    python embryo3d/reperes_video.py embryo3d/queue-VO.png --grille 4x2 --stades CS13,CS14,,CS17,CS15,CS16,CS19,CS20
    python embryo3d/reperes_video.py --verifier          # apres topographie_axe.py : ecarts aux labels, images de controle, manifest

Vérification (par stade) : écart entre la marque « vésicule optique » et le label yeux (dans le plan sagittal, l'œil étant
latéral), écart entre la marque « queue » et l'extrémité caudale de l'enveloppe (bout géodésique) et du tube
neural ; image de contrôle (vidéo + projection sagittale avec labels, axe prolongé, occipital 1 et Co-4) ; ligne de tableau
dans manifest_topographie.json (controle.reperes_utilisateur) avec verdict (ok si les deux écarts <= 0,3 mm).
"""
import sys, os, json, argparse
import numpy as np, cv2
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from video360 import trouver_video, charger_images, silhouette, tour_complet, recalage_2d
from marques_video import image_la_plus_proche
from topographie_axe import bout_geodesique

DOSSIERS = {'CS13': 'CS13.f4v', 'CS14': 'CS14_f4v', 'CS15': 'CS15_f4v', 'CS16': 'CS16_f4v', 'CS17': 'CS17_f4v', 'CS19': 'CS19_f4v', 'CS20': 'CS20_f4v'}
SEUIL_MM = 0.3


def blobs(im):
    """Marques colorées avec leur taille (pixels) : le cercle sur la vésicule otique est le plus gros, le point de queue le plus petit."""
    b, g, r = [c.astype(int) for c in cv2.split(im)]
    sat = ((b - np.maximum(r, g)) > 60) | ((g - np.maximum(r, b)) > 60) | ((r - np.maximum(g, b)) > 60)
    lbl, n = ndi.label(sat); t = np.bincount(lbl.ravel()); t[0] = 0
    c = ndi.center_of_mass(sat, lbl, range(1, n + 1))
    return [((x, y), int(t[i + 1])) for i, (y, x) in enumerate(c) if t[i + 1] >= 3]


def lab_de(z, shape, k):
    return np.unpackbits(z[k])[:int(np.prod(shape))].reshape(shape).astype(bool) if k in z.files else None


def situer(dossier, cell):
    """Recale une case (image vidéo entière) sur sa vidéo : essais direct et miroir sur la projection sagittale, puis, si l'erreur
    dépasse 15 px, projections tournées autour de l'axe crânio-caudal (pas 15°). Renvoie les fonctions de conversion et le bilan."""
    man = json.load(open(os.path.join(dossier, 'out', 'manifest.json'), encoding='utf-8')); vox = man['mm_per_voxel']
    F = charger_images(trouver_video(dossier)); n_tour = tour_complet(F)
    g = cv2.cvtColor(cell, cv2.COLOR_BGR2GRAY)
    corr, k, (x0, y0, x1, y1), (a0, b0, a1, b1) = image_la_plus_proche(g, F)
    if corr < 0.8:
        raise SystemExit(f'{dossier} : case non reconnue comme image video entiere (correlation {corr:.2f})')
    sx = (a1 - a0 + 1) / (x1 - x0 + 1); sy = (b1 - b0 + 1) / (y1 - y0 + 1)
    vers_video = lambda P: np.stack([(P[:, 0] - x0) * sx + a0, (P[:, 1] - y0) * sy + b0], 1)
    sil = silhouette(F[k]); W = sil.shape[1]
    z = np.load(os.path.join(dossier, 'work', 'labels.npz')); shape = tuple(int(q) for q in z['shape'])
    env = lab_de(z, shape, 'enveloppe')
    miroir = lambda A: np.hstack([A[:, :2] @ np.array([[-1.0, 0.0], [0.0, 1.0]]), (A[:, 2] + A[:, :2] @ np.array([W - 1.0, 0.0]))[:, None]])
    def essais(E):
        out = []
        r = recalage_2d(sil, E, rot_max_deg=30)
        if r is not None:
            out.append((r[1], r[0], 'direct'))
        r = recalage_2d(sil[:, ::-1], E, rot_max_deg=30)
        if r is not None:
            out.append((r[1], miroir(r[0]), 'miroir'))
        return out
    E0 = env.any(axis=0); cands = [(e, A, m, 0.0, E0) for e, A, m in essais(E0)]
    best = min(cands, key=lambda c: c[0])
    if best[0] > 15:
        dsp = max(1, int(round(0.08 / vox))); env_p = env[::dsp, ::dsp, ::dsp]
        for th in range(15, 180, 15):
            rot = ndi.rotate(env_p.astype(np.uint8), th, axes=(0, 2), reshape=True, order=0) > 0
            E_th = rot.any(axis=0)
            E_th = cv2.resize(E_th.astype(np.uint8), (E_th.shape[1] * dsp, E_th.shape[0] * dsp), interpolation=cv2.INTER_NEAREST) > 0
            for e, A, m in essais(E_th):
                if e < best[0]:
                    best = (e, A, m, float(th), E_th)
    err, A, mode, theta, E = best
    x_med = float(np.median(np.nonzero(env)[0]))
    if theta:
        rot_full = ndi.rotate(env[::2, ::2, ::2].astype(np.uint8), theta, axes=(0, 2), reshape=True, order=0) > 0
        x_med_r = float(np.median(np.nonzero(rot_full)[0])) * 2
        c_rot = np.array(rot_full.shape, float) * 2 / 2; c_env = np.array(shape, float) / 2
    def vers_3d(xy):
        p2 = A[:, :2] @ xy + A[:, 2]
        if theta:
            th = np.radians(-theta); xr, zr = x_med_r - c_rot[0], p2[0] - c_rot[2]
            x_p = c_env[0] + xr * np.cos(th) - zr * np.sin(th); z_p = c_env[2] + xr * np.sin(th) + zr * np.cos(th)
            return np.array([x_p, p2[1], z_p]) * vox
        return np.array([x_med, p2[1], p2[0]]) * vox
    vue = f'{mode}' + (f', projection tournee de {theta:.0f} deg' if theta else '')
    return vers_video, vers_3d, int(k), float(corr), float(err), vue, F[k], man['stage'], round(360.0 * k / n_tour, 1)


def construire(image, grille, stades, marge=6):
    im = cv2.imread(image); H, W = im.shape[:2]; nc, nl = grille
    cw, ch = W // nc, H // nl; resultats = {}
    for idx, cs in enumerate(stades):
        if not cs:
            continue
        i, j = idx // nc, idx % nc
        cell = im[i * ch + marge:(i + 1) * ch - marge, j * cw + marge:(j + 1) * cw - marge]
        pts = blobs(cell)
        if len(pts) < 2:
            print(f'{cs} : {len(pts)} marque(s) seulement, case ignoree'); continue
        pts.sort(key=lambda p: -p[1]); (xo, to), (xq, tq) = pts[0], pts[-1]
        dossier = DOSSIERS[cs]
        vers_video, vers_3d, k, corr, err, vue, frame, stage, angle = situer(dossier, cell)
        Pv = vers_video(np.array([xo, xq], float))
        p_ot = vers_3d(Pv[0]); p_q = vers_3d(Pv[1])
        out = os.path.join(dossier, 'out', 'topographie', 'video360'); os.makedirs(out, exist_ok=True)
        rep = {'stage': stage, 'source': os.path.basename(image), 'case': [i, j], 'image_video': k, 'angle_deg': angle, 'correlation': round(corr, 3),
               'recalage_erreur_px': round(err, 1), 'vue': vue,
               'vesicule_optique': {'xy_video_px': [round(float(v), 1) for v in Pv[0]], 'xyz_mm_pipeline': [round(float(v), 3) for v in p_ot], 'taille_px': to,
                                   'note': 'cercle vert de l utilisateur = vesicule optique (oeil) ; position dans le plan median (l oeil est lateral : seule l abscisse le long de l axe compte) ; repere cranien, pas l ancre otique'},
               'pointe_queue': {'xy_video_px': [round(float(v), 1) for v in Pv[1]], 'xyz_mm_pipeline': [round(float(v), 3) for v in p_q], 'taille_px': tq,
                                'note': 'point vert de l utilisateur = pointe de la queue, plan median'},
               'confiance': 'utilisateur'}
        fj = os.path.join(out, 'reperes_utilisateur.json'); json.dump(rep, open(fj + '.tmp', 'w', encoding='utf-8'), indent=1); os.replace(fj + '.tmp', fj)
        cv2.imwrite(os.path.join(out, 'reperes_utilisateur_image.png'), frame)
        print(f'{cs} : image video {k} ({angle} deg, corr {corr:.2f}, recalage {err:.0f} px, {vue}) ; VO {p_ot.round(2)} mm, queue {p_q.round(2)} mm -> {fj}')
        resultats[cs] = rep
    return resultats


def trianguler(image, grille, stades, angle_video, marge=6, marge_deg=30.0):
    """Vue supplementaire de l'oeil (montage a un autre angle du viewer 360, un point vert par case). La case est recalee sur la
    projection de l'enveloppe tournee de theta autour de l'axe cranio-caudal (theta cherche autour de +-angle annonce et de
    +-angle de l'image video, marge +-marge_deg, pas 5 deg, direct et miroir). Chaque vue donne une equation lineaire
    -sin(-th)(x-c0) + cos(-th)(z-c2) = z'_obs dans le plan transversal ; la vue de profil (theta 0) donne z. Les vues sont
    accumulees dans reperes_utilisateur.json (vesicule_optique.vues) puis resolues aux moindres carres par resoudre()."""
    im = cv2.imread(image); H, W = im.shape[:2]; nc, nl = grille; cw, ch = W // nc, H // nl
    for idx, cs in enumerate(stades):
        if not cs:
            continue
        i, j = idx // nc, idx % nc; cell = im[i * ch + marge:(i + 1) * ch - marge, j * cw + marge:(j + 1) * cw - marge]
        pts = blobs(cell)
        if not pts:
            print(f'{cs} : pas de marque, case ignoree'); continue
        pts.sort(key=lambda q: -q[1]); (xo, to) = pts[0]
        dossier = DOSSIERS[cs]; fr = os.path.join(dossier, 'out', 'topographie', 'video360', 'reperes_utilisateur.json')
        if not os.path.exists(fr):
            print(f'{cs} : pas de vue de profil (reperes_utilisateur.json)'); continue
        rep = json.load(open(fr, encoding='utf-8'))
        man = json.load(open(os.path.join(dossier, 'out', 'manifest.json'), encoding='utf-8')); vox = man['mm_per_voxel']
        F = charger_images(trouver_video(dossier)); n_tour = tour_complet(F)
        g = cv2.cvtColor(cell, cv2.COLOR_BGR2GRAY); corr, k, (x0, y0, x1, y1), (a0, b0, a1, b1) = image_la_plus_proche(g, F)
        if corr < 0.8:
            print(f'{cs} : case non reconnue (correlation {corr:.2f})'); continue
        sx = (a1 - a0 + 1) / (x1 - x0 + 1); sy = (b1 - b0 + 1) / (y1 - y0 + 1); xy = np.array([(xo[0] - x0) * sx + a0, (xo[1] - y0) * sy + b0])
        sil = silhouette(F[k]); Wv = sil.shape[1]
        z = np.load(os.path.join(dossier, 'work', 'labels.npz')); shape = tuple(int(q) for q in z['shape']); env = lab_de(z, shape, 'enveloppe')
        miroir = lambda A: np.hstack([A[:, :2] @ np.array([[-1.0, 0.0], [0.0, 1.0]]), (A[:, 2] + A[:, :2] @ np.array([Wv - 1.0, 0.0]))[:, None]])
        dsp = max(1, int(round(0.08 / vox))); env_p = env[::dsp, ::dsp, ::dsp]; best = None
        angle_k = 360.0 * k / n_tour
        # convention : l'angle du viewer et theta tournent dans le meme sens (verifie sur la vue 244 : theta retenu 235-260 sur les 7 stades) ;
        # on ne cherche donc qu'autour de +angle (annonce et image video) : autour de -angle la silhouette miroir est ambigue (CS14/19/20 a 302)
        cands = sorted({int(round(a / 5) * 5) % 360 for base in (angle_k, angle_video) for a in np.arange(base - marge_deg, base + marge_deg + 1, 5)})
        for th in cands:
            rot = ndi.rotate(env_p.astype(np.uint8), th, axes=(0, 2), reshape=True, order=0) > 0; E_th = rot.any(axis=0)
            E_th = cv2.resize(E_th.astype(np.uint8), (E_th.shape[1] * dsp, E_th.shape[0] * dsp), interpolation=cv2.INTER_NEAREST) > 0
            for s_, mode in ((sil, 'direct'), (sil[:, ::-1], 'miroir')):
                r = recalage_2d(s_, E_th, rot_max_deg=30)
                if r is None:
                    continue
                A = miroir(r[0]) if mode == 'miroir' else r[0]
                if best is None or r[1] < best[0]:
                    best = (r[1], A, mode, float(th))
        if best is None:
            print(f'{cs} : recalage impossible'); continue
        err, A, mode, theta = best
        p2 = A[:, :2] @ xy + A[:, 2]
        rot_full = ndi.rotate(env[::2, ::2, ::2].astype(np.uint8), theta, axes=(0, 2), reshape=True, order=0) > 0
        c_rot = np.array(rot_full.shape, float); c_env = np.array(shape, float) / 2
        vue = {'source': os.path.basename(image), 'angle_annonce_deg': angle_video, 'image_video': int(k), 'angle_video_deg': round(angle_k, 1), 'theta_projection_deg': theta,
               'mode': mode, 'recalage_erreur_px': round(err, 1), 'xy_video_px': [round(float(v), 1) for v in xy], 'xy_marque_case_px': [round(float(v), 1) for v in xo],
               'zr_obs_vox': round(float(p2[0] - c_rot[2]), 2), 'ax1_vox': round(float(p2[1]), 2), 'c_env_vox': [float(c_env[0]), float(c_env[2])],
               'utilisable': bool(err <= 20 and abs(np.sin(np.radians(theta))) >= 0.2)}
        vues = [v for v in rep['vesicule_optique'].get('vues', []) if v['source'] != vue['source']] + [vue]
        rep['vesicule_optique']['vues'] = vues
        json.dump(rep, open(fr + '.tmp', 'w', encoding='utf-8'), indent=1); os.replace(fr + '.tmp', fr)
        print(f"{cs} : {os.path.basename(image)} -> image video {k} ({angle_k:.0f} deg), theta retenu {theta:.0f} deg ({mode}), recalage {err:.0f} px, annonce {angle_video} deg{'' if vue['utilisable'] else ' : NON utilisable'}")


def resoudre(stades=None):
    """Moindres carres sur toutes les vues de l'oeil (profil = equation z, vues tournees = equations lineaires) -> xyz 3D, residu,
    ecart 3D au label yeux ; ecrit vesicule_optique.xyz_mm_3d / residu_mm / ecart_3d_label_yeux_mm."""
    for cs, dossier in DOSSIERS.items():
        if stades and cs not in stades:
            continue
        fr = os.path.join(dossier, 'out', 'topographie', 'video360', 'reperes_utilisateur.json')
        if not os.path.exists(fr):
            continue
        rep = json.load(open(fr, encoding='utf-8')); vo = rep['vesicule_optique']; vues = [v for v in vo.get('vues', []) if v.get('utilisable')]
        if not vues:
            print(f'{cs} : aucune vue supplementaire utilisable'); continue
        man = json.load(open(os.path.join(dossier, 'out', 'manifest.json'), encoding='utf-8')); vox = man['mm_per_voxel']
        z = np.load(os.path.join(dossier, 'work', 'labels.npz')); shape = tuple(int(q) for q in z['shape'])
        c0, c2 = vues[0]['c_env_vox']; p1 = np.array(vo['xyz_mm_pipeline']) / vox
        rows = [[0.0, 1.0]]; rhs = [p1[2] - c2]; ax1s = [p1[1]]                                # profil : z - c2 = z1 - c2
        for v in vues:
            thr = np.radians(-v['theta_projection_deg']); rows.append([-np.sin(thr), np.cos(thr)]); rhs.append(v['zr_obs_vox']); ax1s.append(v['ax1_vox'])
        Am = np.array(rows); b = np.array(rhs); sol, res, rk, sv = np.linalg.lstsq(Am, b, rcond=None)
        resid = float(np.sqrt(np.mean((Am @ sol - b) ** 2))) * vox
        P3 = np.array([c0 + sol[0], float(np.mean(ax1s)), c2 + sol[1]]) * vox
        yeux = lab_de(z, shape, 'yeux'); env = lab_de(z, shape, 'enveloppe'); ec3 = ec3c = None
        if yeux is not None and yeux.any():
            V = np.stack(np.nonzero(yeux[::2, ::2, ::2]), 1) * 2 * vox; ec3 = float(np.min(np.linalg.norm(V - P3, axis=1)))
            lbl, n = ndi.label(yeux); tt = np.bincount(lbl.ravel()); tt[0] = 0
            cents = [np.array(ndi.center_of_mass(lbl == k_)) * vox for k_ in np.argsort(tt)[::-1][:2] if tt[k_] > 20]
            ec3c = float(min(np.linalg.norm(c_ - P3) for c_ in cents)) if cents else None
        dedans = bool(env[tuple(np.clip(np.round(P3 / vox).astype(int), 0, np.array(shape) - 1))])
        vo['xyz_mm_3d'] = [round(float(q), 3) for q in P3]; vo['n_vues'] = len(vues) + 1; vo['residu_mm'] = round(resid, 3)
        vo['ecart_ax1_entre_vues_mm'] = round(float(np.ptp(ax1s)) * vox, 3); vo['dans_enveloppe'] = dedans
        vo['ecart_3d_label_yeux_mm'] = None if ec3 is None else round(ec3, 3); vo['ecart_3d_centroide_yeux_mm'] = None if ec3c is None else round(ec3c, 3)
        vo['angles_retenus'] = [{'source': v['source'], 'annonce': v['angle_annonce_deg'], 'video': v['angle_video_deg'], 'theta': v['theta_projection_deg'], 'mode': v['mode'], 'err_px': v['recalage_erreur_px']} for v in vues]
        json.dump(rep, open(fr + '.tmp', 'w', encoding='utf-8'), indent=1); os.replace(fr + '.tmp', fr)
        print(f"{cs} : {len(vues) + 1} vues, oeil 3D {P3.round(2)} mm, residu {resid:.3f} mm, ecart ax1 {np.ptp(ax1s) * vox:.2f} mm, ecart 3D label yeux {'-' if ec3 is None else f'{ec3:.2f}'} mm (centroide {'-' if ec3c is None else f'{ec3c:.2f}'}), dans enveloppe {dedans} ; angles : " + ', '.join(f"{v['angle_annonce_deg']}->theta {v['theta_projection_deg']:.0f} (video {v['angle_video_deg']:.0f}, {v['recalage_erreur_px']:.0f} px)" for v in vues))


def verifier(stades=None):
    """Après topographie_axe.py : écarts marques / labels, image de contrôle, ligne de tableau dans le manifest. Retourne le tableau."""
    tableau = []
    for cs, dossier in DOSSIERS.items():
        if stades and cs not in stades:
            continue
        out = os.path.join(dossier, 'out', 'topographie'); fr = os.path.join(out, 'video360', 'reperes_utilisateur.json')
        if not os.path.exists(fr):
            continue
        rep = json.load(open(fr, encoding='utf-8')); man = json.load(open(os.path.join(dossier, 'out', 'manifest.json'), encoding='utf-8')); vox = man['mm_per_voxel']
        z = np.load(os.path.join(dossier, 'work', 'labels.npz')); shape = tuple(int(q) for q in z['shape'])
        env = lab_de(z, shape, 'enveloppe'); otiq = lab_de(z, shape, 'vesicules_otiques'); snc = lab_de(z, shape, 'snc')
        if 'vesicule_otique' in rep and 'vesicule_optique' not in rep:      # anciens fichiers : le cercle etait l oeil (erratum utilisateur 17:20)
            rep['vesicule_optique'] = rep.pop('vesicule_otique'); rep['vesicule_optique']['note'] = 'cercle vert = vesicule optique (oeil), repere cranien'
        po = np.array(rep['vesicule_optique']['xyz_mm_pipeline']); pq = np.array(rep['pointe_queue']['xyz_mm_pipeline'])
        yeux = lab_de(z, shape, 'yeux')
        # --- vesicule optique : ecart dans le plan sagittal (ax1, ax2) au label yeux (0 si la marque tombe dans sa projection) et au centroide
        ec_ot = ec_ot_c = None
        if yeux is not None and yeux.any():
            proj = yeux.any(axis=0); ys, zs = np.nonzero(proj); P2 = np.stack([ys, zs], 1) * vox
            ec_ot = float(np.min(np.linalg.norm(P2 - po[1:], axis=1)))
            lbl, n = ndi.label(yeux); t = np.bincount(lbl.ravel()); t[0] = 0
            cents = [np.array(ndi.center_of_mass(lbl == k_)) * vox for k_ in np.argsort(t)[::-1][:2] if t[k_] > 20]
            ec_ot_c = float(min(np.linalg.norm(c_[1:] - po[1:]) for c_ in cents)) if cents else None
        # --- pointe de la queue : bout geodesique de l'enveloppe et extremite du tube neural la plus proche de la marque
        dse = max(1, int(round(0.05 / vox))); env_d = env[::dse, ::dse, ::dse]; ie = np.nonzero(env_d); top = int(np.argmin(ie[1]))
        tip_d, _ = bout_geodesique(env_d, (ie[0][top], ie[1][top], ie[2][top]), pont=2); tip = np.array(tip_d, float) * dse * vox
        ec_q_env = float(np.linalg.norm(tip[1:] - pq[1:]))                      # plan sagittal (la marque est dans le plan median)
        # distance de la marque a l'enveloppe elle-meme (la pointe doit etre dans le corps, a son bord)
        ie2 = np.stack(np.nonzero(env_d), 1) * dse * vox; ec_q_bord = float(np.min(np.linalg.norm(ie2[:, 1:] - pq[1:], axis=1)))
        ec_q_snc = None
        if snc is not None and snc.any():
            s2 = np.stack(np.nonzero(snc[::2, ::2, ::2]), 1) * 2 * vox; ec_q_snc = float(np.min(np.linalg.norm(s2[:, 1:] - pq[1:], axis=1)))
        # --- niveaux imposes et axe prolonge
        niv = json.load(open(os.path.join(out, 'niveaux_imposes.json'), encoding='utf-8')) if os.path.exists(os.path.join(out, 'niveaux_imposes.json')) else None
        axp = json.load(open(os.path.join(out, 'axe_prolonge.json'), encoding='utf-8'))['axe'] if os.path.exists(os.path.join(out, 'axe_prolonge.json')) else None
        s_ot = niv['reperes'].get('vesicule_otique_bord_caudal_s_mm') if niv else None; s_tip = niv['reperes'].get('pointe_queue_s_mm') if niv else None
        src_ot = niv['reperes'].get('source_vesicule_otique') if niv else None; src_tip = niv['reperes'].get('source_pointe_queue') if niv else None
        occ1 = next((n_ for n_ in niv['niveaux'] if n_['nom'] == 'occipital 1'), None) if niv else None
        co4 = next((n_ for n_ in niv['niveaux'] if n_['nom'] == 'Co-4'), None) if niv else None
        # indices : autres labels dont la projection sagittale contient la marque (oeil / cristallin / membres / cordon) : ambiguite de la vue unique
        def prox(k, pt):
            m = lab_de(z, shape, k)
            if m is None or not m.any():
                return None
            ys, zs = np.nonzero(m.any(axis=0)); return round(float(np.min(np.linalg.norm(np.stack([ys, zs], 1) * vox - pt[1:], axis=1))), 3)
        indices_ot = {k: prox(k, po) for k in ('vesicules_otiques', 'cristallins', 'snc') if prox(k, po) is not None}
        indices_q = {k: prox(k, pq) for k in ('membre_inf_gauche', 'membre_inf_droit', 'membres', 'cordon_ombilical') if prox(k, pq) is not None}
        ot_valide = ec_ot is not None and ec_ot <= SEUIL_MM
        ec3 = rep['vesicule_optique'].get('ecart_3d_label_yeux_mm'); res3 = rep['vesicule_optique'].get('residu_mm')
        if ec3 is not None and res3 is not None and res3 <= SEUIL_MM:   # triangulation coherente (residu <= seuil) : l'ecart 3D fait foi
            ot_valide = ec3 <= SEUIL_MM
        # sinon (residu > seuil : videos en perspective, vues incoherentes) l'ecart sagittal reste le critere, la 3D est informative
        q_valide = min(ec_q_env, ec_q_snc if ec_q_snc is not None else 9) <= SEUIL_MM
        # validation manuelle (champ 'validation' pose a la main apres controle visuel) : conservee
        if rep['pointe_queue'].get('validation'):
            q_valide = bool(rep['pointe_queue'].get('valide', True))
        if rep['vesicule_optique'].get('validation'):
            ot_valide = bool(rep['vesicule_optique'].get('valide', True))
        rep['vesicule_optique']['valide'] = bool(ot_valide); rep['pointe_queue']['valide'] = bool(q_valide)
        rep['vesicule_optique']['ecart_label_yeux_mm'] = None if ec_ot is None else round(ec_ot, 3); rep['pointe_queue']['ecart_tube_neural_mm'] = None if ec_q_snc is None else round(ec_q_snc, 3)
        rep['vesicule_optique']['indices_autres_labels_mm'] = indices_ot; rep['pointe_queue']['indices_autres_labels_mm'] = indices_q
        json.dump(rep, open(fr + '.tmp', 'w', encoding='utf-8'), indent=1); os.replace(fr + '.tmp', fr)
        verdict = 'ok'; motifs = []
        if ec_ot is not None and ec_ot > SEUIL_MM:
            verdict = 'a verifier'; motifs.append(f'marque vesicule optique (oeil) a {ec_ot:.2f} mm du label yeux' + (f' (3D : {ec3:.2f} mm, residu {res3:.2f})' if ec3 is not None else ''))
        if ec3 is not None and res3 is not None and res3 > SEUIL_MM:
            motifs.append(f'triangulation oeil incoherente entre vues (residu {res3:.2f} mm > {SEUIL_MM}) : videos en perspective, 3D informative seulement')
        if ec_ot is None:
            motifs.append('pas de label yeux (pas de controle possible)')
        if not q_valide:
            verdict = 'a verifier'; motifs.append(f"marque queue a {ec_q_env:.2f} mm du bout geodesique de l enveloppe et a {'-' if ec_q_snc is None else f'{ec_q_snc:.2f}'} mm du tube neural")
        elif rep['pointe_queue'].get('validation'):
            motifs.append('queue : ' + rep['pointe_queue']['validation'])
        if not ot_valide and ec_ot is None:
            verdict = 'a verifier'
        s_oeil = niv['reperes'].get('oeil_utilisateur_s_mm') if niv else None
        if rep['recalage_erreur_px'] > 20:
            verdict = 'a verifier'; motifs.append(f"recalage {rep['recalage_erreur_px']:.0f} px")
        if src_ot and 'utilisateur' in src_ot:
            motifs.append('ATTENTION : ancre otique prise sur la marque oeil (ancienne interpretation) : topographie a relancer')
        ligne = {'stade': cs, 'image_video': rep['image_video'], 'angle_deg': rep.get('angle_deg'), 'vue': rep['vue'], 'recalage_erreur_px': rep['recalage_erreur_px'],
                 'ecart_oeil_label_yeux_mm': None if ec_ot is None else round(ec_ot, 3), 'ecart_oeil_centroide_mm': None if ec_ot_c is None else round(ec_ot_c, 3), 's_oeil_mm': s_oeil,
                 'ecart_queue_bout_enveloppe_mm': round(ec_q_env, 3), 'ecart_queue_bord_enveloppe_mm': round(ec_q_bord, 3), 'ecart_queue_tube_neural_mm': None if ec_q_snc is None else round(ec_q_snc, 3),
                 's_otique_mm': s_ot, 's_pointe_mm': s_tip, 'source_otique': src_ot, 'source_pointe': src_tip, 'verdict': verdict, 'motifs': motifs,
                 'oeil_xyz_mm_3d': rep['vesicule_optique'].get('xyz_mm_3d'), 'ecart_oeil_3d_label_yeux_mm': rep['vesicule_optique'].get('ecart_3d_label_yeux_mm'),
                 'ecart_oeil_3d_centroide_yeux_mm': rep['vesicule_optique'].get('ecart_3d_centroide_yeux_mm'), 'oeil_n_vues': rep['vesicule_optique'].get('n_vues'),
                 'oeil_residu_mm': rep['vesicule_optique'].get('residu_mm'), 'oeil_angles_retenus': rep['vesicule_optique'].get('angles_retenus'), 'oeil_dans_enveloppe': rep['vesicule_optique'].get('dans_enveloppe'),
                 'oeil_valide': bool(ot_valide), 'queue_valide': bool(q_valide), 'indices_oeil_autres_labels_mm': indices_ot, 'indices_queue_autres_labels_mm': indices_q,
                 'note': 'ecarts mesures dans le plan sagittal (vue unique : la profondeur laterale est inconnue, une marque peut coincider en projection avec un membre ou l oeil)',
                 'seuil_mm': SEUIL_MM, 'images': ['video360/controle_reperes_utilisateur.png']}
        # --- image de controle : video (marques) + projection sagittale (enveloppe, labels, axe prolonge, occipital 1 / Co-4, bout geodesique)
        fimg = os.path.join(out, 'video360', 'reperes_utilisateur_image.png')
        vid = cv2.imread(fimg, cv2.IMREAD_GRAYSCALE) if os.path.exists(fimg) else None
        E0 = env.any(axis=0); ctrl = cv2.cvtColor(E0.astype(np.uint8) * 70, cv2.COLOR_GRAY2BGR)
        if snc is not None:
            ctrl[snc.any(axis=0)] = (110, 110, 110)
        if otiq is not None:
            ctrl[otiq.any(axis=0)] = (255, 160, 60)
        if yeux is not None:
            ctrl[yeux.any(axis=0)] = (200, 60, 200)
        if axp:
            pts = [(int(a['xyz_mm'][2] / vox), int(a['xyz_mm'][1] / vox)) for a in axp]
            for p1, p2 in zip(pts[:-1], pts[1:]):
                cv2.line(ctrl, p1, p2, (0, 0, 255), 2)
        cv2.circle(ctrl, (int(po[2] / vox), int(po[1] / vox)), 7, (0, 255, 0), 2); cv2.circle(ctrl, (int(pq[2] / vox), int(pq[1] / vox)), 6, (0, 255, 0), -1)
        cv2.drawMarker(ctrl, (int(tip[2] / vox), int(tip[1] / vox)), (255, 255, 0), cv2.MARKER_CROSS, 16, 2)
        for lev, col in ((occ1, (255, 0, 255)), (co4, (255, 0, 255))):
            if lev:
                c_ = lev['xyz_centre_mm']; cv2.circle(ctrl, (int(c_[2] / vox), int(c_[1] / vox)), 5, col, 2)
                cv2.putText(ctrl, lev['nom'], (int(c_[2] / vox) + 8, int(c_[1] / vox) + 4), cv2.FONT_HERSHEY_SIMPLEX, 0.5, col, 1)
        if vid is not None:
            v = cv2.cvtColor(vid, cv2.COLOR_GRAY2BGR)
            for key, col in (('vesicule_optique', (0, 255, 0)), ('pointe_queue', (0, 200, 255))):
                xy = rep[key]['xy_video_px']; cv2.circle(v, (int(xy[0]), int(xy[1])), 7, col, 2)
            sc = ctrl.shape[0] / v.shape[0]; v = cv2.resize(v, (int(v.shape[1] * sc), ctrl.shape[0])); ctrl = np.hstack([v, ctrl])
        cv2.putText(ctrl, f"{cs} reperes utilisateur : image {rep['image_video']} ({rep['vue']}, recalage {rep['recalage_erreur_px']:.0f} px) ; vert = marques, orange = vesicules otiques, violet = yeux, rouge = axe prolonge, magenta = occipital 1 / Co-4, croix = bout geodesique ; "
                    f"ecarts : oeil {'-' if ec_ot is None else f'{ec_ot:.2f}'} mm, queue {ec_q_env:.2f} mm (bord {ec_q_bord:.2f}) -> {verdict}", (6, 16), cv2.FONT_HERSHEY_SIMPLEX, 0.42, (0, 255, 255), 1)
        cv2.imwrite(os.path.join(out, 'video360', 'controle_reperes_utilisateur.png'), ctrl)
        fm = os.path.join(out, 'manifest_topographie.json')
        if os.path.exists(fm):
            mj = json.load(open(fm, encoding='utf-8')); mj['controle']['reperes_utilisateur'] = ligne
            mj['controle']['images'] = list(dict.fromkeys(mj['controle'].get('images', []) + ['video360/controle_reperes_utilisateur.png']))
            tmp = fm + '.tmp'; json.dump(mj, open(tmp, 'w', encoding='utf-8'), ensure_ascii=False, indent=1); os.replace(tmp, fm)
        print(f"{cs} : recalage {ligne['recalage_erreur_px']} px ; oeil {ligne['ecart_oeil_label_yeux_mm']} mm (centroide {ligne['ecart_oeil_centroide_mm']}) ; queue {ligne['ecart_queue_bout_enveloppe_mm']} mm (bord {ligne['ecart_queue_bord_enveloppe_mm']}, tube {ligne['ecart_queue_tube_neural_mm']}) ; s_ot {s_ot} s_tip {s_tip} -> {verdict} {motifs}")
        tableau.append(ligne)
    json.dump(tableau, open(os.path.join('embryons_3D', 'reperes_utilisateur_tableau.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    return tableau


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('image', nargs='?'); ap.add_argument('--grille', default='4x2'); ap.add_argument('--stades', default=None, help='liste par case, ligne par ligne, vide = case sans stade')
    ap.add_argument('--verifier', action='store_true')
    ap.add_argument('--vue2', action='store_true', help='vue supplementaire de l oeil (un point par case) : ajoutee aux vues a trianguler')
    ap.add_argument('--angle', type=float, default=302.0, help='angle annonce (viewer) de cette vue')
    ap.add_argument('--resoudre', action='store_true', help='moindres carres sur toutes les vues de l oeil')
    a = ap.parse_args()
    if a.vue2:
        nc, nl = [int(v) for v in a.grille.lower().split('x')]
        trianguler(a.image, (nc, nl), [s.strip() for s in a.stades.split(',')], a.angle)
    elif a.resoudre:
        resoudre([s.strip() for s in a.stades.split(',')] if a.stades else None)
    elif a.verifier:
        verifier([s.strip() for s in a.stades.split(',')] if a.stades else None)
    else:
        nc, nl = [int(v) for v in a.grille.lower().split('x')]
        construire(a.image, (nc, nl), [s.strip() for s in a.stades.split(',')])
