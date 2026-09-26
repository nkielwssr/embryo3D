"""
Etages vertebraux marques a la main par l'utilisateur sur une image de la video 360 (points de couleur saturee,
un par etage, sur une vue de profil). L'image peut etre redimensionnee ou recadree par rapport a la video.
  1. detection des marques (pixels satures : une composante couleur domine les deux autres de > 60)
  2. identification de l'image de la video la plus proche (silhouettes recalees par leurs boites englobantes)
  3. recalage 2D de cette image sur la projection sagittale de l'enveloppe du pipeline (video360.recalage_2d)
  4. marques -> plan median du pipeline (mm) ; ordre = du cou vers la queue (la marque la plus proche du sommet
     de la tete projete est la premiere) ; frontieres = milieux entre marques successives
Sorties : <dossier>/out/topographie/video360/etages_manuels.json, <CS>_etages_manuels.ply (billes), controle_marques.png ;
bloc controle.etages_manuels + segment etages_manuels dans manifest_topographie.json.
Usage : python embryo3d/marques_video.py CS15_f4v "CS15_f4v/point verts.png"
"""
import sys, os, json, argparse
import numpy as np, cv2
from scipy import ndimage as ndi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from video360 import trouver_video, charger_images, silhouette, tour_complet, recalage_2d
from topographie_axe import bout_geodesique


def marques(im):
    b, g, r = [c.astype(int) for c in cv2.split(im)]
    sat = ((b - np.maximum(r, g)) > 60) | ((g - np.maximum(r, b)) > 60) | ((r - np.maximum(g, b)) > 60)
    lbl, n = ndi.label(sat); t = np.bincount(lbl.ravel()); t[0] = 0
    c = ndi.center_of_mass(sat, lbl, range(1, n + 1))
    return np.array([(x, y) for i, (y, x) in enumerate(c) if t[i + 1] >= 3], float)


def boite(sil):
    ys, xs = np.nonzero(sil); return xs.min(), ys.min(), xs.max(), ys.max()


def image_la_plus_proche(im_gray, F):
    """Silhouette de l'image annotee ramenee a la boite englobante de chaque image video ; correlation."""
    s0 = silhouette(im_gray); x0, y0, x1, y1 = boite(s0)
    crop = im_gray[y0:y1 + 1, x0:x1 + 1].astype(float)
    best = None
    for k, f in enumerate(F):
        s = silhouette(f); a0, b0, a1, b1 = boite(s)
        w, h = a1 - a0 + 1, b1 - b0 + 1
        c2 = cv2.resize(crop, (w, h)); ref = f[b0:b1 + 1, a0:a1 + 1].astype(float)
        c2 = c2 - c2.mean(); ref = ref - ref.mean()
        corr = float((c2 * ref).sum() / (np.sqrt((c2 ** 2).sum() * (ref ** 2).sum()) + 1e-9))
        if best is None or corr > best[0]:
            best = (corr, k, (x0, y0, x1, y1), (a0, b0, a1, b1))
    return best


def situer_zoom(crop_gray, frame_gray, marques_xy, echelles=np.arange(0.25, 1.01, 0.025)):
    """Un zoom (recadrage agrandi) d'une image video : recherche de motif multi-echelle (marques masquees).
    Retourne (marques en coordonnees de l'image video, echelle, score)."""
    best = None
    masque = np.ones(crop_gray.shape, np.uint8)
    for x, y in marques_xy:
        cv2.circle(masque, (int(x), int(y)), 6, 0, -1)
    for e in echelles:
        tpl = cv2.resize(crop_gray, None, fx=e, fy=e, interpolation=cv2.INTER_AREA); mk = cv2.resize(masque, (tpl.shape[1], tpl.shape[0]), interpolation=cv2.INTER_NEAREST)
        if tpl.shape[0] >= frame_gray.shape[0] or tpl.shape[1] >= frame_gray.shape[1] or tpl.shape[0] < 20:
            continue
        # marques effacees par inpainting plutot que masque (le masque rend TM_CCOEFF_NORMED instable / NaN)
        tpl_i = cv2.inpaint(tpl, (mk == 0).astype(np.uint8), 3, cv2.INPAINT_TELEA)
        res = cv2.matchTemplate(frame_gray, tpl_i, cv2.TM_CCOEFF_NORMED)
        _, mx, _, loc = cv2.minMaxLoc(res)
        if not np.isfinite(mx):
            continue
        if best is None or mx > best[0]:
            best = (float(mx), float(e), loc)
    sc, e, (lx, ly) = best
    Pv = np.stack([marques_xy[:, 0] * e + lx, marques_xy[:, 1] * e + ly], 1)
    return Pv, e, sc


def construire(dossier, image, nom_serie='etages_manuels', zooms=(), prolonger=0, reference=False):
    out = os.path.join(dossier, 'out', 'topographie', 'video360'); os.makedirs(out, exist_ok=True)
    man = json.load(open(os.path.join(dossier, 'out', 'manifest.json'), encoding='utf-8')); cs = man['stage']; vox = man['mm_per_voxel']
    im = cv2.imread(image); P = marques(im)
    print(f'{cs} : {len(P)} marques dans {os.path.basename(image)} ({im.shape[1]}x{im.shape[0]})')
    v = trouver_video(dossier); F = charger_images(v); n_tour = tour_complet(F)
    im_gray = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY)
    corr, k, (x0, y0, x1, y1), (a0, b0, a1, b1) = image_la_plus_proche(im_gray, F)
    print(f'  image video la plus proche : {k} ({360.0 * k / n_tour:.0f} deg), correlation {corr:.2f}')
    if corr >= 0.8:
        # image entiere (redimensionnee) : marques -> coordonnees video par les boites englobantes
        sx = (a1 - a0 + 1) / (x1 - x0 + 1); sy = (b1 - b0 + 1) / (y1 - y0 + 1)
        Pv = np.stack([(P[:, 0] - x0) * sx + a0, (P[:, 1] - y0) * sy + b0], 1)
    else:
        # l'image annotee est un recadrage agrandi : recherche de motif multi-echelle dans les images video (1 sur 4)
        best = None
        for kk in range(0, n_tour, 4):
            Pk, e, sc = situer_zoom(im_gray, F[kk], P, echelles=np.arange(0.3, 1.0, 0.05))
            if best is None or sc > best[0]:
                best = (sc, kk, Pk, e)
        sc, k, Pv, e = best
        print(f'  image annotee = recadrage : situee dans l image video {k} ({360.0 * k / n_tour:.0f} deg) a l echelle {e:.2f}, score {sc:.2f}')
    # zooms complementaires (recadrages agrandis de la meme vue) : situes par recherche de motif dans l'image k
    for zimg in zooms:
        zi = cv2.imread(zimg); Pz = marques(zi); zg = cv2.cvtColor(zi, cv2.COLOR_BGR2GRAY)
        Pzv, e, sc = situer_zoom(zg, F[k], Pz)
        # fusion : on ecarte les marques du zoom a moins de 6 px d'une marque deja presente
        nouv = [q for q in Pzv if np.linalg.norm(Pv - q, axis=1).min() > 6]
        print(f'  zoom {os.path.basename(zimg)} : {len(Pz)} marques, situe a l echelle {e:.2f} (score {sc:.2f}), {len(nouv)} nouvelles')
        if nouv:
            Pv = np.vstack([Pv, np.array(nouv)])
    sil = silhouette(F[k])
    # marques de la tete (reperes hors chaine : oeil, oreille, arcs) : dans le disque du point le plus epais de la silhouette
    dist_v = ndi.distance_transform_edt(sil); ty, tx = np.unravel_index(int(np.argmax(dist_v)), dist_v.shape); tete_v = np.array([tx, ty], float)
    # les marques de la chaine sont pres de la surface dorsale (profondeur faible) ; les reperes de tete (oeil, oreille,
    # arcs) sont profonds : on ecarte les marques dont la profondeur depasse 0,55 x le rayon de la tete
    r_tete = float(dist_v.max())
    prof_m = np.array([dist_v[int(np.clip(y, 0, sil.shape[0] - 1)), int(np.clip(x, 0, sil.shape[1] - 1))] for x, y in Pv])
    # type de vue : la vue de profil est la plus large ; a ~90 deg on voit le dos ou le ventre -> projection frontale
    ws = [np.ptp(np.nonzero(silhouette(F[kk]).any(axis=0))[0]) for kk in range(0, n_tour, 2)]
    kL = 2 * int(np.argmax(ws)); d_ang = abs(((k - kL) * 360.0 / n_tour + 180) % 360 - 180)
    dorsale = 50 <= d_ang <= 130
    print(f'  vue : {"dos / ventre" if dorsale else "profil"} (a {d_ang:.0f} deg de la vue de profil {kL})')
    # en vue de dos la chaine est sur la ligne mediane (profonde dans la silhouette) : pas de filtre de profondeur
    dans_tete = (prof_m > 0.55 * r_tete) & (not dorsale)
    if dans_tete.any():
        print(f'  {int(dans_tete.sum())} marque(s) profonde(s) ecartee(s) (reperes de tete hors chaine)')
        Pv = Pv[~dans_tete]
    # recalage 2D sur la projection de l'enveloppe : projection sagittale directe, puis, si elle recale mal, recherche de
    # l'angle de projection (rotation autour de l'axe cranio-caudal) qui reproduit la silhouette de la vue video
    z = np.load(os.path.join(dossier, 'work', 'labels.npz')); shape = tuple(int(q) for q in z['shape'])
    env = np.unpackbits(z['enveloppe'])[:int(np.prod(shape))].reshape(shape).astype(bool)
    theta = 0.0; E = env.any(axis=2).T if dorsale else env.any(axis=0)      # frontale (lignes ax1, colonnes ax0) ou sagittale (ax1, ax2)
    # vue du cote oppose au profil de reference (> 130 deg) : l'ICP force une rotation propre et part tete-beche ; on recale la
    # silhouette retournee (miroir horizontal) puis on compose le miroir dans A, qui s'applique ensuite aux coordonnees d'origine
    oppose = (not dorsale) and d_ang > 130
    if oppose:
        W = sil.shape[1]; r = recalage_2d(sil[:, ::-1], E, rot_max_deg=30)
        if r is None:
            r = recalage_2d(sil[:, ::-1], E)
        A, err = r; A = np.hstack([A[:, :2] @ np.array([[-1.0, 0.0], [0.0, 1.0]]), (A[:, 2] + A[:, :2] @ np.array([W - 1.0, 0.0]))[:, None]])
        print(f'  vue du cote oppose au profil de reference : recalage par miroir, {err:.0f} px')
    else:
        r = recalage_2d(sil, E, rot_max_deg=30)
        if r is None:
            r = recalage_2d(sil, E); print('  recalage a rotation libre (aucune similitude a moins de 30 deg)')
        A, err = r
    if err > 15 and not dorsale and not oppose:
        dsp = max(1, int(round(0.08 / vox))); env_p = env[::dsp, ::dsp, ::dsp]
        meilleur = (err, 0.0, E, A)
        for th in range(15, 180, 15):
            rot = ndi.rotate(env_p.astype(np.uint8), th, axes=(0, 2), reshape=True, order=0) > 0
            E_th = rot.any(axis=0)
            E_th = cv2.resize(E_th.astype(np.uint8), (E_th.shape[1] * dsp, E_th.shape[0] * dsp), interpolation=cv2.INTER_NEAREST) > 0
            A_th, e_th = recalage_2d(sil, E_th, rot_max_deg=30)
            if e_th < meilleur[0]:
                meilleur = (e_th, float(th), E_th, A_th)
        err, theta, E, A = meilleur
        print(f'  projection de l enveloppe tournee de {theta:.0f} deg autour de l axe cranio-caudal : recalage {err:.0f} px')
    mm_px = float(np.sqrt(abs(np.linalg.det(A[:, :2]))) * vox)
    # plan median dans le repere tourne : x' = mediane de l'enveloppe tournee ; retour au repere pipeline par la rotation inverse
    if theta:
        rot_full = ndi.rotate(env[::2, ::2, ::2].astype(np.uint8), theta, axes=(0, 2), reshape=True, order=0) > 0
        x_med_r = float(np.median(np.nonzero(rot_full)[0])) * 2
        c_rot = np.array(rot_full.shape, float) * 2 / 2; c_env = np.array(shape, float) / 2
    x_med = float(np.median(np.nonzero(env)[0]))
    Ainv = np.linalg.inv(np.vstack([A, [0, 0, 1]]))
    # sommet de la tete (pipeline) projete : ordre des marques du cou vers la queue
    dse = max(1, int(round(0.06 / vox))); env_d = env[::dse, ::dse, ::dse]; ie = np.nonzero(env_d); top = int(np.argmin(ie[1]))
    start = np.array([ie[0][top], ie[1][top], ie[2][top]], float) * dse
    xt = Ainv[:2, :2] @ np.array([start[2], start[1]]) + Ainv[:2, 2]
    # ordonner les marques le long de la chaine (plus proche voisin) ; depart = extremite la plus proche de la TETE de la
    # silhouette video = point le plus epais (maximum de la distance au bord), independant du recalage
    def chaine_depuis(i0):
        reste = list(range(len(Pv))); ordre = [i0]; reste.remove(i0)
        while reste:
            d = np.linalg.norm(Pv[reste] - Pv[ordre[-1]], axis=1); j = reste[int(np.argmin(d))]; ordre.append(j); reste.remove(j)
        return ordre
    # extremites de la chaine = les deux points les plus eloignes l'un de l'autre le long du parcours
    o0 = chaine_depuis(int(np.argmax(np.linalg.norm(Pv - Pv.mean(0), axis=1)))); ext_a, ext_b = o0[0], o0[-1]
    i_dep = ext_a if np.linalg.norm(Pv[ext_a] - tete_v) < np.linalg.norm(Pv[ext_b] - tete_v) else ext_b
    Pv = Pv[chaine_depuis(i_dep)]
    if len(Pv) > 8:
        sauts = np.linalg.norm(np.diff(Pv, axis=0), axis=1); med = float(np.median(sauts))
        coupe = [i for i in range(8, len(sauts)) if sauts[i] > 3.0 * med]
        if coupe:
            print(f'  chaine coupee au saut {coupe[0] + 1} ({sauts[coupe[0]]:.0f} px > 3 x {med:.0f}) : {len(Pv) - coupe[0] - 1} marque(s) ecartee(s)')
            Pv = Pv[:coupe[0] + 1]
    # etages extrapoles au bout de la chaine (cote queue) : pas median des marques en px, direction des 3 dernieres marques,
    # ramenes dans la silhouette video si necessaire
    n_reel = len(Pv); extrapole = np.zeros(n_reel, bool)
    if prolonger > 0 and n_reel >= 3:
        pas_px = float(np.median(np.linalg.norm(np.diff(Pv, axis=0), axis=1)))
        dirn = Pv[-1] - Pv[-3]; dirn /= np.linalg.norm(dirn) + 1e-9
        idx_sil = np.stack(np.nonzero(sil), 1)[:, ::-1].astype(float)          # (x, y) des pixels de la silhouette
        for m_ in range(1, prolonger + 1):
            q = Pv[-1] + dirn * pas_px
            if not sil[int(np.clip(q[1], 0, sil.shape[0] - 1)), int(np.clip(q[0], 0, sil.shape[1] - 1))]:
                q = idx_sil[int(np.argmin(np.linalg.norm(idx_sil - q, axis=1)))]  # rappel dans la silhouette
            Pv = np.vstack([Pv, q]); extrapole = np.r_[extrapole, True]
        print(f'  {prolonger} etage(s) extrapole(s) au bout de la chaine (pas {pas_px:.1f} px)')
    pts = []
    fa = os.path.join(dossier, 'out', 'topographie', 'axe_vertebral.json')
    axe_auto = json.load(open(fa, encoding='utf-8'))['axe'] if os.path.exists(fa) else None
    for i, xy in enumerate(Pv):
        p2 = A[:, :2] @ xy + A[:, 2]
        if dorsale:
            ax0, ax1 = p2[0], p2[1]
            if axe_auto:
                Qa_ = np.array([a.get('xyz_corps_mm', a['xyz_mm']) for a in axe_auto]) / vox
                j = int(np.argmin(np.abs(Qa_[:, 1] - ax1))); ax2 = Qa_[j, 2]
            else:
                row = env[:, int(np.clip(ax1, 0, shape[1] - 1)), :]; ax2 = float(np.median(np.nonzero(row)[1])) if row.any() else shape[2] / 2
            p_vox = np.array([ax0, ax1, ax2])
        elif theta:
            # point dans le repere tourne (x' median, ax1 = p2[1], z' = p2[0]) -> rotation inverse autour de ax1, centres des volumes
            th = np.radians(-theta); xr, zr = x_med_r - c_rot[0], p2[0] - c_rot[2]
            x_p = c_env[0] + xr * np.cos(th) - zr * np.sin(th); z_p = c_env[2] + xr * np.sin(th) + zr * np.cos(th)
            p_vox = np.array([x_p, p2[1], z_p])
        else:
            p_vox = np.array([x_med, p2[1], p2[0]])
        pts.append(p_vox * vox)
    pts = np.array(pts); d_arc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]
    pas = float(np.median(np.diff(d_arc))) if len(pts) > 1 else 0
    etages = [{'n': i + 1, 'xy_video_px': [round(float(xy[0]), 1), round(float(xy[1]), 1)], 'xyz_mm_pipeline': [round(float(q), 3) for q in pts[i]],
               's_marques_mm': round(float(d_arc[i]), 3), 'extrapole': bool(extrapole[i])} for i, xy in enumerate(Pv)]
    json.dump({'stage': cs, 'image': os.path.basename(image), 'image_video': int(k), 'angle_deg': round(360.0 * k / n_tour, 1), 'correlation': round(corr, 3), 'projection_theta_deg': theta, 'vue': 'dos/ventre' if dorsale else ('profil oppose (miroir)' if oppose else 'profil'),
               'recalage_erreur_px': round(err, 1), 'mm_par_px': round(mm_px, 5), 'n_marques': n_reel, 'n_extrapoles': int(extrapole.sum()), 'pas_median_mm': round(pas, 3),
               'vue_reference_calage': bool(reference),
               'note': 'marques utilisateur = centres d etages, ordonnees du cou vers la queue ; positions dans le plan median (X = mediane de l enveloppe)',
               'etages': etages}, open(os.path.join(out, f'{nom_serie}.json'), 'w', encoding='utf-8'), indent=1)
    # billes PLY (repere des PLY du pipeline)
    import trimesh
    c0 = np.array(man['center_voxel']); cc = np.array([c0[0], c0[2], -c0[1]]) * vox
    r = float(np.clip(0.012 * man.get('greatest_length_mm_assumed', 10), 0.06, 0.2)); ms = []
    for p in pts:
        Pb = np.array([p[0], p[2], -p[1]]) - cc; sp = trimesh.creation.icosphere(subdivisions=1, radius=r); sp.apply_translation(Pb); ms.append(sp)
    m = trimesh.util.concatenate(ms); fn = f'{cs}_{nom_serie}.ply'; m.export(os.path.join(out, fn))
    # controle : marques numerotees sur l'image video + projection sur l'enveloppe
    img = cv2.cvtColor(F[k], cv2.COLOR_GRAY2BGR)
    for i, xy in enumerate(Pv):
        cv2.circle(img, (int(xy[0]), int(xy[1])), 4, (0, 255, 0), 1 if extrapole[i] else -1); cv2.putText(img, str(i + 1), (int(xy[0]) + 5, int(xy[1]) + 4), cv2.FONT_HERSHEY_SIMPLEX, 0.35, (0, 255, 0), 1)
    cv2.drawMarker(img, (int(xt[0]), int(xt[1])), (255, 0, 255), cv2.MARKER_CROSS, 14, 1)
    E0 = env.any(axis=0)                                    # projection sagittale du pipeline (repere des points)
    ctrl = cv2.cvtColor(E0.astype(np.uint8) * 90, cv2.COLOR_GRAY2BGR)
    for p in pts:
        q = p / vox; cv2.circle(ctrl, (int(q[2]), int(q[1])), 3, (0, 255, 0), -1)
    if axe_auto:                                            # etages automatiques (bleu) pour comparaison
        fe = os.path.join(dossier, 'out', 'topographie', 'etages_vertebraux.json')
        if os.path.exists(fe):
            ej = json.load(open(fe, encoding='utf-8'))['etages']; Sa_ = np.array([a['s_mm'] for a in axe_auto]); Qa_ = np.array([a.get('xyz_corps_mm', a['xyz_mm']) for a in axe_auto])
            for nom, e in ej.items():
                j = int(np.argmin(np.abs(Sa_ - e['s_debut_mm']))); q = Qa_[j] / vox; cv2.circle(ctrl, (int(q[2]), int(q[1])), 2, (255, 120, 0), 1)
    # et la projection frontale (le long de l'axe AP) pour verifier la position laterale des points
    E1 = env.any(axis=2).T; ctrl2 = cv2.cvtColor(E1.astype(np.uint8) * 90, cv2.COLOR_GRAY2BGR)   # lignes = ax1, colonnes = ax0
    for p in pts:
        q = p / vox; cv2.circle(ctrl2, (int(q[0]), int(q[1])), 3, (0, 255, 0), -1)
    ctrl = np.hstack([ctrl, ctrl2])
    sc = img.shape[0] / ctrl.shape[0]; ctrl = cv2.resize(ctrl, (int(ctrl.shape[1] * sc), img.shape[0]))
    plan = np.hstack([img, ctrl])
    cv2.putText(plan, f'{cs} : {len(Pv)} marques utilisateur (vert, numerotees du cou a la queue) sur l image video {k} ; a droite : enveloppe pipeline, sagittale puis frontale ; recalage {err:.0f} px (theta {theta:.0f} deg), pas median {pas:.2f} mm',
                (6, 16), cv2.FONT_HERSHEY_SIMPLEX, 0.42, (0, 255, 255), 1)
    cv2.imwrite(os.path.join(out, f'controle_{nom_serie}.png'), plan)
    fm = os.path.join(dossier, 'out', 'topographie', 'manifest_topographie.json')
    if os.path.exists(fm):
        mj = json.load(open(fm, encoding='utf-8'))
        mj['controle'][nom_serie] = {'image': os.path.basename(image), 'image_video': int(k), 'correlation': round(corr, 3), 'recalage_erreur_px': round(err, 1),
                                     'n_marques': n_reel, 'n_extrapoles': int(extrapole.sum()), 'pas_median_mm': round(pas, 3), 'vue': 'dos/ventre' if dorsale else 'profil',
                                     'vue_reference_calage': bool(reference), 'images': [f'video360/controle_{nom_serie}.png']}
        mj['controle']['images'] = list(dict.fromkeys(mj['controle'].get('images', []) + [f'video360/controle_{nom_serie}.png']))
        mj['segments'] = [sg for sg in mj['segments'] if sg['name'] != nom_serie] + [
            {'name': nom_serie, 'file': f'video360/{fn}', 'color': [0.2, 0.9, 0.2], 'alpha': 1.0, 'faces': int(len(ms) * 80), 'volume_mm3': 0.0,
             'confiance': 'utilisateur', 'note': f'{n_reel} etages marques par l utilisateur sur la video 360 (image {k}), plan median'
                     + (f' + {int(extrapole.sum())} extrapole(s) au pas median (bille creuse sur le controle)' if extrapole.any() else '')
                     + (' ; vue de reference pour le reperage / calage' if reference else '')}]
        tmp = fm + '.tmp'; json.dump(mj, open(tmp, 'w', encoding='utf-8'), ensure_ascii=False, indent=1); os.replace(tmp, fm)
    print(f'  {len(Pv)} etages, pas median {pas:.2f} mm, recalage {err:.0f} px ({mm_px * 1000:.1f} um/px) -> {out}')


def fusionner(dossier, principale, complements, nom_serie='etages_manuels'):
    """Fusion 3D de series de marques prises sur des vues differentes : la serie principale (du cou vers la queue, sans ses
    etages extrapoles) est prolongee par chaque complement, oriente dans le sens qui s'eloigne du bout de la chaine ; les points
    du complement deja couverts (avant le point le plus proche du dernier etage) sont ignores."""
    out = os.path.join(dossier, 'out', 'topographie', 'video360')
    man = json.load(open(os.path.join(dossier, 'out', 'manifest.json'), encoding='utf-8')); cs = man['stage']; vox = man['mm_per_voxel']
    lire = lambda s: json.load(open(os.path.join(out, f'{s}.json'), encoding='utf-8'))
    J0 = lire(principale); pts = [np.array(e['xyz_mm_pipeline']) for e in J0['etages'] if not e.get('extrapole')]
    src = [principale] * len(pts); vues = {principale: {'image': J0['image'], 'image_video': J0['image_video'], 'n': len(pts)}}
    pas = J0['pas_median_mm']
    for comp in complements:
        Jc = lire(comp); Pc = np.array([e['xyz_mm_pipeline'] for e in Jc['etages'] if not e.get('extrapole')])
        fin = pts[-1]; d = np.linalg.norm(Pc - fin, axis=1); j = int(np.argmin(d))
        # sens de prolongement : celui ou l'on s'eloigne du bout de la chaine principale
        avant = d[:j].mean() if j > 0 else -1; apres = d[j + 1:].mean() if j < len(Pc) - 1 else -1
        suite = Pc[j + 1:] if apres >= avant else Pc[:j][::-1]
        suite = [q for q in suite if np.linalg.norm(q - fin) > 0.5 * pas]
        pts += list(suite); src += [comp] * len(suite)
        vues[comp] = {'image': Jc['image'], 'image_video': Jc['image_video'], 'n': len(suite), 'raccord_mm': round(float(d[j]), 3)}
        print(f'  {comp} : raccord a {d[j]:.2f} mm du dernier etage de {principale}, {len(suite)} etage(s) ajoute(s) vers la queue')
    pts = np.array(pts); d_arc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]
    pas = float(np.median(np.diff(d_arc)))
    etages = [{'n': i + 1, 'xyz_mm_pipeline': [round(float(q), 3) for q in pts[i]], 's_marques_mm': round(float(d_arc[i]), 3), 'source': src[i]}
              for i in range(len(pts))]
    fj = os.path.join(out, f'{nom_serie}.json')
    json.dump({'stage': cs, 'fusion': vues, 'n_marques': len(pts), 'pas_median_mm': round(pas, 3), 'vue_reference_calage': bool(J0.get('vue_reference_calage')),
               'image': J0['image'], 'image_video': J0['image_video'],
               'note': 'fusion 3D de marques utilisateur prises sur plusieurs vues de la video 360, ordonnees du cou vers la queue', 'etages': etages},
              open(fj + '.tmp', 'w', encoding='utf-8'), indent=1); os.replace(fj + '.tmp', fj)
    import trimesh
    c0 = np.array(man['center_voxel']); cc = np.array([c0[0], c0[2], -c0[1]]) * vox
    r = float(np.clip(0.012 * man.get('greatest_length_mm_assumed', 10), 0.06, 0.2)); ms = []
    for p in pts:
        sp = trimesh.creation.icosphere(subdivisions=1, radius=r); sp.apply_translation(np.array([p[0], p[2], -p[1]]) - cc); ms.append(sp)
    fn = f'{cs}_{nom_serie}.ply'; trimesh.util.concatenate(ms).export(os.path.join(out, fn))
    # controle : projections sagittale et frontale de l'enveloppe, une couleur par serie
    z = np.load(os.path.join(dossier, 'work', 'labels.npz')); shape = tuple(int(q) for q in z['shape'])
    env = np.unpackbits(z['enveloppe'])[:int(np.prod(shape))].reshape(shape).astype(bool)
    couleurs = [(0, 255, 0), (0, 200, 255), (255, 0, 255), (0, 255, 255)]; series = list(vues)
    E0 = env.any(axis=0); ctrl = cv2.cvtColor(E0.astype(np.uint8) * 90, cv2.COLOR_GRAY2BGR)
    E1 = env.any(axis=2).T; ctrl2 = cv2.cvtColor(E1.astype(np.uint8) * 90, cv2.COLOR_GRAY2BGR)
    for i, p in enumerate(pts):
        q = p / vox; col = couleurs[series.index(src[i]) % len(couleurs)]
        cv2.circle(ctrl, (int(q[2]), int(q[1])), 3, col, -1); cv2.circle(ctrl2, (int(q[0]), int(q[1])), 3, col, -1)
        if i % 5 == 0 or i == len(pts) - 1:
            cv2.putText(ctrl, str(i + 1), (int(q[2]) + 5, int(q[1]) + 4), cv2.FONT_HERSHEY_SIMPLEX, 0.4, col, 1)
    plan = np.hstack([ctrl, ctrl2])
    legende = ' ; '.join(f'{s} ({v["n"]}, image {v["image_video"]})' for s, v in vues.items())
    cv2.putText(plan, f'{cs} : {len(pts)} etages fusionnes ' + legende + f' ; pas median {pas:.2f} mm', (6, 16), cv2.FONT_HERSHEY_SIMPLEX, 0.42, (0, 255, 255), 1)
    cv2.imwrite(os.path.join(out, f'controle_{nom_serie}.png'), plan)
    fm = os.path.join(dossier, 'out', 'topographie', 'manifest_topographie.json')
    if os.path.exists(fm):
        mj = json.load(open(fm, encoding='utf-8'))
        mj['controle'][nom_serie] = {'fusion': vues, 'n_marques': len(pts), 'n_extrapoles': 0, 'pas_median_mm': round(pas, 3),
                                     'vue_reference_calage': bool(J0.get('vue_reference_calage')), 'images': [f'video360/controle_{nom_serie}.png']}
        for s in list(complements) + [principale]:
            if s != nom_serie:
                mj['controle'].pop(s, None)
        autres = [s for s in list(complements) + [principale] if s != nom_serie]
        mj['controle']['images'] = [im for im in mj['controle'].get('images', []) if not any(im.endswith(f'controle_{s}.png') for s in autres)]
        mj['controle']['images'] = list(dict.fromkeys(mj['controle']['images'] + [f'video360/controle_{nom_serie}.png']))
        mj['segments'] = [sg for sg in mj['segments'] if sg['name'] not in [nom_serie] + autres] + [
            {'name': nom_serie, 'file': f'video360/{fn}', 'color': [0.2, 0.9, 0.2], 'alpha': 1.0, 'faces': int(len(ms) * 80), 'volume_mm3': 0.0,
             'confiance': 'utilisateur', 'note': f'{len(pts)} etages marques par l utilisateur sur la video 360, fusion de {len(vues)} vues (images '
                                                 + ', '.join(str(v['image_video']) for v in vues.values()) + '), plan median'}]
        tmp = fm + '.tmp'; json.dump(mj, open(tmp, 'w', encoding='utf-8'), ensure_ascii=False, indent=1); os.replace(tmp, fm)
    print(f'  {len(pts)} etages fusionnes, pas median {pas:.2f} mm -> {fj}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('dossier'); ap.add_argument('image', nargs='?'); ap.add_argument('--nom', default='etages_manuels')
    ap.add_argument('--fusion', nargs='+', help='fusion 3D : serie principale puis series complementaires (noms des JSON de video360/)')
    ap.add_argument('--zoom', action='append', default=[], help='image(s) zoomee(s) de la meme vue avec des marques complementaires')
    ap.add_argument('--prolonger', type=int, default=0, help='nombre d etages extrapoles au bout de la chaine (cote queue), au pas median')
    ap.add_argument('--reference', action='store_true', help='marquer cette image comme vue de reference de reperage / calage')
    a = ap.parse_args()
    if a.fusion:
        fusionner(a.dossier, a.fusion[0], a.fusion[1:], a.nom)
    else:
        construire(a.dossier, a.image, a.nom, a.zoom, a.prolonger, a.reference)
