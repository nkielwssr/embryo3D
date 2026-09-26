"""
Extraction depuis la video de rotation 360 deg d'un stade (rendu volumique semi-transparent, fond noir) :
  1. images de la rotation (N ~ 200 = 360 deg, axe de rotation vertical de l'image)
  2. vues de profil (largeur maximale de la silhouette, et vue opposee) : la chaine segmentaire (somites / corps
     vertebraux) se lit le long du contour dorsal ; profil d'intensite dans une bande interne au contour, periode par
     autocorrelation glissante, frontieres = minima du profil detendu
  3. recalage 2D (similitude + miroir) de la silhouette de profil sur la projection sagittale de l'enveloppe du
     pipeline (labels.npz) : donne l'echelle mm/px et le repere ; les frontieres du contour dorsal sont posees dans le
     plan median du pipeline (X = mediane de l'enveloppe, Y/Z par la similitude)
  4. rapprochement avec l'axe vertebral (out/topographie/axe_vertebral.json) : s de chaque frontiere video projetee
     sur l'axe, comparaison avec les etages ; les frontieres au-dela de l'axe (queue) sont conservees telles quelles
Sorties : <dossier>/out/topographie/video360/ : geometrie.json, recalage2d.json, somites_video.json (frontieres, mm,
repere pipeline), <CS>_somites_video.ply (billes aux frontieres, repere des PLY du pipeline), profil_somites.png,
controle_recalage.png ; bloc controle.video360 dans manifest_topographie.json.
Usage : python embryo3d/video360.py CS14_f4v
"""
import sys, os, json, glob, argparse
import numpy as np, cv2
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def trouver_video(dossier):
    for v in sorted(glob.glob(os.path.join(dossier, '*.mp4')) + glob.glob(os.path.join(dossier, '*.f4v'))):
        c = cv2.VideoCapture(v); ok, f = c.read(); c.release()
        if ok and f.mean() < 60:
            return v
    return None


def charger_images(v):
    c = cv2.VideoCapture(v); F = []
    while True:
        ok, f = c.read()
        if not ok:
            break
        F.append(cv2.cvtColor(f, cv2.COLOR_BGR2GRAY))
    return np.array(F)


def silhouette(f, seuil=18):
    s = ndi.binary_fill_holes(f > seuil); s = ndi.binary_opening(s, iterations=2)
    lbl, n = ndi.label(s); t = np.bincount(lbl.ravel()); t[0] = 0
    return (lbl == np.argmax(t)) if n else s


def tour_complet(F):
    def corr(a, b):
        a = a - a.mean(); b = b - b.mean(); return float((a * b).sum() / (np.sqrt((a * a).sum() * (b * b).sum()) + 1e-9))
    return len(F) - 1 if corr(F[0].astype(float), F[-1].astype(float)) > corr(F[0].astype(float), F[-2].astype(float)) else len(F)


def contour_ordonne(sil):
    cnts, _ = cv2.findContours(sil.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    c = max(cnts, key=cv2.contourArea)[:, 0, :].astype(float)
    d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(c, axis=0), axis=1))]
    return c, d


def profil_bande(f, sil, c, r0=4, r1=16):
    dist = ndi.distance_transform_edt(sil); band = (dist > r0) & (dist < r1)
    ys, xs = np.nonzero(band); _, j = cKDTree(c).query(np.stack([xs, ys], 1))
    acc = np.bincount(j, weights=f[ys, xs].astype(float), minlength=len(c)); cnt = np.bincount(j, minlength=len(c))
    return ndi.gaussian_filter1d(acc / np.maximum(cnt, 1), 2, mode='wrap')


def chaine_contour(prof, p_min=12, p_max=60, fen=500):
    w_min = max(3, int(0.3 * (p_min + p_max) / 2))
    det = prof - ndi.gaussian_filter1d(prof, 40, mode='wrap')
    n = len(det); force = np.zeros(n); per = np.zeros(n)
    for i in range(0, n, 10):
        seg = np.take(det, np.arange(i - fen // 2, i + fen // 2), mode='wrap'); seg = seg - seg.mean()
        ac = np.correlate(seg, seg, 'full')[len(seg) - 1:]; ac /= ac[0] + 1e-9
        cand = [k for k in range(p_min, min(p_max, len(ac) - 1)) if ac[k] > ac[k - 1] and ac[k] >= ac[k + 1]]
        if cand:
            k = max(cand, key=lambda q: ac[q]); force[i:i + 10] = ac[k]; per[i:i + 10] = k
    w = w_min
    mins = [i for i in range(n) if det[i] <= np.take(det, np.arange(i - w, i + w + 1), mode='wrap').min()]
    return det, force, per, mins


def recalage_2d(sil, E, A_init=None, rot_max_deg=None):
    """Similitude 2D (avec miroir) amenant la silhouette (x, y image) sur la projection sagittale E (lignes = ax1, colonnes = ax2)
    du pipeline. Init par moments (centroide, orientation, echelle par l'aire), 4 combinaisons de miroir, ICP sur contours.
    Retourne (A 2x3 : [ax2, ax1] = A @ [x, y, 1], erreur mediane px)."""
    def cont(m):
        cs, _ = cv2.findContours(m.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        return max(cs, key=cv2.contourArea)[:, 0, :].astype(float)
    Ps = cont(sil); Pe = cont(E)                     # (x, y) et (col = ax2, row = ax1)
    def moments(P, m):
        M = cv2.moments(m.astype(np.uint8)); cx, cy = M['m10'] / M['m00'], M['m01'] / M['m00']
        cov = np.array([[M['mu20'], M['mu11']], [M['mu11'], M['mu02']]]) / M['m00']
        w, V = np.linalg.eigh(cov); return np.array([cx, cy]), V[:, ::-1], np.sqrt(M['m00'])
    cs_, Vs, ss = moments(Ps, sil); ce, Ve, se = moments(Pe, E)
    arbre = cKDTree(Pe); meilleur = None
    inits = []
    if A_init is not None:
        inits.append(A_init.copy())
    elif rot_max_deg is not None:
        # video et pipeline ont la tete en haut : rotation limitee ; init = identite (et miroir horizontal), echelle par les aires
        s = se / ss
        for f1 in (1, -1):
            R = np.array([[f1, 0.0], [0.0, 1.0]])
            inits.append(np.hstack([s * R, (ce - s * R @ cs_)[:, None]]))
    else:
        for f1 in (1, -1):
            for f2 in (1, -1):
                V2 = Vs.copy(); V2[:, 0] *= f1; V2[:, 1] *= f2
                R = Ve @ V2.T; s = se / ss
                inits.append(np.hstack([s * R, (ce - s * R @ cs_)[:, None]]))
    for A in inits:
            s = float(np.sqrt(abs(np.linalg.det(A[:, :2]))))
            for it in range(40):
                Q = (A[:, :2] @ Ps.T).T + A[:, 2]
                if not np.all(np.isfinite(Q)):
                    break
                d, j = arbre.query(Q); sel = d < np.percentile(d, 85)
                if sel.sum() < 20:
                    break
                X = Ps[sel]; Y = Pe[j[sel]]
                mx, my = X.mean(0), Y.mean(0); Xc, Yc = X - mx, Y - my
                U, D, Vt = np.linalg.svd(Yc.T @ Xc); Sg = np.eye(2); Sg[1, 1] = np.sign(np.linalg.det(U @ Vt))
                R2 = U @ Sg @ Vt; s2 = s                    # echelle fixee par le rapport des aires (l'ICP a echelle libre s'effondre)
                A = np.hstack([s2 * R2, (my - s2 * R2 @ mx)[:, None]])
            if not np.all(np.isfinite(A)):
                continue
            if rot_max_deg is not None:
                Rn = A[:, :2] / np.sqrt(abs(np.linalg.det(A[:, :2]))); ang = abs(np.degrees(np.arctan2(Rn[1, 0], Rn[0, 0])))
                if min(ang, 180 - ang) > rot_max_deg and ang > rot_max_deg:
                    continue
            Q = (A[:, :2] @ Ps.T).T + A[:, 2]; err = float(np.median(arbre.query(Q)[0]))
            if meilleur is None or err < meilleur[1]:
                meilleur = (A, err)
    return meilleur


def construire(dossier):
    out = os.path.join(dossier, 'out', 'topographie', 'video360'); os.makedirs(out, exist_ok=True)
    man = json.load(open(os.path.join(dossier, 'out', 'manifest.json'), encoding='utf-8')); cs = man['stage']; vox = man['mm_per_voxel']
    v = trouver_video(dossier)
    if v is None:
        raise SystemExit(f'{cs} : pas de video de rotation (fond noir) dans {dossier}')
    F = charger_images(v); n_tour = tour_complet(F)
    print(f'{cs} : {len(F)} images {F.shape[2]}x{F.shape[1]}, tour complet sur {n_tour} images')
    z = np.load(os.path.join(dossier, 'work', 'labels.npz')); shape = tuple(int(q) for q in z['shape'])
    env = np.unpackbits(z['enveloppe'])[:int(np.prod(shape))].reshape(shape).astype(bool)
    E = env.any(axis=0); x_med = float(np.median(np.nonzero(env)[0]))          # projection sagittale (ax1 lignes, ax2 colonnes)
    # bout de la queue = voxel de l'enveloppe le plus loin (geodesique) du sommet de la tete
    from topographie_axe import bout_geodesique
    dse = max(1, int(round(0.06 / vox))); env_d = env[::dse, ::dse, ::dse]
    gl = man.get('greatest_length_mm_assumed', 10.0)
    if gl >= 14 and 'membres' in z.files:                 # a partir de CS17 le point le plus loin de la tete serait un pied : on retire les membres
        mem = np.unpackbits(z['membres'])[:int(np.prod(shape))].reshape(shape).astype(bool)[::dse, ::dse, ::dse]
        env_d = env_d & ~ndi.binary_dilation(mem, iterations=2)
        lbl_e, n_e = ndi.label(env_d); t_e = np.bincount(lbl_e.ravel()); t_e[0] = 0; env_d = lbl_e == np.argmax(t_e)
    ie = np.nonzero(env_d); top = int(np.argmin(ie[1])); start = (ie[0][top], ie[1][top], ie[2][top])
    tip_d, _ = bout_geodesique(env_d, start, pont=2); tip_vox = np.array(tip_d, float) * dse; start_vox = np.array(start, float) * dse
    print(f'  bout de la queue (pipeline) : voxel {tip_vox.round(0)} = {(tip_vox * vox).round(2)} mm')
    # vues de profil : largeur maximale et vue opposee
    S = [silhouette(F[k]) for k in range(n_tour)]
    ws = [np.ptp(np.nonzero(s.any(axis=0))[0]) for s in S]; kL = int(np.argmax(ws)); kO = (kL + n_tour // 2) % n_tour
    fa = os.path.join(dossier, 'out', 'topographie', 'axe_vertebral.json')
    axe = json.load(open(fa, encoding='utf-8')) if os.path.exists(fa) else None
    Qa = np.array([a.get('xyz_corps_mm', a['xyz_mm']) for a in axe['axe']]) if axe else None; Sa = np.array([a['s_mm'] for a in axe['axe']]) if axe else None
    ej = json.load(open(os.path.join(dossier, 'out', 'topographie', 'etages_vertebraux.json'), encoding='utf-8')) if axe else {}
    et = ej.get('etages', {}); fiab = ej.get('fiabilite', {})
    vues = []; planches = []; billes = []
    # recalage independant des deux vues, puis coherence : la vue opposee est l'image miroir (autour de l'axe de rotation x0)
    # de la premiere ; on ancre sur la vue la mieux recalee et on re-derive l'autre depuis le miroir (evite le tete-beche
    # des silhouettes en C, presque symetriques)
    x0 = float(np.mean([np.nonzero(s_.any(axis=0))[0].mean() for s_ in S]))
    reg = {k: recalage_2d(S[k], E) for k in (kL, kO)}
    k_anc = min(reg, key=lambda q: reg[q][1]); k_aut = kO if k_anc == kL else kL
    M = np.array([[-1.0, 0.0, 2 * x0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    A_der = (np.vstack([reg[k_anc][0], [0, 0, 1]]) @ M)[:2]
    reg_der = recalage_2d(S[k_aut], E, A_init=A_der)
    print(f'  recalage : vue {k_anc} ancre ({reg[k_anc][1]:.1f} px) ; vue {k_aut} independante {reg[k_aut][1]:.1f} px, derivee du miroir {reg_der[1]:.1f} px -> derivee retenue')
    reg[k_aut] = reg_der
    for k in (kL, kO):
        f = F[k].astype(float); sil = S[k]
        A, err = reg[k]; mm_px = float(np.sqrt(abs(np.linalg.det(A[:, :2]))) * vox)
        Ainv = np.linalg.inv(np.vstack([A, [0, 0, 1]]))
        def vers_video(p_vox):
            return Ainv[:2, :2] @ np.array([p_vox[2], p_vox[1]]) + Ainv[:2, 2]
        c, d = contour_ordonne(sil); prof = profil_bande(f, sil, c)
        # fenetre de periode autour de la periode connue de l'axe (etages_vertebraux.json), sinon large
        det, force, per, mins = chaine_contour(prof, p_min=8, p_max=60)      # fenetre large : la video a sa propre periode
        bons = [i for i in mins if force[i] >= 0.3 and per[i] > 0]
        fr = []
        for i in bons:
            p = A[:, :2] @ c[i] + A[:, 2]                 # (ax2, ax1) en voxels
            p_vox = np.array([x_med, p[1], p[0]]); p_mm = p_vox * vox
            s_axe = None
            if Qa is not None:
                j = int(np.argmin(np.linalg.norm(Qa - p_mm, axis=1))); dist_axe = float(np.linalg.norm(Qa[j] - p_mm))
                s_axe = {'s_mm': round(float(Sa[j]), 2), 'dist_mm': round(dist_axe, 2), 'au_dela': bool(j == len(Qa) - 1 and dist_axe > 0.3)}
            fr.append({'arc_px': round(float(d[i]), 1), 'xy_px': [round(float(c[i][0]), 1), round(float(c[i][1]), 1)], 'force': round(float(force[i]), 2),
                       'periode_px': float(per[i]), 'xyz_mm_pipeline': [round(float(q), 3) for q in p_mm], 'axe': s_axe})
        # ---- prolongement vers la queue : etages de l'axe projetes sur le contour, puis grille au pas de l'axe (en px)
        # jusqu'au bout de la queue, calee sur les minima video (une frontiere video sur deux : la video a une periode moitie)
        prolong = []
        recalage_ok = err <= 25
        if Qa is not None and et:
            # controle du sens du recalage (une silhouette en C se recale parfois tete-beche) : le sommet de la tete projete
            # doit etre plus pres du premier etage que du dernier
            # test independant du recalage : dans la VIDEO, la tete est epaisse et la queue fine ; l'epaisseur de la silhouette
            # (distance au bord) au point ou se projette le sommet de la tete doit depasser celle au bout de la queue
            dist_v = ndi.distance_transform_edt(sil)
            def ep(xy):
                yi, xi = int(np.clip(xy[1], 0, sil.shape[0] - 1)), int(np.clip(xy[0], 0, sil.shape[1] - 1)); return float(dist_v[yi, xi])
            e_tete = ep(vers_video(start_vox + np.array([0, 0.6 / vox, 0]))); e_queue = ep(vers_video(tip_vox))
            print(f'  image {k} : epaisseur silhouette a la tete {e_tete:.0f} px (projetee en {vers_video(start_vox).round(0)}), a la queue {e_queue:.0f} px (projetee en {vers_video(tip_vox).round(0)})')
            if e_queue > e_tete:
                print(f'  image {k} : recalage tete-beche (epaisseur au bout de queue {e_queue:.0f} px > tete {e_tete:.0f} px) : vue ecartee'); recalage_ok = False
        if Qa is not None and et and recalage_ok:                        # seulement sur une vue bien recalee et dans le bon sens
            arbre_c = cKDTree(c)
            idx_et = []
            for nom, e in et.items():
                jj = int(np.argmin(np.abs(Sa - e['s_debut_mm']))); xy = vers_video(Qa[jj] / vox)
                dd, ii = arbre_c.query(xy); idx_et.append((ii, dd))
            i_tip = int(arbre_c.query(vers_video(tip_vox))[1])
            n_c = len(c)
            if len(idx_et) >= 3:
                # sens de progression le long du contour (indices croissants ou decroissants) et dernier etage
                i0, i1 = idx_et[0][0], idx_et[-1][0]
                sens_c = 1 if ((i1 - i0) % n_c) < ((i0 - i1) % n_c) else -1
                i_last = idx_et[-1][0]
                if fiab.get('fiable') and fiab.get('periode_mm'):
                    pas_px = fiab['periode_mm'] / mm_px; origine_pas = 'axe'
                else:                                          # axe non fiable : la video a une periodicite moitie du pas (2 traits par bloc)
                    pas_px = float(np.median(per[bons]) * 2) if bons else 0.0; origine_pas = 'video x2 (approximatif)'
                # longueur d'arc du dernier etage au bout de la queue dans le sens de progression
                arc_tot = ((i_tip - i_last) * sens_c) % n_c
                # profondeur du dernier etage sous le contour : les etages de queue sont rentres d'autant (normale interieure)
                dist_t = ndi.distance_transform_edt(sil); gy, gx = np.gradient(dist_t)
                jj_last = int(np.argmin(np.abs(Sa - list(et.values())[-1]['s_debut_mm']))); xy_last = vers_video(Qa[jj_last] / vox)
                prof_last = float(dist_t[int(np.clip(xy_last[1], 0, sil.shape[0] - 1)), int(np.clip(xy_last[0], 0, sil.shape[1] - 1))])
                prof_last = float(np.clip(prof_last, 4, 40))
                def rentrer(i_c):
                    x, y = c[i_c]; yi, xi = int(np.clip(y, 0, sil.shape[0] - 1)), int(np.clip(x, 0, sil.shape[1] - 1))
                    g = np.array([gx[yi, xi], gy[yi, xi]]); n_ = np.linalg.norm(g)
                    return c[i_c] + (g / n_ * prof_last if n_ > 1e-6 else 0)
                if arc_tot > n_c / 2:
                    print(f'  image {k} : bout de la queue derriere le dernier etage (arc {arc_tot} > demi-contour) : pas de prolongement')
                elif pas_px > 4 and arc_tot > pas_px:
                    mins_set = np.array(sorted(mins))
                    n_new = min(int(arc_tot // pas_px), 12)      # au plus S1-S5 + Co1-4 (+ marge)
                    for m_ in range(1, n_new + 1):
                        i_pred = int((i_last + sens_c * m_ * pas_px) % n_c)
                        # calage sur le minimum video le plus proche (± 0,3 pas)
                        dmin = np.abs(((mins_set - i_pred + n_c // 2) % n_c) - n_c // 2)
                        jm = int(np.argmin(dmin)); i_use = int(mins_set[jm]) if dmin[jm] <= 0.3 * pas_px else i_pred
                        xy_in = rentrer(i_use); p2 = A[:, :2] @ xy_in + A[:, 2]; p_vox = np.array([x_med, p2[1], p2[0]])
                        prolong.append({'i_contour': i_use, 'cale': bool(dmin[jm] <= 0.3 * pas_px), 'xy_px': [round(float(xy_in[0]), 1), round(float(xy_in[1]), 1)],
                                        'xyz_mm_pipeline': [round(float(q), 3) for q in p_vox * vox]})
                        billes.append(p_vox)
                    print(f"  image {k} : {n_new} etage(s) supplementaire(s) vers la queue (pas {pas_px:.1f} px [{origine_pas}], arc {arc_tot:.0f} px), {sum(pp['cale'] for pp in prolong)} cale(s) sur la video")
        vues.append({'image': int(k), 'angle_deg': round(360.0 * k / n_tour, 1), 'recalage_erreur_px': round(err, 1), 'recalage_ok': bool(recalage_ok), 'mm_par_px': round(mm_px, 5),
                     'etages_supplementaires_queue': prolong, 'n_supplementaires': len(prolong),
                     'n_frontieres': len(bons), 'force_moyenne': round(float(np.mean(force[bons])), 2) if bons else 0.0,
                     'periode_px_mediane': round(float(np.median(per[bons])), 1) if bons else None,
                     'periode_mm': round(float(np.median(per[bons])) * mm_px, 3) if bons else None, 'frontieres': fr})
        img = cv2.cvtColor(F[k], cv2.COLOR_GRAY2BGR)
        if Qa is not None:
            for nom, e in et.items():
                jj = int(np.argmin(np.abs(Sa - e['s_debut_mm']))); q = Qa[jj] / vox      # voxels (ax0, ax1, ax2)
                xy = Ainv[:2, :2] @ np.array([q[2], q[1]]) + Ainv[:2, 2]
                cv2.circle(img, (int(xy[0]), int(xy[1])), 3, (255, 120, 0), 1)
        for i in mins:
            ok = force[i] >= 0.3 and per[i] > 0
            cv2.circle(img, (int(c[i][0]), int(c[i][1])), 3 if ok else 1, (0, 140, 255) if ok else (90, 90, 90), -1)
        for pp in prolong:
            cv2.circle(img, (int(pp['xy_px'][0]), int(pp['xy_px'][1])), 5, (0, 255, 0) if pp['cale'] else (0, 160, 0), 1)
        if Qa is not None:
            xt = vers_video(tip_vox); cv2.drawMarker(img, (int(xt[0]), int(xt[1])), (255, 0, 255), cv2.MARKER_CROSS, 14, 1)
        cv2.putText(img, f'{cs} image {k} ({360.0 * k / n_tour:.0f} deg) : {len(bons)} frontieres video (orange), etages axe (bleu), {len(prolong)} etages queue (vert), bout (magenta) ; periode {vues[-1]["periode_px_mediane"]} px = {vues[-1]["periode_mm"]} mm ; recalage {err:.0f} px, {mm_px * 1000:.1f} um/px',
                    (6, 16), cv2.FONT_HERSHEY_SIMPLEX, 0.4, (0, 255, 255), 1)
        planches.append(img)
        # controle du recalage : contour de la silhouette transforme sur E
        if k == kL:
            Q = (A[:, :2] @ c.T).T + A[:, 2]
            ctrl = cv2.cvtColor((E.astype(np.uint8) * 90), cv2.COLOR_GRAY2BGR)
            for q in Q[::3]:
                yy, xx = int(round(q[1])), int(round(q[0]))
                if 0 <= yy < ctrl.shape[0] and 0 <= xx < ctrl.shape[1]:
                    ctrl[yy, xx] = (0, 0, 255)
            for i in bons:
                q = A[:, :2] @ c[i] + A[:, 2]; cv2.circle(ctrl, (int(q[0]), int(q[1])), 3, (0, 255, 255), -1)
            if Qa is not None:
                for q in Qa[::3] / vox:
                    cv2.circle(ctrl, (int(q[2]), int(q[1])), 1, (255, 120, 0), -1)
            sc = min(1.0, 900 / ctrl.shape[0]); ctrl = cv2.resize(ctrl, None, fx=sc, fy=sc)
            cv2.putText(ctrl, f'{cs} : enveloppe pipeline (gris), silhouette video recalee (rouge), frontieres video (jaune), axe (bleu) ; erreur {err:.0f} px', (6, 16), cv2.FONT_HERSHEY_SIMPLEX, 0.4, (0, 255, 255), 1)
            cv2.imwrite(os.path.join(out, 'controle_recalage.png'), ctrl)
    cv2.imwrite(os.path.join(out, 'profil_somites.png'), np.hstack(planches))
    json.dump({'stage': cs, 'video': os.path.basename(v), 'n_images': int(len(F)), 'n_tour': int(n_tour), 'vues': vues},
              open(os.path.join(out, 'somites_video.json'), 'w', encoding='utf-8'), indent=1)
    # billes PLY dans le repere des PLY du pipeline
    if billes:
        import trimesh
        c0 = np.array(man['center_voxel']); cc = np.array([c0[0], c0[2], -c0[1]]) * vox
        r = float(np.clip(0.012 * man.get('greatest_length_mm_assumed', 10), 0.06, 0.2))
        ms = []
        for p in billes:
            P = np.array(p) * vox; Pb = np.array([P[0], P[2], -P[1]]) - cc
            sp = trimesh.creation.icosphere(subdivisions=1, radius=r); sp.apply_translation(Pb); ms.append(sp)
        m = trimesh.util.concatenate(ms); fn = f'{cs}_somites_video.ply'; m.export(os.path.join(out, fn))
    fm = os.path.join(dossier, 'out', 'topographie', 'manifest_topographie.json')
    if os.path.exists(fm):
        m = json.load(open(fm, encoding='utf-8'))
        m['controle']['video360'] = {'video': os.path.basename(v), 'n_images': int(len(F)),
                                     'vues_profil': [{kk: vv for kk, vv in sv.items() if kk != 'frontieres'} for sv in vues],
                                     'images': ['video360/profil_somites.png', 'video360/controle_recalage.png']}
        m['controle']['images'] = list(dict.fromkeys(m['controle'].get('images', []) + ['video360/profil_somites.png', 'video360/controle_recalage.png']))
        if billes:
            m['segments'] = [sg for sg in m['segments'] if sg['name'] != 'somites_video'] + [
                {'name': 'somites_video', 'file': f'video360/{cs}_somites_video.ply', 'color': [1.0, 0.6, 0.1], 'alpha': 1.0, 'faces': int(len(ms) * 80), 'volume_mm3': 0.0,
                 'confiance': 'video', 'note': f"etages supplementaires vers la queue, au pas de l axe, cales sur la video 360 (vues {[sv['image'] for sv in vues]}), plan median"}]
        json.dump(m, open(fm, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    for sv in vues:
        print(f"  image {sv['image']} ({sv['angle_deg']} deg) : recalage {sv['recalage_erreur_px']} px, {sv['mm_par_px'] * 1000:.1f} um/px ; {sv['n_frontieres']} frontieres, force {sv['force_moyenne']}, periode {sv['periode_px_mediane']} px = {sv['periode_mm']} mm")
    print('->', out)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('dossier'); a = ap.parse_args(); construire(a.dossier)
