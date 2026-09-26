"""Candidats du pharynx (intestin pharyngien) et des poches pharyngiennes, à valider avant de les copier dans digestif_points/<CS>.json.

Principe : le pharynx est la lumière pâle qui prolonge l'œsophage vers le haut jusqu'à la cavité buccale. On part du début (crânial) de
l'œsophage tracé et on prend la composante connexe de lumière pâle qui le contient, au-dessus de lui (seuil calé sur la lumière de
l'œsophage ; exclus : hors enveloppe, cœur, cavité péricardique, vaisseaux, ventricules cérébraux, yeux, vésicules otiques). Le pharynx est
le chemin de moindre coût dans cette lumière (ligne centrale favorisée) entre le début de l'œsophage et le point le plus éloigné proche de
la ligne médiane, c'est-à-dire le fond de la cavité buccale. Il est écrit de haut en bas (cavité buccale → œsophage), comme dans
digestif_points, et finit sur le premier point de l'œsophage.
Poches : le long du pharynx, de chaque côté, l'extension latérale de la lumière passe par un maximum à chaque poche (entre deux poches,
l'arc pharyngien la resserre). Chaque poche va de l'axe du pharynx (début) à son fond latéral (fin), par la lumière. Numérotation
crânio-caudale d'après vaisseaux_points/calendrier_arcs.json quand le nombre de poches d'un côté égale le nombre attendu au stade
(présentes + en formation + en régression, sinon présentes seules), sinon poche_candidat_<i>_<côté>. Côté : celui de l'aorte dorsale
gauche tracée (x plus petit dans les tracés actuels), sinon x plus petit = gauche.

Entrées (<work> = dossier work du stade) : dens.npy, digestif/chemins.json (œsophage recalé ; sinon points de digestif_points/<CS>.json) ;
facultatifs : labels.npz (enveloppe, cœur, cavité péricardique, ventricules, yeux, vésicules otiques), cardio/coeur.npz, cardio/vaisseaux.npz,
cardio/vaisseaux_chemins.json (haut de la boîte, côté gauche).
Sorties : embryo3d/digestif_points/<CS>_pharynx_proposes.json (format de digestif_points/<CS>.json + « qualite », « a_valider ») et
<work>/digestif/pharynx_candidats.png (face et profil médian en projection minimale, contour de la lumière retenue, pharynx, poches, œsophage).

usage : python embryo3d/pharynx_candidats.py <work> <CS> [--seuil N] [--marge 60] [--mediane 15] [--prominence 6] [--pas 1] [--sans-poches]
                                            [--sortie fichier.json] [--planche fichier.png]
        python embryo3d/pharynx_candidats.py fusionner <CS>_pharynx_proposes.json [digestif_points/<CS>.json]
          → copie pharynx et poche_pharyngienne_* (pas les poche_candidat_* non renommées) dans le fichier de points du stade, en remplaçant
            ceux de même nom, sans les blocs « qualite »."""
import numpy as np, json, os, sys, argparse
from scipy import ndimage as ndi
from scipy.signal import find_peaks

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from arcs_candidats import masque, union_npz, aortes, attendus, echantillonne, fusionner, COTES  # noqa: E402
POCHES = (1, 2, 3, 4)

def oesophage(work, cs):
    """centres (voxels, début = crânial), rayons, source"""
    p = os.path.join(work, 'digestif', 'chemins.json')
    if os.path.exists(p):
        ch = json.load(open(p, encoding='utf-8')).get('oesophage', {})
        if 'centres' in ch:
            C, R = np.array(ch['centres'], float), np.array(ch['rayons'], float)
            return (C[::-1], R[::-1], 'digestif/chemins.json') if C[0, 1] > C[-1, 1] else (C, R, 'digestif/chemins.json')
    p = os.path.join(HERE, 'digestif_points', f'{cs}.json')
    if os.path.exists(p):
        seg = next((s for s in json.load(open(p, encoding='utf-8'))['segments'] if s['nom'] == 'oesophage'), None)
        if seg:
            C = np.array([q[:3] for q in seg['points']], float); C = C[::-1] if C[0, 1] > C[-1, 1] else C
            return C, np.full(len(C), seg.get('rmax', 8) / 2.0), f'digestif_points/{cs}.json (points non recalés)'
    return None, None, None

def proposer(work, cs, seuil=None, marge=60, mediane=15, prominence=6.0, pas=1, poches=True):
    from skimage.graph import MCP_Geometric
    dens = np.load(os.path.join(work, 'dens.npy'), mmap_mode='r'); SH = dens.shape
    Coe, Roe, src_oe = oesophage(work, cs)
    if Coe is None: sys.exit('œsophage non tracé (digestif/chemins.json ou digestif_points/<CS>.json) : le tracer d\'abord')
    E0 = Coe[0]; print('début de l\'œsophage', np.round(E0).astype(int).tolist(), '(' + src_oe + ')')
    try: aor = aortes(work)
    except (OSError, ValueError): aor = {}
    pl = os.path.join(work, 'labels.npz'); Lz = np.load(pl) if os.path.exists(pl) else None
    env = masque(Lz, 'enveloppe', SH) if Lz is not None and 'enveloppe' in Lz.files else None
    # boîte : du haut des aortes (ou 3 marges au-dessus de l'œsophage) à une demi-marge sous le début de l'œsophage ; en x et y, l'enveloppe
    s_top = int(E0[1]) - 3 * marge
    if aor: s_top = min(s_top, int(min(C[0, 1] for C, _ in aor.values())) - marge)
    s_top, s_bot = max(s_top, 0), min(int(E0[1]) + marge // 2, SH[1])
    if env is not None and env[:, s_top:s_bot, :].any():
        e = env[:, s_top:s_bot, :]; xs = np.where(e.any(axis=(1, 2)))[0]; ys = np.where(e.any(axis=(0, 1)))[0]
        lo, hi = np.array([xs.min(), s_top, ys.min()]), np.array([xs.max() + 1, s_bot, ys.max() + 1])
    else:
        print('ATTENTION : pas d\'enveloppe dans labels.npz, l\'extérieur de l\'embryon (pâle) n\'est pas exclu')
        c = np.round(E0).astype(int); lo = np.maximum(c - 3 * marge, 0); hi = np.minimum(c + 3 * marge + 1, SH); lo[1], hi[1] = s_top, s_bot
    sub = np.asarray(dens[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]).astype(np.float32)
    D = ndi.gaussian_filter(sub, max(1.0, 0.6 * pas))[::pas, ::pas, ::pas]; del sub
    loc = lambda X: np.clip(np.round((np.asarray(X, float) - lo) / pas).astype(int), 0, np.array(D.shape) - 1)
    glob = lambda I: np.asarray(I) * pas + lo
    print('boîte (x, s, y)', lo.tolist(), '→', hi.tolist(), 'pas', pas, 'voxels', D.size)
    # exclusions
    def reduit(m): return m[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]][::pas, ::pas, ::pas]
    interdit = np.zeros(D.shape, bool); notes = []
    if env is not None: interdit |= ~ndi.binary_erosion(reduit(env), iterations=1); notes.append('hors enveloppe')
    for k, dil in (('cavite_pericardique', 0), ('coeur_detoure', 1), ('coeur', 1), ('cavites_cardiaques', 1), ('ventricules', 1), ('yeux', 1),
                   ('vesicules_otiques', 1), ('cristallins', 0)):
        if Lz is not None and k in Lz.files:
            m = reduit(masque(Lz, k, SH)); interdit |= ndi.binary_dilation(m, iterations=dil) if dil else m; notes.append(k)
    for f, garder in (('cardio/coeur.npz', lambda k: k == 'coeur_plein'), ('cardio/vaisseaux.npz', lambda k: True)):
        m, cles = union_npz(os.path.join(work, f), garder, SH)
        if m is not None: interdit |= ndi.binary_dilation(reduit(m), iterations=1); notes.append(f + ':' + ','.join(cles))
    # seuil calé sur la lumière de l'œsophage (centres recalés dans la boîte)
    Cb = Coe[(Coe[:, 1] >= lo[1]) & (Coe[:, 1] < hi[1])]
    lum = float(np.median(D[tuple(loc(Cb if len(Cb) else E0[None]).T)])); dedans = ~interdit
    tissu = float(np.percentile(D[dedans], 75)) if dedans.any() else float(np.percentile(D, 75))
    if seuil is None: seuil = lum + 0.35 * (tissu - lum)
    print(f'densité : lumière de l\'œsophage {lum:.0f}, tissu (75e centile) {tissu:.0f}, seuil de lumière {seuil:.0f}')
    if tissu - lum < 15: print('ATTENTION : faible contraste lumière / tissu (œsophage collabé ?), donner --seuil')
    pale = (D <= seuil) & ~interdit
    g = tuple(loc(E0))
    if not pale[g]:                                                        # graine dans la lumière la plus proche (8 voxels au plus)
        cand = np.argwhere(pale); dd = np.linalg.norm((cand - np.array(g)) * pas, axis=1)
        if not len(cand) or dd.min() > 8: sys.exit('pas de lumière pâle au début de l\'œsophage : vérifier le tracé ou --seuil')
        g = tuple(cand[np.argmin(dd)])
    lab, _ = ndi.label(pale); comp = lab == lab[g]; del lab
    nvox = int(comp.sum()) * pas ** 3; bords = [comp[0].any(), comp[-1].any(), comp[:, 0].any(), comp[:, :, 0].any(), comp[:, :, -1].any()]
    print(f'lumière retenue : {nvox} voxels' + (' ; ATTENTION : elle touche le bord de la boîte (fuite ? augmenter --marge ou baisser --seuil)' if any(bords) else ''))
    # plan médian x(s) : milieu des aortes dorsales paires à chaque hauteur (embryon symétrique), sinon x du début de l'œsophage ;
    # le pharynx est aplati (large en x, mince en y) : la centralité ne le centre pas en x, d'où une pénalité d'écart au plan médian
    s_g = lo[1] + pas * np.arange(D.shape[1])
    if {'gauche', 'droite'} <= set(aor):
        (G, _), (Dr, _) = aor['gauche'], aor['droite']
        G, Dr = G[np.argsort(G[:, 1])], Dr[np.argsort(Dr[:, 1])]
        x_mid = (np.interp(s_g, G[:, 1], G[:, 0]) + np.interp(s_g, Dr[:, 1], Dr[:, 0])) / 2; src_mid = 'milieu des aortes dorsales'
    else: x_mid = np.full(len(s_g), E0[0]); src_mid = 'x du début de l\'œsophage'
    dx_mid = (lo[0] + pas * np.arange(D.shape[0]))[:, None, None] - x_mid[None, :, None]    # (x, s, 1) en voxels de dens.npy
    print(f'plan médian : {src_mid} (x de {x_mid.min():.0f} à {x_mid.max():.0f})')
    edt = ndi.distance_transform_edt(comp)
    cout_base = np.where(comp, 1.0 + 3.0 / np.maximum(edt, 0.5), np.inf)
    cout = cout_base * (1.0 + (dx_mid / 8.0) ** 2)
    mcp = MCP_Geometric(cout, fully_connected=True); A, _ = mcp.find_costs([g])
    for tol in sorted({3.0, 6.0, float(mediane)}):                       # fond de la cavité buccale : point le plus éloigné, bande médiane stricte d'abord
        med = comp & np.isfinite(A) & (np.abs(dx_mid) <= tol)
        if med.any(): break
    if not med.any(): sys.exit('aucune lumière médiane au-dessus de l\'œsophage (plan médian ? --mediane)')
    e = np.unravel_index(np.argmax(np.where(med, A, -1.0)), A.shape)
    Ph = np.array(mcp.traceback(e))[::-1]                                   # du fond de la cavité buccale (haut) au début de l'œsophage
    Pg = glob(Ph).astype(float)
    if np.linalg.norm(Pg[-1] - E0) > 1: Pg = np.vstack([Pg, E0])
    rl = float(edt[tuple(Ph.T)].max()) * pas
    segments = [{'nom': 'pharynx', 'rmax': int(np.clip(round(2.5 * rl), 8, 25)), 'rayon_recalage': 4, 'points': echantillonne(Pg, 8.0),
                 'qualite': {'longueur_vox': round(float(np.linalg.norm(np.diff(Pg, axis=0), axis=1).sum()), 1), 'demi_largeur_lumiere_max_vox': round(rl, 1),
                             'fond_cavite_buccale': glob(e).astype(int).tolist(), 'lumiere_vox': nvox, 'touche_bord_de_boite': bool(any(bords))}}]
    q = segments[0]['qualite']
    print(f"pharynx : longueur {q['longueur_vox']:.0f} vox, demi-largeur de lumière max {rl:.1f} vox, fond de la cavité buccale {q['fond_cavite_buccale']}")
    # poches : maxima de l'extension latérale de la lumière le long de l'axe du pharynx, de chaque côté
    if poches and not any(attendus(cs, c, 'poche_pharyngienne', POCHES)[0] for c in COTES):
        print('aucune poche attendue à ce stade (calendrier) : pas de recherche'); poches = False
    if poches:
        axe_id = np.full(D.shape, -1, np.int32); axe_id[tuple(Ph.T)] = np.arange(len(Ph))
        ind = np.zeros((3,) + D.shape, np.int32)
        ndi.distance_transform_edt(axe_id < 0, return_distances=False, return_indices=True, indices=ind)
        vox = np.argwhere(comp); t = axe_id[tuple(ind[:, vox[:, 0], vox[:, 1], vox[:, 2]])]; del ind
        dx = dx_mid[vox[:, 0], vox[:, 1], 0]                                # écart latéral au plan médian (voxels de dens.npy)
        sgn_g = np.sign(np.mean(aor['gauche'][0][:, 0]) - np.mean(aor['droite'][0][:, 0])) if {'gauche', 'droite'} <= set(aor) else -1.0
        m2 = MCP_Geometric(cout_base, fully_connected=True); m2.find_costs([tuple(p) for p in Ph])
        w = max(1, int(round(3 / pas)))
        for c in COTES:
            sgn = sgn_g if c == 'gauche' else -sgn_g; sel = dx * sgn > 0
            ext = np.zeros(len(Ph)); np.maximum.at(ext, t[sel], np.abs(dx[sel]))
            ext = ndi.uniform_filter1d(ndi.maximum_filter1d(ext, w), 5)      # comble les indices d'axe sans voxel latéral, lisse
            # une poche dépasse nettement la demi-largeur du corps du pharynx (médiane du profil) et ses voisines
            pics, props = find_peaks(ext, height=np.median(ext) + prominence, prominence=prominence, distance=max(2, int(10 / pas)))
            tous_c, presents = attendus(cs, c, 'poche_pharyngienne', POCHES)
            if len(pics) == len(tous_c): noms, regle = [f'poche_pharyngienne_{k}_{c}' for k in tous_c], 'nombre = poches attendues'
            elif len(pics) == len(presents): noms, regle = [f'poche_pharyngienne_{k}_{c}' for k in presents], 'nombre = poches « présentes » du calendrier'
            else: noms, regle = [f'poche_candidat_{i + 1}_{c}' for i in range(len(pics))], f'{len(pics)} poches pour {len(tous_c)} attendues ({tous_c}) : à nommer'
            print(f'{c} : {len(pics)} poche(s) — {regle}')
            for nom, p, prom in zip(noms, pics, props['prominences']):
                zone = sel & (np.abs(t - p) <= w); iz = np.where(zone)[0]
                pointe = tuple(vox[iz[np.argmax(np.abs(dx[iz]))]])
                ch = np.array(m2.traceback(pointe)); Qg = glob(ch).astype(float)
                rp = float(edt[tuple(ch.T)].max()) * pas
                qq = {'hauteur_s_axe': int(glob(Ph[p])[1]), 'extension_laterale_vox': round(float(ext[p]), 1), 'prominence': round(float(prom), 1),
                      'fond': glob(pointe).astype(int).tolist(), 'longueur_vox': round(float(np.linalg.norm(np.diff(Qg, axis=0), axis=1).sum()), 1)}
                segments.append({'nom': nom, 'rmax': int(np.clip(round(2.2 * rp), 3, 8)), 'rayon_recalage': 2, 'points': echantillonne(Qg, 5.0), 'qualite': qq})
                print(f"  {nom:26s} axe s={qq['hauteur_s_axe']:4d}  extension {qq['extension_laterale_vox']:5.1f} vox  prominence {qq['prominence']:5.1f}  fond {qq['fond']}")
    res = {'stade': cs, 'a_valider': True, 'propose_par': 'pharynx_candidats.py',
           'validation': 'vérifier le pharynx (du fond de la cavité buccale au début de l\'œsophage) et chaque poche sur la planche, renommer les '
                         'poche_candidat_* (poche_pharyngienne_<k>_<côté>) ou les supprimer, corriger les points au besoin (x, s, y de dens.npy), puis '
                         '« pharynx_candidats.py fusionner » et digestif_build.py <work> embryo3d/digestif_points/<CS>.json',
           'parametres': {'debut_oesophage': np.round(E0).astype(int).tolist(), 'source_oesophage': src_oe, 'seuil': round(seuil, 1), 'lumiere_oesophage': round(lum, 1),
                          'tissu': round(tissu, 1), 'plan_median': src_mid, 'marge': marge, 'mediane': mediane, 'prominence': prominence, 'pas': pas,
                          'boite': [lo.tolist(), hi.tolist()], 'exclusions': notes},
           'segments': segments}
    return res, dict(D=D, lo=lo, pas=pas, comp=comp, oe=Coe, xm=g[0])

# ---------------------------------------------------------------- planche de contrôle
def planche(res, ctx, chemin):
    import cv2
    try:
        from tubes_morph import couleur
    except Exception:
        couleur = lambda n: (0.8, 0.3, 0.3) if n == 'pharynx' else (0.6, 0.6, 0.9)
    D, lo, pas, comp, xm = ctx['D'], ctx['lo'], ctx['pas'], ctx['comp'], ctx['xm']
    z = max(1, int(round(2 / pas)))
    bgr = lambda n: tuple(int(255 * v) for v in couleur(n)[::-1])
    def vue(proj, masq, ax_l, ax_c, titre):
        im = cv2.cvtColor(np.clip(255 - proj, 0, 255).astype(np.uint8), cv2.COLOR_GRAY2BGR)
        im = cv2.resize(im, None, fx=z, fy=z, interpolation=cv2.INTER_NEAREST)
        mm = cv2.resize(masq.astype(np.uint8), (im.shape[1], im.shape[0]), interpolation=cv2.INTER_NEAREST)
        cnt, _ = cv2.findContours(mm, cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE); cv2.drawContours(im, cnt, -1, (0, 140, 255), 1)
        pt = lambda X: (int((X[ax_c] - lo[ax_c]) / pas * z), int((X[ax_l] - lo[ax_l]) / pas * z))
        cv2.polylines(im, [np.array([pt(x) for x in ctx['oe']], np.int32)], False, bgr('oesophage'), 2, cv2.LINE_AA)
        for s in res['segments']:
            P = np.array(s['points']); cc = bgr('pharynx') if s['nom'] == 'pharynx' else (160, 40, 160)   # poches en violet (la légende, pâle, ne se lit pas ici)
            if s['nom'].startswith('poche_candidat'): cc = (0, 90, 230)
            cv2.polylines(im, [np.array([pt(x) for x in P], np.int32)], False, cc, 2, cv2.LINE_AA)
            for x in P: cv2.circle(im, pt(x), 2, cc, -1)
            lab = s['nom'].replace('poche_pharyngienne_', 'P').replace('poche_candidat_', 'P? ').replace('_gauche', ' G').replace('_droite', ' D')
            cv2.putText(im, lab, (pt(P[-1])[0] + 4, pt(P[-1])[1] + 4), 0, 0.4, cc, 1, cv2.LINE_AA)
        h, w = im.shape[:2]
        for ax, n, horiz in ((ax_l, h, False), (ax_c, w, True)):          # graduations tous les 50 voxels de dens.npy
            for v in range(int(np.ceil(lo[ax] / 50.0)) * 50, int(lo[ax] + n / z * pas) + 1, 50):
                p = int((v - lo[ax]) / pas * z)
                if horiz: cv2.line(im, (p, 0), (p, 6), (0, 150, 0), 1); cv2.putText(im, str(v), (p + 2, 16), 0, 0.35, (0, 120, 0), 1)
                else: cv2.line(im, (0, p), (6, p), (0, 150, 0), 1); cv2.putText(im, str(v), (8, p + 4), 0, 0.35, (0, 120, 0), 1)
        cv2.putText(im, titre, (8, h - 8), 0, 0.45, (200, 0, 0), 1, cv2.LINE_AA)
        return im
    ys = np.where(comp.any(axis=(0, 1)))[0]
    face = vue(D[:, :, ys.min():ys.max() + 1].min(axis=2).T, comp.any(axis=2).T, 1, 0, 'face (projection minimale sur la profondeur de la lumiere)')
    a, b = int(max(xm - 20 / pas, 0)), int(min(xm + 20 / pas, D.shape[0] - 1)) + 1
    prof = vue(D[a:b].min(axis=0), comp[a:b].any(axis=0), 1, 2, f'profil median (x {int(a * pas + lo[0])}-{int(b * pas + lo[0])}, projection minimale)')
    H = max(face.shape[0], prof.shape[0])
    tu = [cv2.copyMakeBorder(t, 0, H - t.shape[0], 0, 6, cv2.BORDER_CONSTANT, value=(90, 90, 90)) for t in (face, prof)]
    haut = np.hstack(tu)
    lignes = [f"{res['stade']} - seuil de lumiere {res['parametres']['seuil']:.0f} - contour orange = lumiere retenue, rose = oesophage trace, violet = poches - A VALIDER (x, s, y de dens.npy)"]
    for s in res['segments']:
        q = s['qualite']
        lignes.append(f"{s['nom']:26s} " + (f"longueur {q['longueur_vox']:.0f} vox, fond de la cavite buccale {q['fond_cavite_buccale']}" if s['nom'] == 'pharynx'
                      else f"axe s={q['hauteur_s_axe']}, extension {q['extension_laterale_vox']:.0f} vox, prominence {q['prominence']:.0f}, fond {q['fond']}"))
    g = np.full((16 * len(lignes) + 10, haut.shape[1], 3), 255, np.uint8)
    for i, t in enumerate(lignes): cv2.putText(g, t, (8, 16 + 16 * i), 0, 0.42, (40, 40, 40) if i else (160, 0, 0), 1, cv2.LINE_AA)
    os.makedirs(os.path.dirname(os.path.abspath(chemin)), exist_ok=True)
    cv2.imwrite(chemin, np.vstack([haut, g])); print('planche', chemin)

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'fusionner':
        fusionner(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None, ('pharynx', 'poche_pharyngienne_'), 'digestif_points'); sys.exit()
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('work'); ap.add_argument('cs', help='stade (CS13, CS14…), pour le calendrier et le nom de sortie')
    ap.add_argument('--seuil', type=float, help='densité maximale d\'une lumière ; défaut : calé sur l\'œsophage')
    ap.add_argument('--marge', type=int, default=60, help='marge de la boîte (voxels) au-dessus des aortes / sous le début de l\'œsophage')
    ap.add_argument('--mediane', type=float, default=15, help='écart maximal (voxels, en x) du fond de la cavité buccale à la ligne médiane')
    ap.add_argument('--prominence', type=float, default=6.0, help='saillie latérale minimale (voxels) d\'une poche par rapport à ses voisines')
    ap.add_argument('--pas', type=int, default=1, help='sous-échantillonnage du volume (2 = plus rapide)')
    ap.add_argument('--sans-poches', action='store_true'); ap.add_argument('--sortie'); ap.add_argument('--planche')
    a = ap.parse_args()
    res, ctx = proposer(a.work, a.cs, a.seuil, a.marge, a.mediane, a.prominence, max(1, a.pas), not a.sans_poches)
    sortie = a.sortie or os.path.join(HERE, 'digestif_points', f'{a.cs}_pharynx_proposes.json')
    json.dump(res, open(sortie, 'w', encoding='utf-8'), indent=1, ensure_ascii=False); print('propositions', sortie)
    planche(res, ctx, a.planche or os.path.join(a.work, 'digestif', 'pharynx_candidats.png'))
