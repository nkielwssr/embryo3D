"""Planche de contrôle 2D des tubes morphables (embryons_3D/tubes_morph.npz + .json), sans Blender.
Pour chaque valeur de « stage » demandée (entière = stade, fractionnaire = interpolation linéaire entre les deux stades
voisins, comme les shape keys), une vignette avec la vue de profil (colonnes = ventral → dorsal, lignes = haut → bas) et
la vue de face (colonnes = gauche → droite), toutes à la même échelle (cadre commun à tous les stades, barre de 1 mm).
Tubes : anneaux + génératrices (pharynx rouge sombre, digestif jaune/orangé, aortes rouges) ; mésos : nappes remplies
translucides (roses), lignes réduites (largeur nulle) invisibles. Les structures absentes au stade sont listées sous le titre.
usage : python embryo3d/tubes_morph_planche.py sortie.png [0 0.5 1 2 ... 6] [--px-mm 30] [--npz embryons_3D/tubes_morph.npz]"""
import numpy as np, json, os, sys, cv2

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
args = sys.argv[1:]
def opt(nom, defaut, conv=float):
    if nom in args:
        i = args.index(nom); v = conv(args[i + 1]); del args[i:i + 2]; return v
    return defaut
px_mm = opt('--px-mm', 30.0); npz = opt('--npz', os.path.join(ROOT, 'embryons_3D', 'tubes_morph.npz'), str)
out = args[0]; valeurs = [float(v) for v in args[1:]] or [0, 0.5, 1, 2, 3, 3.5, 4, 5, 5.5, 6]
Z = np.load(npz); meta = json.load(open(os.path.splitext(npz)[0] + '.json', encoding='utf-8'))
ST = meta['stades']; structs = meta.get('structures', {})
def info(nom): return structs.get(nom, {'type': 'tube', 'famille': 'aorte' if nom.startswith('aorte') else 'digestif'})
COUL = {'pharynx': (60, 60, 200), 'oesophage': (180, 100, 215), 'estomac': (50, 150, 240), 'duodenum': (90, 190, 80),
        'intestin_moyen': (60, 215, 240), 'intestin_posterieur': (240, 140, 80),
        'aorte_dorsale_gauche': (25, 25, 230), 'aorte_dorsale_droite': (50, 90, 240), 'aorte_commune': (25, 15, 190)}
COUL_MESO = (200, 170, 245)

def forme(nom, v):
    """sommets (N, M|K, 3) à la valeur de stage v (interpolation linéaire entre stades voisins)"""
    A = Z[nom]; i = int(np.clip(np.floor(v), 0, len(ST) - 1)); f = float(v) - i
    if f < 1e-9 or i >= len(ST) - 1: return A[min(i, len(ST) - 1)]
    return (1 - f) * A[i] + f * A[i + 1]

# cadre commun (mm) sur tous les stades et toutes les structures
tout = np.concatenate([Z[n].reshape(-1, 3) for n in Z.files])
lo, hi = tout.min(0) - 0.5, tout.max(0) + 0.5
W_prof = int((hi[1] - lo[1]) * px_mm) + 1; W_face = int((hi[0] - lo[0]) * px_mm) + 1; H = int((hi[2] - lo[2]) * px_mm) + 1
def proj(P, vue):
    """(…, 3) mm -> (…, 2) pixels ; profil : colonnes = Y (ventral → dorsal), face : colonnes = X (gauche → droite) ; lignes = Z décroissant"""
    c = (P[..., 1] - lo[1]) if vue == 'profil' else (P[..., 0] - lo[0])
    return np.stack([c * px_mm, (hi[2] - P[..., 2]) * px_mm], -1)

def vignette(v):
    tuiles = []
    for vue, W in (('profil', W_prof), ('face', W_face)):
        im = np.full((H, W, 3), 255, np.uint8); ov = im.copy(); nappe_dessinee = False
        for nom in Z.files:                                     # nappes d'abord (dessous), remplies en translucide
            if info(nom)['type'] != 'nappe': continue
            P = proj(forme(nom, v), vue)
            for i in range(len(P) - 1):
                q = np.array([P[i, 0], P[i, -1], P[i + 1, -1], P[i + 1, 0]], np.int32)
                if cv2.contourArea(q.astype(np.float32)) < 1.0: continue
                cv2.fillPoly(ov, [q], COUL_MESO); nappe_dessinee = True
                cv2.line(ov, tuple(P[i, -1].astype(int)), tuple(P[i + 1, -1].astype(int)), (150, 90, 210), 1)
        if nappe_dessinee: im = cv2.addWeighted(ov, 0.55, im, 0.45, 0)
        for nom in Z.files:
            if info(nom)['type'] != 'tube': continue
            A = forme(nom, v); P = proj(A, vue).astype(np.int32); c = COUL.get(nom, (120, 120, 120))
            if np.ptp(A.reshape(-1, 3), axis=0).max() < 1e-6:
                cv2.circle(im, tuple(P[0, 0]), 2, c, -1); continue     # réduit à un point (absent au stade)
            Nn, Mn = P.shape[:2]
            for i in range(0, Nn, 2): cv2.polylines(im, [np.ascontiguousarray(P[i])], True, c, 1, cv2.LINE_AA)
            for j in range(0, Mn, max(1, Mn // 4)): cv2.polylines(im, [np.ascontiguousarray(P[:, j])], False, c, 1, cv2.LINE_AA)
        cv2.rectangle(im, (0, 0), (W - 1, H - 1), (90, 90, 90), 1)
        cv2.putText(im, vue, (6, H - 8), 0, 0.45, (60, 60, 60), 1, cv2.LINE_AA)
        tuiles.append(im)
    im = np.hstack([tuiles[0], np.full((H, 4, 3), 90, np.uint8), tuiles[1]])
    b = int(px_mm); cv2.line(im, (10, 14), (10 + b, 14), (0, 0, 0), 2); cv2.putText(im, '1 mm', (14 + b, 19), 0, 0.4, (0, 0, 0), 1, cv2.LINE_AA)
    i0 = int(np.floor(v)); titre = f'stage {v:g}' + (f' = {ST[i0]}' if abs(v - i0) < 1e-9 and i0 < len(ST) else f' ({ST[i0]} -> {ST[min(i0 + 1, len(ST) - 1)]})')
    absents = [n for n in Z.files if abs(v - i0) < 1e-9 and i0 < len(ST) and not meta['presence'][n][ST[i0]]]
    bandeau = np.full((44, im.shape[1], 3), 235, np.uint8)
    cv2.putText(bandeau, titre, (8, 18), 0, 0.55, (20, 20, 20), 1, cv2.LINE_AA)
    if absents: cv2.putText(bandeau, 'absents : ' + ', '.join(absents), (8, 38), 0, 0.36, (120, 40, 40), 1, cv2.LINE_AA)
    return np.vstack([bandeau, im])

vs = [vignette(v) for v in valeurs]
cols = min(5, len(vs)); H0, W0 = vs[0].shape[:2]
vs = [cv2.copyMakeBorder(x, 0, H0 - x.shape[0], 0, W0 - x.shape[1] + 6, cv2.BORDER_CONSTANT, value=(255, 255, 255)) for x in vs]
while len(vs) % cols: vs.append(np.full_like(vs[0], 255))
grille = np.vstack([np.hstack(vs[k:k + cols]) for k in range(0, len(vs), cols)])
os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True); cv2.imwrite(out, grille); print('PLANCHE', out, grille.shape)
