"""Cœur par stade : tissu contenu dans le sac péricardique rempli, puis cavités (lumières pâles internes).
usage : python cardio_coeur.py <work> <params.json>
params : {"seuil_tissu": 12, "graines": [[x,s,y],...], "boite": [[x0,x1],[s0,s1],[y0,y1]] (optionnel),
          "fermeture_sac": 6, "ouverture": 1, "sac": "label" | "boite"}
  sac = "label" : sac = cavite_pericardique du pipeline, fermé puis rempli ; "boite" : pas de sac, seulement la boîte
  (le cœur est alors la composante connexe des graines, coupée des ponts fins par ouverture).
Sorties : <work>/cardio/coeur.npz (packbits : coeur_plein, myocarde, cavites_cardiaques), controle_coeur.png"""
import numpy as np, json, os, sys, cv2
from scipy import ndimage as ndi

work, pfile = sys.argv[1], sys.argv[2]
P = json.load(open(pfile, encoding='utf-8'))['coeur']
dens = np.load(os.path.join(work, 'dens.npy'), mmap_mode='r'); SH = dens.shape
od = os.path.join(work, 'cardio'); os.makedirs(od, exist_ok=True)
L = np.load(os.path.join(work, 'labels.npz')); n = int(np.prod(tuple(L['shape'])))
def lab(k): return np.unpackbits(L[k])[:n].reshape(SH).astype(bool) if k in L.files else np.zeros(SH, bool)
box = P.get('boite')
if box is None:
    peri = lab('cavite_pericardique'); idx = [np.where(peri.any(axis=tuple(j for j in range(3) if j != a)))[0] for a in range(3)]
    box = [[max(int(i[0]) - 15, 0), min(int(i[-1]) + 15, SH[a])] for a, i in enumerate(idx)]
lo = np.array([b[0] for b in box]); hi = np.array([b[1] for b in box])
sl = tuple(slice(lo[a], hi[a]) for a in range(3))
D = ndi.gaussian_filter(np.asarray(dens[sl]).astype(np.float32), 1.0)
tissu = D > P.get('seuil_tissu', 12)
ball = lambda r: np.linalg.norm(np.indices((2 * r + 1,) * 3) - r, axis=0) <= r
if P.get('sac', 'label') == 'label':
    peri = lab('cavite_pericardique')[sl]
    r = P.get('fermeture_sac', 6)
    sac = ndi.binary_fill_holes(ndi.binary_closing(np.pad(peri, r), ball(r))[r:-r, r:-r, r:-r])
    sac = ndi.binary_fill_holes(sac) & ~peri
    cand = tissu & sac
elif P.get('sac') in ('convexe', 'pericarde_auto'):   # enveloppe convexe de la cavité péricardique, coupe par coupe
    if P['sac'] == 'convexe':
        peri = lab('cavite_pericardique')[sl]
    else:                                  # cavité péricardique recalculée : plus grande lumière pâle interne de la boîte
        env = ndi.binary_erosion(lab('enveloppe')[sl], iterations=4)
        pale = (D < P.get('seuil_pale', 6)) & env & ~lab('ventricules')[sl] & ~ndi.binary_dilation(lab('snc')[sl], iterations=2)
        pale = ndi.binary_opening(pale, ball(2)); lp, npl = ndi.label(pale)
        szp = ndi.sum(pale, lp, range(1, npl + 1)); peri = lp == (int(np.argmax(szp)) + 1)
        print('pericarde auto', int(peri.sum()))
    hull = np.zeros_like(peri)
    for k in range(peri.shape[1]):
        pts = np.argwhere(peri[:, k, :])
        if len(pts) < 20: continue
        h = cv2.convexHull(pts[:, ::-1].astype(np.int32)); img = np.zeros((peri.shape[2], peri.shape[0]), np.uint8)
        cv2.fillPoly(img, [h], 1); hull[:, k, :] = img.T.astype(bool)
    hull = ndi.binary_erosion(hull, iterations=P.get('erosion_sac', 2))
    cand = tissu & hull & ~peri
elif P.get('sac') == 'pipeline':          # enveloppe convexe du cœur du pipeline + cavité péricardique
    peri = lab('cavite_pericardique')[sl]
    cand = tissu & (lab('coeur')[sl] | peri) & ~peri
else:
    cand = tissu
for k in P.get('exclure', ['foie']):
    cand &= ~ndi.binary_dilation(lab(k)[sl], iterations=2)
if P.get('ouverture', 1): cand = ndi.binary_opening(cand, ball(P.get('ouverture', 1)))
lb, nlab = ndi.label(cand)
if 'graines' not in P:                        # défaut : tissu dense le plus central du cœur du pipeline
    ch = lab('coeur')[sl] & tissu
    dt = ndi.distance_transform_edt(ch); top = np.argwhere(dt >= np.percentile(dt[ch], 97))
    P['graines'] = (top[np.linspace(0, len(top) - 1, 6).astype(int)] + lo).tolist()
    print('graines auto', P['graines'])
seeds = [tuple(np.array(g) - lo) for g in P['graines']]
ids = set()
if P.get('plus_grande', False):
    szl = ndi.sum(cand, lb, range(1, nlab + 1)); ids.add(int(np.argmax(szl)) + 1)
for g in ([] if P.get('plus_grande', False) else seeds):
    g = tuple(int(v) for v in g)
    if lb[g] > 0: ids.add(lb[g])
    else:                                                    # graine hors tissu : composante la plus proche
        sub = lb[max(g[0]-4,0):g[0]+5, max(g[1]-4,0):g[1]+5, max(g[2]-4,0):g[2]+5]; v = sub[sub > 0]
        if v.size: ids.add(int(np.bincount(v).argmax()))
myo = np.isin(lb, list(ids))
if P.get('sac') == 'label_corrige':           # cœur corrigé à la main (3D Slicer) : pris tel quel
    myo = lab('coeur')[sl]; print('coeur corrige (Slicer)', int(myo.sum()))
if P.get('sac') == 'bassin_densite':          # partage des eaux sur la densité : coupure dans les espaces pâles
    from skimage.segmentation import watershed
    el = -ndi.gaussian_filter(D, P.get('lissage', 2.0))
    mk = np.zeros(D.shape, np.int32)
    for g in seeds:
        g = tuple(int(v) for v in g); r = P.get('r_graine', 5)
        mk[max(g[0]-r,0):g[0]+r, max(g[1]-r,0):g[1]+r, max(g[2]-r,0):g[2]+r] = 1
    faces = P.get('faces', ['x0', 'x1', 's0', 's1', 'y0', 'y1'])
    fm = {'x0': (0, 0), 'x1': (0, -1), 's0': (1, 0), 's1': (1, -1), 'y0': (2, 0), 'y1': (2, -1)}
    for f in faces:
        a, i = fm[f]; idx = [slice(None)] * 3; idx[a] = i; sub = mk[tuple(idx)]; sub[sub == 0] = 2
    for k in P.get('exclure', []): mk[lab(k)[sl] & (mk == 0)] = 2
    ws = watershed(el, mk)
    myo = (ws == 1) & tissu
    lb2, n2 = ndi.label(myo); 
    if n2:
        sz2 = ndi.sum(myo, lb2, range(1, n2 + 1)); myo = np.isin(lb2, 1 + np.where(sz2 >= 0.05 * sz2.max())[0])
    print('bassin densite : coeur', int(myo.sum()))
if P.get('sac') == 'bassin':                  # ligne de partage des eaux : cœur vs reste, coupure aux cols étroits
    from skimage.segmentation import watershed
    t = ndi.binary_opening(tissu, ball(1))
    for k in P.get('exclure', []): t &= ~lab(k)[sl]
    dt = ndi.distance_transform_edt(t)
    mk = np.zeros(t.shape, np.int32)
    for g in seeds:
        g = tuple(int(v) for v in g); r = 4
        mk[max(g[0]-r,0):g[0]+r, max(g[1]-r,0):g[1]+r, max(g[2]-r,0):g[2]+r] = 1
    brd = np.zeros(t.shape, bool); brd[[0, -1], :, :] = True; brd[:, [0, -1], :] = True; brd[:, :, [0, -1]] = True
    for g in P.get('graines_exterieur', []):
        g = tuple(int(v) for v in np.array(g) - lo); brd[max(g[0]-4,0):g[0]+4, max(g[1]-4,0):g[1]+4, max(g[2]-4,0):g[2]+4] = True
    for k in P.get('exterieur', ['foie', 'snc', 'enveloppe_bord']):
        if k == 'enveloppe_bord':
            env = lab('enveloppe')[sl]; brd |= env & ~ndi.binary_erosion(env, iterations=P.get('paroi', 8))
        else: brd |= lab(k)[sl]
    mk[brd & t & (mk == 0)] = 2
    ws = watershed(-ndi.gaussian_filter(dt, 1.0), mk, mask=t)
    myo = ws == 1
    print('bassin : coeur', int(myo.sum()))
rb = P.get('ouverture', 1)
if rb >= 2 and P.get('sac') not in ('bassin', 'bassin_densite', 'label_corrige'):                                   # rendre l'épaisseur retirée par l'ouverture, sans refranchir les ponts
    base = tissu.copy()
    for k in P.get('exclure', ['foie']): base &= ~ndi.binary_dilation(lab(k)[sl], iterations=2)
    myo = ndi.binary_dilation(myo, ball(rb)) & base
plein = ndi.binary_fill_holes(ndi.binary_closing(myo, ball(2)))
# cavités : lumières pâles fermées dans le cœur (remplissage coupe par coupe dans les 3 directions, intersection)
cl = ndi.binary_closing(myo, ball(3))
f = np.ones_like(cl)
for ax in range(3):
    ff = np.stack([ndi.binary_fill_holes(s) for s in np.moveaxis(cl, ax, 0)], 0); f &= np.moveaxis(ff, 0, ax)
plein |= f
cav = plein & ~myo & (D < P.get('seuil_tissu', 12) + 8)
cav = ndi.binary_opening(cav, ball(1))
lc, nc = ndi.label(cav); sz = ndi.sum(cav, lc, range(1, nc + 1)); cav = np.isin(lc, 1 + np.where(sz >= P.get('cavite_min', 200))[0])
def full(m):
    F = np.zeros(SH, bool); F[sl] = m; return F
np.savez_compressed(os.path.join(od, 'coeur.npz'), coeur_plein=np.packbits(full(plein)), myocarde=np.packbits(full(myo)),
                    cavites_cardiaques=np.packbits(full(cav)), shape=np.array(SH))
print(f'coeur plein {plein.sum()} vox, myocarde {myo.sum()}, cavites {cav.sum()} ({len(np.unique(lc[cav]))} composantes)')
# contrôle : 8 coupes transversales dans la boîte
tiles = []
for s in np.linspace(0, D.shape[1] - 1, 10).astype(int)[1:-1]:
    im = cv2.cvtColor((255 - np.clip(D[:, s, :].T, 0, 255)).astype(np.uint8), cv2.COLOR_GRAY2BGR)
    for m, col in ((plein, (0, 140, 255)), (cav, (255, 0, 0))):
        c, _ = cv2.findContours(np.ascontiguousarray(m[:, s, :].T.astype(np.uint8)), cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)
        cv2.drawContours(im, c, -1, col, 1)
    cv2.putText(im, str(s + lo[1]), (5, 15), 0, 0.5, (0, 0, 200), 1); tiles.append(im)
cv2.imwrite(os.path.join(od, 'controle_coeur.png'), np.vstack([np.hstack(tiles[:4]), np.hstack(tiles[4:8])]))
