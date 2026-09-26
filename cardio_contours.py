"""Cœur détouré à la main : polygones sur quelques coupes transversales clés (s), interpolés entre elles par
distance signée, puis tissu / cavités dans le volume détouré.
usage : python cardio_contours.py <work> <params.json>
params["contours"] = {"<s>": [[x, y], ...], ...}  (x = gauche-droite, y = ventral-dorsal, coordonnées de dens.npy)
Au-delà de la première et de la dernière coupe clé, le contour se referme sur "fermeture" coupes (défaut 4).
Sortie : <work>/cardio/coeur.npz (coeur_plein, myocarde, cavites_cardiaques), controle_coeur.png"""
import numpy as np, json, os, sys, cv2
from scipy import ndimage as ndi
work, pfile = sys.argv[1], sys.argv[2]
P = json.load(open(pfile, encoding='utf-8'))['coeur']
dens = np.load(os.path.join(work, 'dens.npy'), mmap_mode='r'); SH = dens.shape
od = os.path.join(work, 'cardio'); os.makedirs(od, exist_ok=True)
C = {int(k): np.array(v, np.int32) for k, v in P['contours'].items()}
ks = sorted(C)
def sdf(poly):
    m = np.zeros((SH[2], SH[0]), np.uint8); cv2.fillPoly(m, [poly[:, ::-1][:, ::-1]], 1)   # image (y lignes, x colonnes)
    m = m.astype(bool); return ndi.distance_transform_edt(~m) - ndi.distance_transform_edt(m)
F = {k: sdf(C[k]) for k in ks}
fe = P.get('fermeture', 4)
plein = np.zeros(SH, bool)
for s in range(ks[0] - fe, ks[-1] + fe + 1):
    if s < 0 or s >= SH[1]: continue
    if s < ks[0]: f = F[ks[0]] + (ks[0] - s) * (np.abs(F[ks[0]].min()) / fe)
    elif s > ks[-1]: f = F[ks[-1]] + (s - ks[-1]) * (np.abs(F[ks[-1]].min()) / fe)
    else:
        i = np.searchsorted(ks, s); 
        if ks[min(i, len(ks)-1)] == s: f = F[s]
        else:
            a, b = ks[i - 1], ks[i]; t = (s - a) / (b - a); f = (1 - t) * F[a] + t * F[b]
    plein[:, s, :] = (f < 0).T
D = np.asarray(dens).astype(np.float32)
thr = P.get('seuil_tissu', 12)
myo = plein & (ndi.gaussian_filter(D, 1.0) > thr)
myo = ndi.binary_opening(myo, iterations=1)
lm, nm = ndi.label(myo); szm = ndi.sum(myo, lm, range(1, nm + 1)); myo = np.isin(lm, 1 + np.where(szm >= 0.02 * szm.max())[0])
cl = ndi.binary_closing(myo, iterations=3)
enc = np.stack([ndi.binary_fill_holes(cl[:, k, :]) for k in range(SH[1])], 1)   # lumières fermées par la paroi
cav = enc & ~myo & plein; cav = ndi.binary_opening(cav, iterations=1)
plein = (myo | cav | (enc & plein))
lc, nc = ndi.label(cav); sz = ndi.sum(cav, lc, range(1, nc + 1)); cav = np.isin(lc, 1 + np.where(sz >= P.get('cavite_min', 150))[0])
np.savez_compressed(os.path.join(od, 'coeur.npz'), coeur_plein=np.packbits(plein), myocarde=np.packbits(myo),
                    cavites_cardiaques=np.packbits(cav), shape=np.array(SH))
print(f'coeur plein {plein.sum()} vox, myocarde {myo.sum()}, cavites {cav.sum()}')
# contrôle : 8 coupes entre les bornes, zoom sur la boîte
xs, ys = np.concatenate([C[k][:, 0] for k in ks]), np.concatenate([C[k][:, 1] for k in ks])
x0, x1, y0, y1 = max(xs.min() - 25, 0), min(xs.max() + 25, SH[0]), max(ys.min() - 25, 0), min(ys.max() + 25, SH[2])
tiles = []
for s in np.linspace(ks[0] - fe + 1, ks[-1] + fe - 1, 8).astype(int):
    im = cv2.cvtColor((255 - np.clip(D[x0:x1, s, y0:y1].T, 0, 255)).astype(np.uint8), cv2.COLOR_GRAY2BGR)
    c, _ = cv2.findContours(np.ascontiguousarray(plein[x0:x1, s, y0:y1].T.astype(np.uint8)), cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)
    cv2.drawContours(im, c, -1, (0, 140, 255), 1)
    if s in C: cv2.polylines(im, [C[s] - [x0, y0]], True, (0, 200, 0), 1)
    cv2.putText(im, str(s), (4, 14), 0, 0.45, (0, 0, 200), 1); tiles.append(im)
cv2.imwrite(os.path.join(od, 'controle_coeur.png'), np.vstack([np.hstack(tiles[:4]), np.hstack(tiles[4:])]))
