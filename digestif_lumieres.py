"""Détection automatique des lumières tubulaires (voxel pâle entouré d'épithélium dense dans presque toutes les
directions), pour repérer tube digestif, bronches, tubules mésonéphriques, etc. Les composantes sont ensuite
choisies à la main par stade (digestif_choix).
usage : python digestif_lumieres.py <work> [R=10] [seuil_encl=0.8]
Sorties : <work>/digestif/lumieres_lab.npy (labels à résolution /2), lumieres.json (taille, centre, boîte),
          lumieres_planche.png (projections face et profil, numéros des composantes)."""
import numpy as np, json, os, sys, cv2
from scipy import ndimage as ndi

work = sys.argv[1]; R = int(sys.argv[2]) if len(sys.argv) > 2 else 10; TH = float(sys.argv[3]) if len(sys.argv) > 3 else 0.8
od = os.path.join(work, 'digestif'); os.makedirs(od, exist_ok=True)
d = np.load(os.path.join(work, 'dens.npy'), mmap_mode='r')
D = ndi.zoom(np.asarray(d).astype(np.float32), 0.5, order=1)            # résolution /2
D = ndi.gaussian_filter(D, 0.6)
L = np.load(os.path.join(work, 'labels.npz')); SH = tuple(L['shape'])
def lab(k):
    if k not in L.files: return np.zeros(D.shape, bool)
    m = np.unpackbits(L[k])[:int(np.prod(SH))].reshape(SH)
    return ndi.zoom(m, 0.5, order=0)[:D.shape[0], :D.shape[1], :D.shape[2]].astype(bool)
env = lab('enveloppe')
excl = lab('ventricules') | lab('cavite_pericardique') | lab('yeux') | lab('vesicules_otiques')
pale = (D < 22) & ~excl
dense = D > 50
Rh = max(3, R // 2)
dirs = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1) if (a, b, c) != (0, 0, 0)]
count = np.zeros(D.shape, np.uint8)
for (a, b, c) in dirs:
    hit = np.zeros(D.shape, bool)
    for k in range(1, Rh + 1):
        sh = np.roll(dense, (-a * k, -b * k, -c * k), axis=(0, 1, 2))
        hit |= sh
    count += hit
encl = count / len(dirs)
lum = pale & (encl >= TH) & ndi.binary_erosion(env, iterations=2)
lum = ndi.binary_opening(lum, iterations=1)
labs, n = ndi.label(lum)
sizes = ndi.sum(lum, labs, range(1, n + 1))
keep = [i + 1 for i, s in enumerate(sizes) if 20 <= s <= 400000]
labs = np.where(np.isin(labs, keep), labs, 0)
info = []
objs = ndi.find_objects(labs)
for i in keep:
    sl = objs[i - 1]
    if sl is None: continue
    m = labs[sl] == i; c = np.array(ndi.center_of_mass(m)) + [s.start for s in sl]
    info.append({'id': int(i), 'taille': int(m.sum()), 'centre_x_s_y_pleine_res': (c * 2).round(0).astype(int).tolist(),
                 'boite': [[s.start * 2, s.stop * 2] for s in sl]})
info.sort(key=lambda r: -r['taille'])
np.save(os.path.join(od, 'lumieres_lab.npy'), labs.astype(np.int32))
json.dump(info, open(os.path.join(od, 'lumieres.json'), 'w'), indent=0)
print('composantes', len(info), 'plus grosses', [(r['id'], r['taille']) for r in info[:15]])
# planche : projection densité + composantes colorées et numérotées (les 60 plus grosses)
rng = np.random.default_rng(1)
def planche(ax, flipT):
    base = 255 - D.max(axis=ax); im = cv2.cvtColor(base.astype(np.uint8), cv2.COLOR_GRAY2BGR)
    if flipT: im = im.transpose(1, 0, 2).copy()
    im = cv2.resize(im, None, fx=2, fy=2, interpolation=cv2.INTER_NEAREST)
    for r in info[:60]:
        m = (labs == r['id']).any(axis=ax)
        if flipT: m = m.T
        m = cv2.resize(m.astype(np.uint8), (im.shape[1], im.shape[0]), interpolation=cv2.INTER_NEAREST).astype(bool)
        col = tuple(int(v) for v in rng.integers(40, 230, 3)); im[m] = col
        ys, xs = np.nonzero(m); cv2.putText(im, str(r['id']), (int(xs.mean()), int(ys.mean())), 0, 0.4, (0, 0, 0), 1)
    return im
a = planche(2, True)     # face : lignes haut-bas, colonnes gauche-droite
b = planche(0, False)    # profil : lignes haut-bas, colonnes ventral-dorsal
h = max(a.shape[0], b.shape[0])
cv2.imwrite(os.path.join(od, 'lumieres_planche.png'), np.hstack([a, np.full((h, 6, 3), 90, np.uint8), b]))
print('OK')
