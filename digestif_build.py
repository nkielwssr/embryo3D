"""Tube digestif guidé par points de passage (semi-automatique), sans toucher à labels.npz.

usage : python digestif_build.py <work_dir> <points.json> [--no-mesh]

points.json : {"segments": [{"nom": "oesophage", "points": [[x, s, y], ...], "rmax": 14}, ...]}
  coordonnées dans dens.npy : x = gauche-droite (axe 0), s = haut-bas (axe 1), y = ventral-dorsal (axe 2).
  Les points peuvent être approximatifs : chacun est recalé sur le centre de lumière le plus proche
  (voxel pâle entouré d'un anneau d'épithélium dense), sauf si le point porte un 4e élément 0 (point imposé).
Pour chaque segment : spline par les points recalés, puis rayon de paroi mesuré dans le plan perpendiculaire
(16 rayons : lumière pâle → pic d'épithélium → bord externe à mi-hauteur), médiane glissante.
Sorties dans <work>/digestif/ : digestif.npz (masques packbits par segment + union), chemins.json
(centres, rayons, recalage), planche de contrôle controle_digestif.png, PLY par segment si trimesh est présent."""
import numpy as np, json, os, sys
from scipy import ndimage as ndi
from scipy.interpolate import CubicSpline

work, pts_file = sys.argv[1], sys.argv[2]
do_mesh = '--no-mesh' not in sys.argv
dens = np.load(os.path.join(work, 'dens.npy'), mmap_mode='r')
SH = dens.shape
out_dir = os.path.join(work, sys.argv[sys.argv.index('--sortie') + 1] if '--sortie' in sys.argv else 'digestif'); os.makedirs(out_dir, exist_ok=True)
NPZ = cfg_name = None
cfg = json.load(open(pts_file, encoding='utf-8'))

def crop_box(P, pad):
    lo = np.maximum(np.floor(P.min(0)) - pad, 0).astype(int); hi = np.minimum(np.ceil(P.max(0)) + pad + 1, SH).astype(int)
    return lo, hi

def rays_profile(D, c, n1, n2, rlen, ndir=16, step=0.5):
    ang = np.linspace(0, 2 * np.pi, ndir, endpoint=False)
    r = np.arange(0, rlen, step)
    dirs = np.cos(ang)[:, None] * n1[None] + np.sin(ang)[:, None] * n2[None]            # ndir x 3
    pts = c[None, None, :] + dirs[:, None, :] * r[None, :, None]                          # ndir x nr x 3
    v = ndi.map_coordinates(D, pts.reshape(-1, 3).T, order=1, mode='nearest').reshape(len(ang), len(r))
    return r, v

DELTA = 25    # contraste minimal paroi - lumière (relatif : les stades n'ont pas la même densité de fond)

def enclosure(D, c, n1, n2, rlen=14, delta=DELTA):
    r, v = rays_profile(D, c, n1, n2, rlen)
    thr = v[:, 0].mean() + delta
    hit = (v > thr)
    first = np.where(hit.any(1), hit.argmax(1), -1)
    return (first > 0).mean(), np.where(first > 0, r[np.maximum(first, 0)], rlen)

def plane(t):
    t = t / (np.linalg.norm(t) + 1e-9); a = np.array([1., 0, 0]) if abs(t[0]) < 0.9 else np.array([0, 0, 1.])
    n1 = np.cross(t, a); n1 /= np.linalg.norm(n1); n2 = np.cross(t, n1); return n1, n2

def snap(D, p, t, rad=6, pale_max=None, lum_min=-1):
    """centre de lumière le plus proche dans le plan perpendiculaire à t : minimum local de densité
    entouré d'un anneau plus dense (contraste relatif)"""
    n1, n2 = plane(t); best = (-1e9, p, 0)
    for a in np.arange(-rad, rad + 0.1, 1.0):
        for b in np.arange(-rad, rad + 0.1, 1.0):
            if a * a + b * b > rad * rad: continue
            q = p + a * n1 + b * n2
            e, dist = enclosure(D, q, n1, n2)
            v0 = ndi.map_coordinates(D, q[:, None], order=1)[0]
            if v0 < lum_min: continue
            score = e + 0.02 * np.median(dist) - 0.01 * np.hypot(a, b) - 0.002 * v0
            if score > best[0]: best = (score, q, e)
    return best[1], best[2]

def wall_radius(D, c, t, rmax, delta=DELTA):
    n1, n2 = plane(t); r, v = rays_profile(D, c, n1, n2, rmax + 6, ndir=24)
    c0 = v[:, 0].mean(); rads = []
    for prof in v:
        above = prof > c0 + delta
        if not above.any(): continue
        k = int(np.argmax(above))
        if k <= 0: continue
        pk = k + int(np.argmax(prof[k:k + int(8 / 0.5)]))
        edge = prof[pk] - 0.5 * (prof[pk] - c0)
        after = np.where(prof[pk:] < edge)[0]
        if len(after) == 0: continue
        rads.append(r[pk + after[0]])
    if len(rads) < 6: return None
    return float(np.clip(np.median(rads), 1.5, rmax))

def sac(seg):
    """poche (estomac…) : lumière remplie coupe par coupe (transversales) à l'intérieur de l'anneau épithélial
    fermé, propagée depuis les coupes des graines vers le haut et le bas (recouvrement avec la coupe voisine),
    puis paroi = dilatation de l'épaisseur donnée, limitée au tissu"""
    P = np.array([p[:3] for p in seg['points']], int); box = seg['boite']          # [[x0,x1],[s0,s1],[y0,y1]]
    lo = np.array([b[0] for b in box]); hi = np.array([b[1] for b in box])
    D = ndi.gaussian_filter(np.asarray(dens[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]).astype(np.float32), 1.0)
    thr = seg.get('seuil_paroi', 40); minov = seg.get('recouvrement', 0.3); amin = seg.get('aire_min', 30)
    n = D.shape[1]; lum = np.zeros(D.shape, bool)
    def cands(k):
        w = ndi.binary_closing(D[:, k, :] > thr, iterations=2)
        f = ndi.binary_fill_holes(w) & ~w
        f = ndi.binary_opening(f, iterations=1)
        lab, nn = ndi.label(f); return lab, nn
    seeds_by_slice = {}
    for p in P - lo: seeds_by_slice.setdefault(int(p[1]), []).append((int(p[0]), int(p[2])))
    for k0, pts in seeds_by_slice.items():
        lab, nn = cands(k0)
        ids = {lab[x, y] for x, y in pts if lab[x, y] > 0}
        cur = np.isin(lab, list(ids)); lum[:, k0, :] |= cur
        for direction in (1, -1):
            prev = cur; k = k0 + direction
            while 0 <= k < n and prev.sum() >= amin:
                lab, nn = cands(k)
                if nn == 0: break
                ov = ndi.sum(prev, lab, range(1, nn + 1)); sz = ndi.sum(np.ones_like(prev), lab, range(1, nn + 1))
                keep = [i + 1 for i in range(nn) if ov[i] >= minov * min(sz[i], prev.sum()) and sz[i] < 4 * prev.sum() + 400]
                if not keep: break
                cur = np.isin(lab, keep); lum[:, k, :] |= cur; prev = cur; k += direction
    lum = ndi.binary_closing(lum, iterations=2)
    ep = seg.get('epaisseur', 6)
    m = ndi.binary_dilation(lum, iterations=ep)
    m = ndi.binary_fill_holes(ndi.binary_closing(m, iterations=2))
    full = np.zeros(SH, bool); full[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]] = m
    return full, int(lum.sum())

result, masks = {}, {}
for seg in cfg['segments']:
    if seg.get('type') == 'sac':
        masks[seg['nom']], nl = sac(seg)
        result[seg['nom']] = {'type': 'sac', 'voxels_lumiere': nl, 'voxels': int(masks[seg['nom']].sum())}
        print(f"{seg['nom']}: poche, lumière {nl} vox, total {masks[seg['nom']].sum()} vox"); continue
    P = np.array([p[:3] for p in seg['points']], float)
    fixed = [len(p) > 3 and p[3] == 0 for p in seg['points']]
    rmax = seg.get('rmax', 12)
    lo, hi = crop_box(P, rmax + 20)
    D = ndi.gaussian_filter(np.asarray(dens[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]).astype(np.float32), 0.8)
    Pl = P - lo
    # tangentes approximatives puis recalage
    T = np.gradient(Pl, axis=0) if len(Pl) > 1 else np.array([[0, 1., 0]])
    snapped, encl = [], []
    for i, p in enumerate(Pl):
        if fixed[i]: snapped.append(p); encl.append(None); continue
        q, e = snap(D, p, T[i], rad=seg.get('rayon_recalage', 6), lum_min=seg.get('lumiere_min', -1)); snapped.append(q); encl.append(round(float(e), 2))
    S_ = np.array(snapped)
    # spline paramétrée par longueur
    dl = np.r_[0, np.cumsum(np.linalg.norm(np.diff(S_, axis=0), axis=1))]
    keep = np.r_[True, np.diff(dl) > 0.5]; S_, dl = S_[keep], dl[keep]
    if len(S_) >= 2:
        cs = CubicSpline(dl, S_, bc_type='natural'); u = np.arange(0, dl[-1] + 1e-6, 1.0)
        C = cs(u); Tg = cs(u, 1)
    else:
        C, Tg = S_, np.array([[0, 1., 0]])
    R = []
    for c, t in zip(C, Tg):
        R.append(wall_radius(D, c, t, rmax))
    R = np.array([np.nan if r is None else r for r in R], float)
    if np.isnan(R).all(): R[:] = seg.get('r_defaut', 4.0)
    idx = np.arange(len(R)); good = ~np.isnan(R)
    R = np.interp(idx, idx[good], R[good])
    R = ndi.median_filter(R, size=min(9, len(R) | 1), mode='nearest'); R = ndi.uniform_filter1d(R, 5, mode='nearest')
    # voxelisation : union de boules
    m = np.zeros(hi - lo, bool)
    zz = np.indices((2 * int(rmax) + 3,) * 3) - (int(rmax) + 1)
    for c, r in zip(C, R):
        ci = np.round(c).astype(int); rr = int(np.ceil(r)) + 1
        sl = tuple(slice(max(ci[k] - rr, 0), min(ci[k] + rr + 1, m.shape[k])) for k in range(3))
        g = np.ogrid[sl[0], sl[1], sl[2]]
        dist2 = (g[0] - c[0]) ** 2 + (g[1] - c[1]) ** 2 + (g[2] - c[2]) ** 2
        m[sl] |= dist2 <= r * r
    full = np.zeros(SH, bool); full[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]] = m
    masks[seg['nom']] = full
    result[seg['nom']] = {'centres': (C + lo).round(2).tolist(), 'rayons': R.round(2).tolist(),
                          'points_recales': (S_ + lo).round(1).tolist(), 'enclosure': encl,
                          'longueur_vox': float(dl[-1]) if len(dl) else 0, 'rayon_median': float(np.median(R))}
    print(f"{seg['nom']}: {len(C)} pts, longueur {dl[-1]:.0f} vox, rayon médian {np.median(R):.1f}, recalage {encl}")

union = np.zeros(SH, bool)
for m in masks.values(): union |= m
np.savez_compressed(os.path.join(out_dir, cfg.get('fichier', 'digestif.npz')), **{k: np.packbits(v) for k, v in masks.items()},
                    union=np.packbits(union), shape=np.array(SH))
json.dump(result, open(os.path.join(out_dir, cfg.get('fichier', 'digestif.npz').replace('.npz', '_chemins.json') if 'fichier' in cfg else 'chemins.json'), 'w', encoding='utf-8'), ensure_ascii=False)

# planche de contrôle : projections (ventrale, latérale) de la densité avec le tube superposé
import cv2
def proj(ax):
    base = np.asarray(dens).max(axis=ax).astype(np.float32)
    base = 255 - base; im = cv2.cvtColor(base.astype(np.uint8), cv2.COLOR_GRAY2BGR)
    colors = [(255, 0, 255), (0, 160, 255), (0, 200, 0), (255, 120, 0), (0, 0, 255), (200, 200, 0), (120, 0, 200)]
    for k, (nom, m) in enumerate(masks.items()):
        mm = m.any(axis=ax)
        ov = im.copy(); ov[mm] = colors[k % len(colors)]; im = cv2.addWeighted(ov, 0.55, im, 0.45, 0)
    return im
a = proj(2).transpose(1, 0, 2)      # vue de face : lignes = haut-bas, colonnes = gauche-droite
b = proj(0)                          # vue latérale : lignes = haut-bas, colonnes = ventral-dorsal
h = max(a.shape[0], b.shape[0])
cv2.imwrite(os.path.join(out_dir, cfg.get('controle', 'controle_digestif.png')), np.hstack([a, np.full((h, 6, 3), 90, np.uint8), b]))

if do_mesh:
    try:
        from skimage import measure; import trimesh
        meta = json.load(open(os.path.join(work, 'fused_meta.json')))
        for nom, m in masks.items():
            sm = ndi.gaussian_filter(m.astype(np.float32), 1.0)
            v, f, _, _ = measure.marching_cubes(sm, 0.5)
            trimesh.Trimesh(v, f).export(os.path.join(out_dir, f'{nom}.ply'))
        print('PLY ecrits (coordonnees voxel)')
    except Exception as e:
        print('maillage non fait :', e)
print('OK', out_dir)
