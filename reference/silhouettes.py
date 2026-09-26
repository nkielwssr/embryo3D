"""
silhouettes.py — exploite les captures de capture_voka.py (référence externe VOKA) :
  1. silhouettes (masque modèle / fond uniforme) pour chaque vue ;
  2. mesures de proportions sur la vue latérale (az 0, el 0) : plus grande longueur (≈ CRL / greatest
     length), hauteur/largeur, part de la tête ; conversion en mm avec le CRL de la littérature
     (voka_stages.json) -> mm/px ;
  3. enveloppe approchée par « visual hull » (intersection des cônes de silhouettes, projection
     orthographique) -> captures/Jxx/Jxx_hull.ply, à usage de référence privée seulement
     (gabarit de proportions pour les stades dont nous n'avons pas les coupes, jamais un livrable).

Usage : python silhouettes.py [J28 J56 ...] [--grid 160] [--no-hull]
Sorties : captures/Jxx/silhouettes/*.png, captures/Jxx/mesures.json, reference/proportions.json (tous stades)
"""
import argparse, json, os, math
import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
CAP = os.path.join(HERE, 'captures')
STAGES = json.load(open(os.path.join(HERE, 'voka_stages.json'), encoding='utf-8'))['stades']
BY_ID = {s['id']: s for s in STAGES}


def mask_of(path, tol=18):
    a = np.asarray(Image.open(path).convert('RGB')).astype(np.int16)
    border = np.concatenate([a[0], a[-1], a[:, 0], a[:, -1]])
    bg = np.median(border, axis=0)
    m = np.abs(a - bg).max(axis=2) > tol
    # plus grande composante connexe (élimine icônes résiduelles)
    from scipy import ndimage
    lab, n = ndimage.label(m)
    if n > 1:
        sizes = ndimage.sum(m, lab, range(1, n + 1))
        m = lab == (1 + int(np.argmax(sizes)))
    return ndimage.binary_fill_holes(m)


def greatest_length(mask):
    """Diamètre de Feret max du contour (px) et direction."""
    ys, xs = np.nonzero(mask)
    if len(xs) < 10:
        return 0, 0
    pts = np.stack([xs, ys], 1).astype(float)
    # enveloppe convexe
    from scipy.spatial import ConvexHull
    h = ConvexHull(pts)
    p = pts[h.vertices]
    d = ((p[:, None, :] - p[None, :, :]) ** 2).sum(-1)
    i, j = np.unravel_index(np.argmax(d), d.shape)
    v = p[j] - p[i]
    return float(math.sqrt(d[i, j])), float(math.degrees(math.atan2(v[1], v[0])))


def measures(sid, d):
    meta = json.load(open(os.path.join(d, 'meta.json'), encoding='utf-8'))
    st = BY_ID[sid]
    out = dict(stage=sid, carnegie=st['carnegie'], crl_mm_litterature=st['crl_mm'], vues={})
    sd = os.path.join(d, 'silhouettes')
    os.makedirs(sd, exist_ok=True)
    lat = None
    for v in meta['views']:
        m = mask_of(os.path.join(d, v['file']))
        Image.fromarray((m * 255).astype(np.uint8)).save(os.path.join(sd, v['file']))
        ys, xs = np.nonzero(m)
        gl, ang = greatest_length(m)
        rec = dict(az=v['az'], el=v['el'], aire_px=int(m.sum()),
                   bbox=[int(xs.min()), int(ys.min()), int(xs.max()), int(ys.max())] if len(xs) else None,
                   greatest_length_px=gl, gl_angle_deg=ang)
        out['vues'][v['file']] = rec
        if v['az'] == 0 and v['el'] == 0:
            lat = rec
    if lat and lat['greatest_length_px'] > 0:
        mm_px = st['crl_mm'] / lat['greatest_length_px']
        x0, y0, x1, y1 = lat['bbox']
        out['mm_par_px'] = mm_px
        out['lateral_mm'] = dict(greatest_length=st['crl_mm'], largeur_bbox=(x1 - x0) * mm_px,
                                 hauteur_bbox=(y1 - y0) * mm_px,
                                 ratio_hauteur_largeur=(y1 - y0) / max(1, x1 - x0))
        # part de la tête : fraction de l'aire au-dessus du 1er creux (nuque) le long de l'axe GL —
        # approximation : moitié supérieure de la bbox selon l'axe principal
        m = mask_of(os.path.join(d, [k for k in out['vues'] if 'az000_el+00' in k][0]))
        ys, xs = np.nonzero(m)
        c = np.cov(np.stack([xs, ys]))
        w, vec = np.linalg.eigh(c)
        ax = vec[:, 1]
        proj = (xs - xs.mean()) * ax[0] + (ys - ys.mean()) * ax[1]
        hist, edges = np.histogram(proj, 40)
        out['profil_largeur_axe_principal'] = hist.tolist()
    json.dump(out, open(os.path.join(d, 'mesures.json'), 'w', encoding='utf-8'), indent=1)
    return out


def visual_hull(sid, d, N=160):
    """Intersection des silhouettes (projection orthographique, caméra à azimut az autour de Z, élévation el).
    Repère : X gauche→droite de l'image latérale, Z vers le haut de l'image, Y profondeur ; centré."""
    meta = json.load(open(os.path.join(d, 'meta.json'), encoding='utf-8'))
    views = meta['views']
    masks = {v['file']: mask_of(os.path.join(d, v['file'])) for v in views}
    # normalisation : on suppose que Unity garde le même zoom pour toutes les vues d'un stade (vue réinitialisée)
    H, W = next(iter(masks.values())).shape
    # pivot de rotation : Unity tourne autour du pivot du modèle, dont la position écran est fixe.
    # Estimation : moyenne des centres de boîte de toutes les vues à el=0 (x) et de toutes les vues (y).
    lat = [v for v in views if v['az'] == 0 and v['el'] == 0][0]
    from scipy import ndimage
    cxs, cys = [], []
    for v in views:
        ys, xs = np.nonzero(masks[v['file']])
        if v['el'] == 0:
            cxs.append((xs.min() + xs.max()) / 2)
        cys.append((ys.min() + ys.max()) / 2)
    cx, cy = float(np.mean(cxs)), float(np.mean(cys))
    ys, xs = np.nonzero(masks[lat['file']])
    half = max(xs.max() - xs.min(), ys.max() - ys.min()) * 0.6
    # tolérance d'alignement : dilatation de 3 px des silhouettes (évite le sur-creusement)
    masks = {k: ndimage.binary_dilation(m, iterations=3) for k, m in masks.items()}
    def grid(n):
        lin = np.linspace(-1, 1, n)
        X, Y, Z = np.meshgrid(lin, lin, lin, indexing='ij')
        return np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)

    def carve(P, sign_az, sign_el, only_el0=False, k_az=1.0, k_el=1.0, dy=0.0):
        inside = np.ones(len(P), bool)
        for v in views:
            if only_el0 and v['el'] != 0:
                continue
            az, el = sign_az * k_az * math.radians(v['az']), sign_el * k_el * math.radians(v['el'])
            ca, sa = math.cos(az), math.sin(az)
            x = ca * P[:, 0] - sa * P[:, 1]
            y = sa * P[:, 0] + ca * P[:, 1]
            z = P[:, 2]
            ce, se = math.cos(el), math.sin(el)
            z2 = se * y + ce * z
            u = np.round(cx + x * half).astype(int)
            w = np.round(cy + dy - z2 * half).astype(int)
            ok = (u >= 0) & (u < W) & (w >= 0) & (w < H)
            m = masks[v['file']]
            hit = np.zeros(len(P), bool)
            hit[ok] = m[w[ok], u[ok]]
            inside &= hit
        return inside

    # Le sens de rotation, l'échelle exacte des angles et le pivot vertical ne sont pas connus a priori :
    # une hypothèse fausse met les silhouettes en désaccord et sur-creuse l'intersection. On garde donc
    # la combinaison qui conserve le plus de voxels (recherche sur grille grossière, carve final à pleine résolution).
    Pc = grid(max(48, N // 3))
    best = None
    for sa_, se_ in [(1, 1), (-1, 1)]:
        n_in = int(carve(Pc, sa_, se_).sum())
        if best is None or n_in > best[0]:
            best = (n_in, sa_, se_, 1.0, 1.0, 0.0)
    sa_, se_ = best[1], best[2]
    for k_az in (0.9, 1.0, 1.1):
        for k_el in (0.75, 1.0, 1.25, 1.5, 1.75):
            for dy in (-80.0, -50.0, -25.0, 0.0, 25.0):
                n_in = int(carve(Pc, sa_, se_, k_az=k_az, k_el=k_el, dy=dy).sum())
                if n_in > best[0]:
                    best = (n_in, sa_, se_, k_az, k_el, dy)
    n_in, sa_, se_, k_az, k_el, dy = best
    P = grid(N)
    inside = carve(P, sa_, se_, k_az=k_az, k_el=k_el, dy=dy)
    n_el0 = int(carve(P, sa_, se_, True, k_az=k_az).sum())
    n_in = int(inside.sum())
    print(f'  hull : az{sa_:+d} x{k_az} el{se_:+d} x{k_el} dy={dy:+.0f}px -> {n_in} voxels (az seul : {n_el0})')
    hull_params = dict(sign_az=sa_, sign_el=se_, k_az=k_az, k_el=k_el, pivot_dy_px=dy, voxels=n_in, voxels_az_seul=n_el0)
    json.dump(hull_params, open(os.path.join(d, f'{sid}_hull_params.json'), 'w'), indent=1)
    vol = inside.reshape(N, N, N)
    from skimage import measure
    vol = ndimage.binary_opening(vol, iterations=1)
    # aperçu : 3 projections (latérale = -X, face = +Y, dessus = +Z)
    prev = np.concatenate([np.flipud(vol.max(axis=0).T), np.flipud(vol.max(axis=1).T), vol.max(axis=2).T], axis=1)
    Image.fromarray((prev * 255).astype(np.uint8)).resize((prev.shape[1] * 2, prev.shape[0] * 2)).save(os.path.join(d, f'{sid}_hull_apercu.png'))
    if vol.sum() == 0:
        print('  hull vide (vérifier calibrage/orientation)')
        return
    # lissage du volume binaire (marche d'escalier des voxels) avant extraction de la surface
    smooth = ndimage.gaussian_filter(vol.astype(np.float32), sigma=1.6)
    verts, faces, _, _ = measure.marching_cubes(smooth, 0.5, spacing=(2 / N,) * 3)
    verts -= 1.0
    mm_px = BY_ID[sid]['crl_mm'] / max(1e-6, greatest_length(masks[lat['file']])[0])
    verts *= half * mm_px   # en mm
    # repère du viewer/pipeline : X gauche-droite, Y ventral→dorsal, Z bas→haut. La caméra « az 0 » du viewer
    # est en -X (écran-droite = -Y) ; l'horizontale de l'image latérale VOKA (X_hull) doit donc devenir -Y.
    import trimesh
    mesh = trimesh.Trimesh(np.stack([verts[:, 1], -verts[:, 0], verts[:, 2]], 1), faces, process=True)
    mesh.export(os.path.join(d, f'{sid}_hull.ply'))
    print(f'  {sid}_hull.ply : {len(verts)} sommets, extension mm {mesh.extents.round(1)}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('stages', nargs='*')
    ap.add_argument('--grid', type=int, default=160)
    ap.add_argument('--no-hull', action='store_true')
    a = ap.parse_args()
    allm = {}
    for sid in sorted(os.listdir(CAP)) if os.path.isdir(CAP) else []:
        if a.stages and sid not in a.stages:
            continue
        d = os.path.join(CAP, sid)
        if not os.path.exists(os.path.join(d, 'meta.json')):
            continue
        print('==', sid)
        allm[sid] = measures(sid, d)
        if not a.no_hull:
            visual_hull(sid, d, a.grid)
    json.dump(allm, open(os.path.join(HERE, 'proportions.json'), 'w', encoding='utf-8'), indent=1)
    print('-> proportions.json')


if __name__ == '__main__':
    main()
