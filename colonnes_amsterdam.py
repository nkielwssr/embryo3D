"""
Colonnes d'Amsterdam calées sur nos embryons : un PSD à calques + un PNG composite par stade (3 vues : profil gauche, dos, face).
usage : python embryo3d/colonnes_amsterdam.py <dossier_transformations> [stades...]
Les transformations atlas -> pipeline (T_<CS>_*.npy : échelle, R 3x3, t) viennent des essais de calage du 25/09/2026
(peau + pointages queue/œil de l'utilisateur ; CS13 affiné par la courbure dorsale). Rien n'est écrit dans les modèles.
Sorties : embryons_3D/colonnes_amsterdam/<CS>_colonne.psd, <CS>_colonne.png, calques/<CS>_<calque>.png
Vues orientées sur les axes ANATOMIQUES d'Amsterdam (haut = crânial, dos = dorsal), donc justes même si notre stade est retourné (CS15).
"""
import sys, os, json, numpy as np, trimesh
from PIL import Image, ImageDraw, ImageFont

TR = sys.argv[1]
STADES = sys.argv[2:] or ['CS13', 'CS15', 'CS16', 'CS17', 'CS20']
DOS = {'CS13': 'CS13.f4v', 'CS15': 'CS15_f4v', 'CS16': 'CS16_f4v', 'CS17': 'CS17_f4v', 'CS20': 'CS20_F4V'}
OUT = 'embryons_3D/colonnes_amsterdam'
H = 1100                      # hauteur d'une vue (px)
FOND = (250, 250, 247)

# calques : (nom, source, [(structure, couleur, alpha)])
COL_A = {'vertebrae': (236, 230, 205), 'skeleton': (236, 230, 205), 'intervertebral_discs': (120, 190, 235), 'ribs': (215, 200, 160),
         'notochord': (250, 215, 60), 'somites': (150, 110, 205)}
COL_N = {'corps_vertebraux': (235, 110, 30), 'arcs_neuraux': (190, 60, 20), 'cotes': (225, 160, 70), 'squelette_axial_cartilage': (245, 150, 60),
         'notochorde': (205, 25, 25), 'somites': (250, 190, 90)}   # nous = tons chauds ; Amsterdam = tons froids / ivoire

def charge_A(cs, nom, T):
    p = f'embryons_3D/modeles/{cs}/{nom}.glb'
    if not os.path.exists(p): return None
    m = trimesh.load(p, force='mesh'); s, R, t = T
    m.vertices = (s * (R @ m.vertices.T)).T + t
    return m

def charge_N(cs, nom):
    p = f'{DOS[cs]}/out/{cs}_{nom}.ply'
    return trimesh.load(p, force='mesh') if os.path.exists(p) else None

def splat(meshes, basis, c, scale, W, Hh, n_par_mm2=9000, nmax=900000):
    """rendu par points ombrés + tampon de profondeur ; renvoie une image RGBA"""
    u, v, w = basis
    img = np.zeros((Hh, W, 4), np.uint8); zb = np.full((Hh, W), -np.inf)
    for m, col, alpha in meshes:
        n = int(min(nmax, max(20000, m.area * n_par_mm2)))
        P, fi = trimesh.sample.sample_surface(m, n, seed=1)
        N = m.face_normals[fi]
        sh = 0.30 + 0.70 * np.abs(N @ w)
        x = ((P - c) @ u) / scale + W / 2; y = Hh / 2 - ((P - c) @ v) / scale; z = (P - c) @ w
        rgb = np.clip(np.array(col)[None, :] * sh[:, None], 0, 255).astype(np.uint8)
        o = np.argsort(z)
        for dx in (0, 1):
            for dy in (0, 1):
                xi = (x[o] + dx).astype(int); yi = (y[o] + dy).astype(int); ok = (xi >= 0) & (xi < W) & (yi >= 0) & (yi < Hh)
                xi, yi, zz, cc = xi[ok], yi[ok], z[o][ok], rgb[o][ok]
                keep = zz >= zb[yi, xi]
                img[yi[keep], xi[keep], :3] = cc[keep]; img[yi[keep], xi[keep], 3] = int(255 * alpha); zb[yi[keep], xi[keep]] = zz[keep]
    return img

def calque_video(cs, rep, man, vue, c, scale, W, H, canvas):
    """Image de la vidéo 360° sur laquelle l'utilisateur a pointé (reperes_utilisateur.json : image_video, vue directe seulement),
    replacée dans la vue de profil : même recalage silhouette -> projection sagittale de l'enveloppe que reperes_video.py
    (similitude 2D), pixel vidéo -> (ax2, ax1) voxels -> point du plan médian -> pixel de la vue."""
    if not rep.get('vue', '').startswith('direct') or 'tourn' in rep.get('vue', ''):
        print(cs, ': vue vidéo non directe, calque vidéo non fait'); return None
    import cv2
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from video360 import trouver_video, silhouette, recalage_2d
    dossier = DOS[cs]; k = int(rep['image_video'])
    cap = cv2.VideoCapture(trouver_video(dossier)); cap.set(cv2.CAP_PROP_POS_FRAMES, k); ok, frame = cap.read(); cap.release()
    if not ok: return None
    z = np.load(os.path.join(dossier, 'work', 'labels.npz')); shape = tuple(int(q) for q in z['shape'])
    env = np.unpackbits(z['enveloppe'])[:int(np.prod(shape))].reshape(shape).astype(bool); z.close()
    A, err = recalage_2d(silhouette(cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)), env.any(axis=0), rot_max_deg=30)
    x_med = float(np.median(np.nonzero(env.any(axis=(1, 2)))[0]))
    vox = man['mm_per_voxel']; c0 = np.array(man['center_voxel']); cc0 = np.array([c0[0], c0[2], -c0[1]]) * vox
    _, u, v, _ = vue
    def vers_vue(xy):
        a2, a1 = A[:, :2] @ np.asarray(xy, float) + A[:, 2]
        P = np.array([x_med, a2, -a1]) * vox - cc0                        # repère des maillages : X = ax0, Y = ax2, Z = -ax1
        return [W / 2 + ((P - c) @ u) / scale, 60 + H / 2 - ((P - c) @ v) / scale]
    hv, wv = frame.shape[:2]; src = np.float32([[0, 0], [wv, 0], [0, hv]]); dst = np.float32([vers_vue(p) for p in src])
    M = cv2.getAffineTransform(src, dst)
    rgba = cv2.cvtColor(frame, cv2.COLOR_BGR2RGBA)
    out = cv2.warpAffine(rgba, M, canvas, flags=cv2.INTER_LINEAR, borderMode=cv2.BORDER_CONSTANT, borderValue=(0, 0, 0, 0))
    out[:, W:] = 0                                                         # seulement sous la vue de profil
    print(f"{cs} : image vidéo {k} replacée (recalage silhouette {err:.1f} px, {np.hypot(*M[:, 0]):.2f} px de vue par px vidéo)")
    return Image.fromarray(out)

def main():
    os.makedirs(f'{OUT}/calques', exist_ok=True)
    try: font = ImageFont.truetype('arial.ttf', 22); fpt = ImageFont.truetype('arial.ttf', 17)
    except Exception: font = fpt = ImageFont.load_default()
    for cs in STADES:
        f = f'{TR}/T_CS13_dos.npy' if cs == 'CS13' else f'{TR}/T_{cs}_queue.npy'
        T = np.load(f); T = (T[0], T[1:10].reshape(3, 3), T[10:])
        R = T[1]; lr, dors, cran = R[:, 0], R[:, 1], R[:, 2]          # axes anatomiques d'Amsterdam dans notre repère
        env = charge_N(cs, 'enveloppe'); skin = charge_A(cs, 'skin', T)
        calA = [(k, charge_A(cs, k, T), c) for k, c in COL_A.items()]; calA = [(k, m, c) for k, m, c in calA if m is not None]
        calN = [(k, charge_N(cs, k), c) for k, c in COL_N.items()]; calN = [(k, m, c) for k, m, c in calN if m is not None]
        # repères utilisateur
        man = json.load(open(f'{DOS[cs]}/out/manifest.json', encoding='utf-8')); vox = man['mm_per_voxel']; c0 = np.array(man['center_voxel'])
        cc0 = np.array([c0[0], c0[2], -c0[1]]) * vox; conv = lambda p: np.array([p[0], p[2], -p[1]]) - cc0
        rep = json.load(open(f'{DOS[cs]}/out/topographie/video360/reperes_utilisateur.json', encoding='utf-8'))
        pts = [('queue (pointage)', conv(rep['pointe_queue']['xyz_mm_pipeline'])), ('œil (pointage)', conv(rep['vesicule_optique']['xyz_mm_pipeline']))]
        allv = np.vstack([env.vertices, skin.vertices]); c = (allv.min(0) + allv.max(0)) / 2
        # 3 vues : profil (regard depuis la gauche anatomique), dos (depuis dorsal), face (depuis ventral) ; haut = crânial
        vues = [('profil', np.cross(cran, -lr), cran, -lr), ('dos', lr, cran, dors), ('face', -lr, cran, -dors)]
        vues = [(nm, u / np.linalg.norm(u), v - (v @ w) * w, w) for nm, u, v, w in vues]
        vues = [(nm, u, v / np.linalg.norm(v), w) for nm, u, v, w in vues]
        ext = max(max(np.ptp(allv @ u), np.ptp(allv @ v)) for _, u, v, _ in vues)
        scale = ext * 1.08 / H; W = H                                     # mm par pixel, identique pour les 3 vues
        couches = [('Fond', None)] + [('Nous - enveloppe', [(env, (225, 70, 60), .28)])] + \
                  [(f'Nous - {k}', [(m, c_, 1.0)]) for k, m, c_ in calN] + \
                  [('Amsterdam - peau', [(skin, (60, 130, 230), .28)])] + \
                  [(f'Amsterdam - {k}', [(m, c_, 1.0)]) for k, m, c_ in calA] + [('Repères et légende', None)]
        canvas = (3 * W, H + 60)
        layers = []
        for nom, meshes in couches:
            L = Image.new('RGBA', canvas, (0, 0, 0, 0))
            if nom == 'Fond': L = Image.new('RGBA', canvas, FOND + (255,))
            elif meshes is not None:
                for k, (vn, u, v, w) in enumerate(vues):
                    L.paste(Image.fromarray(splat(meshes, (u, v, w), c, scale, W, H)), (k * W, 60))
            else:
                d = ImageDraw.Draw(L)
                for k, (vn, u, v, w) in enumerate(vues):
                    d.text((k * W + 16, 14), f'{cs} — {vn}', fill=(20, 20, 20, 255), font=font)
                    if vn == 'profil':
                        for lab, p in pts:
                            x = k * W + W / 2 + ((p - c) @ u) / scale; y = 60 + H / 2 - ((p - c) @ v) / scale
                            d.ellipse([x - 9, y - 9, x + 9, y + 9], outline=(0, 150, 0, 255), width=3); d.text((x + 12, y - 10), lab, fill=(0, 110, 0, 255), font=fpt)
                    if k == 2:                                        # légende des couleurs
                        leg = [('Nous', None)] + [(nm_, c_) for nm_, _, c_ in calN] + [('Amsterdam', None)] + [(nm_, c_) for nm_, _, c_ in calA]
                        yy = 70
                        for nm_, c_ in leg:
                            if c_ is None: d.text((k * W + W - 290, yy), nm_, fill=(20, 20, 20, 255), font=font); yy += 30; continue
                            d.rectangle([k * W + W - 285, yy + 3, k * W + W - 265, yy + 19], fill=tuple(c_) + (255,), outline=(60, 60, 60, 255))
                            d.text((k * W + W - 255, yy), nm_, fill=(30, 30, 30, 255), font=fpt); yy += 24
                    L1 = 1.0 / scale; x0 = k * W + 30; y0 = 60 + H - 30
                    d.line([x0, y0, x0 + L1, y0], fill=(20, 20, 20, 255), width=4); d.text((x0, y0 - 26), '1 mm', fill=(20, 20, 20, 255), font=fpt)
            layers.append((nom, L))
            L.save(f'{OUT}/calques/{cs}_{nom.replace(" - ", "_").replace(" ", "_")}.png')
        comp = Image.new('RGBA', canvas)
        for _, L in layers: comp = Image.alpha_composite(comp, L)
        comp.convert('RGB').save(f'{OUT}/{cs}_colonne.png')
        LV = calque_video(cs, rep, man, vues[0], c, scale, W, H, canvas)       # image originale de la vidéo 360°, sous la vue de profil
        if LV is not None:
            layers.insert(1, (f"Vidéo 360° originale (image {rep['image_video']})", LV))
            LV.save(f'{OUT}/calques/{cs}_video_originale.png')
            comp = Image.new('RGBA', canvas)
            for _, L in layers: comp = Image.alpha_composite(comp, L)
            comp.convert('RGB').save(f'{OUT}/{cs}_colonne_video.png')
        from psd_tools import PSDImage
        from psd_tools.api.layers import PixelLayer
        psd = PSDImage.new('RGB', canvas)
        for nom, L in layers: psd.append(PixelLayer.frompil(L, psd, nom))
        psd.save(f'{OUT}/{cs}_colonne.psd')
        print(cs, 'calques :', [n for n, _ in layers], flush=True)

main()
