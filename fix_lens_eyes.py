"""Yeux par le cristallin : boules denses compactes, symétriques, latérales, dans la tête -> 'cristallins' ; oeil = coque dense autour (<= r) -> 'yeux'.
usage : python fix_lens_eyes.py <dossier_stade> [--diag]"""
import numpy as np, os, sys, re, shutil, subprocess, cv2

class _Npz(dict):
    """contenu d'un .npz chargé en mémoire et FERMÉ (np.load garde le fichier ouvert, ce qui bloque os.replace sous Windows)"""
    files = property(lambda self: list(self.keys()))
def _load_npz(path):
    import numpy as _np
    with _np.load(path) as f: return _Npz({k: f[k] for k in f.files})


def _savez_atomic(path, **kw):
    """écriture atomique (fichier temporaire puis remplacement) : d'autres sessions lisent labels.npz en parallèle"""
    import numpy as _np, os as _os; tmp = path + '.tmp.npz'; _np.savez_compressed(tmp, **kw); _os.replace(tmp, path)

from scipy import ndimage as ndi
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from segment import dilate_r, open_r, close_r, comps_info, _crop, geodesic
stage_dir = os.path.abspath(sys.argv[1]); diag = '--diag' in sys.argv
work = os.path.join(stage_dir, 'work'); out = os.path.join(stage_dir, 'out'); stage = 'CS%s' % re.search(r'CS\s*(\d+)', os.path.basename(stage_dir), re.I).group(1)
z = _load_npz(os.path.join(work, 'labels.npz')); shape = tuple(z['shape']); n = int(np.prod(shape))
labels = {k: np.unpackbits(z[k])[:n].reshape(shape).astype(bool) for k in z.files if k != 'shape'}
ds = np.load(os.path.join(work, 'ds1.npy')).astype(np.float32); env = labels['enveloppe']; snc = labels.get('snc', np.zeros(shape, bool))
lr = np.where(env.any(axis=(1,2)))[0]; MID = 0.5*(lr.min()+lr.max()); LRe = lr.max()-lr.min()
si = np.where(env.any(axis=(0,2)))[0]; si0, si1 = si.min(), si.max(); SIe = si1-si0
head = np.zeros(shape, bool); head[:, :int(si0 + 0.45*SIe), :] = True
lateral = np.zeros(shape, bool); lateral[np.abs(np.arange(shape[0]) - MID) > 0.08*LRe, :, :] = True
r_open = max(3, int(0.006 * SIe))
balls = open_r((ds > 100) & env & head & lateral & ~snc, r_open)
lab, info = comps_info(balls, 500); info = [(i, s, c) for i, s, c in info if s < 0.002 * env.sum()]
cands = []
for i, s, c in info:
    m = lab == i; sl = _crop(m, 1); zz, yy, xx = np.where(m[sl]); P = np.stack([zz, yy, xx], 1).astype(np.float32); P -= P.mean(0)
    ev = np.sort(np.linalg.eigvalsh(np.cov(P.T)))[::-1]; el = float(np.sqrt(ev[0] / max(ev[2], 1e-6)))
    if el < 2.4: cands.append((i, s, c, el))
for i, s, c, el in sorted(cands, key=lambda x: -x[1])[:12]: print('  boule %d : %d vox centroide %s elong %.1f' % (i, s, c.astype(int), el))
best = None
for a in cands:
    for b in cands:
        if a[0] >= b[0]: continue
        ca, cb = a[2], b[2]
        if abs((ca[0]-MID) + (cb[0]-MID)) < 0.08*LRe and abs(ca[0]-cb[0]) > 0.16*LRe and abs(ca[1]-cb[1]) < 0.06*SIe and abs(ca[2]-cb[2]) < 0.08*shape[2]:
            score = a[1] + b[1]
            if best is None or score > best[0]: best = (score, a, b)
if best is None: raise SystemExit('%s : aucune paire de cristallins trouvée (%d boules candidates)' % (stage, len(cands)))
_, a, b = best; lens = (lab == a[0]) | (lab == b[0])
print('%s : cristallins labs %d/%d, %d + %d vox, centroïdes %s %s' % (stage, a[0], b[0], a[1], b[1], a[2].astype(int), b[2].astype(int)))
r_eye = max(8, int(0.02 * SIe))
eye = geodesic(dilate_r(lens, 3) & (ds > 55) & env & ~snc, (ds > 55) & env & ~snc & dilate_r(lens, r_eye), 1000)
eye = ndi.binary_fill_holes(close_r(eye | lens, 3))
lab2, info2 = comps_info(eye, 500); info2.sort(key=lambda x: -x[1]); eye = np.isin(lab2, [i[0] for i in info2[:2]])
print('  yeux %d vox' % eye.sum())
d = np.load(os.path.join(work, 'dens.npy'), mmap_mode='r'); i = int(a[2][0]); j = int(a[2][1])
img = cv2.cvtColor(255 - np.asarray(d[i]), cv2.COLOR_GRAY2BGR); img[eye[i]] = (255, 0, 0); img[lens[i]] = (0, 255, 255)
img2 = cv2.cvtColor(255 - np.asarray(d[:, j, :]), cv2.COLOR_GRAY2BGR); img2[eye[:, j, :]] = (255, 0, 0); img2[lens[:, j, :]] = (0, 255, 255)
h = img.shape[0]; img2 = cv2.resize(img2, (int(img2.shape[1]*h/img2.shape[0]), h)); cv2.imwrite(os.path.join(work, 'check_fix_lens.png'), np.hstack([img, img2]))
if not diag:
    labels['cristallins'] = lens; labels['yeux'] = eye
    shutil.copy(os.path.join(work, 'labels.npz'), os.path.join(work, 'labels_avant_fix_lens.npz'))
    _savez_atomic(os.path.join(work, 'labels.npz'), **{kk: np.packbits(v) for kk, v in labels.items()}, shape=np.array(shape))
    import meshexport, export_nrrd; meshexport.run(work, out, stage); export_nrrd.run(stage_dir)
    B = r'C:/Program Files/Blender Foundation/Blender 5.2/blender.exe'; blend = os.path.join(out, f'{stage}_embryon.blend')
    subprocess.run([B, '-b', '-P', os.path.join(HERE, 'blender_build_scene.py'), '--', os.path.join(out, 'manifest.json'), blend], capture_output=True)
    subprocess.run([B, '-b', blend, '-P', os.path.join(HERE, 'blender_render_scene.py'), '--', os.path.join(out, f'{stage}_render_organs.png'), '90', '0'], capture_output=True)
print('done', stage)
