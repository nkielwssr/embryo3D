"""Ajoute des structures paires de la tête à partir de labels de cavités (cavo_lab) identifiés visuellement.
usage : python fix_pairs.py <dossier_stade> yeux <lab1> <lab2>              (cupule optique + cristallin autour de l'espace intra-oculaire)
        python fix_pairs.py <dossier_stade> vesicules_otiques <lab1> <lab2>  (cavité + paroi épithéliale)
Régénère maillages, NRRD, scène Blender et rendu."""
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
from segment import dilate_r, geodesic, open_r, close_r, comps_info, _crop
stage_dir = os.path.abspath(sys.argv[1]); kind = sys.argv[2]; labs = [int(x) for x in sys.argv[3:]]
work = os.path.join(stage_dir, 'work'); out = os.path.join(stage_dir, 'out'); stage = 'CS%s' % re.search(r'CS\s*(\d+)', os.path.basename(stage_dir), re.I).group(1)
z = _load_npz(os.path.join(work, 'labels.npz')); shape = tuple(z['shape']); n = int(np.prod(shape))
labels = {k: np.unpackbits(z[k])[:n].reshape(shape).astype(bool) for k in z.files if k != 'shape'}
ds = np.load(os.path.join(work, 'ds1.npy')).astype(np.float32); env = labels['enveloppe']; cavo = np.load(os.path.join(work, 'cavo_lab.npy'))
cns = labels.get('snc', np.zeros(shape, bool)); cav = env & (ds < 3.0)
seed = geodesic(np.isin(cavo, labs), cav, 4)
if kind == 'yeux':
    eye = geodesic(dilate_r(seed, 4) & (ds > 50) & ~cns, (ds > 50) & env & ~cns, 14)
    eye = ndi.binary_fill_holes(close_r(eye | seed, 3))
    lab, info = comps_info(eye, 500); info.sort(key=lambda x: -x[1]); eye = np.isin(lab, [i[0] for i in info[:2]])
    lensc = open_r(eye & (ds > 110), 2); lab, info = comps_info(lensc, 300); info.sort(key=lambda x: -x[1])
    lens = np.isin(lab, [i[0] for i in info[:2]]) if info else lensc
    labels['yeux'] = eye; labels['cristallins'] = lens
    print('%s yeux %.3f mm3-eq %d vox, cristallins %d vox' % (stage, 0, eye.sum(), lens.sum()))
    new = eye
elif kind == 'vesicules_otiques':
    otic = seed | (dilate_r(seed, 6) & (ds > 60) & env & ~cns)
    otic = close_r(otic, 2); labels['vesicules_otiques'] = otic; new = otic
    print('%s vesicules_otiques %d vox' % (stage, otic.sum()))
else:
    raise SystemExit('kind inconnu')
shutil.copy(os.path.join(work, 'labels.npz'), os.path.join(work, 'labels_avant_fix_%s.npz' % kind))
_savez_atomic(os.path.join(work, 'labels.npz'), **{kk: np.packbits(v) for kk, v in labels.items()}, shape=np.array(shape))
import meshexport, export_nrrd; meshexport.run(work, out, stage); export_nrrd.run(stage_dir)
B = r'C:/Program Files/Blender Foundation/Blender 5.2/blender.exe'; blend = os.path.join(out, f'{stage}_embryon.blend')
subprocess.run([B, '-b', '-P', os.path.join(HERE, 'blender_build_scene.py'), '--', os.path.join(out, 'manifest.json'), blend], capture_output=True)
subprocess.run([B, '-b', blend, '-P', os.path.join(HERE, 'blender_render_scene.py'), '--', os.path.join(out, f'{stage}_render_organs.png'), '90', '0'], capture_output=True)
d = np.load(os.path.join(work, 'dens.npy'), mmap_mode='r'); zz, yy, xx = np.where(new); i = int(np.median(zz[zz < np.median(zz)])) if len(zz) else shape[0]//2
img = cv2.cvtColor(255 - np.asarray(d[i]), cv2.COLOR_GRAY2BGR); img[new[i]] = (255, 0, 0)
if kind == 'yeux': img[labels['cristallins'][i]] = (0, 255, 255)
j = int(np.median(yy)); img2 = cv2.cvtColor(255 - np.asarray(d[:, j, :]), cv2.COLOR_GRAY2BGR); img2[new[:, j, :]] = (255, 0, 0)
h = img.shape[0]; img2 = cv2.resize(img2, (int(img2.shape[1]*h/img2.shape[0]), h))
cv2.imwrite(os.path.join(work, 'check_fix_%s.png' % kind), np.hstack([img, img2])); print('done', stage, kind)
