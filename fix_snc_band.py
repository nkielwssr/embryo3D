"""Retire du SNC ce qui déborde latéralement sous la tête (somites, peau) : sous la limite de la tête, le SNC est confiné à une
bande médiane |LR - milieu| <= frac * largeur, puis on garde la plus grande composante. Régénère maillages, NRRD, scène, rendu.
usage : python fix_snc_band.py <dossier_stade> [frac=0.13] [tete=0.33]"""
import numpy as np, os, sys, re, shutil, subprocess

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
stage_dir = os.path.abspath(sys.argv[1]); frac = float(sys.argv[2]) if len(sys.argv) > 2 else 0.13; tete = float(sys.argv[3]) if len(sys.argv) > 3 else 0.33
work = os.path.join(stage_dir, 'work'); out = os.path.join(stage_dir, 'out'); stage = 'CS%s' % re.search(r'CS\s*(\d+)', os.path.basename(stage_dir), re.I).group(1)
z = _load_npz(os.path.join(work, 'labels.npz')); shape = tuple(z['shape']); n = int(np.prod(shape))
labels = {k: np.unpackbits(z[k])[:n].reshape(shape).astype(bool) for k in z.files if k != 'shape'}
snc, env = labels['snc'], labels['enveloppe']
lr = np.where(env.any(axis=(1,2)))[0]; MID = 0.5*(lr.min()+lr.max()); LRe = lr.max()-lr.min()
si = np.where(env.any(axis=(0,2)))[0]; si0, si1 = si.min(), si.max(); lim = int(si0 + tete*(si1-si0))
band = np.abs(np.arange(shape[0]) - MID) <= frac*LRe
keep = snc.copy(); keep[~band, lim:, :] = False
lab, k = ndi.label(keep); sizes = np.bincount(lab.ravel())[1:]; keep = lab == (np.argmax(sizes)+1)
print('%s : snc %.2fM -> %.2fM (retiré %.2fM latéral, %d composantes -> 1)' % (stage, snc.sum()/1e6, keep.sum()/1e6, (snc & ~keep).sum()/1e6, k))
labels['snc'] = keep
shutil.copy(os.path.join(work, 'labels.npz'), os.path.join(work, 'labels_avant_fix_snc.npz'))
_savez_atomic(os.path.join(work, 'labels.npz'), **{kk: np.packbits(v) for kk, v in labels.items()}, shape=np.array(shape))
import meshexport, export_nrrd; meshexport.run(work, out, stage); export_nrrd.run(stage_dir)
B = r'C:/Program Files/Blender Foundation/Blender 5.2/blender.exe'
blend = os.path.join(out, f'{stage}_embryon.blend')
subprocess.run([B, '-b', '-P', os.path.join(HERE, 'blender_build_scene.py'), '--', os.path.join(out, 'manifest.json'), blend], capture_output=True)
subprocess.run([B, '-b', blend, '-P', os.path.join(HERE, 'blender_render_scene.py'), '--', os.path.join(out, f'{stage}_render_organs.png'), '90', '0'], capture_output=True)
# planche de contrôle mi-sagittale
import cv2
d = np.load(os.path.join(work, 'dens.npy'), mmap_mode='r'); i = int(MID)
img = cv2.cvtColor(255 - np.asarray(d[i]), cv2.COLOR_GRAY2BGR); img[keep[i]] = (0, 140, 255); img[(snc & ~keep)[i]] = (0, 0, 255)
j = int(MID + 0.2*LRe); img2 = cv2.cvtColor(255 - np.asarray(d[j]), cv2.COLOR_GRAY2BGR); img2[keep[j]] = (0, 140, 255); img2[(snc & ~keep)[j]] = (0, 0, 255)
cv2.imwrite(os.path.join(work, 'check_fix_snc.png'), np.hstack([img, img2])); print('done', stage)
