"""Prolonge le SNC le long du tube neural : noyau dense (>= seuil) érodé, dans la bande médiane sous la tête, connecté au SNC existant.
usage : python fix_snc_tube.py <dossier_stade> [seuil=100] [erosion=2] [frac=0.13] [tete=0.33] [queue=0.80]  (rien n'est ajouté sous la fraction queue)"""
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
from segment import dilate_r, erode_r, geodesic
stage_dir = os.path.abspath(sys.argv[1]); a = sys.argv[2:] + [None]*5
seuil = float(a[0] or 100); ero = int(a[1] or 2); frac = float(a[2] or 0.13); tete = float(a[3] or 0.33); queue = float(a[4] or 0.80)
work = os.path.join(stage_dir, 'work'); out = os.path.join(stage_dir, 'out'); stage = 'CS%s' % re.search(r'CS\s*(\d+)', os.path.basename(stage_dir), re.I).group(1)
z = _load_npz(os.path.join(work, 'labels.npz')); shape = tuple(z['shape']); n = int(np.prod(shape))
labels = {k: np.unpackbits(z[k])[:n].reshape(shape).astype(bool) for k in z.files if k != 'shape'}
ds = np.load(os.path.join(work, 'ds1.npy')).astype(np.float32); env = labels['enveloppe']; snc = labels['snc']
lr = np.where(env.any(axis=(1,2)))[0]; MID = 0.5*(lr.min()+lr.max()); LRe = lr.max()-lr.min()
si = np.where(env.any(axis=(0,2)))[0]; si0, si1 = si.min(), si.max(); lim = int(si0 + tete*(si1-si0))
band = np.zeros(shape, bool); band[np.abs(np.arange(shape[0]) - MID) <= frac*LRe, :, :] = True; band[:, :lim, :] = True
band[:, int(si0 + queue*(si1-si0)):, :] = False        # pas d'ajout dans la queue (tissu uniformément dense)
deep = erode_r(env, 4)
core = erode_r((ds >= seuil) & deep & band, ero)
tube = geodesic(dilate_r(snc, 8) & core, core, 1000)
add = geodesic(dilate_r(tube, ero + 1) & (ds >= 70) & env & band, (ds >= 70) & env & band, 2)
new = snc | add
lab, k = ndi.label(new); sizes = np.bincount(lab.ravel())[1:]; new = lab == (np.argmax(sizes)+1)
for other in ('foie', 'coeur', 'yeux', 'membres', 'cordon_ombilical'):
    if other in labels: new &= ~labels[other]
print('%s : snc %.2fM -> %.2fM (tube ajouté %.2fM)' % (stage, snc.sum()/1e6, new.sum()/1e6, (new & ~snc).sum()/1e6))
labels['snc'] = new
shutil.copy(os.path.join(work, 'labels.npz'), os.path.join(work, 'labels_avant_fix_tube.npz'))
_savez_atomic(os.path.join(work, 'labels.npz'), **{kk: np.packbits(v) for kk, v in labels.items()}, shape=np.array(shape))
import meshexport, export_nrrd; meshexport.run(work, out, stage); export_nrrd.run(stage_dir)
B = r'C:/Program Files/Blender Foundation/Blender 5.2/blender.exe'; blend = os.path.join(out, f'{stage}_embryon.blend')
subprocess.run([B, '-b', '-P', os.path.join(HERE, 'blender_build_scene.py'), '--', os.path.join(out, 'manifest.json'), blend], capture_output=True)
subprocess.run([B, '-b', blend, '-P', os.path.join(HERE, 'blender_render_scene.py'), '--', os.path.join(out, f'{stage}_render_organs.png'), '90', '0'], capture_output=True)
d = np.load(os.path.join(work, 'dens.npy'), mmap_mode='r'); i = int(MID)
img = cv2.cvtColor(255 - np.asarray(d[i]), cv2.COLOR_GRAY2BGR); img[new[i]] = (0, 140, 255); img[(new & ~snc)[i]] = (0, 200, 0)
cv2.imwrite(os.path.join(work, 'check_fix_tube.png'), img); print('done', stage)
