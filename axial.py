"""Squelette axial en développement : colonne paraxiale le long du tube neural.
 - stades jeunes (CS13-16) : 'somites' = blocs denses métamériques flanquant le tube (séparés aux fentes intersomitiques) ;
 - stades âgés (CS17-20) : 'squelette_axial_cartilage' complété par les arcs neuraux / côtes (blobs pâles compacts du couloir paraxial) ;
 - tous : 'notochorde' = cordon dense médian juste ventral au tube.
usage : python axial.py <dossier_stade> [--diag]     (--diag : planches seulement, pas d'écriture)"""
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
from segment import dilate_r, erode_r, open_r, close_r, comps_info, _crop, _dist
stage_dir = os.path.abspath(sys.argv[1]); diag = '--diag' in sys.argv
SEUIL = float(next((a.split('=')[1] for a in sys.argv if a.startswith('--seuil=')), 95))
DOUT = float(next((a.split('=')[1] for a in sys.argv if a.startswith('--dout=')), 0))          # fraction de LRe (0 = défaut)
RMEMB = int(next((a.split('=')[1] for a in sys.argv if a.startswith('--rmembres=')), 3))      # marge autour des membres (voxels)
work = os.path.join(stage_dir, 'work'); out = os.path.join(stage_dir, 'out'); stage = 'CS%s' % re.search(r'CS\s*(\d+)', os.path.basename(stage_dir), re.I).group(1)
num = int(stage[2:]); young = num <= 16
z = _load_npz(os.path.join(work, 'labels.npz')); shape = tuple(z['shape']); n = int(np.prod(shape))
labels = {k: np.unpackbits(z[k])[:n].reshape(shape).astype(bool) for k in z.files if k != 'shape'}
ds = np.load(os.path.join(work, 'ds1.npy')).astype(np.float32); env = labels['enveloppe']; snc = labels['snc']
lr = np.where(env.any(axis=(1,2)))[0]; MID = 0.5*(lr.min()+lr.max()); LRe = lr.max()-lr.min()
si = np.where(env.any(axis=(0,2)))[0]; si0, si1 = si.min(), si.max(); SIe = si1-si0; lim = int(si0 + 0.33*SIe)
# tube neural sous la tête (bande médiane)
mid = np.abs(np.arange(shape[0]) - MID)
cord = snc.copy(); cord[:, :lim, :] = False; cord[mid > 0.13*LRe, :, :] = False
lab, k = ndi.label(cord);
if k: cord = lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)
dcord = _dist(~cord)                                   # distance au tube (voxels)
# demi-largeur du tube -> rayons du couloir
half = float(np.median(_dist(cord)[cord])) * 2 if cord.any() else 8.0
d_in, d_out = half * 0.4, (DOUT if DOUT > 0 else (0.22 if young else 0.20)) * LRe
lateral = np.zeros(shape, bool); lateral[mid > 0.035*LRe, :, :] = True
corridor = (dcord > d_in) & (dcord < d_out) & env & ~snc & lateral; corridor[:, :lim, :] = False
excl = np.zeros(shape, bool)
for kx in ('foie', 'coeur', 'coeur_detoure', 'membres', 'cordon_ombilical', 'yeux', 'vesicules_otiques', 'ventricules'):
    if kx in labels: excl |= dilate_r(labels[kx], RMEMB if kx == 'membres' else 3)
corridor &= ~excl
print('%s : tube %.2fM vox, demi-largeur %.1f, couloir %.1f..%.1f vox, %.2fM vox' % (stage, cord.sum()/1e6, half, d_in, d_out, corridor.sum()/1e6))
new = {}
if young:
    som = open_r((ds > SEUIL) & corridor, 2)
    lab, info = comps_info(som, 1500); info.sort(key=lambda x: -x[1])
    keep = [i for i, s, c in info if s < 0.02 * env.sum()]
    som = np.isin(lab, keep); new['somites'] = som
    print('  somites : %d blocs, %.2fM vox' % (len(keep), som.sum()/1e6))
else:
    pale = open_r((ds < 8) & corridor, 3)                      # cartilage : quasi non coloré, compact
    lab, info = comps_info(pale, 800); info.sort(key=lambda x: -x[1]); keep = []
    for i, s, c in info:
        if s > 0.005 * env.sum(): continue
        m = lab == i; sl = _crop(m, 1); zz, yy, xx = np.where(m[sl]); P = np.stack([zz, yy, xx], 1).astype(np.float32); P -= P.mean(0)
        ev = np.sort(np.linalg.eigvalsh(np.cov(P.T)))[::-1]
        if np.sqrt(ev[0] / max(ev[2], 1e-6)) < 8: keep.append(i)
    arcs = np.isin(lab, keep)
    base = labels.get('squelette_axial_cartilage', np.zeros(shape, bool))
    pav = os.path.join(work, 'labels_avant_axial.npz')
    if os.path.exists(pav):                                   # relance : repartir des corps vertébraux d'origine
        za = np.load(pav)
        if 'squelette_axial_cartilage' in za.files: base = np.unpackbits(za['squelette_axial_cartilage'])[:n].reshape(shape).astype(bool)
    new['squelette_axial_cartilage'] = base | arcs
    print('  arcs/côtes : %d blocs, %.2fM vox (+ corps vertébraux %.2fM)' % (len(keep), arcs.sum()/1e6, base.sum()/1e6))
# notochorde : médian, ventral au tube, dense, fin et long
midline = np.zeros(shape, bool); midline[mid <= 0.02*LRe, :, :] = True
ventral_band = (dcord > 1) & (dcord < max(half*1.2, 12)) & midline & env & ~snc; ventral_band[:, :lim, :] = False
noto = open_r((ds > 90) & ventral_band, 1)
lab, info = comps_info(noto, 2000); info.sort(key=lambda x: -x[1])
if info:
    i, s, c = info[0]; noto = lab == i
    # allongé ?
    sl = _crop(noto, 1); zz, yy, xx = np.where(noto[sl]); P = np.stack([zz, yy, xx], 1).astype(np.float32); P -= P.mean(0); ev = np.sort(np.linalg.eigvalsh(np.cov(P.T)))[::-1]
    el = float(np.sqrt(ev[0] / max(ev[2], 1e-6))); print('  notochorde candidate %d vox, élongation %.1f' % (s, el))
    if el > 4: new['notochorde'] = noto
# planches
d = np.load(os.path.join(work, 'dens.npy'), mmap_mode='r')
def paint(i, transpose=False):
    sl = (slice(None), slice(None), i) if transpose else (i,)
    img = cv2.cvtColor(255 - np.asarray(d[sl]), cv2.COLOR_GRAY2BGR)
    if transpose: img = np.ascontiguousarray(img.transpose(1, 0, 2))
    for kx, col in (('somites', (255, 160, 60)), ('squelette_axial_cartilage', (200, 200, 255)), ('notochorde', (0, 220, 220))):
        if kx in new:
            m = new[kx][sl]; m = m.T if transpose else m; img[m] = col
    return img
i0 = int(MID); i1 = int(MID + 0.08*LRe); i2 = int(MID + 0.15*LRe)
ap = int(np.mean(np.where(cord.any(axis=(0, 1)))[0])) if cord.any() else shape[2]//2
cv2.imwrite(os.path.join(work, 'check_axial.png'), np.hstack([paint(i0), paint(i1), paint(i2), paint(ap, True)]))
if not diag:
    labels.update(new)
    shutil.copy(os.path.join(work, 'labels.npz'), os.path.join(work, 'labels_avant_axial.npz'))
    _savez_atomic(os.path.join(work, 'labels.npz'), **{kk: np.packbits(v) for kk, v in labels.items()}, shape=np.array(shape))
    import meshexport, export_nrrd; meshexport.run(work, out, stage); export_nrrd.run(stage_dir)
    B = r'C:/Program Files/Blender Foundation/Blender 5.2/blender.exe'; blend = os.path.join(out, f'{stage}_embryon.blend')
    subprocess.run([B, '-b', '-P', os.path.join(HERE, 'blender_build_scene.py'), '--', os.path.join(out, 'manifest.json'), blend], capture_output=True)
    subprocess.run([B, '-b', blend, '-P', os.path.join(HERE, 'blender_render_scene.py'), '--', os.path.join(out, f'{stage}_render_organs.png'), '90', '0'], capture_output=True)
print('done', stage)
