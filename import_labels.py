"""Réinjecte une carte de labels corrigée dans 3D Slicer (<CS>_labels_corrige.nrrd + .json) dans work/labels.npz,
puis régénère les maillages, la scène Blender du stade et un rendu.   usage : python import_labels.py <dossier_stade> [fichier.nrrd]"""
import numpy as np, nrrd, json, os, sys, re, subprocess, shutil

class _Npz(dict):
    """contenu d'un .npz chargé en mémoire et FERMÉ (np.load garde le fichier ouvert, ce qui bloque os.replace sous Windows)"""
    files = property(lambda self: list(self.keys()))
def _load_npz(path):
    import numpy as _np
    with _np.load(path) as f: return _Npz({k: f[k] for k in f.files})


def _savez_atomic(path, **kw):
    """écriture atomique (fichier temporaire puis remplacement) : d'autres sessions lisent labels.npz en parallèle"""
    import numpy as _np, os as _os; tmp = path + '.tmp.npz'; _np.savez_compressed(tmp, **kw); _os.replace(tmp, path)

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
stage_dir = os.path.abspath(sys.argv[1]); work = os.path.join(stage_dir, 'work'); out = os.path.join(stage_dir, 'out')
stage = 'CS%s' % re.search(r'CS\s*(\d+)', os.path.basename(stage_dir), re.I).group(1)
path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(out, f'{stage}_labels_corrige.nrrd')
mapping_path = os.path.splitext(path)[0] + '.json'
mapping = json.load(open(mapping_path)) if os.path.exists(mapping_path) else json.load(open(os.path.join(out, f'{stage}_label_ids.json')))
arr, hdr = nrrd.read(path)
z = _load_npz(os.path.join(work, 'labels.npz')); shape = tuple(z['shape'])
# Slicer peut réordonner / retourner les axes à l'export : on réaligne sur la géométrie du NRRD de référence via 'space directions'
ref_path = os.path.join(out, f'{stage}_labels.nrrd')
if os.path.exists(ref_path) and 'space directions' in hdr:
    _, href = nrrd.read_header(ref_path) if False else (None, nrrd.read_header(ref_path))
    Dref = np.array(href['space directions'], float); Dcor = np.array(hdr['space directions'], float)
    perm = []; flips = []
    for j in range(3):                                   # axe j de la référence <- quel axe i du fichier corrigé ?
        dots = [float(np.dot(Dcor[i], Dref[j])) / (np.linalg.norm(Dcor[i]) * np.linalg.norm(Dref[j]) + 1e-12) for i in range(3)]
        i = int(np.argmax(np.abs(dots))); perm.append(i); flips.append(dots[i] < 0)
    arr = np.transpose(arr, perm)
    for ax, f in enumerate(flips):
        if f: arr = np.flip(arr, axis=ax)
    arr = np.ascontiguousarray(arr); print('réalignement axes : permutation', perm, 'retournements', flips)
if arr.shape != shape: raise SystemExit('forme %s incompatible avec %s' % (arr.shape, shape))
old = {k: np.unpackbits(z[k])[:int(np.prod(shape))].reshape(shape).astype(bool) for k in z.files if k != 'shape'}
labels = {'enveloppe': old['enveloppe']}
from collections import Counter
dupl = {v for v, c in Counter(mapping.values()).items() if c > 1 and v != 0}
for name, val in mapping.items():
    if val == 0 or val in dupl:
        print('  %-28s ignoré (valeur %d %s)' % (name, val, 'vide' if val == 0 else 'partagée par plusieurs structures')); continue
    m = arr == val
    if m.any(): labels[name] = m
    print('  %-28s %8d voxels' % (name, m.sum()))
# une carte de labels n'a qu'une valeur par voxel : on restaure les inclusions connues (le cristallin est dans l'oeil)
if 'yeux' in labels and 'cristallins' in labels: labels['yeux'] = labels['yeux'] | labels['cristallins']
shutil.copy(os.path.join(work, 'labels.npz'), os.path.join(work, 'labels_avant_correction.npz'))
_savez_atomic(os.path.join(work, 'labels.npz'), **{k: np.packbits(v) for k, v in labels.items()}, shape=np.array(shape))
print('labels.npz mis à jour ; export des maillages...')
import meshexport, export_nrrd; meshexport.run(work, out, stage); export_nrrd.run(stage_dir)
B = r'C:/Program Files/Blender Foundation/Blender 5.2/blender.exe'
if os.path.exists(B):
    blend = os.path.join(out, f'{stage}_embryon.blend')
    subprocess.run([B, '-b', '-P', os.path.join(HERE, 'blender_build_scene.py'), '--', os.path.join(out, 'manifest.json'), blend], capture_output=True)
    subprocess.run([B, '-b', blend, '-P', os.path.join(HERE, 'blender_render_scene.py'), '--', os.path.join(out, f'{stage}_render_organs.png'), '90', '0'], capture_output=True)
    print('scène Blender et rendu régénérés :', blend)
print('Pour la scène maître à 7 stades : bash embryo3d/build_all.sh')
