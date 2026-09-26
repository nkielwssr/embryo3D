"""(Ré)écrit les NRRD d'un stade pour 3D Slicer : densité + carte de labels, orientation anatomique (tête en haut, face devant),
+ table de couleurs .ctbl et label_ids.json.   usage : python export_nrrd.py <dossier_stade> [...]"""
import numpy as np, nrrd, json, os, sys, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from meshexport import CRL_MM, STYLE

def run(stage_dir):
    work = os.path.join(stage_dir, 'work'); out = os.path.join(stage_dir, 'out')
    stage = 'CS%s' % re.search(r'CS\s*(\d+)', os.path.basename(stage_dir), re.I).group(1)
    z = np.load(os.path.join(work, 'labels.npz')); shape = tuple(z['shape'])
    names = [k for k in z.files if k != 'shape']
    labels = {k: np.unpackbits(z[k])[:int(np.prod(shape))].reshape(shape).astype(bool) for k in names}
    env = labels['enveloppe']; si = np.where(env.any(axis=(0,2)))[0]; scale = CRL_MM.get(stage, 20.0) / (np.ptp(si) + 1)
    c = np.array([np.mean(np.where(env.any(axis=(1,2)))[0]), si.mean(), np.mean(np.where(env.any(axis=(0,1)))[0])])
    d = np.load(os.path.join(work, 'dens.npy'))
    lm = np.zeros(shape, np.uint8); ids = {}
    for i, name in enumerate([n for n in names if n != 'enveloppe']):
        lm[labels[name]] = i + 1; ids[name] = i + 1
    # tableau (LR, SI, AP) : axe0 = +droite, axe1 = +bas, axe2 = +dorsal.  Espace LPS : L = -axe0, P = +axe2, S = -axe1
    dirs = np.array([[-scale, 0, 0], [0, 0, -scale], [0, scale, 0]])
    origin = -(dirs.T @ c)                                    # le centre du corps à l'origine
    hdr = {'space': 'left-posterior-superior', 'space directions': dirs, 'space origin': origin, 'kinds': ['domain', 'domain', 'domain']}
    nrrd.write(os.path.join(out, f'{stage}_density.nrrd'), np.ascontiguousarray(d), dict(hdr), compression_level=2)
    nrrd.write(os.path.join(out, f'{stage}_labels.nrrd'), np.ascontiguousarray(lm), dict(hdr), compression_level=2)
    json.dump(ids, open(os.path.join(out, f'{stage}_label_ids.json'), 'w'), indent=1)
    with open(os.path.join(out, f'{stage}_labels.ctbl'), 'w') as f:
        f.write('0 fond 0 0 0 0\n')
        for name, v in ids.items():
            col = STYLE.get(name, (None, (0.7, 0.7, 0.7)))[1]
            f.write('%d %s %d %d %d 255\n' % (v, name, *[int(255*x) for x in col]))
    print(stage, 'NRRD écrits :', shape, 'mm/voxel %.4f' % scale, len(ids), 'labels')

if __name__ == '__main__':
    for d in sys.argv[1:]: run(os.path.abspath(d))
