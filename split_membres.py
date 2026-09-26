"""Sépare le label 'membres' en 4 bourgeons nommés (membre_sup_gauche/droit, membre_inf_gauche/droit) pour un morphing pièce par pièce.
Les 4 plus grosses composantes connexes sont classées : moitié gauche/droite par rapport au plan médian LR de l'enveloppe, puis la plus
crâniale (SI plus petit) de chaque moitié = supérieur. Écrit dans labels.npz (les PLY sont régénérés par meshexport au prochain passage).
usage : python split_membres.py <dossier_stade> [...]"""
import numpy as np, scipy.ndimage as ndi, os, sys, shutil
NAMES = ('membre_sup_gauche', 'membre_sup_droit', 'membre_inf_gauche', 'membre_inf_droit')
def run(stage_dir):
    work = os.path.join(stage_dir, 'work'); p = os.path.join(work, 'labels.npz')
    with np.load(p) as z:
        shape = tuple(z['shape']); n = int(np.prod(shape)); L = {k: z[k] for k in z.files if k != 'shape'}
    if 'membres' not in L: print(stage_dir, ': pas de membres'); return
    m = np.unpackbits(L['membres'])[:n].reshape(shape).astype(bool)
    env = np.unpackbits(L['enveloppe'])[:n].reshape(shape).astype(bool); lr = np.where(env.any(axis=(1, 2)))[0]; mid = 0.5 * (lr.min() + lr.max())
    lab, nc = ndi.label(m); sizes = np.bincount(lab.ravel())[1:]; order = np.argsort(sizes)[::-1][:12] + 1
    cms = ndi.center_of_mass(m, lab, order); comps = [(i, c[0] < mid, c[1], int(sizes[i - 1])) for i, c in zip(order, cms)]     # (id, gauche ?, SI, taille)
    seuil = 0.05 * comps[0][3]
    out = {}
    for gauche in (True, False):
        side = sorted([c for c in comps if c[1] == gauche and c[3] >= seuil], key=lambda c: -c[3])[:2]   # les 2 plus grosses de ce côté
        side = sorted(side, key=lambda c: c[2])                                                          # crânial d'abord
        if len(side) >= 1: out['membre_sup_' + ('gauche' if gauche else 'droit')] = side[0][0]
        if len(side) >= 2: out['membre_inf_' + ('gauche' if gauche else 'droit')] = side[1][0]
    for k in NAMES:
        L[k] = np.packbits(lab == out[k]) if k in out else np.packbits(np.zeros(shape, bool))
    shutil.copy(p, os.path.join(work, 'labels_avant_membres.npz')); tmp = p + '.tmp.npz'
    np.savez_compressed(tmp, **L, shape=np.array(shape)); os.replace(tmp, p)
    print(os.path.basename(stage_dir), 'composantes', nc, '->', {k: (int(sizes[out[k] - 1]) if k in out else 0) for k in NAMES})
if __name__ == '__main__':
    for d in sys.argv[1:]: run(os.path.abspath(d))
