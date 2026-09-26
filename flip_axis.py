"""Retourne manuellement un stade canonique le long d'un axe (0=LR, 1=SI, 2=AP) : usage python flip_axis.py <work_dir> <axis>"""
import numpy as np, sys, os, json
work, ax = sys.argv[1], int(sys.argv[2])
for f in ['dens.npy', 'ds1.npy', 'env.npy', 'cavo_lab.npy', 'liver_core.npy']:
    p = os.path.join(work, f)
    if os.path.exists(p):
        a = np.load(p); np.save(p, np.ascontiguousarray(np.flip(a, axis=ax))); n = a.shape[ax]; print('flipped', f)
rows = np.load(os.path.join(work, 'cavo_rows.npy'))
rows[:, 2+ax] = (n - 1) - rows[:, 2+ax]; np.save(os.path.join(work, 'cavo_rows.npy'), rows)
axes = json.load(open(os.path.join(work, 'axes.json'))); axes.setdefault('manual_flips', []).append(ax); axes['flips'][ax] = not axes['flips'][ax]
json.dump(axes, open(os.path.join(work, 'axes.json'), 'w'), indent=1); print('axes.json updated', axes['flips'], 'manual', axes['manual_flips'])
