"""Planche de contrôle du morphing : rend les stades demandés (blender_render_morph.py, cadrage fixe sur CS20) puis les assemble en une image.
usage : python planche_morph.py sortie.png [0 0.5 1 ...]"""
import subprocess, sys, os, cv2, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
env_opt = ['--enveloppe'] if '--enveloppe' in sys.argv else []; sys.argv = [a for a in sys.argv if a != '--enveloppe']
out = sys.argv[1]; stades = sys.argv[2:] or ['0', '0.5', '1', '2', '3', '3.5', '4', '5', '5.5', '6']
B = r'C:/Program Files/Blender Foundation/Blender 5.2/blender.exe'; blend = os.path.join(ROOT, 'embryons_3D', 'embryons_CS13-CS20_morph.blend')
tmp = os.path.join(ROOT, 'embryons_3D', 'rendus', 'planche_tmp'); os.makedirs(tmp, exist_ok=True)
r = subprocess.run([B, '-b', blend, '-P', os.path.join(HERE, 'blender_render_morph.py'), '--'] + env_opt + [os.path.join(tmp, 'p')] + stades, capture_output=True, text=True, errors='replace')
if 'Traceback' in r.stdout + r.stderr: print(r.stdout[-1500:], r.stderr[-800:])
ims = []
for st in stades:
    f = os.path.join(tmp, 'p_stage%.2f.png' % float(st)); im = cv2.imread(f)
    if im is None: print('manque', f); continue
    im = cv2.resize(im, (450, 550)); cv2.putText(im, 'stage %s' % st, (10, 30), cv2.FONT_HERSHEY_SIMPLEX, 0.9, (30, 30, 30), 2); ims.append(im)
cols = 5; rows = (len(ims) + cols - 1) // cols
while len(ims) % cols: ims.append(np.full_like(ims[0], 255))
grid = np.vstack([np.hstack(ims[i*cols:(i+1)*cols]) for i in range(rows)]); cv2.imwrite(out, grid); print('PLANCHE', out, grid.shape)
