"""Recalcule vent_lab (système ventriculaire) : plus grande cavité (>=60k vox) dont le centroïde est dans le tiers supérieur du corps,
en préférant une paroi dense (score) ; usage: python fix_vent.py <work_dir> [...]"""
import numpy as np, sys, os, json, cv2
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from segment import _crop, _dist
tiles = []
for work in sys.argv[1:]:
    rows = np.load(os.path.join(work, 'cavo_rows.npy')); env = np.load(os.path.join(work, 'env.npy'), mmap_mode='r')
    ds = np.load(os.path.join(work, 'ds1.npy'), mmap_mode='r'); cavo = np.load(os.path.join(work, 'cavo_lab.npy'), mmap_mode='r')
    si = np.where(np.asarray(env).any(axis=(0,2)))[0]; a, b = si.min(), si.max(); head_lim = a + 0.38*(b-a)
    cands = [r for r in rows if r[1] >= 20000 and r[3] < head_lim]
    best = None; good = []
    for r in sorted(cands, key=lambda r: -r[1])[:10]:
        m = np.asarray(cavo) == int(r[0]); sl = _crop(m, 16); mm = m[sl]; dsc = np.asarray(ds[sl]).astype(np.float32); dist = _dist(~mm)
        prof = [float(np.median(dsc[(dist > k-1) & (dist <= k)])) for k in range(4, 15)]; score = float(np.mean(prof))
        val = np.log10(r[1]) * float(np.clip(score/60.0, 0.3, 1.6))
        print('%s lab %d size %dk SI %d score %.0f val %.2f' % (os.path.basename(os.path.dirname(work)), r[0], r[1]/1000, r[3], score, val))
        if best is None or val > best[0]: best = (val, int(r[0]))
        if score >= 45 or (r[1] >= 1000000 and score >= 25): good.append(int(r[0]))    # grande cavité de la tête à paroi un peu pâle (CS16)
    ax = json.load(open(os.path.join(work, 'axes.json'))); ax['vent_lab'] = best[1]; ax['vent_labs'] = sorted(set(good + [best[1]])); json.dump(ax, open(os.path.join(work, 'axes.json'), 'w'), indent=1)
    print('  -> vent_lab', best[1], 'vent_labs', ax['vent_labs'])
    dens = np.load(os.path.join(work, 'dens.npy'), mmap_mode='r'); LR = dens.shape[0]; i = LR//2
    img = cv2.cvtColor(255-np.asarray(dens[i]), cv2.COLOR_GRAY2BGR); img[np.isin(np.asarray(cavo[i]), ax['vent_labs'])] = (0, 200, 255)
    h = 420; img = cv2.resize(img, (int(img.shape[1]*h/img.shape[0]), h)); cv2.putText(img, os.path.basename(os.path.dirname(work)), (10,25), cv2.FONT_HERSHEY_SIMPLEX, 0.8, (0,0,0), 2); tiles.append(img)
w = max(t.shape[1] for t in tiles); tiles = [cv2.copyMakeBorder(t,0,0,0,w-t.shape[1],cv2.BORDER_CONSTANT,value=(128,128,128)) for t in tiles]
while len(tiles) % 4: tiles.append(np.full_like(tiles[0], 128))
cv2.imwrite('check_vent.png', np.vstack([np.hstack(tiles[i:i+4]) for i in range(0, len(tiles), 4)])); print('montage check_vent.png')
