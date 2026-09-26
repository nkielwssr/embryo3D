"""Candidats de veines cardinales : lumières pâles de part et d'autre des aortes tracées (work/cardio/vaisseaux.npz), coupe par coupe.
usage : python cardinales_candidats.py <work> [pas=15] [seuil=25] [dy=35] [dx=80]"""
import numpy as np, sys, os
from scipy import ndimage as ndi
w = sys.argv[1]; st = int(sys.argv[2]) if len(sys.argv) > 2 else 15; thr = float(sys.argv[3]) if len(sys.argv) > 3 else 25
DY = int(sys.argv[4]) if len(sys.argv) > 4 else 35; DX = int(sys.argv[5]) if len(sys.argv) > 5 else 80
d = np.load(os.path.join(w, 'dens.npy'), mmap_mode='r')
z = np.load(os.path.join(w, 'cardio', 'vaisseaux.npz')); SH = tuple(z['shape'])
ao = np.zeros(SH, bool)
for k in z.files:
    if k.startswith('aorte'): ao |= np.unpackbits(z[k])[:int(np.prod(SH))].reshape(SH).astype(bool)
ss = np.where(ao.any(axis=(0, 2)))[0]
for s in range(ss.min(), ss.max() + 1, st):
    a = ao[:, s, :]
    if not a.any(): continue
    sl = ndi.gaussian_filter(np.asarray(d[:, s, :]).astype(float), 1.0)
    pale = (sl < thr) & ~ndi.binary_dilation(a, iterations=3)
    lab, n = ndi.label(pale); ax = np.argwhere(a); ya = ax[:, 1].mean(); xl, xr = ax[:, 0].min(), ax[:, 0].max()
    out = []
    for i in range(1, n + 1):
        m = lab == i; sz = int(m.sum())
        if not (15 <= sz <= 600): continue
        cx, cy = np.argwhere(m).mean(0)
        if abs(cy - ya) < DY and (xl - DX < cx < xl - 3 or xr + 3 < cx < xr + DX): out.append((int(cx), int(cy), sz))
    print(s, 'aortes x', xl, xr, 'y', int(ya), 'fond', int(np.median(sl[max(xl-20,0):xr+20, int(ya)-20:int(ya)+20])), out)
