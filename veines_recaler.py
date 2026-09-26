"""Recale des points de veines approximatifs (transfert EHD) sur les lumières pâles du volume.
usage : python veines_recaler.py <work> <transfert.json> <veine> [rayon=25] [seuil=25]
Affiche, pour chaque point, la lumière retenue (x, s, y, taille, distance) ou 'aucune'."""
import numpy as np, json, sys, os
from scipy import ndimage as ndi
w, f, v = sys.argv[1], sys.argv[2], sys.argv[3]
R = float(sys.argv[4]) if len(sys.argv) > 4 else 25; thr = float(sys.argv[5]) if len(sys.argv) > 5 else 25
d = np.load(os.path.join(w, 'dens.npy'), mmap_mode='r'); SH = d.shape
L = np.load(os.path.join(w, 'labels.npz')); n = int(np.prod(SH))
env = np.unpackbits(L['enveloppe'])[:n].reshape(SH).astype(bool)
J = json.load(open(f, encoding='utf-8'))['veines'][v]
if isinstance(J, list):                     # format v4 : liste de points avec drapeau 'fiable'
    PTS = [((q['x'], q['s'], q['y']), q['section']) for q in J if q.get('fiable')]
else:
    PTS = list(zip(J['points_xsy'], J['sections_EHD']))
res = []
for (x, s, y), sec in PTS:
    s = int(round(s))
    if sec < 55 or not (0 <= s < SH[1]): continue
    sl = ndi.gaussian_filter(np.asarray(d[:, s, :]).astype(float), 1.0)
    e = ndi.binary_erosion(env[:, s, :], iterations=6)
    lab, k = ndi.label((sl < thr) & e); best = None
    for i in range(1, k + 1):
        m = lab == i; sz = int(m.sum())
        if not (12 <= sz <= 900): continue
        cx, cy = np.argwhere(m).mean(0); dist = np.hypot(cx - x, cy - y)
        if dist <= R and (best is None or dist < best[4]): best = (int(cx), s, int(cy), sz, round(float(dist), 1))
    res.append((sec, best))
    print(sec, s, (round(x), round(y)), '->', best if best else 'aucune')
