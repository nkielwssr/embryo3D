"""Coupes transversales zoomées avec contour du tube digestif reconstruit (contrôle).
usage : python digestif_verif.py <work> s0 s1 pas y0 y1 x0 x1 zoom sortie.png [colonnes]"""
import numpy as np, cv2, sys, os
work = sys.argv[1]; s0, s1, st, r0, r1, c0, c1 = map(int, sys.argv[2:9]); zoom = float(sys.argv[9]); out = sys.argv[10]
ncol = int(sys.argv[11]) if len(sys.argv) > 11 else 5
d = np.load(os.path.join(work, 'dens.npy'), mmap_mode='r')
z = np.load(os.path.join(work, sys.argv[12] if len(sys.argv) > 12 else 'digestif/digestif.npz')); SH = tuple(z['shape'])
segs = {k: np.unpackbits(z[k])[:int(np.prod(SH))].reshape(SH) for k in z.files if k not in ('shape', 'union')}
cols = [(255, 0, 255), (0, 160, 255), (0, 200, 0), (255, 120, 0), (0, 0, 255), (200, 200, 0), (120, 0, 200)]
tiles = []
for s in range(s0, s1 + 1, st):
    im = 255 - np.ascontiguousarray(d[:, s, :].T[r0:r1, c0:c1])
    im = cv2.cvtColor(cv2.resize(im, None, fx=zoom, fy=zoom, interpolation=cv2.INTER_CUBIC), cv2.COLOR_GRAY2BGR)
    for k, (nom, m) in enumerate(segs.items()):
        mm = np.ascontiguousarray(m[:, s, :].T[r0:r1, c0:c1]).astype(np.uint8)
        mm = cv2.resize(mm, (im.shape[1], im.shape[0]), interpolation=cv2.INTER_NEAREST)
        cnt, _ = cv2.findContours(mm, cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE); cv2.drawContours(im, cnt, -1, cols[k % 7], 1)
    for y in range((r0 // 20 + 1) * 20, r1, 20):
        yy = int((y - r0) * zoom); cv2.line(im, (0, yy), (5, yy), (0, 150, 0), 1)
        if y % 40 == 0: cv2.putText(im, str(y), (8, yy + 4), 0, 0.35, (0, 120, 0), 1)
    for x in range((c0 // 20 + 1) * 20, c1, 20):
        xx = int((x - c0) * zoom); cv2.line(im, (xx, 0), (xx, 5), (0, 150, 0), 1)
        if x % 40 == 0: cv2.putText(im, str(x), (xx + 2, 18), 0, 0.35, (0, 120, 0), 1)
    cv2.putText(im, str(s), (im.shape[1] - 45, im.shape[0] - 8), 0, 0.6, (200, 0, 0), 2)
    tiles.append(cv2.copyMakeBorder(im, 0, 3, 0, 3, cv2.BORDER_CONSTANT, value=(90, 90, 90)))
rows = [np.hstack(tiles[k:k + ncol] + [np.full_like(tiles[0], 90)] * (ncol - len(tiles[k:k + ncol]))) for k in range(0, len(tiles), ncol)]
cv2.imwrite(out, np.vstack(rows)); print(out)
