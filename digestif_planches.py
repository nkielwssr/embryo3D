"""Planches de coupes pour placer les points de passage du tube digestif.
usage : python digestif_planches.py <work_dir> <axe 0|1|2> <debut> <fin> <pas> [sortie.png]
Axes du volume dens.npy : 0 = gauche-droite, 1 = haut-bas (0 = haut), 2 = ventral-dorsal (0 = ventral).
Chaque vignette porte son indice de coupe ; une graduation tous les 50 voxels donne les coordonnées.
Contours : foie (rouge), tube digestif actuel (magenta), coeur (orange)."""
import numpy as np, cv2, sys, os

work, ax, a0, a1, step = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
out = sys.argv[6] if len(sys.argv) > 6 else os.path.join(work, f'digestif_ax{ax}_{a0}_{a1}.png')
d = np.load(os.path.join(work, 'dens.npy'), mmap_mode='r')
L = np.load(os.path.join(work, 'labels.npz'))
cols = {'foie': (0, 0, 255), 'tube_digestif': (255, 0, 255), 'coeur': (0, 140, 255)}
shape = tuple(L['shape'])
masks = {k: np.unpackbits(L[k])[:int(np.prod(shape))].reshape(shape) for k in cols if k in L.files}

def take(a, i):
    s = np.take(a, i, axis=ax)
    # affichage : ax0 -> (haut-bas en lignes, ventral-dorsal en colonnes) ; ax1 -> (ventral-dorsal en lignes, gauche-droite en colonnes)
    # ax2 -> (haut-bas en lignes, gauche-droite en colonnes)
    return s.T if ax == 1 else (s if ax == 0 else s.T.T)

tiles = []
for i in range(a0, a1 + 1, step):
    s = 255 - np.asarray(take(d, i))
    if ax == 2: s = s  # (LR, SI) -> transposer pour haut-bas en lignes
    if ax == 2: s = s.T
    im = cv2.cvtColor(np.ascontiguousarray(s.astype(np.uint8)), cv2.COLOR_GRAY2BGR)
    for k, m in masks.items():
        mm = np.asarray(take(m, i)) > 0
        if ax == 2: mm = mm.T
        cnt, _ = cv2.findContours(np.ascontiguousarray(mm.astype(np.uint8)), cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)
        cv2.drawContours(im, cnt, -1, cols[k], 1)
    h, w = im.shape[:2]
    for y in range(0, h, 50):
        cv2.line(im, (0, y), (6, y), (0, 160, 0), 1); cv2.putText(im, str(y), (8, y + 4), 0, 0.33, (0, 130, 0), 1)
    for x in range(0, w, 50):
        cv2.line(im, (x, 0), (x, 6), (0, 160, 0), 1); cv2.putText(im, str(x), (x + 1, 16), 0, 0.33, (0, 130, 0), 1)
    cv2.putText(im, f'{ax}:{i}', (w - 70, h - 8), 0, 0.6, (200, 0, 0), 2)
    tiles.append(im)
H = max(t.shape[0] for t in tiles); W = max(t.shape[1] for t in tiles)
tiles = [cv2.copyMakeBorder(t, 0, H - t.shape[0], 0, W - t.shape[1] + 4, cv2.BORDER_CONSTANT, value=(90, 90, 90)) for t in tiles]
n = len(tiles); c = int(np.ceil(np.sqrt(n * H / W * 1.3)))
rows = [np.hstack(tiles[k:k + c] + [np.full_like(tiles[0], 90)] * (c - len(tiles[k:k + c]))) for k in range(0, n, c)]
cv2.imwrite(out, np.vstack(rows)); print(out)
