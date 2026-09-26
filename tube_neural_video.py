"""
Reperage du tube neural (moelle) et du rhombencephale sur UNE image de la video de rotation 360 deg (rendu volumique
semi-transparent), sans passer par les labels du pipeline (le recalage 2D des labels sur une vue oblique est inutilisable).
Methode (vue de profil, embryon en C, tube neural sous le contour dorsal) :
  - bord dorsal du tube = contour de la silhouette rentre de 3 px (ectoderme), de la pointe de la queue a la jonction cervicale ;
  - bord ventral = queue : 45 % de l epaisseur locale (tube = moitie dorsale d une queue fine) ; tronc : les stries somitiques
    (periode 12-40 px le long du contour) commencent a d_debut et finissent a d_fin sous la peau, tube = d_debut + 0,55 (d_fin - d_debut)
    (le tube neural occupe ~ les 2/3 dorsaux de la hauteur des somites, la notochorde a son pole ventral) ;
  - rhombencephale : toit = contour de la jonction cervicale a l isthme, plancher = polyligne donnee (estimee : dorsal a la
    vesicule otique et au ganglion trijumeau) ; l isthme est ESTIME (limite mesencephale).
Les reperes (pointe de queue, jonction cervicale, isthme, plancher) sont propres a chaque stade et a chaque image : les lire sur
l image (grille) avant de lancer. Sorties : out/topographie/video360/tube_neural_video_<CS>.{png,json} (px de l image 736x736).
Usage : python embryo3d/tube_neural_video.py CS14_f4v --image 125 --queue 455,500,355.5 --cou 252,205 --isthme 415,93           --plancher "415,150;395,175;370,190;345,200;320,208;300,215;278,225" --mid 505,500
  --queue xmin,xmax,ymin : fenetre ou la pointe de la queue est le point le plus haut du contour ; --mid : point du bord dorsal
  de la queue (fixe le sens de parcours). Sans work/rotation/frames.npy, la video est relue (video360.charger_images).
"""
import sys, os, argparse, numpy as np, cv2, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from video360 import silhouette, trouver_video, charger_images
from scipy import ndimage as ndi

ap = argparse.ArgumentParser(); ap.add_argument('dossier'); ap.add_argument('--image', type=int, required=True)
ap.add_argument('--queue', required=True); ap.add_argument('--cou', required=True); ap.add_argument('--isthme', required=True)
ap.add_argument('--plancher', required=True); ap.add_argument('--mid', required=True)
ap.add_argument('--crop', default='', help='x0,y0,x1,y1 (px image) pour une vignette recadree')
args = ap.parse_args()
def xy(t): return tuple(float(v) for v in t.split(','))
dossier = args.dossier; K = args.image
man = json.load(open(os.path.join(dossier, 'out', 'manifest.json'), encoding='utf-8')); CS = man['stage']
SC = os.path.join(dossier, 'out', 'topographie', 'video360'); os.makedirs(SC, exist_ok=True)
fr = os.path.join(dossier, 'work', 'rotation', 'frames.npy')
if os.path.exists(fr):
    F = np.load(fr, mmap_mode='r')
else:
    v = trouver_video(dossier); F = charger_images(v)
f = np.array(F[K]); ff = f.astype(float); H, W = f.shape
sil = silhouette(f)
cs, _ = cv2.findContours(sil.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
c = max(cs, key=cv2.contourArea)[:, 0, :].astype(float)
c = np.stack([ndi.gaussian_filter1d(c[:, 0], 6, mode='wrap'), ndi.gaussian_filter1d(c[:, 1], 6, mode='wrap')], 1); n = len(c)
def near(p, cond=None):
    d = np.linalg.norm(c - np.array(p), axis=1)
    if cond is not None: d[~cond] = 1e9
    return int(np.argmin(d))
# pointe de la queue : point du contour le plus haut dans la fenetre 455<x<500, y>352
qx0, qx1, qy = xy(args.queue); cand = (c[:, 0] > qx0) & (c[:, 0] < qx1) & (c[:, 1] > qy); i_tip = int(np.argmin(np.where(cand, c[:, 1], 1e9)))
i_neck = near(xy(args.cou)); i_mid = near(xy(args.mid)); i_isth = near(xy(args.isthme))
print('pointe queue', c[i_tip].round(0), 'cou', c[i_neck].round(0), 'isthme', c[i_isth].round(0))
sens = 1 if (i_mid - i_tip) % n < (i_tip - i_mid) % n else -1
L = (i_neck - i_tip) * sens % n; idx = (i_tip + sens * np.arange(L + 1)) % n; P = c[idx]
L2 = (i_isth - i_neck) * sens % n; idx2 = (i_neck + sens * np.arange(L2 + 1)) % n; P2 = c[idx2]     # toit du rhombencephale
s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
def normale(P):
    t = np.gradient(P, axis=0); t /= np.linalg.norm(t, axis=1, keepdims=True) + 1e-9
    nr = np.stack([-t[:, 1], t[:, 0]], 1); q = P + 6 * nr
    ins = sil[np.clip(q[:, 1].astype(int), 0, H - 1), np.clip(q[:, 0].astype(int), 0, W - 1)]
    return nr if ins.mean() >= 0.5 else -nr
nrm = normale(P); nrm2 = normale(P2)
depths = np.arange(0, 60); I = np.zeros((60, len(P)))
for j, d in enumerate(depths):
    q = P + d * nrm[:, :]; I[j] = ndi.map_coordinates(ff, [q[:, 1], q[:, 0]], order=1, mode='nearest')
nP = len(P)
ep = np.full(nP, 59.0)
for i in range(nP):
    k = np.nonzero(I[5:, i] < 18)[0]
    if len(k): ep[i] = 5 + k[0]
ep = ndi.median_filter(ep, 41)
E = np.zeros_like(I)
for j in range(60):
    bp = ndi.gaussian_filter1d(I[j], 3) - ndi.gaussian_filter1d(I[j], 12); E[j] = np.sqrt(ndi.uniform_filter1d(bp ** 2, 120))
d_start = np.zeros(nP); d_end = np.zeros(nP)
for i in range(nP):
    col = E[:, i]; on = np.nonzero(col[3:] > 0.4 * col[3:].max())[0] + 3
    d_start[i], d_end[i] = (on[0], on[-1]) if len(on) else (12, 50)
d_start = ndi.median_filter(d_start, 81); d_end = ndi.median_filter(d_end, 81)
T = np.where(ep < 58, 0.45 * ep, d_start + 0.55 * (d_end - d_start)); T = ndi.gaussian_filter1d(np.clip(T, 6, 40), 25)
for a in range(0, nP, 150):
    sel = slice(a, min(a + 150, nP)); print(f's {a}-{a+150}: epaisseur {ep[sel].mean():.0f}  stries {d_start[sel].mean():.0f}-{d_end[sel].mean():.0f}  T {T[sel].mean():.0f}')
d0 = 3.0
dors = P + d0 * nrm; vent = P + T[:, None] * nrm
# rhombencephale : toit = contour (offset 3) du cou a l'isthme ; plancher = polyligne estimee (dorsal au ganglion trijumeau et a la vesicule otique)
toit = P2 + d0 * nrm2
plancher = np.array([xy(t) for t in args.plancher.split(';')], float)
plancher[-1] = vent[-1]                     # raccord exact au bord ventral du tube
rh = np.vstack([toit, plancher])
# ---------- figure
sc = 2; big = cv2.resize(cv2.cvtColor(f, cv2.COLOR_GRAY2BGR), None, fx=sc, fy=sc, interpolation=cv2.INTER_CUBIC)
over = big.copy()
cv2.fillPoly(over, [(np.vstack([dors, vent[::-1]]) * sc).astype(np.int32)], (60, 220, 60))
cv2.fillPoly(over, [(rh * sc).astype(np.int32)], (0, 170, 255))
big = cv2.addWeighted(over, 0.32, big, 0.68, 0)
cv2.polylines(big, [(dors * sc).astype(np.int32)], False, (60, 255, 60), 2); cv2.polylines(big, [(vent * sc).astype(np.int32)], False, (60, 255, 60), 2)
cv2.polylines(big, [(toit * sc).astype(np.int32)], False, (0, 200, 255), 2)
# plancher en pointille (estime)
for i in range(len(plancher) - 1):
    a, b = plancher[i] * sc, plancher[i + 1] * sc; m = int(np.linalg.norm(b - a) / 8)
    for k in range(0, max(m, 1), 2):
        p1 = a + (b - a) * k / max(m, 1); p2 = a + (b - a) * min(k + 1, m) / max(m, 1); cv2.line(big, tuple(p1.astype(int)), tuple(p2.astype(int)), (0, 200, 255), 2)
# isthme et jonction cervicale
cv2.line(big, tuple((toit[-1] * sc).astype(int)), tuple((plancher[0] * sc).astype(int)), (0, 0, 255), 2)
cv2.line(big, tuple((dors[-1] * sc).astype(int)), tuple((vent[-1] * sc).astype(int)), (255, 255, 255), 2)
cv2.circle(big, tuple((P[0] * sc).astype(int)), 6, (60, 255, 60), 2)
def txt(t, p, col=(255, 255, 255)):
    cv2.putText(big, t, p, cv2.FONT_HERSHEY_SIMPLEX, 0.6, (0, 0, 0), 3); cv2.putText(big, t, p, cv2.FONT_HERSHEY_SIMPLEX, 0.6, col, 1)
txt('tube neural (moelle) - trace depuis l image', (40, 1380), (60, 255, 60))
txt('rhombencephale (toit sur le contour, plancher estime)', (40, 1410), (0, 200, 255))
txt('isthme (limite mesencephale) - estime', (40, 1440), (0, 0, 255))
txt('jonction cervicale', (40, 1470), (255, 255, 255))
txt(f'{CS} - video 360 image {K}', (40, 40))
txt('pointe de la queue', (int(P[0][0] * sc) + 10, int(P[0][1] * sc) - 10), (60, 255, 60))
txt('isthme', (int(toit[-1][0] * sc) + 10, int(toit[-1][1] * sc) - 6), (0, 0, 255))
cv2.imwrite(os.path.join(SC, f'tube_neural_video_{CS}.png'), big)
# cropped version like the user's image
if args.crop:
    x0, y0, x1, y1 = (int(v) for v in args.crop.split(',')); cv2.imwrite(os.path.join(SC, f'tube_neural_video_{CS}_crop.png'), big[y0 * sc:y1 * sc, x0 * sc:x1 * sc])
out = {'stage': CS, 'image': K, 'taille_image': [W, H], 'reperes': vars(args), 'note': 'coordonnees en px de l image 736x736 ; bord dorsal = contour de la silhouette a 3 px ; bord ventral = queue 45 % de l epaisseur locale, tronc d_debut + 0,55 (d_fin - d_debut) des stries somitiques',
       'tube_neural': {'bord_dorsal_px': dors.round(1).tolist(), 'bord_ventral_px': vent.round(1).tolist(), 's_px': s.round(1).tolist(), 'longueur_px': float(s[-1])},
       'rhombencephale': {'toit_px': toit.round(1).tolist(), 'plancher_estime_px': plancher.round(1).tolist(), 'isthme_px': [toit[-1].round(1).tolist(), plancher[0].tolist()]},
       'jonction_cervicale_px': [dors[-1].round(1).tolist(), vent[-1].round(1).tolist()]}
json.dump(out, open(os.path.join(SC, f'tube_neural_video_{CS}.json'), 'w', encoding='utf-8'), indent=1)
print('longueur du tube (queue -> jonction cervicale)', round(s[-1]), 'px ; toit rhombencephale', round(float(np.sum(np.linalg.norm(np.diff(toit, axis=0), axis=1)))), 'px')
