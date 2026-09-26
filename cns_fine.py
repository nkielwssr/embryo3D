"""Segmentation fine du système nerveux central à partir des labels 'snc' et 'ventricules'.
Le système ventriculaire est paramétré par la distance géodésique depuis son extrémité antérieure (télencéphale) ; le profil
de section le long de cet axe présente des étranglements (diencéphale/mésencéphale, isthme, entrée du canal central) qui servent de coupures :
  cavités : 'ventricule_prosencephale', 'ventricule_mesencephale', 'ventricule_rhombencephale', 'canal_central'
  parois  : 'prosencephale', 'mesencephale', 'rhombencephale' (plus proche cavité), 'moelle' (tube sous la tête)
  enveloppes : 'meninges_mesenchyme_cranien' (tissu entre paroi neurale et épiderme, hors organes), 'epiderme_cranien'
usage : python cns_fine.py <dossier_stade> [--diag] [--cuts=a,b,c] (coupures en voxels géodésiques, sinon automatiques)
"""
import numpy as np, os, sys, re, shutil, subprocess, cv2

class _Npz(dict):
    """contenu d'un .npz chargé en mémoire et FERMÉ (np.load garde le fichier ouvert, ce qui bloque os.replace sous Windows)"""
    files = property(lambda self: list(self.keys()))
def _load_npz(path):
    import numpy as _np
    with _np.load(path) as f: return _Npz({k: f[k] for k in f.files})


def _savez_atomic(path, **kw):
    """écriture atomique (fichier temporaire puis remplacement) : d'autres sessions lisent labels.npz en parallèle"""
    import numpy as _np, os as _os; tmp = path + '.tmp.npz'; _np.savez_compressed(tmp, **kw); _os.replace(tmp, path)

from scipy import ndimage as ndi
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from segment import dilate_r, erode_r, open_r, close_r, comps_info, _crop, _dist, geodesic
stage_dir = os.path.abspath(sys.argv[1]); diag = '--diag' in sys.argv
CUTS = next((a.split('=')[1] for a in sys.argv if a.startswith('--cuts=')), None)
work = os.path.join(stage_dir, 'work'); out = os.path.join(stage_dir, 'out'); stage = 'CS%s' % re.search(r'CS\s*(\d+)', os.path.basename(stage_dir), re.I).group(1)
z = _load_npz(os.path.join(work, 'labels.npz')); shape = tuple(z['shape']); n = int(np.prod(shape))
L = {k: np.unpackbits(z[k])[:n].reshape(shape).astype(bool) for k in z.files if k != 'shape'}
ds = np.load(os.path.join(work, 'ds1.npy')).astype(np.float32); env = L['enveloppe']; snc = L['snc']; vent = L['ventricules']
lr = np.where(env.any(axis=(1,2)))[0]; MID = 0.5*(lr.min()+lr.max()); LRe = lr.max()-lr.min()
si = np.where(env.any(axis=(0,2)))[0]; si0, si1 = si.min(), si.max(); SIe = si1-si0
F = 3                                                       # facteur de sous-échantillonnage pour le graphe
# ------------------------------------------------------------------ 1. distance géodésique dans les ventricules (grille /F)
# domaine du graphe = lumière + paroi neurale (la paroi est continue là où la lumière est pincée ou mal détectée)
dom_full = vent | snc
dom_full[:, int(si0 + 0.50*SIe):, :] = False                    # tête + cou seulement : la queue recourbée ne perturbe pas l'axe
dom = ndi.maximum_filter(dom_full, size=F)[::F, ::F, ::F]
lab, k = ndi.label(dom); dom = lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)
vsub = dom; vent_sub = ndi.maximum_filter(vent, size=F)[::F, ::F, ::F] & dom
idx = np.flatnonzero(vsub.ravel()); pos = np.array(np.unravel_index(idx, vsub.shape)).T
lut = -np.ones(vsub.size, np.int64); lut[idx] = np.arange(len(idx))
rows, cols, wts = [], [], []
import itertools
for off in itertools.product((-1, 0, 1), repeat=3):                # 26-connexité, poids euclidiens (évite la métrique de Manhattan)
    if off == (0, 0, 0) or off < (0, 0, 0): continue
    nb = pos + np.array(off)
    ok = np.all((nb >= 0) & (nb < np.array(vsub.shape)), axis=1)
    j = lut[np.ravel_multi_index(nb[ok].T, vsub.shape)]; good = j >= 0
    rows.append(np.flatnonzero(ok)[good]); cols.append(j[good]); wts.append(np.full(good.sum(), float(np.sqrt(sum(o*o for o in off)))))
rows = np.concatenate(rows); cols = np.concatenate(cols); wts = np.concatenate(wts)
G = coo_matrix((wts, (rows, cols)), shape=(len(idx), len(idx))).tocsr(); G = G + G.T
# extrémité antérieure : voxel de la cavité le plus ventral (AP min) dans la moitié haute ; départ = son plus lointain géodésique n'est pas nécessaire
in_vent = vent_sub.ravel()[idx]
head_mask = (pos[:, 1] < (si0 + 0.45*SIe) / F) & in_vent
# extrémités de l'axe neural = diamètre géodésique du graphe (2 passes) ; rostrale = celle entourée du plus de lumière ventriculaire
def _far(src):
    dd = dijkstra(G, directed=False, indices=src); dd[~np.isfinite(dd)] = -1; return dd
# extrémité rostrale (télencéphale) = voxel ventriculaire le plus antérieur (AP min) dans la zone de recherche --tipzone=a,b
# (fractions de la hauteur du corps ; défaut 0..0.45 = toute la tête ; 0.2..0.45 quand le toit du mésencéphale déborde en avant : CS15, CS16)
TZ = [float(x) for x in next((a.split('=')[1] for a in sys.argv if a.startswith('--tipzone=')), '0,0.45').split(',')]
zone = in_vent & (pos[:, 1] >= (si0 + TZ[0]*SIe) / F) & (pos[:, 1] < (si0 + TZ[1]*SIe) / F)
if '--tip=largest' in sys.argv:                                   # tip = point le plus antérieur de la PLUS GRANDE cavité ventriculaire (CS15 : télencéphale ventral)
    labv, kv = ndi.label(vent_sub); big = labv == (np.argmax(np.bincount(labv.ravel())[1:]) + 1); zone = big.ravel()[idx] & in_vent
cand = np.flatnonzero(zone) if zone.any() else np.flatnonzero(in_vent)
if '--tip=bottom' in sys.argv: tip = int(cand[np.argmax(pos[cand, 1] * 1.0 - 0.3 * pos[cand, 2])])   # pôle ventral du prosencéphale (le plus bas dans la tête)
else: tip = int(cand[np.argmin(pos[cand, 2] * 1.0 + 0.3 * pos[cand, 1])])
dA = _far(tip); A = tip; B = tip; dB = dA
print('  tip rostral (SI %d, AP %d) zone %s' % (pos[tip, 1]*F, pos[tip, 2]*F, TZ))
dist = (dA if tip == A else dB).copy(); dist[dist < 0] = dist.max()               # en pas de F voxels
dmap_sub = np.full(vsub.shape, -1.0); dmap_sub.ravel()[idx] = dist
dmap = np.repeat(np.repeat(np.repeat(dmap_sub, F, 0), F, 1), F, 2)[:shape[0], :shape[1], :shape[2]]
pad = [(0, shape[i] - dmap.shape[i]) for i in range(3)]; dmap = np.pad(dmap, pad, constant_values=-1)
# voxels de vent hors composante principale : distance du plus proche voxel connu
dmap[~vent] = -1
# ------------------------------------------------------------------ 2. profil de section et coupures
dv = dmap[vent]; dv = dv[dv >= 0]; nb_bins = int(dv.max()) + 1
prof = np.bincount(dv.astype(int), minlength=nb_bins).astype(float) * 1.0          # voxels par pas géodésique (= section x F)
sm = ndi.uniform_filter1d(prof, 7)
MODE = next((a.split('=')[1] for a in sys.argv if a.startswith('--mode=')), 'apex')
if CUTS:
    cuts = [int(x) for x in CUTS.split(',')]
elif MODE == 'apex':
    # repère anatomique : le mésencéphale est au sommet de la flexure céphalique = point le plus HAUT (SI min) de l'axe ventriculaire.
    tail = np.median(sm[int(0.92 * nb_bins):]) if nb_bins > 30 else 0
    above = np.flatnonzero(sm > 4 * tail + 1e-6); end_brain = int(above[-1]) + 1 if len(above) else nb_bins
    # axe médian (SI, AP) par pas géodésique -> courbure ; l'apex de la flexure céphalique = max de courbure dans la 1re moitié
    dvi = dmap[vent]; ok = dvi >= 0; dbin = dvi[ok].astype(int); wv = np.where(vent); si_v = wv[1][ok]; ap_v = wv[2][ok]
    cnt = np.maximum(np.bincount(dbin, minlength=nb_bins), 1)
    si_mean = ndi.uniform_filter1d(np.bincount(dbin, weights=si_v, minlength=nb_bins) / cnt, 9)
    ap_mean = ndi.uniform_filter1d(np.bincount(dbin, weights=ap_v, minlength=nb_bins) / cnt, 9)
    step = max(4, end_brain // 12)
    ang = np.zeros(nb_bins)
    for b in range(step, nb_bins - step):
        v1 = np.array([si_mean[b] - si_mean[b - step], ap_mean[b] - ap_mean[b - step]]); v2 = np.array([si_mean[b + step] - si_mean[b], ap_mean[b + step] - ap_mean[b]])
        c = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-9); ang[b] = np.degrees(np.arccos(np.clip(c, -1, 1)))
    lo, hi = int(0.12 * end_brain), int(0.65 * end_brain)
    apex = lo + int(np.argmax(ndi.uniform_filter1d(ang, 5)[lo:hi]))
    print('  courbure max %.0f deg au pas %d ; SI min au pas %d' % (ang[apex], apex, lo + int(np.argmin(si_mean[lo:hi]))))
    w = max(6, int(0.09 * end_brain))
    cuts = [max(5, apex - w), min(end_brain - 5, apex + w), end_brain]
    print('  apex (mesencephale) au pas %d, demi-largeur %d, fin encephale %d' % (apex, w, end_brain))
else:
    # minima locaux de la section (prominence) avant la chute vers le canal central
    from scipy.signal import find_peaks
    inv = sm.max() - sm
    peaks, props = find_peaks(inv, prominence=0.12 * sm.max(), distance=max(5, nb_bins // 25))
    # niveau du canal : section médiane du dernier tiers ; fin de la partie encéphalique = premier pas où sm passe sous 2x ce niveau durablement
    tail = np.median(sm[int(0.92 * nb_bins):]) if nb_bins > 30 else 0          # niveau du canal central (queue du profil)
    above = np.flatnonzero(sm > 4 * tail + 1e-6)
    end_brain = int(above[-1]) + 1 if len(above) else nb_bins                      # dernier pas encore « encéphalique »
    mins = [int(p) for p in peaks if p < end_brain]
    mins = sorted(sorted(mins, key=lambda p: -props['prominences'][list(peaks).index(p)])[:2])
    cuts = mins + [end_brain]
print('%s : ventricules %.2fM vox, axe %d pas (x%d vox), coupures (pas) %s' % (stage, vent.sum()/1e6, nb_bins, F, cuts))
np.save(os.path.join(work, 'cns_profile.npy'), sm)
bounds = [-1] + cuts + [10**9]
names = ['ventricule_prosencephale', 'ventricule_mesencephale', 'ventricule_rhombencephale', 'canal_central'][:len(bounds) - 1]
full = np.zeros(shape, np.int32)
for k, nm in enumerate(names):
    full[(dmap > bounds[k]) & (dmap <= bounds[k + 1]) & vent] = k + 1
# petites composantes ventriculaires hors graphe -> plus proche
if ((full == 0) & vent).any():
    ind = ndi.distance_transform_edt(full == 0, return_distances=False, return_indices=True); full = np.where(vent, full[tuple(ind)], 0)
cav = {nm: full == (k + 1) for k, nm in enumerate(names)}
for nm, m in cav.items(): print('  %-28s %8d vox' % (nm, m.sum()))
# ------------------------------------------------------------------ 3. parois : plus proche cavité ; moelle sous la tête
ind = ndi.distance_transform_edt(full == 0, return_distances=False, return_indices=True)
wall = full[tuple(ind)] * snc
walls = {'prosencephale': wall == 1, 'mesencephale': wall == 2, 'rhombencephale': wall == 3, 'moelle': (wall == 4) if len(names) > 3 else np.zeros(shape, bool)}
head_lim = int(si0 + 0.42*SIe)
below = np.zeros(shape, bool); below[:, head_lim:, :] = True
walls['moelle'] |= snc & below & (np.abs(np.arange(shape[0]) - MID)[:, None, None] <= 0.13*LRe)
for k in ('prosencephale', 'mesencephale', 'rhombencephale'): walls[k] &= ~walls['moelle']
for nm, m in walls.items(): print('  %-28s %8d vox' % (nm, m.sum()))
# ------------------------------------------------------------------ 4. méninges / mésenchyme crânien et épiderme
head = np.zeros(shape, bool); head[:, :head_lim, :] = True
skin = env & ~erode_r(env, 3)
excl = snc | vent
for k in ('yeux', 'cristallins', 'vesicules_otiques', 'chondrocrane', 'coeur', 'coeur_detoure', 'foie', 'cavite_pericardique'):
    if k in L: excl |= L[k]
mes = open_r(head & env & ~skin & ~excl & (ds > 8), 1)
lab2, info2 = comps_info(mes, 20000); mes = np.isin(lab2, [i for i, s, c in info2])
epid = head & skin & ~excl
new = {**cav, **walls, 'meninges_mesenchyme_cranien': mes, 'epiderme_cranien': epid}
print('  meninges/mesenchyme %.2fM, epiderme cranien %.2fM' % (mes.sum()/1e6, epid.sum()/1e6))
# ------------------------------------------------------------------ planche (2 sagittales + 1 frontale) + profil
d = np.load(os.path.join(work, 'dens.npy'), mmap_mode='r'); i = int(MID); j = int(MID + 0.12*LRe)
COL = {'ventricule_prosencephale': (255, 200, 80), 'ventricule_mesencephale': (255, 120, 200), 'ventricule_rhombencephale': (120, 220, 255), 'canal_central': (200, 255, 200),
       'prosencephale': (200, 120, 0), 'mesencephale': (180, 0, 120), 'rhombencephale': (0, 120, 200), 'moelle': (0, 160, 0),
       'meninges_mesenchyme_cranien': (235, 235, 210), 'epiderme_cranien': (150, 150, 150)}
def paint(sl, transpose=False):
    img = cv2.cvtColor(255 - np.asarray(d[sl]), cv2.COLOR_GRAY2BGR)
    if transpose: img = np.ascontiguousarray(img.transpose(1, 0, 2))
    for nm, col in COL.items():
        if nm in new:
            m = new[nm][sl]; m = m.T if transpose else m; img[m] = col[::-1]
    return img
ap = int(np.mean(np.where(vent.any(axis=(0, 1)))[0]))
plot = np.full((shape[1], 400, 3), 255, np.uint8); pm = sm / max(sm.max(), 1)
for x in range(min(nb_bins, 400)):
    cv2.line(plot, (x, shape[1] - 1), (x, shape[1] - 1 - int(pm[x] * (shape[1] - 20))), (80, 80, 80), 1)
for c in cuts:
    if c < 400: cv2.line(plot, (c, 0), (c, shape[1] - 1), (0, 0, 255), 1)
cv2.putText(plot, 'section vs distance geodesique', (5, 20), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 0, 0), 1)
dm = np.asarray(dmap[i]); dimg = cv2.cvtColor(255 - np.asarray(d[i]), cv2.COLOR_GRAY2BGR); mm = dm >= 0
if mm.any(): dimg[mm] = cv2.applyColorMap(np.clip(dm[mm] / max(cuts[-1] * 1.15, 1) * 255, 0, 255).astype(np.uint8), cv2.COLORMAP_JET)[:, 0, :]
cv2.putText(dimg, 'distance geodesique (bleu=tip)', (5, 20), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 0, 0), 1)
cv2.imwrite(os.path.join(work, 'check_cns_fine.png'), np.hstack([paint(i), dimg, paint((slice(None), slice(None), ap), True), plot]))
if not diag:
    for nm, m in new.items(): L[nm] = m
    shutil.copy(os.path.join(work, 'labels.npz'), os.path.join(work, 'labels_avant_cns_fine.npz'))
    _savez_atomic(os.path.join(work, 'labels.npz'), **{kk: np.packbits(v) for kk, v in L.items()}, shape=np.array(shape))
    import meshexport, export_nrrd; meshexport.run(work, out, stage); export_nrrd.run(stage_dir)
print('done', stage)
