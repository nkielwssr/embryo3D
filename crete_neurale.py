"""Crête neurale (le « quatrième feuillet ») : origine, trajets de migration et dérivés, déduits des labels existants.
Lecture seule de work/labels.npz (comme topographie_axe.py) ; option --fusion pour écrire les labels crete_* dans labels.npz.

Biologie représentée (détail dans le README, section « Crête neurale ») : les cellules se détachent du toit du tube neural
(transition épithélio-mésenchymateuse) puis migrent dans le mésenchyme vers leurs cibles, par niveau d'origine :
  crâniale  (diencéphale → rhombencéphale)  → massif facial et arcs pharyngiens (squelette et cartilages de la face)
  cardiaque (post-otique, somites 1-3)      → pôle artériel du cœur : septation du bulbe, arcs aortiques
  vagale    (somites 1-7) et sacrée (> S28) → paroi du tube digestif : système nerveux entérique
  tronc     (segmentaire)                   → ganglions rachidiens (voie ventro-latérale) puis chaîne sympathique le long de l'aorte dorsale
Structures produites (work/crete_neurale/crete_neurale.npz, PLY + manifest_crete_neurale.json dans out/crete_neurale/) :
  crete_origine         bande dorsale du tube neural d'où émigrent les cellules (surface du tube la plus proche de l'ectoderme)
  crete_craniale        trajets vers la face (3 bandes rostro-caudales par côté) + trajets vers les ganglions crâniens
  crete_cardiaque       trajets vers le pôle artériel du cœur (par les arcs pharyngiens)
  crete_vagale          trajets vers l'intestin antérieur
  crete_enterique       gaine autour du tube digestif, limitée au front de colonisation du stade
  crete_tronc           trajets segmentaires vers chaque ganglion rachidien
  crete_sympathique     prolongements des ganglions vers l'aorte dorsale
  crete_sacree          trajets vers l'intestin postérieur
  crete_bulbe_cardiaque zone de septation : cœur au contact de l'aorte commune
Les trajets sont des plus courts chemins pondérés (Dijkstra) dans le mésenchyme — hors tube neural, organes, membres ; densité élevée
pénalisée — rasterisés en tubes fins : ce sont des SCHÉMAS de routes, pas des cellules observées (confiance 'faible' ; origine,
gaine entérique et bulbe 'moyenne' car calés sur des labels existants). Tout ce qui manque (pas de ganglions, pas d'aorte…) est sauté.
usage : python embryo3d/crete_neurale.py <dossier_stade> [--F 4] [--rayon-mm 0.06] [--fusion] [--diag] [--sans-ply]
  --fusion : écrit aussi les labels crete_* dans work/labels.npz (sauvegarde labels_avant_crete.npz) puis relance meshexport + NRRD.
  Ne pas lancer en parallèle d'un autre script qui écrit labels.npz (voir README : course → labels perdus)."""
import numpy as np, os, sys, re, json, time, itertools, shutil, argparse
from scipy import ndimage as ndi
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from segment import dilate_r, erode_r, close_r, comps_info, _crop
from meshexport import mask_to_mesh, CRL_MM

ap = argparse.ArgumentParser(); ap.add_argument('dossier')
ap.add_argument('--F', type=int, default=4, help='sous-échantillonnage du graphe de migration (voxels)')
ap.add_argument('--rayon-mm', type=float, default=0.06, help='rayon des tubes de trajet')
ap.add_argument('--fusion', action='store_true'); ap.add_argument('--diag', action='store_true', help='accepté pour homogénéité ; la planche de contrôle est toujours écrite')
ap.add_argument('--sans-ply', action='store_true')
ap.add_argument('--confiance-trajets', default='faible', choices=['bonne', 'moyenne', 'faible'])
args = ap.parse_args()
T0 = time.time()
def log(*a): print('[%5.0fs] crete' % (time.time() - T0), *a); sys.stdout.flush()

stage_dir = os.path.abspath(args.dossier); work = os.path.join(stage_dir, 'work'); out = os.path.join(stage_dir, 'out')
stage = 'CS%s' % re.search(r'CS\s*(\d+)', os.path.basename(stage_dir), re.I).group(1); num = int(stage[2:])
F = max(2, args.F)
# ------------------------------------------------------------------ labels
with np.load(os.path.join(work, 'labels.npz')) as z:
    shape = tuple(int(v) for v in z['shape']); n = int(np.prod(shape))
    L = {k: np.unpackbits(z[k])[:n].reshape(shape).astype(bool) for k in z.files if k != 'shape'}
env = L['enveloppe']
def U(*keys):
    m = np.zeros(shape, bool)
    for k in keys:
        if k in L: m |= L[k]
    return m
NT = U('snc', 'moelle', 'moelle_rachidienne', 'prosencephale', 'mesencephale', 'rhombencephale', 'ventricules',
       'ventricule_prosencephale', 'ventricule_mesencephale', 'ventricule_rhombencephale', 'canal_central')
if not NT.any(): sys.exit('%s : aucun label du tube neural (snc / moelle / encéphale) — rien à faire' % stage)
sl = _crop(NT, 3); NTf = np.zeros(shape, bool); NTf[sl] = ndi.binary_fill_holes(close_r(NT[sl], 1))
coeur = U('coeur_detoure', 'coeur', 'myocarde', 'cavites_cardiaques')
SEGMENTS = ['digestif_oesophage', 'digestif_estomac', 'digestif_duodenum', 'digestif_intestin_moyen', 'digestif_intestin_posterieur']
gut_segs = {k: L[k] for k in SEGMENTS if k in L and L[k].any()}
gut = U(*gut_segs) if gut_segs else U('tube_digestif')
foie = U('foie'); yeux = U('yeux', 'cristallins'); otic = U('vesicules_otiques'); aorte = U('vaisseaux_aorte'); gangl = U('ganglions')
cartilage = U('chondrocrane', 'cartilage_autre', 'squelette_axial_cartilage', 'corps_vertebraux', 'arcs_neuraux', 'cotes', 'notochorde', 'somites')
organes = coeur | gut | foie | yeux | otic | aorte | gangl | cartilage | U('membres', 'membre_sup_gauche', 'membre_sup_droit', 'membre_inf_gauche', 'membre_inf_droit', 'cordon_ombilical', 'vaisseaux', 'vaisseaux_cardinales', 'vaisseaux_ombilicaux', 'vaisseaux_vitellins')
# échelle et centre (manifest du pipeline, sinon recalcul comme meshexport)
pman = os.path.join(out, 'manifest.json')
if os.path.exists(pman):
    man = json.load(open(pman, encoding='utf-8')); vox = float(man['mm_per_voxel']); center = np.array(man['center_voxel'])
else:
    si_ext = np.ptp(np.where(env.any(axis=(0, 2)))[0]) + 1; vox = CRL_MM.get(stage, 20.0) / si_ext; center = np.array(ndi.center_of_mass(env))
mm = lambda v: v * vox; vx = lambda d_mm: d_mm / vox
r_tube = max(2, int(round(vx(args.rayon_mm))))
lr = np.where(env.any(axis=(1, 2)))[0]; MID = 0.5 * (lr.min() + lr.max())
log('%s : grille %s, %.4f mm/vox, F=%d, tube r=%d vox ; labels : %s' % (stage, shape, vox, F, r_tube, ', '.join(sorted(L))))
# ------------------------------------------------------------------ grille sous-échantillonnée
def down(m):
    pad = [(0, (-s) % F) for s in m.shape]; mp = np.pad(m, pad)
    a, b, c = (s // F for s in mp.shape)
    return mp.reshape(a, F, b, F, c, F).any(axis=(1, 3, 5))
def up(ms):
    full = np.repeat(np.repeat(np.repeat(ms, F, 0), F, 1), F, 2); return full[:shape[0], :shape[1], :shape[2]]
sub_shape = down(env).shape
envs = down(env); NTs = down(NTf)
lab, k = ndi.label(NTs); NTs = lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)       # composante principale du tube
Cbody = np.array(ndi.center_of_mass(envs))
dens = None
for f in ('ds1.npy', 'dens.npy'):
    p = os.path.join(work, f)
    if os.path.exists(p):
        d = np.load(p, mmap_mode='r')
        if tuple(d.shape) == shape: dens = np.asarray(d[::F, ::F, ::F]).astype(np.float32); break
OFFS = [o for o in itertools.product((-1, 0, 1), repeat=3) if o != (0, 0, 0) and o > (0, 0, 0)]
def graphe(mask, cost=None):
    """graphe 26-connexe des voxels de mask (grille /F) ; poids = distance euclidienne × coût moyen des deux extrémités"""
    idx = np.flatnonzero(mask.ravel()); pos = np.array(np.unravel_index(idx, mask.shape)).T
    lut = -np.ones(mask.size, np.int64); lut[idx] = np.arange(len(idx))
    c = np.ones(len(idx), np.float32) if cost is None else cost.ravel()[idx]
    rows, cols, wts = [], [], []
    for off in OFFS:
        nb = pos + np.array(off); ok = np.all((nb >= 0) & (nb < np.array(mask.shape)), axis=1)
        j = lut[np.ravel_multi_index(nb[ok].T, mask.shape)]; good = j >= 0; i = np.flatnonzero(ok)[good]; j = j[good]
        rows.append(i); cols.append(j); wts.append(np.sqrt(sum(o * o for o in off)) * 0.5 * (c[i] + c[j]))
    rows = np.concatenate(rows); cols = np.concatenate(cols); wts = np.concatenate(wts).astype(np.float32)
    G = coo_matrix((wts, (rows, cols)), shape=(len(idx), len(idx))).tocsr(); G = G + G.T
    return G, idx, pos, lut
# ------------------------------------------------------------------ 1. abscisse curviligne s le long du tube neural (graphe du tube)
Gn, idx_n, pos_n, lut_n = graphe(NTs)
if any(k in L for k in ('ventricule_prosencephale', 'prosencephale')):
    ct = np.array(ndi.center_of_mass(U('ventricule_prosencephale', 'prosencephale'))) / F
    tip = int(np.argmin(np.linalg.norm(pos_n - ct, axis=1)))
else:                                                       # pôle ventro-rostral de la tête (comme cns_fine : AP min, pondéré par la hauteur)
    si = np.where(envs.any(axis=(0, 2)))[0]; head = pos_n[:, 1] < si.min() + 0.45 * np.ptp(si)
    cand = np.flatnonzero(head) if head.any() else np.arange(len(pos_n))
    tip = int(cand[np.argmin(pos_n[cand, 2] * 1.0 + 0.3 * pos_n[cand, 1])])
s_n = dijkstra(Gn, directed=False, indices=tip); fin = np.isfinite(s_n); s_n[~fin] = s_n[fin].max()
tree_n = cKDTree(pos_n)
def s_of(p_full):
    """s (voxels /F) du point du tube le plus proche d'un point en coordonnées pleines"""
    return float(s_n[tree_n.query(np.asarray(p_full, float) / F)[1]])
Ls = float(s_n.max())
s_otic = s_of(ndi.center_of_mass(otic)) if otic.any() else 0.22 * Ls
if not otic.any(): log('  pas de vésicules otiques : niveau otique estimé à 22 % de l arc du tube')
pas = max((Ls - s_otic) / 40.0, 1.0)                       # ≈ un somite (otique → bout de la queue ≈ 40 segments)
log('  tube neural : arc %.2f mm, niveau otique s=%.2f mm, pas somitique ≈ %.3f mm' % (mm(Ls * F), mm(s_otic * F), mm(pas * F)))
# ------------------------------------------------------------------ 2. origine : surface du tube la plus proche de l'ectoderme
dist_skin = ndi.distance_transform_edt(envs)
surf_s = NTs & ~ndi.binary_erosion(NTs)
keep_s = np.zeros(sub_shape, bool)
for kk in np.flatnonzero(NTs.any(axis=(0, 2))):
    l2, n2 = ndi.label(NTs[:, kk, :])
    for i in range(1, n2 + 1):
        comp = l2 == i; sf = comp & surf_s[:, kk, :]
        if not sf.any(): continue
        d = np.where(sf, dist_skin[:, kk, :], np.inf); tol = max(1.0, 0.2 * np.sqrt(comp.sum() / np.pi))
        # moitié dorsale de la section seulement : côté opposé au centre du corps (l'embryon en C enroule son ventre autour du cœur)
        cc = np.array(ndi.center_of_mass(comp)); vdir = np.array([Cbody[0], Cbody[2]]) - cc
        g0, g2 = np.ogrid[:comp.shape[0], :comp.shape[1]]; dorsal = (g0 - cc[0]) * vdir[0] + (g2 - cc[1]) * vdir[1] <= 0
        keep_s[:, kk, :] |= sf & dorsal & (d <= d.min() + tol)
s_keep = np.full(sub_shape, -1.0); s_keep.ravel()[idx_n] = s_n
keep_s &= (s_keep >= 0.25 * s_otic)                         # pas de crête sur le télencéphale (crête neurale antérieure exclue)
surface = NTf & ~erode_r(NTf, 1)
origine = dilate_r(surface & up(ndi.binary_dilation(keep_s)), 1) & (NTf | dilate_r(NTf, 1))
log('  origine : %d voxels sous-grille, %.3f mm3' % (keep_s.sum(), origine.sum() * vox ** 3))
# ------------------------------------------------------------------ 3. graphe du mésenchyme (hors tube neural, organes, membres)
Ms = envs & ~ndi.binary_dilation(NTs) & ~down(organes)
cost = None
if dens is not None:
    dm = dens[tuple(slice(0, s) for s in sub_shape)] if dens.shape != sub_shape else dens
    dm = np.pad(dm, [(0, a - b) for a, b in zip(sub_shape, dm.shape)]) if dm.shape != sub_shape else dm
    p90 = np.percentile(dm[Ms], 90) if Ms.any() else 1.0
    cost = (1.0 + 4.0 * np.clip(dm / max(p90, 1e-6), 0, 1.5) ** 2).astype(np.float32)      # épithéliums denses pénalisés
Gm, idx_m, pos_m, lut_m = graphe(Ms, cost)
d_nt, j_nt = tree_n.query(pos_m); s_m = s_n[j_nt]            # chaque nœud : distance au tube et s du point de tube le plus proche
vent = ((pos_m - pos_n[j_nt]) * (Cbody - pos_n[j_nt])).sum(1) > 0     # côté ventral = vers le centre du corps (embryon en C)
side_m = np.sign(pos_m[:, 0] - MID / F)
log('  mésenchyme : %d nœuds, %d arêtes' % (len(idx_m), Gm.nnz // 2))
def noeuds(mask_s, it=1):
    """nœuds du graphe mésenchyme dans ou au contact (26-voisinage, it couches) d'un masque sous-grille"""
    j = lut_m[np.flatnonzero((ndi.binary_dilation(mask_s, iterations=it) & Ms).ravel())]; return j[j >= 0]
orig_s_nodes = np.flatnonzero(keep_s.ravel())              # indices plats sous-grille des voxels d'origine
def sources(s_lo, s_hi):
    sel = keep_s.copy(); sel.ravel()[orig_s_nodes] = (s_keep.ravel()[orig_s_nodes] >= s_lo) & (s_keep.ravel()[orig_s_nodes] <= s_hi)
    return noeuds(sel, 2)                                   # 2 couches : le mésenchyme commence à 2 voxels /F de la surface du tube
def trajet(src, tgt, lim_factor=4.0):
    """plus court chemin pondéré d'un ensemble de sources vers l'ensemble de cibles ; liste de positions sous-grille ou None"""
    src = np.unique(src); tgt = np.unique(tgt)
    if len(src) == 0 or len(tgt) == 0: return None
    d0 = np.linalg.norm(pos_m[src].mean(0) - pos_m[tgt].mean(0)); lim = lim_factor * d0 + 30
    dist, pred, _ = dijkstra(Gm, directed=False, indices=src, min_only=True, return_predecessors=True, limit=lim)
    dt = dist[tgt]
    if not np.isfinite(dt).any(): return None
    node = int(tgt[np.argmin(dt)]); path = [node]
    while pred[node] >= 0: node = int(pred[node]); path.append(node)
    return pos_m[path[::-1]]
chemins = []                                                 # (structure, côté, label cible, positions sous-grille)
def ajouter(nom, cote, cible, P):
    if P is not None and len(P) >= 2: chemins.append((nom, cote, cible, P)); return True
    return False
SIDES = (('gauche', -1), ('droit', 1))
# --- crête crâniale → face (3 bandes rostro-caudales : fronto-nasale, maxillaire, mandibulaire) par côté
d_face = max(2.0, vx(0.25) / F)
face = (s_m < s_otic) & (d_nt >= d_face) & vent
if coeur.any(): face &= ~np.isin(np.arange(len(idx_m)), noeuds(down(dilate_r(coeur, int(vx(0.3))))))
if (yeux | otic).any(): face &= ~np.isin(np.arange(len(idx_m)), noeuds(down(yeux | otic)))
if face.any(): face &= d_nt >= np.percentile(d_nt[face], 60)                  # cœur du mésenchyme facial (le plus profond sous l'encéphale)
bandes = np.linspace(0.25 * s_otic, s_otic, 4)
for b in range(3):
    src = sources(bandes[b], bandes[b + 1])
    for cote, sg in SIDES:
        ajouter('crete_craniale', cote, 'face (bande %d)' % (b + 1), trajet(src, np.flatnonzero(face & (side_m == sg))))
# --- crête cardiaque → pôle artériel (bulbe au contact de l'aorte commune, sinon sommet du cœur)
bulbe = np.zeros(shape, bool)
if coeur.any():
    if aorte.any(): bulbe = coeur & dilate_r(aorte, int(vx(0.25)))
    if not bulbe.any():
        si_c = np.where(coeur.any(axis=(0, 2)))[0]; bulbe = coeur.copy(); bulbe[:, int(si_c.min() + 0.15 * np.ptp(si_c)):, :] = False
    tgt = noeuds(down(bulbe)); src = sources(s_otic, s_otic + 3 * pas)
    for cote, sg in SIDES: ajouter('crete_cardiaque', cote, 'pole arteriel du coeur', trajet(src, tgt[side_m[tgt] == sg]))
# --- crête vagale → intestin antérieur ; crête sacrée → intestin postérieur
def partie_gut(avant=True):
    if gut_segs:
        for k in (SEGMENTS if avant else SEGMENTS[::-1]):
            if k in gut_segs: return gut_segs[k], k
    if not gut.any(): return None, None
    si_g = np.where(gut.any(axis=(0, 2)))[0]; m = gut.copy()
    if avant: m[:, int(si_g.min() + 0.25 * np.ptp(si_g)):, :] = False
    else: m[:, :int(si_g.max() - 0.25 * np.ptp(si_g)), :] = False
    return m, 'tube_digestif (%s)' % ('quart cranial' if avant else 'quart caudal')
m_av, nom_av = partie_gut(True); m_ar, nom_ar = partie_gut(False)
if m_av is not None:
    tgt = noeuds(down(m_av)); src = sources(s_otic + 1 * pas, s_otic + 7 * pas)
    for cote, sg in SIDES: ajouter('crete_vagale', cote, nom_av, trajet(src, tgt[side_m[tgt] == sg] if (side_m[tgt] == sg).any() else tgt))
if m_ar is not None and nom_ar != nom_av:
    tgt = noeuds(down(m_ar)); src = sources(s_otic + 28 * pas, Ls)
    for cote, sg in SIDES: ajouter('crete_sacree', cote, nom_ar, trajet(src, tgt[side_m[tgt] == sg] if (side_m[tgt] == sg).any() else tgt, 6.0))
# --- gaine entérique : front de colonisation selon le stade (vagale cranio-caudale ; sacrée rejoint le colon en fin de période)
COLONISES = {13: 2, 14: 3, 15: 4, 16: 4, 17: 5, 18: 5, 19: 5, 20: 5}        # nombre de segments colonisés (oesophage → intestin postérieur)
FRACTION = {13: 0.35, 14: 0.50, 15: 0.65, 16: 0.80}                         # si le tube n'est pas segmenté : fraction crâniale colonisée
enterique = np.zeros(shape, bool); front = ''
if gut.any():
    if gut_segs:
        pres = [k for k in SEGMENTS if k in gut_segs]; nb = min(COLONISES.get(num, 5), len(pres)); couvert = U(*pres[:nb]); front = pres[nb - 1]
    else:
        f = FRACTION.get(num, 1.0); si_g = np.where(gut.any(axis=(0, 2)))[0]; couvert = gut.copy()
        couvert[:, int(si_g.min() + f * np.ptp(si_g)):, :] = False; front = 'tube_digestif %.0f %%' % (100 * f)
    enterique = dilate_r(couvert, max(2, int(round(vx(0.05))))) & ~gut & env & ~NTf & ~coeur & ~foie & ~aorte
# --- crête du tronc → ganglions rachidiens (segmentaire) ; chaîne sympathique : ganglion → aorte dorsale
n_tronc = n_cran = n_symp = 0
if gangl.any():
    labg, info = comps_info(gangl, max(20, int(vx(0.02) ** 3)))
    info.sort(key=lambda x: -x[1])
    tgt_ao = noeuds(down(aorte)) if aorte.any() else np.array([], int)
    for gi, gsize, gc in info[:120]:
        sg = s_of(gc); comp_s = down(labg == gi); tgt = noeuds(comp_s)
        src = sources(sg - 0.7 * pas, sg + 0.7 * pas)
        if len(src) == 0: src = sources(sg - 1.5 * pas, sg + 1.5 * pas)
        cote = 'gauche' if gc[0] < MID else 'droit'
        if sg < s_otic + 3 * pas: n_cran += ajouter('crete_craniale', cote, 'ganglion cranien', trajet(src, tgt))
        else:
            n_tronc += ajouter('crete_tronc', cote, 'ganglion rachidien', trajet(src, tgt))
            if len(tgt_ao):
                ta = tgt_ao[np.abs(s_m[tgt_ao] - sg) <= 1.5 * pas]
                n_symp += ajouter('crete_sympathique', cote, 'aorte dorsale', trajet(tgt, ta, 6.0))
log('  trajets : %d (dont tronc %d, ganglions craniens %d, sympathique %d)' % (len(chemins), n_tronc, n_cran, n_symp))
# ------------------------------------------------------------------ 4. rasterisation des trajets en tubes
def tube(P_sub, r):
    P = P_sub * F + (F - 1) / 2.0
    d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    u = np.arange(0, d[-1] + 1e-6, 1.0); C = np.stack([np.interp(u, d, P[:, i]) for i in range(3)], 1)
    lo = np.maximum(np.floor(C.min(0)) - r - 1, 0).astype(int); hi = np.minimum(np.ceil(C.max(0)) + r + 2, shape).astype(int)
    m = np.zeros(hi - lo, bool)
    for c in C - lo:
        ci = np.round(c).astype(int); s_ = tuple(slice(max(ci[k] - r - 1, 0), min(ci[k] + r + 2, m.shape[k])) for k in range(3))
        g = np.ogrid[s_[0], s_[1], s_[2]]; m[s_] |= (g[0] - c[0]) ** 2 + (g[1] - c[1]) ** 2 + (g[2] - c[2]) ** 2 <= r * r
    full = np.zeros(shape, bool); full[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]] = m; return full
masques = {'crete_origine': origine}
for nom, cote, cible, P in chemins:
    masques[nom] = masques.get(nom, np.zeros(shape, bool)) | tube(P, r_tube)
if enterique.any(): masques['crete_enterique'] = enterique
if bulbe.any(): masques['crete_bulbe_cardiaque'] = bulbe
for k in [k for k, m in masques.items() if not m.any()]: del masques[k]
for k in masques: masques[k] &= env
# ------------------------------------------------------------------ 5. sorties
wd = os.path.join(work, 'crete_neurale'); od = os.path.join(out, 'crete_neurale'); os.makedirs(wd, exist_ok=True); os.makedirs(od, exist_ok=True)
np.savez_compressed(os.path.join(wd, 'crete_neurale.npz'), **{k: np.packbits(v) for k, v in masques.items()}, shape=np.array(shape))
def to_mm(P_full):
    P = np.asarray(P_full, float); c = center
    return np.stack([(P[:, 0] - c[0]) * vox, (P[:, 2] - c[2]) * vox, -(P[:, 1] - c[1]) * vox], 1)      # repère Blender du pipeline
PHASE = {13: 'migration active : crête crâniale dans les arcs, cardiaque en route, entérique dans l intestin antérieur, ganglions rachidiens en condensation',
         14: 'migration active : crête cardiaque atteint le pôle artériel, colonisation entérique jusqu au duodénum',
         15: 'crête crâniale en place (ectomésenchyme facial), colonisation entérique de l intestin moyen',
         16: 'début de la septation du tronc artériel, intestin moyen colonisé',
         17: 'migration essentiellement terminée ; dérivés en différenciation (cartilages faciaux, ganglions, plexus entériques jusqu au colon)',
         19: 'dérivés en différenciation, chaîne sympathique et plexus entériques complets',
         20: 'dérivés en différenciation, chaîne sympathique et plexus entériques complets'}
DERIVES = {'ganglions': 'ganglions rachidiens et crâniens (neurones sensitifs, cellules gliales)', 'chondrocrane': 'cartilages de la face et du neurocrâne antérieur (ectomésenchyme)',
           'cartilage_autre': 'cartilages des arcs pharyngiens (partie)', 'vaisseaux_aorte': 'paroi musculo-conjonctive des arcs aortiques (cellules murales)',
           'coeur_detoure': 'septum aortico-pulmonaire du bulbe', 'coeur': 'septum aortico-pulmonaire du bulbe', 'meninges_mesenchyme_cranien': 'méninges antérieures (partie)',
           'epiderme_cranien': 'mélanocytes (voie dorso-latérale, non modélisée)'}
rec = {'stage': stage, 'units': 'mm', 'mm_per_voxel': vox, 'center_voxel': center.tolist(), 'F': F, 'rayon_tube_vox': r_tube,
       'reperes': {'arc_tube_neural_mm': mm(Ls * F), 'niveau_otique_s_mm': mm(s_otic * F), 'pas_somitique_mm': mm(pas * F),
                   'otique_mesure': bool(otic.any()), 'front_enterique': front},
       'phase_du_stade': PHASE.get(num, ''), 'derives_existants': {k: v for k, v in DERIVES.items() if k in L},
       'trajets': [{'structure': nom, 'cote': cote, 'cible': cible, 's_depart_mm': mm(float(s_m[lut_m[np.ravel_multi_index(tuple(P[0]), sub_shape)]]) * F),
                    'longueur_mm': float(np.sum(np.linalg.norm(np.diff(P * F, axis=0), axis=1))) * vox,
                    'xyz_vox': (P * F + (F - 1) / 2.0).round(1).tolist(), 'xyz_mm_blender': to_mm(P * F + (F - 1) / 2.0).round(3).tolist()} for nom, cote, cible, P in chemins],
       'note': 'xyz_vox en coordonnees pleines de labels.npz (LR, SI, AP) ; trajets = plus courts chemins ponderes dans le mesenchyme, schematiques'}
json.dump(rec, open(os.path.join(wd, 'chemins.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
# couleurs / confiance : mêmes valeurs que meshexport.STYLE (import tardif pour rester lisible ici)
from meshexport import STYLE, CONFIANCE
pc = os.path.join(work, 'confiance.json'); conf_stage = json.load(open(pc, encoding='utf-8')) if os.path.exists(pc) else {}
manifest = {'stage': stage, 'units': 'mm', 'mm_per_voxel': vox, 'center_voxel': center.tolist(), 'source': 'embryo3d/crete_neurale.py (labels existants, trajets Dijkstra)',
            'reperes': rec['reperes'], 'phase_du_stade': rec['phase_du_stade'], 'derives_existants': rec['derives_existants'], 'structures': []}
if not args.sans_ply:
    cc = np.array([center[0], center[2], -center[1]]) * vox
    for k, m in masques.items():
        coll, col, alpha, tf, sig = STYLE.get(k, ('Crete neurale', (0.6, 0.3, 0.8), 1.0, 40000, 0.8))
        mesh = mask_to_mesh(m, vox, sigma=sig, target_faces=tf, min_faces_component=20)
        if mesh is None: continue
        mesh.apply_translation(-cc); fn = f'{stage}_{k}.ply'; mesh.export(os.path.join(od, fn))
        conf = conf_stage.get(k, args.confiance_trajets if k not in ('crete_origine', 'crete_enterique', 'crete_bulbe_cardiaque') else CONFIANCE.get(k, 'moyenne'))
        manifest['structures'].append({'name': k, 'file': fn, 'collection': coll, 'color': list(col), 'alpha': alpha, 'faces': int(len(mesh.faces)),
                                       'volume_mm3': float(m.sum()) * vox ** 3, 'confiance': conf, 'trajets': sum(1 for c in chemins if c[0] == k)})
        log('  %-24s %7d faces %8.4f mm3  confiance %s' % (k, len(mesh.faces), m.sum() * vox ** 3, conf))
json.dump(manifest, open(os.path.join(od, 'manifest_crete_neurale.json'), 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
# planche de contrôle : projections sagittale et frontale
try:
    import cv2
    COUL = {k: tuple(int(255 * c) for c in STYLE[k][1][::-1]) for k in masques if k in STYLE}
    def planche(axis):
        proj = lambda m: m.any(axis=0) if axis == 0 else m.any(axis=2).T        # lignes = SI ; colonnes = AP (sagittale) ou LR (frontale)
        base = np.zeros(proj(env).shape + (3,), np.uint8); base[proj(env)] = (45, 45, 45); base[proj(organes)] = (85, 85, 85); base[proj(NTf)] = (120, 120, 120)
        for k, m in masques.items(): base[proj(m)] = COUL.get(k, (200, 100, 255))
        return base
    im = [planche(0), planche(2)]                            # sagittale (lignes SI, colonnes AP) ; frontale (lignes SI, colonnes LR)
    h = max(i.shape[0] for i in im); sc = max(1, 900 // h)
    im = [cv2.resize(i, None, fx=sc, fy=sc, interpolation=cv2.INTER_NEAREST) for i in im]
    pano = np.zeros((max(i.shape[0] for i in im) + 30 + 22 * (len(masques) + 1), sum(i.shape[1] for i in im) + 20, 3), np.uint8); x = 0
    for i in im: pano[30:30 + i.shape[0], x:x + i.shape[1]] = i; x += i.shape[1] + 20
    cv2.putText(pano, '%s crete neurale - sagittale (gauche) / frontale (droite) ; %d trajets' % (stage, len(chemins)), (6, 20), cv2.FONT_HERSHEY_SIMPLEX, 0.55, (255, 255, 255), 1)
    y = 30 + max(i.shape[0] for i in im) + 18
    for k in masques:
        cv2.rectangle(pano, (8, y - 10), (24, y + 2), COUL.get(k, (200, 100, 255)), -1); cv2.putText(pano, k, (32, y), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (230, 230, 230), 1); y += 22
    cv2.imwrite(os.path.join(od, 'controle_crete_neurale.png'), pano)
except ImportError:
    log('  (cv2 absent : pas de planche de contrôle)')
# ------------------------------------------------------------------ 6. fusion optionnelle dans labels.npz
if args.fusion and masques:
    src_npz = os.path.join(work, 'labels.npz')
    with np.load(src_npz) as z: labels = {k: z[k] for k in z.files if k != 'shape'}; shape2 = tuple(z['shape'])
    if shape2 != shape: sys.exit('labels.npz a changé de forme pendant le calcul : fusion annulée')
    for k, m in masques.items(): labels[k] = np.packbits(m)
    shutil.copy(src_npz, os.path.join(work, 'labels_avant_crete.npz'))
    tmp = os.path.join(work, 'labels.tmp.npz'); np.savez_compressed(tmp, **labels, shape=np.array(shape)); os.replace(tmp, src_npz)
    log('  fusion : %s ajoutés à labels.npz (sauvegarde labels_avant_crete.npz)' % ', '.join(masques))
    import meshexport; meshexport.run(work, out, stage)
    try:
        import export_nrrd; export_nrrd.run(stage_dir)
    except Exception as e: log('  NRRD non réécrit :', e)
log('terminé :', ', '.join('%s %.3f mm3' % (k, m.sum() * vox ** 3) for k, m in masques.items()))
