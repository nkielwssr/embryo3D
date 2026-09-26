# -*- coding: utf-8 -*-
"""Somites / dermomyotomes et masses musculaires des membres, en reconstruction séparée (out/muscles/), alignée sur le pipeline.

    python embryo3d/muscles_build.py <dossier_stade> [--dry] [--t-somite 70] [--t-muscle 75] [--band 0.15]

Somites (CS13-CS17) : le tube neural (label snc) est squelettisé ; chaque point du squelette porte sa distance le long du
tube depuis la tête et le rayon local du tube. Candidats = tissu dense (ds > t_somite) de l'enveloppe, hors autres labels,
à moins de `band` mm de la surface du tube et décalé latéralement (gauche-droite) d'au moins le rayon local + 2 voxels par
rapport au point de squelette le plus proche ; on ne garde que la partie « tronc » (au-delà de la tête, repérée par la
chute de section du SNC). Une colonne par côté. Le nombre de blocs = maxima du profil de volume le long du tube (info :
la colonne n'est pas découpée). Muscles des membres (label « membres », CRL ≥ 11 mm) : tissu dense (ds > t_muscle) au
cœur des bourgeons (érodés de 3 voxels), hors cartilages, ouverture r=2, composantes > 2000 voxels.
Sorties : out/muscles/<CS>_somites_gauche|droite.ply, <CS>_muscles_membres.ply, manifest_muscles.json (format des
manifest_digestif), planche out/muscles/<CS>_controle.png (projections), masques work/muscles/muscles.npz.
Confiance déclarée : somites « moyenne » si le nombre de blocs est plausible pour le stade, sinon « faible » ;
muscles des membres « faible » (automatique, non vérifié).
"""
import json
import os
import re
import sys
import time

import numpy as np
from scipy import ndimage as ndi
from scipy.signal import find_peaks
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import meshexport  # noqa: E402

stage_dir = os.path.abspath(sys.argv[1])
def opt(nom, defaut):
    return type(defaut)(sys.argv[sys.argv.index(nom) + 1]) if nom in sys.argv else defaut
T_SOM, T_MUS, BAND, DRY = opt("--t-somite", 70.0), opt("--t-muscle", 75.0), opt("--band", 0.15), "--dry" in sys.argv
SANS_SOM = "--sans-somites" in sys.argv   # les somites sont segmentés par une autre session (label « somites » du pipeline)
work, out = os.path.join(stage_dir, "work"), os.path.join(stage_dir, "out")
stage = "CS%s" % re.search(r"CS\s*(\d+)", os.path.basename(stage_dir), re.I).group(1)
num = int(stage[2:])
man = json.load(open(os.path.join(out, "manifest.json")))
scale = float(man["mm_per_voxel"]); center = np.array(man["center_voxel"], float)
mm = lambda x: max(1, int(round(x / scale)))
t0 = time.time()
log = lambda *a: print(f"[{time.time()-t0:5.0f}s]", *a, flush=True)

z = np.load(os.path.join(work, "labels.npz"))
shape = tuple(int(x) for x in z["shape"]); n = int(np.prod(shape))
lab = {k: np.unpackbits(z[k])[:n].reshape(shape).astype(bool) for k in z.files if k != "shape"}
env = lab["enveloppe"]; cns = lab.get("snc", np.zeros(shape, bool))
ds = np.load(os.path.join(work, "ds1.npy"), mmap_mode="r")
S = ndi.generate_binary_structure(3, 1)
def dil(m, r): return ndi.binary_dilation(m, S, iterations=r) if r > 0 else m
def ero(m, r): return ndi.binary_erosion(m, S, iterations=r) if r > 0 else m
log(stage, "forme", shape, "mm/voxel %.4f" % scale, "labels", sorted(lab))
others = np.zeros(shape, bool)
for k, v in lab.items():
    if k != "enveloppe":
        others |= v
segments = {}

# ------------------------------------------------------------------ somites / dermomyotomes
if 13 <= num <= 17 and cns.any() and not SANS_SOM:
    from skimage.morphology import skeletonize
    area = cns.sum(axis=(0, 2)).astype(float); rows = np.flatnonzero(area > 0); amax = area.max()
    sm = ndi.uniform_filter1d(area, size=max(3, mm(0.15)))
    head_end = int(rows[0])
    for r in rows:
        if sm[r] > 0.35 * amax: head_end = int(r)
        elif r > head_end + mm(0.4): break
    trunk0 = head_end + mm(0.2)
    log("SNC : section max %d, fin de la tête rangée %d, tronc à partir de %d" % (amax, head_end, trunk0))
    # sous-volume autour du SNC dilaté de la bande
    pad = mm(BAND) + 4
    idx = np.argwhere(cns); lo = np.maximum(idx.min(0) - pad, 0); hi = np.minimum(idx.max(0) + pad + 1, shape)
    sl = tuple(slice(int(a), int(b)) for a, b in zip(lo, hi))
    c = cns[sl]; e = env[sl]; o = others[sl]; d = np.asarray(ds[sl]).astype(np.float32)
    # squelette du tube + distance le long du tube depuis la tête (BFS 26-connexe)
    sk = skeletonize(c)
    pts = np.argwhere(sk)
    log("squelette : %d points" % len(pts))
    tree = cKDTree(pts)
    pairs = tree.query_pairs(r=1.75)
    import collections
    adj = collections.defaultdict(list)
    for a, b in pairs:
        adj[a].append(b); adj[b].append(a)
    start = int(np.argmin(pts[:, 1]))                 # point le plus haut (tête)
    along = np.full(len(pts), -1.0); along[start] = 0.0
    q = collections.deque([start])
    while q:
        u = q.popleft()
        for v in adj[u]:
            if along[v] < 0:
                along[v] = along[u] + float(np.linalg.norm(pts[v] - pts[u])); q.append(v)
    orphan = along < 0
    if orphan.any():                                  # branches non connectées : distance approx. par la rangée
        along[orphan] = (pts[orphan, 1] - pts[start, 1]).clip(0)
    inside = ndi.distance_transform_edt(c)            # rayon local du tube aux points du squelette
    rt = inside[tuple(pts.T)]
    # encéphale = plus grande composante du SNC érodé (seul le cerveau survit à une érosion de 1,6 × rayon médian) + marge 0,25 mm
    r_med = float(np.median(rt)); k_er = max(2, int(round(1.6 * r_med)))
    core = ero(c, k_er); lc, kc = ndi.label(core, S)
    brain = np.zeros_like(c)
    if kc:
        cnt = np.bincount(lc.ravel()); cnt[0] = 0
        brain |= dil(lc == int(cnt.argmax()), k_er + mm(0.25))
    # + ventricules cérébraux (épais : survivent à une érosion de 0,04 mm, le canal central non) dilatés de 0,4 mm
    if "ventricules" in lab:
        vth = ero(lab["ventricules"][sl], mm(0.04))
        if vth.any(): brain |= dil(vth, mm(0.04) + mm(0.4))
    log("tête exclue : %.2f M voxels (SNC érodé %d voxels : %.2f M ; ventricules)" % (brain.sum()/1e6, k_er, (lc == int(cnt.argmax())).sum()/1e6 if kc else 0))
    sens = np.zeros_like(c)
    for k in ("yeux", "vesicules_otiques", "chondrocrane", "cristallins"):
        if k in lab: sens |= lab[k][sl]
    sens = dil(sens, mm(0.12)) if sens.any() else sens
    outside = ndi.distance_transform_edt(~c)          # distance à la surface du tube
    cand = e & ~c & (outside <= mm(BAND)) & (d > T_SOM) & ~dil(o, 2)
    cidx = np.argwhere(cand)
    if len(cidx):
        _, near = tree.query(cidx, k=1)
        xs = pts[near, 0]; lat = np.abs(cidx[:, 0] - xs)
        ok = (lat >= rt[near] + 2) & (cidx[:, 1] + sl[1].start >= trunk0) & ~brain[tuple(cidx.T)] & ~sens[tuple(cidx.T)]
        cand2 = np.zeros_like(cand); cand2[tuple(cidx[ok].T)] = True
        cand2 = dil(ero(cand2, 1), 1)
        cand2 = ndi.binary_fill_holes(cand2)
        side_left = np.zeros_like(cand); side_left[tuple(cidx[ok].T)] = cidx[ok, 0] < xs[ok]
        attendu = {13: (22, 40), 14: (26, 44), 15: (28, 48), 16: (28, 50), 17: (28, 50)}[num]
        for side, sel in (("gauche", side_left), ("droite", ~side_left)):
            m = cand2 & sel
            l, k = ndi.label(m, S)
            if k == 0: continue
            cnt = np.bincount(l.ravel()); cnt[0] = 0
            keep = [i for i in range(1, k + 1) if cnt[i] >= max(300, cnt.max() / 30)]
            m = np.isin(l, keep)
            vidx = np.argwhere(m); _, nv = tree.query(vidx, k=1)
            a_mm = along[nv] * scale
            bins = np.arange(0, a_mm.max() + 0.02, 0.02)
            prof, _ = np.histogram(a_mm, bins=bins)
            prof = ndi.gaussian_filter1d(prof.astype(float), 1.0)
            pk, _ = find_peaks(prof, distance=max(2, int(0.06 / 0.02)), prominence=0.15 * prof.max() if prof.max() > 0 else 1)
            full = np.zeros(shape, bool); full[sl] = m
            okc = attendu[0] <= len(pk) <= attendu[1]
            segments[f"somites_{side}"] = {"mask": full, "n_blocs": int(len(pk)), "confiance": "moyenne" if okc else "faible",
                                           "note": f"{len(pk)} blocs détectés le long du tube (attendu {attendu[0]}-{attendu[1]} paires) ; "
                                                   f"colonne non découpée ; longueur suivie {a_mm.max():.2f} mm",
                                           "color": [0.75, 0.25, 0.35] if side == "gauche" else [0.85, 0.35, 0.45]}
            log(f"somites {side} : {m.sum()/1e3:.0f} k voxels, {len(pk)} blocs, longueur {a_mm.max():.2f} mm")
    else:
        log("aucun candidat somite")

# ------------------------------------------------------------------ muscles des membres
if "membres" in lab and meshexport.CRL_MM.get(stage, 0) >= 11.0:
    limbs = lab["membres"]
    core = ero(limbs, 3)
    cart = np.zeros(shape, bool)
    for k in ("squelette_axial_cartilage", "chondrocrane", "cartilage_autre"):
        if k in lab: cart |= lab[k]
    m = core & (np.asarray(ds) > T_MUS) & ~dil(cart, 2) & ~cns
    m = dil(ero(m, 2), 2)
    l, k = ndi.label(m, S)
    if k:
        cnt = np.bincount(l.ravel()); cnt[0] = 0
        m = np.isin(l, np.flatnonzero(cnt >= 2000))
    if m.any():
        segments["muscles_membres"] = {"mask": m, "color": [0.70, 0.30, 0.30], "confiance": "faible",
                                       "note": "tissu dense au cœur des bourgeons, hors cartilage ; automatique non vérifié"}
        log(f"muscles des membres : {m.sum()/1e3:.0f} k voxels")
    else:
        log("muscles des membres : rien de dense trouvé dans les bourgeons")

if not segments:
    raise SystemExit("aucun segment produit")

# ------------------------------------------------------------------ planche de contrôle (projections)
od = os.path.join(out, "muscles"); os.makedirs(od, exist_ok=True)
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(15, 6))
    views = [("face (LR × SI)", (2,)), ("profil (AP × SI)", (0,)), ("dessus (LR × AP)", (1,))]
    for ax, (ttl, axis) in zip(axs, views):
        base = env.any(axis=axis[0]).astype(float) * 0.25 + cns.any(axis=axis[0]).astype(float) * 0.5
        img = np.stack([base, base, base], -1)
        for i, (name, seg) in enumerate(segments.items()):
            pm = seg["mask"].any(axis=axis[0])
            col = np.array(seg["color"]); img[pm] = col
        if axis[0] == 2: img = img.transpose(1, 0, 2)          # SI en vertical
        elif axis[0] == 0: img = img.transpose(1, 0, 2)
        ax.imshow(img, origin="upper"); ax.set_title(f"{stage} — {ttl}"); ax.axis("off")
    fig.suptitle(", ".join(f"{k} : {v['confiance']}" + (f" ({v['n_blocs']} blocs)" if "n_blocs" in v else "") for k, v in segments.items()))
    fig.tight_layout(); fig.savefig(os.path.join(od, f"{stage}_controle.png"), dpi=90); plt.close(fig)
    log("planche :", os.path.join(od, f"{stage}_controle.png"))
except Exception as ex:
    log("planche impossible :", ex)
if DRY:
    log("--dry : rien d'autre écrit"); sys.exit(0)

# ------------------------------------------------------------------ export
wd = os.path.join(work, "muscles"); os.makedirs(wd, exist_ok=True)
manifest = {"stage": stage, "units": "mm", "mm_per_voxel": scale, "center_voxel": center.tolist(),
            "source": "embryo3d/muscles_build.py (automatique)", "segments": []}
cc = np.array([center[0], center[2], -center[1]]) * scale
for name, seg in segments.items():
    m = meshexport.mask_to_mesh(seg["mask"], scale, sigma=1.0, target_faces=80000)
    if m is None:
        log("  (vide)", name); continue
    m.apply_translation(-cc)
    fn = f"{stage}_{name}.ply"; m.export(os.path.join(od, fn))
    manifest["segments"].append({"name": name, "file": fn, "color": seg["color"], "alpha": 1.0, "faces": int(len(m.faces)),
                                 "volume_mm3": float(seg["mask"].sum()) * scale ** 3, "confiance": seg["confiance"],
                                 "note": seg.get("note", ""), **({"n_blocs": seg["n_blocs"]} if "n_blocs" in seg else {})})
    log(f"  {name:20s} {len(m.faces):7d} faces  {seg['mask'].sum()*scale**3:.3f} mm3  {seg['confiance']}")
json.dump(manifest, open(os.path.join(od, "manifest_muscles.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
np.savez_compressed(os.path.join(wd, "muscles.npz"), **{k: np.packbits(v["mask"]) for k, v in segments.items()}, shape=np.array(shape))
log("OK", stage, "->", od)
