# -*- coding: utf-8 -*-
"""Recalcule la cavité péricardique d'un stade à partir du cœur corrigé (labels.npz), sans toucher aux autres labels.

    python embryo3d/fix_pericarde.py <dossier_stade> [--thr 12] [--seed 6] [--steps 25] [--open 1] [--dry]

Règle : la cavité = espace de faible densité (ds < thr, dans l'enveloppe) qui entoure le cœur, obtenu par croissance
géodésique depuis la coque du cœur (rayon `seed`), en restant hors des autres structures, hors du voisinage du foie
(8 voxels) et du SNC (4 voxels), et à moins de `steps` voxels du cœur. Puis ouverture r=`open` et plus grande composante.
Après correction manuelle du cœur (Slicer, « cœur plein »), l'espace péricardique est une lame mince (CS16 : ds < 12,
à 4-8 voxels de la surface), d'où les valeurs par défaut plus douces que dans segment.py (ds < 3, ouverture 3).
Écrit labels.npz (sauvegarde labels_avant_pericarde.npz), régénère les maillages et les NRRD du stade.
"""
import json
import os
import re
import shutil
import sys

import numpy as np
from scipy import ndimage as ndi

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

stage_dir = os.path.abspath(sys.argv[1])
def opt(nom, defaut):
    return type(defaut)(sys.argv[sys.argv.index(nom) + 1]) if nom in sys.argv else defaut
steps, thr, seed_r, open_r = opt("--steps", 25), opt("--thr", 12.0), opt("--seed", 6), opt("--open", 1)
dry = "--dry" in sys.argv
work, out = os.path.join(stage_dir, "work"), os.path.join(stage_dir, "out")
stage = "CS%s" % re.search(r"CS\s*(\d+)", os.path.basename(stage_dir), re.I).group(1)

z = np.load(os.path.join(work, "labels.npz"))
shape = tuple(int(x) for x in z["shape"]); n = int(np.prod(shape))
lab = {k: np.unpackbits(z[k])[:n].reshape(shape).astype(bool) for k in z.files if k != "shape"}
if "coeur" not in lab or not lab["coeur"].any():
    raise SystemExit("pas de label coeur dans labels.npz : rien à faire")
heart = lab["coeur"]
print(f"{stage} : forme {shape}, coeur {heart.sum()/1e6:.3f} M voxels, labels {sorted(lab)}")

# --- sous-volume autour du cœur (marge = steps + 10) pour rester léger
idx = np.argwhere(heart)
lo = np.maximum(idx.min(0) - (steps + 10), 0); hi = np.minimum(idx.max(0) + (steps + 11), shape)
sl = tuple(slice(int(a), int(b)) for a, b in zip(lo, hi))
print("sous-volume", [s_.start for s_ in sl], "->", [s_.stop for s_ in sl])

ds = np.load(os.path.join(work, "ds1.npy"), mmap_mode="r")[sl].astype(np.float32)
env = np.load(os.path.join(work, "env.npy"), mmap_mode="r")[sl].astype(bool)
h = heart[sl]
others = np.zeros_like(h)
for k, v in lab.items():
    if k not in ("enveloppe", "coeur", "cavite_pericardique"):
        others |= v[sl]
S = ndi.generate_binary_structure(3, 1)
def dil(m, r): return ndi.binary_dilation(m, S, iterations=r) if r > 0 else m
def ero(m, r): return ndi.binary_erosion(m, S, iterations=r) if r > 0 else m
def largest(m):
    l, k = ndi.label(m, S)
    if k == 0: return m
    c = np.bincount(l.ravel()); c[0] = 0
    return l == int(c.argmax())

cav = env & (ds < thr) & ~others & ~h
excl = np.zeros_like(h)
if "foie" in lab: excl |= dil(lab["foie"][sl], 8)
if "snc" in lab: excl |= dil(lab["snc"][sl], 4)
domaine = cav & ~excl
seed = domaine & dil(h, seed_r)
print(f"cavité candidate {cav.sum()/1e6:.3f} M, graines autour du coeur {seed.sum()/1e3:.0f} k")
peri = ndi.binary_dilation(seed, S, iterations=steps, mask=domaine)
peri = dil(ero(peri, open_r), open_r)  # ouverture
peri = largest(peri)
# bouche les petits trous et retire ce qui recouvrirait le cœur
peri = ndi.binary_fill_holes(peri) & ~h & ~others
ratio = peri.sum() / max(1, h.sum())
print(f"cavité péricardique : {peri.sum()/1e6:.3f} M voxels (rapport au coeur {ratio:.2f})")
if not (0.05 <= ratio <= 4.0):
    raise SystemExit("rapport cavité/coeur hors plage plausible [0.05, 4] : vérifier avant d'écrire (--dry)")
if dry:
    print("--dry : rien écrit"); sys.exit(0)

full = np.zeros(shape, bool); full[sl] = peri
lab["cavite_pericardique"] = full
shutil.copy(os.path.join(work, "labels.npz"), os.path.join(work, "labels_avant_pericarde.npz"))
np.savez_compressed(os.path.join(work, "labels.npz"), **{k: np.packbits(v) for k, v in lab.items()}, shape=np.array(shape))
st = os.path.join(work, "labels_stats.json")
try:
    stats = json.load(open(st)); stats["cavite_pericardique"] = int(full.sum()); stats["coeur"] = int(heart.sum())
    json.dump(stats, open(st, "w"), indent=1)
except Exception:
    pass
print("labels.npz mis à jour ; export des maillages et NRRD...")
import meshexport, export_nrrd  # noqa: E402
meshexport.run(work, out, stage)
export_nrrd.run(stage_dir)
print("OK", stage, ": cavite_pericardique rétablie ; relancer le blend du stade et bash embryo3d/build_all.sh (ou make_master_manifest + scène maître)")
