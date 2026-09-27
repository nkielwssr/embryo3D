# Plan de calcul pour la machine RTX 3080 (i9 9900, 32 Go) — ce qui ferait vraiment progresser le modèle

Rédigé le 26/09/2026 par la session cloud « Modèle d'embryon complet stades 10-23 », **corrigé le 27/09** d'après la session locale
« Agrégateur modèles 3D et temps embryon » (correction validée par l'utilisateur).

## Correctif du 27/09 : les sources, et ce qui est déjà fait

- Les **vidéos 360° `CSxx_f4v` viennent du HDBR** (tomographie optique ; CS13 = spécimen N475). Ce ne sont **pas** les embryons Carnegie.
- **Hikspoors et Amsterdam** sont les embryons Carnegie de nos coupes ehd / du modèle VHE (3709, 6330, 6344, 8943, 836, 6502, 721, 6517,
  6520, 4430, 462, 9226). Le calage CS13 à 41 µm est fait **sur VHE**, pas sur le volume vidéo.
- **P1 (segmentation apprise) est FAIT localement, sur VHE, avec la RTX 3080** : nnU-Net 19 classes (espacement en pixels), jeu 502
  (CS13, CS16, CS20), Dice EMA 0,875 ; jeu 503 prêt (+ CS23) ; contrôle par les légendes ehd (`valider_ia.py`). Labels VHE de référence
  après les correctifs : CS13 0,91, CS16 0,891, CS20 0,872, CS23 0,855. **CS10 reste à 0,46** (le gris d'Amsterdam ne s'accorde pas avec
  nos coupes). **CS9 et les stades jeunes sont hors du domaine appris.**
- Conséquences : P1 devient « étendre l'IA aux stades jeunes et à CS23 avec Hikspoors calé sur VHE comme exemples » ; P2 (atlas 4D) vient
  après ; P3 et P4 (vidéos HDBR : autres spécimens) passent en priorité basse.

## Diagnostic mis à jour : où le modèle est faible

1. **Stades jeunes (CS9–CS12)** : hors domaine de l'IA, Amsterdam en désaccord avec nos coupes (CS10 Dice 0,46). C'est là que les maillages
   Hikspoors calés sur VHE (12 stades, planches justes, calage json par stade) apportent des labels experts du cœur, des vaisseaux, du tube
   neural, des somites, de l'intestin.
2. **Sources non recalées entre elles** hors VHE↔Hikspoors : vidéos HDBR (autre spécimen, échelle CRL typique), hulls VOKA.
3. **Le passage d'un stade à l'autre est un Shrinkwrap** : correct pour des blobs, pas pour le cerveau, les tubes, la colonne ; CS19, CS21,
   CS22 ne sont qu'interpolés (pas de PDF Hikspoors ; CS19 vient de la vidéo HDBR).
4. **Enveloppe et aspect** : occlusion sur densité, pas de texture réelle (sauf ce que VHE et les vidéos HDBR peuvent donner).

## Priorités, dans l'ordre

### P1' — Étendre la segmentation apprise (VHE) aux stades jeunes et à CS23, avec Hikspoors comme exemples

Pourquoi : le réseau 502/503 est bon de CS13 à CS23 mais ne connaît pas CS9–CS12, et CS10 est faux. Les maillages Hikspoors, **calés sur
VHE** (`embryons_3D/modeles/CSxx_hikspoors/`, calage json de la banque, repère modeles), se voxelisent dans la grille VHE de chaque stade
et donnent des labels experts là où Amsterdam se trompe.
Recette : (a) voxeliser les structures Hikspoors dans la grille VHE de CS9, CS10, CS11, CS12 (et CS14–CS18 pour densifier) avec la
correspondance nomenclature → 19 classes (à écrire : `hikspoors_vers_classes.json`) ; (b) là où Hikspoors et Amsterdam se contredisent
(CS10), Hikspoors fait foi ; le reste des voxels reste « ignore » ; (c) jeu 504 = 503 + stades jeunes ; entraînement (une nuit) ;
(d) validation par `valider_ia.py` (légendes ehd) et par un stade tenu à l'écart ; (e) export dans `modeles/CSxx_recon/` avec Dice et
confiance dans le manifest.
Attendu : Dice > 0,8 sur CS10–CS12, modèles `_recon` pour les stades jeunes, donc une chaîne CS9→CS23 entièrement « nos modèles »
(ordre validé : `_recon` > `_hikspoors` > atlas). Risque : Hikspoors est centré sur le cœur ; peau, membres, yeux restent Amsterdam/coupes
à ces stades.

### P2 — Atlas spatio-temporel CS9→CS23 par recalage déformable (remplace le Shrinkwrap)

Sur les volumes VHE et les cartes de labels IA (P1') : déformation CS9→CS10→…→CS23 (SyN d'ANTs sur CPU 8 cœurs, ou VoxelMorph/TransMorph
via MONAI sur GPU), régularisée dans le temps → morphing par champ, maillages à topologie commune (transport du maillage CS23 vers chaque
stade), vraie interpolation de CS19 (à partir de CS18/CS20), CS21, CS22. Les vertèbres par niveau, le tube neural, les tubes suivent le
champ sans casser. Sorties : `phi_CSxx_CSyy.nii.gz`, shape keys générées depuis le champ (`morph_par_champ.py` à écrire).
Temps : 1 à 3 h par paire sur GPU, une nuit pour la chaîne. Risque : grandes déformations CS12→CS13 (courbure) ; initialiser par les
repères pointés (yeux, otique, cœur, queue) et faire du multi-résolution. Lecture utile : « 4D Human Embryonic Brain Atlas » (arXiv 2503.07177).

### P5 — Rendu GPU (à installer dès le premier jour)

Cycles ou EEVEE sur la 3080 : les séquences de `blender_render_anim.py` passent de heures à minutes ; planches et vidéos de comparaison
VHE / Hikspoors / IA pour tout contrôler visuellement.

### Priorité basse : P3 (splatting gaussien des vidéos 360° HDBR) et P4 (fusion/restauration GPU des volumes vidéo)

Les vidéos HDBR sont d'autres spécimens : elles ne servent ni à P1' ni à P2. Elles restent utiles plus tard pour l'enveloppe texturée et pour
CS19 (seul stade où elles sont la source), avec le même outillage (Postshot/nerfstudio, CuPy, Noise2Void). À ne reprendre qu'une fois P1' et P2
livrés.

### À vérifier : les modèles OPT du HDBR atlas

Les vidéos 360° étant issues du HDBR, les **modèles OPT annotés du HDBR atlas (CS12→CS23, licence CC BY-NC-SA)** sont peut-être les mêmes
spécimens que nos volumes vidéo : leurs annotations labelliseraient directement ces volumes (CS19 compris). Le site est bloqué depuis le cloud :
à regarder depuis le PC (`hdbratlas.org/carnegie_stages.html`).

## Préparation de la clé USB

Logiciels (Windows) : pilote NVIDIA récent + CUDA 12.x ; Python 3.11 portable (WinPython ou conda-pack) avec `torch` (cu12), `monai`,
`nnunetv2`, `SimpleITK`, `antspyx`, `cupy-cuda12x`, `scikit-image`, `trimesh`, `fast_simplification`, `pynrrd`, `opencv-python`, `pillow` ;
Blender 4.x/5.x ; 3D Slicer (contrôle des labels) ; git. Postshot/nerfstudio + COLMAP seulement pour P3.
Données : le dépôt (branche `claude/trusting-clarke-3imqlz`) ; les volumes et labels VHE par stade et les jeux nnU-Net 502/503 (+ le modèle
entraîné) ; `embryons_3D/modeles/` (dont `CSxx_hikspoors/` calés) ; la banque `hikspoors_maastricht/` (npz, calage json, textes) ;
`coupes embryos 9-23/` pour `valider_ia.py` ; les vidéos HDBR seulement si P3/P4. Compter 20 à 30 Go.
Vérifications au premier démarrage : `python -c "import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))"`,
`python embryo3d/test_hikspoors_nomenclature.py`, `python embryo3d/test_hikspoors_synthetique.py`, un rendu Blender GPU d'une scène existante.

## Séquence proposée pour la première semaine

Jour 1 : installation, rendu GPU (P5), voxelisation Hikspoors → grille VHE pour CS13 (stade calé et validé) et contrôle Dice contre les
labels VHE de référence (mesure de la qualité du calage en voxels).
Jours 2–3 : voxelisation CS9–CS12 et CS14–CS18, table nomenclature → 19 classes, jeu 504.
Jour 4 : entraînement (nuit).
Jour 5 : validation (`valider_ia.py`, stade tenu à l'écart), export `modeles/CSxx_recon/` des stades jeunes, agrégateur.
Jours 6–7 : premières paires de recalage déformable (CS13→CS14, CS18→CS20), prototype de morphing par champ, comparaison avec le Shrinkwrap.

Livrables à remonter dans l'agrégateur : par stade, un bloc `ia` dans le manifest (Dice, confiance, jeu) ; par paire de stades, le champ de
déformation et le résidu.
