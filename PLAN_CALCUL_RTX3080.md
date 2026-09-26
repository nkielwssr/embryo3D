# Plan de calcul pour la machine RTX 3080 (i9 9900, 32 Go) — ce qui ferait vraiment progresser le modèle

Rédigé le 26/09/2026 par la session cloud « Modèle d'embryon complet stades 10-23 », à partir de l'état du dépôt embryo3D, des banques
locales (vidéos ehd, coupes légendées ehd, atlas d'Amsterdam, PDF Hikspoors, VOKA) et de ce que le GPU rend possible.

## Diagnostic : où le modèle est faible aujourd'hui

1. **La segmentation est heuristique** (seuils de densité, formes, positions) : SNC incomplet à CS16/CS17, tube digestif et vaisseaux
   « partiels », cartilage confondu avec les cavités, corrections à la main dans Slicer. Le morphing du cerveau a dû être retiré.
2. **Les sources ne sont pas recalées entre elles** : volumes vidéo (échelle = CRL typique), maillages Hikspoors (unités PDF, même spécimen),
   coupes légendées 2D (pas d'échelle ni de recalage entre sections), hulls VOKA. Chaque source vit dans son repère.
3. **Le passage d'un stade à l'autre est un Shrinkwrap** : ça marche pour des blobs, pas pour le cerveau, les tubes, la colonne ; les stades
   manquants (CS18, CS21, CS22 ; CS19 sans Hikspoors) ne sont qu'interpolés linéairement.
4. **L'enveloppe et l'aspect** viennent d'un critère d'occlusion sur la densité, sans couleur ni texture réelle ; VOKA n'est qu'un regard.

Un GPU ne change rien au point 2 (c'est du travail de géométrie, CPU suffit), mais il change tout aux points 1, 3 et 4.

## Priorités, dans l'ordre où je les ferais

### P1 — Transformer les maillages Hikspoors en vérité terrain, puis apprendre une segmentation 3D (le plus gros gain)

Pourquoi : Hikspoors a reconstruit **les mêmes spécimens Carnegie** que nos vidéos (CS13 #836, CS20 #462, CS23 #9226…). Une fois leurs
maillages recalés sur nos volumes fusionnés (déjà fait à 41 µm près pour CS13), on les voxelise et on obtient, pour 7 stades, des labels
experts du cœur (cavités, myocarde, coussins, valves), des veines, des artères, du tube neural, des somites, de l'intestin, du foie, des
poumons. C'est exactement ce qui manque à notre segmentation. Un réseau 3D entraîné dessus généralise ensuite aux stades sans PDF (CS19)
et remplace les heuristiques partout, avec une carte de confiance.

Entrées : `CSxx_f4v/work/dens.npy` (volumes isotropes ≈ 25 µm), maillages `modeles/CSxx_hikspoors/*.ply` recalés (`--calage`), nos labels
actuels comme labels faibles (peau, membres, yeux, cartilage), coupes légendées ehd rasterisées en labels épars là où elles existent.
Outils : nnU-Net v2 (ou MONAI) 3d_fullres, patches 128³, batch 2 : tient dans 10 Go de VRAM ; labels partiels gérés par une classe « ignore ».
Recette : (a) recalage rigide + affine + déformable léger Hikspoors → volume par stade (SimpleITK/ANTs, CPU, contrôle sur les coupes) ;
(b) voxelisation des maillages (trimesh `contains` ou rasterisation par tranches) → `labels_hikspoors.npz` ; (c) entraînement 5 folds
ou 1 fold + validation croisée par stade (laisser un stade de côté pour mesurer la généralisation) ; (d) inférence sur les 7 volumes +
CS19 ; (e) export dans `labels.npz` avec un préfixe `ia_` et confiance, puis `mesh` et comparaison Dice avec les heuristiques.
Temps : préparation 2 jours (CPU), entraînement 12 à 24 h (GPU), inférence minutes. Risque : la fusion vidéo a des artefacts de
compression ; prévoir une augmentation « flou + bruit de bloc ». Résultat attendu : cœur, vaisseaux, tube neural, intestin fiables à tous
les stades vidéo, donc morphing possible sur ces structures.

### P2 — Atlas spatio-temporel CS10→CS23 par recalage déformable (remplace le Shrinkwrap)

Pourquoi : un morphing correct est un **champ de déformation** entre stades consécutifs, pas une projection de surface. Avec des volumes
et des cartes de labels cohérents (P1), on estime une déformation CS13→CS14→…→CS23 (SyN d'ANTs sur CPU 8 cœurs, ou VoxelMorph/TransMorph
sur GPU via MONAI), on la régularise dans le temps, et on obtient un atlas continu : n'importe quel « jour » entre deux stades est une
interpolation géodésique, y compris CS18, CS21, CS22. Les vertèbres par niveau, le tube neural et le cerveau suivent le champ sans casser.
Entrées : volumes + labels de P1, maillages Hikspoors pour les stades sans vidéo (CS9–12, CS23 : voxelisés).
Sorties : champs de déformation `phi_CSxx_CSyy.nii.gz`, maillages à topologie commune (on transporte le maillage CS23 vers chaque stade),
shape keys Blender générées depuis le champ (un script `morph_par_champ.py` à écrire).
Temps : 1 à 3 h par paire sur GPU, une nuit pour la chaîne. Risque : grandes déformations entre CS12 et CS13 (courbure) ; faire du multi-
résolution et initialiser par les repères (yeux, otique, cœur, queue) déjà pointés.

### P3 — Enveloppe et texture réelles par splatting gaussien sur les vidéos 360°

Pourquoi : chaque stade a une vidéo de rendu volumique en rotation (perspective, fond noir). Le Gaussian Splatting (ou un NeRF) reconstruit
une surface colorée de très bonne qualité à partir de ces vues : on obtient l'enveloppe **avec l'aspect** (pigment rétinien, arcs, bourgeons)
au lieu de l'occlusion sur densité et du matériau VOKA procédural. Ça sert aussi de vérité pour l'orientation et la silhouette.
Outils : Postshot (Windows, simple) ou gsplat/nerfstudio ; COLMAP pour les poses (la rotation est régulière, on peut aussi imposer les poses).
Entrées : images extraites des mp4 (`video360.py` sait les lire), masque du fond noir.
Sorties : nuage gaussien → maillage par extraction de surface (SuGaR ou marching cubes sur l'opacité) → `enveloppe_gs.ply` texturée, à
recaler sur le volume par la silhouette (le code de `video360.py` fait déjà ce recalage 2D).
Temps : 20 à 40 min par stade sur GPU. Risque : perspective inconnue exactement ; l'échelle vient toujours de la CRL.

### P4 — Fusion et restauration des volumes vidéo sur GPU

Pourquoi : la fusion de Fourier des trois piles est aujourd'hui limitée par la RAM et le CPU (≈ 790 voxels de hauteur, 8 Go par stade).
Avec CuPy (FFT GPU) et 32 Go de RAM, on peut monter à 1000–1200 voxels, ajouter une déconvolution (PSF anisotrope de chaque pile) et un
débruitage auto-supervisé 3D (Noise2Void) qui efface les artefacts de compression vidéo. Meilleurs volumes = meilleure P1.
Temps : quelques minutes par stade une fois codé. Gain plus modeste que P1–P3, mais il conditionne leur qualité : à faire tôt si P1 montre
que les artefacts gênent.

### P5 — Rendu : animations 4K du morphing sur GPU

Cycles ou EEVEE sur la 3080 : les séquences de `blender_render_anim.py` passent de heures à minutes. À mettre en place dès le premier jour,
c'est le moyen de contrôler visuellement tout le reste (planches, vidéos de comparaison Hikspoors / vidéo / IA).

### À vérifier avant de s'engager (deux pistes de données)

- Les **modèles OPT du HDBR atlas** (tomographie optique, CS12→CS23, annotés, licence CC BY-NC-SA) : si les volumes ou les maillages se
  téléchargent, ce sont de vrais volumes 3D isotropes couvrant CS19, CS21, CS22 — la meilleure source possible pour P1 et P2 sur les stades
  manquants. Le site est bloqué depuis le cloud : à regarder depuis le PC (`hdbratlas.org/carnegie_stages.html`).
- L'article « 4D Human Embryonic Brain Atlas » (arXiv 2503.07177) décrit exactement la méthode de P2 pour le cerveau : lire et réutiliser.

## Préparation de la clé USB

Logiciels (Windows) : pilote NVIDIA récent + CUDA 12.x ; Python 3.11 portable (WinPython ou conda-pack) avec `torch` (cu12), `monai`,
`nnunetv2`, `SimpleITK`, `antspyx`, `cupy-cuda12x`, `scikit-image`, `trimesh`, `fast_simplification`, `pynrrd`, `opencv-python`, `pillow` ;
Blender 4.x/5.x ; 3D Slicer (contrôle des labels) ; Postshot ou nerfstudio + COLMAP ; MeshLab / CloudCompare ; git.
Données : le dépôt (branche `claude/trusting-clarke-3imqlz`) ; par stade `work/dens.npy`, `work/labels.npz`, `work/axes.json`, `out/manifest.json`
et PLY (≈ 1 Go par stade) ; les mp4 de rotation ; la banque `hikspoors_maastricht/` (npz, calage, textes) ; `embryons_3D/modeles/` ;
`coupes embryos 9-23/` si P1 utilise les coupes légendées. Compter 20 à 30 Go.
Vérifications au premier démarrage : `python -c "import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))"`,
`python embryo3d/test_hikspoors_synthetique.py`, un rendu Blender GPU d'une scène existante.

## Séquence proposée pour la première semaine

Jour 1 : installation, rendu GPU (P5), extraction des images 360° et premier splatting sur CS20 (P3, stade le plus propre).
Jours 2–3 : recalage Hikspoors → volumes pour CS13, CS15, CS17, CS20 (contrôles sur coupes), voxelisation, jeu d'entraînement (P1a–b).
Jour 4 : entraînement nnU-Net (nuit), pendant ce temps P4 si les volumes le demandent.
Jour 5 : inférence, Dice contre les heuristiques et contre Hikspoors sur un stade tenu à l'écart, export `ia_*` dans `labels.npz`, maillages.
Jours 6–7 : premières paires de recalage déformable (CS13→CS14, CS19→CS20), prototype de morphing par champ, comparaison avec le Shrinkwrap.

Livrables à remonter dans l'agrégateur : par stade, un bloc `ia` dans le manifest (Dice, confiance) et un statut « appris » distinct de
« heuristique » et « externe » ; par paire de stades, le champ de déformation et le résidu ; par stade vidéo, l'enveloppe texturée.
