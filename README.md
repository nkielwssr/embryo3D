# embryo3d — reconstruction 3D segmentée d'embryons humains (CS13 → CS20, sept stades : 13, 14, 15, 16, 17, 19, 20) pour Blender / 3ds Max

Source : les 4 vidéos par stade des dossiers `CSxx_f4v` (3 piles de coupes histologiques orthogonales
transversale / frontale / sagittale + 1 vidéo de contrôle en rendu volumique).

## Principe

1. **Extraction** (`step_extract`) : chaque vidéo est lue image par image, les frames dupliquées (conversion de cadence)
   sont supprimées, le fond est normalisé et l'image inversée en *densité tissulaire* (0 = fond, 255 = très dense).
   La vidéo à fond noir est reconnue comme vidéo de contrôle et ignorée.
2. **Orientation et recalage** (`step_fuse`) : les 48 permutations/miroirs d'axes sont testées par corrélation pour
   aligner les piles, l'échelle pixel de chaque vidéo est déduite des emprises communes, puis un recalage affine
   (SimpleITK, métrique de corrélation) est fait. La pile de référence retenue est celle vis-à-vis de laquelle les deux
   autres présentent le moins de cisaillement (les vidéos sont légèrement obliques).
3. **Fusion super-résolution** : chaque pile est fine le long de son axe de coupe (700 à 1700 coupes) et grossière
   dans le plan (~424 px). Les trois volumes sont rééchantillonnés sur une grille isotrope (≈ 790 voxels sur la
   hauteur) et fusionnés dans le domaine de Fourier, chaque pile pesant davantage aux fréquences qu'elle résout.
4. **Normalisation et enveloppe** (`step_prep`) : les densités varient fortement d'un spécimen à l'autre (coloration) ;
   les percentiles du tissu (25/50/75/90/99) sont ramenés sur ceux de CS20 par une table linéaire par morceaux
   (`dens_raw.npy` garde l'original). La surface du corps est ensuite obtenue par un critère d'occlusion directionnelle
   (un voxel est « dedans » si presque toutes les directions rencontrent du tissu), robuste aux cavités et aux replis.
   Les **cavités internes** (densité ≈ 0) sont séparées par ouverture morphologique et caractérisées
   (taille, position, élongation, fraction de paroi sombre).
5. **Axes anatomiques** (`step_axes`) : gauche-droite = plus petite emprise de l'enveloppe, haut-bas = plus grande ;
   tête en haut = côté des ventricules cérébraux (cavité à paroi épaisse et dense, région de la tête) ; dorsal = côté
   de la moelle. Ces indices sont fiables aux stades âgés mais fragiles sur les embryons très recourbés (CS13-CS16) :
   chaque stade est **vérifié visuellement** (`check_volume.py`, planches `check_canonical*.png`) et corrigé au besoin par
   `flip_axis.py <work> <axe>` (0 = gauche-droite, 1 = haut-bas, 2 = avant-arrière ; noté dans `axes.json`).
   `fix_vent.py` recalcule les graines ventriculaires (cavités de la tête à paroi dense). Tous les stades sont ainsi dans
   le même repère : X = gauche→droite, Y = ventral→dorsal, Z = bas→haut (dans Blender).
6. **Segmentation** (`segment.py`) — heuristiques anatomiques sur la densité :
   - `ventricules` : système ventriculaire (cavité bordée de neuroépithélium) ;
   - `snc` : cerveau (croissance dans le tissu dense depuis la paroi des ventricules) + moelle (noyau très dense érodé,
     connecté au cerveau, ce qui coupe les ponts fins vers la peau et les membres) ;
   - `foie` : plus gros noyau très dense du tronc ;
   - `yeux`, `cristallins`, `vesicules_otiques` : paires symétriques de petites cavités (la plus antérieure = espace
     intra-oculaire → globe et cristallin ; la plus postérieure = labyrinthe otique) ;
   - `squelette_axial_cartilage`, `chondrocrane`, `cartilage_autre` : ébauches cartilagineuses = blobs pâles compacts
     (le cartilage est presque non coloré), classés par position (près du tube neural = corps vertébraux / arcs / côtes) ;
   - `cavite_pericardique` et `coeur` : cavité épaisse du thorax (hors voisinage du foie et de la moelle) et tissu
     contenu dans son enveloppe convexe ;
   - `membres`, `cordon_ombilical` : protubérances de l'enveloppe ;
   - `tube_digestif` : lumières du tronc bordées d'épithélium sombre + leur paroi ;
   - `ganglions` : amas denses compacts flanquant la moelle ;
   - `vaisseaux` : lumières tubulaires (vesselness de Frangi) hors cavités déjà attribuées.
7. **Maillages** (`meshexport.py`) : marching cubes sur masque lissé, décimation quadrique, export PLY par structure,
   GLB combiné, `manifest.json` (couleurs, collections, volumes) et NRRD (densité + labels) pour 3D Slicer.
   Échelle absolue : la hauteur de l'enveloppe est mise à la *greatest length* typique du stade (CS13 4,5 mm … CS17 11 mm, CS19 17 mm, CS20 20 mm) — à ajuster dans `CRL_MM` si l'on connaît la vraie taille du spécimen.

## Dépendances

`pip install numpy scipy scikit-image opencv-python SimpleITK trimesh fast_simplification pynrrd edt` (Python 3.11, ffmpeg non requis).
Mémoire : ~8 Go par stade ; ne pas lancer plus de 3 stades en parallèle sur 32 Go.

## Lancer

```bash
python embryo3d/pipeline.py CS20_F4V            # tout : extract,fuse,prep,axes,segment,mesh
python embryo3d/check_volume.py CS20_F4V/work   # coupes médianes du volume fusionné (contrôle)
python embryo3d/flip_axis.py CS15_f4v/work 2    # si la face regarde du mauvais côté (2 = avant-arrière)
python embryo3d/fix_vent.py CS15_f4v/work       # graines ventriculaires + planche check_vent.png
python embryo3d/pipeline.py CS19_f4v --steps segment,mesh   # rejouer une partie (segment/mesh se relancent toujours)
bash embryo3d/run_queue.sh segment,mesh CS13.f4v CS14_f4v   # plusieurs stades à la suite
python embryo3d/make_master_manifest.py master.json CS13.f4v/out CS14_f4v/out CS15_f4v/out CS16_f4v/out CS17_f4v/out CS19_f4v/out CS20_F4V/out
"C:/Program Files/Blender Foundation/Blender 5.2/blender.exe" -b -P embryo3d/blender_build_scene.py -- master.json embryon_CS13-20.blend
```

## Dans Blender

- Une collection par stade, sous-collections `Systeme nerveux`, `Organes`, `Squelette`, `Vaisseaux`, `Cavites`, `Enveloppe`.
- Matériaux Principled BSDF (enveloppe et cavités semi-transparentes).
- **Morphing** : le stade le plus âgé porte, pour chaque structure présente aux stades précédents, une shape key par stade
  (sa topologie est projetée par *Shrinkwrap + Corrective Smooth* sur le maillage du stade cible). Un slider
  `stage` (propriété de scène, 0 = CS13 … 6 = CS20) pilote les shape keys par drivers : animer `stage` fait
  grandir l'embryon. Les stades cibles sont conservés masqués comme référence.
- Pour 3ds Max : importer les PLY/GLB ; le morphing se refait avec le modificateur *Morpher* à partir des maillages
  à topologie commune exportés depuis Blender (Fichier → Exporter → FBX avec shape keys).

## Limites et pistes

- La segmentation est automatique et *heuristique* : les frontières sont des frontières de densité, pas de
  vérification histologique. Le cartilage et les cavités liquidiennes ont la même densité (≈ 0) et ne sont
  distingués que par la forme et la position ; les vaisseaux sont les lumières détectées comme tubes.
- Tout se corrige à la main dans 3D Slicer (module Segment Editor) sur `*_labels.nrrd`, puis on relance `mesh`.
- Résolution : ≈ 25 µm/voxel à CS20 (grille ≈ 350 × 800 × 550 voxels), ≈ 800 voxels sur la hauteur à tous les stades.
- Les structures « partielles » à ce jour : `tube_digestif` (lumières trop fines, parois de densité moyenne) et
  `vaisseaux` (quelques segments d'aorte / veines détectés comme tubes). Le cerveau des stades CS16 et CS17 peut être
  incomplet (ventricules mal séparés du réseau de mésenchyme lâche) : à revoir dans Slicer si besoin.

## Retouches in-place (sans re-segmentation) et fusion des autres systèmes

Chaque script réécrit `work/labels.npz` (copie de sauvegarde `labels_avant_*.npz`), puis régénère PLY/GLB/manifest/NRRD/scène/rendu du stade.

```bash
python embryo3d/fix_snc_band.py CS13.f4v                 # SNC confiné à la bande médiane sous la tête (retire somites/peau)
python embryo3d/fix_snc_tube.py CS15_f4v 110 2 0.10 0.33 0.72   # prolonge le SNC le long du tube neural (seuil, érosion, bande, tête, queue)
python embryo3d/fix_pairs.py CS13.f4v yeux 3 67          # yeux depuis 2 labels de cavités (cavo_lab) repérés visuellement
python embryo3d/fix_pairs.py CS15_f4v vesicules_otiques 4 286
python embryo3d/fix_lens_eyes.py CS16_f4v [--diag]       # yeux par le cristallin (boules denses symétriques de la tête)
python embryo3d/axial.py CS15_f4v --seuil=70 [--diag]    # squelette axial : somites (CS13-16) / arcs-côtes (CS17-20) + notochorde
python embryo3d/fusion_systemes.py CS20_F4V [--sans-scene]  # importe cœur/digestif/aorte de la session « VHE » (work/cardio, work/digestif)
python embryo3d/drop_label.py CS13.f4v cristallins       # retire une structure
bash embryo3d/build_all.sh                               # scènes par stade + scène maître + rendus
```

Structures ajoutées par ces étapes : `somites`, `notochorde`, `coeur_detoure`, `myocarde`, `cavites_cardiaques`,
`digestif_pharynx|oesophage|estomac|duodenum|intestin_moyen|intestin_posterieur`, `vaisseaux_aorte` (+ `vaisseaux_cardinales|ombilicaux|vitellins` quand tracés).
Coordination multi-sessions : seule cette chaîne écrit `labels.npz`/`manifest.json` ; la session « Reconstruction 3D depuis vidéos VHE » produit
`work/cardio/*.npz`, `work/digestif/*.npz` ; la session « Agrégateur » lit les manifestes et publie `embryons_3D/agregateur.html` (états dans `taches_etat.json`).

## Topographie par région (topographie_axe.py, 24/09/2026 — session « Extraction squelette »)
`python embryo3d/topographie_axe.py <dossier> [--guide auto|colonne|snc] [--lisser 15] [--demi-mm 1.5]` — lecture seule de
work/labels.npz + work/dens.npy, sorties dans `<dossier>/out/topographie/` (manifest_topographie.json au format lire_recon + bloc
`controle`, axe_vertebral.json, etages_vertebraux.json/.png, zones_axe.json, PLY axe + anneaux d'étages).
Axe = chemin de moindre coût dans le masque colonne (cartilage axial ∪ notochorde ∪ somites) ou, si ces labels couvrent < 60 % du
tube neural ou stade somitique, dans la moelle (snc ∪ ventricules rempli, petits blobs compacts par coupe = moelle, cerveau exclu) ;
spline lissée, abscisse curviligne s. Étages = chaîne segmentaire la plus périodique autour de l'axe (corps vertébraux, ganglions ou
somites), période par autocorrélation, frontières = minima. Zones cou/thorax/lombaire/pelvis par étages si fiable, sinon fractions
du CS23 (marquées « approximatif »). Bilan 24/09 : CS20 et CS17 fiables, CS14 à confirmer, CS13/15/16/19 non fiables (labels
incomplets). Ne jamais afficher le label snc (décision utilisateur) : usage interne uniquement.
Repère : les PLY de topographie sont translatés de −center_voxel comme ceux du pipeline. L'axe affiché est décalé ventralement
sur les corps vertébraux (le guide moelle/arcs est conservé dans axe_vertebral.json) ; le nommage des étages est ancré sur la
racine du membre supérieur (= C-7 / somite 11), les étages en amont de C-1 sont nommés « occipital k ».

## SNC fin, squelette axial recruté, coordination (24/09)

- `cns_fine.py <stade> [--diag] [--tipzone=a,b | --tip=bottom | --cuts=a,b,c]` : découpe 'snc'/'ventricules' en ventricule_prosencephale|mesencephale|rhombencephale,
  canal_central, prosencephale, mesencephale, rhombencephale, moelle, meninges_mesenchyme_cranien, epiderme_cranien (axe géodésique du système
  ventriculaire ; mésencéphale = courbure max de l'axe ; fin d'encéphale = chute de section). Réglages retenus : CS13/14/17/19/20 défaut ;
  CS16 --tipzone=0.2,0.45 ; CS15 --tip=bottom (ventricules fragmentés).
- `axial_recrutement.py <stade> [--diag]` : lit out/topographie/{axe_vertebral,etages_vertebraux,tube_neural}.json (session « Extraction squelette »)
  et écrit corps_vertebraux, arcs_neuraux, cotes (CS17+), notochorde affinée (vidée si douteuse), moelle_rachidienne (sous occipital-1/C-1),
  squelette_axial_cartilage = union ; instances par étage dans work/vertebres_niveaux.npz(+json) et PLY par niveau dans out/vertebres/.
  Stades jeunes : tissu = label 'somites'. Ne pas lancer deux scripts qui écrivent le même labels.npz en parallèle (course → labels perdus).
- Manifests : champ 'confiance' (bonne/moyenne/faible) par structure (défauts dans meshexport.CONFIANCE, surcharge work/confiance.json) ;
  l'agrégateur masque les pièces 'faible' en 3D. Brouillons masqués dans Blender dès qu'un stade a son SNC fin : snc, ganglions, cavite_pericardique.
- Mémoire : ≤ 3 scripts lourds en parallèle ; `queue_seq.sh` (N=…) enchaîne des commandes quand il reste ≤ N python actifs.
Vidéo 360° (`python embryo3d/video360.py <dossier>`, 24/09/2026) : lit le rendu volumique en rotation (mp4 à fond noir), recale la
vue de profil (2D, similitude à échelle fixée par les aires, miroir) sur la projection sagittale de l'enveloppe, projette les étages
de l'axe sur la vidéo et prolonge la chaîne jusqu'au bout de la queue (bout = géodésique de l'enveloppe depuis la tête, membres
retirés à partir de CS17) au pas de l'axe, calé sur la vidéo, à la profondeur du dernier étage. Sorties out/topographie/video360/
(somites_video.json, <CS>_somites_video.ply, profil_somites.png, controle_recalage.png) + bloc controle.video360 du manifeste.
Piège : le contour vidéo a une périodicité moitié du pas des somites.
Marques manuelles (`python embryo3d/marques_video.py <dossier> <image>`) : image de la vidéo 360° annotée par l'utilisateur
(un point de couleur saturée par étage, vue de profil ; image redimensionnée acceptée) → image vidéo retrouvée par corrélation,
recalage 2D sur l'enveloppe, étages en 3D dans le plan médian (etages_manuels.json, billes PLY, contrôle) + bloc/segment
`etages_manuels` dans le manifeste. Fait pour CS15 (24 étages) le 24/09/2026.
Marques manuelles, compléments : `--zoom image2` fusionne un recadrage agrandi (recherche de motif) ; recadrages acceptés pour l'image
principale ; l'angle de projection du volume est cherché automatiquement quand le recalage direct est mauvais (CS13 : 165°).


Marques manuelles, vue de dos (24/09) : le type de vue est déduit de l'angle de l'image vidéo par rapport à la vue de profil (la plus large) ; à 50-130° on est de dos ou de ventre et le recalage se fait sur la projection frontale de l'enveloppe (`env.any(axis=2).T`, marques → ax0/ax1), la profondeur AP étant prise sur l'axe automatique à la même hauteur (sinon médiane de l'enveloppe). Pas de filtre de profondeur en vue de dos (la chaîne est médiane, donc « profonde » dans la silhouette). Exemple : `CS14_f4v/Sans titre-4.png` (19 somites de dos, image vidéo 40) → série `etages_manuels` (nom par défaut), à 0,05-0,36 mm de l'axe automatique.

Marques manuelles, options (24/09) : `--prolonger N` ajoute N étages extrapolés au bout de la chaîne (côté queue) au pas médian des marques, dans la direction des 3 dernières, ramenés dans la silhouette vidéo ; ils portent `extrapole: true` dans le JSON et une bille creuse sur le contrôle. `--reference` marque l'image comme vue de référence de repérage/calage (`vue_reference_calage` dans le bloc de contrôle). CS15 : `point verts.png` (image vidéo 20, 24 marques confirmées par l'utilisateur) + 1 extrapolé = 25 étages, vue de référence.

Marques manuelles, vues multiples et fusion (24/09) : une image prise du côté opposé au profil de référence (> 130° de la vue de profil) est recalée par miroir horizontal de la silhouette (l'ICP force une rotation propre et partait tête-bêche), le miroir étant composé dans A. Fusion 3D de plusieurs séries : `python embryo3d/marques_video.py <dossier> --fusion serie_principale serie2 ... --nom etages_manuels` : la principale (sans ses extrapolés) est prolongée par chaque complément dans le sens qui s'éloigne de son dernier étage, les points déjà couverts sont ignorés ; JSON avec `source` par étage, contrôle une couleur par série, manifest mis à jour (les séries partielles restent sur disque, hors manifest). CS15 : point verts.png (image 20, 24) + Sans titre-5.png (recadrage queue, image 108, 16 nouveaux, raccord 0,21 mm) = 40 étages, les 8 derniers serrés sur la pointe de la queue (coccyx).

Topographie sur les étages manuels (`python embryo3d/topographie_axe.py <dossier> --guide niveaux`, 24/09) : si out/topographie/video360/etages_manuels.json existe, l'axe vertébral est la spline passant par les étages manuels, les étages sont ces points (frontières à mi-chemin), nommage par ancrage membre supérieur, zones par région ; work/vertebres_niveaux.npz (instances par niveau de la session base) ne sert qu'au contrôle d'appariement ('niveaux_corps' lu en priorité quand la base le fournit ; 'niveaux' = corps ∪ arcs ∪ côtes, identifiants non ordonnés : ne pas l'utiliser tel quel comme axe). Si la série manuelle vient d'une projection tournée (projection_theta_deg ≠ 0, CS13/CS16), le X de chaque étage est recentré sur la médiane de l'enveloppe à sa hauteur (fenêtre AP ± 0,5 mm) : sans cela l'axe sort du corps (écart latéral 2,7 mm à CS16). Sans série manuelle : chaîne plus proche voisin des centroïdes. Le manifest réécrit conserve les segments et blocs de contrôle externes (etages_manuels, somites_video, reper…). Résultats : CS19 31 étages (occipital 2 → S-5, un corps par étage sur la bande frontale), CS15 40, CS16 42, CS13 32 ; nommage somitique incertain (±2 à 6, ancrage tiré par une composante crâniale à CS16). Le chemin de moindre coût sur le label corps_vertebraux (--guide colonne) donne des axes tronqués (instances séparées).

Niveaux à nombre imposé (24/09, consigne utilisateur « on sait combien il y a de vertèbres ») : topographie_axe.py impose désormais 37 niveaux identiques à tous les stades (occipital 4→1, C-1…C-7, T-1…T-12, L-1…L-5, S-1…S-5, Co-1…Co-4) et ne compte plus les étages détectés (conservés dans etages_detectes.json / etages_detectes.png comme contrôle). Repères le long de l'axe : centre de C-7 = racine du membre supérieur (label membres, composantes les plus crâniales, celles qui se projettent sur une extrémité de l'axe sont rejetées), centre de S-1 = racine du membre inférieur (composante la plus proche de C-7 + 18 pas détectés, sinon estimée), pointe de la queue = bout géodésique de l'enveloppe (distance droite depuis la fin de l'axe ; ne pas utiliser le prolongement vidéo, aberrant à CS19). Pas cervical/occipital = pas détecté (somites formés), pas thoraco-lombaire = (S-1 − C-7)/18, pas sacro-coccygien = reste jusqu'à la pointe ; puis calage de chaque frontière sur la frontière détectée la plus proche (< 0,3 pas) sans changer le nombre ; niveaux hors axe (au-dessus du début ou au-delà de la fin) = centres extrapolés le long de la tangente, `hors_axe: true`. Sorties : niveaux_imposes.json, etages_vertebraux.json (= imposés), video360/etages_imposes.json (format lu par le recrutement de la base : 37 noms courts occ4…Co4, xyz_mm_pipeline des centres), anneaux PLY, zones par région, manifest (controle.etages = imposés, controle.etages_detectes, controle.niveaux_imposes). Aux stades jeunes (CS13-14) les niveaux caudaux sont prédits dans la queue non segmentée (pas comprimé) : c'est voulu, même nombre partout pour le morphing niveau par niveau.

Extrémités ancrées et occipitaux (24/09, 2e retour utilisateur « pas de corps vertébraux au niveau occipital ») : occipital 1→4 numérotés du crâne vers le cou, `corps: false` (somites intégrés à la base du crâne, à ne pas recruter ni morpher en vertèbres) ; 33 niveaux avec corps. L'axe est prolongé aux deux bouts (`prolonger_axe`) : en crânial, marche dans l'enveloppe tirée vers le centroïde du snc complet (max 0,35 arc) ; en caudal, chemin de moindre coût dans l'enveloppe (coût 1 dans le snc) jusqu'au bout géodésique de la queue, seulement si ce bout prolonge l'axe (point d'axe le plus proche à ≥ 0,85 arc) — une queue enroulée qui touche le corps court-circuite la géodésique (CS15) et la marche libre rentre dans le corps. Ancre occipitale = bord caudal des vésicules otiques par projection curviligne ((V − Q(s))·T(s) = 0 sur la partie crâniale, centre + rayon équivalent), acceptée entre 2 et 8 pas cervicaux au-dessus de C-1 (sinon pas cervical ; échoue à CS17/CS19 où la tête fléchie donne un croisement au niveau du cou). Co-4 finit à la fin de l'axe prolongé, jamais au-delà (pas caudal comprimé) ; membre inférieur écarté s'il ne laisse pas 0,4 pas derrière lui (CS13 : bout de queue étiqueté membres). `hors_axe` = false partout, `sur_prolongement` = centre sur le prolongement. Réserve : à CS15 les deux ancres laissent 2,2 mm aux 4 occipitaux (2 × le pas cervical) → écart ~4 somites entre otique et membre supérieur, à arbitrer.

Bandes de la colonne déroulée (24/09) : le repère (dorsal / latéral) est transporté parallèlement le long de l'axe dans `profils` ; l'ancien repère global orthogonalisé à T devenait singulier quand la tangente passait par la direction AP (queue enroulée : cassure horizontale des deux panneaux à 17,1 mm sur CS20). La planche des niveaux imposés est tracée sur l'axe prolongé (0 → bout de la queue, rééchantillonné à 0,02 mm) ; la planche des étages détectés reste sur l'axe d'origine.

Repères d'extrémités pointés par l'utilisateur (24/09, `embryo3d/reperes_video.py <montage.png> --grille 4x2 --stades CS13,CS14,,CS17,CS15,CS16,CS19,CS20`) : montage de captures vidéo 360°, une case par stade, cercle vert = vésicule otique, point vert = pointe de la queue (le plus gros blob = vésicule). Chaque case est recalée sur sa vidéo comme les étages manuels (profil / profil opposé miroir / dos), les deux marques sont ramenées dans le plan médian → video360/reperes_utilisateur.json + controle_reperes_utilisateur.png. topographie_axe.py les lit en priorité : vésicule otique = point d'axe prolongé le plus proche (partie crâniale) → ancre d'occipital 1 ; pointe de queue = fin de l'axe (prolongement droit jusqu'à la marque si elle dépasse, troncature si elle est en deçà) → fin de Co-4. Sources notées dans controle.niveaux_imposes.reperes (source_vesicule_otique / source_pointe_queue = 'utilisateur (queue-VO)'). Piège du 24/09 16:42 : la base a réécrit CS20 labels.npz (moelle_rachidienne 7,6-18 mm, corps fragmentés) et le guide auto/snc ne redonnait plus l'axe validé ; CS20 restauré depuis les anneaux publiés (30 centres → video360/etages_manuels.json, --guide niveaux) ; `--label-tube` force le label du guide ; le guide tube choisit le premier label couvrant ≥ 85 % de l'étendue max.

Vérification des repères utilisateur (`reperes_video.py --verifier`, après la topographie) : écart marque VO ↔ label vesicules_otiques (plan sagittal + centroïde), marque queue ↔ bout géodésique de l'enveloppe / bord / tube neural, indices d'autres labels (yeux, cristallins, membres, cordon), champs `valide` par repère (seuil 0,3 mm ; une validation manuelle `validation: ...` posée dans le JSON est conservée), ligne controle.reperes_utilisateur dans le manifest, image controle_reperes_utilisateur.png, tableau embryons_3D/reperes_utilisateur_tableau.json. topographie_axe.py n'utilise que les repères valides ; la pointe utilisateur devient la cible du chemin caudal (`cible_caud`), sans troncature. Résultat du 24/09 : CS13 les deux, CS14/CS15 la queue ; CS17/19/20 cercles sur l'œil et points sur les membres inférieurs (non pris). Une vue unique ne donne pas la profondeur : demander une 2e vue à 90°.

Erratum 24/09 17:20 : le cercle vert de queue-VO.png est la vésicule OPTIQUE (œil), pas otique. reperes_video.py compare donc la marque au label `yeux` (clé `vesicule_optique`), et topographie_axe.py ne s'en sert que comme repère crânien (`oeil_utilisateur_s_mm`), jamais comme ancre des occipitaux ; seule la pointe de queue reste une ancre (fin de Co-4). Écarts œil ↔ label yeux : CS13 0,01 mm, CS20 0,05, CS15 0,77, CS16 0,81, CS14 1,47, CS19 2,11, CS17 3,14.

Triangulation de l'œil sur plusieurs vues (`reperes_video.py <montage> --vue2 --angle A` par montage, puis `--resoudre`, puis `--verifier`) : chaque vue supplémentaire est recalée sur la projection de l'enveloppe tournée de theta (cherché autour de +A et de +angle de l'image vidéo, ±30°, pas 5°, direct/miroir ; chercher seulement autour de +A, la silhouette miroir à −A est ambiguë) et donne une équation linéaire −sin(−θ)(x−c0)+cos(−θ)(z−c2)=z'obs ; la vue de profil donne z ; moindres carrés → xyz, résidu, écart d'ax1 entre vues, écart 3D au label yeux. Convention : les degrés du viewer de l'utilisateur = theta (même sens, offset 0-15°). Limite : les vidéos 360° sont en perspective → résidus 0,35-1,75 mm (24/09, vues 302° et 244°) ; la 3D n'est prise pour le verdict que si le résidu ≤ 0,3 mm, sinon l'écart sagittal reste le critère.

## Animation de la scène maître
`embryons_3D/embryons_CS13-CS20_morph.blend` : la propriété de scène `stage` (0 = CS13 … 6 = CS20) est animée par `blender_build_scene.py` : palier de 12 images sur chaque stade puis transition linéaire de 48 images, 24 i/s, images 1 à 372. Appuyer sur Lecture fait défiler les 7 stades ; réglage manuel dans Propriétés > Scène > Propriétés personnalisées > stage, ou dans l'éditeur graphique. `blender_render_morph.py` efface l'animation en mémoire pour ses rendus fixes.

## Vertèbres animées (nommage imposé)
Si `out/topographie/video360/etages_imposes.json` existe (session Extraction squelette : même liste de noms de niveaux pour les 7 stades, ex. C1…Co3), `axial_recrutement.py` le lit en priorité (mode manuel, pas de niveaux de queue) et exporte par niveau un PLY complet (`vertebre_<nom>`) et un PLY corps seul (`corps_<nom>`) dans `out/vertebres/`. `blender_build_scene.py` morphe alors les corps niveau par niveau dans la scène maître (collection « Vertèbres (morph) »), uniquement si les noms sont identiques sur tous les stades ; sinon les vertèbres restent masquées comme toute structure non animée.

## Tubes morphables : pharynx, tube digestif, aortes et mésos (26/09)
`python embryo3d/tubes_morph.py [--portee 3.0] [--colonnes 4] [--sans-mesos]` lit les lignes centrales tracées par points de passage
(work/digestif/chemins.json, work/cardio/vaisseaux_chemins.json ; estomac « sac » : ligne et rayon équivalent tirés du masque) et écrit
`embryons_3D/tubes_morph.npz` + `tubes_morph.json` : une forme par stade et par structure, à topologie commune, en mm dans le repère du pipeline.
- **Tubes** (64 anneaux × 16 sommets) : `pharynx`, `oesophage`, `estomac`, `duodenum`, `intestin_moyen`, `intestin_posterieur`,
  `aorte_dorsale_gauche|droite`, `aorte_commune`. Un segment absent à un stade est réduit à un point sur son raccord (début/fin du voisin) et
  « pousse » pendant le morphing. Le pharynx est un segment digestif comme les autres : il apparaît dès qu'un segment `pharynx` (alias
  `intestin_pharyngien`) est tracé dans `digestif_points/<CS>.json` (points relevés sur `digestif_planches.py`, du fond de la cavité buccale à
  l'origine de l'œsophage derrière la voie de sortie du cœur ; en amont de l'œsophage, sens haut → bas), puis `digestif_build.py`,
  `digestif_export.py` (PLY `<CS>_pharynx.ply`, confiance moyenne) et `fusion_systemes.py` (label `digestif_pharynx`). Aucun pharynx n'est
  tracé à ce jour : tant qu'il manque, il reste réduit au début de l'œsophage.
- **Mésos dorsaux** (nappes 64 lignes × 4 colonnes) : `meso_oesophage`, `mesogastre_dorsal`, `mesoduodenum`, `mesentere`, `mesocolon_dorsal`,
  tendus entre le bord dorsal de chaque anneau digestif et le bord ventral de l'axe aortique du stade (milieu des aortes dorsales paires là où
  elles se font face à < 1 mm, puis aorte commune). L'attache est la projection au plus proche, rendue monotone le long de l'aorte par régression
  isotonique sur toute la chaîne digestive (pas de croisement ; l'anse de l'intestin moyen donne un éventail depuis sa racine). Une ligne dont
  l'anneau est à plus de `--portee` mm de l'aorte, ou dont le pied tombe au-delà d'une extrémité tracée, est réduite à son bord digestif (largeur
  nulle) ; sans aorte au stade le méso reste plaqué sur le tube et s'ouvre pendant le morphing. Conséquence des tracés actuels : CS16 (aortes
  paires au niveau de l'œsophage) ne tend que le méso-œsophage, CS17 (tronc commun court) surtout le mésentère et le mésoduodénum.
  `tubes_morph.json` → `mesos.attaches` donne par méso et par stade le nombre de lignes tendues, l'intervalle d'abscisse aortique et la largeur
  médiane. Le pharynx n'a pas de méso ; les mésos ventraux (mésogastre ventral, ligament falciforme vers le foie) ne sont pas construits.
- **Contrôle 2D sans Blender** : `python embryo3d/tubes_morph_planche.py embryons_3D/rendus/tubes_planche.png [0 0.5 1 … 6] [--px-mm 30]`
  (profil + face à échelle commune, valeurs fractionnaires = interpolation entre stades, structures absentes listées).
- **Scène Blender** : `blender -b -P embryo3d/tubes_morph_blender.py -- embryons_3D/tubes_morph.blend [embryons_3D/rendus/tubes] [--epaisseur 0.05]`
  → collection « Tubes morphables » (sous-collections Tube digestif / Aortes / Mésos dorsaux), une shape key par stade pilotée par `stage`, mêmes
  keyframes que la scène maître ; les nappes portent un Solidify (0,05 mm) et un matériau translucide double face. La scène maître garde son bloc
  tubes désactivé (`TUBES_VHE = False`, décision du 24/09) ; la collection se lie dans une autre scène par `bpy.data.libraries.load`.
- Le site (`viewer_3dh.py`) lit les formes par structure (`liste` dans les métadonnées) : tubes fermés, nappes ouvertes translucides sous le
  système « Tube digestif » ; couleurs lues dans `tubes_morph.json` (`structures.<nom>.couleur`, table `COULEURS` de `tubes_morph.py`).

## Arcs aortiques, sac aortique et poches pharyngiennes (26/09, dynamique de Rana et al. 2014)
La figure CS11 → CS14 de Rana, Sizarov, Christoffels & Moorman (2014, *Am J Med Genet A* 164A:1372) montre ce qui manque à nos jeunes stades :
les artères des arcs pharyngiens entre le sac aortique et les aortes dorsales, leur apparition et leur régression, et les poches pharyngiennes qui
les séparent. La chaîne `tubes_morph.py` les prend en charge comme tubes supplémentaires, dès qu'ils sont tracés :
- **Structures** : `arc_aortique_{1,2,3,4,6}_{gauche,droite}` (1 mandibulaire, 2 hyoïdien, 3 carotidien, 4 aortique, 6 pulmonaire), `sac_aortique`,
  `tronc_arteriel` (lumière de la voie de sortie du cœur), `poche_pharyngienne_{1..4}_{gauche,droite}`. Couleurs de la légende de la figure
  (beige, jaune, vert, cyan, magenta ; sac orange, voie de sortie bleu-violet, aorte dorsale rouge) dans `tubes_morph.py` → `COULEURS`.
- **Tracé** : arcs, sac et tronc dans `vaisseaux_points/<CS>.json` puis `digestif_build.py <work> embryo3d/vaisseaux_points/<CS>.json --sortie cardio`
  (comme les aortes ; `vaisseaux_export.py` exporte les PLY, `fusion_systemes.py` les réunit dans le label `vaisseaux_arcs_aortiques`) ; poches dans
  `digestif_points/<CS>.json` (label fusionné `digestif_poches_pharyngiennes`). Sens : arcs du sac aortique (ventral, début) vers l'aorte dorsale
  (dorsal, fin) ; sac et tronc du cœur vers les arcs ; poches de la lumière du pharynx vers leur fond latéral. Repères dans les volumes : les arcs
  sont des lumières pâles courtes et obliques dans le mésenchyme des arcs pharyngiens, latérales au pharynx, entre deux poches ; leur origine
  ventrale est le sac aortique, juste en avant de la voie de sortie du cœur, leur terminaison dorsale rejoint l'aorte dorsale du même côté ; les
  poches sont les prolongements latéraux pâles de la lumière pharyngienne (`digestif_lumieres.py` les détecte comme composantes distinctes).
- **Relevé semi-automatique des arcs et du sac** (une fois les aortes dorsales tracées) :
  `python embryo3d/arcs_candidats.py CS13.f4v/work CS13 [--sac x,s,y] [--pas 2]` cherche les départs de branches pâles dans une coque autour de la
  partie crâniale de chaque aorte dorsale (côté ventral ou latéral ; intersegmentaires dorsales écartées), puis, aorte bloquée, le chemin de moindre
  coût (lumière pâle, favorisant la ligne centrale) depuis une graine dans le sac aortique (défaut : haut des cavités cardiaques) jusqu'à chaque
  départ : ce chemin ne peut passer que par la branche elle-même, c'est l'arc ; le tronc commun des chemins donne le sac. Seuil de lumière calé sur
  la densité des aortes tracées ; exclus : hors enveloppe, cavité péricardique, veines et tube digestif déjà tracés. Sorties :
  `vaisseaux_points/<CS>_arcs_proposes.json` (format de `<CS>.json`, bloc `qualite` par segment : fraction pâle, plus longue traversée de paroi,
  jonction aortique) et `<work>/cardio/arcs_candidats.png` (profil de chaque côté, face, bandeau récapitulatif, graduations en voxels). Arcs nommés
  d'après le calendrier quand le nombre de branches propres d'un côté égale le nombre attendu, sinon `arc_candidat_<i>_<côté>` ; un chemin qui
  force une paroi (> 3 voxels) devient `arc_candidat_douteux_<i>_<côté>`. **Validation obligatoire** sur la planche : renommer ou supprimer les
  candidats, couper le sac au bout de la voie de sortie (le reste = `tronc_arteriel`), corriger des points au besoin, puis
  `python embryo3d/arcs_candidats.py fusionner embryo3d/vaisseaux_points/CS13_arcs_proposes.json` (recopie `arc_aortique_*`, `sac_aortique`,
  `tronc_arteriel` dans `vaisseaux_points/CS13.json`, sans les blocs `qualite`) et `digestif_build.py … --sortie cardio`.
- **Relevé semi-automatique du pharynx et des poches** (une fois l'œsophage tracé ; à faire avant les arcs, dont il améliore les exclusions) :
  `python embryo3d/pharynx_candidats.py CS13.f4v/work CS13 [--pas 2]` part du début (crânial) de l'œsophage recalé (`digestif/chemins.json`) et
  prend la lumière pâle connexe au-dessus (seuil calé sur l'œsophage ; exclus : hors enveloppe, cœur, cavité péricardique, vaisseaux, ventricules,
  yeux, vésicules otiques). Pharynx = chemin de moindre coût dans cette lumière, centré (distance aux parois + écart au plan médian, pris au milieu
  des aortes dorsales, car le pharynx est aplati) jusqu'au point le plus éloigné sur la ligne médiane (fond de la cavité buccale), écrit de haut en
  bas et fini sur le premier point de l'œsophage. Poches = maxima, de chaque côté, de l'extension latérale de la lumière le long du pharynx
  (au-dessus de la demi-largeur médiane du pharynx), chacune de l'axe du pharynx à son fond latéral ; numérotées d'après le calendrier comme les
  arcs, sinon `poche_candidat_<i>_<côté>` ; pas de recherche aux stades sans poche attendue. Sorties : `digestif_points/<CS>_pharynx_proposes.json`
  et `<work>/digestif/pharynx_candidats.png` (face, profil médian, contour de la lumière retenue). Validation sur la planche, puis
  `python embryo3d/pharynx_candidats.py fusionner embryo3d/digestif_points/CS13_pharynx_proposes.json` et
  `digestif_build.py <work> embryo3d/digestif_points/CS13.json`.
- **Tests sans données** : `python embryo3d/test_tubes_morph_synthetique.py` (chaîne `tubes_morph` + planche : présences, raccords des absents,
  calendrier, mésos) et `python embryo3d/test_arcs_candidats_synthetique.py` (volume synthétique CS13 avec leurres : pharynx, poches, veines
  cardinales, cavité péricardique, intersegmentaire, kyste) et `python embryo3d/test_pharynx_candidats_synthetique.py` (pharynx aplati et courbé,
  4 poches à gauche et 3 à droite, bourgeon pulmonaire, ventricule cérébral) ; chacun finit par « OK : … réussi ». Dépendances : numpy, scipy,
  opencv, scikit-image.
- **Morphing** : un arc absent (pas encore formé ou régressé) est réduit sur la fin du sac aortique, sinon sur le début de l'arc voisin le plus
  proche du même côté, sinon sur l'aorte dorsale : il pousse ou se résorbe depuis son origine ventrale pendant la transition. Une poche absente est
  réduite sur le pharynx à sa hauteur attendue (1/5 … 4/5 de sa longueur). Sous-collections Blender « Arcs aortiques, sac aortique, tronc artériel »
  et « Poches pharyngiennes » ; site : arcs sous « Artères », poches sous « Tube digestif ».
- **Calendrier attendu** (`vaisseaux_points/calendrier_arcs.json`, d'après Rana 2014, Graham 2023 et Congdon 1922) : `tubes_morph.py` compare la présence
  tracée à ce calendrier et affiche, sans bloquer, les structures « attendues mais non tracées » et « tracées mais attendues absentes » (bloc
  `calendrier` de `tubes_morph.json`). Clés latéralisées prioritaires (ex. `arc_aortique_6_droite`).

  | | CS13 | CS14 | CS15 | CS16 | CS17 | CS19 | CS20 |
  |---|---|---|---|---|---|---|---|
  | arc 1 (mandibulaire) | régression | absent | absent | absent | absent | absent | absent |
  | arc 2 (hyoïdien) | régression | absent | absent | absent | absent | absent | absent |
  | arc 3 (carotidien) | présent | présent | présent | présent | présent | présent | présent |
  | arc 4 (aortique) | présent | présent | présent | présent | présent | présent | présent |
  | arc 6 (pulmonaire) | formation | formation (CS14-CS15 selon l'embryon) | présent | présent | présent (droit : régression) | gauche seul | gauche seul |
  | sac aortique, tronc artériel | présents | présents | présents (cloisonnement) | présents | présents | présents | présents |
  | poches 1-4 | 1-3 présentes, 4 formation | présentes | présentes | présentes | régression | absentes | absentes |

- **Non modélisé** : 5e arc (transitoire) ; régression des segments intermédiaires des aortes dorsales (canal carotidien entre arcs 3 et 4 vers
  CS17, aorte dorsale droite caudale vers CS19) — les aortes dorsales restent des lignes uniques, à découper en segments si l'on veut les faire
  disparaître un à un ; artères pulmonaires, carotides externes, intersegmentaires. Aucun arc ni poche n'est tracé à ce jour : la chaîne est prête,
  les points restent à relever (CS13 et CS14 d'abord, où la dynamique est la plus riche), avec `pharynx_candidats.py` pour le pharynx et les
  poches, puis `arcs_candidats.py` pour les arcs et le sac.
  Calendrier aligné le 26/09 sur le résumé de Rana 2014 (arcs 1 et 2 déjà en régression à CS13, disparus ou réduits à des reliquats à CS14) et
  sur Graham 2023 (6e arc formé entre CS14 et CS15 selon l'embryon).

## Membres séparés
`split_membres.py <stade>` sépare le label `membres` en `membre_sup_gauche`, `membre_sup_droit`, `membre_inf_gauche`, `membre_inf_droit` (2 plus grosses composantes de chaque côté du plan médian, la plus crâniale = supérieur) dans `labels.npz` ; `meshexport` les exporte (collection Membres) et la scène maître les morphe pièce par pièce. Contrôle du morphing : `planche_morph.py sortie.png [stades]` ; vidéo : `blender_render_anim.py` (séquence PNG, caméra fixe cadrée sur CS20, option `--suivre`) puis `encode_frames.py`.
