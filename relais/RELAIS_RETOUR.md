# Retour du relais local — 26/09, sur 7c7ea13

Exécuté dans un worktree frère (`Documents/Claude/_embryo3d_cloud`, HEAD détachée sur `origin/claude/trusting-clarke-3imqlz`) pour ne pas
changer de branche sous les autres sessions locales. Sorties dans le scratchpad du relais, **pas** dans `embryons_3D/modeles/`
(l'agrégateur lit ce dossier) : c'est pour cela que `verif_maillages` n'a pas tourné (ta condition `rel` hors de `embryons_3D/modeles`).

Fichiers joints dans `relais/` : `inventaire_noms.txt`, `cs13_A.log` / `cs13_B.log` / `cs13_C.log`, `rapport_CS13_auto.json`,
`rapport_CS13_calage_modeles.json`, `CS13_comparaison_A_B_C.jpg` (les trois planches controle.png l'une sous l'autre).

## 1. Test synthétique

`OK : test synthétique Hikspoors réussi` (écart d'orientation 2.97°, hauteur 4.75 mm, 24 structures). Dépendances toutes présentes.

## 2. Inventaire des noms réels

Liste complète (nom, parent, couleur `|c`, nb de sommets, teintes `|fc`) pour 9, 13, 18_NEWvalves, 20_NEWvalves, 23_NEWvalves :
`relais/inventaire_noms.txt`. Constats :

- `|p` est **identique au nom du nœud** partout : le repli « via parent » n'apporte rien.
- `|c` est en **0-1** (float32) ; `|fc` est en **0-255** (uint8).
- Original contre NEWvalves : CS18 NEWvalves **perd `gut`** (et `CCS` devient `CCS_1`) ; CS20 identiques (hors cube) ; CS23 `Asc_Ao_wall`
  devient `asc_Ao_wall`. Choisir NEWvalves par défaut fait donc perdre l'intestin à CS18.
- Pas de `liver` ni de `transverse_septum` à CS18, CS20, CS23 (ni dans les originaux).
- Noms à couvrir dans la nomenclature (vus dans ces 5 stades) : `arterial_plexus`, `venous_plexus`, `coelom`, `cardiac_jelly`,
  `neural_plate`, `somites_somitomeres` (CS9) ; `PAAs_2/3/4/6` ; `R_/L_HCC` ; `SAN` ; `DMP`, `myocardial_DMP` ; `CCS`, `CCS_1`, `CSS`
  (coquille des auteurs à CS20) ; `NCCs` ; `musc_ventricular_septum`, `myocard_outlet_septum`, `primary_/secondary_atrial_septum` ;
  `Ao_/pulm_{non_adj,parietal,septal}_leaflet`, `mitral_valve`, `tricuspid_valve`, `venous_valves` ; `asc_Ao_wall`, `pulm_trunk_wall`,
  `ascending_aorta`, `pulmonary_trunk`, `pulmonary_arteries`, `left_coronary_artery`, `coronary_arteries`, `lumen_coronary_sinus` ;
  `lumen_OFT_subaortic_part`, `lumen_OFT_subpulm_part` ; `lumen_pulmonary_vein`, `myocard_pulmonary_veins` ; `cardinal_veins`,
  `common_cardinal_veins`, `azygos_venous_system`, `cranial_veins`, `sup_caval_vein`, `inf_caval_vein` ; `L_/R_spinal_ganglia` ; `lungs`,
  `lung_lobes`, `main_branches_lungs` ; `gut` (CS20/23, 10-12 teintes).
- `_scene.json` : `nodes` est une liste de 87 nœuds `{name, kind, resource, parents:[{name, matrix_column_major}], visibility}` ; les trois
  premiers sont `Carnegie Stage 13` (groupe racine), `neural_tube`, `somites`, tous avec une **matrice identité**.
- LISEZMOI de la banque : sa table des spécimens (CS9 3709, CS10 6330, CS11 6344, CS12 8943, CS13 836, CS14 6502, CS15 721, CS16 6517,
  CS17 6520, CS18 4430, CS20 462, CS23 9226) et le format npz sont ceux que tu as déjà. Le reste : section « Chaîne »
  (`telecharger.py`, `extraire.py`, `vues_pdf.py`, `textes.py`, `controle_tous.py`).

## 3. CS13 : trois variantes (la troisième ajoutée par le relais)

| | commande | indices | planche (vue par le relais) |
|---|---|---|---|
| A | auto (`--separer-couleurs`) | score 3.0 : OK crânial ×2, OK dorsal, NON encéphale, NON gauche | **FAUX** : dos à gauche sur le profil, vue de face oblique (rotation autour de Z) |
| B | `--calage T_CS13_hikspoors_vers_vhe_mm.npy --calage-repere pipeline` | -3.0 tel quel, 1.0 avec miroir ; ATTENTION miroir | **FAUX** : tête en bas et dos à gauche, c.-à-d. rotation de 180° autour de X, pas un miroir X |
| C | `--calage M.npy --calage-repere modeles`, M = `matrice_pdf_vers_mm` de `calage/Carnegie_Stage_13.json` | 3.0 : OK crânial ×2, OK gauche, NON dorsal, NON encéphale | **JUSTE** : tête en haut, cœur ventral (à gauche du profil), somites dorsaux, vue de face étroite et symétrique |

Étendues : A 1.82 × 2.44 × 4.04 mm (rapport_hauteur_crl 0.899) ; B 1.21 × 2.79 × 4.06 ; C 1.19 × 2.88 × 4.14 (échelle du json : 1.101 µm/u,
contre 1.078 dans `hikspoors_echelles.json`).

**Les indices anatomiques ne sont pas fiables sur un embryon en C** : ils valident A (faux) et refusent « dorsal » à C (juste). Causes
probables : le tube neural fait tout le tour du C (le prosencéphale est ventral contre le cœur, la queue remonte), donc son barycentre
n'est pas dorsal et « plus large en haut » échoue ; « coeur » contient les veines cardinales (dorsales), et `^veine_cardinale` des DORSAUX
ne trouve jamais rien (voir nomenclature). Le message « ATTENTION : essayer le miroir X » de B propose la mauvaise correction : le repère de
`T_CS13…` n'est pas le repère « pipeline » supposé (Z n'y est pas crânial).

### Nomenclature : erreurs du journal CS13

- `R_/L_cardinal_vein`, `R_/L_common_card_vein` → **« coeur »** (la règle cardiaque attrape « card ») : attendu `veine_cardinale_{droite,gauche}`,
  `veine_cardinale_commune_{droite,gauche}`.
- `epicard`, `pericard`, `pericardial_reflection`, `sup_/inf_endocardial_cushion` → **« coeur »** : attendu `epicarde`, `pericarde`,
  `reflexion_pericardique`, `coussin_endocardique_{superieur,inferieur}`. Conséquence : la structure « coeur » (et son union) est surtout le
  péricarde avec les veines cardinales ; le gros ballon rouge des planches est le péricarde.
- `R_Vit_vein` + `R_Umb_vein` → une seule « veine_r_droit », `L_Vit_vein` + `L_Umb_vein` → « veine_l_gauche » (côté doublé, vitellines et
  ombilicales mélangées) : attendu `veine_vitelline_{d,g}`, `veine_ombilicale_{d,g}`. Le contrôle « gauche à +X (2 paires) » repose sur ces
  mélanges.
- `lumen_venous_sinus` → « veine_lumen » : attendu `cavite_sinus_veineux` ; `myocard_venous_sinus` → « myocarde_coeur » : attendu
  `myocarde_sinus_veineux`.
- Non reconnus : `PAAs_2/3/4` = artères des arcs pharyngiens → `arc_aortique_2/3/4` (l'indice crânial « ^arc_aortique » ne tient
  aujourd'hui que par `sac_aortique`) ; `R_/L_HCC` → `canal_hepatocardiaque_{d,g}` ; `SAN` → nœud sinusal ; `DMP` → protrusion
  mésenchymateuse dorsale.
- `loop_wire_heart_tube` (55 teintes, 4 morceaux « coeur_tube_cardiaque ») ressemble à un fil de repère de la boucle cardiaque, pas à un
  tissu : à ignorer ou à garder à part du morphing.
- `OFT_cardiac_jelly` → « coeur_voie_efferente » : plutôt `gelee_cardiaque_voie_efferente`.

## Calages déjà faits sur le PC (réponse à ta question 3)

La session locale « Carnegie Stage 13 PDF 3D » a calé presque tous les PDF sur **nos modèles de stade** :
`hikspoors_maastricht/calage/Carnegie_Stage_<n>[_NEWvalves].json`, clé `matrice_pdf_vers_mm` (4×4, unités PDF → mm dans le repère de
`embryons_3D/modeles/<CS>/` : Z crânial, Y dorsal, +X gauche, direct). ICP apparié par classe (tube neural, aortes dorsales, lumière
cardiaque, foie), 48 départs, contrôle gauche/droite par les R_/L_. C'est ce calage qui donne la variante C.

| PDF | modèle | spécimen du modèle | µm/unité | réflexion | G/D (mm) | écarts médians (µm) |
|---|---|---|---|---|---|---|
| 9 | CS09 | H712 | 0.986 | non | — | tube 43, aorte 59 |
| 10 | CS10 | 6330 | 0.798 | non | -0.48 | tube 15, aorte 13 |
| 11 | CS11 | 6784 | 0.758 | non | -0.30 | 27 / 25 / cœur 29 |
| 12 | CS12 | 8505A | **4.864** | **oui** | +0.43 | 71 / 49 / 65 |
| 13 | CS13 | 836 | 1.101 | non | -0.54 | 59 / 33 / 28 / foie 27 |
| 14 | **CS15** (pas de modèle CS14) | 3512 | 6.368 (`_14.json`) ou 5.790 (`_14_vers_CS15.json`) | non | -0.47 / -0.43 | 147 / 104 / 94 / 106 |
| 15 | CS15 | 3512 | 0.790 | non | -0.46 | 128 / 140 / 66 |
| 16 | CS16 | 6517 | 1.101 | non | -0.93 | 177 / 94 / 92 |
| 17 | CS17 | 6521 | 0.770 | non | -0.81 | 197 / 151 / 93 |
| 18 (+NEWvalves) | CS18 | 6524 | 0.727 | non | -0.81 | 211 / 161 / 126 |
| 20 (+NEWvalves) | CS20 | 462 | 0.970 | non | -0.81 | 201 / 181 / 100 |
| 23 | — | — | pas encore de calage | | | |

Travail **en cours** de cette session (ses journaux montrent une reprise, et un échec de renommage pour CS14 vers CS15) : CS12 (réflexion
et échelle ×5) et CS14 (autre unité, calé sur CS15) sont à confirmer avec elle. Les modèles de stade n'ont pas toujours le même spécimen
que Hikspoors (CS9, 11, 12, 14, 15, 17, 18) : pour ceux-là c'est un calage inter-spécimens.

**Proposition** : que `hikspoors_modele.py` prenne par défaut `calage/<nom>.json` (clé `matrice_pdf_vers_mm`, repère `modeles`) quand il
existe, et n'utilise l'orientation automatique qu'en repli, avec la planche comme seul juge.

## Autres remarques

- `publiable` vaut `true` par défaut. Règle de l'utilisateur pour harcelon.fr/3dht : seulement nos modèles achevés, jamais d'atlas retravaillé
  présenté comme nôtre. Un modèle issu de Hikspoors devrait donc être `publiable: false` par défaut.
- Point 4 **non lancé** : CS13 n'est pas bon en automatique (orientation et nomenclature). J'attends ta correction (nomenclature + calage json),
  puis je relance CS13 et, si la planche est juste, les onze autres stades.
