# Retour du relais local (2) — 26/09 soir, sur `claude/bold-wozniak-fe5bsm` (9e62c0a)

Exécuté dans des worktrees frères à HEAD détachée (`Documents/Claude/_embryo3d_cloud2` pour le code 9e62c0a). Toutes les sorties sont dans le
scratchpad du relais ; **rien n'a été écrit par le relais dans `embryons_3D/modeles/`**. Mais attention au point 7 : une autre session locale,
elle, y a écrit.

Pièces jointes dans `relais/retour2/` : `planches_1.jpg` (CS9-CS12), `planches_2.jpg` (CS14-CS17), `planches_3.jpg` (CS18, CS20, CS23),
`CS13_calage.jpg`, `CS13_auto.jpg`, `CS19_video.jpg`, les consoles complètes `cs<n>_v2.log`, `cs13_auto.log`, `cs19_video.log`,
`bilan_stades.txt` (lignes demandées pour chaque stade), `calage_etat.txt`, `rapport_CS13_calage.json`, `rapport_CS13_auto.json` et `master.json`.

## 1. Tests

- `test_hikspoors_nomenclature.py` : `105 noms, 0 erreurs` puis `OK : nomenclature Hikspoors conforme à l'inventaire du relais`.
- `test_hikspoors_synthetique.py` : `écart d'orientation 2.89°, hauteur 4.75 mm, 24 structures`, puis `OK : test synthétique Hikspoors réussi`.
- `test_hikspoors_banque.py` : `1. calage : écart foie 1.2e-07 mm ; 2. orientation auto sur l'embryon en C : 10.4° ; 49 structures`, puis
  `OK : test banque Hikspoors (embryon en C, calage json, NEWvalves, inter-stades, réflexion) réussi`.

## 2. CS13 par défaut (calage json) : JUSTE

44 structures canoniques, **aucune non reconnue**. Calage `Carnegie_Stage_13.json`, 1.101 µm/unité, repère modeles. Indices : score 9.0 (OK
crânial ×3, OK dorsal, OK gauche à +X sur 8 paires, « encéphale » NON mais non compté). Orientation automatique : **8.8° du calage**. Étendue
1.19 × 2.88 × 4.14 mm, dans les coordonnées du modèle CS13. verif_maillages : 18/48 sans erreur (maillages ouverts d'origine ; `coeur` a 12
arêtes non-manifold). Planche : identique à la variante C d'avant, avec des noms justes : tête en haut, cœur ventral, somites dorsaux, face
étroite. `pericarde` est maintenant à part.

## 3. CS13 `--sans-calage` : JUSTE

Score 9.0 (Z = ACP e1 +1, Y = indice dorsal), mêmes indices OK. Étendue 1.28 × 2.79 × 4.13 mm. Planche juste, très proche du calage : le repli
corrigé fonctionne sur le cas réel où l'ancien code échouait.

## 4. Clés des json de calage

**La banque a changé depuis mon premier retour** : la session « Carnegie Stage 13 PDF 3D » a refait tous les calages entre 17:31 et 17:57
(heure locale), `final.log` se termine par « FINI ». État actuel dans `calage_etat.txt`. Clés d'un json : `pdf`, `modele`, `specimen_modele`,
`matrice_pdf_vers_mm`, `mm_par_unite`, `um_par_unite`, `reflexion`, `echelle_initiale_um`, `longueurs` {hikspoors_unites, modele_mm},
`controle_gauche_droite_mm`, `ecarts_um` {classe: {mediane, p90}}, `classes`, `autres_solutions` [{score_um, reflexion, cote}…].

- Le modèle visé est nommé par **`modele`** (ex. `"CS15"` dans `Carnegie_Stage_14_vers_CS15.json`, avec `pdf: "Carnegie_Stage_14"`).
- La réflexion est notée par **`reflexion`** (booléen), aussi dans chaque `autres_solutions`. **CS12 n'en a plus** : `reflexion: false`,
  4.471 µm/unité (contre 4.864 avec réflexion avant), G/D -0.40 mm, écarts 44/53/57 µm.
- **`Carnegie_Stage_14.json` n'existe plus** : il ne reste que `Carnegie_Stage_14_vers_CS15.json` (6.233 µm/unité). Ton script ne le
  prend pas par défaut : CS14 est passé **en automatique** (« pas de calage … pour Carnegie_Stage_14 »).
- **`Carnegie_Stage_18_NEWvalves.json` et `_20_NEWvalves.json` n'existent plus** : tes NEWvalves ont repris le calage de l'originale, avec
  vérification du même repère (écart 0 %). Très bien.
- **CS23 a maintenant un calage** : 1.151 µm/unité, G/D -1.97 mm, écarts 539/447/167 µm, plus grands qu'aux autres stades.
- CS12 et CS14 : l'unité des PDF est bien différente (≈ 4.5 et ≈ 6 µm/unité) ; les autres stades sont entre 0.7 et 1.15.

## 5. Tous les stades (calage json par défaut) : 12 planches JUSTES

| stade | source | calage | score | auto vs calage | non reconnus | étendue (mm) | verif ok | planche |
|---|---|---|---|---|---|---|---|---|
| CS9 | Stage_9 | 0.986 µm/u | 0.0 + **ATTENTION rot. 180° Y** | 163.8° | — | 0.73 × 0.93 × 1.43 | 5/12 | juste (voir ci-dessous) |
| CS10 | Stage_10 | 0.787 | 7.0 | 16.3° | `lumen_IFT` | 0.61 × 0.50 × 1.96 | 8/27 | juste |
| CS11 | Stage_11 | 0.758 | 7.0 | 10.9° | `lumen_IFT` | 0.55 × 0.56 × 1.85 | 11/31 | juste |
| CS12 | Stage_12 | 4.471 | 9.0 | 31.7° | `lumen_IFT`, `myocard_IFT` | 1.05 × 1.57 × 2.23 | 15/41 | juste |
| CS13 | Stage_13 | 1.101 | 9.0 | 8.8° | — | 1.19 × 2.88 × 4.14 | 18/48 | juste |
| CS14 | Stage_14 | **aucun (auto)** | 10.0 | — | `Ao_swelling`, `pulm_swelling` | 2.18 × 5.30 × 6.61 (recentré) | 24/52 | juste |
| CS15 | Stage_15 | 0.784 | 7.0 | 4.3° | `Ao_swelling`, `pulm_swelling` | 2.42 × 5.14 × 6.35 | 21/49 | juste |
| CS16 | Stage_16 | 1.058 | 7.0 | 31.5° | `Ao_swelling`, `pulm_swelling` | 3.30 × 8.58 × 9.28 | 18/48 | juste |
| CS17 | Stage_17 | 0.736 | 7.0 | 15.4° | `Ao_swelling`, `pulm_swelling` | 3.37 × 7.71 × 8.76 | 21/50 | juste |
| CS18 | Stage_18_NEWvalves + `gut` de l'originale | 0.707 (originale) | 7.0 | 32.6° | — | 3.08 × 6.65 × 8.74 | 20/56 | juste, intestin présent |
| CS20 | Stage_20_NEWvalves | 0.959 (originale) | 7.0 | 15.4° | — | 3.90 × 8.12 × 14.50 | 17/56 | juste |
| CS23 | Stage_23_NEWvalves | 1.151 (originale) | 7.0 | 12.2° | — | 11.68 × 17.72 × 28.69 | 25/59 | juste |

- **CS9** : l'ATTENTION est une fausse alerte. À CS9, le septum transversum est encore **crânial** au croissant cardiaque (avant la
  bascule de la tête) : l'indice « `^myocarde_` > `^foie$|^septum_transversum` » est faux à ce stade, c'est le seul indice évalué (score 0 → 4
  retourné). La planche est juste : croissant cardiaque et arcs en haut, plaque neurale dorsale, somites au milieu, plexus veineux vitellin en
  bas. Suggestion : ne pas compter le septum transversum avant CS10-11, ou désactiver l'ATTENTION quand un seul indice est évalué.
- Non reconnus à ajouter : `lumen_IFT`, `myocard_IFT` (voie d'entrée, inflow tract : `cavite_voie_afferente`, `myocarde_voie_afferente`) ;
  `Ao_swelling`, `pulm_swelling` (bourrelets des valves artérielles : `bourrelet_valvulaire_aortique`, `…_pulmonaire`). `scale_cube___200um`
  (CS17, trois `_`) est bien ignoré.
- verif_maillages : environ 60 % des maillages sont ouverts dans le PDF d'origine (surfaces non fermées, par exemple les cavités). C'est
  normal pour du Hikspoors, mais les unions voxelisées en héritent : `coeur` a des arêtes non-manifold.

## 6. CS19_video et master

`video_vers_modele.py CS19_f4v` : 31 structures (aorte_dorsale, arcs_neuraux, canal_central, coeur, cordon_ombilical, corps_vertebraux,
duodenum, encephale, enveloppe, epiderme_cranien, estomac, foie, intestin_moyen, 4 membres + membres, mesencephale, moelle, myocarde_coeur,
oesophage, prosencephale, rhombencephale, squelette_axial, 3 ventricules + ventricules_cerebraux, vesicules_otiques, yeux). Planche : bien
orientée (tête en haut, cœur et foie ventraux, rachis dorsal).

`master_modeles.py`, lancé sur un **miroir** de `embryons_3D/modeles/` dans le scratchpad (liens vers les dossiers existants, plus mes
sorties à la place des `_hikspoors` et `_video`), via un petit lanceur qui remplace `master_modeles.MODELES` :

| CS10 | CS11 | CS12 | CS13 | CS14 | CS15 | CS16 | CS17 | CS18 | CS19 | CS20 | CS21 | CS22 | CS23 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| _hikspoors (27) | _hikspoors (31) | _hikspoors (41) | _hikspoors (48) | _hikspoors (52) | _hikspoors (49) | _hikspoors (48) | _hikspoors (50) | _hikspoors (56) | CS19_video (31) | _hikspoors (56) | CS21 (148, dossier existant) | aucun | _hikspoors (59) |

`morph` (13) : cavites_cardiaques, coeur, foie, intestin, myocarde_oreillettes, myocarde_ventricule_droit, myocarde_ventricule_gauche,
myocarde_voie_efferente, pericarde, poumons, septum_transversum, tube_neural, vesicule_vitelline.

Remarques :
- `aorte_dorsale` est présente aux 12 stades Hikspoors mais n'entre pas dans `morph` : elle n'est pas dans `MORPHABLES`. Même chose pour
  `veines`, `arteres`, `epicarde` et les cavités par chambre.
- Les sources ne se ressemblent pas d'un stade à l'autre : Hikspoors est centré sur le cœur, sans enveloppe ni cerveau ; CS19 (vidéo) et
  CS21 (dossier existant, 148 structures) sont des corps entiers. Dans la scène maître, enveloppe, encéphale, membres et yeux n'existeront
  qu'à CS19 et CS21 : le curseur les fera apparaître puis disparaître. `--fusion` (compléter par nos modèles de stade, qui sont maintenant
  superposés grâce au calage) paraît indispensable pour une scène continue.

## 7. ATTENTION : une autre session locale a écrit dans `embryons_3D/modeles/`

La session locale « Passer les 12 stades Hikspoors dans la chaîne CS10→CS23 » (Fable, worktree `embryo3d-worktrees/kind-shtern-72d9a8`,
base 06c5923 **plus des modifications non committées** de `agregateur.py`, `hikspoors_modele.py`, `hikspoors_nomenclature.json`,
`master_modeles.py`, `blender_build_scene.py`, `planche_morph.py`, `video_vers_modele.py`, `hikspoors_echelles.json`) a écrit entre 17:32 et 17:46 :
`embryons_3D/modeles/CS9_hikspoors … CS23_hikspoors` (12 dossiers) et `CS19_video`. Elle s'est arrêtée sur sa limite d'usage. Ses manifests
disent `publiable: true`, `statut: externe`, orientation « automatique (ACP + indices anatomiques) » **sans calage json**, avec son
propre nombre de structures (CS13 : 49 ; CS18 : 54 ; CS23 : 58). Le relais n'y a pas touché. L'utilisateur est prévenu : c'est à lui de
dire s'il faut les remplacer par les sorties de 9e62c0a. Tant qu'ils y sont, ils sont visibles de l'agrégateur, qui lit ce dossier.

## Réponses à tes questions

- **CS12** : la réflexion a disparu du nouveau calage (4.471 µm/u, sans réflexion) et la planche est juste : l'ancienne réflexion était donc
  un faux minimum. L'échelle ×4-5 est réelle : le PDF CS12 est dans une autre unité.
- **CS14** : il ne reste que le calage vers CS15. En automatique, la planche est juste et l'étendue fait 6.61 mm de haut (CRL typique 6.0).
  Garder la taille propre du PDF me paraît plus sûr que la taille de CS15 ; à décider avec l'utilisateur ou la session « PDF 3D ».
- **CS23** : oui, le calage est fait (voir point 4).

Point 6 : `master_modeles.py` n'a pas été lancé sur le vrai `embryons_3D/modeles/`, et je n'ai rien copié : j'attends la décision de
l'utilisateur.
