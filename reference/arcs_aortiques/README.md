# reference/arcs_aortiques — modèles 3D de référence du système des arcs aortiques (Rana et al. 2014, AMC Amsterdam)

Récupère les reconstructions 3D publiées par l'équipe d'Antoon Moorman (AMC Amsterdam), la même que pour le cœur
humain embryonnaire, et les range par stade de Carnegie, dans le format du pipeline : maillages PLY, GLB, manifeste lisible
par `blender_build_scene.py`, rapport, aperçu.

> Rana MS, Sizarov A, Christoffels VM, Moorman AFM. *Development of the human aortic arch system captured in an
> interactive three-dimensional reference model.* Am J Med Genet A 2014;164A:1372-1383. doi:10.1002/ajmg.a.35881, PMID 23613216.

Les modèles viennent de coupes histologiques sériées d'embryons (collection Carnegie et collections néerlandaises). La lumière
des vaisseaux et du pharynx y est segmentée à la main, stade par stade. Il y a un modèle par stade, de la formation des premiers
arcs jusqu'à la disposition définitive. Ils sont livrés en **PDF 3D interactifs** dans le matériel supplémentaire. Chaque structure
est colorée selon la légende de la figure : tube digestif gris, voie de sortie bleu-violet, aorte dorsale rouge, sac aortique orange,
arcs 1/2/3/4/6 beige/jaune/vert/cyan/magenta, tube neural beige clair, otocyste brun.

## 1. Récupérer les PDF (à faire sur le PC)

Les sites des éditeurs refusent les téléchargements scriptés. Dans cette session cloud, `onlinelibrary.wiley.com`, `hdbratlas.org`
et `3dembryoatlas.com` sont de plus bloqués par la politique réseau. Il faut donc passer par le navigateur :

1. Ouvrir <https://onlinelibrary.wiley.com/doi/10.1002/ajmg.a.35881>, onglet **Supporting Information**, puis télécharger chaque fichier.
2. Les déposer **sans les renommer** dans `embryo3d/reference/arcs_aortiques/pdf/` : le stade est lu entre autres dans le nom du fichier.
3. Compléments facultatifs, lus par les mêmes outils : les PDF 3D de Hikspoors et al. (Maastricht, <https://hdbratlas.org/hikspoors-pdf/JH_CS13.html>…)
   et le 3D Atlas of Human Embryology (<https://www.3dembryoatlas.com/>). Voir `sources.json`.

## 2. Extraire

```bash
python embryo3d/reference/arcs_aortiques/extraire.py                      # tous les PDF de pdf/ -> extraits/
python embryo3d/reference/arcs_aortiques/extraire.py x.pdf --lister       # inventaire des flux 3D (format, pages, vues) sans décoder
python embryo3d/reference/arcs_aortiques/extraire.py x.pdf --stade CS14   # stade imposé
"C:/Program Files/Blender Foundation/Blender 5.2/blender.exe" -b -P embryo3d/blender_build_scene.py -- embryo3d/reference/arcs_aortiques/extraits/CS14/manifest.json arcs_CS14.blend
"C:/Program Files/Blender Foundation/Blender 5.2/blender.exe" -b -P embryo3d/reference/arcs_aortiques/blender_arcs.py -- embryo3d/reference/arcs_aortiques/extraits/index.json arcs_aortiques.blend
```

Dépendances : `numpy`. Sont conseillés `pypdf` (vues, JavaScript, texte des pages, PDF chiffrés ; sans lui, un balayage brut
trouve quand même les flux 3D), `cryptography` (PDF chiffrés en AES), `trimesh` (GLB) et `opencv-python` (aperçu).

Sorties dans `extraits/` (non versionnées, comme les autres données lourdes) :

| fichier | contenu |
|---|---|
| `rapport.md` | pour chaque flux 3D : stade retenu et indice utilisé, échelle, vues, tableau objet → structure (méthode d'identification) |
| `index.json` | inventaire (PDF → flux → dossier) ; utilisé par `blender_arcs.py` |
| `<CS>/manifest.json` | format `blender_build_scene.py` (`stage`, `units`, `structures[name, file, collection, color, alpha]`), avec en plus la provenance, le nom d'origine, la méthode d'identification, la confiance, l'aire… |
| `<CS>/<CS>_<structure>[_gauche/_droite].ply` | un maillage par objet du PDF, en mm, couleurs par sommet ; les jeux de lignes (fils) sont des arêtes PLY |
| `<CS>/<CS>_arcs_aortiques.glb` | tous les maillages du stade, colorés |
| `<CS>/apercu.png` | 3 projections (face, profil, dessus) |
| `<CS>/brut/*.u3d / *.prc` | flux 3D bruts extraits du PDF (réutilisables dans d'autres outils) |
| `<CS>/vues.json`, `javascript.js`, `texte_page.txt` | vues prédéfinies du PDF (caméra, nœuds visibles), scripts, texte de la page (légende) |
| `<CS>/correspondance_proposee.json` | nom d'origine → structure proposée, à corriger puis repasser avec `--correspondance` |

## 3. Identification des structures

`legende.json` définit les structures (libellé, couleur de la légende, collection Blender, motifs de noms en anglais,
néerlandais et français). Chaque objet du PDF est identifié dans cet ordre :

1. la correspondance fournie (`--correspondance fichier.json`, format `{"nom d'origine": "arc_4"}` ou `{"…": {"structure": "arc_4", "cote": "droite"}}`) ;
2. une artère d'arc, repérée par un marqueur (PAA, arch artery, boogslagader, arc…) et un numéro (1-6, I-VI, *third*, *derde*…). Un « 5 » est signalé comme inhabituel ;
3. les motifs de noms (`ordre_de_test`) : nerfs V/VII/IX/X, canal carotidien, carotides, subclavière, sac aortique, voie de sortie, aortes, tube neural, otocyste…
4. à défaut, la couleur de la légende la plus proche (écart RGB ≤ 0,22). La confiance est alors notée « faible » et le rapport le signale.

Le côté (gauche/droite, links/rechts, L/R…) est lu à part. Il s'ajoute au nom de l'objet : `arc_3_gauche`, `aorte_dorsale_droite`.

## 4. Stade, unités, repère

- **Stade** : dans l'ordre, titre de l'annotation 3D, texte de la page, noms des vues, noms des nœuds 3D, nom du fichier (`CS14`, `stage 14`, `Carnegie stage 14`, `stadium 14`).
  Si c'est ambigu (plusieurs stades sur une même page), le dossier s'appelle `inconnu_<empreinte>` : relancer avec `--stade` ou `--stades stades.json`
  (`{"<sha1 du flux>" | "<fichier>#<k>" | "<fichier>": "CS14"}`, les empreintes figurent dans `index.json`).
- **Unités** : les modèles Amira sont en général en µm. Le facteur vers le mm est déduit de l'étendue du modèle (région pharyngienne attendue entre 0,2 et 30 mm)
  ou de l'en-tête U3D. Il est noté dans le manifeste et peut être imposé avec `--mm-par-unite`.
- **Repère** : c'est celui du fichier d'origine, translaté pour centrer la boîte englobante (`centre_retire_mm`). Il n'est **pas** recalé sur nos
  embryons : le recalage (otocystes, tube neural, aortes dorsales ↔ `vesicules_otiques`, `snc`, `vaisseaux_aorte`) reste à faire une fois les vrais modèles en main.

## 5. Formats lus (sans logiciel externe, Python pur)

- **U3D** (ECMA-363), dans `u3d.py`. C'est un port fidèle du décodeur de référence d'Intel (bibliothèque U3D, Apache 2.0, vendue avec MeshLab) :
  décodeur arithmétique, histogrammes adaptatifs, maillages CLOD (maillage de base et raffinement progressif par division de sommets),
  jeux de lignes, nœuds et groupes avec leurs matrices, modificateurs d'ombrage, shaders, matériaux, couleurs par sommet, métadonnées.
  Les vues, lumières, textures et animations sont sautées.
- **PRC** (ISO 14739-1), dans `prc.py` : conteneur, sections zlib, globales (couleurs, matériaux, styles), arbre (occurrences, pièces,
  éléments de représentation, transformations) et tessellations 3D régulières (triangles, éventails, bandes). La table des flottants fréquents
  est dans `prc_table.py`. **Non pris en charge** : la tessellation « hautement compressée » (`PRC_TYPE_TESS_3D_Compressed`) et la
  géométrie exacte (B-rep/NURBS). Si les PDF de Rana en contiennent, le rapport l'indique ; on garde les flux bruts `brut/*.prc`.
- **PDF**, dans `pdf3d.py` : annotations `/3D` et `/RichMedia`, vues `/VA` (caméra, visibilité et opacité par nœud), JavaScript
  `/OnInstantiate` et JavaScript de niveau document, texte des pages, et balayage brut de secours.

On ne connaît pas encore le format réel des PDF de Rana et al. Le groupe produisait ses PDF avec Amira puis Adobe Acrobat 9 Pro Extended
(de Boer et al., Development 2011), qui écrit du U3D ou du PRC. Les deux sont couverts ; le cas non couvert est le PRC fortement compressé.

## 6. Données de référence

- `chronologie.json` : chronologie des artères des arcs par stade. Arcs 1 et 2 présents et en régression à CS13, disparus ou réduits à des reliquats
  à CS14. Arcs 3, 4 et 6 persistants, avec le 6e formé entre CS14 et CS15 selon les embryons. Disparition du canal carotidien et de la partie dorsale du 6e droit
  vers CS17, régression plus tardive de l'aorte dorsale droite. S'y ajoutent la correspondance arc/poche/nerf (V, VII, IX, X), le principe
  plexus → vaisseau unique et le devenir définitif de chaque segment. Chaque affirmation porte sa source ; `a_verifier` marque ce qu'il faut contrôler sur les modèles.
- `legende.json` : structures, couleurs, motifs de noms. `sources.json` : références, URL, compléments (Hikspoors 2022, 3D Atlas 2016, de Boer 2011).

## 7. Tests

```bash
python -m pytest embryo3d/reference/arcs_aortiques/tests     # ou : python embryo3d/reference/arcs_aortiques/tests/test_arcs.py
```

Les fichiers de `tests/donnees/` ont été produits par les **écrivains de référence**. Pour le U3D, c'est le convertisseur IDTF→U3D d'Intel,
qui écrit des maillages 100 % progressifs, des couleurs par sommet, deux matériaux, plus de 16 383 sommets et des jeux de lignes. Pour le PRC,
c'est l'écrivain d'Asymptote. `tests/generer_donnees.py` sert à les régénérer. Les tests vérifient le nombre de sommets et de faces, l'aire, le volume
signé (orientation), le caractère fermé et orienté de la surface (une désynchronisation du flux le casse), les couleurs, les transformations imbriquées,
l'extraction PDF (avec ou sans pypdf), l'identification, les stades et unités, et l'export complet au format de `blender_build_scene.py`.
`tests/fabrique.py` écrit des U3D non compressés (maillage de base) et des PDF 3D minimaux.
