# Modèle de croissance céphalique de Fenart (référentiel vestibulaire) — synthèse pour embryo3d

Objet : disposer, pour prolonger « l'embryon en croissance par système » au-delà de CS23, d'un modèle de
croissance de la tête qui soit **géométrique** (points de repère, trajets, vitesses), **daté** (stades du fœtus à
l'adulte) et **superposable** d'un âge à l'autre grâce à un repère anatomique stable. C'est exactement ce que fournit
la base craniométrique de Raphaël Fenart (école de craniologie de Lille, méthode vestibulaire) telle qu'elle est
relue par Captier et Boë.

Les trois sources lues pour cette synthèse sont des documents fournis par l'utilisateur (01/10/2026) et **ne sont
pas versés dans le dépôt** (article et chapitre sous copyright, manuscrit non publié) :

| Source | Nature | Ce qu'on en tire |
|---|---|---|
| Pellerin P. *Creating an Atlas of Growing Skull Templates in Vestibular Orientation for Analysis of Craniofacial Malformations.* J Craniofac Surg 2026;37:16-19 | article, 4 p. | définition opérationnelle de l'orientation vestibulaire, protocole 3D Slicer, atlas 3D de 19 gabarits (181 scanners, 0 à ~18 ans) |
| Captier G., Boë L.-J. *Croissance générale et céphalique.* Chapitre de *Chirurgie plastique de l'enfant et de l'adolescent* (p. 75-93) | chapitre, 19 p. | modèle à deux unités (crâne cérébral / crâne facial), figures 13-18 des neuf âges de Fenart, profils de vitesse, organisation modulaire du crâne |
| Captier G., Subsol G., Laporte M., Boë L.-J. (+ Barbier G.) *3D ontogenetic trajectories of the human skull from prenatal to adulthood: heterochrony.* Manuscrit (titre provisoire) | brouillon .docx | matériel et méthode de la base Fenart (9 stades, 142 points), tableau des points, distances/âges à 90 %, ACP, discussion Hox/non-Hox |

Le manuscrit est un brouillon de travail : plusieurs chiffres y sont marqués « à vérifier » par les auteurs
(notamment les pourcentages de l'ACP et certaines colonnes du tableau 2). Ils sont repris ici tels quels, avec la
mention.

---

## 1. Le référentiel vestibulaire

### 1.1 Définition

Proposée par Broca (1873), formalisée par Perez (1922) puis Delattre et Fenart (1953, 1960), l'orientation
vestibulaire repose sur l'idée que la position de la tête dans l'espace est asservie à la gravité par le vestibule :
la régulation posturale maintient en moyenne le **canal semi-circulaire latéral** à l'horizontale.

- **Plan vestibulaire (horizontal)** : plan des deux canaux semi-circulaires latéraux.
- **Axe de Perez** : droite joignant les deux canaux semi-circulaires latéraux (gauche-droite).
- **Plan frontal vestibulaire** : orthogonal au plan vestibulaire, passant par l'axe de Perez.
- **Plan sagittal médian** : orthogonal aux deux précédents, par le milieu de l'axe de Perez.

L'origine du repère est donc dans le rocher, au milieu de l'axe de Perez. Avantage décisif pour la croissance et
les malformations : le labyrinthe osseux a sa taille adulte très tôt et n'est pas déplacé par la croissance de la
voûte ni de la face. C'est l'un des rares repères qui existe **du stade embryonnaire à l'adulte** (voir § 7.1).

### 1.2 Mise en œuvre

- Jusqu'en 1993, uniquement sur crâne sec, par dissection du rocher (atlas de Fenart = projections 2D
  orthogonales, huit à neuf stades).
- Depuis, *in vivo* par scanner : dissection virtuelle du labyrinthe, d'abord sous Osirix (superpositions 2D
  orthogonales), puis en vrai 3D sous 3D Slicer (Pellerin). Tutoriel vidéo de Pellerin : https://youtu.be/Bs4GWLrFQfA
- Pellerin place les **points de repère à l'intersection des trois plans vestibulaires avec la surface externe du
  crâne**, toujours dans le même ordre et en même nombre pour tous les crânes (figure 1 de l'article), ce qui permet
  ensuite l'analyse de correspondance dense (DeCA) et la superposition par transformation rigide.

### 1.3 Repère de Francfort (variante clinique)

Captier et Boë, pour la représentation clinique, ont **tourné** les données de Fenart du repère vestibulaire vers le
repère de Francfort : plan horizontal par les deux porions et les deux points orbitaires inférieurs, plan frontal
orthogonal passant par les porions, plan sagittal par les points impairs médians ; **origine au porion** (ou à
l'intersection des plans). Les deux repères ne diffèrent que d'une rotation rigide ; les trajets sont les mêmes à
cette rotation près. Le manuscrit le note comme un biais assumé (« choix d'un repère anatomique plus simple pour
l'interprétation clinique »).

---

## 2. La base craniométrique de Fenart

### 2.1 Échantillon

- Laboratoire de craniologie comparée de la Faculté libre de médecine de Lille (Delattre et Fenart 1960, Cousin
  1969, Fenart et Biecq-Sellier 2004, Fenart 2006).
- Crânes secs **entiers avec mandibule**, étude **transversale**, **sans distinction de sexe**. Chaque stade est une
  moyenne (le manuscrit indique « d'une trentaine d'individus », à préciser).
- **Neuf stades ontogéniques**, datés par la dentition et l'édification radiculaire :

| Code | Stade | Âge | Repère sur notre axe (jours post-fécondation, approx.) | Denture |
|---|---|---|---|---|
| A | 5 mois fœtal | ~20 sem. | ≈ J140 | — |
| B | 7 mois fœtal | ~28 sem. | ≈ J196 | — |
| C | naissance | 0 j | ≈ J266 | — |
| D | 8 mois ½ | 0,7 an | ≈ J525 | denture temporaire en constitution |
| E | 2 ans | 2 a | ≈ J996 | denture temporaire stable |
| F | 4 ans | 4 a | ≈ J1 726 | denture temporaire stable |
| G | 8 ans ½ | 8,5 a | ≈ J3 370 | denture mixte |
| H | 14 ans | 14 a | ≈ J5 380 | denture adolescente / M2 |
| Ad | adulte | — | — | denture adulte ; **crânes de ~65 ans** (remodelage sénile de l'orbite, manuscrit) |

Les jours sont calculés avec « mois fœtal » = mois depuis la fécondation (comme les mois lunaires de la figure 1 de
Scammon reprise par Captier) et une naissance à 38 semaines post-fécondation (J266) ; ils servent à placer les stades
sur l'axe J0-J91 du viewer, pas à dater finement.

### 2.2 Points de repère

142 points par crâne : 32 impairs (médians) et 55 paires ; en traitant les deux hémicrânes comme identiques, le
manuscrit travaille sur **87 points**. Points anatomiques ou construits (« vestibiens » : définis par le repère,
ex. frontal antérieur vestibien, vertex vestibien, pré-vestibion). Tableau 1 du manuscrit, regroupé par unité :

| Unité | Point | Abr. | Position |
|---|---|---|---|
| Calvaria | Glabelle | Gl | médian |
| Calvaria | Pré-vestibion | PV | médian |
| Calvaria | Frontal antérieur vestibien | Ant | médian (construit) |
| Calvaria | Bregma exocrânien | Br | médian |
| Calvaria | Vertex vestibien | Ve | médian (construit) |
| Calvaria | Lambda exocrânien | La | médian |
| Calvaria | Opisthocrânion | Opc | médian |
| Calvaria | Euryon | Eu | bilatéral |
| Calvaria | (non identifié dans le brouillon) | Cr | bilatéral |
| Base du crâne | Nasion | Na | médian |
| Base du crâne | Centre de la selle turcique | S (noté « SI » dans le brouillon, homonyme de la symphyse) | médian |
| Base du crâne | Basion | Ba | médian |
| Base du crâne | Inion | In | médian |
| Base du crâne | Opisthion | Opt | médian |
| Base du crâne | Porion | Po | bilatéral (origine du repère de Francfort) |
| Base du crâne | Processus mastoïde | Mas | bilatéral |
| Base du crâne | Foramen lacerum (trou déchiré antérieur) | tda | bilatéral |
| Base du crâne | Foramen jugulaire (trou déchiré postérieur) | tdp | bilatéral |
| Base du crâne | Méat acoustique interne | VIII | bilatéral |
| Massif facial | Épine nasale antérieure | ENA | médian |
| Massif facial | Épine nasale postérieure | — | médian |
| Massif facial | Rhinion | Rh | médian |
| Massif facial | Maxillo-zygomatique | MMI | bilatéral |
| Massif facial | Naso-maxillaire | NM | bilatéral |
| Massif facial | Narinaire | Nae | bilatéral |
| Massif facial | Largeur maxillaire maximale | Mac | bilatéral |
| Massif facial | Largeur maxillaire canine-incisive | Mi | bilatéral |
| Massif facial | Prosthion | Pr | médian |
| Orbite | Orbitaire externe | Oex | bilatéral |
| Orbite | Orbitaire inférieur | Oif | bilatéral (plan de Francfort) |
| Orbite | Orbitaire médian | Om | bilatéral |
| Orbite | Orbitaire interne (dacryon) | Oit | bilatéral |
| Orbite | Orbitaire supérieur | Os | bilatéral |
| Mandibule | (non identifié dans le brouillon) | De | médian |
| Mandibule | Infradentale | Id | médian |
| Mandibule | Pogonion | Pg | médian |
| Mandibule | Menton | Me | médian |
| Mandibule | Symphyse | SI | médian |
| Mandibule | Condyle externe | Ce / Co | bilatéral |
| Mandibule | Gonion | Go | bilatéral |

Les figures 1 et 2 du manuscrit (profil et face) montrent en outre « ta » (près de la mastoïde, non défini) et les
trajets de tous les points sur les neuf stades.

Le chapitre numérote différemment 20 points principaux (figure 13) : 1 condylion, 2 gonion, 3 symphyse, 9 épine
nasale antérieure, 12 nasion, 13 glabelle, 14 bregma, 17 lambda, 19 inion, 20 opisthion ; **points 1-12 = unité
faciale F, points 12-20 = unité cérébrale EC**.

### 2.3 Ce que les documents ne contiennent pas

Les **coordonnées (x, y, z) des 142 points aux neuf stades** ne figurent dans aucun des trois documents : elles sont
« publiées et disponibles » (Fenart et Biecq-Sellier 2004, *International Orthodontics* 2(4):265-278 ; Fenart 2006,
*Ontogenèse craniométrique vestibulaire*, chez l'auteur / Reproflash, Lille) et dans les fichiers de Captier et Boë
(figures Matlab). Sans elles, le modèle ci-dessous reste qualitatif + quelques distances. Voir § 7.4.

---

## 3. Le modèle de croissance à deux unités (Captier et Boë)

### 3.1 Les deux unités de part et d'autre de la base du crâne

| Unité | Contenu | Moteur | Profil |
|---|---|---|---|
| **EC — crâne cérébral** (calvaria + encéphale) | voûte, région occipitale ; points 12-20 | croissance encéphalique (deux sous-unités : compartiment supra-tentoriel et infra-tentoriel) | **rapide et précoce** : proche de la taille adulte très tôt |
| **F — crâne facial** (massif facial + mandibule) | de la glabelle à la symphyse ; points 1-12 | alimentation, ventilation, parole ; maintien de l'occlusion (modèle cybernétique de Petrovic) | **lent et prolongé**, jusqu'à la fin de l'adolescence |

La **base du crâne** est le socle commun : croissance lente, à dominante génétique (synchondroses, sensibles à la GH).

### 3.2 De profil (figures 13-15 du chapitre, figure 1 du manuscrit)

1. **Expansion radiale depuis le porion** pour les deux unités : les trajets des points rayonnent autour de Po
   (Po est par construction immobile : déplacement 0 dans le tableau 2).
2. **Unité F : expansion radiale quasi exclusive**, « projection de la face en bas et en avant ». Trajets
   globalement **rectilignes**. La région est divisée par la cavité nasale : partie supérieure fronto-nasale (Ant →
   rhinion) dirigée en avant et vers le haut, partie inférieure maxillo-mandibulaire dirigée en avant et vers le bas
   (directions divergentes autour de l'orifice piriforme).
3. **Unité EC : expansion radiale + rotation**. Trajets **curvilignes**. Rotation vers l'arrière pour les points du
   sommet (15, 16, 17 : vertex, lambda), vers le bas pour la région occipitale (18, 19, 20 : opisthocrânion, inion,
   opisthion). Pas de composante de rotation au front. Sur les profils normalisés en taille (figure 14) : de la
   naissance à 4 ans, **horizontalisation de l'inion**, stable jusque vers 14 ans, puis nouvelle rotation vers le bas
   jusqu'à l'adulte (poussée pubertaire de l'écaille et de la base de l'occipital).
4. **Gonion** : légèrement en bas et en avant jusqu'à la naissance, puis trajet vertical rectiligne jusqu'à l'adulte.
5. **Proportions** : la hauteur de la face sous le plan de Francfort passe de **25 % à la naissance à 40 % chez
   l'adulte** ; **à 4 ans, la moitié du trajet** jusqu'à l'adulte est parcourue ; l'angle mandibulaire devient
   prédominant à partir de l'adolescence.

### 3.3 De face (figures 16-17 du chapitre, figure 2 du manuscrit)

- Expansion radiale des deux unités, avec rotation vers le bas de la partie latérale d'EC (euryon), nette après 2 ans.
- Entre 4 et 14 ans, EC ne croît que **latéralement** ; le sommet de la calvaria ne bouge pas.
- **L'écart inter-orbitaire ne change pratiquement pas** de la naissance à l'adulte ; la région inter-orbitaire et
  l'orifice piriforme sont très peu modifiés transversalement.
- Orbite : croissance vers le haut et latéralement, rapide jusqu'à 4 ans puis lente ; hauteur quasi fixe après 2 ans ;
  net élargissement latéral après 14 ans (l'adulte à 65 ans ajoute un remodelage latéro-inférieur).
- Déplacements verticaux du viscérocrâne > neurocrâne ; élargissement de la mandibule (gonion) très net après 14 ans ;
  la mastoïde descend et s'écarte jusqu'à 14 ans puis se médialise.
- Résultat global : passage d'une tête fœtale **trapézoïde à base supérieure** à un ovoïde allongé adulte ; le
  viscérocrâne s'élargit plus que le neurocrâne, d'où l'impression de rétrécissement de celui-ci.

### 3.4 Analyse en composantes principales (manuscrit, chiffres « à vérifier » par les auteurs)

- Sur les 87 points, données brutes non normées : **deux composantes expliquent > 98,5 % de la variance**, en 3D
  comme dans chaque plan.
- Profil : **PC1 ≈ 84 % = expansion radiale ; PC2 ≈ 12 % = rotation postérieure (dorsale et caudale) du
  neurocrâne**, effet modeste sur la face sauf autour de la cavité nasale.
- Face : PC1 = composante verticale (plus marquée au maxillo-facial qu'au neurocrâne) ; PC2 = composante transversale
  (élargissement homogène, inter-orbitaire et cavité nasale inchangés).

Autrement dit, pour un modèle numérique, **une homothétie centrée sur le porion plus une rotation du seul
neurocrâne postérieur** rendent compte de l'essentiel de la forme à chaque âge.

---

## 4. Vitesses et hétérochronie

### 4.1 Trois périodes (figure 18 du chapitre, figure 6 du manuscrit)

| Période | Stades | Comportement |
|---|---|---|
| Synchrone fœtale | A → B (5 → 7 mois fœtal) | toutes les vitesses **ralentissent** jusqu'à 7 mois, puis remontent vers la naissance |
| Asynchrone | B → G (7 mois fœtal → 8,5 ans) | accélération jusque vers 7-8 mois postnatals (vertex et lambda les plus rapides), puis décroissance jusqu'à 8,5 ans ; petit regain à 4 ans (symphyse, nasion) ; gonion presque constant |
| Synchrone pubertaire | G → Ad | forte accélération après 14 ans (symphyse très au-dessus des autres, puis inion, nasion, lambda, vertex, ant), puis chute jusqu'à la taille adulte |

Le manuscrit écrit « à part l'inion » pour le pic pubertaire, mais la courbe de l'inion (figure 6) montre bien un pic
après 14 ans : incohérence du brouillon, à trancher avec les auteurs.

Repères généraux du chapitre : croissance cérébrale très rapide chez le nourrisson, qui **ralentit fortement à 2 ans ½**
et atteint **95 % à 7 ans** (profil neural de Scammon, sans pic pubertaire) ; la stature suit le profil somatique
(pic pubertaire) ; calvaria rapide (encéphale), base lente (génétique/GH), face lente et prolongée (occlusion).

### 4.2 Distances parcourues de 5 mois fœtal à l'adulte et âge où 90 % du trajet est fait (tableau 2 du manuscrit)

| Point | Déplacement total (mm) | Âge à 90 % (ans) |
|---|---|---|
| Glabelle (Gl) | 104 | 15,2 |
| Bregma (Br) | 108 | **1,75** |
| Vertex (Ve) | 109 | **1,75** |
| Lambda (La) | 97 | **1,9** |
| Opisthion (Opt) | 96 | 7,75 |
| Nasion (Na) | 98 | 11,75 |
| Inion (In) | 80 | 14,9 |
| Opisthocrânion (Opc) | 45 | 15,6 |
| Porion (Po) | 0 | — (origine) |
| Mastoïde (Mas) | 37 | 20,6 |
| Prosthion (Pr) | 111 | 17,6 |
| Infradentale (Id) | 118 | 18,25 |
| Pogonion (Pg) | 132 | 18,75 |
| Symphyse (SI) | 128 | 19,25 |
| Condyle (Co) | 22 | 18,1 |
| Orbitaire inférieur (Oif) | 86 | 15,4 |

Lecture : la **voûte** (Br, Ve, La) a fait 90 % de son chemin **avant 2 ans** ; la **mandibule** (Pg, SI, Id) a les plus
grands déplacements et finit vers 19 ans ; le condyle a le plus petit déplacement transversal. Les autres colonnes du
tableau (vitesses max 1 et 2) sont vides dans le brouillon.

### 4.3 Proportions corporelles (rappels du chapitre, utiles pour l'axe du viewer)

- Tête + cou = moitié de la longueur du corps chez l'embryon de 3 semaines ; hauteur de la tête = 1/4 de la taille
  à la naissance et jusqu'à 1 an, 1/3 (sic, tel quel dans le chapitre) à 3 ans, 1/6 à 5 ans, 1/7 à 10 ans, 1/8 adulte.
- Période embryonnaire (sem. 5-8) = croissance **allométrique** (chaque organe à sa vitesse, forme extérieure très
  changeante) ; période fœtale à partir du 5e mois lunaire = croissance globalement **isométrique**, tête plus ronde,
  front réduit par rapport à la face.

---

## 5. Organisation modulaire du crâne (lecture embryologique)

| Module | Origine / contrôle | Ossification | Dans embryo3d |
|---|---|---|---|
| **Neurocrâne préchordal** (en avant de la fosse hypophysaire) | **Hox négatif**, crêtes neurales ; capsules **optique** et **nasale** | enchondrale (chondrocrâne antérieur) puis sutures | `yeux`, `cristallins`, partie antérieure de `chondrocrane` |
| **Neurocrâne parachordal** (en arrière) | **Hox positif**, comme la colonne ; capsule **otique** dans le rocher | enchondrale ; synchondrose sphéno-occipitale active jusqu'à la puberté (direction bas-arrière) | `vesicules_otiques`, partie postérieure de `chondrocrane`, continuité avec `squelette_axial_cartilage` / `notochorde` |
| **Splanchnocrâne** | arcs pharyngiens ; 1er arc (mandibule) Hox−, arcs suivants Hox+ (hyoïde, styloïde, stapes, larynx) | enchondrale (cartilage de Meckel régresse) | non segmenté à ce jour |
| **Dermatocrâne** | mésenchyme péri-encéphalique et des bourgeons faciaux | membraneuse (sutures = « joints de dilatation ») | `enveloppe` seulement ; pas de label osseux |
| Encéphale | moteur d'EC ; compartiments supra/infra-tentoriel déjà différenciés chez le fœtus | — | `snc`, `ventricules`, segmentation fine `prosencephale`/`mesencephale`/`rhombencephale` |

Points d'embryologie directement utiles à nos stades :

- Le **chondrocrâne est mature et continu à la fin de la période embryonnaire (CS23)** (Müller et O'Rahilly 1980) ;
  sa croissance prénatale est multidirectionnelle à partir d'un **centre géométrique à la synchondrose
  intra-sphénoïdale**, et c'est lui qui fixe largeur et longueur de la base, donc la position des os de membrane.
- La limite Hox+/Hox− passe au **tuberculum sellae** ; le champ facial (non Hox) et le champ vertébral (Hox) sont
  deux moitiés du crâne de part et d'autre de la selle (figure 12 du chapitre).
- Deux grands types de mouvement : **ventral linéaire** (champ cranio-facial antérieur, crêtes neurales) et
  **postérieur avec expansion-rotation horaire** (neurocrâne postérieur), déjà présents chez le fœtus.

---

## 6. L'atlas 3D de Pellerin (2026) : complément moderne de l'atlas de Fenart

- Motivation : l'atlas de Fenart est en 2D (projections orthogonales), sur crânes secs (densité osseuse altérée),
  huit stades seulement ; comparer un patient à un seul crâne normal est un biais.
- Matériel : 184 scanners d'enfants « normaux » (Lille neurochirurgie pédiatrique, Lille imagerie, New Mexico Decedent
  Image Database), 181 exploitables, triés en **19 groupes d'âge** → 19 gabarits « Skull Atlas ###-###d » (âges en
  jours du plus jeune et du plus âgé du groupe ; exemple cité : 645-765 j).
- Méthode (3D Slicer / SlicerMorph) : segmentation os à 250 UH, modèle décimé à 0,8, points de repère aux
  intersections plans vestibulaires × surface externe, **DeCA** (Dense Correspondence Analysis, Rolfe et Maga 2023)
  pour le modèle moyen, lissage (Surface Toolbox), vérification statistique par **GPA** (ellipses équiprobables).
- Livrable : fichiers **.mrb** (bundle 3D Slicer) contenant le modèle moyen + ses repères + les trois plans
  vestibulaires. Utilisation : **FastModelAlign** (Procuste, forme seule) ou **superposition vestibulaire** par une
  transformée appliquée au dossier atlas pour faire coïncider les plans (forme + position), atlas en fil de fer.
- Disponibilité : sur demande à l'auteur (philippepellerin@me.com).

Intérêt pour nous : notre pipeline exporte déjà densité + labels en **NRRD pour 3D Slicer** ; les gabarits .mrb
seraient directement superposables à un volume fœtal ou postnatal placé dans le même repère vestibulaire.

---

## 7. Application à embryo3d

### 7.1 Un repère commun de CS13 à l'adulte : l'axe de Perez porté par les otocystes

Les plans de Francfort (porion, orbitaire) n'existent pas chez l'embryon : pas d'os. En revanche l'**otocyste**
(label `vesicules_otiques`, fermé dès CS13) devient le labyrinthe membraneux puis la capsule otique du chondrocrâne ;
le canal semi-circulaire latéral en est issu. La droite joignant les centres des deux vésicules otiques est donc, à
une petite translation près, **l'axe de Perez embryonnaire**, et le plan sagittal médian est déjà défini par la
notochorde / le tube neural. Il manque la troisième direction (le plan horizontal vestibulaire) : chez l'embryon
recourbé on peut l'approcher par le plan contenant l'axe de Perez et un point médian stable de la base
(hypophyse / selle, à la limite pré/parachordale, c'est-à-dire près de l'extrémité rostrale de la `notochorde`).

Proposition : calculer dans `topographie_axe.py` ou un script dédié, pour chaque stade, le repère
`{origine = milieu des otocystes, X = axe de Perez, plan sagittal = notochorde + ventricules}` et l'écrire dans
`out/topographie/repere_vestibulaire.json`. Il servirait à (a) superposer les stades CS13-CS20 autrement que par la
boîte englobante, (b) raccorder la tête embryonnaire aux stades fœtaux de Fenart (A = 5 mois) et à l'atlas de Pellerin,
(c) exprimer la tête dans le repère qu'utilisent cliniciens et craniologues.

### 7.2 Prolonger l'axe du temps du viewer

Aujourd'hui l'axe du viewer et de l'agrégateur va de J15 (CS7) à J91 (VOKA 13 semaines, référence externe). La base de
Fenart apporte neuf jalons de **J140 à l'adulte** pour la tête seule (voir tableau § 2.1). Entre J91 et J140 il n'y a
rien dans nos sources : la biométrie fœtale (Guihard-Costa et Larroche 1995) et la céphalométrie fœtale de Captier
(2009, 2011 ; Herlin et al. 2011 pour la base) sont les pistes citées par le chapitre.

### 7.3 Règles à encoder pour une tête « qui grandit » après CS23

Avec les seules données disponibles (qualitatif + tableau 2), un premier modèle paramétrique de la tête peut être :

1. Repère vestibulaire ou de Francfort, origine Po (ou milieu de l'axe de Perez).
2. Unité F (points en avant de la verticale de Po, sous Na) : **homothétie centrée Po**, facteur s(t) tiré d'une
   double logistique (pics vers 8 mois postnatals et vers 15-16 ans), déplacements totaux du tableau 2 pour caler
   l'amplitude (Pg 132 mm, SI 128 mm, Pr 111 mm, Oif 86 mm) ; directions rectilignes.
3. Unité EC postérieure (de Br à Opt) : homothétie **plus rotation** d'angle θ(t) autour de Po, θ positif vers
   l'arrière-bas, le gros de la rotation entre la naissance et 4 ans puis après 14 ans ; amplitudes du tableau 2
   (Br 108, Ve 109, La 97, Opt 96, In 80, Opc 45 mm) avec 90 % du chemin atteint à 1,75-1,9 ans pour la voûte.
4. Front (Ant, Gl, Na) : homothétie seule (pas de rotation), mais tardive pour Gl (sinus frontal, 90 % à 15 ans).
5. De face : largeur inter-orbitaire constante ; élargissement homogène ; euryon rotation vers le bas après 2 ans ;
   mandibule élargie après 14 ans.
6. Contrainte globale : hauteur faciale sous Francfort 25 % → 40 % de la hauteur totale ; mi-chemin à 4 ans.

Ces règles sont reprises sous forme structurée dans `docs/modele_croissance_fenart.json` (stades, points, tableau 2,
phases de vitesse, règles par unité) pour qu'un script puisse les lire.

### 7.4 Ce qu'il faut obtenir pour passer au quantitatif

- Les **coordonnées des 142 points × 9 stades** : article Fenart et Biecq-Sellier 2004 (*International
  Orthodontics*), ouvrage Fenart 2006, ou fichiers de Captier/Boë (les figures 1-2 et 7-8 du manuscrit sont tracées
  à partir de ces fichiers).
- Les **gabarits .mrb de Pellerin** (19 âges, 0 à ~18 ans) pour une surface 3D moyenne du crâne par âge, dans le même
  repère.
- Pour la période J60-J140 : céphalométrie fœtale de Captier (2009, 2011) et modèles géométriques encéphale /
  base du crâne de Captier et al. 2013 (*Morphologie* 97:38-47).

### 7.5 Limites et biais à garder en tête

- Crânes **secs**, étude **transversale** (moyennes par stade, pas de trajets individuels), **sans sexe** : acceptable
  jusqu'à ~14 ans (pas de différence de forme 0-4 ans ; croissance des mâchoires similaire avant 14 ans), le
  dimorphisme apparaît ensuite (filles plus précoces, moins intenses, plus horizontales).
- **Adultes = crânes de ~65 ans** : remodelage sénile de l'orbite et perte dentaire possible ; l'extrémité des trajets
  n'est pas celle d'un adulte jeune.
- Os seulement : aucune information sur l'encéphale, les parties molles ni les cavités ; c'est le squelette qui
  « enregistre » la croissance des tissus mous (matrices fonctionnelles de Moss, conformateurs de Couly).
- La grande variabilité individuelle de l'adolescence (début, durée, intensité du pic) rend toute prédiction
  individuelle illusoire ; le modèle vaut pour une tête moyenne.

---

## Bibliographie (telle que citée par les trois sources)

- Broca P. Nouvelles recherches sur le plan horizontal de la tête et sur le degré d'inclinaison des divers plans
  crâniens. BMSAP 1873;8:542-566.
- Perez F. Crâniologie vestibienne, ethnique et zoologique. BMSAP 1922;3:16-32.
- Delattre A., Fenart R. Méthode vestibulaire et craniométrie ; détermination des axes vestibulaires et coordonnées
  vestibulaires. BMSAP 1953;4(5):543-549.
- Delattre A., Fenart R. *L'hominisation du crâne, étudiée par la méthode vestibulaire.* Paris, CNRS, 1960, 418 p.
- Cousin R.-P. *Étude en projection sagittale de crânes d'enfants orientés dans les axes vestibulaires.* 1969.
- Fenart R., Biecq-Sellier M. Trajets ontogéniques des points du squelette céphalique : distances parcourues et
  vitesses de déplacement dans le référentiel vestibulaire. International Orthodontics 2004;2(4):265-278.
- Fenart R. *Ontogenèse craniométrique vestibulaire. Analyse morphométrique positionnelle.* Lille, Reproflash /
  chez l'auteur, 2006-2007.
- Fenart R., Landouzy J. Valuation of facial size by vestibular coordinates, ontogenesis-phylogenesis. 2009.
- Pellerin P. Creating an Atlas of Growing Skull Templates in Vestibular Orientation for Analysis of Craniofacial
  Malformations. J Craniofac Surg 2026;37(1/2):16-19. doi:10.1097/SCS.0000000000012132
- Piral T., Pertuzon B., Pellerin P. et al. Rev Stomatol Chir Maxillofac 1994;95:151-153 et 1995;96:207-209
  (première application clinique et principes géométriques de l'orientation vestibulaire au scanner).
- Vinchon M., Pellerin P., Pertuzon B. et al. Vestibular orientation for craniofacial surgery. Childs Nerv Syst
  2007;23:1403-1409.
- Rolfe S.M., Maga A.M. DeCA: a dense correspondence analysis toolkit for shape analysis. Shape in Medical Imaging,
  Springer 2023:259-270.
- Captier G., Boë L.-J. Croissance générale et céphalique. In *Chirurgie plastique de l'enfant et de l'adolescent*,
  p. 75-93.
- Captier G., Boë L.-J., Badin P., Guihard-Costa A.-M., Canovas F., Larroche J.-C. Modèles géométriques de croissance
  du cerveau, cervelet, tronc cérébral et modification des angles de la base du crâne au cours de la période fœtale.
  Morphologie 2013;97(317):38-47.
- Captier G. *Analyse céphalométrique de la croissance faciale fœtale.* Éditions Universitaires Européennes, 2011.
- Herlin C., Largey A., de Matteï C., Daurès J.-P., Bigorre M., Captier G. Modeling of the human fetal skull base
  growth. Early Hum Dev 2011;87(4):239-245.
- Müller F., O'Rahilly R. The human chondrocranium at the end of the embryonic period, proper. Am J Anat
  1980;159:33-58.
- Guihard-Costa A.-M., Ramirez-Rossi F. Growth of the human brain and skull slows down at about 2.5 years. 2004.
- Scammon R.E. The measurement of the body in children. 1930 ; Scammon R., Calkins L. 1929 (dimensions externes du
  fœtus).
- Enlow D.H., Hans M.G. *Essentials of facial growth.* Saunders, 1996.
- Moss M.L., Salentijn L. The primary role of functional matrices in facial growth. Am J Orthod 1969;55:566-577.
- Couly G. La dynamique de croissance céphalique. Le principe de conformation organo-fonctionnelle. Act Odontol
  1977;118:63-96.
- Petrovic A. Control of postnatal growth of secondary cartilages of the mandible by mechanisms regulating occlusion.
  Trans Eur Orthod Soc 1974:69-75.
- Liang C., Profico A., Buzi C., Khonsari R.H. et al. Normal human craniofacial growth and development from 0 to 4
  years. Sci Rep 2023;13:9641.
- Ursi W.J., Trotman C.A., McNamara J.A., Behrents R.G. Sexual dimorphism in normal craniofacial growth. Angle
  Orthod 1993;63:47-56.
