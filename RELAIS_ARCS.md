# RELAIS_ARCS.md — arcs aortiques, sac aortique et poches pharyngiennes : échanges entre la session cloud et la session relais locale

Suite de la tâche Fable « Morphing aorte et tube digestif » (branche `claude/sleepy-noether-goh2of`, reprise sur `claude/eloquent-edison-s471eq`).
La chaîne `tubes_morph.py` accepte les arcs, le sac, le tronc artériel et les poches depuis le 26/09. Ce qui manque, c'est le relevé des points sur
les volumes CS13 et CS14, qui ne sont pas dans le dépôt. La session cloud ne peut pas écrire aux autres sessions : elle répond ici et pousse sur
`claude/eloquent-edison-s471eq`. Le relais renvoie les sorties texte par message, ou sur une branche `relais/arcs` (même principe que `relais/retours`).

## 26/09 — état de la branche

- `test_tubes_morph_synthetique.py` : test sans données de `tubes_morph.py` et de `tubes_morph_planche.py` (présences, raccords des absents,
  calendrier, mésos, planche). Le code de Fable passe sans modification.
- `arcs_candidats.py` : propose les arcs et le sac à partir des volumes réels (détail dans le README, section « Arcs aortiques… », point
  « Relevé semi-automatique »), plus un mode `fusionner` qui recopie les segments validés dans `vaisseaux_points/<CS>.json`.
  `test_arcs_candidats_synthetique.py` le vérifie sur un volume synthétique avec leurres : 7 arcs sur 7 retrouvés à ±2 voxels, bien numérotés,
  un kyste marqué « douteux ».
- `vaisseaux_points/calendrier_arcs.json` : aligné sur le résumé de Rana 2014 fourni par l'utilisateur (arcs 1 et 2 déjà en régression à CS13,
  disparus ou réduits à des reliquats à CS14) et sur Graham 2023 (6e arc formé entre CS14 et CS15 selon l'embryon).
- `tubes_morph_planche.py` : la liste des absents passe à la ligne (elle était tronquée) et les vues restent alignées d'une vignette à l'autre.

## 26/09 (suite) — pharynx et poches d'abord (décision de l'utilisateur)

- `pharynx_candidats.py` : propose le pharynx (du fond de la cavité buccale au début de l'œsophage, par la lumière pâle au-dessus de l'œsophage
  tracé, centré sur le plan médian) et les poches pharyngiennes (maxima de l'extension latérale de cette lumière), plus un mode `fusionner` vers
  `digestif_points/<CS>.json`. `test_pharynx_candidats_synthetique.py` : pharynx aplati et courbé, 4 poches à gauche et 3 à droite retrouvées à
  ±3 voxels et bien numérotées, leurres écartés (bourgeon pulmonaire, ventricule cérébral).
- Ordre : pharynx et poches, puis arcs. Une fois le pharynx construit (`digestif.npz`), `arcs_candidats.py` l'exclut des chemins des arcs.

## Demandes au relais (dans l'ordre)

Toujours depuis la racine `Documents\Claude` (parent du dépôt `embryo3d`). Pour ne pas déranger le dépôt principal, travailler dans un
worktree voisin : les scripts trouvent les dossiers de stades (`CS13.f4v`, `CS14_f4v`…) dans le dossier parent, comme depuis `embryo3d`.

1. Préparer le worktree et lancer les tests :
   ```
   git -C embryo3d fetch origin claude/eloquent-edison-s471eq
   git -C embryo3d worktree add ../embryo3d_arcs origin/claude/eloquent-edison-s471eq
   python embryo3d_arcs/test_tubes_morph_synthetique.py
   python embryo3d_arcs/test_pharynx_candidats_synthetique.py
   python embryo3d_arcs/test_arcs_candidats_synthetique.py
   ```
   (Si le worktree existe déjà : `git -C embryo3d fetch origin claude/eloquent-edison-s471eq` puis
   `git -C embryo3d_arcs checkout --detach origin/claude/eloquent-edison-s471eq`.)
   → renvoyer la dernière ligne de chaque test (« OK : test synthétique … réussi »). Dépendances : numpy scipy opencv-python scikit-image.

2. Vérifier les entrées de CS13 et CS14 :
   ```
   python -c "import numpy as np,json,os,sys; w=sys.argv[1]; print('labels', np.load(w+'/labels.npz').files); p=w+'/cardio/coeur.npz'; print('coeur', np.load(p).files if os.path.exists(p) else 'absent'); print('vaisseaux', list(json.load(open(w+'/cardio/vaisseaux_chemins.json')))); p=w+'/digestif/chemins.json'; print('digestif', list(json.load(open(p))) if os.path.exists(p) else 'absent'); print('dens', np.load(w+'/dens.npy', mmap_mode='r').shape)" CS13.f4v/work
   ```
   (idem avec `CS14_f4v/work`) → renvoyer les cinq lignes. Il faut au moins `oesophage` dans `digestif` ; `aorte_dorsale_gauche` et
   `aorte_dorsale_droite` dans `vaisseaux` ; `enveloppe` et `cavite_pericardique` dans `labels` (et si possible `ventricules`, `yeux`,
   `vesicules_otiques`, exclus pour le pharynx) ; un masque cardiaque (`coeur.npz` ou `coeur_detoure` / `coeur` dans `labels`) pour la graine du sac.

3. Proposer le pharynx et les poches, CS13 puis CS14 (lecture seule : rien n'est écrit dans `labels.npz`, `manifest.json` ni
   `digestif_points/<CS>.json`) :
   ```
   python embryo3d_arcs/pharynx_candidats.py CS13.f4v/work CS13
   python embryo3d_arcs/pharynx_candidats.py CS14_f4v/work CS14
   ```
   Si c'est trop lent ou si la mémoire manque (la boîte couvre toute la région pharyngienne), ajouter `--pas 2`. → renvoyer la sortie console
   complète et pousser sur `relais/arcs` : `CS13.f4v/work/digestif/pharynx_candidats.png`, `CS14_f4v/work/digestif/pharynx_candidats.png`,
   `embryo3d_arcs/digestif_points/CS13_pharynx_proposes.json`, `embryo3d_arcs/digestif_points/CS14_pharynx_proposes.json`.
   Si la console signale que la lumière touche le bord de la boîte, le contour orange de la planche montre où elle fuit : relancer avec un
   `--seuil` plus bas (la valeur affichée moins 5 à 10).

4. Proposer les arcs, CS13 puis CS14 (lecture seule : rien n'est écrit dans `labels.npz`, `manifest.json` ni `vaisseaux_points/<CS>.json`) :
   ```
   python embryo3d_arcs/arcs_candidats.py CS13.f4v/work CS13
   python embryo3d_arcs/arcs_candidats.py CS14_f4v/work CS14
   ```
   Si c'est trop lent ou si la mémoire manque, ajouter `--pas 2`. → renvoyer la sortie console complète et pousser sur `relais/arcs` :
   `CS13.f4v/work/cardio/arcs_candidats.png`, `CS14_f4v/work/cardio/arcs_candidats.png`, `embryo3d_arcs/vaisseaux_points/CS13_arcs_proposes.json`,
   `embryo3d_arcs/vaisseaux_points/CS14_arcs_proposes.json`.
   Si la console affiche « graine déplacée » ou si la graine (point orange de la planche) n'est pas au bout crânial de la voie de sortie du cœur,
   relever sur `digestif_planches.py` un point dans la lumière du sac aortique et relancer avec `--sac x,s,y`.
   Les arcs ne sont cherchés que le long des aortes dorsales **tracées** : si leur tracé s'arrête sous la région des arcs 1-2 (tête), ces arcs
   ne peuvent pas sortir ; le dire, il faudra d'abord prolonger `aorte_dorsale_*` vers le haut dans `vaisseaux_points/<CS>.json`.

5. **Rien n'est fusionné avant validation** des planches (par l'utilisateur, ou par la session cloud sur les fichiers poussés). Après
   validation, la session cloud recopie pharynx, poches et arcs validés dans `digestif_points/<CS>.json` et `vaisseaux_points/<CS>.json` sur la
   branche (`pharynx_candidats.py fusionner`, `arcs_candidats.py fusionner`) et le signale ici. Puis, une fois la branche intégrée au dépôt
   `embryo3d` (les scripts d'aval lisent `digestif_points/` et `vaisseaux_points/` dans leur propre dossier), et avec l'accord de l'utilisateur,
   car ces étapes écrivent `digestif.npz`, `vaisseaux.npz`, `labels.npz` et les PLY du stade :
   ```
   python embryo3d/digestif_build.py CS13.f4v/work embryo3d/digestif_points/CS13.json
   python embryo3d/digestif_export.py CS13.f4v
   python embryo3d/digestif_build.py CS13.f4v/work embryo3d/vaisseaux_points/CS13.json --sortie cardio
   python embryo3d/vaisseaux_export.py CS13.f4v
   python embryo3d/fusion_systemes.py CS13.f4v --sans-scene
   ```
   (idem CS14), puis `python embryo3d/tubes_morph.py` et `python embryo3d/tubes_morph_planche.py embryons_3D/rendus/tubes_planche.png 0 0.5 1 2`
   → renvoyer la sortie de `tubes_morph.py` (tableau des présences, écarts au calendrier) et la planche.

## Questions ouvertes

- (Réglé le 26/09 : le pharynx est tracé en premier, avec les poches, par `pharynx_candidats.py`.)
- Les PDF 3D de Rana et al. (PR « Arcs aortiques : extraction des modèles 3D ») donneraient une référence indépendante à superposer aux
  propositions, une fois récupérés via le navigateur et recalés.
