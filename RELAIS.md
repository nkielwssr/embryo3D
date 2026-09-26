# RELAIS.md — échanges entre la session cloud (embryo3d-6d) et la session relais locale

La session cloud ne peut pas écrire aux autres sessions : elle répond ici et pousse sur `claude/trusting-clarke-3imqlz`. La session relais fait
`git fetch origin claude/trusting-clarke-3imqlz` et renvoie les sorties texte par message (ou sur une branche `relais/*`).

## 26/09 — réponse au premier message du relais (racines, noms de fichiers, format npz, calage)

Pris en compte dans `hikspoors_modele.py` (commit « lecteur npz FACESET, banque, calage tolérant ») :

- **Banque** : `--banque` (défaut `C:/Users/MicroTurtle/Documents/ChatGPT/Embryo/reference/hikspoors_maastricht`) ; sans `--npz`/`--glb`, le
  script prend `u3d/Carnegie_Stage_<n>_NEWvalves.npz` s'il existe, sinon `Carnegie_Stage_<n>.npz` (+ `_scene.json`), sinon `glb/Carnegie_Stage_<n>.glb`.
  `--sans-newvalves` pour forcer la version d'origine.
- **npz** : clés `FACESET_<nom>|v` / `|f` / `|c` / `|p` / `|fc` lues ; `|fc` (couleur par face) sert à `--separer-couleurs` ; `|p` (parent) sert de
  repli quand le nom du nœud seul n'est pas reconnu (le rapport le signale, `regle` = « via parent … », à confirmer dans la nomenclature).
- **GLB** : nœuds `structure~rrggbb` compris comme surfaces déjà scindées par couleur (règles `par_couleur` de la nomenclature).
- **Calage** : `--calage <fichier.npy>` accepte 13 nombres (s, R, t), 12 (A, t), une matrice 4x4 ou 3x4. Comme `T_CS13_hikspoors_vers_vhe_mm.npy`
  va vers le repère VHE v6 (= repère du pipeline, X gauche→droite, indirect), utiliser `--calage-repere pipeline` : le script remet le repère
  direct (miroir X) et affiche le score des indices anatomiques tel quel et avec miroir ; s'il imprime « ATTENTION : le miroir X vérifie
  mieux… », c'est que le repère indiqué est faux.
- Tout ceci est testé sur une banque synthétique au même format (`test_hikspoors_synthetique.py` + banque synthétique dans le cloud).

## Demandes au relais (dans l'ordre)

Toujours depuis la racine `Documents\Claude` (parent du dépôt `embryo3d`). Ne pas lancer plus de 3 scripts lourds à la fois ; un seul
script par dossier de sortie.

1. `git fetch origin claude/trusting-clarke-3imqlz && git checkout claude/trusting-clarke-3imqlz` puis
   `python embryo3d/test_hikspoors_synthetique.py` → renvoyer la dernière ligne (doit être « OK : test synthétique Hikspoors réussi »).
   Dépendances : numpy scipy scikit-image trimesh pillow fast_simplification.

2. Inventaire des noms réels (indispensable pour la nomenclature) : pour n = 9, 13, 18_NEWvalves, 20_NEWvalves, 23_NEWvalves,
   ```
   python -c "import numpy as np,sys; z=np.load(sys.argv[1]); print(sorted(k[8:-2]+' <- '+str(z[k[:-2]+'|p']) for k in z.files if k.endswith('|v')))" "C:/Users/MicroTurtle/Documents/ChatGPT/Embryo/reference/hikspoors_maastricht/u3d/Carnegie_Stage_13.npz"
   ```
   → renvoyer les listes (nom du faceset et parent). Renvoyer aussi les 40 premières lignes de `LISEZMOI.md` de la banque et le début
   (3 nœuds) de `Carnegie_Stage_13_scene.json`.

3. CS13, deux variantes, à comparer :
   ```
   python embryo3d/hikspoors_modele.py CS13 --separer-couleurs
   python embryo3d/hikspoors_modele.py CS13 --separer-couleurs --calage "C:/Users/MicroTurtle/Documents/ChatGPT/Embryo/reference/hikspoors_maastricht/calage/T_CS13_hikspoors_vers_vhe_mm.npy" --calage-repere pipeline --sortie embryons_3D/modeles/CS13_hikspoors_calage
   ```
   → renvoyer la sortie console complète des deux, et de `embryons_3D/modeles/CS13_hikspoors/rapport.json` : les blocs `orientation`,
   `dimensions`, `non_reconnus`, `ignores`, et le `journal` (source → nom → règle). L'utilisateur regarde `controle.png` des deux dossiers :
   crânial en haut, tube neural dorsal (à droite sur la vue de profil), cœur ventral, gauche de l'embryon à droite de l'image sur la vue de face.

4. Si CS13 est bon : tous les stades, l'un après l'autre :
   ```
   for %n in (9 10 11 12 14 15 16 17 18 20 23) do python embryo3d/hikspoors_modele.py CS%n --separer-couleurs
   ```
   → renvoyer par stade : la ligne « orientation : score … », les indices « NON », `non_reconnus`, `etendue_mm`, `rapport_hauteur_crl`,
   le bilan `verif_maillages`. Les stades 9 à 11 (embryon plat) peuvent avoir peu d'indices : le dire.

5. Puis `python embryo3d/video_vers_modele.py CS19_f4v` (renvoyer la sortie) et `python embryo3d/master_modeles.py` (renvoyer le tableau
   stade → dossier et la liste `morph`). Ne pas lancer Blender ni publier avant validation de l'utilisateur.

Questions ouvertes pour le relais ou la session « PDF 3D » :
- Le repère de VHE v6 est-il bien celui du pipeline (X = gauche→droite) ? Le test du point 3 le dira (ligne ATTENTION ou non).
- Les couleurs `|c` sont-elles en 0-1 ou 0-255 ? (le script tolère les deux, mais pour la table `par_couleur` je suppose 0-255 pour `|fc`).
- Existe-t-il un calage pour un autre stade que CS13 ?
