# reference/ — VOKA J28→J91 comme référence externe (regard seulement)

VOKA (catalog.voka.io, section Embryologie / 1er trimestre) fournit 10 modèles **externes** de l'embryon
puis du fœtus, du jour 28 au jour 91 post-fécondation. Ils n'ont **aucune structure interne** (le panneau
« Éléments » du viewer ne liste que « Embryon, N jours »), mais l'enveloppe est très bien faite : forme,
proportions, peau translucide avec réseau vasculaire.

Ce module sert à trois choses, sans jamais extraire de maillage du viewer (conditions d'utilisation VOKA) :

1. **Confirmer nos emplacements** : comparer côte à côte une vue VOKA et notre GLB sous le même angle
   (`viewer.html`), pour vérifier que cœur, foie, yeux, membres, cordon… tombent au bon endroit sous la peau.
2. **Aspect de l'enveloppe** : les rendus VOKA servent de référence de matériau ; `blender_skin_material.py`
   reproduit cet aspect en procédural (couleur, sous-surface, vaisseaux Voronoi, brillance humide).
3. **Dimensions et proportions des stades qui nous manquent** (CS21–23 et fœtus 9–13 semaines) :
   captures multi-angles → silhouettes → mesures et « visual hull » (gabarit privé de l'enveloppe).

## Fichiers

| Fichier | Rôle |
|---|---|
| `voka_stages.json` | Les 10 stades : jour, stade Carnegie correspondant, semaine, CRL de la littérature (mm), repères anatomiques, URL page + URL catalogue (id du modèle), image, quelles coupes nous avons. |
| `voka_pages/voka_Jxx.webp` | Image publique de chaque stade (vue latérale, fond transparent) — base pour la peau. |
| `capture_voka.py` | Pilote le viewer VOKA (Playwright, Chromium) : réinitialise la vue, tourne par glissements progressifs, enregistre `captures/Jxx/Jxx_azAAA_elSEE.png`, une planche et `meta.json`. Calibre d'abord les pixels par tour (`capture_config.json`). |
| `silhouettes.py` | Silhouettes des captures, mesures (plus grande longueur en px → mm/px via le CRL), `proportions.json`, et `Jxx_hull.ply` (visual hull orthographique, repère du pipeline, en mm). |
| `viewer.html` + `serve.py` | Viewer WebGL local (three.js) : VOKA à gauche (image ou capture choisie), notre GLB/PLY à droite, bouton « caler » pour prendre le même azimut/élévation, case « miroir X » si le GLB est dans l'autre chiralité, enveloppe semi-transparente et une case par structure du GLB pour la masquer (les matériaux glTF sont forcés à metalness 0, sinon three.js les rend noirs sans carte d'environnement). Entrée « embryo-ref » de `.claude/launch.json` (port 8769). |
| `blender_skin_material.py` | Matériau `VOKA_peau` pour les objets « enveloppe » (ou une collection). |
| `probe_voka.py` | Diagnostic du viewer (iframes, boutons, rotation pas à pas) ; à relancer si VOKA change son interface. |
| `captures/Jxx/` | Vues PNG (`az` 0-315° × `el` -60/0/+60), `planche.png`, `meta.json`, `silhouettes/`, `mesures.json`, `Jxx_hull.ply` + aperçu + paramètres. |

## Lancer

```bash
python embryo3d/reference/capture_voka.py --headed       # 1re fois : se connecter dans la fenêtre qui s'ouvre
python embryo3d/reference/capture_voka.py --stages J56 J63 J70 J77 J84 J91 --az-step 30 --elev 0 60 -60
python embryo3d/reference/silhouettes.py                 # mesures + hulls
python embryo3d/reference/serve.py                       # ouvre http://127.0.0.1:8769/embryo3d/reference/viewer.html
"C:/Program Files/Blender Foundation/Blender 5.2/blender.exe" -b embryon.blend -P embryo3d/reference/blender_skin_material.py -- --collection Enveloppe --save
```

## Correspondance jours VOKA ↔ nos stades

| VOKA | Carnegie | CRL litt. (mm) | Nos coupes |
|---|---|---|---|
| J28 | CS13 | 4–6 | CS13.f4v |
| J35 | CS15 (14–15) | 7–9 | CS15_f4v |
| J42 | CS17 (16–17) | 11–14 | dossier CS18_f4v (= CS17) |
| J49 | CS19 (19–20) | 16–22 | CS19_f4v / CS20_F4V |
| J56 | CS23 (22–23) | 27–31 | — (manquant) |
| J63 | fœtus 9 sem. | 40–50 | — |
| J70 | 10 sem. | 55–61 | — |
| J77 | 11 sem. | 68–75 | — |
| J84 | 12 sem. | 80–87 | — |
| J91 | 13 sem. | 95–105 | — |

Les CRL sont des valeurs de la littérature (O'Rahilly & Müller 2010 pour l'embryon, Moore & Persaud
pour le fœtus, âges post-fécondation), pas des mesures VOKA : le viewer n'a pas de règle. Les captures
donnent donc des **proportions** ; l'échelle absolue vient de `crl_mm` dans `voka_stages.json` (à ajuster
si l'on connaît mieux). Note : le modèle « J28 » de VOKA montre déjà des palettes de main, ce qui ressemble
plutôt à CS15–16 qu'à CS13 ; se fier aux repères anatomiques plus qu'à l'étiquette de jour.

## Pièges connus (chèrement acquis)

- **Cadence d'images** : en headless, Chromium rend par défaut avec SwiftShader (~10 i/s) ; Unity lisse l'axe
  souris selon la cadence, la sensibilité devient fausse et la vue dérive. `capture_voka.py` impose
  `--use-gl=angle --use-angle=d3d11` : 60 i/s sur le GPU même sans fenêtre.
- **Saut de curseur** : poser le curseur au point de départ 150 ms *avant* `mouse.down`, sinon Unity compte
  le déplacement vers ce point dans le geste (chaque drag ajoutait une rotation parasite ; 5 × 20 px ≠ 100 px).
  Avec ce délai les drags s'additionnent exactement, sans inertie ni dépendance à la vitesse.
- Un glissement trop rapide = clic de sélection (surbrillance verte). `drag()` bouge par pas de 20 px / 30 ms.
- Sensibilité ≈ 2300 px par tour (calibrée par retour de silhouette, `capture_config.json`) ; les grands angles
  sont découpés en segments qui restent dans la zone 3D (`drag_total`).
- Structure : page Angular > iframe `unityIframe` (iframe.voka.io/model : barre d'outils droite, boutons
  réinitialiser / annuler / rétablir = `.button` en haut à gauche, bulle « Astuces ») > iframe
  `unity/viewer/index.html` (canvas Unity). Cross-origin : on n'agit qu'à la souris, mais Playwright voit les DOM.
- Chargement long (20-60 s) : on attend que le canvas Unity soit dimensionné (largeur > 300) et que le texte
  « % / Finalisation » ait disparu, puis que la silhouette soit stable.
- Le bouton « réinitialiser la vue » remet toujours la même orientation et le même zoom : toutes les vues d'un
  stade sont comparables entre elles (l'échelle entre stades, elle, est normalisée par le viewer → CRL littérature).
- **Visual hull** : le sens de rotation du viewer est inversé par rapport au repère supposé, et l'élévation réelle
  diffère de la consigne ; `silhouettes.py` teste signes, échelles d'angle et pivot vertical, et garde l'intersection
  la plus volumineuse (`Jxx_hull_params.json`). Approximation orthographique d'une caméra perspective : bon pour
  proportions et placement, pas pour les détails. Reste dans `captures/` (privé), jamais dans `out/`.
- Le profil Chromium de connexion est dans `.voka_profile/` (ignoré par git) ; 1re connexion avec `--headed`.
