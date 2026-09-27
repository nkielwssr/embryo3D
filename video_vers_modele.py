# -*- coding: utf-8 -*-
"""Nos reconstructions vidéo (pipeline.py, dossier CSxx_f4v/out) -> modèle au format de l'agrégateur, dans le même repère et la même
nomenclature que les modèles Hikspoors, pour entrer dans la chaîne CS10→CS23 (master_modeles.py).

    python embryo3d/video_vers_modele.py CS19_f4v [--sortie embryons_3D/modeles/CS19_video] [--sans-miroir] [--garder-brouillons]

Repère : le pipeline écrit X = gauche→droite, Y = ventral→dorsal, Z = bas→haut, ce qui est un repère INDIRECT (miroir) ; la convention
des modèles est directe (X = Y×Z = gauche anatomique). Le passage est donc X -> -X (faces réorientées), sauf --sans-miroir.
Noms : table « video » de hikspoors_nomenclature.json (snc -> tube_neural, coeur_detoure -> coeur, vaisseaux_aorte -> aorte_dorsale…) ;
un label absent de la table garde son nom, système déduit de la collection. Deux labels vers le même nom canonique : le plus fin l'emporte
(coeur_detoure > coeur, moelle_rachidienne > moelle…). Les labels de confiance « faible » sont écartés sauf --garder-brouillons.
"""
import argparse
import json
import os
import sys

import numpy as np
import trimesh

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ICI = os.path.dirname(os.path.abspath(__file__))
CWD0 = os.getcwd()
os.chdir(RACINE)
sys.path.insert(0, ICI)
from hikspoors_modele import NOMENCLATURE, stade_norme, planche_controle, CRL_MM


def abs0(p):
    return None if p is None else (p if os.path.isabs(p) else os.path.normpath(os.path.join(CWD0, p)))

COLLECTION_SYSTEME = {"Systeme nerveux": "nerveux", "SNC": "nerveux", "Organes des sens": "nerveux", "Organes": "digestif", "Digestif": "digestif",
                      "Squelette": "os_cartilages", "Vaisseaux": "vasculaire", "Cardio": "vasculaire", "Cavites": "annexes", "Enveloppe": "derme",
                      "Membres": "appendiculaire", "Muscles": "muscles", "Divers": "autre"}
PRIORITE = {"coeur_detoure": 2, "coeur": 1, "moelle_rachidienne": 2, "moelle": 1, "moelle_morph": 0, "encephale_morph": 1, "snc": 1}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dossier", help="dossier du stade (CS19_f4v) contenant out/manifest.json")
    ap.add_argument("--sortie"); ap.add_argument("--sans-miroir", action="store_true"); ap.add_argument("--garder-brouillons", action="store_true")
    ap.add_argument("--nomenclature", default=NOMENCLATURE)
    a = ap.parse_args()
    if not os.path.isdir(a.dossier):
        a.dossier = abs0(a.dossier)
    a.sortie = abs0(a.sortie)
    man_in = json.load(open(os.path.join(a.dossier, "out", "manifest.json"), encoding="utf-8"))
    cs = stade_norme(man_in.get("stage") or os.path.basename(a.dossier.rstrip("/\\")).split("_")[0].split(".")[0])
    sortie = a.sortie or os.path.join("embryons_3D", "modeles", cs + "_video")
    os.makedirs(sortie, exist_ok=True)
    table = json.load(open(a.nomenclature, encoding="utf-8")).get("video", {})
    choix = {}
    for s in man_in["structures"]:
        lab = s["name"]
        if s.get("confiance") == "faible" and not a.garder_brouillons:
            continue
        nom, syst = table.get(lab, [lab, COLLECTION_SYSTEME.get(s.get("collection", ""), "autre")])
        pr = PRIORITE.get(lab, 1)
        if nom in choix and choix[nom][0] >= pr:
            continue
        choix[nom] = (pr, lab, syst, s)
    # centre : enveloppe si présente, sinon ensemble
    meshes = {}
    for nom, (pr, lab, syst, s) in choix.items():
        f = os.path.join(a.dossier, "out", s["file"])
        if not os.path.exists(f):
            print("  absent :", f); continue
        m = trimesh.load(f, force="mesh", process=False)
        V = np.asarray(m.vertices, float).copy(); F = np.asarray(m.faces).copy()
        if not a.sans_miroir:
            V[:, 0] = -V[:, 0]; F = F[:, [0, 2, 1]]
        meshes[nom] = (trimesh.Trimesh(V, F, process=False), lab, syst, s)
    if not meshes:
        raise SystemExit("aucun maillage lu")
    ref = meshes.get("enveloppe", next(iter(meshes.values())))[0]
    centre = (ref.bounds[0] + ref.bounds[1]) / 2
    manifest = {"stade": cs, "source": "reconstruction 3D depuis les vidéos de coupes ehd.org (pipeline.py, session du 24/09), spécimen typique du stade",
                "session": "video_vers_modele", "unites": "mm",
                "repere": "Z crânial, Y dorsal, X = Y×Z = GAUCHE anatomique (repère direct" + (", miroir X appliqué au repère du pipeline)" if not a.sans_miroir else ", sans miroir)"),
                "statut": "recon", "publiable": True, "licence": "", "specimen": "",
                "dimensions": {"longueur_atlas_mm": man_in.get("greatest_length_mm_assumed"), "crl_typique_mm": CRL_MM.get(cs)},
                "notes": "Échelle = plus grande longueur typique du stade (CRL_MM de meshexport.py), pas la taille du spécimen. Segmentation heuristique sur la densité.",
                "structures": []}
    structs_mm = {}
    for nom, (m, lab, syst, s) in sorted(meshes.items(), key=lambda kv: (kv[1][2], kv[0])):
        m.apply_translation(-centre)
        fn = nom + ".ply"; m.export(os.path.join(sortie, fn))
        col = s.get("color", [0.7, 0.7, 0.7])
        manifest["structures"].append({"nom": nom, "nom_fr": nom.replace("_", " "), "nom_source": [lab], "systeme": syst, "fichier": fn, "couleur": col,
                                       "confiance": s.get("confiance", "moyenne"), "volume_mm3": s.get("volume_mm3"), "faces": int(len(m.faces)),
                                       "note": "label pipeline « %s »" % lab})
        structs_mm[nom] = {"meshes": [m], "couleur": col}
    json.dump(manifest, open(os.path.join(sortie, "manifest.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    planche_controle(structs_mm, os.path.join(sortie, "controle.png"), "%s — reconstruction vidéo — %d structures" % (cs, len(structs_mm)))
    try:
        import verif_maillages
        rel = os.path.relpath(sortie, os.path.join("embryons_3D", "modeles"))
        if not rel.startswith(".."):
            verif_maillages.main(rel)
    except Exception as e:
        print("   verif_maillages :", e)
    print("-> %s : %d structures (%s)" % (sortie, len(manifest["structures"]), ", ".join(sorted(structs_mm))))


if __name__ == "__main__":
    main()
