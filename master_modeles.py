# -*- coding: utf-8 -*-
"""Chaîne CS10→CS23 : assemble les dossiers embryons_3D/modeles/<CS>*/ (Hikspoors, vidéo, reconstructions, brouillons atlas) en un manifest
multi-stades pour blender_build_scene.py (shape keys entre stades sur les structures de même nom, slider « stage »).

    python embryo3d/master_modeles.py [--sortie embryons_3D/master_CS10-CS23.json] [--stades CS10,CS11,...,CS23]
                                      [--priorite hikspoors,video,recon,brouillon] [--fusion] [--sans-brouillons] [--publiables-seuls]
    "C:/Program Files/Blender Foundation/Blender 5.2/blender.exe" -b -P embryo3d/blender_build_scene.py -- embryons_3D/master_CS10-CS23.json embryons_3D/embryon_CS10-CS23.blend

Par stade : les dossiers dont le manifest porte ce stade sont classés par --priorite (suffixe du dossier : _hikspoors, _video, _recon, sans
suffixe = brouillon atlas ou livraison directe). Sans --fusion, seul le premier est pris ; avec --fusion, les structures des dossiers suivants
dont le nom canonique manque encore sont ajoutées (ex. enveloppe, membres, yeux de la vidéo à côté du cœur Hikspoors). Tous les dossiers
doivent être dans le repère des modèles (mm, Z crânial, Y dorsal, X gauche) : l'assemblage ne recale pas, il signale les écarts d'étendue.
Les dossiers « publiable": false » (Hikspoors par défaut : atlas d'auteurs retravaillé) sont pris : la scène Blender reste locale ;
--publiables-seuls les écarte. La publication sur le site est filtrée par agregateur.py, pas ici.
Sortie : {"stages":[{"stage","units","source","structures":[{"name","file","collection","color","alpha","confiance"}]}], "morph":[noms],
"stades_manquants":[...]} ; les stades sans modèle sont sautés (le morphing interpole alors entre les stades voisins).
"""
import argparse
import json
import os
import sys

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ICI = os.path.dirname(os.path.abspath(__file__))
CWD0 = os.getcwd()
os.chdir(RACINE)
sys.path.insert(0, ICI)
from agregateur import AXE, SYSTEMES, VUE_SEULE, MODELES, norm

LIBELLE = {i: l for i, l, _ in SYSTEMES + VUE_SEULE}
ALPHA = {"enveloppe": 0.25, "amnios": 0.2, "pericarde": 0.3, "coelome": 0.3, "vesicule_vitelline": 0.5, "cavites_cardiaques": 0.6}
# structures « blob » raisonnables pour des shape keys (Shrinkwrap) ; les tubes fins et les pièces multiples restent affichées par stade
MORPHABLES = {"enveloppe", "foie", "coeur", "cavites_cardiaques", "tube_neural", "encephale", "intestin", "estomac", "poumons", "vesicule_vitelline",
              "septum_transversum", "yeux", "vesicules_otiques", "pericarde", "myocarde_coeur", "myocarde_ventricule_gauche", "myocarde_ventricule_droit",
              "myocarde_oreillette_gauche", "myocarde_oreillette_droite", "myocarde_voie_efferente", "myocarde_ventricules", "myocarde_oreillettes",
              "myocarde_tube_cardiaque", "notochorde", "intestin_anterieur", "intestin_moyen", "intestin_posterieur", "oesophage", "duodenum"}


def rang(nom_dossier, priorite):
    suf = nom_dossier.split("_", 1)[1].lower() if "_" in nom_dossier else "brouillon"
    for i, p in enumerate(priorite):
        if p == suf or (p == "brouillon" and "_" not in nom_dossier):
            return i
    return len(priorite)


def dossiers(cs, priorite, sans_brouillons, publiables_seuls=False):
    out = []
    if not os.path.isdir(MODELES):
        return out
    for n in sorted(os.listdir(MODELES)):
        f = os.path.join(MODELES, n, "manifest.json")
        if not os.path.exists(f) or norm(n.split("_")[0]) != norm(cs):
            continue
        try:
            m = json.load(open(f, encoding="utf-8"))
        except Exception:
            continue
        if norm(m.get("stade", cs)) != norm(cs) or (publiables_seuls and m.get("publiable", True) is False):
            continue
        if sans_brouillons and m.get("statut") == "brouillon":
            continue
        out.append((rang(n, priorite), n, m))
    return [(n, m) for _, n, m in sorted(out, key=lambda t: (t[0], t[1]))]


def etendue(dossier, m):
    try:
        import numpy as np, trimesh
        lo, hi = [], []
        for x in m.get("structures", [])[:40]:
            f = os.path.join(MODELES, dossier, x.get("fichier", ""))
            if os.path.exists(f):
                b = trimesh.load(f, force="mesh", process=False).bounds
                lo.append(b[0]); hi.append(b[1])
        if lo:
            return (np.max(hi, 0) - np.min(lo, 0)).round(2).tolist()
    except Exception:
        pass
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sortie", default=os.path.join("embryons_3D", "master_CS10-CS23.json"))
    ap.add_argument("--stades", default=",".join(cs for cs, *_ in AXE if 10 <= norm(cs) <= 23))
    ap.add_argument("--priorite", default="hikspoors,video,recon,brouillon")
    ap.add_argument("--fusion", action="store_true"); ap.add_argument("--sans-brouillons", action="store_true")
    ap.add_argument("--publiables-seuls", action="store_true", help="écarter les dossiers publiable=false (Hikspoors par défaut)")
    a = ap.parse_args()
    if not os.path.isabs(a.sortie) and not a.sortie.startswith("embryons_3D"):
        a.sortie = os.path.normpath(os.path.join(CWD0, a.sortie))
    priorite = [p.strip().lower() for p in a.priorite.split(",") if p.strip()]
    base_out = os.path.dirname(os.path.abspath(a.sortie)) or "."
    master = {"stages": [], "stades_manquants": [], "morph": [], "priorite": priorite, "fusion": a.fusion,
              "repere": "modèles : mm, Z crânial, Y dorsal, X gauche (direct)"}
    presence = {}
    for cs in [c.strip() for c in a.stades.split(",") if c.strip()]:
        ds = dossiers(cs, priorite, a.sans_brouillons, a.publiables_seuls)
        if not ds:
            master["stades_manquants"].append(cs); print("%-5s : aucun modèle" % cs); continue
        st = {"stage": cs, "units": "mm", "source": [], "structures": []}
        noms = set()
        for k, (n, m) in enumerate(ds):
            if k > 0 and not a.fusion:
                break
            ajout = 0
            for x in m.get("structures", []):
                nom = x.get("nom")
                if not nom or nom in noms or x.get("confiance") == "faible":
                    continue
                f = os.path.join(MODELES, n, x.get("fichier", ""))
                if not os.path.exists(f):
                    continue
                noms.add(nom); ajout += 1
                st["structures"].append({"name": nom, "file": os.path.relpath(f, base_out).replace(os.sep, "/"),
                                         "collection": LIBELLE.get(x.get("systeme"), x.get("systeme", "Divers")), "color": x.get("couleur") or [0.7, 0.7, 0.7],
                                         "alpha": ALPHA.get(nom, 1.0), "confiance": x.get("confiance", "moyenne"), "dossier": n})
            st["source"].append({"dossier": n, "statut": m.get("statut", ""), "publiable": m.get("publiable", True) is not False,
                                 "structures": ajout, "etendue_mm": etendue(n, m)})
            if a.fusion and k > 0 and ajout:
                e0, e1 = st["source"][0]["etendue_mm"], st["source"][-1]["etendue_mm"]
                if e0 and e1 and max(e0) and abs(max(e1) / max(e0) - 1) > 0.25:
                    print("   ATTENTION %s : étendue %s (%s) vs %s (%s) — sources non recalées entre elles" % (cs, e0, ds[0][0], e1, n))
        master["stages"].append(st)
        for nom in noms:
            presence.setdefault(nom, []).append(cs)
        print("%-5s : %s" % (cs, " + ".join("%s (%d)" % (s["dossier"], s["structures"]) for s in st["source"])))
    master["morph"] = sorted(n for n, cs_ in presence.items() if n in MORPHABLES and len(cs_) >= 2)
    master["presence"] = {n: cs_ for n, cs_ in sorted(presence.items())}
    os.makedirs(base_out, exist_ok=True)
    json.dump(master, open(a.sortie, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print("-> %s : %d stades, %d structures morphées (%s), stades manquants : %s" % (a.sortie, len(master["stages"]), len(master["morph"]),
                                                                                   ", ".join(master["morph"]) or "aucune", ", ".join(master["stades_manquants"]) or "aucun"))


if __name__ == "__main__":
    main()
