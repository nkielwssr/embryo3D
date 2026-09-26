# -*- coding: utf-8 -*-
"""Contrôle d'erreurs des maillages d'un modèle livré (embryons_3D/modeles/<dossier>/).

    python embryo3d/verif_maillages.py CS10            # rapport console + modeles/CS10/verif_maillages.json

Par structure : faces, sommets, composantes, étanche (watertight), orientation cohérente, arêtes non-manifold,
faces dégénérées (aire nulle), faces dupliquées, sommets isolés, volume signé. Verdict : ok / a_corriger.
"""
import json
import os
import sys

import numpy as np
import trimesh

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(RACINE)


def charger(f):
    m = trimesh.load(f, force="mesh", process=False)
    return m


def verifier(m):
    faces = m.faces
    aires = m.area_faces
    arete = np.sort(m.edges, axis=1)
    _, cnt = np.unique(arete, axis=0, return_counts=True)
    nm = int((cnt > 2).sum())                       # arêtes partagées par plus de 2 faces
    bord = int((cnt == 1).sum())                    # arêtes de bord (trous)
    dup = len(faces) - len(np.unique(np.sort(faces, axis=1), axis=0))
    utilises = np.zeros(len(m.vertices), bool); utilises[faces.ravel()] = True
    comp = len(m.split(only_watertight=False)) if len(faces) < 400000 else -1
    r = {"faces": int(len(faces)), "sommets": int(len(m.vertices)), "composantes": comp,
         "etanche": bool(m.is_watertight), "orientation_coherente": bool(m.is_winding_consistent),
         "aretes_non_manifold": nm, "aretes_de_bord": bord, "faces_degenerees": int((aires < 1e-14).sum()),
         "faces_dupliquees": int(dup), "sommets_isoles": int((~utilises).sum()),
         "volume_mm3": round(float(m.volume), 5) if m.is_watertight else None}
    pb = []
    if not r["etanche"]: pb.append(f"non étanche ({bord} arêtes de bord)")
    if not r["orientation_coherente"]: pb.append("orientation incohérente")
    if nm: pb.append(f"{nm} arêtes non-manifold")
    if r["faces_degenerees"]: pb.append(f"{r['faces_degenerees']} faces dégénérées")
    if dup: pb.append(f"{dup} faces dupliquées")
    if r["sommets_isoles"]: pb.append(f"{r['sommets_isoles']} sommets isolés")
    r["problemes"] = pb
    r["verdict"] = "ok" if not pb else "a_corriger"
    return r


def main(dossier):
    d = os.path.join("embryons_3D", "modeles", dossier)
    man = json.load(open(os.path.join(d, "manifest.json"), encoding="utf-8"))
    out = {}
    for x in man.get("structures", []):
        f = os.path.join(d, x.get("fichier", ""))
        if not os.path.exists(f):
            out[x["nom"]] = {"verdict": "absent", "problemes": ["fichier absent"]}
            continue
        out[x["nom"]] = verifier(charger(f))
    ok = sum(1 for v in out.values() if v["verdict"] == "ok")
    res = {"dossier": dossier, "structures": out, "bilan": {"ok": ok, "total": len(out)}}
    json.dump(res, open(os.path.join(d, "verif_maillages.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(f"{dossier} : {ok}/{len(out)} maillages sans erreur")
    for n, v in out.items():
        if v["verdict"] != "ok":
            print(f"  {n:32s} {v.get('faces','-'):>7} faces  {v.get('composantes','-'):>4} comp.  " + " ; ".join(v["problemes"]))
    return res


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "CS10")
