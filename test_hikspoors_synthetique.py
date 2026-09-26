# -*- coding: utf-8 -*-
"""Test sans données : embryon synthétique nommé comme dans les PDF Hikspoors (Neural_tube, Myocardium_LV, gut_yolk_sac bicolore, scale_cube…),
placé dans un repère tourné au hasard en unités de 1,078 µm, puis hikspoors_modele.py doit retrouver le repère anatomique, l'échelle, les noms.

    python embryo3d/test_hikspoors_synthetique.py [dossier_de_travail]     (défaut : dossier temporaire, effacé à la fin)
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import trimesh

ICI = os.path.dirname(os.path.abspath(__file__))


def ell(c, r, sub=3):
    m = trimesh.creation.icosphere(subdivisions=sub, radius=1.0); m.apply_scale(r); m.apply_translation(c); return m


def tube(pts, r, n=24):
    parts = [trimesh.creation.cylinder(radius=r, segment=[a, b], sections=n) for a, b in zip(pts[:-1], pts[1:])]
    parts += [ell(pts[0], r), ell(pts[-1], r)]
    return trimesh.util.concatenate(parts)


def generer(out, graine=7, echelle=1.0):
    """anatomie CS13-like en mm (x gauche, y dorsal, z crânial) -> GLB nommés à la Hikspoors, tournés, en unités de 1,078 µm"""
    os.makedirs(out, exist_ok=True)
    rs = np.random.RandomState(graine)
    z = np.linspace(-2.0, 2.0, 30)
    neural = trimesh.util.concatenate([tube(np.c_[np.zeros_like(z), 0.9 - 0.15 * z ** 2, z], 0.18), ell([0, 0.75, 2.05], [0.45, 0.42, 0.5]), ell([0, 0.85, 1.55], [0.32, 0.30, 0.35])])
    gut = tube(np.c_[np.zeros(20), 0.35 - 0.05 * np.linspace(-1.8, 1.4, 20) ** 2, np.linspace(-1.8, 1.4, 20)], 0.12)
    parts = {"Neural_tube": neural, "Myocardium_LV": ell([0.35, -0.25, 0.45], [0.42, 0.40, 0.42]), "Myocardium_RV": ell([-0.35, -0.25, 0.45], [0.38, 0.36, 0.40]),
             "Myocardium_LA": ell([0.28, 0.05, 0.95], [0.25, 0.22, 0.22]), "Myocardium_RA": ell([-0.28, 0.05, 0.95], [0.25, 0.22, 0.22]),
             "Lumen_LV": ell([0.35, -0.25, 0.45], [0.25, 0.24, 0.25]), "Lumen_RV": ell([-0.35, -0.25, 0.45], [0.22, 0.21, 0.23]),
             "OFT_myocardium": tube(np.array([[0.0, -0.3, 0.8], [0.0, -0.1, 1.3], [0.0, 0.2, 1.6]]), 0.16), "Sinus_venosus": ell([0.0, 0.1, 0.0], [0.3, 0.2, 0.2]),
             "Liver": ell([0.05, -0.15, -0.75], [0.55, 0.45, 0.5]), "gut_yolk_sac": trimesh.util.concatenate([gut, ell([0.0, -1.0, -1.2], [0.6, 0.5, 0.6])]),
             "Septum_transversum": ell([0.0, -0.25, -0.25], [0.5, 0.15, 0.12]),
             "Umbilical_vein_L": tube(np.array([[0.7, -0.6, -1.4], [0.5, -0.4, -0.8], [0.3, -0.2, -0.3]]), 0.06),
             "Umbilical_vein_R": tube(np.array([[-0.7, -0.6, -1.4], [-0.5, -0.4, -0.8], [-0.3, -0.2, -0.3]]), 0.06),
             "Dorsal_aorta_left": tube(np.c_[0.15 * np.ones(20), 0.55 * np.ones(20), np.linspace(-1.9, 1.2, 20)], 0.05),
             "Dorsal_aorta_right": tube(np.c_[-0.15 * np.ones(20), 0.55 * np.ones(20), np.linspace(-1.9, 1.2, 20)], 0.05),
             "Aortic_arch_3_left": tube(np.array([[0.0, 0.0, 1.6], [0.2, 0.3, 1.75], [0.15, 0.55, 1.6]]), 0.05),
             "scale_cube_200um": trimesh.creation.box(extents=[0.2, 0.2, 0.2]).apply_translation([2.0, 2.0, 2.0]),
             "Somites": trimesh.util.concatenate([ell([0.28 * s, 0.75, zz], [0.09, 0.09, 0.06]) for zz in np.linspace(-1.6, 1.0, 12) for s in (1, -1)]),
             "Mystery_blob": ell([0.9, 0.9, -1.5], [0.15, 0.15, 0.15])}
    q = rs.normal(size=4); q /= np.linalg.norm(q)
    R0 = trimesh.transformations.quaternion_matrix(q)[:3, :3]
    t0 = np.array([12345.0, -777.0, 4200.0]); K = echelle * 1000.0 / 1.078
    for nom, m in parts.items():
        V = (R0 @ (np.asarray(m.vertices) * K).T).T + t0
        mm = trimesh.Trimesh(V, m.faces, process=False)
        if nom == "gut_yolk_sac":
            fc = np.zeros((len(mm.faces), 4), np.uint8) + 255
            cy = (R0 @ (np.array([0.0, -1.0, -1.2]) * K)) + t0
            d = np.linalg.norm(mm.triangles_center - cy, axis=1) < 0.62 * K
            fc[:, :3] = (174, 174, 174); fc[d, :3] = (120, 120, 120); mm.visual.face_colors = fc
        mm.export(os.path.join(out, nom + ".glb"))
    return R0


def main():
    tmp = sys.argv[1] if len(sys.argv) > 1 else tempfile.mkdtemp(prefix="hikspoors_test_")
    glb, sortie = os.path.join(tmp, "glb13"), os.path.join(tmp, "out_CS13")
    R0 = generer(glb)
    cmd = [sys.executable, os.path.join(ICI, "hikspoors_modele.py"), "CS13", "--glb", glb, "--separer-couleurs", "--um-par-unite", "1.078", "--sortie", sortie]
    print(" ".join(cmd)); subprocess.run(cmd, check=True)
    man = json.load(open(os.path.join(sortie, "manifest.json"), encoding="utf-8"))
    rap = json.load(open(os.path.join(sortie, "rapport.json"), encoding="utf-8"))
    R = np.array(rap["orientation"]["R_lignes_XYZ_dans_repere_source"])
    angle = float(np.degrees(np.arccos(np.clip((np.trace(R @ R0) - 1) / 2, -1, 1))))
    noms = {x["nom"] for x in man["structures"]}
    erreurs = []
    if angle > 10:
        erreurs.append("orientation : %.1f° d'écart avec le repère vrai" % angle)
    for n in ("tube_neural", "myocarde_ventricule_gauche", "cavite_ventricule_droit", "foie", "intestin", "vesicule_vitelline", "veine_ombilicale_gauche",
              "aorte_dorsale_droit", "arc_aortique_3_gauche", "somites", "septum_transversum", "coeur_sinus_veineux", "coeur", "arteres", "veines", "cavites_cardiaques"):
        if n not in noms:
            erreurs.append("structure attendue absente : " + n)
    if "scale_cube_200um" in noms or "scale_cube" in " ".join(noms):
        erreurs.append("le cube d'échelle n'a pas été ignoré")
    if rap["non_reconnus"] != ["Mystery_blob"]:
        erreurs.append("non reconnus attendus ['Mystery_blob'], obtenu %s" % rap["non_reconnus"])
    h = rap["dimensions"]["hauteur_mm"]
    if not 4.3 < h < 5.2:
        erreurs.append("hauteur %.2f mm hors de [4.3, 5.2] : échelle fausse" % h)
    if not all(d["ok"] for d in rap["orientation"]["indices"] if d["poids"] > 0):
        erreurs.append("indices d'orientation non tous vérifiés : %s" % [d["indice"] for d in rap["orientation"]["indices"] if not d["ok"] and d["poids"] > 0])
    vg = trimesh.load(os.path.join(sortie, "myocarde_ventricule_gauche.ply")); vd = trimesh.load(os.path.join(sortie, "myocarde_ventricule_droit.ply"))
    if vg.centroid[0] <= vd.centroid[0]:
        erreurs.append("ventricule gauche à X <= ventricule droit : miroir")
    if man["statut"] != "externe" or "BY-NC-SA" not in man["licence"]:
        erreurs.append("manifest : statut/licence inattendus")
    print("\nécart d'orientation %.2f°, hauteur %.2f mm, %d structures" % (angle, h, len(noms)))
    if len(sys.argv) < 2:
        shutil.rmtree(tmp, ignore_errors=True)
    if erreurs:
        print("ÉCHEC :\n  " + "\n  ".join(erreurs)); sys.exit(1)
    print("OK : test synthétique Hikspoors réussi")


if __name__ == "__main__":
    main()
