# -*- coding: utf-8 -*-
"""Test sans données, au format réel de la banque hikspoors_maastricht : embryon synthétique EN C façon CS13 (tube neural qui fait le
tour du cœur, queue relevée), noms réels des PDF (inventaire du relais : R_cardinal_vein, pericard, PAAs_3, L_Vit_vein, loop_wire…),
npz FACESET_<nom>|v/|f/|c/|p/|fc en unités PDF (1,101 µm/unité, repère tourné), calage/<nom>.json (matrice_pdf_vers_mm) vers le repère
du modèle. Vérifie :
  1. calage de la banque pris par défaut : coordonnées du modèle exactes, pas de recentrage, noms, unions sans veines cardinales ni
     péricarde dans « coeur », publiable = false ;
  2. --sans-calage : orientation automatique juste sur l'embryon en C (écart < 15°, gauche/droite, indices pondérés tous vérifiés) ;
  3. CS18 _NEWvalves : « gut » repris de la version d'origine (même repère), CCS / CCS_1 non doublés ;
  4. CS14 calé sur le modèle CS15 : signalé « inter_stades » ; CS12 calage avec réflexion : faces réorientées (volumes positifs).

    python embryo3d/test_hikspoors_banque.py [dossier_de_travail]     (défaut : dossier temporaire, effacé à la fin)
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
sys.path.insert(0, ICI)
from test_hikspoors_synthetique import ell, tube          # noqa: E402

UM = 1.101                                                  # µm par unité PDF (CS13, calage de la session « PDF 3D »)


def arc(phis_deg, facteur=1.0, x=0.0, ay=1.35, az=1.95):
    """points de la courbe en C (plan YZ) : dorsal = +Y à mi-hauteur, tête en haut qui revient en avant du cœur, queue relevée en bas"""
    p = np.radians(np.asarray(phis_deg, float))
    return np.c_[np.full(len(p), x), facteur * ay * np.cos(p), facteur * az * np.sin(p)]


def anatomie():
    """{nom réel: (maillage en mm dans le repère du modèle, couleur |c 0-1, couleurs par face 0-255 ou None)}"""
    P = {}
    def ajoute(nom, m, c, fc=None):
        P[nom] = (m, np.asarray(c, float), fc)
    spine = tube(arc(np.linspace(-140, 70, 40)), 0.11)
    cerveau = [ell(arc([a])[0], r) for a, r in ((85, [0.28, 0.26, 0.26]), (115, [0.33, 0.3, 0.3]), (145, [0.35, 0.3, 0.32]))]
    ajoute("neural_tube", trimesh.util.concatenate([spine] + cerveau), [1.0, 0.958, 0.64])
    som = [ell(arc([a], 0.97, s * 0.22)[0], [0.08, 0.08, 0.06], sub=2) for a in np.linspace(-120, 62, 22) for s in (1, -1)]
    m = trimesh.util.concatenate(som)
    fc = np.zeros((len(m.faces), 3), np.uint8) + (255, 224, 218)
    fc[: len(fc) // 5] = (0, 0, 0)                          # surface texturée : deux teintes (000000, ffe0da) comme dans le PDF
    ajoute("somites", m, [0.5, 0.5, 0.5], fc)
    ajoute("dorsal_aorta", trimesh.util.concatenate([tube(arc(np.linspace(-110, 60, 30), 0.80, s * 0.12), 0.06) for s in (1, -1)]), [1, 0, 0])
    for c, s in (("L", 1), ("R", -1)):
        ajoute(c + "_cardinal_vein", tube(arc(np.linspace(-70, 35, 20), 0.83, s * 0.42), 0.05), [0, 0, 0.447])
        ajoute(c + "_common_card_vein", tube(np.array([arc([10], 0.83, s * 0.42)[0], [s * 0.35, 0.3, -0.05], [s * 0.25, -0.2, -0.15]]), 0.05), [0, 0, 0.447])
        ajoute(c + "_Vit_vein", tube(np.array([[s * 0.15, -0.9, -1.3], [s * 0.2, -0.5, -0.7], [s * 0.2, -0.2, -0.2]]), 0.05), [0.704, 0.945, 0.945])
        ajoute(c + "_Umb_vein", tube(np.array([[s * 0.55, -1.0, -1.1], [s * 0.45, -0.6, -0.6], [s * 0.3, -0.25, -0.2]]), 0.05), [1, 0, 1])
        ajoute(c + "_HCC", tube(np.array([[s * 0.2, -0.4, -0.5], [s * 0.2, -0.25, -0.25]]), 0.04), [0.69, 0.424, 0.933])
    for c, x, rgb in (("LV", 0.28, (66, 219, 86)), ("RV", -0.28, (165, 255, 140))):
        ajoute("lumen_" + c, ell([x, -0.75, 0.0], [0.26, 0.24, 0.26]), [0, 0.83, 0.91])
        m = ell([x, -0.75, 0.0], [0.32, 0.3, 0.32])
        fc = np.zeros((len(m.faces), 3), np.uint8) + rgb
        fc[np.abs(m.triangles_center[:, 0]) < 0.12] = (127, 173, 87)          # courbure interne
        ajoute("myocard_" + c, m, [0.5, 0.5, 0.5], fc)
    for c, x in (("LA", 0.25), ("RA", -0.25)):
        ajoute("lumen_" + c, ell([x, -0.25, 0.45], [0.18, 0.15, 0.15]), [0.28, 0.79, 1.0])
    m = trimesh.util.concatenate([ell([x, -0.25, 0.45], [0.22, 0.19, 0.19]) for x in (0.25, -0.25)])
    fc = np.zeros((len(m.faces), 3), np.uint8) + (0, 115, 8); fc[: len(fc) // 6] = (99, 115, 94)
    ajoute("myocard_atriums", m, [0.5, 0.5, 0.5], fc)
    oft = np.array([[0, -0.8, 0.3], [0, -0.6, 0.8], [0, -0.45, 1.05]])
    ajoute("lumen_OFT", tube(oft, 0.07), [0.737, 1, 1]); ajoute("myocard_OFT", tube(oft, 0.12), [0.682, 1, 0.698])
    ajoute("OFT_cardiac_jelly", tube(oft, 0.1), [0.918, 0.651, 1])
    ajoute("aortic_sac", ell([0, -0.45, 1.12], [0.12, 0.1, 0.1]), [1, 0.584, 0])
    for k, n in enumerate((2, 3, 4)):
        ajoute("PAAs_%d" % n, trimesh.util.concatenate([tube(np.array([[s * 0.08, -0.45, 1.12 + 0.06 * k], [s * 0.3, 0.1, 1.25 + 0.06 * k],
                                                                       [s * 0.12, 0.69, 1.2 + 0.05 * k]]), 0.04) for s in (1, -1)]), [0.66, 0.385, 0])
    ajoute("myocard_AV_canal", ell([0, -0.45, 0.2], [0.15, 0.12, 0.1]), [1, 1, 0])
    ajoute("lumen_venous_sinus", ell([0, -0.2, -0.15], [0.3, 0.15, 0.12]), [0, 0, 1])
    ajoute("myocard_venous_sinus", ell([0, -0.2, -0.15], [0.34, 0.18, 0.15]), [0.486, 0.64, 0.535])
    ajoute("SAN", ell([-0.2, -0.15, -0.05], [0.05, 0.05, 0.05], sub=2), [1, 1, 0.639])
    ajoute("DMP", ell([0, -0.1, 0.2], [0.06, 0.06, 0.06], sub=2), [1, 0, 0.8])
    ajoute("sup_endocardial_cushion", ell([0, -0.5, 0.25], [0.08, 0.06, 0.05], sub=2), [1, 0.655, 0.1])
    ajoute("inf_endocardial_cushion", ell([0, -0.58, 0.12], [0.08, 0.06, 0.05], sub=2), [1, 0.5, 0])
    ajoute("epicard", ell([0, -0.6, 0.3], [0.68, 0.5, 0.66]), [0.76, 0.485, 0.301])
    ajoute("pericard", ell([0, -0.6, 0.3], [0.8, 0.64, 0.8]), [0.82, 0.316, 0.082])
    ajoute("pericardial_reflection", ell([0, -0.2, 0.72], [0.2, 0.08, 0.06], sub=2), [1, 0.808, 0.675])
    m = tube(np.array([[0, -0.2, -0.15], [0.25, -0.75, 0.0], [-0.28, -0.75, 0.05], [0, -0.45, 1.05]]), 0.02)
    fc = np.zeros((len(m.faces), 3), np.uint8); fc[:] = (188, 255, 255); fc[len(fc) // 3:] = (61, 252, 255); fc[2 * len(fc) // 3:] = (0, 232, 255)
    ajoute("loop_wire_heart_tube", m, [0.5, 0.5, 0.5], fc)
    ajoute("pulmonary_vein", tube(np.array([[0, 0.05, 0.35], [0, -0.1, 0.42]]), 0.03), [0.044, 0.772, 0.687])
    ajoute("liver", ell([0, -0.55, -0.75], [0.55, 0.4, 0.4]), [0.866, 0.99, 0.495])
    ajoute("transverse_septum", ell([0, -0.45, -0.35], [0.55, 0.35, 0.12]), [0.99, 0.518, 0.412])
    gut = tube(arc(np.linspace(-100, 70, 30), 0.55), 0.09)
    yolk = ell([0, -1.35, -0.9], [0.45, 0.4, 0.4])
    m = trimesh.util.concatenate([gut, yolk])
    fc = np.zeros((len(m.faces), 3), np.uint8) + (174, 174, 174); fc[len(gut.faces):] = (120, 119, 119)
    ajoute("gut_yolk_sac", m, [0.5, 0.5, 0.5], fc)
    ajoute("scale_cube_200um", trimesh.creation.box(extents=[0.2, 0.2, 0.2]).apply_translation([1.6, 1.6, -2.2]), [0.6, 0.6, 0.6])
    return P


def rotation(graine):
    q = np.random.RandomState(graine).normal(size=4); q /= np.linalg.norm(q)
    return trimesh.transformations.quaternion_matrix(q)[:3, :3]


def ecrire_npz(chemin, P, R0, t0, um=UM, reflexion=False, omettre=(), renommer=None):
    """x_pdf = R0ᵀ (x_mm - t0) / s (miroir X d'abord si reflexion) ; renvoie la matrice 4x4 pdf -> mm du calage"""
    s = um / 1000.0
    D = np.diag([-1.0, 1, 1]) if reflexion else np.eye(3)
    arr = {}
    for nom, (m, c, fc) in P.items():
        if nom in omettre:
            continue
        n = (renommer or {}).get(nom, nom)
        V = ((np.asarray(m.vertices) - t0) @ R0) / s @ D.T
        k = "FACESET_" + n
        F = np.asarray(m.faces, np.int32)
        if reflexion:
            F = F[:, ::-1]                                    # maillage du PDF bien orienté dans son propre repère (normales sortantes)
        arr[k + "|v"] = V.astype(np.float32); arr[k + "|f"] = F; arr[k + "|c"] = np.asarray(c, np.float32)
        arr[k + "|p"] = np.array(n)
        if fc is not None:
            arr[k + "|fc"] = np.asarray(fc, np.uint8)
    np.savez(chemin, **arr)
    M = np.eye(4); M[:3, :3] = s * R0 @ D; M[:3, 3] = t0
    return M


def lancer(*args):
    cmd = [sys.executable, os.path.join(ICI, "hikspoors_modele.py")] + list(args)
    print("\n$ " + " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8")
    print(r.stdout[-3500:])
    if r.returncode:
        print(r.stderr[-3000:]); raise SystemExit("échec de hikspoors_modele.py")
    return r.stdout


def lire(sortie):
    return (json.load(open(os.path.join(sortie, "manifest.json"), encoding="utf-8")),
            json.load(open(os.path.join(sortie, "rapport.json"), encoding="utf-8")))


def main():
    tmp = sys.argv[1] if len(sys.argv) > 1 else tempfile.mkdtemp(prefix="hikspoors_banque_")
    banque = os.path.join(tmp, "hikspoors_maastricht")
    os.makedirs(os.path.join(banque, "u3d"), exist_ok=True); os.makedirs(os.path.join(banque, "calage"), exist_ok=True)
    P = anatomie()
    R0, t0 = rotation(11), np.array([0.31, -0.42, 0.17])          # origine du modèle hors du centre de la boîte : pas de recentrage attendu
    M = ecrire_npz(os.path.join(banque, "u3d", "Carnegie_Stage_13.npz"), P, R0, t0)
    json.dump({"matrice_pdf_vers_mm": M.tolist(), "modele": "CS13", "um_par_unite": UM}, open(os.path.join(banque, "calage", "Carnegie_Stage_13.json"), "w"))
    erreurs = []

    # 1. calage de la banque par défaut
    s1 = os.path.join(tmp, "CS13_calage")
    lancer("CS13", "--banque", banque, "--separer-couleurs", "--sortie", s1)
    man, rap = lire(s1)
    noms = {x["nom"] for x in man["structures"]}
    if not rap.get("calage") or os.path.basename(rap["calage"]["fichier"]) != "Carnegie_Stage_13.json":
        erreurs.append("1. calage de la banque non pris")
    foie = trimesh.load(os.path.join(s1, "foie.ply"), process=False)
    err = float(np.abs(np.asarray(foie.vertices) - np.asarray(P["liver"][0].vertices)).max())
    if err > 1e-3:
        erreurs.append("1. foie à %.3g mm de sa place dans le repère du modèle (recentrage ou calage mal appliqué)" % err)
    if rap["dimensions"]["recentre"]:
        erreurs.append("1. recentré malgré le calage vers le modèle")
    if abs(rap["echelle"]["um_par_unite"] - UM) > 1e-3:
        erreurs.append("1. échelle %.4g µm/unité, attendu %.4g" % (rap["echelle"]["um_par_unite"], UM))
    for n in ("veine_cardinale_gauche", "veine_cardinale_droit", "veine_cardinale_commune_gauche", "veine_vitelline_droit", "veine_ombilicale_gauche",
              "canal_hepatocardiaque_gauche", "arc_aortique_3", "pericarde", "epicarde", "reflexion_pericardique", "coussin_endocardique_inferieur",
              "gelee_cardiaque_voie_efferente", "cavite_sinus_veineux", "myocarde_sinus_veineux", "noeud_sinusal", "protrusion_mesenchymateuse_dorsale",
              "axe_tube_cardiaque", "intestin", "vesicule_vitelline", "courbure_interne_lv", "coeur", "veines", "arteres", "cavites_cardiaques"):
        if n not in noms:
            erreurs.append("1. structure attendue absente : " + n)
    if rap["non_reconnus"]:
        erreurs.append("1. non reconnus : %s" % rap["non_reconnus"])
    coeur = next(x for x in man["structures"] if x["nom"] == "coeur")
    intrus = [c for c in coeur["composantes"] if c.startswith(("veine", "pericarde", "reflexion", "axe_"))]
    if intrus:
        erreurs.append("1. union « coeur » contient %s" % intrus)
    if man["publiable"] is not False:
        erreurs.append("1. publiable devrait être false par défaut")
    if not os.path.exists(os.path.join(s1, "verif_maillages.json")) or not os.path.exists(os.path.join(s1, "controle.png")):
        erreurs.append("1. verif_maillages.json ou controle.png absent")
    if not all(d["ok"] for d in rap["orientation"]["indices"] if d["poids"] > 0):
        erreurs.append("1. indices non vérifiés avec le bon calage : %s" % [d["indice"] for d in rap["orientation"]["indices"] if not d["ok"] and d["poids"] > 0])
    if rap["orientation"].get("controle_auto", {}).get("ecart_deg", 99) > 15:
        erreurs.append("1. contrôle de l'orientation automatique contre le calage : %s" % rap["orientation"].get("controle_auto"))

    # 2. orientation automatique sur l'embryon en C
    s2 = os.path.join(tmp, "CS13_auto")
    lancer("CS13", "--banque", banque, "--separer-couleurs", "--sans-calage", "--sans-unions", "--um-par-unite", str(UM), "--sortie", s2)
    man2, rap2 = lire(s2)
    R = np.array(rap2["orientation"]["R_lignes_XYZ_dans_repere_source"])
    angle = float(np.degrees(np.arccos(np.clip((np.trace(R @ R0.T) - 1) / 2, -1, 1))))  # x_mm = s R0 x_pdf + t0 : R attendu = R0
    if angle > 15:
        erreurs.append("2. orientation automatique à %.1f° du repère vrai (candidat %s / %s)" % (angle, rap2["orientation"]["z"], rap2["orientation"]["y"]))
    if not all(d["ok"] for d in rap2["orientation"]["indices"] if d["poids"] > 0):
        erreurs.append("2. indices pondérés non tous vérifiés : %s" % [d["indice"] for d in rap2["orientation"]["indices"] if not d["ok"] and d["poids"] > 0])
    g = trimesh.load(os.path.join(s2, "veine_cardinale_gauche.ply")); d = trimesh.load(os.path.join(s2, "veine_cardinale_droit.ply"))
    if g.centroid[0] <= d.centroid[0]:
        erreurs.append("2. veine cardinale gauche à X <= droite : miroir")

    # 3. CS18 NEWvalves sans « gut », CCS -> CCS_1 ; l'originale a les deux
    M18 = ecrire_npz(os.path.join(banque, "u3d", "Carnegie_Stage_18_NEWvalves.npz"), P, R0, t0, omettre=("gut_yolk_sac",), renommer={"SAN": "CCS_1"})
    ecrire_npz(os.path.join(banque, "u3d", "Carnegie_Stage_18.npz"), P, R0, t0, renommer={"SAN": "CCS", "gut_yolk_sac": "gut"})
    json.dump({"matrice_pdf_vers_mm": M18.tolist()}, open(os.path.join(banque, "calage", "Carnegie_Stage_18_NEWvalves.json"), "w"))
    s3 = os.path.join(tmp, "CS18")
    lancer("CS18", "--banque", banque, "--sans-unions", "--sortie", s3)
    man3, rap3 = lire(s3)
    noms3 = {x["nom"]: x for x in man3["structures"]}
    if "intestin" not in noms3:
        erreurs.append("3. « gut » de la version d'origine non repris")
    if "systeme_conduction" not in noms3 or noms3["systeme_conduction"]["nom_source"] != ["CCS_1"]:
        erreurs.append("3. système de conduction : %s (attendu CCS_1 seul)" % (noms3.get("systeme_conduction") or {}).get("nom_source"))
    if rap3.get("version_origine", {}).get("ajoutees") != ["gut"]:
        erreurs.append("3. ajoutées depuis l'originale : %s, attendu ['gut']" % rap3.get("version_origine", {}).get("ajoutees"))

    # 4. CS14 calé sur le modèle CS15 ; CS12 avec réflexion
    M14 = ecrire_npz(os.path.join(banque, "u3d", "Carnegie_Stage_14.npz"), P, R0, t0, um=5.79)
    json.dump({"matrice_pdf_vers_mm": M14.tolist(), "modele": "CS15"}, open(os.path.join(banque, "calage", "Carnegie_Stage_14.json"), "w"))
    s4 = os.path.join(tmp, "CS14")
    lancer("CS14", "--banque", banque, "--sans-unions", "--sortie", s4)
    _, rap4 = lire(s4)
    if rap4["calage"].get("inter_stades") != "CS15":
        erreurs.append("4. calage CS14 -> CS15 non signalé (inter_stades = %s)" % rap4["calage"].get("inter_stades"))
    M12 = ecrire_npz(os.path.join(banque, "u3d", "Carnegie_Stage_12.npz"), P, R0, t0, um=4.864, reflexion=True)
    json.dump({"matrice_pdf_vers_mm": M12.T.tolist()}, open(os.path.join(banque, "calage", "Carnegie_Stage_12.json"), "w"))  # rangée par colonnes
    s5 = os.path.join(tmp, "CS12")
    lancer("CS12", "--banque", banque, "--sans-unions", "--sortie", s5)
    _, rap5 = lire(s5)
    foie12 = trimesh.load(os.path.join(s5, "foie.ply"), process=False)
    if not rap5["calage"]["reflexion"] or not rap5["calage"].get("transposee"):
        erreurs.append("4. CS12 : réflexion ou matrice rangée par colonnes non détectée")
    if float(np.abs(np.asarray(foie12.vertices) - np.asarray(P["liver"][0].vertices)).max()) > 1e-3 or foie12.volume <= 0:   # normales sortantes
        erreurs.append("4. CS12 : foie mal placé ou faces retournées (volume %.4g)" % foie12.volume)

    print("\n1. calage : écart foie %.2g mm ; 2. orientation auto sur l'embryon en C : %.1f° ; %d structures" % (err, angle, len(noms)))
    if len(sys.argv) < 2:
        shutil.rmtree(tmp, ignore_errors=True)
    if erreurs:
        print("ÉCHEC :\n  " + "\n  ".join(erreurs)); sys.exit(1)
    print("OK : test banque Hikspoors (embryon en C, calage json, NEWvalves, inter-stades, réflexion) réussi")


if __name__ == "__main__":
    main()
