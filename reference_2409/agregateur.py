# -*- coding: utf-8 -*-
"""
Agrégateur des trois sources « embryon 3D + temps » :
  A. Reconstruction 3D segmentée depuis nos vidéos de coupes  (CSxx_f4v/out, embryons_3D/master_manifest.json)
  B. Coupes histologiques légendées ehd.org                   (coupes embryos 9-23/embryo_stage_S/legendes)
  C. Référence externe VOKA J28→J91                           (embryo3d/reference/captures, proportions.json)

Scanne le disque, place tout sur un axe unique (jours post-fécondation / stade Carnegie / CRL),
calcule ce qui existe et ce qu'il reste à faire, et écrit embryons_3D/agregateur.html (autonome).

    python embryo3d/agregateur.py            # écrit embryons_3D/agregateur.html
    python embryo3d/agregateur.py --json     # écrit aussi embryons_3D/agregateur.json
    python embryo3d/agregateur.py --site D   # paquet autonome publiable dans D/ (index.html + glb/ rendus/ overlays/ voka/)
                                             # les hulls VOKA (captures/, privés) ne sont jamais copiés
"""
import shutil
import hashlib
import glob
import json
import os
import re
import sys
from datetime import date, datetime

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(RACINE)

SORTIE_HTML = os.path.join("embryons_3D", "agregateur.html")
SORTIE_JSON = os.path.join("embryons_3D", "agregateur.json")

# --------------------------------------------------------------------------------------
# Axe de temps commun : stades Carnegie (O'Rahilly & Müller 2010) puis semaines fœtales
# jours = post-fécondation ; CRL = plus grande longueur (mm)
# --------------------------------------------------------------------------------------
AXE = [
    # id, jours min, jours max, CRL min, CRL max, repères
    ("CS2", 2, 3, 0.1, 0.2, "Segmentation, morula"),
    ("CS3", 4, 5, 0.1, 0.2, "Blastocyste libre"),
    ("CS4", 6, 6, 0.1, 0.2, "Blastocyste attaché"),
    ("CS5", 7, 12, 0.1, 0.2, "Implantation, disque bilaminaire"),
    ("CS6", 13, 15, 0.2, 0.2, "Villosités choriales, ligne primitive"),
    ("CS7", 15, 17, 0.4, 0.4, "Notochorde"),
    ("CS8", 17, 19, 1.0, 1.5, "Gouttière neurale, fossette primitive"),
    ("CS9", 19, 21, 1.5, 2.5, "1-3 somites, plis neuraux"),
    ("CS10", 22, 23, 2.0, 3.5, "4-12 somites, fermeture neurale"),
    ("CS11", 23, 26, 2.5, 4.5, "13-20 somites, neuropore rostral fermé"),
    ("CS12", 26, 30, 3.0, 5.0, "21-29 somites, bourgeons MS"),
    ("CS13", 28, 32, 4.0, 6.0, "4 arcs, otocyste fermé, bourgeons MI"),
    ("CS14", 31, 35, 5.0, 7.0, "Cupule optique, fosse nasale"),
    ("CS15", 35, 38, 7.0, 9.0, "Palette de la main, cristallin"),
    ("CS16", 37, 42, 8.0, 11.0, "Pigment rétinien, palette du pied"),
    ("CS17", 42, 44, 11.0, 14.0, "Rayons digitaux, bourrelets auriculaires"),
    ("CS18", 44, 48, 13.0, 17.0, "Ossification débute, coudes"),
    ("CS19", 48, 51, 16.0, 18.0, "Tronc s'allonge et se redresse"),
    ("CS20", 51, 53, 18.0, 22.0, "Doigts séparés, MS fléchis"),
    ("CS21", 53, 54, 22.0, 24.0, "Doigts allongés, mains se rejoignent"),
    ("CS22", 54, 56, 23.0, 28.0, "Paupières, pavillon de l'oreille"),
    ("CS23", 56, 60, 27.0, 31.0, "Fin de l'embryon : tête arrondie, queue disparue"),
    ("F9", 61, 67, 40.0, 50.0, "Fœtus 9 sem. : yeux fermés, hernie ombilicale"),
    ("F10", 68, 74, 55.0, 61.0, "Fœtus 10 sem. : intestin rentré, visage humain"),
    ("F11", 75, 81, 68.0, 75.0, "Fœtus 11 sem. : sexe externe, cou"),
    ("F12", 82, 88, 80.0, 87.0, "Fœtus 12 sem. : ossification étendue"),
    ("F13", 89, 95, 95.0, 105.0, "Fœtus 13 sem. : fin du 1er trimestre"),
    ("F14-15", 96, 109, 120.0, 135.0, "2e trimestre : lanugo, mouvements coordonnés"),
    ("F16-17", 110, 123, 140.0, 155.0, "Ossification visible, sexe externe identifiable"),
    ("F18-19", 124, 137, 160.0, 175.0, "Vernix, mouvements perçus par la mère"),
    ("F20-21", 138, 151, 190.0, 200.0, "Sourcils, cheveux ; poids ≈ 460 g"),
    ("F22-23", 152, 165, 210.0, 220.0, "Peau ridée, alvéoles primitives"),
    ("F24-25", 166, 179, 230.0, 240.0, "Ongles ; surfactant débute"),
    ("F26-27", 180, 193, 250.0, 260.0, "3e trimestre : yeux ouverts"),
    ("F28-29", 194, 207, 270.0, 275.0, "Graisse sous-cutanée, pupilles réactives"),
    ("F30-31", 208, 221, 280.0, 290.0, "Peau lisse et rose"),
    ("F32-33", 222, 235, 300.0, 315.0, "Ongles aux orteils"),
    ("F34-35", 236, 249, 320.0, 335.0, "Membres potelés, testicules descendus"),
    ("F36-38", 250, 266, 340.0, 360.0, "À terme : naissance ≈ J266 (38 sem. p.f. = 40 sem. d'aménorrhée)"),
]
AXE_IDX = {a[0]: i for i, a in enumerate(AXE)}

# --------------------------------------------------------------------------------------
# Source A : reconstruction 3D (nos vidéos)
# --------------------------------------------------------------------------------------
DOSSIERS_A = {"CS13": "CS13.f4v", "CS14": "CS14_f4v", "CS15": "CS15_f4v", "CS16": "CS16_f4v",
              "CS17": "CS17_f4v", "CS19": "CS19_f4v", "CS20": "CS20_F4V"}
# Évaluation visuelle (embryons_3D/LISEZMOI.txt, 24/09/2026) — à mettre à jour quand on corrige
QUALITE_A = {
    "CS20": ("tres_bon", "SNC complet, cœur, foie, yeux, oreilles, corps vertébraux, membres, cordon"),
    "CS19": ("tres_bon", "très bon ; moelle basse manquante"),
    "CS17": ("bon", "cœur et foie corrects, cerveau approximatif"),
    "CS14": ("bon", "cœur et foie corrects, cerveau approximatif"),
    "CS15": ("moyen", "SNC déborde sur les somites, cerveau partiel"),
    "CS13": ("moyen", "SNC déborde sur les somites"),
    "CS16": ("bon", "labels corrigés le 24/09/2026 par une autre session (CS16_labels_corrige) : cœur, foie, ventricules et cerveau présents sur le rendu"),
}
NOTE_QUALITE = {"tres_bon": 4, "bon": 3, "moyen": 2, "faible": 1}
STRUCTURES_REF = ["enveloppe", "snc", "ventricules", "ganglions", "yeux", "cristallins", "vesicules_otiques",
                  "coeur", "cavite_pericardique", "foie", "tube_digestif", "squelette_axial_cartilage",
                  "chondrocrane", "cartilage_autre", "membres", "cordon_ombilical", "vaisseaux"]
STRUCTURES_PARTIELLES = {"tube_digestif", "vaisseaux"}  # partielles partout (LISEZMOI)
SEGMENTS_DIGESTIF = ["oesophage", "estomac", "duodenum", "intestin_moyen", "intestin_posterieur"]
# Hiérarchie du viewer : système → sous-groupe → structure (les noms sont ceux des maillages du pipeline)
UTILISER_B3D = False   # 24/09 : « reconstructions à partir de coupes (squelette, muscles) catastrophiques, oublier »
# … mais « les nouvelles prod on les affiche » : seules ces structures des stades issus des coupes vont dans la vue 3D
B3_VUE_AUTORISEES = ["axe_vertebral", "corps_vertebraux_predits"]
SYSTEMES = ["Système digestif", "Système nerveux", "Système vasculaire", "Os + cartilages", "Muscles",
            "Appendiculaire", "Derme (peau + enveloppes)", "Autres"]
STRUCT_SYS = {
    "tube_digestif": ["Système digestif", "Tube digestif (pipeline)"],
    "foie": ["Système digestif", "Glandes annexes"],
    "oesophage": ["Système digestif", "Tube digestif (reconstruction séparée)"],
    "estomac": ["Système digestif", "Tube digestif (reconstruction séparée)"],
    "duodenum": ["Système digestif", "Tube digestif (reconstruction séparée)"],
    "intestin_moyen": ["Système digestif", "Tube digestif (reconstruction séparée)"],
    "intestin_posterieur": ["Système digestif", "Tube digestif (reconstruction séparée)"],
    "snc": ["Système nerveux", "Système nerveux central"],
    "prosencephale": ["Système nerveux", "Vésicules cérébrales"], "mesencephale": ["Système nerveux", "Vésicules cérébrales"],
    "rhombencephale": ["Système nerveux", "Vésicules cérébrales"],
    "ventricule_prosencephale": ["Système nerveux", "Ventricules"], "ventricule_mesencephale": ["Système nerveux", "Ventricules"],
    "ventricule_rhombencephale": ["Système nerveux", "Ventricules"], "canal_central": ["Système nerveux", "Moelle et canal central"],
    "moelle": ["Système nerveux", "Moelle et canal central"], "moelle_rachidienne": ["Système nerveux", "Moelle et canal central"], "meninges_mesenchyme_cranien": ["Système nerveux", "Membranes"],
    "epiderme_cranien": ["Derme (peau + enveloppes)", "Épiderme crânien"],
    "ventricules": ["Système nerveux", "Système nerveux central"],
    "ganglions": ["Système nerveux", "Système nerveux périphérique"],
    "yeux": ["Système nerveux", "Organes des sens"],
    "cristallins": ["Système nerveux", "Organes des sens"],
    "vesicules_otiques": ["Système nerveux", "Organes des sens"],
    "coeur": ["Système vasculaire", "Cœur"],
    "cavite_pericardique": ["Système vasculaire", "Cœur"],
    "vaisseaux": ["Système vasculaire", "Vaisseaux"],
    "cavites_cardiaques": ["Système vasculaire", "Cœur"],
    "aorte_dorsale_gauche": ["Système vasculaire", "Aortes dorsales"],
    "aorte_dorsale_droite": ["Système vasculaire", "Aortes dorsales"],
    "aorte_commune": ["Système vasculaire", "Aortes dorsales"],
    "aorte_dorsale": ["Système vasculaire", "Aortes dorsales"],
    "somites_gauche": ["Muscles", "Somites et myotomes"],
    "somites_droite": ["Muscles", "Somites et myotomes"],
    "somites_myotomes": ["Muscles", "Somites et myotomes"],
    "muscles_membres": ["Appendiculaire", "Muscles des membres"],
    "muscles_dos": ["Muscles", "Muscles du dos"],
    "muscles_paroi": ["Muscles", "Muscles de la paroi"],
    "muscles_cou_tete": ["Muscles", "Muscles du cou et de la tête"],
    "muscles_membre_superieur": ["Appendiculaire", "Muscles des membres"],
    "muscles_membre_inferieur": ["Appendiculaire", "Muscles des membres"],
    "somites": ["Muscles", "Somites et myotomes"],
    "myotomes": ["Muscles", "Somites et myotomes"],
    "notochorde": ["Os + cartilages", "Notochorde"],
    "axe_vertebral": ["Os + cartilages", "Colonne : axe et prédiction"],
    "corps_vertebraux_predits": ["Os + cartilages", "Colonne : axe et prédiction"],
    "etages_vertebraux": ["Os + cartilages", "Colonne : axe et prédiction"], "somites_video": ["Os + cartilages", "Colonne : axe et prédiction"], "etages_manuels": ["Os + cartilages", "Colonne : axe et prédiction"], "etages_manuels_dorsal": ["Os + cartilages", "Colonne : axe et prédiction"], "reper1": ["Os + cartilages", "Colonne : axe et prédiction"], "reper2": ["Os + cartilages", "Colonne : axe et prédiction"],
    "corps_vertebraux": ["Os + cartilages", "Colonne : corps et arcs recrutés"], "arcs_neuraux": ["Os + cartilages", "Colonne : corps et arcs recrutés"],
    "cotes": ["Os + cartilages", "Côtes et sternum"], "sternum": ["Os + cartilages", "Côtes et sternum"],
    "ossifications": ["Os + cartilages", "Ossifications"],
    "membre_superieur": ["Appendiculaire", "Squelette des membres"], "membre_inferieur": ["Appendiculaire", "Squelette des membres"],
    "coeur_detoure": ["Système vasculaire", "Cœur"],
    "myocarde": ["Système vasculaire", "Cœur"],
    "digestif_oesophage": ["Système digestif", "Tube digestif (fusionné)"],
    "digestif_estomac": ["Système digestif", "Tube digestif (fusionné)"],
    "digestif_duodenum": ["Système digestif", "Tube digestif (fusionné)"],
    "digestif_intestin_moyen": ["Système digestif", "Tube digestif (fusionné)"],
    "digestif_intestin_posterieur": ["Système digestif", "Tube digestif (fusionné)"],
    "vaisseaux_aorte": ["Système vasculaire", "Aortes dorsales"],
    "vaisseaux_veines": ["Système vasculaire", "Veines"],
    "vaisseaux_cardinales": ["Système vasculaire", "Veines"],
    "vaisseaux_ombilicaux": ["Système vasculaire", "Veines"],
    "vaisseaux_vitellins": ["Système vasculaire", "Veines"],
    "vaisseaux_ombilicales": ["Système vasculaire", "Veines"],
    "vaisseaux_vitellines": ["Système vasculaire", "Veines"],
    "cardinale_commune_gauche": ["Système vasculaire", "Veines"],
    "cardinale_commune_droite": ["Système vasculaire", "Veines"],
    "ombilicale_gauche": ["Système vasculaire", "Veines"],
    "ombilicale_droite": ["Système vasculaire", "Veines"],
    "ombilicale_commune": ["Système vasculaire", "Veines"],
    "vitelline_gauche": ["Système vasculaire", "Veines"],
    "vitelline_droite": ["Système vasculaire", "Veines"],
    "cardinale_posterieure_gauche": ["Système vasculaire", "Veines"],
    "cardinale_posterieure_droite": ["Système vasculaire", "Veines"],
    "cardinale_anterieure_gauche": ["Système vasculaire", "Veines"],
    "cardinale_anterieure_droite": ["Système vasculaire", "Veines"],
    "veine_ombilicale_gauche": ["Système vasculaire", "Veines"],
    "veine_ombilicale_droite": ["Système vasculaire", "Veines"],
    "veine_vitelline_gauche": ["Système vasculaire", "Veines"],
    "veine_vitelline_droite": ["Système vasculaire", "Veines"],
    "squelette_axial_cartilage": ["Os + cartilages", "Squelette axial"],
    "sternum": ["Os + cartilages", "Squelette axial"],
    "membre_inferieur": ["Os + cartilages", "Membres"],
    "chondrocrane": ["Os + cartilages", "Chondrocrâne"],
    "cartilage_autre": ["Os + cartilages", "Autres cartilages"],
    "enveloppe": ["Derme (peau + enveloppes)", "Peau"],
    "membres": ["Appendiculaire", "Bourgeons des membres"],
    "cordon_ombilical": ["Derme (peau + enveloppes)", "Cordon ombilical"],
}


def taille_dossier(p):
    t = 0
    for r, _, fs in os.walk(p):
        for f in fs:
            try:
                t += os.path.getsize(os.path.join(r, f))
            except OSError:
                pass
    return t


def source_A():
    man = {}
    mm = os.path.join("embryons_3D", "master_manifest.json")
    if os.path.exists(mm):
        for s in json.load(open(mm, encoding="utf-8"))["stages"]:
            man[s["stage"]] = s
    res = {}
    for cs, d in DOSSIERS_A.items():
        out = os.path.join(d, "out")
        work = os.path.join(d, "work")
        videos = [f for f in os.listdir(d) if f.lower().endswith((".mp4", ".f4v"))] if os.path.isdir(d) else []
        fs = os.listdir(out) if os.path.isdir(out) else []
        s = man.get(cs)
        # le manifeste du stade (out/manifest.json) est la source la plus fraîche : les fusions d'autres sessions y arrivent d'abord
        ms = os.path.join(out, "manifest.json")
        master_retard = False
        if os.path.exists(ms):
            s_stage = json.load(open(ms, encoding="utf-8"))
            for x in s_stage.get("structures", []):
                x.setdefault("file", f"../{d}/out/{x.get('file', '')}") if False else None
            master_retard = os.path.exists(mm) and os.path.getmtime(ms) > os.path.getmtime(mm) + 60
            s = {**(s or {}), **s_stage, "structures": [dict(x, file=f"../{d}/out/{x['file']}") for x in s_stage["structures"]]}
        structs = []
        if s:
            for x in s["structures"]:
                structs.append({"nom": x["name"], "collection": x["collection"], "faces": x["faces"], "color": x.get("color", [0.7, 0.7, 0.7]),
                                "volume_mm3": round(x["volume_mm3"], 4), "confiance": x.get("confiance"),
                                "fichier_ok": os.path.exists(os.path.normpath(os.path.join("embryons_3D", x["file"])))})
        noms = {x["nom"] for x in structs}
        manquantes = [n for n in STRUCTURES_REF if n not in noms]
        q, qtxt = QUALITE_A.get(cs, ("inconnu", ""))
        etapes = {
            "extract": os.path.exists(os.path.join(work, "stack0.npy")),
            "fuse": os.path.exists(os.path.join(work, "dens_raw.npy")) or os.path.exists(os.path.join(work, "fused_meta.json")),
            "prep": os.path.exists(os.path.join(work, "dens.npy")),
            "axes": os.path.exists(os.path.join(work, "axes.json")),
            "segment": os.path.exists(os.path.join(work, "labels.npz")),
            "mesh": f"{cs}_all.glb" in fs,
        }
        # une étape aval faite prouve les étapes amont (CS20 a l'ancienne organisation de work/)
        cles = list(etapes)
        for i, k in enumerate(cles):
            if any(etapes[k2] for k2 in cles[i + 1:]):
                etapes[k] = True
        # reconstructions séparées, alignées sur le pipeline : out/<sous-dossier>/manifest_<nom>.json
        def lire_recon(nom, attendus, npz_nom, cle_confiance=None, dossier=None):
            dossier = dossier or nom
            md = os.path.join(out, dossier, f"manifest_{nom}.json")
            if not os.path.exists(md):
                return None
            m = json.load(open(md, encoding="utf-8"))
            import unicodedata
            norm = lambda t: unicodedata.normalize("NFKD", t).encode("ascii", "ignore").decode().lower()
            conf_txt = (m.get(cle_confiance) or "") if cle_confiance else ""
            conf_glob = conf_txt.split(" (")[0]
            # « aortes : bonne ; cardinales postérieures : moyenne » → confiance par famille de segments
            familles = []
            for part in conf_txt.split(";"):
                if ":" in part:
                    lib, cf = part.split(":", 1)
                    familles.append((lib.strip().lower(), cf.strip().split(" (")[0]))
            def conf_seg(nom):
                n = norm(nom).replace("_", " ")
                best, score = None, 0
                for lib, cf in familles:
                    toks = {w[:6] for w in re.findall(r"[a-z]+", norm(lib)) if len(w) > 3}
                    # les mots de côté / de regroupement pèsent peu : c'est le nom anatomique qui décide
                    sc = sum((0.25 if w in ("gauche", "droite", "commun", "cordon") else 1.0) for w in toks if w in n)
                    if sc > score:
                        best, score = cf, sc
                return best or conf_glob or "?"
            segs = []
            for x in m.get("segments", []):
                segs.append({"nom": x["name"], "confiance": x.get("confiance") or conf_seg(x["name"]), "faces": x.get("faces"),
                             "volume_mm3": round(x.get("volume_mm3", 0), 4), "color": x.get("color", [0.9, 0.5, 0.3]),
                             "url": f"../{d}/out/{dossier}/{x['file']}",
                             "fichier_ok": os.path.exists(os.path.join(out, dossier, x["file"]))})
            aligne = (s is not None and abs(m.get("mm_per_voxel", 0) - s["mm_per_voxel"]) < 1e-9
                      and all(abs(a - b) < 1e-6 for a, b in zip(m.get("center_voxel", []), s["center_voxel"])))
            presents = {x["nom"] for x in segs}
            manquants = [n for n in attendus if n not in presents]
            # absences déclarées non réalisables : champ « absents » (texte libre) ou famille « non traçable » dans la confiance
            na = {}
            abs_txt = m.get("absents")
            if isinstance(abs_txt, dict):
                for k, v in abs_txt.items():
                    na[k] = str(v)
            elif abs_txt:
                t = norm(str(abs_txt))
                for n_ in manquants:
                    mot = norm(n_.split("_")[-1])[:6]
                    if mot in t:
                        na[n_] = str(abs_txt)
            for lib, cf in familles:
                if "non tra" in norm(cf):
                    na[lib.split(" (")[0]] = cf
            return {"segments": segs, "aligne": aligne, "npz": os.path.exists(os.path.join(work, nom, npz_nom)),
                    "manquants": [n for n in manquants if n not in na], "na": na,
                    "confiance_globale": m.get(cle_confiance) if cle_confiance else None}
        dig = lire_recon("digestif", SEGMENTS_DIGESTIF, "digestif.npz")
        cardio = lire_recon("cardio", ["coeur", "cavites_cardiaques"], "coeur.npz", "confiance_coeur")
        vais = lire_recon("vaisseaux", [], "vaisseaux.npz", "confiance", dossier="cardio")
        musc = lire_recon("muscles", [], "muscles.npz")
        topo = lire_recon("topographie", [], "topographie.npz")
        if topo:
            mt = json.load(open(os.path.join(out, "topographie", "manifest_topographie.json"), encoding="utf-8"))
            ct = mt.get("controle") or {}
            fiab = ((ct.get("etages") or {}).get("fiabilite") or {})
            nv = (ct.get("etages") or {}).get("niveaux") or {}
            hors = [k for k, v in (nv.items() if isinstance(nv, dict) else []) if isinstance(v, dict) and v.get("hors_axe")]
            cal = sum(1 for v in (nv.values() if isinstance(nv, dict) else []) if isinstance(v, dict) and v.get("frontiere_calee"))
            sans_corps = [k for k, v in (nv.items() if isinstance(nv, dict) else []) if isinstance(v, dict) and v.get("corps") is False]
            prol = [k for k, v in (nv.items() if isinstance(nv, dict) else []) if isinstance(v, dict) and v.get("sur_prolongement")]
            topo["controle"] = {"axe_vertebral": ct.get("axe_vertebral"), "fiabilite": fiab, "zones": ct.get("zones"),
                                "niveaux_imposes": ({"composition": (ct.get("niveaux_imposes") or {}).get("composition"), "n": len(nv), "cales": cal,
                                                     "hors_axe": hors, "sans_corps": sans_corps, "sur_prolongement": prol,
                                                     "reperes": (fiab.get("reperes") or {}).get("source")} if fiab.get("mode") == "nombre impose" else None),
                                "avertissement": ct.get("avertissement"), "images": ct.get("images") or [],
                                "images_url": [("../" + os.path.join(out, "topographie", im)).replace("\\", "/") for im in (ct.get("images") or [])],
                                "video360": ct.get("video360"), "reperes_utilisateur": ct.get("reperes_utilisateur"), "etages_manuels": ct.get("etages_manuels"), "etages_manuels_dorsal": ct.get("etages_manuels_dorsal")}
            topo["fiable"] = bool(fiab.get("fiable")) and not ct.get("avertissement")
            for x in topo["segments"]:
                if x.get("confiance") == "utilisateur":                    # marques posées à la main par l'utilisateur : toujours en vue 3D
                    continue
                x["confiance"] = "bonne" if topo["fiable"] else "faible"
                x["fichier_ok"] = x["fichier_ok"] and topo["fiable"]      # non fiable → pas en vue 3D (mais lu dans le tableau)
        res[cs] = {
            "dossier": d, "videos": len(videos), "n_ply": len([f for f in fs if f.endswith(".ply")]),
            "digestif": dig, "cardio": cardio, "vais": vais, "musc": musc, "topo": topo,
            "glb": f"{cs}_all.glb" in fs, "blend": f"{cs}_embryon.blend" in fs,
            "nrrd_densite": f"{cs}_density.nrrd" in fs, "nrrd_labels": f"{cs}_labels.nrrd" in fs,
            "labels_corriges": any("corrige" in f for f in fs),
            "rendu_organes": os.path.exists(os.path.join("embryons_3D", "rendus", f"{cs}_organes.png")),
            "rendu_peau": os.path.exists(os.path.join("embryons_3D", "rendus", f"{cs}_peau.png")),
            "rendu_organes_url": f"rendus/{cs}_organes.png", "rendu_peau_url": f"rendus/{cs}_peau.png",
            "mm_per_voxel": round(s["mm_per_voxel"], 5) if s else None,
            "gl_mm": s.get("greatest_length_mm_assumed") if s else None,
            "structures": structs, "manquantes": manquantes,
            "qualite": q, "qualite_note": NOTE_QUALITE.get(q, 0), "qualite_txt": qtxt,
            "etapes": etapes, "taille_out_mo": round(taille_dossier(out) / 1e6) if os.path.isdir(out) else 0,
            "dans_morph": True, "master_retard": master_retard,
        }
    return res


# --------------------------------------------------------------------------------------
# Source B : coupes ehd.org
# --------------------------------------------------------------------------------------
BASE_B = "coupes embryos 9-23"


def source_B_3d():
    """Stades ou les coupes legendees (source B) ont produit un volume 3D (coupes embryos 9-23/
    embryo_stage_S/out/manifest.json, ecrit par squelette_3d.py). Rendu comme une entree de type A
    (meme format de manifeste) pour que la vue 3D et le tableau le lisent tels quels."""
    res = {}
    for md in glob.glob(os.path.join(BASE_B, "embryo_stage_*", "out", "manifest.json")):
        m = json.load(open(md, encoding="utf-8"))
        cs = m["stage"]
        d = os.path.dirname(os.path.dirname(md)).replace("\\", "/")
        out = os.path.join(d, "out")
        structs = [{"nom": x["name"], "collection": x.get("collection", "Os + cartilages"), "faces": x["faces"],
                    "color": x.get("color", [0.8, 0.8, 0.7]), "volume_mm3": round(x.get("volume_mm3", 0), 4),
                    "confiance": x.get("confiance", "?"), "n_masques": x.get("n_masques"),
                    "fichier_ok": os.path.exists(os.path.join(out, x["file"]))} for x in m["structures"]]
        q = m.get("qualite", "faible")
        def lire_recon_B(nom, attendus):
            md2 = os.path.join(out, nom, f"manifest_{nom}.json")
            if not os.path.exists(md2):
                return None
            m2 = json.load(open(md2, encoding="utf-8"))
            segs = [{"nom": x["name"], "confiance": x.get("confiance", "?"), "faces": x.get("faces"),
                     "volume_mm3": round(x.get("volume_mm3", 0), 4), "color": x.get("color", [0.9, 0.5, 0.3]),
                     "url": f"../{d}/out/{nom}/{x['file']}", "fichier_ok": os.path.exists(os.path.join(out, nom, x["file"])),
                     "repere": bool(x.get("repere"))} for x in m2.get("segments", [])]
            aligne = bool(m.get("reperage_seul")) or (abs(m2.get("mm_per_voxel", 0) - m["mm_per_voxel"]) < 1e-9
                      and all(abs(a - b) < 1e-6 for a, b in zip(m2.get("center_voxel", []), m["center_voxel"])))
            return {"segments": segs, "aligne": aligne, "npz": False,
                    "manquants": [n for n in attendus if n not in {x["nom"] for x in segs}], "na": {},
                    "confiance_globale": None, "source": m2.get("source", "coupes légendées ehd (source B)"),
                    "evaluation_visuelle": m2.get("evaluation_visuelle"), "confiance_auto": m2.get("confiance_auto"),
                    "certitude_moyenne": m2.get("certitude_moyenne"), "n_coupes": m2.get("sections"),
                    "controle": ({"bilan_zones": m2["controle"].get("bilan_zones"), "inventaire_incomplet": m2["controle"].get("inventaire_incomplet")}
                                 if isinstance(m2.get("controle"), dict) else None)}
        res[cs] = {
            "dossier": d, "origine": "B", "videos": 0, "n_ply": len(structs), "reperage_seul": bool(m.get("reperage_seul")),
            "echelle_source": m.get("echelle_source"),
            "digestif": lire_recon_B("digestif", SEGMENTS_DIGESTIF), "cardio": lire_recon_B("cardio", []),
            "vais": lire_recon_B("vaisseaux", []), "musc": lire_recon_B("muscles", []),
            "glb": os.path.exists(os.path.join(out, f"{cs}_all.glb")), "blend": False,
            "nrrd_densite": False, "nrrd_labels": False, "labels_corriges": False,
            "rendu_organes": False, "rendu_peau": False, "rendu_organes_url": None, "rendu_peau_url": None,
            "mm_per_voxel": round(m["mm_per_voxel"], 5), "gl_mm": m.get("greatest_length_mm_assumed"),
            "structures": structs, "manquantes": [],
            "qualite": q, "qualite_note": NOTE_QUALITE.get(q, 1),
            "qualite_txt": f"squelette seul, reconstruit depuis les {m.get('sections', '?')} coupes légendées ehd (source B) "
                           f"segmentées par MobileSAM ; "
                           + (f"échelle relevée : {m['echelle_source']}" if m.get("echelle_source")
                              else f"échelle supposée (GL {m.get('greatest_length_mm_assumed')} mm, pixel {m.get('px_mm_assumed')} mm)")
                           + " ; recalage par translation seulement",
            "controle_B": ({"rapport": f"../{d}/out/{m['controle'].get('rapport')}" if m["controle"].get("rapport") else None,
                            "n_masques": m["controle"].get("n_masques"), "n_rejets": m["controle"].get("n_rejets"),
                            "n_chaines": m["controle"].get("n_chaines"), "n_coupes_interpolees": m["controle"].get("n_coupes_interpolees"),
                            "bilan_zones": m["controle"].get("bilan_zones"), "inventaire_incomplet": m["controle"].get("inventaire_incomplet"),
                            "symetries": m["controle"].get("symetries"), "axe_vertebral": m["controle"].get("axe_vertebral"),
                            "etages": m["controle"].get("etages"),
                            "images": [f"../{d}/out/{im}" for im in dict.fromkeys(list(m["controle"].get("images") or []) + ["coupes_medianes.png", "axe_vertebral.png", "etages_vertebraux.png",
                                                                       f"{cs}_controle_colonne.png", f"{cs}_controle_axial.png", "symetrie_controle.png"])
                                       if os.path.exists(os.path.join(out, im))],
                            "rapports": [f"../{d}/out/{os.path.basename(f)}" for f in sorted(glob.glob(os.path.join(out, "controle_*.html")))]}
                           if isinstance(m.get("controle"), dict) else None),
            "etapes": {"extract": True, "fuse": False, "prep": False, "axes": False, "segment": True, "mesh": True},
            "taille_out_mo": round(taille_dossier(out) / 1e6) if os.path.isdir(out) else 0,
            "dans_morph": False,
        }
    return res


def source_B():
    res = {}
    for sd in glob.glob(os.path.join(BASE_B, "embryo_stage_*")):
        s = int(sd.rsplit("_", 1)[-1])
        cs = f"CS{s}"
        n_img = len(glob.glob(os.path.join(sd, "images_fond", "*.jpg")))
        n_gif = len(glob.glob(os.path.join(sd, "calques_gif", "*.gif")))
        n_txt = len(glob.glob(os.path.join(sd, "textes", "*.txt")))
        leg = os.path.join(sd, "legendes")
        f_struct = glob.glob(os.path.join(leg, "*_structures.json"))
        f_leg = glob.glob(os.path.join(leg, "*_legendes.json"))
        f_rap = glob.glob(os.path.join(leg, "*_rapport.txt"))
        n_struct = None
        n_points = 0
        top = []
        if f_struct:
            try:
                d = json.load(open(f_struct[0], encoding="utf-8"))
                n_struct = len(d)
                n_points = sum(len(v) for v in d.values())
                top = sorted(d.items(), key=lambda kv: -len(kv[1]))[:12]
                top = [{"nom": k, "points": len(v)} for k, v in top]
            except Exception as e:  # noqa
                n_struct = -1
        n_labels = n_traits = n_anom = 0
        if f_rap:
            txt = open(f_rap[0], encoding="utf-8", errors="replace").read()
            for m in re.finditer(r"(\d+) labels, (\d+) blocs, (\d+) traits apparies, (\d+) anomalies", txt):
                n_labels += int(m.group(1)); n_traits += int(m.group(3)); n_anom += int(m.group(4))
        n_overlays = len(glob.glob(os.path.join(leg, "overlays", "*.png")))
        # densité : ehd.org publie des dizaines/centaines de sections par stade ; ≤3 = scraping tronqué
        if n_img >= 50:
            dens = "dense"
        elif n_img >= 10:
            dens = "moyen"
        else:
            dens = "squelette"
        res[cs] = {
            "sections": n_img, "gifs": n_gif, "textes": n_txt, "legendes_json": bool(f_leg),
            "structures": n_struct, "points": n_points, "labels": n_labels, "traits": n_traits,
            "anomalies": n_anom, "overlays": n_overlays,
            "viewer": bool(glob.glob(os.path.join(leg, "*_viewer.html"))),
            "densite": dens, "top": top,
            "exemple_overlay": (sorted(glob.glob(os.path.join(leg, "overlays", "*.png")))[len(glob.glob(os.path.join(leg, "overlays", "*.png"))) // 2]
                                if n_overlays else None),
        }
        if res[cs]["exemple_overlay"]:
            res[cs]["overlay_url"] = "../" + res[cs]["exemple_overlay"].replace("\\", "/")
    return res


# --------------------------------------------------------------------------------------
# Source C : VOKA
# --------------------------------------------------------------------------------------
REF = os.path.join("embryo3d", "reference")


def source_C():
    res = {}
    vs = json.load(open(os.path.join(REF, "voka_stages.json"), encoding="utf-8"))
    props = {}
    pp = os.path.join(REF, "proportions.json")
    if os.path.exists(pp):
        props = json.load(open(pp, encoding="utf-8"))
    for st in vs["stades"]:
        J = st["id"]
        cap = os.path.join(REF, "captures", J)
        fs = os.listdir(cap) if os.path.isdir(cap) else []
        pngs = [f for f in fs if f.endswith(".png") and "_az" in f]
        p = props.get(J, {})
        lat = p.get("lateral_mm", {})
        res[J] = {
            "jour": st["jour"], "semaine": st["semaine_pf"], "carnegie": st["carnegie"], "plage": st["carnegie_plage"],
            "crl_min": st["crl_mm_min"], "crl_max": st["crl_mm_max"], "crl": st["crl_mm"],
            "reperes": st["reperes"], "nos_coupes": st["nos_coupes"], "page": st["page_voka"],
            "image": st["image_voka"], "image_ok": os.path.exists(os.path.join(REF, st["image_voka"])),
            "captures": len(pngs), "planche": "planche.png" in fs, "meta": "meta.json" in fs,
            "silhouettes": "silhouettes" in fs, "mesures": "mesures.json" in fs,
            "hull": any(f.endswith("_hull.ply") for f in fs),
            "image_url": f"../embryo3d/reference/{st['image_voka']}" if os.path.exists(os.path.join(REF, st["image_voka"])) else None,
            "mm_par_px": p.get("mm_par_px"), "ratio_h_l": lat.get("ratio_hauteur_largeur"),
            "largeur_mm": lat.get("largeur_bbox"), "hauteur_mm": lat.get("hauteur_bbox"),
        }
    return res


# --------------------------------------------------------------------------------------
# Fusion sur l'axe + reste à faire
# --------------------------------------------------------------------------------------
def fusion(A, B, C):
    voka_par_cs = {}
    for J, v in C.items():
        voka_par_cs.setdefault(v["carnegie"], []).append(J)
    # bornes contiguës pour la frise : les plages Carnegie se chevauchent, on coupe au milieu du chevauchement
    bornes = []
    for i, a in enumerate(AXE):
        deb = a[1] if i == 0 else (AXE[i - 1][2] + 1 + a[1]) / 2
        fin = a[2] + 1 if i == len(AXE) - 1 else (a[2] + 1 + AXE[i + 1][1]) / 2
        bornes.append((deb, max(fin, deb + 0.6)))
    # J VOKA dont la plage Carnegie couvre chaque stade (pour le hull et l'image de comparaison)
    def voka_couvrant(cs):
        n = int(cs[2:]) if cs.startswith("CS") else None
        for J, v in sorted(C.items(), key=lambda kv: kv[1]["jour"]):
            m = re.match(r"CS(\d+)(?:-(\d+))?$", v["plage"])
            if n is not None and m and int(m.group(1)) <= n <= int(m.group(2) or m.group(1)):
                return J
        return None
    B3 = source_B_3d()   # 3D depuis les coupes : gardée dans le tableau (repérage vertébral, topographie par zones) mais écartée de la vue 3D si UTILISER_B3D est faux
    lignes = []
    for i, (cs, j0, j1, c0, c1, rep) in enumerate(AXE):
        a = A.get(cs) or B3.get(cs); b = B.get(cs); cj = voka_par_cs.get(cs, [])
        if a and a.get("reperage_seul"):
            a["ecarte_3d"] = True; a["qualite"] = "faible"; a["qualite_note"] = 0
            a["qualite_txt"] = "repérage seul (axe vertébral, étages, zones) sans maillage 3D" + (f" ; échelle : {a['echelle_source']}" if a.get("echelle_source") else "") + (" ; tube digestif depuis les coupes (hors vue 3D)" if a.get("digestif") else "")
            a["b3_autorisees"] = []
            a["vue3d"] = {"glb": None, "J": voka_couvrant(cs), "hull": None, "voka_img": None}
        elif a and a.get("origine") == "B" and not UTILISER_B3D:
            a["ecarte_3d"] = True
            a["qualite"] = "faible"; a["qualite_note"] = 0
            a["qualite_txt"] = ("3D écartée de la vue (maillages depuis les coupes jugés inutilisables le 24/09) sauf les nouvelles productions : "
                                + ", ".join(B3_VUE_AUTORISEES) + " ; reconstructions séparées (digestif, muscles) gardées dans le tableau ; "
                                "le repérage vertébral et la topographie par zones sont la méthode retenue")
            noms = {x["nom"] for x in a["structures"]}
            a["b3_autorisees"] = [n for n in B3_VUE_AUTORISEES if n in noms]
            J = voka_couvrant(cs)
            a["vue3d"] = ({"glb": f"../{a['dossier']}/out/{cs}_all.glb", "J": J,
                           "hull": f"../embryo3d/reference/captures/{J}/{J}_hull.ply" if J and C[J]["hull"] else None,
                           "voka_img": f"../embryo3d/reference/{C[J]['image']}" if J and C[J]["image_ok"] else None}
                          if a["glb"] and a["b3_autorisees"] else {"glb": None, "J": J, "hull": None, "voka_img": None})
        elif a:
            J = voka_couvrant(cs)
            a["vue3d"] = {"glb": f"../{a['dossier']}/out/{cs}_all.glb" if a["glb"] else None, "J": J,
                          "hull": f"../embryo3d/reference/captures/{J}/{J}_hull.ply" if J and C[J]["hull"] else None,
                          "voka_img": f"../embryo3d/reference/{C[J]['image']}" if J and C[J]["image_ok"] else None}
        # score de couverture 0..3 (interne 3D / coupes 2D / enveloppe externe)
        interne3d = (a["qualite_note"] if a and not a.get("ecarte_3d") else 0)
        coupes = {"dense": 3, "moyen": 2, "squelette": 1}[b["densite"]] if b else 0
        env = 3 if any(C[J]["hull"] for J in cj) else (1 if cj else 0)
        if interne3d >= 3:
            statut = "3D interne"
        elif interne3d > 0:
            statut = "3D à corriger"
        elif a and a.get("ecarte_3d") and coupes >= 3:
            statut = "coupes denses + repérage par zones (3D écartée)"
        elif a and a.get("reperage_seul"):
            statut = "repérage vertébral seul (coupes tronquées)"
        elif coupes >= 3 and env:
            statut = "coupes denses + enveloppe (pas de 3D)"
        elif coupes >= 3:
            statut = "coupes denses (pas de 3D)"
        elif env >= 3:
            statut = "enveloppe VOKA seule"
        elif coupes == 2:
            statut = "coupes moyennes (pas de 3D)"
        elif coupes:
            statut = "coupes squelettiques"
        elif cs.startswith("F") and int(cs[1:3]) >= 14:
            statut = "rien (2e/3e trimestre : aucune source)"
        else:
            statut = "rien"
        lignes.append({"stade": cs, "jours": [j0, j1], "frise": list(bornes[i]), "crl": [c0, c1], "reperes": rep,
                       "A": a, "B": b, "C": cj, "score": {"interne3d": interne3d, "coupes": coupes, "enveloppe": env},
                       "statut": statut})
    return lignes


# structures qu'on VEUT détourer, même si aucun stade ne les a encore (objectif) — par système
OBJECTIF_SUPPL = {
    "Muscles": ["somites_myotomes", "muscles_membres"],
    "Système vasculaire": ["veines_cardinales_ant", "veines_cardinales_post", "veines_ombilicales", "veines_vitellines"],
}
NOTE_CONF = {"bonne": 3, "bon": 3, "tres_bon": 3, "moyenne": 2, "moyen": 2, "approximative": 1, "faible": 1}
# confiance déclarée par structure et par stade pour les labels du pipeline (sinon : qualité globale du stade)
# source : session « Reconstruction embryons CS13-CS20 (base) », message du 24/09/2026 08:00 (axial.py)
def conf_de(cs, st):
    """Confiance d'une structure du pipeline : le champ `confiance` du manifest (session base, meshexport.py) prime sur la table CONF_STRUCT."""
    c = (st or {}).get("confiance")
    if c in ("bonne", "moyenne", "faible"):
        return c
    return CONF_STRUCT.get((cs, (st or {}).get("nom")))


CONF_STRUCT = {
    ("CS13", "somites"): "faible", ("CS14", "somites"): "moyenne", ("CS15", "somites"): "bonne", ("CS16", "somites"): "moyenne",
    ("CS17", "squelette_axial_cartilage"): "bonne", ("CS19", "squelette_axial_cartilage"): "faible", ("CS20", "squelette_axial_cartilage"): "moyenne",
    # notochorde : faible partout (trop épaisse, probablement colonne périchordale / aorte) ; veines fusionnées : faible (session VHE)
    **{(cs, "notochorde"): "faible" for cs in ("CS13", "CS14", "CS15", "CS16", "CS17", "CS19", "CS20")},
    # SNC fin CS20 (session base, 24/09) : distance géodésique le long de l'axe neural
    **{("CS20", n): "bonne" for n in ("ventricule_prosencephale", "ventricule_mesencephale", "ventricule_rhombencephale", "canal_central",
                                       "prosencephale", "mesencephale", "rhombencephale", "moelle")},
    ("CS20", "meninges_mesenchyme_cranien"): "faible", ("CS20", "epiderme_cranien"): "moyenne",
    **{(cs, n): "faible" for cs in ("CS13", "CS14", "CS15") for n in ("vaisseaux_cardinales", "vaisseaux_ombilicaux", "vaisseaux_vitellins")},
}
# variantes de nommage regroupées sur une seule ligne de progression
ALIAS_PROG = {"aorte_dorsale_gauche": "aortes_dorsales", "aorte_dorsale_droite": "aortes_dorsales",
              "aorte_commune": "aortes_dorsales", "aorte_dorsale": "aortes_dorsales",
              "somites_gauche": "somites_myotomes", "somites_droite": "somites_myotomes", "somites": "somites_myotomes", "myotomes": "somites_myotomes",
              "coeur_detoure": "coeur", "vaisseaux_aorte": "aortes_dorsales", "vaisseaux_veines": "veines_cardinales_ant",
              "cardinale_posterieure_gauche": "veines_cardinales_post", "cardinale_posterieure_droite": "veines_cardinales_post",
              "vaisseaux_cardinales": "veines_cardinales_ant", "vaisseaux_ombilicales": "veines_ombilicales", "vaisseaux_vitellines": "veines_vitellines",
              "vaisseaux_ombilicaux": "veines_ombilicales", "vaisseaux_vitellins": "veines_vitellines",
              "cardinale_anterieure_gauche": "veines_cardinales_ant", "cardinale_anterieure_droite": "veines_cardinales_ant",
              "cardinale_commune_gauche": "veines_cardinales_ant", "cardinale_commune_droite": "veines_cardinales_ant",
              "ombilicale_gauche": "veines_ombilicales", "ombilicale_droite": "veines_ombilicales", "ombilicale_commune": "veines_ombilicales",
              "vitelline_gauche": "veines_vitellines", "vitelline_droite": "veines_vitellines",
              "veine_ombilicale_gauche": "veines_ombilicales", "veine_ombilicale_droite": "veines_ombilicales",
              "veine_vitelline_gauche": "veines_vitellines", "veine_vitelline_droite": "veines_vitellines",
              "digestif_oesophage": "oesophage", "digestif_estomac": "estomac", "digestif_duodenum": "duodenum",
              "digestif_intestin_moyen": "intestin_moyen", "digestif_intestin_posterieur": "intestin_posterieur"}


def norm_txt(t):
    import unicodedata
    return unicodedata.normalize("NFKD", str(t)).encode("ascii", "ignore").decode().lower()


def progression(A):
    """Par système : structures attendues × stades A, avec source (pipeline / reconstruction séparée) et confiance."""
    stades = sorted(A, key=lambda k: AXE_IDX[k])
    # structures observées
    vues = {}  # nom -> systeme
    al = lambda n: ALIAS_PROG.get(n, n)
    for a in A.values():
        for x in a["structures"]:
            vues.setdefault(al(x["nom"]), STRUCT_SYS.get(x["nom"], ["Autres"])[0])
        for rec in ("digestif", "cardio", "vais", "musc"):
            if a.get(rec):
                for x in a[rec]["segments"]:
                    vues.setdefault(al(x["nom"]), STRUCT_SYS.get(x["nom"], ["Autres"])[0])
    for sy, names in OBJECTIF_SUPPL.items():
        for n in names:
            vues.setdefault(n, sy)
    systemes = []
    tot_ok = tot_att = 0
    for sy in SYSTEMES:
        noms = [n for n, s2 in vues.items() if s2 == sy]
        if not noms:
            continue
        noms.sort(key=lambda n: (STRUCT_SYS.get(n, ["", "zz"])[1], n))
        lignes_ = []
        ok_sy = 0
        for n in noms:
            cells = {}
            for cs in stades:
                a = A[cs]
                c = None
                st = next((x for x in a["structures"] if al(x["nom"]) == n), None)
                if st and n not in STRUCTURES_PARTIELLES:
                    cf = conf_de(cs, st)
                    c = ({"src": "pipeline", "note": NOTE_CONF.get(cf, 2), "txt": cf + (" (manifest)" if st.get("confiance") in ("bonne", "moyenne", "faible") else " (déclarée)")} if cf
                         else {"src": "pipeline", "note": NOTE_CONF.get(a["qualite"], 2), "txt": QUALITE_A.get(cs, ("", ""))[0]})
                elif st:
                    c = {"src": "pipeline", "note": 1, "txt": "partiel"}
                for rec, lib in (("digestif", "digestif"), ("cardio", "cœur détouré"), ("vais", "aortes"), ("musc", "muscles auto")):
                    if a.get(rec):
                        seg = next((x for x in a[rec]["segments"] if al(x["nom"]) == n), None)
                        if seg:
                            conf = (seg.get("confiance") or "").split(" (")[0]
                            c2 = {"src": lib, "note": NOTE_CONF.get(conf, 2), "txt": conf}
                            if c is None or c2["note"] > c["note"]:
                                c = c2
                if c is None:
                    for rec in ("digestif", "cardio", "vais", "musc"):
                        for k, why in ((a.get(rec) or {}).get("na", {}) or {}).items():
                            kk = al(k.replace(" ", "_"))
                            kn = norm_txt(k)
                            if kk == n or (al(k) == n) or (n == "veines_cardinales_post" and "cardinale" in kn and "poster" in kn) or (n == "veines_cardinales_ant" and "cardinale" in kn and "anter" in kn):
                                c = {"na": True, "src": "déclaré non réalisable", "note": 0, "txt": why}
                if c and not c.get("na"):
                    ok_sy += 1
                cells[cs] = c
            sg = STRUCT_SYS.get(n, STRUCT_SYS.get(next((k for k, v in ALIAS_PROG.items() if v == n), n), ["", "Veines" if n.startswith("veines") else "—"]))[1]
            lignes_.append({"nom": n, "sous_groupe": sg, "cells": cells,
                            "n_ok": sum(1 for v in cells.values() if v and not v.get("na")),
                            "n_na": sum(1 for v in cells.values() if v and v.get("na"))})
        att = len(noms) * len(stades) - sum(x["n_na"] for x in lignes_)
        tot_ok += ok_sy; tot_att += att
        systemes.append({"systeme": sy, "structures": lignes_, "ok": ok_sy, "attendu": att,
                         "pct": round(100 * ok_sy / att) if att else 0})
    return {"stades": stades, "systemes": systemes, "ok": tot_ok, "attendu": tot_att,
            "pct": round(100 * tot_ok / tot_att) if tot_att else 0}


def reste_a_faire(A, B, C, lignes):
    T = []

    def add(prio, source, stade, action, detail="", cmd=""):
        T.append({"prio": prio, "source": source, "stade": stade, "action": action, "detail": detail, "cmd": cmd})

    # --- A : reconstruction 3D
    for cs, a in sorted(A.items(), key=lambda kv: AXE_IDX[kv[0]]):
        if a["qualite"] == "faible":
            add("haute", "A", cs, "Reprendre la segmentation (qualité faible)", a["qualite_txt"],
                f"python embryo3d/fix_vent.py {a['dossier']}  puis  python embryo3d/pipeline.py {a['dossier']} --steps segment,mesh")
        elif a["qualite"] == "moyen":
            add("moyenne", "A", cs, "Corriger le SNC qui déborde (qualité moyenne)", a["qualite_txt"],
                f"3D Slicer : {a['dossier']}/out/{cs}_density.nrrd + {cs}_labels.nrrd → réinjecter dans work/labels.npz → --steps mesh")
        manq = [m for m in a["manquantes"]]
        if a.get("cardio") and manq:
            manq = [m for m in manq if m not in ("coeur",)]
        if manq:
            add("moyenne" if len(manq) >= 3 else "basse", "A", cs, "Structures non détectées : " + ", ".join(manq),
                "présentes à CS20 (référence) mais absentes de ce stade", "")
        for e, ok in a["etapes"].items():
            if not ok:
                add("haute", "A", cs, f"Étape pipeline manquante : {e}", "", f"python embryo3d/pipeline.py {a['dossier']} --steps {e}")
        for x in a["structures"]:
            if not x["fichier_ok"]:
                add("haute", "A", cs, f"PLY manquant : {x['nom']}", "référencé par master_manifest.json", "bash embryo3d/build_all.sh")
        if not a["rendu_organes"] or not a["rendu_peau"]:
            add("basse", "A", cs, "Rendu de contrôle manquant", "", "bash embryo3d/build_all.sh")
    add("haute", "A", "CS13", "Orientation différente de CS20 dans la vue 3D (dos vers la droite sous l'angle de référence)",
        "l'étape axes est fragile sur les embryons très recourbés : vérifier check_canonical*.png, corriger avec flip_axis.py, puis --steps mesh ; "
        "en attendant, décalages az/el/roll de CS13 dans la trajectoire de la vue 3D", "python embryo3d/check_volume.py CS13.f4v")
    add("moyenne", "A", "tous", "Vaisseaux partiels partout",
        "lumières tubulaires seulement ; envisager un suivi de tube (tracking) sur les coupes ou une correction Slicer", "")
    n_dig = sum(1 for a in A.values() if a.get("digestif"))
    if n_dig:
        add("moyenne", "A", "tous", f"Fusionner le tube digestif séparé ({n_dig} stades) dans labels.npz",
            "work/digestif/digestif.npz (packbits, même format) → label tube_digestif propre, puis --steps mesh et build_all.sh", "")
        for cs, a in sorted(A.items(), key=lambda kv: AXE_IDX[kv[0]]):
            dg = a.get("digestif")
            if not dg:
                continue
            approx = [x["nom"] for x in dg["segments"] if x["confiance"] in ("approximative", "faible")]
            if approx:
                add("basse", "A", cs, "Segments digestifs approximatifs : " + ", ".join(approx),
                    "anse herniée difficile à suivre ; points modifiables dans embryo3d/digestif_points/" + cs + ".json",
                    f"python embryo3d/digestif_build.py {a['dossier']}")
            if dg["manquants"]:
                add("basse", "A", cs, "Segments digestifs absents : " + ", ".join(dg["manquants"]), "", "")
            for k, why in (dg.get("na") or {}).items():
                add("basse", "A", cs, f"Déclaré non réalisable : {k}", why, "")
            if not dg["aligne"]:
                add("haute", "A", cs, "Tube digestif non aligné sur le pipeline (mm/voxel ou centre différents)", "", "")
    else:
        add("moyenne", "A", "tous", "Tube digestif partiel partout", "reconstruction séparée à lancer (digestif_build.py)", "")
    n_vai = sum(1 for a in A.values() if a.get("vais"))
    if n_vai:
        for cs, a in sorted(A.items(), key=lambda kv: AXE_IDX[kv[0]]):
            va = a.get("vais")
            if not va:
                continue
            conf = (va.get("confiance_globale") or "").split(" (")[0]
            if conf.startswith("faible"):
                add("moyenne", "A", cs, "Aortes dorsales tracées avec confiance faible", va.get("confiance_globale", ""), "")
            if not va["aligne"]:
                add("haute", "A", cs, "Aortes dorsales non alignées sur le pipeline", "", "")
        add("moyenne", "A", "tous", "Veines (cardinales, ombilicales, vitellines) non tracées",
            "seules les aortes dorsales sont faites (vaisseaux_points/<CS>.json) ; même méthode par points", "")
    n_topo = [cs for cs, a in A.items() if a.get("topo")]
    if n_topo:
        nf = [cs for cs in n_topo if not A[cs]["topo"]["fiable"]]
        if nf:
            add("haute", "A", ", ".join(sorted(nf, key=lambda k: AXE_IDX[k])), "Topographie par région : axe vertébral non fiable (labels colonne trop courts)",
                "piste proposée par la session Extraction squelette : guider l'axe par le tube neural (label snc, complet partout) — en attente de votre décision", "")
        fi = [cs for cs in n_topo if A[cs]["topo"]["fiable"]]
        if fi:
            st = ", ".join(sorted(fi, key=lambda k: AXE_IDX[k]))
            add("haute", "A", st, "Colonne : caler l'axe et les étages modélisés dans l'embryon (adaptations restantes)",
                "l'axe vertébral et les anneaux d'étages doivent rentrer dans le volume du stade (position, échelle, courbure) avant toute segmentation", "")
            add("haute", "A", st, "Colonne : recruter les volumes vertébraux dans le volume du stade",
                "une fois l'axe calé : segmenter corps et arcs dans work/labels.npz autour de chaque étage (guide = axe + étages + notochorde) → PLY par vertèbre, fusion", "")
    add("haute", "A", "tous", "Morphing web : exporter les 7 clés par structure (topologie CS20 projetée, lissée) pour le viewer 3dht",
        "embryons_3D/morph_web/<structure>.npz {positions (7,N,3), faces} + morph_web.json ; session base ; puis interpolation linéaire par le curseur temps dans viewer_3dh (agrégateur)", "")
    add("moyenne", "A", "tous", "Vidéo du morphing (mp4 1080p, caméra fixe, 7 stades) pour la version mobile de 3dht", "rendu Blender de la scène maître validée (planche 16:25)", "")
    n_car = sum(1 for a in A.values() if a.get("cardio"))
    if n_car:
        add("moyenne", "A", "tous", f"Fusionner le cœur détouré ({n_car} stades) dans labels.npz",
            "work/cardio/coeur.npz (coeur_plein, myocarde, cavites_cardiaques) → labels coeur propres, puis --steps mesh et build_all.sh", "")
        for cs, a in sorted(A.items(), key=lambda kv: AXE_IDX[kv[0]]):
            ca = a.get("cardio")
            if not ca:
                continue
            conf = (ca.get("confiance_globale") or "").split(" (")[0]
            if conf.startswith("faible"):
                add("moyenne", "A", cs, "Cœur détouré avec confiance faible", ca.get("confiance_globale", ""),
                    f"python embryo3d/cardio_contours.py {a['dossier']}")
            elif conf.startswith("moyenne"):
                add("basse", "A", cs, "Cœur détouré avec confiance moyenne", ca.get("confiance_globale", ""), "")
            if not ca["aligne"]:
                add("haute", "A", cs, "Cœur détouré non aligné sur le pipeline", "", "")
    add("haute", "A", "tous", "SNC : segmenter finement (vésicules cérébrales, ventricules, membranes / os du crâne, moelle précise) — snc et ventricules automatiques ne sont plus affichés",
        "consigne utilisateur 24/09 : rien du système nerveux automatique tant que ce n'est pas segmenté proprement (« ventricules » jugé inutilisable)", "")
    add("moyenne", "A", "CS19", "Moelle basse manquante", "LISEZMOI : moelle épinière basse non retenue à CS19", "")
    n_mus = sum(1 for a in A.values() if a.get("musc"))
    if not n_mus:
        add("moyenne", "A", "tous", "Système « Muscles » vide : aucune structure musculaire segmentée",
            "somites / myotomes (CS13-17) puis masses musculaires des membres (CS19-20)", "python embryo3d/muscles_build.py <dossier>")
    else:
        add("moyenne", "A", "tous", f"Somites / muscles : première passe automatique sur {n_mus} stades, à vérifier et découper",
            "colonnes paraxiales non découpées en blocs ; planches out/muscles/<CS>_controle.png ; affiner (seuils --t-somite/--band) ou repasser par points", "")
        for cs, a in sorted(A.items(), key=lambda kv: AXE_IDX[kv[0]]):
            for x in a["musc"]["segments"] if a.get("musc") else []:
                if x["nom"].startswith("somites") and x.get("confiance") == "faible":
                    add("basse", "A", cs, f"{x['nom']} : confiance faible", "", "")
    add("basse", "A", "CS16", "Vérifier CS16_labels_corrige.json : cavite_pericardique porte le label 1, comme ventricules",
        "fusion voulue ou erreur de table ? le PLY cavite_pericardique existe pourtant dans out/", "")
    add("basse", "A", "CS18", "Aucune vidéo de coupes : le morphing saute CS17→CS19",
        "interpoler (shape keys à 0.5) ou chercher une source CS18 (ehd n'a que 2 sections)", "")
    add("haute", "A", "CS21-23", "Aucune reconstruction interne après CS20",
        "candidat : CS23 = 209 coupes ehd (images_fond) + hull VOKA J56 comme enveloppe → nouveau volume", "")
    retard = [cs for cs, a in A.items() if a.get("master_retard")]
    if retard:
        add("haute", "A", ", ".join(sorted(retard, key=lambda k: AXE_IDX[k])), "Scène maître et GLB en retard sur les manifestes de stade",
            "des structures ont été ajoutées ou refaites au niveau du stade (fusions, somites…) : relancer l'assemblage", "bash embryo3d/build_all.sh")
    add("haute", "A", "tous", "Nouveau viewer « embryon en croissance par système » sur harcelon.fr/3dh (nouveau thread)",
        "point de départ : embryo3d/viewer_3dh.py + embryons_3D/site_3dh/PASSATION_3DH.md ; à tester, corriger, publier dans site/3dh/", "python embryo3d/viewer_3dh.py --local")
    add("moyenne", "A", "tous", "Rendu final sur la 2e machine (i9 9900 + RTX 3080)",
        "copier embryons_3D/ + CSxx_*/out/ (chemins relatifs) ; Blender 5.2", "")

    # --- B : coupes ehd
    for cs, b in sorted(B.items(), key=lambda kv: AXE_IDX[kv[0]]):
        if b["densite"] == "squelette":
            add("haute" if AXE_IDX[cs] >= AXE_IDX["CS13"] else "basse", "B", cs,
                f"Scraping tronqué : {b['sections']} section(s) seulement", "relancer le scraper ehd.org pour ce stade puis extraire_legendes.py",
                f'python "{BASE_B}/extraire_legendes.py" {cs[2:]}  &&  python "{BASE_B}/agreger_structures.py" {cs[2:]}')
        if b["structures"] == -1:
            add("haute", "B", cs, "stage_S_structures.json illisible", "encodage ?", "")
        if not b["legendes_json"]:
            add("haute", "B", cs, "Légendes non extraites", "", f'python "{BASE_B}/extraire_legendes.py" {cs[2:]}')
        if b["anomalies"] and b["sections"] >= 10:
            add("basse", "B", cs, f"{b['anomalies']} anomalies OCR/traits à relire", "voir overlays/ et stage_S_rapport.txt", "")
    add("haute", "B→A", "CS13", "Recaler le nuage de points ehd (247 sections, 363 structures) sur le volume 3D CS13",
        "seul stade où les deux sources sont denses : sert à valider / nommer nos labels automatiques", "")
    add("moyenne", "B", "CS23", "Squelette : compléter par symétrie les pièces légendées d'un seul côté",
        "côtes 1-2 et 7-10, tibia/fibula/talus/calcanéum, patella, Reichert, maxillaire/palatin", 'python "coupes embryos 9-23/squelette_3d.py" 23')
    add("moyenne", "B", "CS23", "Squelette : attribuer un niveau aux vertèbres sans niveau et recréer L1-L5 (jamais légendées)",
        "au pas vertébral mesuré (0,64 mm)", 'python "coupes embryos 9-23/controle_squelette.py" 23')
    add("moyenne", "B", "CS23", "Segmenter directement l'autre membre inférieur (tibia, fibula, talus, calcanéum, patella) et les fragments droits radius/ulna/fémur, non légendés",
        "les membres sont exclus de la symétrie (pose propre à chaque membre) ; script de symétrie : coupes embryos 9-23/symetrie.py", "")
    add("moyenne", "B", "CS23", "Faire valider par l'expert attendu_squelette.json / attendu_zones.json",
        "longueurs hypothétiques par zone (tête, cou, thorax, abdomen, pelvis, membres) ; principe d'extraction prédictive par zone à généraliser aux organes et muscles", "")
    add("moyenne", "B→A", "CS23", "Exploiter les 209 coupes ehd CS23 (1309 structures) pour un volume 3D",
        "il faut un recalage section à section (pas de vidéo) ; enveloppe = hull VOKA J56", "")
    add("basse", "B", "CS2-CS12", "Stades pré-somitiques : hors périmètre 3D actuel", "garder les nuages de points comme documentation", "")

    # --- C : VOKA
    for J, c in C.items():
        if c["captures"] < 24:
            add("moyenne", "C", J, f"Captures incomplètes ({c['captures']}/24)", "", f"python embryo3d/reference/capture_voka.py --stages {J}")
        if not c["hull"]:
            add("moyenne", "C", J, "Visual hull manquant", "", "python embryo3d/reference/silhouettes.py")
    add("moyenne", "C→A", "CS13,15,17,19,20", "Comparer côte à côte VOKA / notre GLB (emplacements sous la peau)",
        "viewer.html : J28↔CS13, J35↔CS15, J42↔CS17, J49↔CS19/20 ; non tracé comme fait", "python embryo3d/reference/serve.py")
    add("moyenne", "C→A", "tous", "Appliquer le matériau VOKA_peau au blend master",
        "blender_skin_material.py sur la collection Enveloppe de embryons_CS13-CS20_morph.blend", "")
    add("haute", "—", "F14-F38", "2e et 3e trimestres : aucune source (ni coupes, ni vidéos, ni modèle externe) jusqu'à la naissance",
        "l'objectif est de couvrir toute la grossesse : chercher une source fœtale (échographies 3D, IRM fœtale, atlas) ou extrapoler par mise à l'échelle CRL", "")
    add("haute", "C", "F9-F13", "Fœtus 9-13 sem. : enveloppe VOKA seule, aucun intérieur",
        "décider : arrêter la chaîne à CS23, ou extrapoler l'intérieur CS20/23 dans les hulls fœtaux", "")
    add("basse", "C", "J42", "voka_stages.json cite encore « CS18_f4v (mal nommé) »", "le dossier s'appelle maintenant CS17_f4v", "")
    add("basse", "C", "J28", "Le modèle VOKA « J28 » a des palettes de main (≈CS15-16)", "se fier aux repères anatomiques, pas au jour", "")

    ordre = {"haute": 0, "moyenne": 1, "basse": 2}
    T.sort(key=lambda t: (ordre[t["prio"]], t["source"], AXE_IDX.get(t["stade"], 99)))
    etat = {}
    fe = os.path.join("embryons_3D", "taches_etat.json")
    if os.path.exists(fe):
        try:
            etat = json.load(open(fe, encoding="utf-8"))
        except Exception:
            etat = {}
    for t in T:
        t["id"] = hashlib.md5(f"{t['source']}|{t['stade']}|{t['action']}".encode("utf-8")).hexdigest()[:8]
        e = etat.get(t["id"], {})
        t["instance"] = e.get("instance", "")
        t["statut"] = e.get("statut", "a_faire")
        t["maj"] = e.get("maj", "")
    return T


# --------------------------------------------------------------------------------------
# HTML
# --------------------------------------------------------------------------------------
HTML = r"""<!DOCTYPE html>
<html lang="fr">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
__ROBOTS__
<title>Embryon 3D · agrégateur</title>
<style>
:root{--bg:#f6f4ef;--panel:#fff;--ink:#1d1d1b;--mut:#6b6a66;--line:#e2ded4;
 --a:#c6552e;--b:#2f6f8f;--c:#5f8f3e;--q4:#3f8f4f;--q3:#8fb04a;--q2:#e2b13c;--q1:#d9583b;--q0:#d8d4cb;
 --haute:#d9583b;--moy:#e2b13c;--basse:#9aa3a8}
@media (prefers-color-scheme: dark){:root:not([data-theme=light]){--bg:#17181a;--panel:#212326;--ink:#ecebe6;--mut:#a5a39c;--line:#34373b;--q0:#3a3d42}}
:root[data-theme=dark]{--bg:#17181a;--panel:#212326;--ink:#ecebe6;--mut:#a5a39c;--line:#34373b;--q0:#3a3d42}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.45 system-ui,Segoe UI,Roboto,sans-serif}
.wrap{max-width:1280px;margin:0 auto;padding:20px 16px 60px}
h1{font-size:24px;margin:0 0 4px}h2{font-size:18px;margin:34px 0 10px;border-bottom:2px solid var(--line);padding-bottom:4px}
.sub{color:var(--mut);margin:0 0 18px}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:12px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:14px 16px;border-top:4px solid var(--line)}
.card.A{border-top-color:var(--a)}.card.B{border-top-color:var(--b)}.card.C{border-top-color:var(--c)}
#cards .card{cursor:pointer;position:relative;user-select:none;padding-bottom:26px}#cards .card::after{content:'actif · cliquer pour masquer';position:absolute;right:12px;bottom:8px;font-size:11px;color:var(--mut)}
body.off-A .card.A,body.off-B .card.B,body.off-C .card.C{opacity:.42}body.off-A #cards .card.A::after,body.off-B #cards .card.B::after,body.off-C #cards .card.C::after{content:'masqué · cliquer pour réafficher';color:var(--q1)}
body.off-A [data-src="A"],body.off-B [data-src="B"],body.off-C [data-src="C"],body.off-A [data-lane="A"],body.off-B [data-lane="B"],body.off-C [data-lane="C"]{display:none}
body.off-C #vokaLbl,body.off-C #vokaImg,body.off-C #hullLbl{display:none}
.card h3{margin:0 0 6px;font-size:16px}.card .big{font-size:28px;font-weight:700;line-height:1.1}
.card p{margin:4px 0;color:var(--mut);font-size:13.5px}
.kv{display:grid;grid-template-columns:auto 1fr;gap:2px 10px;font-size:13.5px;margin-top:8px}.kv b{font-weight:600}
svg{width:100%;height:auto;display:block}
.tl{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:10px 12px;overflow-x:auto}
.leg{display:flex;flex-wrap:wrap;gap:14px;font-size:13px;color:var(--mut);margin:8px 0 0}
.leg span i{display:inline-block;width:12px;height:12px;border-radius:3px;vertical-align:-1px;margin-right:5px}
table{width:100%;border-collapse:collapse;background:var(--panel);border:1px solid var(--line);border-radius:10px;overflow:hidden;font-size:13.5px}
th,td{padding:7px 9px;border-bottom:1px solid var(--line);vertical-align:top;text-align:left}
th{background:var(--bg);font-weight:600;position:sticky;top:0;z-index:1}
tr.st{cursor:pointer}tr.st:hover td{background:rgba(127,127,127,.06)}
tr.det td{background:rgba(127,127,127,.04);padding:10px 14px}
.pill{display:inline-block;padding:1px 8px;border-radius:999px;font-size:12px;font-weight:600;color:#fff;white-space:nowrap}
.q4{background:var(--q4)}.q3{background:var(--q3)}.q2{background:var(--q2);color:#1d1d1b}.q1{background:var(--q1)}.q0{background:var(--q0);color:var(--mut)}
.src{font-weight:700}.src.A{color:var(--a)}.src.B{color:var(--b)}.src.C{color:var(--c)}
.mut{color:var(--mut)}.small{font-size:12.5px}
.det-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:14px}
.det-grid h4{margin:0 0 6px;font-size:14px}.det-grid ul{margin:0;padding-left:18px}.det-grid li{margin:1px 0}
.thumbs img{max-height:150px;max-width:100%;border-radius:6px;border:1px solid var(--line);background:#fff;margin:4px 6px 0 0}
.todo{display:grid;grid-template-columns:auto auto auto 1fr 180px auto;gap:6px 12px;align-items:start;background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:8px 12px}
.todo .st{font-size:12px;font-weight:600;padding:1px 7px;border-radius:999px;display:inline-block;color:#fff}.st.a_faire{background:var(--q0);color:var(--mut)}.st.en_cours{background:var(--q2);color:#1d1d1b}.st.termine{background:var(--q4)}.st.echec{background:var(--q1)}
.todo .inst{font-size:12.5px}.todo .act button{font-size:12px;padding:3px 8px;border:1px solid var(--line);border-radius:6px;background:var(--bg);color:var(--ink);cursor:pointer;margin:1px 2px 1px 0}.todo .act button:disabled{opacity:.45;cursor:not-allowed}
#lanceurEtat{font-size:12.5px}#lanceurEtat.on{color:var(--c)}#lanceurEtat.off{color:var(--mut)}
.layers{max-height:300px;overflow:auto;border:1px solid var(--line);border-radius:6px;padding:4px 6px;font-size:12.5px}
.layers details{margin:2px 0}.layers summary{list-style:none;cursor:pointer;display:flex;align-items:center;gap:5px;font-weight:700;font-size:12.5px}.layers summary::before{content:'▸';color:var(--mut);font-size:11px}.layers details[open]>summary::before{content:'▾'}
.layers .sg{margin:2px 0 0 14px}.layers .sg>label{font-weight:600;font-size:12px;margin-left:0}.layers label{display:flex;align-items:center;gap:5px;margin:1px 0 1px 26px;font-size:12px}.layers .vide{margin-left:14px;color:var(--mut);font-size:11.5px}
.layers i{display:inline-block;width:11px;height:11px;border-radius:3px;border:1px solid rgba(0,0,0,.25)}
.v3d .buf{color:var(--q2);font-weight:600}
.todo>div{padding:6px 0;border-bottom:1px solid var(--line)}
.todo .p{font-weight:700;text-transform:uppercase;font-size:11.5px;letter-spacing:.04em}
.p.haute{color:var(--haute)}.p.moyenne{color:var(--moy)}.p.basse{color:var(--basse)}
code{font-family:ui-monospace,Consolas,monospace;font-size:12px;background:rgba(127,127,127,.12);padding:1px 5px;border-radius:4px;word-break:break-all}
.filters{display:flex;flex-wrap:wrap;gap:8px;margin:0 0 10px;align-items:center}
.filters label{font-size:13px;display:inline-flex;gap:4px;align-items:center;cursor:pointer}
input[type=search]{padding:6px 10px;border:1px solid var(--line);border-radius:8px;background:var(--panel);color:var(--ink);min-width:220px}
.bar{height:8px;border-radius:4px;background:var(--q0);overflow:hidden}.bar i{display:block;height:100%}
#prog td.c{text-align:center;padding:4px 6px;font-size:11.5px;white-space:nowrap}#prog td.c span{display:inline-block;min-width:54px;padding:2px 6px;border-radius:5px;color:#fff;font-weight:600}
#prog .cna{background:repeating-linear-gradient(45deg,var(--q0),var(--q0) 4px,transparent 4px,transparent 8px);color:var(--mut)!important}#prog .c3{background:var(--q4)}#prog .c2{background:var(--q2);color:#1d1d1b!important}#prog .c1{background:var(--q1)}#prog .c0{background:var(--q0);color:var(--mut)!important}
#prog tr.sys td{background:rgba(127,127,127,.09);font-weight:700}#prog tr.sys .bar{width:160px;display:inline-block;vertical-align:middle;margin-left:10px}
#prog td.nom{white-space:nowrap}#prog td.nom .mut{font-weight:400}
.v3d{display:grid;grid-template-columns:1fr 250px;gap:12px;background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:10px}
.v3d canvas{width:100%;height:540px;display:block;border-radius:8px;background:#2b2d31;cursor:grab}
.v3d .side{font-size:13px}.v3d .side img{width:100%;border-radius:6px;background:#000;margin:4px 0}
.v3d .side h4{margin:8px 0 4px;font-size:13px}.v3d label{display:block;margin:3px 0}
.v3d button{padding:5px 10px;border:1px solid var(--line);border-radius:7px;background:var(--bg);color:var(--ink);cursor:pointer;margin:2px 2px 2px 0}
.v3d button.on{background:var(--a);color:#fff;border-color:var(--a)}.v3d button.big{font-weight:700;padding:6px 12px;border-width:2px}.v3d button.big.on{background:var(--b);border-color:var(--b)}
#dayBtns button{padding:2px 5px;font-size:11px;min-width:34px;border-radius:5px;line-height:1.2}#dayBtns button b{display:block;font-size:9.5px;font-weight:600;opacity:.8}
#dayBtns button.q4{background:var(--q4);color:#fff;border-color:transparent}#dayBtns button.q3{background:var(--q3);color:#fff;border-color:transparent}#dayBtns button.q2{background:var(--q2);color:#1d1d1b;border-color:transparent}#dayBtns button.q1{background:var(--q1);color:#fff;border-color:transparent}#dayBtns button.q0{background:var(--q0);color:var(--mut)}
#dayBtns button.on{outline:2px solid var(--a);outline-offset:1px}
.v3d input[type=range]{width:100%}
.rng{position:relative;height:34px}.rng input[type=range]{position:absolute;left:0;top:10px;width:100%;margin:0;-webkit-appearance:none;appearance:none;background:transparent;pointer-events:none;height:16px}
.rng input[type=range]::-webkit-slider-runnable-track{height:6px;background:transparent}.rng #day::-webkit-slider-runnable-track{background:var(--q0);border-radius:3px}
.rng input[type=range]::-webkit-slider-thumb{-webkit-appearance:none;pointer-events:auto;width:16px;height:16px;border-radius:50%;background:var(--ink);border:2px solid #fff;margin-top:-5px;cursor:grab}
.rng input.bar::-webkit-slider-thumb{width:8px;height:22px;border-radius:3px;background:var(--q1);border:1px solid #fff;margin-top:-8px;cursor:ew-resize}
.rng .band{position:absolute;top:13px;height:6px;background:var(--b);opacity:.45;border-radius:3px;pointer-events:none}
#sysMenu,#dayBtns,#stageBtns,#stageTicks,#shown,#playRecon,.v3d .tl,.v3d .tl + .leg,#v3d-msg,.v3d .stage > div[style*="margin-top:10px"],#hullLbl,#hullInfo,#vokaLbl,#vokaImg,#skinT,#wire,#digLbl,#cardioLbl,#vaisLbl,#muscLbl,.v3d details.traj,.v3d .side > div:first-child + .small.mut,#stageBtns + .small.mut,#dayBtns + .small.mut,.v3d .side h4{display:none!important}
.v3d .side label:has(#skinT),.v3d .side label:has(#wire){display:none!important}
.v3d .side .small.mut{display:none}.v3d .side .layers{display:block!important;max-height:60vh}
.v3d .stage > .small.mut{display:none}
.topbar{display:flex;align-items:baseline;justify-content:space-between;gap:12px;margin:0 0 10px}.topbar h1{margin:0;font-size:22px}.topbar a{color:var(--mut)}
#reste{margin-top:48px}#reste .sub{margin-top:6px}
.status{margin-top:8px;border:1px solid var(--line);border-left:5px solid var(--a);border-radius:8px;padding:8px 12px;background:var(--panel)}
.status .l1{font-size:16px;font-weight:700}.status .l2{font-size:13px;margin-top:2px}.status .pbar{height:6px;background:var(--q0);border-radius:3px;margin-top:6px;overflow:hidden}.status .pbar i{display:block;height:100%;background:var(--a);width:0}.v3d .mode{font-weight:700}.v3d .mode.auto{color:var(--c)}.v3d .mode.manuel{color:var(--a)}
.v3d textarea{width:100%;font:11.5px ui-monospace,Consolas,monospace;height:80px;background:var(--bg);color:var(--ink);border:1px solid var(--line);border-radius:6px}
#v3d-msg{position:absolute;left:12px;top:12px;color:#fff;background:rgba(0,0,0,.55);padding:4px 8px;border-radius:6px;font-size:12.5px;pointer-events:none}
.v3d .stage{position:relative}
.sysmenu{display:flex;flex-direction:column;gap:4px;margin:0 0 10px;font-size:13px}
.sysmenu .sys{display:flex;flex-direction:row;align-items:center;gap:8px;border:0;border-radius:8px;padding:1px 4px;background:transparent;min-height:28px}
.sysmenu .sys>button{font-weight:700;font-size:12.5px;padding:3px 10px;width:190px;border:1px solid var(--line);border-radius:8px;background:var(--panel);color:var(--ink);cursor:pointer;text-align:left;flex:0 0 auto}
.sysmenu .sub{display:flex;flex-wrap:wrap;gap:4px;flex:1 1 auto}.sysmenu .sub button{font-size:12px;padding:2px 9px;border:1px solid var(--line);border-radius:999px;background:var(--panel);color:var(--ink);cursor:pointer}
.sysmenu button.off{opacity:.4;text-decoration:line-through}.sysmenu button.mix{border-style:dashed}
.sysmenu .hint{color:var(--mut);font-size:12px;padding:0 6px}
.playrecon{margin-top:8px;border:2px solid var(--b);border-radius:8px;padding:6px 8px;background:var(--bg)}
.playrecon .t{font-weight:700;font-size:12.5px;margin-bottom:4px}.playrecon .fam{display:flex;flex-wrap:wrap;gap:4px 10px;align-items:center;font-size:12px}
.playrecon label{display:inline-flex;align-items:center;gap:4px;border:1px solid var(--line);border-radius:999px;padding:2px 8px;background:var(--panel);cursor:pointer}
.playrecon label.absent{opacity:.45}.playrecon label b{font-weight:600}.playrecon label .iv{color:var(--mut);font-size:11px}
.playrecon .act{display:flex;flex-wrap:wrap;gap:6px;margin-top:6px}.playrecon .act button{font-size:12px}
#prog td.c span{cursor:pointer}#prog td.nom{cursor:pointer}#prog td.nom:hover,#prog td.c span:hover{outline:2px solid var(--a);outline-offset:1px}
.v3d .stage.blink{animation:blink .9s 2}@keyframes blink{0%{box-shadow:0 0 0 0 var(--a)}50%{box-shadow:0 0 0 4px var(--a)}100%{box-shadow:0 0 0 0 var(--a)}}
@media (max-width:900px){.v3d{grid-template-columns:1fr}.v3d canvas{height:380px}}
@media (max-width:720px){.todo{grid-template-columns:auto 1fr}.todo>div:nth-child(6n+2),.todo>div:nth-child(6n+3),.todo>div:nth-child(6n+5),.todo>div:nth-child(6n+6){display:none}
 th:nth-child(n+6),td:nth-child(n+6){display:none}}
/* ===================== design ===================== */
:root{--font:'Inter',system-ui,Segoe UI,Roboto,sans-serif;--r:12px;--sh:0 1px 2px rgba(0,0,0,.06),0 8px 24px -12px rgba(0,0,0,.18);--dark:#15171a;--dark2:#1e2125;--dark-line:#2d3137;--dark-ink:#e8e6e1;--dark-mut:#9a9890}
body{font-family:var(--font);font-size:14.5px;letter-spacing:-.005em}
.wrap{max-width:1360px;padding:24px 20px 80px}
h1{font-size:26px;letter-spacing:-.02em;font-weight:800}h1 small{font-weight:500;color:var(--mut);font-size:14px;margin-left:8px}
h2{font-size:15px;letter-spacing:.04em;text-transform:uppercase;color:var(--mut);border:0;padding:0;margin:40px 0 12px;display:flex;align-items:baseline;gap:10px}
h2::before{content:'';display:inline-block;width:22px;height:3px;background:var(--a);border-radius:2px;transform:translateY(-3px)}
h2 .small{text-transform:none;letter-spacing:0;font-weight:400}
.card{border-radius:var(--r);box-shadow:var(--sh);border-color:transparent}
/* bande compacte des sources */
#cards{grid-template-columns:repeat(3,1fr);gap:10px}
#cards .card{padding:10px 14px 30px;border-top-width:5px}#cards .card h3{font-size:14px;margin-bottom:2px}#cards .card .big{font-size:22px}#cards .card p{font-size:12.5px;margin:2px 0}
#cards .card .kv{display:none;font-size:12.5px}#cards .card.open .kv{display:grid}
#cards .card .more{position:absolute;left:14px;bottom:8px;font-size:11px;color:var(--mut);cursor:pointer;text-decoration:underline dotted}
/* vue 3D = pièce maîtresse, cadre sombre */
.v3d{background:var(--dark);color:var(--dark-ink);border:0;border-radius:var(--r);box-shadow:var(--sh);padding:14px;grid-template-columns:1fr 270px}
.v3d .mut,.v3d .small.mut{color:var(--dark-mut)}.v3d .side h4{color:var(--dark-mut);text-transform:uppercase;letter-spacing:.06em;font-size:11px;margin:12px 0 6px}
.v3d canvas{border-radius:10px;background:#0f1113;height:72vh;min-height:420px}
.v3d .stage{display:flex;flex-direction:column}.v3d canvas{order:-3}.v3d #v3d-msg{order:-3}.v3d #status{order:-2}.v3d #sysMenu{order:6;margin-top:10px}
.v3d button{background:var(--dark2);color:var(--dark-ink);border-color:var(--dark-line)}.v3d button:hover{border-color:var(--dark-mut)}
.v3d button.on{background:var(--a);border-color:var(--a);color:#fff}.v3d button.big.on{background:var(--b);border-color:var(--b)}
.sysmenu .sys{background:transparent}.sysmenu .sys>button,.sysmenu .sub button{background:var(--dark);border-color:var(--dark-line);color:var(--dark-ink)}
.sysmenu .hint{color:var(--dark-mut)}.sysmenu .hint button{background:var(--dark2);color:var(--dark-ink);border-color:var(--dark-line)}
.status{background:var(--dark2);border-color:var(--dark-line);border-left-color:var(--a)}.status .l2{color:var(--dark-mut)}.status .pbar{background:var(--dark-line)}
.playrecon{background:var(--dark2);border-color:var(--b)}.playrecon label{background:var(--dark);border-color:var(--dark-line);color:var(--dark-ink)}.playrecon .t{color:var(--dark-ink)}
#dayBtns button.q0{background:var(--dark2);color:var(--dark-mut);border-color:var(--dark-line)}#stageBtns button{background:var(--dark2)}
.layers{background:var(--dark2);border-color:var(--dark-line)}.layers summary::before{color:var(--dark-mut)}.layers .vide{color:var(--dark-mut)}
.v3d .tl{background:var(--dark2);border-color:var(--dark-line)}.v3d .tl svg text{fill:var(--dark-ink)}.v3d .leg{color:var(--dark-mut)}
.v3d textarea{background:var(--dark);color:var(--dark-ink);border-color:var(--dark-line)}.v3d details summary{color:var(--dark-mut)}
.v3d .side img{background:#000;border:1px solid var(--dark-line)}
.rng #day::-webkit-slider-runnable-track{background:var(--dark-line)}.rng input[type=range]::-webkit-slider-thumb{background:#fff;border-color:var(--dark)}
#v3d-msg{background:rgba(0,0,0,.7);font-weight:600}
/* tableaux */
table{border-radius:var(--r);box-shadow:var(--sh);border-color:transparent}th{background:var(--panel);border-bottom:2px solid var(--line)}
/* tâches */
.todo{border-radius:var(--r);box-shadow:var(--sh);border-color:transparent}
@media (max-width:900px){#cards{grid-template-columns:1fr}.v3d{grid-template-columns:1fr;padding:8px}.v3d .stage{display:flex;flex-direction:column}.v3d canvas{order:-2;height:52vh;min-height:260px}.v3d #status{order:-1}.v3d #sysMenu{order:5}.wrap{padding:10px 10px 60px}.topbar{flex-direction:column;gap:2px}}
</style>
<link rel="preconnect" href="https://fonts.googleapis.com"><link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&display=swap" rel="stylesheet">
</head>
<body><div class="wrap">
<!-- ouverture : le viewer en lecture, le reste plus bas -->


<div class="topbar"><h1>Embryon 3D <small>croissance par système · agrégateur des sources</small></h1><a href="#reste" class="small">sources, progression, tâches ▾</a></div>
<div class="v3d-embed"><iframe id="viewer3dh" src="__VIEWER3DH__" title="Embryon humain · croissance par système" loading="eager" allow="fullscreen" allowfullscreen></iframe></div>
<style>.v3d-embed{border-radius:12px;overflow:hidden;box-shadow:var(--sh);background:#0e1013}.v3d-embed iframe{width:100%;height:86vh;min-height:520px;border:0;display:block}
.v3d-embed.fs{position:fixed;inset:0;z-index:1000;border-radius:0;height:100dvh}.v3d-embed.fs iframe{height:100%!important;min-height:0}
@media (max-width:900px){.v3d-embed{border-radius:8px}.v3d-embed iframe{height:calc(100dvh - 48px);min-height:360px}}</style>
<script>(function(){const f=document.getElementById('viewer3dh');function taille(){if(!window.matchMedia('(max-width:900px)').matches){f.style.height='';return;}const top=f.getBoundingClientRect().top+window.scrollY;f.style.height=Math.max(360,window.innerHeight-top-6)+'px';}
taille();window.addEventListener('resize',taille);window.addEventListener('load',taille);setTimeout(taille,500);})();
window.addEventListener('message',e=>{if(e.data&&e.data.v3d==='fullscreen'){document.querySelector('.v3d-embed').classList.toggle('fs',!!e.data.on);document.documentElement.style.overflow=e.data.on?'hidden':'';}});</script>
<div id="reste"><h2>Sources <span class="small mut">(cliquer un cadre pour masquer / réafficher sa source partout)</span></h2>

<p class="sub">Généré le <b id="date"></b> par <code>embryo3d/agregateur.py</code> — axe commun : jours post-fécondation · stade Carnegie · CRL (O'Rahilly &amp; Müller 2010, Moore &amp; Persaud), de la fécondation à la naissance (J266). Cliquer un cadre ci-dessous pour masquer / réafficher sa source partout. Relancer le script pour rafraîchir.</p>

<div class="cards" id="cards"></div>
<h2>Progression du détourage par système <span class="mut small" id="progTot"></span></h2>
<div class="small mut">Ce qui est détouré (pipeline ou reconstruction séparée) contre ce qui ne l'est pas, pour chaque structure attendue et chaque stade reconstruit. Vert = bonne confiance, jaune = moyenne, rouge = faible ou partiel, gris = absent, hachuré « n/r » = déclaré non réalisable dans ce volume (hors dénominateur). Les structures « objectif » (muscles, veines) n'existent encore nulle part.</div>
<div id="progCards" class="cards" style="margin:10px 0"></div>
<table id="prog"><thead></thead><tbody></tbody></table>

<h2>Contenu par stade <span class="mut small">(cliquer une ligne pour le détail)</span></h2>
<div class="filters"><input type="search" id="q" placeholder="filtrer : stade, structure, statut…">
<label><input type="checkbox" id="onlyGaps"> seulement les stades sans 3D interne</label></div>
<table id="tab"><thead><tr><th>Stade</th><th>Jours</th><th>CRL mm</th><th>Statut</th><th data-src="A"><span class="src A">A</span> ehd vidéos → 3D</th><th data-src="B"><span class="src B">B</span> ehd coupes</th><th data-src="C"><span class="src C">C</span> VOKA (enveloppe)</th></tr></thead><tbody></tbody></table>

<h2>Reste à faire <span class="mut small" id="ntodo"></span> <span id="lanceurEtat" class="off"></span></h2>
<div class="small mut">Colonne « instance » : qui fait le travail (session Claude nommée, lanceur local, personne). « Lancer » exécute la commande de la tâche sur ce PC via <code>python embryo3d/lanceur.py</code> (port 8791) ; « Confier » attribue la tâche à une instance et copie un prompt prêt à coller dans cette session Claude. L'état est partagé dans <code>embryons_3D/taches_etat.json</code>.</div>
<div id="enCours" class="cards" style="margin:8px 0 12px"></div>
<div class="filters" id="todoF"><label><input type="checkbox" value="haute" checked> haute</label><label><input type="checkbox" value="moyenne" checked> moyenne</label><label><input type="checkbox" value="basse" checked> basse</label>
<span class="mut small">· source :</span><label><input type="checkbox" value="A" checked> A</label><label><input type="checkbox" value="B" checked> B</label><label><input type="checkbox" value="C" checked> C</label>
<span class="mut small">· statut :</span><label><input type="checkbox" value="s:a_faire" checked> à faire</label><label><input type="checkbox" value="s:en_cours" checked> en cours</label><label><input type="checkbox" value="s:termine"> terminées</label><label><input type="checkbox" value="s:echec" checked> échec</label></div>
<div class="todo" id="todo"></div>

<h2>Comment lire les trois sources</h2>
<div class="cards">
<div class="card A"><h3>A · ehd.org, vidéos → reconstruction 3D</h3><p>ehd.org fournit deux sources : les vidéos de coupes (A) et les coupes légendées (B). A = volumes fusionnés à partir des 3 piles vidéo par stade, segmentés automatiquement (pipeline.py). Précision : voxel 6-25 µm, structures internes nommées, PLY/GLB/NRRD/Blend. Échelle = plus grande longueur typique du stade (CRL_MM dans meshexport.py). Intervalle : 7 stades ponctuels CS13→CS20, sans CS18.</p></div>
<div class="card B"><h3>B · ehd.org, coupes légendées</h3><p>Deuxième source ehd.org : coupes histologiques légendées (Virtual Human Embryo). Précision : 2D par section, points (section, x, y) par structure nommée, pas d'échelle mm ni de recalage entre sections. Intervalle : 22 stades CS2→CS23, mais seuls 8 stades ont &gt; 10 sections (scraping tronqué ailleurs).</p></div>
<div class="card C"><h3>C · VOKA</h3><p>10 modèles externes J28→J91, hebdomadaires, jusqu'à la fin du 1er trimestre (la grossesse va jusqu'à J266). Précision : enveloppe seulement, captures 24 vues → silhouettes → visual hull orthographique (proportions, placement), échelle absolue = CRL littérature. Les enveloppes (visual hull) se superposent au modèle dans la vue 3D.</p></div>
</div>
</div>
<script>window.__errs=[];addEventListener('error',e=>{window.__errs.push((e.error&&e.error.stack)||e.message)});</script>
<script id="data" type="application/json">__DATA__</script>
<script>
const D=JSON.parse(document.getElementById('data').textContent);
document.getElementById('date').textContent=D.date;
// rafraîchissement automatique : la page se recharge quand une nouvelle génération est publiée
if(D.version){setInterval(async()=>{try{const r=await fetch('version.txt?t='+Date.now(),{cache:'no-store'});if(r.ok){const v=(await r.text()).trim();if(v&&v!==D.version)location.reload();}}catch(e){}},20000);}
const QC={4:'q4',3:'q3',2:'q2',1:'q1',0:'q0'};const QT={tres_bon:'très bon',bon:'bon',moyen:'moyen',faible:'faible'};
const fmt=n=>n==null?'—':(+n).toLocaleString('fr-FR',{maximumFractionDigits:3});
// ---- cartes
const A=D.A,B=D.B,C=D.C;
const nA=Object.keys(A).length,nStruct=Object.values(A).reduce((s,a)=>s+a.structures.length,0);
const nB=Object.keys(B).length,nSec=Object.values(B).reduce((s,b)=>s+b.sections,0),nPts=Object.values(B).reduce((s,b)=>s+b.points,0),nBd=Object.values(B).filter(b=>b.densite!='squelette').length;
const nC=Object.keys(C).length,nCap=Object.values(C).reduce((s,c)=>s+c.captures,0),nHull=Object.values(C).filter(c=>c.hull).length;
const jA=D.lignes.filter(l=>l.A).map(l=>l.jours),jB=D.lignes.filter(l=>l.B).map(l=>l.jours),jC=Object.values(C).map(c=>c.jour);
document.getElementById('cards').innerHTML=`
<div class="card A"><h3><span class="src A">A</span> ehd.org · vidéos → reconstruction 3D</h3><div class="big">${nA} stades</div>
<p>Vidéos de coupes ehd.org (3D Atlas), fusionnées et segmentées par nous · CS13→CS20 (sans CS18) · jours ${jA[0][0]}→${jA[jA.length-1][1]} · ${nStruct} maillages</p>
<div class="kv"><b>Très bon</b><span>${Object.keys(A).filter(k=>A[k].qualite=='tres_bon').join(', ')||'—'}</span><b>Bon</b><span>${Object.keys(A).filter(k=>A[k].qualite=='bon').join(', ')||'—'}</span><b>Moyen</b><span>${Object.keys(A).filter(k=>A[k].qualite=='moyen').join(', ')||'—'}</span><b>Faible</b><span>${Object.keys(A).filter(k=>A[k].qualite=='faible').join(', ')||'—'}</span><b>Livrable</b><span>embryons_CS13-CS20_morph.blend (slider stage 0..6)</span><b>Tube digestif séparé</b><span>${Object.values(A).filter(a=>a.digestif).length} stades, ${Object.values(A).reduce((n,a)=>n+(a.digestif?a.digestif.segments.length:0),0)} segments (digestif_CS13-CS20.blend)</span><b>Cœur détouré</b><span>${Object.values(A).filter(a=>a.cardio).length} stades (coeur + cavités cardiaques, out/cardio/)</span><b>Aortes dorsales</b><span>${Object.values(A).filter(a=>a.vais).length} stades, ${Object.values(A).reduce((n,a)=>n+(a.vais?a.vais.segments.length:0),0)} segments (tracé par points)</span><b>Somites / muscles (auto)</b><span>${Object.values(A).filter(a=>a.musc).length} stades (out/muscles/, première passe)</span></div></div>
<div class="card B"><h3><span class="src B">B</span> ehd.org · coupes légendées</h3><div class="big">${nB} stades</div>
<p>Coupes histologiques légendées (Virtual Human Embryo) · CS2→CS23 · jours ${jB[0][0]}→${jB[jB.length-1][1]} · ${nSec} sections · ${fmt(nPts)} points légendés</p>
<div class="kv"><b>Denses</b><span>${Object.keys(B).filter(k=>B[k].densite=='dense').sort((a,b)=>+a.slice(2)-+b.slice(2)).join(', ')}</span><b>Moyens</b><span>${Object.keys(B).filter(k=>B[k].densite=='moyen').sort((a,b)=>+a.slice(2)-+b.slice(2)).join(', ')||'—'}</span><b>Tronqués (≤3 sections)</b><span>${nB-nBd} stades</span><b>Légendes extraites</b><span>${Object.values(B).filter(b=>b.legendes_json).length}/${nB}</span></div></div>
<div class="card C"><h3><span class="src C">C</span> VOKA — enveloppe externe</h3><div class="big">${nC} modèles</div>
<p>J28→J91 (hebdo) · CS13→fœtus 13 sem. · ${nCap} captures · ${nHull} visual hulls</p>
<div class="kv"><b>Recouvre nos stades</b><span>J28/CS13, J35/CS15, J42/CS17, J49/CS19-20</span><b>Prolonge après nos coupes</b><span>J56/CS23 puis J63→J91 (fœtus 9-13 sem.), s'arrête à la fin du 1er trimestre</span><b>Reste à couvrir</b><span>semaines 14→38 (naissance) : aucune source pour l'instant</span><b>Usage</b><span>enveloppe de référence : forme, proportions, placement des organes sous la peau, aspect de la peau</span></div></div>`;
// ---- cadres = interrupteurs
const srcOff=new Set();window.srcOff=srcOff;function toggleSrc(k){if(srcOff.has(k))srcOff.delete(k);else srcOff.add(k);document.body.classList.toggle('off-'+k,srcOff.has(k));
 const cb=document.querySelector('#todoF input[value="'+k+'"]');if(cb){cb.checked=!srcOff.has(k);}if(typeof rendTodo==='function')rendTodo();if(window.updateShown)setTimeout(window.updateShown,50);}
['A','B','C'].forEach(k=>{const c=document.querySelector('#cards .card.'+k);if(!c)return;const m=document.createElement('span');m.className='more';m.textContent='détails ▾';m.onclick=e=>{e.stopPropagation();c.classList.toggle('open');m.textContent=c.classList.contains('open')?'réduire ▴':'détails ▾';};c.appendChild(m);c.onclick=()=>toggleSrc(k);});
// ---- timeline
(function(){
 const W=1180,L=90,R=20,x0=L,x1=W-R,y0=34,rowH=36,rows=['A · ehd vidéos→3D','B · ehd coupes','C · VOKA','Bilan'];
 const J1=96,JMAX=267,F1=0.62;   // semaines 1-13 sur 62 % de la largeur, semaines 14-38 sur le reste
 const X=j=>j<=J1?x0+(x1-x0)*F1*j/J1:x0+(x1-x0)*(F1+(1-F1)*(j-J1)/(JMAX-J1));
 const H=y0+rows.length*rowH;
 let s=`<svg viewBox="0 0 ${W} ${H+34}" xmlns="http://www.w3.org/2000/svg" font-family="system-ui,sans-serif" font-size="11">`;
 for(let j=0;j<=JMAX;j+=7){const w=j/7;const show=j<=J1||w%2===0;if(!show)continue;s+=`<line x1="${X(j)}" y1="${y0-6}" x2="${X(j)}" y2="${H}" stroke="var(--line)"/><text x="${X(j)}" y="${y0-10}" text-anchor="middle" fill="var(--mut)" font-size="${j<=J1?10:9}">J${j}</text>`;}
 for(let w=1;w<=38;w++){const a=X((w-1)*7),b=X(Math.min(w*7,JMAX));if(w<=13||w%2===0)s+=`<text x="${(a+b)/2}" y="${H+16}" text-anchor="middle" fill="var(--mut)" font-size="${w<=13?10:9}">${w<=13?'sem. '+w:w}</text>`;}
 s+=`<line x1="${X(J1)}" y1="${y0-24}" x2="${X(J1)}" y2="${H}" stroke="var(--mut)" stroke-dasharray="3 3"/><text x="${X(J1)+4}" y="${y0-26}" fill="var(--mut)" font-size="10">fin du 1er trimestre · échelle compressée →</text>`;
 s+=`<line x1="${X(266)}" y1="${y0-24}" x2="${X(266)}" y2="${H}" stroke="var(--q1)" stroke-dasharray="3 3"/><text x="${X(266)-4}" y="${y0-26}" text-anchor="end" fill="var(--q1)" font-size="10">naissance J266</text>`;
 rows.forEach((r,i)=>{s+=`<text x="${L-8}" y="${y0+i*rowH+rowH/2+4}" text-anchor="end" font-weight="600" fill="var(--ink)" data-lane="${'ABC'[i]||''}">${r}</text>`;});
 const col=q=>`var(--q${q})`;const lane={A:'',B:'',C:'',S:''};
 D.lignes.forEach(l=>{const xa=X(l.frise[0]),xb=Math.max(X(l.frise[1]),xa+3);const sc=l.score;
  const lab=(y,q,t,title)=>`<rect x="${xa}" y="${y+5}" width="${xb-xa-1}" height="${rowH-10}" rx="3" fill="${col(q)}"><title>${title}</title></rect>`+((t&&(xb-xa)>=14)?`<text x="${(xa+xb)/2}" y="${y+rowH/2+4}" text-anchor="middle" fill="${q>=3||q==1?'#fff':'var(--ink)'}" font-size="${(xb-xa)>30?10:8.5}">${t}</text>`:'');
  const court=l.stade.replace('CS','').replace(/^F(\d+)-\d+$/,'$1');
  lane.A+=lab(y0,l.A?sc.interne3d:0,l.A?((xb-xa)>30?l.stade:court):'',l.A?`${l.stade} : ${QT[l.A.qualite]} — ${l.A.structures.length} structures`:`${l.stade} : pas de reconstruction`);
  lane.B+=lab(y0+rowH,l.B?sc.coupes*1.33|0:0,l.B?l.B.sections+' s.':'',l.B?`${l.stade} : ${l.B.sections} sections, ${l.B.structures} structures`:`${l.stade} : pas de coupes`);
  const bq=sc.interne3d>=3?4:sc.interne3d>0?2:(sc.coupes>=3?3:sc.enveloppe>=3?2:sc.coupes?1:0);
  lane.S+=lab(y0+3*rowH,bq,(xb-xa)>30?l.stade:court,`${l.stade} : ${l.statut}`);
 });
 Object.values(C).forEach(c=>{const x=X(c.jour),y=y0+2*rowH+rowH/2;lane.C+=`<circle cx="${x}" cy="${y}" r="9" fill="${c.hull?'var(--q4)':'var(--q2)'}"><title>J${c.jour} · ${c.plage} · CRL ${c.crl_min}-${c.crl_max} mm · ${c.captures} captures${c.hull?' · hull':''}</title></circle><text x="${x}" y="${y+16}" text-anchor="middle" font-size="9.5" fill="var(--mut)">J${c.jour}</text>`;});
 s+=`<g data-lane="A">${lane.A}</g><g data-lane="B">${lane.B}</g><g data-lane="C">${lane.C}</g><g>${lane.S}</g><line id="tlCursor" x1="${X(28)}" y1="${y0-8}" x2="${X(28)}" y2="${H}" stroke="var(--a)" stroke-width="2"/></svg>`;const tlEl=document.getElementById('timeline');if(tlEl)tlEl.innerHTML=s;
 window.tlX=X;window.setDayCursor=d=>{const c=document.getElementById('tlCursor');if(c){c.setAttribute('x1',X(d));c.setAttribute('x2',X(d));}};
})();
// ---- progression du détourage
(function(){const P=D.progression;if(!P)return;const col=n=>n>=3?'c3':n==2?'c2':n==1?'c1':'c0';const lib={3:'bon',2:'moyen',1:'faible',0:'—'};
 document.getElementById('progTot').textContent=`(${P.ok} / ${P.attendu} cases structure × stade détourées, ${P.pct} %)`;
 document.getElementById('progCards').innerHTML=P.systemes.map(sy=>{const q=sy.pct>=75?'q4':sy.pct>=40?'q2':sy.pct>0?'q1':'q0';return `<div class="card"><h3>${sy.systeme}</h3><div class="big">${sy.pct} %</div><p>${sy.ok} / ${sy.attendu} · ${sy.structures.length} structure${sy.structures.length>1?'s':''} × ${P.stades.length} stades</p><div class="bar"><i class="${q}" style="width:${sy.pct}%;background:var(--${q})"></i></div><p class="small">manque : ${sy.structures.filter(x=>x.n_ok+x.n_na<P.stades.length).map(x=>x.nom+(x.n_ok?' ('+(P.stades.length-x.n_ok-x.n_na)+' st.)':'')).join(', ')||'rien'}${sy.structures.some(x=>x.n_na)?' · <span class="mut">n/r = déclaré non réalisable : '+sy.structures.filter(x=>x.n_na).map(x=>x.nom+' ('+x.n_na+')').join(', ')+'</span>':''}</p></div>`;}).join('');
 document.querySelector('#prog thead').innerHTML='<tr><th>Système · structure</th>'+P.stades.map(cs=>`<th>${cs}</th>`).join('')+'<th>détouré</th></tr>';
 document.querySelector('#prog tbody').innerHTML=P.systemes.map(sy=>`<tr class="sys"><td>${sy.systeme}<div class="bar"><i style="width:${sy.pct}%;background:var(--${sy.pct>=75?'q4':sy.pct>=40?'q2':'q1'})"></i></div></td>`+P.stades.map(cs=>{const n=sy.structures.filter(x=>x.cells[cs]).length;return `<td class="c">${n}/${sy.structures.length}</td>`;}).join('')+`<td class="c">${sy.ok}/${sy.attendu}</td></tr>`+sy.structures.map(x=>`<tr><td class="nom" data-nom="${x.nom}" title="cliquer : lire cette structure sur ses stades">▶ ${x.nom} <span class="mut small">· ${x.sous_groupe}</span></td>`+P.stades.map(cs=>{const c=x.cells[cs];return `<td class="c"><span class="${c&&c.na?'cna':col(c?c.note:0)}" data-stage="${cs}" data-nom="${x.nom}" title="${c?c.src+' · '+c.txt+' · cliquer : voir dans la vue 3D':'non détouré'}">${c?(c.na?'n/r':(c.src==='pipeline'?lib[c.note]:c.src)):'—'}</span></td>`;}).join('')+`<td class="c">${x.n_ok}/${P.stades.length-x.n_na}${x.n_na?' <span class="mut">(+'+x.n_na+' n/r)</span>':''}</td></tr>`).join('')).join('');
 document.querySelectorAll('#prog td.c span[data-nom]').forEach(sp=>sp.onclick=()=>{if(!window.viewerShow){alert('vue 3D non disponible');return;}window.viewerShow(sp.dataset.stage,sp.dataset.nom,false);});
 document.querySelectorAll('#prog td.nom[data-nom]').forEach(td=>td.onclick=()=>{if(!window.viewerShow){alert('vue 3D non disponible');return;}window.viewerShow(null,td.dataset.nom,true);});
})();
// ---- tableau
const tb=document.querySelector('#tab tbody');
function cellA(a){if(!a)return '<span class="pill q0">aucune</span>';const q=a.qualite_note;
 return `<span class="pill ${QC[q]}">${QT[a.qualite]||'?'}</span> ${a.structures.length} structures · ${fmt(a.mm_per_voxel*1000)} µm/voxel · GL ${a.gl_mm} mm<br><span class="small mut">${a.glb?'GLB':'<s>GLB</s>'} · ${a.blend?'Blend':'<s>Blend</s>'} · ${a.nrrd_labels?'NRRD':'<s>NRRD</s>'} · ${a.rendu_organes?'rendus':'<s>rendus</s>'}${a.labels_corriges?' · labels corrigés':''}${a.manquantes.length?' · manque : '+a.manquantes.join(', '):''}${a.digestif?' · <b>digestif</b> '+a.digestif.segments.length+' seg.':''}${a.cardio?' · <b>cœur</b> '+(a.cardio.confiance_globale||'').split(' (')[0]:''}${a.vais?' · <b>aortes</b> '+(a.vais.confiance_globale||'').split(' (')[0]:''}${a.musc?' · <b>muscles</b> '+a.musc.segments.length+' seg.':''}</span>`;}
function cellB(b){if(!b)return '<span class="pill q0">aucune</span>';const q={dense:4,moyen:2,squelette:1}[b.densite];
 return `<span class="pill ${QC[q]}">${b.densite}</span> ${b.sections} sections · ${b.structures} structures · ${fmt(b.points)} points<br><span class="small mut">${b.labels} labels, ${b.traits} traits, ${b.anomalies} anomalies · ${b.legendes_json?'légendes ok':'<b>légendes à extraire</b>'}${b.viewer?' · viewer':''}</span>`;}
function cellC(js){if(!js.length)return '<span class="pill q0">aucune</span>';return js.map(J=>{const c=C[J];return `<span class="pill ${c.hull?'q4':'q2'}">${J}</span> ${c.plage} · CRL ${c.crl_min}-${c.crl_max} mm<br><span class="small mut">${c.captures} captures${c.hull?' · hull':''}${c.ratio_h_l?' · H/L '+fmt(c.ratio_h_l):''}</span>`;}).join('<br>');}
function detail(l){const a=l.A,b=l.B;let h='<div class="det-grid">';
 h+=`<div><h4>Repères</h4><div class="small">${l.reperes}</div><div class="small mut">Statut : ${l.statut}</div></div>`;
 if(a){const byCol={};a.structures.forEach(x=>{const sy=(D.struct_sys[x.nom]||['Autres'])[0];(byCol[sy]=byCol[sy]||[]).push(x)});if(a.digestif)a.digestif.segments.forEach(x=>{(byCol['Système digestif']=byCol['Système digestif']||[]).push({nom:x.nom+' (séparé, '+x.confiance+')',volume_mm3:x.volume_mm3})});if(a.cardio)a.cardio.segments.forEach(x=>{(byCol['Système vasculaire']=byCol['Système vasculaire']||[]).push({nom:x.nom+' (détouré, '+x.confiance+')',volume_mm3:x.volume_mm3})});if(a.vais)a.vais.segments.forEach(x=>{(byCol['Système vasculaire']=byCol['Système vasculaire']||[]).push({nom:x.nom+' (tracé, '+x.confiance+')',volume_mm3:x.volume_mm3})});if(a.musc)a.musc.segments.forEach(x=>{(byCol['Muscles']=byCol['Muscles']||[]).push({nom:x.nom+' ('+(a.origine==='B'?'coupes':'auto')+', '+x.confiance+')',volume_mm3:x.volume_mm3})});D.systemes.forEach(sy=>{if(!byCol[sy]&&sy!=='Autres')byCol[sy]=[]});
  h+=`<div><h4><span class="src A">A</span> ${a.dossier}/out — ${a.taille_out_mo} Mo</h4><div class="small mut">${a.qualite_txt}</div><ul class="small">`+Object.entries(byCol).map(([c,xs])=>`<li><b>${c}</b> : `+(xs.length?xs.map(x=>`${x.nom} <span class="mut">${fmt(x.volume_mm3)} mm³</span>`).join(', '):'<span class="mut">— aucune structure</span>')+'</li>').sort((u,v)=>D.systemes.findIndex(z=>u.includes('<b>'+z))-D.systemes.findIndex(z=>v.includes('<b>'+z))).join('')+'</ul>';
  h+=`<div class="small mut">Étapes : ${Object.entries(a.etapes).map(([k,v])=>v?k:'<s>'+k+'</s>').join(' › ')}</div>`;
  if(a.controle_B){const cb=a.controle_B;h+=`<div class="small" style="margin-top:6px">${a.ecarte_3d?'<span class="pill q1">3D écartée</span> ':''}<b>Topographie par zones et repérage vertébral (méthode retenue)</b> : `+Object.entries(cb.bilan_zones||{}).map(([z,v])=>{const tot=Object.values(v).reduce((x,y)=>x+y,0);const ok=v.complet||0;return `${z} <span class="pill ${ok===tot?'q4':ok?'q2':'q1'}" style="font-size:11px">${ok}/${tot}</span>`;}).join(' ')+`<div class="mut">${cb.n_masques} masques, ${cb.n_rejets} rejets, ${cb.n_chaines} chaînes, ${cb.n_coupes_interpolees} coupes interpolées${(cb.rapports&&cb.rapports.length?cb.rapports:(cb.rapport?[cb.rapport]:[])).map(r=>` · <a href="${r}" target="_blank" rel="noopener">${r.split('/').pop().replace('.html','').replace('controle_','rapport ')}</a>`).join('')}</div>`+(cb.axe_vertebral?`<div class="mut">axe vertébral courbe : arc ${fmt(cb.axe_vertebral.longueur_arc_mm)} mm (${fmt(cb.axe_vertebral.longueur_z_mm)} mm en z), obliquité des coupes ${typeof cb.axe_vertebral.obliquite_deg==='object'?`${cb.axe_vertebral.obliquite_deg.min}–${cb.axe_vertebral.obliquite_deg.max} ° (moy. ${cb.axe_vertebral.obliquite_deg.moyenne})`:cb.axe_vertebral.obliquite_deg+' °'}, hauteur vertébrale ${fmt(cb.axe_vertebral.hauteur_vertebrale_mm)} mm${cb.axe_vertebral.niveaux_manquants&&cb.axe_vertebral.niveaux_manquants.length?' · niveaux manquants : '+cb.axe_vertebral.niveaux_manquants.join(', '):''}</div>`:'')+(cb.etages&&cb.etages.niveaux?`<div class="mut">étages : ${Object.keys(cb.etages.niveaux).length} niveaux${cb.etages.pas_local_mm?` · pas local ${typeof cb.etages.pas_local_mm==='number'?fmt(cb.etages.pas_local_mm)+' mm':''}`:''} — ${Object.entries(cb.etages.niveaux).slice(0,40).map(([k,v])=>k+(v.origine&&v.origine!=='legende'?'*':'')).join(' ')}${Object.keys(cb.etages.niveaux).length>40?'…':''} <span class="small">(* = prédit)</span></div>`:'')+(cb.images&&cb.images.length?`<div class="thumbs">${cb.images.map(im=>`<a href="${im}" target="_blank" rel="noopener"><img loading="lazy" src="${im}" alt=""></a>`).join('')}</div>`:'')+(cb.symetries&&cb.symetries.length?`<div class="mut">complétées par symétrie (${cb.symetries.length}) : ${cb.symetries.map(x=>x.piece||x).join(', ')}</div>`:'')+(cb.inventaire_incomplet&&cb.inventaire_incomplet.length?`<div class="mut">inventaire incomplet (${cb.inventaire_incomplet.length}) : ${cb.inventaire_incomplet.slice(0,12).map(x=>typeof x==='string'?x:(x.nom||x.piece||JSON.stringify(x))).join(', ')}${cb.inventaire_incomplet.length>12?'…':''}</div>`:'')+`</div>`;}
  if(a.digestif)h+=`<div class="small" style="margin-top:6px"><b>Tube digestif (${a.origine==='B'?'depuis les coupes, hors vue 3D':'reconstruction séparée'}${a.digestif.aligne?', aligné':', <span style="color:var(--q1)">NON aligné</span>'})</b>${a.digestif.evaluation_visuelle?` — <i>${a.digestif.evaluation_visuelle}</i>`:''} : `+a.digestif.segments.map(x=>`${x.nom} <span class="mut">${x.confiance}, ${fmt(x.volume_mm3)} mm³</span>`).join(', ')+(a.digestif.manquants.length?` · <span class="mut">absents : ${a.digestif.manquants.join(', ')}</span>`:'')+`</div>`;
  if(a.topo){const c=a.topo.controle||{};const ax=c.axe_vertebral||{};const f=c.fiabilite||{};h+=`<div class="small" style="margin-top:4px"><b>Topographie par région (axe vertébral porté sur le volume)</b> <span class="pill ${a.topo.fiable?'q4':'q1'}" style="font-size:11px">${a.topo.fiable?'fiable':'non fiable'}</span> : arc ${fmt(ax.longueur_arc_mm)} mm (${fmt(ax.longueur_z_mm)} mm en z, écart latéral ${fmt(ax.ecart_lateral_mm)} mm)${c.niveaux_imposes?` · <b>${c.niveaux_imposes.sans_corps&&c.niveaux_imposes.sans_corps.length?`${c.niveaux_imposes.sans_corps.length} somites occipitaux + ${c.niveaux_imposes.n-c.niveaux_imposes.sans_corps.length} niveaux vertébraux`:`${c.niveaux_imposes.n} niveaux imposés`}</b> (${c.niveaux_imposes.composition||''}), ${c.niveaux_imposes.cales} frontières calées${c.niveaux_imposes.sur_prolongement&&c.niveaux_imposes.sur_prolongement.length?`, sur prolongement de l'axe : ${c.niveaux_imposes.sur_prolongement[0]} → ${c.niveaux_imposes.sur_prolongement[c.niveaux_imposes.sur_prolongement.length-1]}`:''}${c.niveaux_imposes.hors_axe&&c.niveaux_imposes.hors_axe.length?`, extrapolés hors axe : ${c.niveaux_imposes.hors_axe[0]} → ${c.niveaux_imposes.hors_axe[c.niveaux_imposes.hors_axe.length-1]}`:''}${c.niveaux_imposes.reperes?` · repères : ${c.niveaux_imposes.reperes}`:''}`:(f.n_etages?` · ${f.n_etages} étages, période ${fmt(f.periode_mm)} mm, cv ${fmt(f.cv_hauteurs)}`:'')}${c.zones&&c.zones.methode?` · zones : ${c.zones.methode}`:''}${c.video360?` · vidéo 360° : ${c.video360.n_images||'?'} images, ${(c.video360.vues_profil||[]).filter(v=>v.recalage_ok).length}/${(c.video360.vues_profil||[]).length} vues de profil recalées, +${(c.video360.vues_profil||[]).reduce((t,v)=>t+(v.n_supplementaires||0),0)} étages de queue`:''}${c.etages_manuels?` · <b>étages marqués à la main par l'utilisateur</b> : ${c.etages_manuels.n_marques||'?'} marques${c.etages_manuels.n_extrapoles?` + ${c.etages_manuels.n_extrapoles} extrapolé(s) côté queue`:''}, pas ${fmt(c.etages_manuels.pas_median_mm)} mm, recalage ${fmt(c.etages_manuels.recalage_erreur_px)} px${c.etages_manuels.vue_reference_calage?' · vue de référence de calage':''}`:''}${c.etages_manuels_dorsal?` · <b>étages marqués à la main (vue de dos)</b> : ${c.etages_manuels_dorsal.n_marques||'?'} marques, pas ${fmt(c.etages_manuels_dorsal.pas_median_mm)} mm, recalage ${fmt(c.etages_manuels_dorsal.recalage_erreur_px)} px`:''}${c.avertissement?`<div class="mut">⚠ ${c.avertissement}</div>`:''}${c.reperes_utilisateur?(()=>{const r=c.reperes_utilisateur;const pil=(ok,t)=>`<span class="pill ${ok?'q4':'q1'}" style="font-size:11px">${t}</span>`;const mm=x=>x==null?'—':fmt(x)+' mm';
 return `<div class="small" style="margin-top:2px"><b>Repères pointés par l'utilisateur (queue-VO)</b> ${pil(r.oeil_valide,'œil '+(r.oeil_valide?'ok':'écarté'))} ${pil(r.queue_valide,'queue '+(r.queue_valide?'ok':'écartée'))} : recalage ${fmt(r.recalage_erreur_px)} px (${r.vue||''}) · œil ↔ label yeux ${mm(r.ecart_oeil_label_yeux_mm)} (sagittal)${r.ecart_oeil_3d_label_yeux_mm!=null?`, ${mm(r.ecart_oeil_3d_label_yeux_mm)} en 3D (${r.oeil_n_vues||'?'} vues, résidu ${mm(r.oeil_residu_mm)})`:''} · queue ↔ tube neural ${mm(r.ecart_queue_tube_neural_mm)}, ↔ bout de l'enveloppe ${mm(r.ecart_queue_bout_enveloppe_mm)}${r.motifs&&r.motifs.length?`<div class="mut">${[].concat(r.motifs).join(' ; ')}</div>`:''}</div>`;})():''}${(c.images_url||[]).length?`<div class="mut">contrôles : ${c.images_url.map(u=>`<a href="${u}" target="_blank">${u.split('/').pop().replace(/^CS\d+_/,'')}</a>`).join(' · ')}</div>`:''}</div>`;}
  if(a.musc)h+=`<div class="small" style="margin-top:4px"><b>Muscles (${a.origine==='B'?'depuis les coupes':'automatique'}${a.musc.aligne?', aligné':', <span style="color:var(--q1)">NON aligné</span>'})</b> : `+a.musc.segments.map(x=>`${x.nom} <span class="mut">${x.confiance}${x.n_blocs!=null?', '+x.n_blocs+' blocs':''}, ${fmt(x.volume_mm3)} mm³</span>`).join(', ')+`</div>`;
  ['digestif','cardio','vais','musc'].forEach(k=>{const r=a[k];if(r&&r.controle&&r.controle.bilan_zones)h+=`<div class="small"><b>Contrôle par zones (${k})</b> : `+Object.entries(r.controle.bilan_zones).map(([z,v])=>{const tot=Object.values(v).reduce((x,y)=>x+y,0);const ok=v.complet||0;return `${z} <span class="pill ${ok===tot?'q4':ok?'q2':'q1'}" style="font-size:11px">${ok}/${tot}</span>`;}).join(' ')+(r.controle.inventaire_incomplet&&r.controle.inventaire_incomplet.length?`<div class="mut">inventaire incomplet : ${r.controle.inventaire_incomplet.slice(0,12).map(x=>typeof x==='string'?x:(x.nom||x.piece||JSON.stringify(x))).join(', ')}</div>`:'')+`</div>`;});
  if(a.vais)h+=`<div class="small" style="margin-top:4px"><b>Aortes dorsales (tracé séparé${a.vais.aligne?', aligné':', <span style="color:var(--q1)">NON aligné</span>'})</b> : `+a.vais.segments.map(x=>`${x.nom} <span class="mut">${fmt(x.volume_mm3)} mm³</span>`).join(', ')+` · <span class="mut">confiance ${a.vais.confiance_globale||'?'}</span></div>`;
  if(a.cardio)h+=`<div class="small" style="margin-top:4px"><b>Cœur détouré (séparé${a.cardio.aligne?', aligné':', <span style="color:var(--q1)">NON aligné</span>'})</b> : `+a.cardio.segments.map(x=>`${x.nom} <span class="mut">${fmt(x.volume_mm3)} mm³</span>`).join(', ')+` · <span class="mut">confiance ${a.cardio.confiance_globale||'?'}</span></div>`;
  h+=`<div class="thumbs">${a.rendu_organes?`<img loading="lazy" src="${a.rendu_organes_url}" alt="">`:''}${a.rendu_peau?`<img loading="lazy" src="${a.rendu_peau_url}" alt="">`:''}</div></div>`;}
 if(b){h+=`<div><h4><span class="src B">B</span> coupes embryos 9-23/embryo_stage_${l.stade.slice(2)}</h4><div class="small">Structures les plus documentées :</div><ul class="small">`+b.top.map(t=>`<li>${t.nom} <span class="mut">${t.points} pts</span></li>`).join('')+'</ul>';
  if(b.overlay_url)h+=`<div class="thumbs"><img loading="lazy" src="${b.overlay_url}" alt=""></div>`;h+='</div>';}
 if(l.C.length){h+=`<div><h4><span class="src C">C</span> VOKA</h4>`+l.C.map(J=>{const c=C[J];return `<div class="small"><b>${J}</b> (${c.semaine} sem.) — ${c.reperes}<br><span class="mut">largeur ${fmt(c.largeur_mm)} × hauteur ${fmt(c.hauteur_mm)} mm (vue latérale, échelle CRL litt.) · <a href="${c.page}" target="_blank" rel="noopener">page</a></span></div>`+(c.image_url?`<div class="thumbs"><img loading="lazy" src="${c.image_url}" alt=""></div>`:'');}).join('')+'</div>';}
 return h+'</div>';}
function rendTab(){const q=(document.getElementById('q').value||'').toLowerCase();const gaps=document.getElementById('onlyGaps').checked;tb.innerHTML='';
 D.lignes.forEach(l=>{if(gaps&&l.score.interne3d>=3)return;
  const txt=JSON.stringify(l).toLowerCase();if(q&&!txt.includes(q))return;
  const tr=document.createElement('tr');tr.className='st';
  tr.innerHTML=`<td><b>${l.stade}</b></td><td>${l.jours[0]}–${l.jours[1]}</td><td>${l.crl[0]}–${l.crl[1]}</td><td>${l.statut}</td><td data-src="A">${cellA(l.A)}</td><td data-src="B">${cellB(l.B)}</td><td data-src="C">${cellC(l.C)}</td>`;
  const td=document.createElement('tr');td.className='det';td.style.display='none';td.innerHTML=`<td colspan="7">${detail(l)}</td>`;
  tr.onclick=()=>{td.style.display=td.style.display=='none'?'':'none'};tb.appendChild(tr);tb.appendChild(td);});}
document.getElementById('q').oninput=rendTab;document.getElementById('onlyGaps').onchange=rendTab;rendTab();
// ---- todo
const LANCEUR='http://127.0.0.1:8791';const LOCAL=/^(localhost|127\.0\.0\.1)$/.test(location.hostname);let lanceurOk=false;const STL={a_faire:'à faire',en_cours:'en cours',termine:'terminé',echec:'échec'};
function promptTache(t){return `Tâche de l'agrégateur embryon 3D (id ${t.id}, priorité ${t.prio}, source ${t.source}, stade ${t.stade}) :\n${t.action}\n${t.detail||''}\n${t.cmd?'Commande suggérée : '+t.cmd+'\n':''}Contexte : C:\\Users\\MicroTurtle\\Documents\\Claude (embryo3d/README.md, embryons_3D/LISEZMOI.txt, page embryons_3D/agregateur.html). Quand c'est fait, mettre le statut « termine » pour l'id ${t.id} dans embryons_3D/taches_etat.json et relancer python embryo3d/agregateur.py.`;}
async function api(path,body){const r=await fetch(LANCEUR+path,{method:body?'POST':'GET',headers:{'Content-Type':'application/json'},body:body?JSON.stringify(body):undefined});const j=await r.json();if(!r.ok)throw new Error(j.erreur||r.status);return j;}
async function majEtat(){if(!LOCAL){try{const r=await fetch('taches_etat.json?t='+Date.now(),{cache:'no-store'});if(r.ok){const e=await r.json();D.todo.forEach(t=>{const x=e[t.id];if(x){t.instance=x.instance||'';t.statut=x.statut||'a_faire';t.maj=x.maj||'';}});document.getElementById('lanceurEtat').textContent='· état publié le '+(e._publie||'?')+' (relu toutes les 30 s)';document.getElementById('lanceurEtat').className='on';}}catch(e){}rendTodo();return;}try{const j=await api('/etat');lanceurOk=true;D.todo.forEach(t=>{const e=j.etat[t.id];if(e){t.instance=e.instance||'';t.statut=e.statut||'a_faire';t.maj=e.maj||'';t.log=e.log;}});
 document.getElementById('lanceurEtat').textContent='· lanceur local connecté';document.getElementById('lanceurEtat').className='on';}
 catch(e){lanceurOk=false;document.getElementById('lanceurEtat').textContent='· lanceur local absent (python embryo3d/lanceur.py)';document.getElementById('lanceurEtat').className='off';}rendTodo();}
async function lancer(t){try{await api('/lancer',{id:t.id});}catch(e){alert('Lancement impossible : '+e.message);}majEtat();}
async function confier(t){const inst=prompt('Confier la tâche à quelle instance ? (nom de la session Claude, ou personne)',t.instance||'');if(inst===null)return;
 try{await api('/confier',{id:t.id,instance:inst});}catch(e){alert(e.message);}try{await navigator.clipboard.writeText(promptTache(t));}catch(e){}majEtat();}
async function statut(t,st){try{await api('/statut',{id:t.id,statut:st});}catch(e){alert(e.message);}majEtat();}
function rendEnCours(){const ec=D.todo.filter(t=>t.statut==='en_cours'),ok=D.todo.filter(t=>t.statut==='termine'),ko=D.todo.filter(t=>t.statut==='echec');
 const byInst={};ec.forEach(t=>{const k=t.instance||'(sans instance)';(byInst[k]=byInst[k]||[]).push(t)});
 document.getElementById('enCours').innerHTML=`<div class="card"><h3>En cours : ${ec.length}</h3>${Object.entries(byInst).map(([i,ts])=>`<p><b>${i}</b> — ${ts.map(t=>`${t.stade} · ${t.action}`).join(' ; ')}</p>`).join('')||'<p>aucune tâche en cours</p>'}</div><div class="card"><h3>Terminées : ${ok.length} · échecs : ${ko.length}</h3><p>${ok.slice(-5).map(t=>`${t.stade} · ${t.action}${t.instance?' <span class="mut">('+t.instance+')</span>':''}`).join('<br>')||'—'}</p></div>`;}
function rendTodo(){rendEnCours();const on=[...document.querySelectorAll('#todoF input')].filter(i=>i.checked).map(i=>i.value);
 const el=document.getElementById('todo');el.innerHTML='<div class="p mut">prio</div><div class="mut small">source</div><div class="mut small">stade</div><div class="mut small">action</div><div class="mut small">instance · statut</div><div class="mut small">lancer</div>';let n=0;
 D.todo.forEach(t=>{const st=t.statut||'a_faire';const srcK='ABC'.includes(t.source[0])?t.source[0]:null;if(!on.includes(t.prio)||(srcK&&!on.includes(srcK))||!on.includes('s:'+st))return;n++;
  const d=document.createElement('div');d.className='act';
  const bL=document.createElement('button');bL.textContent='▶ Lancer';bL.title=t.cmd?'exécuter la commande sur ce PC':'pas de commande : à confier';bL.disabled=!(LOCAL&&lanceurOk&&t.cmd&&st!=='en_cours');bL.onclick=()=>lancer(t);
  const bC=document.createElement('button');bC.textContent='Confier…';bC.disabled=!(LOCAL&&lanceurOk);bC.onclick=()=>confier(t);
  const bF=document.createElement('button');bF.textContent=st==='termine'?'Rouvrir':'Fait';bF.disabled=!(LOCAL&&lanceurOk);bF.onclick=()=>statut(t,st==='termine'?'a_faire':'termine');
  const bP=document.createElement('button');bP.textContent='📋';bP.title='copier le prompt pour une session Claude';bP.onclick=async()=>{try{await navigator.clipboard.writeText(promptTache(t));bP.textContent='✓';setTimeout(()=>bP.textContent='📋',1200);}catch(e){alert(promptTache(t));}};
  d.append(bL,bC,bF,bP);
  el.insertAdjacentHTML('beforeend',`<div class="p ${t.prio}">${t.prio}</div><div><span class="src ${'ABC'.includes(t.source[0])?t.source[0]:''}">${t.source}</span></div><div><b>${t.stade}</b></div><div>${t.action}${t.detail?`<div class="small mut">${t.detail}</div>`:''}${t.cmd?`<div><code>${t.cmd}</code></div>`:''}<div class="small mut">id ${t.id}</div></div><div class="inst"><span class="st ${st}">${STL[st]||st}</span> ${t.instance?'<b>'+t.instance+'</b>':'<span class="mut">—</span>'}${t.maj?`<div class="small mut">${t.maj}</div>`:''}${t.log&&LOCAL?`<div class="small"><a href="${LANCEUR}/log/${t.id}" target="_blank" rel="noopener">journal</a></div>`:''}</div>`);
  el.appendChild(d);});
 document.getElementById('ntodo').textContent=`(${n} sur ${D.todo.length})`;}
document.querySelectorAll('#todoF input').forEach(i=>i.onchange=rendTodo);rendTodo();majEtat();setInterval(majEtat,LOCAL?5000:30000);
window.__D=D;
</script>
<script type="importmap">{"imports":{"three":"https://cdn.jsdelivr.net/npm/three@0.170.0/build/three.module.js","three/addons/":"https://cdn.jsdelivr.net/npm/three@0.170.0/examples/jsm/","three-mesh-bvh":"https://cdn.jsdelivr.net/npm/three-mesh-bvh@0.7.8/build/index.module.js"}}</script>
<script type="module">
const $=id=>document.getElementById(id);if(!$('gl')){throw new Error('viewer intégré remplacé par 3dh (iframe) — module inactif');}const msg=t=>{$('v3d-msg').textContent=t;$('v3d-msg').style.display=t?'':'none'};
const D=window.__D;const ST=D.lignes.filter(l=>l.A&&l.A.vue3d&&l.A.vue3d.glb).map(l=>({id:l.stade,j0:l.jours[0],j1:l.jours[1],mid:(l.jours[0]+l.jours[1])/2,dig:(l.A.digestif&&!l.A.ecarte_3d)?l.A.digestif.segments.filter(x=>x.fichier_ok):[],cardio:(l.A.cardio&&!l.A.ecarte_3d)?l.A.cardio.segments.filter(x=>x.fichier_ok):[],vais:(l.A.vais&&!l.A.ecarte_3d)?l.A.vais.segments.filter(x=>x.fichier_ok):[],musc:(l.A.musc&&!l.A.ecarte_3d)?l.A.musc.segments.filter(x=>x.fichier_ok):[],topo:(l.A.topo&&!l.A.ecarte_3d)?l.A.topo.segments.filter(x=>x.fichier_ok):[],b3:l.A.b3_autorisees||null,...l.A.vue3d}));
const DMAX=266;$('day').max=DMAX;
$('stageTicks').textContent='Stades reconstruits : '+ST.map(s=>`${s.id} J${s.j0}-${s.j1}`).join(' · ');
$('stageBtns').innerHTML=ST.map(s=>`<button data-stage="${s.id}" title="J${s.j0}-${s.j1}">${s.id}</button>`).join('');
const bqL=l=>{const sc=l.score;return sc.interne3d>=3?4:sc.interne3d>0?2:(sc.coupes>=3?3:sc.enveloppe>=3?2:sc.coupes?1:0);};
$('dayBtns').innerHTML=D.lignes.map(l=>`<button class="q${bqL(l)}" data-day="${l.frise[0]}" data-stade="${l.stade}" title="${l.stade} · J${l.jours[0]}-${l.jours[1]} · CRL ${l.crl[0]}-${l.crl[1]} mm · ${l.statut}">J${l.jours[0]}<b>${l.stade}</b></button>`).join('');if(!ST.some(s=>s.hull))$('hullLbl').style.display='none';
if(location.protocol==='file:'){msg('Le GLB ne peut pas être lu en file:// — ouvrir la page par un serveur HTTP (localhost:8766/embryons_3D/agregateur.html ou harcelon.fr/3dht/).');}
else try{
const THREE=await import('three');const {OrbitControls}=await import('three/addons/controls/OrbitControls.js');
const {GLTFLoader}=await import('three/addons/loaders/GLTFLoader.js');const {PLYLoader}=await import('three/addons/loaders/PLYLoader.js');
let MeshBVH=null;try{({MeshBVH}=await import('three-mesh-bvh'));}catch(e){console.warn('three-mesh-bvh indisponible : pas d\'interpolation',e);}
const canvas=$('gl');const renderer=new THREE.WebGLRenderer({canvas,antialias:true,preserveDrawingBuffer:true});renderer.setPixelRatio(Math.min(devicePixelRatio,2));
const scene=new THREE.Scene();scene.background=new THREE.Color(0x2b2d31);const V3SCENE=scene;
scene.add(new THREE.HemisphereLight(0xffffff,0x665544,3.2));const key=new THREE.DirectionalLight(0xffffff,3.2);key.position.set(-3,2,4);scene.add(key);const fill=new THREE.DirectionalLight(0xffffff,0.9);fill.position.set(3,-2,-1);scene.add(fill);
const cam=new THREE.OrthographicCamera(-1,1,1,-1,0.1,5000);let aspect=1;
function resize(){const w=canvas.clientWidth,h=canvas.clientHeight;renderer.setSize(w,h,false);aspect=w/h;cam.left=-aspect;cam.right=aspect;cam.top=1;cam.bottom=-1;cam.updateProjectionMatrix();}
new ResizeObserver(resize).observe(canvas);resize();
// ---- trajectoire
let TR=JSON.parse($('traj').value);$('trajApply').onclick=()=>{try{TR=JSON.parse($('traj').value);halfPrev=0;}catch(e){alert('JSON invalide : '+e.message)}};
const SD=new THREE.Vector2(-0.6,0.8).normalize();
function pose(azOff,elOff,roll){ // repère pipeline : X gauche→droite, Y ventral→dorsal, Z haut. Vue côté gauche : on regarde vers +X.
 const t=THREE.MathUtils.degToRad(TR.tilt||15);const Dd=new THREE.Vector3(0,1,0);
 let v=new THREE.Vector3(Math.cos(t),-Math.sin(t),0).normalize();
 v.applyAxisAngle(new THREE.Vector3(0,0,1),THREE.MathUtils.degToRad(azOff||0));
 const ax=new THREE.Vector3().crossVectors(new THREE.Vector3(0,0,1),v).normalize();v.applyAxisAngle(ax,THREE.MathUtils.degToRad(elOff||0)).normalize();
 const Dp=Dd.clone().sub(v.clone().multiplyScalar(Dd.dot(v))).normalize();const w=new THREE.Vector3().crossVectors(v,Dp);
 const up=Dp.clone().multiplyScalar(SD.y).add(w.clone().multiplyScalar(-SD.x)).normalize();if(roll)up.applyAxisAngle(v,THREE.MathUtils.degToRad(roll));const right=new THREE.Vector3().crossVectors(v,up).normalize();
 const q=new THREE.Quaternion().setFromRotationMatrix(new THREE.Matrix4().makeBasis(right,up,v.clone().negate()));return {v,up,right,q};}
// ---- modèles
const loaded={},hulls={};let cur=null,root=null,hullMesh=null,skin=null,box=new THREE.Box3(),center=new THREE.Vector3(),radius=10,halfPrev=0;
function setSkin(){if(skin){if($('skinContour')&&$('skinContour').checked){if(!skin.userData.contour){skin.userData.contour=true;skin.material=new THREE.ShaderMaterial({transparent:true,depthWrite:false,side:THREE.DoubleSide,uniforms:{col:{value:new THREE.Color(0xf0d0c0)},pw:{value:2.2},mx:{value:0.6}},
  vertexShader:'varying vec3 vN;varying vec3 vV;void main(){vN=normalize(normalMatrix*normal);vec4 mv=modelViewMatrix*vec4(position,1.0);vV=normalize(-mv.xyz);gl_Position=projectionMatrix*mv;}',
  fragmentShader:'uniform vec3 col;uniform float pw;uniform float mx;varying vec3 vN;varying vec3 vV;void main(){float f=pow(1.0-abs(dot(normalize(vN),normalize(vV))),pw);gl_FragColor=vec4(col,f*mx);}'});}
  skin.userData.recol=true;updateShown();return;}skin.material.transparent=$('skinT').checked;skin.material.opacity=$('skinT').checked?0.22:1;skin.material.depthWrite=!$('skinT').checked;skin.material.vertexColors=false;skin.material.color=new THREE.Color(0xe8b8a0);skin.material.roughness=0.6;}updateShown();}
window.updateShown=updateShown;function updateShown(){const s=ST.find(x=>x.id===cur);if(!s){$('shown').textContent='';return;}if(true){$('shown').innerHTML='';return;}const parts=[reconOnly?`<b>${s.id}</b> · <b>reconstructions seules</b> (modèle automatique du pipeline et derme masqués)`:(window.srcOff&&window.srcOff.has('A'))?`<b>${s.id}</b> · pipeline masqué (cadre A désactivé)`:`<b>${s.id}</b> · <span class="src A">A</span> pipeline (GLB, couleurs des structures ; <b>derme = « enveloppe »</b>, ${$('skinT').checked?'translucide rosée':'opaque'})`];
 Object.entries(RECONS).forEach(([k,R])=>{if($(R.cb).checked&&s[R.champ].length)parts.push(`<span class="src A">A</span> ${R.sg.toLowerCase()}`);});
 if($('hull').checked&&s.hull)parts.push('<span class="src C">C</span> hull VOKA '+s.J+' (translucide beige, centré)');
 $('shown').innerHTML='Affiché : '+parts.join(' · ');
 document.querySelectorAll('#stageBtns button').forEach(b=>b.classList.toggle('on',b.dataset.stage===cur));}
$('skinT').onchange=setSkin;$('skinContour').onchange=()=>{if(skin){skin.userData.contour=false;skin.userData.recol=false;skin.material=new THREE.MeshStandardMaterial({color:0xe8b8a0});}setSkin();applyLayers();};$('wire').onchange=()=>root&&root.traverse(o=>{if(o.isMesh)o.material.wireframe=$('wire').checked});
const loading={};function preload(s){if(!s)return;Object.values(RECONS).forEach(R=>s[R.champ].forEach(seg=>{if(!recCache[seg.url]&&!recCache['@'+seg.url]){recCache['@'+seg.url]=true;new PLYLoader().load(seg.url,geo=>{geo.computeVertexNormals();recCache[seg.url]=geo;},undefined,()=>{delete recCache['@'+seg.url];});}}));
 const k=keyOf(s);if(loaded[k]||loading[k])return;loading[k]=true;new GLTFLoader().load(glbOf(s),g=>{loaded[k]=g.scene;loading[k]=false;},undefined,()=>{loading[k]=false;});}
// ---- interpolation : chaque maillage de reconstruction du stade courant reçoit comme cible morph sa projection sur l'homologue (même nom) du stade suivant
const morphDone=new Set();const morphQueue=[];
function targetsOf(next){const out={};const sc=loaded[keyOf(next)];if(sc)sc.traverse(o=>{if(o.isMesh)out[o.name]=o.geometry;});Object.entries(RECONS).forEach(([k,R])=>next[R.champ].forEach(seg=>{const g=recCache[seg.url];if(g&&g.attributes)out[k+':'+seg.nom]=g;}));return out;}
function displayedRecon(){const list=[];if(root)root.traverse(o=>{if(o.isMesh&&RECON_RE.test(o.name))list.push([o.name,o]);});Object.entries(recGroup).forEach(([k,g])=>{if(g)g.children.forEach(m=>list.push([k+':'+m.name,m]));});return list;}
function ensureMorph(cur_,next){if(!MeshBVH||!next||!loaded[keyOf(next)])return;const tg=targetsOf(next);
 displayedRecon().forEach(([key,m])=>{const id=cur_.id+'|'+next.id+'|'+key;if(morphDone.has(id))return;const t=tg[key];if(!t){morphDone.add(id);return;}morphDone.add(id);morphQueue.push({m,t,key,next:next.id});});}
const MORPH_CHUNK=15000;   // sommets projetés par image (≈ 8 ms) : la lecture reste fluide
function morphStep(){if(!MORPH_ACTIF)return;const job=morphQueue[0];if(!job)return;const {m,t}=job;try{if(!t.boundsTree)t.boundsTree=new MeshBVH(t);const pos=m.geometry.attributes.position;
  if(!job.out){job.out=new Float32Array(pos.count*3);job.i=0;}
  const p=new THREE.Vector3(),tgt={point:new THREE.Vector3()};const end=Math.min(pos.count,job.i+MORPH_CHUNK);
  for(let i=job.i;i<end;i++){p.fromBufferAttribute(pos,i);t.boundsTree.closestPointToPoint(p,tgt);job.out[3*i]=tgt.point.x;job.out[3*i+1]=tgt.point.y;job.out[3*i+2]=tgt.point.z;}
  job.i=end;if(end<pos.count)return;
  m.geometry.morphAttributes.position=[new THREE.Float32BufferAttribute(job.out,3)];m.geometry.morphTargetsRelative=false;m.updateMorphTargets();m.material.needsUpdate=true;m.userData.morphTo=job.next;}catch(e){console.warn('morph',job.key,e);}
 morphQueue.shift();}
const MORPH_ACTIF=false;   // 24/09 « plus jamais ça » : la projection par point le plus proche crée des pointes entre pièces non correspondantes → désactivé
function applyMorph(){if(!MORPH_ACTIF)return;const s=stageAt(day);const next=s?nextOf(s):null;if(!s||!next){return;}ensureMorph(s,next);if(!morphDone.size)return;
 const t=Math.max(0,Math.min(1,(day-s.mid)/Math.max(0.01,next.mid-s.mid)));
 displayedRecon().forEach(([key,m])=>{if(m.morphTargetInfluences&&m.morphTargetInfluences.length&&m.userData.morphTo===next.id)m.morphTargetInfluences[0]=t;});}
const glbOf=s=>(reconOnly||!s.glb_full)?s.glb:s.glb_full;const keyOf=s=>s.id+(reconOnly||!s.glb_full?'':'|full');
function showStage(s){if(!s){if(cur!==null){cur=null;if(root){scene.remove(root);root=null;}if(hullMesh){scene.remove(hullMesh);hullMesh=null;}Object.keys(recGroup).forEach(k=>{if(recGroup[k]){scene.remove(recGroup[k]);recGroup[k]=null;}});$('vokaImg').style.display='none';$('vokaLbl').textContent='VOKA correspondant : —';buildLayers();}
  const l=D.lignes.slice().reverse().find(x=>day>=x.frise[0])||D.lignes[0];msg('');$('shown').innerHTML=`<b>${l.stade}</b> (J${l.jours[0]}-${l.jours[1]}) · pas de reconstruction 3D · données : ${l.statut}`;return;}
 if(cur===keyOf(s))return;cur=keyOf(s);
 if(root){scene.remove(root);root=null;}if(hullMesh){scene.remove(hullMesh);hullMesh=null;}Object.keys(recGroup).forEach(k=>{if(recGroup[k]){scene.remove(recGroup[k]);recGroup[k]=null;}});
 const place=g=>{if(cur!==keyOf(s))return;root=g;scene.add(root);if(soloFilter){Object.keys(layerOn).forEach(k=>delete layerOn[k]);}box.setFromObject(root);box.getCenter(center);const sz=box.getSize(new THREE.Vector3());radius=Math.max(sz.x,sz.y,sz.z);
  let best=null,bv=0;root.traverse(o=>{if(o.isMesh){o.material.side=THREE.DoubleSide;if(o.name==='enveloppe')best=o;const b=new THREE.Box3().setFromObject(o).getSize(new THREE.Vector3());const vol=b.x*b.y*b.z;if(!best&&vol>bv){bv=vol;}}});skin=best /* peau = maillage nomme 'enveloppe' seulement (CS23 issu des coupes n'en a pas) */;setSkin();
  msg(`${s.id} — ${sz.x.toFixed(1)} × ${sz.y.toFixed(1)} × ${sz.z.toFixed(1)} mm`);setTimeout(()=>{if(cur===s.id)msg('')},2500);showHull();showDig();applyLayers();buildLayers();};
 const k=keyOf(s);if(loaded[k]){place(loaded[k]);}else{msg(`chargement ${s.id}…`);new GLTFLoader().load(glbOf(s),g=>{loaded[k]=g.scene;place(g.scene);},x=>{if(x.total)msg(`chargement ${s.id} : ${Math.round(100*x.loaded/x.total)} %`)},e=>msg('erreur GLB : '+e.message));}
 $('vokaLbl').textContent='VOKA correspondant : '+(s.J||'—');if(s.voka_img){$('vokaImg').src=s.voka_img;$('vokaImg').style.display='';}else $('vokaImg').style.display='none';}
// matériau « contours » : opaque sur les bords (normale ⟂ vue), transparent de face → silhouette légère
function hullMaterial(mode){if(mode==='fil')return new THREE.MeshBasicMaterial({color:0xf0d9c8,wireframe:true,transparent:true,opacity:0.18,depthWrite:false});
 if(mode==='voile')return new THREE.MeshStandardMaterial({color:0xe8c9bd,transparent:true,opacity:0.12,depthWrite:false,side:THREE.DoubleSide,roughness:0.9});
 return new THREE.ShaderMaterial({transparent:true,depthWrite:false,side:THREE.DoubleSide,uniforms:{col:{value:new THREE.Color(0xffd9c2)},pw:{value:3.0},mx:{value:0.45}},
  vertexShader:'varying vec3 vN;varying vec3 vV;void main(){vN=normalize(normalMatrix*normal);vec4 mv=modelViewMatrix*vec4(position,1.0);vV=normalize(-mv.xyz);gl_Position=projectionMatrix*mv;}',
  fragmentShader:'uniform vec3 col;uniform float pw;uniform float mx;varying vec3 vN;varying vec3 vV;void main(){float f=pow(1.0-abs(dot(normalize(vN),normalize(vV))),pw);gl_FragColor=vec4(col,f*mx);}'});}
const hullFit={};   // par stade : {matrix, residu}
function showHull(){const s=ST.find(x=>x.id===cur);if(hullMesh){scene.remove(hullMesh);hullMesh=null;}$('hullInfo').textContent='';if(true)return; /* VOKA retiré du viewer (24/09) */ if(!$('hull').checked||!s||!s.hull)return;
 const put=geo=>{if(cur!==s.id||!$('hull').checked)return;const m=new THREE.Mesh(geo,hullMaterial($('hullMode').value));m.name='hull_voka';
  const c=new THREE.Box3().setFromBufferAttribute(geo.attributes.position).getCenter(new THREE.Vector3());m.position.copy(center).sub(c);hullMesh=m;scene.add(m);
  if(hullFit[s.id]){m.matrixAutoUpdate=false;m.matrix.copy(hullFit[s.id].matrix);$('hullInfo').textContent=`calage 3D : écart moyen ${hullFit[s.id].residu.toFixed(3)} mm (${hullFit[s.id].init})`;}
  else fitHull(m,s);};
 if(hulls[s.J])put(hulls[s.J]);else new PLYLoader().load(s.hull,geo=>{geo.computeVertexNormals();hulls[s.J]=geo;put(geo);});}
$('hullMode').onchange=showHull;
// ---- calage 3D de l'enveloppe VOKA sur la peau du modèle : ICP avec échelle (Umeyama), 4 orientations initiales, la meilleure gagne
function svd3(H){ // Jacobi sur H^T H → V, valeurs ; U = H V / σ
 const A=[[0,0,0],[0,0,0],[0,0,0]];for(let i=0;i<3;i++)for(let j=0;j<3;j++)for(let k=0;k<3;k++)A[i][j]+=H[k][i]*H[k][j];
 const V=[[1,0,0],[0,1,0],[0,0,1]];for(let it=0;it<30;it++){let off=0;for(let p_=0;p_<3;p_++)for(let q=p_+1;q<3;q++){off+=A[p_][q]*A[p_][q];if(Math.abs(A[p_][q])<1e-12)continue;const th=0.5*Math.atan2(2*A[p_][q],A[q][q]-A[p_][p_]);const c=Math.cos(th),s_=Math.sin(th);
  for(let k=0;k<3;k++){const akp=A[k][p_],akq=A[k][q];A[k][p_]=c*akp-s_*akq;A[k][q]=s_*akp+c*akq;}for(let k=0;k<3;k++){const apk=A[p_][k],aqk=A[q][k];A[p_][k]=c*apk-s_*aqk;A[q][k]=s_*apk+c*aqk;}for(let k=0;k<3;k++){const vkp=V[k][p_],vkq=V[k][q];V[k][p_]=c*vkp-s_*vkq;V[k][q]=s_*vkp+c*vkq;}}if(off<1e-18)break;}
 const sig=[0,1,2].map(i=>Math.sqrt(Math.max(0,A[i][i])));const U=[[0,0,0],[0,0,0],[0,0,0]];for(let j=0;j<3;j++){for(let i=0;i<3;i++){let v=0;for(let k=0;k<3;k++)v+=H[i][k]*V[k][j];U[i][j]=sig[j]>1e-12?v/sig[j]:0;}}
 return {U,V,sig};}
function umeyama(P,Q){ // P (source) → Q (cible), retourne {s,R,t}
 const n=P.length;const mp=[0,0,0],mq=[0,0,0];P.forEach(p_=>{mp[0]+=p_[0];mp[1]+=p_[1];mp[2]+=p_[2];});Q.forEach(q=>{mq[0]+=q[0];mq[1]+=q[1];mq[2]+=q[2];});for(let i=0;i<3;i++){mp[i]/=n;mq[i]/=n;}
 const H=[[0,0,0],[0,0,0],[0,0,0]];let vp=0;for(let k=0;k<n;k++){const p_=[P[k][0]-mp[0],P[k][1]-mp[1],P[k][2]-mp[2]],q=[Q[k][0]-mq[0],Q[k][1]-mq[1],Q[k][2]-mq[2]];vp+=p_[0]*p_[0]+p_[1]*p_[1]+p_[2]*p_[2];for(let i=0;i<3;i++)for(let j=0;j<3;j++)H[i][j]+=q[i]*p_[j];}
 const {U,V,sig}=svd3(H);const R=[[0,0,0],[0,0,0],[0,0,0]];for(let i=0;i<3;i++)for(let j=0;j<3;j++)for(let k=0;k<3;k++)R[i][j]+=U[i][k]*V[j][k];
 const det=R[0][0]*(R[1][1]*R[2][2]-R[1][2]*R[2][1])-R[0][1]*(R[1][0]*R[2][2]-R[1][2]*R[2][0])+R[0][2]*(R[1][0]*R[2][1]-R[1][1]*R[2][0]);let d3=1;if(det<0){d3=-1;for(let i=0;i<3;i++)for(let j=0;j<3;j++){R[i][j]=0;for(let k=0;k<3;k++)R[i][j]+=U[i][k]*(k===2?-1:1)*V[j][k];}}
 const sc=(sig[0]+sig[1]+d3*sig[2])/vp;const t=[0,1,2].map(i=>mq[i]-sc*(R[i][0]*mp[0]+R[i][1]*mp[1]+R[i][2]*mp[2]));return {s:sc,R,t};}
function fitHull(m,s){if(!MeshBVH||!root)return;let skin=null;root.traverse(o=>{if(o.isMesh&&o.name==='enveloppe')skin=o;});if(!skin){$('hullInfo').textContent='pas de peau du pipeline : enveloppe VOKA simplement centrée';return;}
 const tg=skin.geometry;if(!tg.boundsTree)tg.boundsTree=new MeshBVH(tg);
 const pos=m.geometry.attributes.position;const N=Math.min(2500,pos.count);const step=Math.max(1,Math.floor(pos.count/N));const src=[];for(let i=0;i<pos.count;i+=step)src.push([pos.getX(i)+m.position.x,pos.getY(i)+m.position.y,pos.getZ(i)+m.position.z]);
 const tgt={point:new THREE.Vector3()};const p=new THREE.Vector3();const c0=[center.x,center.y,center.z];
 const inits=[['identité',[[1,0,0],[0,1,0],[0,0,1]]],['180° X',[[1,0,0],[0,-1,0],[0,0,-1]]],['180° Y',[[-1,0,0],[0,1,0],[0,0,-1]]],['180° Z',[[-1,0,0],[0,-1,0],[0,0,1]]]];
 let best=null;
 for(const [nom,R0] of inits){let T={s:1,R:R0,t:[0,1,2].map(i=>c0[i]-(R0[i][0]*c0[0]+R0[i][1]*c0[1]+R0[i][2]*c0[2]))};let res=1e9;
  for(let it=0;it<12;it++){const P=[],Q=[];let acc=0;for(const v of src){const x=[0,1,2].map(i=>T.s*(T.R[i][0]*v[0]+T.R[i][1]*v[1]+T.R[i][2]*v[2])+T.t[i]);p.set(x[0],x[1],x[2]);tg.boundsTree.closestPointToPoint(p,tgt);P.push(v);Q.push([tgt.point.x,tgt.point.y,tgt.point.z]);acc+=p.distanceTo(tgt.point);}
   res=acc/src.length;const U=umeyama(P,Q);T=U;}
  if(!best||res<best.res)best={res,T,nom};}
 const T=best.T;const M=new THREE.Matrix4().set(T.s*T.R[0][0],T.s*T.R[0][1],T.s*T.R[0][2],T.t[0],T.s*T.R[1][0],T.s*T.R[1][1],T.s*T.R[1][2],T.t[1],T.s*T.R[2][0],T.s*T.R[2][1],T.s*T.R[2][2],T.t[2],0,0,0,1);
 // la matrice s'applique aux coordonnées déjà décalées par m.position : on l'intègre
 const P0=new THREE.Matrix4().makeTranslation(m.position.x,m.position.y,m.position.z);const full=new THREE.Matrix4().multiplyMatrices(M,P0);
 m.matrixAutoUpdate=false;m.matrix.copy(full);hullFit[s.id]={matrix:full,residu:best.res,init:best.nom,scale:T.s};
 $('hullInfo').textContent=`calage 3D sur la peau : écart moyen ${best.res.toFixed(3)} mm, échelle ×${T.s.toFixed(2)} (départ ${best.nom})`;}
$('hull').onchange=()=>{showHull();updateShown();};
// tube digestif séparé : PLY dans le même repère (mm, même centre), pas de recentrage
const RECONS={dig:{cb:'dig',lbl:'digLbl',champ:'dig',sg:'Tube digestif (reconstruction séparée)'},cardio:{cb:'cardio',lbl:'cardioLbl',champ:'cardio',sg:'Cœur (détourage séparé)'},vais:{cb:'vais',lbl:'vaisLbl',champ:'vais',sg:'Aortes dorsales (tracé séparé)'},musc:{cb:'musc',lbl:'muscLbl',champ:'musc',sg:'Muscles (reconstruction séparée)',banni:true},topo:{cb:'topo',lbl:'topoLbl',champ:'topo',sg:'Colonne : axe et étages (topographie)'}};
const recCache={},recGroup={};let digGroup=null;
function showRecon(nom){const R=RECONS[nom];const s=ST.find(x=>x.id===cur);if(recGroup[nom]){scene.remove(recGroup[nom]);recGroup[nom]=null;}if(R.banni||!$(R.cb).checked||!s||!s[R.champ].length)return;
 const g=new THREE.Group();g.name='recon:'+nom;recGroup[nom]=g;scene.add(g);
 const inGlb=new Set();if(root)root.traverse(o=>{if(o.isMesh)inGlb.add(o.name);});
 const fused=n=>inGlb.has(n)||inGlb.has('digestif_'+n)||(n==='coeur'&&inGlb.has('coeur_detoure'))||(n.startsWith('aorte')&&inGlb.has('vaisseaux_aorte'))||(n.startsWith('cardinale')&&inGlb.has('vaisseaux_cardinales'))||(n.startsWith('ombilicale')&&(inGlb.has('vaisseaux_ombilicales')||inGlb.has('vaisseaux_ombilicaux')))||(n.startsWith('vitelline')&&(inGlb.has('vaisseaux_vitellines')||inGlb.has('vaisseaux_vitellins')));
 s[R.champ].forEach(seg=>{if(fused(seg.nom))return;const put=geo=>{if(recGroup[nom]!==g)return;const col=new THREE.Color(...seg.color);const m=new THREE.Mesh(geo,new THREE.MeshStandardMaterial({color:col,roughness:0.45,side:THREE.DoubleSide}));m.name=seg.nom;m.userData.recon=nom;m.userData.conf=seg.confiance;g.add(m);applyLayers();buildLayers();};
  if(recCache[seg.url])put(recCache[seg.url]);else new PLYLoader().load(seg.url,geo=>{geo.computeVertexNormals();recCache[seg.url]=geo;put(geo);});});}
function showDig(){Object.keys(RECONS).forEach(showRecon);}
Object.values(RECONS).forEach(R=>{$(R.cb).onchange=()=>{showRecon(R.champ);buildLayers();updateShown();};if(!ST.some(s=>s[R.champ].length))$(R.lbl).style.display='none';});
// ---- couches : état global (par nom de structure), appliqué à chaque stade chargé
const layerOn={};const COLL_DIG='Tube digestif (séparé)';const openSys=new Set(); // systèmes dépliés (tous repliés à l'ouverture)
// brouillons du pipeline automatique remplacés par une reconstruction : jamais affichés ni listés
const BANNIS=new Set(['cavite_pericardique','muscles_membres']);   // brouillons de ma main ; le travail du pipeline s'affiche
const SNC_GROSSIER=new Set(['snc','ventricules','ganglions']);   // 24/09 : rien du système nerveux automatique (snc, ventricules, ganglions) tant qu'il n'est pas segmenté proprement
function isScratch(name){const s=ST.find(x=>x.id===cur||x.id+'|full'===cur);if(BANNIS.has(name))return true;if(!s)return false;
 if(SNC_GROSSIER.has(name)){let fin=false;if(root)root.traverse(o=>{if(o.isMesh&&/^(prosencephale|rhombencephale|moelle)$/.test(o.name))fin=true;});return fin;}   // SNC automatique remplacé par la segmentation fine (CS20)
 if(s.b3&&s.b3.length)return !s.b3.includes(name);const names=new Set();if(root)root.traverse(o=>{if(o.isMesh)names.add(o.name);});
 if(name==='tube_digestif')return s.dig.length>0||[...names].some(n=>n.startsWith('digestif_'));
 if(name==='vaisseaux')return s.vais.length>0||names.has('vaisseaux_aorte');
 if(name==='coeur')return s.cardio.length>0||names.has('coeur_detoure');
 return false;}
var soloFilter=null,playEnd=266;
function structInfo(name){const l=D.lignes.find(x=>x.stade===cur);const st=l&&l.A?l.A.structures.find(x=>x.nom===name):null;return st||{nom:name,collection:'Autres',color:[.7,.7,.7]};}
// code couleur des reconstructions (demande utilisateur) : digestif jaune, veines bleues, artères rouges
const PALETTE=[[/^(dig:|digestif_|oesophage|estomac|duodenum|intestin_)/,0xf2c14e],[/^(vais:aorte|vaisseaux_aorte|aorte)/,0xd62828],
 [/^(vais:(cardinale|ombilicale|vitelline|veine)|vaisseaux_(cardinales|ombilicales|vitellines|ombilicaux|vitellins|veines)|cardinale_|veine_|ombilicale_|vitelline_)/,0x2f6fd6],
 [/^(cardio:coeur|coeur_detoure$|myocarde$)/,0xb0202a],[/^(cardio:cavites|cavites_cardiaques$)/,0xf08a8a],[/^(somites|musc:somites)/,0x9b59b6],[/^notochorde$/,0xf4f1e8],[/^(musc:|muscles_)/,0xa0522d],[/^ventricules$/,0x8fd3ff],[/^snc$/,0x4f7fc9],[/^ganglions$/,0x9fc5e8],[/^prosencephale$/,0x3f6fc4],[/^mesencephale$/,0x5b8be0],[/^rhombencephale$/,0x7fa8f0],[/^ventricule_/,0xbfe3ff],[/^canal_central$/,0xdff2ff],[/^moelle$/,0x2f5fb0],[/^meninges_mesenchyme_cranien$/,0xa9b8c8],[/^epiderme_cranien$/,0xe8c8b8],[/^(topo:axe_vertebral|axe_vertebral$)/,0xff9f43],[/^(topo:etages|etages_vertebraux$|corps_vertebraux_predits$)/,0xf1ead6],[/^squelette_axial_cartilage$/,0xf1ead6],[/^chondrocrane$/,0xe8e0c8],[/^cartilage_autre$/,0xd9d0b8]];
function recolor(o,key){if(o.userData.recol)return;for(const [re,c] of PALETTE){if(re.test(key)){o.material=o.material.clone();o.material.vertexColors=false;o.material.color=new THREE.Color(c);o.material.roughness=0.55;o.material.metalness=0;o.material.needsUpdate=true;break;}}o.userData.recol=true;}
let reconOnly=false;   // défaut : modèle complet
const RECON_RE=/^(somites$|notochorde$|squelette_axial_cartilage$|prosencephale$|mesencephale$|rhombencephale$|ventricule_(pro|mes|rhomb)encephale$|canal_central$|moelle$|meninges_mesenchyme_cranien$|epiderme_cranien$|topo:|axe_vertebral$|corps_vertebraux_predits$|etages_vertebraux$|coeur_detoure|myocarde|cavites_cardiaques|digestif_|vaisseaux_aorte|vaisseaux_veines|vaisseaux_cardinales|vaisseaux_ombilicales|vaisseaux_vitellines|vaisseaux_ombilicaux|vaisseaux_vitellins|cardinale_|veine_|ombilicale_|vitelline_|aorte|oesophage|estomac|duodenum|intestin_|muscles_membres)/;
function applyLayers(){if(root){root.visible=!(window.srcOff&&window.srcOff.has('A'));root.traverse(o=>{if(o.isMesh){if(o.name==='enveloppe'&&$('skinContour')&&$('skinContour').checked){o.visible=true;return;}recolor(o,o.name);if(isScratch(o.name)){o.visible=false;return;}if(soloFilter&&layerOn[o.name]===undefined)layerOn[o.name]=soloFilter(o.name);o.visible=layerOn[o.name]!==false&&(!reconOnly||RECON_RE.test(o.name));}});}Object.entries(recGroup).forEach(([k,g])=>{if(g)g.children.forEach(m=>{const kk=k+':'+m.name;recolor(m,kk);if(soloFilter&&layerOn[kk]===undefined)layerOn[kk]=soloFilter(kk);m.visible=layerOn[kk]!==false;});});}
function buildLayers(){const el=$('layers');if(!root){el.innerHTML='<span class="mut">chargement…</span>';return;}
 // système → sous-groupe → structures ; état = layerOn[key]
 const tree={};D.systemes.forEach(sy=>tree[sy]={});
 const add=(key,nom,color,sgForce)=>{const base=nom.replace(/ \(.*$/,'');const m=D.struct_sys[base];const sy=(m||['Autres'])[0],sg=(m&&m[1])||sgForce||'—';(tree[sy][sg]=tree[sy][sg]||[]).push({key,nom,color});};
 root.traverse(o=>{if(o.isMesh&&!isScratch(o.name))add(o.name,o.name,structInfo(o.name).color);});
 const s=ST.find(x=>x.id===cur);Object.entries(RECONS).forEach(([k,R])=>{if(s&&s[R.champ].length&&$(R.cb).checked)s[R.champ].forEach(x=>add(k+':'+x.nom,x.nom+' ('+x.confiance+')',x.color,R.sg));});
 const hex=c=>'#'+(c||[.7,.7,.7]).map(v=>Math.round(v*255).toString(16).padStart(2,'0')).join('');
 const st=keys=>{const on=keys.filter(k=>layerOn[k]!==false).length;return on===0?'off':on===keys.length?'on':'mix';};
 const cb=(attr,val,keys)=>{const e=st(keys);return `<input type="checkbox" data-${attr}="${val}" ${e==='on'?'checked':''} ${e==='mix'?'data-mix="1"':''}>`;};
 let h='';D.systemes.forEach(sy=>{const sgs=tree[sy];const keysSy=Object.values(sgs).flat().map(x=>x.key);if(!keysSy.length&&sy==='Autres')return;
  h+=`<details data-sys="${sy}" ${openSys.has(sy)?'open':''}><summary>${cb('sys',sy,keysSy)} ${sy} <span class="mut" style="font-weight:400">${keysSy.length}</span></summary>`;
  if(!keysSy.length)h+=`<div class="vide">aucune structure segmentée à ce stade</div>`;
  Object.entries(sgs).forEach(([sg,xs])=>{const ks=xs.map(x=>x.key);h+=`<div class="sg"><label>${cb('sg',sy+'|'+sg,ks)} ${sg}</label>`+xs.map(x=>`<label><input type="checkbox" data-key="${x.key}" ${layerOn[x.key]!==false?'checked':''}><i style="background:${hex(x.color)}"></i>${x.nom}</label>`).join('')+'</div>';});
  h+='</details>';});
 el.innerHTML=h;el.querySelectorAll('input[data-mix]').forEach(i=>i.indeterminate=true);
 // menu horizontal : un bouton par système, un sous-bouton par groupe (clic = bascule, double-clic = ce groupe seul)
 const M=$('sysMenu');const allKeys=Object.values(tree).flatMap(sgs=>Object.values(sgs).flat().map(x=>x.key));
 M.innerHTML=`<div class="hint">Isoler un système ou un groupe pour voir les incohérences : clic = afficher / masquer · double-clic = ce groupe seul · <button data-all="1" style="font-size:11px;padding:1px 6px">tout réafficher</button></div>`+D.systemes.map(sy=>{const sgs=tree[sy];const keysSy=Object.values(sgs).flat().map(x=>x.key);if(!keysSy.length)return '';
  return `<div class="sys"><button data-msys="${sy}" class="${st(keysSy)==='on'?'':st(keysSy)==='off'?'off':'mix'}">${sy} <span class="mut" style="font-weight:400">${keysSy.length}</span></button><div class="sub">`+Object.entries(sgs).map(([sg,xs])=>{const ks=xs.map(x=>x.key);return `<button data-msg="${sy}|${sg}" class="${st(ks)==='on'?'':st(ks)==='off'?'off':'mix'}" title="${xs.map(x=>x.nom).join(', ')}">${sg}</button>`;}).join('')+`</div></div>`;}).join('');
 const solo=keys=>{allKeys.forEach(k=>layerOn[k]=keys.includes(k));applyLayers();buildLayers();};
 M.querySelectorAll('button[data-msys]').forEach(b=>{const ks=Object.values(tree[b.dataset.msys]).flat().map(x=>x.key);b.onclick=()=>setAll(ks,st(ks)!=='on');b.ondblclick=e=>{e.preventDefault();solo(ks);};});
 M.querySelectorAll('button[data-msg]').forEach(b=>{const [sy,sg]=b.dataset.msg.split('|');const ks=tree[sy][sg].map(x=>x.key);b.onclick=()=>setAll(ks,st(ks)!=='on');b.ondblclick=e=>{e.preventDefault();solo(ks);};});
 const ba=M.querySelector('button[data-all]');if(ba)ba.onclick=()=>setAll(allKeys,true);
 el.querySelectorAll('details').forEach(d=>d.addEventListener('toggle',()=>{if(d.open)openSys.add(d.dataset.sys);else openSys.delete(d.dataset.sys);}));
 const setAll=(keys,v)=>{keys.forEach(k=>layerOn[k]=v);applyLayers();buildLayers();};
 el.querySelectorAll('input[data-key]').forEach(i=>i.onchange=()=>setAll([i.dataset.key],i.checked));
 el.querySelectorAll('input[data-sys]').forEach(i=>{i.onclick=e=>e.stopPropagation();i.onchange=()=>setAll(Object.values(tree[i.dataset.sys]).flat().map(x=>x.key),i.checked);});
 el.querySelectorAll('input[data-sg]').forEach(i=>i.onchange=()=>{const [sy,sg]=i.dataset.sg.split('|');setAll(tree[sy][sg].map(x=>x.key),i.checked);});}
window.__v3d={scene,cam,get root(){return root},get box(){return box},get center(){return center},get radius(){return radius},get mode(){return mode},get tPlay(){return tPlay}};
// ---- temps
let playing=false,tPlay=0,day=28,last=performance.now();
try{const sv=JSON.parse(sessionStorage.getItem('v3d')||'null');if(sv&&typeof sv.day==='number'){day=sv.day;}}catch(e){}
setInterval(()=>{try{sessionStorage.setItem('v3d',JSON.stringify({day}));}catch(e){}},2000);
const nearest=d=>ST.reduce((b,s)=>Math.abs(s.mid-d)<Math.abs(b.mid-d)?s:b,ST[0]);
const lastReached=d=>{let r=ST[0];for(const s of ST){if(s.mid<=d)r=s;}return r;};
const stageAt=d=>{const s=lastReached(d);return (d>=ST[0].j0-4&&d<=ST[ST.length-1].j1+1)?s:null;};   // hors de l'intervalle des modèles : rien
const nextOf=s=>{const i=ST.indexOf(s);return i>=0&&i<ST.length-1?ST[i+1]:null;};
const glide=()=>1;   // 24/09 : plus aucune animation de caméra (glissement retiré)
$('play').onclick=()=>{playing=!playing;$('play').textContent=playing?'❚❚ Pause':'▶ Lecture';$('play').classList.toggle('on',playing);if(playing&&(day<dayLo-4||day>=dayHi)){day=Math.max(0,dayLo-4);$('day').value=day;}};
let dayLo=+$('dayLo').value,dayHi=+$('dayHi').value;
function rngDraw(){const w=$('rngBand').parentElement.clientWidth||300;$('rngBand').style.left=(dayLo/DMAX*w)+'px';$('rngBand').style.width=(Math.max(0,dayHi-dayLo)/DMAX*w)+'px';$('rngLbl').textContent=`lecture bornée J${dayLo}–J${dayHi}`;}
function setRange(lo,hi){dayLo=Math.max(0,Math.min(DMAX,Math.round(lo)));dayHi=Math.max(dayLo,Math.min(DMAX,Math.round(hi)));$('dayLo').value=dayLo;$('dayHi').value=dayHi;rngDraw();}
$('dayLo').oninput=()=>{const v=+$('dayLo').value;if(v>dayHi)$('dayLo').value=dayHi;setRange(+$('dayLo').value,dayHi);};
$('dayHi').oninput=()=>{const v=+$('dayHi').value;if(v<dayLo)$('dayHi').value=dayLo;setRange(dayLo,+$('dayHi').value);};
new ResizeObserver(rngDraw).observe($('rngBand').parentElement);setRange(ST[0].j0,ST[ST.length-1].j1);
$('day').oninput=()=>{day=+$('day').value;tPlay=(day-28)/2;};
function setDay(d){playing=false;$('play').textContent='▶ Lecture';$('play').classList.remove('on');day=Math.max(0,Math.min(DMAX,d));$('day').value=day;tPlay=(day-28)/2;}
document.querySelectorAll('#dayBtns button').forEach(b=>b.onclick=()=>setDay(+b.dataset.day));
document.querySelectorAll('#stageBtns button').forEach(b=>b.onclick=()=>{const s=ST.find(x=>x.id===b.dataset.stage);setDay(s.mid);});
// familles de reconstructions : test sur le nom d'une couche (maillage du GLB ou clé recon « xxx:segment »)
const FAM=[
 {id:'coeur',lbl:'Cœur détouré',re:/^(cardio:|coeur_detoure$|myocarde$|cavites_cardiaques$)/},
 {id:'dig',lbl:'Tube digestif',re:/^(dig:|digestif_)/},
 {id:'aortes',lbl:'Aortes dorsales',re:/^(vais:aorte|vaisseaux_aorte$)/},
 {id:'veines',lbl:'Veines',re:/^(vais:(cardinale|ombilicale|vitelline|veine)|vaisseaux_(cardinales|ombilicales|vitellines|ombilicaux|vitellins|veines)$)/},
 {id:'somites',lbl:'Somites / myotomes',re:/^(somites|musc:somites)/},
 {id:'notochorde',lbl:'Notochorde',re:/^notochorde$/},
 {id:'musc',lbl:'Muscles',re:/^(musc:|muscles_)/},
 {id:'colonne',lbl:'Vertèbres, arcs et côtes (axe, étages, corps prédits, cartilage axial fiable)',re:/^(topo:|axe_vertebral$|corps_vertebraux_predits$|etages_vertebraux$|squelette_axial_cartilage$)/},
 {id:'sncfin',lbl:'Système nerveux, segmentation fine (vésicules, ventricules, moelle, méninges)',re:/^(prosencephale|mesencephale|rhombencephale|ventricule_|canal_central|moelle$|meninges_mesenchyme_cranien)/},
];
// clés potentielles d'un stade (avant chargement) : maillages du manifeste + segments des reconstructions séparées
const keysOf=l=>{const a=l.A;if(!a)return [];const k=(a.b3_autorisees?a.structures.filter(x=>a.b3_autorisees.includes(x.nom)):a.structures).map(x=>x.nom);[['digestif','dig'],['cardio','cardio'],['vais','vais'],['musc','musc']].forEach(([f,pfx])=>{if(a[f])a[f].segments.forEach(x=>k.push(pfx+':'+x.nom));});return k;};
FAM.forEach(f=>{const st=D.lignes.filter(l=>l.A&&l.A.vue3d&&l.A.vue3d.glb&&keysOf(l).some(k=>f.re.test(k)));f.stades=st.map(l=>l.stade);f.j0=st.length?Math.min(...st.map(l=>l.jours[0])):null;f.j1=st.length?Math.max(...st.map(l=>l.jours[1])):null;});
$('playRecon').innerHTML=`<div class="t">Familles à lire (chacune sur ses stades)</div><div class="small" style="margin:2px 0 6px"><i style="display:inline-block;width:11px;height:11px;border-radius:3px;background:#f2c14e"></i> digestif · <i style="display:inline-block;width:11px;height:11px;border-radius:3px;background:#d62828"></i> artères · <i style="display:inline-block;width:11px;height:11px;border-radius:3px;background:#2f6fd6"></i> veines · <i style="display:inline-block;width:11px;height:11px;border-radius:3px;background:#b0202a"></i> cœur · <i style="display:inline-block;width:11px;height:11px;border-radius:3px;background:#f08a8a"></i> cavités cardiaques · <i style="display:inline-block;width:11px;height:11px;border-radius:3px;background:#9b59b6"></i> somites · <i style="display:inline-block;width:11px;height:11px;border-radius:3px;background:#f4f1e8"></i> notochorde · <i style="display:inline-block;width:11px;height:11px;border-radius:3px;background:#a0522d"></i> muscles · <i style="display:inline-block;width:11px;height:11px;border-radius:3px;background:#3f6fc4"></i>/<i style="display:inline-block;width:11px;height:11px;border-radius:3px;background:#7fa8f0"></i> vésicules cérébrales · <i style="display:inline-block;width:11px;height:11px;border-radius:3px;background:#bfe3ff"></i> ventricules · <i style="display:inline-block;width:11px;height:11px;border-radius:3px;background:#2f5fb0"></i> moelle · <i style="display:inline-block;width:11px;height:11px;border-radius:3px;background:#f1ead6"></i> vertèbres, arcs, côtes</div><div class="fam">`+FAM.map(f=>`<label class="${f.stades.length?'':'absent'}" title="${f.stades.join(', ')||'aucun stade'}"><input type="checkbox" data-fam="${f.id}" ${f.stades.length?'':'disabled'}><b>${f.lbl}</b><span class="iv">${f.stades.length?`J${f.j0}–J${f.j1} · ${f.stades.length} st.`:'aucun stade'}</span></label>`).join('')+`</div><div class="act"><button id="playSel">▶ Lire la sélection</button><button id="playAll">▶ Lire toutes les reconstructions</button><button id="showSel">◉ Afficher la sélection (sans lecture)</button><button id="playReset">↺ Tout réafficher</button></div>`;
// soloFilter : fonction(clé) → visible ?, appliquée après chaque chargement (déclarés plus haut)
function famSel(){return FAM.filter(f=>{const c=document.querySelector('input[data-fam="'+f.id+'"]');return c&&c.checked;});}
function setSolo(fams){soloFilter=fams.length?(k=>fams.some(f=>f.re.test(k))):null;halfPrev=0;applyLayersSolo();}
function applyLayersSolo(){if(!soloFilter)return;Object.keys(layerOn).forEach(k=>layerOn[k]=soloFilter(k));if(root)root.traverse(o=>{if(o.isMesh)layerOn[o.name]=soloFilter(o.name);});Object.entries(recGroup).forEach(([k,g])=>{if(g)g.children.forEach(m=>{layerOn[k+':'+m.name]=soloFilter(k+':'+m.name);});});applyLayers();}
let famActives=[];function playFams(fams,lecture){if(!fams.length){alert('Cocher au moins une famille');return;}const j0=Math.min(...fams.map(f=>f.j0)),j1=Math.max(...fams.map(f=>f.j1));setRange(j0,j1);famActives=fams;setSolo(fams);if(lecture){setDay(Math.max(0,j0-4));playEnd=j1+1;playing=true;$('play').textContent='❚❚ Pause';$('play').classList.add('on');}else{setDay(Math.min(...fams.map(f=>f.j0))+1);}document.querySelector('.v3d .stage').classList.remove('blink');void document.querySelector('.v3d .stage').offsetWidth;document.querySelector('.v3d .stage').classList.add('blink');}
$('playSel').onclick=()=>playFams(famSel(),true);$('playAll').onclick=()=>playFams(FAM.filter(f=>f.stades.length),true);
// ouverture : toutes les reconstructions cochées, mode « reconstructions seules », lecture lancée
setTimeout(()=>{const d0=day;setRange(ST[0].j0,ST[ST.length-1].j1);playing=true;$('play').textContent='❚❚ Pause';$('play').classList.add('on');if(!(d0>28&&d0<60)){day=Math.max(0,ST[0].j0-4);$('day').value=day;}},300);$('showSel').onclick=()=>playFams(famSel(),false);
$('playReset').onclick=()=>{soloFilter=null;playEnd=266;famActives=[];setRange(ST[0].j0,ST[ST.length-1].j1);document.querySelectorAll('input[data-fam]').forEach(c=>c.checked=false);Object.keys(layerOn).forEach(k=>layerOn[k]=true);applyLayers();buildLayers();};
// depuis la matrice de progression : cellule = ce stade + cette structure seule ; nom = lecture de la structure sur ses stades
const AL=D.alias||{};const baseKey=k=>{const b=k.replace(/^[a-z]+:/,'').replace(/^digestif_/,'');return [b,AL[b]||b];};
window.viewerShow=(stade,nom,lecture)=>{const st=D.lignes.filter(l=>l.A&&keysOf(l).some(k=>baseKey(k).includes(nom)));if(!st.length){alert('Aucune couche « '+nom+' » chargeable');return;}
 soloFilter=k=>baseKey(k).includes(nom);playEnd=266;famActives=[{lbl:nom}];setRange(Math.min(...st.map(l=>l.jours[0])),Math.max(...st.map(l=>l.jours[1])));
 if(lecture){const j0=Math.min(...st.map(l=>l.jours[0])),j1=Math.max(...st.map(l=>l.jours[1]));setDay(Math.max(0,j0-4));playEnd=j1+1;playing=true;$('play').textContent='❚❚ Pause';$('play').classList.add('on');}
 else{const l=st.find(x=>x.stade===stade)||st[0];setDay((l.jours[0]+l.jours[1])/2);}
 applyLayersSolo();document.querySelector('.v3d').scrollIntoView({behavior:'smooth',block:'start'});const sg=document.querySelector('.v3d .stage');sg.classList.remove('blink');void sg.offsetWidth;sg.classList.add('blink');};

$('reconOnly').onclick=()=>{reconOnly=!reconOnly;$('reconOnly').classList.toggle('on',reconOnly);$('reconOnly').textContent=reconOnly?'◉ Tout afficher':'◐ Reconstructions seules';cur=null;applyLayers();updateShown();};
$('snap').onclick=()=>{const a=document.createElement('a');a.download=`agregateur_${cur}_J${day.toFixed(0)}.png`;a.href=canvas.toDataURL('image/png');a.click();};
// ---- caméra
let mode='auto',controls=null,ret=null;
function goalPose(){const s=stageAt(day)||nearest(day);const off=TR[s.id]||{az:0,el:0};let P=pose(off.az,off.el,off.roll);
 const gk=glide();if(gk<1&&TR.depart){const P0=pose((off.az||0)+(TR.depart.az||0),(off.el||0)+(TR.depart.el||0),off.roll);const k=gk;const e=k*k*(3-2*k);
  const q=P0.q.clone().slerp(P.q,e);const m=new THREE.Matrix4().makeRotationFromQuaternion(q);const right=new THREE.Vector3().setFromMatrixColumn(m,0),up=new THREE.Vector3().setFromMatrixColumn(m,1),v=new THREE.Vector3().setFromMatrixColumn(m,2).negate();P={v,up,right,q};}
 // cadrage orthographique sur l'emprise, échelle jamais décroissante en auto
 let ex=1,ey=1;if(root){const vb=new THREE.Box3();let any=false;const tmp=new THREE.Box3();
  const consider=o=>{if(o.isMesh&&o.visible&&o.parent&&o.parent.visible&&o.name!=='hull_voka'&&o.geometry&&o.geometry.boundingBox!==undefined){if(!o.geometry.boundingBox)o.geometry.computeBoundingBox();tmp.copy(o.geometry.boundingBox).applyMatrix4(o.matrixWorld);vb.union(tmp);any=true;}};
  root.traverse(consider);Object.values(recGroup).forEach(g=>{if(g)g.traverse(consider);});
  const bb=any?vb:box;const pts=[];for(let i=0;i<8;i++)pts.push(new THREE.Vector3(i&1?bb.max.x:bb.min.x,i&2?bb.max.y:bb.min.y,i&4?bb.max.z:bb.min.z));
  if(any)bb.getCenter(center);
  const xs=pts.map(p=>p.dot(P.right)),ys=pts.map(p=>p.dot(P.up));ex=Math.max(...xs)-Math.min(...xs);ey=Math.max(...ys)-Math.min(...ys);}
 let half=Math.max(ey/2,ex/2/aspect)*1.15;if(playing&&day>=28){half=Math.max(half,halfPrev);halfPrev=half;}else halfPrev=half;
 return {q:P.q,up:P.up,pos:center.clone().sub(P.v.clone().multiplyScalar(radius*4+5)),zoom:1/half};}
function applyPose(g,k){cam.position.lerp(g.pos,k);cam.quaternion.slerp(g.q,k);cam.up.lerp(g.up,k).normalize();cam.zoom+=(g.zoom-cam.zoom)*k;cam.updateProjectionMatrix();}
function toManual(){if(mode==='manuel')return;mode='manuel';$('mode').textContent='manuel';$('mode').className='mode manuel';$('backAuto').style.display='';
 controls=new OrbitControls(cam,canvas);controls.target.copy(center);controls.enableDamping=true;controls.update();}
canvas.addEventListener('pointerdown',toManual);canvas.addEventListener('wheel',toManual,{passive:true});
$('backAuto').onclick=()=>{if(controls){controls.dispose();controls=null;}mode='retour';$('backAuto').style.display='none';$('mode').textContent='retour…';ret=performance.now();};
function step(now){const dt=Math.max(0,Math.min(0.1,(now-last)/1000));last=Math.max(last,now);
 let buffering=false;
 if(playing){const next=stageAt(day+3);if(next)preload(next);
  const cible=stageAt(day+dt*2);
  if(!cible||loaded[keyOf(cible)]||(keyOf(cible)===cur&&root)){day=Math.min(DMAX,day+dt*2);$('day').value=day;if(day>=Math.min(DMAX,playEnd,dayHi+1)){if($('loop').checked){day=Math.max(0,dayLo-4);$('day').value=day;halfPrev=0;}else{playing=false;playEnd=266;$('play').textContent='▶ Lecture';$('play').classList.remove('on');}}}
  else{buffering=true;preload(cible);}}
 tPlay=(day-28)/2;const gk=glide();
 const sAt=stageAt(day);const lig=D.lignes.slice().reverse().find(x=>day>=x.frise[0])||D.lignes[0];
 $('dayLbl').innerHTML=`${playing?'▶ Lecture':'⏸ Arrêt'} · J${day.toFixed(1)} · ${sAt?sAt.id:lig.stade+' (pas de modèle)'}`+(gk>0&&gk<1?` · glissement caméra ${Math.round(gk*100)} %`:'')+(buffering?' · <span class="buf">chargement du stade suivant…</span>':'');
 let nvis=0;V3SCENE.traverse(o=>{if(o.isMesh&&o.visible&&o.parent&&o.parent.visible)nvis++;});
 const sA=stageAt(day),nA=sA?nextOf(sA):null;const interp=(MORPH_ACTIF&&sA&&nA&&day>sA.mid)?` · interpolation ${sA.id}→${nA.id} ${Math.round(100*Math.max(0,Math.min(1,(day-sA.mid)/Math.max(0.01,nA.mid-sA.mid))))} %${morphQueue.length?' (calcul…)':''}`:'';
 $('statusMode').innerHTML=(famActives.length?`Lecture : <b>${famActives.map(f=>f.lbl).join(' + ')}</b> · `:'')+interp.replace(/^ · /,'')+(interp?' · ':'')+(reconOnly?'<b>reconstructions seules</b> · ':'')+(mode==='manuel'?'caméra manuelle · ':'caméra auto · ')+`bornes J${dayLo}–J${dayHi} · ${nvis} couche${nvis>1?'s':''} affichée${nvis>1?'s':''}`+'';
 $('pfill').style.width=(Math.max(0,Math.min(1,(day-dayLo)/Math.max(1,dayHi-dayLo)))*100)+'%';
 if(window.setDayCursor)window.setDayCursor(day);document.querySelectorAll('#dayBtns button').forEach(b=>b.classList.toggle('on',+b.dataset.day<=day&&(b.nextElementSibling?+b.nextElementSibling.dataset.day>day:true)));
 showStage(stageAt(day));if(root)root.visible=!(window.srcOff&&window.srcOff.has('A'));
 {const s_=stageAt(day);const n_=s_?nextOf(s_):null;if(n_)preload(n_);morphStep();applyMorph();}
 if(mode==='auto'){const g=goalPose();applyPose(g,1);}
 else if(mode==='retour'){const g=goalPose();applyPose(g,1);mode='auto';$('mode').textContent='automatique';$('mode').className='mode auto';}
 else if(controls)controls.update();
 renderer.render(scene,cam);}
(function loop(now){requestAnimationFrame(loop);step(now);})(performance.now());
window.__v3d.step=step;window.__v3d.dbg=()=>({playing,day:+day.toFixed(2),playEnd,dayLo,dayHi,fam:famActives.length,cur,loaded:Object.keys(loaded),loading:Object.keys(loading).filter(k=>loading[k]),mq:morphQueue.length,md:morphDone.size,bvh:!!MeshBVH,mq0:morphQueue[0]?{key:morphQueue[0].key,i:morphQueue[0].i||0,n:morphQueue[0].m.geometry.attributes.position.count,to:morphQueue[0].next}:null});
}catch(e){msg('three.js non chargé (CDN inaccessible ?) : '+e.message);}
</script>
</body></html>
"""


# structures jamais affichées (brouillons du pipeline, système nerveux automatique) : absentes des GLB allégés
BANNIS_VUE = {"snc", "ventricules", "ganglions", "cavite_pericardique", "coeur", "tube_digestif", "vaisseaux", "muscles_membres",
              "squelette_axial_cartilage", "chondrocrane", "cartilage_autre"}   # cartilages automatiques en morceaux : « vire tes horreurs » (24/09)
DECIMER = {"enveloppe": 60000, "membres": 30000, "cordon_ombilical": 15000}
DECIMER_DEFAUT = 80000    # plafond de faces pour toute pièce du GLB allégé (SNC fin exporté à 200-250 k faces par pièce)


def glb_allege(a, cs, dossier_out, dest_glb, autorisees=None):
    """Assemble un GLB avec les seules structures affichables (PLY du stade), l'enveloppe décimée ; renvoie le chemin ou None."""
    import trimesh
    import numpy as np
    try:
        import fast_simplification
    except Exception:
        fast_simplification = None
    src_files = []
    cart_ok = conf_de(cs, next((x for x in a["structures"] if x["nom"] == "squelette_axial_cartilage"), None)) in ("bonne", "moyenne")
    for x in a["structures"]:
        n = x["nom"]
        if (n in BANNIS_VUE and not (n == "squelette_axial_cartilage" and cart_ok)) or (autorisees is not None and n not in autorisees):
            continue
        if x.get("confiance") == "faible":      # déclaré faible par la session base dans le manifest : tableau seulement, hors vue 3D
            continue
        f = os.path.join(dossier_out, f"{cs}_{n}.ply")
        if os.path.exists(f):
            src_files.append((n, f, x.get("color", [0.7, 0.7, 0.7])))
    if not src_files:
        return None
    latest = max(os.path.getmtime(f) for _, f, _ in src_files)
    if os.path.exists(dest_glb) and os.path.getmtime(dest_glb) >= latest:
        return dest_glb
    sc = trimesh.Scene()
    for n, f, col in src_files:
        m = trimesh.load(f, force="mesh")
        cible = DECIMER.get(n, DECIMER_DEFAUT)
        if cible and fast_simplification is not None and len(m.faces) > cible * 1.1:
            v, fc = fast_simplification.simplify(m.vertices.astype(np.float32), m.faces.astype(np.int32),
                                                 target_reduction=float(1 - cible / len(m.faces)))
            m = trimesh.Trimesh(v, fc, process=True)
        m.visual = trimesh.visual.ColorVisuals(m, face_colors=np.tile((np.array(list(col) + [1.0]) * 255).astype(np.uint8), (len(m.faces), 1)))
        sc.add_geometry(m, node_name=n, geom_name=n)
    sc.export(dest_glb)
    return dest_glb


def paquet_site(data, dest):
    """Copie les ressources dans dest/ et réécrit les chemins ; les hulls VOKA (privés) sont exclus."""
    for sub in ("glb", "rendus", "overlays", "voka", "hull", "digestif", "cardio", "vais", "musc", "topo", "rapports"):
        os.makedirs(os.path.join(dest, sub), exist_ok=True)

    def cp(src, sub, prefixe=""):
        if not src or not os.path.exists(src):
            return None
        nom = os.path.basename(src)
        if prefixe and not nom.startswith(prefixe):
            nom = prefixe + nom
        dst = os.path.join(dest, sub, nom)
        if (not os.path.exists(dst) or os.path.getsize(dst) != os.path.getsize(src)
                or int(os.path.getmtime(src)) > int(os.path.getmtime(dst))):
            shutil.copy2(src, dst)
        return f"{sub}/{nom}"

    for l in data["lignes"]:
        if l["A"] and l["A"].get("origine") == "B" and l["stade"] not in data["A"]:
            data["A"][l["stade"]] = l["A"]
    entrees = list(data["A"].items())
    vus = {id(a) for _, a in entrees}
    for l in data["lignes"]:
        if l.get("A") and id(l["A"]) not in vus:
            entrees.append((l["stade"], l["A"])); vus.add(id(l["A"]))
    for cs, a in entrees:
        v = a.get("vue3d") or {}
        v["glb_full"] = cp(os.path.join(a["dossier"], "out", f"{cs}_all.glb"), "glb") if a["glb"] else None
        allege = glb_allege(a, cs, os.path.join(a["dossier"], "out"), os.path.join(dest, "glb", f"{cs}_recon.glb"),
                            autorisees=a.get("b3_autorisees") if a.get("origine") == "B" else None) if a["glb"] else None
        v["glb"] = f"glb/{cs}_recon.glb" if allege else v["glb_full"]
        if a.get("controle_B"):
            if a["controle_B"].get("rapport"):
                a["controle_B"]["rapport"] = cp(a["controle_B"]["rapport"].replace("../", "", 1), "rapports", cs + "_")
            a["controle_B"]["images"] = [x for x in (cp(im.replace("../", "", 1), "rapports", cs + "_") for im in a["controle_B"].get("images", [])) if x]
            a["controle_B"]["rapports"] = [x for x in (cp(r.replace("../", "", 1), "rapports", cs + "_") for r in a["controle_B"].get("rapports", [])) if x]
        v["hull"] = cp(v["hull"].replace("../", "", 1), "hull") if v.get("hull") else None
        v["voka_img"] = cp(v["voka_img"].replace("../", "", 1), "voka") if v.get("voka_img") else None
        a["vue3d"] = v
        for nom in ("digestif", "cardio", "vais", "musc", "topo"):
            if a.get(nom):
                for x in a[nom]["segments"]:
                    x["url"] = cp(x["url"].replace("../", "", 1), nom) or x["url"]
        if a.get("topo") and a["topo"].get("controle"):
            a["topo"]["controle"]["images_url"] = [x for x in (cp(u.replace("../", "", 1), "rapports", cs + "_") for u in a["topo"]["controle"].get("images_url", [])) if x]
            # sauvegarde versionnée des JSON de topographie dans le paquet (incident CS20 du 24/09 : seuls les PLY publiés avaient survécu)
            for jn in ("manifest_topographie.json", "niveaux_imposes.json", "axe_vertebral.json", "etages_vertebraux.json", "zones_axe.json"):
                cp(os.path.join(a["dossier"], "out", "topographie", jn), "topo", cs + "_")
        a["rendu_organes_url"] = cp(os.path.join("embryons_3D", "rendus", f"{cs}_organes.png"), "rendus")
        a["rendu_peau_url"] = cp(os.path.join("embryons_3D", "rendus", f"{cs}_peau.png"), "rendus")
        a["rendu_organes"] = bool(a["rendu_organes_url"]); a["rendu_peau"] = bool(a["rendu_peau_url"])
    for b in data["B"].values():
        b["overlay_url"] = cp(b.get("exemple_overlay"), "overlays")
    for c in data["C"].values():
        c["image_url"] = cp(os.path.join(REF, c["image"]), "voka") if c.get("image_ok") else None
    data["angle_img"] = cp(os.path.join(REF, "voka_pages", "angle_camera_des_CS13.webp"), "voka")
    fe = os.path.join("embryons_3D", "taches_etat.json")
    e = json.load(open(fe, encoding="utf-8")) if os.path.exists(fe) else {}
    e["_publie"] = date.today().strftime("%d/%m/%Y") + " (génération)"
    with open(os.path.join(dest, "taches_etat.json"), "w", encoding="utf-8") as f:
        json.dump(e, f, ensure_ascii=False, indent=1)
    data["site"] = True
    return data


def main():
    A, B, C = source_A(), source_B(), source_C()
    lignes = fusion(A, B, C)
    todo = reste_a_faire(A, B, C, lignes)
    data = {"date": date.today().strftime("%d/%m/%Y"), "A": A, "B": B, "C": C, "lignes": lignes, "todo": todo,
            "systemes": SYSTEMES, "struct_sys": STRUCT_SYS, "progression": progression(A), "alias": ALIAS_PROG,
            "cartilage_ok": sorted(cs for cs, a in A.items() if conf_de(cs, next((x for x in a["structures"] if x["nom"] == "squelette_axial_cartilage"), None)) in ("bonne", "moyenne")),
            "angle_img": "../embryo3d/reference/voka_pages/angle_camera_des_CS13.webp"}
    os.makedirs("embryons_3D", exist_ok=True)
    if "--site" in sys.argv:
        dest = sys.argv[sys.argv.index("--site") + 1]
        data = paquet_site(data, dest)
        data["version"] = datetime.now().strftime("%Y%m%d-%H%M%S")
        with open(os.path.join(dest, "version.txt"), "w") as f:
            f.write(data["version"])
        js = json.dumps(data, ensure_ascii=False).replace("</", "<\\/")
        with open(os.path.join(dest, "index.html"), "w", encoding="utf-8", newline="") as f:
            f.write(HTML.replace("__DATA__", js).replace("__ROBOTS__", '<meta name="robots" content="noindex">').replace("__VIEWER3DH__", "viewer.html"))
        tot = sum(os.path.getsize(os.path.join(r, f)) for r, _, fs in os.walk(dest) for f in fs)
        print(f"paquet site : {dest} ({tot/1e6:.0f} Mo) — publier index.html ET version.txt (rafraîchissement auto des pages ouvertes)")
        return
    js = json.dumps(data, ensure_ascii=False).replace("</", "<\\/")
    data["version"] = datetime.now().strftime("%Y%m%d-%H%M%S")
    os.makedirs(os.path.join("embryons_3D", "glb_allege"), exist_ok=True)
    for l in data["lignes"]:
        a = l.get("A")
        if a and a.get("vue3d") and a["vue3d"].get("glb"):
            cs = l["stade"]
            a["vue3d"]["glb_full"] = a["vue3d"]["glb"]
            dst = os.path.join("embryons_3D", "glb_allege", f"{cs}_recon.glb")
            if glb_allege(a, cs, os.path.join(a["dossier"], "out"), dst, autorisees=a.get("b3_autorisees") if a.get("origine") == "B" else None):
                a["vue3d"]["glb"] = f"glb_allege/{cs}_recon.glb"
    js = json.dumps(data, ensure_ascii=False).replace("</", "<\\/")
    with open(SORTIE_HTML, "w", encoding="utf-8") as f:
        f.write(HTML.replace("__DATA__", js).replace("__ROBOTS__", "").replace("__VIEWER3DH__", "site_3dh/index_local.html"))
    with open(os.path.join("embryons_3D", "version.txt"), "w") as f:
        f.write(data["version"])
    with open(SORTIE_JSON, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=1)
    n_h = sum(t["prio"] == "haute" for t in todo)
    print(f"OK {SORTIE_HTML} : {len(A)} stades 3D, {len(B)} stades ehd, {len(C)} VOKA, {len(todo)} tâches ({n_h} hautes)")
    for l in lignes:
        print(f"  {l['stade']:5s} J{l['jours'][0]:>2}-{l['jours'][1]:<2}  {l['statut']}")


if __name__ == "__main__":
    main()
