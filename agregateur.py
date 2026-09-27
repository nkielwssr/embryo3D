# -*- coding: utf-8 -*-
"""Agrégateur « nouveaux modèles 3D à partir de CS10 » (remis à zéro le 25/09/2026).

    python embryo3d/agregateur.py                       # embryons_3D/agregateur.html + agregateur.json (tâches, lu par lanceur.py)
    python embryo3d/agregateur.py --site DIR            # paquet autonome à publier (harcelon.fr/3dht) : index.html, modeles/, version.txt

Le travail du 24/09 (reconstruction CS13-CS20 depuis les vidéos, coupes, colonne, morphing…) est gardé en référence :
code dans embryo3d/reference_2409/, page dans embryons_3D/reference_2409/, en ligne sur harcelon.fr/3dht/reference-2409/.

CONVENTION DES NOUVEAUX MODÈLES (une livraison = un dossier par stade) :
    embryons_3D/modeles/<CS>/manifest.json
    {
      "stade": "CS10", "source": "texte libre (coupes ehd, vidéo 360°, …)", "session": "nom de la session qui livre",
      "unites": "mm", "repere": "Z crânial, Y dorsal, X = Y×Z = GAUCHE anatomique (repère direct, sans miroir), origine = centre du corps",
      "structures": [
        {"nom": "tube_neural", "systeme": "nerveux", "fichier": "tube_neural.ply", "couleur": [0.25, 0.45, 0.8],
         "confiance": "bonne|moyenne|faible", "volume_mm3": 0.12, "note": "…"}
      ],
      "notes": "…"
    }
    systeme ∈ digestif, nerveux, vasculaire, os_cartilages, muscles, derme, appendiculaire, autre.
    fichier = PLY ou GLB dans le même dossier. confiance « faible » : listée, pas affichée en 3D.
"""
import hashlib
import json
import os
import shutil
import sys
from datetime import datetime

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(RACINE)
MODELES = os.path.join("embryons_3D", "modeles")
SORTIE_HTML = os.path.join("embryons_3D", "agregateur.html")
SORTIE_JSON = os.path.join("embryons_3D", "agregateur.json")
ETAT = os.path.join("embryons_3D", "taches_etat.json")
COUPES = "coupes embryos 9-23"
REF_URL = "reference-2409/"          # relatif à la page publiée

# stade, jours min, jours max, CRL min, CRL max (mm), repères (O'Rahilly & Müller 2010)
AXE = [
    ("CS7", 15, 17, 0.4, 0.4, "Notochorde, début de la gastrulation"),
    ("CS8", 17, 19, 1.0, 1.5, "Gouttière neurale, fossette primitive"),
    ("CS9", 19, 21, 1.5, 2.5, "1-3 paires de somites, plis neuraux"),
    ("CS10", 22, 23, 2.0, 3.5, "4-12 somites, fermeture neurale"),
    ("CS11", 23, 26, 2.5, 4.5, "13-20 somites, neuropore rostral fermé"),
    ("CS12", 26, 30, 3.0, 5.0, "21-29 somites, bourgeons des membres supérieurs"),
    ("CS13", 28, 32, 4.0, 6.0, "4 arcs, otocyste fermé, bourgeons des membres inférieurs"),
    ("CS14", 31, 35, 5.0, 7.0, "Cupule optique, fosse nasale"),
    ("CS15", 35, 38, 7.0, 9.0, "Palette de la main, cristallin"),
    ("CS16", 37, 42, 8.0, 11.0, "Pigment rétinien, palette du pied"),
    ("CS17", 42, 44, 11.0, 14.0, "Rayons digitaux, bourrelets auriculaires"),
    ("CS18", 44, 48, 13.0, 17.0, "Début d'ossification, coudes"),
    ("CS19", 48, 51, 16.0, 18.0, "Tronc s'allonge et se redresse"),
    ("CS20", 51, 53, 18.0, 22.0, "Doigts séparés, membres supérieurs fléchis"),
    ("CS21", 53, 54, 22.0, 24.0, "Doigts allongés, les mains se rejoignent"),
    ("CS22", 54, 56, 23.0, 28.0, "Paupières, pavillon de l'oreille"),
    ("CS23", 56, 60, 27.0, 31.0, "Fin de la période embryonnaire, queue disparue"),
]
# hiérarchie voulue par l'utilisateur (couleurs : digestif jaune, artères rouge, veines bleu, cerveau en bleus)
SYSTEMES = [
    ("nerveux", "Nerveux", "#4f7fc9"),
    ("digestif", "Digestif", "#e0b43a"),
    ("vasculaire", "Vasculaire", "#c8323a"),
    ("os_cartilages", "Os + cartilages", "#d9ceb0"),
    ("muscles", "Muscles", "#9b59b6"),
    ("derme", "Derme (peau + enveloppes)", "#e3a98f"),
    ("appendiculaire", "Appendiculaire", "#c98f6e"),
]
# systèmes supplémentaires : colonnes de la matrice (organes) et catégories de vue seulement (annexes, tissus, autre)
SYSTEMES += [
    ("respiratoire", "Respiratoire", "#6cc3c1"),
    ("urogenital", "Uro-génital", "#b98bd6"),
    ("endocrine", "Glandes endocrines", "#e39b3c"),
]
VUE_SEULE = [
    ("annexes", "Annexes et cavités", "#aeb6bf"),
    ("tissus", "Tissus embryonnaires (mésoderme, mésenchyme)", "#8a9098"),
    ("autre", "Autres (oreille moyenne, sillons pharyngés)", "#999999"),
]
SYS_IDS = [s[0] for s in SYSTEMES]
TOUS_IDS = SYS_IDS + [s[0] for s in VUE_SEULE]
RECLASSEMENT = [   # (motif sur le nom, système) — première correspondance
    (r"^(skin|peau|nipples|ectoderme_epidermique|epidermal_ectoderm)", "derme"),
    (r"(amnion|amnios|amniotic|yolk_sac|vesicule_vitelline|allanto|connecting_stalk|pedicule_fixation|exocoelom|extra_embryonic|trophoblast)", "annexes"),
    (r"(pericardial_cavity|cavite_pericardique|pericardioperitoneal|pericardo_peritoneaux|peritoneal_cavity|pleural_cavity|intra_embryonic_coelom|vesicules_coelomiques)", "annexes"),
    (r"(^mesoderm|^mesoderme|intermediate_and_lateral_plate|^epiblast|^hypoblast|primitive_node|primitive_streak)", "tissus"),
    (r"(respiratory|lung|trachea|bronch|pneumato_enteric)", "respiratoire"),
    (r"(mesonephr|metanephr|ureter|renal_pelvis|urinary_bladder|urethr|urogeni|urachus|gonad|genital|paramesonephric)", "urogenital"),
    (r"(thyroid_gland|thyroid_diverticulum|parathyroid|thymus|adrenal|adenohypophysis)", "endocrine"),
    (r"gut_surrounding_mesenchyme", "digestif"),
]


def reclasser(nom, systeme):
    import re
    for motif, sid in RECLASSEMENT:
        if re.search(motif, nom):
            return sid
    return systeme if systeme in TOUS_IDS else "autre"
# dossiers du travail du 24/09 (référence)
REF_VIDEO = {"CS13": "CS13.f4v", "CS14": "CS14_f4v", "CS15": "CS15_f4v", "CS16": "CS16_f4v", "CS17": "CS17_f4v", "CS19": "CS19_f4v", "CS20": "CS20_F4V"}


def compter(dossier, ext):
    if not os.path.isdir(dossier):
        return 0
    return sum(1 for f in os.listdir(dossier) if f.lower().endswith(ext))


def sources(cs):
    n = cs[2:]
    d = os.path.join(COUPES, f"embryo_stage_{n}")
    coupes = compter(os.path.join(d, "images_fond"), (".jpg", ".png"))
    legendes = len(os.listdir(os.path.join(d, "legendes"))) if os.path.isdir(os.path.join(d, "legendes")) else 0
    videos = compter(REF_VIDEO[cs], ".mp4") if cs in REF_VIDEO else 0
    jours = next(a for a in AXE if a[0] == cs)[1]
    ref = cs in REF_VIDEO or cs in ("CS21", "CS22", "CS23")
    return {"coupes": coupes, "legendes": legendes, "videos360": videos, "voka": jours >= 28, "reference_2409": ref}


CACHE_ORIENT = os.path.join("embryons_3D", "orientation_cache.json")
RE_NERF = r"(^spinal_cord$|neural_tube|tube_neural|moelle)"            # tube neural / moelle : dorsal à tous les stades
RE_COEUR = r"(myocard|tube_cardiaque|heart_tube|coeur|^heart$|^(left|right)_(ventricle|atrium)$|^heart_lumen$)"


def controle_orientation(dossier, structures):
    """{'rot_z_180': bool, 'dy_coeur_nerf_mm': float|None} : cœur dorsal au tube neural ⇒ demi-tour autour de Z."""
    import re
    f = os.path.join(MODELES, dossier, "manifest.json")
    cle = f"{dossier}|{int(os.path.getmtime(f))}"
    try:
        cache = json.load(open(CACHE_ORIENT, encoding="utf-8"))
    except Exception:
        cache = {}
    if cle in cache:
        return cache[cle]
    import numpy as np
    import trimesh

    def centre(motif):
        pts = []
        for x in structures:
            if x["ok"] and re.search(motif, x["nom"]) and not re.search("lumen|jelly|gelee|lumiere", x["nom"]):
                try:
                    g = trimesh.load(os.path.join(MODELES, dossier, x["fichier"]), force="mesh")
                    pts.append(g.vertices.mean(0))
                except Exception:
                    pass
        return np.mean(pts, 0) if pts else None
    n, c = centre(RE_NERF), centre(RE_COEUR)
    dy = None if n is None or c is None else round(float(c[1] - n[1]), 3)
    vg, vd = centre(r"^(left_ventricle|myocarde_ventricule_gauche|cavite_ventricule_gauche|coeur_ventricule_gauche)$"), centre(r"^(right_ventricle|myocarde_ventricule_droit|cavite_ventricule_droit|coeur_ventricule_droit)$")
    dx = None if vg is None or vd is None else round(float(vg[0] - vd[0]), 3)     # attendu > 0 : +X = gauche anatomique
    rot = bool(dy is not None and dy > 0.2)                                          # marge 0,2 mm contre les cas limites
    res = {"rot_z_180": rot, "dy_coeur_nerf_mm": dy, "dx_vg_vd_mm": dx,
           "miroir_suspect": bool(dx is not None and ((dx < 0) != rot) and abs(dx) > 0.1)}
    cache = {k: v for k, v in cache.items() if not k.startswith(dossier + "|")}
    cache[cle] = res
    json.dump(cache, open(CACHE_ORIENT, "w", encoding="utf-8"), indent=1)
    return res


def norm(cs):
    """CS07 ≡ CS7."""
    try:
        return int("".join(ch for ch in cs if ch.isdigit()))
    except ValueError:
        return cs


def dossiers_modele(cs):
    """Dossiers de modeles/ dont le manifest porte ce stade (CS10, CS10_recon…), publiables d'abord."""
    out = []
    if os.path.isdir(MODELES):
        for n in sorted(os.listdir(MODELES)):
            f = os.path.join(MODELES, n, "manifest.json")
            if os.path.exists(f) and norm(n.split("_")[0]) == norm(cs):
                try:
                    m = json.load(open(f, encoding="utf-8"))
                except Exception:
                    continue
                if norm(m.get("stade", cs)) == norm(cs):
                    out.append((m.get("publiable", True) is not False, m.get("statut") != "brouillon", n))
    # publiable d'abord, puis reconstruit (non brouillon) avant brouillon
    return [n for _, _, n in sorted(out, key=lambda t: (not t[0], not t[1], t[2]))]


def lire_modele(cs, publiable_seul=False):
    for n in dossiers_modele(cs):
        mod = lire_dossier(cs, n)
        if mod and (mod["publiable"] or not publiable_seul):
            return mod
    return None


def coupes_ams(c):
    """Repérage dans les coupes légendées de l'atlas d'Amsterdam (champ coupes_amsterdam écrit par « coupes amsterdam/integrer_manifests.py ») :
    noms et positions seulement, pas d'images (licence ND)."""
    if not c or not c.get("coupes"):
        return None
    return {"l": c.get("label_fr") or c.get("label"), "g": c.get("methode") == "géométrie", "a": c["coupes"][0], "b": c["coupes"][-1],
            "c": c.get("coupe_max"), "d": c.get("accord_geometrie") is False}


def lire_dossier(cs, nom_dossier):
    d = os.path.join(MODELES, nom_dossier)
    f = os.path.join(d, "manifest.json")
    if not os.path.exists(f):
        return None
    m = json.load(open(f, encoding="utf-8"))
    fv = os.path.join(d, "verif_maillages.json")
    verif = json.load(open(fv, encoding="utf-8")) if os.path.exists(fv) else None
    vs = (verif or {}).get("structures", {})
    structs = []
    for x in m.get("structures", []):
        fich = os.path.join(d, x.get("fichier", ""))
        structs.append({"nom": x.get("nom", "?"), "nom_fr": x.get("nom_fr") or x.get("nom", "?").replace("_", " "), "legende": x.get("legende", ""),
                        "dimensions_um": x.get("dimensions_um"), "couleur_atlas": x.get("couleur_atlas"),
                        "systeme": reclasser(x.get("nom", ""), x.get("systeme", "autre")), "systeme_source": x.get("systeme"),
                        "fichier": x.get("fichier"), "ok": os.path.exists(fich), "couleur": x.get("couleur"),
                        "confiance": x.get("confiance", "?"), "volume_mm3": x.get("volume_mm3"), "note": x.get("note", ""),
                        "verif": (vs.get(x.get("nom")) or {}).get("verdict"), "verif_pb": (vs.get(x.get("nom")) or {}).get("problemes", []),
                        "ca": coupes_ams(x.get("coupes_amsterdam"))})
    return {"source": m.get("source", ""), "session": m.get("session", ""), "notes": m.get("notes", ""), "structures": structs,
            "dossier": nom_dossier, "publiable": m.get("publiable", True) is not False, "usage": m.get("usage", ""),
            "brouillon": m.get("statut") == "brouillon", "externe": m.get("statut") == "externe", "attribution": m.get("attribution", ""),
            "echelle": m.get("echelle") or {}, "orientation_auto": m.get("orientation") or {},
            "verif": (verif or {}).get("bilan"), "coupes_ams": m.get("coupes_amsterdam"),
            "orientation": controle_orientation(nom_dossier, structs),
            "specimen": m.get("specimen", ""), "licence": m.get("licence", ""), "temps": m.get("temps") or {}, "dimensions": m.get("dimensions") or {},
            "date": datetime.fromtimestamp(os.path.getmtime(f)).strftime("%d/%m %H:%M")}


NON_FORME = ({(c, "appendiculaire") for c in ("CS7", "CS8", "CS9", "CS10", "CS11")} | {(c, s) for c in ("CS7", "CS8") for s in ("muscles", "digestif")}
             | {(c, s) for c in ("CS7", "CS8", "CS9", "CS10", "CS11") for s in ("respiratoire", "urogenital")}   # bourgeon pulmonaire, mésonéphros : CS12
             | {(c, "endocrine") for c in ("CS7", "CS8", "CS9", "CS10")})                                     # diverticule thyroïdien : CS11   # bourgeons des membres à CS12


def statut_systeme(mod, sid, cs=""):
    if (cs, sid) in NON_FORME:
        return "na"
    if not mod:
        return "a_faire"
    xs = [x for x in mod["structures"] if x["systeme"] == sid and x["ok"]]
    if not xs:
        return "a_faire"
    if any(x["confiance"] in ("bonne", "moyenne") for x in xs):
        if mod.get("brouillon"):
            return "brouillon"
        return "externe" if mod.get("externe") else "fait"
    return "partiel"


def tid(*parts):
    return hashlib.md5("|".join(parts).encode("utf-8")).hexdigest()[:8]


def taches(lignes):
    T = []
    def add(prio, stade, action, detail=""):
        T.append({"id": tid("v2", stade, action), "prio": prio, "source": "nouveaux modèles", "stade": stade, "action": action, "detail": detail, "cmd": ""})
    add("haute", "tous", "Définir la méthode des nouveaux modèles (source par stade, repère, échelle)",
        "convention de livraison : embryons_3D/modeles/CSxx/manifest.json (voir en-tête d'embryo3d/agregateur.py)")
    for l in lignes:
        manque = [dict((s[0], s[1]) for s in SYSTEMES)[k] for k, v in l["systemes"].items() if v not in ("fait", "externe", "na")]
        if not l["modele"]:
            src = l["sources"]
            dispo = ", ".join(x for x, ok in (("coupes ehd", src["coupes"]), ("vidéo 360°", src["videos360"]), ("VOKA (regard)", src["voka"])) if ok) or "aucune source locale"
            add("haute" if l["stade"] in ("CS10", "CS11", "CS12") else "moyenne", l["stade"], "Nouveau modèle 3D", f"sources : {dispo}")
        elif l["modele"].get("brouillon"):
            add("haute", l["stade"], "Reconstruire par notre code (remplacer le brouillon tiré de l'atlas)",
                "livraison dans embryons_3D/modeles/" + l["stade"] + "_recon/ (statut ≠ brouillon)")
        elif manque:
            add("moyenne", l["stade"], "Compléter le modèle : " + ", ".join(manque))
    return T


RE_VIT = r"^(yolk_sac_wall|yolk_sac|vesicule_vitelline|yolk_sac_cavity)$"
CACHE_VIT = os.path.join("embryons_3D", "vitellus_cache.json")


def vitellus_regulier(lignes):
    """Ajoute modele['vitellus'] = {centre, rayon, source} aux stades CS8-CS12."""
    import re
    import numpy as np
    import trimesh
    try:
        cache = json.load(open(CACHE_VIT, encoding="utf-8"))
    except Exception:
        cache = {}

    def charge(mod, motif):
        out = []
        for x in mod["structures"]:
            if x["ok"] and re.search(motif, x["nom"]):
                try:
                    out.append(trimesh.load(os.path.join(MODELES, mod["dossier"], x["fichier"]), force="mesh"))
                except Exception:
                    pass
        return out

    mesures = {}
    for l in lignes:
        mod, n = l["modele"], norm(l["stade"])
        if not mod or not 8 <= n <= 11:
            continue
        cle = f'{mod["dossier"]}|{int(os.path.getmtime(os.path.join(MODELES, mod["dossier"], "manifest.json")))}'
        if cle not in cache:
            noms = {x["nom"] for x in mod["structures"]}
            motif = RE_VIT
            parois = [x["nom"] for x in mod["structures"] if re.search(RE_VIT, x["nom"]) and not (x["nom"].endswith("_cavity") and x["nom"].replace("_cavity", "_wall") in noms)]
            ms = charge(mod, "^(" + "|".join(map(re.escape, parois)) + ")$") if parois else []
            if not ms:
                cache[cle] = None
            else:
                P = np.vstack([m.vertices for m in ms])
                A = np.c_[2 * P, np.ones(len(P))]; b = (P ** 2).sum(1)
                x, *_ = np.linalg.lstsq(A, b, rcond=None); c = x[:3]; r = float(np.sqrt(x[3] + c @ c))
                mid = charge(mod, r"(midgut|intestin_moyen)")
                cache[cle] = {"centre": [round(float(v), 4) for v in c], "rayon_mesure": round(r, 4),
                              "intestin": [round(float(v), 4) for v in mid[0].vertices.mean(0)] if mid else None}
        if cache[cle]:
            mesures[n] = cache[cle]
    json.dump(cache, open(CACHE_VIT, "w", encoding="utf-8"), indent=1)
    if len(mesures) < 2:
        return
    xs = np.array(sorted(mesures)); ys = np.array([mesures[k]["rayon_mesure"] for k in xs])
    pente, orig = np.polyfit(xs, ys, 1)
    rayon = lambda n: round(float(orig + pente * n), 4)
    # distance centre ↔ intestin moyen, rapportée au rayon, aux stades mesurés
    rel = [np.linalg.norm(np.array(v["centre"]) - np.array(v["intestin"])) / v["rayon_mesure"] for v in mesures.values() if v["intestin"]]
    k = float(np.median(rel)) if rel else 1.3
    for l in lignes:
        mod, n = l["modele"], norm(l["stade"])
        if not mod:
            continue
        if n in mesures:
            mod["vitellus"] = {"centre": mesures[n]["centre"], "rayon": rayon(n),
                               "source": f"sphère régulière : centre ajusté sur l'atlas, rayon lissé {rayon(n):.2f} mm (mesuré {mesures[n]['rayon_mesure']:.2f})"}
        elif n == 12:
            mid = charge(mod, r"(midgut|intestin_moyen)")
            if mid:
                ic = mid[0].vertices.mean(0); r = rayon(12)
                # côté ventral : −Y (Y dorsal) ; posée à la distance relative médiane des stades mesurés
                c = ic + np.array([0.0, -1.0, 0.0]) * k * r
                mod["vitellus"] = {"centre": [round(float(v), 4) for v in c], "rayon": r,
                                   "source": f"sphère régulière extrapolée (absente de l'atlas à CS12) : rayon {r:.2f} mm, devant l'intestin moyen"}


K_AMNIOS = 1.10
CACHE_AMN = os.path.join("embryons_3D", "amnios_cache.json")


def amnios_regulier(lignes):
    """Ajoute modele['amnios'] = {centre, demi_axes, source} aux stades CS9-CS12 (ellipsoïde englobant le corps)."""
    import re
    import numpy as np
    import trimesh
    try:
        cache = json.load(open(CACHE_AMN, encoding="utf-8"))
    except Exception:
        cache = {}
    for l in lignes:
        mod, n = l["modele"], norm(l["stade"])
        if True:                       # CS9 : amnios des auteurs (PDF 3D) remis, comme dans le PDF ; CS10-CS12 : retiré (consignes 25/09)
            continue
        if not mod or n != 9:
            continue
        cle = f'{mod["dossier"]}|{int(os.path.getmtime(os.path.join(MODELES, mod["dossier"], "manifest.json")))}'
        if cle not in cache:
            corps = [x for x in mod["structures"] if x["ok"] and re.search(r"^(skin|peau|ectoderme_epidermique|epidermal_ectoderm)", x["nom"])]
            if not corps:
                cache[cle] = None
            else:
                P = trimesh.load(os.path.join(MODELES, mod["dossier"], corps[0]["fichier"]), force="mesh").vertices
                lo, hi = np.percentile(P, 2, 0), np.percentile(P, 98, 0)
                cache[cle] = {"centre": [round(float(v), 4) for v in (lo + hi) / 2], "demi_corps": [round(float(v), 4) for v in (hi - lo) / 2]}
        c = cache[cle]
        if c:
            mod["amnios"] = {"centre": c["centre"], "demi_axes": [round(v * K_AMNIOS, 4) for v in c["demi_corps"]],
                             "source": f"ellipsoïde régulier englobant le corps (×{K_AMNIOS:.1f})" + (" — absent de l'atlas à ce stade" if n == 12 else "")}
    json.dump(cache, open(CACHE_AMN, "w", encoding="utf-8"), indent=1)


R_AMNIOS_CS8 = 0.60      # mm : amnios sphérique « beaucoup plus gros » que le dôme de l'atlas (demi-longueur du disque ≈ 0,38)


def disque_cs8(lignes):
    """CS8 : repère du disque par ACP de l'endoderme (normale = petit axe, orientée vers l'amnios de l'atlas = dorsal).
    Vésicule vitelline : ellipsoïde ventral dont la section par le plan du disque est exactement l'ellipse de l'endoderme
    (proportions de la vésicule de l'atlas). Amnios : sphère dorsale de rayon R_AMNIOS_CS8 dont le cercle d'attache
    au plan du disque a le rayon moyen du bord de l'épiblaste."""
    import re
    import numpy as np
    import trimesh
    for l in lignes:
        mod = l["modele"]
        if not mod or norm(l["stade"]) != 8:
            continue
        def pts(motif):
            xs = [x for x in mod["structures"] if x["ok"] and re.search(motif, x["nom"])]
            return np.vstack([trimesh.load(os.path.join(MODELES, mod["dossier"], x["fichier"]), force="mesh").vertices for x in xs]) if xs else None
        E, A, V, B = pts(r"^endoderm$"), pts(r"^(amnion_wall|amniotic_cavity)$"), pts(r"^yolk_sac_wall$"), pts(r"^(epiblast_or_ectoderm|epiblast)$")
        if E is None:
            return
        c = E.mean(0); w, vec = np.linalg.eigh(np.cov((E - c).T))
        n, u, t = vec[:, 0], vec[:, 1], vec[:, 2]
        if A is not None and (A.mean(0) - c) @ n < 0:
            n = -n
        t = np.cross(u, n)                                   # base directe (rotation propre pour le viewer)
        demi = lambda P, ax: float((np.percentile((P - c) @ ax, 98) - np.percentile((P - c) @ ax, 2)) / 2)
        au, at = demi(E, u), demi(E, t)
        # vésicule : proportions de l'atlas (demi-étendues le long de u, n, t)
        if V is not None:
            vu, vn, vt = demi(V, u), demi(V, n), demi(V, t)
        else:
            vu, vn, vt = au * 1.5, 0.46, at * 1.4
        k = float(np.clip((au / vu + at / vt) / 2, 0.3, 0.95))      # section relative au plan du disque
        h = vn * np.sqrt(1 - k * k)
        centre_v = c - n * h
        base = [[round(float(x), 5) for x in v] for v in (u, n, t)]
        mod["vitellus"] = {"centre": [round(float(x), 4) for x in centre_v], "demi_axes": [round(au / k, 4), round(float(vn), 4), round(at / k, 4)],
                           "base": base, "source": "ellipsoïde ventral dont le bord épouse exactement le contour de l'endoderme du disque"}
        # amnios : sphère dorsale sur le bord de l'épiblaste
        P = B if B is not None else E
        a = (demi(P, u) + demi(P, t)) / 2
        R = max(R_AMNIOS_CS8, a * 1.05)
        centre_a = c + n * np.sqrt(R * R - a * a)
        mod["amnios"] = {"centre": [round(float(x), 4) for x in centre_a], "demi_axes": [round(R, 4)] * 3,
                         "source": f"sphère dorsale de {R:.2f} mm de rayon, attachée au bord du disque"}


def complements(lignes):
    """Amnios / vésicule tracés par l'utilisateur (prioritaires)."""
    for l in lignes:
        mod = l["modele"]
        if not mod:
            continue
        d = mod["dossier"] + "_complements"
        f = os.path.join(MODELES, d, "complements.json")
        if not os.path.exists(f):
            continue
        c = json.load(open(f, encoding="utf-8"))
        mod["complements"] = d
        for cle in ("amnios", "vitellus"):
            if c.get(cle) and os.path.exists(os.path.join(MODELES, d, c[cle]["fichier"])):
                mod[cle] = {"fichier": f"modeles/{d}/{c[cle]['fichier']}", "source": c[cle].get("source", "tracé de l'utilisateur"),
                            "etendue_mm": c[cle].get("etendue_mm")}


def construire(publiable_seul=False):
    lignes = []
    for cs, j0, j1, c0, c1, rep in AXE:
        mod = lire_modele(cs, publiable_seul)
        lignes.append({"stade": cs, "j": [j0, j1], "crl": [c0, c1], "reperes": rep, "sources": sources(cs), "modele": mod,
                       "systemes": {s: statut_systeme(mod, s, cs) for s in SYS_IDS}})
    cellules = len(lignes) * len(SYS_IDS)
    try:
        vitellus_regulier(lignes)
    except Exception as e:
        print("vésicule vitelline régulière : échec", e)
    try:
        amnios_regulier(lignes)
    except Exception as e:
        print("amnios régulier : échec", e)
    try:
        disque_cs8(lignes)
    except Exception as e:
        print("disque CS8 : échec", e)
    complements(lignes)
    cellules = sum(1 for l in lignes for v in l["systemes"].values() if v != "na")
    faites = sum(1 for l in lignes for v in l["systemes"].values() if v in ("fait", "externe"))
    etat = json.load(open(ETAT, encoding="utf-8")) if os.path.exists(ETAT) else {}
    todo = taches(lignes)
    for t in todo:
        t["etat"] = etat.get(t["id"], {})
    return {"date": datetime.now().strftime("%d/%m/%Y %H:%M"), "version": datetime.now().strftime("%Y%m%d-%H%M%S"),
            "systemes": [{"id": i, "lbl": l, "col": c} for i, l, c in SYSTEMES],
            "vue_seule": [{"id": i, "lbl": l, "col": c} for i, l, c in VUE_SEULE], "lignes": lignes, "todo": todo,
            "progression": {"faites": faites, "cellules": cellules}, "ref_url": REF_URL}


HTML = r"""<!DOCTYPE html>
<html lang="fr"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="robots" content="noindex">
<title>Embryon 3D · nouveaux modèles</title>
<link rel="preconnect" href="https://fonts.googleapis.com"><link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
<style>
:root{--bg:#f6f4ef;--panel:#fff;--ink:#1d1e20;--mut:#6b6d70;--line:#e3e0d8;--ok:#2f9e6e;--mid:#d9a53a;--no:#e9e6de;--acc:#c8323a}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){--bg:#141517;--panel:#1d1f22;--ink:#ecebe6;--mut:#a3a19b;--line:#2f3236;--no:#2a2c30}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.45 Inter,system-ui,sans-serif}
main{max-width:1180px;margin:0 auto;padding:22px 16px 60px}
h1{font-size:24px;margin:0}h1 small{font-weight:400;color:var(--mut);font-size:14px;margin-left:8px}
h2{font-size:13px;letter-spacing:.06em;text-transform:uppercase;color:var(--mut);margin:30px 0 10px}
.top{display:flex;flex-wrap:wrap;align-items:baseline;gap:10px 18px}.top a{color:var(--ink)}
.card{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:14px 16px}
.kpi{display:flex;gap:26px;flex-wrap:wrap;margin-top:14px}.kpi b{display:block;font-size:22px}.kpi span{color:var(--mut);font-size:12px}
.wrap{overflow-x:auto;border-radius:12px;border:1px solid var(--line);background:var(--panel)}
table{border-collapse:collapse;width:100%;min-width:900px}th,td{padding:7px 8px;border-bottom:1px solid var(--line);text-align:left;vertical-align:middle}
th{font-size:11px;color:var(--mut);font-weight:600;background:var(--panel);position:sticky;top:0}
td.c{text-align:center;width:74px}.dot{display:inline-block;width:14px;height:14px;border-radius:4px;background:var(--no)}
.dot.fait{background:var(--ok)}.dot.brouillon{background:repeating-linear-gradient(45deg,#7fb3e6 0 3px,#cfe3f7 3px 6px)}.dot.externe{background:#2f6fb5}
.badge{font-size:11px;font-weight:600;padding:1px 7px;border-radius:9px;background:#dbe9f8;color:#1d4d80;margin-left:6px}.dot.na{background:transparent;border:1px dashed var(--mut)}.dot.partiel{background:var(--mid)}.src{font-size:12px;color:var(--mut)}.src b{color:var(--ink);font-weight:600}
tr.has td:first-child{box-shadow:inset 3px 0 0 var(--ok)}.st{font-weight:700}.mut{color:var(--mut);font-size:12px}
.todo td{font-size:13px}.pill{font-size:11px;padding:1px 7px;border-radius:9px;background:var(--no)}.pill.haute{background:#f3d3d5;color:#7a1d22}
.pill.en_cours{background:#d8ecff;color:#194a7a}.pill.termine{background:#d5f0e3;color:#1c5c40}
#v3d{display:grid;grid-template-columns:1fr;gap:10px}#gl{width:100%;height:62vh;min-height:360px;display:block;border-radius:12px;background:#0e1013}
#scene{position:relative;display:flex;flex-direction:column;gap:8px}#vue{position:relative}#tiroir{order:-1;display:flex;flex-direction:column;gap:6px}
#ovl{position:absolute;left:12px;top:10px;color:#eceae4;font-size:12px;line-height:1.35;pointer-events:none;text-shadow:0 1px 3px #000;max-width:70%}
#ovl b{display:block;font-size:20px}#ovkpi{position:absolute;right:12px;top:10px;color:#a9aeb4;font-size:11px;text-align:right;pointer-events:none;text-shadow:0 1px 3px #000}
#barre{position:absolute;left:10px;bottom:10px;display:flex;gap:8px;align-items:center}#barre button{font:inherit;font-size:15px;min-width:42px;height:38px;border-radius:10px;border:1px solid #3a3d42;background:rgba(20,21,23,.8);color:#eceae4;cursor:pointer}
#barre #barst{color:#eceae4;font-size:13px;text-shadow:0 1px 3px #000}#tir{display:none}
.stbtn{display:flex;flex-wrap:wrap;gap:6px}.stbtn button{font:inherit;font-size:13px;padding:5px 11px;border-radius:8px;border:1px solid var(--line);background:var(--panel);color:var(--ink);cursor:pointer}
.stbtn button.on{background:var(--ink);color:var(--panel)}.stbtn select{font:inherit;font-size:13px;padding:4px 6px;border-radius:8px;border:1px solid var(--line);background:var(--panel);color:var(--ink)}
#ctl{align-items:center}#play{min-width:40px}.empty{padding:34px 16px;text-align:center;color:var(--mut)}
@media (max-width:800px){
 main{padding:8px 8px 40px}.top{gap:2px 10px}h1{font-size:17px}h1 small{display:none}.top a{font-size:12px}#kpi,#h2v3d{display:none}
 h2{margin:18px 0 8px}#v3d{padding:0;border:0;background:transparent;margin:0 -8px}#v3d>.mut,#fiche{padding:0 10px}
 #gl{height:calc(100dvh - 56px);min-height:420px;border-radius:0}
 #tir{display:inline-block}#barre{left:0;right:0;bottom:0;padding:8px;background:linear-gradient(transparent,rgba(0,0,0,.55))}#barre #barst{flex:1}
 #tiroir{order:0;position:absolute;left:0;right:0;bottom:54px;max-height:58%;overflow:auto;padding:10px;background:rgba(24,25,28,.94);border-top:1px solid #34373b;
  transform:translateY(calc(100% + 60px));transition:transform .22s ease;z-index:3}
 #scene.ouvert #tiroir{transform:none}#tiroir .stbtn button,#tiroir select{background:#26282c;color:#eceae4;border-color:#3a3d42}#tiroir .stbtn button.on{background:#eceae4;color:#111}
 #ovl b{font-size:22px}}
</style></head><body><main>
<div class="top"><h1>Embryon 3D <small>nouveaux modèles à partir de CS10 · remis à zéro le 25/09</small></h1>
 <a href="__REF__">Référence : travail du 24/09 →</a> · <a href="vue/?cs=CS13">Vue détaillée CS11–CS15 (original, modèle, cœur animé) →</a></div>
<div class="kpi" id="kpi"></div>

<h2 id="h2v3d">Modèles 3D</h2>
<div id="v3d" class="card"></div>

<h2>Contenu et reste à faire, par stade et par système</h2>
<div class="wrap"><table id="mat"></table></div>
<p class="mut">Vert : système modélisé par notre code · bleu plein : maillages d'auteurs intégrés (Hikspoors et al. 2022, HDBR atlas, CC BY-NC-SA) · hachuré bleu : brouillon tiré de l'atlas d'Amsterdam, en attente de notre reconstruction · jaune : partiel ou faible · gris : à faire · pointillé : pas encore formé à ce stade. Sources : coupes légendées ehd (nombre d'images), vidéo 360° ehd, VOKA (regard seulement, à partir de J28), référence du 24/09.</p>

<h2>Tâches</h2>
<div class="wrap"><table class="todo" id="todo"></table></div>
<p class="mut" id="gen"></p>
</main>
<script id="data" type="application/json">__DATA__</script>
<script>
const D=JSON.parse(document.getElementById('data').textContent);const $=id=>document.getElementById(id);
const SYS=D.systemes;const avec=D.lignes.filter(l=>l.modele);
const nb=avec.filter(l=>l.modele.brouillon).length;const nx=avec.filter(l=>l.modele.externe).length;
$('kpi').innerHTML=`<div><b>${avec.length} / ${D.lignes.length}</b><span>stades avec un modèle${nb?` (dont ${nb} brouillon${nb>1?'s':''})`:''}${nx?` (dont ${nx} externe${nx>1?'s':''} Hikspoors)`:''}</span></div><div><b>${D.progression.faites} / ${D.progression.cellules}</b><span>cases stade × système modélisées</span></div><div><b>${D.todo.length}</b><span>tâches ouvertes</span></div>`;
// matrice
let h='<tr><th>Stade</th><th>Jours · CRL</th><th>Sources disponibles</th>'+SYS.map(s=>`<th style="text-align:center">${s.lbl}</th>`).join('')+'</tr>';
for(const l of D.lignes){const s=l.sources;const src=[s.coupes?`<b>${s.coupes}</b> coupe${s.coupes>1?'s':''} légendée${s.coupes>1?'s':''}`:'',s.videos360?'<b>vidéo 360°</b>':'',s.voka?'VOKA':'',s.reference_2409?`<a href="${D.ref_url}">réf. 24/09</a>`:''].filter(Boolean).join(' · ')||'—';
 h+=`<tr class="${l.modele?'has':''}"><td><span class="st">${l.stade}</span><div class="mut">${l.reperes}</div></td><td class="src">J${l.j[0]}–${l.j[1]}<br>${l.crl[0]}–${l.crl[1]} mm</td><td class="src">${src}${l.modele?`<div class="mut">${l.modele.brouillon?'<span class="badge">brouillon</span> ':''}${l.modele.externe?'<span class="badge">externe</span> ':''}modèle : ${(l.modele.source||'?').split(' — ')[0]} · ${l.modele.date}</div>`:''}</td>`+
 SYS.map(sy=>`<td class="c" title="${sy.lbl} : ${({fait:'modélisé',externe:'maillages d\'auteurs (Hikspoors)',brouillon:'brouillon (atlas)',partiel:'partiel',a_faire:'à faire',na:'pas encore formé'})[l.systemes[sy.id]]}"><span class="dot ${l.systemes[sy.id]}"></span></td>`).join('')+'</tr>';}
$('mat').innerHTML=h;
// tâches (état partagé publié par le lanceur, rafraîchi toutes les 20 s)
function rendreTaches(etat){let t='<tr><th>Priorité</th><th>Stade</th><th>Action</th><th>Instance</th><th>Statut</th></tr>';
 for(const x of D.todo){const e=(etat&&etat[x.id])||x.etat||{};const st=e.statut||'a_faire';
  t+=`<tr><td><span class="pill ${x.prio}">${x.prio}</span></td><td>${x.stade}</td><td>${x.action}${x.detail?`<div class="mut">${x.detail}</div>`:''}</td><td class="mut">${e.instance||'—'}</td><td><span class="pill ${st}">${st.replace('_',' ')}</span></td></tr>`;}
 $('todo').innerHTML=t;}
rendreTaches(null);
async function majEtat(){try{const r=await fetch('taches_etat.json?t='+Date.now(),{cache:'no-store'});if(r.ok)rendreTaches(await r.json());}catch(e){}}
majEtat();setInterval(majEtat,20000);
$('gen').textContent=`Page générée le ${D.date} par embryo3d/agregateur.py · version ${D.version}`;
setInterval(async()=>{try{const r=await fetch('version.txt?t='+Date.now(),{cache:'no-store'});if(r.ok){const v=(await r.text()).trim();if(v&&v!==D.version)location.reload();}}catch(e){}},20000);
if(!avec.length){$('v3d').innerHTML='<div class="empty">Aucun nouveau modèle pour l\'instant. Le premier stade livré apparaîtra ici, avec un bouton par stade.</div>';}
</script>
<script type="importmap">{"imports":{"three":"https://cdn.jsdelivr.net/npm/three@0.170.0/build/three.module.js","three/addons/":"https://cdn.jsdelivr.net/npm/three@0.170.0/examples/jsm/"}}</script>
<script type="module">
const D=JSON.parse(document.getElementById('data').textContent);const avec=D.lignes.filter(l=>l.modele&&l.modele.structures.some(x=>x.ok&&x.confiance!=='faible'));
if(avec.length){
 const THREE=await import('three');const {OrbitControls}=await import('three/addons/controls/OrbitControls.js');const {PLYLoader}=await import('three/addons/loaders/PLYLoader.js');const {GLTFLoader}=await import('three/addons/loaders/GLTFLoader.js');
 const box=document.getElementById('v3d');box.innerHTML='<div id="scene"><div id="vue"><canvas id="gl"></canvas><div id="ovl"><b id="ovst"></b><span id="ovinfo"></span></div><div id="ovkpi"></div>'+
  '<div id="barre"><button id="play" class="on" title="lecture / pause">❚❚</button><span id="barst"></span><button id="tir" title="stades et systèmes">☰</button></div></div>'+
  '<div id="tiroir"><div class="stbtn" id="ctl"><select id="dur" title="durée par stade"><option value="3">3 s / stade</option><option value="5" selected>5 s / stade</option><option value="8">8 s / stade</option></select><span id="stb" class="stbtn"></span></div><div class="stbtn" id="syb"></div></div>'+
  '</div><div class="mut" id="leg"></div><div id="fiche"></div>';
 const tir=document.getElementById('tir');document.getElementById('gl').addEventListener('pointerdown',()=>document.getElementById('scene').classList.remove('ouvert'));tir.onclick=()=>document.getElementById('scene').classList.toggle('ouvert');
 const canvas=document.getElementById('gl');const R=new THREE.WebGLRenderer({canvas,antialias:true});R.setPixelRatio(Math.min(devicePixelRatio,2));
 const scene=new THREE.Scene();scene.background=new THREE.Color(0x0e1013);scene.add(new THREE.HemisphereLight(0xffffff,0x334455,2.4));const k=new THREE.DirectionalLight(0xffffff,2.4);k.position.set(-3,2,4);scene.add(k);
 const cam=new THREE.PerspectiveCamera(32,1,0.01,500);cam.up.set(0,0,1);const ctl=new OrbitControls(cam,canvas);ctl.enableDamping=true;
 function size(){const w=canvas.clientWidth,h=canvas.clientHeight;R.setSize(w,h,false);cam.aspect=w/h;cam.updateProjectionMatrix();if(window.__pret)cadrer();}new ResizeObserver(size).observe(canvas);size();
 const COL=Object.fromEntries(D.systemes.concat(D.vue_seule||[]).map(s=>[s.id,s.col]));let grp=null;const cache={};const on={};
 const TRANSL=/^(amnios|amnion|vesicule_vitelline|yolk_sac|chorion|extra_embryonic|exocoelom)/;
 // peau : toujours très transparente, on voit dessous (consigne utilisateur 25/09 : « transparente tout court »)
 const PEAU=/^(skin|peau|enveloppe|ectoderme_epidermique|epidermal_ectoderm)/;const peauCachee=true;
 const PEAU_OPACITE=0.18;
 const SOMITES=/^(somites$|somite |plaque segmentaire|plaque_segmentaire)/i;
 function vis(){if(grp)grp.children.forEach(o=>{o.visible=on[o.userData.sys]!==false;
  if(o.userData.peau)o.traverse(m=>{if(m.isMesh)m.renderOrder=2;});});}
 // cadrage sur l'embryon lui-même : les annexes (cavité amniotique, vésicule vitelline…) sont exclues de la boîte (CS7-CS9 : disque embryonnaire)
 function cadrer(){grp.updateMatrixWorld(true);const b=new THREE.Box3();grp.children.forEach(o=>{if(o.userData.sys!=='annexes'||o.userData.cadre)b.expandByObject(o);});/* l'amnios régulier entre dans le cadre */if(b.isEmpty())b.setFromObject(grp);const c=b.getCenter(new THREE.Vector3());
  const dir=new THREE.Vector3(-1,0,0.15).normalize();const fwd=dir.clone().negate(),right=new THREE.Vector3().crossVectors(fwd,cam.up).normalize(),upv=new THREE.Vector3().crossVectors(right,fwd).normalize();
  const tv=Math.tan(THREE.MathUtils.degToRad(cam.fov/2)),th=tv*Math.max(0.3,cam.aspect||1);let need=0;
  for(let i=0;i<8;i++){const q=new THREE.Vector3(i&1?b.max.x:b.min.x,i&2?b.max.y:b.min.y,i&4?b.max.z:b.min.z).sub(c);const z=q.dot(dir);need=Math.max(need,Math.abs(q.dot(right))/th+z,Math.abs(q.dot(upv))/tv+z);}
  const d=need*(window.innerWidth<=800?1.04:1.12);cam.position.copy(c).addScaledVector(dir,d);ctl.target.copy(c);cam.near=d/200;cam.far=d*20;cam.updateProjectionMatrix();ctl.update();}
 // amnios et vésicule vitelline : voiles colorés transparents jusqu'à CS12 (l'atlas s'arrête à CS11) ; paroi préférée à la cavité si les deux existent
 const AMNIOS=/^(amnion_wall|amnios|amniotic_cavity)$/,VITELLUS=/^(yolk_sac_wall|yolk_sac|vesicule_vitelline|yolk_sac_cavity)$/;
 function voile(c){return new THREE.ShaderMaterial({transparent:true,depthWrite:false,side:THREE.DoubleSide,uniforms:{col:{value:new THREE.Color(c)}},
  vertexShader:'varying vec3 vN;varying vec3 vV;void main(){vN=normalize(normalMatrix*normal);vec4 mv=modelViewMatrix*vec4(position,1.0);vV=normalize(-mv.xyz);gl_Position=projectionMatrix*mv;}',
  fragmentShader:'uniform vec3 col;varying vec3 vN;varying vec3 vV;void main(){float f=pow(1.0-abs(dot(normalize(vN),normalize(vV))),2.0);gl_FragColor=vec4(col,0.06+0.5*f);}'});}
 async function charger(l,afficher){const base=`modeles/${l.modele.dossier}/`;const g=new THREE.Group();const nst=parseInt(l.stade.replace(/\D/g,''));
  const noms=new Set(l.modele.structures.map(x=>x.nom));
  const xs=l.modele.structures.filter(x=>x.ok&&x.confiance!=='faible').filter(x=>{
   if(VITELLUS.test(x.nom)&&l.modele.vitellus)return false;   // remplacée par la sphère régulière
   if(AMNIOS.test(x.nom)&&(l.modele.amnios||(nst>=10&&nst<=12)))return false;   /* remplacé (CS8, CS9) ou retiré (CS10-CS12) */
   if(/^(trophoblastic_lacunae|exocoelomic_cyst|extra_embryonic_coelom|extra_embryonic_(splanchnopleuric|somatopleuric)_mesoderm)$/.test(x.nom)||(nst<=8&&x.nom==='connecting_stalk'))return false;   // grandes cavités extra-embryonnaires : masquaient le disque (CS7)
   if(AMNIOS.test(x.nom)||VITELLUS.test(x.nom)){if(nst>12)return false;if(/_cavity$/.test(x.nom)&&noms.has(x.nom.replace('_cavity','_wall')))return false;if(x.nom==='amniotic_cavity'&&noms.has('amnion_wall'))return false;}
   return true;});let fait=0;
  const un=async x=>{const url=base+x.fichier;
   const col=(x.couleur&&!['annexes','tissus'].includes(x.systeme))?new THREE.Color(...x.couleur):new THREE.Color(COL[x.systeme]||'#999');   // annexes et tissus : gris neutre (la vésicule vitelline jaune se confondait avec le digestif)
   const mat=AMNIOS.test(x.nom)?voile('#9fd0ff'):VITELLUS.test(x.nom)?voile('#f3c872'):PEAU.test(x.nom)?new THREE.ShaderMaterial({transparent:true,depthWrite:false,side:THREE.DoubleSide,uniforms:{col:{value:col.clone().lerp(new THREE.Color(1,1,1),0.25)}},
     vertexShader:'varying vec3 vN;varying vec3 vV;void main(){vN=normalize(normalMatrix*normal);vec4 mv=modelViewMatrix*vec4(position,1.0);vV=normalize(-mv.xyz);gl_Position=projectionMatrix*mv;}',
     fragmentShader:'uniform vec3 col;varying vec3 vN;varying vec3 vV;void main(){float f=pow(1.0-abs(dot(normalize(vN),normalize(vV))),2.2);gl_FragColor=vec4(col,0.07+0.6*f);}'})
    :new THREE.MeshStandardMaterial({color:col,roughness:0.55,side:THREE.DoubleSide});   // peau : contour lumineux, centre transparent
   if(TRANSL.test(x.nom)||x.systeme==='tissus'||x.systeme==='annexes'){mat.transparent=true;mat.opacity=x.systeme==='tissus'?0.22:0.14;mat.depthWrite=false;}   // annexes, cavités et tissus : translucides
   else if(SOMITES.test(x.nom)){mat.transparent=true;mat.opacity=0.35;mat.depthWrite=false;}   // somites (et leurs parties) : translucides, on voit le tube neural et les aortes
   try{if(/\.glb$/i.test(x.fichier)){const s=(await new GLTFLoader().loadAsync(url)).scene;s.traverse(o=>{if(o.isMesh){if(!o.geometry.attributes.normal)o.geometry.computeVertexNormals();o.material=mat;}});s.userData.sys=x.systeme;s.userData.nom=x.nom;s.userData.peau=PEAU.test(x.nom);g.add(s);}
    else{const geo=await new PLYLoader().loadAsync(url);geo.computeVertexNormals();const me=new THREE.Mesh(geo,mat);me.userData.sys=x.systeme;me.userData.peau=PEAU.test(x.nom);g.add(me);}}catch(e){console.warn(url,e);}
   fait++;if(afficher)document.getElementById('leg').textContent=`chargement ${l.stade} : ${fait} / ${xs.length} structures…`;};
  for(let i=0;i<xs.length;i+=8)await Promise.all(xs.slice(i,i+8).map(un));   // 8 fichiers à la fois
  const ellipsoide=(e,c)=>{const me=new THREE.Mesh(new THREE.SphereGeometry(1,72,48),voile(c));const d=e.demi_axes||[e.rayon,e.rayon,e.rayon];me.scale.set(...d);me.position.set(...e.centre);
   if(e.base){const m=new THREE.Matrix4().makeBasis(new THREE.Vector3(...e.base[0]),new THREE.Vector3(...e.base[1]),new THREE.Vector3(...e.base[2]));me.quaternion.setFromRotationMatrix(m);}me.userData.sys='annexes';return me;};
  const enveloppe=async(e,c)=>{if(!e.fichier)return ellipsoide(e,c);const geo=await new PLYLoader().loadAsync(e.fichier);geo.computeVertexNormals();const me=new THREE.Mesh(geo,voile(c));me.userData.sys='annexes';return me;};
  if(l.modele.amnios&&nst<=12){const me=await enveloppe(l.modele.amnios,'#9fd0ff');me.userData.cadre=true;g.add(me);}
  if(l.modele.vitellus&&nst<=12)g.add(await enveloppe(l.modele.vitellus,'#f3c872'));
  if(l.modele.orientation&&l.modele.orientation.rot_z_180)g.rotation.z=Math.PI;   // demi-tour autour de l'axe crânio-caudal
  return g;}
 let courant=null;async function montrer(l){courant=l;document.querySelectorAll('#stb button').forEach(b=>b.classList.toggle('on',b.dataset.s===l.stade));
  if(!cache[l.stade])cache[l.stade]=await charger(l,true);if(courant!==l)return;if(grp)scene.remove(grp);grp=cache[l.stade];scene.add(grp);vis();cadrer();window.__pret=true;window.__grp=grp;
  const pres=[...new Set(l.modele.structures.filter(x=>x.ok&&x.confiance!=='faible').map(x=>x.systeme))];const LBL=Object.fromEntries(D.systemes.concat(D.vue_seule||[]).map(s=>[s.id,s.lbl]));
  document.getElementById('syb').innerHTML=pres.map(sid=>`<button data-y="${sid}" class="${on[sid]!==false?'on':''}"><span style="display:inline-block;width:10px;height:10px;border-radius:3px;background:${COL[sid]||'#999'};margin-right:6px"></span>${LBL[sid]||sid}</button>`).join('');
  document.querySelectorAll('#syb button').forEach(b=>b.onclick=()=>{const y=b.dataset.y;on[y]=on[y]===false;b.classList.toggle('on',on[y]!==false);vis();});
  const m=l.modele,t=m.temps||{},dm=m.dimensions||{};const par={};m.structures.forEach(x=>(par[x.systeme]=par[x.systeme]||[]).push(x));
  document.getElementById('fiche').innerHTML=`<div class="mut" style="margin:8px 0">${t.jours_post_fecondation?`J${t.jours_post_fecondation.join('–')} · `:''}${dm.longueur_atlas_mm?`${dm.longueur_atlas_mm} mm · `:''}${dm.paires_somites?`${dm.paires_somites} paires de somites · `:''}${m.specimen?`spécimen ${m.specimen} · `:''}${m.vitellus?`vésicule vitelline : ${m.vitellus.source} · `:''}${m.amnios?`amnios : ${m.amnios.source} · `:''}${m.orientation&&m.orientation.rot_z_180?'<b>orientation corrigée à l\'affichage (demi-tour, le cœur était dorsal)</b> · ':''}${t.reperes_stade||''}</div>`+
   Object.entries(par).map(([sid,xs])=>`<details><summary><b>${LBL[sid]||sid}</b> <span class="mut">(${xs.length})</span></summary><ul style="margin:4px 0 8px;padding-left:18px">${xs.map(x=>`<li><b>${x.nom_fr}</b>${x.confiance==='faible'?' <span class="pill">faible, non affiché</span>':''}${x.verif==='a_corriger'?` <span class="pill haute" title="${(x.verif_pb||[]).join(' ; ')}">maillage à corriger</span>`:''}${x.legende?` — <span class="mut">${x.legende}</span>`:''}${x.ca?`<br><span class="mut" style="font-size:.92em">coupes Amsterdam : ${x.ca.l} (coupes ${x.ca.a}–${x.ca.b}, max en ${x.ca.c})${x.ca.g?' <i>via la géométrie</i>':''}${x.ca.d?' <span class="pill">nom et géométrie en désaccord</span>':''}</span>`:''}</li>`).join('')}</ul></details>`).join('')+
   (m.verif?`<p class="mut" style="margin:8px 0 0">Contrôle des maillages : <b>${m.verif.ok} / ${m.verif.total}</b> sans erreur (étanchéité, arêtes, faces dégénérées ou en double).</p>`:'')+(m.brouillon?`<p style="margin:10px 0 0"><span class="badge">BROUILLON</span> <span class="mut">Modèle provisoire tiré de l'atlas, en attendant la reconstruction par notre code.</span></p>`:'')+(!m.publiable?`<p class="mut" style="margin-top:10px"><b>Référence interne, non publiable</b>${m.usage?' ('+m.usage+')':''} : visible seulement en local.</p>`:'')+(m.licence?`<p class="mut" style="margin-top:10px">Source : ${m.source}. Licence ${m.licence} — <a href="https://creativecommons.org/licenses/${/by-nc-sa/i.test(m.licence)?'by-nc-sa':/by-sa/i.test(m.licence)?'by-sa':/by-nc-nd/i.test(m.licence)?'by-nc-nd':'by'}/4.0/deed.fr" target="_blank">conditions</a>. ${m.attribution?`Attribution : ${m.attribution}. `:''}Géométrie originale des auteurs, seulement mise à l'échelle et positionnée ; couleurs par système.${m.echelle&&m.echelle.um_par_unite?` Échelle ${m.echelle.um_par_unite} µm/unité.`:''}${m.orientation_auto&&m.orientation_auto.score!==undefined?` Orientation automatique (score ${m.orientation_auto.score}) : vérifier controle.png.`:''}</p>`:'');
  const nS=l.modele.structures.filter(x=>x.ok&&x.confiance!=='faible').length,tt=l.modele.temps||{},dd=l.modele.dimensions||{};
  document.getElementById('leg').textContent=`${l.stade} · ${nS} structures · source : ${l.modele.source||'?'}${l.modele.session?' · '+l.modele.session:''}`;
  document.getElementById('ovst').textContent=l.stade+(l.modele.brouillon?' · brouillon':'');
  document.getElementById('ovinfo').innerHTML=`J${(tt.jours_post_fecondation||l.j).join('–')} · ${dd.longueur_atlas_mm?dd.longueur_atlas_mm+' mm':l.crl.join('–')+' mm'} · ${nS} structures<br>${l.reperes}`;
  document.getElementById('barst').textContent=l.stade;
  const nb=avec.length,nbb=avec.filter(x=>x.modele.brouillon).length;document.getElementById('ovkpi').innerHTML=`${nb} stades · ${nbb} brouillons<br>${D.progression.faites}/${D.progression.cellules} cases modélisées`;}
 document.getElementById('stb').innerHTML=avec.map(l=>`<button data-s="${l.stade}">${l.stade}${l.modele.brouillon?' · brouillon':''}</button>`).join('');
 // diaporama en boucle : stade suivant toutes les N secondes (préchargé), clic sur un stade = pause
 let lecture=true,tSuiv=0;const btn=document.getElementById('play');
 const majBtn=()=>{btn.textContent=lecture?'❚❚':'▶';btn.classList.toggle('on',lecture);};
 btn.onclick=()=>{lecture=!lecture;tSuiv=performance.now()+(+document.getElementById('dur').value)*1000;majBtn();};
 document.querySelectorAll('#stb button').forEach(b=>b.onclick=()=>{lecture=false;majBtn();montrer(avec.find(l=>l.stade===b.dataset.s));});
 const suivant=()=>avec[(avec.indexOf(courant)+1)%avec.length];
 async function precharger(){const n=suivant();if(n&&!cache[n.stade])cache[n.stade]=await charger(n);}
 await montrer(avec[0]);precharger();tSuiv=performance.now()+(+document.getElementById('dur').value)*1000;
 (function loop(now){requestAnimationFrame(loop);
  ctl.update();R.render(scene,cam);})(performance.now());
 setInterval(()=>{const now=performance.now();if(lecture&&avec.length>1&&now>=tSuiv){const n=suivant();if(cache[n.stade]){tSuiv=now+(+document.getElementById('dur').value)*1000;montrer(n).then(precharger);}else precharger();}},250);
}
</script></body></html>
"""


def page(data):
    js = json.dumps(data, ensure_ascii=False).replace("</", "<\\/")
    return HTML.replace("__DATA__", js).replace("__REF__", data["ref_url"])


def main():
    data = construire()
    with open(SORTIE_JSON, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=1)
    loc = dict(data, ref_url="reference_2409/site_3dht/index.html")
    with open(SORTIE_HTML, "w", encoding="utf-8") as f:
        f.write(page(loc))
    with open(os.path.join("embryons_3D", "version.txt"), "w") as f:      # la page locale compare sa version à ce fichier
        f.write(data["version"])
    print(f"OK {SORTIE_HTML} · {sum(1 for l in data['lignes'] if l['modele'])} stades avec modèle · {len(data['todo'])} tâches")
    if "--site" in sys.argv:
        dest = sys.argv[sys.argv.index("--site") + 1]
        os.makedirs(dest, exist_ok=True)
        data = construire(publiable_seul=True)          # jamais de modèle publiable=false en ligne
        data["version"] = datetime.now().strftime("%Y%m%d-%H%M%S")
        garde = {l["modele"]["dossier"] for l in data["lignes"] if l["modele"]} | {l["modele"]["complements"] for l in data["lignes"] if l["modele"] and l["modele"].get("complements")}
        if os.path.isdir(os.path.join(dest, "modeles")):
            for n in os.listdir(os.path.join(dest, "modeles")):
                if n not in garde:
                    shutil.rmtree(os.path.join(dest, "modeles", n))
        for l in data["lignes"]:
            if l["modele"]:
                src = os.path.join(MODELES, l["modele"]["dossier"]); dst = os.path.join(dest, "modeles", l["modele"]["dossier"])
                os.makedirs(dst, exist_ok=True)
                for x in l["modele"]["structures"]:
                    if x["ok"]:
                        shutil.copy2(os.path.join(src, x["fichier"]), os.path.join(dst, x["fichier"]))
                if l["modele"].get("complements"):
                    cs_, cd_ = os.path.join(MODELES, l["modele"]["complements"]), os.path.join(dest, "modeles", l["modele"]["complements"])
                    os.makedirs(cd_, exist_ok=True)
                    for f_ in os.listdir(cs_):
                        if f_.endswith((".ply", ".json")):
                            shutil.copy2(os.path.join(cs_, f_), os.path.join(cd_, f_))
        with open(os.path.join(dest, "index.html"), "w", encoding="utf-8") as f:
            f.write(page(data))
        with open(os.path.join(dest, "version.txt"), "w") as f:
            f.write(data["version"])
        print(f"paquet site : {dest} → publier modeles/ puis index.html et version.txt")


if __name__ == "__main__":
    main()
