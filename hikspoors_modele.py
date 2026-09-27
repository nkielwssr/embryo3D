# -*- coding: utf-8 -*-
"""Maillages des PDF 3D Hikspoors et al. 2022 (Maastricht, hdbratlas.org/hikspoors-pdf, licence HDBR CC BY-NC-SA 4.0), déjà extraits
des PDF par la session « Carnegie Stage 13 PDF 3D » (bloc U3D 0x100), -> modèle au format de l'agrégateur :
embryons_3D/modeles/<CS>_hikspoors/{<structure>.ply, manifest.json, rapport.json, controle.png, verif_maillages.json}

    python embryo3d/hikspoors_modele.py CS13 --banque "C:/Users/MicroTurtle/Documents/ChatGPT/Embryo/reference/hikspoors_maastricht"
        -> u3d/Carnegie_Stage_13[_NEWvalves].npz (+ _scene.json) de la banque de la session « PDF 3D » ; --banque est la valeur par défaut
    python embryo3d/hikspoors_modele.py CS13 --npz .../u3d/Carnegie_Stage_13.npz [--scene .../Carnegie_Stage_13_scene.json]
    python embryo3d/hikspoors_modele.py CS13 --glb .../glb/Carnegie_Stage_13.glb   # scène GLB (nœuds « structure~rrggbb » déjà scindés par couleur)
    options : --calage C.json|T.npy [--calage-repere modeles|pipeline]   similitude connue unités PDF -> mm au lieu de l'orientation auto :
                                            .json (clé matrice_pdf_vers_mm, 4x4), ou 13 nombres s,R,t, 12 nombres A,t, matrice 4x4 / 3x4.
                                            DÉFAUT : calage/<nom du npz>.json de la banque s'il existe (calage ICP de la session « PDF 3D » sur
                                            embryons_3D/modeles/<CS>/, repère « modeles ») ; les coordonnées restent alors celles du modèle du
                                            stade (pas de recentrage) : les structures Hikspoors se superposent à notre modèle.
              --sans-calage                 ignorer le calage de la banque (orientation automatique)
              --um-par-unite 1.078          échelle sans calage (µm par unité U3D ; défaut : hikspoors_echelles.json puis 1.0)
              --retourner z|y               inverse un axe après l'orientation automatique (x suit pour rester direct)
              --separer-couleurs            une partie portant plusieurs couleurs de faces (texture) est scindée par couleur
              --sans-completer              version _NEWvalves : ne pas reprendre de l'originale les structures absentes (ex. gut à CS18)
              --sans-unions                 ne pas écrire les unions coeur / cavites_cardiaques / arteres / veines / intestin
              --publiable                   autoriser la publication (défaut : référence interne, publiable=false : atlas d'auteurs retravaillé)
              --sortie DOSSIER              défaut embryons_3D/modeles/<CS>_hikspoors

Format npz de la banque (decode_u3d_local.py) : clés « FACESET_<nom>|v » (sommets float32, unités PDF), « |f » (faces int32), « |c » (couleur
matériau, 3 flottants 0-1), « |p » (structure parente, chaîne), « |fc » (couleur par face uint8 0-255 pour les surfaces texturées : deux
structures sur une surface, scindées avec --separer-couleurs). Les stades disponibles : 9…18, 20, 23 (+ 18/20/23 _NEWvalves, préférées ; les
structures qui n'existent que dans l'originale y sont ajoutées si les deux versions sont dans le même repère).

Repère de sortie (convention agrégateur) : mm, Z crânial, Y dorsal, X = Y×Z = gauche anatomique (repère direct).
Orientation automatique (repli, sans calage) : candidats Z = ±axes de l'ACP et ±axe anatomique (foie -> cœur -> arcs aortiques) ; dorsal =
aorte dorsale, somites, veines cardinales contre le cœur, mesurés dans la tranche crânio-caudale du cœur (un embryon en C fait mentir les
barycentres : le tube neural fait le tour du cœur) ; choix par vote (cœur au-dessus du foie, arcs au-dessus du cœur, voie efférente au-dessus du
sinus veineux, structures gauches à +X…). Tout est consigné dans rapport.json et controle.png : VÉRIFIER la planche, corriger avec --retourner.
"""
import argparse
import json
import os
import re
import sys
import unicodedata

import numpy as np
import trimesh

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ICI = os.path.dirname(os.path.abspath(__file__))
CWD0 = os.getcwd()
os.chdir(RACINE)


def abs0(p):
    """chemin donné sur la ligne de commande : relatif au dossier de lancement (les scripts travaillent depuis RACINE)"""
    return None if p is None else (p if os.path.isabs(p) else os.path.normpath(os.path.join(CWD0, p)))
sys.path.insert(0, ICI)
try:
    from agregateur import AXE, SYSTEMES, VUE_SEULE, reclasser
except Exception:                       # agrégateur absent : valeurs minimales
    AXE = [("CS%d" % i, 0, 0, 0, 0, "") for i in range(7, 24)]
    SYSTEMES, VUE_SEULE = [], []
    def reclasser(nom, systeme):
        return systeme

# CRL typique (mm) pour le contrôle d'échelle (O'Rahilly & Müller) — identique à meshexport.CRL_MM, étendu aux stades jeunes
CRL_MM = {'CS9': 1.8, 'CS10': 2.5, 'CS11': 3.5, 'CS12': 4.0, 'CS13': 4.5, 'CS14': 6.0, 'CS15': 7.5, 'CS16': 9.0, 'CS17': 11.0, 'CS18': 14.0,
          'CS19': 17.0, 'CS20': 20.0, 'CS21': 22.0, 'CS22': 25.0, 'CS23': 28.0}
SPECIMENS = {'CS9': '3709', 'CS10': '6330', 'CS11': '6344', 'CS12': '8943', 'CS13': '836', 'CS14': '6502', 'CS15': '721', 'CS16': '6517',
             'CS17': '6520', 'CS18': '4430', 'CS20': '462', 'CS23': '9226'}          # collection Carnegie (Hikspoors 2022, table suppl. 2)
SOURCE = "Hikspoors et al. 2022, Commun Biol (doi 10.1038/s42003-022-03153-x) — PDF 3D hdbratlas.org/hikspoors-pdf"
LICENCE = "CC BY-NC-SA 4.0 (HDBR atlas, hdbratlas.org/copyright.html)"
BANQUE = "C:/Users/MicroTurtle/Documents/ChatGPT/Embryo/reference/hikspoors_maastricht"     # banque de la session « Carnegie Stage 13 PDF 3D »
NOMENCLATURE = os.path.join(ICI, "hikspoors_nomenclature.json")
ECHELLES = os.path.join(ICI, "hikspoors_echelles.json")
FEMININ = {"veine", "veines", "artère", "artères", "valve", "valvules", "carotide", "vésicule", "vésicules", "crête", "crêtes", "cavité",
           "courbure", "cellules", "glande", "trachée", "rate", "voies"}


# ----------------------------------------------------------------------------------------------------------------------- utilitaires
def normaliser(nom):
    s = unicodedata.normalize("NFKD", str(nom)).encode("ascii", "ignore").decode("ascii").lower()
    s = re.sub(r"[^a-z0-9]+", "_", s).strip("_")
    return re.sub(r"_+", "_", s)


def slug(nom):
    return normaliser(nom) or "sans_nom"


def stade_norme(cs):
    n = "".join(ch for ch in str(cs) if ch.isdigit())
    return "CS%d" % int(n) if n else str(cs)


def couleur_255(c):
    if c is None:
        return None
    c = np.asarray(c, float).ravel()[:3]
    if c.max() <= 1.0:
        c = c * 255.0
    return [int(round(float(v))) for v in c]


def couleur_01(c):
    c = couleur_255(c)
    return None if c is None else [round(v / 255.0, 4) for v in c]


# ----------------------------------------------------------------------------------------------------------------------- lecture
def trouver_dans_banque(banque, cs, prefere_newvalves=True):
    """u3d/Carnegie_Stage_<n>[_NEWvalves].npz (+ _scene.json) ou glb/Carnegie_Stage_<n>.glb dans la banque.
    Renvoie {npz, scene, glb, base, original} (original : {npz, scene, base} de la version d'origine quand la _NEWvalves est prise) ou None."""
    n = str(norm_int(cs))
    cands = (["Carnegie_Stage_%s_NEWvalves" % n] if prefere_newvalves else []) + ["Carnegie_Stage_%s" % n]

    def npz_de(base):
        npz = os.path.join(banque, "u3d", base + ".npz")
        if os.path.exists(npz):
            sc = os.path.join(banque, "u3d", base + "_scene.json")
            return {"npz": npz, "scene": sc if os.path.exists(sc) else None, "glb": None, "base": base}
        return None
    for base in cands:
        r = npz_de(base)
        if r:
            r["original"] = npz_de(base[:-len("_NEWvalves")]) if base.endswith("_NEWvalves") else None
            return r
    for base in cands:
        glb = os.path.join(banque, "glb", base + ".glb")
        if os.path.exists(glb):
            return {"npz": None, "scene": None, "glb": glb, "base": base, "original": None}
    return None


def trouver_calage(banque, base):
    """calage/<base>.json de la banque (session « PDF 3D ») ; pour une version _NEWvalves sans calage propre, celui de l'originale (signalé).
    Renvoie (chemin | None, note)."""
    if not banque or not base:
        return None, ""
    for b, note in ((base, ""), (re.sub(r"_NEWvalves$", "", base), "calage de la version d'origine")):
        f = os.path.join(banque, "calage", b + ".json")
        if os.path.exists(f):
            return f, note
    return None, ""


CLES_MATRICE = ("matrice_pdf_vers_mm", "matrice", "matrix", "M", "T", "transformation", "transform")


def lire_calage(chemin):
    """.json (clé matrice_pdf_vers_mm : 4x4, unités PDF -> mm) ; .npy/.txt : 13 nombres (s, R, t), 12 (A, t), 4x4 ou 3x4.
    Renvoie (A 3x3, t, infos) avec x_mm = A x_pdf + t ; infos : clé, clés du json, stade visé s'il est indiqué (fichier « _vers_CS15 » ou champ)."""
    infos = {"fichier": chemin}
    if chemin.lower().endswith(".json"):
        j = json.load(open(chemin, encoding="utf-8"))
        if not isinstance(j, dict):
            raise SystemExit("calage %s : objet JSON attendu" % chemin)
        infos["cles_json"] = sorted(j.keys())
        cle = next((k for k in CLES_MATRICE if k in j), None)
        if cle is None:
            raise SystemExit("calage %s : aucune clé %s (clés : %s)" % (chemin, " / ".join(CLES_MATRICE), infos["cles_json"]))
        infos["cle"] = cle
        for k in ("modele", "modele_cible", "stade_modele", "cible", "reference", "vers", "stade"):
            if isinstance(j.get(k), str) and re.search(r"CS\s*\d+", j[k], re.I):
                infos["cible"] = j[k]
                break
        T = np.asarray(j[cle], float)
    elif chemin.lower().endswith(".txt"):
        T = np.loadtxt(chemin)
    else:
        T = np.load(chemin, allow_pickle=True)
    T = np.asarray(T, float)
    if T.shape == (4, 4) and not np.allclose(T[3], [0, 0, 0, 1]) and np.allclose(T[:, 3], [0, 0, 0, 1]):
        T = T.T; infos["transposee"] = True                  # matrice rangée par colonnes (translation en bas)
    if T.shape in ((4, 4), (3, 4)):
        A, t = T[:3, :3], T[:3, 3]
    elif T.size == 13:
        T = T.ravel(); A, t = T[0] * T[1:10].reshape(3, 3), T[10:13]
    elif T.size == 12:
        T = T.ravel(); A, t = T[:9].reshape(3, 3), T[9:12]
    else:
        raise SystemExit("calage : forme %s non reconnue (json matrice_pdf_vers_mm ; 13 nombres s,R,t ; 12 nombres A,t ; 4x4 ou 3x4)" % (T.shape,))
    m = re.search(r"_vers_(CS\d+)", os.path.basename(chemin), re.I)
    if m and "cible" not in infos:
        infos["cible"] = m.group(1)
    return np.asarray(A, float), np.asarray(t, float), infos


def meme_repere(pa, pb, tol=0.02):
    """Deux versions d'un même PDF (_NEWvalves / originale) sont-elles dans le même repère ? Boîtes des parties de même nom comparées à
    l'étendue totale. Renvoie (True | False | None, texte)."""
    da = {normaliser(p["nom"]): p["mesh"] for p in pa}
    db = {normaliser(p["nom"]): p["mesh"] for p in pb}
    communs = [k for k in da if k in db]
    if not communs:
        return None, "aucune partie de même nom"
    B = np.vstack([da[k].bounds for k in communs])
    ext = max(float((B.max(0) - B.min(0)).max()), 1e-9)
    ecarts = [float(np.abs(da[k].bounds - db[k].bounds).max()) / ext for k in communs]
    med = float(np.median(ecarts))
    return med < tol, "%d parties de même nom, écart médian des boîtes %.2g %% de l'étendue" % (len(communs), 100 * med)


def completer_depuis_original(parties, originales, nomen):
    """Parties de la version d'origine dont le nom canonique manque dans la version _NEWvalves (ex. « gut » perdu à CS18 NEWvalves) ;
    CCS / CCS_1 ou Asc_Ao_wall / asc_Ao_wall ont le même nom canonique et ne sont pas doublés."""
    def canon(p):
        r = nomen.chercher(p["nom"])
        return None if r.get("ignorer") else r["nom"]
    presents = {canon(p) for p in parties}
    return [dict(p, complement=os.path.basename(p["fichier"])) for p in originales if canon(p) is not None and canon(p) not in presents]


def norm_int(cs):
    n = "".join(ch for ch in str(cs) if ch.isdigit())
    return int(n) if n else 0


def _nom_noeud(scene, noeud, geom):
    n = str(noeud)
    if re.match(r"^(geometry|mesh|node|world)_?\d*$", n, re.I) or n in ("world", ""):
        n = str(geom)
    return n


def _couleur_visuel(m):
    """couleur dominante d'un maillage trimesh (matériau ou couleurs de faces/sommets), en 0-255, et couleurs par face si présentes"""
    try:
        v = m.visual
        if v.kind == "face" and len(v.face_colors):
            fc = np.asarray(v.face_colors)[:, :3]
            return couleur_255(np.median(fc, axis=0)), fc
        if v.kind == "vertex" and len(v.vertex_colors):
            vc = np.asarray(v.vertex_colors)[:, :3]
            fc = vc[np.asarray(m.faces)].mean(axis=1)
            return couleur_255(np.median(vc, axis=0)), fc
        if hasattr(v, "material") and v.material is not None:
            mat = v.material
            if getattr(mat, "baseColorFactor", None) is not None:
                return couleur_255(mat.baseColorFactor), None
            if getattr(mat, "main_color", None) is not None:
                return couleur_255(mat.main_color), None
            if getattr(mat, "diffuse", None) is not None:
                return couleur_255(mat.diffuse), None
    except Exception:
        pass
    return None, None


def charger_fichiers(chemin, cs):
    """Dossier (ou fichier) de maillages : GLB/GLTF (scène ou objet), PLY, OBJ, STL. Renvoie [{nom, mesh, couleur, couleurs_faces, fichier}]."""
    fichiers = []
    if os.path.isdir(chemin):
        sous = [os.path.join(chemin, d) for d in (cs, cs.lower(), stade_norme(cs), stade_norme(cs).lower())]
        for d in sous:
            if os.path.isdir(d):
                chemin = d
                break
        for f in sorted(os.listdir(chemin)):
            if f.lower().endswith((".glb", ".gltf", ".ply", ".obj", ".stl")):
                fichiers.append(os.path.join(chemin, f))
        if not fichiers:                       # fichiers nommés <CS>_* ou *<CS>* dans le dossier
            for f in sorted(os.listdir(chemin)):
                if stade_norme(cs).lower() in f.lower() and f.lower().endswith((".glb", ".gltf", ".ply", ".obj", ".stl")):
                    fichiers.append(os.path.join(chemin, f))
    else:
        fichiers = [chemin]
    parties = []
    for f in fichiers:
        try:
            obj = trimesh.load(f, force=None, process=False)
        except Exception as e:
            print("  lecture impossible :", f, e)
            continue
        base = os.path.splitext(os.path.basename(f))[0]
        if isinstance(obj, trimesh.Scene):
            for noeud in obj.graph.nodes_geometry:
                T, geom = obj.graph[noeud]
                m = obj.geometry[geom]
                if not isinstance(m, trimesh.Trimesh) or len(m.faces) == 0:
                    continue
                m = m.copy()
                if T is not None:
                    m.apply_transform(T)
                nom = _nom_noeud(obj, noeud, geom)
                if len(obj.graph.nodes_geometry) == 1 and len(fichiers) > 1:
                    nom = base                 # un GLB par structure : le nom du fichier fait foi
                col, fc = _couleur_visuel(m)
                part = {"nom": nom, "mesh": m, "couleur": col, "couleurs_faces": fc, "fichier": f}
                mh = re.match(r"^(.*?)~([0-9a-fA-F]{6})(\.\d+)?$", nom)      # banque : « structure~rrggbb » = surface déjà scindée par couleur
                if mh:
                    rgb = [int(mh.group(2)[i:i + 2], 16) for i in (0, 2, 4)]
                    part.update({"nom": mh.group(1) + "__c" + mh.group(2).lower(), "source_partie": mh.group(1), "couleur_scindee": rgb, "couleur": rgb})
                parties.append(part)
        elif isinstance(obj, trimesh.Trimesh):
            if len(obj.faces) == 0:
                continue
            col, fc = _couleur_visuel(obj)
            parties.append({"nom": base, "mesh": obj, "couleur": col, "couleurs_faces": fc, "fichier": f})
    return parties


_CLES = {"sommets": ["positions", "vertices", "sommets", "points", "verts", "xyz", "v", "pos"],
         "faces": ["faces", "triangles", "indices", "tris", "f", "tri"],
         "partie": ["face_part", "face_parts", "part", "parts", "partie", "parties", "face_partie", "objet", "objets", "object", "objects",
                    "face_object", "labels", "label", "ids", "id", "part_id", "face_id", "groupe", "group", "face_label"],
         "noms": ["names", "noms", "part_names", "parts_names", "objects_names", "noms_parties", "labels_noms", "nom", "name", "parties_noms"],
         "couleurs": ["face_colors", "face_colours", "couleurs_faces", "couleurs", "colors", "colours", "couleur_face", "face_color", "couleur"],
         "sommet_partie": ["vertex_part", "vertex_parts", "sommet_partie", "vertex_object", "vertex_label", "vertex_id"]}


def charger_npz(chemin, scene_json=None):
    """npz de la session « PDF 3D » : clés reconnues au mieux (voir _CLES). Formats acceptés :
    (a) sommets (N,3) + faces (M,3) + partie par face (M,) [+ noms (P,)] [+ couleurs par face (M,3)] ;
    (b) sommets + faces + partie par sommet (N,) ;
    (c) un jeu de tableaux par partie : <nom>_vertices/<nom>_faces (ou <nom>/vertices…), ou une liste d'objets (dtype object).
    _scene.json (facultatif) : noms, couleurs, visibilité par nœud, sous forme {nom: {...}} ou [{"nom"|"name":…, "couleur"|"color":…, "visible":…}]."""
    z = np.load(chemin, allow_pickle=True)
    cles = list(z.files)
    lc = {k.lower(): k for k in cles}
    # format de la banque hikspoors_maastricht (decode_u3d_local.py) : FACESET_<nom>|v, |f, |c, |p, |fc
    faceset = [k for k in cles if k.endswith("|v") and "|" in k]
    if faceset:
        parties = []
        for kv in faceset:
            pref = kv[:-2]
            nom = pref[len("FACESET_"):] if pref.startswith("FACESET_") else pref
            kf = pref + "|f"
            if kf not in z:
                continue
            V = np.asarray(z[kv], float).reshape(-1, 3); F = np.asarray(z[kf]).reshape(-1, 3).astype(np.int64)
            if len(F) == 0 or len(V) == 0:
                continue
            c = z[pref + "|c"] if pref + "|c" in z else None
            fc = np.asarray(z[pref + "|fc"]) if pref + "|fc" in z else None
            if fc is not None and (fc.ndim != 2 or len(fc) != len(F)):
                fc = None
            par = z[pref + "|p"] if pref + "|p" in z else None
            par = str(np.asarray(par).ravel()[0]) if par is not None and np.asarray(par).size else ""
            col = couleur_255(np.median(fc[:, :3], axis=0)) if fc is not None else (couleur_255(np.asarray(c, float).ravel()[:3]) if c is not None else None)
            parties.append({"nom": nom, "mesh": trimesh.Trimesh(V, F, process=False), "couleur": col, "couleurs_faces": fc[:, :3] if fc is not None else None,
                            "fichier": chemin, "parent": par, "couleur_materiau": couleur_255(np.asarray(c, float).ravel()[:3]) if c is not None else None})
        if parties:
            _completer_scene(parties, scene_json)
            return parties

    def trouve(cat):
        for k in _CLES[cat]:
            if k in lc:
                return lc[k]
        return None
    parties = []
    ks, kf = trouve("sommets"), trouve("faces")
    if ks and kf:
        V = np.asarray(z[ks], float).reshape(-1, 3)
        F = np.asarray(z[kf]).reshape(-1, 3).astype(np.int64)
        kp, kn, kc, kvp = trouve("partie"), trouve("noms"), trouve("couleurs"), trouve("sommet_partie")
        noms = [str(x) for x in np.asarray(z[kn]).ravel()] if kn else None
        FC = np.asarray(z[kc]) if kc else None
        if FC is not None and FC.ndim == 2 and len(FC) == len(V) and len(V) != len(F):
            FC = FC[F].mean(axis=1)
        if kp and len(np.asarray(z[kp]).ravel()) == len(F):
            P = np.asarray(z[kp]).ravel()
        elif kvp and len(np.asarray(z[kvp]).ravel()) == len(V):
            P = np.asarray(z[kvp]).ravel()[F[:, 0]]
        elif kp and len(np.asarray(z[kp]).ravel()) == len(V):
            P = np.asarray(z[kp]).ravel()[F[:, 0]]
        else:
            P = np.zeros(len(F), int)
        for i, p in enumerate(np.unique(P)):
            sel = P == p
            fs = F[sel]
            idx, inv = np.unique(fs.ravel(), return_inverse=True)
            m = trimesh.Trimesh(V[idx], inv.reshape(-1, 3), process=False)
            nom = noms[int(p)] if noms is not None and isinstance(p, (int, np.integer)) and int(p) < len(noms) else str(p)
            col = couleur_255(np.median(FC[sel][:, :3], axis=0)) if FC is not None and FC.ndim == 2 else None
            parties.append({"nom": nom, "mesh": m, "couleur": col, "couleurs_faces": FC[sel][:, :3] if FC is not None and FC.ndim == 2 else None,
                            "fichier": chemin})
    else:
        groupes = {}
        for k in cles:
            m = re.match(r"^(.*?)[/_](vertices|verts|positions|sommets|faces|triangles|indices|colors|couleurs|face_colors)$", k, re.I)
            if m:
                groupes.setdefault(m.group(1), {})[m.group(2).lower()] = k
        for nom, g in groupes.items():
            ks2 = next((g[a] for a in ("vertices", "verts", "positions", "sommets") if a in g), None)
            kf2 = next((g[a] for a in ("faces", "triangles", "indices") if a in g), None)
            if not ks2 or not kf2:
                continue
            V = np.asarray(z[ks2], float).reshape(-1, 3); F = np.asarray(z[kf2]).reshape(-1, 3).astype(np.int64)
            kc2 = next((g[a] for a in ("face_colors", "colors", "couleurs") if a in g), None)
            FC = np.asarray(z[kc2]) if kc2 else None
            parties.append({"nom": nom, "mesh": trimesh.Trimesh(V, F, process=False), "couleur": couleur_255(np.median(FC[:, :3], 0)) if FC is not None else None,
                            "couleurs_faces": FC[:, :3] if FC is not None else None, "fichier": chemin})
        if not parties:
            for k in cles:                     # liste d'objets (dtype object) : dicts {nom, vertices, faces…}
                arr = z[k]
                if arr.dtype == object:
                    for o in np.atleast_1d(arr).ravel():
                        if isinstance(o, dict):
                            V = o.get("vertices", o.get("positions", o.get("sommets")))
                            F = o.get("faces", o.get("triangles"))
                            if V is None or F is None:
                                continue
                            parties.append({"nom": str(o.get("nom", o.get("name", k))), "mesh": trimesh.Trimesh(np.asarray(V, float).reshape(-1, 3), np.asarray(F).reshape(-1, 3), process=False),
                                            "couleur": couleur_255(o.get("couleur", o.get("color"))), "couleurs_faces": None, "fichier": chemin})
    if not parties:
        raise SystemExit("npz non reconnu ; clés présentes : %s\nAttendu : sommets (N,3) + faces (M,3) + partie par face (M,) + noms (P,), voir charger_npz()." % cles)
    _completer_scene(parties, scene_json)
    return parties


def _infos_scene(scene_json):
    """_scene.json (facultatif) : {nom: {...}} ou [{"nom"|"name":…, "couleur"|"color":…, "visible":…}] ou {"noeuds"|"nodes"|…: …}"""
    infos = {}
    if scene_json and os.path.exists(scene_json):
        try:
            sc = json.load(open(scene_json, encoding="utf-8"))
        except Exception as e:
            print("   _scene.json illisible :", e); return infos
        items = sc
        if isinstance(sc, dict):
            for k in ("noeuds", "nodes", "parties", "parts", "objets", "structures", "facesets", "meshes"):
                if isinstance(sc.get(k), (dict, list)):
                    items = sc[k]; break
        if isinstance(items, dict):
            for k, v in items.items():
                if isinstance(v, dict):
                    infos[normaliser(k)] = v
        elif isinstance(items, list):
            for v in items:
                if isinstance(v, dict):
                    n = v.get("nom") or v.get("name") or v.get("node") or v.get("noeud")
                    if n:
                        infos[normaliser(n)] = v
    return infos


def _completer_scene(parties, scene_json):
    infos = _infos_scene(scene_json)
    for p in parties:
        inf = infos.get(normaliser(p["nom"])) or infos.get(normaliser("FACESET_" + p["nom"]))
        if inf:
            if p.get("couleur") is None:
                c = inf.get("couleur") or inf.get("color") or inf.get("diffuse") or inf.get("rgb")
                p["couleur"] = couleur_255(c) if c is not None else None
            p["visible"] = inf.get("visible", inf.get("V", True))
            if not p.get("parent") and (inf.get("parent") or inf.get("p")):
                p["parent"] = str(inf.get("parent") or inf.get("p"))


# ----------------------------------------------------------------------------------------------------------------------- nomenclature
class Nomenclature:
    def __init__(self, chemin=NOMENCLATURE):
        self.t = json.load(open(chemin, encoding="utf-8"))
        self.ignorer = [re.compile(m, re.I) for m in self.t.get("ignorer", [])]
        self.cote_g = re.compile(self.t["cotes"]["gauche"], re.I)
        self.cote_d = re.compile(self.t["cotes"]["droit"], re.I)

        def regle(r, prefixe=""):
            return {"motif": re.compile(r["motif"], re.I), "exige": None, "nom": r["nom"], "nom_fr": r.get("nom_fr", r["nom"]),
                    "systeme": r.get("systeme", "autre"), "groupe": r.get("groupe"), "couleur": r.get("couleur"),
                    "confiance": r.get("confiance", "moyenne"), "capture": re.compile(r["capture"], re.I) if r.get("capture") else None,
                    "regle": prefixe + r["motif"]}
        # noms réels qui contiennent « card », « atri », « ventric »… (cardinal_vein, pericard, primary_atrial_septum) : avant les règles cardiaques
        self.regles = [regle(r, "prioritaire : ") for r in self.t.get("prioritaires", [])]
        card = self.t.get("cardiaque", {})
        for ty in card.get("types", []):
            for ch in card.get("chambres", []):
                self.regles.append({"motif": re.compile(ty["motif"], re.I), "exige": re.compile(ch["motif"], re.I),
                                    "nom": ty["prefixe"] + "_" + ch["suffixe"], "nom_fr": ty["nom_fr"] + " " + ch["nom_fr"], "systeme": "vasculaire",
                                    "groupe": ty.get("groupe"), "couleur": ty.get("couleur"), "confiance": ty.get("confiance", "bonne"), "capture": None,
                                    "regle": "cardiaque:%s×%s" % (ty["prefixe"], ch["suffixe"])})
        cs_ = card.get("chambre_seule")
        if cs_:
            for ch in card.get("chambres", []):
                self.regles.append({"motif": re.compile(ch["motif"], re.I), "exige": None, "nom": cs_["prefixe"] + "_" + ch["suffixe"],
                                    "nom_fr": cs_["nom_fr"] + " " + ch["nom_fr"], "systeme": "vasculaire", "groupe": cs_.get("groupe"), "couleur": cs_.get("couleur"),
                                    "confiance": cs_.get("confiance", "bonne"), "capture": None, "regle": "cardiaque:chambre seule %s" % ch["suffixe"]})
        for r in self.t.get("regles", []):
            self.regles.append(regle(r))
        self.par_couleur = [dict(r, _re=re.compile(r.get("source", "."), re.I), _cap=re.compile(r["capture"], re.I) if r.get("capture") else None)
                            for r in self.t.get("par_couleur", []) if "rgb" in r]
        self.groupes = self.t.get("groupes", {})

    def cote(self, n):
        g, d = bool(self.cote_g.search(n)), bool(self.cote_d.search(n))
        if g and not d:
            return "gauche"
        if d and not g:
            return "droit"
        return ""

    @staticmethod
    def _remplir(txt, cote, cap, fr=False):
        if fr:
            mot = txt.split(" ")[0].lower() if txt else ""
            cf = ("droite" if mot in FEMININ else "droit") if cote == "droit" else cote
            txt = txt.replace("{cote_fr}", cf)
        txt = txt.replace("{cote}", cote).replace("{1}", cap)
        txt = re.sub(r"\s+", " ", txt).strip() if fr else re.sub(r"_+", "_", txt).strip("_")
        if not fr:
            txt = re.sub(r"^([a-z0-9]+)_\1$", r"\1", txt)          # « coeur_coeur » (chambre seule « heart ») -> « coeur »
        return txt

    def chercher(self, nom_source):
        n = normaliser(nom_source)
        for ig in self.ignorer:
            if ig.search(n):
                return {"ignorer": True}
        cote = self.cote(n)
        for r in self.regles:
            if r["motif"].search(n) and (r["exige"] is None or r["exige"].search(n)):
                cap = ""
                if r["capture"] is not None:
                    m = r["capture"].search(n)
                    cap = normaliser(m.group(1)) if m and m.groups() and m.group(1) else ""
                nom = self._remplir(r["nom"], cote, cap) or slug(nom_source)
                return {"nom": nom, "nom_fr": self._remplir(r["nom_fr"], cote, cap.replace("_", " "), fr=True), "systeme": r["systeme"],
                        "groupe": r["groupe"], "couleur": r["couleur"], "confiance": r["confiance"], "regle": r["regle"], "non_reconnu": False}
        nom = slug(nom_source)
        return {"nom": nom, "nom_fr": str(nom_source), "systeme": reclasser(nom, "autre"), "groupe": None, "couleur": None,
                "confiance": "moyenne", "regle": None, "non_reconnu": True}

    def par_couleur_regle(self, nom_source, rgb):
        """règle « par couleur » (structure scindée par couleur de texture) ; renvoie un dict prêt pour le manifest ou None"""
        n = normaliser(nom_source)
        for r in self.par_couleur:
            if r["_re"].search(n) and rgb is not None and np.abs(np.asarray(rgb, float) - np.asarray(r["rgb"], float)).max() <= r.get("tol", 30):
                cap = ""
                if r["_cap"] is not None:
                    m = r["_cap"].search(n)
                    cap = normaliser(m.group(1)) if m else ""
                return {"nom": self._remplir(r["nom"], self.cote(n), cap), "nom_fr": self._remplir(r.get("nom_fr", r["nom"]), self.cote(n), cap.upper(), fr=True),
                        "systeme": r.get("systeme", "autre"), "groupe": r.get("groupe"), "couleur": r.get("couleur"), "confiance": r.get("confiance", "moyenne"),
                        "regle": "par_couleur %s" % r["rgb"], "non_reconnu": False}
        return None


def scinder_par_couleur(p, seuil=0.05):
    """Une partie portant plusieurs couleurs de faces (texture des auteurs : deux structures sur une surface) -> sous-parties."""
    fc = p.get("couleurs_faces")
    if fc is None or len(fc) != len(p["mesh"].faces):
        return [p]
    q = (np.asarray(fc, float)[:, :3] // 24).astype(int)
    cles, inv, cnt = np.unique(q, axis=0, return_inverse=True, return_counts=True)
    ordre = np.argsort(-cnt)
    gros = [k for k in ordre if cnt[k] >= seuil * len(fc)]
    if len(gros) < 2:
        return [p]
    F = np.asarray(p["mesh"].faces); V = np.asarray(p["mesh"].vertices)
    reste = np.ones(len(F), bool)
    out = []
    for k in gros:
        sel = inv == k
        reste &= ~sel
        idx, inv2 = np.unique(F[sel].ravel(), return_inverse=True)
        m = trimesh.Trimesh(V[idx], inv2.reshape(-1, 3), process=False)
        col = couleur_255(np.median(np.asarray(fc)[sel][:, :3], axis=0))
        out.append({"nom": "%s__c%02x%02x%02x" % (p["nom"], *col), "mesh": m, "couleur": col, "couleurs_faces": None, "fichier": p["fichier"],
                    "source_partie": p["nom"], "couleur_scindee": col, "fraction_faces": round(float(sel.mean()), 4),
                    "parent": p.get("parent", ""), "complement": p.get("complement")})
    if reste.any():                            # petites couleurs restantes : rattachées à la plus grosse sous-partie
        idx, inv2 = np.unique(F[reste].ravel(), return_inverse=True)
        out[0]["mesh"] = trimesh.util.concatenate([out[0]["mesh"], trimesh.Trimesh(V[idx], inv2.reshape(-1, 3), process=False)])
    return out


# ----------------------------------------------------------------------------------------------------------------------- orientation
def _echantillon(m, n=4000):
    V = np.asarray(m.vertices)
    if len(V) <= n:
        return V
    return V[np.random.RandomState(0).choice(len(V), n, replace=False)]


def _centre(structs, motifs):
    pts = []
    for nom, s in structs.items():
        if any(re.search(mo, nom) for mo in motifs):
            pts.append(np.vstack([_echantillon(m) for m in s["meshes"]]).mean(0))
    return np.mean(pts, 0) if pts else None


STADE_COURANT = None      # numéro du stade en cours (fixé par main) : certains indices ne valent qu'à partir d'un stade
CRANIAL = [   # (A crânial à B) : motifs sur les noms canoniques, poids, stade minimal où l'indice est valable
    # cœur au-dessus du foie / septum transversum : FAUX avant la bascule de la tête (CS9-CS10 : le septum transversum est encore crânial au
    # croissant cardiaque ; fausse alerte CS9 du relais, 26/09) -> compté à partir de CS11
    (["^myocarde_", "^cavite_", "^coeur$", "^tube_cardiaque", "^myocarde_coeur"], ["^foie$", "^septum_transversum", "^intestin_moyen", "^intestin_posterieur", "^veine_ombilicale", "^artere_ombilicale", "^cordon"], 2.0, 11),
    (["^arc_aortique", "^arcs_aortiques", "^sac_aortique", "^carotide", "^pharynx", "^encephale", "^yeux", "^vesicules_otiques"], ["^myocarde_", "^cavite_", "^coeur$", "^foie$", "^intestin"], 2.0, 0),
    (["^myocarde_voie_efferente", "^cavite_voie_efferente", "^myocarde_ventricule", "^cavite_ventricule"], ["^myocarde_sinus_veineux", "^cavite_sinus_veineux", "^veine_vitelline", "^canal_hepatocardiaque", "^veine_cave_inferieure"], 1.0, 0),
    (["^intestin_anterieur", "^oesophage", "^estomac", "^poumons"], ["^intestin_posterieur", "^veine_ombilicale", "^artere_ombilicale", "^cordon"], 1.0, 0),
]
# Dorsal : structures de la paroi dorsale. Le tube neural n'est qu'un repli : dans un embryon en C (CS12-CS17) il fait le tour du cœur
# (prosencéphale ventral contre le cœur, queue relevée), son barycentre n'est pas dorsal (retour du relais, CS13).
DORSAUX = ["^aorte_dorsale", "^somites", "^notochorde", "^ganglions_spinaux", "^veine_cardinale(?!_commune)"]
DORSAUX_REPLI = ["^tube_neural", "^encephale"]
CARDIAQUES = ["^myocarde_", "^cavite_", "^coeur$", "^tube_cardiaque", "^pericarde", "^epicarde", "^gelee_cardiaque"]
VENTRAUX = CARDIAQUES + ["^foie$", "^septum_transversum", "^vesicule_vitelline"]


def _points_par_structure(structs, motifs, n=3000):
    return [np.vstack([_echantillon(m, n) for m in s["meshes"]]) for nom, s in structs.items() if any(re.search(mo, nom) for mo in motifs)]


def _dorsal(structs, Z):
    """Vecteur ventral -> dorsal (⊥ Z, non normé) : barycentre des structures dorsales, chacune restreinte à la tranche crânio-caudale du
    cœur, moins celui des structures cardiaques (une voix par structure) ; barycentres entiers si la tranche est vide.
    Renvoie (vecteur | None, description)."""
    Gd, lab = _points_par_structure(structs, DORSAUX), "aorte/somites/cardinales"
    if not Gd:
        Gd, lab = _points_par_structure(structs, DORSAUX_REPLI), "tube neural (repli)"
    Gv = _points_par_structure(structs, CARDIAQUES) or _points_par_structure(structs, VENTRAUX)
    if not Gd or not Gv:
        return None, "pas d'indice dorsal"
    lo, hi = np.percentile(np.vstack(Gv) @ Z, [5, 95]); marge = 0.15 * (hi - lo)
    cd = []
    for P in Gd:
        z = P @ Z
        sel = (z >= lo - marge) & (z <= hi + marge)
        if sel.sum() >= 10:
            cd.append(P[sel].mean(0))
    if cd:
        lab += ", tranche du cœur"
    else:
        cd = [P.mean(0) for P in Gd]; lab += ", barycentres"
    d = np.mean(cd, 0) - np.mean([P.mean(0) for P in Gv], 0)
    d = d - (d @ Z) * Z
    return (d if np.linalg.norm(d) > 1e-12 else None), lab


def _axe_anatomique(structs):
    """axe caudal -> crânial d'après les paires CRANIAL (foie -> cœur -> arcs aortiques…), ou None"""
    v = np.zeros(3)
    for a, b, w, smin in CRANIAL:
        if STADE_COURANT is not None and STADE_COURANT < smin:
            continue
        ca, cb = _centre(structs, a), _centre(structs, b)
        if ca is not None and cb is not None and np.linalg.norm(ca - cb) > 1e-12:
            v += w * (ca - cb) / np.linalg.norm(ca - cb)
    return v / np.linalg.norm(v) if np.linalg.norm(v) > 1e-9 else None


def _score(structs, R):
    """score d'un repère R (lignes = X gauche, Y dorsal, Z crânial dans le repère source) : somme des indices vérifiés, détail par indice"""
    def c(motifs):
        v = _centre(structs, motifs)
        return None if v is None else R @ v
    detail = []; total = 0.0
    for a, b, w, smin in CRANIAL:
        if STADE_COURANT is not None and STADE_COURANT < smin:
            continue
        ca, cb = c(a), c(b)
        if ca is None or cb is None:
            continue
        ok = ca[2] > cb[2]
        total += w if ok else -w
        detail.append({"indice": "crânial : %s > %s" % (a[0], b[0]), "ok": bool(ok), "ecart_mm": None, "poids": w})
    d, lab = _dorsal(structs, R[2])
    if d is not None:
        ok = float(R[1] @ d) > 0
        total += 2.0 if ok else -2.0
        detail.append({"indice": "dorsal : %s au-dessus (Y) du cœur" % lab, "ok": bool(ok), "poids": 2.0})
    # tube neural plus large du côté crânial : information seulement (faux sur l'embryon en C de CS13, retour du relais)
    if "tube_neural" in structs:
        P = R @ np.vstack([_echantillon(m, 20000) for m in structs["tube_neural"]["meshes"]]).T
        z = P[2]; lo, hi = np.percentile(z, [15, 85])
        def rayon(sel):
            Q = P[:2, sel]
            return float(np.sqrt(((Q - Q.mean(1, keepdims=True)) ** 2).sum(0)).mean()) if sel.sum() > 20 else None
        r_haut, r_bas = rayon(z >= hi), rayon(z <= lo)
        if r_haut and r_bas and abs(r_haut - r_bas) > 0.15 * max(r_haut, r_bas):
            detail.append({"indice": "tube neural plus large en haut (encéphale ; information, non compté)", "ok": bool(r_haut > r_bas), "poids": 0.0,
                           "rayon_haut": r_haut, "rayon_bas": r_bas})
    # gauche/droite : structures appariées (…_gauche / …_droit ou …_droite)
    dxs = []
    for nom in structs:
        if not nom.endswith("_gauche"):
            continue
        dr = next((nom[:-7] + s for s in ("_droit", "_droite") if nom[:-7] + s in structs), None)
        if dr:
            g = R @ np.vstack([_echantillon(m) for m in structs[nom]["meshes"]]).mean(0)
            d_ = R @ np.vstack([_echantillon(m) for m in structs[dr]["meshes"]]).mean(0)
            dxs.append(float(g[0] - d_[0]))
    if dxs:
        ok = np.median(dxs) > 0
        total += 2.0 if ok else -2.0
        detail.append({"indice": "gauche à +X (%d paires)" % len(dxs), "ok": bool(ok), "poids": 2.0, "dx_median": float(np.median(dxs))})
    return total, detail


def orienter(structs):
    """Renvoie le meilleur candidat {R (3x3, lignes X/Y/Z dans le repère source), score, detail, z, y} et tous les candidats.
    Z = ±e1, ±e2, ±e3 (ACP) et ±axe anatomique ; Y = indice dorsal ⊥ Z (sinon ±axe ACP) ; X = Y×Z. Meilleur score, e1 d'abord à égalité."""
    P = np.vstack([_echantillon(m) for s in structs.values() for m in s["meshes"]])
    _, _, Vt = np.linalg.svd(P - P.mean(0), full_matrices=False)
    axes = [("ACP e1", Vt[0]), ("ACP e2", Vt[1]), ("ACP e3", Vt[2])]
    a = _axe_anatomique(structs)
    if a is not None:
        axes.append(("axe anatomique", a))
    cand = []
    for lab_z, e in axes:
        for sz in (1, -1):
            Z = sz * e / np.linalg.norm(e)
            d, lab_d = _dorsal(structs, Z)
            ys = [("indice dorsal (%s)" % lab_d, d / np.linalg.norm(d))] if d is not None else []
            if not ys:
                autre = Vt[1] if lab_z == "ACP e1" else Vt[0]
                for sy in (1, -1):
                    y = sy * autre - ((sy * autre) @ Z) * Z
                    if np.linalg.norm(y) > 1e-9:
                        ys.append(("ACP (%+d)" % sy, y / np.linalg.norm(y)))
            for lab, Y in ys:
                R = np.vstack([np.cross(Y, Z), Y, Z])
                sc, det = _score(structs, R)
                cand.append({"R": R, "score": sc, "detail": det, "z": "%s (%+d)" % (lab_z, sz), "y": lab})
    cand.sort(key=lambda c: -c["score"])          # tri stable : à score égal, l'ordre ci-dessus (e1 d'abord)
    return cand[0], cand


VARIANTES = [("tel quel", np.eye(3)), ("rotation 180° autour de X", np.diag([1., -1, -1])), ("rotation 180° autour de Y", np.diag([-1., 1, -1])),
             ("rotation 180° autour de Z", np.diag([-1., -1, 1])), ("miroir X", np.diag([-1., 1, 1]))]


# ----------------------------------------------------------------------------------------------------------------------- unions, rendu
def union_groupe(meshes, faces_max=60000):
    """union « blob » de plusieurs maillages : concaténation, puis re-maillage par voxels (surface fermée unique) si possible"""
    cat = trimesh.util.concatenate([m.copy() for m in meshes])
    try:
        diag = float(np.linalg.norm(cat.bounds[1] - cat.bounds[0]))
        pitch = max(diag / 160.0, 1e-6)
        vg = cat.voxelized(pitch).fill()
        u = vg.marching_cubes
        u.apply_transform(vg.transform)          # marching_cubes de trimesh rend des coordonnées voxel : retour en mm
        if len(u.faces) > 0:
            trimesh.smoothing.filter_laplacian(u, iterations=4)
            if len(u.faces) > faces_max:
                try:
                    u = u.simplify_quadric_decimation(faces_max)
                except Exception:
                    pass
            return u, "voxels (pas %.4g)" % pitch
    except Exception as e:
        print("   union par voxels impossible :", e)
    return cat, "concaténation"


def planche_controle(structs_mm, chemin, titre):
    """3 vues orthographiques (profil gauche, face, dos) par splat de points avec tampon de profondeur ; légende par structure."""
    try:
        from PIL import Image, ImageDraw
    except Exception:
        print("   PIL absent : pas de planche de contrôle"); return
    H = 640; marge = 30
    pts = []
    for nom, s in structs_mm.items():
        col = s.get("couleur") or [0.6, 0.6, 0.6]
        for m in s["meshes"]:
            n = int(min(60000, max(3000, m.area * 2500)))
            Pp, fi = trimesh.sample.sample_surface(m, n, seed=1)
            pts.append((Pp, m.face_normals[fi], np.array(col)))
    if not pts:
        return
    allP = np.vstack([p[0] for p in pts]); c = (allP.max(0) + allP.min(0)) / 2; ext = float((allP.max(0) - allP.min(0)).max())
    sc = ext / (H - 2 * marge)
    vues = [("profil gauche (ventral ← → dorsal)", np.array([0, 1, 0.]), np.array([0, 0, 1.]), np.array([-1, 0, 0.])),
            ("face (gauche de l'embryon à droite de l'image)", np.array([-1, 0, 0.]), np.array([0, 0, 1.]), np.array([0, -1, 0.])),
            ("dos (gauche de l'embryon à gauche de l'image)", np.array([1, 0, 0.]), np.array([0, 0, 1.]), np.array([0, 1, 0.]))]
    W = H
    img = Image.new("RGB", (3 * W + 2 * 10, H + 170), (250, 250, 247)); dr = ImageDraw.Draw(img)
    police = None
    for f in ("arial.ttf", "DejaVuSans.ttf", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", "C:/Windows/Fonts/arial.ttf"):
        try:
            from PIL import ImageFont
            police = ImageFont.truetype(f, 13); break
        except Exception:
            continue
    if police is None:      # police bitmap par défaut : pas d'accents
        vues = [(normaliser(l).replace("_", " "), u, v, w) for l, u, v, w in vues]; titre = normaliser(titre).replace("_", " ")
    def texte(xy, t, fill):
        dr.text(xy, t, fill=fill, font=police)
    for k, (lab, u, v, w) in enumerate(vues):
        buf = np.zeros((H, W, 3), np.uint8) + np.array([250, 250, 247], np.uint8); zb = np.full((H, W), -np.inf)
        for Pp, N, col in pts:
            x = ((Pp - c) @ u) / sc + W / 2; y = H / 2 - ((Pp - c) @ v) / sc; z = (Pp - c) @ w
            sh = 0.35 + 0.65 * np.abs(N @ w)
            rgb = np.clip(col[None, :] * 255 * sh[:, None], 0, 255).astype(np.uint8)
            o = np.argsort(z)
            xi = x[o].astype(int); yi = y[o].astype(int); ok = (xi >= 0) & (xi < W) & (yi >= 0) & (yi < H)
            xi, yi, zz, cc = xi[ok], yi[ok], z[o][ok], rgb[o][ok]
            keep = zz >= zb[yi, xi]
            buf[yi[keep], xi[keep]] = cc[keep]; zb[yi[keep], xi[keep]] = zz[keep]
        img.paste(Image.fromarray(buf), (k * (W + 10), 0))
        texte((k * (W + 10) + 8, 6), lab, (30, 30, 30))
        texte((k * (W + 10) + 8, H - 18), "cranial en haut - barre = 1 mm", (60, 60, 60))
        dr.line([(k * (W + 10) + 8, H - 26), (k * (W + 10) + 8 + int(1.0 / sc), H - 26)], fill=(0, 0, 0), width=3)
    texte((8, H + 8), titre, (20, 20, 20))
    x, y = 8, H + 28
    for nom, s in sorted(structs_mm.items()):
        col = tuple(int(255 * v) for v in (s.get("couleur") or [0.6, 0.6, 0.6]))
        dr.rectangle([x, y + 2, x + 10, y + 12], fill=col); texte((x + 14, y), nom[:34], (30, 30, 30))
        y += 14
        if y > H + 150:
            y = H + 28; x += 250
            if x > 3 * W - 240:
                break
    img.save(chemin)


def verifier_sortie(sortie, manifest):
    """verif_maillages.json dans le dossier de sortie, où qu'il soit (l'agrégateur le lit pour embryons_3D/modeles/<dossier>/)"""
    try:
        import verif_maillages
    except Exception as e:
        print("   verif_maillages :", e); return None
    out = {}
    for x in manifest["structures"]:
        f = os.path.join(sortie, x["fichier"])
        try:
            out[x["nom"]] = verif_maillages.verifier(verif_maillages.charger(f)) if os.path.exists(f) else {"verdict": "absent", "problemes": ["fichier absent"]}
        except Exception as e:
            out[x["nom"]] = {"verdict": "a_corriger", "problemes": ["lecture : %s" % e]}
    ok = sum(1 for v in out.values() if v["verdict"] == "ok")
    res = {"dossier": os.path.basename(os.path.normpath(sortie)), "structures": out, "bilan": {"ok": ok, "total": len(out)}}
    json.dump(res, open(os.path.join(sortie, "verif_maillages.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("   verif_maillages : %d/%d maillages sans erreur" % (ok, len(out)))
    for n, v in out.items():
        if v["verdict"] != "ok":
            print("      %-40s %s" % (n, " ; ".join(v["problemes"])))
    return res["bilan"]


# ----------------------------------------------------------------------------------------------------------------------- principal
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stade")
    ap.add_argument("--glb", help="dossier (ou fichier) de maillages GLB/PLY/OBJ/STL")
    ap.add_argument("--npz"); ap.add_argument("--scene")
    ap.add_argument("--banque", default=BANQUE, help="banque hikspoors_maastricht (u3d/Carnegie_Stage_<n>[_NEWvalves].npz, calage/*.json)")
    ap.add_argument("--sans-newvalves", action="store_true", help="prendre Carnegie_Stage_<n>.npz même si la version _NEWvalves existe")
    ap.add_argument("--sans-completer", action="store_true", help="_NEWvalves : ne pas reprendre de l'originale les structures absentes")
    ap.add_argument("--um-par-unite", type=float, default=None)
    ap.add_argument("--echelles", default=ECHELLES)
    ap.add_argument("--calage", help="calage explicite (.json matrice_pdf_vers_mm, .npy) ; défaut : calage/<nom du npz>.json de la banque")
    ap.add_argument("--calage-repere", default="modeles", choices=["modeles", "pipeline"],
                    help="repère d'arrivée du calage : modeles (défaut, calages de la banque) ; pipeline = X gauche→droite (indirect), remis direct")
    ap.add_argument("--sans-calage", action="store_true", help="ignorer le calage de la banque : orientation automatique")
    ap.add_argument("--recentrer", action="store_true", help="avec un calage vers les modèles : recentrer quand même (perd la superposition)")
    ap.add_argument("--retourner", default="", help="axes à inverser après l'orientation auto : z, y ou zy")
    ap.add_argument("--separer-couleurs", action="store_true")
    ap.add_argument("--sans-unions", action="store_true")
    ap.add_argument("--publiable", action="store_true", help="autoriser la publication sur le site (défaut : référence interne)")
    ap.add_argument("--non-publiable", action="store_true", help=argparse.SUPPRESS)       # ancien interrupteur : c'est maintenant le défaut
    ap.add_argument("--nomenclature", default=NOMENCLATURE)
    ap.add_argument("--sortie"); ap.add_argument("--session", default="Hikspoors → modèles (cloud, 26/09)")
    ap.add_argument("--faces-max", type=int, default=0, help="décimation des structures au-delà de N faces (0 = aucune)")
    a = ap.parse_args()
    for k in ("glb", "npz", "scene", "calage", "sortie", "banque"):
        setattr(a, k, abs0(getattr(a, k)))
    if not os.path.exists(a.echelles):
        a.echelles = abs0(a.echelles)
    if not os.path.exists(a.nomenclature):
        a.nomenclature = abs0(a.nomenclature)
    cs = stade_norme(a.stade)
    global STADE_COURANT
    STADE_COURANT = norm_int(cs)
    sortie = a.sortie or os.path.join("embryons_3D", "modeles", cs + "_hikspoors")
    os.makedirs(sortie, exist_ok=True)
    nomen = Nomenclature(a.nomenclature)
    rapport = {"stade": cs, "source": SOURCE}

    # 1. lecture (+ structures de la version d'origine absentes de la _NEWvalves)
    original = None
    if not a.npz and not a.glb:
        tr = trouver_dans_banque(a.banque, cs, not a.sans_newvalves)
        if tr is None:
            ap.error("stade %s introuvable dans la banque %s (u3d/Carnegie_Stage_%d[_NEWvalves].npz ou glb/…) ; sinon --npz ou --glb" % (cs, a.banque, norm_int(cs)))
        a.npz, a.glb, a.scene, original = tr["npz"], tr["glb"], a.scene or tr["scene"], tr["original"]
        print("   banque : %s" % os.path.basename(a.npz or a.glb))
    fichier_source = a.npz or a.glb
    base = re.sub(r"_scene$", "", os.path.splitext(os.path.basename(fichier_source))[0]) if os.path.isfile(fichier_source) else ""
    rapport["fichier_source"] = fichier_source
    parties = charger_npz(a.npz, a.scene) if a.npz else charger_fichiers(a.glb, cs)
    print("%s : %d parties lues" % (cs, len(parties)))
    meme_rep = None
    if original:
        parties_o = charger_npz(original["npz"], original["scene"])
        meme_rep, txt = meme_repere(parties, parties_o)
        rapport["version_origine"] = {"fichier": original["npz"], "meme_repere": meme_rep, "controle": txt}
        if a.sans_completer:
            pass
        elif meme_rep:
            ajout = completer_depuis_original(parties, parties_o, nomen)
            parties += ajout
            rapport["version_origine"]["ajoutees"] = [p["nom"] for p in ajout]
            print("   version d'origine (%s) : %s ; ajoutées : %s" % (os.path.basename(original["npz"]), txt, [p["nom"] for p in ajout] or "aucune"))
        else:
            print("   ATTENTION version d'origine %s : repère différent ou non vérifiable (%s) : rien n'est repris" % (os.path.basename(original["npz"]), txt))
    if a.separer_couleurs:
        parties = [q for p in parties for q in scinder_par_couleur(p)]
        print("   après scission par couleur : %d parties" % len(parties))

    # 2. nomenclature -> structures canoniques (fusion des parties de même nom)
    structs, ignores, non_reconnus, journal = {}, [], [], []
    for p in parties:
        r = nomen.par_couleur_regle(p.get("source_partie", p["nom"]), p.get("couleur_scindee")) if p.get("couleur_scindee") is not None else None
        if r is None:
            r = nomen.chercher(p.get("source_partie", p["nom"]))
        if r.get("non_reconnu") and p.get("parent") and normaliser(p["parent"]) != normaliser(p.get("source_partie", p["nom"])):
            r2 = nomen.chercher(p["parent"] + "_" + p.get("source_partie", p["nom"]))    # nœud inconnu sous un parent connu
            if not r2.get("non_reconnu") and not r2.get("ignorer"):
                r2["regle"] = "via parent « %s » : %s" % (p["parent"], r2["regle"]); r2["non_reconnu"] = True; r = r2   # à confirmer dans la table
            else:
                r3 = nomen.chercher(p["parent"])
                if not r3.get("non_reconnu") and not r3.get("ignorer"):     # système et groupe du parent, nom propre, à compléter dans la table
                    r = {"nom": r3["nom"] + "_" + slug(p.get("source_partie", p["nom"])), "nom_fr": r3["nom_fr"] + " : " + str(p["nom"]),
                         "systeme": r3["systeme"], "groupe": r3["groupe"], "couleur": r3["couleur"], "confiance": "moyenne",
                         "regle": "parent « %s » reconnu (%s), nom à compléter" % (p["parent"], r3["regle"]), "non_reconnu": True}
        if r.get("ignorer"):
            ignores.append(p["nom"]); continue
        if r["non_reconnu"]:
            non_reconnus.append(p["nom"])
        s = structs.setdefault(r["nom"], {"meshes": [], "sources": [], "meta": r, "couleurs_source": []})
        s["meshes"].append(p["mesh"]); s["sources"].append(p["nom"])
        if p.get("couleur") is not None:
            s["couleurs_source"].append(p["couleur"])
        j = {"source": p["nom"], "parent": p.get("parent", ""), "nom": r["nom"], "regle": r["regle"], "faces": int(len(p["mesh"].faces)), "couleur_source": p.get("couleur")}
        if p.get("complement"):
            j["complement"] = p["complement"]
        journal.append(j)
    if not structs:
        raise SystemExit("aucune structure retenue")
    print("   %d structures canoniques ; non reconnues : %s ; ignorées : %s" % (len(structs), non_reconnus or "aucune", ignores or "aucune"))
    rapport.update({"non_reconnus": non_reconnus, "ignores": ignores, "journal": journal})

    # 3. calage (banque par défaut) ou échelle + orientation automatique
    calage, note_calage, repere_calage = a.calage, "", a.calage_repere
    if not calage and not a.sans_calage:
        calage, note_calage = trouver_calage(a.banque, base)
        repere_calage = "modeles"
        if calage and note_calage and meme_rep is not True:
            print("   calage %s ignoré : %s, mais le même repère n'est pas vérifié" % (os.path.basename(calage), note_calage))
            calage = None
    reflexion = False
    if calage:
        A, t_, infos = lire_calage(calage)
        det_A = float(np.linalg.det(A))
        if det_A == 0:
            raise SystemExit("calage %s : matrice singulière" % calage)
        s_ = float(np.cbrt(abs(det_A)))
        miroir = repere_calage == "pipeline"
        M = np.diag([-1.0, 1, 1]) @ A if miroir else A
        tt = np.diag([-1.0, 1, 1]) @ t_ if miroir else t_
        reflexion = float(np.linalg.det(M)) < 0
        um = s_ * 1000.0
        k_mm = 1.0
        inter = None
        mc = re.search(r"CS\s*0*(\d+)", str(infos.get("cible", "")), re.I)
        if mc and int(mc.group(1)) != norm_int(cs):
            inter = "CS%d" % int(mc.group(1))
        def transformer(V):
            return V @ M.T + tt
        Rc = M / s_
        variantes = []
        for lab, D in VARIANTES:
            sc, det = _score(structs, D @ Rc)
            variantes.append({"variante": lab, "score": sc, "indices": det})
        sc_tel, det_tel = variantes[0]["score"], variantes[0]["indices"]
        meilleure = max(variantes, key=lambda v: v["score"])
        rapport["echelle"] = {"um_par_unite": um, "source": "calage %s" % os.path.basename(calage)}
        rapport["calage"] = dict(infos, repere=repere_calage, note=note_calage, echelle_um_par_unite=um, reflexion=reflexion, inter_stades=inter,
                                 recentre=bool(a.recentrer or repere_calage == "pipeline"))
        rapport["orientation"] = {"methode": "calage %s (repère %s%s)" % (os.path.basename(calage), repere_calage, ", " + note_calage if note_calage else ""),
                                  "score_indices": sc_tel, "indices": det_tel, "variantes": [{"variante": v["variante"], "score": v["score"]} for v in variantes]}
        print("   calage : %s (%s%s) ; %.4g µm/unité ; repère %s%s" % (os.path.basename(calage), infos.get("cle", "matrice"), ", " + note_calage if note_calage else "",
                                                                   um, repere_calage, " ; RÉFLEXION (det < 0) : faces réorientées, à confirmer" if reflexion else ""))
        if inter:
            print("   ATTENTION : calage vers le modèle %s, pas %s : l'échelle est celle de %s (à confirmer avec la session « PDF 3D »)" % (inter, cs, inter))
        print("   indices anatomiques (contrôle, le calage fait foi) : score %.1f" % sc_tel)
        for d in det_tel:
            print("      %s %s" % ("OK " if d["ok"] else "NON", d["indice"]))
        if meilleure["score"] >= sc_tel + 2 and len(det_tel) >= 2:
            print("   ATTENTION : « %s » vérifie mieux les indices (%.1f contre %.1f) : regarder controle.png" % (meilleure["variante"], meilleure["score"], sc_tel))
        elif len(det_tel) < 2:
            print("   (%d indice évalué seulement : pas de verdict, le calage fait foi ; regarder controle.png)" % len(det_tel))
        if np.linalg.det(Rc) > 0:          # l'orientation automatique (repli des stades sans calage, CS23) jugée sur ce stade calé
            auto, _ = orienter(structs)
            ecart = float(np.degrees(np.arccos(np.clip((np.trace(auto["R"] @ Rc.T) - 1) / 2, -1, 1))))
            rapport["orientation"]["controle_auto"] = {"ecart_deg": round(ecart, 1), "score": auto["score"], "z": auto["z"], "y": auto["y"]}
            print("   orientation automatique (repli sans calage) : %.1f° du calage (score %.1f, Z %s)" % (ecart, auto["score"], auto["z"]))
    else:
        um = a.um_par_unite
        if um is None and os.path.exists(a.echelles):
            try:
                um = float(json.load(open(a.echelles, encoding="utf-8")).get(cs, {}).get("um_par_unite"))
            except Exception:
                um = None
        if um is None:
            um = 1.0
            print("   échelle inconnue : 1 unité = 1 µm supposé (option --um-par-unite ou hikspoors_echelles.json)")
        k_mm = um / 1000.0
        rapport["echelle"] = {"um_par_unite": um, "source": "--um-par-unite" if a.um_par_unite else "hikspoors_echelles.json"}
        if not a.sans_calage and not a.calage:
            print("   pas de calage dans %s pour %s : orientation automatique" % (os.path.join(a.banque or "", "calage"), base or cs))
        meilleur, cand = orienter(structs)
        R = meilleur["R"]
        if "z" in a.retourner.lower():
            R = np.vstack([-R[0], R[1], -R[2]])
        if "y" in a.retourner.lower():
            R = np.vstack([-R[0], -R[1], R[2]])
        sc_final, det_final = _score(structs, R)
        def transformer(V):
            return (R @ V.T).T * k_mm
        rapport["orientation"] = {"methode": "automatique (ACP + indices anatomiques)" + (" puis --retourner " + a.retourner if a.retourner else ""),
                                  "R_lignes_XYZ_dans_repere_source": R.round(6).tolist(), "score": sc_final, "indices": det_final,
                                  "z": meilleur["z"], "y": meilleur["y"], "candidats": [{"z": c["z"], "y": c["y"], "score": c["score"]} for c in cand]}
        print("   orientation : score %.1f (Z %s, Y %s)" % (sc_final, meilleur["z"], meilleur["y"]))
        for d in det_final:
            print("      %s %s" % ("OK " if d["ok"] else "NON", d["indice"]))
        if not [d for d in det_final if d["poids"] > 0]:
            print("      AUCUN indice anatomique disponible : orientation à vérifier sur controle.png (--retourner z / y)")
    for s in structs.values():
        s["meshes_mm"] = []
        for m in s["meshes"]:
            mm = trimesh.Trimesh(transformer(np.asarray(m.vertices, float)), np.asarray(m.faces), process=False)
            if reflexion:
                mm.invert()                    # une réflexion retourne les faces : normales remises vers l'extérieur
            s["meshes_mm"].append(mm)
    allV = np.vstack([np.asarray(m.vertices) for s in structs.values() for m in s["meshes_mm"]])
    centre = (allV.max(0) + allV.min(0)) / 2
    etendue = allV.max(0) - allV.min(0)
    recentrer = not calage or repere_calage == "pipeline" or a.recentrer
    if recentrer:                              # origine = centre de la boîte (le repère source n'a pas d'origine utile)
        for s in structs.values():
            for m in s["meshes_mm"]:
                m.apply_translation(-centre)
    crl = CRL_MM.get(cs)
    rapport["dimensions"] = {"etendue_mm": etendue.round(3).tolist(), "hauteur_mm": round(float(etendue[2]), 3), "crl_typique_mm": crl,
                             "rapport_hauteur_crl": round(float(etendue[2]) / crl, 3) if crl else None,
                             "centre_boite_mm": centre.round(4).tolist(), "recentre": bool(recentrer),
                             "note": "les PDF Hikspoors sont centrés sur le cœur : l'étendue n'est pas la longueur de l'embryon"}
    print("   étendue (mm) X %.2f  Y %.2f  Z %.2f  (CRL typique %s mm) ; %s" % (*etendue, crl, "recentré sur la boîte" if recentrer else
                                                                              "coordonnées du modèle %s (centre de la boîte %s)" % (cs, centre.round(3).tolist())))

    # 4. écriture des structures
    publiable = bool(a.publiable)
    origine = "centre de la boîte des structures" if recentrer else "celle de embryons_3D/modeles/%s/ (calage)" % cs
    manifest = {"stade": cs, "source": SOURCE, "session": a.session, "unites": "mm",
                "repere": "Z crânial, Y dorsal, X = Y×Z = GAUCHE anatomique (repère direct, sans miroir), origine = " + origine,
                "specimen": "Carnegie #%s" % SPECIMENS[cs] if cs in SPECIMENS else "", "licence": LICENCE,
                "attribution": "Hikspoors JPJM, Lamers WH et al., Maastricht University ; HDBR atlas (hdbratlas.org) ; adaptation : embryo3D (mise à l'échelle, repère, nomenclature)",
                "statut": "externe", "publiable": publiable,
                "usage": "" if publiable else "référence interne : atlas d'auteurs retravaillé, jamais présenté comme notre modèle (harcelon.fr/3dht)",
                "temps": {"jours_post_fecondation": next(([x[1], x[2]] for x in AXE if x[0] == cs), None)},
                "dimensions": {"longueur_atlas_mm": None, "etendue_mm": etendue.round(3).tolist()},
                "echelle": rapport["echelle"], "orientation": {k: v for k, v in rapport["orientation"].items() if k not in ("candidats", "variantes")},
                "structures": [], "notes": "Maillages originaux des auteurs (reconstruction Amira des coupes de la collection Carnegie), seulement "
                                          "mis à l'échelle, orientés et renommés. Unions (coeur, cavites_cardiaques, arteres, veines, intestin) dérivées pour le morphing."}
    if calage:
        manifest["calage"] = {k: rapport["calage"].get(k) for k in ("fichier", "cle", "repere", "echelle_um_par_unite", "reflexion", "inter_stades", "note")}
    structs_mm = {}
    for nom in sorted(structs, key=lambda n: (structs[n]["meta"]["systeme"], n)):
        s = structs[nom]; r = s["meta"]
        m = trimesh.util.concatenate(s["meshes_mm"]) if len(s["meshes_mm"]) > 1 else s["meshes_mm"][0]
        if a.faces_max and len(m.faces) > a.faces_max:
            try:
                m = m.simplify_quadric_decimation(a.faces_max)
            except Exception:
                pass
        fn = nom + ".ply"
        m.export(os.path.join(sortie, fn))
        col_src = couleur_255(np.median(np.array(s["couleurs_source"], float), 0)) if s["couleurs_source"] else None
        col = r["couleur"] or (couleur_01(col_src) if col_src else [0.7, 0.7, 0.7])
        vol = float(m.volume) if m.is_watertight else None
        x = {"nom": nom, "nom_fr": r["nom_fr"], "nom_source": s["sources"], "systeme": r["systeme"], "fichier": fn, "couleur": col,
             "couleur_atlas": col_src, "confiance": r["confiance"], "volume_mm3": round(vol, 5) if vol is not None else None,
             "faces": int(len(m.faces)), "legende": "", "note": ("nom non reconnu par la nomenclature" if r["non_reconnu"] else "")}
        if r.get("groupe"):
            x["groupe"] = r["groupe"]
        manifest["structures"].append(x)
        structs_mm[nom] = {"meshes": [m], "couleur": col}
    # unions par groupe
    if not a.sans_unions:
        for g, gi in nomen.groupes.items():
            membres = [n for n in structs if structs[n]["meta"].get("groupe") == g]
            if len(membres) < 2 or g in structs:
                continue
            u, methode = union_groupe([m for n in membres for m in structs[n]["meshes_mm"]])
            fn = g + ".ply"; u.export(os.path.join(sortie, fn))
            manifest["structures"].append({"nom": g, "nom_fr": gi.get("nom_fr", g), "systeme": gi.get("systeme", "autre"), "fichier": fn,
                                           "couleur": gi.get("couleur", [0.7, 0.7, 0.7]), "confiance": gi.get("confiance", "moyenne"),
                                           "volume_mm3": round(float(u.volume), 5) if u.is_watertight else None, "faces": int(len(u.faces)),
                                           "derive": True, "composantes": membres, "note": "union dérivée (%s) pour le morphing entre stades" % methode})
            print("   union %-18s <- %d parties (%s, %d faces)" % (g, len(membres), methode, len(u.faces)))
    json.dump(manifest, open(os.path.join(sortie, "manifest.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    score = rapport["orientation"].get("score")
    planche_controle(structs_mm, os.path.join(sortie, "controle.png"),
                     "%s — Hikspoors 2022 — %d structures — %.4g µm/unité — %s" % (cs, len(structs_mm), um, "orientation auto, score %s" % score if score is not None
                                                                                 else "calage %s" % os.path.basename(calage)))
    rapport["verif_maillages"] = verifier_sortie(sortie, manifest)       # contrôle des maillages (lu par l'agrégateur)
    json.dump(rapport, open(os.path.join(sortie, "rapport.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False, default=str)
    print("-> %s : %d structures, manifest.json, rapport.json, controle.png%s" % (sortie, len(manifest["structures"]), "" if publiable else " (non publiable)"))


if __name__ == "__main__":
    main()
