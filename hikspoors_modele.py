# -*- coding: utf-8 -*-
"""Maillages des PDF 3D Hikspoors et al. 2022 (Maastricht, hdbratlas.org/hikspoors-pdf, licence HDBR CC BY-NC-SA 4.0), déjà extraits
des PDF par la session « Carnegie Stage 13 PDF 3D » (bloc U3D 0x100), -> modèle au format de l'agrégateur :
embryons_3D/modeles/<CS>_hikspoors/{<structure>.ply, manifest.json, rapport.json, controle.png, verif_maillages.json}

    python embryo3d/hikspoors_modele.py CS13 --glb reference/hikspoors_maastricht/glb/CS13     # un GLB/PLY/OBJ par structure, ou une scène GLB
    python embryo3d/hikspoors_modele.py CS13 --npz reference/hikspoors_maastricht/u3d/CS13.npz [--scene reference/.../CS13_scene.json]
    options : --um-par-unite 1.078          échelle (µm par unité U3D ; défaut : hikspoors_echelles.json puis 1.0)
              --calage T.npy [--calage-repere modeles|pipeline]   similitude connue (13 nombres : s, R 3x3, t) au lieu de l'orientation auto
              --retourner z|y               inverse un axe après l'orientation automatique (x suit pour rester direct)
              --separer-couleurs            une partie portant plusieurs couleurs de faces (texture) est scindée par couleur
              --sans-unions                 ne pas écrire les unions coeur / cavites_cardiaques / arteres / veines / intestin
              --non-publiable               modèle gardé en local (champ publiable=false)
              --sortie DOSSIER              défaut embryons_3D/modeles/<CS>_hikspoors

Repère de sortie (convention agrégateur) : mm, Z crânial, Y dorsal, X = Y×Z = gauche anatomique (repère direct), origine = centre du corps.
Orientation automatique : ACP de l'ensemble des sommets (grand axe = crânio-caudal) ; dorsal = du cœur (ou intestin/foie) vers le tube neural
(ou somites/ganglions/notochorde/aorte dorsale) ; le signe crânial est choisi par vote (cœur au-dessus du foie et de l'intestin, arcs aortiques
au-dessus du cœur, tube neural plus large du côté crânial, structures gauches à +X…). Tout est consigné dans rapport.json et controle.png :
VÉRIFIER la planche avant publication, corriger avec --retourner si besoin.
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
                parties.append({"nom": nom, "mesh": m, "couleur": col, "couleurs_faces": fc, "fichier": f})
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

    def trouve(cat):
        for k in _CLES[cat]:
            if k in lc:
                return lc[k]
        return None
    infos = {}
    if scene_json and os.path.exists(scene_json):
        sc = json.load(open(scene_json, encoding="utf-8"))
        items = sc.get("noeuds") or sc.get("nodes") or sc.get("parties") or sc.get("parts") or sc.get("objets") or sc.get("structures") or sc
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
    for p in parties:                          # compléments du _scene.json
        inf = infos.get(normaliser(p["nom"]))
        if inf:
            if p["couleur"] is None:
                c = inf.get("couleur") or inf.get("color") or inf.get("diffuse") or inf.get("rgb")
                p["couleur"] = couleur_255(c) if c is not None else None
            p["visible"] = inf.get("visible", inf.get("V", True))
            if inf.get("nom_affiche") or inf.get("label") or inf.get("libelle"):
                p["nom"] = inf.get("nom_affiche") or inf.get("label") or inf.get("libelle")
    return parties


# ----------------------------------------------------------------------------------------------------------------------- nomenclature
class Nomenclature:
    def __init__(self, chemin=NOMENCLATURE):
        self.t = json.load(open(chemin, encoding="utf-8"))
        self.ignorer = [re.compile(m, re.I) for m in self.t.get("ignorer", [])]
        self.cote_g = re.compile(self.t["cotes"]["gauche"], re.I)
        self.cote_d = re.compile(self.t["cotes"]["droit"], re.I)
        self.regles = []
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
            self.regles.append({"motif": re.compile(r["motif"], re.I), "exige": None, "nom": r["nom"], "nom_fr": r.get("nom_fr", r["nom"]),
                                "systeme": r.get("systeme", "autre"), "groupe": r.get("groupe"), "couleur": r.get("couleur"),
                                "confiance": r.get("confiance", "moyenne"), "capture": re.compile(r["capture"], re.I) if r.get("capture") else None,
                                "regle": r["motif"]})
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
                    "source_partie": p["nom"], "couleur_scindee": col, "fraction_faces": round(float(sel.mean()), 4)})
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


CRANIAL = [   # (A crânial à B) : motifs sur les noms canoniques, poids
    (["^myocarde_", "^cavite_", "^coeur$", "^tube_cardiaque", "^myocarde_coeur"], ["^foie$", "^septum_transversum", "^intestin_moyen", "^intestin_posterieur", "^veine_ombilicale", "^artere_ombilicale", "^cordon"], 2.0),
    (["^arc_aortique", "^arcs_aortiques", "^sac_aortique", "^carotide", "^pharynx", "^encephale", "^yeux", "^vesicules_otiques"], ["^myocarde_", "^cavite_", "^coeur$", "^foie$", "^intestin"], 2.0),
    (["^myocarde_voie_efferente", "^cavite_voie_efferente", "^myocarde_ventricule", "^cavite_ventricule"], ["^myocarde_sinus_veineux", "^cavite_sinus_veineux", "^veine_vitelline", "^canal_hepatocardiaque", "^veine_cave_inferieure"], 1.0),
    (["^intestin_anterieur", "^oesophage", "^estomac", "^poumons"], ["^intestin_posterieur", "^veine_ombilicale", "^artere_ombilicale", "^cordon"], 1.0),
]
DORSAUX = ["^tube_neural", "^somites", "^ganglions_spinaux", "^notochorde", "^aorte_dorsale", "^encephale", "^veine_cardinale"]
VENTRAUX = ["^myocarde_", "^cavite_", "^coeur$", "^tube_cardiaque", "^foie$", "^intestin", "^vesicule_vitelline", "^estomac", "^pericarde", "^epicarde", "^septum_transversum"]


def _score(structs, R):
    """score d'un repère R (lignes = X gauche, Y dorsal, Z crânial dans le repère source) : somme des indices vérifiés, détail par indice"""
    def c(motifs):
        v = _centre(structs, motifs)
        return None if v is None else R @ v
    detail = []; total = 0.0
    for a, b, w in CRANIAL:
        ca, cb = c(a), c(b)
        if ca is None or cb is None:
            continue
        ok = ca[2] > cb[2]
        total += w if ok else -w
        detail.append({"indice": "crânial : %s > %s" % (a[0], b[0]), "ok": bool(ok), "ecart_mm": None, "poids": w})
    cd, cv = c(DORSAUX), c(VENTRAUX)
    if cd is not None and cv is not None:
        ok = cd[1] > cv[1]
        total += 2.0 if ok else -2.0
        detail.append({"indice": "dorsal : tube neural/somites au-dessus (Y) du cœur/intestin", "ok": bool(ok), "poids": 2.0})
    # tube neural plus large du côté crânial (vésicules cérébrales)
    if "tube_neural" in structs:
        P = R @ np.vstack([_echantillon(m, 20000) for m in structs["tube_neural"]["meshes"]]).T
        z = P[2]; lo, hi = np.percentile(z, [15, 85])
        def rayon(sel):
            Q = P[:2, sel]
            return float(np.sqrt(((Q - Q.mean(1, keepdims=True)) ** 2).sum(0)).mean()) if sel.sum() > 20 else None
        r_haut, r_bas = rayon(z >= hi), rayon(z <= lo)
        if r_haut and r_bas and abs(r_haut - r_bas) > 0.15 * max(r_haut, r_bas):
            ok = r_haut > r_bas
            total += 1.0 if ok else -1.0
            detail.append({"indice": "tube neural plus large en haut (encéphale)", "ok": bool(ok), "poids": 1.0, "rayon_haut": r_haut, "rayon_bas": r_bas})
    # gauche/droite : structures appariées
    paires = 0; dxs = []
    for nom in structs:
        if nom.endswith("_gauche") and nom[:-7] + "_droit" in structs:
            g = R @ np.vstack([_echantillon(m) for m in structs[nom]["meshes"]]).mean(0)
            d = R @ np.vstack([_echantillon(m) for m in structs[nom[:-7] + "_droit"]["meshes"]]).mean(0)
            dxs.append(float(g[0] - d[0])); paires += 1
    if paires:
        ok = np.median(dxs) > 0
        total += 2.0 if ok else -2.0
        detail.append({"indice": "gauche à +X (%d paires)" % paires, "ok": bool(ok), "poids": 2.0, "dx_median": float(np.median(dxs))})
    return total, detail


def orienter(structs):
    """Renvoie R (3x3, lignes X/Y/Z dans le repère source), le score et le détail. Candidats : Z = ±e1 (ACP), Y = indice dorsal ⊥ Z (sinon ±e2)."""
    P = np.vstack([_echantillon(m) for s in structs.values() for m in s["meshes"]])
    c0 = P.mean(0)
    _, _, Vt = np.linalg.svd(P - c0, full_matrices=False)
    e1, e2 = Vt[0], Vt[1]
    cd, cv = _centre(structs, DORSAUX), _centre(structs, VENTRAUX)
    cand = []
    for sz in (1, -1):
        Z = sz * e1
        ys = []
        if cd is not None and cv is not None:
            d = cd - cv; d = d - (d @ Z) * Z
            if np.linalg.norm(d) > 1e-9:
                ys.append(("indice dorsal", d / np.linalg.norm(d)))
        if not ys:
            for sy in (1, -1):
                y = sy * e2 - ((sy * e2) @ Z) * Z
                ys.append(("ACP e2 (%+d)" % sy, y / np.linalg.norm(y)))
        for lab, Y in ys:
            X = np.cross(Y, Z)
            R = np.vstack([X, Y, Z])
            sc, det = _score(structs, R)
            cand.append({"R": R, "score": sc, "detail": det, "z": "ACP e1 (%+d)" % sz, "y": lab})
    cand.sort(key=lambda c: -c["score"])
    return cand[0], cand


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


# ----------------------------------------------------------------------------------------------------------------------- principal
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stade")
    ap.add_argument("--glb", help="dossier (ou fichier) de maillages GLB/PLY/OBJ/STL")
    ap.add_argument("--npz"); ap.add_argument("--scene")
    ap.add_argument("--um-par-unite", type=float, default=None)
    ap.add_argument("--echelles", default=ECHELLES)
    ap.add_argument("--calage"); ap.add_argument("--calage-repere", default="modeles", choices=["modeles", "pipeline"])
    ap.add_argument("--retourner", default="", help="axes à inverser après l'orientation auto : z, y ou zy")
    ap.add_argument("--separer-couleurs", action="store_true")
    ap.add_argument("--sans-unions", action="store_true")
    ap.add_argument("--non-publiable", action="store_true")
    ap.add_argument("--nomenclature", default=NOMENCLATURE)
    ap.add_argument("--sortie"); ap.add_argument("--session", default="Hikspoors → modèles (cloud, 26/09)")
    ap.add_argument("--faces-max", type=int, default=0, help="décimation des structures au-delà de N faces (0 = aucune)")
    a = ap.parse_args()
    for k in ("glb", "npz", "scene", "calage", "sortie"):
        setattr(a, k, abs0(getattr(a, k)))
    if not os.path.exists(a.echelles):
        a.echelles = abs0(a.echelles)
    if not os.path.exists(a.nomenclature):
        a.nomenclature = abs0(a.nomenclature)
    cs = stade_norme(a.stade)
    sortie = a.sortie or os.path.join("embryons_3D", "modeles", cs + "_hikspoors")
    os.makedirs(sortie, exist_ok=True)
    nomen = Nomenclature(a.nomenclature)

    # 1. lecture
    if a.npz:
        parties = charger_npz(a.npz, a.scene)
    elif a.glb:
        parties = charger_fichiers(a.glb, cs)
    else:
        ap.error("--glb ou --npz requis")
    print("%s : %d parties lues" % (cs, len(parties)))
    if a.separer_couleurs:
        parties = [q for p in parties for q in scinder_par_couleur(p)]
        print("   après scission par couleur : %d parties" % len(parties))

    # 2. nomenclature -> structures canoniques (fusion des parties de même nom)
    structs, ignores, non_reconnus, journal = {}, [], [], []
    for p in parties:
        r = nomen.par_couleur_regle(p.get("source_partie", p["nom"]), p.get("couleur_scindee")) if p.get("couleur_scindee") is not None else None
        if r is None:
            r = nomen.chercher(p.get("source_partie", p["nom"]))
        if r.get("ignorer"):
            ignores.append(p["nom"]); continue
        if r["non_reconnu"]:
            non_reconnus.append(p["nom"])
        s = structs.setdefault(r["nom"], {"meshes": [], "sources": [], "meta": r, "couleurs_source": []})
        s["meshes"].append(p["mesh"]); s["sources"].append(p["nom"])
        if p.get("couleur") is not None:
            s["couleurs_source"].append(p["couleur"])
        journal.append({"source": p["nom"], "nom": r["nom"], "regle": r["regle"], "faces": int(len(p["mesh"].faces)), "couleur_source": p.get("couleur")})
    if not structs:
        raise SystemExit("aucune structure retenue")
    print("   %d structures canoniques ; non reconnues : %s ; ignorées : %s" % (len(structs), non_reconnus or "aucune", ignores or "aucune"))

    # 3. échelle
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

    # 4. repère
    rapport = {"stade": cs, "source": SOURCE, "echelle": {"um_par_unite": um}, "non_reconnus": non_reconnus, "ignores": ignores, "journal": journal}
    if a.calage:
        T = np.load(a.calage).ravel()
        s_, R_, t_ = float(T[0]), T[1:10].reshape(3, 3), T[10:13]
        def transformer(V):
            W = (s_ * (R_ @ V.T)).T + t_
            if a.calage_repere == "pipeline":
                W[:, 0] = -W[:, 0]
            return W
        rapport["orientation"] = {"methode": "calage fourni (%s, repère %s)" % (a.calage, a.calage_repere)}
        print("   calage fourni : échelle %.4g, repère %s" % (s_, a.calage_repere))
    else:
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
                                  "candidats": [{"z": c["z"], "y": c["y"], "score": c["score"]} for c in cand]}
        print("   orientation : score %.1f (%s, %s)" % (sc_final, meilleur["z"], meilleur["y"]))
        for d in det_final:
            print("      %s %s" % ("OK " if d["ok"] else "NON", d["indice"]))
        if not det_final:
            print("      AUCUN indice anatomique disponible : orientation à vérifier sur controle.png (--retourner z / y)")
    # centrage sur le centre du corps (boîte englobante de l'ensemble)
    for s in structs.values():
        s["meshes_mm"] = []
        for m in s["meshes"]:
            s["meshes_mm"].append(trimesh.Trimesh(transformer(np.asarray(m.vertices, float)), np.asarray(m.faces), process=False))
    if not a.calage and np.linalg.det(R) < 0:      # ne devrait pas arriver (X = Y×Z), sécurité
        for s in structs.values():
            for m in s["meshes_mm"]:
                m.invert()
    allV = np.vstack([np.asarray(m.vertices) for s in structs.values() for m in s["meshes_mm"]])
    centre = (allV.max(0) + allV.min(0)) / 2
    etendue = allV.max(0) - allV.min(0)
    for s in structs.values():
        for m in s["meshes_mm"]:
            m.apply_translation(-centre)
    crl = CRL_MM.get(cs)
    rapport["dimensions"] = {"etendue_mm": etendue.round(3).tolist(), "hauteur_mm": round(float(etendue[2]), 3), "crl_typique_mm": crl,
                             "rapport_hauteur_crl": round(float(etendue[2]) / crl, 3) if crl else None,
                             "note": "les PDF Hikspoors sont centrés sur le cœur : l'étendue n'est pas la longueur de l'embryon"}
    print("   étendue (mm) X %.2f  Y %.2f  Z %.2f  (CRL typique %s mm)" % (*etendue, crl))

    # 5. écriture des structures
    manifest = {"stade": cs, "source": SOURCE, "session": a.session, "unites": "mm",
                "repere": "Z crânial, Y dorsal, X = Y×Z = GAUCHE anatomique (repère direct, sans miroir), origine = centre du corps",
                "specimen": "Carnegie #%s" % SPECIMENS[cs] if cs in SPECIMENS else "", "licence": LICENCE,
                "attribution": "Hikspoors JPJM, Lamers WH et al., Maastricht University ; HDBR atlas (hdbratlas.org) ; adaptation : embryo3D (mise à l'échelle, repère, nomenclature)",
                "statut": "externe", "publiable": not a.non_publiable, "usage": "" if not a.non_publiable else "référence interne",
                "temps": {"jours_post_fecondation": next(([x[1], x[2]] for x in AXE if x[0] == cs), None)},
                "dimensions": {"longueur_atlas_mm": None, "etendue_mm": etendue.round(3).tolist()},
                "echelle": {"um_par_unite": um}, "orientation": {k: v for k, v in rapport["orientation"].items() if k != "candidats"},
                "structures": [], "notes": "Maillages originaux des auteurs (reconstruction Amira des coupes de la collection Carnegie), seulement "
                                          "mis à l'échelle, orientés et renommés. Unions (coeur, cavites_cardiaques, arteres, veines, intestin) dérivées pour le morphing."}
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
    json.dump(rapport, open(os.path.join(sortie, "rapport.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False, default=str)
    planche_controle(structs_mm, os.path.join(sortie, "controle.png"),
                     "%s — Hikspoors 2022 — %d structures — échelle %.4g µm/unité — score orientation %s" % (cs, len(structs_mm), um, rapport["orientation"].get("score", "calage")))
    # contrôle des maillages (verif_maillages.json lu par l'agrégateur)
    try:
        import verif_maillages
        rel = os.path.relpath(sortie, os.path.join("embryons_3D", "modeles"))
        if not rel.startswith(".."):
            verif_maillages.main(rel)
    except Exception as e:
        print("   verif_maillages :", e)
    print("-> %s : %d structures, manifest.json, rapport.json, controle.png" % (sortie, len(manifest["structures"])))


if __name__ == "__main__":
    main()
