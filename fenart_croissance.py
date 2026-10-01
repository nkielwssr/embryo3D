#!/usr/bin/env python3
"""Modèle de croissance céphalique de Fenart : lecture des coordonnées vestibulaires (docs/coordonnees_fenart.csv),
interpolation en âge, trajets, vitesses et ajustement « expansion radiale + rotation » par unité de croissance.

Bibliothèque standard seulement (pas de numpy) : utilisable tel quel sur n'importe quel poste.

    python embryo3d/fenart_croissance.py                       # contrôle du fichier + tableaux (déplacements, ajustements)
    python embryo3d/fenart_croissance.py --age 6.5             # coordonnées interpolées à 6,5 ans depuis la conception
    python embryo3d/fenart_croissance.py --derive docs/fenart_derive.json --svg docs   # exports (JSON dérivé, planches SVG)

Repère (vestibulaire, cf. docs/modele_croissance_fenart.md § 1 et § 2.4) : origine au milieu de l'axe de Perez ;
x = antéro-postérieur (négatif = avant, face), y = vertical (positif = haut), z = transversal (demi-largeur, un seul
hémicrâne, z = 0 pour les points médians). Unités : mm, âges en années depuis la conception.
"""
import argparse
import csv
import json
import math
import os
import sys

ICI = os.path.dirname(os.path.abspath(__file__))
CSV_DEFAUT = os.path.join(ICI, "docs", "coordonnees_fenart.csv")

# Stades du fichier → stade du manuscrit Captier et al. / chapitre Captier & Boë, âge (ans depuis la conception), libellé.
# La naissance est prise à 0,756 an (≈ 276 j) dans le fichier ; l'« adulte » est posé à 30 ans.
STADES = [
    ("A", "A", 0.416, "5 mois fœtal"), ("B", "B", 0.586, "7 mois fœtal"), ("C", "C", 0.756, "naissance"),
    ("D1", "D", 1.456, "8 mois 1/2"), ("D2", "E", 2.756, "2 ans"), ("E", "F", 4.756, "4 ans"),
    ("F", "G", 9.256, "8 ans 1/2"), ("G", "H", 14.756, "14 ans"), ("Ad", "Ad", 30.756, "adulte"),
]
NAISSANCE = 0.756
REGIONS = {"F": "face (massif facial, orbites)", "U": "voûte (calvaria)", "M": "mandibule", "B": "base du crâne"}
# unités de croissance testées (§ 3 du document) : F + M = crâne facial, voûte antérieure / postérieure (x < 0 / x ≥ 0 à 5 mois)
UNITES = {
    "crane_facial (F+M)": lambda p: p["region"] in ("F", "M"),
    "face (F)": lambda p: p["region"] == "F",
    "mandibule (M)": lambda p: p["region"] == "M",
    "voute (U)": lambda p: p["region"] == "U",
    "voute_anterieure (U, x<0)": lambda p: p["region"] == "U" and p["xyz"][0][0] < 0,
    "voute_posterieure (U, x>=0)": lambda p: p["region"] == "U" and p["xyz"][0][0] >= 0,
    "base (B)": lambda p: p["region"] == "B",
}


def charger(chemin=CSV_DEFAUT):
    """→ dict point → {region, pair, fiable, xyz: [(x,y,z)] dans l'ordre des STADES}. Vérifie la structure."""
    with open(chemin, encoding="utf-8") as f:
        lignes = list(csv.DictReader(f, delimiter=";"))
    codes = [s[0] for s in STADES]
    pts = {}
    for l in lignes:
        p = pts.setdefault(l["point"], {"region": l["region"], "pair": l["pair"] == "1", "fiable": l["fiable"] == "1", "stades": [], "xyz": []})
        p["stades"].append(l["stade"]); p["xyz"].append((float(l["x_mm"]), float(l["y_mm"]), float(l["z_mm"])))
        if abs(float(l["ans_depuis_conception"]) - dict((s[0], s[2]) for s in STADES)[l["stade"]]) > 1e-6:
            raise ValueError(f"{l['point']} {l['stade']} : âge {l['ans_depuis_conception']} inattendu")
    erreurs = []
    for nom, p in pts.items():
        if p["stades"] != codes: erreurs.append(f"{nom} : stades {p['stades']}")
        if not p["pair"] and any(z != 0 for _, _, z in p["xyz"]): erreurs.append(f"{nom} : point médian avec z ≠ 0")
        if p["pair"] and any(z < 0 for _, _, z in p["xyz"]): erreurs.append(f"{nom} : point pair avec z < 0")
        if p["region"] not in REGIONS: erreurs.append(f"{nom} : région {p['region']} inconnue")
    if erreurs:
        raise ValueError("fichier incohérent :\n  " + "\n  ".join(erreurs))
    return pts


def ages():
    return [s[2] for s in STADES]


def interpoler(pts, age, echelle="log", miroir=False):
    """Coordonnées de tous les points à un âge donné (ans depuis la conception), interpolation linéaire par morceaux
    en log(âge) (défaut, la croissance est plus régulière ainsi) ou en âge ; borné aux stades extrêmes.
    miroir=True : ajoute les points pairs du côté opposé (suffixe « _g », z négatif)."""
    a = ages()
    t = max(a[0], min(a[-1], age))
    f = (lambda v: math.log(v)) if echelle == "log" else (lambda v: v)
    i = max(k for k in range(len(a) - 1) if a[k] <= t) if t < a[-1] else len(a) - 2
    w = (f(t) - f(a[i])) / (f(a[i + 1]) - f(a[i]))
    out = {}
    for nom, p in pts.items():
        p0, p1 = p["xyz"][i], p["xyz"][i + 1]
        q = tuple(u + w * (v - u) for u, v in zip(p0, p1))
        out[nom] = q
        if miroir and p["pair"]:
            out[nom + "_g"] = (q[0], q[1], -q[2])
    return out


def _d(a, b):
    return math.sqrt(sum((u - v) ** 2 for u, v in zip(a, b)))


def deplacements(pts):
    """Par point : distance 5 mois → adulte, longueur du chemin, âge (depuis la conception et postnatal) où 50 % et 90 %
    du chemin sont parcourus (interpolation linéaire entre stades), direction du déplacement total dans le plan sagittal
    (angle en degrés, 0 = vers l'avant, 90 = vers le haut) et transversal (dz)."""
    a = ages()
    res = {}
    for nom, p in pts.items():
        P = p["xyz"]
        seg = [_d(P[k], P[k + 1]) for k in range(len(P) - 1)]
        chemin = sum(seg)
        cum = [0.0]
        for s in seg: cum.append(cum[-1] + s)
        def age_frac(fr):
            if chemin == 0: return None
            cible = fr * chemin
            for k in range(len(seg)):
                if cum[k + 1] >= cible:
                    w = (cible - cum[k]) / seg[k] if seg[k] else 0
                    return a[k] + w * (a[k + 1] - a[k])
            return a[-1]
        dx, dy, dz = (P[-1][i] - P[0][i] for i in range(3))
        ang = math.degrees(math.atan2(dy, -dx)) if (dx or dy) else None   # -dx : l'avant est en x négatif
        a50, a90 = age_frac(0.5), age_frac(0.9)
        res[nom] = {"region": p["region"], "pair": p["pair"], "fiable": p["fiable"],
                    "distance_mm": round(_d(P[0], P[-1]), 1), "chemin_mm": round(chemin, 1),
                    "dx_mm": round(dx, 1), "dy_mm": round(dy, 1), "dz_mm": round(dz, 1),
                    "angle_sagittal_deg": None if ang is None else round(ang, 1),
                    "age_50pct_conception": None if a50 is None else round(a50, 2), "age_90pct_conception": None if a90 is None else round(a90, 2),
                    "age_50pct_postnatal": None if a50 is None else round(a50 - NAISSANCE, 2), "age_90pct_postnatal": None if a90 is None else round(a90 - NAISSANCE, 2),
                    "vitesse_mm_par_an": [round(s / (a[k + 1] - a[k]), 1) for k, s in enumerate(seg)]}
    return res


def _similitude(A, B, translation=False):
    """Meilleure similitude plane (échelle k, rotation θ, ± translation) envoyant les points A sur B (moindres carrés,
    nombres complexes). Sans translation : rotation/homothétie autour de l'origine du repère. → (k, θ_deg, tx, ty, rms)."""
    za = [complex(x, y) for x, y in A]; zb = [complex(x, y) for x, y in B]
    if translation:
        ma, mb = sum(za) / len(za), sum(zb) / len(zb)
        za = [z - ma for z in za]; zb = [z - mb for z in zb]
    c = sum(a.conjugate() * b for a, b in zip(za, zb)) / sum(abs(a) ** 2 for a in za)
    rms = math.sqrt(sum(abs(c * a - b) ** 2 for a, b in zip(za, zb)) / len(za))
    t = (mb - c * ma) if translation else 0j
    return abs(c), math.degrees(math.atan2(c.imag, c.real)), t.real, t.imag, rms


def ajustements(pts):
    """Pour chaque unité et chaque stade : similitude plane (sagittal x,y) depuis le stade A, autour de l'origine
    vestibulaire (hypothèse « expansion radiale ± rotation ») et avec translation libre (contrôle) ; en frontal (z,y)
    pour les points pairs : facteurs d'échelle largeur (z) et hauteur (y) séparés."""
    out = {}
    for nom_u, filtre in UNITES.items():
        sel = [p for p in pts.values() if filtre(p)]
        if len(sel) < 3: continue
        A = [(p["xyz"][0][0], p["xyz"][0][1]) for p in sel]
        pairs = [p for p in sel if p["pair"]]
        lignes = []
        for i, (code, _, age, lib) in enumerate(STADES):
            B = [(p["xyz"][i][0], p["xyz"][i][1]) for p in sel]
            k, th, _, _, rms = _similitude(A, B)
            k2, th2, tx, ty, rms2 = _similitude(A, B, translation=True)
            l = {"stade": code, "age": age, "libelle": lib, "echelle": round(k, 3), "rotation_deg": round(th, 1), "rms_mm": round(rms, 1),
                 "avec_translation": {"echelle": round(k2, 3), "rotation_deg": round(th2, 1), "tx_mm": round(tx, 1), "ty_mm": round(ty, 1), "rms_mm": round(rms2, 1)}}
            if len(pairs) >= 2:
                z0 = [p["xyz"][0][2] for p in pairs]; z1 = [p["xyz"][i][2] for p in pairs]
                y0 = [p["xyz"][0][1] for p in pairs]; y1 = [p["xyz"][i][1] for p in pairs]
                kz = sum(a * b for a, b in zip(z0, z1)) / sum(a * a for a in z0)
                ky = sum(a * b for a, b in zip(y0, y1)) / sum(a * a for a in y0) if any(y0) else None
                l["frontal"] = {"echelle_largeur_z": round(kz, 3), "echelle_hauteur_y": None if ky is None else round(ky, 3)}
            lignes.append(l)
        out[nom_u] = {"n_points": len(sel), "points": sorted(pts_nom for pts_nom, p in pts.items() if filtre(p)), "stades": lignes}
    return out


def _svg_trajets(pts, plan, chemin):
    """Planche SVG des trajets 5 mois → adulte : plan 'profil' (x,y) ou 'face' (z,y, avec miroir). Couleur par région."""
    coul = {"F": "#d62828", "U": "#2f6fd6", "M": "#2a9d3f", "B": "#8c5a2b"}
    ech = 3.0
    def proj(q, signe=1):
        return (q[0], q[1]) if plan == "profil" else (signe * q[2], q[1])
    xs, ys = [], []
    for p in pts.values():
        for q in p["xyz"]:
            for s in ((1, -1) if (plan == "face" and p["pair"]) else (1,)):
                u, v = proj(q, s); xs.append(u); ys.append(v)
    x0, x1, y0, y1 = min(xs) - 10, max(xs) + 10, min(ys) - 10, max(ys) + 10
    W, H = max(620.0, (x1 - x0) * ech), (y1 - y0) * ech + 20
    def px(u, v): return (u - x0) * ech, (y1 - v) * ech + 20
    s = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.0f}" height="{H:.0f}" viewBox="0 0 {W:.0f} {H:.0f}" font-family="sans-serif" font-size="9">',
         '<rect width="100%" height="100%" fill="white"/>']
    ox, oy = px(0, 0)
    s.append(f'<line x1="0" y1="{oy:.1f}" x2="{W:.0f}" y2="{oy:.1f}" stroke="#bbb" stroke-dasharray="4 3"/>')
    s.append(f'<line x1="{ox:.1f}" y1="0" x2="{ox:.1f}" y2="{H:.0f}" stroke="#bbb" stroke-dasharray="4 3"/>')
    for nom, p in sorted(pts.items()):
        for signe in ((1, -1) if (plan == "face" and p["pair"]) else (1,)):
            if plan == "face" and not p["pair"] and signe == -1: continue
            pl = " ".join("%.1f,%.1f" % px(*proj(q, signe)) for q in p["xyz"])
            c = coul[p["region"]]; op = "1" if p["fiable"] else "0.45"
            s.append(f'<polyline points="{pl}" fill="none" stroke="{c}" stroke-width="1.2" opacity="{op}"/>')
            for q in p["xyz"]:
                cx, cy = px(*proj(q, signe)); s.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="1.6" fill="{c}" opacity="{op}"/>')
            if signe == 1:
                tx, ty = px(*proj(p["xyz"][-1], signe)); s.append(f'<text x="{tx + 3:.1f}" y="{ty - 2:.1f}" fill="{c}">{nom}</text>')
    titre = ("Profil (x antéro-postérieur, y vertical), avant à gauche" if plan == "profil" else "Face (z transversal, y vertical)") + \
            " · trajets 5 mois fœtal → adulte · repère vestibulaire (Fenart)"
    s.append(f'<text x="8" y="14" font-size="11" fill="#333">{titre}</text>')
    s.append('<text x="8" y="28" font-size="10" fill="#555">rouge = face · bleu = voûte · vert = mandibule · brun = base · pâle = point peu fiable · un point par stade (9)</text>')
    s.append(f'<line x1="{W - 70:.0f}" y1="{H - 12:.0f}" x2="{W - 70 + 20 * ech:.0f}" y2="{H - 12:.0f}" stroke="#333" stroke-width="2"/><text x="{W - 70:.0f}" y="{H - 16:.0f}" fill="#333">20 mm</text>')
    s.append("</svg>")
    with open(chemin, "w", encoding="utf-8") as f:
        f.write("\n".join(s))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--csv", default=CSV_DEFAUT)
    ap.add_argument("--age", type=float, help="âge (ans depuis la conception) : imprime les coordonnées interpolées")
    ap.add_argument("--postnatal", action="store_true", help="--age est un âge postnatal (ans après la naissance)")
    ap.add_argument("--derive", help="écrit le JSON dérivé (déplacements, ajustements) à ce chemin")
    ap.add_argument("--svg", help="dossier où écrire fenart_trajets_profil.svg et fenart_trajets_face.svg")
    ap.add_argument("--json", action="store_true", help="sortie JSON au lieu des tableaux texte")
    a = ap.parse_args(argv)
    pts = charger(a.csv)
    n_med = sum(1 for p in pts.values() if not p["pair"]); n_pair = len(pts) - n_med
    if a.age is not None:
        age = a.age + (NAISSANCE if a.postnatal else 0)
        q = interpoler(pts, age, miroir=True)
        if a.json: print(json.dumps({"age_conception": age, "points": q}, ensure_ascii=False, indent=1))
        else:
            for nom in sorted(q): print(f"{nom:7s} {q[nom][0]:7.1f} {q[nom][1]:7.1f} {q[nom][2]:7.1f}")
        return 0
    dep, aj = deplacements(pts), ajustements(pts)
    if a.derive:
        with open(a.derive, "w", encoding="utf-8") as f:
            json.dump({"_source": os.path.relpath(a.csv, ICI), "_description": "Dérivés du fichier de coordonnées de Fenart (fenart_croissance.py) : déplacements et âges à 50/90 % par point, similitudes par unité et stade.",
                       "stades": [{"code_csv": c, "code_manuscrit": m, "age_conception": ag, "libelle": l} for c, m, ag, l in STADES],
                       "n_points": len(pts), "n_medians": n_med, "n_pairs": n_pair, "deplacements": dep, "ajustements": aj}, f, ensure_ascii=False, indent=1)
        print(f"→ {a.derive}")
    if a.svg:
        for plan in ("profil", "face"):
            ch = os.path.join(a.svg, f"fenart_trajets_{plan}.svg"); _svg_trajets(pts, plan, ch); print(f"→ {ch}")
    if a.json:
        print(json.dumps({"deplacements": dep, "ajustements": aj}, ensure_ascii=False, indent=1)); return 0
    print(f"{len(pts)} points ({n_med} médians, {n_pair} pairs), {len(STADES)} stades, fichier cohérent.\n")
    print("Déplacements 5 mois fœtal → adulte (mm) et âge postnatal (ans) à 50 % / 90 % du chemin, par région :")
    for r, lib in REGIONS.items():
        print(f"\n[{r}] {lib}")
        print(f"{'point':7s} {'dist':>6s} {'chemin':>7s} {'dx':>6s} {'dy':>6s} {'dz':>6s} {'angle':>6s} {'a50':>6s} {'a90':>6s}  fiable")
        for nom in sorted(k for k, v in dep.items() if v["region"] == r):
            d = dep[nom]
            f = lambda v: "   —" if v is None else f"{v:6.1f}"
            print(f"{nom:7s} {d['distance_mm']:6.1f} {d['chemin_mm']:7.1f} {d['dx_mm']:6.1f} {d['dy_mm']:6.1f} {d['dz_mm']:6.1f} {f(d['angle_sagittal_deg'])} {f(d['age_50pct_postnatal'])} {f(d['age_90pct_postnatal'])}  {'oui' if d['fiable'] else 'non'}")
    print("\nSimilitude plane depuis 5 mois fœtal (plan sagittal, autour de l'origine vestibulaire) : échelle k, rotation θ (°, sens trigonométrique"
          " avec l'avant à gauche : θ < 0 = le haut bascule vers l'ARRIÈRE et l'occiput vers le bas, θ > 0 = le bas bascule vers l'arrière), résidu rms ;"
          " colonnes _t : même ajustement avec translation libre (|t| mm) ; frontal : échelles largeur (z) / hauteur (y) des points pairs"
          " (k_y sans sens pour les groupes proches de y = 0).")
    for u, v in aj.items():
        print(f"\n{u} — {v['n_points']} points")
        print(f"{'stade':6s} {'âge':>6s} {'k':>6s} {'θ°':>6s} {'rms':>5s}   {'k_t':>6s} {'θ_t°':>6s} {'|t|':>5s} {'rms_t':>5s}   {'k_z':>6s} {'k_y':>6s}")
        for l in v["stades"]:
            t = l["avec_translation"]; fr = l.get("frontal", {})
            kz = fr.get("echelle_largeur_z"); ky = fr.get("echelle_hauteur_y")
            print(f"{l['stade']:6s} {l['age']:6.2f} {l['echelle']:6.3f} {l['rotation_deg']:6.1f} {l['rms_mm']:5.1f}   {t['echelle']:6.3f} {t['rotation_deg']:6.1f} {math.hypot(t['tx_mm'], t['ty_mm']):5.1f} {t['rms_mm']:5.1f}   "
                  f"{'' if kz is None else f'{kz:6.3f}'} {'' if ky is None else f'{ky:6.3f}'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
