# -*- coding: utf-8 -*-
"""Amnios et vésicule ombilicale construits à partir d'une figure au trait (reconstruction médiane + vue dorsale alignées),
par exemple ehd.org, Carnegie stage 7, figure 1 (No. 7802) enregistrée par l'utilisateur dans CS7.png.

    python embryo3d/annexes_figure.py CS7.png CS07

Lecture de la figure (tête à droite, dos en haut) :
  - reconstruction médiane : A.C. = région blanche fermée par la paroi noire (amnios), U.V. = région blanche fermée par la ligne
    hachurée (vésicule) ; les hachures sont scellées par une fermeture morphologique ;
  - vue dorsale au-dessus, alignée : contour du disque et de l'amnios → demi-largeur latérale de l'amnios le long de l'axe ;
  - axe et échelle : extrémités du disque dans la coupe médiane (--queue X,Y --tete X,Y), calées sur la longueur de l'épiblaste
    de notre modèle (les proportions de la figure sont conservées, pas sa barre d'échelle).
Seuls des contours (mesures) sont tirés de la figure ; la géométrie est générée par notre code (loft de tranches elliptiques).
Sorties : embryons_3D/modeles/<dossier>_complements/{amnios_trace.ply, vesicule_trace.ply, complements.json, controle_trace.png}.
"""
import json
import os
import sys

import numpy as np
import trimesh
from PIL import Image, ImageDraw
from scipy import ndimage

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(RACINE)

# réglages de la figure CS7 (pixels de CS7.png, 1022×970)
ZONE_MEDIANE = (315, 790, 95, 510)          # y0, y1, x0, x1
ZONE_DORSALE = (92, 330, 100, 488)   # sans la règle du haut ni les pointillés latéraux
GRAINE_AC, GRAINE_UV = (400, 490), (330, 640)   # (x, y) dans les régions « A.C. » et « U.V. »
QUEUE, TETE = (185, 585), (378, 555)        # extrémités du disque dans la coupe médiane
ECHELLE_FIGURE = (665, 735)                 # barre verticale de 0,1 mm (y0, y1) : contrôle seulement


def composante(m, graine):
    lab, _ = ndimage.label(m)
    return lab == lab[graine[1], graine[0]]


def main(figure, dossier, opts=None):
    opts = opts or {}
    img = Image.open(figure).convert("RGB")
    gris = np.array(img.convert("L")).astype(float)
    noir = gris < 110

    y0, y1, x0, x1 = ZONE_MEDIANE
    zone = np.zeros_like(noir); zone[y0:y1, x0:x1] = True
    ferme = ndimage.binary_closing(noir & zone, iterations=3)
    blanc = ~ferme & zone
    ac = ndimage.binary_fill_holes(composante(blanc, GRAINE_AC))          # lettres « A.C. » comblées
    uv = ndimage.binary_fill_holes(composante(blanc, GRAINE_UV))
    amn = ndimage.binary_dilation(ac, iterations=6)                       # + épaisseur de la paroi
    vit = ndimage.binary_dilation(uv, iterations=3)
    amn &= ~ndimage.binary_erosion(vit, iterations=4)                     # pas de chevauchement au niveau du disque

    y0, y1, x0, x1 = ZONE_DORSALE
    zd = np.zeros_like(noir); zd[y0:y1, x0:x1] = True
    trait = (gris < 200) & zd                                             # contour fin (anticrénelé), avec de petites lacunes
    dorsal = np.zeros_like(noir)                                          # colonne par colonne : du trait le plus haut au plus bas
    for x in range(x0, x1):
        ys_ = np.nonzero(trait[:, x])[0]
        if len(ys_) > 1:
            dorsal[ys_.min():ys_.max() + 1, x] = True

    tete, queue = np.array(opts.get("tete", TETE), float), np.array(opts.get("queue", QUEUE), float)
    ax = (tete - queue) / np.linalg.norm(tete - queue)
    L_px = float(np.linalg.norm(tete - queue)); c_img = (tete + queue) / 2
    perp = np.array([ax[1], -ax[0]])
    if perp[1] > 0:
        perp = -perp                                                      # vers le haut de l'image = dorsal

    d = os.path.join("embryons_3D", "modeles", dossier)
    man = json.load(open(os.path.join(d, "manifest.json"), encoding="utf-8"))
    corps = [x for x in man["structures"] if x["nom"] in ("epiblast", "epiblast_or_ectoderm")]
    V = np.vstack([trimesh.load(os.path.join(d, x["fichier"]), force="mesh").vertices for x in corps])
    # repère du disque : normale = petit axe (ACP) orientée vers la cavité amniotique, crânial = Z projeté dans le plan du disque
    c0 = V.mean(0); _, vec = np.linalg.eigh(np.cov((V - c0).T)); n = vec[:, 0]
    A = [x for x in man["structures"] if x["nom"] in ("amniotic_cavity", "amnion_wall")]
    if A:
        cA = np.vstack([trimesh.load(os.path.join(d, x["fichier"]), force="mesh").vertices for x in A]).mean(0)
        n = n if (cA - c0) @ n > 0 else -n
    elif n[1] < 0:
        n = -n
    tz = np.array([0, 0, 1.0]) - n[2] * n; tz /= np.linalg.norm(tz)
    lat = np.cross(n, tz)                                                  # base (lat, n, tz) directe comme (X, Y, Z)
    pt, pn, pl = (V - c0) @ tz, (V - c0) @ n, (V - c0) @ lat
    L_mod = float(np.percentile(pt, 99) - np.percentile(pt, 1))
    c_mod = c0 + tz * (np.percentile(pt, 99) + np.percentile(pt, 1)) / 2 + lat * (np.percentile(pl, 99) + np.percentile(pl, 1)) / 2         + n * (np.percentile(pn, 99) + np.percentile(pn, 1)) / 2
    s = L_mod / L_px

    # demi-largeur latérale (vue dorsale alignée : même abscisse que la coupe médiane)
    dys, dxs = np.nonzero(dorsal)
    def demi_largeur_px(x):
        sel = np.abs(dxs - x) <= 1
        return 0.0 if not sel.any() else (dys[sel].max() - dys[sel].min()) / 2

    def loft(masque, largeur, n_z=90, n_t=56):
        yy, xx = np.nonzero(masque); Q = np.c_[xx, yy].astype(float)
        u = (Q - c_img) @ ax; vv = (Q - c_img) @ perp
        u0, u1 = u.min(), u.max(); anneaux = []
        for k in range(n_z):
            a = u0 + (u1 - u0) * (k + 0.5) / n_z
            sel = np.abs(u - a) < (u1 - u0) / n_z
            if sel.sum() < 3:
                continue
            haut, bas = vv[sel].max(), vv[sel].min()
            cy = (haut + bas) / 2 * s; sy = max((haut - bas) / 2 * s, 1e-3)
            x_img = (c_img + a * ax)[0]
            sx = largeur(x_img, sy)
            t = np.linspace(0, 2 * np.pi, n_t, endpoint=False)
            anneaux.append(c_mod + np.outer(sx * np.cos(t), lat) + np.outer(cy + sy * np.sin(t), n) + a * s * tz)
        R_ = np.array(anneaux); nz = len(R_)
        verts = list(R_.reshape(-1, 3)); faces = []
        for i in range(nz - 1):
            for j in range(n_t):
                a0, a1 = i * n_t + j, i * n_t + (j + 1) % n_t
                b0, b1 = a0 + n_t, a1 + n_t
                faces += [[a0, b0, a1], [a1, b0, b1]]
        for bout, i in ((0, 0), (1, nz - 1)):
            verts.append(R_[i].mean(0)); cidx = len(verts) - 1
            for j in range(n_t):
                a0, a1 = i * n_t + j, i * n_t + (j + 1) % n_t
                faces.append([cidx, a1, a0] if bout == 0 else [cidx, a0, a1])
        m = trimesh.Trimesh(np.array(verts), np.array(faces), process=True)
        m.fix_normals()
        return trimesh.smoothing.filter_taubin(m, iterations=10) or m

    larg_amn = lambda x, sy: max(demi_largeur_px(x) * s, sy * 0.5)
    m_amn = loft(amn, larg_amn)
    m_vit = loft(vit, lambda x, sy: sy)                                   # section ronde

    out = os.path.join("embryons_3D", "modeles", dossier + "_complements")
    os.makedirs(out, exist_ok=True)
    m_amn.export(os.path.join(out, "amnios_trace.ply"))
    m_vit.export(os.path.join(out, "vesicule_trace.ply"))
    s_fig = 0.1 / (ECHELLE_FIGURE[1] - ECHELLE_FIGURE[0])
    src = f"figure au trait {os.path.basename(figure)} (ehd.org, Carnegie No. 7802, reconstruction médiane + vue dorsale), contours seulement"
    info = {"repere_disque": {"normale": [round(float(x), 4) for x in n], "cranial": [round(float(x), 4) for x in tz],
                              "inclinaison_deg": round(float(np.degrees(np.arccos(abs(n[1])))), 1)},
            "stade": man.get("stade"), "source_trace": os.path.basename(figure), "echelle_mm_par_px": round(float(s), 6),
            "echelle_figure_mm_par_px": round(s_fig, 6), "rapport_modele_sur_figure": round(float(s / s_fig), 3),
            "longueur_disque_mm": round(float(L_mod), 3), "longueur_disque_figure_mm": round(L_px * s_fig, 3),
            "amnios": {"fichier": "amnios_trace.ply", "etendue_mm": [round(float(v), 3) for v in m_amn.extents],
                       "source": src + " ; largeur latérale = vue dorsale"},
            "vitellus": {"fichier": "vesicule_trace.ply", "etendue_mm": [round(float(v), 3) for v in m_vit.extents],
                         "source": src + " ; section ronde"}}
    json.dump(info, open(os.path.join(out, "complements.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)

    # contrôle : régions lues + maillages reprojetés sur la figure
    ctrl = img.copy(); dr = ImageDraw.Draw(ctrl, "RGBA")
    for mk, col in ((amn, (80, 170, 255, 70)), (vit, (255, 200, 90, 70)), (dorsal, (120, 255, 120, 60))):
        yy, xx = np.nonzero(mk)
        for x_, y_ in list(zip(xx, yy))[::3]:
            dr.point((int(x_), int(y_)), fill=col)
    vers_img = lambda Pm: c_img + np.outer((Pm - c_mod) @ tz / s, ax) + np.outer((Pm - c_mod) @ n / s, perp)
    for m, col in ((m_amn, (30, 90, 255)), (m_vit, (230, 140, 0))):
        for p in vers_img(m.vertices[::4]):
            dr.ellipse((p[0] - 1, p[1] - 1, p[0] + 1, p[1] + 1), fill=col)
    dr.line([tuple(queue), tuple(tete)], fill=(255, 0, 255), width=3)
    ctrl.save(os.path.join(out, "controle_trace.png"))
    print(json.dumps(info, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    a = sys.argv[3:]; o = {}
    for i, k in enumerate(a):
        if k in ("--tete", "--queue"):
            o[k[2:]] = [float(v) for v in a[i + 1].split(",")]
    main(sys.argv[1], sys.argv[2], o)
