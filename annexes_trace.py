# -*- coding: utf-8 -*-
"""Amnios et vésicule vitelline construits à partir des tracés de l'utilisateur (PSD), pour un stade.

    python embryo3d/annexes_trace.py CS9.psd CS09

PSD attendu (calque de tracé au-dessus de la photo) :
  - à gauche : vue DORSALE, contour vert de l'embryon (tête en haut) → largeur de l'amnios le long de l'axe ;
  - à droite : vue de PROFIL, tête à droite, dos en haut : contour vert = amnios, contour rouge = vésicule vitelline.
Échelle et position : l'embryon clair de la photo de profil (seuillé à l'intérieur de l'amnios) est calé sur la longueur crânio-caudale
et le centre de notre modèle (peau / ectoderme du stade). Repère du modèle : Z crânial, Y dorsal, X gauche.
Chaque enveloppe est un « loft » : tranches le long de Z, section elliptique (amnios : demi-largeur = demi-largeur dorsale × 1,1 ;
vésicule : section ronde). Sorties : embryons_3D/modeles/<dossier>_complements/{amnios_trace.ply, vesicule_trace.ply,
complements.json, controle_trace.png} (géométrie générée par notre code à partir des tracés).
"""
import json
import os
import sys

import numpy as np
import trimesh
from psd_tools import PSDImage
from scipy import ndimage

RACINE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(RACINE)
MARGE_LARGEUR = 1.10


def masques(psd_path):
    """(trace, photo) : masque des pixels de tracé (couleur saturée, dernier calque pixel) et photo en niveaux de gris sans le tracé."""
    psd = PSDImage.open(psd_path)
    calques = list(psd.descendants())
    trace = next(l for l in calques[::-1] if l.kind == "pixel")
    ox, oy = trace.bbox[:2]
    im = np.array(trace.composite().convert("RGBA")).astype(int)
    rgb, A = im[..., :3], im[..., 3]
    m = ((rgb.max(-1) - rgb.min(-1)) > 90) & (A > 100)
    # teinte (0-360) pour séparer les tracés de couleurs différentes (contours qui se touchent)
    r_, g_, b_ = [rgb[..., i].astype(float) for i in range(3)]
    mx, mn = rgb.max(-1).astype(float), rgb.min(-1).astype(float); d_ = np.maximum(mx - mn, 1)
    h = np.where(mx == r_, ((g_ - b_) / d_) % 6, np.where(mx == g_, (b_ - r_) / d_ + 2, (r_ - g_) / d_ + 4)) * 60
    trace.visible = False
    photo = psd.composite(force=True).convert("RGB")
    fond = np.array(photo.convert("L")).astype(float)
    Hh, Ww = fond.shape
    def plein(mm):
        o = np.zeros((Hh, Ww), bool); o[oy:oy + mm.shape[0], ox:ox + mm.shape[1]] = mm[:Hh - oy, :Ww - ox]; return o
    couleurs = []
    for c0 in range(0, 360, 30):                                   # une classe par tranche de teinte bien remplie
        mm = m & (((h - c0) % 360) < 30)
        if mm.sum() > 1500:
            couleurs.append(plein(mm))
    return plein(m), couleurs, fond, photo


def remplir(contour):
    return ndimage.binary_fill_holes(ndimage.binary_closing(contour, iterations=4))


def composantes(m):
    lab, n = ndimage.label(m)
    return sorted([(lab == i) for i in range(1, n + 1)], key=lambda c: -c.sum())


def main(psd_path, dossier, opts=None):
    opts = opts or {}
    trace, couleurs, fond, photo = masques(psd_path)
    comps = [c for mc in couleurs for c in composantes(ndimage.binary_closing(mc, iterations=3)) if c.sum() > 1500]
    def forme(c):
        ys_, xs_ = np.nonzero(c); return (ys_.max() - ys_.min()) / max(1, xs_.max() - xs_.min())
    if opts.get("profil_seul"):
        debout, couches = [], comps
        dorsal = None
    else:
        debout = [c for c in comps if forme(c) > 1.4]; couches = [c for c in comps if forme(c) <= 1.4]
        dorsal = remplir(max(debout, key=lambda c: c.sum()))
    pleins = [remplir(c) for c in couches]
    pleins.sort(key=lambda f: np.nonzero(f)[0].mean())         # amnios : le contour de profil le plus haut (dos en haut), vésicule en dessous
    amn = pleins[0]
    vit = pleins[1] if len(pleins) > 1 else None
    if opts.get("role") == "vitellus" and vit is None:
        amn, vit = None, pleins[0]
    vert = trace

    # --- embryon clair dans l'amnios (vue de profil) → axe, longueur, centre
    DOS = {"haut": (0, -1), "bas": (0, 1), "gauche": (-1, 0), "droite": (1, 0)}
    if opts.get("tete") and opts.get("queue"):                # extrémités données par l'utilisateur (photos courbées ou retournées)
        tete, queue = np.array(opts["tete"], float), np.array(opts["queue"], float)
        ax = (tete - queue) / np.linalg.norm(tete - queue); L_px = float(np.linalg.norm(tete - queue)); c_img = (tete + queue) / 2
        dd = np.array(DOS[opts.get("dos", "haut")], float); perp = dd - (dd @ ax) * ax; perp /= np.linalg.norm(perp)
        P = np.array([tete, queue]); proj = (P - c_img) @ ax
    else:
        seuil = np.percentile(fond[amn], 70)
        cl = composantes(ndimage.binary_opening(amn & (fond > seuil) & ~ndimage.binary_dilation(trace, iterations=6), iterations=2))
        emb = np.any([c for c in cl if c.sum() > 0.08 * cl[0].sum()], axis=0)      # tête et queue peuvent être séparées au cou
        ys, xs = np.nonzero(emb); P = np.c_[xs, ys].astype(float); c_img = P.mean(0)
        w, v = np.linalg.eigh(np.cov((P - c_img).T)); ax = v[:, 1]
        if ax[0] < 0:
            ax = -ax                                               # tête à droite
        perp = np.array([ax[1], -ax[0]])                           # vers le haut de l'image = dorsal
        if perp[1] > 0:
            perp = -perp
        proj = (P - c_img) @ ax
        L_px = np.percentile(proj, 99.5) - np.percentile(proj, 0.5)

    # --- modèle : longueur et centre du corps
    d = os.path.join("embryons_3D", "modeles", dossier)
    man = json.load(open(os.path.join(d, "manifest.json"), encoding="utf-8"))
    corps = [x for x in man["structures"] if any(k in x["nom"] for k in ("skin", "epidermal_ectoderm", "ectoderme_epidermique", "epiblast"))]
    V = np.vstack([trimesh.load(os.path.join(d, x["fichier"]), force="mesh").vertices for x in corps])
    lo, hi = np.percentile(V, 1, 0), np.percentile(V, 99, 0)
    L_mod = hi[2] - lo[2]; c_mod = (lo + hi) / 2
    s = L_mod / L_px                                           # mm par pixel
    img2mod = lambda q: np.c_[c_mod[0] + 0 * q[:, 0], c_mod[1] + ((q - c_img) @ perp) * s, c_mod[2] + ((q - c_img) @ ax) * s]

    # --- largeur dorsale le long de l'axe (tête en haut dans la vue dorsale), normalisée sur la longueur
    dys, dxs = np.nonzero(dorsal) if dorsal is not None else (np.array([0, 1]), np.array([0, 0]))
    y0, y1 = dys.min(), dys.max()
    def demi_largeur(t):                                       # t ∈ [0,1] : 0 = queue, 1 = tête
        if dorsal is None:
            return 0.0
        row = int(round(y1 - t * (y1 - y0)))
        cols = dxs[dys == row]
        return 0.0 if len(cols) == 0 else (cols.max() - cols.min()) / 2 / (y1 - y0) * L_mod

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
            z = c_mod[2] + a * s
            sx = largeur(z, sy)
            t = np.linspace(0, 2 * np.pi, n_t, endpoint=False)
            anneaux.append(np.c_[c_mod[0] + sx * np.cos(t), c_mod[1] + cy + sy * np.sin(t), np.full(n_t, z)])
        R_ = np.array(anneaux); nz = len(R_)
        verts = list(R_.reshape(-1, 3)); faces = []
        for i in range(nz - 1):
            for j in range(n_t):
                a0, a1 = i * n_t + j, i * n_t + (j + 1) % n_t
                b0, b1 = a0 + n_t, a1 + n_t
                faces += [[a0, b0, a1], [a1, b0, b1]]
        for bout, i in ((0, 0), (1, nz - 1)):                  # bouchons
            verts.append(R_[i].mean(0)); cidx = len(verts) - 1
            for j in range(n_t):
                a0, a1 = i * n_t + j, i * n_t + (j + 1) % n_t
                faces.append([cidx, a1, a0] if bout == 0 else [cidx, a0, a1])
        m = trimesh.Trimesh(np.array(verts), np.array(faces), process=True)
        m.fix_normals()
        return trimesh.smoothing.filter_taubin(m, iterations=10) or m

    def larg_amnios(z, sy):
        t = np.clip((z - lo[2]) / L_mod, 0, 1)
        if dorsal is None:
            return sy                                             # profil seul : section ronde
        return max(demi_largeur(t) * MARGE_LARGEUR, sy * 0.5, 0.05)

    m_amn = loft(amn, larg_amnios) if amn is not None else None
    m_vit = loft(vit, lambda z, sy: sy) if vit is not None else None   # section ronde
    out = os.path.join("embryons_3D", "modeles", dossier + "_complements")
    os.makedirs(out, exist_ok=True)
    if m_amn is not None:
        m_amn.export(os.path.join(out, "amnios_trace.ply"))
    if m_vit is not None:
        m_vit.export(os.path.join(out, "vesicule_trace.ply"))
    elif os.path.exists(os.path.join(out, "vesicule_trace.ply")):
        os.remove(os.path.join(out, "vesicule_trace.ply"))
    info = {"stade": man.get("stade"), "source_trace": os.path.basename(psd_path), "echelle_mm_par_px": round(float(s), 6),
            "longueur_embryon_mm": round(float(L_mod), 3),
            "amnios": None if m_amn is None else {"fichier": "amnios_trace.ply", "etendue_mm": [round(float(v), 3) for v in m_amn.extents],
                       "source": "tracé de l'utilisateur (" + ("profil seul, section ronde" if dorsal is None else "profil + largeur dorsale ×1,1") + "), calé sur la longueur de l'embryon"},
            "vitellus": None if m_vit is None else {"fichier": "vesicule_trace.ply", "etendue_mm": [round(float(v), 3) for v in m_vit.extents],
                         "source": "tracé de l'utilisateur (profil, section ronde), calé sur la longueur de l'embryon"}}
    json.dump(info, open(os.path.join(out, "complements.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)

    # --- contrôle : contours des maillages reprojetés sur la photo de profil
    try:
        from PIL import Image, ImageDraw
        img = photo.copy(); dr = ImageDraw.Draw(img)
        for mk, col in ((trace, (0, 255, 0)),):                                                   # tracés de l'utilisateur, fins
            yy, xx = np.nonzero(mk[::1, ::1])
            for x_, y_ in list(zip(xx, yy))[::9]:
                dr.point((int(x_), int(y_)), fill=col)
        def vers_img(Pm):
            a = (Pm[:, 2] - c_mod[2]) / s; b = (Pm[:, 1] - c_mod[1]) / s
            return c_img + np.outer(a, ax) + np.outer(b, perp)
        for m, col in ([(m_amn, (80, 170, 255))] if m_amn is not None else []) + ([(m_vit, (255, 200, 90))] if m_vit is not None else []):
            for p in vers_img(m.vertices[::3]):
                dr.ellipse((float(p[0]) - 1.5, float(p[1]) - 1.5, float(p[0]) + 1.5, float(p[1]) + 1.5), fill=col)
        for q in P[::40]:                                                                        # embryon détecté (échelle)
            dr.point((float(q[0]), float(q[1])), fill=(255, 0, 255))
        img.save(os.path.join(out, "controle_trace.png"))
    except Exception as e:
        print("image de contrôle impossible :", e)
    print(json.dumps(info, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    a = sys.argv[3:]; o = {}
    for i, k in enumerate(a):
        if k == "--profil-seul":
            o["profil_seul"] = True
        elif k in ("--tete", "--queue"):
            o[k[2:]] = [float(v) for v in a[i + 1].split(",")]
        elif k in ("--dos", "--role"):
            o[k[2:]] = a[i + 1]
    main(sys.argv[1], sys.argv[2], o)
