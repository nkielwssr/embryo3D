"""
Topographie par region d'un stade du pipeline video (volume canonique work/labels.npz + work/dens.npy) :
axe vertebral courbe, colonne deroulee, etages vertebraux par periodicite, zones en coordonnees curvilignes.
Portage de la methode mise au point sur les coupes ehd (coupes embryos 9-23/axe_vertebral.py,
etages_vertebraux.py, coupes_medianes.py) ; ici aucune legende : tout vient des masques et de la densite.

  1. axe : chemin de moindre cout a l'interieur du masque colonne (squelette_axial_cartilage U notochorde),
     cout = 1 / (1 + distance au bord) -> le chemin suit l'axe medial ; de l'extremite craniale a l'extremite
     caudale ; spline fortement lissee, abscisse curviligne s (0 = extremite craniale = atlas)
  2. colonne deroulee : bandes sagittale et frontale courbes echantillonnees dans dens.npy (± demi_mm)
  3. etages : profil le long de s de la section du masque colonne (corps larges / disques etroits) et de
     la densite au centre ; periode par autocorrelation, disques = minima ; comptage depuis l'atlas ->
     C-1..Co ; fiabilite = nombre de disques trouves vs 24-32 et regularite du pas
  4. zones : cou / thorax / lombaire / pelvis / queue par les etages si fiables, sinon fractions de s
     du CS23 (marque « approximatif ») ; tete = au-dessus de l'extremite craniale de l'axe

Sorties : <dossier>/out/topographie/ : axe_vertebral.json, etages_vertebraux.json, zones_axe.json,
coupes_medianes.png, etages_vertebraux.png, <CS>_axe_vertebral.ply, <CS>_etages.ply (anneaux aux
frontieres), manifest_topographie.json (format lire_recon + bloc controle). Ne touche ni labels.npz
ni out/manifest.json.
Usage : python embryo3d/topographie_axe.py CS20_F4V [--demi-mm 1.5] [--lisser 15]
"""
import sys, os, re, json, argparse, collections
import numpy as np
from scipy import ndimage as ndi
from scipy.interpolate import UnivariateSpline
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from meshexport import mask_to_mesh
import trimesh

NIVEAUX = [f'C-{i}' for i in range(1, 8)] + [f'T-{i}' for i in range(1, 13)] + [f'L-{i}' for i in range(1, 6)] + \
          [f'S-{i}' for i in range(1, 6)] + [f'Co-{i}' for i in range(1, 5)]
FRACTIONS_CS23 = {'cou': (0.0, 0.224), 'thorax': (0.224, 0.645), 'lombaire': (0.645, 0.824), 'pelvis': (0.824, 1.0)}   # s / arc, CS23
AP_SIGNE = 1.0      # +1 : +ax2 = dorsal (u positif = dorsal) ; fixe par direction_ventrale() a chaque stade
COULEURS = {'tete': (255, 200, 120), 'cou': (120, 220, 255), 'thorax': (140, 255, 140), 'lombaire': (255, 255, 120), 'pelvis': (255, 160, 200), 'queue': (200, 200, 200)}


def charger(dossier):
    work = os.path.join(dossier, 'work'); out = os.path.join(dossier, 'out')
    z = np.load(os.path.join(work, 'labels.npz')); shape = tuple(int(v) for v in z['shape'])
    lab = lambda k: np.unpackbits(z[k])[:int(np.prod(shape))].reshape(shape).astype(bool) if k in z.files else None
    man = json.load(open(os.path.join(out, 'manifest.json'), encoding='utf-8'))
    dens = np.load(os.path.join(work, 'dens.npy'), mmap_mode='r')
    return shape, lab, man, dens


def bout_geodesique(mask, start, pont=6):
    """Voxel du masque le plus eloigne de `start` au sens geodesique : l'extremite caudale d'une queue enroulee n'est
    pas le point le plus bas en z. La distance est calculee dans le masque dilate de `pont` voxels (les masques ont des
    trous que le chemin sait franchir), mais l'extremite est choisie parmi les voxels du masque d'origine."""
    from skimage.graph import MCP_Geometric
    pontage = ndi.binary_dilation(mask, iterations=pont)
    cost = np.where(pontage, 1.0, np.inf)
    mcp = MCP_Geometric(cost, fully_connected=True)
    dist, _ = mcp.find_costs([tuple(int(v) for v in start)])
    dist = np.where(np.isfinite(dist) & mask, dist, -1)
    return np.unravel_index(int(np.argmax(dist)), dist.shape), float(dist.max())


def axe_par_chemin(colonne, ds=2):
    """Chemin medial dans le masque colonne, de l'extremite craniale (ax1 min) a l'extremite caudale (ax1 max)."""
    from skimage.graph import route_through_array
    idx = np.nonzero(colonne)
    lo = np.array([i.min() for i in idx]) - 6; hi = np.array([i.max() for i in idx]) + 7
    lo = np.maximum(lo, 0); hi = np.minimum(hi, colonne.shape)
    sub = colonne[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]
    # sous-echantillonnage (max par bloc) + fermeture pour relier les corps
    sh = tuple((np.array(sub.shape) + ds - 1) // ds * ds)
    pad = np.zeros(sh, bool); pad[:sub.shape[0], :sub.shape[1], :sub.shape[2]] = sub
    small = pad.reshape(sh[0] // ds, ds, sh[1] // ds, ds, sh[2] // ds, ds).any(axis=(1, 3, 5))
    small = ndi.binary_closing(small, iterations=2)
    d = ndi.distance_transform_edt(small); dout = ndi.distance_transform_edt(~small)
    cost = np.where(small, 1.0 / (1.0 + d), 3.0 + 0.3 * dout)   # hors masque : franchissable, de plus en plus cher en s'eloignant
    xs, ys, zs = np.nonzero(small)
    xmid = np.median(xs)
    # penalite laterale : les corps vertebraux sont medians, les arcs neuraux / cotes sont lateraux (le chemin
    # prenait la chaine des arcs, a 1,3 mm de la colonne des corps a CS20)
    gx = np.abs(np.arange(small.shape[0]) - xmid)[:, None, None] * ds
    cost = cost + 0.04 * gx
    cran = np.argmin(ys + 0.3 * np.abs(xs - xmid))
    start = (xs[cran], ys[cran], zs[cran]); end, _ = bout_geodesique(small, start)
    path, c = route_through_array(cost, start, end, fully_connected=True, geometric=True)
    P = np.array(path, float) * ds + lo + ds / 2
    return P, c


def queue_depuis_membres(membres, tube_full, vox, gl):
    """CS13-16 : la queue n'est pas dans le snc segmente mais dans le label 'membres' (appendiculaire). On en extrait les
    composantes MEDIANES (les membres sont lateraux) situees dans la moitie caudale du tube, pour prolonger le guide."""
    if membres is None or not membres.any():
        return None
    xs_t = np.nonzero(tube_full.any(axis=(1, 2)))[0]; xmid = float(np.median(np.nonzero(tube_full)[0])) if tube_full.any() else membres.shape[0] / 2
    zs_t = np.nonzero(tube_full.any(axis=(0, 2)))[0]; z_mid = float(np.median(zs_t)) if len(zs_t) else 0
    # la queue est fusionnee avec les bourgeons de membres dans le label : on prend la tranche MEDIANE du label
    # (|x - xmid| < 0,2 mm : les bourgeons sont lateraux, la queue est sur la ligne mediane), partie caudale seulement
    x = np.arange(membres.shape[0]); tranche = np.abs(x - xmid) * vox < 0.2
    med = membres & tranche[:, None, None]
    med[:, :int(z_mid), :] = False
    lbl, n = ndi.label(med); tailles = np.bincount(lbl.ravel()); tailles[0] = 0
    garde = tailles > 2000; garde[0] = False
    if not garde.any():
        return None
    q = garde[lbl]
    return ndi.binary_dilation(q, iterations=int(round(0.1 / vox)))    # re-epaissir (~0,1 mm) pour rejoindre la moelle


def axe_par_tube(snc, vent, vox, gl, ds=None, queue=None):
    """Axe guide par le tube neural (label snc, complet sur tous les stades) : on remplit la lumiere (ventricules) pour
    epaissir le cerveau, on ne garde que les parties fines (moelle : rayon < r_max), puis chemin medial de l'extremite
    caudale a l'extremite craniale de la partie fine (jonction bulbo-cervicale ~ atlas). Le cerveau est ainsi exclu."""
    from skimage.graph import route_through_array
    if ds is None:
        ds = max(1, int(round(0.05 / vox)))                      # ~50 um par voxel reduit
    filled = snc | vent if vent is not None else snc.copy()
    if queue is not None:
        filled = filled | queue
    idx = np.nonzero(filled)
    lo = np.array([i.min() for i in idx]) - 4; hi = np.array([i.max() for i in idx]) + 5
    lo = np.maximum(lo, 0); hi = np.minimum(hi, filled.shape)
    sub = filled[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]
    sh = tuple((np.array(sub.shape) + ds - 1) // ds * ds)
    pad = np.zeros(sh, bool); pad[:sub.shape[0], :sub.shape[1], :sub.shape[2]] = sub
    small = pad.reshape(sh[0] // ds, ds, sh[1] // ds, ds, sh[2] // ds, ds).any(axis=(1, 3, 5))
    small = ndi.binary_closing(small, iterations=2); small = ndi.binary_fill_holes(small)
    r_max = float(np.clip(0.035 * gl, 0.22, 0.55))
    # moelle = dans chaque coupe transversale (axe cranio-caudal) un blob compact d'aire < pi (2 r_max)^2 ;
    # le cerveau (rempli) et ses parois forment de grandes composantes 2D -> exclus. Une paroi fine mais etendue
    # (nappe) ne passe pas ce critere, contrairement au critere de rayon local.
    a_max = np.pi * (2.0 * r_max) ** 2 / (vox * ds) ** 2
    tube = np.zeros_like(small)
    for y in range(small.shape[1]):
        sl = small[:, y, :]
        if not sl.any():
            continue
        l2, n2 = ndi.label(sl); ar = np.bincount(l2.ravel()); ar[0] = 0
        ok = np.flatnonzero((ar > 3) & (ar <= a_max))
        if len(ok):
            tube[:, y, :] = np.isin(l2, ok)
    if queue is not None:                                     # la queue (epaisse) rejoint le tube fin : on la garde entiere
        q = queue[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]
        qp = np.zeros(sh, bool); qp[:q.shape[0], :q.shape[1], :q.shape[2]] = q
        tube |= qp.reshape(sh[0] // ds, ds, sh[1] // ds, ds, sh[2] // ds, ds).any(axis=(1, 3, 5))
    tube = ndi.binary_closing(tube, iterations=4 if queue is not None else 1)   # relie moelle fine et queue epaisse
    lbl, n = ndi.label(tube); tailles = np.bincount(lbl.ravel()); tailles[0] = 0
    grandes = tailles >= 0.1 * tailles.max(); grandes[0] = False
    print(f'    composantes du guide : {sorted(tailles[grandes].tolist(), reverse=True)} (voxels reduits), gardees')
    tube = grandes[lbl]
    dt = ndi.distance_transform_edt(tube); dout = ndi.distance_transform_edt(~tube)
    cost = np.where(tube, 1.0 / (1.0 + dt), 3.0 + 0.3 * dout)
    xs, ys, zs = np.nonzero(tube); xmid = np.median(xs)
    cost = cost + 0.04 * (np.abs(np.arange(tube.shape[0]) - xmid)[:, None, None] * ds)
    cran = np.argmin(ys + 0.3 * np.abs(xs - xmid))
    start = (xs[cran], ys[cran], zs[cran]); end, lg = bout_geodesique(tube, start)
    path, c = route_through_array(cost, start, end, fully_connected=True, geometric=True)
    P = np.array(path, float) * ds + lo + ds / 2
    return P, c, r_max, int(tube.sum())


def direction_ventrale(Q, T, coeur, vox):
    """Signe de l'axe AP (+ax2) : le coeur est ventral ; on regarde de quel cote de l'axe (point le plus proche) est le
    centroide du coeur. Retourne +1 si +ax2 = dorsal (coeur du cote -ax2), -1 sinon, et la distance signee."""
    if coeur is None or not coeur.any():
        return 1.0, None
    c = np.array(ndi.center_of_mass(coeur)) * vox
    j = int(np.argmin(np.linalg.norm(Q - c, axis=1)))
    ap = np.array([0, 0, 1.0]); u = ap - (ap @ T[j]) * T[j]; u /= np.linalg.norm(u) + 1e-9
    d = float((c - Q[j]) @ u)
    return (1.0 if d < 0 else -1.0), d


def lisser_axe(P, vox, lisser=15.0, pas_mm=0.02):
    """Spline lissee sur le chemin (parametre = longueur cumulee), reechantillonnee au pas pas_mm. Retourne points (mm), s, tangentes."""
    Pm = P * vox
    d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(Pm, axis=0), axis=1))]
    keep = np.r_[True, np.diff(d) > 1e-6]; Pm, d = Pm[keep], d[keep]
    spl = [UnivariateSpline(d, Pm[:, k], s=lisser * len(d) * 0.01, k=3) for k in range(3)]
    t = np.linspace(0, d[-1], 4000)
    Q = np.stack([f(t) for f in spl], 1)
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(Q, axis=0), axis=1))]
    ss = np.arange(0, s[-1], pas_mm)
    Q = np.stack([np.interp(ss, s, Q[:, k]) for k in range(3)], 1)
    T = np.gradient(Q, ss, axis=0); T /= np.linalg.norm(T, axis=1, keepdims=True)
    return Q, ss, T


def recentrer_sur_corps(Q, ss, T, dens, vox, u_min=-2.2, u_max=0.4, v_max=1.8, pas=0.05, sigma_mm=0.2):
    """Le masque colonne ne contient pas toujours les corps vertebraux (le chemin suit alors les arcs) : dans le plan
    perpendiculaire a chaque point, on cherche le blob le plus dense (densite lissee) dans une fenetre ventro-mediane
    et on deplace l'axe dessus ; decalages filtres (mediane sur ±0,5 mm de s)."""
    from scipy.ndimage import map_coordinates
    ap = np.array([0, 0, AP_SIGNE]); lr = np.array([1.0, 0, 0])
    uu = np.arange(u_min, u_max + 1e-9, pas); vv = np.arange(-v_max, v_max + 1e-9, pas)
    U, V = np.meshgrid(uu, vv, indexing='ij')
    du = np.zeros(len(ss)); dv = np.zeros(len(ss)); score = np.zeros(len(ss))
    G = []; cartes = []
    for i in range(0, len(ss)):
        u = ap - (ap @ T[i]) * T[i]; u /= np.linalg.norm(u) + 1e-9
        v = lr - (lr @ T[i]) * T[i]; v /= np.linalg.norm(v) + 1e-9
        Pts = (Q[i][None, :] + U.ravel()[:, None] * u[None, :] + V.ravel()[:, None] * v[None, :]) / vox
        g = map_coordinates(dens, Pts.T, order=1, mode='constant', cval=0).reshape(U.shape)
        cartes.append(ndi.gaussian_filter(g, sigma_mm / pas))
    # 1) decalage global : la chaine des corps est un blob brillant a offset a peu pres constant de l'axe des arcs
    Gm = np.mean(cartes, axis=0); k0 = np.unravel_index(np.argmax(Gm), Gm.shape); du0, dv0 = uu[k0[0]], vv[k0[1]]
    # 2) affinage local dans ± 0,35 mm autour du decalage global
    fen = (np.abs(U - du0) <= 0.35) & (np.abs(V - dv0) <= 0.35)
    for i, g in enumerate(cartes):
        gg = np.where(fen, g, -1); k = np.unravel_index(np.argmax(gg), gg.shape)
        du[i], dv[i], score[i] = uu[k[0]], vv[k[1]], g[k]
    w = int(0.5 / (ss[1] - ss[0]))
    du = ndi.median_filter(du, size=2 * w + 1, mode='nearest'); dv = ndi.median_filter(dv, size=2 * w + 1, mode='nearest')
    Q2 = Q.copy()
    for i in range(len(ss)):
        u = ap - (ap @ T[i]) * T[i]; u /= np.linalg.norm(u) + 1e-9
        v = lr - (lr @ T[i]) * T[i]; v /= np.linalg.norm(v) + 1e-9
        Q2[i] = Q[i] + du[i] * u + dv[i] * v
    return Q2, float(du0), float(dv0)


def profils(Q, ss, T, colonne, dens, vox, demi_mm, rayon_prof=0.15):
    """Section du masque colonne par tranche de s (voxels projetes sur l'axe) ; densite au centre ; bandes courbes."""
    from scipy.ndimage import map_coordinates
    arbre = cKDTree(Q)
    idx = np.nonzero(colonne); V = np.stack(idx, 1) * vox
    dist, j = arbre.query(V, distance_upper_bound=2.0)
    ok = np.isfinite(dist)
    aire = np.bincount(j[ok], minlength=len(ss)).astype(float) * vox ** 2 / 0.02 * vox   # ~ mm2 (voxels par bin / epaisseur du bin)
    # bandes courbes : u_sag = axe AP (ax2) orthogonalise a T ; v_fro = axe LR (ax0) orthogonalise
    ap = np.array([0, 0, AP_SIGNE]); lr = np.array([1.0, 0, 0])
    tt = np.arange(-demi_mm, demi_mm + 1e-9, vox)
    sag = np.zeros((len(ss), len(tt)), np.float32); fro = np.zeros_like(sag); centre = np.zeros(len(ss))
    # repere transporte parallelement le long de l'axe (initialise sur AP / LR) : un repere global orthogonalise a T devient
    # singulier quand la tangente passe par l'axe AP (queue enroulee), ce qui cassait les bandes
    u = ap - (ap @ T[0]) * T[0]; u /= np.linalg.norm(u) + 1e-9
    for i in range(len(ss)):
        u = u - (u @ T[i]) * T[i]; u /= np.linalg.norm(u) + 1e-9
        v = np.cross(T[i], u); v /= np.linalg.norm(v) + 1e-9
        if i == 0 and (v @ lr) < 0:
            v = -v
        if i == 0:
            signe_v = 1.0 if (np.cross(T[0], u) @ v) > 0 else -1.0
        v = signe_v * np.cross(T[i], u); v /= np.linalg.norm(v) + 1e-9
        Ps = (Q[i][None, :] + tt[:, None] * u[None, :]) / vox; Pf = (Q[i][None, :] + tt[:, None] * v[None, :]) / vox
        sag[i] = map_coordinates(dens, Ps.T, order=1, mode='constant', cval=0)
        fro[i] = map_coordinates(dens, Pf.T, order=1, mode='constant', cval=0)
        c = np.abs(tt) <= rayon_prof
        centre[i] = 0.5 * (sag[i][c].mean() + fro[i][c].mean())
    return aire, centre, sag, fro, tt


def chaine_segmentaire(Q, ss, T, dens, vox, rayon=2.5, pas=0.25, pas_attendu=0.56):
    """Cherche autour de l'axe (plan perpendiculaire, grille ± rayon mm) la position dont le profil de densite le long
    de s est le plus periodique (corps vertebraux / ganglions rachidiens / myotomes : une unite par etage).
    Retourne (offset (u, v), periode_mm, profil detendu, force)."""
    from scipy.ndimage import map_coordinates
    ap = np.array([0, 0, AP_SIGNE]); lr = np.array([1.0, 0, 0])
    uu = np.arange(-rayon, rayon + 1e-9, pas); vv = np.arange(-rayon, rayon + 1e-9, pas)
    U, V = np.meshgrid(uu, vv, indexing='ij')
    C = np.zeros((len(ss), len(uu), len(vv)), np.float32)
    for i in range(len(ss)):
        u = ap - (ap @ T[i]) * T[i]; u /= np.linalg.norm(u) + 1e-9
        v = lr - (lr @ T[i]) * T[i]; v /= np.linalg.norm(v) + 1e-9
        Pts = (Q[i][None, :] + U.ravel()[:, None] * u[None, :] + V.ravel()[:, None] * v[None, :]) / vox
        C[i] = ndi.gaussian_filter(map_coordinates(dens, Pts.T, order=1, mode='constant', cval=0).reshape(U.shape), 1.0)
    # fenetre de periode serree autour du pas attendu (sinon une harmonique ou une autre chaine l'emporte)
    ds_ = ss[1] - ss[0]; lag0, lag1 = max(2, int(0.7 * pas_attendu / ds_)), int(1.5 * pas_attendu / ds_)
    w_lis = max(1, int(0.12 * pas_attendu / ds_)); w_det = max(5, int(2.5 * pas_attendu / ds_))
    n_m = (5 * len(ss) // 6) - (len(ss) // 6); lag1 = min(lag1, n_m - 2)
    if lag1 <= lag0 + 2:
        return (0.0, 0.0), 0.0, np.zeros(len(ss)), 0.0, (np.zeros(U.shape), uu, vv)
    force = np.zeros(U.shape); per = np.zeros(U.shape)
    m = slice(len(ss) // 6, 5 * len(ss) // 6)
    for i in range(len(uu)):
        for j in range(len(vv)):
            x = ndi.gaussian_filter1d(C[:, i, j], w_lis) - ndi.gaussian_filter1d(C[:, i, j], w_det)
            x = x[m] - x[m].mean()
            if x.std() < 1e-6:
                continue
            ac = np.correlate(x, x, 'full')[len(x) - 1:]; ac /= ac[0]
            k = lag0 + int(np.argmax(ac[lag0:lag1]))
            if ac[k] > ac[k - 1] and ac[k] >= ac[k + 1]:
                force[i, j] = ac[k]; per[i, j] = k * ds_
    k = np.unravel_index(np.argmax(force), force.shape)
    prof = ndi.gaussian_filter1d(C[:, k[0], k[1]], w_lis) - ndi.gaussian_filter1d(C[:, k[0], k[1]], w_det)
    return (float(uu[k[0]]), float(vv[k[1]])), float(per[k]), prof, float(force[k]), (force, uu, vv)


def etages_par_periodicite(ss, aire, centre, arc, pas_attendu, prof=None, periode=None, force=0.0, arc_attendu=0.0, somitique=False):
    """Disques = minima de la section (lissee) ; periode par autocorrelation ; comptage depuis l'atlas."""
    # signal = profil de la chaine segmentaire la plus periodique (chaine_segmentaire) ; une unite (maximum) par etage,
    # frontieres = minima entre deux maxima
    a = prof if prof is not None else (ndi.gaussian_filter1d(centre, 3) - ndi.gaussian_filter1d(centre, 40))
    if not periode:
        periode = pas_attendu
    w = max(2, int(0.25 * periode / (ss[1] - ss[0])))
    mins = [i for i in range(w, len(a) - w) if a[i] <= a[i - w:i + w + 1].min()]
    disques = []
    for i in mins:
        if not disques or ss[i] - ss[disques[-1]] >= 0.55 * periode:
            disques.append(i)
    # calage predictif : phase de la grille (periode) estimee sur les minima detectes (mediane circulaire), grille etendue
    # a toute la bande, chaque frontiere predite calee sur le minimum detecte le plus proche (± 0,25 periode) s'il existe
    det = np.array([float(ss[i]) for i in disques])
    if len(det) >= 3 and periode > 0:
        # periode variable le long de s (les etages retrecissent vers la queue) : periode locale = mediane des ecarts
        # entre minima detectes dans une fenetre glissante de 3 mm, bornee a [0,6 ; 1,3] x periode globale, lissee
        ecarts = np.diff(det); mil = 0.5 * (det[1:] + det[:-1])
        def per_loc(s):
            sel = np.abs(mil - s) <= 1.5
            sel &= (ecarts >= 0.55 * periode) & (ecarts <= 1.6 * periode)
            return float(np.clip(np.median(ecarts[sel]), 0.6 * periode, 1.3 * periode)) if sel.sum() >= 2 else periode
        # phase sur la partie centrale (periodicite la plus sure), puis integration de la grille vers les deux bouts
        cen = det[(det >= 0.3 * arc) & (det <= 0.7 * arc)]
        cen = cen if len(cen) >= 2 else det
        ang = 2 * np.pi * (cen % periode) / periode
        s0 = float(cen[len(cen) // 2])
        g_bas = [s0]
        while g_bas[-1] + per_loc(g_bas[-1]) < arc - 0.3 * periode:
            g_bas.append(g_bas[-1] + per_loc(g_bas[-1]))
        g_haut = []
        while (g_haut[-1] if g_haut else s0) - per_loc((g_haut[-1] if g_haut else s0)) > 0.3 * periode:
            g_haut.append((g_haut[-1] if g_haut else s0) - per_loc((g_haut[-1] if g_haut else s0)))
        grille = np.array(sorted(g_haut + g_bas))
        b_int = []; origines_b = []
        for g in grille:
            j = int(np.argmin(np.abs(det - g))) if len(det) else -1
            if j >= 0 and abs(det[j] - g) <= 0.25 * per_loc(g):
                b_int.append(float(det[j])); origines_b.append('cale')
            else:
                b_int.append(float(g)); origines_b.append('predit')
        b = [0.0] + b_int + [float(arc)]; orig = ['debut'] + origines_b + ['fin']
    else:
        b = [0.0] + [float(ss[i]) for i in disques] + [float(arc)]; orig = ['debut'] + ['cale'] * len(disques) + ['fin']
    keep = [k for k, v in enumerate(b) if k == 0 or v - b[k - 1] >= 0.4 * periode]
    b = [b[k] for k in keep]; orig = [orig[k] for k in keep]
    n_cale = sum(o == 'cale' for o in orig); n_pred = sum(o == 'predit' for o in orig)
    etages = {}
    noms = [f'somite {k + 1}' for k in range(len(b))] if somitique else NIVEAUX      # somites : 1-4 occipitaux, 5-12 cervicaux...
    for k in range(min(len(b) - 1, len(noms))):
        etages[noms[k]] = {'s_debut_mm': round(b[k], 2), 's_fin_mm': round(b[k + 1], 2), 'hauteur_mm': round(b[k + 1] - b[k], 2),
                           'origine': f'{orig[k]} / {orig[k + 1]}'}
    n = len(b) - 1
    hs = np.diff(b[1:-1]) if len(b) > 3 else np.diff(b); cv = float(np.std(hs) / np.mean(hs)) if len(hs) > 1 else 1.0
    part_cale = n_cale / max(1, n_cale + n_pred)
    if somitique:   # somites : chaine bien visible en HREM ; periodicite forte, >= 12 somites ; frontieres calees a >= 50 %,
                    # ou >= 35 % si la chaine est tres reguliere (cv < 0,2) et tres periodique (force >= 0,6)
        fiable = n >= 12 and cv < 0.35 and force >= 0.4 and abs(periode - pas_attendu) / pas_attendu < 0.6 and arc >= 0.4 * arc_attendu             and (part_cale >= 0.5 or (part_cale >= 0.35 and cv < 0.2 and force >= 0.6))
    else:
        fiable = 20 <= n <= 41 and cv < 0.35 and abs(periode - pas_attendu) / pas_attendu < 0.5 and force >= 0.25 and arc >= 0.4 * arc_attendu and part_cale >= 0.5
    return etages, {'periode_mm': round(float(periode), 3), 'n_etages': n, 'cv_hauteurs': round(cv, 2), 'pas_attendu_mm': pas_attendu,
                    'force_periodicite': round(force, 2), 'frontieres_calees': n_cale, 'frontieres_predites': n_pred, 'fiable': bool(fiable)}


def racine_membre_superieur(membres, Q, S, vox):
    """s de la racine du membre superieur : composantes 'membres' les plus craniales (gauche + droite), s median des 5 %
    de voxels les plus proches de l'axe. Retourne (s, details) ou (None, ...)."""
    if membres is None or not membres.any():
        return None, None, 'pas de label membres'
    lbl, n = ndi.label(membres); tailles = np.bincount(lbl.ravel()); tailles[0] = 0
    comps = [k for k in np.argsort(tailles)[::-1][:6] if tailles[k] > 2000]
    arbre = cKDTree(Q); racines = []
    for k in comps:
        idx = np.nonzero(lbl == k); V = np.stack(idx, 1) * vox
        if len(V) > 150000:
            V = V[np.random.RandomState(0).choice(len(V), 150000, replace=False)]
        dist, j = arbre.query(V); q = np.percentile(dist, 5); sel = dist <= q
        s_r = float(np.median(S[j[sel]]))
        if s_r < 0.02 * S[-1] or s_r > 0.98 * S[-1]:       # composante qui se projette sur une extremite de l'axe (hors de son etendue) : pas une racine
            continue
        racines.append((s_r, float(q), int(tailles[k])))
    if not racines:
        return None, None, 'aucune composante'
    racines.sort()
    # membres superieurs = les plus craniaux ; on prend les composantes a moins de 25 % de l'arc de la plus craniale ;
    # membres inferieurs = les autres (plus caudales), s moyen
    s0 = racines[0][0]; sup = [r for r in racines if r[0] - s0 < 0.25 * S[-1]]; inf = [r for r in racines if r[0] - s0 >= 0.25 * S[-1]]
    s_sup = float(np.median([r[0] for r in sup])); s_inf = float(np.median([r[0] for r in inf])) if inf else None   # mediane : robuste a une composante parasite
    return s_sup, s_inf, (f'{len(sup)} composante(s) craniales, s = {[round(r[0], 2) for r in sup]}, distances {[round(r[1], 2) for r in sup]}'
                          + (f' ; membre inferieur : {len(inf)} composante(s), s = {[round(r[0], 2) for r in inf]}' if inf else ' ; pas de membre inferieur'))


def renommer_par_ancre(etages, s_ancre, somitique):
    """Recale le nommage : l'etage contenant la racine du membre superieur devient C-7 (vertebres) ou somite 11 (somites :
    membre superieur = somites 9-13, cervicaux = somites 5-12). Les etages en amont de C-1 deviennent 'occipital k'."""
    cles = list(etages)
    idx_ancre = next((i for i, k in enumerate(cles) if etages[k]['s_debut_mm'] <= s_ancre < etages[k]['s_fin_mm']), None)
    if idx_ancre is None:
        return etages, 0
    cible = 10 if somitique else 6                  # index (base 0) de somite 11 / de C-7 dans la liste des noms
    decalage = cible - idx_ancre                     # >0 : nos etages etaient nommes trop cranialement ; <0 : trop caudalement
    noms = [f'somite {k + 1}' for k in range(60)] if somitique else NIVEAUX
    nouveau = collections.OrderedDict()
    for i, k in enumerate(cles):
        j = i + decalage
        if j < 0:
            nom = f'occipital {abs(j)}' if abs(j) <= 4 else f'crane {abs(j) - 4}'   # au-dessus de C-1 : 4 sclerotomes occipitaux, puis base du crane
        elif j < len(noms):
            nom = noms[j]
        else:
            nom = f'caudal {j - len(noms) + 1}'
        e = dict(etages[k]); e['nom_avant_ancre'] = k
        nouveau[nom] = e
    return nouveau, decalage


NOMS_IMPOSES = [f'occipital {k}' for k in (1, 2, 3, 4)] + NIVEAUX[:29]      # 4 somites occipitaux (1 = le plus cranial, sous la vesicule otique, sans corps), C-1..C-7, T-1..T-12, L-1..L-5, S-1..S-5 (+ Co-1..Co-n)


def prolonger_axe(Q, ss, T, env, snc, vox, gl, arc, max_cran, max_caud, pas=0.04, tissu=None, exclure=None, cible_caud=None):
    """Prolonge l'axe aux deux bouts, a l'interieur de l'enveloppe, guide par le tissu nerveux (snc) quand il existe : a chaque pas
    le point avance le long de la tangente puis est tire vers le centroide des voxels snc (sinon enveloppe) dans un rayon
    proportionnel au stade. S'arrete en sortant de l'enveloppe ou a la longueur max. Retourne (Qx, ssx, Tx, n_cran, n_caud) ;
    ssx garde s = 0 au debut de l'axe d'origine (prolongement cranial en s negatif)."""
    r = float(np.clip(0.06 * gl, 0.15, 0.6))
    dse = max(1, int(round(0.05 / vox))); env_d = env[::dse, ::dse, ::dse]
    Ve = np.stack(np.nonzero(env_d), 1) * dse * vox; arbre_e = cKDTree(Ve)
    Vs = np.stack(np.nonzero(snc[::2, ::2, ::2]), 1) * 2 * vox if snc is not None and snc.any() else None
    arbre_s = cKDTree(Vs) if Vs is not None else None
    def dedans(pt):
        i = np.round(pt / vox).astype(int)
        return bool(np.all(i >= 0) and np.all(i < np.array(env.shape)) and env[tuple(i)])
    def marcher(p0, d0, max_len):
        pts = []; p = p0.copy(); d = d0.copy(); L = 0.0
        while L < max_len:
            q = p + pas * d
            if not dedans(q):
                break
            idx = arbre_s.query_ball_point(q, r) if arbre_s is not None else []
            if len(idx) >= 5:
                c = Vs[idx].mean(0)
            else:
                idx = arbre_e.query_ball_point(q, r); c = Ve[idx].mean(0) if len(idx) >= 5 else q
            q = 0.5 * (q + c)                                   # tire vers le tissu, sans sauter
            if not dedans(q):
                break
            dn = q - p; n = np.linalg.norm(dn)
            if n < 1e-6:
                break
            d = 0.7 * d + 0.3 * dn / n; d /= np.linalg.norm(d)
            L += n; p = q; pts.append(p.copy())
        return pts
    cran = marcher(Q[0], -T[0], max_cran)
    # caudal : chemin de moindre cout (enveloppe sous-echantillonnee, cout 1 dans le snc, 4 ailleurs) de la fin de l'axe au bout
    # geodesique de la queue (voxel de l'enveloppe le plus loin du sommet de la tete), retenu s'il prolonge bien l'axe
    caud = []
    if max_caud > 0:
        from skimage.graph import route_through_array
        if cible_caud is not None:                          # pointe pointee par l'utilisateur : chemin jusqu'a elle, pas de troncature
            tip_mm = np.array(cible_caud, float); tip_d = tuple(np.clip(np.round(tip_mm / vox / dse).astype(int), 0, np.array(env_d.shape) - 1))
            j_tip = len(Q) - 1; tissu = None; exclure = None; max_caud = 1e9
        else:
            ie = np.nonzero(env_d); top = int(np.argmin(ie[1]))
            tip_d, _ = bout_geodesique(env_d, (ie[0][top], ie[1][top], ie[2][top]), pont=2); tip_mm = np.array(tip_d, float) * dse * vox
            j_tip = int(np.argmin(np.linalg.norm(Q - tip_mm, axis=1)))
        if ss[j_tip] >= 0.85 * arc and np.linalg.norm(tip_mm - Q[-1]) > 2 * dse * vox:
            env_p = ndi.binary_dilation(env_d, iterations=1)
            snc_d = snc[::dse, ::dse, ::dse] if snc is not None else None
            cost = np.where(env_p, 4.0, 1e4);
            if snc_d is not None:
                cost[snc_d] = 1.0
            a = tuple(np.clip(np.round(Q[-1] / vox / dse).astype(int), 0, np.array(env_d.shape) - 1)); b_ = tuple(int(v) for v in tip_d)
            path, c_ = route_through_array(cost, a, b_, fully_connected=True, geometric=True)
            P_ = np.array(path, float) * dse * vox
            # le chemin ne vaut que tant qu'il reste dans le tissu axial (tube neural / notochorde / queue) : au-dela, il file vers
            # le membre inferieur ou le cordon (queue regressee a CS20). Arret a la premiere sortie de plus de 2 pas.
            # 1) jamais dans les membres ni le cordon (queue regressee : le bout geodesique est ailleurs) ; 2) si le tissu axial est
            #    etiquete au-dela de la fin de l'axe (premiers pas dans le tissu), on s'arrete a sa fin ; sinon (queue non etiquetee,
            #    stades jeunes) le chemin vaut jusqu'au bout de l'enveloppe
            n_ok = len(path)
            if exclure is not None and exclure.any():
                exc_d = ndi.binary_dilation(exclure[::dse, ::dse, ::dse], iterations=1)
                for k_, pt in enumerate(path):
                    if k_ > 2 and exc_d[tuple(int(v) for v in pt)]:
                        n_ok = min(n_ok, k_); break
            if tissu is not None and tissu.any():
                tis_d = ndi.binary_dilation(tissu[::dse, ::dse, ::dse], iterations=2)
                dans = [bool(tis_d[tuple(int(v) for v in pt)]) for pt in path]
                if sum(dans[1:16]) >= 10:                           # tissu axial etiquete nettement au-dela de la fin de l'axe (>= 10 pas)
                    hors = 0; n_t = 0
                    for k_, ok_ in enumerate(dans):
                        if ok_:
                            hors = 0; n_t = k_ + 1
                        else:
                            hors += 1
                            if hors > 2:
                                break
                    n_ok = min(n_ok, n_t)
            P_ = P_[:n_ok]
            if len(P_) > 3:
                P_ = ndi.gaussian_filter1d(P_, 2, axis=0)
            caud = [q for q in P_[1:] if np.linalg.norm(q - Q[-1]) > 1e-6]
            caud = caud[:int(max_caud / (dse * vox)) + 1]
    Qx = np.vstack([np.array(cran[::-1]).reshape(-1, 3), Q, np.array(caud).reshape(-1, 3)])
    d_ = np.r_[0, np.cumsum(np.linalg.norm(np.diff(Qx, axis=0), axis=1))]; ssx = d_ - d_[len(cran)]
    Tx = np.gradient(Qx, ssx, axis=0) if len(Qx) > 2 else np.tile(T[0], (len(Qx), 1)); Tx /= np.linalg.norm(Tx, axis=1, keepdims=True) + 1e-9
    return Qx, ssx, Tx, len(cran), len(caud)
REGION_DE = lambda nom: 'tete' if nom.startswith('occipital') else {'C': 'cou', 'T': 'thorax', 'L': 'lombaire', 'S': 'pelvis'}.get(nom[0], 'queue')


def imposer_niveaux(ss, Q, T, arc, s_sup, s_inf, s_tip, etages_det=None, n_co=4, pas_det=None, s_ot=None, ot_utilisateur=False):
    """Nombre de niveaux impose (consigne utilisateur) : 4 occipitaux, 7 cervicaux, 12 thoraciques, 5 lombaires, 5 sacres, n_co
    coccygiens, memes noms a tous les stades. Reperes le long de l'axe : centre de C-7 = racine du membre superieur (bourgeon C5-T1),
    centre de S-1 = racine du membre inferieur (L2-S2), pointe de la queue = fin de l'axe (ou bout video). Entre deux reperes, les
    niveaux sont repartis a abscisse curviligne reguliere ; chaque frontiere est ensuite calee sur la frontiere detectee la plus
    proche (marques utilisateur ou periodicite) si elle est a moins de 0,3 pas, sans changer le nombre. Retourne (etages, info, zones)."""
    N_TL = 18                                   # intervalles entre le centre de C-7 et le centre de S-1
    n_q = 4 + n_co                              # niveaux sous S-1 : S-2..S-5 + coccygiens
    s_tip = float(arc) if s_tip is None else float(s_tip)
    if s_sup is None and s_inf is None:
        p = s_tip / (6.5 + N_TL + 0.5 + n_q); s_sup = 6.5 * p; s_inf = s_sup + N_TL * p
        source = 'aucun repere de membre : C-1 au debut de l axe, pas uniforme jusqu a la pointe de la queue'
    elif s_inf is None:
        p = (s_tip - s_sup) / (N_TL + 0.5 + n_q); s_inf = s_sup + N_TL * p
        source = 'membre superieur (C-7) + pointe de la queue ; membre inferieur estime au pas uniforme'
    elif s_sup is None:
        p_q0 = (s_tip - s_inf) / (0.5 + n_q); s_sup = s_inf - N_TL * p_q0
        source = 'membre inferieur (S-1) + pointe de la queue ; membre superieur estime au pas uniforme'
    else:
        source = 'membre superieur (C-7), membre inferieur (S-1), pointe de la queue'
    p_tl = (s_inf - s_sup) / N_TL
    # membre inferieur implausible (queue trop courte derriere lui : pas caudal < 0,4 pas thoraco-lombaire, cas d'un bout de queue
    # etiquete 'membres' aux stades jeunes) : on l'ignore et on repartit uniformement de C-7 a la pointe
    if (s_tip - (s_inf - 0.5 * p_tl)) / (1 + n_q) < 0.4 * p_tl:
        p = (s_tip - s_sup) / (N_TL + 0.5 + n_q); s_inf = s_sup + N_TL * p; p_tl = p
        source += ' ; membre inferieur ecarte (queue trop courte derriere lui), S-1 estime au pas uniforme'
    # pas cervical / occipital : le pas detecte (somites ou corps deja formes) s'il existe et reste plausible, sinon le pas
    # thoraco-lombaire ; aux stades jeunes (queue courte, tous les niveaux pas encore segmentes) cela evite d'ecraser le cou
    p_c = float(pas_det) if pas_det and 0.5 * p_tl <= pas_det <= 2.0 * p_tl else p_tl
    c1 = s_sup - 6.5 * p_c                      # debut de C-1 (frontiere occipito-cervicale)
    s1 = s_inf - 0.5 * p_tl                     # debut de S-1
    p_q = (s_tip - s1) / (1 + n_q)              # pas sacro-coccygien : S-1 .. Co-n comprimes jusqu a la fin reelle de l'axe (jamais au-dela)
    # somites occipitaux : du bord caudal de la vesicule otique (s_ot) a C-1, repartis regulierement ; sans ancre otique, pas cervical
    if s_ot is not None and ((c1 - 8 * p_c) <= s_ot <= (c1 - 2 * p_c) or (ot_utilisateur and s_ot < c1 - 0.5 * p_c)):
        p_o = (c1 - s_ot) / 4.0; source_occ = ('vesicule otique pointee par l utilisateur' if ot_utilisateur else 'bord caudal des vesicules otiques') + ' -> C-1'
        if ot_utilisateur and not ((c1 - 8 * p_c) <= s_ot <= (c1 - 2 * p_c)):
            source_occ += f' (hors plage 2-8 pas cervicaux : {(c1 - s_ot) / p_c:.1f} pas)'
    else:
        p_o = p_c; source_occ = 'pas cervical (ancre otique absente ou hors plage 2-8 pas au-dessus de C-1' + (f' : s_ot = {s_ot:.2f}, C-1 = {c1:.2f})' if s_ot is not None else ')')
    noms = NOMS_IMPOSES + [f'Co-{k}' for k in range(1, n_co + 1)]
    b = [c1 - (4 - k) * p_o for k in range(4)] + [c1 + k * p_c for k in range(7)] + [s_sup + 0.5 * p_tl + k * p_tl for k in range(17)] \
        + [s1 + k * p_q for k in range(1 + n_q)] + [s_tip]
    b = np.array(b, float); fixe = np.zeros(len(b), bool)
    # calage sur les frontieres detectees (sans changer le nombre ni l ordre)
    calees = 0
    if etages_det:
        det = np.array(sorted({e['s_debut_mm'] for e in etages_det.values()} | {max(e['s_fin_mm'] for e in etages_det.values())}))
        for i in range(1, len(b) - 1):
            pas_loc = min(b[i] - b[i - 1], b[i + 1] - b[i]); j = int(np.argmin(np.abs(det - b[i])))
            if abs(det[j] - b[i]) <= 0.3 * pas_loc and b[i - 1] < det[j] < b[i + 1]:
                b[i] = det[j]; fixe[i] = True; calees += 1
    etages = collections.OrderedDict()
    for i, nom in enumerate(noms):
        sc = 0.5 * (b[i] + b[i + 1]); hors = not (ss[0] <= sc <= ss[-1]); prolong = not (0 <= sc <= arc)
        if sc < ss[0]:
            xyz = Q[0] + T[0] * (sc - ss[0])
        elif sc > ss[-1]:
            xyz = Q[-1] + T[-1] * (sc - ss[-1])
        else:
            xyz = Q[int(np.argmin(np.abs(ss - sc)))]
        origine = ('repere membre superieur' if b[i] <= s_sup < b[i + 1] else 'repere membre inferieur' if b[i] <= s_inf < b[i + 1]
                   else 'repere vesicule otique' if (nom == 'occipital 1' and s_ot is not None and p_o != p_c) else 'repere pointe' if nom == noms[-1] else 'reparti')
        etages[nom] = {'s_debut_mm': round(float(b[i]), 3), 's_fin_mm': round(float(b[i + 1]), 3), 'hauteur_mm': round(float(b[i + 1] - b[i]), 3),
                       's_centre_mm': round(float(sc), 3), 'xyz_centre_mm': [round(float(v), 3) for v in xyz], 'region': REGION_DE(nom),
                       'corps': not nom.startswith('occipital'), 'origine': origine, 'frontiere_calee': bool(fixe[i]), 'hors_axe': bool(hors),
                       'sur_prolongement': bool(prolong)}
    d = lambda k: etages[k]['s_debut_mm']
    zones = collections.OrderedDict([('tete', [None, d('C-1')]), ('cou', [d('C-1'), d('T-1')]), ('thorax', [d('T-1'), d('L-1')]),
                                     ('lombaire', [d('L-1'), d('S-1')]), ('pelvis', [d('S-1'), d('Co-1')]), ('queue', [d('Co-1'), round(float(s_tip), 3)])])
    info = {'mode': 'nombre impose', 'n_niveaux': len(etages), 'composition': f'4 occipitaux + 7 C + 12 T + 5 L + 5 S + {n_co} Co',
            'reperes': {'membre_superieur_s_mm': round(float(s_sup), 3), 'membre_inferieur_s_mm': round(float(s_inf), 3), 'pointe_queue_s_mm': round(float(s_tip), 3),
                        'vesicule_otique_bord_caudal_s_mm': round(float(s_ot), 3) if s_ot is not None else None, 'occipitaux': source_occ,
                        'debut_C1_s_mm': round(float(c1), 3), 'arc_axe_mm': round(float(arc), 3), 'axe_prolonge_s_mm': [round(float(ss[0]), 3), round(float(ss[-1]), 3)], 'source': source},
            'pas_occipital_mm': round(float(p_o), 3), 'pas_cervical_mm': round(float(p_c), 3), 'pas_thoraco_lombaire_mm': round(float(p_tl), 3), 'pas_sacro_coccygien_mm': round(float(p_q), 3), 'frontieres_calees': int(calees),
            'n_hors_axe': int(sum(e['hors_axe'] for e in etages.values())), 'n_sur_prolongement': int(sum(e['sur_prolongement'] for e in etages.values())),
            'n_corps': int(sum(e['corps'] for e in etages.values())),
            'principe': 'nombre de niveaux fixe et noms identiques a tous les stades : 4 somites occipitaux sans corps vertebral (occipital 1 sous la '
                        'vesicule otique) + 33 niveaux avec corps ; C-7 sur la racine du membre superieur, S-1 sur celle du membre inferieur, Co-4 finit '
                        'a la fin reelle de l axe (prolonge dans le tissu nerveux jusqu au bout de la queue, jamais au-dela) ; repartition reguliere entre '
                        'reperes puis calage des frontieres sur les frontieres detectees (< 0,3 pas) ; sur_prolongement = centre sur le prolongement de l axe'}
    return etages, info, zones


def zones_depuis(etages, arc, fiab, somitique=False):
    if fiab['fiable'] and somitique:
        # convention somitique : somites 1-4 occipitaux (tete), 5-12 cervicaux, 13-24 thoraciques, 25-29 lombaires,
        # 30-34 sacres, 35+ caudaux ; l'axe (moelle) commence a la jonction bulbo-cervicale ~ somite 5
        d = lambda k: etages[f'somite {k}']['s_debut_mm'] if f'somite {k}' in etages else arc
        b = collections.OrderedDict([('tete', [None, 0.0]), ('cou', [0.0, d(9)]), ('thorax', [d(9), d(21)]), ('lombaire', [d(21), d(26)]),
                                     ('pelvis', [d(26), d(31)]), ('queue', [d(31), arc])])
        return b, 'somites comptes depuis la jonction bulbo-cervicale (somite 5 suppose) : cou = 8 somites, thorax 12, lombaire 5, sacre 5'
    if fiab['fiable'] and all(k in etages for k in ('T-1', 'L-1')):
        # l'axe (masque colonne) peut s'arreter avant le sacrum : les zones absentes sont vides (bornes = fin d'axe)
        d = lambda k: etages[k]['s_debut_mm'] if k in etages else arc
        b = collections.OrderedDict([('tete', [None, 0.0]), ('cou', [0.0, d('T-1')]), ('thorax', [d('T-1'), d('L-1')]),
                                     ('lombaire', [d('L-1'), d('S-1')]), ('pelvis', [d('S-1'), d('Co-1')]), ('queue', [d('Co-1'), arc])])
        return b, 'etages' + ('' if 'S-1' in etages else ' (axe arrete avant le sacrum : pelvis / queue hors axe)')
    b = collections.OrderedDict([('tete', [None, 0.0])] + [(k, [round(f0 * arc, 2), round(f1 * arc, 2)]) for k, (f0, f1) in FRACTIONS_CS23.items()] + [('queue', [arc, arc])])
    return b, 'fractions CS23 (approximatif)'


def planches(out, cs, ss, tt, sag, fro, centre, aire, etages, fiab, zones, vox):
    import cv2
    def norm(b):
        b = np.clip(b, np.percentile(b, 1), np.percentile(b, 99)); b = (255 * (b - b.min()) / (b.max() - b.min() + 1e-9)).astype(np.uint8)
        return cv2.cvtColor(b, cv2.COLOR_GRAY2BGR)
    k = max(1, int(round(0.04 / vox)))                              # 1 px = 1 voxel ; hauteur d'une ligne = 0,02 mm
    S = norm(sag); F = norm(fro)
    band = np.hstack([S, np.full((len(ss), 4, 3), 255, np.uint8), F])
    band = cv2.resize(band, (band.shape[1] * 2, int(len(ss) * 0.02 / vox * 2)))          # isotrope, x2
    Hp, Wb = band.shape[:2]; mg, lp = 130, 220
    P = np.full((Hp + 40, mg + Wb + lp + 20, 3), 255, np.uint8); P[20:20 + Hp, mg:mg + Wb] = band
    ypx = lambda s: int(20 + s / vox * 2)
    for mm in range(0, int(ss[-1]) + 1):
        y = ypx(mm); cv2.line(P, (mg - 6, y), (mg, y), (0, 0, 0), 1)
        if mm % 2 == 0:
            cv2.putText(P, f'{mm} mm', (mg - 60, y + 4), cv2.FONT_HERSHEY_SIMPLEX, 0.4, (0, 0, 0), 1)
    for nom, e in etages.items():
        y = ypx(e['s_debut_mm'])
        if e.get('origine', '').startswith('predit'):
            for x in range(mg, mg + Wb, 8):
                cv2.line(P, (x, y), (x + 4, y), (0, 150, 0), 1)
        else:
            cv2.line(P, (mg, y), (mg + Wb, y), (0, 120, 220), 1)
        cv2.putText(P, nom, (mg + Wb + 4, ypx((e['s_debut_mm'] + e['s_fin_mm']) / 2) + 4), cv2.FONT_HERSHEY_SIMPLEX, 0.4, (0, 120, 220), 1)
    a = centre; a = (a - a.min()) / (a.max() - a.min() + 1e-9)
    x0 = mg + Wb + 60
    for i in range(1, len(ss)):
        cv2.line(P, (x0 + int(a[i - 1] * (lp - 80)), ypx(ss[i - 1])), (x0 + int(a[i] * (lp - 80)), ypx(ss[i])), (60, 60, 60), 1)
    cv2.putText(P, 'profil chaine segmentaire', (x0, 14), cv2.FONT_HERSHEY_SIMPLEX, 0.4, (60, 60, 60), 1)
    cv2.putText(P, 'sagittale courbe (ventral | dorsal)', (mg + 4, 14), cv2.FONT_HERSHEY_SIMPLEX, 0.4, (0, 0, 0), 1)
    cv2.putText(P, 'frontale courbe', (mg + S.shape[1] * 2 + 12, 14), cv2.FONT_HERSHEY_SIMPLEX, 0.4, (0, 0, 0), 1)
    cv2.putText(P, f"{cs} colonne deroulee (volume pipeline) ; periode {fiab['periode_mm']} mm, {fiab['n_etages']} etages, cv {fiab['cv_hauteurs']}, {'FIABLE' if fiab['fiable'] else 'NON FIABLE'} ; zones : {zones[1]}",
                (8, Hp + 34), cv2.FONT_HERSHEY_SIMPLEX, 0.42, (0, 0, 0), 1)
    cv2.imwrite(os.path.join(out, 'etages_vertebraux.png'), P)


def construire(dossier, demi_mm=1.5, lisser=15.0, guide='auto', label_tube=None):
    shape, lab, man, dens = charger(dossier)
    cs = man['stage']; vox = man['mm_per_voxel']; out = os.path.join(dossier, 'out', 'topographie'); os.makedirs(out, exist_ok=True)
    # masque « colonne » = cartilage axial U notochorde U somites (stades precoces : le cartilage axial est court ou absent,
    # la notochorde et les somites courent sur tout le tronc) ; la penalite laterale garde le chemin entre les somites
    # priorite au label 'corps_vertebraux' (livre par la session base) : chemin direct dans les corps, sans arcs ni cotes
    parts = {k: lab(k) for k in (('corps_vertebraux', 'notochorde') if lab('corps_vertebraux') is not None else ('squelette_axial_cartilage', 'notochorde', 'somites'))}
    parts = {k: v for k, v in parts.items() if v is not None}
    if not parts:
        raise SystemExit(f'{cs} : ni squelette_axial_cartilage, ni notochorde, ni somites dans labels.npz')
    colonne = np.zeros(shape, bool)
    for k, v in parts.items():
        colonne |= v
    print(f'{cs} : masque colonne = ' + ' U '.join(f'{k} ({int(v.sum())} vox)' for k, v in parts.items()))
    # la colonne segmentee est morcelee (corps separes, trous entre blocs) : fermeture, suppression des miettes,
    # et le chemin est autorise a franchir les trous (cout croissant avec la distance au masque)
    colonne = ndi.binary_closing(colonne, iterations=4)
    lbl, n = ndi.label(colonne); tailles = np.bincount(lbl.ravel()); tailles[0] = 0
    keep = tailles >= 500; colonne = keep[lbl]          # bool (l'indexation par tailles[lbl] creait un int64 de 1 Go)
    print(f'{cs} : masque colonne {colonne.sum()} voxels ({n} composantes, miettes < 500 voxels retirees), voxel {vox * 1000:.1f} um')
    gl = man.get('greatest_length_mm_assumed') or 20.0
    somitique = gl < 10.5                                       # CS13-16 (somites) ; CS17 et plus : corps vertebraux
    # guide tube neural : 'moelle_rachidienne' (moelle amputee au-dessus de C-1, session base) > 'moelle' > 'snc'
    # guide tube neural : le premier label (moelle amputee > moelle > snc) dont l'etendue cranio-caudale atteint 85 % de la plus
    # grande (la base peut amputer 'moelle_rachidienne' bien plus bas que C-1 : CS20 a 16:42, 7,6-18 mm au lieu de 2,7-18)
    cand = {k: lab(k) for k in ('moelle_rachidienne', 'moelle', 'snc') if lab(k) is not None}
    if label_tube and lab(label_tube) is not None:
        snc = lab(label_tube); print(f'  guide tube neural force sur le label {label_tube}')
    elif cand:
        ext0 = lambda m: float(np.ptp(np.nonzero(m.any(axis=(0, 2)))[0])) * vox if m.any() else 0.0
        exts = {k: ext0(m) for k, m in cand.items()}; e_max = max(exts.values())
        k_sel = next(k for k in cand if exts[k] >= 0.85 * e_max); snc = cand[k_sel]
        print(f'  guide tube neural : label {k_sel} (etendues ' + ', '.join(f'{k} {v:.1f} mm' for k, v in exts.items()) + ')')
    else:
        snc = None
    ext = lambda m: float(np.ptp(np.nonzero(m.any(axis=(0, 2)))[0])) * vox if m is not None and m.any() else 0.0
    guide_tube = (guide == 'snc') or (guide == 'auto' and snc is not None and (somitique or ext(colonne) < 0.6 * ext(snc)))
    cents = None
    if guide == 'niveaux':
        fnv = os.path.join(dossier, 'work', 'vertebres_niveaux.npz')
        if not os.path.exists(fnv):
            raise SystemExit(f'{cs} : pas de work/vertebres_niveaux.npz (recrutement par niveau de la session base)')
        zn = np.load(fnv); N = zn['niveaux_corps'] if 'niveaux_corps' in zn.files else zn['niveaux']   # corps seuls si la base les fournit
        ids = [int(k) for k in np.unique(N) if k > 0]
        cents = np.array(ndi.center_of_mass(N > 0, N, ids), float)          # voxels, dans l'ordre des identifiants
        # les identifiants ne sont pas forcement ordonnes et le fichier peut contenir des instances parasites (chaine dorsale des
        # arcs, miettes) : on apparie chaque etage manuel (video 360) au niveau le plus proche, dans l'ordre de la chaine manuelle
        fman = os.path.join(out, 'video360', 'etages_manuels.json')
        if os.path.exists(fman):
            jm = json.load(open(fman, encoding='utf-8')); em = jm['etages']
            Pm = np.array([e['xyz_mm_pipeline'] for e in em]) / vox
            # X approche (vue video recalee sur une projection tournee, CS13/CS16) : recentrage de chaque point sur le plan
            # median local = mediane en ax0 des voxels de l'enveloppe de la tranche ax1, dans une fenetre AP de +-0,5 mm
            theta_m = jm.get('projection_theta_deg', 0) or any(v.get('projection_theta_deg', 0) for v in jm.get('fusion', {}).values() if isinstance(v, dict))
            env = lab('enveloppe')
            if theta_m and env is not None:
                w = max(2, int(round(0.5 / vox))); x_avant = Pm[:, 0].copy()
                for i_, pm in enumerate(Pm):
                    y0 = int(np.clip(round(pm[1]), 0, shape[1] - 1)); z0, z1 = int(max(0, pm[2] - w)), int(min(shape[2], pm[2] + w + 1))
                    tr = env[:, max(0, y0 - 2):y0 + 3, z0:z1]
                    if tr.any():
                        Pm[i_, 0] = float(np.median(np.nonzero(tr)[0]))
                print(f'  X des etages manuels recentre sur le plan median de l enveloppe (projection tournee) : correction mediane {np.median(np.abs(Pm[:, 0] - x_avant)) * vox:.2f} mm, max {np.abs(Pm[:, 0] - x_avant).max() * vox:.2f} mm')
            pas_m = float(np.median(np.linalg.norm(np.diff(Pm, axis=0), axis=1)))
            choisis = []
            for pm in Pm:
                dm = np.linalg.norm(cents - pm, axis=1); j = int(np.argmin(dm))
                if dm[j] <= 0.75 * pas_m and j not in choisis:
                    choisis.append(j)
            print(f'  niveaux apparies aux etages manuels : {len(choisis)} / {len(ids)} instances a moins de 0,75 pas ({len(em)} etages manuels)')
            # l'axe passe par les etages manuels eux-memes (reference de l'utilisateur, chaine complete jusqu'a la queue) ;
            # les instances recrutees ne servent qu'au controle d'appariement
            cents = Pm.copy(); ids = list(range(1, len(em) + 1)); n_app = len(choisis)
        else:
            ordre = [int(np.argmin(cents[:, 1]))]                             # chaine plus proche voisin depuis le plus cranial
            while len(ordre) < len(ids):
                dm = np.linalg.norm(cents - cents[ordre[-1]], axis=1); dm[ordre] = np.inf; ordre.append(int(np.argmin(dm)))
            cents = cents[ordre]; ids = [ids[j] for j in ordre]; n_app = None
        P = []
        for a_, b_ in zip(cents[:-1], cents[1:]):                             # chemin dense (pas ~ 1 voxel) entre centroides
            P += list(np.linspace(a_, b_, max(2, int(np.linalg.norm(b_ - a_))), endpoint=False))
        P.append(cents[-1]); P = np.array(P); cout = 0.0
        source_axe = (f'etages manuels de la video 360 ({len(ids)} points, {n_app} niveaux recrutes apparies)' if n_app is not None
                      else f'centroides des {len(ids)} niveaux vertebraux recrutes (work/vertebres_niveaux.npz, session base)')
        print(f'  guide = niveaux recrutes : {len(ids)} niveaux, de ax1 {cents[0, 1] * vox:.2f} a {cents[-1, 1] * vox:.2f} mm')
    if cents is not None:
        pass
    elif guide_tube:
        queue = queue_depuis_membres(lab('membres'), snc, vox, gl) if gl < 14 else None     # CS13-16 : queue dans 'membres'
        if queue is not None:
            print(f'  queue ajoutee au guide depuis le label membres : {int(queue.sum())} voxels')
        P, cout, r_max, n_tube = axe_par_tube(snc, lab('ventricules'), vox, gl, queue=queue)
        print(f'  guide = tube neural (colonne {ext(colonne):.1f} mm vs snc {ext(snc):.1f} mm) : parties fines r < {r_max:.2f} mm, {n_tube} voxels reduits')
        source_axe = f'chemin medial dans le tube neural (snc rempli par les ventricules, rayon < {r_max:.2f} mm = moelle, cerveau exclu)'
    else:
        P, cout = axe_par_chemin(colonne); source_axe = 'chemin medial dans squelette_axial_cartilage U notochorde U somites'
    Q, ss, T = lisser_axe(P, vox, lisser)
    # (pas de recentrage sur la densite : l'axe medial du masque colonne suffit pour la topographie ;
    #  les etages viennent de la chaine segmentaire la plus periodique autour de l'axe, cf. chaine_segmentaire)
    arc = float(ss[-1]); z_ext = float(Q[:, 1].max() - Q[:, 1].min())
    global AP_SIGNE
    coeur = lab('coeur_detoure') if lab('coeur_detoure') is not None else lab('coeur')
    AP_SIGNE, d_coeur = direction_ventrale(Q, T, coeur, vox)
    print(f"  orientation AP : coeur a {d_coeur:+.2f} mm de l'axe le long de +ax2 -> {'+ax2 = dorsal' if AP_SIGNE > 0 else '+ax2 = ventral (u inverse)'}" if d_coeur is not None else '  orientation AP : pas de coeur, +ax2 = dorsal suppose')
    print(f'  axe : {len(P)} voxels de chemin, arc {arc:.2f} mm pour {z_ext:.2f} mm cranio-caudal, ecart lateral {np.ptp(Q[:, 0]):.2f} mm')
    demi_mm = float(min(demi_mm, max(0.4, 0.075 * gl)))       # bandes ± 1,5 mm a CS20, ± 0,4 mm a CS13
    aire, centre, sag, fro, tt = profils(Q, ss, T, colonne, dens, vox, demi_mm)
    pas_attendu = round(0.76 * gl / 31.0 * 1.15, 3)          # hauteur T8 CS23 x GL, +15 % (CS23 est petit pour sa taille en z)
    if somitique:                                               # stades somitiques : longueur d'un somite mesuree ~ 0,047 x GL (CS14 : 0,28 mm)
        pas_attendu = round(0.042 * gl, 3)
    rayon_ch = float(np.clip(12 * pas_attendu, 0.6, 2.5)); pas_ch = float(max(0.05, min(0.25, pas_attendu / 2)))
    (ou, ov), periode, prof, force, carte = chaine_segmentaire(Q, ss, T, dens, vox, rayon=rayon_ch, pas=pas_ch, pas_attendu=pas_attendu)
    print(f'  chaine segmentaire la plus periodique : decalage ventral {ou:+.2f} mm, lateral {ov:+.2f} mm, periode {periode:.2f} mm, force {force:.2f}')
    arc_attendu = 1.2 * gl if somitique else 20.2 * gl / 31.0   # tube neural en C (somitique) ; sinon arc CS23 x GL
    etages, fiab = etages_par_periodicite(ss, aire, centre, arc, pas_attendu, prof, periode, force, arc_attendu, somitique)
    if cents is not None:
        arb_q = cKDTree(Q); Sc = np.array([float(ss[arb_q.query(c_ * vox)[1]]) for c_ in cents])
        ordre = np.argsort(Sc); Sc = Sc[ordre]
        bornes = [0.0] + [0.5 * (Sc[i] + Sc[i + 1]) for i in range(len(Sc) - 1)] + [float(arc)]
        noms = [f'somite {k + 1}' for k in range(80)] if somitique else NIVEAUX + [f'caudal {k}' for k in range(1, 40)]
        etages = collections.OrderedDict((noms[i], {'s_debut_mm': round(bornes[i], 2), 's_fin_mm': round(bornes[i + 1], 2),
                                                    'hauteur_mm': round(bornes[i + 1] - bornes[i], 2), 'origine': 'niveau recrute (etages manuels)'})
                                         for i in range(len(Sc)))
        h = np.diff(Sc)
        fiab.update({'fiable': True, 'n_etages': int(len(Sc)), 'periode_mm': round(float(np.median(h)), 3), 'cv_hauteurs': round(float(h.std() / h.mean()), 2),
                     'frontieres_calees': int(len(Sc)), 'frontieres_predites': 0,
                     'source_etages': source_axe})
    fiab['guide'] = 'niveaux recrutes' if cents is not None else ('tube neural' if guide_tube else 'colonne'); fiab['somitique'] = somitique
    fiab['arc_mm'] = round(arc, 2); fiab['arc_attendu_mm'] = round(arc_attendu, 1)
    if arc < 0.4 * arc_attendu:
        print(f'  AVERTISSEMENT : axe partiel ({arc:.1f} mm pour ~{arc_attendu:.1f} attendus) : masques colonne incomplets')
    fiab['chaine_offset_mm'] = [round(ou, 2), round(ov, 2)]
    fiab['nommage'] = 'compte depuis l extremite craniale du masque colonne = C-1 suppose : incertitude +-1 a 2 niveaux'
    centre = prof
    # ancrage anatomique du nommage : racine du membre superieur = C-7 (ou somite 11)
    s_ancre, s_inf, det_ancre = racine_membre_superieur(lab('membres'), Q, ss, vox)
    decalage = 0
    if s_ancre is not None and fiab['fiable']:
        etages, decalage = renommer_par_ancre(etages, s_ancre, somitique)
        print(f"  ancrage : racine du membre superieur a s = {s_ancre:.2f} mm ({det_ancre}) -> decalage de nommage {decalage:+d} niveau(x)")
    fiab['ancre_membre_sup_s_mm'] = round(s_ancre, 2) if s_ancre is not None else None
    fiab['decalage_nommage'] = decalage
    fiab['nommage'] = ('racine du membre superieur = C-7 (somite 11 pour les stades somitiques) ; ' + ('etages en amont = occipitaux' if decalage > 0 else 'premier etage = C-1 + ' + str(-decalage))) if s_ancre is not None else fiab['nommage']
    zones = zones_depuis(etages, arc, fiab, somitique)
    print(f"  etages : periode {fiab['periode_mm']} mm (attendu ~{pas_attendu}), {fiab['n_etages']} etages, cv {fiab['cv_hauteurs']} -> {'fiable' if fiab['fiable'] else 'NON fiable'} ; zones par {zones[1]}")
    planches(out, cs, ss, tt, sag, fro, centre, aire, etages, fiab, zones, vox)
    os.replace(os.path.join(out, 'etages_vertebraux.png'), os.path.join(out, 'etages_detectes.png'))
    etages_det, fiab_det, zones_det = etages, dict(fiab), zones
    # --- nombre de niveaux impose (consigne utilisateur) : reperes = racines des membres + pointe de la queue
    # prolongement de l'axe aux deux bouts dans le tissu nerveux (jusqu'a la vesicule otique en haut, jusqu'au bout de la queue en bas)
    env_ = lab('enveloppe'); Qx, ssx, Tx = Q, ss, T; s_ot = None
    rep_u = None; ot_ok = q_ok = False
    fru = os.path.join(out, 'video360', 'reperes_utilisateur.json')
    if os.path.exists(fru):
        rep_u = json.load(open(fru, encoding='utf-8'))
        if 'vesicule_otique' in rep_u and 'vesicule_optique' not in rep_u:      # erratum utilisateur : le cercle est l oeil
            rep_u['vesicule_optique'] = rep_u.pop('vesicule_otique')
        ot_ok = False                                                         # l'oeil n'est pas l ancre otique : repere cranien seulement
        oeil_ok = bool(rep_u['vesicule_optique'].get('valide', False)); q_ok = bool(rep_u['pointe_queue'].get('valide', False))
        if not (oeil_ok or q_ok):
            print('  reperes utilisateur presents mais non valides (reperes_video.py --verifier) : ancres automatiques conservees'); rep_u = None
    if env_ is not None:
        snc_plein = lab('snc') if lab('snc') is not None else snc      # snc complet (la moelle amputee s'arrete a C-1)
        tissu_ax = np.zeros(shape, bool)                                 # tissu axial pour borner le prolongement caudal
        for k_t in ('snc', 'moelle_rachidienne', 'notochorde', 'squelette_axial_cartilage'):
            if lab(k_t) is not None:
                tissu_ax |= lab(k_t)
        if gl < 14 and lab('membres') is not None and snc is not None:
            q_ = queue_depuis_membres(lab('membres'), snc, vox, gl)
            if q_ is not None:
                tissu_ax |= q_
        excl = np.zeros(shape, bool)
        for k_e in ('membres', 'cordon_ombilical'):
            if lab(k_e) is not None:
                excl |= lab(k_e)
        cible = np.array(rep_u['pointe_queue']['xyz_mm_pipeline']) if (rep_u is not None and q_ok) else None
        Qx, ssx, Tx, n_cr, n_cd = prolonger_axe(Q, ss, T, env_, snc_plein, vox, gl, arc, max_cran=0.35 * arc, max_caud=0.6 * arc, tissu=tissu_ax, exclure=excl, cible_caud=cible)
        if cible is not None:
            print(f'  prolongement caudal dirige vers la pointe de queue pointee par l utilisateur ({n_cd} pts, fin a {np.linalg.norm(Qx[-1] - cible):.2f} mm de la marque)')
        print(f'  axe prolonge : {ssx[0]:+.2f} mm en cranial ({n_cr} pts), +{ssx[-1] - arc:.2f} mm en caudal ({n_cd} pts) -> s de {ssx[0]:.2f} a {ssx[-1]:.2f} mm')
    otiq = lab('vesicules_otiques')
    if otiq is not None and otiq.any():                     # ancre occipitale : bord caudal des vesicules otiques projete sur l'axe
        # projection curviligne : s ou la vesicule est dans le plan perpendiculaire a l'axe ((V - Q(s)) . T(s) = 0), partie craniale
        # seulement (une tete flechie rapproche les vesicules du cou en distance euclidienne) ; bord caudal = centre + rayon equivalent
        lo, no = ndi.label(otiq[::2, ::2, ::2]); tailles_o = np.bincount(lo.ravel()); tailles_o[0] = 0
        sel_c = np.nonzero(ssx <= 0.25 * arc)[0]; s_cands = []
        for k_o in np.argsort(tailles_o)[::-1][:2]:            # les deux plus grosses composantes = gauche / droite
            if tailles_o[k_o] < 20:
                continue
            co = np.array(ndi.center_of_mass(lo == k_o)) * 2 * vox; r_o = (3 * tailles_o[k_o] * (2 * vox) ** 3 / (4 * np.pi)) ** (1 / 3)
            dd = np.array([(co - Qx[i]) @ Tx[i] for i in sel_c]); dist = np.linalg.norm(Qx[sel_c] - co, axis=1)
            z_ = np.nonzero(np.diff(np.sign(dd)) != 0)[0]
            if len(z_):
                i_best = z_[int(np.argmin(dist[z_]))]; s_cands.append(float(ssx[sel_c[i_best]]) + r_o)
        s_ot = float(np.mean(s_cands)) if s_cands else None
        if s_ot is not None:
            print(f'  vesicules otiques : bord caudal (projection curviligne) a s = {s_ot:.2f} mm ({len(s_cands)} vesicule(s))')
        else:
            print('  vesicules otiques : pas de projection curviligne sur la partie craniale de l axe')
    s_tip = float(ssx[-1]); src_ot = 'label vesicules_otiques'; src_tip = 'fin de l axe prolonge'
    s_oeil = None
    if rep_u is not None and oeil_ok:
        # oeil (vesicule optique) : repere cranien = abscisse du point d'axe prolonge le plus proche (partie craniale)
        po = np.array(rep_u['vesicule_optique']['xyz_mm_pipeline']); sel_c = np.nonzero(ssx <= 0.25 * arc)[0]
        jo = sel_c[int(np.argmin(np.linalg.norm(Qx[sel_c] - po, axis=1)))]; s_oeil = float(ssx[jo])
    if rep_u is not None and q_ok:
        src_tip = 'utilisateur (queue-VO)'                    # l'axe prolonge se termine sur la marque (cible du chemin caudal)
    if rep_u is not None:
        print(f"  reperes utilisateur : oeil {'a s = %.2f mm (repere cranien)' % s_oeil if oeil_ok else 'NON valide'}, pointe de la queue {'a s = %.2f mm' % s_tip if q_ok else 'NON valide (fin d axe conservee)'}")
    etages, info_imp, zones_b = imposer_niveaux(ssx, Qx, Tx, arc, s_ancre, s_inf, s_tip, etages_det, n_co=4, pas_det=fiab_det.get('periode_mm') if fiab_det.get('fiable') else None,
                                                s_ot=s_ot, ot_utilisateur=False)
    info_imp['reperes']['source_vesicule_otique'] = src_ot; info_imp['reperes']['source_pointe_queue'] = src_tip
    info_imp['reperes']['oeil_utilisateur_s_mm'] = round(s_oeil, 3) if s_oeil is not None else None
    if rep_u is not None:
        info_imp['reperes']['utilisateur'] = {'fichier': 'video360/reperes_utilisateur.json', 'image_video': rep_u.get('image_video'), 'vue': rep_u.get('vue'),
                                              'vesicule_optique_xyz_mm': rep_u['vesicule_optique']['xyz_mm_pipeline'], 'pointe_queue_xyz_mm': rep_u['pointe_queue']['xyz_mm_pipeline'],
                                              'note': 'oeil = repere cranien (pas l ancre otique) ; pointe de queue = fin de Co-4'}
    zones = (zones_b, 'niveaux imposes (C-7 = membre superieur, S-1 = membre inferieur, Co-4 = pointe de la queue)')
    fiab = dict(fiab_det); fiab.update({'fiable': True, 'mode': 'nombre impose', 'n_etages': info_imp['n_niveaux'], 'periode_mm': info_imp['pas_thoraco_lombaire_mm'],
                                        'cv_hauteurs': 0.0, 'frontieres_calees': info_imp['frontieres_calees'], 'frontieres_predites': info_imp['n_niveaux'] + 1 - info_imp['frontieres_calees'],
                                        'reperes': info_imp['reperes'], 'nommage': info_imp['principe'], 'source_etages': 'niveaux imposes ; detection conservee dans etages_detectes.json'})
    print(f"  niveaux imposes : {info_imp['n_niveaux']} ({info_imp['composition']}), otique {info_imp['reperes']['vesicule_otique_bord_caudal_s_mm']}, C-1 a s = {info_imp['reperes']['debut_C1_s_mm']}, C-7 a s = {info_imp['reperes']['membre_superieur_s_mm']}, S-1 a s = {info_imp['reperes']['membre_inferieur_s_mm']}, "
          f"pointe {info_imp['reperes']['pointe_queue_s_mm']} mm ; pas occ {info_imp['pas_occipital_mm']}, C {info_imp['pas_cervical_mm']}, TL {info_imp['pas_thoraco_lombaire_mm']}, SCo {info_imp['pas_sacro_coccygien_mm']} mm ; "
          f"{info_imp['frontieres_calees']} frontiere(s) calee(s), {info_imp['n_hors_axe']} hors axe, {info_imp['n_sur_prolongement']} sur prolongement ; reperes : {info_imp['reperes']['source']} ; occipitaux : {info_imp['reperes']['occipitaux']}")
    if ssx[-1] > arc + 0.05:                                  # bandes prolongees jusqu'au bout de la queue pour la planche des niveaux imposes
        ssr = np.arange(0, ssx[-1], ss[1] - ss[0]); Qr = np.stack([np.interp(ssr, ssx, Qx[:, k_]) for k_ in range(3)], 1)
        Tr = np.gradient(Qr, ssr, axis=0); Tr /= np.linalg.norm(Tr, axis=1, keepdims=True) + 1e-9
        aire_r, centre_r, sag_r, fro_r, tt_r = profils(Qr, ssr, Tr, colonne, dens, vox, demi_mm)
        planches(out, cs, ssr, tt_r, sag_r, fro_r, centre_r, aire_r, collections.OrderedDict((k, e) for k, e in etages.items() if 0 <= e['s_centre_mm'] <= ssr[-1]), fiab, zones, vox)
    else:
        planches(out, cs, ss, tt, sag, fro, centre, aire, collections.OrderedDict((k, e) for k, e in etages.items() if 0 <= e['s_centre_mm'] <= arc), fiab, zones, vox)
    # --- decalage ventral de l'axe exporte : le guide (moelle ou arcs) est dorsal aux corps vertebraux ; on place l'axe
    # affiche sur la chaine des corps (offset ventral de la chaine segmentaire si elle est mediane) ou, a defaut, sur la
    # notochorde (offset median des voxels notochorde par rapport au guide), en restant sur la ligne mediane (lateral = 0)
    ap = np.array([0, 0, AP_SIGNE]); noto = parts.get('notochorde')
    dec_u = None; source_dec = ''
    if cents is not None or ('corps_vertebraux' in parts and not guide_tube):
        dec_u = 0.0; source_dec = 'axe deja dans les corps vertebraux (label corps_vertebraux)'
    elif ou < 0 and abs(ov) <= 0.6:
        dec_u = float(ou); source_dec = 'chaine segmentaire mediane (corps vertebraux)'
    elif noto is not None and noto.any():
        idxn = np.nonzero(noto); Vn = np.stack(idxn, 1) * vox
        arb = cKDTree(Q); dn, jn = arb.query(Vn, distance_upper_bound=2.0); okn = np.isfinite(dn)
        if okn.sum() > 100:
            us = []  # offset notochorde ; s'il ressort dorsal (>= 0) le label est douteux (plancher du tube) -> repli
            for i_, v_ in zip(jn[okn], Vn[okn]):
                u = ap - (ap @ T[i_]) * T[i_]; u /= np.linalg.norm(u) + 1e-9
                us.append(float((v_ - Q[i_]) @ u))
            dec_u = float(np.median(us)); source_dec = 'notochorde (offset median)'
            if dec_u >= 0:
                dec_u = None; source_dec = ''
    if dec_u is None:
        dec_u = -float(np.clip(0.035 * gl, 0.22, 0.55)) - 0.05; source_dec = 'valeur par defaut (rayon de la moelle + 0,05 mm)'
    Qv = Q.copy()
    for i in range(len(ss)):
        u = ap - (ap @ T[i]) * T[i]; u /= np.linalg.norm(u) + 1e-9
        Qv[i] = Q[i] + dec_u * u
    print(f'  axe exporte decale de {dec_u:+.2f} mm (ventral) : {source_dec}')
    # maillages : tube de l'axe + anneaux aux frontieres d'etages (si fiables)
    vol = np.zeros(shape, bool)
    R = int(round(0.10 / vox)); g = np.arange(-R, R + 1); ball = (g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2) <= R ** 2
    def stamp(c):
        i0 = np.round(c / vox).astype(int); lo = i0 - R; hi = i0 + R + 1
        blo = np.maximum(0, -lo); bhi = ball.shape[0] - np.maximum(0, hi - np.array(shape)); lo = np.maximum(lo, 0); hi = np.minimum(hi, shape)
        if np.all(hi > lo):
            vol[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]] |= ball[blo[0]:bhi[0], blo[1]:bhi[1], blo[2]:bhi[2]]
    for q in Qv[::2]:
        stamp(q)
    m_axe = mask_to_mesh(vol, vox, sigma=0.8, target_faces=30000, min_faces_component=20)
    # meme repere que les PLY du pipeline (meshexport.run) : centre de l'enveloppe (center_voxel) a l'origine
    c = np.array(man['center_voxel']); cc = np.array([c[0], c[2], -c[1]]) * vox
    m_axe.apply_translation(-cc)
    segs = []
    fa = f'{cs}_axe_vertebral.ply'; m_axe.export(os.path.join(out, fa))
    segs.append({'name': 'axe_vertebral', 'file': fa, 'color': [0.9, 0.1, 0.1], 'alpha': 1.0, 'faces': int(len(m_axe.faces)), 'volume_mm3': float(abs(m_axe.volume)),
                 'confiance': 'bonne' if cout < 1e4 else 'moyenne', 'note': f'axe {source_axe.split(" (")[0]}, decale de {dec_u:+.2f} mm sur les corps vertebraux ({source_dec}), arc {arc:.2f} mm'})
    if fiab['fiable']:
        vol[:] = False
        for nom, e in etages.items():
            if e.get('hors_axe') or not (0 <= e['s_debut_mm'] <= arc):
                continue
            i = int(np.argmin(np.abs(ss - e['s_debut_mm'])))
            u = np.array([0, 0, 1.0]); u -= (u @ T[i]) * T[i]; u /= np.linalg.norm(u); v = np.cross(T[i], u)
            r_ann = float(np.clip(0.03 * gl, 0.15, 0.5))          # rayon des anneaux a l'echelle du stade
            for a_ in np.linspace(0, 2 * np.pi, 24, endpoint=False):
                stamp(Qv[i] + r_ann * (np.cos(a_) * u + np.sin(a_) * v))
        m_et = mask_to_mesh(vol, vox, sigma=0.6, target_faces=60000, min_faces_component=20)
        m_et.apply_translation(-cc)
        fe = f'{cs}_etages.ply'; m_et.export(os.path.join(out, fe))
        segs.append({'name': 'etages_vertebraux', 'file': fe, 'color': [0.95, 0.75, 0.3], 'alpha': 1.0, 'faces': int(len(m_et.faces)), 'volume_mm3': float(abs(m_et.volume)),
                     'confiance': 'moyenne', 'note': f"anneaux aux frontieres des {fiab['n_etages']} niveaux imposes (occipital 4 -> Co-4), reperes membres + pointe de la queue"})
    # fichiers
    def ecrire_json(obj, chemin, **kw):            # ecriture atomique : les autres sessions lisent ces fichiers en continu
        tmp = chemin + '.tmp'
        with open(tmp, 'w', encoding='utf-8') as fh:
            json.dump(obj, fh, **kw)
        os.replace(tmp, chemin)
    par_s = [{'s_mm': round(float(ss[i]), 3), 'xyz_mm': [round(float(v), 3) for v in Q[i]], 'xyz_corps_mm': [round(float(v), 3) for v in Qv[i]],
              'tangente': [round(float(v), 4) for v in T[i]]} for i in range(0, len(ss), 5)]
    ecrire_json({'stage': cs, 'mm_per_voxel': vox, 'longueur_arc_mm': round(arc, 2), 'longueur_z_mm': round(z_ext, 2), 'ecart_lateral_mm': round(float(np.ptp(Q[:, 0])), 2),
               'source': source_axe + ', spline lissee', 'cout_chemin': float(cout),
               'repere': 'coordonnees volume (ax0 = LR, ax1 = cranio-caudal, ax2 = AP) en mm ; maillage : X = ax0, Y = ax2, Z = -ax1',
               'decalage_ventral_mm': round(dec_u, 3), 'decalage_source': source_dec,
               'note': 'xyz_mm = axe guide (moelle ou arcs) ; xyz_corps_mm = axe decale sur les corps vertebraux (celui des maillages)',
               'axe': par_s}, os.path.join(out, 'axe_vertebral.json'), indent=1)
    ecrire_json({'stage': cs, 'fiabilite': fiab, 'etages': etages}, os.path.join(out, 'etages_vertebraux.json'), indent=1)
    ecrire_json({'stage': cs, 'fiabilite': fiab_det, 'etages': etages_det, 'zones': {'methode': zones_det[1], 'bornes_s_mm': zones_det[0]}}, os.path.join(out, 'etages_detectes.json'), indent=1)
    ecrire_json({'stage': cs, 'note': 'axe prolonge aux deux bouts (s < 0 : prolongement cranial ; s > arc : caudal), pas irregulier',
                 'arc_mm': round(float(arc), 3), 'axe': [{'s_mm': round(float(ssx[i]), 3), 'xyz_mm': [round(float(v), 3) for v in Qx[i]]} for i in range(0, len(ssx), 3)]},
                os.path.join(out, 'axe_prolonge.json'), indent=1)
    ecrire_json({'stage': cs, 'mm_per_voxel': vox, 'repere': 'coordonnees volume (ax0 = LR, ax1 = cranio-caudal, ax2 = AP) en mm', **info_imp,
                 'niveaux': [{'nom': k, **e} for k, e in etages.items()]}, os.path.join(out, 'niveaux_imposes.json'), ensure_ascii=False, indent=1)
    # format lu par le recrutement de la session base (meme forme que etages_manuels.json, noms sans espace, memes 37 noms a tous les stades)
    os.makedirs(os.path.join(out, 'video360'), exist_ok=True)
    court = lambda k: k.replace('occipital ', 'occ').replace('-', '')
    ecrire_json({'stage': cs, 'source': 'niveaux imposes (topographie_axe.py)', 'n_etages': len(etages), 'noms': [court(k) for k in etages],
                 'note': 'xyz_mm_pipeline = centre du niveau sur l axe (repere pipeline, mm) ; hors_axe = centre extrapole au-dela de l axe (tete ou queue)',
                 'note_occipitaux': 'occ1..occ4 = somites occipitaux integres a la base du crane, sans corps vertebral (corps = false) : ne pas les recruter ni les morpher en vertebres',
                 'etages': [{'nom': court(k), 'n': i + 1, 'xyz_mm_pipeline': e['xyz_centre_mm'], 's_debut_mm': e['s_debut_mm'], 's_fin_mm': e['s_fin_mm'],
                             'region': 'occipital' if e['region'] == 'tete' else e['region'], 'corps': e['corps'], 'hors_axe': e['hors_axe'],
                             'sur_prolongement': e['sur_prolongement'], 'origine': e['origine']} for i, (k, e) in enumerate(etages.items())]},
                os.path.join(out, 'video360', 'etages_imposes.json'), ensure_ascii=False, indent=1)
    json.dump({'stage': cs, 'principe': 'zone = s du point de l axe le plus proche ; tete = au-dessus de l extremite craniale de l axe', 'methode': zones[1], 'bornes_s_mm': zones[0]},
              open(os.path.join(out, 'zones_axe.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    manifest = {'stage': cs, 'units': 'mm', 'mm_per_voxel': vox, 'center_voxel': man['center_voxel'], 'source': 'embryo3d/topographie_axe.py (lecture seule de labels.npz)',
                'segments': segs,
                'controle': {'axe_vertebral': {'longueur_arc_mm': round(arc, 2), 'longueur_z_mm': round(z_ext, 2), 'ecart_lateral_mm': round(float(np.ptp(Q[:, 0])), 2), 'source': fiab['guide']},
                             'etages': {'fiabilite': fiab, 'niveaux': etages}, 'etages_detectes': {'fiabilite': fiab_det, 'niveaux': etages_det},
                             'niveaux_imposes': info_imp, 'zones': {'methode': zones[1], 'bornes_s_mm': zones[0]},
                             'images': [f for f in ('etages_vertebraux.png', 'etages_detectes.png', 'controle_axe_3d.png', 'controle_extremites.png', 'controle_recrutement.png', 'controle_tube_neural.png')
                                        if os.path.exists(os.path.join(out, f))],
                             'avertissement': None if fiab['fiable'] else f"etages non fiables ({fiab['n_etages']} etages, periode {fiab['periode_mm']} mm, cv {fiab['cv_hauteurs']}) : zones approximatives"}}
    # conservation des segments et blocs de controle ajoutes par video360.py / marques_video.py (somites_video, etages_manuels, ...)
    fm = os.path.join(out, 'manifest_topographie.json')
    if os.path.exists(fm):
        try:
            ancien = json.load(open(fm, encoding='utf-8'))
            noms = {sg['name'] for sg in manifest['segments']}
            manifest['segments'] += [sg for sg in ancien.get('segments', []) if sg['name'] not in noms]
            for k, v in ancien.get('controle', {}).items():
                if k not in manifest['controle']:
                    manifest['controle'][k] = v
                elif k == 'images':
                    manifest['controle']['images'] = list(dict.fromkeys(manifest['controle']['images'] + [im for im in v if '/' in im]))
        except (ValueError, OSError) as e:
            print('  ancien manifest illisible, segments externes non conserves :', e)
    ecrire_json(manifest, fm, ensure_ascii=False, indent=1)
    print('->', out)
    return manifest


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('dossier'); ap.add_argument('--demi-mm', type=float, default=1.5); ap.add_argument('--lisser', type=float, default=15.0)
    ap.add_argument('--label-tube', default=None, help='forcer le label du guide tube neural (moelle_rachidienne, moelle, snc)')
    ap.add_argument('--guide', default='auto', choices=['auto', 'colonne', 'snc', 'niveaux'], help='auto : tube neural si les labels colonne couvrent < 60 %% du snc')
    a = ap.parse_args(); construire(a.dossier, a.demi_mm, a.lisser, a.guide, label_tube=a.label_tube)
