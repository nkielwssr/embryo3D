"""Tubes à topologie commune (pharynx + tube digestif + aortes) et nappes de mésos pour le morphing CS13 → CS20.

Chaque structure est reconstruite, à chaque stade, comme un tube de N anneaux × M sommets le long de sa ligne
centrale (work/digestif/chemins.json, work/cardio/vaisseaux_chemins.json ; pour une poche « sac » comme l'estomac,
ligne centrale et rayon équivalent tirés du masque, coupe par coupe). Coordonnées : mm, repère Blender du pipeline
(mêmes mm par voxel et même centre que meshexport.py). Un segment absent à un stade est réduit à un point posé sur
son raccord (fin ou début du segment voisin du même stade) : il « pousse » pendant le morphing.
Le pharynx (intestin pharyngien) est un segment digestif comme les autres : il est lu dans chemins.json s'il a été tracé
(segment « pharynx » de embryo3d/digestif_points/<CS>.json), sinon réduit au début de l'œsophage.

Mésos dorsaux (méso-œsophage, mésogastre dorsal, mésoduodénum, mésentère, mésocôlon dorsal) : pour chaque segment
digestif, nappe de N lignes × K colonnes tendue entre le bord dorsal de chaque anneau et le bord ventral de l'axe
aortique du stade (milieu des aortes dorsales paires quand elles se font face, puis aorte commune). L'attache glisse de
façon monotone le long de l'aorte (projection au plus proche puis régression isotonique sur toute la chaîne digestive :
pas de croisement, l'anse de l'intestin moyen donne un éventail depuis sa racine). Les lignes dont l'anneau est trop
loin de l'aorte (> --portee mm) ou dont le pied tombe au-delà d'une extrémité tracée sont réduites à leur bord digestif
(largeur nulle) ; un méso sans aorte au stade reste plaqué sur le tube (largeur nulle) et s'ouvre pendant le morphing.
Le pharynx n'a pas de méso (il est plaqué sous la notochorde, entre les aortes).

Arcs aortiques (artères des arcs pharyngiens 1, 2, 3, 4, 6, gauche et droite), sac aortique, tronc artériel et poches pharyngiennes 1-4
(dynamique CS11-CS16 de Rana et al. 2014) : tubes comme les autres, lus dans vaisseaux_chemins.json (arcs, sac, tronc) et chemins.json
(poches) dès qu'ils sont tracés ; un arc absent est réduit sur le sac aortique (sinon l'arc voisin, sinon l'aorte dorsale), une poche absente
sur le pharynx à sa hauteur attendue. vaisseaux_points/calendrier_arcs.json donne la présence attendue par stade ; les écarts sont signalés.

Sortie : embryons_3D/tubes_morph.npz (clé '<structure>' -> tableau (7, N, M, 3) pour un tube, (7, N, K, 3) pour un méso)
+ tubes_morph.json (stades, N, M, K, présence par stade, type/famille/forme par structure, attaches des mésos).
Lecteurs : tubes_morph_blender.py (scène Blender), tubes_morph_planche.py (contrôle 2D), viewer_3dh.py (site), blender_build_scene.py.
usage : python embryo3d/tubes_morph.py [--portee 3.0] [--colonnes 4] [--sans-mesos] [--sortie embryons_3D]"""
import numpy as np, json, os, argparse
from scipy import ndimage as ndi
from scipy.interpolate import interp1d

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STADES = [('CS13', 'CS13.f4v'), ('CS14', 'CS14_f4v'), ('CS15', 'CS15_f4v'), ('CS16', 'CS16_f4v'),
          ('CS17', 'CS17_f4v'), ('CS19', 'CS19_f4v'), ('CS20', 'CS20_F4V')]
N, M = 64, 16
# structure -> (fichier de chemins, clés acceptées par ordre de préférence, raccords si absent : liste de (structure, 'debut'|'fin'|fraction 0-1)
#               essayés dans l'ordre puis, s'ils manquent aussi, leurs propres raccords ; une fraction = abscisse relative le long de la structure)
STRUCT = {
    'pharynx':              ('digestif', ['pharynx', 'intestin_pharyngien'], [('oesophage', 'debut')]),
    'oesophage':            ('digestif', ['oesophage'], [('pharynx', 'fin')]),
    'estomac':              ('digestif', ['estomac'], [('oesophage', 'fin')]),
    'duodenum':             ('digestif', ['duodenum'], [('estomac', 'fin')]),
    'intestin_moyen':       ('digestif', ['intestin_moyen'], [('duodenum', 'fin')]),
    'intestin_posterieur':  ('digestif', ['intestin_posterieur'], [('intestin_moyen', 'fin')]),
    'aorte_dorsale_gauche': ('cardio', ['aorte_dorsale_gauche', 'aorte_dorsale'], [('aorte_commune', 'debut')]),
    'aorte_dorsale_droite': ('cardio', ['aorte_dorsale_droite'], [('aorte_dorsale_gauche', 'debut')]),
    'aorte_commune':        ('cardio', ['aorte_commune'], [('aorte_dorsale_gauche', 'fin')]),
    # voie de sortie du cœur (lumière du tronc artériel) et sac aortique, tracés du cœur (début) vers les arcs (fin)
    'tronc_arteriel':       ('cardio', ['tronc_arteriel', 'voie_de_sortie'], [('sac_aortique', 'debut'), ('aorte_dorsale_gauche', 'debut')]),
    'sac_aortique':         ('cardio', ['sac_aortique'], [('tronc_arteriel', 'fin')] + [(f'arc_aortique_{k}_gauche', 'debut') for k in (4, 3, 6, 2, 1)]
                                                          + [('aorte_dorsale_gauche', 'debut')]),
}
# artères des arcs pharyngiens 1 (mandibulaire), 2 (hyoïdien), 3 (carotidien), 4 (aortique), 6 (pulmonaire) — dynamique CS11-CS16 de
# Rana et al. 2014 —, tracées du sac aortique (ventral, début) vers l'aorte dorsale (dorsal, fin). Un arc absent (pas encore formé ou
# régressé) est réduit sur la fin du sac aortique, sinon sur le début de l'arc voisin le plus proche du même côté, sinon sur l'aorte dorsale :
# il pousse ou se résorbe depuis son origine ventrale pendant la transition.
ARCS, COTES = (1, 2, 3, 4, 6), ('gauche', 'droite')
for k in ARCS:
    for c in COTES:
        voisins = sorted((j for j in ARCS if j != k), key=lambda j: (abs(j - k), j))
        STRUCT[f'arc_aortique_{k}_{c}'] = ('cardio', [f'arc_aortique_{k}_{c}', f'arc_{k}_{c}'],
                                           [('sac_aortique', 'fin')] + [(f'arc_aortique_{j}_{c}', 'debut') for j in voisins]
                                           + [(f'aorte_dorsale_{c}', 'debut'), ('aorte_dorsale_gauche', 'debut')])
# poches pharyngiennes 1 à 4 (diverticules latéraux du pharynx entre les arcs), tracées de la lumière du pharynx (début) vers leur fond (fin) ;
# une poche absente est réduite sur le pharynx à la hauteur attendue (1/5, 2/5, 3/5, 4/5 de sa longueur), sinon au début de l'œsophage
for k in (1, 2, 3, 4):
    for c in COTES:
        STRUCT[f'poche_pharyngienne_{k}_{c}'] = ('digestif', [f'poche_pharyngienne_{k}_{c}', f'poche_{k}_{c}'], [('pharynx', k / 5), ('oesophage', 'debut')])

def famille(nom):
    if nom.startswith(('arc_aortique', 'sac_aortique', 'tronc_arteriel')): return 'arc'
    if nom.startswith('poche_pharyngienne'): return 'poche'
    if nom.startswith('aorte'): return 'aorte'
    if nom.startswith(('meso', 'mesentere')): return 'meso'
    return 'digestif'
FAMILLE = {n: famille(n) for n in STRUCT}
# couleurs RGB 0-1, partagées par la scène Blender, la planche et le site (écrites dans tubes_morph.json) ; arcs, sac, voie de sortie et
# aorte dorsale reprennent la légende de Rana et al. 2014 (mandibulaire beige, hyoïdien jaune, carotidien vert, aortique cyan, pulmonaire
# magenta, sac orange, voie de sortie bleu-violet, aorte rouge) ; le préfixe le plus long l'emporte
COULEURS = {'pharynx': (0.80, 0.25, 0.25), 'oesophage': (0.85, 0.40, 0.70), 'estomac': (0.95, 0.60, 0.20), 'duodenum': (0.30, 0.75, 0.35),
            'intestin_moyen': (0.95, 0.85, 0.25), 'intestin_posterieur': (0.30, 0.55, 0.95),
            'aorte_dorsale_gauche': (0.90, 0.10, 0.10), 'aorte_dorsale_droite': (0.95, 0.35, 0.20), 'aorte_commune': (0.75, 0.05, 0.10),
            'sac_aortique': (0.95, 0.55, 0.15), 'tronc_arteriel': (0.40, 0.30, 0.85),
            'arc_aortique_1': (0.93, 0.87, 0.70), 'arc_aortique_2': (0.95, 0.80, 0.20), 'arc_aortique_3': (0.20, 0.70, 0.30),
            'arc_aortique_4': (0.20, 0.80, 0.85), 'arc_aortique_6': (0.80, 0.20, 0.80), 'poche_pharyngienne': (0.78, 0.78, 0.86),
            'meso_oesophage': (0.96, 0.80, 0.72), 'mesogastre_dorsal': (0.96, 0.76, 0.66), 'mesoduodenum': (0.94, 0.78, 0.70),
            'mesentere': (0.97, 0.82, 0.74), 'mesocolon_dorsal': (0.93, 0.74, 0.68)}
def couleur(nom):
    for k in sorted(COULEURS, key=len, reverse=True):
        if nom == k or nom.startswith(k): return COULEURS[k]
    return (0.7, 0.7, 0.7)
# segments digestifs dans l'ordre crânio-caudal (chaîne pour la régression isotonique des attaches)
CHAINE = ['oesophage', 'estomac', 'duodenum', 'intestin_moyen', 'intestin_posterieur']
# segment digestif -> méso dorsal qui le relie à l'aorte
MESOS = {'oesophage': 'meso_oesophage', 'estomac': 'mesogastre_dorsal', 'duodenum': 'mesoduodenum',
         'intestin_moyen': 'mesentere', 'intestin_posterieur': 'mesocolon_dorsal'}
CRANIO_CAUDAL = ('pharynx', 'oesophage', 'estomac', 'aorte')     # segments dont le sens est forcé de haut en bas
DORSAL = np.array([0.0, 1.0, 0.0])                              # repère Blender du pipeline : +Y = dorsal

def to_mm(P, man):
    m = man['mm_per_voxel']; c = np.array(man['center_voxel'])
    P = np.asarray(P, float)
    return np.stack([(P[:, 0] - c[0]) * m, (P[:, 2] - c[2]) * m, -(P[:, 1] - c[1]) * m], 1)

def sac_ligne(work, nom):
    z = np.load(os.path.join(work, 'digestif', 'digestif.npz')); SH = tuple(z['shape'])
    m = np.unpackbits(z[nom])[:int(np.prod(SH))].reshape(SH).astype(bool)
    C, R = [], []
    for s in np.where(m.any(axis=(0, 2)))[0]:
        sl = m[:, s, :]; a = sl.sum()
        if a < 10: continue
        cx, cy = np.argwhere(sl).mean(0); C.append([cx, s, cy]); R.append(np.sqrt(a / np.pi))
    C, R = np.array(C), np.array(R)
    k = 5; C = ndi.uniform_filter1d(C, k, axis=0, mode='nearest'); R = ndi.uniform_filter1d(R, k, mode='nearest')
    return C, R

def reechantillonne(C, R, n=N):
    d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(C, axis=0), axis=1))]
    keep = np.r_[True, np.diff(d) > 1e-6]; C, R, d = C[keep], R[keep], d[keep]
    if len(C) < 2: return np.repeat(C[:1], n, 0), np.repeat(R[:1], n)
    u = np.linspace(0, d[-1], n)
    return interp1d(d, C, axis=0)(u), interp1d(d, R)(u)

def anneaux(C, R, ref=np.array([1.0, 0, 0])):
    """anneaux par transport parallèle, angle de départ aligné sur la direction de référence (moins de torsion entre stades)"""
    T = np.gradient(C, axis=0); T /= np.linalg.norm(T, axis=1, keepdims=True) + 1e-12
    V = np.zeros((len(C), M, 3)); ang = np.linspace(0, 2 * np.pi, M, endpoint=False)
    n1 = ref - T[0] * T[0].dot(ref)
    if np.linalg.norm(n1) < 1e-6: n1 = np.cross(T[0], [0, 0, 1.0])
    n1 /= np.linalg.norm(n1)
    for i in range(len(C)):
        if i > 0:
            n1 = n1 - T[i] * T[i].dot(n1); nn = np.linalg.norm(n1)
            n1 = n1 / nn if nn > 1e-9 else np.cross(T[i], [0, 0, 1.0])
        n2 = np.cross(T[i], n1)
        V[i] = C[i] + R[i] * (np.cos(ang)[:, None] * n1 + np.sin(ang)[:, None] * n2)
    return V

# ---------------------------------------------------------------- mésos : axe aortique, projection monotone, nappe
def axe_aortique(lg, st, seuil_paire=1.0):
    """ligne médiane des aortes d'un stade (centres mm, rayons mm) : moyenne des aortes dorsales paires là où elles se font
    face (point le plus proche à < seuil_paire mm ; la gauche fait foi ailleurs), puis aorte commune ; morceaux enchaînés
    du plus crânial au plus caudal (Z moyen décroissant). None si aucune aorte tracée au stade."""
    G, D, Cm = (lg.get((st, k)) for k in ('aorte_dorsale_gauche', 'aorte_dorsale_droite', 'aorte_commune'))
    morceaux = []
    if G is not None and D is not None:
        (CG, RG), (CD, RD) = G, D
        d = np.linalg.norm(CG[:, None] - CD[None], axis=2); j = d.argmin(1)
        proche = d[np.arange(len(CG)), j] < seuil_paire
        morceaux.append((np.where(proche[:, None], (CG + CD[j]) / 2, CG), np.where(proche, (RG + RD[j]) / 2, RG)))
    elif G is not None: morceaux.append(G)
    elif D is not None: morceaux.append(D)
    if Cm is not None: morceaux.append(Cm)
    if not morceaux: return None
    morceaux.sort(key=lambda cr: -cr[0][:, 2].mean())
    C = np.concatenate([c for c, _ in morceaux]); R = np.concatenate([r for _, r in morceaux])
    keep = np.r_[True, np.linalg.norm(np.diff(C, axis=0), axis=1) > 1e-6]
    return C[keep], R[keep]

def projeter(A, P):
    """projection de points P (n,3) sur la polyligne A (m,3) : abscisse curviligne (mm) du point le plus proche"""
    S, V = A[:-1], np.diff(A, axis=0); L = np.linalg.norm(V, axis=1) + 1e-12
    t = np.clip(((P[:, None, :] - S[None]) * V[None]).sum(2) / L[None] ** 2, 0, 1)      # n x (m-1)
    Q = S[None] + t[..., None] * V[None]
    k = np.linalg.norm(P[:, None, :] - Q, axis=2).argmin(1); n = np.arange(len(P))
    return np.r_[0, np.cumsum(L)][k] + t[n, k] * L[k]

def isotone(y):
    """régression isotonique (non décroissante, moindres carrés) par « pool adjacent violators »"""
    blocs = []
    for v in y:
        blocs.append([float(v), 1])
        while len(blocs) > 1 and blocs[-2][0] / blocs[-2][1] > blocs[-1][0] / blocs[-1][1]:
            s, n = blocs.pop(); blocs[-1][0] += s; blocs[-1][1] += n
    return np.concatenate([[s / n] * n for s, n in blocs])

def attaches(axe, C, portee, tol=0.3):
    """pour des centres d'anneaux C (n,3) en ordre crânio-caudal : point d'attache sur l'axe aortique (n,3), rayon aortique (n),
    abscisse (n) et masque des lignes actives (pied sur la ligne tracée, distance ≤ portee)"""
    A, RA = axe
    cum = np.r_[0, np.cumsum(np.linalg.norm(np.diff(A, axis=0), axis=1))]
    if len(A) < 2 or cum[-1] < 1e-6:
        Pa = np.repeat(A[:1], len(C), 0); return Pa, np.repeat(RA[:1], len(C)), np.zeros(len(C)), np.linalg.norm(C - Pa, axis=1) <= portee
    u = np.clip(isotone(projeter(A, C)), 0, cum[-1])
    Pa = interp1d(cum, A, axis=0)(u); Ra = interp1d(cum, RA)(u)
    T0 = A[1] - A[0]; T0 /= np.linalg.norm(T0); T1 = A[-1] - A[-2]; T1 /= np.linalg.norm(T1)
    dist = np.linalg.norm(C - Pa, axis=1)
    avant = (u <= 1e-6) & ((C - A[0]) @ T0 < -tol)          # pied avant le début tracé de l'aorte
    apres = (u >= cum[-1] - 1e-6) & ((C - A[-1]) @ T1 > tol)  # pied après sa fin tracée
    return Pa, Ra, u, (dist <= portee) & ~avant & ~apres

def nappe(C, R, Pa, Ra, actif, K):
    """nappe (n, K, 3) du bord dorsal des anneaux (centres C, rayons R) au bord ventral de l'aorte (attaches Pa, rayons Ra) ;
    lignes inactives ou aortes accolées au tube : largeur nulle sur le bord du tube"""
    n = Pa - C; dist = np.linalg.norm(n, axis=1)
    n = np.where(dist[:, None] > 1e-9, n / (dist[:, None] + 1e-12), DORSAL)
    t0 = R.copy(); t1 = np.where(actif, dist - Ra, t0); t1 = np.where(t1 > t0, t1, t0)
    lam = np.linspace(0, 1, K)
    return C[:, None, :] + n[:, None, :] * (t0[:, None] + (t1 - t0)[:, None] * lam[None, :])[:, :, None]

# ---------------------------------------------------------------- lecture des lignes centrales
def charger(stades=STADES):
    """(stade, structure) -> (centres mm (N,3), rayons mm (N,)) pour chaque segment tracé"""
    lignes = {}
    for st, dossier in stades:
        work, out = os.path.join(ROOT, dossier, 'work'), os.path.join(ROOT, dossier, 'out')
        pm = os.path.join(out, 'manifest.json')
        if not os.path.exists(pm): print(st, ': manifest absent, stade ignoré'); continue
        man = json.load(open(pm, encoding='utf-8'))
        ch = {}
        for f, key in (('digestif', os.path.join(work, 'digestif', 'chemins.json')), ('cardio', os.path.join(work, 'cardio', 'vaisseaux_chemins.json'))):
            ch[f] = json.load(open(key, encoding='utf-8')) if os.path.exists(key) else {}
        for nom, (f, cles, _) in STRUCT.items():
            k = next((c for c in cles if c in ch[f]), None)
            if k is None: continue
            e = ch[f][k]
            if e.get('type') == 'sac': C, R = sac_ligne(work, k)
            else: C, R = np.array(e['centres'], float), np.array(e['rayons'], float)
            if len(C) < 2: continue
            if C[0, 1] > C[-1, 1] and nom.startswith(CRANIO_CAUDAL): C, R = C[::-1], R[::-1]   # sens crânio-caudal (axe 1 : 0 = haut)
            Cm = to_mm(C, man); Rm = R * man['mm_per_voxel']
            lignes[(st, nom)] = reechantillonne(Cm, Rm)
    return lignes

def point_raccord(lignes, st, nom):
    """point (mm) où une structure absente est réduite : premier raccord tracé au stade (début, fin ou fraction de la longueur d'une
    structure voisine), en élargissant aux raccords des voisines absentes ; sinon centre des structures tracées"""
    file, vus = list(STRUCT[nom][2]), {nom}
    while file:
        r, ou = file.pop(0)
        if (st, r) in lignes:
            C, _ = lignes[(st, r)]
            return C[-1] if ou == 'fin' else C[0] if ou == 'debut' else C[int(round(float(ou) * (len(C) - 1)))]
        if r not in vus: vus.add(r); file += list(STRUCT[r][2])
    pres = [lignes[(st, n)][0].mean(0) for n in STRUCT if (st, n) in lignes]
    return np.mean(pres, axis=0) if pres else np.zeros(3)

HERE = os.path.dirname(os.path.abspath(__file__))
def verifier_calendrier(presence, chemin=os.path.join(HERE, 'vaisseaux_points', 'calendrier_arcs.json')):
    """compare la présence tracée des arcs, du sac, du tronc et des poches au calendrier attendu (littérature) :
    {'a_tracer': {stade: [structures attendues mais absentes]}, 'inattendus': {stade: [tracées mais attendues absentes]}}"""
    if not os.path.exists(chemin): return None
    cal = json.load(open(chemin, encoding='utf-8')).get('stades', {}); ecarts = {'a_tracer': {}, 'inattendus': {}}
    for nom, pres in presence.items():
        base = nom.rsplit('_', 1)[0] if nom.endswith(('_gauche', '_droite')) else nom
        for st, ok in pres.items():
            attendu = cal.get(st, {}).get(nom, cal.get(st, {}).get(base))     # clé latéralisée prioritaire (ex. arc_aortique_6_droite)
            if attendu is None: continue
            if ok and attendu == 'absent': ecarts['inattendus'].setdefault(st, []).append(nom)
            if not ok and attendu in ('present', 'formation', 'regression'): ecarts['a_tracer'].setdefault(st, []).append(f'{nom} ({attendu})')
    return ecarts

def construire(lignes, stades=STADES, portee=3.0, K=4, mesos=True):
    """tableaux (n_stades, N, M|K, 3) par structure + métadonnées (présence, type/famille/forme, attaches des mésos)"""
    noms_st = [s for s, _ in stades]
    res, presence, structures = {}, {}, {}
    for nom in STRUCT:
        arr = np.zeros((len(stades), N, M, 3)); pres = []
        for i, st in enumerate(noms_st):
            if (st, nom) in lignes:
                C, R = lignes[(st, nom)]; arr[i] = anneaux(C, R); pres.append(True)
            else:
                arr[i] = np.broadcast_to(point_raccord(lignes, st, nom), (N, M, 3)); pres.append(False)
        res[nom] = arr; presence[nom] = dict(zip(noms_st, pres))
        structures[nom] = {'type': 'tube', 'famille': FAMILLE[nom], 'forme': [N, M], 'couleur': list(couleur(nom))}
    attach = {}
    if mesos:
        # attaches calculées sur toute la chaîne digestive du stade (monotones le long de l'aorte), puis découpées par segment
        par_stade = {}
        for i, st in enumerate(noms_st):
            axe = axe_aortique(lignes, st); segs = [s for s in CHAINE if (st, s) in lignes]
            if axe is None or not segs: par_stade[st] = None; continue
            Ctot = np.concatenate([lignes[(st, s)][0] for s in segs])
            Pa, Ra, u, actif = attaches(axe, Ctot, portee)
            par_stade[st] = {s: (Pa[k * N:(k + 1) * N], Ra[k * N:(k + 1) * N], u[k * N:(k + 1) * N], actif[k * N:(k + 1) * N]) for k, s in enumerate(segs)}
        for seg, nom in MESOS.items():
            arr = np.zeros((len(stades), N, K, 3)); pres = []; attach[nom] = {}
            for i, st in enumerate(noms_st):
                if (st, seg) not in lignes:                       # segment absent : même point que le tube réduit
                    arr[i] = np.broadcast_to(res[seg][i][0, 0], (N, K, 3)); pres.append(False); continue
                C, R = lignes[(st, seg)]
                a = (par_stade.get(st) or {}).get(seg)
                if a is None:                                     # pas d'aorte au stade : plaqué sur le bord dorsal du tube
                    arr[i] = nappe(C, R, C + DORSAL, np.zeros(N), np.zeros(N, bool), K); pres.append(False); continue
                Pa, Ra, u, actif = a
                arr[i] = nappe(C, R, Pa, Ra, actif, K); pres.append(bool(actif.any()))
                larg = np.linalg.norm(arr[i][:, -1] - arr[i][:, 0], axis=1)
                attach[nom][st] = {'lignes_actives': int(actif.sum()), 'abscisse_aorte_mm': [round(float(u[actif].min()), 2), round(float(u[actif].max()), 2)] if actif.any() else None,
                                   'largeur_mediane_mm': round(float(np.median(larg[actif])), 3) if actif.any() else 0.0}
            res[nom] = arr; presence[nom] = dict(zip(noms_st, pres))
            structures[nom] = {'type': 'nappe', 'famille': 'meso', 'forme': [N, K], 'couleur': list(couleur(nom)), 'segment': seg,
                               'attache': 'axe aortique (aortes dorsales / aorte commune)'}
    meta = {'stades': noms_st, 'N': N, 'M': M, 'K': K, 'presence': presence, 'structures': structures,
            'mesos': {'portee_max_mm': portee, 'colonnes': K, 'attaches': attach} if mesos else None,
            'calendrier': verifier_calendrier(presence)}
    return res, meta

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--portee', type=float, default=3.0, help='distance maximale anneau digestif ↔ aorte pour tendre le méso (mm)')
    ap.add_argument('--colonnes', type=int, default=4, help='colonnes de la nappe (bord digestif → bord aortique)')
    ap.add_argument('--sans-mesos', action='store_true', help='tubes seulement (pas de nappes)')
    ap.add_argument('--sortie', default=os.path.join(ROOT, 'embryons_3D'), help='dossier de tubes_morph.npz / .json')
    args = ap.parse_args()
    lignes = charger()
    res, meta = construire(lignes, portee=args.portee, K=max(2, args.colonnes), mesos=not args.sans_mesos)
    for nom in res:
        pres = meta['presence'][nom]; typ = meta['structures'][nom]['type']
        print(f'{nom:22s} {typ:5s}', ' '.join(('X' if pres[s] else '.') for s in meta['stades']),
              '' if typ == 'tube' else '  ' + ' '.join(f"{s}:{a['lignes_actives']}" for s, a in meta['mesos']['attaches'][nom].items()))
    if meta['calendrier']:
        for cle, titre in (('a_tracer', 'attendus (calendrier) mais non tracés'), ('inattendus', 'tracés mais attendus absents')):
            for st in meta['stades']:
                if st in meta['calendrier'][cle]: print(f'  {st} {titre} : ' + ', '.join(meta['calendrier'][cle][st]))
    os.makedirs(args.sortie, exist_ok=True)
    np.savez_compressed(os.path.join(args.sortie, 'tubes_morph.npz'), **res)
    json.dump(meta, open(os.path.join(args.sortie, 'tubes_morph.json'), 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    print('OK', os.path.join(args.sortie, 'tubes_morph.npz'))
