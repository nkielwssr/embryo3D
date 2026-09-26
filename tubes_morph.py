"""Tubes à topologie commune (tube digestif + aortes) pour le morphing CS13 → CS20.

Chaque structure est reconstruite, à chaque stade, comme un tube de N anneaux × M sommets le long de sa ligne
centrale (work/digestif/chemins.json, work/cardio/vaisseaux_chemins.json ; pour une poche « sac » comme l'estomac,
ligne centrale et rayon équivalent tirés du masque, coupe par coupe). Coordonnées : mm, repère Blender du pipeline
(mêmes mm par voxel et même centre que meshexport.py). Un segment absent à un stade est réduit à un point posé sur
son raccord (fin ou début du segment voisin du même stade) : il « pousse » pendant le morphing.
Sortie : embryons_3D/tubes_morph.npz  (clé '<structure>' -> tableau (7, N, M, 3)), + tubes_morph.json (stades, présence).
usage : python embryo3d/tubes_morph.py"""
import numpy as np, json, os
from scipy import ndimage as ndi
from scipy.interpolate import interp1d

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STADES = [('CS13', 'CS13.f4v'), ('CS14', 'CS14_f4v'), ('CS15', 'CS15_f4v'), ('CS16', 'CS16_f4v'),
          ('CS17', 'CS17_f4v'), ('CS19', 'CS19_f4v'), ('CS20', 'CS20_F4V')]
N, M = 64, 16
# structure -> (fichier de chemins, clés acceptées par ordre de préférence, raccord si absent : (structure, 'fin'|'debut'))
STRUCT = {
    'oesophage':            ('digestif', ['oesophage'], None),
    'estomac':              ('digestif', ['estomac'], ('oesophage', 'fin')),
    'duodenum':             ('digestif', ['duodenum'], ('estomac', 'fin')),
    'intestin_moyen':       ('digestif', ['intestin_moyen'], ('duodenum', 'fin')),
    'intestin_posterieur':  ('digestif', ['intestin_posterieur'], ('intestin_moyen', 'fin')),
    'aorte_dorsale_gauche': ('cardio', ['aorte_dorsale_gauche', 'aorte_dorsale'], ('aorte_commune', 'debut')),
    'aorte_dorsale_droite': ('cardio', ['aorte_dorsale_droite'], ('aorte_dorsale_gauche', 'debut')),
    'aorte_commune':        ('cardio', ['aorte_commune'], ('aorte_dorsale_gauche', 'fin')),
}

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

lignes = {}      # (stade, structure) -> (C mm, R mm)
for st, dossier in STADES:
    work, out = os.path.join(ROOT, dossier, 'work'), os.path.join(ROOT, dossier, 'out')
    man = json.load(open(os.path.join(out, 'manifest.json'), encoding='utf-8'))
    ch = {}
    for f, key in (('digestif', os.path.join(work, 'digestif', 'chemins.json')), ('cardio', os.path.join(work, 'cardio', 'vaisseaux_chemins.json'))):
        ch[f] = json.load(open(key, encoding='utf-8')) if os.path.exists(key) else {}
    for nom, (f, cles, _) in STRUCT.items():
        k = next((c for c in cles if c in ch[f]), None)
        if k is None: continue
        e = ch[f][k]
        if e.get('type') == 'sac': C, R = sac_ligne(work, k)
        else: C, R = np.array(e['centres']), np.array(e['rayons'])
        if len(C) < 2: continue
        if C[0, 1] > C[-1, 1] and nom.startswith(('oesophage', 'estomac', 'aorte')): C, R = C[::-1], R[::-1]   # sens crânio-caudal
        Cm = to_mm(C, man); Rm = R * man['mm_per_voxel']
        lignes[(st, nom)] = reechantillonne(Cm, Rm)

res, presence = {}, {}
for nom, (_, _, raccord) in STRUCT.items():
    arr = np.zeros((len(STADES), N, M, 3)); pres = []
    for i, (st, _) in enumerate(STADES):
        if (st, nom) in lignes:
            C, R = lignes[(st, nom)]; arr[i] = anneaux(C, R); pres.append(True); continue
        pres.append(False)
        p = None; r = raccord; vus = {nom}
        while r is not None and p is None and r[0] not in vus:   # remonte la chaîne des raccords
            vus.add(r[0])
            if (st, r[0]) in lignes:
                C, _ = lignes[(st, r[0])]; p = C[-1] if r[1] == 'fin' else C[0]
            else: r = STRUCT[r[0]][2]
        if p is None:
            p = np.mean([lignes[(st, n)][0].mean(0) for n in STRUCT if (st, n) in lignes], axis=0)
        arr[i] = np.broadcast_to(p, (N, M, 3))
    res[nom] = arr; presence[nom] = dict(zip([s for s, _ in STADES], pres))
    print(f'{nom:22s}', ' '.join(('X' if v else '.') for v in pres))
os.makedirs(os.path.join(ROOT, 'embryons_3D'), exist_ok=True)
np.savez_compressed(os.path.join(ROOT, 'embryons_3D', 'tubes_morph.npz'), **res)
json.dump({'stades': [s for s, _ in STADES], 'N': N, 'M': M, 'presence': presence},
          open(os.path.join(ROOT, 'embryons_3D', 'tubes_morph.json'), 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
print('OK')
