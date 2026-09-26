# -*- coding: utf-8 -*-
"""Test sans données de arcs_candidats.py : volume synthétique de type CS13 (dens.npy uint8, tissu dense, lumières pâles bruitées)
avec cœur + voie de sortie + sac aortique, aortes dorsales, arcs 1-2-3-4-6 à gauche et 3-4 à droite, et des leurres pâles
(pharynx, poches pharyngiennes, veines cardinales, cavité péricardique, extérieur de l'embryon, artère intersegmentaire dorsale,
kyste pâle collé à l'aorte). L'outil doit retrouver chaque arc à sa hauteur, du bon côté, le numéroter d'après le calendrier CS13
sans se laisser décaler par le kyste (candidat « douteux »), proposer un sac aortique, écrire la planche, puis « fusionner » doit
recopier les arcs validés dans un fichier de points de stade.

    python embryo3d/test_arcs_candidats_synthetique.py [dossier_de_travail]     (défaut : dossier temporaire, effacé à la fin)
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

ICI = os.path.dirname(os.path.abspath(__file__))
SH = (240, 320, 240)                                         # x gauche-droite, s haut-bas, y ventral-dorsal
XM = 120                                                     # ligne médiane
ARCS_G = {1: 70, 2: 90, 3: 110, 4: 128, 6: 146}              # hauteur s de la jonction avec l'aorte dorsale
ARCS_D = {3: 110, 4: 128}
SAC_FIN = (XM, 140, 85)


def peindre(m, pts, r):
    """ajoute au masque m un tube de rayon r (voxels) le long de la polyligne pts"""
    P = np.array(pts, float)
    for a, b in zip(P[:-1], P[1:]):
        lo = np.maximum(np.floor(np.minimum(a, b) - r - 1).astype(int), 0); hi = np.minimum(np.ceil(np.maximum(a, b) + r + 2).astype(int), m.shape)
        g = np.stack(np.meshgrid(*[np.arange(lo[k], hi[k]) for k in range(3)], indexing='ij'), -1).astype(float)
        v = b - a; t = np.clip(((g - a) @ v) / max(v @ v, 1e-9), 0, 1)
        m[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]] |= np.linalg.norm(g - (a + t[..., None] * v), axis=-1) <= r


def ellipsoide(c, r):
    g = np.ogrid[:SH[0], :SH[1], :SH[2]]
    return sum(((g[k] - c[k]) / r[k]) ** 2 for k in range(3)) <= 1


def arc(s, cote):
    sx = -1 if cote == 'gauche' else 1
    return [SAC_FIN, (XM + sx * 15, s, 95), (XM + sx * 42, s, 130), (XM + sx * 30, s, 170)]


def generer(work):
    os.makedirs(os.path.join(work, 'cardio'), exist_ok=True)
    rs = np.random.RandomState(3)
    env = ellipsoide((XM, 170, 125), (105, 160, 110))
    coeur = ellipsoide((XM, 205, 60), (32, 26, 26)); peindre(coeur, [(XM, 205, 60), (XM, 160, 75)], 7)
    cav = ellipsoide((XM, 205, 60), (24, 18, 18)); peindre(cav, [(XM, 200, 62), (XM, 150, 80)], 5)   # cavités + voie de sortie
    peri = ellipsoide((XM, 205, 60), (42, 34, 34)) & ~ellipsoide((XM, 205, 60), (36, 29, 29)) & env
    lum = cav.copy()
    peindre(lum, [(XM, 150, 80), SAC_FIN], 5)                                    # sac aortique
    ao = {c: np.zeros(SH, bool) for c in ('gauche', 'droite')}
    for c, x in (('gauche', XM - 30), ('droite', XM + 30)): peindre(ao[c], [(x, 40, 170), (x, 300, 170)], 4)
    for c, arcs in (('gauche', ARCS_G), ('droite', ARCS_D)):
        for k, s in arcs.items(): peindre(lum, arc(s, c), 2.5)
    leurres = np.zeros(SH, bool)
    peindre(leurres, [(XM, 50, 130), (XM, 200, 132)], 10)                        # pharynx
    for s in (80, 100, 119, 137):                                                # poches, entre les arcs
        for sx in (-1, 1): peindre(leurres, [(XM, s, 130), (XM + sx * 32, s, 130)], 3)
    peindre(leurres, [(XM - 30, 100, 170), (XM - 30, 100, 200)], 2)            # artère intersegmentaire (dorsale, reliée à l'aorte gauche)
    peindre(leurres, [(XM - 40, 98, 168), (XM - 40, 98, 168)], 3)              # kyste pâle latéral à l'aorte gauche, non relié au sac
    card = np.zeros(SH, bool)
    for x in (XM - 50, XM + 50): peindre(card, [(x, 60, 176), (x, 260, 176)], 3)
    d = np.where(env, 110.0, 0.0) + rs.normal(0, 8, SH)
    for m in (lum, ao['gauche'], ao['droite'], leurres, card, peri): d[m] = 14 + rs.normal(0, 6, int(m.sum()))
    np.save(os.path.join(work, 'dens.npy'), np.clip(d, 0, 255).astype(np.uint8))
    pk = lambda m: np.packbits(m)
    np.savez_compressed(os.path.join(work, 'labels.npz'), enveloppe=pk(env), coeur=pk(coeur), cavite_pericardique=pk(peri), shape=np.array(SH))
    np.savez_compressed(os.path.join(work, 'cardio', 'vaisseaux.npz'), aorte_dorsale_gauche=pk(ao['gauche']), aorte_dorsale_droite=pk(ao['droite']),
                        cardinale_anterieure_gauche=pk(card & (np.arange(SH[0]) < XM)[:, None, None]),
                        cardinale_anterieure_droite=pk(card & (np.arange(SH[0]) >= XM)[:, None, None]), shape=np.array(SH))
    # aortes tracées du bas vers le haut (l'outil doit remettre le début en haut)
    ch = {f'aorte_dorsale_{c}': {'centres': [[x, s, 170] for s in range(300, 39, -1)], 'rayons': [4.0] * 261}
          for c, x in (('gauche', XM - 30), ('droite', XM + 30))}
    json.dump(ch, open(os.path.join(work, 'cardio', 'vaisseaux_chemins.json'), 'w', encoding='utf-8'))


def main():
    racine = sys.argv[1] if len(sys.argv) > 1 else tempfile.mkdtemp(prefix='arcs_candidats_')
    garder = len(sys.argv) > 1
    try:
        work = os.path.join(racine, 'CS13.f4v', 'work'); generer(work)
        sortie, png = os.path.join(racine, 'CS13_arcs_proposes.json'), os.path.join(racine, 'arcs_candidats.png')
        r = subprocess.run([sys.executable, os.path.join(ICI, 'arcs_candidats.py'), work, 'CS13', '--sortie', sortie, '--planche', png],
                           capture_output=True, text=True)
        print(r.stdout.rstrip())
        if r.returncode != 0: print(r.stderr.rstrip()); sys.exit('ÉCHEC : arcs_candidats.py a échoué')
        res = json.load(open(sortie, encoding='utf-8')); segs = {s['nom']: s for s in res['segments']}
        erreurs = []
        def ok(cond, msg):
            if not cond: erreurs.append(msg)
        attendu = {f'arc_aortique_{k}_gauche': s for k, s in ARCS_G.items()}; attendu.update({f'arc_aortique_{k}_droite': s for k, s in ARCS_D.items()})
        arcs = {n for n in segs if n.startswith('arc_aortique_')}; autres = sorted(n for n in segs if n.startswith('arc_candidat'))
        ok(arcs == set(attendu), f'arcs trouvés {sorted(arcs)} ≠ attendus {sorted(attendu)}')
        ok(autres == ['arc_candidat_douteux_1_gauche'], f'candidats hors arcs : {autres} (attendu : le kyste seul, douteux)')
        for n, s in attendu.items():
            if n not in segs: continue
            q = segs[n]['qualite']; j = q['jonction_aortique']
            ok(abs(j[1] - s) <= 4, f'{n} : jonction s={j[1]}, attendue {s}')
            ok((j[0] < XM) == n.endswith('gauche'), f'{n} : jonction du mauvais côté (x={j[0]})')
            ok(q['fraction_pale'] >= 0.9 and q['traversee_dense_max_vox'] <= 2, f'{n} : chemin hors lumière {q}')
            P = np.array(segs[n]['points']); ok(P[0][2] < P[-1][2], f'{n} : sens sac → aorte dorsale inversé')
        ok('sac_aortique' in segs, 'pas de proposition de sac aortique')
        if 'sac_aortique' in segs:
            P = np.array(segs['sac_aortique']['points']); ok(np.linalg.norm(P[-1] - np.array(SAC_FIN)) <= 12, f'sac aortique : fin {P[-1].tolist()} loin de {SAC_FIN}')
        ok(os.path.exists(png), 'planche absente')
        # fusion dans un fichier de points de stade (copie) : arcs + sac ajoutés, aortes conservées, pas de bloc « qualite »
        cible = os.path.join(racine, 'CS13.json')
        json.dump({'stade': 'CS13', 'fichier': 'vaisseaux.npz', 'segments': [{'nom': 'aorte_dorsale_gauche', 'points': [[90, 40, 170], [90, 300, 170]]},
                                                                             {'nom': 'arc_aortique_3_gauche', 'points': [[0, 0, 0], [1, 1, 1]]}]},
                  open(cible, 'w', encoding='utf-8'))
        r = subprocess.run([sys.executable, os.path.join(ICI, 'arcs_candidats.py'), 'fusionner', sortie, cible], capture_output=True, text=True)
        print(r.stdout.rstrip())
        T = json.load(open(cible, encoding='utf-8')); noms = [s['nom'] for s in T['segments']]
        ok(r.returncode == 0 and 'aorte_dorsale_gauche' in noms and set(attendu) <= set(noms) and 'sac_aortique' in noms, f'fusion : {noms}')
        ok(noms.count('arc_aortique_3_gauche') == 1 and all('qualite' not in s for s in T['segments']), 'fusion : doublon ou bloc qualite restant')
        ok(next(s for s in T['segments'] if s['nom'] == 'arc_aortique_3_gauche')['points'] != [[0, 0, 0], [1, 1, 1]], 'fusion : arc existant non remplacé')
        for e in erreurs: print('ÉCHEC :', e)
        if erreurs: sys.exit(1)
        print('OK : test synthétique arcs_candidats réussi')
    finally:
        if not garder: shutil.rmtree(racine, ignore_errors=True)


if __name__ == '__main__':
    main()
