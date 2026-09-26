# -*- coding: utf-8 -*-
"""Test sans données de pharynx_candidats.py : volume synthétique de type CS13 (dens.npy uint8, tissu dense, lumières pâles bruitées)
avec un pharynx aplati (large en x, mince en y) qui monte depuis le début de l'œsophage puis se courbe vers la bouche (ventrale),
4 poches à gauche et 3 à droite, et des leurres : bourgeon pulmonaire relié à la jonction pharynx-œsophage (caudal), ventricule
cérébral, cœur, aortes dorsales (étiquetés). L'outil doit tracer le pharynx du fond de la cavité buccale au début de l'œsophage, trouver
chaque poche à sa hauteur et du bon côté, les numéroter d'après le calendrier CS13, écrire la planche, puis « fusionner » doit recopier
pharynx et poches dans un fichier de points de stade.

    python embryo3d/test_pharynx_candidats_synthetique.py [dossier_de_travail]     (défaut : dossier temporaire, effacé à la fin)
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

ICI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ICI)
from test_arcs_candidats_synthetique import peindre, ellipsoide, SH, XM  # noqa: E402

E0 = (XM, 210, 150)                                          # début de l'œsophage
AXE = [(XM, 75, 100), (XM, 100, 130), (XM, 150, 145), E0]    # axe du pharynx, du fond de la cavité buccale à l'œsophage
POCHES_G = {1: 110, 2: 135, 3: 160, 4: 182}                  # hauteur s de chaque poche
POCHES_D = {1: 110, 2: 135, 3: 160}


def generer(work):
    os.makedirs(os.path.join(work, 'cardio'), exist_ok=True); os.makedirs(os.path.join(work, 'digestif'), exist_ok=True)
    rs = np.random.RandomState(5)
    env = ellipsoide((XM, 170, 125), (105, 160, 110))
    lum = np.zeros(SH, bool)
    for dx in range(-14, 15, 2): peindre(lum, [(x + dx, s, y) for x, s, y in AXE], 5)     # pharynx aplati
    for c, poches in (('gauche', POCHES_G), ('droite', POCHES_D)):
        sx = -1 if c == 'gauche' else 1
        for k, s in poches.items(): peindre(lum, [(XM + sx * 15, s, 142), (XM + sx * 45, s, 140)], 4)
    peindre(lum, [E0, (XM, 290, 160)], 4)                                                 # œsophage
    peindre(lum, [(XM, 214, 140), (XM, 285, 128)], 3)                                     # bourgeon pulmonaire (caudal)
    vent = ellipsoide((XM, 55, 185), (25, 20, 18)) & env                                  # ventricule cérébral
    coeur = ellipsoide((XM, 200, 60), (30, 24, 24))
    ao = np.zeros(SH, bool)
    for x in (XM - 30, XM + 30): peindre(ao, [(x, 60, 175), (x, 300, 175)], 4)
    d = np.where(env, 110.0, 0.0) + rs.normal(0, 8, SH)
    for m in (lum, vent, coeur, ao): d[m] = 14 + rs.normal(0, 6, int(m.sum()))
    np.save(os.path.join(work, 'dens.npy'), np.clip(d, 0, 255).astype(np.uint8))
    pk = np.packbits
    np.savez_compressed(os.path.join(work, 'labels.npz'), enveloppe=pk(env), coeur=pk(coeur), ventricules=pk(vent), shape=np.array(SH))
    np.savez_compressed(os.path.join(work, 'cardio', 'vaisseaux.npz'), aorte_dorsale_gauche=pk(ao & (np.arange(SH[0]) < XM)[:, None, None]),
                        aorte_dorsale_droite=pk(ao & (np.arange(SH[0]) >= XM)[:, None, None]), shape=np.array(SH))
    json.dump({f'aorte_dorsale_{c}': {'centres': [[x, s, 175] for s in range(60, 301)], 'rayons': [4.0] * 241} for c, x in (('gauche', XM - 30), ('droite', XM + 30))},
              open(os.path.join(work, 'cardio', 'vaisseaux_chemins.json'), 'w', encoding='utf-8'))
    # œsophage tracé du bas vers le haut (l'outil doit remettre le début en haut)
    json.dump({'oesophage': {'centres': [[XM, s, 150 + (s - 210) / 8] for s in range(290, 209, -1)], 'rayons': [4.0] * 81}},
              open(os.path.join(work, 'digestif', 'chemins.json'), 'w', encoding='utf-8'))


def dist_axe(P):
    """distance de chaque point à la polyligne AXE"""
    A = np.array(AXE, float); out = []
    for p in np.asarray(P, float):
        best = np.inf
        for a, b in zip(A[:-1], A[1:]):
            v = b - a; t = np.clip((p - a) @ v / (v @ v), 0, 1); best = min(best, np.linalg.norm(p - (a + t * v)))
        out.append(best)
    return np.array(out)


def main():
    racine = sys.argv[1] if len(sys.argv) > 1 else tempfile.mkdtemp(prefix='pharynx_candidats_')
    garder = len(sys.argv) > 1
    try:
        work = os.path.join(racine, 'CS13.f4v', 'work'); generer(work)
        sortie, png = os.path.join(racine, 'CS13_pharynx_proposes.json'), os.path.join(racine, 'pharynx_candidats.png')
        r = subprocess.run([sys.executable, os.path.join(ICI, 'pharynx_candidats.py'), work, 'CS13', '--sortie', sortie, '--planche', png],
                           capture_output=True, text=True)
        print(r.stdout.rstrip())
        if r.returncode != 0: print(r.stderr.rstrip()); sys.exit('ÉCHEC : pharynx_candidats.py a échoué')
        res = json.load(open(sortie, encoding='utf-8')); segs = {s['nom']: s for s in res['segments']}
        erreurs = []
        def ok(cond, msg):
            if not cond: erreurs.append(msg)
        ok('pharynx' in segs, 'pas de pharynx proposé')
        if 'pharynx' in segs:
            P = np.array(segs['pharynx']['points'], float)
            ok(np.linalg.norm(P[0] - np.array(AXE[0])) <= 10, f'pharynx : début {P[0].tolist()} loin du fond de la cavité buccale {AXE[0]}')
            ok(np.linalg.norm(P[-1] - np.array(E0)) <= 2, f'pharynx : fin {P[-1].tolist()} ≠ début de l\'œsophage {E0}')
            ok(dist_axe(P).max() <= 6, f'pharynx : écart à l\'axe {dist_axe(P).max():.1f} vox (lumière : rayon 5 autour de l\'axe)')
        attendu = {f'poche_pharyngienne_{k}_gauche': s for k, s in POCHES_G.items()}
        attendu.update({f'poche_pharyngienne_{k}_droite': s for k, s in POCHES_D.items()})
        poches = {n for n in segs if n.startswith('poche')}
        ok(poches == set(attendu), f'poches trouvées {sorted(poches)} ≠ attendues {sorted(attendu)}')
        for n, s in attendu.items():
            if n not in segs: continue
            P = np.array(segs[n]['points']); fond = P[-1]
            ok(abs(fond[1] - s) <= 4, f'{n} : fond à s={fond[1]}, attendu {s}')
            ok((fond[0] < XM - 30) if n.endswith('gauche') else (fond[0] > XM + 30), f'{n} : fond {fond.tolist()} pas assez latéral ou du mauvais côté')
            ok(dist_axe(P[:1]).max() <= 16, f'{n} : ne part pas de la lumière du pharynx ({P[0].tolist()})')
        ok(os.path.exists(png), 'planche absente')
        cible = os.path.join(racine, 'CS13.json')
        json.dump({'stade': 'CS13', 'segments': [{'nom': 'oesophage', 'points': [[120, 210, 150], [120, 290, 160]]}]}, open(cible, 'w', encoding='utf-8'))
        r = subprocess.run([sys.executable, os.path.join(ICI, 'pharynx_candidats.py'), 'fusionner', sortie, cible], capture_output=True, text=True)
        print(r.stdout.rstrip())
        noms = [s['nom'] for s in json.load(open(cible, encoding='utf-8'))['segments']]
        ok(r.returncode == 0 and noms[0] == 'oesophage' and 'pharynx' in noms and set(attendu) <= set(noms), f'fusion : {noms}')
        for e in erreurs: print('ÉCHEC :', e)
        if erreurs: sys.exit(1)
        print('OK : test synthétique pharynx_candidats réussi')
    finally:
        if not garder: shutil.rmtree(racine, ignore_errors=True)


if __name__ == '__main__':
    main()
