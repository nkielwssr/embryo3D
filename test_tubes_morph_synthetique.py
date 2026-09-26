# -*- coding: utf-8 -*-
"""Test sans données de tubes_morph.py (+ tubes_morph_planche.py) : deux stades synthétiques au format du pipeline
(out/manifest.json, work/digestif/chemins.json, work/cardio/vaisseaux_chemins.json) avec pharynx, tube digestif, aortes,
arcs aortiques, sac aortique, tronc artériel et poches pharyngiennes, puis vérification de la présence, des raccords des
structures absentes, du calendrier, des nappes de mésos, de l'écriture npz/json et de la planche 2D.

    python embryo3d/test_tubes_morph_synthetique.py [dossier_de_travail]     (défaut : dossier temporaire, effacé à la fin)
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
import tubes_morph as tm  # noqa: E402

MM = 0.01                       # mm par voxel
CENTRE = [200, 200, 200]        # voxel (x gauche-droite, s haut-bas, y ventral-dorsal)


def ligne(*pts, n=12, r=4.0):
    """chemin {centres, rayons} en voxels par segments droits entre les points donnés"""
    P = np.array(pts, float); C = [P[0]]
    for a, b in zip(P[:-1], P[1:]):
        C += [a + (b - a) * t for t in np.linspace(0, 1, n)[1:]]
    return {'centres': np.round(C, 2).tolist(), 'rayons': [r] * len(C)}


def arc(s, cote, y_sac=150, y_ao=260):
    """arc aortique au niveau s : du sac aortique (ventral, médian) à l'aorte dorsale du côté donné, en contournant le pharynx"""
    x = 170 if cote == 'gauche' else 230; xl = 160 if cote == 'gauche' else 240
    return ligne([200, s, y_sac], [xl, s, 205], [x, s, y_ao], r=2.5)


def poche(s, cote):
    return ligne([200, s, 220], [165 if cote == 'gauche' else 235, s, 215], r=3.0)


def stade(racine, dossier, digestif, cardio):
    work, out = os.path.join(racine, dossier, 'work'), os.path.join(racine, dossier, 'out')
    for d in (os.path.join(work, 'digestif'), os.path.join(work, 'cardio'), out): os.makedirs(d, exist_ok=True)
    json.dump({'stage': dossier, 'units': 'mm', 'mm_per_voxel': MM, 'center_voxel': CENTRE, 'structures': []},
              open(os.path.join(out, 'manifest.json'), 'w', encoding='utf-8'))
    json.dump(digestif, open(os.path.join(work, 'digestif', 'chemins.json'), 'w', encoding='utf-8'))
    json.dump(cardio, open(os.path.join(work, 'cardio', 'vaisseaux_chemins.json'), 'w', encoding='utf-8'))


def generer(racine):
    base_dig = {'oesophage': ligne([200, 180, 225], [200, 260, 228]), 'estomac': ligne([200, 260, 228], [205, 300, 225], r=8),
                'duodenum': ligne([205, 300, 225], [205, 320, 222]), 'intestin_moyen': ligne([205, 320, 222], [205, 350, 180], [205, 370, 222]),
                'intestin_posterieur': ligne([205, 370, 222], [200, 390, 235])}
    # aortes dorsales tracées du bas vers le haut (le sens crânio-caudal doit être rétabli), fusion en aorte commune à s = 300
    base_ao = {'aorte_dorsale_gauche': ligne([185, 300, 255], [170, 100, 260]), 'aorte_dorsale_droite': ligne([215, 300, 255], [230, 100, 260]),
               'aorte_commune': ligne([200, 300, 255], [200, 380, 252], r=5)}
    # CS13 : pharynx, arcs 2-4 des deux côtés, arc 6 gauche seul (formation), arc 1 absent (régression), poches 1-3 + 4 gauche
    dig13 = dict(base_dig, pharynx=ligne([200, 100, 220], [200, 180, 225], r=6))
    for k, s in ((1, 118), (2, 136), (3, 154)):
        for c in ('gauche', 'droite'): dig13[f'poche_pharyngienne_{k}_{c}'] = poche(s, c)
    dig13['poche_pharyngienne_4_gauche'] = poche(170, 'gauche')
    ao13 = dict(base_ao, sac_aortique=ligne([200, 200, 150], [200, 170, 150], r=5), tronc_arteriel=ligne([200, 240, 140], [200, 200, 150], r=6))
    for k, s in ((2, 128), (3, 146), (4, 162)):
        for c in ('gauche', 'droite'): ao13[f'arc_aortique_{k}_{c}'] = arc(s, c)
    ao13['arc_aortique_6_gauche'] = arc(176, 'gauche')
    stade(racine, 'CS13.f4v', dig13, ao13)
    # CS14 : pas de pharynx ni de poches tracés, arcs 3, 4, 6 des deux côtés ; arc 1 tracé alors qu'il est attendu absent
    ao14 = dict(base_ao, sac_aortique=ligne([200, 200, 150], [200, 170, 150], r=5))
    for k, s in ((3, 146), (4, 162), (6, 176), (1, 112)):
        for c in ('gauche', 'droite'): ao14[f'arc_aortique_{k}_{c}'] = arc(s, c)
    stade(racine, 'CS14_f4v', base_dig, ao14)
    # CS15 : manifest seul (rien de tracé) ; CS16… : absents (stades ignorés)
    stade(racine, 'CS15_f4v', {}, {})


def verifier(racine):
    tm.ROOT = racine
    lignes = tm.charger()
    res, meta = tm.construire(lignes)
    ST, P = meta['stades'], meta['presence']
    i13, i14 = ST.index('CS13'), ST.index('CS14')
    erreurs = []
    def ok(cond, msg):
        if not cond: erreurs.append(msg)
    # présence lue
    for n in ('pharynx', 'sac_aortique', 'tronc_arteriel', 'arc_aortique_2_gauche', 'arc_aortique_4_droite', 'arc_aortique_6_gauche', 'poche_pharyngienne_4_gauche'):
        ok(P[n]['CS13'], f'{n} devrait être présent à CS13')
    for n in ('arc_aortique_1_gauche', 'arc_aortique_6_droite', 'poche_pharyngienne_4_droite'):
        ok(not P[n]['CS13'], f'{n} devrait être absent à CS13')
    ok(P['arc_aortique_6_droite']['CS14'] and not P['tronc_arteriel']['CS14'], 'présences CS14')
    ok(not any(P[n]['CS15'] for n in tm.STRUCT), 'CS15 : rien ne doit être présent')
    # sens crânio-caudal rétabli pour les aortes (début = haut = Z maximal)
    C = lignes[('CS13', 'aorte_dorsale_gauche')][0]; ok(C[0, 2] > C[-1, 2], 'aorte dorsale gauche : sens crânio-caudal non rétabli')
    # formes
    for n, a in res.items():
        ok(a.shape[0] == len(ST) and a.shape[1] == tm.N and np.isfinite(a).all(), f'{n} : forme ou valeurs invalides {a.shape}')
    # raccords des absents : arc 1 (CS13) -> fin du sac aortique ; arc 6 droit (CS13) -> fin du sac ; poche 4 droite -> pharynx à 4/5
    fin_sac = lignes[('CS13', 'sac_aortique')][0][-1]
    for n in ('arc_aortique_1_gauche', 'arc_aortique_6_droite'):
        a = res[n][i13]; ok(np.ptp(a.reshape(-1, 3), 0).max() < 1e-9 and np.allclose(a[0, 0], fin_sac), f'{n} CS13 : non réduit sur la fin du sac aortique')
    Cph = lignes[('CS13', 'pharynx')][0]
    ok(np.allclose(res['poche_pharyngienne_4_droite'][i13][0, 0], Cph[int(round(0.8 * (len(Cph) - 1)))]), 'poche 4 droite CS13 : raccord pharynx 4/5 faux')
    # CS14 : pas de pharynx -> poches réduites au début de l'œsophage (raccord du pharynx absent)
    deb_oe = lignes[('CS14', 'oesophage')][0][0]
    ok(np.allclose(res['poche_pharyngienne_1_gauche'][i14][0, 0], deb_oe), 'poche 1 gauche CS14 : devrait être réduite au début de l\'œsophage')
    ok(np.allclose(res['pharynx'][i14][0, 0], deb_oe), 'pharynx CS14 : devrait être réduit au début de l\'œsophage')
    # CS14 : tronc artériel absent -> début du sac aortique
    ok(np.allclose(res['tronc_arteriel'][i14][0, 0], lignes[('CS14', 'sac_aortique')][0][0]), 'tronc artériel CS14 : raccord faux')
    # arcs tracés : du sac (ventral, Y petit) vers l'aorte dorsale (dorsal, Y grand)
    A = res['arc_aortique_3_gauche'][i13].mean(1); ok(A[0, 1] < A[-1, 1], 'arc 3 gauche : sens sac -> aorte dorsale inversé')
    # calendrier
    cal = meta['calendrier'] or {}
    a13 = ' '.join(cal.get('a_tracer', {}).get('CS13', []))
    ok('arc_aortique_1_gauche (regression)' in a13 and 'arc_aortique_6_droite (formation)' in a13 and 'poche_pharyngienne_4_droite (formation)' in a13,
       f'calendrier CS13 : écarts attendus manquants ({a13})')
    ok('arc_aortique_2_gauche' not in a13, 'calendrier CS13 : arc 2 gauche tracé signalé à tort')
    ok(set(cal.get('inattendus', {}).get('CS14', [])) >= {'arc_aortique_1_gauche', 'arc_aortique_1_droite'}, 'calendrier CS14 : arc 1 tracé non signalé')
    # mésos : tendus à CS13 (aortes proches du tube), largeur nulle à CS15 (rien de tracé)
    at = meta['mesos']['attaches']
    ok(at['meso_oesophage'].get('CS13', {}).get('lignes_actives', 0) > 0, f"méso-œsophage CS13 non tendu : {at['meso_oesophage'].get('CS13')}")
    ok(res['meso_oesophage'].shape[2] == meta['K'], 'méso : nombre de colonnes')
    for n in tm.MESOS.values(): ok(not P[n]['CS15'], f'{n} ne doit pas être présent à CS15')
    # familles / couleurs
    S = meta['structures']
    ok(S['arc_aortique_4_gauche']['famille'] == 'arc' and S['poche_pharyngienne_2_droite']['famille'] == 'poche' and S['sac_aortique']['famille'] == 'arc', 'familles')
    ok(tuple(S['arc_aortique_6_droite']['couleur']) == tm.COULEURS['arc_aortique_6'], 'couleur arc 6 (préfixe)')
    return res, meta, erreurs


def main():
    racine = sys.argv[1] if len(sys.argv) > 1 else tempfile.mkdtemp(prefix='tubes_morph_')
    garder = len(sys.argv) > 1
    try:
        generer(racine)
        res, meta, erreurs = verifier(racine)
        sortie = os.path.join(racine, 'embryons_3D'); os.makedirs(sortie, exist_ok=True)
        np.savez_compressed(os.path.join(sortie, 'tubes_morph.npz'), **res)
        json.dump(meta, open(os.path.join(sortie, 'tubes_morph.json'), 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
        png = os.path.join(sortie, 'tubes_planche.png')
        r = subprocess.run([sys.executable, os.path.join(ICI, 'tubes_morph_planche.py'), png, '0', '0.5', '1', '2', '--px-mm', '60',
                            '--npz', os.path.join(sortie, 'tubes_morph.npz')], capture_output=True, text=True)
        if r.returncode != 0 or not os.path.exists(png): erreurs.append('planche : ' + (r.stderr.strip().splitlines() or ['échec'])[-1])
        for e in erreurs: print('ÉCHEC :', e)
        if erreurs: sys.exit(1)
        print('planche :', png if garder else '(dossier temporaire)')
        print('OK : test synthétique tubes_morph réussi')
    finally:
        if not garder: shutil.rmtree(racine, ignore_errors=True)


if __name__ == '__main__':
    main()
