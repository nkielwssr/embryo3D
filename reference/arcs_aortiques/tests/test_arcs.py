"""Tests de l'extraction des modèles 3D de référence (U3D, PRC, PDF, identification, export).
Lancer : python -m pytest embryo3d/reference/arcs_aortiques/tests   (ou python .../tests/test_arcs.py)
Les fichiers de donnees/ ont été produits par les écrivains de référence (voir generer_donnees.py)."""
import json
import os
import sys
import tempfile

import numpy as np

ICI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(ICI))
sys.path.insert(0, ICI)

import fabrique  # noqa: E402
from extraire import candidats_stade, charger_legende, facteur_mm, identifier, traiter  # noqa: E402
from pdf3d import extraire_pdf  # noqa: E402
from prc import lire_prc  # noqa: E402
from u3d import lire_u3d  # noqa: E402

DON = os.path.join(ICI, 'donnees')
ATT = json.load(open(os.path.join(DON, 'attendu.json')))


def _lire(nom):
    return open(os.path.join(DON, nom), 'rb').read()


def _aire_volume(V, F):
    a, b, c = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    cr = np.cross(b - a, c - a)
    return np.linalg.norm(cr, axis=1).sum() / 2, np.einsum('ij,ij->i', a, np.cross(b, c)).sum() / 6


def _ferme_et_oriente(F):
    """Chaque arête orientée apparaît une fois et son opposée une fois : surface fermée, orientation cohérente.
    Une désynchronisation du flux casserait immédiatement cette propriété."""
    aretes = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    cles = set(map(tuple, aretes.tolist()))
    return len(cles) == len(aretes) and all((b, a) in cles for a, b in cles)


def _verifier(objets, attendu, tol=1e-3, tol_aire=None, tol_pos=None):
    par_nom = {o['noeud']: o for o in objets}
    assert set(attendu) <= set(par_nom), (set(attendu), set(par_nom))
    for nom, e in attendu.items():
        o = par_nom[nom]
        assert len(o['sommets']) == e['sommets'] and len(o['faces']) == e['faces'], nom
        assert _ferme_et_oriente(o['faces']), nom
        aire, vol = _aire_volume(o['sommets'], o['faces'])
        ta = tol if tol_aire is None else tol_aire
        assert abs(aire - e['aire']) <= ta * max(1, e['aire']), (nom, aire, e['aire'])
        assert abs(vol - e['volume']) <= ta * max(1, abs(e['volume'])), (nom, vol, e['volume'])   # signe : orientation
        tp = tol * 10 if tol_pos is None else tol_pos
        assert np.allclose(o['sommets'].min(0), e['min'], atol=tp) and np.allclose(o['sommets'].max(0), e['max'], atol=tp), nom
        if 'couleur' in e:
            assert np.allclose(o['couleur'][:3], e['couleur'], atol=1e-3), (nom, o['couleur'])
        if 'opacite' in e:
            assert abs(o['opacite'] - e['opacite']) < 0.01, (nom, o['opacite'])
        if 'groupes' in e:
            assert o['groupes'] == e['groupes'], (nom, o['groupes'])


def test_u3d_progressif_groupes_materiaux():
    sc = lire_u3d(_lire('progressif.u3d'))
    assert not sc.erreurs, sc.erreurs
    _verifier(sc.objets(), ATT['progressif.u3d'])


def test_u3d_couleurs_par_sommet():
    o = lire_u3d(_lire('couleurs.u3d')).objets()[0]
    ref = np.loadtxt(os.path.join(DON, 'couleurs_attendues.txt'))
    d = np.linalg.norm(o['sommets'][:, None, :] - ref[None, :, :3], axis=2)
    i = d.argmin(1)
    assert d.min(1).max() < 1e-4
    assert np.abs(o['couleurs_sommets'][:, :3] - ref[i, 3:]).max() < 1e-3


def test_u3d_deux_materiaux():
    sc = lire_u3d(_lire('deux_materiaux.u3d'))
    m = next(iter(sc.maillages.values()))
    assert np.bincount(m.f_mat[:m.nb_faces]).tolist() == ATT['deux_materiaux.u3d']['repartition']


def test_u3d_plus_de_16383_sommets():
    sc = lire_u3d(_lire('grand.u3d'))
    assert not sc.erreurs, sc.erreurs
    pas = next(iter(sc.maillages.values())).iq_pos           # qualité 100 : quantification grossière (pas ~0,018)
    _verifier(sc.objets(), ATT['grand.u3d'], tol_aire=0.05, tol_pos=pas)


def test_u3d_maillage_de_base_non_compresse():
    s = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    F = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]
    M = [[0, -1, 0, 1], [1, 0, 0, 2], [0, 0, 1, 3], [0, 0, 0, 1]]
    vc = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0]]
    buf = fabrique.u3d_base([dict(nom='tetra', sommets=s.tolist(), faces=F, couleur=(0.9, 0.5, 0.1), parent='G', M=M,
                                  couleurs_sommets=vc)], groupes=[dict(nom='G')])
    o = lire_u3d(buf).objets()[0]
    attendu = s @ np.array(M)[:3, :3].T + np.array(M)[:3, 3]
    assert np.allclose(o['sommets'], attendu, atol=1e-6)
    assert o['faces'].tolist() == F
    assert np.allclose(o['couleur'], (0.9, 0.5, 0.1), atol=1e-6)
    assert np.allclose(o['couleurs_sommets'][:, :3], vc, atol=1e-6)
    assert o['groupes'] == ['G']


def test_u3d_jeu_de_lignes():
    o = lire_u3d(_lire('lignes.u3d')).objets()[0]
    tt = np.linspace(0, 4 * np.pi, 60)
    P = np.stack([np.cos(tt) + 2, np.sin(tt), tt / 5], 1)            # hélice translatée de +2 en x par son nœud
    d = np.linalg.norm(o['sommets'][:, None] - P[None], axis=2)
    i = d.argmin(1)
    assert d.min(1).max() < 1e-3
    seg = {tuple(sorted((int(i[a]), int(i[b])))) for a, b in o['segments']}
    assert seg == {(k, k + 1) for k in range(59)} | {(0, 30)}
    assert o['noeud'] == 'Nervus_vagus' and len(o['faces']) == 0


def test_prc_arbre_tessellation():
    sc = lire_prc(_lire('arbre.prc'))
    assert not sc.erreurs, sc.erreurs
    _verifier(sc.objets(), ATT['arbre.prc'])


def test_pdf_u3d_et_prc():
    with tempfile.TemporaryDirectory() as d:
        p1 = os.path.join(d, 'a.pdf')
        open(p1, 'wb').write(fabrique.pdf_3d(_lire('progressif.u3d'), b'U3D', vues=[('Ventral', {'PAA3_left': True, 'boite': False})],
                                             titre='Carnegie stage 14'))
        p2 = os.path.join(d, 'b.pdf')
        open(p2, 'wb').write(fabrique.pdf_3d(_lire('arbre.prc'), b'PRC', compresser=False))
        r1, r2 = extraire_pdf(p1), extraire_pdf(p2)
        assert [f['format'] for f in r1['flux3d']] == ['U3D'] and r1['flux3d'][0]['donnees'] == _lire('progressif.u3d')
        v = r1['flux3d'][0]['vues'][0]
        assert v['nom'] == 'Ventral' and [(n['nom'], n['visible']) for n in v['noeuds']] == [('PAA3_left', True), ('boite', False)]
        assert [f['format'] for f in r2['flux3d']] == ['PRC'] and r2['flux3d'][0]['donnees'] == _lire('arbre.prc')


def test_pdf_balayage_sans_pypdf():
    import pdf3d
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, 'a.pdf')
        open(p, 'wb').write(fabrique.pdf_3d(_lire('couleurs.u3d'), b'U3D'))
        res = {'flux3d': []}
        pdf3d._par_balayage(open(p, 'rb').read(), res)
        assert len(res['flux3d']) == 1 and res['flux3d'][0]['donnees'] == _lire('couleurs.u3d')


def test_identification_par_nom_et_couleur():
    L = charger_legende()
    cas = {'PAA3 left': ('arc_3', 'gauche'), 'right 4th pharyngeal arch artery': ('arc_4', 'droite'),
           'Arch artery VI L': ('arc_6', 'gauche'), 'linker zesde boogslagader': ('arc_6', 'gauche'),
           '2e arc aortique droit': ('arc_2', 'droite'), 'Aortic_sac': ('sac_aortique', None),
           'Dorsal aorta R': ('aorte_dorsale', 'droite'), 'OFT lumen': ('voie_de_sortie', None),
           'Pharynx': ('tube_digestif', None), 'Neural tube': ('tube_neural', None), 'otocyst_left': ('otocyste', 'gauche'),
           'N. IX': ('nerf_IX', None), 'vagus nerve': ('nerf_X', None), 'Carotid duct left': ('canal_carotidien', 'gauche'),
           'Ductus arteriosus': ('canal_arteriel', None), 'Right subclavian artery': ('subclaviere', 'droite')}
    for nom, (structure, cote) in cas.items():
        r = identifier(nom, None, L)
        assert (r['structure'], r['cote']) == (structure, cote), (nom, r)
    r = identifier('Material #27', (0.21, 0.73, 0.24), L)
    assert r['structure'] == 'arc_3' and r['methode'] == 'couleur'
    assert identifier('Material #28', (0.0, 0.0, 0.0), L)['structure'] is None
    r = identifier('Material #28', None, L, corr={'Material #28': {'structure': 'arc_4', 'cote': 'droite'}})
    assert (r['structure'], r['cote'], r['methode']) == ('arc_4', 'droite', 'correspondance')


def test_stades_et_unites():
    assert candidats_stade('Rana_suppl_stage17.pdf') == ['CS17']
    assert candidats_stade('Carnegie stage 13', 'JH_CS15.pdf') == ['CS13', 'CS15']
    assert candidats_stade('focus 14 CSV') == []
    V = np.array([[0, 0, 0], [2400, 900, 1500]], float)
    assert facteur_mm([{'sommets': V}])[0] == 1e-3
    assert facteur_mm([{'sommets': V / 1000}])[0] == 1.0


def test_extraction_complete():
    L = charger_legende()
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, 'Supplementary_CS14.pdf')
        open(p, 'wb').write(fabrique.pdf_3d(_lire('progressif.u3d'), b'U3D'))
        idx = traiter([p], os.path.join(d, 'sortie'), L, avec_apercu=False)
        fl = idx['fichiers'][0]['flux'][0]
        assert fl['stade'] == 'CS14' and fl['dossier'] == 'CS14'
        man = json.load(open(os.path.join(d, 'sortie', 'CS14', 'manifest.json'), encoding='utf-8'))
        assert man['stage'] == 'CS14' and man['units'] == 'mm'
        noms = {s['name']: s for s in man['structures']}
        assert {'arc_3_gauche', 'sac_aortique'} <= set(noms)
        for s in man['structures']:
            assert {'name', 'file', 'collection', 'color', 'alpha'} <= set(s)          # format de blender_build_scene.py
            assert os.path.exists(os.path.join(d, 'sortie', 'CS14', s['file']))
        assert noms['sac_aortique']['alpha'] == 0.5
        assert os.path.exists(os.path.join(d, 'sortie', 'rapport.md'))


if __name__ == '__main__':
    n = 0
    for nom, f in sorted(globals().items()):
        if nom.startswith('test_') and callable(f):
            f()
            n += 1
            print('ok', nom)
    print(n, 'tests réussis')
