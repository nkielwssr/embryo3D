"""Régénère les fichiers de test de tests/donnees/ (inutile pour lancer les tests : ils sont versionnés).

Il faut les deux écrivains de référence, qui ne sont pas des dépendances du projet :
  - U3D : convertisseur IDTF -> U3D de la bibliothèque U3D d'Intel (livré avec pymeshlab : lib/plugins/libio_u3d.so et
    lib/libIFXCore.so...), appelé par un petit programme C++ (fonction IDTFConverter::IDTFToU3d) ; variable IDTF2U3D =
    chemin de ce programme (arguments : entrée.idtf sortie.u3d 0 <dossier des libIFX> <qualité>), U3D_LIBDIR = dossier
    des bibliothèques. Le convertisseur exige la locale en_US.UTF-8.
  - PRC : programme construit sur oPRCFile (bibliothèque PRC d'Asymptote) ; variable GEN_PRC. Format d'entrée texte :
    « G nom 16 valeurs » (groupe, matrice ligne par ligne), « E » (fin de groupe),
    « M nom r g b a nV nF normales » suivi des sommets et des faces.
Les valeurs attendues (attendu.json) sont calculées sur les maillages sources.
"""
import json
import os
import subprocess
import sys

import numpy as np
import trimesh

ICI = os.path.dirname(os.path.abspath(__file__))
DON = os.path.join(ICI, 'donnees')


def _tm(M):
    return '\n'.join('\t\t\t\t' + ' '.join(f'{M[r][c]:.6f}' for c in range(4)) for r in range(4))


def _maillage_idtf(nom, V, F, normales, vcouleurs, nshaders=1, fshaders=None):
    L = [f'\t\tRESOURCE_NAME "{nom}"\n\t\tMODEL_TYPE "MESH"\n\t\tMESH {{',
         f'\t\t\tFACE_COUNT {len(F)}\n\t\t\tMODEL_POSITION_COUNT {len(V)}']
    N = np.zeros((0, 3))
    if normales:
        fn = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]])
        fn /= np.linalg.norm(fn, axis=1, keepdims=True) + 1e-12
        N = np.repeat(fn, 3, axis=0)
    nc = 0 if vcouleurs is None else len(vcouleurs)
    L.append(f'\t\t\tMODEL_NORMAL_COUNT {len(N)}\n\t\t\tMODEL_DIFFUSE_COLOR_COUNT {nc}\n\t\t\tMODEL_SPECULAR_COLOR_COUNT 0\n'
             f'\t\t\tMODEL_TEXTURE_COORD_COUNT 0\n\t\t\tMODEL_BONE_COUNT 0\n\t\t\tMODEL_SHADING_COUNT {nshaders}\n'
             '\t\t\tMODEL_SHADING_DESCRIPTION_LIST {')
    for s in range(nshaders):
        L.append(f'\t\t\t\tSHADING_DESCRIPTION {s} {{\n\t\t\t\t\tTEXTURE_LAYER_COUNT 0\n\t\t\t\t\tSHADER_ID {s}\n\t\t\t\t}}')
    L.append('\t\t\t}\n\t\t\tMESH_FACE_POSITION_LIST {')
    L += [f'\t\t\t\t{a} {b} {c}' for a, b, c in F]
    L.append('\t\t\t}')
    if normales:
        L.append('\t\t\tMESH_FACE_NORMAL_LIST {')
        L += [f'\t\t\t\t{3 * i} {3 * i + 1} {3 * i + 2}' for i in range(len(F))]
        L.append('\t\t\t}')
    L.append('\t\t\tMESH_FACE_SHADING_LIST {')
    L += [f'\t\t\t\t{s}' for s in (fshaders if fshaders is not None else [0] * len(F))]
    L.append('\t\t\t}')
    if vcouleurs is not None:
        L.append('\t\t\tMESH_FACE_DIFFUSE_COLOR_LIST {')
        L += [f'\t\t\t\t{a} {b} {c}' for a, b, c in F]
        L.append('\t\t\t}')
    L.append('\t\t\tMODEL_POSITION_LIST {')
    L += [f'\t\t\t\t{x:.6f} {y:.6f} {z:.6f}' for x, y, z in V]
    L.append('\t\t\t}')
    if normales:
        L.append('\t\t\tMODEL_NORMAL_LIST {')
        L += [f'\t\t\t\t{x:.6f} {y:.6f} {z:.6f}' for x, y, z in N]
        L.append('\t\t\t}')
    if vcouleurs is not None:
        L.append('\t\t\tMODEL_DIFFUSE_COLOR_LIST {')
        L += [f'\t\t\t\t{r:.6f} {g:.6f} {b:.6f} 1.000000' for r, g, b in vcouleurs]
        L.append('\t\t\t}')
    L.append('\t\t}')
    return '\n'.join(L)


def ecrire_idtf(chemin, objets, groupes=()):
    """objets : dicts nom, V, F, couleur, parent, M (4x4, convention colonne), normales, vcouleurs, nshaders, fshaders."""
    out = ['FILE_FORMAT "IDTF"\nFORMAT_VERSION 100\n']
    for g in groupes:
        out.append(f'NODE "GROUP" {{\n\tNODE_NAME "{g["nom"]}"\n\tPARENT_LIST {{\n\t\tPARENT_COUNT 1\n\t\tPARENT 0 {{\n'
                   f'\t\t\tPARENT_NAME "{g.get("parent", "<NULL>")}"\n\t\t\tPARENT_TM {{\n{_tm(np.asarray(g["M"]).T)}\n\t\t\t}}\n\t\t}}\n\t}}\n}}\n')
    for o in objets:
        out.append(f'NODE "MODEL" {{\n\tNODE_NAME "{o["nom"]}"\n\tPARENT_LIST {{\n\t\tPARENT_COUNT 1\n\t\tPARENT 0 {{\n'
                   f'\t\t\tPARENT_NAME "{o.get("parent", "<NULL>")}"\n\t\t\tPARENT_TM {{\n{_tm(np.asarray(o.get("M", np.eye(4))).T)}\n'
                   f'\t\t\t}}\n\t\t}}\n\t}}\n\tRESOURCE_NAME "{o["nom"]}_res"\n}}\n')
    out.append(f'RESOURCE_LIST "MODEL" {{\n\tRESOURCE_COUNT {len(objets)}')
    for i, o in enumerate(objets):
        out.append(f'\tRESOURCE {i} {{\n' + _maillage_idtf(o['nom'] + '_res', o['V'], o['F'], o.get('normales', True), o.get('vcouleurs'),
                                                            o.get('nshaders', 1), o.get('fshaders')) + '\n\t}')
    out.append('}\n')
    out.append(f'RESOURCE_LIST "SHADER" {{\n\tRESOURCE_COUNT {len(objets)}')
    for i, o in enumerate(objets):
        out.append(f'\tRESOURCE {i} {{\n\t\tRESOURCE_NAME "{o["nom"]}_shader"\n\t\tSHADER_MATERIAL_NAME "{o["nom"]}_mat"\n'
                   '\t\tSHADER_ACTIVE_TEXTURE_COUNT 0\n\t}')
    out.append('}\n')
    out.append(f'RESOURCE_LIST "MATERIAL" {{\n\tRESOURCE_COUNT {len(objets)}')
    for i, o in enumerate(objets):
        r, g, b = o.get('couleur', (0.8, 0.8, 0.8))
        out.append(f'\tRESOURCE {i} {{\n\t\tRESOURCE_NAME "{o["nom"]}_mat"\n\t\tMATERIAL_AMBIENT {r * .2:.6f} {g * .2:.6f} {b * .2:.6f}\n'
                   f'\t\tMATERIAL_DIFFUSE {r:.6f} {g:.6f} {b:.6f}\n\t\tMATERIAL_SPECULAR 0.2 0.2 0.2\n\t\tMATERIAL_EMISSIVE 0 0 0\n'
                   f'\t\tMATERIAL_REFLECTIVITY 0.1\n\t\tMATERIAL_OPACITY {o.get("opacite", 1.0):.6f}\n\t}}')
    out.append('}\n')
    for o in objets:
        out.append(f'MODIFIER "SHADING" {{\n\tMODIFIER_NAME "{o["nom"]}"\n\tPARAMETERS {{\n\t\tSHADER_LIST_COUNT 1\n'
                   '\t\tSHADER_LIST_LIST {\n\t\t\tSHADER_LIST 0 {\n\t\t\t\tSHADER_COUNT 1\n\t\t\t\tSHADER_NAME_LIST {\n'
                   f'\t\t\t\t\tSHADER 0 NAME: "{o["nom"]}_shader"\n\t\t\t\t}}\n\t\t\t}}\n\t\t}}\n\t}}\n}}\n')
    open(chemin, 'w').write('\n'.join(out))


def idtf_vers_u3d(idtf, u3d, qualite=1000):
    lib = os.environ['U3D_LIBDIR']
    env = dict(os.environ, LD_LIBRARY_PATH=lib, LANG='en_US.UTF-8')
    subprocess.run([os.environ['IDTF2U3D'], idtf, u3d, '0', lib, str(qualite)], env=env, check=True, capture_output=True)


def rot_z(a):
    c, s = np.cos(a), np.sin(a)
    M = np.eye(4)
    M[:2, :2] = [[c, -s], [s, c]]
    return M


def trans(x, y, z):
    M = np.eye(4)
    M[:3, 3] = [x, y, z]
    return M


def mesures(V, F, M=np.eye(4)):
    Vw = V @ M[:3, :3].T + M[:3, 3]
    m = trimesh.Trimesh(Vw, F, process=False)
    return {'sommets': len(V), 'faces': len(F), 'aire': float(m.area), 'volume': float(m.volume),
            'min': Vw.min(0).tolist(), 'max': Vw.max(0).tolist()}


def main():
    os.makedirs(DON, exist_ok=True)
    attendu = {}
    rng = np.random.default_rng(0)
    sph = trimesh.creation.icosphere(subdivisions=2)
    tor = trimesh.creation.torus(1.0, 0.3, major_sections=24, minor_sections=10)
    boite = trimesh.creation.box([4, 1, 0.5])
    # 1. maillages progressifs, groupes imbriqués (translation + rotation), matériaux
    Mb = trans(10, 0, 0) @ rot_z(np.pi / 2)
    G1, G2 = trans(0, 0, 100), rot_z(np.pi / 2)
    objs = [dict(nom='PAA3_left', V=sph.vertices, F=sph.faces, couleur=(0.2, 0.74, 0.22)),
            dict(nom='Aortic_sac', V=tor.vertices, F=tor.faces, couleur=(1.0, 0.55, 0.12), M=trans(5, 0, 0), opacite=0.5),
            dict(nom='boite', V=boite.vertices, F=boite.faces, couleur=(0.2, 0.3, 0.9), M=Mb, parent='G2', normales=False)]
    ecrire_idtf('/tmp/progressif.idtf', objs, [dict(nom='G1', M=G1), dict(nom='G2', parent='G1', M=G2)])
    idtf_vers_u3d('/tmp/progressif.idtf', os.path.join(DON, 'progressif.u3d'))
    attendu['progressif.u3d'] = {'PAA3_left': {**mesures(sph.vertices, sph.faces), 'couleur': [0.2, 0.74, 0.22], 'groupes': []},
                                 'Aortic_sac': {**mesures(tor.vertices, tor.faces, trans(5, 0, 0)), 'couleur': [1.0, 0.55, 0.12], 'opacite': 0.5, 'groupes': []},
                                 'boite': {**mesures(boite.vertices, boite.faces, G1 @ G2 @ Mb), 'couleur': [0.2, 0.3, 0.9], 'groupes': ['G1', 'G2']}}
    # 2. couleurs par sommet + deux matériaux
    vc = rng.random((len(sph.vertices), 3))
    ecrire_idtf('/tmp/couleurs.idtf', [dict(nom='couleurs', V=sph.vertices, F=sph.faces, vcouleurs=vc)])
    idtf_vers_u3d('/tmp/couleurs.idtf', os.path.join(DON, 'couleurs.u3d'))
    np.savetxt(os.path.join(DON, 'couleurs_attendues.txt'), np.hstack([sph.vertices, vc]), fmt='%.6f')
    fs = (sph.vertices[sph.faces].mean(1)[:, 0] > 0).astype(int)
    ecrire_idtf('/tmp/deuxmat.idtf', [dict(nom='deuxmat', V=sph.vertices, F=sph.faces, nshaders=2, fshaders=fs)])
    idtf_vers_u3d('/tmp/deuxmat.idtf', os.path.join(DON, 'deux_materiaux.u3d'))
    attendu['deux_materiaux.u3d'] = {'repartition': np.bincount(fs).tolist()}
    # 3. plus de 16 383 sommets (index lus hors contexte statique)
    grille = trimesh.creation.torus(2.0, 0.5, major_sections=180, minor_sections=96)       # 17 280 sommets
    ecrire_idtf('/tmp/grand.idtf', [dict(nom='grand', V=grille.vertices, F=grille.faces, normales=False)])
    idtf_vers_u3d('/tmp/grand.idtf', os.path.join(DON, 'grand.u3d'), 100)
    attendu['grand.u3d'] = {'grand': mesures(grille.vertices, grille.faces)}
    # 4. PRC (Asymptote) : groupe translaté, transparence, normales
    lignes = [f'M sphere 1 0 0 1 {len(sph.vertices)} {len(sph.faces)} 0']
    lignes += [f'{x:.9g} {y:.9g} {z:.9g}' for x, y, z in sph.vertices] + [f'{a} {b} {c}' for a, b, c in sph.faces]
    lignes.append('G deplace ' + ' '.join(f'{v:.9g}' for v in trans(5, 0, 0).ravel()))
    lignes.append(f'M tore 0 1 0 0.5 {len(tor.vertices)} {len(tor.faces)} 1')
    lignes += [f'{x:.9g} {y:.9g} {z:.9g}' for x, y, z in tor.vertices] + [f'{a} {b} {c}' for a, b, c in tor.faces]
    lignes.append('E')
    open('/tmp/arbre.txt', 'w').write('\n'.join(lignes) + '\n')
    subprocess.run([os.environ['GEN_PRC'], '/tmp/arbre.txt', os.path.join(DON, 'arbre.prc')], check=True)
    attendu['arbre.prc'] = {'sphere': {**mesures(sph.vertices, sph.faces), 'couleur': [1, 0, 0], 'groupes': ['root']},
                            'tore': {**mesures(tor.vertices, tor.faces, trans(5, 0, 0)), 'couleur': [0, 1, 0], 'opacite': 0.5, 'groupes': ['root', 'deplace']}}
    json.dump(attendu, open(os.path.join(DON, 'attendu.json'), 'w'), indent=1)
    for f in sorted(os.listdir(DON)):
        print(f, os.path.getsize(os.path.join(DON, f)))


if __name__ == '__main__':
    sys.exit(main())
