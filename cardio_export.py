"""Exporte le cœur (work/cardio/coeur.npz) en PLY alignés sur les maillages du pipeline : <out>/cardio/<CS>_coeur.ply,
<CS>_cavites_cardiaques.ply + manifest_cardio.json (la confiance vient de embryo3d/cardio_points/<CS>.json).
Le réseau vasculaire reste celui du pipeline (<out>/<CS>_vaisseaux.ply, détection automatique), référencé tel quel.
usage : python cardio_export.py <dossier_stade> [...]"""
import numpy as np, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from meshexport import mask_to_mesh
HERE = os.path.dirname(os.path.abspath(__file__))
for d in sys.argv[1:]:
    work, out = os.path.join(d, 'work'), os.path.join(d, 'out')
    man = json.load(open(os.path.join(out, 'manifest.json'), encoding='utf-8'))
    st, s_mm, c = man['stage'], man['mm_per_voxel'], np.array(man['center_voxel'])
    pts = json.load(open(os.path.join(HERE, 'cardio_points', f'{st}.json'), encoding='utf-8'))
    z = np.load(os.path.join(work, 'cardio', 'coeur.npz')); SH = tuple(z['shape'])
    od = os.path.join(out, 'cardio'); os.makedirs(od, exist_ok=True)
    res = {'stage': st, 'units': 'mm', 'mm_per_voxel': s_mm, 'center_voxel': c.tolist(),
           'source_parametres': f'embryo3d/cardio_points/{st}.json', 'confiance_coeur': pts.get('confiance_coeur', ''), 'segments': []}
    for key, name, col, alpha in (('coeur_plein', 'coeur', (0.80, 0.15, 0.20), 1.0), ('cavites_cardiaques', 'cavites_cardiaques', (0.95, 0.55, 0.55), 1.0)):
        m = np.unpackbits(z[key])[:int(np.prod(SH))].reshape(SH).astype(bool)
        mesh = mask_to_mesh(m, s_mm, sigma=1.0, target_faces=80000 if name == 'coeur' else 30000)
        if mesh is None: print(st, name, 'vide'); continue
        mesh.apply_translation(-np.array([c[0], c[2], -c[1]]) * s_mm)
        fn = f'{st}_{name}.ply'; mesh.export(os.path.join(od, fn))
        res['segments'].append({'name': name, 'file': fn, 'color': list(col), 'alpha': alpha, 'faces': int(len(mesh.faces)),
                                'volume_mm3': float(m.sum()) * s_mm ** 3})
        print(f'{st} {name:20s} {len(mesh.faces):6d} faces {m.sum() * s_mm ** 3:8.3f} mm3')
    if os.path.exists(os.path.join(out, f'{st}_vaisseaux.ply')):
        res['vaisseaux_pipeline'] = {'file': f'../{st}_vaisseaux.ply', 'confiance': 'détection automatique, fragmentaire'}
    json.dump(res, open(os.path.join(od, 'manifest_cardio.json'), 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
