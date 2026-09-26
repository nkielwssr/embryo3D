"""Exporte les aortes tracées (work/cardio/vaisseaux.npz) en PLY alignés sur le pipeline : <out>/cardio/<CS>_<segment>.ply
+ <out>/cardio/manifest_vaisseaux.json. usage : python vaisseaux_export.py <dossier_stade> [...]"""
import numpy as np, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from meshexport import mask_to_mesh
from tubes_morph import couleur as couleur_tube      # arcs aortiques, sac aortique, tronc artériel : légende commune
HERE = os.path.dirname(os.path.abspath(__file__))
COUL = {'aorte_dorsale_gauche': (0.90, 0.10, 0.10), 'aorte_dorsale_droite': (0.95, 0.35, 0.20), 'aorte_commune': (0.75, 0.05, 0.10),
        'aorte_dorsale': (0.85, 0.10, 0.10),
        'cardinale_anterieure_gauche': (0.35, 0.55, 0.95), 'cardinale_anterieure_droite': (0.45, 0.65, 1.0),
        'cardinale_commune_gauche': (0.10, 0.20, 0.70), 'ombilicale_gauche': (0.55, 0.20, 0.75), 'ombilicale_commune': (0.45, 0.15, 0.65), 'vitelline_gauche': (0.20, 0.60, 0.55),
        'cardinale_posterieure_gauche': (0.15, 0.30, 0.90), 'cardinale_posterieure_droite': (0.25, 0.45, 0.95)}
for d in sys.argv[1:]:
    work, out = os.path.join(d, 'work'), os.path.join(d, 'out')
    man = json.load(open(os.path.join(out, 'manifest.json'), encoding='utf-8'))
    st, s_mm, c = man['stage'], man['mm_per_voxel'], np.array(man['center_voxel'])
    pts = json.load(open(os.path.join(HERE, 'vaisseaux_points', f'{st}.json'), encoding='utf-8'))
    z = np.load(os.path.join(work, 'cardio', 'vaisseaux.npz')); SH = tuple(z['shape'])
    od = os.path.join(out, 'cardio'); os.makedirs(od, exist_ok=True)
    res = {'stage': st, 'units': 'mm', 'mm_per_voxel': s_mm, 'center_voxel': c.tolist(),
           'source_points': f'embryo3d/vaisseaux_points/{st}.json', 'confiance': pts.get('confiance', ''), 'segments': []}
    for k in z.files:
        if k in ('shape', 'union'): continue
        m = np.unpackbits(z[k])[:int(np.prod(SH))].reshape(SH).astype(bool)
        mesh = mask_to_mesh(m, s_mm, sigma=1.0, target_faces=20000)
        if mesh is None: continue
        mesh.apply_translation(-np.array([c[0], c[2], -c[1]]) * s_mm)
        fn = f'{st}_{k}.ply'; mesh.export(os.path.join(od, fn))
        res['segments'].append({'name': k, 'file': fn, 'color': list(COUL.get(k) or couleur_tube(k)), 'faces': int(len(mesh.faces)),
                                'volume_mm3': float(m.sum()) * s_mm ** 3})
        print(f'{st} {k:22s} {len(mesh.faces):6d} faces')
    json.dump(res, open(os.path.join(od, 'manifest_vaisseaux.json'), 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
