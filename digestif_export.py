"""Exporte le tube digestif (work/digestif/digestif.npz) en PLY alignés sur les maillages du pipeline
(mêmes mm par voxel, même centre, mêmes axes Blender que meshexport.py) : <out>/digestif/<CS>_<segment>.ply
+ <out>/digestif/manifest_digestif.json.
usage : python digestif_export.py <dossier_stade> [<dossier_stade> ...]      (ex. CS20_F4V CS19_f4v)"""
import numpy as np, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from meshexport import mask_to_mesh
from tubes_morph import couleur as couleur_tube      # poches pharyngiennes : couleur commune

COULEURS = {'pharynx': (0.80, 0.25, 0.25), 'oesophage': (0.85, 0.40, 0.70), 'estomac': (0.95, 0.60, 0.20), 'duodenum': (0.30, 0.75, 0.35),
            'intestin_moyen': (0.95, 0.85, 0.25), 'intestin_posterieur': (0.30, 0.55, 0.95)}
CONFIANCE = {'pharynx': 'moyenne', 'oesophage': 'bonne', 'estomac': 'bonne', 'duodenum': 'moyenne',
             'intestin_moyen': 'approximative', 'intestin_posterieur': 'approximative'}

for d in sys.argv[1:]:
    work, out = os.path.join(d, 'work'), os.path.join(d, 'out')
    man = json.load(open(os.path.join(out, 'manifest.json'), encoding='utf-8'))
    stage, s_mm, c = man['stage'], man['mm_per_voxel'], np.array(man['center_voxel'])
    z = np.load(os.path.join(work, 'digestif', 'digestif.npz')); SH = tuple(z['shape'])
    od = os.path.join(out, 'digestif'); os.makedirs(od, exist_ok=True)
    res = {'stage': stage, 'units': 'mm', 'mm_per_voxel': s_mm, 'center_voxel': c.tolist(),
           'source_points': f'embryo3d/digestif_points/{stage}.json', 'segments': []}
    for k in z.files:
        if k in ('shape', 'union'): continue
        m = np.unpackbits(z[k])[:int(np.prod(SH))].reshape(SH).astype(bool)
        mesh = mask_to_mesh(m, s_mm, sigma=1.0, target_faces=40000)
        if mesh is None: print(stage, k, 'vide'); continue
        mesh.apply_translation(-np.array([c[0], c[2], -c[1]]) * s_mm)
        fn = f'{stage}_{k}.ply'; mesh.export(os.path.join(od, fn))
        res['segments'].append({'name': k, 'file': fn, 'color': list(COULEURS.get(k) or couleur_tube(k)),
                                'confiance': CONFIANCE.get(k, 'approximative'), 'faces': int(len(mesh.faces)),
                                'volume_mm3': float(m.sum()) * s_mm ** 3})
        print(f'{stage} {k:22s} {len(mesh.faces):6d} faces {m.sum() * s_mm ** 3:8.3f} mm3')
    json.dump(res, open(os.path.join(od, 'manifest_digestif.json'), 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
