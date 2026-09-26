"""Assemble les manifests de plusieurs stades (du plus jeune au plus âgé) en un manifest multi-stades pour blender_build_scene.py.
usage : python make_master_manifest.py out_master.json CS13.f4v/out CS14_f4v/out ... CS20_F4V/out"""
import json, sys, os
out = sys.argv[1]; dirs = sys.argv[2:]
master = {'stages': []}
for d in dirs:
    m = json.load(open(os.path.join(d, 'manifest.json'), encoding='utf-8'))
    rel = os.path.relpath(os.path.abspath(d), os.path.dirname(os.path.abspath(out)))
    for s in m['structures']:
        s['file'] = os.path.join(rel, s['file']).replace(os.sep, '/')
    master['stages'].append(m)
json.dump(master, open(out, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
print('stages:', [m['stage'] for m in master['stages']], '->', out)
