"""Scène Blender des tubes morphables (tube digestif + aortes), CS13 → CS20.
Lit embryons_3D/tubes_morph.npz : un objet par structure (N anneaux × M sommets, extrémités fermées), une shape key par stade,
pilotée par la propriété de scène « stage » (0 = CS13 … 6 = CS20), comme la scène maître ; mêmes keyframes
(palier de 12 images puis transition de 48 images) : Lecture fait défiler les stades.
usage : blender -b -P tubes_morph_blender.py -- <sortie.blend> [<dossier_rendus>]
Réutilisable par la scène maître : bpy.data.libraries.load(<sortie.blend>) → collection « Tubes morphables »."""
import bpy, sys, os, json, numpy as np, mathutils
argv = sys.argv[sys.argv.index('--') + 1:]
out_blend = argv[0]; rendus = argv[1] if len(argv) > 1 else None
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
Z = np.load(os.path.join(ROOT, 'embryons_3D', 'tubes_morph.npz'))
meta = json.load(open(os.path.join(ROOT, 'embryons_3D', 'tubes_morph.json'), encoding='utf-8'))
ST = meta['stades']; N, M = meta['N'], meta['M']
COUL = {'oesophage': (0.85, 0.40, 0.70), 'estomac': (0.95, 0.60, 0.20), 'duodenum': (0.30, 0.75, 0.35),
        'intestin_moyen': (0.95, 0.85, 0.25), 'intestin_posterieur': (0.30, 0.55, 0.95),
        'aorte_dorsale_gauche': (0.90, 0.10, 0.10), 'aorte_dorsale_droite': (0.95, 0.35, 0.20), 'aorte_commune': (0.75, 0.05, 0.10)}

bpy.ops.wm.read_factory_settings(use_empty=True)
sc = bpy.context.scene
sc.unit_settings.system = 'METRIC'; sc.unit_settings.scale_length = 0.001; sc.unit_settings.length_unit = 'MILLIMETERS'
coll = bpy.data.collections.new('Tubes morphables'); sc.collection.children.link(coll)
sc['stage'] = 0.0
ui = sc.id_properties_ui('stage'); ui.update(min=0.0, max=float(len(ST) - 1), soft_min=0.0, soft_max=float(len(ST) - 1),
                                              description='Stade : ' + ', '.join(f'{i} = {s}' for i, s in enumerate(ST)))

faces = []
for i in range(N - 1):
    for j in range(M):
        a, b = i * M + j, i * M + (j + 1) % M; faces.append((a, b, b + M, a + M))
cap0, cap1 = N * M, N * M + 1
for j in range(M): faces.append((cap0, (j + 1) % M, j)); faces.append((cap1, (N - 1) * M + j, (N - 1) * M + (j + 1) % M))

def verts(A):  # A (N, M, 3) -> sommets + deux centres de bouchons
    v = A.reshape(-1, 3); return np.vstack([v, A[0].mean(0), A[-1].mean(0)])

def drv(id_data, path, expr):
    d = id_data.driver_add(path).driver; d.type = 'SCRIPTED'
    var = d.variables.new(); var.name = 'st'; var.type = 'SINGLE_PROP'
    var.targets[0].id_type = 'SCENE'; var.targets[0].id = sc; var.targets[0].data_path = '["stage"]'; d.expression = expr

for nom in Z.files:
    A = Z[nom]                                          # (7, N, M, 3)
    me = bpy.data.meshes.new(nom); me.from_pydata(verts(A[-1]).tolist(), [], faces); me.update()
    o = bpy.data.objects.new(f'{nom} (morph)', me); coll.objects.link(o)
    for p in me.polygons: p.use_smooth = True
    m = bpy.data.materials.new(nom); m.use_nodes = True
    m.node_tree.nodes['Principled BSDF'].inputs['Base Color'].default_value = (*COUL.get(nom, (0.7, 0.7, 0.7)), 1)
    me.materials.append(m)
    o.shape_key_add(name='Basis', from_mix=False)
    for i, st in enumerate(ST):
        k = o.shape_key_add(name=st, from_mix=False)
        k.data.foreach_set('co', verts(A[i]).astype(np.float32).ravel())
        drv(k, 'value', f'max(0, 1 - abs(st - {i}))')
    # la Basis vaut la forme CS20 : on la ramène à zéro pour que les clés soient absolues
    base = verts(np.zeros_like(A[-1])).astype(np.float32).ravel()
    o.data.shape_keys.key_blocks['Basis'].data.foreach_set('co', base)
    for v, c in zip(me.vertices, base.reshape(-1, 3)): v.co = c
    o['presence'] = json.dumps(meta['presence'][nom])

# keyframes identiques à la scène maître : palier 12 images, transition 48 images, 24 i/s
sc.render.fps = 24; f = 1; P, T = 12, 48
for i in range(len(ST)):
    sc['stage'] = float(i); sc.keyframe_insert(data_path='["stage"]', frame=f); f += P
    sc.keyframe_insert(data_path='["stage"]', frame=f)
    if i < len(ST) - 1: f += T
act = sc.animation_data.action; fcs = []
try: fcs = list(act.fcurves)
except Exception:
    for layer in act.layers:
        for strip in layer.strips:
            for cb in strip.channelbags: fcs += list(cb.fcurves)
for fc in fcs:
    for kp in fc.keyframe_points: kp.interpolation = 'LINEAR'
sc.frame_start, sc.frame_end = 1, f; sc.frame_set(1)

cam = bpy.data.objects.new('Caméra', bpy.data.cameras.new('cam')); sc.collection.objects.link(cam); sc.camera = cam
cam.data.type = 'ORTHO'; cam.data.ortho_scale = 22; cam.location = (80, 0, 0); cam.rotation_euler = (1.5708, 0, 1.5708)
for dvec, e in [((1, -0.5, 1), 3.0), ((-1, 0.5, 0.3), 1.0)]:
    L = bpy.data.objects.new('Soleil', bpy.data.lights.new('sun', 'SUN')); L.data.energy = e; sc.collection.objects.link(L)
    L.rotation_euler = mathutils.Vector(dvec).to_track_quat('Z', 'Y').to_euler()
w = bpy.data.worlds.new('Monde'); w.use_nodes = True; w.node_tree.nodes['Background'].inputs[0].default_value = (1, 1, 1, 1); sc.world = w
sc.view_settings.view_transform = 'Standard'; sc.render.engine = 'BLENDER_EEVEE'; sc.render.resolution_x = 700; sc.render.resolution_y = 900
txt = bpy.data.texts.new('LIRE_MOI_tubes')
txt.write('Tubes morphables (tube digestif + aortes) CS13-CS20 : topologie commune 64 anneaux x 16 sommets, une shape key par stade, '
          'pilotées par la propriété de scène stage. Espace = lecture. Un segment absent à un stade est réduit à un point sur son raccord. '
          'Formes simplifiées (tube de rayon mesuré) : les maillages détaillés restent dans cardio/digestif_CS13-CS20.blend.\n')
bpy.ops.wm.save_as_mainfile(filepath=out_blend)
if rendus:
    os.makedirs(rendus, exist_ok=True)
    for fr in (1, 61, 121, 181, 241, 301, f):
        sc.frame_set(fr); sc.render.filepath = os.path.join(rendus, f'tubes_morph_{fr:04d}.png'); bpy.ops.render.render(write_still=True)
    sc.frame_set(1); bpy.ops.wm.save_mainfile()
print('TUBES OK', f)
