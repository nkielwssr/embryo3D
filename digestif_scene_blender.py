"""Scène Blender du tube digestif CS13 → CS20 : une collection par stade (enveloppe translucide et foie du pipeline
pour le repère, segments digestifs colorés), stades côte à côte, rendu de contrôle.
usage : blender -b -P digestif_scene_blender.py -- <sortie.blend> <rendu.png> <dossier_stade> [...]"""
import bpy, sys, os, json, mathutils
argv = sys.argv[sys.argv.index('--') + 1:]
out_blend, out_png, dirs = argv[0], argv[1], argv[2:]
bpy.ops.wm.read_factory_settings(use_empty=True)
sc = bpy.context.scene
sc.unit_settings.system = 'METRIC'; sc.unit_settings.scale_length = 0.001; sc.unit_settings.length_unit = 'MILLIMETERS'

def mat(name, rgb, alpha=1.0):
    m = bpy.data.materials.get(name)
    if m: return m
    m = bpy.data.materials.new(name); m.use_nodes = True
    b = m.node_tree.nodes['Principled BSDF']; b.inputs['Base Color'].default_value = (*rgb, 1)
    b.inputs['Alpha'].default_value = alpha; b.inputs['Roughness'].default_value = 0.45
    if alpha < 1: m.surface_render_method = 'BLENDED'
    return m

def load(path, name, coll, material):
    bpy.ops.wm.ply_import(filepath=path); o = bpy.context.selected_objects[0]; o.name = name
    for c in list(o.users_collection): c.objects.unlink(o)
    coll.objects.link(o); o.data.materials.clear(); o.data.materials.append(material)
    for p in o.data.polygons: p.use_smooth = True
    return o

x = 0.0; spans = []
for d in dirs:
    out = os.path.join(d, 'out'); md = json.load(open(os.path.join(out, 'digestif', 'manifest_digestif.json'), encoding='utf-8'))
    st = md['stage']
    coll = bpy.data.collections.new(f'{st} - tube digestif'); sc.collection.children.link(coll)
    ctx = bpy.data.collections.new(f'{st} - repères (enveloppe, foie)'); coll.children.link(ctx)
    objs = []
    for nom, rgb, a in (('enveloppe', (0.92, 0.80, 0.74), 0.10), ('foie', (0.55, 0.20, 0.15), 0.12)):
        p = os.path.join(out, f'{st}_{nom}.ply')
        if os.path.exists(p): objs.append(load(p, f'{st} {nom}', ctx, mat(f'repere_{nom}', rgb, a)))
    for s in md['segments']:
        objs.append(load(os.path.join(out, 'digestif', s['file']), f"{st} {s['name']}", coll, mat(s['name'], s['color'])))
    xs = [(o.matrix_world @ mathutils.Vector(v))[0] for o in objs for v in o.bound_box]
    w = max(xs) - min(xs); shift = x - min(xs)
    for o in objs: o.location.x += shift
    spans.append((st, x, x + w)); x += w + 3.0

# caméra orthographique de profil droit sur l'ensemble
allo = [o for o in sc.objects if o.type == 'MESH']
bb = [(o.matrix_world @ mathutils.Vector(v)) for o in allo for v in o.bound_box]
lo = mathutils.Vector([min(v[i] for v in bb) for i in range(3)]); hi = mathutils.Vector([max(v[i] for v in bb) for i in range(3)])
ctr = (lo + hi) / 2
cam = bpy.data.objects.new('Caméra', bpy.data.cameras.new('cam')); sc.collection.objects.link(cam); sc.camera = cam
cam.data.type = 'ORTHO'; cam.data.ortho_scale = (hi.x - lo.x) * 1.05; cam.data.clip_end = 1000
cam.location = (ctr.x, ctr.y - 200, ctr.z); cam.rotation_euler = (1.5708, 0, 0)       # vue de face (depuis -Y)
for dvec, e in [((0.4, -1, 0.8), 3.0), ((-0.6, -0.3, 0.2), 1.0)]:
    L = bpy.data.objects.new('Soleil', bpy.data.lights.new('sun', 'SUN')); L.data.energy = e; sc.collection.objects.link(L)
    L.rotation_euler = mathutils.Vector(dvec).to_track_quat('Z', 'Y').to_euler()
w = bpy.data.worlds.new('Monde'); w.use_nodes = True; w.node_tree.nodes['Background'].inputs[0].default_value = (1, 1, 1, 1); sc.world = w
sc.view_settings.view_transform = 'Standard'; sc.render.engine = 'BLENDER_EEVEE'; sc.render.resolution_x = 2000
sc.render.resolution_y = int(2000 * (hi.z - lo.z) / (hi.x - lo.x) * 1.1) + 40
txt = bpy.data.texts.new('LIRE_MOI_digestif')
txt.write('Tube digestif CS13-CS20, reconstruit par points de passage dans les volumes HREM du pipeline embryo3d.\n'
          'Confiance : oesophage et estomac bonnes ; duodénum moyenne ; intestin moyen (anse, hernie) et postérieur approximatifs.\n'
          'Points modifiables : embryo3d/digestif_points/<CS>.json, puis digestif_build.py et digestif_export.py.\n')
bpy.ops.wm.save_as_mainfile(filepath=out_blend)
sc.render.filepath = out_png; bpy.ops.render.render(write_still=True)
print('SCENE OK', spans)
