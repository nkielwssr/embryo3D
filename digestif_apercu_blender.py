"""Aperçu Blender du tube digestif d'un stade : enveloppe translucide + segments colorés, vues de profil droit et de face.
usage : blender -b -P digestif_apercu_blender.py -- <work_dir> <sortie_prefixe>
Coordonnées PLY en voxels (x = gauche-droite, s = haut-bas vers le bas, y = ventral-dorsal) ; converties en
Blender X = x, Y = y (dorsal +), Z = -s."""
import bpy, sys, os, math, glob, mathutils
argv = sys.argv[sys.argv.index('--') + 1:]
work, prefix = argv[0], argv[1]
dd = os.path.join(work, 'digestif')
bpy.ops.wm.read_factory_settings(use_empty=True)
sc = bpy.context.scene
cols = {'oesophage': (0.85, 0.35, 0.75), 'estomac': (0.95, 0.55, 0.15), 'duodenum': (0.2, 0.75, 0.3),
        'intestin_moyen': (0.95, 0.85, 0.2), 'intestin_posterieur': (0.25, 0.55, 0.95), 'rectum': (0.3, 0.3, 0.9),
        'pharynx': (0.9, 0.3, 0.3), 'foie': (0.45, 0.12, 0.1)}
def mat(name, rgb, alpha=1.0):
    m = bpy.data.materials.new(name); m.use_nodes = True
    b = m.node_tree.nodes['Principled BSDF']; b.inputs['Base Color'].default_value = (*rgb, 1); b.inputs['Alpha'].default_value = alpha
    b.inputs['Roughness'].default_value = 0.45
    if alpha < 1: m.surface_render_method = 'BLENDED'
    return m
objs = []
for p in sorted(glob.glob(os.path.join(dd, '*.ply'))):
    nom = os.path.splitext(os.path.basename(p))[0]
    bpy.ops.wm.ply_import(filepath=p); o = bpy.context.selected_objects[0]; o.name = nom
    o.data.transform(mathutils.Matrix(((1, 0, 0, 0), (0, 0, 1, 0), (0, -1, 0, 0), (0, 0, 0, 1))))   # (x,s,y) -> (x, y, -s)
    for poly in o.data.polygons: poly.use_smooth = True
    if nom.startswith('enveloppe'):
        o.data.materials.append(mat('peau', (0.9, 0.8, 0.75), 0.12))
    else:
        o.data.materials.append(mat(nom, cols.get(nom, (0.7, 0.7, 0.7))))
    objs.append(o)
env = [o for o in objs if o.name.startswith('enveloppe')]
ref = env[0] if env else objs[0]
bb = [ref.matrix_world @ mathutils.Vector(c) for c in ref.bound_box]
ctr = sum(bb, mathutils.Vector()) / 8; size = max((max(v[i] for v in bb) - min(v[i] for v in bb)) for i in range(3))
cam = bpy.data.objects.new('cam', bpy.data.cameras.new('cam')); sc.collection.objects.link(cam); sc.camera = cam
cam.data.type = 'ORTHO'; cam.data.ortho_scale = size * 1.1; cam.data.clip_end = size * 10
for d, e in [((1, -1, 1), 3.0), ((-1, 0.5, 0.3), 1.2)]:
    L = bpy.data.objects.new('sun', bpy.data.lights.new('sun', 'SUN')); L.data.energy = e; sc.collection.objects.link(L)
    L.rotation_euler = mathutils.Vector(d).to_track_quat('Z', 'Y').to_euler()
w = bpy.data.worlds.new('w'); w.use_nodes = True; w.node_tree.nodes['Background'].inputs[0].default_value = (1, 1, 1, 1); sc.world = w
sc.render.engine = 'BLENDER_EEVEE'; sc.render.resolution_x = 700; sc.render.resolution_y = 900
vues = {'profil_droit': mathutils.Vector((1, 0, 0)), 'face': mathutils.Vector((0, -1, 0))}
for k, dirv in vues.items():
    cam.location = ctr + dirv * size * 3
    cam.rotation_euler = (-dirv).to_track_quat('-Z', 'Y').to_euler()
    sc.render.filepath = f'{prefix}_{k}.png'; bpy.ops.render.render(write_still=True)
print('APERCU OK')
