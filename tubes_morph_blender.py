"""Scène Blender des tubes morphables (pharynx + tube digestif + aortes) et des nappes de mésos, CS13 → CS20.
Lit embryons_3D/tubes_morph.npz + .json : un objet par structure, une shape key par stade, pilotées par la propriété de
scène « stage » (0 = CS13 … 6 = CS20), comme la scène maître ; mêmes keyframes (palier de 12 images puis transition de
48 images) : Lecture fait défiler les stades.
- type « tube » (N anneaux × M sommets, extrémités fermées) : pharynx, poches pharyngiennes, œsophage, estomac, duodénum, intestins,
  aortes dorsales, arcs aortiques, sac aortique, tronc artériel (couleurs lues dans tubes_morph.json, légende de Rana et al. 2014) ;
- type « nappe » (N lignes × K colonnes, bord digestif → bord aortique) : mésos dorsaux, avec un modificateur Solidify
  (épaisseur --epaisseur mm, 0,05 par défaut) et un matériau translucide double face.
Sous-collections « Tube digestif », « Poches pharyngiennes », « Aortes », « Arcs aortiques… », « Mésos dorsaux » dans « Tubes morphables ».
usage : blender -b -P tubes_morph_blender.py -- <sortie.blend> [<dossier_rendus>] [--epaisseur 0.05] [--npz <tubes_morph.npz>]
Réutilisable par la scène maître : bpy.data.libraries.load(<sortie.blend>) → collection « Tubes morphables »."""
import bpy, sys, os, json, numpy as np, mathutils
argv = sys.argv[sys.argv.index('--') + 1:]
def opt(nom, defaut, conv=float):
    if nom in argv:
        i = argv.index(nom); v = conv(argv[i + 1]); del argv[i:i + 2]; return v
    return defaut
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
epaisseur = opt('--epaisseur', 0.05); npz = opt('--npz', os.path.join(ROOT, 'embryons_3D', 'tubes_morph.npz'), str)
out_blend = argv[0]; rendus = argv[1] if len(argv) > 1 else None
Z = np.load(npz)
meta = json.load(open(os.path.splitext(npz)[0] + '.json', encoding='utf-8'))
ST = meta['stades']; STRUCTS = meta.get('structures', {})
def info(nom): return STRUCTS.get(nom, {'type': 'tube', 'famille': 'aorte' if nom.startswith('aorte') else 'digestif'})
COUL = {'pharynx': (0.80, 0.25, 0.25), 'oesophage': (0.85, 0.40, 0.70), 'estomac': (0.95, 0.60, 0.20), 'duodenum': (0.30, 0.75, 0.35),
        'intestin_moyen': (0.95, 0.85, 0.25), 'intestin_posterieur': (0.30, 0.55, 0.95),
        'aorte_dorsale_gauche': (0.90, 0.10, 0.10), 'aorte_dorsale_droite': (0.95, 0.35, 0.20), 'aorte_commune': (0.75, 0.05, 0.10),
        'meso_oesophage': (0.96, 0.80, 0.72), 'mesogastre_dorsal': (0.96, 0.76, 0.66), 'mesoduodenum': (0.94, 0.78, 0.70),
        'mesentere': (0.97, 0.82, 0.74), 'mesocolon_dorsal': (0.93, 0.74, 0.68)}
SOUS_COLL = {'digestif': 'Tube digestif', 'poche': 'Poches pharyngiennes', 'aorte': 'Aortes', 'arc': 'Arcs aortiques, sac aortique, tronc artériel', 'meso': 'Mésos dorsaux'}

bpy.ops.wm.read_factory_settings(use_empty=True)
sc = bpy.context.scene
sc.unit_settings.system = 'METRIC'; sc.unit_settings.scale_length = 0.001; sc.unit_settings.length_unit = 'MILLIMETERS'
coll = bpy.data.collections.new('Tubes morphables'); sc.collection.children.link(coll)
sous = {}
for fam, nom in SOUS_COLL.items():
    sous[fam] = bpy.data.collections.new(nom); coll.children.link(sous[fam])
sc['stage'] = 0.0
ui = sc.id_properties_ui('stage'); ui.update(min=0.0, max=float(len(ST) - 1), soft_min=0.0, soft_max=float(len(ST) - 1),
                                              description='Stade : ' + ', '.join(f'{i} = {s}' for i, s in enumerate(ST)))

def faces_tube(N, M):
    """quads entre anneaux successifs (anneau fermé) + deux bouchons en éventail sur les centres des extrémités"""
    faces = []
    for i in range(N - 1):
        for j in range(M):
            a, b = i * M + j, i * M + (j + 1) % M; faces.append((a, b, b + M, a + M))
    cap0, cap1 = N * M, N * M + 1
    for j in range(M): faces.append((cap0, (j + 1) % M, j)); faces.append((cap1, (N - 1) * M + j, (N - 1) * M + (j + 1) % M))
    return faces

def faces_nappe(N, K):
    """quads d'une grille ouverte N lignes × K colonnes"""
    return [(i * K + j, i * K + j + 1, (i + 1) * K + j + 1, (i + 1) * K + j) for i in range(N - 1) for j in range(K - 1)]

def verts(A, typ):
    """A (N, M|K, 3) -> sommets ; un tube reçoit en plus les deux centres de bouchons"""
    v = A.reshape(-1, 3)
    return np.vstack([v, A[0].mean(0), A[-1].mean(0)]) if typ == 'tube' else v

def drv(id_data, path, expr):
    d = id_data.driver_add(path).driver; d.type = 'SCRIPTED'
    var = d.variables.new(); var.name = 'st'; var.type = 'SINGLE_PROP'
    var.targets[0].id_type = 'SCENE'; var.targets[0].id = sc; var.targets[0].data_path = '["stage"]'; d.expression = expr

def materiau(nom, rgb, alpha=1.0):
    m = bpy.data.materials.new(nom); m.use_nodes = True
    b = m.node_tree.nodes['Principled BSDF']; b.inputs['Base Color'].default_value = (*rgb, 1); b.inputs['Roughness'].default_value = 0.5
    if alpha < 1:
        b.inputs['Alpha'].default_value = alpha; m.surface_render_method = 'BLENDED'; m.use_backface_culling = False
    return m

n_obj = {'tube': 0, 'nappe': 0}
for nom in Z.files:
    A = Z[nom]; inf = info(nom); typ = inf['type']; fam = inf.get('famille', 'digestif')
    if A.ndim != 4 or A.shape[0] != len(ST): print('structure', nom, ': forme inattendue', A.shape, 'ignorée'); continue
    _, Nn, Mn, _ = A.shape
    faces = faces_tube(Nn, Mn) if typ == 'tube' else faces_nappe(Nn, Mn)
    me = bpy.data.meshes.new(nom); me.from_pydata(verts(A[-1], typ).tolist(), [], faces); me.update()
    o = bpy.data.objects.new(f'{nom} (morph)', me); sous.get(fam, coll).objects.link(o)
    for p in me.polygons: p.use_smooth = True
    me.materials.append(materiau(nom, tuple(inf.get('couleur') or COUL.get(nom, (0.7, 0.7, 0.7))), 0.75 if typ == 'nappe' else 1.0))
    o.shape_key_add(name='Basis', from_mix=False)
    for i, st in enumerate(ST):
        k = o.shape_key_add(name=st, from_mix=False)
        k.data.foreach_set('co', verts(A[i], typ).astype(np.float32).ravel())
        drv(k, 'value', f'max(0, 1 - abs(st - {i}))')
    # la Basis vaut la forme CS20 : on la ramène à zéro pour que les clés soient absolues (mélange = interpolation linéaire)
    base = np.zeros_like(verts(A[-1], typ)).astype(np.float32).ravel()
    o.data.shape_keys.key_blocks['Basis'].data.foreach_set('co', base)
    for v, c in zip(me.vertices, base.reshape(-1, 3)): v.co = c
    if typ == 'nappe':
        mod = o.modifiers.new('Épaisseur', 'SOLIDIFY'); mod.thickness = epaisseur; mod.offset = 0.0; mod.use_even_offset = True
        o['segment'] = inf.get('segment', ''); o['attache'] = inf.get('attache', '')
    o['presence'] = json.dumps(meta['presence'][nom]); o['type'] = typ; o['famille'] = fam
    n_obj[typ] += 1

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
txt.write('Tubes morphables CS13-CS20 : pharynx, tube digestif (oesophage, estomac, duodenum, intestins moyen et posterieur) et aortes\n'
          '(topologie commune %d anneaux x %d sommets), mesos dorsaux (meso-oesophage, mesogastre dorsal, mesoduodenum, mesentere,\n'
          'mesocolon dorsal : nappes %d lignes x %d colonnes du bord dorsal du tube au bord ventral de l\'axe aortique, Solidify %.3f mm).\n'
          'Une shape key par stade, pilotees par la propriete de scene stage (0 = %s ... %d = %s). Espace = lecture.\n'
          'Un segment absent a un stade est reduit a un point sur son raccord ; un meso sans aorte au stade reste plaque sur le tube\n'
          '(largeur nulle) et s\'ouvre pendant le morphing ; les lignes trop loin de l\'aorte (> %s mm) ou au-dela de ses extremites\n'
          'tracees sont reduites. Formes simplifiees (tube de rayon mesure) : les maillages detailles restent dans cardio/digestif_CS13-CS20.blend.\n'
          % (meta['N'], meta['M'], meta['N'], meta.get('K', 0), epaisseur, ST[0], len(ST) - 1, ST[-1],
             (meta.get('mesos') or {}).get('portee_max_mm', '?')))
bpy.ops.wm.save_as_mainfile(filepath=out_blend)
if rendus:
    os.makedirs(rendus, exist_ok=True)
    for fr in (1, 61, 121, 181, 241, 301, f):
        sc.frame_set(fr); sc.render.filepath = os.path.join(rendus, f'tubes_morph_{fr:04d}.png'); bpy.ops.render.render(write_still=True)
    sc.frame_set(1); bpy.ops.wm.save_mainfile()
print('TUBES OK', n_obj['tube'], 'tubes,', n_obj['nappe'], 'nappes, derniere image', f)
