"""Scène Blender de tous les stades extraits (modèles de référence des arcs aortiques).

usage : blender -b -P embryo3d/reference/arcs_aortiques/blender_arcs.py -- extraits/index.json arcs_aortiques.blend [--superposer]
  défaut        : un stade par collection, côte à côte le long de X (du plus jeune au plus âgé), étiquette sous chacun
  --superposer  : stades superposés à l'origine ; la propriété de scène « stade_arcs » (0 = premier stade) n'affiche
                  que le stade le plus proche (pilotes sur la visibilité) : l'animer fait défiler les stades.
Pour un seul stade, blender_build_scene.py lit aussi directement extraits/<CS>/manifest.json.
"""
import json
import os
import re
import sys

import bpy

argv = sys.argv[sys.argv.index('--') + 1:]
index_path, out_blend = argv[0], argv[1]
superposer = '--superposer' in argv
racine = os.path.dirname(os.path.abspath(index_path))
index = json.load(open(index_path, encoding='utf-8'))

manifestes = []
for f in index['fichiers']:
    for fl in f['flux']:
        if fl.get('manifest'):
            manifestes.append(os.path.join(racine, fl['manifest']))


def cle_stade(chemin):
    m = re.match(r'CS(\d+)', json.load(open(chemin, encoding='utf-8'))['stage'])
    return (0, int(m.group(1))) if m else (1, 0)


manifestes.sort(key=cle_stade)

bpy.ops.wm.read_factory_settings(use_empty=True)
scene = bpy.context.scene
scene.unit_settings.system = 'METRIC'
scene.unit_settings.length_unit = 'MILLIMETERS'
scene.unit_settings.scale_length = 0.001


def collection(nom, parent=None):
    parent = parent or scene.collection
    for c in parent.children:
        if c.name == nom:
            return c
    c = bpy.data.collections.new(nom)
    parent.children.link(c)
    return c


def materiau(nom, couleur, alpha):
    m = bpy.data.materials.get(nom)
    if m:
        return m
    m = bpy.data.materials.new(nom)
    m.use_nodes = True
    bsdf = m.node_tree.nodes.get('Principled BSDF')
    bsdf.inputs['Base Color'].default_value = (*couleur, 1.0)
    bsdf.inputs['Roughness'].default_value = 0.45
    if alpha < 1.0:
        bsdf.inputs['Alpha'].default_value = alpha
        if hasattr(m, 'surface_render_method'):
            m.surface_render_method = 'BLENDED'
        elif hasattr(m, 'blend_method'):
            m.blend_method = 'BLEND'
    m.diffuse_color = (*couleur, alpha)
    return m


def pilote_visibilite(obj, i):
    for chemin in ('hide_viewport', 'hide_render'):
        fc = obj.driver_add(chemin)
        d = fc.driver
        d.type = 'SCRIPTED'
        v = d.variables.new()
        v.name = 'st'
        v.targets[0].id_type = 'SCENE'
        v.targets[0].id = scene
        v.targets[0].data_path = '["stade_arcs"]'
        d.expression = f'abs(st - {i}) >= 0.5'


decalage = 0.0
etendues = []
for i, chemin in enumerate(manifestes):
    M = json.load(open(chemin, encoding='utf-8'))
    base = os.path.dirname(chemin)
    col_stade = collection(M['stage'])
    objs = []
    for s in M['structures']:
        p = os.path.join(base, s['file'])
        if not os.path.exists(p):
            print('ABSENT', p)
            continue
        avant = set(bpy.data.objects)
        bpy.ops.wm.ply_import(filepath=p)
        nouveaux = [o for o in set(bpy.data.objects) - avant if o.type == 'MESH']
        if not nouveaux:
            continue
        o = nouveaux[0]
        o.name = f"{M['stage']}_{s['name']}"
        for c in list(o.users_collection):
            c.objects.unlink(o)
        collection(f"{M['stage']} {s.get('collection', 'Divers')}", col_stade).objects.link(o)
        o.data.materials.clear()
        o.data.materials.append(materiau(s['name'] if s.get('structure') else o.name, s.get('color', [0.7, 0.7, 0.7]), s.get('alpha', 1.0)))
        for pg in o.data.polygons:
            pg.use_smooth = True
        o['structure'] = s.get('structure') or ''
        o['nom_original'] = s.get('nom_original', '')
        o['stade'] = M['stage']
        objs.append(o)
    if not objs:
        continue
    xs = [(o.matrix_world @ v.co).x for o in objs for v in o.data.vertices]
    ys = [(o.matrix_world @ v.co).y for o in objs for v in o.data.vertices]
    largeur = max(xs) - min(xs)
    if superposer:
        for o in objs:
            pilote_visibilite(o, i)
    else:
        for o in objs:
            o.location.x += decalage - min(xs)
    txt = bpy.data.curves.new(f"etiquette_{M['stage']}", 'FONT')
    txt.body = M['stage']
    txt.size = max(largeur, 0.5) * 0.12
    t = bpy.data.objects.new(f"etiquette_{M['stage']}", txt)
    t.location = ((0 if superposer else decalage), min(ys) - txt.size * 2, 0)
    col_stade.objects.link(t)
    if superposer:
        pilote_visibilite(t, i)
    etendues.append((M['stage'], largeur))
    decalage += largeur * 1.25 + 0.2

scene['stade_arcs'] = 0.0
ui = scene.id_properties_ui('stade_arcs')
ui.update(min=0.0, max=float(max(len(manifestes) - 1, 0)),
          description='Stade affiché : ' + ', '.join(f'{k} = {json.load(open(c, encoding="utf-8"))["stage"]}' for k, c in enumerate(manifestes)))
bpy.ops.wm.save_as_mainfile(filepath=os.path.abspath(out_blend))
print('SAUVÉ', out_blend, etendues)
