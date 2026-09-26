# blender -b scene.blend -P blender_render_morph.py -- out_prefix stage1 [stage2 ...]   (rendu Workbench sans enveloppe)
import bpy, sys, math, mathutils
argv = sys.argv[sys.argv.index('--')+1:]; garder_env = '--enveloppe' in argv; argv = [a for a in argv if a != '--enveloppe']; prefix = argv[0]; stages = [float(x) for x in argv[1:]]
scene = bpy.context.scene
# la scène est animée (propriété stage + parent « Embryon » normalisant la taille) : on se place à l'IMAGE dont la valeur de stage est la plus
# proche de celle demandée, au lieu de forcer la propriété (sinon le parent resterait à l'échelle de l'image courante)
anim = scene.animation_data is not None and scene.animation_data.action is not None
def frame_for(st):
    if not anim: return None
    best = (1e9, scene.frame_start)
    for f in range(scene.frame_start, scene.frame_end + 1):
        scene.frame_set(f); d = abs(scene['stage'] - st)
        if d < best[0]: best = (d, f)
    return best[1]
def set_stage(st):
    if anim: scene.frame_set(frame_for(st))
    else: scene['stage'] = st; scene.frame_set(scene.frame_current)
if scene.get('stage') is not None: set_stage(float(max(stages)) if stages else 6.0)      # cadrage fixe mesuré sur le plus grand stade
for o in scene.objects:
    if o.type == 'MESH' and o.get('structure') == 'enveloppe' and not garder_env: o.hide_render = True     # --enveloppe : la garder (translucide)
dg = bpy.context.evaluated_depsgraph_get()
objs = [o for o in scene.objects if o.type == 'MESH' and not o.hide_render and not o.hide_viewport and not o.evaluated_get(dg).hide_render]
mins = mathutils.Vector((1e9,)*3); maxs = mathutils.Vector((-1e9,)*3)
for o in objs:
    e = o.evaluated_get(dg)
    for c in e.bound_box:
        w = e.matrix_world @ mathutils.Vector(c); mins = mathutils.Vector(map(min, mins, w)); maxs = mathutils.Vector(map(max, maxs, w))
center = (mins+maxs)/2; size = max(maxs-mins)
cam = scene.camera; cam.data.type = 'ORTHO'; cam.data.ortho_scale = size*1.15
a = math.radians(90); cam.location = center + mathutils.Vector((math.sin(a)*size*3, -math.cos(a)*size*3, size*0.1))
cam.rotation_euler = (center - cam.location).to_track_quat('-Z','Y').to_euler()
scene.render.engine = 'BLENDER_WORKBENCH'; sh = scene.display.shading; sh.light = 'STUDIO'; sh.color_type = 'MATERIAL'; sh.show_cavity = True; sh.show_shadows = True
scene.render.resolution_x = 900; scene.render.resolution_y = 1100
for st in stages:
    set_stage(st)
    bpy.context.view_layer.update()
    for o in objs:
        if o.data.shape_keys:
            for kb in o.data.shape_keys.key_blocks: pass
    scene.frame_set(scene.frame_current)   # force drivers
    scene.render.filepath = '%s_stage%.2f.png' % (prefix, st); bpy.ops.render.render(write_still=True)
