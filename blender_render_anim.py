# blender -b scene.blend -P blender_render_anim.py -- sortie.mp4 [largeur hauteur] [azimut_deg]
# Rendu Workbench de l'animation « stage » de la scène maître : caméra orthographique FIXE cadrée sur le plus grand stade,
# fond clair, étiquette du stade. Séquence PNG (ce Blender n'a pas de sortie FFMPEG),
# à encoder ensuite avec encode_frames.py.
import bpy, sys, os, math, mathutils, shutil
argv = sys.argv[sys.argv.index('--')+1:]; out = argv[0]
W, H = (int(argv[1]), int(argv[2])) if len(argv) > 2 else (1280, 720); az = float(argv[3]) if len(argv) > 3 else 90.0
scene = bpy.context.scene; n = int(round(scene['stage'])) + 1 if scene.get('stage') is not None else 1
# clés de la propriété stage : (image, valeur)
keys = []
ad = scene.animation_data
try: fcs = list(ad.action.fcurves)
except Exception: fcs = list(ad.action.layers[0].strips[0].channelbag(ad.action_slot).fcurves)
for fc in fcs:
    if fc.data_path == '["stage"]': keys = [(int(kp.co[0]), kp.co[1]) for kp in fc.keyframe_points]
n = int(round(max(v for _, v in keys))) + 1
stage_names = {}
for o in scene.objects:
    if o.get('stage') and o.get('stage') != 'morph': stage_names.setdefault(o['stage'], None)
names_sorted = sorted(stage_names, key=lambda s: int(''.join(ch for ch in s if ch.isdigit()) or 0))
def visible_bbox(frame):
    scene.frame_set(frame); dg = bpy.context.evaluated_depsgraph_get()
    mins = mathutils.Vector((1e9,)*3); maxs = mathutils.Vector((-1e9,)*3)
    for o in scene.objects:
        if o.type != 'MESH' or o.hide_render or o.evaluated_get(dg).hide_render: continue
        e = o.evaluated_get(dg)
        for c in e.bound_box:
            w = e.matrix_world @ mathutils.Vector(c); mins = mathutils.Vector(map(min, mins, w)); maxs = mathutils.Vector(map(max, maxs, w))
    return (mins+maxs)/2, max(maxs-mins)
cam = scene.camera; cam.data.type = 'ORTHO'; cam.animation_data_clear(); cam.data.animation_data_clear()
a = math.radians(az); dirv = mathutils.Vector((math.sin(a), -math.cos(a), 0.1)).normalized()
rot = (-dirv).to_track_quat('-Z', 'Y'); cam.rotation_euler = rot.to_euler()
right = rot @ mathutils.Vector((1, 0, 0)); up = rot @ mathutils.Vector((0, 1, 0))
# étiquettes : un objet texte par stade, caché hors de son stade (mêmes pilotes que les structures non morphées)
mat_txt = bpy.data.materials.new('etiquette'); mat_txt.diffuse_color = (0.15, 0.15, 0.18, 1)
labels = []
for si in range(n):
    nm = names_sorted[si] if si < len(names_sorted) else 'stade %d' % si
    cu = bpy.data.curves.new('lbl%d' % si, 'FONT'); cu.body = 'Stade de Carnegie %s' % nm[2:] if nm.startswith('CS') else nm; cu.align_x = 'LEFT'
    t = bpy.data.objects.new('lbl%d' % si, cu); scene.collection.objects.link(t); t.data.materials.append(mat_txt)
    fc = t.driver_add('hide_render'); fc.driver.type = 'SCRIPTED'; v = fc.driver.variables.new(); v.name = 'st'; v.type = 'SINGLE_PROP'
    v.targets[0].id_type = 'SCENE'; v.targets[0].id = scene; v.targets[0].data_path = '["stage"]'; fc.driver.expression = 'abs(st - %d) >= 0.5' % si
    labels.append(t)
# caméra FIXE par défaut (consigne utilisateur : l'embryon grandit dans un cadre immobile) : une seule pose et une seule échelle, calculées sur le plus grand
# stade (dernière image). Option --suivre : échelle et position recadrées en douceur sur chaque stade (taille normalisée), clés aux images des paliers.
aspect = W / H; SUIVRE = '--suivre' in argv
def place(center, size, frame=None):
    osc = max(size * 1.15, 0.5); cam.data.ortho_scale = osc; cam.location = center + dirv * (size * 3 + 5)
    for t in labels:
        t.location = center + right * (-0.47 * osc) + up * (-0.46 * osc / aspect) + dirv * (size * 1.5); t.scale = (osc * 0.04,) * 3; t.rotation_euler = rot.to_euler()
    if frame is not None:
        cam.keyframe_insert('location', frame=frame); cam.data.keyframe_insert('ortho_scale', frame=frame)
        for t in labels: t.keyframe_insert('location', frame=frame); t.keyframe_insert('scale', frame=frame)
    return osc
if SUIVRE:
    for frame, val in keys: place(*visible_bbox(frame), frame=frame)
    for idb in [cam, cam.data] + labels:
        try:
            for fc in idb.animation_data.action.fcurves:
                for kp in fc.keyframe_points: kp.interpolation = 'BEZIER'; kp.handle_left_type = kp.handle_right_type = 'AUTO_CLAMPED'
        except Exception: pass
    print('CAMERA suiveuse (taille normalisée par stade)')
else:
    center, size = visible_bbox(scene.frame_end); osc = place(center, size)
    print('CAMERA fixe : centre %s, ortho_scale %.1f mm' % ([round(v, 1) for v in center], osc))
scene.render.engine = 'BLENDER_WORKBENCH'; sh = scene.display.shading; sh.light = 'STUDIO'; sh.color_type = 'MATERIAL'; sh.show_cavity = True; sh.show_shadows = False
sh.background_type = 'VIEWPORT'; sh.background_color = (0.96, 0.96, 0.95)
scene.render.resolution_x = W; scene.render.resolution_y = H; scene.render.resolution_percentage = 100; scene.render.film_transparent = False
frames_dir = os.path.splitext(out)[0] + '_frames'; shutil.rmtree(frames_dir, ignore_errors=True); os.makedirs(frames_dir, exist_ok=True)
scene.render.image_settings.file_format = 'PNG'; scene.render.image_settings.color_mode = 'RGB'
scene.render.filepath = os.path.join(frames_dir, 'f_'); scene.frame_set(scene.frame_start)
bpy.ops.render.render(animation=True)
print('ANIM_DONE', frames_dir, scene.frame_start, scene.frame_end)
