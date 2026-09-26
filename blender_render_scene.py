# Rendu de contrôle (Workbench, éclairage studio, couleurs matériaux) : blender -b scene.blend -P blender_render_scene.py -- out.png [azimut] [enveloppe: 0=cachée 1=visible]
import bpy, sys, math, mathutils
argv = sys.argv[sys.argv.index('--')+1:]
out = argv[0]; az = float(argv[1]) if len(argv) > 1 else 90.0; show_env = (argv[2] == '1') if len(argv) > 2 else False
scene = bpy.context.scene
for o in scene.objects:
    if o.type == 'MESH' and o.get('structure') == 'enveloppe':
        o.hide_render = not show_env
objs = [o for o in scene.objects if o.type == 'MESH' and not o.hide_render]
mins = mathutils.Vector((1e9,)*3); maxs = mathutils.Vector((-1e9,)*3)
for o in objs:
    for c in o.bound_box:
        w = o.matrix_world @ mathutils.Vector(c); mins = mathutils.Vector(map(min, mins, w)); maxs = mathutils.Vector(map(max, maxs, w))
center = (mins+maxs)/2; size = max(maxs-mins)
cam = scene.camera; cam.data.type = 'ORTHO'; cam.data.ortho_scale = size*1.1
a = math.radians(az); cam.location = center + mathutils.Vector((math.sin(a)*size*3, -math.cos(a)*size*3, size*0.15))
cam.rotation_euler = (center - cam.location).to_track_quat('-Z','Y').to_euler()
scene.render.engine = 'BLENDER_WORKBENCH'
sh = scene.display.shading; sh.light = 'STUDIO'; sh.color_type = 'MATERIAL'; sh.show_cavity = True; sh.show_shadows = True; sh.show_object_outline = True
sh.cavity_type = 'BOTH'
scene.display.render_aa = '8'
scene.render.resolution_x = 1400; scene.render.resolution_y = 1400; scene.render.film_transparent = False
scene.render.filepath = out; bpy.ops.render.render(write_still=True)
