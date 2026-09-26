# rendu de contrôle : squelette axial + SNC translucide + cœur/foie, vue de 3/4 dorsale. blender -b scene.blend -P blender_render_axial.py -- out.png
import bpy, sys, math, mathutils
out = sys.argv[sys.argv.index('--')+1]
scene = bpy.context.scene
KEEP = {'somites', 'squelette_axial_cartilage', 'notochorde', 'snc', 'coeur_detoure', 'foie', 'vaisseaux_aorte', 'vaisseaux_cardinales', 'chondrocrane'}
for o in scene.objects:
    if o.type == 'MESH':
        o.hide_render = o.get('structure') not in KEEP
        if o.get('structure') == 'snc':
            for m in o.data.materials:
                m.diffuse_color = (*m.diffuse_color[:3], 0.25)
objs = [o for o in scene.objects if o.type == 'MESH' and not o.hide_render]
mins = mathutils.Vector((1e9,)*3); maxs = mathutils.Vector((-1e9,)*3)
for o in objs:
    for c in o.bound_box:
        w = o.matrix_world @ mathutils.Vector(c); mins = mathutils.Vector(map(min, mins, w)); maxs = mathutils.Vector(map(max, maxs, w))
center = (mins+maxs)/2; size = max(maxs-mins)
cam = scene.camera; cam.data.type = 'ORTHO'; cam.data.ortho_scale = size*1.1
a = math.radians(140); cam.location = center + mathutils.Vector((math.sin(a)*size*3, -math.cos(a)*size*3, size*0.5))
cam.rotation_euler = (center - cam.location).to_track_quat('-Z','Y').to_euler()
scene.render.engine = 'BLENDER_WORKBENCH'; sh = scene.display.shading; sh.light = 'STUDIO'; sh.color_type = 'MATERIAL'; sh.show_cavity = True; sh.show_shadows = True
sh.show_xray = True; sh.xray_alpha = 0.6
scene.render.resolution_x = 900; scene.render.resolution_y = 1000
scene.render.filepath = out; bpy.ops.render.render(write_still=True)
