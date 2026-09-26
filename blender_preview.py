# usage: blender -b -P blender_preview.py -- in.glb out.png [azimuth_deg]
import bpy, sys, math, os
argv = sys.argv[sys.argv.index('--')+1:]
src, out = argv[0], argv[1]; az = float(argv[2]) if len(argv) > 2 else 30.0
bpy.ops.wm.read_factory_settings(use_empty=True)
ext = os.path.splitext(src)[1].lower()
if ext == '.glb': bpy.ops.import_scene.gltf(filepath=src)
elif ext == '.ply': bpy.ops.wm.ply_import(filepath=src)
elif ext == '.obj': bpy.ops.wm.obj_import(filepath=src)
objs = [o for o in bpy.context.scene.objects if o.type == 'MESH']
for o in objs:
    if not o.data.materials:
        m = bpy.data.materials.new('m'); m.diffuse_color = (0.85, 0.6, 0.5, 1); o.data.materials.append(m)
# bounds
import mathutils
mins = mathutils.Vector((1e9,)*3); maxs = mathutils.Vector((-1e9,)*3)
for o in objs:
    for c in o.bound_box:
        w = o.matrix_world @ mathutils.Vector(c)
        mins = mathutils.Vector(map(min, mins, w)); maxs = mathutils.Vector(map(max, maxs, w))
center = (mins+maxs)/2; size = max(maxs-mins)
scene = bpy.context.scene
cam = bpy.data.cameras.new('cam'); camo = bpy.data.objects.new('cam', cam); scene.collection.objects.link(camo); scene.camera = camo
cam.type = 'ORTHO'; cam.ortho_scale = size*1.15
a = math.radians(az)
camo.location = center + mathutils.Vector((math.sin(a)*size*3, -math.cos(a)*size*3, size*0.3))
camo.rotation_euler = (center - camo.location).to_track_quat('-Z','Y').to_euler()
sun = bpy.data.lights.new('sun','SUN'); sun.energy = 3; so = bpy.data.objects.new('sun', sun); scene.collection.objects.link(so); so.rotation_euler = (math.radians(50), math.radians(20), math.radians(az+40))
scene.render.engine = 'BLENDER_WORKBENCH'
scene.display.shading.light = 'STUDIO'; scene.display.shading.color_type = 'MATERIAL'
scene.display.shading.show_cavity = True
scene.render.resolution_x = 900; scene.render.resolution_y = 1200
scene.render.filepath = out; scene.render.image_settings.file_format = 'PNG'
bpy.ops.render.render(write_still=True)
