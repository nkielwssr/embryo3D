"""
blender_skin_material.py — matériau « peau d'embryon » inspiré de l'aspect VOKA (voir voka_pages/) :
peau translucide rose pâle, diffusion sous-surface, réseau vasculaire superficiel rougeâtre,
légère brillance humide. Procédural (aucune texture image), donc portable entre machines.

Utilisation (Blender 4.x / 5.x) :
  - en batch :  blender -b scene.blend -P blender_skin_material.py -- [--collection Enveloppe] [--save]
  - dans l'éditeur de scripts : Run Script -> crée/actualise le matériau 'VOKA_peau' et l'applique
    aux objets dont le nom contient 'enveloppe' (insensible à la casse) ou à la collection donnée.
"""
import bpy, sys


def skin_material(name='VOKA_peau'):
    mat = bpy.data.materials.get(name) or bpy.data.materials.new(name)
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    n = nt.nodes
    L = nt.links

    out = n.new('ShaderNodeOutputMaterial'); out.location = (900, 0)
    bsdf = n.new('ShaderNodeBsdfPrincipled'); bsdf.location = (600, 0)
    L.new(bsdf.outputs['BSDF'], out.inputs['Surface'])

    # --- couleurs de base (échantillonnées sur les rendus VOKA : peau ~ (0.93,0.80,0.74), vaisseaux ~ (0.80,0.35,0.35))
    peau = n.new('ShaderNodeRGB'); peau.location = (-200, 300); peau.outputs[0].default_value = (0.93, 0.80, 0.74, 1)
    sang = n.new('ShaderNodeRGB'); sang.location = (-200, 100); sang.outputs[0].default_value = (0.78, 0.33, 0.33, 1)
    ombre = n.new('ShaderNodeRGB'); ombre.location = (-200, -100); ombre.outputs[0].default_value = (0.86, 0.66, 0.62, 1)

    # --- réseau vasculaire : Voronoi « distance to edge » à deux échelles, seuillé finement
    tex = n.new('ShaderNodeTexCoord'); tex.location = (-900, 0)
    v1 = n.new('ShaderNodeTexVoronoi'); v1.location = (-700, 200); v1.feature = 'DISTANCE_TO_EDGE'; v1.inputs['Scale'].default_value = 6.0
    v2 = n.new('ShaderNodeTexVoronoi'); v2.location = (-700, -50); v2.feature = 'DISTANCE_TO_EDGE'; v2.inputs['Scale'].default_value = 18.0
    L.new(tex.outputs['Object'], v1.inputs['Vector']); L.new(tex.outputs['Object'], v2.inputs['Vector'])
    r1 = n.new('ShaderNodeMapRange'); r1.location = (-500, 200); r1.inputs['From Min'].default_value = 0.0; r1.inputs['From Max'].default_value = 0.035; r1.inputs['To Min'].default_value = 1; r1.inputs['To Max'].default_value = 0
    r2 = n.new('ShaderNodeMapRange'); r2.location = (-500, -50); r2.inputs['From Min'].default_value = 0.0; r2.inputs['From Max'].default_value = 0.02; r2.inputs['To Min'].default_value = 0.6; r2.inputs['To Max'].default_value = 0
    L.new(v1.outputs['Distance'], r1.inputs['Value']); L.new(v2.outputs['Distance'], r2.inputs['Value'])
    vmax = n.new('ShaderNodeMath'); vmax.location = (-300, 100); vmax.operation = 'MAXIMUM'
    L.new(r1.outputs['Result'], vmax.inputs[0]); L.new(r2.outputs['Result'], vmax.inputs[1])
    # atténuation aléatoire (les vaisseaux ne couvrent pas tout) : bruit basse fréquence
    noise = n.new('ShaderNodeTexNoise'); noise.location = (-700, -300); noise.inputs['Scale'].default_value = 2.5; noise.inputs['Detail'].default_value = 3
    L.new(tex.outputs['Object'], noise.inputs['Vector'])
    gate = n.new('ShaderNodeMapRange'); gate.location = (-500, -300); gate.inputs['From Min'].default_value = 0.35; gate.inputs['From Max'].default_value = 0.65
    L.new(noise.outputs['Fac'], gate.inputs['Value'])
    vess = n.new('ShaderNodeMath'); vess.location = (-100, 0); vess.operation = 'MULTIPLY'
    L.new(vmax.outputs['Value'], vess.inputs[0]); L.new(gate.outputs['Result'], vess.inputs[1])

    # --- marbrure douce (variations de teinte)
    mottle = n.new('ShaderNodeTexNoise'); mottle.location = (-700, 450); mottle.inputs['Scale'].default_value = 9; mottle.inputs['Detail'].default_value = 4
    L.new(tex.outputs['Object'], mottle.inputs['Vector'])
    mix1 = n.new('ShaderNodeMix'); mix1.location = (100, 300); mix1.data_type = 'RGBA'; mix1.inputs['Factor'].default_value = 0.5
    L.new(mottle.outputs['Fac'], mix1.inputs['Factor']); L.new(peau.outputs[0], mix1.inputs[6]); L.new(ombre.outputs[0], mix1.inputs[7])
    mix2 = n.new('ShaderNodeMix'); mix2.location = (300, 200); mix2.data_type = 'RGBA'
    L.new(vess.outputs['Value'], mix2.inputs['Factor']); L.new(mix1.outputs[2], mix2.inputs[6]); L.new(sang.outputs[0], mix2.inputs[7])
    L.new(mix2.outputs[2], bsdf.inputs['Base Color'])

    # --- translucidité / brillance humide
    def set_in(name, val):
        if name in bsdf.inputs:
            bsdf.inputs[name].default_value = val
    set_in('Subsurface Weight', 0.6)
    set_in('Subsurface Radius', (1.0, 0.35, 0.2))
    set_in('Subsurface Scale', 0.15)
    set_in('Roughness', 0.32)
    set_in('IOR', 1.4)
    set_in('Coat Weight', 0.35)
    set_in('Coat Roughness', 0.15)
    set_in('Specular IOR Level', 0.5)
    # relief fin (pores / plis) via bump léger
    bump = n.new('ShaderNodeBump'); bump.location = (400, -300); bump.inputs['Strength'].default_value = 0.08
    L.new(mottle.outputs['Fac'], bump.inputs['Height']); L.new(bump.outputs['Normal'], bsdf.inputs['Normal'])
    mat.blend_method = 'BLEND' if hasattr(mat, 'blend_method') else mat.blend_method
    return mat


def apply(mat, collection=None):
    objs = []
    if collection and collection in bpy.data.collections:
        objs = [o for o in bpy.data.collections[collection].all_objects if o.type == 'MESH']
    else:
        objs = [o for o in bpy.data.objects if o.type == 'MESH' and 'enveloppe' in o.name.lower()]
    for o in objs:
        if o.data.materials:
            o.data.materials[0] = mat
        else:
            o.data.materials.append(mat)
    print(f'matériau {mat.name} appliqué à {len(objs)} objet(s)')


if __name__ == '__main__':
    argv = sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else []
    coll = argv[argv.index('--collection') + 1] if '--collection' in argv else None
    m = skin_material()
    apply(m, coll)
    if '--save' in argv and bpy.data.filepath:
        bpy.ops.wm.save_mainfile()
