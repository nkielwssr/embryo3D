"""
Construit la scène Blender de l'embryon à partir d'un manifest JSON.
usage :  blender -b -P blender_build_scene.py -- manifest.json out.blend
Le manifest liste les structures : {"stage":"CS20", "units":"mm", "structures":[{"name":..., "file":..., "collection":..., "color":[r,g,b], "alpha":1.0}]}
Si le manifest contient "stages": [ {"stage":..,"structures":[..]}, ... ] (du plus jeune au plus âgé), les structures de même nom
reçoivent des shape keys (topologie du dernier stade, projetée par Shrinkwrap sur les stades précédents) et un slider "stage".
"""
import bpy, sys, json, os, math
argv = sys.argv[sys.argv.index('--')+1:]
manifest_path, out_blend = argv[0], argv[1]
M = json.load(open(manifest_path, encoding='utf-8'))
base = os.path.dirname(os.path.abspath(manifest_path))
bpy.ops.wm.read_factory_settings(use_empty=True)
scene = bpy.context.scene
scene.unit_settings.system = 'METRIC'; scene.unit_settings.length_unit = 'MILLIMETERS'; scene.unit_settings.scale_length = 0.001

def get_coll(name, parent=None):
    parent = parent or scene.collection
    for c in parent.children:
        if c.name == name: return c
    c = bpy.data.collections.new(name); parent.children.link(c); return c

def make_material(name, color, alpha=1.0, roughness=0.45, sss=0.0):
    m = bpy.data.materials.new(name); m.use_nodes = True
    bsdf = m.node_tree.nodes.get('Principled BSDF')
    bsdf.inputs['Base Color'].default_value = (*color, 1.0)
    bsdf.inputs['Roughness'].default_value = roughness
    if 'Subsurface Weight' in bsdf.inputs: bsdf.inputs['Subsurface Weight'].default_value = sss
    if alpha < 1.0:
        bsdf.inputs['Alpha'].default_value = alpha
        m.surface_render_method = 'BLENDED' if hasattr(m, 'surface_render_method') else None
        m.blend_method = 'BLEND' if hasattr(m, 'blend_method') else None
        m.show_transparent_back = False
    m.diffuse_color = (*color, alpha)
    return m

def import_mesh(path):
    before = set(bpy.data.objects)
    ext = os.path.splitext(path)[1].lower()
    if ext == '.ply': bpy.ops.wm.ply_import(filepath=path)
    elif ext == '.obj': bpy.ops.wm.obj_import(filepath=path, forward_axis='Y', up_axis='Z')
    elif ext in ('.glb', '.gltf'): bpy.ops.import_scene.gltf(filepath=path)
    elif ext == '.stl': bpy.ops.wm.stl_import(filepath=path)
    new = [o for o in set(bpy.data.objects) - before if o.type == 'MESH']
    return new[0] if new else None

def build_stage(stage):
    objs = {}
    for s in stage['structures']:
        path = os.path.join(base, s['file'])
        if not os.path.exists(path): print('MISSING', path); continue
        o = import_mesh(path)
        if o is None: continue
        o.name = f"{s['name']}"
        for c in list(o.users_collection): c.objects.unlink(o)
        get_coll(s.get('collection', 'Divers'), get_coll(stage['stage'])).objects.link(o)
        mat = make_material(f"{s['name']}", s.get('color', [0.8,0.7,0.6]), s.get('alpha', 1.0), sss=s.get('sss', 0.0))
        o.data.materials.clear(); o.data.materials.append(mat)
        for p in o.data.polygons: p.use_smooth = True
        o['structure'] = s['name']; o['stage'] = stage['stage']
        objs[s['name']] = o
    return objs

stages = M.get('stages') or [M]
built = [build_stage(st) for st in stages]

# ---------------- morphing entre stades
# structures "blob" morphables (shape keys sur la topologie du dernier stade, projetées par Shrinkwrap après pré-alignement),
# les autres (multi-morceaux) sont simplement affichées quand le slider est proche de leur stade.
MORPH = {'enveloppe', 'foie', 'coeur', 'cavite_pericardique', 'yeux', 'coeur_detoure', 'myocarde', 'digestif_estomac', 'digestif_oesophage', 'notochorde'}
# 24/09 16:00 : subdivisions du SNC (prosencephale, mesencephale, rhombencephale, moelle, ventricule_*) et meninges retirées du morphing :
# volumes incohérents d'un stade à l'autre (prosencéphale CS15 0,06 mm³ vs 2,0 à CS16 ; moelle CS16 0,19 vs 13,8 à CS17) → cerveau « cassé » ;
# bloc méningé à coupe plane ; 'ventricules' global retiré aussi (8,9 mm³ à CS16, 2,1 à CS17, 80 à CS19 : sort du cerveau).
# 17:00 : SNC MASQUÉ du morphing jusqu'à la re-segmentation de CS15/CS16 (consigne utilisateur) ; il reviendra en deux pièces encephale_morph /
# moelle_morph (découpe du snc global à C1 par axial_recrutement, labels déjà produits) : les remettre dans MORPH à ce moment. Membres retirés du morphing ;
# cordon_ombilical retiré aussi (17:40 : nappe plate aux jeunes stades, volumes erratiques 0,22 → 8,9 → 5,8 mm³, plaques dans le morph)
# (« tu les enlèves carrément ») ; le découpage fin reste dans les scènes par stade.
def bbox_of(obj):
    import mathutils
    pts = [obj.matrix_world @ mathutils.Vector(c) for c in obj.bound_box]
    mn = mathutils.Vector((min(p.x for p in pts), min(p.y for p in pts), min(p.z for p in pts)))
    mx = mathutils.Vector((max(p.x for p in pts), max(p.y for p in pts), max(p.z for p in pts)))
    return mn, mx
def add_stage_driver(id_data, path, expr):
    fc = id_data.driver_add(path); drv = fc.driver; drv.type = 'SCRIPTED'
    var = drv.variables.new(); var.name = 'st'; var.type = 'SINGLE_PROP'; var.targets[0].id_type = 'SCENE'; var.targets[0].id = scene; var.targets[0].data_path = '["stage"]'
    drv.expression = expr
if len(stages) > 1:
    n = len(stages)
    scene['stage'] = float(n-1)
    ui = scene.id_properties_ui('stage'); ui.update(min=0.0, max=float(n-1), description='Stade : ' + ', '.join('%d = %s' % (i, st['stage']) for i, st in enumerate(stages)))
    def project(src_coords, ref_obj, tgt):
        """coordonnées de ref_obj (topologie de référence) posées sur la surface de tgt, en partant de src_coords (stade voisin déjà calculé)"""
        tmp = ref_obj.copy(); tmp.data = ref_obj.data.copy(); tmp.shape_key_clear(); tmp.animation_data_clear(); scene.collection.objects.link(tmp)
        for v, c in zip(tmp.data.vertices, src_coords): v.co = c
        import mathutils
        smn = mathutils.Vector([min(c[i] for c in src_coords) for i in range(3)]); smx = mathutils.Vector([max(c[i] for c in src_coords) for i in range(3)])   # boîte des coordonnées réelles (bound_box de tmp serait celle de la référence, non rafraîchie)
        tmn, tmx = bbox_of(tgt)
        sc = (smn+smx)/2; tc = (tmn+tmx)/2; se = smx-smn; te = tmx-tmn
        scale = [min(max((te[i] / se[i]) if se[i] > 1e-6 else 1.0, 0.3), 3.0) for i in range(3)]
        for v in tmp.data.vertices:
            for i in range(3): v.co[i] = tc[i] + (v.co[i] - sc[i]) * scale[i]
        tmp.data.update()
        sw = tmp.modifiers.new('sw', 'SHRINKWRAP'); sw.target = tgt; sw.wrap_method = 'NEAREST_SURFACEPOINT'; sw.wrap_mode = 'ON_SURFACE'
        cs = tmp.modifiers.new('cs', 'CORRECTIVE_SMOOTH'); cs.iterations = 30; cs.factor = 0.7; cs.use_only_smooth = True; cs.smooth_type = 'LENGTH_WEIGHTED'   # paramètres du matin (validés par l'utilisateur)
        dg = bpy.context.evaluated_depsgraph_get(); ev = tmp.evaluated_get(dg)
        coords = [v.co.copy() for v in ev.data.vertices]
        bpy.data.objects.remove(tmp, do_unlink=True)
        return lisser_etirements(coords, ref_obj)
    _adj_cache = {}
    def lisser_etirements(coords, ref_obj, seuil=2.5, iters=8):
        """anti-pointes : les arêtes étirées à plus de `seuil` fois leur allongement médian (par rapport à la topologie de référence)
        voient leurs sommets remplacés par la moyenne de leurs voisins (l'utilisateur refuse les pointes / faces en éventail)"""
        import mathutils
        me = ref_obj.data
        if me.name not in _adj_cache:
            edges = [(e.vertices[0], e.vertices[1]) for e in me.edges]; adj = [[] for _ in range(len(me.vertices))]
            for a, b in edges: adj[a].append(b); adj[b].append(a)
            base_len = [(me.vertices[a].co - me.vertices[b].co).length for a, b in edges]
            _adj_cache[me.name] = (edges, adj, base_len)
        edges, adj, base_len = _adj_cache[me.name]
        n_bad = 0; bad = None
        for _ in range(iters):
            if bad is None:                                   # sommets fautifs déterminés UNE fois (sinon le lissage se propage de proche en proche)
                ratios = [((coords[a] - coords[b]).length / bl) if bl > 1e-9 else 1.0 for (a, b), bl in zip(edges, base_len)]
                srt = sorted(ratios); med = srt[len(srt) // 2] or 1.0
                lo_ = mathutils.Vector([min(c[i] for c in coords) for i in range(3)]); hi_ = mathutils.Vector([max(c[i] for c in coords) for i in range(3)]); diag = (hi_ - lo_).length
                bad = set()
                for (a, b), r in zip(edges, ratios):
                    if r > seuil * med or (coords[a] - coords[b]).length > 0.10 * diag: bad.add(a); bad.add(b)     # relatif ou absolu (filament)
            if not bad: break
            n_bad = max(n_bad, len(bad)); new = {}
            for v in bad:
                if adj[v]: new[v] = sum((coords[u] for u in adj[v]), mathutils.Vector()) / len(adj[v])
            for v, c in new.items(): coords[v] = c
        if n_bad: print('  lissage étirements %s : %d sommets' % (ref_obj.name, n_bad))
        return coords
    all_names = sorted({nm for b in built for nm in b}); morph_web = []
    for name in all_names:
        have = [i for i, b in enumerate(built) if name in b]
        ref_idx = have[-1]; obj = built[ref_idx][name]
        animable = name in MORPH and len(have) == n          # consigne utilisateur : animé = présent aux 7 stades ; sinon pas montré dans le morphing
        if animable:
            # décimation pour le web (agrégateur : ≤ 40 k faces, enveloppe ≤ 60 k) : la topologie de référence est décimée AVANT les clés,
            # la scène Blender et l'export morph_web partagent donc exactement les mêmes maillages
            fmax = 60000 if name == 'enveloppe' else 40000
            if len(obj.data.polygons) > fmax:
                dm = obj.modifiers.new('dec', 'DECIMATE'); dm.ratio = fmax / len(obj.data.polygons); dm.use_collapse_triangulate = True
                dg_ = bpy.context.evaluated_depsgraph_get(); me2 = bpy.data.meshes.new_from_object(obj.evaluated_get(dg_)); me2.name = obj.data.name + '_web'
                old_me = obj.data; obj.modifiers.remove(dm); obj.data = me2; bpy.data.meshes.remove(old_me)
                for pg in obj.data.polygons: pg.use_smooth = True
            # projection DIRECTE de la topologie de référence (dernier stade) sur chaque stade, comme le matin (l'utilisateur veut « comme avant » ;
            # la projection en chaîne stade à stade, essayée à 15:25, est abandonnée : CHAINE = True pour la réessayer)
            CHAINE = False
            obj.shape_key_add(name='Basis', from_mix=False)
            base = [v.co.copy() for v in obj.data.vertices]
            coords = {ref_idx: base}
            for order in (range(ref_idx - 1, -1, -1), range(ref_idx + 1, n)):
                prev = base
                for si in order:
                    tgt = built[si].get(name)
                    coords[si] = project(prev if CHAINE else base, obj, tgt) if tgt is not None else prev     # stade absent : forme du voisin (pas de retour à la base)
                    prev = coords[si]
            for si in range(n):
                if si == ref_idx: continue
                sk = obj.shape_key_add(name=stages[si]['stage'], from_mix=False)
                for v, c in zip(sk.data, coords[si]): v.co = c
                add_stage_driver(sk, 'value', 'max(0, 1 - abs(st - %d))' % si)
            morph_web.append((name, obj, [coords[si] for si in range(n)]))
            missing = [si for si in range(n) if si not in have]
            if missing:                                                             # caché quand on est plus près d'un stade où la structure n'existe pas
                expr = '(' + ' + '.join('(abs(st - %d) < 0.5)' % m for m in missing) + ') > 0'
                add_stage_driver(obj, 'hide_render', expr); add_stage_driver(obj, 'hide_viewport', expr)
            for i in have:
                if i != ref_idx: built[i][name].hide_render = True; built[i][name].hide_viewport = True
        else:
            # non animé : pas montré dans la scène de morphing (consigne utilisateur), mais conservé dans le fichier
            for i in have:
                o = built[i][name]; o.hide_render = True; o.hide_viewport = True; o['non_anime'] = True
    # ---------------- vertèbres morphées niveau par niveau (demande utilisateur : nombre de niveaux imposé, identique aux 7 stades)
    # actif seulement si chaque stade a un out/vertebres/manifest_vertebres.json dont les noms de niveaux sont identiques d'un stade à l'autre
    # et ne sont pas des étages manuels « etage k » / « queue k ». Corps seuls (entrées corps_<nom>) si présents, sinon niveau complet (vertebre_<nom>).
    import re as _re
    mdir = os.path.dirname(os.path.abspath(manifest_path)); vert_sets = []
    for st in stages:
        f0 = next((x['file'] for x in st['structures'] if x.get('file')), None)
        vdir = os.path.join(os.path.dirname(os.path.join(mdir, f0)), 'vertebres') if f0 else None
        pm = os.path.join(vdir, 'manifest_vertebres.json') if vdir else None
        if not pm or not os.path.exists(pm): vert_sets.append(None); continue
        vm = json.load(open(pm, encoding='utf-8'))['structures']
        pref = 'corps_' if any(x['name'].startswith('corps_') for x in vm) else 'vertebre_'
        vert_sets.append({x['name'][len(pref):]: os.path.join(vdir, x['file']) for x in vm if x['name'].startswith(pref)})
    names_ok = all(v is not None for v in vert_sets) and len({frozenset(v) for v in vert_sets}) == 1 and vert_sets and \
               not any(_re.match(r'^(etage|queue)[ _]', k) for k in vert_sets[0]) and 15 <= len(vert_sets[0]) <= 45
    if names_ok:
        coll_v = get_coll('Vertèbres (morph)', get_coll(stages[-1]['stage'])); mat_v = make_material('vertebre_morph', (0.9, 0.9, 0.85), 1.0)
        for lvl in sorted(vert_sets[0]):
            if _re.match(r'^occ', lvl): continue                     # consigne utilisateur : pas de corps vertébraux au niveau occipital (33 vertèbres C1–Co4, 37 somites)
            objs_l = []
            for si in range(n):
                o = import_mesh(vert_sets[si][lvl])
                if o is None: objs_l.append(None); continue
                o.name = 'vert_%s_%s' % (lvl, stages[si]['stage']); objs_l.append(o)
            if objs_l[-1] is None or sum(o is not None for o in objs_l) < 2:
                for o in objs_l:
                    if o is not None: bpy.data.objects.remove(o, do_unlink=True)
                continue
            ref = objs_l[-1]; ref.name = 'vertebre_' + lvl; ref['structure'] = 'vertebre_' + lvl; ref['stage'] = 'morph'
            for c in list(ref.users_collection): c.objects.unlink(ref)
            coll_v.objects.link(ref); ref.data.materials.clear(); ref.data.materials.append(mat_v)
            ref.shape_key_add(name='Basis', from_mix=False); base = [v.co.copy() for v in ref.data.vertices]
            prev = base; coords = {}
            for si in range(n - 2, -1, -1):
                coords[si] = project(prev, ref, objs_l[si]) if objs_l[si] is not None else prev; prev = coords[si]
            for si in range(n - 1):
                sk = ref.shape_key_add(name=stages[si]['stage'], from_mix=False)
                for v, c in zip(sk.data, coords[si]): v.co = c
                add_stage_driver(sk, 'value', 'max(0, 1 - abs(st - %d))' % si)
            morph_web.append(('vertebre_' + lvl, ref, [coords[si] if si < n - 1 else base for si in range(n)]))
            missing = [si for si in range(n) if objs_l[si] is None]
            if missing:
                expr = '(' + ' + '.join('(abs(st - %d) < 0.5)' % m for m in missing) + ') > 0'
                add_stage_driver(ref, 'hide_render', expr); add_stage_driver(ref, 'hide_viewport', expr)
            for o in objs_l[:-1]:
                if o is not None: bpy.data.objects.remove(o, do_unlink=True)
        print('VERTEBRES morphées :', len([k for k in vert_sets[0] if not _re.match(r'^occ', k)]), 'niveaux (occipitaux exclus)')
    else:
        print('VERTEBRES non morphées (nommage des niveaux non commun aux stades)')
    # ---------------- tubes morphables (session VHE) : pharynx, tube digestif et aortes à topologie commune, une forme par stade (embryons_3D/tubes_morph.npz)
    tubes_npz = os.path.join(os.path.dirname(os.path.abspath(out_blend)), 'tubes_morph.npz')
    if not os.path.exists(tubes_npz): tubes_npz = os.path.join(os.path.dirname(os.path.abspath(manifest_path)), 'tubes_morph.npz')
    TUBES_VHE = False   # 24/09 : l'utilisateur ne veut pas des tubes morphés (digestif + aortes) ; bloc conservé mais désactivé
    if TUBES_VHE and os.path.exists(tubes_npz):
        import numpy as np
        tz = np.load(tubes_npz); tj_path = os.path.splitext(tubes_npz)[0] + '.json'
        tj = json.load(open(tj_path, encoding='utf-8')) if os.path.exists(tj_path) else {}
        t_stades = tj.get('stades', [st['stage'] for st in stages]); idx_of = {st['stage']: i for i, st in enumerate(stages)}
        coll_t = get_coll('Tubes morphables', get_coll(stages[-1]['stage']))
        col_tube = {'aorte': (0.9, 0.15, 0.15)}
        for key in tz.files:
            A = tz[key]                                                   # (n_stades, N anneaux, M sommets, 3) en mm, repère de meshexport
            if tj.get('structures', {}).get(key, {}).get('type', 'tube') != 'tube': continue   # nappes de mésos (N × K) : voir tubes_morph_blender.py
            if A.ndim != 4 or A.shape[0] != len(t_stades): print('tube %s : forme inattendue %s' % (key, A.shape)); continue
            S_, N_, M_, _ = A.shape
            def tube_coords(k):
                ring = A[k].reshape(-1, 3); caps = np.stack([A[k][0].mean(0), A[k][-1].mean(0)])
                return [tuple(map(float, c)) for c in np.concatenate([ring, caps])]
            faces = []
            for r in range(N_ - 1):
                for m in range(M_):
                    m2 = (m + 1) % M_; faces.append((r*M_ + m, r*M_ + m2, (r+1)*M_ + m2, (r+1)*M_ + m))
            c0, c1 = N_*M_, N_*M_ + 1
            for m in range(M_):
                m2 = (m + 1) % M_; faces.append((c0, m2, m)); faces.append((c1, (N_-1)*M_ + m, (N_-1)*M_ + m2))
            ref_k = S_ - 1
            me = bpy.data.meshes.new('tube_' + key); me.from_pydata(tube_coords(ref_k), [], faces); me.update()
            for pg in me.polygons: pg.use_smooth = True
            o = bpy.data.objects.new('tube_' + key, me); coll_t.objects.link(o)
            fam = 'aorte' if key.startswith('aorte') else 'digestif'
            o.data.materials.append(make_material('tube_' + fam, col_tube.get(fam, (0.95, 0.6, 0.2)), 1.0))
            o['structure'] = 'tube_' + key; o['stage'] = 'morph'; o['source'] = 'VHE tubes_morph.npz'
            o.shape_key_add(name='Basis', from_mix=False)
            for k, stn in enumerate(t_stades):
                if k == ref_k or stn not in idx_of: continue
                sk = o.shape_key_add(name=stn, from_mix=False)
                for v, c in zip(sk.data, tube_coords(k)): v.co = c
                add_stage_driver(sk, 'value', 'max(0, 1 - abs(st - %d))' % idx_of[stn])
        print('TUBES', len(tz.files), 'tubes morphables depuis', tubes_npz)
    # ---------------- export des clés de morphing pour le viewer web (agrégateur) : embryons_3D/morph_web/<structure>.npz + morph_web.json
    try:
        import numpy as np
        wdir = os.path.join(os.path.dirname(os.path.abspath(out_blend)), 'morph_web'); os.makedirs(wdir, exist_ok=True)
        meta = {'stades': [st['stage'] for st in stages], 'unites': 'mm', 'repere': 'pipeline (meshexport : X=LR, Y=AP, Z=-SI)', 'faces_max': 60000, 'structures': []}
        for name, obj, keys_ in morph_web:
            me = obj.data; me.calc_loop_triangles()
            faces = np.array([[t.vertices[0], t.vertices[1], t.vertices[2]] for t in me.loop_triangles], np.int32)
            pos = np.array([[[c.x, c.y, c.z] for c in k] for k in keys_], np.float32)
            np.savez_compressed(os.path.join(wdir, name + '.npz'), positions=pos, faces=faces)
            mat = obj.data.materials[0] if obj.data.materials else None
            col = list(mat.diffuse_color) if mat else [0.7, 0.7, 0.7, 1.0]
            meta['structures'].append({'nom': name, 'fichier': name + '.npz', 'sommets': int(pos.shape[1]), 'faces': int(len(faces)), 'couleur': [round(float(c), 3) for c in col[:3]], 'alpha': round(float(col[3]), 3)})
        json.dump(meta, open(os.path.join(wdir, 'morph_web.json'), 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
        print('MORPH_WEB', len(morph_web), 'structures ->', wdir)
    except Exception as e: print('MORPH_WEB échec :', e)
    # animation : Lecture fait défiler les stades (demande utilisateur). 48 images par transition à 24 i/s (2 s), pause de 12 images sur chaque stade.
    FPS, TRANS, PAUSE = 24, 48, 12
    scene.render.fps = FPS; f = 1
    for si in range(n):
        scene['stage'] = float(si); scene.keyframe_insert(data_path='["stage"]', frame=f)
        scene.keyframe_insert(data_path='["stage"]', frame=f + PAUSE)          # palier
        f += PAUSE + TRANS
    scene.frame_start = 1; scene.frame_end = f - TRANS; scene.frame_current = 1
    # taille NORMALISÉE, caméra fixe (consigne utilisateur) : un parent vide « Embryon » porte toutes les pièces ; son échelle (taille CS20 / taille du stade)
    # et sa position (centre du corps ramené à l'origine) sont animées aux mêmes images que « stage ». Tous les stades remplissent donc le cadre.
    import mathutils
    root = bpy.data.objects.new('Embryon', None); scene.collection.objects.link(root); root.empty_display_size = 1.0
    for o in scene.objects:
        if o.type == 'MESH' and o.parent is None: o.parent = root
    def emprise(frame):
        scene.frame_set(frame); dg_ = bpy.context.evaluated_depsgraph_get()
        mn = mathutils.Vector((1e9,)*3); mx = mathutils.Vector((-1e9,)*3)
        for o in scene.objects:
            if o.type != 'MESH' or o.hide_render or o.evaluated_get(dg_).hide_render: continue
            e = o.evaluated_get(dg_)
            for c in e.bound_box:
                w = e.matrix_world @ mathutils.Vector(c); mn = mathutils.Vector(map(min, mn, w)); mx = mathutils.Vector(map(max, mx, w))
        return (mn + mx) / 2, max(mx - mn)
    root.scale = (1, 1, 1); root.location = (0, 0, 0)
    stage_frames = []; f = 1
    for si in range(n): stage_frames.append(f); f += PAUSE + TRANS
    mesures = [emprise(f0) for f0 in stage_frames]                      # toutes les mesures AVANT la première clé (le parent est encore à l'identité)
    ref_center, ref_size = mesures[n - 1]                               # référence : le dernier stade (CS20), facteur 1
    for si, f0 in enumerate(stage_frames):
        c, sz = mesures[si]; k = (ref_size / sz) if sz > 1e-6 else 1.0
        root.scale = (k, k, k); root.location = -c * k                  # le centre du corps du stade reste à l'origine
        for fr in (f0, f0 + PAUSE): root.keyframe_insert('scale', frame=fr); root.keyframe_insert('location', frame=fr)
        print('  normalisation %s : taille %.1f mm, facteur %.2f' % (stages[si]['stage'], sz, k))
    def fcurves_of(idb):
        ad_ = idb.animation_data
        try: return list(ad_.action.fcurves)
        except Exception: return list(ad_.action.layers[0].strips[0].channelbag(ad_.action_slot).fcurves)
    try:
        for fc in fcurves_of(root):
            for kp in fc.keyframe_points: kp.interpolation = 'BEZIER'; kp.handle_left_type = kp.handle_right_type = 'AUTO_CLAMPED'
    except Exception as e: print('clés parent :', e)
    fcs = []
    ad = scene.animation_data
    try: fcs = list(ad.action.fcurves)
    except Exception:
        try: fcs = list(ad.action.layers[0].strips[0].channelbag(ad.action_slot).fcurves)
        except Exception as e: print('fcurves introuvables :', e)
    for fc in fcs:
        for kp in fc.keyframe_points: kp.interpolation = 'BEZIER'; kp.handle_left_type = kp.handle_right_type = 'AUTO_CLAMPED'   # paliers plats, transitions en S
    scene.frame_set(1)
    print('ANIM stage : %d clés, images 1-%d (%d i/s)' % (2 * n, scene.frame_end, FPS))

# brouillons masqués par défaut (consigne utilisateur : pas de SNC grossier tant que la segmentation fine n'est pas là)
BROUILLONS = {'snc', 'ganglions', 'cavite_pericardique'}
for b in built:
    fine = any(k in b for k in ('prosencephale', 'moelle'))
    for nm, o in b.items():
        if nm in BROUILLONS and fine and not (len(stages) > 1 and nm == 'snc'):      # scène maître : le snc global est la seule pièce SNC morphée
            o.hide_render = True; o.hide_viewport = True; o['brouillon'] = True
# caméra + lumières
cam = bpy.data.cameras.new('Camera'); camo = bpy.data.objects.new('Camera', cam); scene.collection.objects.link(camo); scene.camera = camo
camo.location = (0, -60, 0); camo.rotation_euler = (math.radians(90), 0, 0); cam.lens = 85
for i,(loc, e) in enumerate([((30,-40,40), 2000), ((-40,-20,20), 800), ((0,40,30), 600)]):
    l = bpy.data.lights.new(f'Light{i}', 'AREA'); l.energy = e; l.size = 30
    lo = bpy.data.objects.new(f'Light{i}', l); scene.collection.objects.link(lo); lo.location = loc
    lo.rotation_euler = (bpy.data.objects['Camera'].location*0 - lo.location*0).to_track_quat('Z','Y').to_euler() if False else (0,0,0)
    import mathutils
    lo.rotation_euler = (mathutils.Vector((0,0,0)) - lo.location).to_track_quat('-Z','Y').to_euler()
scene.render.engine = 'CYCLES'; scene.cycles.samples = 128; scene.cycles.device = 'GPU'
scene.render.resolution_x = 1920; scene.render.resolution_y = 1080
bpy.ops.wm.save_as_mainfile(filepath=out_blend)
print('SAVED', out_blend)
