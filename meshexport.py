"""Masques de labels -> maillages lissés/décimés (mm, axes Blender X=gauche-droite, Y=antéro-postérieur, Z=haut) + manifest.json + NRRD."""
import numpy as np, trimesh, json, os, time
from scipy import ndimage as ndi
from skimage import measure
import fast_simplification

# longueur vertex-coccyx (greatest length) typique par stade, mm  (O'Rahilly & Müller) -> échelle absolue approx.
CRL_MM = {'CS13': 4.5, 'CS14': 6.0, 'CS15': 7.5, 'CS16': 9.0, 'CS17': 11.0, 'CS18': 14.0, 'CS19': 17.0, 'CS20': 20.0, 'CS21': 22.0, 'CS22': 25.0, 'CS23': 28.0}

# confiance par défaut par structure (bonne/moyenne/faible) ; surchargeable par stade via work/confiance.json {label: niveau}
CONFIANCE = {'enveloppe': 'bonne', 'snc': 'faible', 'ventricules': 'bonne', 'ganglions': 'faible', 'yeux': 'moyenne', 'cristallins': 'moyenne',
 'vesicules_otiques': 'moyenne', 'coeur': 'moyenne', 'cavite_pericardique': 'faible', 'foie': 'bonne', 'tube_digestif': 'faible',
 'squelette_axial_cartilage': 'moyenne', 'chondrocrane': 'faible', 'cartilage_autre': 'faible', 'membres': 'bonne', 'encephale_morph': 'moyenne', 'moelle_morph': 'moyenne', 'membre_sup_gauche': 'bonne', 'membre_sup_droit': 'bonne', 'membre_inf_gauche': 'moyenne', 'membre_inf_droit': 'moyenne', 'cordon_ombilical': 'moyenne',
 'vaisseaux': 'faible', 'coeur_detoure': 'moyenne', 'myocarde': 'moyenne', 'cavites_cardiaques': 'faible', 'digestif_pharynx': 'moyenne', 'digestif_oesophage': 'moyenne',
 'digestif_estomac': 'moyenne', 'digestif_duodenum': 'moyenne', 'digestif_intestin_moyen': 'moyenne', 'digestif_intestin_posterieur': 'moyenne',
 'vaisseaux_aorte': 'moyenne', 'vaisseaux_arcs_aortiques': 'moyenne', 'digestif_poches_pharyngiennes': 'moyenne', 'vaisseaux_cardinales': 'faible', 'vaisseaux_ombilicaux': 'faible', 'vaisseaux_vitellins': 'faible',
 'somites': 'moyenne', 'notochorde': 'faible', 'corps_vertebraux': 'moyenne', 'arcs_neuraux': 'moyenne', 'cotes': 'faible',
 'ventricule_prosencephale': 'bonne', 'ventricule_mesencephale': 'moyenne', 'ventricule_rhombencephale': 'bonne', 'canal_central': 'moyenne',
 'prosencephale': 'moyenne', 'mesencephale': 'moyenne', 'rhombencephale': 'moyenne', 'moelle': 'moyenne', 'moelle_rachidienne': 'moyenne',
 'meninges_mesenchyme_cranien': 'faible', 'epiderme_cranien': 'moyenne'}
# nom -> (collection, couleur RGB 0-1, alpha, faces cible, sigma lissage)
STYLE = {
 'enveloppe':                 ('Enveloppe',   (0.93, 0.78, 0.70), 0.25, 150000, 1.5),
 'snc':                       ('Systeme nerveux', (1.00, 0.72, 0.30), 1.0, 80000, 1.2),
 'ventricules':               ('Systeme nerveux', (0.55, 0.85, 1.00), 0.6, 80000, 1.2),
 'ganglions':                 ('Systeme nerveux', (1.00, 0.55, 0.55), 1.0, 60000, 1.0),
 'yeux':                      ('Organes des sens', (0.25, 0.35, 0.90), 1.0, 40000, 1.0),
 'cristallins':               ('Organes des sens', (0.95, 0.95, 0.60), 1.0, 10000, 0.8),
 'vesicules_otiques':         ('Organes des sens', (0.30, 0.90, 0.90), 1.0, 20000, 1.0),
 'coeur':                     ('Organes',     (0.80, 0.15, 0.15), 1.0, 150000, 1.2),
 'cavite_pericardique':       ('Cavites',     (0.75, 0.85, 1.00), 0.3, 80000, 1.5),
 'foie':                      ('Organes',     (0.55, 0.20, 0.30), 1.0, 80000, 1.2),
 'tube_digestif':             ('Organes',     (0.95, 0.60, 0.20), 1.0, 100000, 1.0),
 'squelette_axial_cartilage': ('Squelette',   (0.92, 0.92, 0.85), 1.0, 150000, 0.9),
 'chondrocrane':              ('Squelette',   (0.90, 0.90, 0.80), 1.0, 80000, 0.9),
 'cartilage_autre':           ('Squelette',   (0.85, 0.85, 0.80), 1.0, 80000, 0.9),
 'membres':                   ('Enveloppe',   (0.90, 0.75, 0.65), 1.0, 80000, 1.2),
 'encephale_morph':            ('SNC',         (0.85, 0.75, 0.35), 1.0, 60000, 1.2),
 'moelle_morph':               ('SNC',         (0.80, 0.70, 0.30), 1.0, 40000, 1.2),
 'membre_sup_gauche':         ('Membres',     (0.90, 0.75, 0.65), 1.0, 30000, 1.2),
 'membre_sup_droit':          ('Membres',     (0.90, 0.75, 0.65), 1.0, 30000, 1.2),
 'membre_inf_gauche':         ('Membres',     (0.90, 0.75, 0.65), 1.0, 30000, 1.2),
 'membre_inf_droit':          ('Membres',     (0.90, 0.75, 0.65), 1.0, 30000, 1.2),
 'cordon_ombilical':          ('Enveloppe',   (0.80, 0.70, 0.60), 1.0, 30000, 1.2),
 'vaisseaux':                 ('Vaisseaux',   (0.85, 0.10, 0.10), 1.0, 150000, 0.8),
 'coeur_detoure':             ('Cardio',      (0.75, 0.12, 0.12), 1.0, 80000, 1.0),
 'myocarde':                  ('Cardio',      (0.85, 0.30, 0.25), 1.0, 80000, 1.0),
 'cavites_cardiaques':        ('Cardio',      (0.95, 0.60, 0.60), 0.7, 100000, 1.0),
 'digestif_pharynx':          ('Digestif',    (0.90, 0.55, 0.45), 1.0, 60000, 1.0),
 'digestif_oesophage':        ('Digestif',    (0.95, 0.75, 0.35), 1.0, 60000, 1.0),
 'digestif_estomac':          ('Digestif',    (0.95, 0.60, 0.20), 1.0, 80000, 1.0),
 'digestif_duodenum':         ('Digestif',    (0.90, 0.55, 0.25), 1.0, 60000, 1.0),
 'digestif_intestin_moyen':   ('Digestif',    (0.85, 0.50, 0.15), 1.0, 80000, 1.0),
 'digestif_intestin_posterieur': ('Digestif', (0.80, 0.45, 0.20), 1.0, 60000, 1.0),
 'vaisseaux_aorte':           ('Vaisseaux',   (0.90, 0.15, 0.15), 1.0, 120000, 0.8),
 'vaisseaux_arcs_aortiques':  ('Vaisseaux',   (0.25, 0.70, 0.40), 1.0, 120000, 0.8),
 'digestif_poches_pharyngiennes': ('Digestif', (0.78, 0.78, 0.86), 1.0, 60000, 1.0),
 'vaisseaux_cardinales':      ('Vaisseaux',   (0.25, 0.35, 0.85), 1.0, 120000, 0.8),
 'vaisseaux_ombilicaux':      ('Vaisseaux',   (0.60, 0.20, 0.60), 1.0, 80000, 0.8),
 'vaisseaux_vitellins':       ('Vaisseaux',   (0.85, 0.50, 0.20), 1.0, 80000, 0.8),
 'ventricule_prosencephale':  ('Systeme nerveux', (1.00, 0.80, 0.35), 0.6, 80000, 1.2),
 'ventricule_mesencephale':   ('Systeme nerveux', (1.00, 0.50, 0.80), 0.6, 80000, 1.2),
 'ventricule_rhombencephale': ('Systeme nerveux', (0.50, 0.85, 1.00), 0.6, 80000, 1.2),
 'canal_central':             ('Systeme nerveux', (0.75, 1.00, 0.75), 0.7, 60000, 1.0),
 'prosencephale':             ('Systeme nerveux', (0.80, 0.50, 0.10), 1.0, 120000, 1.0),
 'mesencephale':              ('Systeme nerveux', (0.70, 0.10, 0.50), 1.0, 120000, 1.0),
 'rhombencephale':            ('Systeme nerveux', (0.10, 0.45, 0.80), 1.0, 120000, 1.0),
 'moelle':                    ('Systeme nerveux', (0.20, 0.65, 0.25), 1.0, 100000, 1.0),
 'meninges_mesenchyme_cranien': ('Enveloppe',  (0.92, 0.90, 0.80), 0.35, 80000, 1.5),
 'epiderme_cranien':          ('Enveloppe',    (0.85, 0.75, 0.70), 0.4, 60000, 1.5),
 'corps_vertebraux':          ('Squelette',   (0.93, 0.92, 0.84), 1.0, 120000, 0.8),
 'moelle_rachidienne':        ('Systeme nerveux', (0.20, 0.65, 0.25), 1.0, 80000, 1.0),
 'arcs_neuraux':              ('Squelette',   (0.88, 0.88, 0.80), 1.0, 100000, 0.8),
 'cotes':                     ('Squelette',   (0.85, 0.86, 0.78), 1.0, 120000, 0.8),
 'somites':                   ('Squelette',   (0.70, 0.80, 0.95), 1.0, 200000, 0.8),
 'notochorde':                ('Squelette',   (0.95, 0.95, 0.50), 1.0, 60000, 0.9),
 'sclerotomes':               ('Squelette',   (0.80, 0.85, 0.95), 1.0, 150000, 0.8),
}

def mask_to_mesh(mask, scale_mm, sigma=1.2, target_faces=None, min_faces_component=60, pad=3, iso=0.5):
    if not mask.any(): return None
    lo = []; hi = []
    for ax in range(3):
        idx = np.flatnonzero(mask.any(axis=tuple(a for a in range(3) if a != ax)))
        lo.append(max(int(idx[0]) - pad, 0)); hi.append(min(int(idx[-1]) + pad + 1, mask.shape[ax]))
    lo = np.array(lo); hi = np.array(hi)
    sub = mask[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]].astype(np.float32)
    if sigma > 0: sub = ndi.gaussian_filter(sub, sigma)
    if sub.max() <= iso: return None
    verts, faces, _, _ = measure.marching_cubes(sub, level=iso, step_size=1)
    m = trimesh.Trimesh(verts + lo, faces, process=True)
    if min_faces_component > 0:
        parts = [p for p in m.split(only_watertight=False) if len(p.faces) >= min_faces_component]
        if not parts: return None
        m = trimesh.util.concatenate(parts)
    if target_faces is not None and len(m.faces) > target_faces * 1.05:
        v, f = fast_simplification.simplify(m.vertices.astype(np.float32), m.faces.astype(np.int32), target_reduction=float(1.0 - target_faces/len(m.faces)))
        m = trimesh.Trimesh(v, f, process=True)
    P = m.vertices * scale_mm
    Pb = np.stack([P[:,0], P[:,2], -P[:,1]], 1)          # X=LR, Y=AP, Z=-SI (haut)
    m = trimesh.Trimesh(Pb, m.faces, process=False); m.fix_normals()
    return m

def run(work, out, stage):
    meta = json.load(open(os.path.join(work, 'fused_meta.json')))
    z = np.load(os.path.join(work, 'labels.npz')); shape = tuple(z['shape'])
    labels = {k: np.unpackbits(z[k])[:int(np.prod(shape))].reshape(shape).astype(bool) for k in z.files if k != 'shape'}
    env = labels['enveloppe']
    si_extent = np.ptp(np.where(env.any(axis=(0,2)))[0]) + 1          # en voxels
    crl = CRL_MM.get(stage, 20.0)
    scale_mm = crl / si_extent                                          # mm par voxel (axe SI ~ greatest length)
    # centre (milieu de l'enveloppe) placé à l'origine
    c = np.array(ndi.center_of_mass(env))
    manifest = {'stage': stage, 'units': 'mm', 'mm_per_voxel': scale_mm, 'voxel_step_px': meta['step'], 'greatest_length_mm_assumed': crl,
                'center_voxel': c.tolist(), 'structures': []}
    os.makedirs(out, exist_ok=True)
    pc = os.path.join(work, 'confiance.json'); conf_stage = json.load(open(pc, encoding='utf-8')) if os.path.exists(pc) else {}
    for name, mask in labels.items():
        coll, col, alpha, tf, sig = STYLE.get(name, ('Divers', (0.7,0.7,0.7), 1.0, 60000, 1.0))
        t = time.time(); m = mask_to_mesh(mask, scale_mm, sigma=sig, target_faces=tf)
        if m is None: print('  (vide)', name); continue
        cc = np.array([c[0], c[2], -c[1]]) * scale_mm; m.apply_translation(-cc)
        fn = f'{stage}_{name}.ply'; m.export(os.path.join(out, fn))
        conf = conf_stage.get(name, CONFIANCE.get(name, 'moyenne'))
        manifest['structures'].append({'name': name, 'file': fn, 'collection': coll, 'color': list(col), 'alpha': alpha, 'faces': int(len(m.faces)), 'volume_mm3': float(mask.sum()) * scale_mm**3, 'confiance': conf})
        print('  %-28s %8d faces  %7.2f mm3  %.0fs' % (name, len(m.faces), mask.sum()*scale_mm**3, time.time()-t))
    # scène GLB combinée
    sc = trimesh.Scene()
    for s in manifest['structures']:
        m = trimesh.load(os.path.join(out, s['file']), force='mesh')
        m.visual = trimesh.visual.ColorVisuals(m, face_colors=np.tile((np.array(s['color']+[s['alpha']])*255).astype(np.uint8), (len(m.faces),1)))
        sc.add_geometry(m, node_name=s['name'], geom_name=s['name'])
    sc.export(os.path.join(out, f'{stage}_all.glb'))
    json.dump(manifest, open(os.path.join(out, 'manifest.json'), 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    # NRRD pour 3D Slicer : voir export_nrrd.py (appelé par pipeline.py et import_labels.py)
    print('meshes written to', out)
