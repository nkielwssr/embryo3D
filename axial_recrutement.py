"""Recrutement des volumes vertébraux guidé par l'axe de la session « Extraction squelette » (out/topographie/axe_vertebral.json
+ etages_vertebraux.json). Pour chaque voxel proche de l'axe on calcule un repère local (s curviligne, offset dorsal d, offset latéral l),
puis on segmente dans le volume de densité :
  - 'corps_vertebraux' : blob compact médian autour de l'axe des corps (xyz_corps_mm) : pâle (cartilage, CS17+) ou dense (sclérotome, CS13-16)
  - 'arcs_neuraux'     : blobs pâles dorsolatéraux flanquant la moelle (CS17+) ; à CS13-16 = partie dorsale des somites/sclérotomes
  - 'cotes'            : prolongements latéraux au-delà des arcs (|l| grand), CS17+
  - 'notochorde'       : cordon dense fin le plus proche de l'axe des corps (|offset| < r_noto), affiné
  - séparation SNC : tout ce qui est dans moelle / snc / ganglions est exclu ; 'moelle' perd ce qui tombe dans les corps
  - instances par étage : work/vertebres_niveaux.npz (int16 : index d'étage + 1 ; noms dans vertebres_niveaux.json) + PLY par niveau dans out/vertebres/
usage : python axial_recrutement.py <dossier_stade> [--diag] [--ecrire-moelle]
"""
import numpy as np, os, sys, re, json, shutil, subprocess, cv2
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from segment import dilate_r, erode_r, open_r, close_r, comps_info, _crop, _dist
stage_dir = os.path.abspath(sys.argv[1]); diag = '--diag' in sys.argv
work = os.path.join(stage_dir, 'work'); out = os.path.join(stage_dir, 'out'); topo = os.path.join(out, 'topographie')
stage = 'CS%s' % re.search(r'CS\s*(\d+)', os.path.basename(stage_dir), re.I).group(1); num = int(stage[2:]); young = num <= 16
pman = os.path.join(topo, 'video360', 'etages_manuels.json')
prec = os.path.join(topo, 'video360', 'etages_manuels_recentres.json')     # X recentré par hauteur (session squelette) : prioritaire
if os.path.exists(prec): pman = prec; sys.argv = [a for a in sys.argv if a != '--recentrer-lr']
pimp = os.path.join(topo, 'video360', 'etages_imposes.json')                # nommage imposé, identique aux 7 stades (session squelette) : priorité absolue, mode manuel, pas de niveaux de queue
impose = os.path.exists(pimp)
if impose:
    pman = pimp; sys.argv = [a for a in sys.argv if a != '--recentrer-lr'] + (['--manuel'] if '--manuel' not in sys.argv else [])
if '--manuel' in sys.argv or (not os.path.exists(os.path.join(topo, 'axe_vertebral.json')) and os.path.exists(pman)):
    # étages marqués à la main sur la vidéo 360° (centres 3D du cou vers la queue) : l'axe est la polyligne des centres
    man = json.load(open(pman, encoding='utf-8')); vox = float(json.load(open(os.path.join(out, 'manifest.json'), encoding='utf-8'))['mm_per_voxel'])
    C = np.array([e['xyz_mm_pipeline'] for e in man['etages']]) / vox
    if '--recentrer-lr' in sys.argv:                                  # X latéral approché (vue tournée) : on force la ligne médiane du corps
        envm = np.unpackbits(np.load(os.path.join(work, 'labels.npz'))['enveloppe'])[:int(np.prod(np.load(os.path.join(work, 'labels.npz'))['shape']))].reshape(tuple(np.load(os.path.join(work, 'labels.npz'))['shape'])).astype(bool)
        lr_idx = np.where(envm.any(axis=(1, 2)))[0]; C[:, 0] = 0.5 * (lr_idx.min() + lr_idx.max()); print('  centres recentrés sur LR médian = %.0f' % C[0, 0])
    seg = np.linalg.norm(np.diff(C, axis=0), axis=1); s_c = np.concatenate([[0.0], np.cumsum(seg)]) * vox; pas = float(np.median(seg) * vox)
    S = s_c; Pc = C; Pg = C.copy(); T = np.gradient(C, axis=0)
    et = {'fiabilite': {'periode_mm': pas, 'fiable': True, 'guide': 'etages manuels video360'}, 'etages': {}}
    for k in range(len(C)):
        nm = man['etages'][k].get('nom') or 'etage %d' % (k + 1)
        et['etages'][nm] = {'s_debut_mm': s_c[k] - pas / 2, 's_fin_mm': s_c[k] + pas / 2, 'hauteur_mm': pas, 'origine': 'manuel', 'nom_avant_ancre': nm, 'corps': man['etages'][k].get('corps', True)}
    print('  axe = %d étages manuels (%s ; pas %.2f mm, arc %.2f mm)' % (len(C), os.path.basename(pman), pas, s_c[-1]))
else:
    axe = json.load(open(os.path.join(topo, 'axe_vertebral.json'), encoding='utf-8')); et = json.load(open(os.path.join(topo, 'etages_vertebraux.json'), encoding='utf-8'))
    vox = float(axe['mm_per_voxel']); A = axe['axe']
    S = np.array([p['s_mm'] for p in A]); Pc = np.array([p['xyz_corps_mm'] for p in A]) / vox; Pg = np.array([p['xyz_mm'] for p in A]) / vox; T = np.array([p['tangente'] for p in A])
T /= np.linalg.norm(T, axis=1, keepdims=True) + 1e-9
# --- étages supplémentaires de la queue (vidéo 360°, session squelette) : prolongent l'axe au-delà du dernier étage
pvid = os.path.join(topo, 'video360', 'somites_video.json'); extra_names = []
if os.path.exists(pvid) and not impose:
    vues = json.load(open(pvid, encoding='utf-8')).get('vues', [])
    vue = max(vues, key=lambda v: len(v.get('etages_supplementaires_queue', []))) if vues else None
    q = (vue or {}).get('etages_supplementaires_queue', [])
    if q:
        Q = np.array([e['xyz_mm_pipeline'] for e in q]) / vox
        # ordonner depuis l'extrémité de l'axe
        order = np.argsort(np.linalg.norm(Q - Pc[-1], axis=1)); Q = Q[order]
        s_last = S[-1]; newS, newP, newT = [], [], []
        prev = Pc[-1]; s_cur = s_last
        for k, c in enumerate(Q):
            d = c - prev; L_ = np.linalg.norm(d); t = d / (L_ + 1e-9); s_cur += L_ * vox
            newS.append(s_cur); newP.append(c); newT.append(t); prev = c
        S = np.concatenate([S, newS]); Pc = np.vstack([Pc, newP]); Pg = np.vstack([Pg, newP]); T = np.vstack([T, newT])
        # étages : frontières aux milieux entre centres, hauteur = pas médian
        pas = float(np.median(np.diff(newS))) if len(newS) > 1 else float(vue.get('periode_mm', 0.3))
        et['etages'] = dict(et['etages'])
        for k, sc in enumerate(newS):
            nm = 'queue %d' % (k + 1); extra_names.append(nm)
            et['etages'][nm] = {'s_debut_mm': sc - pas / 2, 's_fin_mm': sc + pas / 2, 'hauteur_mm': pas, 'origine': 'video360 %s' % ('cale' if q[order[k]].get('cale') else 'predit'), 'nom_avant_ancre': nm}
        print('  + %d étages de queue (vidéo 360°, vue image %s, pas %.2f mm) -> axe prolongé à %.2f mm' % (len(newS), vue.get('image'), pas, S[-1]))
# repère local : dorsal = +AP (repère canonique) orthogonalisé à la tangente ; latéral = T x D (≈ gauche-droite)
D = np.tile(np.array([[0.0, 0.0, 1.0]]), (len(T), 1)); D -= (D * T).sum(1, keepdims=True) * T; D /= np.linalg.norm(D, axis=1, keepdims=True) + 1e-9
Lat = np.cross(T, D); Lat /= np.linalg.norm(Lat, axis=1, keepdims=True) + 1e-9
z = np.load(os.path.join(work, 'labels.npz')); shape = tuple(z['shape']); n = int(np.prod(shape))
L = {k: np.unpackbits(z[k])[:n].reshape(shape).astype(bool) for k in z.files if k != 'shape'}; z.close()   # fermer : sinon os.replace de labels.npz échoue sous Windows
ds = np.load(os.path.join(work, 'ds1.npy')).astype(np.float32); env = L['enveloppe']
cns = L.get('moelle', np.zeros(shape, bool)) | L.get('snc', np.zeros(shape, bool)) | L.get('ganglions', np.zeros(shape, bool))
for k in ('prosencephale', 'mesencephale', 'rhombencephale', 'ventricules', 'canal_central'):
    if k in L: cns |= L[k]
# ------------------------------------------------------------------ repère local dans un manchon autour de l'axe des corps
period_vox = float(et['fiabilite'].get('periode_mm', 0.4)) / vox
approx_x = '--manuel' in sys.argv and not impose and not (os.path.exists(prec) and pman == prec)   # X approché seulement si les centres n'ont pas été recentrés
R = max(30.0, 3.2 * period_vox, 1.2 / vox if approx_x else 0)   # rayon du manchon (voxels) ; ≥ 1,2 mm avec des centres manuels non recentrés (X approché)
lo = np.maximum(np.floor(Pc.min(0) - R).astype(int), 0); hi = np.minimum(np.ceil(Pc.max(0) + R).astype(int) + 1, shape)
sl = tuple(slice(a, b) for a, b in zip(lo, hi)); sub_shape = tuple(b - a for a, b in zip(lo, hi))
gz, gy, gx = np.mgrid[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]; pts = np.stack([gz.ravel(), gy.ravel(), gx.ravel()], 1).astype(np.float32)
tree = cKDTree(Pc); dist, idx = tree.query(pts, k=1, distance_upper_bound=R)
near = np.isfinite(dist); idx = np.where(near, idx, 0)
rel = pts - Pc[idx]; along = (rel * T[idx]).sum(1); dors = (rel * D[idx]).sum(1); lat = (rel * Lat[idx]).sum(1); dlr = np.abs(rel[:, 0])   # écart gauche-droite brut à l'axe
s_vox = S[idx] / vox + along
def field(v, fill=np.nan):
    f = np.full(shape, fill, np.float32); f[sl] = np.where(near, v, fill).reshape(sub_shape); return f
Fdist = field(dist); Fdors = field(dors); Flat = field(lat); Fs = field(s_vox); Flr = field(dlr)
# distance à l'axe de la moelle pour situer les arcs : tube_neural.json (ligne médiane de la moelle) si présent, sinon l'axe guide
ptube = os.path.join(topo, 'tube_neural.json')
if os.path.exists(ptube):
    tj = json.load(open(ptube, encoding='utf-8')); Pt = np.array([p['xyz_mm'] for p in tj['tube']]) / vox
    print('  guide moelle : tube_neural.json (%d pts, %.0f%% dans snc)' % (len(Pt), 100 * float(tj.get('fraction_points_dans_snc', 0))))
else:
    Pt = Pg
tree_g = cKDTree(Pt); dist_g, _ = tree_g.query(pts, k=1, distance_upper_bound=R * 1.5); Fdg = field(dist_g)
manchon = np.isfinite(Fdist) & env & ~cns
print('%s : axe %d pts, période %.2f mm (%.1f vox), manchon R=%.0f vox -> %.2fM vox hors SNC' % (stage, len(S), period_vox * vox, period_vox, R, manchon.sum() / 1e6))
# ------------------------------------------------------------------ tissus
r_corps = 0.9 * period_vox
if young:
    # stades somitiques : le tissu de référence est le label 'somites' (blocs métamériques) s'il existe, sinon le tissu dense
    # label 'somites' utilisé seulement s'il couvre l'axe (≥ 50 % de la longueur en s), sinon tissu dense
    use_som = False
    if 'somites' in L and L['somites'].any():
        s_som = Fs[L['somites'] & np.isfinite(Fs)]
        couv = (np.ptp(s_som) / max(np.ptp(S) / vox, 1)) if s_som.size else 0.0
        use_som = couv >= 0.5; print('  label somites : couverture de l axe %.0f%% -> %s' % (100 * couv, 'utilisé' if use_som else 'ignoré (tissu dense)'))
    tissu = L['somites'] if use_som else (ds > 70)
    corps = manchon & tissu & (np.abs(Flat) < 1.2 * r_corps) & (Fdors > -1.5 * r_corps) & (Fdors < 0.6 * r_corps)      # sclérotome ventro-médial
    arcs = manchon & tissu & ~corps & (Fdg < 2.0 * r_corps) & (np.abs(Flat) < 2.2 * r_corps) & (Fdors > -1.2 * r_corps)   # reste du somite (dermomyotome / partie dorsale), sans la paroi latérale ni le ventre
    if not use_som: arcs &= (Flr >= 0.5 * r_corps)
    cotes = np.zeros(shape, bool)
else:
    pale = (ds < 12)                                                     # cartilage non coloré
    corps = manchon & pale & (np.abs(Flat) < r_corps) & (np.abs(Fdors) < 1.0 * r_corps)
    arcs = manchon & pale & (Flr >= 0.5 * r_corps) & (np.abs(Flat) < 2.6 * r_corps) & (Fdors > 0.0) & (Fdg < 2.6 * r_corps) & ~corps
    cotes = manchon & pale & (Flr >= 2.0 * r_corps) & (Fdors > -0.8 * r_corps) & (Fdors < 2.5 * r_corps) & ~corps
    cotes &= ~arcs
# nettoyage morphologique + composantes
def clean(m, r, minvox):
    m = open_r(m, r); lab, info = comps_info(m, minvox); return np.isin(lab, [i for i, s_, c in info])
corps = clean(corps, 2 if young else 2, 400); arcs = clean(arcs, 2, 300); cotes = clean(cotes, 2, 300) if not young else cotes
# notochorde affinée : cordon dense à moins de r_noto de l'axe des corps
r_noto = max(3.0, 0.22 * period_vox)
noto = manchon & (ds > 85) & (Fdist < r_noto)
noto = clean(noto, 1, 200); lab, info = comps_info(noto, 1); info.sort(key=lambda x: -x[1]); noto = np.isin(lab, [i for i, s_, c in info[:3]])
corps &= ~noto
# ------------------------------------------------------------------ séparation SNC : moelle perd ce qui tombe dans les corps
moelle_new = None; moelle_rach = None
if 'moelle' in L:
    moelle_new = L['moelle'] & ~corps & ~noto
    # 'moelle_rachidienne' : moelle sous la frontière occipital-1 / C-1 (s_fin_mm de « occipital 1 » si présent, sinon début de l'axe)
    s_c1 = (et['etages']['C1']['s_debut_mm'] if 'C1' in et['etages'] else et['etages'].get('occipital 1', {}).get('s_fin_mm', 0.0)) / vox
    moelle_rach = moelle_new & (Fs >= s_c1)
    print('  moelle_rachidienne : sous s=%.2f mm -> %.2fM vox (moelle %.2fM)' % (s_c1 * vox, moelle_rach.sum() / 1e6, moelle_new.sum() / 1e6))
# ------------------------------------------------------------------ niveaux sans corps vertébral (occipitaux : somites seuls, consigne utilisateur)
sans_corps = [k for k, v in et['etages'].items() if v.get('corps', True) is False]
if sans_corps:
    m_sc = np.zeros(shape, bool)
    for k in sans_corps: m_sc |= np.isfinite(Fs) & (Fs >= et['etages'][k]['s_debut_mm'] / vox) & (Fs < et['etages'][k]['s_fin_mm'] / vox)
    retire = int((corps & m_sc).sum()); corps &= ~m_sc
    print('  niveaux sans corps (%s) : %d vox de corps retirés' % (', '.join(sans_corps), retire))
# ------------------------------------------------------------------ SNC animé en deux pièces (consigne utilisateur 24/09 16:40) : encéphale / moelle,
# découpe du 'snc' global à la frontière occ1 / C1 (nommage imposé) ; les pièces fines ne sont pas touchées
enc_m = moelle_m = None
if impose and 'snc' in L and 'C1' in et['etages']:
    s_c1i = et['etages']['C1']['s_debut_mm'] / vox
    snc_ = L['snc']; moelle_m = snc_ & np.isfinite(Fs) & (Fs >= s_c1i); enc_m = snc_ & ~moelle_m
    lab_e, info_e = comps_info(enc_m, 1); info_e.sort(key=lambda x: -x[1])
    if len(info_e) > 1:                                   # tout ce qui n'est pas la plus grosse composante crâniale (queue hors manchon, miettes) va à la moelle
        keep = info_e[0][0]; moelle_m |= enc_m & (lab_e != keep); enc_m = lab_e == keep
    print('  SNC animé : encephale_morph %.2fM vox, moelle_morph %.2fM vox (frontière C1 à s=%.2f mm)' % (enc_m.sum() / 1e6, moelle_m.sum() / 1e6, s_c1i * vox))
# ------------------------------------------------------------------ instances par étage
names = list(et['etages'].keys()); bounds = [(et['etages'][k]['s_debut_mm'] / vox, et['etages'][k]['s_fin_mm'] / vox) for k in names]
inst = np.zeros(shape, np.int16); squel = corps | arcs | cotes
for k, (a, b) in enumerate(bounds):
    inst[squel & (Fs >= a) & (Fs < b)] = k + 1
for k, nm in enumerate(names):
    v = int((inst == k + 1).sum()); print('  %-14s s %.2f-%.2f mm  %6d vox' % (nm, bounds[k][0] * vox, bounds[k][1] * vox, v))
print('  corps %.2fM  arcs %.2fM  cotes %.2fM  notochorde %.3fM' % (corps.sum() / 1e6, arcs.sum() / 1e6, cotes.sum() / 1e6, noto.sum() / 1e6))
# ------------------------------------------------------------------ planche
d = np.load(os.path.join(work, 'dens.npy'), mmap_mode='r'); lr = np.where(env.any(axis=(1, 2)))[0]; MID = int(np.median(Pc[:, 0])); LRe = lr.max() - lr.min()   # coupe au LR médian de l'axe
def paint(i):
    img = cv2.cvtColor(255 - np.asarray(d[i]), cv2.COLOR_GRAY2BGR)
    for m, col in ((corps, (255, 200, 200)), (arcs, (200, 200, 255)), (cotes, (200, 255, 200)), (noto, (0, 220, 220))):
        img[m[i]] = col[::-1]
    if moelle_new is not None: img[(moelle_new[i]) & ~corps[i]] = (0, 140, 0)[::-1]
    return img
mid_inst = inst[MID]; img_inst = cv2.cvtColor(255 - np.asarray(d[MID]), cv2.COLOR_GRAY2BGR)
rng = np.random.RandomState(3); cols = rng.randint(60, 255, (len(names) + 1, 3)); img_inst[mid_inst > 0] = cols[mid_inst[mid_inst > 0]]
cv2.imwrite(os.path.join(work, 'check_recrutement.png'), np.hstack([paint(MID), paint(int(MID + 0.08 * LRe)), img_inst]))
if not diag:
    L['corps_vertebraux'] = corps; L['arcs_neuraux'] = arcs; L['notochorde'] = noto
    if not young: L['cotes'] = cotes
    L['squelette_axial_cartilage'] = corps | arcs | cotes
    if moelle_new is not None and '--ecrire-moelle' in sys.argv: L['moelle'] = moelle_new
    if moelle_rach is not None and 'moelle_rachidienne' not in L: L['moelle_rachidienne'] = moelle_rach   # jamais réécrite si elle existe : Extraction squelette s'en sert de guide (--rachidienne pour forcer)
    elif moelle_rach is not None and '--rachidienne' in sys.argv: L['moelle_rachidienne'] = moelle_rach
    if enc_m is not None: L['encephale_morph'] = enc_m; L['moelle_morph'] = moelle_m
    shutil.copy(os.path.join(work, 'labels.npz'), os.path.join(work, 'labels_avant_recrutement.npz'))
    tmp = os.path.join(work, 'labels.tmp.npz'); np.savez_compressed(tmp, **{kk: np.packbits(v) for kk, v in L.items()}, shape=np.array(shape)); os.replace(tmp, os.path.join(work, 'labels.npz'))   # atomique
    inst_corps = np.where(corps, inst, 0).astype(np.int16)                      # chaîne des corps seule (sans arcs ni côtes) pour l'appariement d'axe
    tmp = os.path.join(work, 'vertebres_niveaux.tmp.npz'); np.savez_compressed(tmp, niveaux=inst, niveaux_corps=inst_corps, shape=np.array(shape)); os.replace(tmp, os.path.join(work, 'vertebres_niveaux.npz'))
    json.dump({'stage': stage, 'niveaux': {str(k + 1): nm for k, nm in enumerate(names)}, 'source_axe': os.path.relpath(topo, stage_dir)}, open(os.path.join(work, 'vertebres_niveaux.json'), 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    from meshexport import mask_to_mesh, CRL_MM
    si = np.where(env.any(axis=(0, 2)))[0]; scale = CRL_MM.get(stage, 20.0) / (np.ptp(si) + 1); c = np.array(ndi.center_of_mass(env))
    os.makedirs(os.path.join(out, 'vertebres'), exist_ok=True); man = {'stage': stage, 'units': 'mm', 'structures': []}
    for k, nm in enumerate(names):
        nm_s = re.sub(r'[^A-Za-z0-9]+', '_', nm)
        for pref, src, minvox in (('vertebre_', inst, 200), ('corps_', inst_corps, 60)):     # niveau complet (corps+arcs+côtes) et corps seul (blob compact, morphable)
            m = src == k + 1
            if m.sum() < minvox: continue
            mesh = mask_to_mesh(m, scale, sigma=0.8, target_faces=20000 if pref == 'vertebre_' else 6000)
            if mesh is None: continue
            mesh.apply_translation(-np.array([c[0], c[2], -c[1]]) * scale); fn = '%s_%s%s.ply' % (stage, '' if pref == 'vertebre_' else 'corps_', nm_s); mesh.export(os.path.join(out, 'vertebres', fn))
            man['structures'].append({'name': pref + nm_s, 'file': fn, 'collection': 'Squelette', 'color': [0.9, 0.9, 0.85], 'alpha': 1.0, 'faces': int(len(mesh.faces)), 'volume_mm3': float(m.sum()) * scale ** 3})
    json.dump(man, open(os.path.join(out, 'vertebres', 'manifest_vertebres.json'), 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    import meshexport, export_nrrd; meshexport.run(work, out, stage); export_nrrd.run(stage_dir)
print('done', stage)
