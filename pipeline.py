"""
Reconstruction 3D d'un embryon (3D Atlas of Human Embryology, vidéos de coupes) -> volume fusionné -> segmentation -> maillages.
usage :  python pipeline.py <dossier_stade> [--stage CS20] [--steps extract,fuse,prep,axes,segment,mesh] [--force]
Sorties dans <dossier_stade>/work (intermédiaires) et <dossier_stade>/out (maillages + manifest.json + volume NRRD).
"""
import argparse, os, sys, json, time, itertools, glob
import numpy as np, cv2
from scipy import ndimage as ndi
from scipy.ndimage import zoom, gaussian_filter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
T0 = time.time()
def log(*a):
    print('[%6.0fs]' % (time.time()-T0), *a); sys.stdout.flush()

def ball(r):
    z,y,x = np.ogrid[-r:r+1,-r:r+1,-r:r+1]; return (x*x+y*y+z*z) <= r*r
def largest(m):
    lab,n = ndi.label(m)
    if n==0: return m
    s = np.bincount(lab.ravel())[1:]; return lab == (np.argmax(s)+1)
def bbox(m):
    idx = np.where(m); return np.array([[i.min(), i.max()+1] for i in idx])

# =====================================================================================  1. extraction
def step_extract(stage_dir, work):
    vids = sorted(glob.glob(os.path.join(stage_dir, '*.mp4')))
    info = {'stacks': [], 'control': None}
    for vi, fn in enumerate(vids):
        cap = cv2.VideoCapture(fn); frames = []
        while True:
            ok, f = cap.read()
            if not ok: break
            frames.append(cv2.cvtColor(f, cv2.COLOR_BGR2GRAY))
        v = np.stack(frames)
        corner = np.concatenate([v[:, :4, :].reshape(len(v), -1), v[:, -4:, :].reshape(len(v), -1)], axis=1)
        bgm = np.median(corner)
        if bgm < 100:                       # fond noir -> vidéo de contrôle (rendu volumique)
            info['control'] = os.path.basename(fn); cv2.imwrite(os.path.join(work, 'control_frame.png'), v[0]); log('control video:', os.path.basename(fn)); continue
        d = np.abs(v[1:].astype(np.int16) - v[:-1].astype(np.int16)).mean(axis=(1,2))
        keep = np.concatenate([[True], d > 0.05]); v = v[keep]
        bg = np.median(corner[keep], axis=1).astype(np.float32)
        dens = np.clip(255 - v.astype(np.float32) * (255.0/np.maximum(bg,1))[:,None,None], 0, 255).astype(np.uint8)
        name = 'stack%d' % len(info['stacks'])
        np.save(os.path.join(work, name + '.npy'), dens)
        info['stacks'].append({'name': name, 'file': os.path.basename(fn), 'shape': list(dens.shape)})
        log(name, os.path.basename(fn), dens.shape, 'unique slices')
    assert len(info['stacks']) == 3, 'il faut 3 piles de coupes'
    json.dump(info, open(os.path.join(work, 'stacks.json'), 'w'), indent=1)

# =====================================================================================  2. orientation + recalage + fusion
def coarse_cube(v, b, n=48):
    sub = v[b[0,0]:b[0,1], b[1,0]:b[1,1], b[2,0]:b[2,1]].astype(np.float32)
    z = np.array([n,n,n]) / np.array(sub.shape)
    return zoom(gaussian_filter(sub, sigma=[max(0.5/zi, 0.01) for zi in z]), z, order=1)

def find_orient(ref_c, c0, topk=1):
    res = []
    for perm in itertools.permutations(range(3)):
        for flips in itertools.product([False,True], repeat=3):
            c = np.transpose(c0, perm)
            for i,f in enumerate(flips):
                if f: c = np.flip(c, axis=i)
            a = ref_c - ref_c.mean(); bb = c - c.mean()
            ncc = float((a*bb).sum() / np.sqrt((a*a).sum()*(bb*bb).sum()))
            res.append((ncc, perm, flips))
    res.sort(key=lambda r: -r[0])
    return res[0] if topk == 1 else res[:topk]

def quick_register_score(fixed_arr, moving_arr):
    """recalage affine rapide (cubes 48^3, spacing 1) -> métrique de corrélation (négative = mieux)"""
    import SimpleITK as sitk
    f = sitk.GetImageFromArray(fixed_arr.astype(np.float32)); m = sitk.GetImageFromArray(moving_arr.astype(np.float32))
    R = sitk.ImageRegistrationMethod(); R.SetMetricAsCorrelation(); R.SetInterpolator(sitk.sitkLinear)
    tx = sitk.CenteredTransformInitializer(f, m, sitk.AffineTransform(3), sitk.CenteredTransformInitializerFilter.GEOMETRY)
    R.SetInitialTransform(tx, inPlace=True)
    R.SetOptimizerAsRegularStepGradientDescent(learningRate=0.5, minStep=1e-4, numberOfIterations=150, relaxationFactor=0.6)
    R.SetOptimizerScalesFromPhysicalShift(); R.SetShrinkFactorsPerLevel([2,1]); R.SetSmoothingSigmasPerLevel([2,1])
    try:
        R.Execute(f, m); return R.GetMetricValue()
    except Exception as e:
        return 0.0

def apply_orient(v, perm, flips):
    c = np.transpose(v, perm)
    for i,f in enumerate(flips):
        if f: c = np.flip(c, axis=i)
    return np.ascontiguousarray(c)

def step_fuse(work, step=0.5):
    import SimpleITK as sitk
    info = json.load(open(os.path.join(work, 'stacks.json')))
    vols = {s['name']: np.load(os.path.join(work, s['name']+'.npy')) for s in info['stacks']}
    names = list(vols)
    # pile A = la plus grande hauteur d'image (pixel de référence)
    A = max(names, key=lambda k: vols[k].shape[1])
    bb = {k: bbox(vols[k] > 40) for k in names}
    refc = coarse_cube(vols[A], bb[A])
    orient = {A: ((0,1,2), (False,False,False))}
    for k in names:
        if k == A: continue
        ck = coarse_cube(vols[k], bb[k]); cands = find_orient(refc, ck, topk=6)
        scored = []
        for ncc, perm, flips in cands:
            c = np.transpose(ck, perm)
            for i,f in enumerate(flips):
                if f: c = np.flip(c, axis=i)
            scored.append((quick_register_score(refc, np.ascontiguousarray(c)), ncc, perm, flips))
        scored.sort(key=lambda r: r[0]); met, ncc, perm, flips = scored[0]; orient[k] = (perm, flips)
        log('orient', k, '-> axes of', A, perm, flips, 'ncc %.3f reg-metric %.3f' % (ncc, met), '| candidates:', [(round(r[0],3), r[2], r[3]) for r in scored[:3]])
    ov = {k: apply_orient(vols[k], *orient[k]) for k in names}
    ext = {k: (bbox(ov[k] > 40)[:,1] - bbox(ov[k] > 40)[:,0]).astype(float) for k in names}
    # échelle pixel de chaque pile relative à A ; extents physiques (unité = pixel in-plane de A)
    # axe 0 de A = axe de coupe (les tableaux sont (z=coupe, y, x) pour A)
    slice_axis = 0
    pxscale = {A: 1.0}
    others = [k for k in names if k != A]
    phys = np.zeros(3); phys[1] = ext[A][1]; phys[2] = ext[A][2]
    est = []
    for k in others:
        # dans k, quel axe est l'axe de coupe (nb de coupes >> extents) : celui dont l'extent est le plus proche du nb de coupes original
        sa = int(np.argmax([ext[k][a] / (ext[A][a] if a != slice_axis else 1e9) for a in range(3)]))  # l'axe avec le plus grand ratio est l'axe de coupe de k
        inplane = [a for a in range(3) if a != sa]
        common = [a for a in inplane if a != slice_axis]           # axe in-plane commun avec A
        if common:
            pxscale[k] = ext[A][common[0]] / ext[k][common[0]]
            if slice_axis in inplane: est.append(ext[k][slice_axis] * pxscale[k])
        log('stack', k, 'slice axis', sa, 'pixel scale vs A %.3f' % pxscale.get(k, float('nan')))
    phys[0] = np.mean(est) if est else ext[A][0]*0.25
    log('physical extents (A frame, px of A):', np.round(phys,1))

    def to_sitk(k, shrink=1.0):
        v = ov[k]; b = bbox(v > 40); e = (b[:,1]-b[:,0]).astype(float)
        spacing = phys / e
        sub = v[b[0,0]:b[0,1], b[1,0]:b[1,1], b[2,0]:b[2,1]].astype(np.float32)
        if shrink > 1:
            z = np.minimum(1.0, spacing/shrink); sub = zoom(sub, z, order=1); spacing = spacing/z
        img = sitk.GetImageFromArray(np.ascontiguousarray(sub)); img.SetSpacing(tuple(spacing[::-1].tolist()))
        return img, spacing
    def register(fixed, moving):
        R = sitk.ImageRegistrationMethod(); R.SetMetricAsCorrelation()
        R.SetMetricSamplingStrategy(R.RANDOM); R.SetMetricSamplingPercentage(0.25, seed=1)
        R.SetInterpolator(sitk.sitkLinear)
        tx = sitk.CenteredTransformInitializer(fixed, moving, sitk.AffineTransform(3), sitk.CenteredTransformInitializerFilter.GEOMETRY)
        R.SetInitialTransform(tx, inPlace=True)
        R.SetOptimizerAsRegularStepGradientDescent(learningRate=1.0, minStep=1e-5, numberOfIterations=400, relaxationFactor=0.6)
        R.SetOptimizerScalesFromPhysicalShift()
        R.SetShrinkFactorsPerLevel([8,4,2,1]); R.SetSmoothingSigmasPerLevel([6,3,1.5,0.7]); R.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
        out = R.Execute(fixed, moving); return out, R.GetMetricValue()
    # choix de la référence : celle pour laquelle les 2 autres se recalent avec le moins de cisaillement
    best = None
    for ref in names:
        fixed, _ = to_sitk(ref, shrink=1.5); shear = 0; txs = {}; mets = []
        for k in names:
            if k == ref: continue
            mov, _ = to_sitk(k, shrink=1.5); tx, m = register(fixed, mov); txs[k] = tx; mets.append(m)
            Amat = np.array(tx.GetMatrix()).reshape(3,3); shear += np.abs(Amat - np.diag(np.diag(Amat))).sum()
        log('ref candidate', ref, 'shear %.3f' % shear, 'metrics', np.round(mets,3))
        if best is None or shear < best[0]: best = (shear, ref, txs)
    shear, ref, txs = best; log('reference stack:', ref)
    # grille commune : pas normalisé pour ~790 voxels sur la plus grande dimension (comparable entre stades)
    step = float(phys.max() / 790.0); log('grid step %.4f (px of A)' % step)
    size_phys = phys + 6; shape = np.ceil(size_phys/step).astype(int); origin = -3.0
    refimg = sitk.Image([int(s) for s in shape[::-1]], sitk.sitkFloat32); refimg.SetSpacing((step,)*3); refimg.SetOrigin((origin,)*3)
    res = {}; native = {}
    for k in names:
        img, sp = to_sitk(k, 1.0); native[k] = sp
        tx = txs.get(k, sitk.Transform(3, sitk.sitkIdentity))
        r = sitk.Resample(img, refimg, tx, sitk.sitkBSpline, 0.0, sitk.sitkFloat32)
        res[k] = np.clip(sitk.GetArrayFromImage(r), 0, 255); log('resampled', k, res[k].shape, 'native spacing', np.round(sp,3))
    # appariement d'histogramme (tissu seulement) puis remise du fond à 0
    ref_t = res[A][res[A] > 20]
    for k in names:
        if k == A: continue
        a = res[k]; src = np.sort(a[a>20][::97]); dst = np.sort(ref_t[::97]); q = np.linspace(0,1,257)
        lut_x = np.quantile(src,q); lut_y = np.quantile(dst,q); lut_x[0] = 0; lut_y[0] = 0
        res[k] = np.interp(a, lut_x, lut_y).astype(np.float32)
    # fusion de Fourier pondérée par la résolution native de chaque pile
    import scipy.fft as F                      # scipy.fft conserve la simple précision (complex64) -> 2x moins de mémoire que numpy.fft
    num = np.zeros(shape, np.complex64); den = np.zeros(shape, np.float32)
    freqs = [F.fftfreq(int(n), d=step).astype(np.float32) for n in shape]
    for k in names:
        w = np.ones(shape, np.float32)
        for ax in range(3):
            sig = 0.6*native[k][ax]; wa = np.exp(-2*np.pi**2 * sig**2 * freqs[ax]**2).astype(np.float32)
            sh = [1,1,1]; sh[ax] = -1; w = w * wa.reshape(sh)
        Fk = F.fftn(res[k].astype(np.float32), workers=-1); Fk *= w; num += Fk; del Fk
        den += w; del w; res[k] = None; log('fft', k)
    num /= np.maximum(den, 1e-3); del den
    fused = np.clip(np.real(F.ifftn(num, workers=-1)), 0, 255).astype(np.uint8)
    del num, res
    corner = np.concatenate([fused[:8].ravel(), fused[-8:].ravel(), fused[:, :8].ravel(), fused[:, -8:].ravel()])
    bgv = float(np.median(corner)); dens = np.clip((fused.astype(np.float32) - bgv) * (255.0/(255.0-bgv)), 0, 255).astype(np.uint8)
    np.save(os.path.join(work, 'dens.npy'), dens)
    json.dump({'step': step, 'origin': origin, 'shape': shape.tolist(), 'phys_ext': phys.tolist(), 'reference': ref, 'A': A,
               'native_spacing': {k: v.tolist() for k,v in native.items()}, 'orient': {k: [list(o[0]), list(o[1])] for k,o in orient.items()},
               'pixel_unit': 'in-plane pixel of stack %s' % A}, open(os.path.join(work, 'fused_meta.json'), 'w'), indent=1)
    log('fused volume', dens.shape)

# =====================================================================================  3. enveloppe + cavités
def envelope_by_occlusion(tissue, D=90, thresh=24):
    t2 = np.pad(ndi.maximum_filter(tissue, size=2)[::2, ::2, ::2], 1)
    dirs = [np.array(v) for v in itertools.product([-1,0,1], repeat=3) if any(v)]
    cnt = np.zeros(t2.shape, np.uint8)
    for v in dirs:
        b = t2.copy()
        for s in range(D):
            sh = np.roll(b, shift=tuple(v), axis=(0,1,2))
            for ax,c in enumerate(v):
                if c==1: sh[(slice(None),)*ax + (0,)] = False
                elif c==-1: sh[(slice(None),)*ax + (-1,)] = False
            nb = b | sh
            if s > 8 and (nb == b).all(): break
            b = nb
        cnt += b.astype(np.uint8)
    inside2 = (cnt >= thresh)[1:-1,1:-1,1:-1]
    up = np.repeat(np.repeat(np.repeat(inside2, 2, 0), 2, 1), 2, 2)
    inside = np.zeros(tissue.shape, bool); L = [min(a,b) for a,b in zip(up.shape, tissue.shape)]
    inside[:L[0], :L[1], :L[2]] = up[:L[0], :L[1], :L[2]]
    return inside

REF_PCT = [25, 50, 75, 90, 99]            # percentiles de densité du tissu de référence (CS20) ...
REF_VAL = [35., 75., 125., 170., 207.]     # ... et leurs valeurs : les autres stades sont ramenés sur cette échelle
def normalize_density(d):
    """LUT linéaire par morceaux : les percentiles du tissu (d>8) sont alignés sur ceux de la référence ; 0->0, 255->255"""
    sub = d[::3, ::3, ::3]; t = sub[sub > 8].astype(np.float32)
    q = np.percentile(t, REF_PCT)
    xs = np.concatenate([[0.0, 8.0], q, [255.0]]); ys = np.concatenate([[0.0, 8.0], REF_VAL, [255.0]])
    xs = np.maximum.accumulate(xs + np.arange(len(xs))*1e-3)
    lut = np.interp(np.arange(256, dtype=np.float32), xs, ys)
    return np.clip(lut[d], 0, 255).astype(np.uint8), q

def step_prep(work):
    d = np.load(os.path.join(work, 'dens.npy'))
    raw = os.path.join(work, 'dens_raw.npy')
    if not os.path.exists(raw): np.save(raw, d)
    else: d = np.load(raw)
    d, q = normalize_density(d); np.save(os.path.join(work, 'dens.npy'), d); log('density normalized; tissue percentiles were', np.round(q).astype(int))
    ds = gaussian_filter(d.astype(np.float32), 1.0); np.save(os.path.join(work, 'ds1.npy'), ds.astype(np.float16))
    tissue = ds > 4
    env = tissue | envelope_by_occlusion(tissue)
    env = ndi.binary_closing(env, structure=ball(2)); env = ndi.binary_fill_holes(env); env = largest(env)
    env = ndi.binary_opening(env, structure=ball(2)); env = largest(env)
    core = env.copy()
    for _ in range(7): core = ndi.binary_erosion(core, border_value=0)
    env = env & ~((env & ~core) & (ds < 14))
    env = gaussian_filter(env.astype(np.float32), 1.5) > 0.5; env = ndi.binary_fill_holes(env); env = largest(env)
    np.save(os.path.join(work, 'env.npy'), env); log('envelope %.1fM vox' % (env.sum()/1e6))
    # cavités : ouverture r=3 pour couper le réseau de poches du mésenchyme
    cav = env & (ds < 3.0)
    cav_o = ndi.binary_opening(cav, structure=ball(3))
    lab, n = ndi.label(cav_o); sizes = np.bincount(lab.ravel(), minlength=n+1)[1:].astype(np.float64); order = np.argsort(-sizes); objs = ndi.find_objects(lab)
    dark = d > 90; rows = []
    for r,i in enumerate(order[:60]):
        sl = objs[i]; pad=6
        sl2 = tuple(slice(max(s.start-pad,0), min(s.stop+pad, sh)) for s,sh in zip(sl, lab.shape))
        m = lab[sl2] == (i+1)
        shl = ndi.binary_dilation(m, structure=ball(5)) & ~ndi.binary_dilation(m, structure=ball(2))
        fd = float(dark[sl2][shl].mean()) if shl.any() else 0.0
        zz,yy,xx = np.where(m); c = np.array([zz.mean()+sl2[0].start, yy.mean()+sl2[1].start, xx.mean()+sl2[2].start])
        ext = (np.ptp(zz)+1, np.ptp(yy)+1, np.ptp(xx)+1)
        P = np.stack([zz,yy,xx],1).astype(np.float32); P -= P.mean(0)
        ev = np.sort(np.linalg.eigvalsh(np.cov(P.T)))[::-1] if len(P) > 3 else np.ones(3); el = float(np.sqrt(ev[0]/max(ev[2],1e-6)))
        rows.append([int(i+1), int(sizes[i]), *c.tolist(), *[int(e) for e in ext], fd, el])
    np.save(os.path.join(work, 'cavo_lab.npy'), lab.astype(np.int32)); np.save(os.path.join(work, 'cavo_sizes.npy'), sizes); np.save(os.path.join(work, 'cavo_rows.npy'), np.array(rows))
    log('cavities: %d comps, top sizes' % n, [int(r[1]/1000) for r in rows[:8]], 'k')

# =====================================================================================  4. axes anatomiques -> orientation canonique (LR, SI haut->bas, AP avant->arrière)
def step_axes(work):
    d = np.load(os.path.join(work, 'dens.npy')); ds = np.load(os.path.join(work, 'ds1.npy')); env = np.load(os.path.join(work, 'env.npy'))
    cavo = np.load(os.path.join(work, 'cavo_lab.npy')); rows = np.load(os.path.join(work, 'cavo_rows.npy')); sizes = np.load(os.path.join(work, 'cavo_sizes.npy'))
    # LR = axe de meilleure symétrie miroir (sur le volume lissé sous-échantillonné)
    small = zoom(gaussian_filter(d.astype(np.float32), 2.0), 0.25, order=1)
    sym = []
    for ax in range(3):
        f = np.flip(small, axis=ax); a = small - small.mean(); b = f - f.mean()
        sym.append(float((a*b).sum()/np.sqrt((a*a).sum()*(b*b).sum())))
    # règle géométrique : gauche-droite = plus petite emprise de l'enveloppe, haut-bas = plus grande (la symétrie miroir n'est qu'un contrôle)
    ext = [int(np.ptp(np.where(env.any(axis=tuple(b for b in range(3) if b != a)))[0])) for a in range(3)]
    LRax = int(np.argmin(ext)); SIax = int(np.argmax(ext)); APax = [a for a in range(3) if a not in (LRax, SIax)][0]
    log('envelope extents per axis', ext)
    log('mirror symmetry per axis', np.round(sym,3), '-> LR axis', LRax, 'SI axis', SIax, 'AP axis', APax)
    perm = (LRax, SIax, APax)
    def P(v): return np.ascontiguousarray(np.transpose(v, perm))
    d, ds, env, cavo = P(d), P(ds), P(env), P(cavo)
    # ventricules = la grande cavité dont la paroi est la plus épaisse et la plus dense (profil radial de ds à 4..14 voxels)
    from segment import _crop
    top = sorted(rows, key=lambda r: -r[1])[:10]
    best = None
    for r in top:
        if r[1] < 60000: continue
        m = cavo == int(r[0]); sl = _crop(m, 16); mm = m[sl]; dsc = ds[sl].astype(np.float32); from segment import _dist; dist = _dist(~mm)
        prof = [float(np.median(dsc[(dist > k-1) & (dist <= k)])) for k in range(4, 15)]
        score = float(np.mean(prof)); log('  cavity lab %d size %dk wall score %.0f' % (r[0], r[1]/1000, score))
        if best is None or score > best[0]: best = (score, int(r[0]))
    vent_lab = best[1]; vent = cavo == vent_lab
    cv = np.array(ndi.center_of_mass(vent))
    # foie = plus gros noyau très sombre hors SNC (ouverture r=6)
    darkL = (ds > 110) & env
    from segment import open_r, dilate_r
    coreL = open_r(darkL, 6); labL, nL = ndi.label(coreL); sL = np.bincount(labL.ravel(), minlength=nL+1)[1:]
    # exclure les noyaux qui touchent la dilatation des ventricules (paroi neurale)
    nearv = dilate_r(vent, 30)
    cand = []
    for i in range(nL):
        m = labL == (i+1)
        if (m & nearv).sum() < 0.05*sL[i]: cand.append((sL[i], i+1))
    cand.sort(reverse=True); liver_lab = cand[0][1]; cl = np.array(ndi.center_of_mass(labL == liver_lab))
    # SNC (noyau érodé connecté aux ventricules) : la moelle est DORSALE -> définit l'axe AP ; la tête (ventricules) est en HAUT
    from segment import erode_r, geodesic
    core = erode_r((ds > 85) & env, 2)
    cns_core = geodesic(dilate_r(vent, 10) & core, core, 1000)
    ce = np.array(ndi.center_of_mass(env))
    si_idx = np.where(env.any(axis=(0,2)))[0]
    # partie caudale du corps (moitié opposée aux ventricules) : comparer la position AP de la moelle et du corps
    if cv[1] < ce[1]: caudal = slice(int(ce[1]), si_idx.max()+1)
    else: caudal = slice(si_idx.min(), int(ce[1]))
    cns_c = cns_core[:, caudal, :]; env_c = env[:, caudal, :]
    ap_cns = ndi.center_of_mass(cns_c)[2] if cns_c.any() else ce[2]; ap_env = ndi.center_of_mass(env_c)[2]
    log('ventricle centroid', cv.astype(int), 'env centroid', ce.astype(int), 'liver centroid', cl.astype(int), 'caudal AP cord %.0f vs body %.0f' % (ap_cns, ap_env))
    flips = [False, False, False]
    if cv[1] > ce[1]: flips[1] = True                 # tête en haut (SI petit)
    if ap_cns < ap_env: flips[2] = True               # moelle dorsale = AP grand
    if cl[0] < d.shape[0]/2: flips[0] = True          # foie plutôt à droite = LR grand (handedness approximative)
    def Fl(v):
        for ax,f in enumerate(flips):
            if f: v = np.flip(v, axis=ax)
        return np.ascontiguousarray(v)
    d, ds, env, cavo = Fl(d), Fl(ds), Fl(env), Fl(cavo)
    np.save(os.path.join(work, 'dens.npy'), d); np.save(os.path.join(work, 'ds1.npy'), ds); np.save(os.path.join(work, 'env.npy'), env); np.save(os.path.join(work, 'cavo_lab.npy'), cavo)
    # recalcul des lignes (centroïdes) dans le nouveau repère
    objs = ndi.find_objects(cavo); rows2 = []
    for r in rows:
        i = int(r[0]); sl = objs[i-1]
        if sl is None: continue
        m = cavo[sl] == i; zz,yy,xx = np.where(m); c = [zz.mean()+sl[0].start, yy.mean()+sl[1].start, xx.mean()+sl[2].start]
        e = [np.ptp(zz)+1, np.ptp(yy)+1, np.ptp(xx)+1]
        rows2.append([i, int(r[1]), *c, *e, r[8], r[9]])
    np.save(os.path.join(work, 'cavo_rows.npy'), np.array(rows2))
    coreL = Fl(coreL); labL = Fl(labL)
    np.save(os.path.join(work, 'liver_core.npy'), labL == liver_lab)
    json.dump({'perm': list(perm), 'flips': flips, 'symmetry': sym, 'vent_lab': vent_lab, 'shape': list(d.shape)}, open(os.path.join(work, 'axes.json'), 'w'), indent=1)
    log('canonical shape (LR,SI,AP)', d.shape, 'flips', flips)

# =====================================================================================  main
STEPS = ['extract', 'fuse', 'prep', 'axes', 'segment', 'mesh']
if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('stage_dir'); ap.add_argument('--stage', default=None); ap.add_argument('--steps', default=','.join(STEPS)); ap.add_argument('--force', action='store_true')
    args = ap.parse_args()
    stage_dir = os.path.abspath(args.stage_dir); import re as _re
    _m = _re.search(r'CS\s*(\d+)', os.path.basename(stage_dir), _re.I)
    stage = args.stage or ('CS%s' % _m.group(1) if _m else os.path.basename(stage_dir))
    work = os.path.join(stage_dir, 'work'); out = os.path.join(stage_dir, 'out'); os.makedirs(work, exist_ok=True); os.makedirs(out, exist_ok=True)
    steps = args.steps.split(',')
    done = lambda f: os.path.exists(os.path.join(work, f)) and not args.force
    if 'extract' in steps and not done('stacks.json'): step_extract(stage_dir, work)
    if 'fuse' in steps and not done('fused_meta.json'): step_fuse(work)
    if 'prep' in steps and not done('cavo_rows.npy'): step_prep(work)
    if 'axes' in steps and not done('axes.json'): step_axes(work)
    if 'segment' in steps:
        import segment; segment.run(work, stage)
    if 'mesh' in steps:
        import meshexport, export_nrrd; meshexport.run(work, out, stage); export_nrrd.run(stage_dir)
    log('pipeline done for', stage)
