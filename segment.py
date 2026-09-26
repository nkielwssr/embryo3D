"""Segmentation heuristique du volume fusionné (repère canonique LR, SI, AP). Résultat : work/labels.npz (masques packbits) + planche de contrôle."""
import numpy as np, cv2, time, json, sys, os
from scipy import ndimage as ndi
T0 = time.time()
def log(*a): print('[%6.0fs] seg' % (time.time()-T0), *a); sys.stdout.flush()
def ball(r):
    z,y,x = np.ogrid[-r:r+1,-r:r+1,-r:r+1]; return (x*x+y*y+z*z) <= r*r
def largest(m):
    lab,n = ndi.label(m)
    if n==0: return m
    s = np.bincount(lab.ravel())[1:]; return lab == (np.argmax(s)+1)
def _crop(m, pad):
    """boîte englobante (slices) de m élargie de pad, via projections (pas de np.where sur le volume entier)"""
    lo = []; hi = []
    for ax in range(3):
        proj = m.any(axis=tuple(a for a in range(3) if a != ax)); idx = np.flatnonzero(proj)
        if len(idx) == 0: return None
        lo.append(max(int(idx[0]) - pad, 0)); hi.append(min(int(idx[-1]) + pad + 1, m.shape[ax]))
    return tuple(slice(a, b) for a, b in zip(lo, hi))
def geodesic(seed, mask, steps):
    """reconstruction géodésique de seed dans mask. steps >= 100 : illimitée (composantes connexes de mask touchant seed) ;
    sinon limitée : composantes touchant seed, coupées à une distance euclidienne <= steps de seed (approximation rapide)."""
    seed = seed & mask
    if not seed.any(): return seed
    lab, n = ndi.label(mask)
    ids = np.unique(lab[seed]); ids = ids[ids > 0]
    rec = np.isin(lab, ids)
    if steps >= 100: return rec
    sl = _crop(rec, 1)
    out = np.zeros_like(rec); sub = rec[sl]
    out[sl] = sub & (_dist(~seed[sl]) <= steps)
    return out

def _crop(m, pad):
    """boîte englobante (slices) de m élargie de pad, via projections (pas de np.where sur le volume entier)"""
    lo = []; hi = []
    for ax in range(3):
        proj = m.any(axis=tuple(a for a in range(3) if a != ax)); idx = np.flatnonzero(proj)
        if len(idx) == 0: return None
        lo.append(max(int(idx[0]) - pad, 0)); hi.append(min(int(idx[-1]) + pad + 1, m.shape[ax]))
    return tuple(slice(a, b) for a, b in zip(lo, hi))
try:
    import edt as _edt
    def _dist(mask):
        """distance euclidienne à l'intérieur de mask (0 hors mask), float32, multithread"""
        return _edt.edt(np.ascontiguousarray(mask, dtype=np.uint8), parallel=0)
except ImportError:
    def _dist(mask): return ndi.distance_transform_edt(mask)
def dilate_r(m, r):
    sl = _crop(m, r+2)
    if sl is None: return m.copy()
    out = np.zeros_like(m); sub = m[sl]
    out[sl] = _dist(~sub) <= r if r > 1 else ndi.binary_dilation(sub)
    return out
def erode_r(m, r):
    sl = _crop(m, 2)
    if sl is None: return m.copy()
    out = np.zeros_like(m); sub = np.pad(m[sl], 1)
    out[sl] = (_dist(sub) > r)[1:-1, 1:-1, 1:-1] if r > 1 else ndi.binary_erosion(sub, border_value=0)[1:-1,1:-1,1:-1]
    return out
def open_r(m, r):  return dilate_r(erode_r(m, r), r)
def close_r(m, r): return erode_r(dilate_r(m, r), r)

def comps_info(m, minsize=0):
    lab,n = ndi.label(m); s = np.bincount(lab.ravel(), minlength=n+1)[1:]; objs = ndi.find_objects(lab); out=[]
    for i in range(n):
        if s[i] < minsize: continue
        sl = objs[i]; c = np.array(ndi.center_of_mass(lab[sl]==(i+1))) + np.array([x.start for x in sl]); out.append((i+1, int(s[i]), c))
    return lab, out

PALETTE = {'ventricules':(0,200,255),'snc':(255,120,0),'foie':(120,0,160),'yeux':(0,0,255),'cristallins':(255,255,0),'vesicules_otiques':(0,255,255),
 'squelette_axial_cartilage':(255,255,255),'chondrocrane':(230,230,200),'cartilage_autre':(200,200,230),'coeur':(0,0,200),'cavite_pericardique':(180,220,255),
 'membres':(0,160,0),'cordon_ombilical':(160,120,60),'tube_digestif':(255,0,255),'ganglions':(255,150,150),'vaisseaux':(200,0,0)}

def guarded(name, fn):
    try:
        return fn()
    except Exception as e:
        log('!! %s skipped: %r' % (name, e)); return None

def run(work, stage):
    d = np.load(os.path.join(work,'dens.npy')); LR,SI,AP = d.shape
    ds = np.load(os.path.join(work,'ds1.npy')).astype(np.float32); env = np.load(os.path.join(work,'env.npy'))
    cavo = np.load(os.path.join(work,'cavo_lab.npy')); rows = np.load(os.path.join(work,'cavo_rows.npy'))
    axes = json.load(open(os.path.join(work,'axes.json'))); vent_lab = axes['vent_lab']
    liver_core = np.load(os.path.join(work,'liver_core.npy'))
    # bornes de l'embryon (fractions calculées sur l'enveloppe)
    si_idx = np.where(env.any(axis=(0,2)))[0]; si0, si1 = si_idx.min(), si_idx.max(); SIe = si1 - si0
    lr_idx = np.where(env.any(axis=(1,2)))[0]; MID = 0.5*(lr_idx.min()+lr_idx.max()); LRe = lr_idx.max()-lr_idx.min()
    ap_idx = np.where(env.any(axis=(0,1)))[0]; ap0, ap1 = ap_idx.min(), ap_idx.max(); APe = ap1 - ap0
    fS = lambda f: si0 + f*SIe          # fraction de la hauteur (0 = sommet de la tête)
    fA = lambda f: ap0 + f*APe          # fraction antéro-postérieure (0 = ventral)
    log('embryo box SI %d..%d LR mid %.0f (ext %d) AP %d..%d' % (si0, si1, MID, LRe, ap0, ap1))
    labels = {}
    cav = env & (ds < 3.0)
    # ---- ventricules + SNC
    vent_labs = axes.get('vent_labs') or [vent_lab]
    vent = geodesic(np.isin(cavo, vent_labs), cav, 6)
    # canal central / recessus fins : lumières entourées de tissu sombre (les poches du mésenchyme sont bordées de tissu clair)
    enclosed = ndi.uniform_filter((ds > 90).astype(np.float32), size=9) > 0.45
    vent = geodesic(vent, cav & (enclosed | vent), 400)
    labels['ventricules'] = vent; log('ventricles+canal %.2fM' % (vent.sum()/1e6))
    darkT = (ds > 55) & env
    brain = geodesic(dilate_r(vent, 3) & darkT, darkT, 28)                    # cerveau : paroi neurale autour des ventricules
    brain = largest(close_r(brain, 2) & env)
    deep = erode_r(env, 4)                                                     # hors couche superficielle (ectoderme / dermomyotome)
    midband = np.zeros_like(env); midband[int(MID - 0.13*LRe):int(MID + 0.13*LRe) + 1, :, :] = True
    midband[:, :int(fS(0.33)), :] = True                                       # pas de contrainte dans la tête
    core = erode_r((ds > 130) & deep & midband, 3)                             # moelle : noyau très dense érodé, connecté au cerveau, bande médiane
    cord_core = geodesic(dilate_r(brain, 6) & core, core, 1000)
    cord = geodesic(dilate_r(cord_core, 4) & (ds > 70) & env & midband, (ds > 70) & env & midband, 2)
    cns = largest(brain | cord) & ~vent; labels['snc'] = cns
    log('CNS %.2fM (brain %.2fM, cord core %.2fM)' % (cns.sum()/1e6, brain.sum()/1e6, cord_core.sum()/1e6))
    # ---- foie
    # foie : noyau très dense compact, ventral, dans le tronc (rejette la bande somitique / neurale des jeunes stades)
    coreL = open_r((ds > 110) & env & ~dilate_r(cns, 6), 6)
    labL, infoL = comps_info(coreL, 60000); infoL.sort(key=lambda x: -x[1]); liver_core2 = None
    for i, s_, c in infoL[:8]:
        m = labL == i; sl = _crop(m, 1); zz, yy, xx = np.where(m[sl]); P = np.stack([zz, yy, xx], 1).astype(np.float32); P -= P.mean(0)
        ev = np.sort(np.linalg.eigvalsh(np.cov(P.T)))[::-1]; el = float(np.sqrt(ev[0] / max(ev[2], 1e-6)))
        ok = el < 3.5 and c[2] < fA(0.7) and c[1] > fS(0.35) and s_ < 0.14 * env.sum()
        log('  liver candidate %d size %dk centroid %s elong %.1f %s' % (i, s_/1000, c.astype(int), el, 'OK' if ok else 'rejected'))
        if ok: liver_core2 = m; break
    if liver_core2 is None:
        # repli : graine = noyau dense compact dans une boîte ventrale du milieu du tronc, puis croissance bornée
        box = np.zeros_like(env); box[:, int(fS(0.40)):int(fS(0.72)), :int(fA(0.55))] = True
        coreB = open_r((ds > 120) & env & box & ~dilate_r(cns, 6), 5); labB, infoB = comps_info(coreB, 100000); infoB.sort(key=lambda x: -x[1])
        if infoB:
            i, s_, c = infoB[0]; liver_core2 = geodesic(labB == i, (ds > 100) & env & ~dilate_r(cns, 4), 30)
            log('  liver: fallback box seed %dk at %s -> grown %.2fM' % (s_/1000, c.astype(int), liver_core2.sum()/1e6))
        else:
            liver_core2 = np.zeros_like(env); log('  liver: no acceptable candidate -> empty')
    liver = geodesic(liver_core2, (ds > 80) & env & ~cns, 10)
    liver = largest(close_r(liver, 3) & env); labels['foie'] = liver; log('liver %.2fM' % (liver.sum()/1e6))
    # ---- paires symétriques (yeux = paire la plus antérieure, otiques = la plus postérieure)
    def find_pairs(cands):
        pairs=[]
        for a in cands:
            for b in cands:
                if a[0] >= b[0]: continue
                ca, cb = a[2:5], b[2:5]
                if abs((ca[0]-MID)+(cb[0]-MID)) < 0.075*LRe and abs(ca[0]-cb[0]) > 0.17*LRe and abs(ca[1]-cb[1]) < 0.035*SIe and abs(ca[2]-cb[2]) < 0.06*APe:
                    pairs.append((a,b))
        return pairs
    cands = [r for r in rows if 1500 < r[1] < 60000 and r[3] < fS(0.53)]
    pairs = find_pairs(cands)
    for a,b in pairs: log('pair labs %d,%d centroids %s %s sizes %d %d' % (a[0], b[0], a[2:5].astype(int), b[2:5].astype(int), a[1], b[1]))
    biglab = int([r for r in sorted(rows, key=lambda r:-r[1]) if int(r[0]) not in vent_labs][0][0])   # cavité corporelle (coelome + réseau)
    body_labs = [int(r[0]) for r in rows if r[1] >= 120000 and int(r[0]) not in vent_labs]     # cavités corporelles (coelome, réseau de mésenchyme)
    used = set(vent_labs) | {biglab} | set(body_labs); eye = np.zeros_like(env); peri = None; heart = np.zeros_like(env)
    def _pairs():
      nonlocal eye, used
      if True:
        pairs.sort(key=lambda p: (p[0][4]+p[1][4])); eye_p, otic_p = pairs[0], pairs[-1]
        used |= set(int(x) for x in [eye_p[0][0], eye_p[1][0], otic_p[0][0], otic_p[1][0]])
        intraoc = geodesic((cavo==eye_p[0][0])|(cavo==eye_p[1][0]), cav, 4)
        otic = geodesic((cavo==otic_p[0][0])|(cavo==otic_p[1][0]), cav, 4); labels['vesicules_otiques'] = otic
        eye = geodesic(dilate_r(intraoc, 4) & (ds > 50) & ~cns, (ds > 50) & env & ~cns, 14)
        eye = ndi.binary_fill_holes(close_r(eye | intraoc, 3))
        lab, info = comps_info(eye, 500); info.sort(key=lambda x:-x[1]); eye = np.isin(lab, [i[0] for i in info[:2]]); labels['yeux'] = eye
        lensc = open_r(eye & (ds > 110), 2); lab, info = comps_info(lensc, 300); info.sort(key=lambda x:-x[1])
        labels['cristallins'] = np.isin(lab, [i[0] for i in info[:2]]) if info else lensc
        log('eyes %.2fM lens %.3fM otic %.3fM' % (eye.sum()/1e6, labels['cristallins'].sum()/1e6, otic.sum()/1e6))
    if pairs: guarded('pairs', _pairs)
    # ---- cartilages : blobs pâles compacts
    sizes_all = np.load(os.path.join(work,'cavo_sizes.npy'))
    cart_labs = [int(r[0]) for r in rows if int(r[0]) not in used and 1500 < r[1] < 120000 and r[8] < 0.25 and r[9] < 7]
    small = [i+1 for i in range(len(sizes_all)) if 800 <= sizes_all[i] <= 1500 and (i+1) not in used]
    cart = geodesic(np.isin(cavo, cart_labs + small), cav, 3) & ~vent
    near_cns = dilate_r(cns, 45)
    lab, info = comps_info(cart, 1); axial = np.zeros_like(cart); head = np.zeros_like(cart); other = np.zeros_like(cart)
    for i, s, c in info:
        m = lab == i
        if c[1] < fS(0.36): head |= m
        elif (m & near_cns).sum() > 0.05*s: axial |= m
        else: other |= m
    labels['squelette_axial_cartilage'] = axial; labels['chondrocrane'] = head; labels['cartilage_autre'] = other
    log('cartilage axial %.2fM head %.2fM other %.2fM' % (axial.sum()/1e6, head.sum()/1e6, other.sum()/1e6))
    # ---- cavité péricardique + coeur : cavité épaisse touchant la boîte thoracique, hors voisinage du foie ; coeur = tissu dans son enveloppe convexe
    big = np.isin(cavo, body_labs + [biglab]) & ~dilate_r(vent, 6)
    thick = open_r(big, 6)
    liver_ok = liver.sum() > 300000 and (np.where(liver.any(axis=(0,2)))[0].mean() > fS(0.45))
    top_liver = int(np.where(liver.any(axis=(0,2)))[0].min()) if liver_ok else int(fS(0.66))
    if not liver_ok: log('liver doubtful (%.2fM) -> default thoracic box' % (liver.sum()/1e6))
    box = np.zeros_like(env); box[:, int(fS(0.34)):top_liver+10, :int(fA(0.62))] = True
    lab, info = comps_info(thick & box, 20000); info.sort(key=lambda x:-x[1])
    for i,s_,c in info[:6]: log('thoracic thick cavity comp %d size %dk centroid %s' % (i, s_/1000, c.astype(int)))
    if info:
        region = np.zeros_like(env); region[:, :int(top_liver + 0.15*SIe), :int(fA(0.82))] = True
        peri = geodesic(lab == info[0][0], big & region & ~dilate_r(liver, 8) & ~dilate_r(cns, 4), 80)
        peri = largest(open_r(peri, 3))
        try:
            from scipy.spatial import Delaunay
            sl = _crop(peri, 4); sub = peri[sl]; surf = sub & ~ndi.binary_erosion(sub)
            pts = np.argwhere(surf)[::11].astype(np.float64); tri = Delaunay(pts)
            g2 = np.indices(tuple((x+1)//2 for x in sub.shape)).reshape(3,-1).T.astype(np.float64) * 2.0
            inside2 = (tri.find_simplex(g2) >= 0).reshape(tuple((x+1)//2 for x in sub.shape))
            hull = np.zeros_like(env); hull[sl] = np.repeat(np.repeat(np.repeat(inside2, 2, 0), 2, 1), 2, 2)[:sub.shape[0], :sub.shape[1], :sub.shape[2]]
        except Exception as e:
            log('convex hull failed, fallback closing:', e); hull = ndi.binary_fill_holes(close_r(peri, 12))
        heart = largest(open_r(hull & ~peri & (ds > 25) & ~liver & ~cns & env, 3))
        heart = geodesic(heart, hull & ~peri & (ds > 20) & env, 4)
        heart = ndi.binary_fill_holes(close_r(heart, 5)) & hull & ~peri & ~liver & ~cns & env     # coeur plein (trabécules comblées)
        labels['coeur'] = heart; labels['cavite_pericardique'] = peri
        log('heart %.2fM pericard %.2fM' % (heart.sum()/1e6, peri.sum()/1e6))
    # ---- membres + cordon : protubérances de l'enveloppe (ouverture r=45), latérales = membres, médianes ventrales basses = cordon
    core = open_r(env, 45); prot = open_r(env & ~core, 3)
    lab, info = comps_info(prot, 30000); info.sort(key=lambda x:-x[1]); limbs = np.zeros_like(env); cord = np.zeros_like(env)
    lr_off = np.abs(np.arange(LR) - MID)[:, None, None] > 0.13*LRe
    for i,s_,c in info[:12]:
        m = lab == i; log('protrusion %d size %dk centroid %s' % (i, s_/1000, c.astype(int)))
        if c[1] < fS(0.38): continue                         # tête
        lat = m & lr_off; med = m & ~lr_off
        if lat.sum() > 15000: limbs |= lat
        if med.sum() > 15000 and c[2] < fA(0.5) and c[1] > fS(0.55): cord |= med
    limbs = open_r(limbs, 2); cord = open_r(cord, 2)
    labels['membres'] = limbs; labels['cordon_ombilical'] = cord; log('limbs %.2fM cord %.2fM' % (limbs.sum()/1e6, cord.sum()/1e6))
    # ---- tube digestif
    excl = cns | liver | vent | cart | heart
    if peri is not None: excl |= peri
    lum_cands = [int(r[0]) for r in rows if r[8] > 0.25 and r[3] > fS(0.44) and int(r[0]) not in used]
    lum = np.isin(cavo, lum_cands) & ~excl
    gut = close_r((dilate_r(lum, 6) & (ds > 70) & env & ~excl) | lum, 2)
    lab, info = comps_info(gut, 3000); gut = np.isin(lab, [i[0] for i in info]); labels['tube_digestif'] = gut; log('gut %.2fM' % (gut.sum()/1e6))
    # ---- ganglions
    darkG = (ds > 100) & env & ~cns & ~liver & ~eye & ~gut & ~heart
    gcore = open_r(darkG, 2)
    nearsc = dilate_r(cns, 18) & ~dilate_r(cns, 2)
    lab, info = comps_info(gcore & nearsc, 400); gang = np.zeros_like(env)
    for i,s,c in info:
        if s <= 40000 and c[1] > fS(0.30): gang |= (lab == i)
    labels['ganglions'] = geodesic(gang, darkG, 3); log('ganglia %.3fM' % (labels['ganglions'].sum()/1e6))
    # ---- vaisseaux : tubes creux (vesselness de Frangi sur l'indicateur de lumière, grille 1 px)
    try:
        from skimage.filters import frangi
        excl2 = vent | cart | eye | gut
        if peri is not None: excl2 |= peri
        if 'vesicules_otiques' in labels: excl2 |= labels['vesicules_otiques']
        lum_ind = ((ds < 8) & env & ~excl2).astype(np.float32)
        small_v = ndi.zoom(ndi.gaussian_filter(lum_ind, 1.0), 0.5, order=1)
        V = frangi(small_v, sigmas=[1.0, 1.6, 2.4], black_ridges=False, alpha=0.5, beta=0.5, gamma=None)
        thr = np.percentile(V[V > 0], 97.5) if (V > 0).any() else 1.0
        vm = ndi.zoom((V > thr).astype(np.float32), 2.0, order=1) > 0.5
        vm2 = np.zeros_like(env); L = [min(a,b) for a,b in zip(vm.shape, env.shape)]; vm2[:L[0],:L[1],:L[2]] = vm[:L[0],:L[1],:L[2]]
        ves = vm2 & env & ~excl2 & ~cns & ~liver & ~heart
        ves = geodesic(ves, (ds < 10) & env & ~excl2 & ~cns, 2)
        lab, info = comps_info(ves, 300); ves = np.isin(lab, [i[0] for i in info if i[1] >= 300])
        # filtrage : composantes >= 3000 vox et allongées
        lab, info = comps_info(ves, 3000); keep = []
        for i,s_,c in info:
            zz,yy,xx = np.where(lab == i); P = np.stack([zz,yy,xx],1).astype(np.float32); P -= P.mean(0)
            ev = np.sort(np.linalg.eigvalsh(np.cov(P.T)))[::-1]
            if np.sqrt(ev[0]/max(ev[2],1e-6)) > 3.5: keep.append(i)
        ves = np.isin(lab, keep)
        # lumières à paroi fine de densité moyenne (vaisseaux ombilicaux, aorte, veines cardinales) : profil radial autour des cavités
        for r_ in rows:
            i = int(r_[0])
            if i in used or r_[3] < fS(0.3) or r_[1] < 800: continue
            m = cavo == i; sl = _crop(m, 14)
            if sl is None: continue
            mm = m[sl]; dsc = ds[sl]; dist = ndi.distance_transform_edt(~mm)
            prof = [float(np.median(dsc[(dist > k-1) & (dist <= k)])) if ((dist > k-1) & (dist <= k)).any() else 0.0 for k in range(1, 12)]
            wall = max(prof[2:6]); far = min(prof[7:11])
            if 20 <= wall <= 75 and far < 12 and (cart[sl] & mm).sum() == 0:
                ves[sl] |= mm | (dilate_r(mm, 4)[...] & (dsc > 15) & (dsc < 80) & ~cart[sl])
                log('  thin-walled lumen -> vessel: lab %d size %d wall %.0f' % (i, r_[1], wall))
        labels['vaisseaux'] = ves; log('vessels %.3fM' % (ves.sum()/1e6))
    except Exception as e:
        log('vessels skipped:', e)
    labels['enveloppe'] = env
    np.savez_compressed(os.path.join(work,'labels.npz'), **{k: np.packbits(v) for k,v in labels.items()}, shape=np.array(env.shape))
    json.dump({k: int(v.sum()) for k,v in labels.items()}, open(os.path.join(work,'labels_stats.json'),'w'), indent=1)
    def paint(sl):
        img = cv2.cvtColor(255-d[sl], cv2.COLOR_GRAY2BGR)
        for k,col in PALETTE.items():
            if k in labels: img[labels[k][sl]] = col[::-1]
        return img
    i = int(MID); j=int(MID-0.15*LRe); k=int(fA(0.45)); j2=int(MID+0.15*LRe)
    cv2.imwrite(os.path.join(work,'check_labels.png'), np.hstack([paint((i,)), paint((j,)), paint((j2,)), np.transpose(paint((slice(None),slice(None),k)), (1,0,2))]))
    log('segmentation done')
