"""Candidats d'artères des arcs pharyngiens et du sac aortique, à valider avant de les copier dans vaisseaux_points/<CS>.json.

Principe : une artère d'arc est une lumière pâle qui relie le sac aortique (ventral, au bout crânial de la voie de sortie du cœur) à
l'aorte dorsale du même côté, à travers le mésenchyme d'un arc pharyngien. On calcule le chemin de moindre coût depuis une graine dans
le sac aortique (coût faible dans les lumières pâles, fort dans le tissu dense ; interdit hors de l'enveloppe, dans la cavité
péricardique, les veines cardinales et le tube digestif déjà tracés). Les départs d'arcs sont cherchés dans une coque de quelques
voxels autour de la partie crâniale de chaque aorte dorsale tracée : composantes pâles du côté ventral ou latéral (les artères
intersegmentaires, dorsales, et le prolongement de l'aorte au-delà de ses extrémités tracées sont écartés). Puis, l'aorte et le reste
de la coque étant bloqués, le chemin de moindre coût depuis la graine jusqu'à chaque départ ne peut passer que par sa propre branche :
c'est l'arc. Le tronc commun à tous les chemins donne le sac aortique (avec le bout de la voie de sortie si la graine est dans le cœur).
(Un minimum du coût d'arrivée le long de l'aorte ne marche pas : sac, arcs et aorte forment une échelle, le coût y croît régulièrement.)
Le seuil de lumière est calé sur la densité mesurée au centre des aortes tracées (même type de vaisseau, même coloration).

Entrées (<work> = dossier work du stade) : dens.npy, cardio/vaisseaux_chemins.json (aortes dorsales tracées) ; facultatifs :
cardio/coeur.npz ou labels.npz (cavités cardiaques pour la graine ; enveloppe, cavité péricardique), cardio/vaisseaux.npz (veines
exclues), digestif/digestif.npz (tube digestif exclu).
Sorties : embryo3d/vaisseaux_points/<CS>_arcs_proposes.json (format de vaisseaux_points/<CS>.json + « qualite » par segment,
« a_valider ») et <work>/cardio/arcs_candidats.png (par côté : vue de profil en projection minimale de la tranche entre la ligne
médiane et l'aorte ; vue de face ; bandeau récapitulatif ; graduations en voxels de dens.npy).
Numérotation : si le nombre d'arcs propres (sans traversée de paroi) d'un côté égale le nombre d'arcs attendus au stade
(vaisseaux_points/calendrier_arcs.json : présents + en formation + en régression, sinon présents seuls), ils sont nommés
arc_aortique_<k>_<côté> dans l'ordre crânio-caudal ; sinon arc_candidat_<i>_<côté> ; ceux qui forcent une paroi sont nommés
arc_candidat_douteux_<i>_<côté>. Les arc_candidat_* sont à renommer (ou supprimer) à la validation.
Qualité : fraction pâle du chemin et plus longue traversée de tissu dense (voxels) ; une traversée > 3 voxels signale un chemin qui
force une paroi (arc interrompu, ou faux chemin par une poche ou une veine) : à vérifier sur la planche.

usage : python embryo3d/arcs_candidats.py <work> <CS> [--sac x,s,y] [--seuil N] [--marge 40] [--coque 4] [--taille-min 3] [--pas 1]
                                          [--sortie fichier.json] [--planche fichier.png]
        python embryo3d/arcs_candidats.py fusionner <CS>_arcs_proposes.json [vaisseaux_points/<CS>.json]
          → copie les segments validés (arc_aortique_*, sac_aortique, tronc_arteriel ; pas les arc_candidat_* non renommés) dans le
            fichier de points du stade, en remplaçant ceux de même nom, sans les blocs « qualite »."""
import numpy as np, json, os, sys, argparse
from scipy import ndimage as ndi

HERE = os.path.dirname(os.path.abspath(__file__))
COTES, ARCS = ('gauche', 'droite'), (1, 2, 3, 4, 6)
PALETTE = [(0, 160, 255), (255, 0, 200), (0, 200, 120), (200, 120, 0), (120, 0, 255), (0, 90, 200)]

def masque(z, k, SH):
    return np.unpackbits(z[k])[:int(np.prod(SH))].reshape(SH).astype(bool)

def union_npz(path, garder, SH):
    """union des masques d'un npz packbits dont la clé passe le filtre garder(clé) ; None si fichier absent ou aucune clé"""
    if not os.path.exists(path): return None, []
    z = np.load(path); cles = [k for k in z.files if k not in ('shape', 'union') and garder(k)]
    if not cles: return None, []
    m = np.zeros(SH, bool)
    for k in cles: m |= masque(z, k, SH)
    return m, cles

def aortes(work):
    ch = json.load(open(os.path.join(work, 'cardio', 'vaisseaux_chemins.json'), encoding='utf-8')); res = {}
    for c in COTES:
        k = next((k for k in [f'aorte_dorsale_{c}'] + (['aorte_dorsale'] if c == 'gauche' else []) if k in ch and 'centres' in ch[k]), None)
        if k is None: continue
        C, R = np.array(ch[k]['centres'], float), np.array(ch[k]['rayons'], float)
        if C[0, 1] > C[-1, 1]: C, R = C[::-1], R[::-1]                       # début = crânial (axe 1 : 0 = haut)
        res[c] = (C, R)
    return res

def graine_auto(work, SH, dens):
    """voxel le plus pâle près du centre des 8 coupes les plus crâniales des cavités cardiaques (bout de la voie de sortie)"""
    m, src = None, None
    for f, cles in (('cardio/coeur.npz', ('cavites_cardiaques', 'coeur_plein')), ('labels.npz', ('cavites_cardiaques', 'coeur_detoure', 'coeur'))):
        p = os.path.join(work, f)
        if not os.path.exists(p): continue
        z = np.load(p); k = next((k for k in cles if k in z.files), None)
        if k is None: continue
        m = masque(z, k, SH); src = f'{f}:{k}'
        if m.any(): break
    if m is None or not m.any(): return None, 'aucun masque cardiaque (donner --sac x,s,y)'
    ss = np.where(m.any(axis=(0, 2)))[0]; s0 = int(ss.min()); haut = m[:, s0:s0 + 8, :]
    c = np.round(np.argwhere(haut).mean(0) + [0, s0, 0]).astype(int)
    lo = np.maximum(c - 12, 0); hi = np.minimum(c + 13, SH)                    # lissage large (σ = 2) : minimum au centre de la lumière
    D = ndi.gaussian_filter(np.asarray(dens[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]).astype(np.float32), 2.0)
    fen = np.full(D.shape, False); w = [slice(max(c[k] - lo[k] - 4, 0), c[k] - lo[k] + 5) for k in range(3)]; fen[tuple(w)] = True
    D[~(fen & ndi.binary_dilation(m[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]], iterations=2))] = np.inf
    return np.unravel_index(np.argmin(D), D.shape) + lo, f'{src}, coupe crâniale s={s0}'

def attendus(cs, cote, chemin=os.path.join(HERE, 'vaisseaux_points', 'calendrier_arcs.json')):
    """(tous les arcs attendus au stade, arcs « présents » seuls) d'après le calendrier ; clé latéralisée prioritaire"""
    if not os.path.exists(chemin): return list(ARCS), list(ARCS)
    cal = json.load(open(chemin, encoding='utf-8')).get('stades', {}).get(cs, {})
    st = {k: cal.get(f'arc_aortique_{k}_{cote}', cal.get(f'arc_aortique_{k}')) for k in ARCS}
    return [k for k in ARCS if st[k] in ('present', 'formation', 'regression')], [k for k in ARCS if st[k] == 'present']

def prefixe(a, b, tol=1.5):
    """longueur du début du chemin a qui reste à moins de tol voxels du chemin b"""
    loin = np.linalg.norm(a[:, None].astype(float) - b[None].astype(float), axis=2).min(1) > tol
    return int(np.argmax(loin)) if loin.any() else len(a)

def echantillonne(P, ecart=6.0):
    """points de passage tous les ~ecart voxels le long du chemin (extrémités comprises), entiers"""
    d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    n = max(2, int(round(d[-1] / ecart)) + 1); u = np.linspace(0, d[-1], n)
    return np.round(np.stack([np.interp(u, d, P[:, k]) for k in range(3)], 1)).astype(int).tolist()

def qualite(vals, seuil, pas):
    dense = vals > seuil; run = best = 0
    for v in dense:
        run = run + 1 if v else 0; best = max(best, run)
    return {'fraction_pale': round(float(1 - dense.mean()), 3), 'traversee_dense_max_vox': int(best * pas),
            'densite_moyenne': round(float(vals.mean()), 1)}

def proposer(work, cs, sac=None, seuil=None, marge=40, coque_vox=4.0, taille_min=3, pas=1):
    from skimage.graph import MCP_Geometric
    dens = np.load(os.path.join(work, 'dens.npy'), mmap_mode='r'); SH = dens.shape
    aor = aortes(work)
    if not aor: sys.exit('aucune aorte dorsale dans cardio/vaisseaux_chemins.json : tracer les aortes d\'abord')
    if sac is None:
        sac, src = graine_auto(work, SH, dens)
        if sac is None: sys.exit(src)
    else: src = 'imposée (--sac)'
    sac = np.array(sac, int); print('graine du sac aortique', sac.tolist(), '(' + src + ')')
    # boîte : du haut des aortes (moins une demi-marge) à la graine (plus la marge), des aortes à la graine en x et y
    s_top = int(min(C[0, 1] for C, _ in aor.values())); s_bot = int(sac[1] + marge)
    if len(aor) == 2:                                           # arrêt au-dessus de la fusion des aortes dorsales (chemins croisés)
        (CG, RG), (CD, RD) = aor['gauche'], aor['droite']
        dist = np.linalg.norm(CG[:, None] - CD[None], axis=2); j = dist.argmin(1)
        fus = np.where(dist[np.arange(len(CG)), j] < RG + RD[j] + 2)[0]
        if len(fus): s_bot = min(s_bot, int(CG[fus[0], 1]) - 5)
    P = np.vstack([C[C[:, 1] <= s_bot] for C, _ in aor.values()] + [sac[None]])
    lo = np.maximum(np.floor(P.min(0)).astype(int) - marge, 0); lo[1] = max(s_top - marge // 2, 0)
    hi = np.minimum(np.ceil(P.max(0)).astype(int) + marge + 1, SH); hi[1] = min(max(s_bot, int(sac[1])) + 6, SH[1])
    sub = np.asarray(dens[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]).astype(np.float32)
    D = ndi.gaussian_filter(sub, max(1.0, 0.6 * pas))[::pas, ::pas, ::pas]
    loc = lambda X: np.clip(np.round((np.asarray(X, float) - lo) / pas).astype(int), 0, np.array(D.shape) - 1)
    glob = lambda I: np.asarray(I) * pas + lo
    print('boîte (x, s, y)', lo.tolist(), '→', hi.tolist(), 'pas', pas, 'voxels', D.size)
    # exclusions
    interdit = np.zeros(D.shape, bool); notes = []
    def reduit(m): return m[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]][::pas, ::pas, ::pas]
    pl = os.path.join(work, 'labels.npz')
    Lz = np.load(pl) if os.path.exists(pl) else None
    if Lz is not None and 'enveloppe' in Lz.files: interdit |= ~ndi.binary_erosion(reduit(masque(Lz, 'enveloppe', SH)), iterations=1); notes.append('hors enveloppe')
    else: print('ATTENTION : pas d\'enveloppe dans labels.npz, l\'extérieur de l\'embryon (pâle) n\'est pas exclu')
    if Lz is not None and 'cavite_pericardique' in Lz.files: interdit |= reduit(masque(Lz, 'cavite_pericardique', SH)); notes.append('cavité péricardique')
    m, cles = union_npz(os.path.join(work, 'cardio', 'vaisseaux.npz'), lambda k: not k.startswith(('aorte', 'arc', 'sac', 'tronc')), SH)
    if m is not None: interdit |= ndi.binary_dilation(reduit(m), iterations=1); notes.append('vaisseaux ' + ','.join(cles))
    m, cles = union_npz(os.path.join(work, 'digestif', 'digestif.npz'), lambda k: True, SH)
    if m is not None: interdit |= ndi.binary_dilation(reduit(m), iterations=1); notes.append('digestif ' + ','.join(cles))
    # calage du seuil sur la lumière des aortes
    centres = np.vstack([loc(C[(C[:, 1] >= lo[1]) & (C[:, 1] < hi[1])]) for C, _ in aor.values()])
    lum = float(np.median(D[tuple(centres.T)])); dedans = ~interdit
    tissu = float(np.percentile(D[dedans], 75)) if dedans.any() else float(np.percentile(D, 75))
    if seuil is None: seuil = lum + 0.35 * (tissu - lum)
    print(f'densité : lumière des aortes {lum:.0f}, tissu (75e centile) {tissu:.0f}, seuil de lumière {seuil:.0f}')
    if tissu - lum < 15: print('ATTENTION : faible contraste lumière / tissu, candidats peu fiables')
    cout = 1.0 + 400.0 * np.clip((D - seuil) / max(tissu - seuil, 1.0), 0, 1) ** 2
    # centralité : coût multiplié par 1 + 3 / (distance à la paroi), maximal (×7) hors lumière → chemins sur la ligne centrale, confondus
    # tant qu'ils partagent un tronc (sinon, dans une lumière large, chacun coupe au plus court et le tronc commun disparaît)
    cout *= 1.0 + 3.0 / np.maximum(ndi.distance_transform_edt((D <= seuil) & ~interdit), 0.5)
    cout[interdit] = np.inf
    # tube de chaque aorte dorsale (partie crâniale) : axe rastérisé, distance à l'axe, point d'axe le plus proche
    Qg, Rl, Tl, cote_id = [], [], [], []
    for c, (C, R) in aor.items():
        garde = (C[:, 1] >= lo[1]) & (C[:, 1] <= min(s_bot, hi[1] - 1)); C, R = C[garde], R[garde]
        if len(C) < 3: print(c, ': aorte trop courte dans la boîte'); continue
        d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(C, axis=0), axis=1))]; u = np.arange(0, d[-1], 0.5 * pas)
        Q = np.stack([np.interp(u, d, C[:, k]) for k in range(3)], 1); T = np.gradient(Q, axis=0)
        Qg.append(Q); Rl.append(np.interp(u, d, R) / pas); Tl.append(T / (np.linalg.norm(T, axis=1, keepdims=True) + 1e-9)); cote_id += [c] * len(Q)
    if not Qg: sys.exit('aucune aorte dorsale utilisable dans la boîte')
    Qg, Rl, Tl = np.concatenate(Qg), np.concatenate(Rl), np.concatenate(Tl); Ql = (Qg - lo) / pas
    axe_id = np.full(D.shape, -1, np.int32); axe_id[tuple(loc(Qg).T)] = np.arange(len(Qg))
    ind = np.zeros((3,) + D.shape, np.int32)
    dist = ndi.distance_transform_edt(axe_id < 0, return_indices=True, indices=ind)
    proche = axe_id[ind[0], ind[1], ind[2]]; del ind
    Rn = Rl[proche]; lumiere_ao = dist <= Rn + 1
    coque = (dist > Rn + 1) & (dist <= Rn + 1 + coque_vox / pas)
    # départs de branches : composantes pâles de la coque, côté ventral ou latéral (pas les intersegmentaires dorsales, pas le
    # prolongement de l'aorte au-delà de ses extrémités tracées)
    lab, n = ndi.label(coque & (D <= seuil) & ~interdit, structure=np.ones((3, 3, 3)))
    departs, gardees = [], np.zeros(D.shape, bool)
    for i, sl in enumerate(ndi.find_objects(lab), 1):
        if sl is None: continue
        vox = np.argwhere(lab[sl] == i) + [x.start for x in sl]
        if len(vox) * pas ** 3 < taille_min: continue
        j = int(proche[tuple(vox[np.argmin(dist[tuple(vox.T)])])]); o = vox.mean(0) - Ql[j]; T = Tl[j]
        le_long = abs(o @ T) / (np.linalg.norm(o) + 1e-9); op = o - (o @ T) * T; ny = op[2] / (np.linalg.norm(op) + 1e-9)
        if le_long > 0.7 or ny > 0.25: continue
        departs.append({'cote': cote_id[j], 'j': j, 'vox': vox, 'ventral': round(float(-ny), 2)}); gardees[tuple(vox.T)] = True
    if not departs: sys.exit('aucun départ de branche pâle ventral le long des aortes : vérifier le seuil (--seuil), la coque (--coque) et la boîte (--marge)')
    # chemins depuis le sac, aortes et reste de la coque bloqués : chaque départ n'est joignable que par sa propre branche
    cout[lumiere_ao | (coque & ~gardees)] = np.inf
    g = tuple(loc(sac))
    if not np.isfinite(cout[g]) or D[g] > seuil:                            # graine hors lumière ou dans une zone exclue : lumière la plus proche
        libre = np.argwhere(np.isfinite(cout) & (D <= seuil))
        if not len(libre): sys.exit('aucune lumière accessible dans la boîte')
        g = tuple(libre[np.argmin(np.linalg.norm(libre - np.array(g), axis=1))]); sac = glob(g).astype(int)
        print('graine déplacée dans la lumière la plus proche :', sac.tolist(), '(vérifier sur la planche, sinon --sac)')
    mcp = MCP_Geometric(cout, fully_connected=True); A, _ = mcp.find_costs([g])
    cands = []
    for k in departs:
        a = A[tuple(k['vox'].T)]
        if not np.isfinite(a).any(): print(f"  départ {k['cote']} s={int(Qg[k['j'], 1])} injoignable depuis le sac (ignoré)"); continue
        e = tuple(k['vox'][np.argmin(a)]); k.update(cout=float(A[e]), chemin=np.array(mcp.traceback(e)), s=float(Qg[k['j'], 1]))
        cands.append(k)
    # doublons (plusieurs départs d'une même branche) : même côté, jonctions à moins de max(6, 2 rayons) voxels et chemins confondus
    # sur plus de 60 % → on garde le moins coûteux
    retenus = []
    for k in sorted(cands, key=lambda k: k['cout']):
        dup = next((r for r in retenus if r['cote'] == k['cote']
                    and np.linalg.norm(Qg[k['j']] - Qg[r['j']]) <= max(6.0, 2 * Rl[k['j']] * pas)
                    and np.mean(np.min(np.linalg.norm(k['chemin'][len(k['chemin']) // 3:, None] - r['chemin'][None], axis=2), axis=1) <= 2) > 0.6), None)
        if dup is None: k['fusionnes'] = 0; retenus.append(k)
        else: dup['fusionnes'] += 1
    # tronc commun (sac aortique), calculé sur les seuls chemins propres (un chemin qui force une paroi le raccourcirait), puis arcs
    for k in retenus: k['propre'] = qualite(D[tuple(k['chemin'].T)], seuil, pas)['traversee_dense_max_vox'] <= 3
    propres_tous = [k for k in retenus if k['propre']]
    for k in retenus:                                                          # sortie du tronc commun, chemin par chemin
        k['k_sac'] = min((prefixe(k['chemin'], r['chemin']) for r in propres_tous if r is not k), default=0)
    ref = propres_tous[0] if len(propres_tous) > 1 else None; k_sac = ref['k_sac'] if ref else 0
    segments = []
    Rao = float(np.median(np.concatenate([R for _, R in aor.values()])))
    if k_sac >= 3:
        Pg = glob(ref['chemin'][:k_sac])
        segments.append({'nom': 'sac_aortique', 'rmax': int(np.clip(round(1.3 * Rao), 4, 8)), 'rayon_recalage': 3, 'points': echantillonne(Pg),
                         'qualite': qualite(D[tuple(ref['chemin'][:k_sac].T)], seuil, pas), 'note': 'tronc commun des chemins (du cœur vers les arcs) ; '
                         'contient le bout de la voie de sortie si la graine est dans le cœur : couper à la validation'})
    for k in retenus:
        k['seg'] = k['chemin'][max(k['k_sac'] - 1, 0):]                         # départ : bout du tronc commun (sac)
        k['q'] = qualite(D[tuple(k['seg'].T)], seuil, pas)
    for c in COTES:
        ks = sorted((k for k in retenus if k['cote'] == c), key=lambda k: k['s'])
        if not ks: continue
        tous_c, presents = attendus(cs, c); propres = [k for k in ks if k['propre']]
        if len(propres) == len(tous_c): arcs, regle = tous_c, 'nombre = arcs attendus (présents, formation, régression)'
        elif len(propres) == len(presents): arcs, regle = presents, 'nombre = arcs « présents » du calendrier'
        else: arcs, regle = None, f'{len(propres)} branches propres pour {len(tous_c)} arcs attendus ({tous_c}) : à nommer'
        it = iter(arcs or []); nd = 0; noms = []
        for i, k in enumerate(ks):
            if not k['propre']: nd += 1; noms.append(f'arc_candidat_douteux_{nd}_{c}')
            else: noms.append(f'arc_aortique_{next(it)}_{c}' if arcs else f'arc_candidat_{i + 1}_{c}')
        print(f'{c} : {len(ks)} branche(s), {len(propres)} propre(s) — {regle}')
        for nom, k in zip(noms, ks):
            seg, q = k['seg'], k['q']
            Pg = np.vstack([glob(seg).astype(float), Qg[k['j']][None]])          # arrivée : axe de l'aorte dorsale
            q.update({'jonction_aortique': np.round(Qg[k['j']]).astype(int).tolist(), 'cout': round(k['cout'], 1), 'ventral': k['ventral'],
                      'longueur_vox': round(float(np.linalg.norm(np.diff(Pg, axis=0), axis=1).sum()), 1), 'departs_fusionnes': k['fusionnes']})
            segments.append({'nom': nom, 'rmax': int(np.clip(round(Rao), 3, 6)), 'rayon_recalage': 2, 'points': echantillonne(Pg), 'qualite': q})
            alerte = ' ← traversée de paroi' if q['traversee_dense_max_vox'] > 3 else ''
            print(f"  {nom:26s} jonction s={q['jonction_aortique'][1]:4d}  pâle {q['fraction_pale']:.2f}  paroi {q['traversee_dense_max_vox']:2d} vox"
                  f"  longueur {q['longueur_vox']:5.0f}  coût {q['cout']:7.0f}  ventral {q['ventral']:+.2f}{alerte}")
    res = {'stade': cs, 'a_valider': True, 'propose_par': 'arcs_candidats.py',
           'validation': 'vérifier chaque chemin sur la planche, renommer les arc_candidat_* (arc_aortique_<k>_<côté>) ou les supprimer, corriger '
                         'les points au besoin (coordonnées de dens.npy : x, s, y), couper le sac aortique au bout de la voie de sortie, puis '
                         '« arcs_candidats.py fusionner » et digestif_build.py … --sortie cardio',
           'parametres': {'graine_sac': sac.tolist(), 'source_graine': src, 'seuil': round(seuil, 1), 'lumiere_aortes': round(lum, 1), 'tissu': round(tissu, 1),
                          'marge': marge, 'coque': coque_vox, 'taille_min': taille_min, 'pas': pas, 'boite': [lo.tolist(), hi.tolist()], 'exclusions': notes},
           'segments': segments}
    return res, dict(D=D, lo=lo, pas=pas, aor=aor, sac=np.asarray(sac))

# ---------------------------------------------------------------- planche de contrôle
def planche(res, ctx, chemin):
    import cv2
    sys.path.insert(0, HERE)
    try:
        from tubes_morph import couleur
    except Exception:
        couleur = None
    D, lo, pas, aor, sac = ctx['D'], ctx['lo'], ctx['pas'], ctx['aor'], ctx['sac']
    z = max(1, int(round(2 / pas)))                                     # agrandissement
    def col(nom, i):
        base = nom.rsplit('_', 1)[0]
        if couleur is not None and not nom.startswith('arc_candidat'):
            r, g, b = couleur(base if base.startswith('arc_aortique') else nom); return (int(255 * b), int(255 * g), int(255 * r))
        return PALETTE[i % len(PALETTE)]
    def grad(im, ax_l, ax_c):
        """graduations tous les 50 voxels de dens.npy (lignes : axe ax_l, colonnes : axe ax_c)"""
        h, w = im.shape[:2]
        for ax, n, horiz in ((ax_l, h, False), (ax_c, w, True)):
            for v in range(int(np.ceil(lo[ax] / 50.0)) * 50, int(lo[ax] + n / z * pas) + 1, 50):
                p = int((v - lo[ax]) / pas * z)
                if horiz: cv2.line(im, (p, 0), (p, 6), (0, 150, 0), 1); cv2.putText(im, str(v), (p + 2, 16), 0, 0.35, (0, 120, 0), 1)
                else: cv2.line(im, (0, p), (6, p), (0, 150, 0), 1); cv2.putText(im, str(v), (8, p + 4), 0, 0.35, (0, 120, 0), 1)
        return im
    def vue(proj, ax_l, ax_c, titre, cote=None):
        im = cv2.cvtColor(np.clip(255 - proj, 0, 255).astype(np.uint8), cv2.COLOR_GRAY2BGR)
        im = cv2.resize(im, None, fx=z, fy=z, interpolation=cv2.INTER_NEAREST)
        pt = lambda X: (int((X[ax_c] - lo[ax_c]) / pas * z), int((X[ax_l] - lo[ax_l]) / pas * z))
        for c, (C, _) in aor.items(): cv2.polylines(im, [np.array([pt(x) for x in C], np.int32)], False, (40, 40, 220), 1, cv2.LINE_AA)
        for i, s in enumerate(res['segments']):
            if cote and s['nom'].endswith(('_gauche', '_droite')) and not s['nom'].endswith(cote): continue
            P = np.array(s['points']); cc = col(s['nom'], i)
            cv2.polylines(im, [np.array([pt(x) for x in P], np.int32)], False, cc, 2, cv2.LINE_AA)
            for x in P: cv2.circle(im, pt(x), 2, cc, -1)
            lab = s['nom'].replace('arc_aortique_', 'arc ').replace('arc_candidat_douteux_', '?? ').replace('arc_candidat_', '? ').replace('_gauche', ' G').replace('_droite', ' D').replace('sac_aortique', 'sac')
            cv2.putText(im, lab, (pt(P[-1])[0] + 4, pt(P[-1])[1] + 4), 0, 0.4, cc, 1, cv2.LINE_AA)
        cv2.circle(im, pt(sac), 4, (0, 140, 255), -1)
        cv2.putText(im, titre, (8, im.shape[0] - 8), 0, 0.45, (200, 0, 0), 1, cv2.LINE_AA)
        return grad(im, ax_l, ax_c)
    xs, ys = (sac[0] - lo[0]) / pas, (sac[2] - lo[2]) / pas; tuiles = []
    for c, (C, _) in aor.items():                                       # profil : tranche de la ligne médiane à l'aorte du côté (+ marge)
        a, b = sorted((xs, (np.median(C[:, 0]) - lo[0]) / pas)); a, b = int(max(a - 6 / pas, 0)), int(min(b + 6 / pas, D.shape[0] - 1)) + 1
        proj = D[a:b].min(axis=0)                                       # (s, y) : lignes = haut → bas, colonnes = ventral → dorsal
        tuiles.append(vue(proj, 1, 2, f'profil {c} (x {int(a * pas + lo[0])}-{int(b * pas + lo[0])}, projection minimale)', c))
    ya = max((np.median(C[:, 2]) - lo[2]) / pas for C, _ in aor.values())
    proj = D[:, :, int(max(min(ys, ya) - 3, 0)):int(max(ys, ya)) + 4].min(axis=2).T   # (s, x) : lignes = haut → bas, colonnes = gauche → droite
    tuiles.append(vue(proj, 1, 0, 'face (projection minimale, du sac aux aortes)'))
    H = max(t.shape[0] for t in tuiles)
    tuiles = [cv2.copyMakeBorder(t, 0, H - t.shape[0], 0, 6, cv2.BORDER_CONSTANT, value=(90, 90, 90)) for t in tuiles]
    haut = np.hstack(tuiles)
    # bandeau : une ligne par segment proposé (nom, jonction, qualité)
    lignes = [f"{res['stade']} - graine du sac {ctx['sac'].tolist()} - seuil de lumiere {res['parametres']['seuil']:.0f} - A VALIDER (points en voxels de dens.npy : x, s, y)"]
    for sg in res['segments']:
        q = sg['qualite']; j = q.get('jonction_aortique')
        lignes.append(f"{sg['nom']:26s}" + (f" jonction s={j[1]:4d}" if j else ' ' * 16) + f"  pale {q['fraction_pale']:.2f}  paroi {q['traversee_dense_max_vox']} vox"
                      + ('  <- traversee de paroi, a verifier' if q['traversee_dense_max_vox'] > 3 else ''))
    g = np.full((16 * len(lignes) + 10, haut.shape[1], 3), 255, np.uint8)
    for i, t in enumerate(lignes): cv2.putText(g, t, (8, 16 + 16 * i), 0, 0.42, (40, 40, 40) if i else (160, 0, 0), 1, cv2.LINE_AA)
    os.makedirs(os.path.dirname(os.path.abspath(chemin)), exist_ok=True)
    cv2.imwrite(chemin, np.vstack([haut, g])); print('planche', chemin)

# ---------------------------------------------------------------- fusion des propositions validées
def fusionner(proposes, cible=None):
    P = json.load(open(proposes, encoding='utf-8'))
    cible = cible or os.path.join(HERE, 'vaisseaux_points', f"{P['stade']}.json")
    T = json.load(open(cible, encoding='utf-8'))
    ok = [s for s in P['segments'] if s['nom'].startswith(('arc_aortique_', 'sac_aortique', 'tronc_arteriel'))]
    ignores = [s['nom'] for s in P['segments'] if s not in ok]
    noms = {s['nom'] for s in ok}
    T['segments'] = [s for s in T['segments'] if s['nom'] not in noms] + [{k: v for k, v in s.items() if k not in ('qualite', 'note')} for s in ok]
    json.dump(T, open(cible, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    print('fusionnés dans', cible, ':', ', '.join(sorted(noms)) or '(rien)')
    if ignores: print('non fusionnés (à renommer ou supprimer) :', ', '.join(ignores))

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'fusionner':
        fusionner(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None); sys.exit()
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('work'); ap.add_argument('cs', help='stade (CS13, CS14…), pour le calendrier et le nom de sortie')
    ap.add_argument('--sac', help='graine x,s,y (voxels de dens.npy) dans la lumière du sac aortique ; défaut : haut des cavités cardiaques')
    ap.add_argument('--seuil', type=float, help='densité maximale d\'une lumière ; défaut : calé sur les aortes')
    ap.add_argument('--marge', type=int, default=40, help='marge de la boîte autour des aortes et de la graine (voxels)')
    ap.add_argument('--coque', type=float, default=4.0, help='épaisseur (voxels) de la coque autour des aortes où l\'on cherche les départs de branches')
    ap.add_argument('--taille-min', type=int, default=3, help='taille minimale (voxels) d\'un départ de branche dans la coque')
    ap.add_argument('--pas', type=int, default=1, help='sous-échantillonnage du volume (2 = plus rapide)')
    ap.add_argument('--sortie'); ap.add_argument('--planche')
    a = ap.parse_args()
    sac = [int(v) for v in a.sac.split(',')] if a.sac else None
    res, ctx = proposer(a.work, a.cs, sac, a.seuil, a.marge, a.coque, a.taille_min, max(1, a.pas))
    sortie = a.sortie or os.path.join(HERE, 'vaisseaux_points', f'{a.cs}_arcs_proposes.json')
    json.dump(res, open(sortie, 'w', encoding='utf-8'), indent=1, ensure_ascii=False); print('propositions', sortie)
    planche(res, ctx, a.planche or os.path.join(a.work, 'cardio', 'arcs_candidats.png'))
