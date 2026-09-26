"""Récupère les modèles 3D des PDF interactifs du système des arcs aortiques (Rana et al. 2014, équipe Moorman,
AMC Amsterdam ; ou tout autre PDF 3D U3D / PRC, p. ex. Hikspoors, 3D Atlas) et les range par stade :
maillages PLY + GLB, manifeste lisible par blender_build_scene.py, identification des structures (légende),
rapport et aperçu.

usage :
  python embryo3d/reference/arcs_aortiques/extraire.py                         # tous les PDF de reference/arcs_aortiques/pdf/
  python embryo3d/reference/arcs_aortiques/extraire.py mon.pdf autre.pdf --sortie extraits
  python embryo3d/reference/arcs_aortiques/extraire.py x.pdf --stade CS14       # impose le stade (un seul modèle)
  python embryo3d/reference/arcs_aortiques/extraire.py x.pdf --lister           # inventaire des flux 3D sans décoder
  options : --correspondance corr.json (noms d'origine -> structure, à partir de correspondance_proposee.json)
            --stades stades.json ({"<sha1 du flux ou fichier#k>": "CS14"}) --mm-par-unite 0.001 --sans-centrage --sans-apercu
Blender : blender -b -P embryo3d/blender_build_scene.py -- extraits/CS14/manifest.json arcs_CS14.blend
"""
import argparse
import json
import os
import re
import sys
import unicodedata

import numpy as np

ICI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ICI)
from pdf3d import extraire_pdf  # noqa: E402
from prc import ErreurPRC, lire_prc  # noqa: E402
from u3d import ErreurU3D, lire_u3d  # noqa: E402

ROMAINS = {'i': 1, 'ii': 2, 'iii': 3, 'iv': 4, 'v': 5, 'vi': 6}
RE_STADE = re.compile(r'(?<![a-z])(?:cs|carnegie[\s_-]*stage|stage|stadium|stade)[\s_.-]*(\d{1,2})(?!\d)', re.I)
COULEURS_LEGENDE = ['tube_digestif', 'voie_de_sortie', 'sac_aortique', 'aorte_dorsale', 'arc_1', 'arc_2', 'arc_3', 'arc_4',
                    'arc_6', 'tube_neural', 'otocyste']


def sans_accents(s):
    return ''.join(c for c in unicodedata.normalize('NFD', s) if unicodedata.category(c) != 'Mn')


def normaliser(nom):
    s = sans_accents(nom or '').lower()
    s = re.sub(r'[_\-.:/\\]+', ' ', s)
    s = re.sub(r'(\d)([a-z])', r'\1 \2', s)
    s = re.sub(r'([a-z])(\d)', r'\1 \2', s)
    return ' ' + re.sub(r'\s+', ' ', s).strip() + ' '


def charger_legende(chemin=None):
    return json.load(open(chemin or os.path.join(ICI, 'legende.json'), encoding='utf-8'))


def lire_cote(n, L):
    g = any(re.search(m, n) for m in L['cotes']['gauche'])
    d = any(re.search(m, n) for m in L['cotes']['droite'])
    return 'gauche' if g and not d else ('droite' if d and not g else None)


def lire_arc(n, L):
    A = L['motifs_arcs']
    if not any(re.search(m, n) for m in A['marqueurs']):
        return None
    for mot, v in A['ordinaux'].items():
        if re.search(r'\b' + re.escape(sans_accents(mot)) + r'\b', n):
            return v
    m = re.search(r'\b([1-6])\b', n)
    if m:
        return int(m.group(1))
    for r, v in sorted(ROMAINS.items(), key=lambda kv: -len(kv[0])):
        if re.search(r'\b' + r + r'\b', n):
            return v
    return None


def distance_couleur(a, b):
    return float(np.linalg.norm(np.asarray(a[:3], float) - np.asarray(b[:3], float)))


def identifier(nom, couleur, L, corr=None):
    """-> {'structure', 'cote', 'methode', 'detail'} ; methode : correspondance | nom | couleur | inconnu."""
    if corr and nom in corr:
        c = corr[nom]
        c = {'structure': c} if isinstance(c, str) else dict(c)
        c.setdefault('cote', None)
        return {**c, 'methode': 'correspondance', 'detail': 'fichier de correspondance'}
    n = normaliser(nom)
    cote = lire_cote(n, L)
    arc = lire_arc(n, L)
    if arc is not None:
        if arc == 5:
            return {'structure': 'arc_5', 'cote': cote, 'methode': 'nom', 'detail': '5e arc nommé (inhabituel chez l\'homme)'}
        return {'structure': f'arc_{arc}', 'cote': cote, 'methode': 'nom', 'detail': f'artère d\'arc n° {arc}'}
    S = L['structures']
    for cle in L.get('ordre_de_test', list(S)):
        for m in S[cle].get('motifs', []):
            if re.search(m, n):
                return {'structure': cle, 'cote': cote, 'methode': 'nom', 'detail': f'motif « {m} »'}
    if couleur is not None:
        best = min(COULEURS_LEGENDE, key=lambda k: distance_couleur(couleur, S[k]['couleur']))
        d = distance_couleur(couleur, S[best]['couleur'])
        if d <= L.get('tolerance_couleur', 0.22):
            return {'structure': best, 'cote': cote, 'methode': 'couleur', 'detail': f'couleur la plus proche (écart {d:.2f})'}
    return {'structure': None, 'cote': cote, 'methode': 'inconnu', 'detail': ''}


def candidats_stade(*textes):
    out = []
    for t in textes:
        for m in RE_STADE.finditer(t or ''):
            v = int(m.group(1))
            if 7 <= v <= 23:
                out.append(f'CS{v}')
    return out


def choisir_stade(flux, res_pdf, objets, nom_fichier, force=None, table=None, k=0):
    if force:
        return force, 'imposé (--stade)'
    if table:
        for cle in (flux['sha1'], f'{os.path.basename(nom_fichier)}#{k}', os.path.basename(nom_fichier)):
            if cle in table:
                return table[cle], f'table des stades ({cle})'
    essais = [('titre de l\'annotation', [flux.get('titre', '')]),
              ('vues', [v['nom'] or '' for v in flux.get('vues', [])]),
              ('noms des nœuds 3D', [o['nom'] for o in objets] + [g for o in objets for g in o.get('groupes', [])]),
              ('nom du fichier', [os.path.basename(nom_fichier)])]
    if flux.get('page'):
        essais.insert(1, ('texte de la page', [res_pdf['texte_pages'][flux['page'] - 1]] if flux['page'] <= len(res_pdf['texte_pages']) else []))
    for quoi, textes in essais:
        c = candidats_stade(*textes)
        if c:
            uniques = sorted(set(c), key=c.index)
            if len(uniques) == 1:
                return uniques[0], quoi
            maj = max(uniques, key=c.count)
            if c.count(maj) > len(c) / 2:
                return maj, f'{quoi} (majoritaire parmi {uniques})'
    return None, 'aucun indice de stade'


def facteur_mm(objets, unites_u3d=None, force=None):
    if force:
        return float(force), 'imposé (--mm-par-unite)'
    V = np.vstack([o['sommets'] for o in objets]) if objets else np.zeros((1, 3))
    L = float(np.ptp(V, axis=0).max()) if len(V) else 0.0
    if unites_u3d:
        f = unites_u3d * 1000.0
        if 0.2 <= L * f <= 30:
            return f, f'en-tête U3D (1 unité = {unites_u3d} m)'
    for f, quoi in ((1.0, 'mm'), (1e-3, 'µm'), (1e3, 'm'), (10.0, 'cm')):
        if 0.2 <= L * f <= 30:
            return f, f'déduit de l\'étendue ({L:.4g} unités, supposées en {quoi})'
    return 1.0, f'INDÉTERMINÉ : étendue {L:.4g} unités (laissé tel quel, préciser --mm-par-unite)'


def ecrire_ply(chemin, V, F, rgb=None, segments=None):
    V = np.asarray(V, np.float32)
    F = np.asarray(F, np.int32).reshape(-1, 3)
    tete = ['ply', 'format binary_little_endian 1.0', 'comment embryo3d reference/arcs_aortiques',
            f'element vertex {len(V)}', 'property float x', 'property float y', 'property float z']
    if rgb is not None:
        tete += ['property uchar red', 'property uchar green', 'property uchar blue']
    tete += [f'element face {len(F)}', 'property list uchar int vertex_indices']
    if segments is not None and len(segments):
        tete += [f'element edge {len(segments)}', 'property int vertex1', 'property int vertex2']
    tete += ['end_header']
    with open(chemin, 'wb') as fh:
        fh.write(('\n'.join(tete) + '\n').encode('ascii'))
        if rgb is None:
            fh.write(V.tobytes())
        else:
            c = np.clip(np.asarray(rgb, float) * 255 + 0.5, 0, 255).astype(np.uint8)
            dt = np.dtype([('p', '<f4', 3), ('c', 'u1', 3)])
            a = np.empty(len(V), dt)
            a['p'] = V
            a['c'] = c if c.ndim == 2 else np.tile(c, (len(V), 1))
            fh.write(a.tobytes())
        dt = np.dtype([('n', 'u1'), ('i', '<i4', 3)])
        a = np.empty(len(F), dt)
        a['n'] = 3
        a['i'] = F
        fh.write(a.tobytes())
        if segments is not None and len(segments):
            fh.write(np.asarray(segments, '<i4').tobytes())


def ecrire_glb(chemin, objets_sortie):
    try:
        import trimesh
    except ImportError:
        return False
    sc = trimesh.Scene()
    for o in objets_sortie:
        if not len(o['F']):
            continue                                  # jeux de lignes : dans les PLY seulement
        m = trimesh.Trimesh(o['V'], o['F'], process=False)
        c = (np.array(list(o['couleur']) + [o['alpha']]) * 255).astype(np.uint8)
        m.visual = trimesh.visual.ColorVisuals(m, face_colors=np.tile(c, (len(o['F']), 1)))
        sc.add_geometry(m, node_name=o['nom'], geom_name=o['nom'])
    if not sc.geometry:
        return False
    try:
        sc.export(chemin)
    except Exception as e:                            # l'export GLB est un plus : les PLY restent la référence
        print('GLB non écrit :', e)
        return False
    return True


def apercu(chemin, objets_sortie, titre):
    """Planche de 3 projections orthographiques (face, profil, dessus), rendu par peintre (OpenCV)."""
    try:
        import cv2
    except ImportError:
        return False
    if not objets_sortie:
        return False
    V = np.vstack([o['V'] for o in objets_sortie])
    lo, hi = V.min(0), V.max(0)
    taille, marge = 420, 20
    vues = [('face (X, Y)', (0, 1, 2)), ('profil (Z, Y)', (2, 1, 0)), ('dessus (X, Z)', (0, 2, 1))]
    planche = np.full((taille + 60, 3 * taille, 3), 255, np.uint8)
    for k, (nom_vue, (a, b, p)) in enumerate(vues):
        ech = (taille - 2 * marge) / max(hi[a] - lo[a], hi[b] - lo[b], 1e-9)
        img = np.full((taille, taille, 3), 255, np.uint8)
        tris, prof, coul = [], [], []
        for o in objets_sortie:
            F = o['F']
            if not len(F):
                continue
            if len(F) > 60000:
                F = F[np.linspace(0, len(F) - 1, 60000).astype(int)]
            P = o['V'][F]
            n = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
            n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-12
            lum = 0.35 + 0.65 * np.abs(n[:, p])
            c = np.array(o['couleur'])[None, :] * lum[:, None]
            xy = np.stack([(P[:, :, a] - lo[a]) * ech + marge, taille - ((P[:, :, b] - lo[b]) * ech + marge)], -1)
            tris.append(xy)
            prof.append(P[:, :, p].mean(1))
            coul.append(c)
        if tris:
            tris, prof, coul = np.concatenate(tris), np.concatenate(prof), np.concatenate(coul)
            for i in np.argsort(prof):
                cv2.fillConvexPoly(img, np.round(tris[i]).astype(np.int32), tuple(int(255 * x) for x in coul[i][::-1]), lineType=cv2.LINE_AA)
        for o in objets_sortie:
            for s0, s1 in o.get('segments', []):
                pa = [int((o['V'][s0, a] - lo[a]) * ech + marge), int(taille - ((o['V'][s0, b] - lo[b]) * ech + marge))]
                pb = [int((o['V'][s1, a] - lo[a]) * ech + marge), int(taille - ((o['V'][s1, b] - lo[b]) * ech + marge))]
                cv2.line(img, pa, pb, tuple(int(255 * x) for x in o['couleur'][::-1]), 2, cv2.LINE_AA)
        planche[40:40 + taille, k * taille:(k + 1) * taille] = img
        cv2.putText(planche, nom_vue, (k * taille + 10, 30), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (60, 60, 60), 1, cv2.LINE_AA)
    cv2.putText(planche, sans_accents(titre)[:90], (10, taille + 55), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 0, 0), 1, cv2.LINE_AA)
    cv2.imwrite(chemin, planche)
    return True


def decoder(flux):
    if flux['format'] == 'U3D':
        sc = lire_u3d(flux['donnees'])
        return sc.objets(), sc.resume(), sc.unites, sc.erreurs
    if flux['format'] == 'PRC':
        sc = lire_prc(flux['donnees'])
        return sc.objets(), sc.resume(), None, sc.erreurs
    raise ValueError(f'format {flux["format"]} inconnu')


def slug(s):
    s = re.sub(r'[^A-Za-z0-9]+', '_', sans_accents(s or '')).strip('_')
    return s[:40] or 'sans_nom'


def traiter(fichiers, sortie, L, stade_force=None, corr=None, table_stades=None, mm_force=None, centrer=True, avec_apercu=True, lister=False):
    os.makedirs(sortie, exist_ok=True)
    index = {'sortie': os.path.abspath(sortie), 'fichiers': []}
    rapport = ['# Extraction des modèles 3D (arcs aortiques)', '']
    deja = {}
    for fichier in fichiers:
        res = extraire_pdf(fichier)
        ent = {'fichier': os.path.basename(fichier), 'pages': res['pages'], 'erreurs': res['erreurs'], 'flux': []}
        index['fichiers'].append(ent)
        rapport += [f'## {os.path.basename(fichier)}', '', f'{res["pages"]} page(s), {len(res["flux3d"])} flux 3D.', '']
        rapport += [f'- ⚠ {e}' for e in res['erreurs']]
        for k, flux in enumerate(res['flux3d']):
            info = {'index': k, 'format': flux['format'], 'octets': flux['octets'], 'sha1': flux['sha1'], 'page': flux['page'],
                    'source': flux['source'], 'titre': flux.get('titre', ''), 'vues': [v['nom'] for v in flux['vues']]}
            ent['flux'].append(info)
            if lister:
                rapport.append(f'- flux {k} : {flux["format"]}, {flux["octets"]} octets, page {flux["page"]}, vues {info["vues"]}')
                continue
            if flux['sha1'] in deja:
                info['doublon_de'] = deja[flux['sha1']]
                continue
            try:
                objets, resume, unites, erreurs = decoder(flux)
            except (ErreurU3D, ErreurPRC, ValueError) as e:
                info['erreur'] = str(e)
                rapport.append(f'- flux {k} ({flux["format"]}) : **illisible** — {e}')
                continue
            except Exception as e:                     # un flux défectueux ne doit pas arrêter les autres
                info['erreur'] = f'{type(e).__name__} : {e}'
                rapport.append(f'- flux {k} ({flux["format"]}) : **illisible** — {e}')
                continue
            stade, pourquoi = choisir_stade(flux, res, objets, fichier, stade_force, table_stades, k)
            dossier_nom = stade or f'inconnu_{flux["sha1"][:8]}'
            if dossier_nom in set(deja.values()):          # deux modèles du même stade (autre PDF, autre série)
                dossier_nom = f'{dossier_nom}_{flux["sha1"][:6]}'
            deja[flux['sha1']] = dossier_nom
            info.update({'stade': stade, 'stade_indice': pourquoi, 'dossier': dossier_nom, 'resume_decodage': resume, 'erreurs_decodage': erreurs})
            d = os.path.join(sortie, dossier_nom)
            os.makedirs(os.path.join(d, 'brut'), exist_ok=True)
            for ancien in os.listdir(d):                 # maillages d'une extraction précédente (noms éventuellement changés)
                if ancien.startswith(dossier_nom + '_') and ancien.endswith(('.ply', '.glb')):
                    os.remove(os.path.join(d, ancien))
            ext = '.u3d' if flux['format'] == 'U3D' else '.prc'
            open(os.path.join(d, 'brut', f'{slug(os.path.splitext(os.path.basename(fichier))[0])}_{k}{ext}'), 'wb').write(flux['donnees'])
            json.dump(flux['vues'], open(os.path.join(d, 'vues.json'), 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
            if flux.get('javascript') or res['javascript_document']:
                with open(os.path.join(d, 'javascript.js'), 'w', encoding='utf-8') as fh:
                    if flux.get('javascript'):
                        fh.write('// --- OnInstantiate de l\'annotation 3D\n' + flux['javascript'] + '\n')
                    for js in res['javascript_document']:
                        fh.write(f'// --- JavaScript du document : {js["nom"]}\n{js["code"]}\n')
            if flux['page'] and flux['page'] <= len(res['texte_pages']):
                open(os.path.join(d, 'texte_page.txt'), 'w', encoding='utf-8').write(res['texte_pages'][flux['page'] - 1])
            fmm, pourquoi_mm = facteur_mm(objets, unites, mm_force)
            V_all = np.vstack([o['sommets'] for o in objets]) * fmm if objets else np.zeros((1, 3))
            centre = (V_all.min(0) + V_all.max(0)) / 2 if centrer else np.zeros(3)
            prefixe = dossier_nom
            noms_pris, structures, sorties, corr_prop = {}, [], [], {}
            for o in objets:
                ident = identifier(o['noeud'], o['couleur'], L, corr)
                cle = ident['structure']
                base = cle or f'inconnu_{slug(o["noeud"])}'
                if ident.get('cote'):
                    base += '_' + ident['cote']
                n = noms_pris.get(base, 0)
                noms_pris[base] = n + 1
                nom = base if n == 0 else f'{base}_{n + 1}'
                V = o['sommets'] * fmm - centre
                F = o['faces']
                Sleg = L['structures'].get(cle or '', {})
                couleur_leg = Sleg.get('couleur')
                couleur = list(o['couleur']) if o['couleur'] is not None else (couleur_leg or [0.7, 0.7, 0.7])
                alpha = float(o.get('opacite') or 1.0)
                fichier_ply = f'{prefixe}_{nom}.ply'
                seg = o.get('segments')
                ecrire_ply(os.path.join(d, fichier_ply), V, F,
                           o['couleurs_sommets'][:, :3] if o.get('couleurs_sommets') is not None else couleur, seg)
                aire = float(np.linalg.norm(np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]]), axis=1).sum() / 2) if len(F) else 0.0
                structures.append({'name': nom, 'file': fichier_ply, 'collection': Sleg.get('collection', 'Non identifie'),
                                   'color': [round(float(x), 4) for x in couleur[:3]], 'alpha': round(alpha, 3),
                                   'structure': cle, 'cote': ident.get('cote'), 'libelle': Sleg.get('libelle'),
                                   'nom_original': o['noeud'], 'groupes_originaux': o.get('groupes', []),
                                   'identification': {'methode': ident['methode'], 'detail': ident['detail']},
                                   'couleur_originale': None if o['couleur'] is None else [round(float(x), 4) for x in o['couleur'][:3]],
                                   'visible_par_defaut': o.get('visible', True), 'sommets': int(len(V)), 'faces': int(len(F)),
                                   'segments': int(len(seg)) if seg is not None else 0,
                                   'aire_mm2': round(aire, 5), 'confiance': 'bonne' if ident['methode'] in ('nom', 'correspondance') else 'faible',
                                   'avertissements': o.get('avertissements', [])})
                sorties.append({'nom': nom, 'V': V, 'F': F, 'couleur': couleur[:3], 'alpha': alpha,
                                'segments': seg if seg is not None else []})
                corr_prop[o['noeud']] = {'structure': cle, 'cote': ident.get('cote'), 'methode': ident['methode']}
            man = {'stage': prefixe, 'units': 'mm', 'reference': 'Rana et al. 2014 (Am J Med Genet A 164A:1372) — équipe Moorman, AMC',
                   'source': {'pdf': os.path.basename(fichier), 'flux': k, 'format': flux['format'], 'sha1': flux['sha1'], 'page': flux['page']},
                   'stade_indice': pourquoi, 'mm_par_unite': fmm, 'unites_indice': pourquoi_mm,
                   'centre_retire_mm': [round(float(x), 5) for x in centre],
                   'repere': 'repère du fichier d\'origine (non recalé sur nos embryons), translaté pour centrer la boîte englobante',
                   'structures': structures}
            json.dump(man, open(os.path.join(d, 'manifest.json'), 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
            json.dump(corr_prop, open(os.path.join(d, 'correspondance_proposee.json'), 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
            glb = ecrire_glb(os.path.join(d, f'{prefixe}_arcs_aortiques.glb'), sorties)
            ap = apercu(os.path.join(d, 'apercu.png'), sorties, f'{prefixe} - {os.path.basename(fichier)}') if avec_apercu else False
            info.update({'manifest': os.path.join(dossier_nom, 'manifest.json'), 'glb': glb, 'apercu': ap, 'structures': len(structures)})
            rapport += ['', f'### Flux {k} → `{dossier_nom}/` ({flux["format"]}, page {flux["page"]})', '',
                        f'- stade : **{stade or "?"}** ({pourquoi})',
                        f'- échelle : {fmm:g} mm par unité ({pourquoi_mm})',
                        f'- vues du PDF : {", ".join(info["vues"]) or "aucune"}',
                        '- décodage : ' + resume.replace('\n', ' ; '), '',
                        '| objet | nom d\'origine | structure | côté | méthode | faces | couleur d\'origine |',
                        '|---|---|---|---|---|---|---|']
            for s in structures:
                rapport.append(f'| {s["name"]} | {s["nom_original"]} | {s["structure"] or "?"} | {s["cote"] or ""} | '
                               f'{s["identification"]["methode"]} | {s["faces"]} | {s["couleur_originale"]} |')
            inconnus = [s['nom_original'] for s in structures if s['structure'] is None]
            if inconnus:
                rapport += ['', f'Non identifiés ({len(inconnus)}) : compléter `{dossier_nom}/correspondance_proposee.json` '
                                f'puis relancer avec `--correspondance`.']
    index['rapport'] = 'rapport.md'
    json.dump(index, open(os.path.join(sortie, 'index.json'), 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    open(os.path.join(sortie, 'rapport.md'), 'w', encoding='utf-8').write('\n'.join(rapport) + '\n')
    return index


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('entrees', nargs='*', help='PDF, U3D, PRC ou dossiers (défaut : reference/arcs_aortiques/pdf/)')
    ap.add_argument('--sortie', default=os.path.join(ICI, 'extraits'))
    ap.add_argument('--stade')
    ap.add_argument('--stades', help='JSON {sha1 | fichier#k | fichier : stade}')
    ap.add_argument('--correspondance', help='JSON {nom d\'origine : structure ou {structure, cote}}')
    ap.add_argument('--legende', help='autre fichier de légende')
    ap.add_argument('--mm-par-unite', type=float)
    ap.add_argument('--sans-centrage', action='store_true')
    ap.add_argument('--sans-apercu', action='store_true')
    ap.add_argument('--lister', action='store_true', help='inventaire des flux 3D seulement')
    a = ap.parse_args(argv)
    entrees = a.entrees or [os.path.join(ICI, 'pdf')]
    fichiers = []
    for e in entrees:
        if os.path.isdir(e):
            fichiers += sorted(os.path.join(e, f) for f in os.listdir(e) if f.lower().endswith(('.pdf', '.u3d', '.prc')))
        elif os.path.exists(e):
            fichiers.append(e)
        else:
            print('introuvable :', e)
    if not fichiers:
        print('Aucun fichier. Déposer les PDF 3D du matériel supplémentaire de Rana et al. 2014 dans', os.path.join(ICI, 'pdf'))
        return 1
    L = charger_legende(a.legende)
    corr = json.load(open(a.correspondance, encoding='utf-8')) if a.correspondance else None
    if corr:
        corr = {k: v for k, v in corr.items() if v and (not isinstance(v, dict) or v.get('structure'))}
    table = json.load(open(a.stades, encoding='utf-8')) if a.stades else None
    idx = traiter(fichiers, a.sortie, L, a.stade, corr, table, a.mm_par_unite, not a.sans_centrage, not a.sans_apercu, a.lister)
    for f in idx['fichiers']:
        print(f['fichier'], f'— {len(f["flux"])} flux 3D')
        for fl in f['flux']:
            if 'erreur' in fl:
                print(f'   flux {fl["index"]} : ERREUR {fl["erreur"]}')
            elif 'dossier' in fl:
                print(f'   flux {fl["index"]} : {fl["format"]} -> {fl["dossier"]}/ ({fl["structures"]} structures, stade {fl["stade"]} : {fl["stade_indice"]})')
            elif 'doublon_de' in fl:
                print(f'   flux {fl["index"]} : doublon de {fl["doublon_de"]}')
            else:
                print(f'   flux {fl["index"]} : {fl["format"]} {fl["octets"]} octets, page {fl["page"]}')
        for e in f['erreurs']:
            print('   !', e)
    print('rapport :', os.path.join(a.sortie, 'rapport.md'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
