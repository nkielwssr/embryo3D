"""Lecteur PRC (ISO 14739-1, « Product Representation Compact ») en Python pur, limité à ce qu'il faut pour
récupérer des maillages : conteneur, sections compressées (zlib), globales (couleurs, matériaux, styles),
arbre (occurrences de produits, définitions de pièces, éléments de représentation, transformations) et
tessellations 3D « régulières » (triangles, éventails, bandes).

Non pris en charge (signalé dans `PRCScene.erreurs`) : tessellation « hautement compressée »
(PRC_TYPE_TESS_3D_Compressed), géométrie exacte B-rep / NURBS (seule la tessellation est lue), textures.
Les structures de l'arbre suivent l'écrivain PRC d'Asymptote (writePRC.cc, LGPL) et la spécification ;
si une entité inconnue interrompt la lecture de l'arbre, les maillages restent récupérés, nommés « tess_<i> ».

Usage :
    from prc import lire_prc
    sc = lire_prc(open('modele.prc', 'rb').read())
    for o in sc.objets(): print(o['nom'], o['sommets'].shape, o['couleur'])
"""
import struct
import zlib

import numpy as np

from prc_table import TABLE

# types (PRC.h)
T_ATTRIBUTE = 201
T_CARTESIAN = 202
T_GENERAL = 207
T_REF_PRCBASE = 205
RI_BREP, RI_CURVE, RI_DIRECTION, RI_PLANE, RI_POINTSET, RI_POLYBREP, RI_POLYWIRE, RI_SET, RI_CS = range(232, 241)
ASM_MODELFILE, ASM_FS = 301, 302
ASM_GLOBALS, ASM_TREE, ASM_TESS, ASM_GEOM, ASM_EXTRA = 303, 304, 305, 306, 307
ASM_PO, ASM_PART, ASM_FILTER = 310, 311, 320
TESS_3D, TESS_3D_COMP, TESS_FACE, TESS_WIRE, TESS_MARKUP = 172, 173, 174, 175, 176
G_STYLE, G_MATERIAL, G_PICTURE, G_TEXAPP, G_TEXDEF, G_LINEPATTERN = 701, 702, 703, 711, 712, 721
MKP_VIEW = 501

# drapeaux de faces tessellées
F_POLY, F_TRI, F_FAN, F_STRIP = 0x1, 0x2, 0x4, 0x8
NORMAL_SINGLE = 0x40000000
NORMAL_MASK = 0x3FFFFFFF

TYPES_REFERENCABLES = {203, 204, RI_BREP, RI_CURVE, RI_DIRECTION, RI_PLANE, RI_POINTSET, RI_POLYBREP, RI_POLYWIRE,
                       RI_SET, RI_CS, ASM_PO, ASM_PART, ASM_FILTER, 501, 502, 503, 504, 505, 506, G_STYLE, G_MATERIAL,
                       G_TEXAPP, G_TEXDEF, G_LINEPATTERN, 723, 724, 725, 726, 731, 732, 733, 734, 741, 742}

_DOUBLES = {(nb, code): (est_double, bits) for est_double, nb, code, bits in TABLE}
_BITS8 = [format(i, '08b') for i in range(256)]


class ErreurPRC(Exception):
    pass


class Bits:
    """Flux PRC : bits de poids fort d'abord (le flux est converti une fois en chaîne '0'/'1')."""

    def __init__(self, data):
        self.s = ''.join(map(_BITS8.__getitem__, data))
        self.p = 0
        self.n = len(self.s)
        self.nom_courant = ''
        self.graph_courant = (-1, -1, 1)

    def _fin(self, k):
        if self.p + k > self.n:
            raise ErreurPRC(f'fin de section atteinte (bit {self.p} + {k} > {self.n})')

    def bit(self):
        self._fin(1)
        b = self.s[self.p] == '1'
        self.p += 1
        return b

    def bits(self, k):
        self._fin(k)
        v = int(self.s[self.p:self.p + k], 2) if k else 0
        self.p += k
        return v

    def char(self):
        return self.bits(8)

    def uint(self):
        r, k = 0, 0
        while self.bit():
            r |= self.char() << (8 * k)
            k += 1
        return r

    def int(self):
        r, k = 0, 0
        while self.bit():
            r |= self.char() << (8 * k)
            k += 1
        if k and r & (1 << (8 * k - 1)):
            r -= 1 << (8 * k)
        return r

    def chaine(self):
        if not self.bit():
            return ''
        n = self.uint()
        b = bytes(self.char() for _ in range(n))
        try:
            return b.decode('utf-8')
        except UnicodeDecodeError:
            return b.decode('latin-1')

    def double(self):
        code = 0
        ent = None
        for i in range(1, 23):
            code = (code << 1) | (1 if self.bit() else 0)
            ent = _DOUBLES.get((i, code))
            if ent is not None:
                nb = i
                break
        if ent is None:
            raise ErreurPRC(f'code de double inconnu au bit {self.p}')
        est_double, v = ent
        if nb == 2 and code == 1 and est_double:
            return 0.0
        if self.bit():
            v |= 1 << 63
        else:
            v &= ~(1 << 63)
        if not est_double and self.bit():
            v |= self.bits(4) << 48
            B = list(v.to_bytes(8, 'little'))
            cur = 5
            while cur >= 0:
                if self.bit():
                    B[cur] = self.char()
                    cur -= 1
                    continue
                off = self.bits(3)
                if off == 0:
                    prev = B[cur + 1]
                    for k in range(cur, -1, -1):
                        B[k] = prev
                    break
                if off == 6:
                    prev = B[cur + 1]
                    for k in range(cur, 0, -1):
                        B[k] = prev
                    B[0] = self.char()
                    break
                B[cur] = B[cur + off]
                cur -= 1
            v = int.from_bytes(bytes(B), 'little')
        return struct.unpack('<d', v.to_bytes(8, 'little'))[0]

    def vec3(self):
        return (self.double(), self.double(), self.double())

    # --- éléments communs ---
    def code(self, attendu):
        c = self.uint()
        if c != attendu:
            raise ErreurPRC(f'code {c} lu, {attendu} attendu (bit {self.p})')
        return c

    def attributs(self):
        res = []
        for _ in range(self.uint()):
            self.code(T_ATTRIBUTE)
            titre = self.uint() if self.bit() else self.chaine()
            cles = []
            for _ in range(self.uint()):
                cle = self.uint() if self.bit() else self.chaine()
                t = self.uint()
                val = {1: self.int, 2: self.double, 3: self.uint, 4: self.chaine}.get(t, lambda: None)()
                cles.append((cle, val))
            res.append((titre, cles))
        return res

    def nom(self):
        if not self.bit():
            self.nom_courant = self.chaine()
        return self.nom_courant

    def base(self, referencable):
        attrs = self.attributs()
        nom = self.nom()
        if referencable:
            self.uint(); self.uint(); self.uint()
        return nom, attrs

    def graphiques(self):
        if not self.bit():
            self.graph_courant = (self.uint() - 1, self.uint() - 1, self.char() | (self.char() << 8))
        return self.graph_courant

    def donnees_utilisateur(self):
        n = self.uint()
        self._fin(n)
        self.p += n

    def uuid(self):
        return tuple(self.uint() for _ in range(4))

    def reinit(self):
        self.nom_courant = ''
        self.graph_courant = (-1, -1, 1)


# --- transformations ------------------------------------------------------------------------------------------
def _transformation3d(b):
    t = b.uint()
    if t == T_GENERAL:
        vals = [b.double() for _ in range(16)]
        return np.array(vals).reshape(4, 4).T                   # stockage colonne par colonne
    if t != T_CARTESIAN:
        raise ErreurPRC(f'transformation de type {t} inconnue')
    comp = b.char()
    M = np.eye(4)
    O = np.array(b.vec3()) if comp & 0x01 else np.zeros(3)
    if comp & 0x20:
        X, Y, Z = np.array(b.vec3()), np.array(b.vec3()), np.array(b.vec3())
    elif comp & 0x02:
        X, Y = np.array(b.vec3()), np.array(b.vec3())
        Z = np.cross(X, Y)
        if comp & 0x04:
            Z = -Z
    else:
        X, Y, Z = np.eye(3)
        if comp & 0x04:
            Z = -Z
    if comp & 0x10:
        s = np.array(b.vec3())
    elif comp & 0x08:
        s = np.full(3, b.double())
    else:
        s = np.ones(3)
    if comp & 0x40:
        [b.double() for _ in range(4)]
    M[:3, 0], M[:3, 1], M[:3, 2], M[:3, 3] = X * s[0], Y * s[1], Z * s[2], O
    return M


def _balises_vides(b, quoi):
    for _ in range(4):
        if b.uint():
            raise ErreurPRC(f'{quoi} : annotations (markups) non prises en charge')


# --- sections ---------------------------------------------------------------------------------------------------
def _schema(b):
    n = b.uint()
    for _ in range(n):
        b.uint()                                    # type d'entité
        for _ in range(b.uint()):                   # jetons du schéma
            b.uint()
    return n


def _globales(b, sc):
    _schema(b)
    b.code(ASM_GLOBALS)
    b.base(False)
    for _ in range(b.uint()):
        b.uuid()
    b.double(); b.double()
    b.chaine()
    for _ in range(b.uint()):                       # polices
        b.chaine(); b.uint()
        for _ in range(b.uint()):
            b.uint(); b.char()
    sc.couleurs.extend(b.vec3() for _ in range(b.uint()))
    for _ in range(b.uint()):                       # images
        b.code(G_PICTURE)
        b.base(False)
        b.int(); b.uint(); b.uint(); b.uint()
    if b.uint():
        raise ErreurPRC('définitions de textures non prises en charge')
    for _ in range(b.uint()):
        t = b.uint()
        b.base(True)
        if t == G_MATERIAL:
            a, d, e, s = (b.uint() - 1 for _ in range(4))
            sh = b.double()
            alphas = [b.double() for _ in range(4)]
            sc.materiaux.append({'ambiant': a, 'diffus': d, 'emissif': e, 'speculaire': s, 'brillance': sh,
                                 'alpha_diffus': alphas[1]})
        elif t == G_TEXAPP:
            m = b.uint() - 1
            b.uint(); b.uint(); b.uint()
            sc.materiaux.append({'texture_de': m})
        else:
            raise ErreurPRC(f'matériau de type {t} inconnu')
    for _ in range(b.uint()):                       # motifs de trait
        b.code(G_LINEPATTERN)
        b.base(True)
        for _ in range(b.uint()):
            b.double()
        b.double()
        b.bit()
    for _ in range(b.uint()):                       # styles
        b.code(G_STYLE)
        b.base(True)
        largeur = b.double()
        b.bit()
        b.uint()
        est_mat = b.bit()
        idx = b.uint() - 1
        transp = b.char() if b.bit() else None
        for _ in range(3):
            if b.bit():
                b.char()
        sc.styles.append({'materiau': est_mat, 'index': idx, 'transparence': transp, 'largeur': largeur})
    if b.uint():
        raise ErreurPRC('motifs de remplissage non pris en charge')
    for _ in range(b.uint()):                       # systèmes de coordonnées de référence
        b.code(RI_CS)
        ri = _contenu_ri(b)
        ri['matrice'] = _transformation3d(b)
        b.donnees_utilisateur()
        sc.reperes.append(ri['matrice'])
    b.donnees_utilisateur()


def _contenu_ri(b):
    nom, attrs = b.base(True)
    g = b.graphiques()
    return {'nom': nom, 'style': g[1], 'repere': b.uint() - 1, 'tess': b.uint() - 1}


def _ri(b):
    t = b.uint()
    ri = _contenu_ri(b)
    ri['type'] = t
    if t == RI_BREP:
        if b.bit():
            b.uint(); b.uint()
        b.bit()
    elif t == RI_CURVE:
        if b.bit():
            b.uint(); b.uint()
    elif t == RI_POLYBREP:
        b.bit()
    elif t == RI_POLYWIRE:
        pass
    elif t == RI_POINTSET:
        ri['points'] = [b.vec3() for _ in range(b.uint())]
    elif t == RI_DIRECTION:
        if b.bit():
            b.vec3()
        b.vec3()
    elif t == RI_SET:
        ri['elements'] = [_ri(b) for _ in range(b.uint())]
    elif t == RI_CS:
        ri['matrice'] = _transformation3d(b)
    else:
        raise ErreurPRC(f'élément de représentation de type {t} non pris en charge')
    b.donnees_utilisateur()
    return ri


def _vue(b):
    b.code(MKP_VIEW)
    b.base(True)
    b.graphiques()
    if b.uint():
        raise ErreurPRC('vues annotées non prises en charge')
    raise ErreurPRC('vues de pièce/produit non prises en charge')


def _arbre(b, sc):
    b.code(ASM_TREE)
    b.base(False)
    for _ in range(b.uint()):
        b.code(ASM_PART)
        nom, _ = b.base(True)
        g = b.graphiques()
        b.vec3(); b.vec3()
        part = {'nom': nom, 'style': g[1], 'ris': []}
        sc.pieces.append(part)
        for _ in range(b.uint()):
            part['ris'].append(_ri(b))
        _balises_vides(b, f'pièce {nom!r}')
        for _ in range(b.uint()):
            _vue(b)
        b.donnees_utilisateur()
    for _ in range(b.uint()):
        b.code(ASM_PO)
        nom, _ = b.base(True)
        g = b.graphiques()
        po = {'nom': nom, 'style': g[1], 'piece': b.uint() - 1, 'prototype': -1, 'fils': [], 'loc': np.eye(4)}
        sc.occurrences.append(po)
        po['prototype'] = b.uint() - 1
        if po['prototype'] >= 0 and not b.bit():
            b.uuid()
        ext = b.uint() - 1
        if ext >= 0 and not b.bit():
            b.uuid()
        po['fils'] = [b.uint() for _ in range(b.uint())]
        b.char()
        b.bit(); po['unite'] = b.double()
        b.char()
        b.int()
        if b.bit():
            po['loc'] = _transformation3d(b)
        for _ in range(b.uint()):                   # références
            b.code(T_REF_PRCBASE)
            b.uint()
            if not b.bit():
                b.uuid()
            b.uint()
        _balises_vides(b, f'occurrence {nom!r}')
        for _ in range(b.uint()):
            _vue(b)
        if b.bit():
            raise ErreurPRC('filtre d\'entités non pris en charge')
        if b.uint():
            raise ErreurPRC('filtres d\'affichage non pris en charge')
        if b.uint():
            raise ErreurPRC('paramètres d\'affichage de scène non pris en charge')
        b.donnees_utilisateur()
    b.code(ASM_FS)
    b.base(False)
    b.uint()
    sc.racine = b.uint()


def _base_tess(b):
    b.bit()
    n = b.uint()
    return np.array([b.double() for _ in range(n)], dtype=np.float64)


def _tess3d(b):
    coords = _base_tess(b)
    b.bit(); b.bit()
    recalc = b.bit()
    if recalc:
        b.char(); b.double()
    normales = [b.double() for _ in range(b.uint())]
    [b.uint() for _ in range(b.uint())]             # indices d'arêtes
    tri = np.array([b.uint() for _ in range(b.uint())], dtype=np.int64)
    avec_normales = len(normales) > 0 and not recalc
    triangles = []
    tri_styles = []
    couleurs_sommets = []
    for _ in range(b.uint()):
        b.code(TESS_FACE)
        nla = b.uint()
        styles = [b.uint() - 1 for _ in range(nla)]
        b.uint()
        [b.uint() for _ in range(b.uint())]
        drapeaux = b.uint()
        debut = b.uint()
        tailles = [b.uint() for _ in range(b.uint())]
        ntex = b.uint()
        t, nsom = _triangles_face(tri, drapeaux, debut, tailles, ntex, avec_normales)
        if b.bit():                                  # couleurs par sommet (une par sommet de triangle / éventail / bande)
            nbv = 4 if b.bit() else 3
            if b.bit():
                raise ErreurPRC('couleurs de sommets « optimisées » non prises en charge')
            cs = [tuple(b.char() for _ in range(nbv))]
            for _ in range(nsom - 1):
                cs.append(cs[-1] if b.bit() else tuple(b.char() for _ in range(nbv)))
            couleurs_sommets.extend(cs)
        if nla:
            b.uint()
        triangles.extend(t)
        tri_styles.extend((styles[i] if i < len(styles) else styles[-1]) if styles else -1 for i in range(len(t)))
    [b.double() for _ in range(b.uint())]           # coordonnées de texture
    if len(coords) % 3:
        raise ErreurPRC('nombre de coordonnées non multiple de 3')
    V = coords.reshape(-1, 3)
    F = np.array(triangles, dtype=np.int64).reshape(-1, 3) // 3
    return {'sommets': V, 'faces': F, 'styles_faces': tri_styles, 'couleurs_sommets_triangles': couleurs_sommets}


def _triangles_face(idx, drapeaux, p, tailles, ntex, avec_normales):
    """Triangles (indices de coordonnées, multiples de 3) d'une face tessellée et nombre de sommets parcourus.
    Pour chaque drapeau présent, dans l'ordre croissant : triangles (1 taille = nombre de triangles), éventails et
    bandes (1 taille = nombre d'éventails, puis une taille par éventail, bit NORMAL_Single possible).
    Par sommet : [indice de normale] [indices de texture] indice de point."""
    k = 0
    out = []
    nsom = 0

    def sommet(avec_n, nt):
        nonlocal p, nsom
        if avec_n:
            p += 1
        p += nt
        v = int(idx[p])
        p += 1
        nsom += 1
        return v

    for bitv in range(16):
        if not drapeaux & (1 << bitv):
            continue
        genre = 1 << (bitv % 4)
        une_normale = (bitv // 4) in (1, 3)
        nt = max(ntex, 1) if bitv >= 8 else 0
        if genre == F_POLY:
            raise ErreurPRC('faces polygonales (polyface) non prises en charge')
        if genre == F_TRI:
            n = tailles[k]
            k += 1
            for _ in range(n):
                if une_normale and avec_normales:
                    p += 1
                an = avec_normales and not une_normale
                out.append((sommet(an, nt), sommet(an, nt), sommet(an, nt)))
            continue
        nb = tailles[k]
        k += 1
        for _ in range(nb):
            brut = tailles[k]
            k += 1
            nv = brut & NORMAL_MASK
            seule = une_normale or bool(brut & NORMAL_SINGLE)
            if seule and avec_normales:
                p += 1
            vs = [sommet(avec_normales and not seule, nt) for _ in range(nv)]
            for i in range(nv - 2):
                if genre == F_FAN:
                    out.append((vs[0], vs[i + 1], vs[i + 2]))
                elif i % 2 == 0:
                    out.append((vs[i], vs[i + 1], vs[i + 2]))
                else:
                    out.append((vs[i + 1], vs[i], vs[i + 2]))
    return out, nsom


def _tessellations(b, sc):
    b.code(ASM_TESS)
    b.base(False)
    for i in range(b.uint()):
        t = b.uint()
        if t == TESS_3D:
            sc.tess.append(_tess3d(b))
        elif t == TESS_WIRE:
            _base_tess(b)
            [b.uint() for _ in range(b.uint())]
            if b.bit():
                raise ErreurPRC('couleurs de fils non prises en charge')
            sc.tess.append(None)
        elif t == TESS_MARKUP:
            _base_tess(b)
            [b.uint() for _ in range(b.uint())]
            [b.chaine() for _ in range(b.uint())]
            b.chaine(); b.char()
            sc.tess.append(None)
        elif t == TESS_3D_COMP:
            raise ErreurPRC(f'tessellation {i} « hautement compressée » (PRC_TYPE_TESS_3D_Compressed) non prise en charge')
        else:
            raise ErreurPRC(f'tessellation {i} de type {t} inconnu')


# --- scène ------------------------------------------------------------------------------------------------------
class PRCScene:
    def __init__(self):
        self.couleurs, self.materiaux, self.styles, self.reperes = [], [], [], []
        self.pieces, self.occurrences, self.tess = [], [], []
        self.racine = None
        self.erreurs = []
        self.version = None

    def couleur_style(self, s):
        if s is None or s < 0 or s >= len(self.styles):
            return None, 1.0
        st = self.styles[s]
        alpha = 1.0 if st['transparence'] is None else st['transparence'] / 255.0
        if st['materiau']:
            m = self.materiaux[st['index']] if 0 <= st['index'] < len(self.materiaux) else {}
            ci = m.get('diffus', -1)
            if m.get('alpha_diffus') is not None and st['transparence'] is None:
                alpha = m['alpha_diffus']
        else:
            ci = st['index']
        if ci is None or ci < 0:
            return None, alpha
        ci = ci // 3 if ci // 3 < len(self.couleurs) and ci % 3 == 0 else ci
        if 0 <= ci < len(self.couleurs):
            return tuple(self.couleurs[ci]), alpha
        return None, alpha

    def _ris_plats(self, ris, M, chemin, style_herite):
        for ri in ris:
            st = ri['style'] if ri['style'] >= 0 else style_herite
            Mi = M
            if 0 <= ri['repere'] < len(self.reperes):
                Mi = M @ self.reperes[ri['repere']]
            if ri['type'] == RI_SET:
                yield from self._ris_plats(ri.get('elements', []), Mi, chemin + [ri['nom']], st)
            elif ri['tess'] >= 0:
                yield ri, Mi, chemin, st

    def objets(self):
        out = []
        utilises = set()
        if self.occurrences:
            fils = {f for po in self.occurrences for f in po['fils']}
            racines = [i for i in range(len(self.occurrences)) if i not in fils]
            if self.racine is not None and self.racine < len(self.occurrences):
                racines = [self.racine] + [r for r in racines if r != self.racine]
            vus = set()

            def parcourir(i, M, chemin, style, prof=0):
                if i in vus or prof > 64 or i >= len(self.occurrences):
                    return
                vus.add(i)
                po = self.occurrences[i]
                Mi = M @ po['loc']
                st = po['style'] if po['style'] >= 0 else style
                piece, fils_ = po['piece'], po['fils']
                proto = po['prototype']
                while piece < 0 and not fils_ and 0 <= proto < len(self.occurrences) and proto != i:
                    piece, fils_ = self.occurrences[proto]['piece'], self.occurrences[proto]['fils']
                    proto = self.occurrences[proto]['prototype']
                ch = chemin + ([po['nom']] if po['nom'] else [])
                if 0 <= piece < len(self.pieces):
                    pc = self.pieces[piece]
                    stp = pc['style'] if pc['style'] >= 0 else st
                    for ri, Mr, chr_, st_ri in self._ris_plats(pc['ris'], Mi, ch, stp):
                        self._objet(out, ri, Mr, chr_, st_ri)
                        utilises.add(ri['tess'])
                for f in fils_:
                    parcourir(f, Mi, ch, st, prof + 1)
                vus.discard(i)

            for r in racines:
                parcourir(r, np.eye(4), [], -1)
        # tessellations non rattachées à l'arbre (arbre illisible) : on les garde
        for i, t in enumerate(self.tess):
            if t is not None and i not in utilises and len(t['faces']):
                self._objet(out, {'nom': f'tess_{i}', 'tess': i}, np.eye(4), [], -1)
        return out

    def _objet(self, out, ri, M, chemin, style):
        i = ri['tess']
        t = self.tess[i] if 0 <= i < len(self.tess) else None
        if t is None or not len(t['faces']):
            return
        V = t['sommets'] @ M[:3, :3].T + M[:3, 3]
        F = t['faces']
        ok = (F < len(V)).all(1)
        F = F[ok]
        couleur, alpha = self.couleur_style(style)
        if couleur is None:
            sts = [s for s in t['styles_faces'] if s >= 0]
            if sts:
                couleur, alpha = self.couleur_style(max(set(sts), key=sts.count))
        nom = ri['nom'] or (chemin[-1] if chemin else f'tess_{i}')
        out.append({'nom': nom, 'noeud': nom, 'groupes': [c for c in chemin if c and c != nom], 'sommets': V, 'faces': F,
                    'couleurs_sommets': None, 'materiau': None, 'couleur': couleur, 'opacite': alpha, 'visible': True,
                    'avertissements': [], 'tessellation': i})

    def resume(self):
        l = [f'PRC version {self.version} : {len(self.occurrences)} occurrences, {len(self.pieces)} pièces, '
             f'{sum(1 for t in self.tess if t is not None)} tessellations 3D, {len(self.couleurs)} couleurs, {len(self.styles)} styles']
        for e in self.erreurs:
            l.append('ERREUR ' + e)
        return '\n'.join(l)


def _decompresser(buf, off):
    d = zlib.decompressobj()
    return d.decompress(buf[off:]) + d.flush()


def lire_prc(buf):
    buf = bytes(buf)
    if buf[:3] != b'PRC':
        raise ErreurPRC('signature PRC absente')
    sc = PRCScene()
    p = 3
    sc.version, _ = struct.unpack_from('<II', buf, p)
    p += 8 + 32
    nfs, = struct.unpack_from('<I', buf, p)
    p += 4
    structures = []
    for _ in range(nfs):
        p += 16 + 4
        noff, = struct.unpack_from('<I', buf, p)
        p += 4
        structures.append(struct.unpack_from(f'<{noff}I', buf, p))
        p += 4 * noff
    for offs in structures:
        sections = []
        for o in offs[1:]:
            try:
                sections.append(_decompresser(buf, o))
            except zlib.error as e:
                sections.append(None)
                sc.erreurs.append(f'section à l\'octet {o} : zlib {e}')
        base_tess = len(sc.tess)
        base_part = len(sc.pieces)
        base_po = len(sc.occurrences)
        nb_styles, nb_mat, nb_coul, nb_rep = len(sc.styles), len(sc.materiaux), len(sc.couleurs), len(sc.reperes)
        for fonction, idx, quoi in ((_globales, 0, 'globales'), (_tessellations, 2, 'tessellations'), (_arbre, 1, 'arbre')):
            if idx >= len(sections) or sections[idx] is None:
                continue
            b = Bits(sections[idx])
            try:
                fonction(b, sc)
            except (ErreurPRC, IndexError, ValueError, struct.error, OverflowError) as e:
                sc.erreurs.append(f'{quoi} : {e}')
        if nfs > 1:
            # plusieurs structures : leurs index sont locaux -> décalage vers les listes communes de la scène
            def decaler_ris(ris):
                for ri in ris:
                    ri['tess'] += base_tess if ri['tess'] >= 0 else 0
                    ri['style'] += nb_styles if ri['style'] >= 0 else 0
                    ri['repere'] += nb_rep if ri['repere'] >= 0 else 0
                    decaler_ris(ri.get('elements', []))
            for st in sc.styles[nb_styles:]:
                if st['index'] >= 0:
                    st['index'] += nb_mat if st['materiau'] else 3 * nb_coul
            for m in sc.materiaux[nb_mat:]:
                for k in ('ambiant', 'diffus', 'emissif', 'speculaire'):
                    if m.get(k, -1) >= 0:
                        m[k] += 3 * nb_coul
            for pc in sc.pieces[base_part:]:
                pc['style'] += nb_styles if pc['style'] >= 0 else 0
                decaler_ris(pc['ris'])
            for po in sc.occurrences[base_po:]:
                po['piece'] += base_part if po['piece'] >= 0 else 0
                po['prototype'] += base_po if po['prototype'] >= 0 else 0
                po['style'] += nb_styles if po['style'] >= 0 else 0
                po['fils'] = [f + base_po for f in po['fils']]
    return sc


if __name__ == '__main__':
    import sys
    for chemin in sys.argv[1:]:
        s = lire_prc(open(chemin, 'rb').read())
        print(chemin)
        print(s.resume())
        for o in s.objets():
            print(f"  {o['nom']:30s} {len(o['sommets']):7d} sommets {len(o['faces']):7d} faces  couleur={o['couleur']} groupes={o['groupes']}")
