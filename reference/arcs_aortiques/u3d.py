"""Décodeur U3D (ECMA-363) en Python pur, sans dépendance compilée.

Port du décodeur de référence de la bibliothèque U3D d'Intel (Apache 2.0, livrée dans les sources de MeshLab :
src/external/u3d/src/RTL/Component/{BitStream,Importing,CLODAuthor}). Les passages délicats reproduisent le code
C++ ligne à ligne, car le moindre écart (ordre d'un ensemble, contexte d'un symbole) désynchronise le flux :
  - CIFXBitStreamX         : décodeur arithmétique 16 bits (contextes statiques / dynamiques, octets « rapides »)
  - IFXHistogramDynamic    : histogrammes adaptatifs (symbole d'échappement 0, division par 2 au-delà de 0x1FFF)
  - CIFXSetX               : ensembles triés par ordre DÉCROISSANT (l'index local d'un sommet dépend de cet ordre)
  - CIFXAuthorCLODDecoder  : maillage CLOD (déclaration, maillage de base, raffinement progressif par vertex split)

Blocs lus : en-tête, chaînes de modificateurs, nœuds groupe / modèle, maillages CLOD, modificateurs d'ombrage,
shaders « lit texture », matériaux, métadonnées. Les autres blocs (vues, lumières, textures, animations,
jeux de points / lignes) sont sautés grâce à leur taille et signalés dans `U3DScene.ignores`.

Usage :
    from u3d import lire_u3d
    scene = lire_u3d(open('modele.u3d', 'rb').read())
    for obj in scene.objets():        # un objet par (nœud modèle, parent), coordonnées dans le repère monde
        print(obj['nom'], obj['sommets'].shape, obj['faces'].shape, obj['couleur'])
"""
import struct
from bisect import bisect_left

import numpy as np

# --- constantes (IFXACContext.h, IFXBlockTypes.h) -------------------------------------------------------------
STATIC_FULL = 0x400
MAX_RANGE = STATIC_FULL + 0x3FFF
ELEPHANT = 0x1FFF
MAX_SYMBOL_IN_HISTOGRAM = 0xFFFF

B_HEADER = 0x00443355
B_FILE_REFERENCE = 0xFFFFFF12
B_MODIFIER_CHAIN = 0xFFFFFF14
B_PRIORITY = 0xFFFFFF15
B_NEW_OBJECT_TYPE = 0xFFFFFF16
B_GROUP = 0xFFFFFF21
B_MODEL = 0xFFFFFF22
B_LIGHT = 0xFFFFFF23
B_VIEW = 0xFFFFFF24
B_CLOD_DECL = 0xFFFFFF31
B_CLOD_BASE = 0xFFFFFF3B
B_CLOD_PROG = 0xFFFFFF3C
B_POINTSET = 0xFFFFFF36
B_POINTSET_CONT = 0xFFFFFF3E
B_LINESET = 0xFFFFFF37
B_LINESET_CONT = 0xFFFFFF3F
B_SHADING = 0xFFFFFF45
B_LIGHT_RES = 0xFFFFFF51
B_VIEW_RES = 0xFFFFFF52
B_SHADER = 0xFFFFFF53
B_MATERIAL = 0xFFFFFF54
B_TEXTURE = 0xFFFFFF55
B_TEXTURE_CONT = 0xFFFFFF5C
B_MOTION = 0xFFFFFF56

NOMS_BLOCS = {B_HEADER: 'en-tête', B_FILE_REFERENCE: 'référence de fichier', B_MODIFIER_CHAIN: 'chaîne de modificateurs',
              B_PRIORITY: 'priorité', B_NEW_OBJECT_TYPE: 'type d\'objet étendu', B_GROUP: 'nœud groupe',
              B_MODEL: 'nœud modèle', B_LIGHT: 'nœud lumière', B_VIEW: 'nœud vue', B_CLOD_DECL: 'maillage CLOD',
              B_CLOD_BASE: 'maillage de base', B_CLOD_PROG: 'maillage progressif', B_POINTSET: 'jeu de points',
              B_POINTSET_CONT: 'jeu de points (suite)', B_LINESET: 'jeu de lignes', B_LINESET_CONT: 'jeu de lignes (suite)',
              B_SHADING: 'modificateur d\'ombrage', B_LIGHT_RES: 'lumière', B_VIEW_RES: 'vue', B_SHADER: 'shader',
              B_MATERIAL: 'matériau', B_TEXTURE: 'texture', B_TEXTURE_CONT: 'texture (suite)', B_MOTION: 'animation',
              0xFFFFFF41: 'glyphe 2D', 0xFFFFFF42: 'subdivision', 0xFFFFFF43: 'animation (modificateur)',
              0xFFFFFF44: 'poids d\'os', 0xFFFFFF46: 'modificateur CLOD'}

PROFILE_NOCOMPRESSION = 0x4
PROFILE_UNITSSCALE = 0x8

# contextes du maillage progressif
C_BASE_SHADING = 1
C_NUM_NEW_FACES = 1
C_SHADING = 65
C_ORIENTATION = 2
C_THIRD_TYPE = 3
C_LOCAL_3RD = 4
C_STAY_MOVE = 15
C_POS_SIGNS, C_POS_X, C_POS_Y, C_POS_Z = 20, 21, 22, 23
C_TEX_SIGNS_NEW = 103
C_TEX_MAG = (33, 34, 35, 36)
C_TEX_DUP, C_TEX_SPLIT_TYPE, C_TEX_LOCAL, C_TEX_GLOBAL = 39, 29, 121, 122
C_NUM_LOCAL_NORMALS, C_NORMAL_SIGNS, C_NORMAL_X, C_NORMAL_Y, C_NORMAL_Z, C_NORMAL_LOCAL = 40, 41, 42, 43, 44, 45
C_COLOR_DUP, C_COLOR_SPLIT_TYPE = 56, 55
C_COLOR_MAG = (60, 61, 62, 63)
C_NEW_DIFFUSE_COUNT, C_DIFFUSE_SIGN, C_NEW_SPECULAR_COUNT, C_SPECULAR_SIGN = 99, 100, 101, 102
C_DIFF_KEEP, C_DIFF_TYPE, C_DIFF_NEW, C_DIFF_LOCAL, C_DIFF_GLOBAL = 104, 105, 106, 107, 108
C_SPEC_KEEP, C_SPEC_TYPE, C_SPEC_NEW, C_SPEC_LOCAL, C_SPEC_GLOBAL = 109, 110, 111, 112, 113
C_TEX_KEEP, C_TEX_TYPE, C_TEX_NEW, C_TEX_LOCALI, C_TEX_GLOBALI = 114, 115, 116, 117, 118
C_COLOR_LOCAL, C_COLOR_GLOBAL = 119, 120
C_NEW_TEX_COUNT = 123

# inversion de l'ordre des bits sur 8 et 15 bits
_SWAP8 = bytes(int('{:08b}'.format(i)[::-1], 2) for i in range(256))
_REV15 = [int('{:015b}'.format(i)[::-1], 2) for i in range(1 << 15)]


class ErreurU3D(Exception):
    pass


# --- histogramme dynamique (arbre de Fenwick pour les fréquences cumulées) ----------------------------------------
class Histogramme:
    __slots__ = ('n', 'cnt', 'tree', 'total')

    def __init__(self):
        self.n = 128
        self.cnt = [0] * self.n
        self.tree = [0] * (self.n + 1)
        self.total = 0
        self._inc(0, 1)                       # symbole d'échappement, compte initial 1

    def _inc(self, s, d):
        self.cnt[s] += d
        self.total += d
        i = s + 1
        t, n = self.tree, self.n
        while i <= n:
            t[i] += d
            i += i & -i

    def _reconstruire(self):
        n = self.n
        t = [0] * (n + 1)
        for s, c in enumerate(self.cnt):
            if c:
                i = s + 1
                while i <= n:
                    t[i] += c
                    i += i & -i
        self.tree = t
        self.total = sum(self.cnt)

    def cumul(self, s):                        # somme des comptes des symboles < s
        if s >= self.n:
            return self.total
        r, i, t = 0, s, self.tree
        while i > 0:
            r += t[i]
            i -= i & -i
        return r

    def freq(self, s):
        return self.cnt[s] if s < self.n else 0

    def symbole(self, f):                      # symbole s tel que cumul(s) <= f < cumul(s) + freq(s)
        if f >= self.total:
            return 0
        pos, rest, t, n = 0, f, self.tree, self.n
        step = 1 << (n.bit_length() - 1)
        while step:
            nxt = pos + step
            if nxt <= n and t[nxt] <= rest:
                pos = nxt
                rest -= t[nxt]
            step >>= 1
        return pos

    def ajouter(self, s):
        if s > MAX_SYMBOL_IN_HISTOGRAM:
            return
        if self.total >= ELEPHANT:            # division par deux ; l'échappement garde au moins 1
            self.cnt = [c >> 1 for c in self.cnt]
            self.cnt[0] += 1
            self._reconstruire()
        if s >= self.n:
            nn = self.n
            while nn <= s:
                nn *= 2
            self.cnt.extend([0] * (nn - self.n))
            self.n = nn
            self._reconstruire()
        self._inc(s, 1)


# --- flux binaire compressé (CIFXBitStreamX, lecture seule) ------------------------------------------------------
class FluxBits:
    def __init__(self, data, sans_compression=False):
        self.d = bytes(data) + b'\x00' * 16      # lire au-delà de la fin renvoie des zéros (comme la référence)
        self.nbits = len(data) * 8
        self.pos = 0
        self.low = 0
        self.high = 0xFFFF
        self.uf = 0
        self.hist = {}
        self.sans_compression = sans_compression

    # bits : le flux est lu bit de poids faible d'abord, octet par octet (= mots U32 petit-boutistes)
    def _bit(self, p):
        return (self.d[p >> 3] >> (p & 7)) & 1

    def _code(self):
        p = self.pos
        d = self.d
        first = (d[p >> 3] >> (p & 7)) & 1
        q = p + 1 + self.uf
        b = q >> 3
        v = (d[b] | (d[b + 1] << 8) | (d[b + 2] << 16)) >> (q & 7)
        return (first << 15) | _REV15[v & 0x7FFF]

    def _renorm(self, low, high):
        n = 0
        while (low & 0x8000) == (high & 0x8000):
            low = (low & 0x7FFF) << 1
            high = ((high & 0x7FFF) << 1) | 1
            n += 1
        if n:
            n += self.uf
            self.uf = 0
        u = 0
        while (low & 0x4000) and not (high & 0x4000):
            low = (low & 0x3FFF) << 1
            high = ((high & 0x3FFF) << 1) | 1
            u += 1
        self.uf += u
        self.low = low & 0x7FFF
        self.high = high | 0x8000
        self.pos += n

    def _statique(self, ctx):
        total = ctx - STATIC_FULL
        code = self._code()
        low = self.low
        rng = self.high + 1 - low
        cum = (total * (1 + code - low) - 1) // rng
        high = low - 1 + rng * (cum + 1) // total
        low = low + rng * cum // total
        self._renorm(low, high)
        return cum + 1

    def _dynamique(self, ctx):
        h = self.hist.get(ctx)
        if h is None:
            h = self.hist[ctx] = Histogramme()
        code = self._code()
        total = h.total
        low = self.low
        rng = self.high + 1 - low
        f = (total * (1 + code - low) - 1) // rng
        s = h.symbole(f)
        cum = h.cumul(s)
        fr = h.freq(s)
        high = low - 1 + rng * (cum + fr) // total
        low = low + rng * cum // total
        h.ajouter(s)
        self._renorm(low, high)
        return s

    def _symbole(self, ctx):
        if ctx == 0:
            return self._statique(STATIC_FULL + 256)
        if ctx > STATIC_FULL:
            return self._statique(ctx)
        return self._dynamique(ctx)

    # lectures non compressées
    def u8(self):
        if self.high == 0xFFFF and self.low == 0 and self.uf == 0:
            p = self.pos
            b = p >> 3
            v = ((self.d[b] | (self.d[b + 1] << 8)) >> (p & 7)) & 0xFF
            self.pos = p + 8
            return v
        return _SWAP8[self._statique(STATIC_FULL + 256) - 1]

    def u16(self):
        lo = self.u8()
        return lo | (self.u8() << 8)

    def u32(self):
        lo = self.u16()
        return lo | (self.u16() << 16)

    def u64(self):
        lo = self.u32()
        return lo | (self.u32() << 32)

    def i32(self):
        return struct.unpack('<i', struct.pack('<I', self.u32()))[0]

    def f32(self):
        return struct.unpack('<f', struct.pack('<I', self.u32()))[0]

    def f64(self):
        return struct.unpack('<d', struct.pack('<Q', self.u64()))[0]

    def chaine(self):
        n = self.u16()
        b = bytes(self.u8() for _ in range(n))
        try:
            return b.decode('utf-8')
        except UnicodeDecodeError:
            return b.decode('latin-1')

    # lectures compressées
    def _compresse(self, ctx, brut):
        if self.sans_compression:
            return brut()
        if ctx and ctx < MAX_RANGE:
            s = self._symbole(ctx)
            if s:
                return s - 1
            v = brut()
            if ctx <= STATIC_FULL:
                h = self.hist.get(ctx)
                if h is None:
                    h = self.hist[ctx] = Histogramme()
                h.ajouter(v + 1)
            return v
        return brut()

    def cu32(self, ctx):
        return self._compresse(ctx, self.u32)

    def cu16(self, ctx):
        return self._compresse(ctx, self.u16)

    def cu8(self, ctx):
        return self._compresse(ctx, self.u8)

    def depasse(self):
        return self.pos > self.nbits


# --- ensembles triés par ordre décroissant (CIFXSetX) ------------------------------------------------------------
class Ens:
    """Ensemble d'entiers ; membre(j) = j-ième plus grand, comme CIFXSetX (tri décroissant)."""
    __slots__ = ('neg',)

    def __init__(self, it=()):
        self.neg = sorted({-v for v in it})

    def add(self, v):
        i = bisect_left(self.neg, -v)
        if i == len(self.neg) or self.neg[i] != -v:
            self.neg.insert(i, -v)

    def remove(self, v):
        i = bisect_left(self.neg, -v)
        if i < len(self.neg) and self.neg[i] == -v:
            del self.neg[i]

    def __contains__(self, v):
        i = bisect_left(self.neg, -v)
        return i < len(self.neg) and self.neg[i] == -v

    def __len__(self):
        return len(self.neg)

    def membre(self, j):
        if j >= len(self.neg):
            raise ErreurU3D(f'index local {j} hors de l\'ensemble (taille {len(self.neg)})')
        return -self.neg[j]

    def __iter__(self):
        return (-v for v in self.neg)


# --- maillage CLOD ------------------------------------------------------------------------------------------------
class MaillageCLOD:
    def __init__(self, nom):
        self.nom = nom
        self.pret = False
        self.avertissements = []

    # 0xFFFFFF31
    def declaration(self, b):
        b.chaine()
        b.u32()                                   # index de chaîne
        attrs = b.u32()
        self.sans_normales = bool(attrs & 1)
        (self.max_faces, self.max_positions, self.max_normales, self.max_diffus, self.max_speculaires,
         self.max_tex, nmat) = (b.u32() for _ in range(7))
        self.materiaux = []
        for _ in range(nmat):
            a = b.u32()
            ncouches = b.u32()
            dims = [b.u32() for _ in range(ncouches)]
            orig = b.u32()
            self.materiaux.append({'diffus': bool(a & 1), 'speculaire': bool(a & 2), 'couches': ncouches,
                                   'dims': dims, 'id_origine': orig})
        self.res_min = b.u32()
        self.res_max_finale = b.u32()
        self.qualite = [b.u32() for _ in range(3)]
        self.iq_pos, self.iq_norm, self.iq_tex, self.iq_diff, self.iq_spec = (b.f32() for _ in range(5))
        self.param_pli, self.param_maj, self.tol_normale = (b.f32() for _ in range(3))
        self.nb_os = b.u32()
        if self.nb_os:
            self.avertissements.append(f'{self.nb_os} os (squelette) ignorés')
        # état reconstruit
        n = max(self.res_max_finale, self.max_positions)
        self.positions = [None] * n
        self.f_pos = []                           # coins [a, b, c] (listes mutables)
        self.f_mat = []
        self.f_diff = []
        self.f_spec = []
        self.f_tex = [[] for _ in range(8)]
        self.diffus = []
        self.speculaires = []
        self.texcoords = []
        self.nb_faces = 0
        self.nb_diffus = 0
        self.nb_spec = 0
        self.nb_tex = 0
        self.res = 0
        self.adj = {}                             # position -> Ens des faces qui l'utilisent
        self.pret = True

    def _mat(self, m):
        if 0 <= m < len(self.materiaux):
            return self.materiaux[m]
        return {'diffus': False, 'speculaire': False, 'couches': 0}

    def _set_face(self, liste, i, val):
        while len(liste) <= i:
            liste.append(None)
        liste[i] = val

    def _faces_de(self, p):
        e = self.adj.get(p)
        if e is None:
            e = self.adj[p] = Ens()
        return e

    def _positions_de(self, faces):
        e = Ens()
        for f in faces:
            a, b_, c = self.f_pos[f]
            e.add(a)
            e.add(b_)
            e.add(c)
        return e

    def _attr_set(self, p, liste_faces, cle, couche=0):
        """Indices d'attribut (couleur diffuse / spéculaire / coord. texture) utilisés au sommet p."""
        e = Ens()
        for f in self._faces_de(p):
            m = self._mat(self.f_mat[f])
            if cle == 'couches':
                if m['couches'] <= couche:
                    continue
            elif not m[cle]:
                continue
            fp = self.f_pos[f]
            fa = liste_faces[f] if f < len(liste_faces) else None
            if fa is None:
                continue
            if fp[0] == p:
                e.add(fa[0])
            elif fp[1] == p:
                e.add(fa[1])
            elif fp[2] == p:
                e.add(fa[2])
        return e

    # 0xFFFFFF3B
    def base(self, b):
        b.chaine()
        b.u32()
        nf, npos, nnorm, ndiff, nspec, ntex = (b.u32() for _ in range(6))
        for i in range(npos):
            self.positions[i] = (b.f32(), b.f32(), b.f32())
        for _ in range(nnorm):
            b.f32(), b.f32(), b.f32()
        # couleurs gardées dans l'ordre mémoire de la référence (B, G, R, A) : le flux progressif code cet ordre-là
        self.diffus = [(lambda r, g, b_, a: (b_, g, r, a))(b.f32(), b.f32(), b.f32(), b.f32()) for _ in range(ndiff)]
        self.speculaires = [(lambda r, g, b_, a: (b_, g, r, a))(b.f32(), b.f32(), b.f32(), b.f32()) for _ in range(nspec)]
        self.texcoords = [(b.f32(), b.f32(), b.f32(), b.f32()) for _ in range(ntex)]
        for i in range(nf):
            m = b.cu32(C_BASE_SHADING)
            mat = self._mat(m)
            coins = {'pos': [0, 0, 0], 'diff': [0, 0, 0], 'spec': [0, 0, 0], 'tex': [[0, 0, 0] for _ in range(mat['couches'])]}
            for k in range(3):
                coins['pos'][k] = b.cu32(STATIC_FULL + npos)
                if not self.sans_normales:
                    b.cu32(STATIC_FULL + nnorm)
                if mat['diffus']:
                    coins['diff'][k] = b.cu32(STATIC_FULL + ndiff)
                if mat['speculaire']:
                    coins['spec'][k] = b.cu32(STATIC_FULL + nspec)
                for j in range(mat['couches']):
                    coins['tex'][j][k] = b.cu32(STATIC_FULL + ntex)
            self._set_face(self.f_pos, i, coins['pos'])
            self._set_face(self.f_mat, i, m)
            self._set_face(self.f_diff, i, coins['diff'] if mat['diffus'] else None)
            self._set_face(self.f_spec, i, coins['spec'] if mat['speculaire'] else None)
            for j in range(mat['couches']):
                self._set_face(self.f_tex[j], i, coins['tex'][j])
        for i in range(nf):
            for p in self.f_pos[i]:
                self._faces_de(p).add(i)
        self.nb_faces, self.res = nf, npos
        self.nb_diffus, self.nb_spec, self.nb_tex = ndiff, nspec, ntex
        if b.depasse():
            raise ErreurU3D(f'{self.nom} : maillage de base tronqué')

    def _nouveaux_attributs(self, b, n, liste, iq, p_split, cle_liste_faces, cle_mat, ctx_signe, ctx_mag, couche=0):
        pred = [0.0, 0.0, 0.0, 0.0]
        ens = self._attr_set(p_split, cle_liste_faces, cle_mat, couche)
        if len(ens):
            for idx in ens:
                v = liste[idx] if idx < len(liste) else (0.0, 0.0, 0.0, 0.0)
                for k in range(4):
                    pred[k] += v[k]
            pred = [x / len(ens) for x in pred]
        base = len(liste)
        for _ in range(n):
            s = b.cu8(ctx_signe)
            mags = [b.cu32(c) for c in ctx_mag]
            liste.append(tuple(pred[k] + (-iq if s & (1 << k) else iq) * mags[k] for k in range(4)))
        return base

    def _lire_triplet_attr(self, b, ctx_dup, ctx_type, ctx_local, ctx_global, ens_split, ens_third, prev, bug_update_local=False):
        dup = b.cu8(ctx_dup)
        vals = []
        for k, (flag, ens) in enumerate(((1, ens_split), (2, ens_split), (4, ens_third))):
            if dup & flag:
                vals.append(prev[k])
                continue
            t = b.cu8(ctx_type)
            if t == 2:                                  # local
                if k == 1 and bug_update_local:
                    # référence : GetMemberX(0) avant la lecture de l'index (spéculaire uniquement), reproduit tel quel
                    v = ens.membre(0) if len(ens) else 0
                    b.cu32(ctx_local)
                else:
                    li = b.cu32(ctx_local)
                    v = ens.membre(li)
            else:
                v = b.cu32(ctx_global)
            vals.append(v)
        return vals

    # 0xFFFFFF3C
    def progressif(self, b):
        b.chaine()
        b.u32()
        debut = b.u32()
        fin = b.u32()
        if debut != self.res:
            self.avertissements.append(f'bloc progressif {debut}-{fin} alors que la résolution courante est {self.res}')
        prev_diff = [0, 0, 0]
        prev_spec = [0, 0, 0]
        prev_tex = [0, 0, 0]
        for i in range(debut, fin):
            r = i
            split = b.cu32(STATIC_FULL + i)
            faces_split = self._faces_de(split)
            ens_c = self._positions_de(faces_split)
            ens_c.remove(split)

            # nouvelles couleurs / coordonnées de texture
            n = b.cu16(C_NEW_DIFFUSE_COUNT)
            nd_base = self.nb_diffus
            if n:
                self._nouveaux_attributs(b, n, self.diffus, self.iq_diff, split, self.f_diff, 'diffus', C_DIFFUSE_SIGN, C_COLOR_MAG)
            n_new_diff = n
            n = b.cu16(C_NEW_SPECULAR_COUNT)
            ns_base = self.nb_spec
            if n:
                self._nouveaux_attributs(b, n, self.speculaires, self.iq_spec, split, self.f_spec, 'speculaire', C_SPECULAR_SIGN, C_COLOR_MAG)
            n_new_spec = n
            n = b.cu16(C_NEW_TEX_COUNT)
            nt_base = self.nb_tex
            if n:
                self._nouveaux_attributs(b, n, self.texcoords, self.iq_tex, split, self.f_tex[0], 'couches', C_TEX_SIGNS_NEW, C_TEX_MAG, 0)
            n_new_tex = n

            # nouvelles faces
            nnew = b.cu32(C_NUM_NEW_FACES)
            gauche, droite = Ens(), Ens()
            nf = self.nb_faces
            for j in range(nnew):
                sh = b.cu32(C_SHADING)
                self._set_face(self.f_mat, nf + j, sh)
                orient = b.cu8(C_ORIENTATION)
                ttype = b.cu8(C_THIRD_TYPE)
                if ttype == 1:
                    third = ens_c.membre(b.cu32(C_LOCAL_3RD))
                else:
                    third = b.cu32(STATIC_FULL + r)
                    ens_c.add(third)
                if orient == 1:
                    gauche.add(third)
                    face = [split, r, third]
                elif orient == 2:
                    droite.add(third)
                    face = [r, split, third]
                else:
                    raise ErreurU3D(f'{self.nom} : orientation de face invalide ({orient}) à la résolution {r}')
                self._set_face(self.f_pos, nf + j, face)

            # faces qui restent au sommet scindé ou passent au nouveau sommet
            deplacees, pos_dep, pos_reste = Ens(), Ens(), Ens()
            for f in list(faces_split):
                a, b_, c = self.f_pos[f]
                if split == a:
                    x, y = b_, c
                elif split == b_:
                    x, y = c, a
                elif split == c:
                    x, y = a, b_
                else:
                    raise ErreurU3D(f'{self.nom} : face {f} sans le sommet scindé {split}')
                if x in droite:
                    pred = 1
                elif y in droite:
                    pred = 2
                elif x in gauche:
                    pred = 2
                elif y in gauche:
                    pred = 1
                else:
                    pred = 0
                if pred == 0:
                    if a in pos_dep or b_ in pos_dep or c in pos_dep:
                        pred = 3
                    elif a in pos_reste or b_ in pos_reste or c in pos_reste:
                        pred = 4
                sm = b.cu8(C_STAY_MOVE + pred)
                cible = pos_dep if sm == 1 else pos_reste
                if sm == 1:
                    deplacees.add(f)
                cible.add(a)
                cible.add(b_)
                cible.add(c)
                cible.remove(split)

            # mises à jour des faces déplacées (appliquées après la lecture de la position, comme SetResolution)
            maj = []
            for f in deplacees:
                coin = self.f_pos[f].index(split)
                maj.append((self.f_pos, f, coin, r))
                m = self._mat(self.f_mat[f])
                if m['diffus']:
                    if b.cu8(C_DIFF_KEEP) == 1:
                        t = b.cu8(C_DIFF_TYPE)
                        if t == 1:
                            v = b.cu32(C_DIFF_NEW) + nd_base
                        elif t == 2:
                            v = self._attr_set(split, self.f_diff, 'diffus').membre(b.cu32(C_DIFF_LOCAL))
                        else:
                            v = b.cu32(C_DIFF_GLOBAL)
                        maj.append((self.f_diff, f, coin, v))
                if m['speculaire']:
                    if b.cu8(C_SPEC_KEEP) == 1:
                        t = b.cu8(C_SPEC_TYPE)
                        if t == 1:
                            v = b.cu32(C_SPEC_NEW) + ns_base
                        elif t == 2:
                            v = self._attr_set(split, self.f_spec, 'speculaire').membre(b.cu32(C_SPEC_LOCAL))
                        else:
                            v = b.cu32(C_SPEC_GLOBAL)
                        maj.append((self.f_spec, f, coin, v))
                for L in range(m['couches']):
                    if b.cu8(C_TEX_KEEP) == 1:
                        t = b.cu8(C_TEX_TYPE)
                        if t == 1:
                            v = b.cu32(C_TEX_NEW) + nt_base
                        elif t == 2:
                            v = self._attr_set(split, self.f_tex[L], 'couches', L).membre(b.cu32(C_TEX_LOCALI))
                        else:
                            v = b.cu32(C_TEX_GLOBALI)
                        maj.append((self.f_tex[L], f, coin, v))

            # attributs des nouvelles faces
            for j in range(nnew):
                face = self.f_pos[nf + j]
                c_new = face.index(r)
                c_split = face.index(split)
                c_third = 3 - c_new - c_split
                third = face[c_third]
                m = self._mat(self.f_mat[nf + j])
                if m['diffus']:
                    v = self._lire_triplet_attr(b, C_COLOR_DUP, C_COLOR_SPLIT_TYPE, C_COLOR_LOCAL, C_COLOR_GLOBAL,
                                                self._attr_set(split, self.f_diff, 'diffus'),
                                                self._attr_set(third, self.f_diff, 'diffus'), prev_diff)
                    fa = [0, 0, 0]
                    fa[c_split], fa[c_new], fa[c_third] = v
                    self._set_face(self.f_diff, nf + j, fa)
                    prev_diff = v
                else:
                    self._set_face(self.f_diff, nf + j, None)
                if m['speculaire']:
                    v = self._lire_triplet_attr(b, C_COLOR_DUP, C_COLOR_SPLIT_TYPE, C_COLOR_LOCAL, C_COLOR_GLOBAL,
                                                self._attr_set(split, self.f_spec, 'speculaire'),
                                                self._attr_set(third, self.f_spec, 'speculaire'), prev_spec,
                                                bug_update_local=True)
                    fa = [0, 0, 0]
                    fa[c_split], fa[c_new], fa[c_third] = v
                    self._set_face(self.f_spec, nf + j, fa)
                    prev_spec = v
                else:
                    self._set_face(self.f_spec, nf + j, None)
                for L in range(m['couches']):
                    v = self._lire_triplet_attr(b, C_TEX_DUP, C_TEX_SPLIT_TYPE, C_TEX_LOCAL, C_TEX_GLOBAL,
                                                self._attr_set(split, self.f_tex[L], 'couches', L),
                                                self._attr_set(third, self.f_tex[L], 'couches', L), prev_tex)
                    fa = [0, 0, 0]
                    fa[c_split], fa[c_new], fa[c_third] = v
                    self._set_face(self.f_tex[L], nf + j, fa)
                    prev_tex = v

            # adjacence : faces déplacées puis nouvelles faces
            fs = self._faces_de(split)
            fr = self._faces_de(r)
            for f in list(deplacees):
                fs.remove(f)
                fr.add(f)
            for j in range(nnew):
                for p in self.f_pos[nf + j]:
                    self._faces_de(p).add(nf + j)
            ens_j = self._positions_de(self._faces_de(r))   # coins pas encore mis à jour (comme la référence)

            # position du nouveau sommet
            pred = self.positions[split] if (r > 0 and split < len(self.positions) and self.positions[split]) else (0.0, 0.0, 0.0)
            s = b.cu8(C_POS_SIGNS)
            dx, dy, dz = b.cu32(C_POS_X), b.cu32(C_POS_Y), b.cu32(C_POS_Z)
            q = self.iq_pos
            if r >= len(self.positions):
                self.positions.extend([None] * (r + 1 - len(self.positions)))
            self.positions[r] = (pred[0] + (-q if s & 1 else q) * dx,
                                 pred[1] + (-q if s & 2 else q) * dy,
                                 pred[2] + (-q if s & 4 else q) * dz)

            # normales : on lit le flux (nombres de symboles dictés par la topologie) sans reconstruire les valeurs
            if not self.sans_normales:
                for p in ens_j:
                    nloc = b.cu32(C_NUM_LOCAL_NORMALS)
                    for _ in range(nloc):
                        b.cu8(C_NORMAL_SIGNS)
                        b.cu32(C_NORMAL_X)
                        b.cu32(C_NORMAL_Y)
                        b.cu32(C_NORMAL_Z)
                    for _ in range(len(self._faces_de(p))):
                        b.cu32(C_NORMAL_LOCAL)

            # SetResolution(r + 1) : mises à jour des faces et compteurs
            for liste, f, coin, v in maj:
                liste[f][coin] = v
            self.nb_faces += nnew
            self.nb_diffus = nd_base + n_new_diff
            self.nb_spec = ns_base + n_new_spec
            self.nb_tex = nt_base + n_new_tex
            self.res = r + 1
            if b.depasse():
                raise ErreurU3D(f'{self.nom} : flux progressif épuisé à la résolution {r} (désynchronisation)')

    def resultat(self):
        """(sommets Nx3, faces Mx3, couleurs par sommet Nx4 ou None) à la résolution courante."""
        n = self.res
        pos = np.array([p if p is not None else (np.nan,) * 3 for p in self.positions[:n]], dtype=np.float64).reshape(-1, 3)
        faces = np.array(self.f_pos[:self.nb_faces], dtype=np.int64).reshape(-1, 3)
        if len(faces):
            ok = (faces < n).all(1) & (faces[:, 0] != faces[:, 1]) & (faces[:, 1] != faces[:, 2]) & (faces[:, 0] != faces[:, 2])
            faces = faces[ok]
        couleurs = None
        if self.diffus and any(fd is not None for fd in self.f_diff[:self.nb_faces]):
            acc = np.zeros((n, 4))
            cnt = np.zeros(n)
            for fp, fd in zip(self.f_pos[:self.nb_faces], self.f_diff[:self.nb_faces]):
                if fd is None:
                    continue
                for k in range(3):
                    if fp[k] < n and fd[k] < len(self.diffus):
                        acc[fp[k]] += self.diffus[fd[k]]
                        cnt[fp[k]] += 1
            ok = cnt > 0
            acc[ok] /= cnt[ok, None]
            couleurs = acc[:, [2, 1, 0, 3]]           # BGRA -> RGBA
        return pos, faces, couleurs


# --- jeu de lignes (CIFXLineSetDecoder) -----------------------------------------------------------------------
class JeuDeLignes:
    """Segments (p. ex. nerfs ou axes tracés en fil). Positions et segments reconstruits ; normales, couleurs et
    coordonnées de texture lues pour rester synchronisé mais non conservées."""

    def __init__(self, nom):
        self.nom = nom
        self.pret = False
        self.avertissements = []

    def declaration(self, b):
        b.chaine()
        b.u32()
        b.u32()                                   # paramètre réservé
        (self.max_lignes, self.max_positions, self.max_normales, self.max_diffus, self.max_speculaires,
         self.max_tex, nmat) = (b.u32() for _ in range(7))
        self.materiaux = []
        for _ in range(nmat):
            a = b.u32()
            nc = b.u32()
            [b.u32() for _ in range(nc)]
            b.u32()
            self.materiaux.append({'diffus': bool(a & 1), 'speculaire': bool(a & 2), 'couches': nc})
        [b.u32() for _ in range(3)]
        self.iq_pos = b.f32()
        [b.f32() for _ in range(4)]
        [b.u32() for _ in range(3)]
        self.positions = []
        self.segments = []
        self.pret = True

    def continuation(self, b):
        b.chaine()
        b.u32()
        debut, fin = b.u32(), b.u32()
        q = self.iq_pos
        for _ in range(debut, fin):
            cur = len(self.positions)
            split = b.cu32(STATIC_FULL + (cur if cur > 1 else 1))
            s = b.cu8(C_POS_SIGNS)
            dx, dy, dz = b.cu32(C_POS_X), b.cu32(C_POS_Y), b.cu32(C_POS_Z)
            pred = self.positions[split] if (cur > 0 and split < cur) else (0.0, 0.0, 0.0)
            self.positions.append((pred[0] + (-q if s & 1 else q) * dx, pred[1] + (-q if s & 2 else q) * dy,
                                   pred[2] + (-q if s & 4 else q) * dz))
            for _ in range(b.cu32(C_NUM_LOCAL_NORMALS)):
                b.cu8(C_NORMAL_SIGNS); b.cu32(C_NORMAL_X); b.cu32(C_NORMAL_Y); b.cu32(C_NORMAL_Z)
            for _ in range(b.cu32(C_NUM_NEW_FACES)):
                m = b.cu32(1)
                fin_seg = b.cu32(STATIC_FULL + cur)
                mat = self.materiaux[m] if m < len(self.materiaux) else {'diffus': False, 'speculaire': False, 'couches': 0}
                for _ in range(2):
                    b.cu32(C_NORMAL_LOCAL)
                    for cle in ('diffus', 'speculaire'):
                        if mat[cle] and b.cu8(C_COLOR_DUP) == 0:
                            b.cu8(C_NORMAL_SIGNS)
                            for c in C_COLOR_MAG:
                                b.cu32(c)
                    if self.max_tex:
                        for _ in range(mat['couches']):
                            if b.cu8(C_TEX_DUP) == 0:
                                b.cu8(32)
                                for c in C_TEX_MAG:
                                    b.cu32(c)
                self.segments.append((fin_seg, cur))
            if b.depasse():
                raise ErreurU3D(f'{self.nom} : flux du jeu de lignes épuisé ({cur} positions)')

    def resultat(self):
        return np.array(self.positions, dtype=np.float64).reshape(-1, 3), np.array(self.segments, dtype=np.int64).reshape(-1, 2)


# --- lecture des blocs --------------------------------------------------------------------------------------------
def _pad4(n):
    return (n + 3) & ~3


def iter_blocs(buf, debut=0, fin=None):
    fin = len(buf) if fin is None else fin
    p = debut
    while p + 12 <= fin:
        t, ds, ms = struct.unpack_from('<III', buf, p)
        d0 = p + 12
        m0 = d0 + _pad4(ds)
        if d0 + ds > fin:
            raise ErreurU3D(f'bloc 0x{t:08X} tronqué à l\'octet {p}')
        yield t, buf[d0:d0 + ds], buf[m0:m0 + ms], p
        p = m0 + _pad4(ms)


def _lire_metadonnees(meta):
    if not meta:
        return {}
    b = FluxBits(meta)
    res = {}
    for _ in range(b.u32()):
        attr = b.u32()
        cle = b.chaine()
        if attr & 1:
            n = b.u32()
            res[cle] = bytes(b.u8() for _ in range(n))
        else:
            res[cle] = b.chaine()
    return res


def _lire_parents(b):
    parents = []
    for _ in range(b.u32()):
        nom = b.chaine()
        m = [b.f32() for _ in range(16)]
        parents.append((nom, np.array(m, dtype=np.float64).reshape(4, 4).T))   # stockage colonne par colonne
    return parents


class U3DScene:
    def __init__(self):
        self.version = None
        self.profil = 0
        self.unites = None
        self.noeuds = {}          # nom -> {'type', 'parents': [(nom, M)], 'ressource', 'visibilite', 'meta'}
        self.maillages = {}       # nom de ressource -> MaillageCLOD
        self.lignes = {}          # nom de ressource -> JeuDeLignes
        self.materiaux = {}       # nom -> {'ambiant','diffus','speculaire','emissif','reflectivite','opacite'}
        self.shaders = {}         # nom -> {'materiau': nom}
        self.ombrages = {}        # nom de nœud / ressource -> [[shaders de la liste 0], ...]
        self.meta_blocs = []
        self.ignores = {}
        self.erreurs = []

    def _ignorer(self, t):
        k = NOMS_BLOCS.get(t, f'0x{t:08X}')
        self.ignores[k] = self.ignores.get(k, 0) + 1

    def lire(self, buf, debut=0, fin=None, type_chaine=None):
        for t, data, meta, off in iter_blocs(buf, debut, fin):
            md = {}
            try:
                md = _lire_metadonnees(meta)
            except Exception as e:                      # métadonnées illisibles : pas bloquant
                self.erreurs.append(f'métadonnées du bloc à {off} : {e}')
            if md:
                self.meta_blocs.append((NOMS_BLOCS.get(t, hex(t)), md))
            try:
                self._bloc(t, data, md, type_chaine)
            except ErreurU3D as e:
                self.erreurs.append(str(e))
            except (IndexError, ValueError, ZeroDivisionError, struct.error) as e:
                self.erreurs.append(f'bloc {NOMS_BLOCS.get(t, hex(t))} à l\'octet {off} : {type(e).__name__} {e}')

    def _flux(self, data):
        return FluxBits(data, bool(self.profil & PROFILE_NOCOMPRESSION))

    def _bloc(self, t, data, md, type_chaine):
        if t == B_HEADER:
            b = FluxBits(data)
            self.version = b.u32()
            self.profil = b.u32()
            b.u32()
            b.u64()
            b.u32()
            if self.profil & PROFILE_UNITSSCALE:
                self.unites = b.f64()
        elif t == B_MODIFIER_CHAIN:
            b = FluxBits(data)
            nom = b.chaine()
            typ = b.u32()
            attrs = b.u32()
            p = 2 + len(nom.encode('utf-8')) + 8
            if attrs & 1:
                [b.f32() for _ in range(4)]
                p += 16
            if attrs & 2:
                [b.f32() for _ in range(6)]
                p += 24
            pad = (4 - ((2 + len(nom.encode('utf-8'))) & 3)) & 3
            for _ in range(pad):
                b.u8()
            p += pad
            b.u32()                                     # nombre de modificateurs
            p += 4
            self.lire(data, p, len(data), typ)
        elif t in (B_GROUP, B_MODEL, B_LIGHT, B_VIEW):
            b = self._flux(data)
            nom = b.chaine()
            parents = _lire_parents(b)
            n = {'type': {B_GROUP: 'groupe', B_MODEL: 'modele', B_LIGHT: 'lumiere', B_VIEW: 'vue'}[t],
                 'parents': parents, 'meta': md}
            if t in (B_MODEL, B_LIGHT, B_VIEW):
                n['ressource'] = b.chaine()
            if t == B_MODEL:
                n['visibilite'] = b.u32()
            self.noeuds[nom] = n
        elif t == B_CLOD_DECL:
            b = self._flux(data)
            nom = b.chaine()
            m = self.maillages.get(nom) or MaillageCLOD(nom)
            self.maillages[nom] = m
            m.meta = md
            m.declaration(self._flux(data))
        elif t in (B_CLOD_BASE, B_CLOD_PROG):
            nom = self._flux(data).chaine()
            m = self.maillages.get(nom)
            if m is None or not m.pret:
                raise ErreurU3D(f'suite de maillage « {nom} » sans déclaration')
            (m.base if t == B_CLOD_BASE else m.progressif)(self._flux(data))
        elif t == B_LINESET:
            nom = self._flux(data).chaine()
            j = self.lignes.get(nom) or JeuDeLignes(nom)
            self.lignes[nom] = j
            j.declaration(self._flux(data))
        elif t == B_LINESET_CONT:
            nom = self._flux(data).chaine()
            j = self.lignes.get(nom)
            if j is None or not j.pret:
                raise ErreurU3D(f'suite de jeu de lignes « {nom} » sans déclaration')
            j.continuation(self._flux(data))
        elif t == B_MATERIAL:
            b = self._flux(data)
            nom = b.chaine()
            b.u32()
            c = [tuple(b.f32() for _ in range(3)) for _ in range(4)]
            self.materiaux[nom] = {'ambiant': c[0], 'diffus': c[1], 'speculaire': c[2], 'emissif': c[3],
                                   'reflectivite': b.f32(), 'opacite': b.f32(), 'meta': md}
        elif t == B_SHADER:
            b = self._flux(data)
            nom = b.chaine()
            b.u32(); b.f32(); b.u32(); b.u32(); b.u32()
            b.u32()                                     # canaux de texture
            b.u32()
            self.shaders[nom] = {'materiau': b.chaine(), 'meta': md}
        elif t == B_SHADING:
            b = self._flux(data)
            nom = b.chaine()
            b.u32()                                     # position dans la chaîne
            b.u32()                                     # attributs
            listes = []
            for _ in range(b.u32()):
                listes.append([b.chaine() for _ in range(b.u32())])
            self.ombrages[nom] = listes
        elif t == B_PRIORITY:
            pass
        else:
            self._ignorer(t)

    # --- scène aplatie ---
    def matrice_monde(self, nom, parent_idx=0, _prof=0):
        n = self.noeuds.get(nom)
        if n is None or not n['parents'] or _prof > 64:
            return [np.eye(4)]
        res = []
        for pnom, M in n['parents']:
            if pnom in ('', '<NULL>') or pnom not in self.noeuds:
                res.append(M)
            else:
                for Mp in self.matrice_monde(pnom, _prof=_prof + 1):
                    res.append(Mp @ M)
        return res

    def materiau_de(self, noeud, ressource):
        for cle in (noeud, ressource):
            listes = self.ombrages.get(cle)
            if listes:
                for liste in listes:
                    for sh in liste:
                        mat = self.shaders.get(sh, {}).get('materiau')
                        if mat in self.materiaux:
                            return mat, self.materiaux[mat]
        # à défaut : shader ou matériau portant le nom du nœud
        for cle in (noeud, ressource):
            mat = self.shaders.get(cle, {}).get('materiau')
            if mat in self.materiaux:
                return mat, self.materiaux[mat]
            if cle in self.materiaux:
                return cle, self.materiaux[cle]
        return None, None

    def chemin(self, nom, _prof=0):
        """Noms des groupes ancêtres (premier parent), du plus haut au plus bas."""
        n = self.noeuds.get(nom)
        if n is None or not n['parents'] or _prof > 64:
            return []
        p = n['parents'][0][0]
        if p in ('', '<NULL>') or p not in self.noeuds:
            return []
        return self.chemin(p, _prof + 1) + [p]

    def objets(self):
        out = []
        for nom, n in self.noeuds.items():
            if n['type'] != 'modele':
                continue
            m = self.maillages.get(n.get('ressource'))
            if m is None or not m.pret:
                j = self.lignes.get(n.get('ressource'))
                if j is not None and j.pret and j.segments:
                    v, seg = j.resultat()
                    mat_nom, mat = self.materiau_de(nom, n.get('ressource'))
                    for k, M in enumerate(self.matrice_monde(nom)):
                        out.append({'nom': nom if k == 0 else f'{nom}#{k}', 'noeud': nom, 'ressource': n.get('ressource'),
                                    'groupes': self.chemin(nom), 'visible': n.get('visibilite', 3) != 0,
                                    'sommets': v @ M[:3, :3].T + M[:3, 3], 'faces': np.zeros((0, 3), np.int64), 'segments': seg,
                                    'couleurs_sommets': None, 'materiau': mat_nom, 'couleur': couleur_visible(mat),
                                    'opacite': mat['opacite'] if mat else 1.0, 'avertissements': list(j.avertissements)})
                continue
            v, f, col = m.resultat()
            if not len(f):
                continue
            mat_nom, mat = self.materiau_de(nom, n.get('ressource'))
            for k, M in enumerate(self.matrice_monde(nom)):
                vw = v @ M[:3, :3].T + M[:3, 3]
                out.append({'nom': nom if k == 0 else f'{nom}#{k}', 'noeud': nom, 'ressource': n.get('ressource'),
                            'groupes': self.chemin(nom), 'visible': n.get('visibilite', 3) != 0,
                            'sommets': vw, 'faces': f, 'couleurs_sommets': col,
                            'materiau': mat_nom, 'couleur': couleur_visible(mat),
                            'opacite': mat['opacite'] if mat else 1.0,
                            'avertissements': list(m.avertissements),
                            'resolution': (m.res, m.res_max_finale)})
        return out

    def resume(self):
        lignes = [f'U3D version 0x{self.version or 0:08X}, profil 0x{self.profil:X}'
                  + (f', unités {self.unites} m' if self.unites else '')]
        lignes.append(f'{sum(1 for n in self.noeuds.values() if n["type"] == "modele")} nœuds modèle, '
                      f'{sum(1 for n in self.noeuds.values() if n["type"] == "groupe")} groupes, '
                      f'{len(self.maillages)} maillages, {len(self.lignes)} jeux de lignes, {len(self.materiaux)} matériaux')
        for nom, m in self.maillages.items():
            if m.pret:
                lignes.append(f'  maillage {nom!r} : résolution {m.res}/{m.res_max_finale}, {m.nb_faces}/{m.max_faces} faces'
                              + (f' ; {"; ".join(m.avertissements)}' if m.avertissements else ''))
        for nom, j in self.lignes.items():
            if j.pret:
                lignes.append(f'  jeu de lignes {nom!r} : {len(j.positions)}/{j.max_positions} positions, {len(j.segments)}/{j.max_lignes} segments')
        if self.ignores:
            lignes.append('blocs ignorés : ' + ', '.join(f'{k} ×{v}' for k, v in self.ignores.items()))
        for e in self.erreurs:
            lignes.append('ERREUR ' + e)
        return '\n'.join(lignes)


def couleur_visible(mat):
    """Couleur diffuse du matériau ; si elle est noire, la composante émissive ou ambiante (certains exports y
    mettent la couleur)."""
    if not mat:
        return None
    for cle in ('diffus', 'emissif', 'ambiant'):
        if sum(mat[cle]) > 0.05:
            return tuple(mat[cle])
    return tuple(mat['diffus'])


def lire_u3d(buf):
    if buf[:4] != b'U3D\x00':
        raise ErreurU3D('signature U3D absente')
    s = U3DScene()
    s.lire(bytes(buf))
    return s


if __name__ == '__main__':
    import sys
    for chemin in sys.argv[1:]:
        sc = lire_u3d(open(chemin, 'rb').read())
        print(chemin)
        print(sc.resume())
        for o in sc.objets():
            print(f"  {o['nom']:30s} {len(o['sommets']):7d} sommets {len(o['faces']):7d} faces  couleur={o['couleur']}")
