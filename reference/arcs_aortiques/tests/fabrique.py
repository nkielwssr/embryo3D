"""Fabrique de fichiers de test : U3D non compressé (profil 0x4, maillage de base seul) et PDF 3D minimal.

Sert à tester, sans logiciel externe, le chemin « maillage de base » du décodeur U3D et l'extraction des flux 3D
d'un PDF. Les U3D compressés (maillage progressif) des tests viennent du convertisseur IDTF de référence
(fichiers `donnees/*.u3d`, voir test_u3d.py).
"""
import struct
import zlib


def _s(txt):
    b = txt.encode('utf-8')
    return struct.pack('<H', len(b)) + b


def _bloc(t, data, meta=b''):
    pad = lambda b: b + b'\x00' * ((4 - len(b) % 4) % 4)
    return struct.pack('<III', t, len(data), len(meta)) + pad(data) + pad(meta)


def _parents(parent, M):
    """M : matrice 4x4 (liste de lignes, convention colonne : translation en dernière colonne)."""
    cols = [M[r][c] for c in range(4) for r in range(4)]           # stockage colonne par colonne
    return struct.pack('<I', 1) + _s(parent) + struct.pack('<16f', *cols)


def _chaine(nom, typ, blocs):
    tete = _s(nom) + struct.pack('<II', typ, 0)
    pad = (4 - ((2 + len(nom.encode('utf-8'))) & 3)) & 3
    return _bloc(0xFFFFFF14, tete + b'\x00' * pad + struct.pack('<I', len(blocs)) + b''.join(blocs))


def u3d_base(objets, groupes=()):
    """objets : dicts {nom, sommets [(x,y,z)], faces [(a,b,c)], couleur (r,g,b), parent, M, couleurs_sommets}.
    Écrit un U3D au profil « sans compression » : chaque maillage n'a qu'un bloc de base (résolution minimale = max)."""
    I = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
    out = [_bloc(0x00443355, struct.pack('<IIIQI', 0, 0x4, 0, 0, 106))]
    for g in groupes:
        out.append(_chaine(g['nom'], 0, [_bloc(0xFFFFFF21, _s(g['nom']) + _parents(g.get('parent', ''), g.get('M', I)))]))
    for o in objets:
        n = o['nom']
        noeud = _bloc(0xFFFFFF22, _s(n) + _parents(o.get('parent', ''), o.get('M', I)) + _s(n + '_res') + struct.pack('<I', 3))
        omb = _bloc(0xFFFFFF45, _s(n) + struct.pack('<III', 1, 1, 1) + struct.pack('<I', 1) + _s(n + '_shader'))
        out.append(_chaine(n, 0, [noeud, omb]))
        V, F = o['sommets'], o['faces']
        vc = o.get('couleurs_sommets')
        nd = len(vc) if vc is not None else 0
        decl = (_s(n + '_res') + struct.pack('<I', 0) + struct.pack('<I', 1)            # attributs : sans normales
                + struct.pack('<7I', len(F), len(V), 0, nd, 0, 0, 1)
                + struct.pack('<III', 1 if nd else 0, 0, 0)                              # matériau 0
                + struct.pack('<II', len(V), len(V))                                     # résolution min = max
                + struct.pack('<III', 1000, 1000, 1000) + struct.pack('<5f', 1e-3, 1e-3, 1e-3, 1e-3, 1e-3)
                + struct.pack('<3f', 0.25, 0.5, 0.985) + struct.pack('<I', 0))
        base = _s(n + '_res') + struct.pack('<I', 0) + struct.pack('<6I', len(F), len(V), 0, nd, 0, 0)
        base += b''.join(struct.pack('<3f', *v) for v in V)
        base += b''.join(struct.pack('<4f', *c, 1.0) for c in (vc or []))
        for f in F:
            base += struct.pack('<I', 0)                                                   # matériau de la face
            for k in range(3):
                base += struct.pack('<I', f[k])
                if nd:
                    base += struct.pack('<I', f[k])
        out.append(_chaine(n + '_res', 1, [_bloc(0xFFFFFF31, decl), _bloc(0xFFFFFF3B, base)]))
        r, g, b = o.get('couleur', (0.8, 0.8, 0.8))
        out.append(_bloc(0xFFFFFF53, _s(n + '_shader') + struct.pack('<IfIIIII', 0, 0.0, 0, 0, 1, 0, 0) + _s(n + '_mat')))
        out.append(_bloc(0xFFFFFF54, _s(n + '_mat') + struct.pack('<I', 0x3F)
                         + struct.pack('<12f', r * .2, g * .2, b * .2, r, g, b, .2, .2, .2, 0, 0, 0)
                         + struct.pack('<ff', 0.1, o.get('opacite', 1.0))))
    return b''.join(out)


def pdf_3d(flux, sous_type=b'U3D', vues=(), compresser=True, titre='test'):
    """PDF minimal d'une page portant une annotation 3D (flux U3D ou PRC). vues : [(nom, {noeud: visible})]."""
    objs = []

    def ajouter(corps):
        objs.append(corps)
        return len(objs)

    donnees = zlib.compress(flux) if compresser else flux
    filtre = b'/Filter /FlateDecode ' if compresser else b''
    vues_refs = []
    for nom, vis in vues:
        noeuds = b' '.join(b'<< /Type /3DNode /N (' + k.encode() + b') /V ' + (b'true' if v else b'false') + b' >>'
                           for k, v in vis.items())
        vues_refs.append(ajouter(b'<< /Type /3DView /XN (' + nom.encode() + b') /IN (' + nom.encode()
                                 + b') /MS /M /C2W [1 0 0 0 1 0 0 0 1 0 0 -10] /NA [' + noeuds + b'] >>'))
    vn = b' '.join(b'%d 0 R' % r for r in vues_refs)
    i3d = ajouter(b'<< /Type /3D /Subtype /' + sous_type + b' ' + filtre + b'/Length %d ' % len(donnees)
                  + (b'/VA [' + vn + b'] ' if vues_refs else b'') + b'>>\nstream\n' + donnees + b'\nendstream')
    annot = ajouter(b'<< /Type /Annot /Subtype /3D /Rect [0 0 400 400] /3DD %d 0 R /Contents (%s) >>' % (i3d, titre.encode()))
    page = ajouter(b'<< /Type /Page /Parent 5 0 R /MediaBox [0 0 400 400] /Annots [%d 0 R] >>' % annot)
    cat = ajouter(b'<< /Type /Catalog /Pages 5 0 R >>')
    assert len(objs) == len(vues_refs) + 4
    pages_num = ajouter(b'<< /Type /Pages /Kids [%d 0 R] /Count 1 >>' % page)
    # renumérotation : l'objet Pages doit être le 5 dans les références ci-dessus -> on corrige
    corps = [o.replace(b'/Parent 5 0 R', b'/Parent %d 0 R' % pages_num).replace(b'/Pages 5 0 R', b'/Pages %d 0 R' % pages_num)
             for o in objs]
    out = bytearray(b'%PDF-1.6\n%\xe2\xe3\xcf\xd3\n')
    offs = []
    for i, o in enumerate(corps, 1):
        offs.append(len(out))
        out += b'%d 0 obj\n' % i + o + b'\nendobj\n'
    xref = len(out)
    out += b'xref\n0 %d\n0000000000 65535 f \n' % (len(corps) + 1)
    for o in offs:
        out += b'%010d 00000 n \n' % o
    out += b'trailer\n<< /Size %d /Root %d 0 R >>\nstartxref\n%d\n%%%%EOF\n' % (len(corps) + 1, cat, xref)
    return bytes(out)
