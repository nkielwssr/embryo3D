"""Extraction du contenu 3D d'un PDF interactif (« PDF 3D ») : flux U3D / PRC, vues prédéfinies, JavaScript,
texte des pages (légendes).

Deux voies complémentaires :
  1. pypdf (si installé) : annotations /3D et /RichMedia de chaque page, vues (/VA), JavaScript (/OnInstantiate),
     texte des pages ; gère les PDF chiffrés sans mot de passe utilisateur (cas des suppléments d'éditeurs).
  2. balayage brut du fichier : tout flux dont le contenu (décompressé si besoin) commence par « U3D\\0 » ou
     « PRC » ; rattrape les PDF abîmés et les flux que la première voie n'atteint pas.
Les flux identiques (même empreinte) ne sont gardés qu'une fois.
"""
import hashlib
import re
import zlib

MAGIES = {b'U3D\x00': 'U3D', b'PRC': 'PRC'}


def _format(data):
    for m, f in MAGIES.items():
        if data[:len(m)] == m:
            return f
    return None


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _pyobj(o, prof=0):
    """Objet pypdf -> types Python simples (pour JSON)."""
    if prof > 8:
        return None
    try:
        from pypdf.generic import (ArrayObject, BooleanObject, DictionaryObject, FloatObject, IndirectObject,
                                   NameObject, NumberObject, StreamObject, TextStringObject, ByteStringObject)
    except ImportError:
        return str(o)
    if isinstance(o, IndirectObject):
        try:
            o = o.get_object()
        except Exception:
            return None
    if isinstance(o, StreamObject):
        return '<flux>'
    if isinstance(o, DictionaryObject):
        return {str(k): _pyobj(v, prof + 1) for k, v in o.items() if k not in ('/Parent', '/P')}
    if isinstance(o, ArrayObject):
        return [_pyobj(v, prof + 1) for v in o]
    if isinstance(o, BooleanObject):
        return bool(getattr(o, 'value', o == True))    # bool(BooleanObject(False)) vaut True dans pypdf
    if isinstance(o, (FloatObject, NumberObject)):
        return float(o)
    if isinstance(o, NameObject):
        return str(o)
    if isinstance(o, (TextStringObject, ByteStringObject)):
        return str(o)
    return str(o)


def _vues(flux_dict):
    vues = []
    va = flux_dict.get('/VA')
    if va is None:
        return vues
    try:
        va = va.get_object()
    except AttributeError:
        pass
    for v in va or []:
        try:
            v = v.get_object()
        except AttributeError:
            pass
        d = _pyobj(v)
        if not isinstance(d, dict):
            continue
        vue = {'nom': d.get('/XN') or d.get('/IN') or '', 'nom_interne': d.get('/IN'),
               'camera_vers_monde': d.get('/C2W'), 'centre_orbite': d.get('/CO'), 'noeuds': []}
        for n in d.get('/NA') or []:
            if isinstance(n, dict):
                vue['noeuds'].append({'nom': n.get('/N'), 'visible': n.get('/V'), 'opacite': n.get('/O'),
                                      'matrice': n.get('/M'), 'rendu': n.get('/RM')})
        vues.append(vue)
    return vues


def _javascript(flux_dict):
    js = flux_dict.get('/OnInstantiate')
    if js is None:
        return None
    try:
        js = js.get_object()
        data = js.get_data() if hasattr(js, 'get_data') else bytes(str(js), 'latin-1')
        return data.decode('utf-8', 'replace')
    except Exception:
        return None


def _par_pypdf(chemin, res):
    try:
        from pypdf import PdfReader
    except (KeyboardInterrupt, SystemExit):
        raise
    except BaseException as e:          # ImportError, ou module de chiffrement cassé (panique pyo3)
        res['erreurs'].append(f'pypdf indisponible ({type(e).__name__}) : seul le balayage brut est utilisé (pip install pypdf)')
        return
    try:
        r = PdfReader(chemin, strict=False)
        if r.is_encrypted:
            try:
                r.decrypt('')
            except Exception as e:
                res['erreurs'].append(f'PDF chiffré : {e} (installer « cryptography » pour l\'AES)')
                return
    except Exception as e:
        res['erreurs'].append(f'pypdf : {e}')
        return
    res['pages'] = len(r.pages)
    try:
        info = r.metadata or {}
        res['metadonnees'] = {str(k): str(v) for k, v in info.items()}
    except Exception:
        pass
    for ip, page in enumerate(r.pages):
        try:
            res['texte_pages'].append(page.extract_text() or '')
        except Exception as e:
            res['texte_pages'].append('')
            res['erreurs'].append(f'texte page {ip + 1} : {e}')
        try:
            annots = page.get('/Annots') or []
            annots = annots.get_object() if hasattr(annots, 'get_object') else annots
        except Exception:
            annots = []
        for a in annots:
            try:
                a = a.get_object()
                st = a.get('/Subtype')
                if st == '/3D':
                    d3 = a.get('/3DD')
                    d3 = d3.get_object() if d3 is not None else None
                    if d3 is not None and d3.get('/Type') == '/3DRef':
                        d3 = d3['/3DD'].get_object()
                    if d3 is None or not hasattr(d3, 'get_data'):
                        continue
                    data = d3.get_data()
                    res['flux3d'].append({'page': ip + 1, 'source': 'annotation 3D', 'format': _format(data) or str(d3.get('/Subtype')),
                                          'donnees': data, 'vues': _vues(d3), 'javascript': _javascript(d3),
                                          'titre': str(a.get('/Contents') or a.get('/T') or ''),
                                          'vue_par_defaut': _pyobj(d3.get('/DV')) if d3.get('/DV') is not None else None})
                elif st == '/RichMedia':
                    contenu = a['/RichMediaContent'].get_object()
                    noms = contenu['/Assets'].get_object().get('/Names') or []
                    for k in range(0, len(noms) - 1, 2):
                        spec = noms[k + 1].get_object()
                        ef = spec['/EF'].get_object()
                        fl = (ef.get('/F') or ef.get('/UF')).get_object()
                        data = fl.get_data()
                        if _format(data):
                            res['flux3d'].append({'page': ip + 1, 'source': f'RichMedia « {noms[k]} »', 'format': _format(data),
                                                  'donnees': data, 'vues': [], 'javascript': None, 'titre': str(noms[k])})
            except Exception as e:
                res['erreurs'].append(f'annotation page {ip + 1} : {type(e).__name__} {e}')
    # JavaScript de niveau document (souvent : légende, bascule des structures)
    try:
        noms = r.trailer['/Root'].get_object().get('/Names')
        if noms is not None:
            js = noms.get_object().get('/JavaScript')
            if js is not None:
                arbre = js.get_object().get('/Names') or []
                for k in range(0, len(arbre) - 1, 2):
                    act = arbre[k + 1].get_object()
                    code = act.get('/JS')
                    code = code.get_object() if hasattr(code, 'get_object') else code
                    code = code.get_data().decode('utf-8', 'replace') if hasattr(code, 'get_data') else str(code)
                    res['javascript_document'].append({'nom': str(arbre[k]), 'code': code})
    except Exception as e:
        res['erreurs'].append(f'JavaScript du document : {e}')


_RE_STREAM = re.compile(rb'stream\r?\n')


def _par_balayage(buf, res):
    for m in _RE_STREAM.finditer(buf):
        d0 = m.end()
        fin = buf.find(b'endstream', d0)
        if fin < 0:
            continue
        brut = buf[d0:fin].rstrip(b'\r\n')
        tete = buf[max(0, m.start() - 600):m.start()]
        cands = [brut]
        if b'FlateDecode' in tete or not _format(brut):
            try:
                cands.insert(0, zlib.decompressobj().decompress(brut))
            except zlib.error:
                pass
        for data in cands:
            f = _format(data)
            if f:
                res['flux3d'].append({'page': None, 'source': f'balayage (octet {d0})', 'format': f, 'donnees': data,
                                      'vues': [], 'javascript': None, 'titre': ''})
                break


def extraire_pdf(chemin):
    buf = open(chemin, 'rb').read()
    res = {'fichier': chemin, 'pages': None, 'metadonnees': {}, 'texte_pages': [], 'flux3d': [],
           'javascript_document': [], 'erreurs': []}
    if not buf.startswith(b'%PDF'):
        f = _format(buf)
        if f:                                    # fichier U3D / PRC nu
            res['flux3d'].append({'page': None, 'source': 'fichier', 'format': f, 'donnees': buf, 'vues': [],
                                  'javascript': None, 'titre': ''})
            return res
        res['erreurs'].append('ni PDF, ni U3D, ni PRC')
        return res
    _par_pypdf(chemin, res)
    vus = {hashlib.sha1(f['donnees']).hexdigest() for f in res['flux3d']}
    tmp = {'flux3d': []}
    _par_balayage(buf, tmp)
    for f in tmp['flux3d']:
        h = hashlib.sha1(f['donnees']).hexdigest()
        if h not in vus:
            vus.add(h)
            res['flux3d'].append(f)
    for f in res['flux3d']:
        f['sha1'] = hashlib.sha1(f['donnees']).hexdigest()
        f['octets'] = len(f['donnees'])
    return res


if __name__ == '__main__':
    import sys
    for c in sys.argv[1:]:
        r = extraire_pdf(c)
        print(c, f"{r['pages']} pages, {len(r['flux3d'])} flux 3D")
        for f in r['flux3d']:
            print(f"  {f['format']} {f['octets']} octets, page {f['page']}, {f['source']}, {len(f['vues'])} vues")
        for e in r['erreurs']:
            print('  !', e)
