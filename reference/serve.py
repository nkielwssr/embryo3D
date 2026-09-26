"""
serve.py — petit serveur local pour viewer.html (les navigateurs refusent fetch() en file://).
Sert le dossier Documents/Claude (racine des CSxx_f4v et d'embryo3d) sur http://127.0.0.1:8769/
et expose /api/models : liste des GLB/PLY disponibles (nos sorties <dossier>/out/*.glb et les hulls).

Usage : python serve.py   puis ouvrir http://127.0.0.1:8769/embryo3d/reference/viewer.html
"""
import http.server, json, os, sys, glob, webbrowser

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))


def hierarchie():
    """Même classement que l'agrégateur (embryo3d/agregateur.py : SYSTEMES, STRUCT_SYS), lu à chaque requête."""
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location('agregateur', os.path.join(HERE, '..', 'agregateur.py'))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return dict(systemes=list(mod.SYSTEMES), struct_sys=dict(mod.STRUCT_SYS))
    except Exception as e:
        return dict(systemes=['Autres'], struct_sys={}, erreur=str(e))
PORT = int(sys.argv[1]) if len(sys.argv) > 1 else 8769


class H(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *a, **k):
        super().__init__(*a, directory=ROOT, **k)

    def do_GET(self):
        if self.path.startswith('/api/models'):
            items = []
            for p in sorted(glob.glob(os.path.join(ROOT, '*', 'out', '*.glb'))):
                items.append(dict(type='glb', path='/' + os.path.relpath(p, ROOT).replace('\\', '/'),
                                  label=os.path.basename(os.path.dirname(os.path.dirname(p))) + ' / ' + os.path.basename(p)))
            for p in sorted(glob.glob(os.path.join(HERE, 'captures', '*', '*_hull.ply'))):
                items.append(dict(type='ply', path='/' + os.path.relpath(p, ROOT).replace('\\', '/'),
                                  label='VOKA hull ' + os.path.basename(p)))
            caps = {}
            for d in sorted(glob.glob(os.path.join(HERE, 'captures', 'J*'))):
                sid = os.path.basename(d)
                caps[sid] = ['/' + os.path.relpath(p, ROOT).replace('\\', '/') for p in sorted(glob.glob(os.path.join(d, '*.png'))) if 'planche' not in p]
            body = json.dumps(dict(models=items, captures=caps, **hierarchie())).encode()
            self.send_response(200); self.send_header('Content-Type', 'application/json'); self.send_header('Content-Length', str(len(body)))
            self.end_headers(); self.wfile.write(body); return
        return super().do_GET()

    def log_message(self, *a):
        pass


if __name__ == '__main__':
    url = f'http://127.0.0.1:{PORT}/embryo3d/reference/viewer.html'
    print('racine servie :', ROOT); print('viewer :', url)
    if '--open' in sys.argv:
        webbrowser.open(url)
    http.server.ThreadingHTTPServer(('127.0.0.1', PORT), H).serve_forever()
