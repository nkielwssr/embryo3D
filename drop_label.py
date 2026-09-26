"""Vide une structure de labels.npz (sans régénérer). usage: python drop_label.py <dossier_stade> <nom>"""
import numpy as np, os, sys

class _Npz(dict):
    """contenu d'un .npz chargé en mémoire et FERMÉ (np.load garde le fichier ouvert, ce qui bloque os.replace sous Windows)"""
    files = property(lambda self: list(self.keys()))
def _load_npz(path):
    import numpy as _np
    with _np.load(path) as f: return _Npz({k: f[k] for k in f.files})


def _savez_atomic(path, **kw):
    """écriture atomique (fichier temporaire puis remplacement) : d'autres sessions lisent labels.npz en parallèle"""
    import numpy as _np, os as _os; tmp = path + '.tmp.npz'; _np.savez_compressed(tmp, **kw); _os.replace(tmp, path)

work = os.path.join(os.path.abspath(sys.argv[1]), 'work'); nom = sys.argv[2]
z = _load_npz(os.path.join(work, 'labels.npz')); d = {k: z[k] for k in z.files if k not in ('shape', nom)}
_savez_atomic(os.path.join(work, 'labels.npz'), **d, shape=z['shape']); print('retiré', nom, 'de', work)
