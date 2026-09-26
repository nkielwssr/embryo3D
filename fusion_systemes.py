"""Fusionne dans labels.npz les structures produites par la session « cardio / digestif » (work/cardio/coeur.npz, work/digestif/digestif.npz)
sous des noms préfixés, puis régénère maillages, NRRD, scène et rendu du stade. Relançable (remplace les clés fusionnées).
usage : python fusion_systemes.py <dossier_stade> [...]   (option --sans-scene pour ne pas relancer Blender)"""
import numpy as np, os, sys, re, shutil, subprocess

class _Npz(dict):
    """contenu d'un .npz chargé en mémoire et FERMÉ (np.load garde le fichier ouvert, ce qui bloque os.replace sous Windows)"""
    files = property(lambda self: list(self.keys()))
def _load_npz(path):
    import numpy as _np
    with _np.load(path) as f: return _Npz({k: f[k] for k in f.files})

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
MAP = {'cardio/coeur.npz': {'coeur_plein': 'coeur_detoure', 'cavites_cardiaques': 'cavites_cardiaques', 'myocarde': 'myocarde'},
       'digestif/digestif.npz': {'pharynx': 'digestif_pharynx', 'oesophage': 'digestif_oesophage', 'estomac': 'digestif_estomac', 'duodenum': 'digestif_duodenum',
                                 'intestin_moyen': 'digestif_intestin_moyen', 'intestin_posterieur': 'digestif_intestin_posterieur'}}
def run(stage_dir, scene=True):
    work = os.path.join(stage_dir, 'work'); out = os.path.join(stage_dir, 'out'); stage = 'CS%s' % re.search(r'CS\s*(\d+)', os.path.basename(stage_dir), re.I).group(1)
    z = _load_npz(os.path.join(work, 'labels.npz')); shape = tuple(z['shape']); n = int(np.prod(shape))
    labels = {k: z[k] for k in z.files if k != 'shape'}          # packbits tels quels (pas de dépaquetage inutile)
    added = []
    for f, m in MAP.items():
        p = os.path.join(work, f)
        if not os.path.exists(p): continue
        s = np.load(p)
        if 'shape' in s.files and tuple(s['shape']) != shape: print('  %s : shape %s != %s, ignoré' % (f, tuple(s['shape']), shape)); continue
        for src, dst in m.items():
            if src in s.files:
                pb = s[src]
                if pb.dtype != np.uint8 or pb.size != (n + 7)//8: print('  %s/%s : format inattendu (%s, %d), ignoré' % (f, src, pb.dtype, pb.size)); continue
                labels[dst] = pb; added.append('%s(%d vox)' % (dst, int(np.unpackbits(pb)[:n].sum())))
    # vaisseaux : regroupement par famille (union des segments)
    pv = os.path.join(work, 'cardio', 'vaisseaux.npz')
    if os.path.exists(pv):
        s = np.load(pv)
        if 'shape' not in s.files or tuple(s['shape']) == shape:
            fam = {'vaisseaux_aorte': lambda k: k.startswith('aorte'), 'vaisseaux_cardinales': lambda k: k.startswith('cardinale'),
                   'vaisseaux_ombilicaux': lambda k: k.startswith(('ombilical', 'ombilicale')), 'vaisseaux_vitellins': lambda k: k.startswith(('vitellin', 'vitelline'))}
            for dst, test in fam.items():
                keys = [k for k in s.files if k not in ('shape', 'union') and test(k)]
                if not keys: continue
                u = np.zeros((n + 7)//8, np.uint8)
                for k in keys: u |= s[k]
                labels[dst] = u; added.append('%s(%s, %d vox)' % (dst, '+'.join(keys), int(np.unpackbits(u)[:n].sum())))
    print(stage, 'fusionné :', ', '.join(added) if added else 'rien')
    shutil.copy(os.path.join(work, 'labels.npz'), os.path.join(work, 'labels_avant_fusion.npz'))
    tmp = os.path.join(work, 'labels.tmp.npz'); np.savez_compressed(tmp, **labels, shape=np.array(shape)); os.replace(tmp, os.path.join(work, 'labels.npz'))   # atomique : d'autres sessions lisent labels.npz
    import meshexport, export_nrrd; meshexport.run(work, out, stage); export_nrrd.run(stage_dir)
    if scene:
        B = r'C:/Program Files/Blender Foundation/Blender 5.2/blender.exe'; blend = os.path.join(out, f'{stage}_embryon.blend')
        subprocess.run([B, '-b', '-P', os.path.join(HERE, 'blender_build_scene.py'), '--', os.path.join(out, 'manifest.json'), blend], capture_output=True)
        subprocess.run([B, '-b', blend, '-P', os.path.join(HERE, 'blender_render_scene.py'), '--', os.path.join(out, f'{stage}_render_organs.png'), '90', '0'], capture_output=True)
if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    for d in args: run(os.path.abspath(d), scene='--sans-scene' not in sys.argv)
