# -*- coding: utf-8 -*-
"""Fonctions de correction semi-automatiques pour la console Python de 3D Slicer (après ouvrir()).
    exec(open(r"C:\\Users\\MicroTurtle\\Documents\\Claude\\embryo3d\\slicer_correct.py", encoding="utf-8").read())
    vue('sag')                 # grande vue sagittale ('4' = 4 vues, 'cor' frontale, 'ax' transversale, '3d')
    coupe('sag', -1.2)         # position de la coupe en mm
    peindre('foie', (x,y,z), 1.0, seuil=110)     # sphère de rayon 1 mm centrée en RAS, ne garde que la densité >= seuil
    effacer('snc', (x,y,z), 0.8)
    ilots('snc')               # garde le plus gros morceau
    fermer('foie', 3)          # fermeture morphologique (voxels)
    volume_mm3('foie')
"""
import slicer, numpy as np, vtk
try:
    from scipy import ndimage as _ndi
except Exception:
    _ndi = None
LAYOUTS = {'4': 3, 'sag': 7, 'ax': 6, 'cor': 8, '3d': 4}
VIEWS = {'sag': 'Yellow', 'ax': 'Red', 'cor': 'Green'}

def _seg(): return ETAT['seg']
def _vol(): return ETAT['vol']
def vue(k): slicer.app.layoutManager().setLayout(LAYOUTS[k]); slicer.util.resetSliceViews()
def coupe(k, mm):
    slicer.app.layoutManager().sliceWidget(VIEWS[k]).sliceLogic().SetSliceOffset(float(mm))
def _sid(nom, creer=True):
    sg = _seg().GetSegmentation()
    for i in range(sg.GetNumberOfSegments()):
        if sg.GetNthSegment(i).GetName() == nom: return sg.GetNthSegmentID(i)
    if not creer: return None
    sid = sg.AddEmptySegment('', nom); sg.GetSegment(sid).SetColor(*COULEURS.get(nom, (0.7, 0.7, 0.7))); return sid
def _arr(nom):
    sid = _sid(nom)
    a = slicer.util.arrayFromSegmentBinaryLabelmap(_seg(), sid, _vol())
    return (a > 0) if a is not None else np.zeros(slicer.util.arrayFromVolume(_vol()).shape, bool), sid
def _set(nom, m):
    sid = _sid(nom)
    slicer.util.updateSegmentBinaryLabelmapFromArray(m.astype(np.uint8), _seg(), sid, _vol())
    _seg().Modified()
def _sphere(centre_ras, rayon_mm):
    v = _vol(); ras2ijk = vtk.vtkMatrix4x4(); v.GetRASToIJKMatrix(ras2ijk)
    ijk = [ras2ijk.MultiplyPoint(list(centre_ras) + [1.0])[i] for i in range(3)]
    sp = v.GetSpacing(); shape = slicer.util.arrayFromVolume(v).shape          # (K, J, I)
    k, j, i = np.ogrid[:shape[0], :shape[1], :shape[2]]
    d2 = ((i - ijk[0]) * sp[0]) ** 2 + ((j - ijk[1]) * sp[1]) ** 2 + ((k - ijk[2]) * sp[2]) ** 2
    return d2 <= rayon_mm ** 2
def peindre(nom, centre_ras, rayon_mm, seuil=None, seuil_max=None):
    m, _ = _arr(nom); s = _sphere(centre_ras, rayon_mm)
    if seuil is not None or seuil_max is not None:
        d = slicer.util.arrayFromVolume(_vol())
        if seuil is not None: s &= d >= seuil
        if seuil_max is not None: s &= d <= seuil_max
    _set(nom, m | s); print(nom, ': +%d voxels -> %.3f mm3' % (s.sum(), volume_mm3(nom)))
def effacer(nom, centre_ras, rayon_mm):
    m, _ = _arr(nom); s = _sphere(centre_ras, rayon_mm); _set(nom, m & ~s); print(nom, ': -%d voxels' % (m & s).sum())
def retirer_de(nom, autre):
    """retire de `nom` tout ce qui appartient à `autre`"""
    m, _ = _arr(nom); o, _ = _arr(autre); _set(nom, m & ~o); print(nom, ': -%d voxels (chevauchement avec %s)' % ((m & o).sum(), autre))
def ilots(nom, min_vox=None):
    m, _ = _arr(nom)
    lab, n = _ndi.label(m); s = np.bincount(lab.ravel())[1:]
    keep = (lab == (np.argmax(s) + 1)) if min_vox is None else np.isin(lab, np.flatnonzero(s >= min_vox) + 1)
    _set(nom, keep); print(nom, ': %d ilots -> %d gardés, %.3f mm3' % (n, 1 if min_vox is None else int((s >= min_vox).sum()), volume_mm3(nom)))
def fermer(nom, r):
    m, _ = _arr(nom); st = _ndi.generate_binary_structure(3, 1)
    mm = _ndi.binary_closing(m, structure=st, iterations=int(r)); mm = _ndi.binary_fill_holes(mm)
    _set(nom, mm); print(nom, ': fermeture r=%d, %.3f mm3' % (r, volume_mm3(nom)))
def ouvrir_morpho(nom, r):
    m, _ = _arr(nom); st = _ndi.generate_binary_structure(3, 1)
    _set(nom, _ndi.binary_opening(m, structure=st, iterations=int(r))); print(nom, ': ouverture r=%d, %.3f mm3' % (r, volume_mm3(nom)))
def volume_mm3(nom):
    m, _ = _arr(nom); sp = _vol().GetSpacing(); return float(m.sum()) * sp[0] * sp[1] * sp[2]
def liste():
    sg = _seg().GetSegmentation()
    for i in range(sg.GetNumberOfSegments()):
        n = sg.GetNthSegment(i).GetName(); print('%-28s %10.3f mm3' % (n, volume_mm3(n)))
def montrer(*noms):
    """n'affiche que ces structures (toutes si vide)"""
    d = _seg().GetDisplayNode(); sg = _seg().GetSegmentation()
    for i in range(sg.GetNumberOfSegments()):
        sid = sg.GetNthSegmentID(i); d.SetSegmentVisibility(sid, (not noms) or sg.GetNthSegment(i).GetName() in noms)
def probe(ras):
    """densité au point RAS"""
    v = _vol(); m = vtk.vtkMatrix4x4(); v.GetRASToIJKMatrix(m); ijk = [int(round(m.MultiplyPoint(list(ras) + [1.0])[i])) for i in range(3)]
    return int(slicer.util.arrayFromVolume(v)[ijk[2], ijk[1], ijk[0]])
print('slicer_correct chargé : vue, coupe, peindre, effacer, retirer_de, ilots, fermer, ouvrir_morpho, volume_mm3, liste, montrer, probe')

def _sphere_ijk(ijk, rayon_vox):
    shape = slicer.util.arrayFromVolume(_vol()).shape          # (K, J, I)
    k, j, i = np.ogrid[:shape[0], :shape[1], :shape[2]]
    return ((i - ijk[0]) ** 2 + (j - ijk[1]) ** 2 + (k - ijk[2]) ** 2) <= rayon_vox ** 2
def peindre_ijk(nom, ijk, rayon_vox, seuil=None, seuil_max=None, sauf=()):
    """comme peindre mais en indices voxel (LR, SI, AP) du pipeline ; `sauf` = structures à ne pas recouvrir"""
    m, _ = _arr(nom); s = _sphere_ijk(ijk, rayon_vox); d = slicer.util.arrayFromVolume(_vol())
    if seuil is not None: s &= d >= seuil
    if seuil_max is not None: s &= d <= seuil_max
    for autre in sauf:
        o, _ = _arr(autre); s &= ~o
    _set(nom, m | s); print(nom, ': +%d voxels -> %.3f mm3' % (s.sum(), volume_mm3(nom)))
def effacer_ijk(nom, ijk, rayon_vox):
    m, _ = _arr(nom); s = _sphere_ijk(ijk, rayon_vox); _set(nom, m & ~s); print(nom, ': -%d voxels' % (m & s).sum())
def coupe_ijk(k, index):
    """positionne la coupe 'sag'/'cor'/'ax' sur un indice voxel (LR / AP / SI)"""
    v = _vol(); m = vtk.vtkMatrix4x4(); v.GetIJKToRASMatrix(m)
    ijk = {'sag': [index, 0, 0], 'ax': [0, index, 0], 'cor': [0, 0, index]}[k]
    ras = m.MultiplyPoint(ijk + [1.0]); off = {'sag': ras[0], 'ax': ras[2], 'cor': ras[1]}[k]
    slicer.app.layoutManager().sliceWidget(VIEWS[k]).sliceLogic().SetSliceOffset(off)
def recharger():
    slicer.mrmlScene.Clear(0); return ouvrir(ETAT['out'])

def fusionner(src, dst):
    """verse la structure src dans dst (src devient vide)"""
    a, _ = _arr(src); b, _ = _arr(dst); _set(dst, a | b); _set(src, np.zeros_like(a)); print('%s -> %s : %d voxels' % (src, dst, a.sum()))
def paroi(nom, cavite, epaisseur_vox=14, seuil=70, marge_ext=6):
    """ajoute à `nom` le tissu dense (>= seuil) situé à moins de `epaisseur_vox` de `cavite`, en évitant la surface du corps (fond à moins de marge_ext)"""
    d = slicer.util.arrayFromVolume(_vol()); cav, _ = _arr(cavite); m, _ = _arr(nom)
    near = _ndi.distance_transform_edt(~cav) <= epaisseur_vox
    fond = (d < 3) & ~cav
    ext = _ndi.distance_transform_edt(~fond) <= marge_ext
    add = near & (d >= seuil) & ~ext
    _set(nom, m | add); print('%s : +%d voxels (paroi de %s) -> %.3f mm3' % (nom, add.sum(), cavite, volume_mm3(nom)))
