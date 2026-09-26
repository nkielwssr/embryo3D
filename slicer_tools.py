# -*- coding: utf-8 -*-
"""Outils pour 3D Slicer (à exécuter DANS la console Python de Slicer : Affichage > Python console).
    exec(open(r"C:\\Users\\MicroTurtle\\Documents\\Claude\\embryo3d\\slicer_tools.py", encoding="utf-8").read())
    ouvrir(r"C:\\Users\\MicroTurtle\\Documents\\Claude\\CS16_f4v\\out")     # charge densité + structures nommées, ouvre Segment Editor
    ...  corrections avec Paint / Erase / Scissors / Islands ...
    sauver()                                                           # écrit <CS>_labels_corrige.nrrd + <CS>_labels_corrige.json dans out/
Puis, dans un terminal :  python embryo3d/import_labels.py CS16_f4v
"""
import slicer, json, os, glob, qt
COULEURS = {
 'snc': (1.00, 0.72, 0.30), 'ventricules': (0.55, 0.85, 1.00), 'ganglions': (1.00, 0.55, 0.55), 'yeux': (0.25, 0.35, 0.90),
 'cristallins': (0.95, 0.95, 0.60), 'vesicules_otiques': (0.30, 0.90, 0.90), 'coeur': (0.80, 0.15, 0.15), 'cavite_pericardique': (0.75, 0.85, 1.00),
 'foie': (0.55, 0.20, 0.30), 'tube_digestif': (0.95, 0.60, 0.20), 'squelette_axial_cartilage': (0.92, 0.92, 0.85), 'chondrocrane': (0.90, 0.90, 0.80),
 'cartilage_autre': (0.85, 0.85, 0.80), 'membres': (0.90, 0.75, 0.65), 'cordon_ombilical': (0.80, 0.70, 0.60), 'vaisseaux': (0.85, 0.10, 0.10)}
ETAT = {}

def ouvrir(out_dir):
    out_dir = str(out_dir)
    dens = glob.glob(os.path.join(out_dir, 'CS*_density.nrrd'))[0]
    cs = os.path.basename(dens).split('_')[0]
    vol = slicer.util.loadVolume(dens)
    vol.SetName(cs + '_densite')
    lab = slicer.util.loadLabelVolume(os.path.join(out_dir, cs + '_labels.nrrd'))
    ids = json.load(open(os.path.join(out_dir, cs + '_label_ids.json')))
    seg = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentationNode', cs + '_structures')
    seg.SetReferenceImageGeometryParameterFromVolumeNode(vol)
    slicer.modules.segmentations.logic().ImportLabelmapToSegmentationNode(lab, seg)
    slicer.mrmlScene.RemoveNode(lab)
    sg = seg.GetSegmentation()
    for i in range(sg.GetNumberOfSegments()):
        s = sg.GetNthSegment(i)
        val = int(s.GetLabelValue())
        for name, v in ids.items():
            if v == val:
                s.SetName(name); s.SetColor(*COULEURS.get(name, (0.7, 0.7, 0.7)))
    seg.CreateClosedSurfaceRepresentation()
    # affichage : densité en niveaux de gris inversés (tissu sombre comme sur les coupes)
    disp = vol.GetDisplayNode(); disp.SetAutoWindowLevel(False); disp.SetWindowLevel(255, 127)
    slicer.util.setSliceViewerLayers(background=vol)
    slicer.util.selectModule('SegmentEditor')
    w = slicer.modules.segmenteditor.widgetRepresentation().self().editor
    w.setSegmentationNode(seg); w.setSourceVolumeNode(vol) if hasattr(w, 'setSourceVolumeNode') else w.setMasterVolumeNode(vol)
    slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutFourUpView)
    slicer.util.resetSliceViews()
    ETAT.update(cs=cs, out=out_dir, seg=seg, vol=vol, ids=ids)
    print('Chargé', cs, ':', sg.GetNumberOfSegments(), 'structures. Corrigez dans Segment Editor puis appelez sauver().')
    return seg

def sauver():
    seg, vol, cs, out = ETAT['seg'], ETAT['vol'], ETAT['cs'], ETAT['out']
    lm = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLabelMapVolumeNode', cs + '_labels_corrige')
    ok = slicer.modules.segmentations.logic().ExportAllSegmentsToLabelmapNode(seg, lm, slicer.vtkSegmentation.EXTENT_REFERENCE_GEOMETRY)
    if not ok: raise RuntimeError('export labelmap impossible')
    # valeur réellement écrite pour chaque segment (Slicer renumérote à l'export) : on la mesure dans le labelmap exporté
    import numpy as np
    exported = slicer.util.arrayFromVolume(lm)
    mapping = {}
    sg = seg.GetSegmentation()
    for i in range(sg.GetNumberOfSegments()):
        s = sg.GetNthSegment(i); sid = sg.GetNthSegmentID(i)
        try:
            m = slicer.util.arrayFromSegmentBinaryLabelmap(seg, sid, lm) > 0
            vals, counts = np.unique(exported[m], return_counts=True)
            vals = vals[vals > 0]; counts = counts[-len(vals):] if len(vals) else counts
            mapping[s.GetName()] = int(vals[np.argmax(counts)]) if len(vals) else 0     # 0 = segment vide, ignoré à l'import
        except Exception as e:
            mapping[s.GetName()] = int(s.GetLabelValue()); print('  (valeur par défaut pour', s.GetName(), ':', e, ')')
    path = os.path.join(out, cs + '_labels_corrige.nrrd')
    slicer.util.saveNode(lm, path)
    json.dump(mapping, open(os.path.join(out, cs + '_labels_corrige.json'), 'w'), indent=1)
    slicer.mrmlScene.RemoveNode(lm)
    print('Sauvé :', path, '\n  puis lancer :  python embryo3d/import_labels.py', os.path.basename(os.path.dirname(out)))
