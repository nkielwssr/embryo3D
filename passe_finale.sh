#!/bin/bash
# Passe complète : 7 recrutements imposés (séquentiels) -> scène maître (parent normalisant, export morph_web) -> planche (avec enveloppe) -> vidéo 1080p.
cd "C:/Users/MicroTurtle/Documents/Claude"
N=1 bash embryo3d/queue_seq.sh < queue_imposes.txt
bash embryo3d/build_master.sh CS13.f4v CS14_f4v CS15_f4v CS16_f4v CS17_f4v CS19_f4v CS20_F4V
python embryo3d/planche_morph.py embryons_3D/rendus/planche_morph.png --enveloppe 0 0.5 1 2 3 3.5 4 5 5.5 6
"C:/Program Files/Blender Foundation/Blender 5.2/blender.exe" -b embryons_3D/embryons_CS13-CS20_morph.blend -P embryo3d/blender_render_anim.py -- "$PWD/embryons_3D/rendus/morph_animation.mp4" 1920 1080 90 2>&1 | grep -E "ANIM_DONE|CAMERA|Traceback"
python embryo3d/encode_frames.py embryons_3D/rendus/morph_animation_frames embryons_3D/rendus/morph_animation.mp4 24
echo PASSE_FINALE_DONE
