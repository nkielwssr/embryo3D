#!/bin/bash
# Reconstruit les scènes des stades passés en argument (optionnel) puis la scène maître + rendus morph. usage: bash build_master.sh [CS13.f4v ...]
cd "C:/Users/MicroTurtle/Documents/Claude"
B="C:/Program Files/Blender Foundation/Blender 5.2/blender.exe"; OUT=embryons_3D
for d in "$@"; do
  st=$(python -c "import re;print(re.search(r'CS\d+','$d').group(0))")
  "$B" -b -P embryo3d/blender_build_scene.py -- "$PWD/$d/out/manifest.json" "$PWD/$d/out/${st}_embryon.blend" 2>&1 | grep -E "SAVED|Traceback"
  "$B" -b "$PWD/$d/out/${st}_embryon.blend" -P embryo3d/blender_render_scene.py -- "$PWD/$OUT/rendus/${st}_organes.png" 90 0 2>&1 | grep -E "Traceback"
done
python embryo3d/make_master_manifest.py "$OUT/master_manifest.json" CS13.f4v/out CS14_f4v/out CS15_f4v/out CS16_f4v/out CS17_f4v/out CS19_f4v/out CS20_F4V/out
"$B" -b -P embryo3d/blender_build_scene.py -- "$PWD/$OUT/master_manifest.json" "$PWD/$OUT/embryons_CS13-CS20_morph.blend" 2>&1 | grep -E "SAVED|Traceback"
"$B" -b "$PWD/$OUT/embryons_CS13-CS20_morph.blend" -P embryo3d/blender_render_morph.py -- "$PWD/$OUT/rendus/morph" 0 1 2 3 4 5 6 2>&1 | grep -E "Traceback"
echo BUILD_MASTER_DONE
