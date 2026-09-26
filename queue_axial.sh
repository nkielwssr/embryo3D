#!/bin/bash
# file : pour chaque stade -> (corrections) -> axial -> fusion systèmes
cd "C:/Users/MicroTurtle/Documents/Claude"
run() { echo "== $* $(date +%H:%M:%S)"; python -u "$@" 2>&1 | grep -vE "faces|^\s*$"; }
for st in "$@"; do
  case $st in
    CS16_f4v) run embryo3d/fix_lens_eyes.py CS16_f4v; run embryo3d/axial.py CS16_f4v --seuil=70 ;;
    CS13.f4v) run embryo3d/drop_label.py CS13.f4v cristallins; run embryo3d/axial.py CS13.f4v --seuil=70 ;;
    CS14_f4v|CS15_f4v) run embryo3d/axial.py $st --seuil=70 ;;
    *) run embryo3d/axial.py $st ;;
  esac
  run embryo3d/fusion_systemes.py $st --sans-scene
  echo "OK $st"
done
echo QUEUE_DONE
