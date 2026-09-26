#!/bin/bash
# Refusionne chaque stade dont coeur.npz / vaisseaux.npz / digestif.npz est plus récent que labels.npz, puis reconstruit
# les scènes touchées + la scène maître. Boucle jusqu'à ce qu'il n'y ait plus rien à faire (max 3 tours).
cd "C:/Users/MicroTurtle/Documents/Claude"
for tour in 1 2 3; do
  redo=()
  for d in CS13.f4v CS14_f4v CS15_f4v CS16_f4v CS17_f4v CS19_f4v CS20_F4V; do
    for f in cardio/coeur.npz cardio/vaisseaux.npz digestif/digestif.npz; do
      if [ -f "$d/work/$f" ] && [ "$d/work/$f" -nt "$d/work/labels.npz" ]; then redo+=("$d"); break; fi
    done
  done
  if [ ${#redo[@]} -eq 0 ]; then echo "tour $tour : rien à refaire"; break; fi
  echo "tour $tour : refusion [${redo[*]}] $(date +%H:%M:%S)"
  for d in "${redo[@]}"; do
    python embryo3d/fusion_systemes.py "$d" --sans-scene 2>&1 | grep -E "fusionn|Traceback" | grep -oE "^CS[0-9]+|vaisseaux_[a-z]+\([^)]*\)|Traceback.*" | tr '\n' ' '; echo
  done
  bash embryo3d/build_master.sh "${redo[@]}" 2>&1 | tail -1
done
echo SYNC_DONE
