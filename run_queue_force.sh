#!/bin/bash
STEPS="$1"; shift
cd "C:/Users/MicroTurtle/Documents/Claude"
for d in "$@"; do
  echo "===== $d  $(date +%H:%M:%S)"
  python -u embryo3d/pipeline.py "$d" --steps "$STEPS" --force >> "$d/pipeline_log.txt" 2>&1 && echo "OK $d" || echo "FAILED $d"
done
echo QUEUE_DONE
