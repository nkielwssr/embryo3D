#!/bin/bash
# attend que les jobs python lourds soient finis (<= N restants), puis exécute les commandes séquentiellement
cd "C:/Users/MicroTurtle/Documents/Claude"; N=${N:-1}
until [ "$(tasklist | grep -c python3.11)" -le "$N" ]; do sleep 20; done
while IFS= read -r cmd; do [ -z "$cmd" ] && continue; echo "== $cmd  $(date +%H:%M:%S)"; PYTHONIOENCODING=utf-8 eval "$cmd" 2>&1 | grep -E "^CS|corps |moelle_rach|tip|ventricule_pro|ventricule_rho|Traceback|Error|done"; done
echo QUEUE_SEQ_DONE
