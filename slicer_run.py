"""Envoie du code Python à l'instance Slicer en cours (boucle slicer_poll.py) et affiche la sortie. usage: python slicer_run.py "code" [timeout_s]"""
import sys, os, time
HERE = os.path.dirname(os.path.abspath(__file__)); CMD = os.path.join(HERE, 'slicer_cmd.py'); OUT = os.path.join(HERE, 'slicer_cmd_out.txt')
code = sys.argv[1]; timeout = float(sys.argv[2]) if len(sys.argv) > 2 else 600
if os.path.exists(OUT): os.remove(OUT)
open(CMD, 'w', encoding='utf-8').write(code)
t0 = time.time()
while time.time() - t0 < timeout:
    if os.path.exists(OUT):
        s = open(OUT, encoding='utf-8').read()
        if '__FIN__' in s: print(s.replace('__FIN__', '').strip()); sys.exit(0)
    time.sleep(0.5)
print('TIMEOUT'); sys.exit(1)
