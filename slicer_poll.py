# Boucle de commande : Slicer exécute embryo3d/slicer_cmd.py dès qu'il apparaît, écrit la sortie dans slicer_cmd_out.txt puis supprime le fichier.
import qt, os, io, sys, traceback, slicer
_CMD = r"C:\Users\MicroTurtle\Documents\Claude\embryo3d\slicer_cmd.py"
_OUT = r"C:\Users\MicroTurtle\Documents\Claude\embryo3d\slicer_cmd_out.txt"
def _poll():
    if os.path.exists(_CMD):
        code = open(_CMD, encoding='utf-8').read(); os.remove(_CMD)
        buf = io.StringIO(); old = sys.stdout; sys.stdout = buf
        try: exec(code, globals())
        except Exception: buf.write(traceback.format_exc())
        finally: sys.stdout = old
        open(_OUT, 'w', encoding='utf-8').write(buf.getvalue() + '\n__FIN__\n')
    qt.QTimer.singleShot(1000, _poll)
qt.QTimer.singleShot(1000, _poll); print('poll actif')
