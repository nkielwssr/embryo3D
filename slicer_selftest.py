# test automatique : Slicer.exe --python-script slicer_selftest.py  (charge CS14, sauvegarde, quitte)
import slicer, os, traceback
try:
    exec(open(r"C:\Users\MicroTurtle\Documents\Claude\embryo3d\slicer_tools.py", encoding="utf-8").read())
    seg = ouvrir(r"C:\Users\MicroTurtle\Documents\Claude\CS14_f4v\out")
    sauver()
    open(r"C:\Users\MicroTurtle\Documents\Claude\embryo3d\slicer_selftest_result.txt", "w").write("OK %d segments\n" % seg.GetSegmentation().GetNumberOfSegments())
except Exception:
    open(r"C:\Users\MicroTurtle\Documents\Claude\embryo3d\slicer_selftest_result.txt", "w").write("ERREUR\n" + traceback.format_exc())
slicer.util.exit(0)
