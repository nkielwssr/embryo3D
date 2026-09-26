"""Encode une séquence PNG (dossier) en mp4 : python encode_frames.py dossier sortie.mp4 [fps]  (ffmpeg si présent, sinon OpenCV)"""
import sys, os, glob, shutil, subprocess
d, out = sys.argv[1], sys.argv[2]; fps = int(sys.argv[3]) if len(sys.argv) > 3 else 24
files = sorted(glob.glob(os.path.join(d, '*.png')))
if shutil.which('ffmpeg'):
    lst = os.path.join(d, 'liste.txt')
    with open(lst, 'w') as f:
        for fn in files: f.write("file '%s'\n" % os.path.basename(fn))
    subprocess.run(['ffmpeg', '-y', '-r', str(fps), '-f', 'concat', '-safe', '0', '-i', lst, '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-crf', '20', out], check=True, capture_output=True)
else:
    import cv2
    im = cv2.imread(files[0]); h, w = im.shape[:2]
    vw = cv2.VideoWriter(out, cv2.VideoWriter_fourcc(*'mp4v'), fps, (w, h))
    for fn in files: vw.write(cv2.imread(fn))
    vw.release()
print('ENCODE_DONE', out, len(files), 'images')
