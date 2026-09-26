"""planche de contrôle : 3 coupes médianes du volume fusionné (dens.npy) + histogramme; usage: python check_volume.py <work_dir> [out.png]"""
import numpy as np, cv2, sys, os
work = sys.argv[1]; out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(work, 'check_fused.png')
d = np.load(os.path.join(work, 'dens.npy')); a,b,c = d.shape
sl = [255-d[a//2], 255-d[:, b//2, :], 255-d[:, :, c//2]]
h = max(s.shape[0] for s in sl)
sl = [cv2.copyMakeBorder(s, 0, h-s.shape[0], 0, 4, cv2.BORDER_CONSTANT, value=128) for s in sl]
img = np.hstack(sl); cv2.putText(img, os.path.basename(os.path.dirname(work)) + ' ' + str(d.shape), (10, 30), cv2.FONT_HERSHEY_SIMPLEX, 0.9, 0, 2)
cv2.imwrite(out, img); print(out, img.shape)
