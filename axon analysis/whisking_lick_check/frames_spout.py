"""Montage of frames around licks with the whisking ME boxes drawn.  usage: python frames_spout.py <session>"""
import sys, os, json, glob, numpy as np, cv2, scipy.io as sio, h5py

sess = sys.argv[1]
here = os.path.dirname(os.path.abspath(__file__))
vid = glob.glob(r'C:\Users\rbondarenko\DLC\work\%s\camera\*.mp4' % sess)[0]
J = json.load(open(r'C:\Users\rbondarenko\DLC\work\%s\whisker_results\%s_whisker.json' % (sess, sess)))
ax_, ay_ = J['anchor_median_xy_px']
boxes = J['boxes_rel_to_nose_px']
cols = {'fan_front': (0, 255, 0), 'fan_down': (255, 0, 255), 'pad': (0, 165, 255), 'ctrl': (255, 255, 0)}
z = np.load(os.path.join(here, 'lt_%s.npz' % sess))
i0 = z['i0']
rng = np.random.default_rng(3)
pick = i0[rng.choice(len(i0), 4, replace=False)]
cap = cv2.VideoCapture(vid)
tiles = []
for lab, frames in (('lick+0.15s', pick + 8), ('lick+0.35s', pick + 18), ('lick-1.0s', pick - 50)):
    for fr in frames:
        cap.set(cv2.CAP_PROP_POS_FRAMES, int(fr))
        ok, f = cap.read()
        if not ok:
            f = np.zeros((1020, 1856, 3), np.uint8)
        for k, (x0, y0, x1, y1) in boxes.items():
            cv2.rectangle(f, (ax_ + x0, ay_ + y0), (ax_ + x1, ay_ + y1), cols[k], 3)
        cv2.circle(f, (ax_, ay_), 6, (0, 0, 255), -1)
        c = f[max(ay_ - 380, 0): ay_ + 380, max(ax_ - 420, 0): ax_ + 560]
        c = cv2.resize(c, (420, 330))
        cv2.putText(c, '%s #%d' % (lab, fr), (5, 20), cv2.FONT_HERSHEY_SIMPLEX, 0.55, (0, 255, 255), 2)
        tiles.append(c)
rows = [np.hstack(tiles[i:i + 4]) for i in range(0, len(tiles), 4)]
cv2.imwrite(os.path.join(here, 'frames_%s.png' % sess), np.vstack(rows))
print('anchor', ax_, ay_)
