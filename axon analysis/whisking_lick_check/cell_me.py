"""Motion energy of a fine grid of cells around the nose, for every frame (same ME definition as DLC/face_model/face_motion.py: 1 px blur, |frame t - frame t-1|,
both frames cropped at the same absolute, nose-anchored position; keyframe + following frame interpolated).
usage: python cell_me.py <session>     -> cells_<session>.npy  float32 [frames, NY, NX],  cells_<session>.json (grid definition)
Grid: x -320..320, y -160..320 relative to the nose anchor, 20 px cells (NX = 32, NY = 24)."""
import sys, os, glob, json, subprocess, time
import numpy as np, pandas as pd, cv2

sys.path.insert(0, r'C:\Users\rbondarenko\projects\2P_Scripts\DLC\face_model')
import face_motion as fm

sess = sys.argv[1]
here = os.path.dirname(os.path.abspath(__file__))
work = r'C:\Users\rbondarenko\DLC\work\%s' % sess
video = glob.glob(os.path.join(work, 'camera', '*.mp4'))[0]
face = pd.read_csv(glob.glob(os.path.join(work, 'face_samples', '*_face.csv'))[0], index_col=0)
n_frames, width, height = fm.video_frames(video)
ax, ay = fm.nose_anchor(face, n_frames)
print(sess, 'frames', n_frames, 'anchor median', int(np.median(ax)), int(np.median(ay)), 'range', int(ax.max() - ax.min()), int(ay.max() - ay.min()), flush=True)
CS = 20
GX0, GX1, GY0, GY1 = -320, 320, -160, 320
NX, NY = (GX1 - GX0) // CS, (GY1 - GY0) // CS
x0 = max(int(ax.min() + GX0), 0)
x1 = min(int(ax.max() + GX1), width)
y0 = max(int(ay.min() + GY0), 0)
y1 = min(int(ay.max() + GY1), height)
w, h = x1 - x0, y1 - y0
cmd = [fm.find_tool('ffmpeg'), '-v', 'error', '-i', video, '-vf', 'crop=%d:%d:%d:%d:exact=1,format=gray' % (w, h, x0, y0), '-f', 'rawvideo', 'pipe:1']
proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=w * h * 32)
out = np.lib.format.open_memmap(os.path.join(here, 'cells_%s.npy' % sess), mode='w+', dtype=np.float32, shape=(n_frames, NY, NX))
out[:] = np.nan
prev = None
i = 0
t0 = time.time()
while True:
    buf = proc.stdout.read(w * h)
    if len(buf) < w * h:
        break
    g = cv2.GaussianBlur(np.frombuffer(buf, np.uint8).reshape(h, w), (0, 0), fm.BLUR_SIGMA_PX).astype(np.float32)
    if prev is not None and i < n_frames:
        a = min(i, len(ax) - 1)
        c0, r0 = int(ax[a] + GX0 - x0), int(ay[a] + GY0 - y0)
        if c0 >= 0 and r0 >= 0 and c0 + NX * CS <= w and r0 + NY * CS <= h:
            d = np.abs(g[r0:r0 + NY * CS, c0:c0 + NX * CS] - prev[r0:r0 + NY * CS, c0:c0 + NX * CS])
            out[i] = cv2.resize(d, (NX, NY), interpolation=cv2.INTER_AREA)
    prev = g
    i += 1
    if i % 20000 == 0:
        print('frame %d / %d  %.0f s' % (i, n_frames, time.time() - t0), flush=True)
proc.wait()
print('decoded', i, 'frames in %.0f s' % (time.time() - t0), flush=True)
# codec noise: keyframe + following frame interpolated (same rule as face_motion.drop_codec_frames)
is_key = fm.keyframes(video)
n = n_frames
key = is_key[:n]
bad = key.copy() if key.mean() <= 0.5 else np.zeros(n, bool)
if key.mean() <= 0.5:
    bad[1:] |= key[:-1]
bad[0] = True
good = np.flatnonzero(~bad)
print('keyframes %.4f of frames; bad frames %d' % (key.mean(), bad.sum()), flush=True)
badi = np.flatnonzero(bad)
hi = np.searchsorted(good, badi)
lo = np.clip(hi - 1, 0, len(good) - 1)
hi = np.clip(hi, 0, len(good) - 1)
gl, gh = good[lo], good[hi]
wgt = np.where(gh > gl, (badi - gl) / np.maximum(gh - gl, 1), 0.0).astype(np.float32)
for k0 in range(0, len(badi), 20000):
    s = slice(k0, k0 + 20000)
    out[badi[s]] = out[gl[s]] * (1 - wgt[s, None, None]) + out[gh[s]] * wgt[s, None, None]
out.flush()
json.dump(dict(session=sess, frames=n_frames, CS=CS, GX=[GX0, GX1], GY=[GY0, GY1], NX=NX, NY=NY, anchor_median=[int(np.median(ax)), int(np.median(ay))]),
          open(os.path.join(here, 'cells_%s.json' % sess), 'w'))
print('DONE', flush=True)
