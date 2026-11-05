"""Camera-vs-ball timing check: motion energy of the ball / paw region of the raw video against the hardware ball speed.

    python video_ball_lag_check.py <video.mp4> <export.mat> [--box x y w h] [--out <png>]

export.mat (v5, made in MATLAB from an aggregated file) holds tc (camera trigger times, s), bt (ball trigger times), bs
(ball speed per reading), t1 (2P volume start times, TimeCa(1,:)) and off (2P lateral offset per volume).
The ball and the paws are seen by the camera and the ball by its own sensor; both are caused by the same movement, so the
cross-correlation of the two must peak at ~0 lag when the video frames are paired with the right trigger times.
A peak at +-0.02 s is one camera frame. Run in an environment with ffmpeg on PATH (dlc-tf) and scipy.
"""
import argparse
import subprocess

import numpy as np
import scipy.io as sio
from scipy.signal import correlate, correlation_lags

ap = argparse.ArgumentParser()
ap.add_argument("video")
ap.add_argument("export")
ap.add_argument("--box", type=int, nargs=4, default=[700, 780, 1156, 240], metavar=("X", "Y", "W", "H"))
ap.add_argument("--down", type=int, default=4)
a = ap.parse_args()

x, y, w, h = a.box
ow, oh = w // a.down, h // a.down
cmd = ["ffmpeg", "-v", "error", "-i", a.video, "-vf", f"crop={w}:{h}:{x}:{y},scale={ow}:{oh},format=gray", "-f", "rawvideo", "-"]
p = subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=ow * oh * 64)
me, prev, n = [], None, 0
while True:
    buf = p.stdout.read(ow * oh)
    if len(buf) < ow * oh:
        break
    f = np.frombuffer(buf, dtype=np.uint8).astype(np.float32)
    if prev is not None:
        me.append(np.abs(f - prev).mean())
    prev = f
    n += 1
p.wait()
me = np.array([me[0]] + me)               # one value per frame (the first repeats)
print("frames decoded:", n)

D = sio.loadmat(a.export, squeeze_me=True)
tc, bt, bs = D["tc"].astype(float), D["bt"].astype(float), D["bs"].astype(float)
assert abs(len(tc) - n) <= 1, f"{n} decoded frames vs {len(tc)} camera triggers"
m = min(n, len(tc))
tc, me = tc[:m], me[:m]
# codec: keyframes (every 12th frame) and the next frame carry noise; a 1 s window averages it evenly
dt = 0.02
nb = min(len(bt), len(bs))
bt, bs = bt[:nb], bs[:nb]
grid = np.arange(max(tc[0], bt[0]) + 1, min(tc[-1], bt[-1]) - 1, dt)
ve = np.interp(grid, tc, me)
bsp = np.interp(grid, bt, bs)
k = int(1.0 / dt)
sm = lambda v: np.convolve(v, np.ones(k) / k, mode="same")
v1, b1 = sm(ve) - sm(ve).mean(), sm(bsp) - sm(bsp).mean()
# high-pass 20 s so slow arousal-related drifts do not dominate
hp = lambda v: v - np.convolve(v, np.ones(int(20 / dt)) / int(20 / dt), mode="same")
v1, b1 = hp(v1), hp(b1)
cc = correlate(v1, b1, mode="full") / (np.linalg.norm(v1) * np.linalg.norm(b1))
lags = correlation_lags(len(v1), len(b1), mode="full") * dt
sel = np.abs(lags) <= 5
i = np.argmax(cc[sel])
print(f"video ME (ball/paw box) vs ball speed: peak r = {cc[sel][i]:.3f} at lag {lags[sel][i]:+.2f} s (video later when positive)")
rng = np.random.default_rng(0)
null = []
for _ in range(60):
    s = rng.integers(int(60 / dt), len(b1) - int(60 / dt))
    c = correlate(v1, np.roll(b1, s), mode="full") / (np.linalg.norm(v1) * np.linalg.norm(b1))
    null.append(np.abs(c[sel]).max())
print(f"null (circular shifts) 99th pct of the peak: {np.percentile(null, 99):.3f}")
# peak lag at finer detail: parabolic interpolation
j = np.where(sel)[0][i]
if 0 < j < len(cc) - 1:
    d = 0.5 * (cc[j - 1] - cc[j + 1]) / (cc[j - 1] - 2 * cc[j] + cc[j + 1])
    print(f"refined peak lag {lags[j] + d * dt:+.3f} s")
# lag of the peak in the first and second half separately (does the offset change over the session?)
for name, sl in (("first half", slice(0, len(v1) // 2)), ("second half", slice(len(v1) // 2, None))):
    c = correlate(v1[sl], b1[sl], mode="full") / (np.linalg.norm(v1[sl]) * np.linalg.norm(b1[sl]) + 1e-12)
    l = correlation_lags(len(v1[sl]), len(b1[sl]), mode="full") * dt
    s2 = np.abs(l) <= 5
    print(f"  {name}: peak r {c[s2].max():.3f} at {l[s2][np.argmax(c[s2])]:+.2f} s")
