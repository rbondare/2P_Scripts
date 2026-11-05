"""Synthetic 5-plane (+1 flyback) bouton recording with known axial motion and known activity.
usage: python make_synth.py <out_dir> [signal_scale]   (1 = dim; see eval_synth.py for the full test)

Layout per volume in the tif: planes 0..4 then one flyback frame (noise only), like ScanImage fastZ.
Tifs: 6 acquisitions of 250 volumes; a 1-frame stub tif sits between acquisitions 3 and 4, and acquisition 5 ends
2 frames into a new volume (partial volume) - both must be handled without shifting the following tifs.
"""
import sys
from pathlib import Path

import numpy as np
import tifffile

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(1)

NZ, NY, NX = 5, 256, 256
DZ, PX, FS = 3.0, 0.5, 5.0
NT = 1500
N_B = 60
SIG_XY, SIG_Z = 1.5, 2.5          # px, um
zk = (np.arange(NZ) - 2) * DZ      # plane depths (um)

# bouton positions
by = rng.uniform(20, NY - 20, N_B)
bx = rng.uniform(20, NX - 20, N_B)
bz = rng.uniform(-4.5, 4.5, N_B)
AMP_SCALE = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
amp = rng.uniform(10, 24, N_B) * AMP_SCALE     # counts above background at focus

# calcium activity: sparse events, exp decay tau = 1 s
dff = np.zeros((N_B, NT))
decay = np.exp(-1.0 / (FS * 1.0))
for i in range(N_B):
    ev = rng.random(NT) < 0.012
    a = ev * rng.uniform(1.0, 3.0, NT)
    s = 0.0
    for t in range(NT):
        s = s * decay + a[t]
        dff[i, t] = s

# axial motion (um): running-like excursions with relaxation + slow drift + jitter
z = np.zeros(NT)
state = 0.0
onsets = np.sort(rng.choice(np.arange(60, NT - 60), 14, replace=False))
target = np.zeros(NT)
for o in onsets:
    target[o:o + rng.integers(20, 70)] = rng.choice([-1, 1]) * rng.uniform(1.0, 3.2)
for t in range(NT):
    state += (target[t] - state) * 0.35
    z[t] = state
z += np.linspace(0, 0.9, NT) + rng.normal(0, 0.15, NT)
# lateral jitter (integer px)
sx = np.clip(np.cumsum(rng.integers(-1, 2, NT)), -4, 4)
sy = np.clip(np.cumsum(rng.integers(-1, 2, NT)), -4, 4)

# spatial templates
yy, xx = np.mgrid[0:NY, 0:NX]
S = np.exp(-((yy[None] - by[:, None, None]) ** 2 + (xx[None] - bx[:, None, None]) ** 2) / (2 * SIG_XY ** 2))
S = S.reshape(N_B, -1).astype(np.float32)

# weights per (t, plane, bouton)
zeff = bz[None, None, :] + z[:, None, None]                         # tissue moves by z
w = amp[None, None, :] * (1 + dff.T[:, None, :]) * np.exp(-(zk[None, :, None] - zeff) ** 2 / (2 * SIG_Z ** 2))
signal = (w.reshape(NT * NZ, N_B).astype(np.float32) @ S).reshape(NT, NZ, NY, NX)

haze = 4 * np.exp(-((yy - 128) ** 2 + (xx - 100) ** 2) / (2 * 60 ** 2)).astype(np.float32)
BG = 20.0
vol = np.empty((NT, NZ + 1, NY, NX), dtype=np.int16)
for t in range(NT):
    mean = BG + haze[None] * (1 + 0.1 * np.sin(t / 50)) + signal[t]
    frame = rng.poisson(mean).astype(np.int16)
    frame = np.roll(frame, (int(sy[t]), int(sx[t])), axis=(1, 2))
    vol[t, :NZ] = frame
    vol[t, NZ] = rng.poisson(BG, (NY, NX))
print("movie", vol.shape, "mean counts", vol[:, :NZ].mean())

np.savez(out / "truth.npz", z_um=z, dff=dff, by=by, bx=bx, bz=bz, amp=amp, sy=sy, sx=sx, fs=FS, dz=DZ)

per = 250
frames = vol.reshape(NT * (NZ + 1), NY, NX)
names = []
for k in range(6):
    chunk = frames[k * per * (NZ + 1):(k + 1) * per * (NZ + 1)]
    if k == 4:   # acquisition ends 2 frames into a new volume
        chunk = np.concatenate([chunk, rng.poisson(BG, (2, NY, NX)).astype(np.int16)])
    name = f"synth_{k + 1:05d}.tif" if k < 3 else f"synth_{k + 2:05d}.tif"
    tifffile.imwrite(out / name, chunk)
    names.append((name, len(chunk)))
tifffile.imwrite(out / "synth_00004.tif", rng.poisson(BG, (1, NY, NX)).astype(np.int16))   # stub acquisition
print(names, "+ synth_00004.tif (1 frame)")
