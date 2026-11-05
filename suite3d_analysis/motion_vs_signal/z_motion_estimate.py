"""Registration-free axial-displacement estimate per volume from the plane profile itself.

    python z_motion_estimate.py <diag_dir> <s3d_job_dir> [--sigma 6] [--margin 12]

Model, in band-passed space (removes haze / neuropil background), for every volume t with M = time-mean volume
(planes, y, x) and D = dM/dplane:
        I_t - M  ~=  dg_t * M  +  d_t * D                  (least squares over all pixels, closed form)
dg_t = shape-preserving gain change (calcium-like brightening, laser/PMT fluctuation); d_t = axial displacement in
PLANES (positive = the sampled stack sits where the tissue profile is higher in plane index). Gain and shift are
separable because the profile gradient relative to brightness differs between objects (some flat, some steep).

Quality checks printed and plotted:
  * split-half reliability: d_t estimated independently from the left / right (and top / bottom) half of the field
    (real axial motion is coherent over the field; estimation noise is not)
  * correlation with Suite3D's rigid z shift (registration_fused offsets), an independent estimator
  * s.e. of d_t from the Gram matrix and the empirical noise
Outputs: <diag_dir>/z_motion.npz (d, dg, se, halves), z_motion.png
"""
import argparse
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

ap = argparse.ArgumentParser()
ap.add_argument("diag")
ap.add_argument("job")
ap.add_argument("--sigma", type=float, default=6.0, help="band-pass: subtract a Gaussian of this sigma (binned px)")
ap.add_argument("--margin", type=int, default=12, help="binned px excluded at the border (registration padding)")
ap.add_argument("--movie", default="reg_bin2.npy")
a = ap.parse_args()

mov = np.load(os.path.join(a.diag, a.movie), mmap_mode="r")     # (T, nz, y, x) int16 sums
T, nz, ny, nx = mov.shape
print("movie", mov.shape)

# time-mean volume
M = np.zeros((nz, ny, nx), dtype=np.float64)
for i in range(0, T, 250):
    M += np.asarray(mov[i:i + 250], dtype=np.float32).sum(0)
M = (M / T).astype(np.float32)


def bandpass(v):          # v (nz, y, x): (1 - G) v, per plane
    return v - gaussian_filter(v, (0, a.sigma, a.sigma), mode="reflect")


D = np.gradient(M, axis=0)                  # per plane step
valid = np.zeros((ny, nx), dtype=bool)
valid[a.margin:-a.margin, a.margin:-a.margin] = True
valid &= (M > 0).all(0)
BM = bandpass(M) * valid
BD = bandpass(D) * valid


def estimate(sel):
    """d_t, dg_t for pixels in sel (bool (y,x)).
    sum_pix [(1-G) I_t] * R  =  sum_pix I_t * [(1-G) R]  for a regressor R that is zero outside sel (G is symmetric)."""
    bm, bd = BM * sel, BD * sel
    p1, p2 = bandpass(bm), bandpass(bd)
    A, B, C = (bm * bm).sum(), (bm * bd).sum(), (bd * bd).sum()
    det = A * C - B * B
    s1 = np.empty(T)
    s2 = np.empty(T)
    for i in range(0, T, 200):
        blk = np.asarray(mov[i:i + 200], dtype=np.float32)
        s1[i:i + 200] = np.tensordot(blk, p1, axes=([1, 2, 3], [0, 1, 2]))
        s2[i:i + 200] = np.tensordot(blk, p2, axes=([1, 2, 3], [0, 1, 2]))
    s1 -= s1.mean()
    s2 -= s2.mean()
    dg = (C * s1 - B * s2) / det
    d = (A * s2 - B * s1) / det
    return d, dg, (A, B, C)


full = np.ones((ny, nx), dtype=bool)
d, dg, (A, B, C) = estimate(full)
print(f"Gram: |corr(M, D)| over pixels = {abs(B) / np.sqrt(A * C):.3f} (1 = gain and shift not separable)")

xx = np.arange(nx)[None, :] * np.ones((ny, 1))
yy = np.arange(ny)[:, None] * np.ones((1, nx))
halves = {}
for name, sel in {"left": xx < nx / 2, "right": xx >= nx / 2, "top": yy < ny / 2, "bottom": yy >= ny / 2}.items():
    halves[name] = estimate(sel)[0]
r_lr = np.corrcoef(halves["left"], halves["right"])[0, 1]
r_tb = np.corrcoef(halves["top"], halves["bottom"])[0, 1]
# the reliability of the full-field estimate from a split-half r (Spearman-Brown)
rel = 2 * np.mean([r_lr, r_tb]) / (1 + np.mean([r_lr, r_tb]))
print(f"d_t: SD {d.std():.3f} planes (range {d.min():+.2f}..{d.max():+.2f}); gain SD {dg.std():.4f}")
print(f"split-half r: left/right {r_lr:+.3f}, top/bottom {r_tb:+.3f}  -> full-field reliability ~ {rel:.2f}")
print(f"estimation noise SD ~ {d.std() * np.sqrt(max(0, 1 - rel)):.3f} planes, true signal SD ~ {d.std() * np.sqrt(max(0, rel)):.3f} planes")

# Suite3D's z estimate for the same volumes
reg = os.path.join(a.job, "registered_fused_data")
offs = [np.load(os.path.join(reg, f), allow_pickle=True).item() for f in sorted(os.listdir(reg)) if f.startswith("offsets")]
sh = np.concatenate([o["sub_pixel_shifts"] for o in offs])
zs, ys, xs = sh[:, 0], sh[:, 1], sh[:, 2]
n = min(len(zs), T)
r_s3d = np.corrcoef(d[:n], zs[:n])[0, 1]
print(f"corr(my d_t, Suite3D rigid z) = {r_s3d:+.3f}   corr(my d_t, Suite3D y/x) = "
      f"{np.corrcoef(d[:n], ys[:n])[0, 1]:+.3f}/{np.corrcoef(d[:n], xs[:n])[0, 1]:+.3f}   corr(gain, d) = {np.corrcoef(dg, d)[0, 1]:+.3f}")
print(f"corr(Suite3D z, Suite3D y/x) = {np.corrcoef(zs[:n], ys[:n])[0, 1]:+.3f}/{np.corrcoef(zs[:n], xs[:n])[0, 1]:+.3f}")

# relative brightness gradient of the objects: how big is the artefact per plane of displacement?
mask = valid & (gaussian_filter(M.mean(0), 2) > np.percentile(M.mean(0), 90))
rel_grad = (np.abs(D).mean(0)[mask] / M.mean(0)[mask])
print(f"relative brightness change per plane of displacement over the brightest 10% of pixels: "
      f"median {np.median(rel_grad):.3f}, 90th pct {np.percentile(rel_grad, 90):.3f}")
print(f"=> expected z-motion artefact (SD of dF/F): median {np.median(rel_grad) * d.std():.4f}, "
      f"90th pct {np.percentile(rel_grad, 90) * d.std():.4f}")

np.savez(os.path.join(a.diag, "z_motion.npz"), d=d, dg=dg, left=halves["left"], right=halves["right"], top=halves["top"],
         bottom=halves["bottom"], s3d_z=zs, s3d_y=ys, s3d_x=xs, rel_grad=rel_grad)

fs = 5.0017
t = np.arange(T) / fs
fig, ax = plt.subplots(4, 2, figsize=(14, 11), constrained_layout=True)
ax[0, 0].plot(t, d, lw=0.5)
ax[0, 0].set_ylabel("axial displacement d_t (planes)")
ax[0, 0].set_title("registration-free estimate from the plane profile")
ax[1, 0].plot(t, dg, lw=0.5, color="C1")
ax[1, 0].set_ylabel("gain change (fraction)")
ax[2, 0].plot(t[:n], zs[:n], lw=0.5, color="C2")
ax[2, 0].set_ylabel("Suite3D rigid z (planes)")
ax[3, 0].plot(t[:n], ys[:n], lw=0.5, label="y")
ax[3, 0].plot(t[:n], xs[:n], lw=0.5, label="x")
ax[3, 0].legend()
ax[3, 0].set_ylabel("Suite3D xy shift (px)")
ax[3, 0].set_xlabel("time (s)")
ax[0, 1].scatter(halves["left"], halves["right"], s=2, alpha=0.3)
ax[0, 1].set_xlabel("d_t, left half")
ax[0, 1].set_ylabel("d_t, right half")
ax[0, 1].set_title(f"split-half r = {r_lr:+.2f}")
ax[1, 1].scatter(halves["top"], halves["bottom"], s=2, alpha=0.3)
ax[1, 1].set_xlabel("d_t, top half")
ax[1, 1].set_ylabel("d_t, bottom half")
ax[1, 1].set_title(f"split-half r = {r_tb:+.2f}")
ax[2, 1].scatter(d[:n], zs[:n], s=2, alpha=0.3, color="C2")
ax[2, 1].set_xlabel("my d_t (planes)")
ax[2, 1].set_ylabel("Suite3D rigid z (planes)")
ax[2, 1].set_title(f"r = {r_s3d:+.2f}")
f, p = np.fft.rfftfreq(T, 1 / fs), np.abs(np.fft.rfft(d - d.mean())) ** 2
ax[3, 1].loglog(f[1:], p[1:], lw=0.7)
ax[3, 1].set_xlabel("frequency (Hz)")
ax[3, 1].set_ylabel("power of d_t")
fig.savefig(os.path.join(a.diag, "z_motion.png"), dpi=90)
print("saved z_motion.png / z_motion.npz")
