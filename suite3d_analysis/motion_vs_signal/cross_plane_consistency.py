"""Do the five planes carry the same activity? For each plane-sum suite2p ROI, the fast-band trace of the same footprint in
each of the 5 registered planes; correlation between planes as a function of plane separation.

    python cross_plane_consistency.py <diag_dir> <suite2p_zsum_plane_dir>

Calcium in a thick section brightens the footprint in every plane together (correlation ~ the same for neighbouring and
distant planes, limited only by noise); a movement of the tissue through distinct depths would make planes diverge
(neighbouring planes correlate more than distant ones, and can anti-correlate at the profile's slopes).
"""
import os
import sys

import numpy as np
from scipy.ndimage import median_filter

diag, plane_dir = sys.argv[1:3]
movie = sys.argv[3] if len(sys.argv) > 3 else "regraw_bin2.npy"
masks_path = sys.argv[4] if len(sys.argv) > 4 else os.path.join(diag, "audit_zsum", "roi_audit", "motion_masks.npz")
FS = 5.0017
mov = np.load(os.path.join(diag, movie), mmap_mode="r")      # (T, 5, 256, 512) sums of 2x2
T, nz, ny, nx = mov.shape
stat = np.load(os.path.join(plane_dir, "stat.npy"), allow_pickle=True)
F = np.load(os.path.join(plane_dir, "F.npy"))[:, :T]
Fn = np.load(os.path.join(plane_dir, "Fneu.npy"))[:, :T]
ok = (np.median(F, 1) > np.median(Fn, 1)) & (np.median(F, 1) > 0) & np.array([len(s["ypix"]) >= 15 for s in stat])
idx = np.where(ok)[0]
masks = np.load(masks_path)
quiet = masks["quiet"]
print(len(idx), "ROIs;", quiet.sum(), "quiet volumes")

w20 = int(20 * FS) | 1
tr = np.zeros((len(idx), nz, T), dtype=np.float32)
for b in range(0, T, 250):
    blk = np.asarray(mov[b:b + 250], dtype=np.float32) / 4.0       # per-pixel mean
    for k, i in enumerate(idx):
        s = stat[i]
        yb, xb = np.clip(np.asarray(s["ypix"]) // 2, 0, ny - 1), np.clip(np.asarray(s["xpix"]) // 2, 0, nx - 1)
        lam = np.asarray(s["lam"], dtype=np.float32)
        tr[k, :, b:b + 250] = (blk[:, :, yb, xb] * lam).sum(-1).T / lam.sum()
hp = tr - median_filter(tr, size=(1, 1, w20), mode="nearest")
hp = hp[:, :, quiet]
hp = (hp - hp.mean(-1, keepdims=True)) / (hp.std(-1, keepdims=True) + 1e-9)
C = np.einsum("kpt,kqt->kpq", hp, hp) / hp.shape[-1]               # (roi, plane, plane)
Cm = C.mean(0)
np.set_printoptions(precision=3, suppress=True)
print("mean correlation between planes (fast band, quiet volumes):\n", Cm)
for sep in range(1, nz):
    v = np.concatenate([np.diagonal(C, offset=sep, axis1=1, axis2=2).ravel()])
    print(f"  plane separation {sep}: median {np.median(v):+.3f}, mean {v.mean():+.3f}, 10-90 pct [{np.percentile(v, 10):+.3f}, {np.percentile(v, 90):+.3f}]")
# per-plane signal brightness of the footprint (raw), as a fraction of the plane-mean
prof = tr.mean(-1)
print("footprint brightness by plane, relative to the ROI's plane-mean: median", np.median(prof / prof.mean(1, keepdims=True), 0).round(3))
# anti-correlation between neighbouring planes (would signal a z-shift through a profile slope)
v1 = np.diagonal(C, offset=1, axis1=1, axis2=2).ravel()
print(f"fraction of neighbouring-plane pairs with correlation < 0: {np.mean(v1 < 0):.3f}")
