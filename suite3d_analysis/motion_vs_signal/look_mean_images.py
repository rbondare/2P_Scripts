"""Quick look: per-plane mean images + Suite3D rigid shift statistics. python look_mean_images.py <diag_dir> <s3d_job_dir> <out_png>"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

diag, job, out = sys.argv[1:4]
mean = np.load(os.path.join(diag, "mean_img.npy"))          # (5, 512, 1024)
nz = mean.shape[0]
print("mean intensity per plane:", mean.mean(axis=(1, 2)).round(2))
print("p99 per plane:", np.percentile(mean, 99, axis=(1, 2)).round(1), " p99.9:", np.percentile(mean, 99.9, axis=(1, 2)).round(1))

fig, ax = plt.subplots(2, nz, figsize=(4 * nz, 5.2), constrained_layout=True)
hi = np.percentile(mean, 99.7)
for p in range(nz):
    ax[0, p].imshow(mean[p], cmap="gray", vmin=0, vmax=hi)
    ax[0, p].set_title(f"plane {p} (mean)")
    hp = gaussian_filter(mean[p], 1) - gaussian_filter(mean[p], 15)
    ax[1, p].imshow(hp, cmap="gray", vmin=-np.percentile(np.abs(hp), 99.5) / 2, vmax=np.percentile(np.abs(hp), 99.5))
    ax[1, p].set_title("band-passed")
    for a in ax[:, p]:
        a.axis("off")
fig.savefig(out, dpi=90)

reg = os.path.join(job, "registered_fused_data")
offs = [np.load(os.path.join(reg, f), allow_pickle=True).item() for f in sorted(os.listdir(reg)) if f.startswith("offsets")]
sh = np.concatenate([o["sub_pixel_shifts"] for o in offs])
print("volumes registered:", len(sh))
for i, n in enumerate("zyx"):
    s = sh[:, i]
    print(f"{n}: mean {s.mean():+.2f} sd {s.std():.2f} range [{s.min():+.2f}, {s.max():+.2f}] "
          f"fraction at |max|: {np.mean(np.isclose(np.abs(s), np.abs(s).max())):.3f}")
z = sh[:, 0]
print("distinct z values:", np.unique(np.round(z, 1))[:30])
