"""How different are the 5 planes really? Axial intensity profile of the brightest local maxima of the time-mean image.
python axial_profile_check.py <diag_dir> [out_png]

If the planes sampled distinct depths (3 um apart, bouton axial extent ~ 2-4 um + PSF), different objects should peak in
different planes and most objects should fall off by >30-50% within 1-2 planes. Flat profiles = the five planes are
(nearly) the same optical section: no usable z information, and axial motion of a few um would barely change intensity.
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, maximum_filter, uniform_filter

diag = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else None
m = np.load(os.path.join(diag, "mean_img.npy"))            # (nz, y, x) full res
nz, ly, lx = m.shape
sm = gaussian_filter(m, (0, 3, 3))                           # ~ bouton scale (0.11 um/px -> 1 um sigma)
bg = uniform_filter(m, (1, 61, 61))
res = sm - gaussian_filter(m, (0, 15, 15))
ref = res.mean(axis=0)
loc = (ref == maximum_filter(ref, 13)) & (ref > np.percentile(ref, 98.5)) & (ref < np.percentile(ref, 99.97))
loc[:30] = loc[-30:] = False
loc[:, :30] = loc[:, -30:] = False
ys, xs = np.where(loc)
print(len(ys), "bright local maxima")
prof = np.array([res[:, y, x] for y, x in zip(ys, xs)])      # (n, nz) band-passed brightness
ok = prof.max(axis=1) > 0
prof = prof[ok]
norm = prof / prof.max(axis=1, keepdims=True)
print("mean normalised profile over objects (plane 0..4):", norm.mean(0).round(3))
print("fraction of objects peaking in each plane:", np.bincount(prof.argmax(1), minlength=nz) / len(prof))
print("median min/max ratio over the 5 planes:", np.median(norm.min(1)).round(3),
      " (1 = flat across 12 um; a bouton with a 3 um FWHM would be <0.1)")
print("median |plane p - plane p+1| / mean:", np.median(np.abs(np.diff(prof, axis=1)).mean(1) / prof.mean(1)).round(3))
# object-by-object correlation of the profile with the object's brightness: do profiles differ between objects?
prof_c = norm - norm.mean(0)
print("SD across objects of the normalised profile (per plane):", norm.std(0).round(3))

fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
ax[0].plot(norm[:60].T, color="k", alpha=0.25, lw=0.8)
ax[0].plot(norm.mean(0), color="r", lw=2)
ax[0].set_xlabel("plane (3 um apart)")
ax[0].set_ylabel("band-passed brightness / max")
ax[0].set_title(f"axial profile of {len(prof)} bright objects")
ax[1].hist(norm.min(1), bins=25)
ax[1].set_xlabel("min/max across the 5 planes")
ax[1].set_ylabel("objects")
if out:
    fig.savefig(out, dpi=100)
