"""Event-triggered averages: does anything move when the population signal fires?

    python event_triggered.py <diag_dir> <suite2p_plane_dir> [--out <png>] [--min-snr ...]

Population trace P(t) = mean dF/F (fraction of raw baseline F, neuropil-corrected, 20 s median removed) of the ROIs of one
suite2p plane that are brighter than their surround. Events = peaks of the 1 s-smoothed P above 4 robust SD that start in
a still stretch (no lateral bout or contrast dip within +-3 s before the event). Triggered averages over -6..+10 s of
  P, lateral offset, contrast gain (registered movie), my axial estimate d_t, Suite3D rigid z
against a null of random still-period trigger times. Positive control: onsets of lateral-motion bouts.
"""
import argparse
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import binary_dilation, median_filter, percentile_filter, uniform_filter1d

ap = argparse.ArgumentParser()
ap.add_argument("diag")
ap.add_argument("plane_dir")
ap.add_argument("--out", default=None)
ap.add_argument("--fs", type=float, default=5.0017)
a = ap.parse_args()
FS = a.fs
masks = np.load(os.path.join(a.diag, "audit_zsum", "roi_audit", "motion_masks.npz"))
zm = np.load(os.path.join(a.diag, "z_motion.npz"))
off, gain_dip, moving, quiet = masks["off"], masks["gain_dip"], masks["moving"], masks["quiet"]
T = len(off)
d_t, dg = zm["d"][:T], zm["dg"][:T]
s3d_z = zm["s3d_z"][:T]

F = np.load(os.path.join(a.plane_dir, "F.npy"))[:, :T].astype(np.float32)
Fn = np.load(os.path.join(a.plane_dir, "Fneu.npy"))[:, :T].astype(np.float32)
stat = np.load(os.path.join(a.plane_dir, "stat.npy"), allow_pickle=True)
ok = (np.median(F, 1) > np.median(Fn, 1)) & (np.median(F, 1) > 0) & np.array([len(s["ypix"]) >= 15 for s in stat])
F, Fn = F[ok], Fn[ok]
w20 = int(20 * FS) | 1
Fc = F - 0.7 * Fn
hp = Fc - median_filter(Fc, size=(1, w20), mode="nearest")
base = percentile_filter(F, 8, size=(1, int(30 * FS)), mode="nearest")
dff = hp / np.maximum(base, 1.0)
P = dff.mean(0)
Ps = uniform_filter1d(P, int(round(1.0 * FS)), mode="nearest")
print(f"{ok.sum()} ROIs; population dF/F SD {P.std():.3f}")

sd = 1.4826 * np.median(np.abs(Ps[quiet] - np.median(Ps[quiet])))
pre, post = int(6 * FS), int(10 * FS)
# events: local maxima of Ps above 4 SD, >= 5 s apart, still in the 3 s before the peak
from scipy.signal import find_peaks
pk_all, _ = find_peaks(Ps, height=3 * sd, distance=int(5 * FS))
pk_all = np.array([p for p in pk_all if pre <= p < T - post])
w2 = int(2 * FS)
is_still = np.array([not moving[p - w2:p + w2 + 1].any() for p in pk_all])
pk = pk_all[is_still]
print(f"{len(pk_all)} population events (peak of the 1 s-smoothed mean > 3 SD = {3 * sd:.3f} dF/F): "
      f"{len(pk)} occur in still periods (no bout / contrast dip within +-2 s), {len(pk_all) - len(pk)} during or next to motion; "
      f"still-period share of the recording = {quiet.mean():.2f}")
print(f"  still-period events: median peak {np.median(Ps[pk]):.3f} dF/F, 90th pct {np.percentile(Ps[pk], 90):.3f}")

# bout onsets (positive control): rising edge of the motion mask after >= 4 s of stillness
edge = np.where(moving[1:] & ~moving[:-1])[0] + 1
on = np.array([p for p in edge if pre <= p < T - post and not moving[p - int(4 * FS):p].any()])
print(f"{len(on)} motion-bout onsets (rising edge of the motion mask after >= 4 s of stillness)")

rng = np.random.default_rng(0)
still_idx = np.where(quiet)[0]
still_idx = still_idx[(still_idx >= pre) & (still_idx < T - post)]


def trig(sig, idx):
    """windows -pre..+post around each trigger, with the mean of -5..-1 s removed from each window"""
    X = np.stack([sig[i - pre:i + post] for i in idx])
    return X - X[:, :int(5 * FS)].mean(1, keepdims=True)


def null_band(sig, n, reps=200):
    """99% band of the triggered mean for n random trigger times in still periods"""
    m = np.stack([trig(sig, rng.choice(still_idx, n, replace=False)).mean(0) for _ in range(reps)])
    return np.percentile(m, [0.5, 99.5], axis=0)


sigs = {"population dF/F": Ps, "lateral offset (px)": off, "image contrast gain (fraction)": dg,
        "my axial estimate d_t (planes)": d_t, "Suite3D rigid z (planes)": s3d_z}
tt = np.arange(-pre, post) / FS
fig, ax = plt.subplots(len(sigs), 2, figsize=(11, 2.6 * len(sigs)), constrained_layout=True, sharex=True)
print("\ntriggered mean at 0..+2 s relative to the pre-event baseline (-5..-1 s), and the 99% band of random still-period triggers")
for r, (name, sig) in enumerate(sigs.items()):
    for c, (idx, title) in enumerate(((pk, f"population events in still periods (n={len(pk)})"),
                                      (on, f"motion-bout onsets (n={len(on)})"))):
        X = trig(sig, idx)
        m, se = X.mean(0), X.std(0) / np.sqrt(len(idx))
        lo, hi = null_band(sig, len(idx))
        ax[r, c].fill_between(tt, lo, hi, color="0.85", label="99% band, random still times")
        ax[r, c].plot(tt, m, lw=1.2)
        ax[r, c].fill_between(tt, m - 2 * se, m + 2 * se, alpha=0.3)
        ax[r, c].axvline(0, color="k", lw=0.5)
        ax[r, c].set_ylabel(name, fontsize=7)
        if r == 0:
            ax[r, c].set_title(title, fontsize=9)
        w = (tt >= 0) & (tt <= 2)
        if c == 0:
            print(f"  events  {name:34s} mean 0..2 s: {m[w].mean():+.4f}   (+-2 SE {2 * se[w].mean():.4f}); "
                  f"signal SD over the recording {sig.std():.4f}")
        else:
            print(f"  bouts   {name:34s} mean 0..2 s: {m[w].mean():+.4f}   (+-2 SE {2 * se[w].mean():.4f})")
ax[-1, 0].set_xlabel("time from trigger (s)")
ax[-1, 1].set_xlabel("time from trigger (s)")
fig.savefig(a.out or os.path.join(a.diag, "event_triggered.png"), dpi=85)
print("saved", a.out or os.path.join(a.diag, "event_triggered.png"))
