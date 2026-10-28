#!/usr/bin/env python3
"""
View single-plane suite2p output as binned, baseline-corrected dF/F
===================================================================
Reads F.npy / Fneu.npy / iscell.npy from a suite2p plane folder (no suite2p
re-run needed), and for the ROIs marked as cells:

  1. subtracts NEUCOEFF * Fneu             (suite2p's default would be 0.7)
  2. bins in time by BIN_FRAMES            (photon shot noise averages down as 1/sqrt(bin))
  3. computes dF/F against a running percentile baseline
     (BASELINE_PCTL over BASELINE_WIN_S), which follows photobleaching drift
  4. ranks ROIs by event SNR and plots them

Outputs (in PLANE_DIR/dff_view/):
  dff_examples.png          top-N ROIs, full recording, shared y-scale
  dff_neuropil_compare.png  best ROI zoomed: raw 30 Hz F vs binned dF/F, NEUCOEFF vs COMPARE_NEUCOEFF
  dff_heatmap.png           all selected ROIs x time
  dff_binned.npz            dff (rois x bins), t_s, roi_idx (row in stat.npy/F.npy), params
"""

import os
import numpy as np
from scipy.ndimage import percentile_filter, uniform_filter1d

# ============================================================================
# SETTINGS
# ============================================================================

PLANE_DIR = r"C:\Users\rbondarenko\suite2p_output\AnimalRB24_261021_1630\suite2p\plane0"

BIN_FRAMES = 5              # 30 Hz / 5 = 6 Hz. Per-bin noise: 1 fr 0.119, 3 fr 0.070, 5 fr 0.055, 10 fr 0.041 (c=0, RB23)
NEUCOEFF = 0.2              # fraction of Fneu subtracted. 0.7 removed ~60% of F for these axons and cost signal-to-noise
COMPARE_NEUCOEFF = 0.7      # only used in the comparison figure
BASELINE_PCTL = 20.0        # running baseline percentile
BASELINE_WIN_S = 60.0       # running baseline window; should be longer than the events of interest
MIN_BASELINE_ADU = 1.0      # floor on the baseline so dF/F cannot blow up when F - c*Fneu is close to 0
CELLS_ONLY = True           # iscell[:, 0] == 1 (set False to include every detected ROI)

N_EXAMPLES = 20              # rows in dff_examples.png
ZOOM_ROI = None             # row in stat.npy / F.npy for dff_neuropil_compare.png; None -> a typical ROI (median event SNR among ROIs with SNR > 5)
ZOOM_LEN_S = 60.0           # window length in dff_neuropil_compare.png (centred on that ROI's strongest event)
SHOW = True                # True: open the figures in windows after saving

# ============================================================================

import matplotlib
if not SHOW:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Reference data-viz palette (light mode)
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984"
BLUE, ORANGE = "#2a78d6", "#eb6834"
SEQ_BLUE = LinearSegmentedColormap.from_list(
    "seq_blue", [SURFACE, "#cde2fb", "#6da7ec", "#2a78d6", "#184f95", "#0d366b"])

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.edgecolor": INK2, "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 11, "axes.titleweight": "regular", "axes.titlelocation": "left",
    "font.size": 10, "legend.frameon": False,
})


def bin_time(x, n):
    """Mean over non-overlapping bins of n frames (remainder dropped). x: rois x frames."""
    nb = x.shape[1] // n
    return x[:, :nb * n].reshape(x.shape[0], nb, n).mean(axis=2)


def delta_f_over_f(F, Fneu, c, n, fs):
    """Binned dF/F against a running-percentile baseline. Returns (dff, baseline_ADU)."""
    x = bin_time(F - c * Fneu, n)
    win = int(round(BASELINE_WIN_S * fs / n)) | 1
    base = percentile_filter(x, BASELINE_PCTL, size=(1, win), mode="nearest")
    return (x - base) / np.maximum(base, MIN_BASELINE_ADU), base


def white_noise_sigma(x):
    """Robust frame-to-frame noise sigma along axis 1 (insensitive to slow signal and sparse events)."""
    return np.median(np.abs(np.diff(x, axis=1)), axis=1) / 0.6745 / np.sqrt(2)


def nice_bar(span):
    for b in (0.05, 0.1, 0.2, 0.25, 0.5, 1, 2, 5):
        if b >= span / 4:
            return b
    return 5


# ----------------------------------------------------------------------------
# Load
# ----------------------------------------------------------------------------
F = np.load(os.path.join(PLANE_DIR, "F.npy"))
Fneu = np.load(os.path.join(PLANE_DIR, "Fneu.npy"))
iscell = np.load(os.path.join(PLANE_DIR, "iscell.npy"))[:, 0] > 0
ops = np.load(os.path.join(PLANE_DIR, "ops.npy"), allow_pickle=True).item()
fs = float(ops["fs"])
out_dir = os.path.join(PLANE_DIR, "dff_view")
os.makedirs(out_dir, exist_ok=True)

roi_idx = np.flatnonzero(iscell) if CELLS_ONLY else np.arange(F.shape[0])
if roi_idx.size == 0:
    raise SystemExit("No ROIs selected - iscell.npy has no cells; set CELLS_ONLY = False.")
print(f"{F.shape[0]} ROIs, {F.shape[1]} frames @ {fs:g} Hz; using {roi_idx.size} ROIs "
      f"({'iscell==1' if CELLS_ONLY else 'all'})")
print(f"bin = {BIN_FRAMES} frames ({fs / BIN_FRAMES:g} Hz), NEUCOEFF = {NEUCOEFF}, "
      f"baseline = {BASELINE_PCTL:g}th pct over {BASELINE_WIN_S:g} s")

# ----------------------------------------------------------------------------
# dF/F + ranking
# ----------------------------------------------------------------------------
dff, base = delta_f_over_f(F[roi_idx], Fneu[roi_idx], NEUCOEFF, BIN_FRAMES, fs)
t_s = (np.arange(dff.shape[1]) + 0.5) * BIN_FRAMES / fs

noise = white_noise_sigma(dff)
score = (np.percentile(dff, 99.5, axis=1) - np.median(dff, axis=1)) / noise    # event SNR
order = np.argsort(-score)
dff, base, roi_idx, score, noise = dff[order], base[order], roi_idx[order], score[order], noise[order]

low = np.median(base, axis=1) < 3
print(f"median binned noise {np.median(noise):.3f} dF/F; event SNR (99.5th pct / noise) median {np.median(score):.1f}, "
      f"top {score[0]:.1f}, ROIs with event SNR > 5: {int((score > 5).sum())}/{len(score)}")
if low.any():
    print(f"note: {int(low.sum())} ROIs have median baseline < 3 ADU after neuropil subtraction - their dF/F is unreliable")
print("top ROIs (row in stat.npy): " + ", ".join(f"{r} ({s:.1f})" for r, s in zip(roi_idx[:N_EXAMPLES], score[:N_EXAMPLES])))

np.savez(os.path.join(out_dir, "dff_binned.npz"), dff=dff.astype(np.float32), t_s=t_s, roi_idx=roi_idx,
         event_snr=score, fs=fs, bin_frames=BIN_FRAMES, neucoeff=NEUCOEFF,
         baseline_pctl=BASELINE_PCTL, baseline_win_s=BASELINE_WIN_S)

# ----------------------------------------------------------------------------
# Figure 1: stacked example traces
# ----------------------------------------------------------------------------
n = min(N_EXAMPLES, len(roi_idx))
sub = dff[:n]
ylo, yhi = np.percentile(sub, 0.05), np.percentile(sub, 99.95)
pad = 0.05 * (yhi - ylo)
t_min = t_s / 60
fig, axes = plt.subplots(n, 1, sharex=True, sharey=True, figsize=(12, 0.95 * n + 1.3), squeeze=False)
axes = axes[:, 0]
for i, ax in enumerate(axes):
    ax.plot(t_min, sub[i], color=BLUE, lw=0.7)
    ax.set_ylim(ylo - pad, yhi + pad)
    ax.set_yticks([])
    ax.text(-0.005, 0.5, f"ROI {roi_idx[i]}", transform=ax.transAxes, ha="right", va="center", fontsize=9, color=INK2)
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.tick_params(bottom=False)
axes[-1].spines["bottom"].set_visible(True)
axes[-1].tick_params(bottom=True)
axes[-1].set_xlabel("Time (min)")
axes[-1].set_xlim(t_min[0], t_min[-1])
bar = nice_bar(yhi - ylo)
x_bar = t_min[-1] + 0.006 * (t_min[-1] - t_min[0])
axes[0].plot([x_bar, x_bar], [0, bar], color=INK2, lw=1.5, clip_on=False)
axes[0].text(x_bar + 0.004 * (t_min[-1] - t_min[0]), bar / 2, f"{bar:g} ΔF/F", va="center", fontsize=9, color=INK2)
fig.suptitle(f"Binned ΔF/F, top {n} of {len(roi_idx)} ROIs by event SNR  ·  {fs / BIN_FRAMES:g} Hz bins, "
             f"F − {NEUCOEFF}·Fneu, baseline {BASELINE_PCTL:g}th pct / {BASELINE_WIN_S:g} s",
             x=0.01, ha="left", fontsize=11)
fig.subplots_adjust(left=0.07, right=0.93, top=0.93, bottom=0.12, hspace=0)
fig.savefig(os.path.join(out_dir, "dff_examples.png"), dpi=150)

# ----------------------------------------------------------------------------
# Figure 2: one ROI, raw 30 Hz F vs binned dF/F, two neuropil coefficients
# ----------------------------------------------------------------------------
good = np.flatnonzero(score > 5)
k = int(good[len(good) // 2]) if good.size else 0           # position in the ranked arrays
if ZOOM_ROI is not None:
    k = int(np.flatnonzero(roi_idx == ZOOM_ROI)[0])
best = roi_idx[k]
peak_bin = int(np.argmax(uniform_filter1d(dff[k], 5)))
half = int(ZOOM_LEN_S * fs / BIN_FRAMES / 2)
b0 = int(np.clip(peak_bin - half, 0, dff.shape[1] - 2 * half))
b1 = b0 + 2 * half
f0, f1 = b0 * BIN_FRAMES, b1 * BIN_FRAMES
dff_cmp, _ = delta_f_over_f(F[[best]], Fneu[[best]], COMPARE_NEUCOEFF, BIN_FRAMES, fs)
t_raw = (np.arange(f0, f1) + 0.5) / fs - f0 / fs
t_bin = t_s[b0:b1] - f0 / fs

fig, (ax0, ax1) = plt.subplots(2, 1, sharex=True, figsize=(12, 5.4), gridspec_kw={"height_ratios": [1, 1.3]})
ax0.plot(t_raw, F[best, f0:f1], color=MUTED, lw=0.6)
ax0.set_ylabel("F (ADU)")
ax0.set_title(f"ROI {best} (rank {k}, event SNR {score[k]:.1f})  ·  raw F, {fs:g} Hz (t = {f0 / fs / 60:.1f}–{f1 / fs / 60:.1f} min)")
ax1.plot(t_bin, dff_cmp[0, b0:b1], color=ORANGE, lw=1.0, label=f"F − {COMPARE_NEUCOEFF}·Fneu (suite2p default)")
ax1.plot(t_bin, dff[k, b0:b1], color=BLUE, lw=1.0, label=f"F − {NEUCOEFF}·Fneu")
ax1.set_ylabel("ΔF/F")
ax1.set_xlabel("Time in window (s)")
ax1.set_title(f"Binned ΔF/F, {fs / BIN_FRAMES:g} Hz bins")
ax1.legend(loc="upper left", ncol=2, fontsize=9)
ax1.set_xlim(0, t_bin[-1])
fig.subplots_adjust(left=0.07, right=0.98, top=0.94, bottom=0.11, hspace=0.35)
fig.align_ylabels([ax0, ax1])
fig.savefig(os.path.join(out_dir, "dff_neuropil_compare.png"), dpi=150)

# ----------------------------------------------------------------------------
# Figure 3: heatmap of all selected ROIs
# ----------------------------------------------------------------------------
vmax = np.percentile(dff, 99.5)
fig, ax = plt.subplots(figsize=(12, max(3.5, 0.11 * len(roi_idx) + 1.6)))
im = ax.imshow(dff, aspect="auto", cmap=SEQ_BLUE, vmin=0, vmax=vmax, interpolation="nearest",
               extent=[t_min[0], t_min[-1], len(roi_idx), 0])
ax.set_xlabel("Time (min)")
ax.set_ylabel("ROI rank (event SNR)")
ax.set_title(f"ΔF/F, {len(roi_idx)} ROIs  ·  {fs / BIN_FRAMES:g} Hz bins, F − {NEUCOEFF}·Fneu  ·  values < 0 shown as 0")
cb = fig.colorbar(im, ax=ax, pad=0.01, fraction=0.025)
cb.set_label("ΔF/F", color=INK2)
cb.outline.set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, "dff_heatmap.png"), dpi=150)

print(f"saved to {out_dir}")
if SHOW:
    plt.show()
