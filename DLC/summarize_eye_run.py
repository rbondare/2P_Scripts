"""Quick-look quality report for one eye-tracking run (output of RB_DLC_analysis.py), optionally vs an older run.

    python summarize_eye_run.py <experiment_folder> <outdir> [--old <older_experiment_folder>] [--thr 0.8]

Writes into <outdir>:  timeseries.png (gaze, pupil size, confidence over time)  and  frames_to_inspect.png
(contact sheet: least-confident frames, biggest frame-to-frame jumps, random frames), and prints a text summary.
"""
import argparse
import glob
import os

import cv2
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("exp")
ap.add_argument("outdir")
ap.add_argument("--old", default=None)
ap.add_argument("--thr", type=float, default=0.8, help="likelihood threshold the pipeline used")
ap.add_argument("--fps", type=float, default=50.0)
ap.add_argument("--batch", type=int, default=100000, help="frames per split eye video (batch_size in experiment_config.yaml)")
args = ap.parse_args()
os.makedirs(args.outdir, exist_ok=True)


def load(exp):
    eye = os.path.join(exp, "eye_results")
    name = os.path.basename(os.path.normpath(exp))
    # long recordings are split into batches (one eye video + one DLC csv each): stitch them back together
    dlc_csvs = sorted(glob.glob(os.path.join(eye, "*DLC_resnet50*.csv")))
    dlc = pd.concat([pd.read_csv(c, header=[1, 2], index_col=0) for c in dlc_csvs], ignore_index=True)
    pp = pd.read_csv(os.path.join(eye, name + "_eye.csv"), skiprows=2)
    vids = sorted(v for v in glob.glob(os.path.join(eye, "*_eye_[0-9]*.mp4")) if "DLC" not in os.path.basename(v))
    return dlc, pp, vids


dlc, pp, vids = load(args.exp)
parts = list(dict.fromkeys(dlc.columns.get_level_values(0)))
pupil = [p for p in parts if p.startswith("pupil")]
n = min(len(dlc), len(pp))
lik = np.stack([dlc[p]["likelihood"].values[:n] for p in parts], 1)
t = np.arange(n) / args.fps
f = lambda c: pp[c].astype(float).values[:n]
az, el = f("pupil_az"), f("pupil_el")
diam = np.sqrt(f("pupil_width") * f("pupil_height"))
min_pupil_lik = lik[:, [parts.index(p) for p in pupil]].min(1)

print(f"frames: {n}  ({n / args.fps / 60:.1f} min @ {args.fps:g} Hz)")
print(f"\nper-point confidence over the whole recording (threshold in use = {args.thr}):")
print(f"  {'point':18s} {'mean':>6s} {'<0.95':>7s} {'<'+str(args.thr):>7s}")
for k, p in enumerate(parts):
    print(f"  {p:18s} {lik[:, k].mean():6.3f} {100 * (lik[:, k] < 0.95).mean():6.1f}% {100 * (lik[:, k] < args.thr).mean():6.1f}%")
print("\nanalysis outcome:")
print(f"  frames with az/el = NaN            : {int(np.isnan(az).sum())} ({100 * np.isnan(az).mean():.2f}%)")
print(f"  frames where eye-outline fit failed: {int((f('eye3D_x') == 0).sum())}")
print(f"  frames where pupil fit failed      : {int((f('pupil_width') == 0).sum())}")
print(f"  frames where every pupil point >= {args.thr}: {100 * (min_pupil_lik >= args.thr).mean():.1f}%   (>= 0.95: {100 * (min_pupil_lik >= 0.95).mean():.1f}%)")
jump_az, jump_el = np.abs(np.diff(az)), np.abs(np.diff(el))
print(f"  frame-to-frame jitter (median |diff|): az {np.nanmedian(jump_az):.3f} deg, el {np.nanmedian(jump_el):.3f} deg")
print(f"  single-frame jumps > 5 deg          : {int((jump_az > 5).sum())} in az, {int((jump_el > 5).sum())} in el   (saccades are real jumps; isolated spikes are suspicious)")
print(f"  pupil diameter: median {np.nanmedian(diam):.1f} px, range {np.nanpercentile(diam, 1):.1f}-{np.nanpercentile(diam, 99):.1f} px (1st-99th pct)")
print(f"  az range {np.nanpercentile(az, 1):.1f}..{np.nanpercentile(az, 99):.1f} deg, el range {np.nanpercentile(el, 1):.1f}..{np.nanpercentile(el, 99):.1f} deg (1st-99th pct)")

oaz = oel = None
if args.old:
    _, ppo, _ = load(args.old)
    m = min(n, len(ppo))
    oaz, oel = ppo["pupil_az"].astype(float).values[:m], ppo["pupil_el"].astype(float).values[:m]
    d_az, d_el = np.abs(az[:m] - oaz), np.abs(el[:m] - oel)
    print(f"\nvs older model on the same frames ({m}):")
    print(f"  |d az| median {np.nanmedian(d_az):.2f}, 95th {np.nanpercentile(d_az, 95):.2f}, 99.9th {np.nanpercentile(d_az, 99.9):.1f} deg")
    print(f"  |d el| median {np.nanmedian(d_el):.2f}, 95th {np.nanpercentile(d_el, 95):.2f}, 99.9th {np.nanpercentile(d_el, 99.9):.1f} deg")
    print(f"  signed mean offset (new - old): az {np.nanmean(az[:m] - oaz):+.2f} deg, el {np.nanmean(el[:m] - oel):+.2f} deg")
    print(f"  jitter old: az {np.nanmedian(np.abs(np.diff(oaz))):.3f}, el {np.nanmedian(np.abs(np.diff(oel))):.3f} deg")

# ---- timeseries figure ----
fig, ax = plt.subplots(4, 1, figsize=(14, 9), sharex=True)
ax[0].plot(t, az, lw=0.4, label="new")
if oaz is not None:
    ax[0].plot(t[:len(oaz)], oaz, lw=0.3, alpha=0.6, label="old")
ax[0].set_ylabel("azimuth (deg)")
ax[0].legend(loc="upper right")
ax[1].plot(t, el, lw=0.4)
if oel is not None:
    ax[1].plot(t[:len(oel)], oel, lw=0.3, alpha=0.6)
ax[1].set_ylabel("elevation (deg)")
ax[2].plot(t, diam, lw=0.4, color="k")
ax[2].set_ylabel("pupil diameter (px)")
ax[3].plot(t, min_pupil_lik, lw=0.4, color="C3")
ax[3].axhline(args.thr, color="gray", ls="--", lw=0.8)
ax[3].set_ylabel("lowest pupil-point\nlikelihood")
ax[3].set_ylim(0, 1.02)
ax[3].set_xlabel("time (s)")
fig.suptitle(os.path.basename(os.path.normpath(args.exp)))
fig.tight_layout()
fig.savefig(os.path.join(args.outdir, "timeseries.png"), dpi=110)

# ---- contact sheet of frames worth inspecting ----
def pick_spaced(order, k, gap=100):
    chosen = []
    for i in order:
        if all(abs(int(i) - c) > gap for c in chosen):
            chosen.append(int(i))
        if len(chosen) == k:
            break
    return chosen


jumps = np.nan_to_num(np.maximum(np.r_[0, jump_az], np.r_[0, jump_el]))
groups = [("least confident", pick_spaced(np.argsort(min_pupil_lik), 6)),
          ("biggest jumps", pick_spaced(np.argsort(-jumps), 6)),
          ("random", sorted(np.random.default_rng(0).choice(n, 6, replace=False).tolist()))]
S = 2
caps = [cv2.VideoCapture(v) for v in vids]
rows = []
for title, idxs in groups:
    tiles = []
    for i in idxs:
        cap = caps[min(i // args.batch, len(caps) - 1)]
        cap.set(cv2.CAP_PROP_POS_FRAMES, i % args.batch)
        ok, fr = cap.read()
        img = cv2.resize(fr, None, fx=S, fy=S, interpolation=cv2.INTER_CUBIC) if ok else np.zeros((150 * S, 150 * S, 3), np.uint8)
        r = pp.iloc[i]
        if float(r["pupil_width"]) > 0:
            c = (int(round(float(r["pupil3d_x"]) * S)), int(round(float(r["pupil3d_y"]) * S)))
            ax_ = (int(round(float(r["pupil_width"]) / 2 * S)), int(round(float(r["pupil_height"]) / 2 * S)))
            cv2.ellipse(img, c, ax_, float(np.degrees(float(r["pupil_theta"]))), 0, 360, (0, 255, 0), 1, cv2.LINE_AA)
        for k, p in enumerate(parts):
            x, y, l = dlc[p]["x"].iloc[i], dlc[p]["y"].iloc[i], dlc[p]["likelihood"].iloc[i]
            col = (255, 220, 0) if p.startswith("pupil") else ((0, 255, 255) if p.startswith("corneal") else (255, 0, 255))
            pt = (int(round(x * S)), int(round(y * S)))
            if l >= args.thr:
                cv2.circle(img, pt, 2, col, -1, cv2.LINE_AA)
            else:
                cv2.circle(img, pt, 4, (0, 0, 255), 1, cv2.LINE_AA)
        cv2.putText(img, f"{title}: f{i} t={i / args.fps:.0f}s", (4, 14), cv2.FONT_HERSHEY_SIMPLEX, 0.4, (255, 255, 255), 1, cv2.LINE_AA)
        cv2.putText(img, f"minL={min_pupil_lik[i]:.2f} az={az[i]:.1f} el={el[i]:.1f}", (4, img.shape[0] - 6), cv2.FONT_HERSHEY_SIMPLEX, 0.38, (255, 255, 255), 1, cv2.LINE_AA)
        tiles.append(img)
    rows.append(np.hstack(tiles))
cv2.imwrite(os.path.join(args.outdir, "frames_to_inspect.png"), np.vstack(rows))
print(f"\nwrote {os.path.join(args.outdir, 'timeseries.png')} and frames_to_inspect.png")
