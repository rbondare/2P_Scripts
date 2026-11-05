"""Is the ROI signal calcium or motion? Audit of one or more ROI sets (Suite3D 3D ROIs, suite2p per-plane ROIs).

    python roi_motion_audit.py <diag_dir> <s3d_job_dir> [--s2p <suite2p root, the folder holding plane0..4>] [--out <dir>]

Motion regressors come from the Suite3D registration (x/y shifts per volume) and from the registered-movie contrast gain
(z_motion.npz, written by z_motion_estimate.py). The tests do not assume that locomotion-locked activity is an artefact
(LC-NE axons are driven by locomotion / arousal): they ask whether the *shape* of the signal looks like calcium.

Per ROI (dF/F = (F - 0.7 Fneu - baseline) / baseline, baseline = 8th pct running 30 s):
  quiet   volumes with the lateral offset < --quiet-px from the still position, |gain dip| small, not near a bout
  moving  the complement (dilated by 2 volumes)
  asym_q  (p99 - median) / (median - p1) of dF/F in quiet volumes. Noise / symmetric artefacts ~ 1, calcium > 1.5
  sd_ratio robust SD (MAD) moving / quiet: an artefact-dominated ROI is much noisier while the brain moves
  R2_fast variance of fast dF/F (4 s running median removed) explained by fast motion regressors (y, x, y^2, x^2, yx,
          gain) vs a circular-shift null (same regressors shifted >= 120 volumes): fraction of ROIs above the null
  art_1pl the dF/F change a ONE-PLANE shift of the sampled stack would cause for this ROI (from its own plane profile)
"""
import argparse
import csv
import os
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter, percentile_filter, binary_dilation

ap = argparse.ArgumentParser()
ap.add_argument("diag")
ap.add_argument("job")
ap.add_argument("--s2p", default=None, help="folder that contains plane0..plane4 (suite2p output)")
ap.add_argument("--out", default=None)
ap.add_argument("--movie", default="regraw_bin2.npy", help="binned movie in <diag> used for the plane profiles")
ap.add_argument("--skip-s3d", action="store_true", help="do not audit the Suite3D ROIs (only --s2p)")
ap.add_argument("--s2p-name", default="suite2p_perplane", help="label of the suite2p ROI set in the outputs")
ap.add_argument("--quiet-px", type=float, default=2.0, help="lateral offset (full-res px) below which the brain counts as still")
ap.add_argument("--min-npix", type=int, default=15)
ap.add_argument("--fs", type=float, default=5.0017)
a = ap.parse_args()
out = Path(a.out or a.diag) / "roi_audit"
out.mkdir(parents=True, exist_ok=True)
FS = a.fs
NPIL = 0.7

zm = np.load(os.path.join(a.diag, "z_motion.npz"))
T = len(zm["d"])
y, x, dg = zm["s3d_y"][:T].astype(float), zm["s3d_x"][:T].astype(float), zm["dg"].astype(float)


def rolling_low(v, win):                      # slow "still position": running 10th percentile
    return percentile_filter(v, 10, size=win, mode="nearest")


win = int(120 * FS)
yo, xo = y - rolling_low(y, win), x - rolling_low(x, win)
off = np.hypot(yo, xo)                        # lateral offset from the still position (full-res px)
gain_dip = np.maximum(0, -(dg - np.median(dg)))
moving = binary_dilation((off > a.quiet_px) | (gain_dip > 0.08), iterations=2)
quiet = ~moving
print(f"{T} volumes: moving {moving.mean() * 100:.1f}%, quiet {quiet.mean() * 100:.1f}%  (lateral offset > {a.quiet_px} px or gain dip > 0.08)")

# fast motion regressors (4 s running median removed so that slow drifts / calcium plateaus do not count)
w4 = int(4 * FS) | 1


def fast(v):
    return v - median_filter(v, size=w4, mode="nearest")


regs = np.stack([fast(v) for v in (yo, xo, yo ** 2, xo ** 2, yo * xo, dg)])
regs = (regs - regs.mean(1, keepdims=True)) / (regs.std(1, keepdims=True) + 1e-9)


def dff_of(F, Fneu):
    """Baseline-removed, neuropil-corrected trace as a FRACTION OF THE ROI's RAW baseline F.

    Fc = F - 0.7 Fneu; the baseline of Fc is a 20 s running median (symmetric noise, unlike an 8th-percentile baseline,
    which makes any ratio trace look positively skewed); the denominator is the 8th-percentile running baseline of the RAW F.
    The usual (Fc - base_c) / base_c explodes here: for 39% of Suite3D ROIs median(F - 0.7 Fneu) <= 0 (ROI dimmer than its
    surround), so base_c is ~0 or negative. All tests below are scale-free except the one-plane-shift amplitude comparison.
    """
    F = F.astype(np.float32)
    Fc = F - NPIL * Fneu.astype(np.float32)
    w20 = max(3, int(round(20 * FS))) | 1
    hp = Fc - median_filter(Fc, size=(1, w20), mode="nearest")
    base_raw = percentile_filter(F, 8, size=(1, max(3, int(round(30 * FS)))), mode="nearest")
    return hp / np.maximum(base_raw, 1.0)


def robust_sd(v, axis=-1):
    return 1.4826 * np.median(np.abs(v - np.median(v, axis=axis, keepdims=True)), axis=axis)


def r2_of(Y):                                  # Y (n, T) fast traces, regs (6, T): R^2 per row
    Yc = Y - Y.mean(1, keepdims=True)
    B = np.linalg.lstsq(regs.T, Yc.T, rcond=None)[0]
    res = Yc.T - regs.T @ B
    return 1 - (res ** 2).sum(0) / ((Yc ** 2).sum(1) + 1e-12)


mean_cache = os.path.join(a.diag, os.path.splitext(a.movie)[0] + "_mean.npy")
M = np.load(mean_cache) if os.path.exists(mean_cache) else None
if M is None:
    mov = np.load(os.path.join(a.diag, a.movie), mmap_mode="r")
    M = np.zeros(mov.shape[1:], dtype=np.float64)
    for i in range(0, mov.shape[0], 250):
        M += np.asarray(mov[i:i + 250], dtype=np.float32).sum(0)
    M = (M / mov.shape[0]).astype(np.float32)
    np.save(mean_cache, M)


def one_plane_artefact(foot):
    """foot: (ys, xs, w) in FULL-res px. dF/F change caused by a one-plane shift of the sampled stack."""
    ys, xs, w = foot
    yb, xb = np.clip(ys // 2, 0, M.shape[1] - 1), np.clip(xs // 2, 0, M.shape[2] - 1)
    prof = np.array([(M[p, yb, xb] * w).sum() / w.sum() for p in range(M.shape[0])])
    slope = (prof[-1] - prof[0]) / (len(prof) - 1)          # mean change per plane across the stack
    return abs(slope) / max(prof.mean(), 1e-6), prof / max(prof.mean(), 1e-6)


def audit(name, F, Fneu, foots, extra=None):
    n = F.shape[0]
    dff = dff_of(F, Fneu)
    ok = np.isfinite(dff).all(1) & (np.median(F, 1) > 0) & (np.median(F, 1) > np.median(Fneu, 1))   # brighter than its surround
    keep = ok & np.array([len(f[0]) >= a.min_npix for f in foots])
    dq, dm = dff[:, quiet], dff[:, moving]
    med_q = np.median(dq, 1)
    p99, p1 = np.percentile(dq, 99, 1), np.percentile(dq, 1, 1)
    asym = (p99 - med_q) / np.maximum(med_q - p1, 1e-4)
    sd_q, sd_m = robust_sd(dq), robust_sd(dm)
    sd_ratio = sd_m / np.maximum(sd_q, 1e-4)
    snr_q = (p99 - med_q) / np.maximum(sd_q, 1e-4)

    Yf = np.stack([fast(r) for r in dff])
    r2 = r2_of(Yf)
    rng = np.random.default_rng(0)
    null = []
    for _ in range(30):
        s = int(rng.integers(120, T - 120))
        keep_regs = regs.copy()
        regs[:] = np.roll(keep_regs, s, axis=1)
        null.append(r2_of(Yf))
        regs[:] = keep_regs
    null = np.stack(null)                                     # (30, n)
    thr = np.percentile(null[:, keep], 99)
    art = np.array([one_plane_artefact(f)[0] for f in foots])

    # activity during moving vs quiet: fraction of large positive events that occur while moving
    thr_ev = 4 * np.maximum(sd_q, 1e-4)[:, None]
    ev = dff > thr_ev
    ev_m = (ev[:, moving].sum(1)) / np.maximum(ev.sum(1), 1)

    k = keep
    print(f"\n=== {name}: {n} ROIs, {k.sum()} analysed (>= {a.min_npix} px, median F > 0 and > median Fneu) ===")
    print(f"  dF/F robust SD in quiet volumes: median {np.median(sd_q[k]):.3f}   moving {np.median(sd_m[k]):.3f}   "
          f"ratio median {np.median(sd_ratio[k]):.2f}")
    print(f"  transient asymmetry in quiet volumes: median {np.median(asym[k]):.2f}; "
          f"{np.mean(asym[k] > 1.5) * 100:.0f}% of ROIs > 1.5 (calcium-like), {np.mean(asym[k] < 1.1) * 100:.0f}% < 1.1 (noise-like)")
    print(f"  quiet-period SNR (p99-median)/SD: median {np.median(snr_q[k]):.1f}; ROIs with SNR > 4: {np.mean(snr_q[k] > 4) * 100:.0f}%")
    print(f"  fast-motion R2: median {np.median(r2[k]):.3f} (null 99th pct {thr:.3f}); "
          f"{np.mean(r2[k] > thr) * 100:.1f}% of ROIs exceed the null (chance 1%); 90th pct R2 {np.percentile(r2[k], 90):.3f}")
    print(f"  fraction of events (> 4 SD_quiet) that fall in moving volumes: median {np.median(ev_m[k]):.2f} "
          f"(moving = {moving.mean():.2f} of the time)")
    print(f"  dF/F change from a ONE-plane stack shift: median {np.median(art[k]):.3f}, 90th pct {np.percentile(art[k], 90):.3f}; "
          f"p99 transient in quiet volumes: median {np.median((p99 - med_q)[k]):.3f}")
    good = k & (asym > 1.5) & (snr_q > 4)
    print(f"  'calcium-like' ROIs (quiet asymmetry > 1.5 and quiet SNR > 4): {good.sum()}")

    # slow (calcium-like) vs white: autocorrelation of the baseline-removed trace inside quiet stretches
    idx_k = np.where(k)[0]
    Zq = dff[idx_k]
    Zq = (Zq - Zq.mean(1, keepdims=True)) / (Zq.std(1, keepdims=True) + 1e-9)
    ac = {}
    for lag in (1, 2, 3, 5):
        both = quiet[:-lag] & quiet[lag:]
        ac[lag] = (Zq[:, :-lag][:, both] * Zq[:, lag:][:, both]).mean(1)
    print("  autocorrelation of dF/F in quiet volumes, median over ROIs: " +
          ", ".join(f"lag {l} ({l / FS:.1f}s) {np.median(v):+.3f}" for l, v in ac.items()) +
          f";  ROIs with lag-1 > 0.15: {np.mean(ac[1] > 0.15) * 100:.0f}%")
    # pair correlation vs distance (does a nearby ROI = same axon co-fluctuate?)
    cen = np.array([[f[0].mean(), f[1].mean()] for f in foots])[idx_k]
    Zs = Zq[:, quiet]
    Zs = (Zs - Zs.mean(1, keepdims=True)) / (Zs.std(1, keepdims=True) + 1e-9)
    Cm = Zs @ Zs.T / Zs.shape[1]
    D = np.hypot(cen[:, None, 0] - cen[None, :, 0], cen[:, None, 1] - cen[None, :, 1])
    iu = np.triu_indices(len(idx_k), 1)
    dd, cc = D[iu], Cm[iu]
    rng2 = np.random.default_rng(1)
    parts = []
    for lo, hi in ((0, 20), (20, 50), (50, 120), (120, 300), (300, 2000)):
        m_ = (dd >= lo) & (dd < hi)
        parts.append(f"{lo}-{hi}px: {cc[m_].mean():+.3f} (n={m_.sum()})" if m_.sum() else f"{lo}-{hi}px: n/a")
    print("  mean pair correlation (quiet volumes) vs centre distance: " + "; ".join(parts))
    # spatial scale: noise floor of a correlation estimated from this many samples
    print(f"  sampling noise of one pair correlation ~ +-{1 / np.sqrt(Zs.shape[1]):.3f}")
    # common mode
    Z = dff[k][:, quiet]
    Z = (Z - Z.mean(1, keepdims=True)) / (Z.std(1, keepdims=True) + 1e-9)
    if Z.shape[0] > 3:
        s = np.linalg.svd(Z[: min(Z.shape[0], 3000)], full_matrices=False, compute_uv=False) ** 2
        print(f"  PC1 of quiet-volume dF/F carries {s[0] / s.sum() * 100:.1f}% of the variance")
    with open(out / f"{name}_roi_stats.csv", "w", newline="") as fh:
        wtr = csv.writer(fh)
        wtr.writerow(["roi", "npix", "sd_quiet", "sd_moving", "asym_quiet", "snr_quiet", "R2_fast_motion", "frac_events_moving",
                      "artefact_one_plane", "calcium_like"] + (list(extra.keys()) if extra else []))
        for i in range(n):
            wtr.writerow([i, len(foots[i][0]), f"{sd_q[i]:.4f}", f"{sd_m[i]:.4f}", f"{asym[i]:.3f}", f"{snr_q[i]:.2f}",
                          f"{r2[i]:.4f}", f"{ev_m[i]:.3f}", f"{art[i]:.4f}", int(good[i])]
                         + ([f"{v[i]}" for v in extra.values()] if extra else []))
    return dict(dff=dff, keep=k, asym=asym, snr=snr_q, r2=r2, thr=thr, sd_q=sd_q, sd_m=sd_m, art=art, p99=p99 - med_q,
                good=good, ev_m=ev_m, name=name)


results = []
# ---------------- Suite3D 3D ROIs ----------------
rois = os.path.join(a.job, "rois")      # job/rois: stats.npy and F.npy are paired (the s3d-results copy is NOT: 9582 vs 9498)
if not a.skip_s3d and os.path.exists(os.path.join(rois, "F.npy")):
    F, Fneu = np.load(os.path.join(rois, "F.npy")), np.load(os.path.join(rois, "Fneu.npy"))
    stats = np.load(os.path.join(rois, "stats.npy"), allow_pickle=True)
    if len(stats) != F.shape[0]:
        ext = np.load(os.path.join(rois, "iscell_extracted.npy"))
        stats = stats[ext.astype(bool)]
    assert len(stats) == F.shape[0], "ROI / trace pairing"
    foots, nplanes = [], []
    for s in stats:
        zz, yy, xx = (np.asarray(c) for c in s["coords"])
        foots.append((yy, xx, np.asarray(s["lam"], dtype=float) + 1e-9))
        nplanes.append(len(np.unique(zz)))
    print("Suite3D:", F.shape)
    results.append(audit("suite3d", F, Fneu, foots, {"n_planes": nplanes}))

# ---------------- suite2p per-plane ROIs ----------------
if a.s2p:
    Fs, Fns, fts, pl = [], [], [], []
    for p in range(5):
        d = os.path.join(a.s2p, f"plane{p}")
        if not os.path.exists(os.path.join(d, "F.npy")):
            continue
        F_, Fn_ = np.load(os.path.join(d, "F.npy")), np.load(os.path.join(d, "Fneu.npy"))
        st = np.load(os.path.join(d, "stat.npy"), allow_pickle=True)
        isc = np.load(os.path.join(d, "iscell.npy"))[:, 0] > 0 if os.path.exists(os.path.join(d, "iscell.npy")) else np.ones(len(st), bool)
        for i, s in enumerate(st):
            Fs.append(F_[i][:T]); Fns.append(Fn_[i][:T])
            fts.append((np.asarray(s["ypix"]), np.asarray(s["xpix"]), np.asarray(s["lam"], dtype=float) + 1e-9))
            pl.append(p)
    if Fs:
        print("suite2p:", len(Fs), "ROIs over planes", sorted(set(pl)))
        results.append(audit(a.s2p_name, np.stack(Fs), np.stack(Fns), fts, {"plane": pl}))

# ---------------- figures ----------------
if results:
    fig, ax = plt.subplots(len(results), 4, figsize=(18, 3.8 * len(results)), constrained_layout=True, squeeze=False)
    for r, row in zip(results, ax):
        k = r["keep"]
        row[0].hist(np.clip(r["asym"][k], 0, 6), bins=60)
        row[0].axvline(1.0, color="k", ls=":")
        row[0].set_xlabel("transient asymmetry, quiet volumes")
        row[0].set_title(f"{r['name']}: n={k.sum()}")
        row[1].hist(r["r2"][k], bins=60, alpha=0.7, label="data")
        row[1].axvline(r["thr"], color="r", ls="--", label="null 99th pct")
        row[1].set_xlabel("R2 of fast dF/F on fast motion regressors")
        row[1].legend()
        row[2].scatter(r["art"][k] * 100, r["p99"][k] * 100, s=3, alpha=0.3)
        lim = max(np.percentile(r["p99"][k] * 100, 99), 10)
        row[2].plot([0, lim], [0, lim], "k:")
        row[2].set_xlim(0, np.percentile(r["art"][k] * 100, 99.5))
        row[2].set_ylim(0, lim)
        row[2].set_xlabel("dF/F % from a one-plane stack shift")
        row[2].set_ylabel("p99 transient in quiet volumes (%)")
        row[3].hist(np.clip(r["sd_m"][k] / np.maximum(r["sd_q"][k], 1e-4), 0, 6), bins=60)
        row[3].set_xlabel("noise SD moving / quiet")
    fig.savefig(out / "roi_audit_summary.png", dpi=90)

    # raster / mean traces with motion
    t = np.arange(T) / FS
    fig, ax = plt.subplots(3, 1, figsize=(15, 8), constrained_layout=True, sharex=True)
    ax[0].plot(t, off, lw=0.5)
    ax[0].plot(t, gain_dip * 10, lw=0.5, color="C3", label="gain dip x10")
    ax[0].fill_between(t, 0, ax[0].get_ylim()[1], where=moving, color="k", alpha=0.08, step="mid")
    ax[0].set_ylabel("lateral offset (px)")
    ax[0].legend(loc="upper right")
    r0 = results[0]
    sel = np.where(r0["good"])[0]
    if len(sel) == 0:
        sel = np.where(r0["keep"])[0][np.argsort(-r0["snr"][r0["keep"]])[:100]]
    mean_tr = r0["dff"][sel].mean(0)
    ax[1].plot(t, mean_tr, lw=0.5)
    ax[1].set_ylabel(f"mean dF/F of {len(sel)} best {r0['name']} ROIs")
    order = sel[np.argsort(-r0["snr"][sel])][:150]
    ax[2].imshow(np.clip(r0["dff"][order], -0.2, 1.0), aspect="auto", cmap="magma", extent=[0, t[-1], len(order), 0])
    ax[2].set_ylabel("ROI (best SNR first)")
    ax[2].set_xlabel("time (s)")
    fig.savefig(out / "roi_raster_vs_motion.png", dpi=90)
    np.savez(out / "motion_masks.npz", moving=moving, quiet=quiet, off=off, gain_dip=gain_dip)
print("\nwritten to", out)
