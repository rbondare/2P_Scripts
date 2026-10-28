#!/usr/bin/env python3
"""
Merge axon segments of one suite2p plane by activity correlation
================================================================
Detection splits one axon into several short segments. This script joins the
segments that belong together, outside the suite2p GUI (whose merge does not
save reliably). For one recording it:

  1. copies <rec>/suite2p -> <rec>/suite2p_unmerged, once. The backup is never
     overwritten, and a re-run always merges starting from it, so you can try other
     thresholds as often as you like. data.bin and reg_tif/ are not copied
     (never modified) unless COPY_BINARY = True.
  2. finds groups of ROIs to merge. Traces are binned dF/F (same recipe as
     view_suite2p_dff.py) with the population-wide component (top PC) projected
     out, because most axons share a common signal. Then greedy agglomeration:
     repeatedly join the two groups whose summed trace correlates best, as long as
     their nearest pixels are <= MAX_GAP_PX apart and r >= CORR_THRESH.
  3. builds one ROI per group (union of the pixels) and re-extracts F / Fneu for it
     from the registered binary with suite2p's own mask and extraction code (the
     neuropil ring is built around the whole union), instead of averaging the
     constituents' traces like the GUI does. A few untouched cells are re-extracted
     too as a check that the binary matches F.npy.
  4. writes stat / F / Fneu / spks / iscell / Fall.mat into <rec>/suite2p/plane<N>
     in the GUI's merge convention: merged ROIs are appended (imerge = rows they
     were built from, inmerge = -1); the constituents keep their rows but get
     inmerge = row of the merged ROI and iscell = 0.

Everything it writes besides those files goes to <plane>/merge_report/:
  merge_log.json    parameters, the merge steps and per-group noise / event SNR
  merge_pairs.png   correlation vs gap for every pair of candidate ROIs (tune CORR_THRESH / MAX_GAP_PX here)
  merge_map.png     where the merged groups are in the field of view
  merge_traces.png  best single segment vs merged ROI for the largest groups

Usage:
    conda activate suite2p_venv
    python merge_suite2p_rois.py                 # uses REC_DIR below
    python merge_suite2p_rois.py "<recording folder>"
    python merge_suite2p_rois.py --dry-run       # or DRY_RUN = True: only write merge_report/, change nothing else
"""

import os
import sys
import json
import time
import shutil
import numpy as np
import torch
from pathlib import PureWindowsPath
from scipy.ndimage import percentile_filter
from scipy.spatial import cKDTree
from scipy.stats import skew

from suite2p.io import BinaryFile
from suite2p.detection.stats import roi_stats
from suite2p.extraction.dcnv import oasis
from suite2p.extraction.extract import extract_traces
from suite2p.extraction.masks import create_cell_pix, create_cell_mask, create_neuropil_masks

from export_suite2p_to_mat import export_plane

# ============================================================================
# SETTINGS
# ============================================================================

REC_DIR = r"Z:\Group Members\Rima\DATA_2P\AnimalRB22_261023_1527"   # recording folder that contains suite2p/ (a command-line argument overrides this)
PLANE = 0
BINARY_PATH = None          # registered binary; None -> <suite2p>/plane<PLANE>/data.bin
COPY_BINARY = False         # True: also copy data.bin / reg_tif into suite2p_unmerged (~63 GB each, unchanged by the merge)
DRY_RUN = False             # True: only write merge_report/ (no backup, no merge)

CELLS_ONLY = True           # only ROIs with iscell == 1 take part (curate in the GUI first)

# traces used for the correlation - same recipe as view_suite2p_dff.py
BIN_FRAMES = 5              # 30 Hz / 5 = 6 Hz
NEUCOEFF = 0.2              # fraction of Fneu subtracted for the correlation only; F / Fneu are saved untouched
BASELINE_PCTL = 20.0        # running baseline percentile
BASELINE_WIN_S = 60.0       # running baseline window
MIN_BASELINE_ADU = 1.0      # floor on the baseline so dF/F cannot blow up

# merge criteria
REMOVE_POP_PCS = 1          # population-wide principal components projected out before correlating (0 = plain correlation)
CORR_THRESH = 0.3           # min correlation (after PC removal) between two groups ...
MAX_GAP_PX = 10             # ... whose nearest pixels are at most this far apart (None = no distance limit)

N_CHECK_ROIS = 3            # untouched cells re-extracted alongside the merged ROIs as a check of the binary
CHECK_TOL = 0.02            # abort without writing if their mean |F_new - F_old| / F_old exceeds this
N_EXAMPLES = 12             # groups shown in merge_traces.png
SHOW = False                # True: open the figures in windows after saving

# ============================================================================

ARGS = [a for a in sys.argv[1:] if not a.startswith("--")]
DRY_RUN = DRY_RUN or "--dry-run" in sys.argv[1:]

PARAMS = ("PLANE", "COPY_BINARY", "CELLS_ONLY", "BIN_FRAMES", "NEUCOEFF", "BASELINE_PCTL", "BASELINE_WIN_S",
          "MIN_BASELINE_ADU", "REMOVE_POP_PCS", "CORR_THRESH", "MAX_GAP_PX", "N_CHECK_ROIS", "CHECK_TOL")

import matplotlib
if not SHOW:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Reference data-viz palette (light mode), categorical slots in fixed order
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984"
BLUE, ORANGE = "#2a78d6", "#eb6834"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.edgecolor": INK2, "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 11, "axes.titleweight": "regular", "axes.titlelocation": "left",
    "font.size": 10, "legend.frameon": False,
})

# Z: is mapped at different levels of the share depending on the session
Z_ROOTS = ("Z:\\group\\joeschgrp\\", "Z:\\joeschgrp\\", "Z:\\")


def resolve_runtime_path(path_str):
    """Resolve a Windows-style path on Windows or WSL/Linux."""
    if os.name == "nt":
        if os.path.exists(path_str):
            return path_str
        for root in Z_ROOTS:
            if path_str.lower().startswith(root.lower()):
                tail = path_str[len(root):]
                for alt in Z_ROOTS:
                    if os.path.exists(alt + tail):
                        return alt + tail
                break
        return path_str

    windows_path = PureWindowsPath(path_str)
    drive = windows_path.drive
    if drive:
        drive_letter = drive.rstrip(":").lower()
        return os.path.join("/mnt", drive_letter, *windows_path.parts[1:])

    return path_str.replace("\\", "/")


# ----------------------------------------------------------------------------
# Traces (same definitions as view_suite2p_dff.py)
# ----------------------------------------------------------------------------
def bin_time(x, n):
    """Mean over non-overlapping bins of n frames (remainder dropped). x: rois x frames."""
    nb = x.shape[1] // n
    return x[:, :nb * n].reshape(x.shape[0], nb, n).mean(axis=2)


def delta_f_over_f(F, Fneu, fs):
    """Binned dF/F against a running-percentile baseline. Returns (dff, baseline_ADU)."""
    x = bin_time(F - NEUCOEFF * Fneu, BIN_FRAMES)
    win = int(round(BASELINE_WIN_S * fs / BIN_FRAMES)) | 1
    base = percentile_filter(x, BASELINE_PCTL, size=(1, win), mode="nearest")
    return (x - base) / np.maximum(base, MIN_BASELINE_ADU), base


def white_noise_sigma(x):
    """Robust frame-to-frame noise sigma along axis 1 (insensitive to slow signal and sparse events)."""
    return np.median(np.abs(np.diff(x, axis=1)), axis=1) / 0.6745 / np.sqrt(2)


def event_snr(dff):
    """(99.5th percentile - median) / noise per ROI, and the noise."""
    noise = white_noise_sigma(dff)
    return (np.percentile(dff, 99.5, axis=1) - np.median(dff, axis=1)) / noise, noise


# ----------------------------------------------------------------------------
# Backup / load
# ----------------------------------------------------------------------------
def backup_unmerged(suite2p_dir, unmerged_dir):
    """Copy suite2p/ to suite2p_unmerged/ (via a temp name, so an interrupted copy is never mistaken for a backup)."""
    tmp_dir = unmerged_dir + "_partial"
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    skip = () if COPY_BINARY else ("data.bin", "reg_tif")
    print(f"Copying {suite2p_dir}\n     -> {unmerged_dir}" + ("" if COPY_BINARY else f"   (without {', '.join(skip)})"))
    t0 = time.time()
    shutil.copytree(suite2p_dir, tmp_dir, ignore=shutil.ignore_patterns(*skip))
    os.rename(tmp_dir, unmerged_dir)
    print(f"  done in {time.time() - t0:.0f} s")


def load_plane(plane_dir):
    def load(fname, **kw):
        return np.load(os.path.join(plane_dir, fname), **kw)
    stat = load("stat.npy", allow_pickle=True)
    F, Fneu, spks, iscell = load("F.npy"), load("Fneu.npy"), load("spks.npy"), load("iscell.npy")
    ops = load("ops.npy", allow_pickle=True).item()
    return stat, F, Fneu, spks, iscell, ops


def find_binary(plane_dirs):
    candidates = ([resolve_runtime_path(BINARY_PATH)] if BINARY_PATH else
                  [os.path.join(d, "data.bin") for d in plane_dirs])
    for c in candidates:
        if os.path.exists(c):
            return c
    raise SystemExit("Registered binary not found. Looked for:\n  " + "\n  ".join(candidates) +
                     "\nSet BINARY_PATH (the merge reads it but never modifies it).")


# ----------------------------------------------------------------------------
# Which ROIs take part, and how far apart they are
# ----------------------------------------------------------------------------
def merge_units(stat, iscell):
    """
    ROIs that may merge: cells that are not hidden constituents of an earlier (GUI) merge.
    An earlier merged ROI takes part as one unit and stands for its imerge rows.
    Returns the stat row of each unit and the original ROI rows it stands for.
    """
    rows, leaves = [], []
    for n, s in enumerate(stat):
        if s.get("inmerge", 0) > 0:
            continue
        if CELLS_ONLY and iscell[n, 0] < 0.5:
            continue
        im = s.get("imerge")
        leaves.append([int(k) for k in im] if im is not None and len(im) > 0 else [n])
        rows.append(n)
    return np.array(rows), leaves


def unit_pixels(stat, rows, leaves):
    """ypix, xpix, lam per unit. lam is rescaled to the summed lam of the original ROIs it stands for,
    so a GUI-merged ROI (whose lam was normalised) keeps the weight of its pieces."""
    out = []
    for r, lv in zip(rows, leaves):
        s = stat[r]
        mass = sum(stat[k]["lam"].sum() for k in lv)
        out.append((s["ypix"], s["xpix"], s["lam"] * (mass / s["lam"].sum())))
    return out


def pixel_gaps(pix):
    """Distance in pixels between the nearest pixels of every pair of units."""
    trees = [cKDTree(np.c_[y, x]) for y, x, _ in pix]
    n = len(pix)
    gap = np.zeros((n, n))
    for a in range(n):
        for b in range(a + 1, n):
            gap[a, b] = gap[b, a] = trees[a].query(np.c_[pix[b][0], pix[b][1]])[0].min()
    return gap


# ----------------------------------------------------------------------------
# Correlation + agglomeration
# ----------------------------------------------------------------------------
def population_pcs(d, k):
    """Top-k temporal principal components of the z-scored traces (k x bins, orthonormal) and the variance they explain."""
    if k <= 0:
        return np.zeros((0, d.shape[1])), np.zeros(0)
    z = d - d.mean(axis=1, keepdims=True)
    z /= z.std(axis=1, keepdims=True) + 1e-12
    _, s, vt = np.linalg.svd(z, full_matrices=False)
    return vt[:k], (s ** 2 / (s ** 2).sum())[:k]


def residual_vectors(d, pcs):
    """Rows of d: mean removed, unit variance, population PCs projected out, unit length.
    The dot product of two rows is their correlation after PC removal."""
    z = d - d.mean(axis=1, keepdims=True)
    z = z / (z.std(axis=1, keepdims=True) + 1e-12)
    if len(pcs):
        z = z - (z @ pcs.T) @ pcs
    return z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-12)


def agglomerate(d, w, gap, pcs):
    """
    Greedy agglomeration. d: unit dF/F (units x bins), w: how much each unit contributes to a summed trace,
    gap: unit-to-unit pixel gap. Each step joins the pair of groups with the highest correlation (between
    their weighted-mean dF/F, after PC removal) among pairs closer than MAX_GAP_PX, until none reaches CORR_THRESH.
    Because the group trace is recomputed after every join, a noisy segment can join a group it only
    correlates with once the group's own noise has averaged down.
    Returns the groups (lists of unit indices), each group's estimated dF/F, and the steps taken.
    """
    n, T = d.shape
    max_gap = np.inf if MAX_GAP_PX is None else MAX_GAP_PX
    vec = np.zeros((2 * n, T))
    vec[:n] = residual_vectors(d, pcs)
    R = np.full((2 * n, 2 * n), -np.inf)
    R[:n, :n] = vec[:n] @ vec[:n].T
    np.fill_diagonal(R, -np.inf)
    GAP = np.full((2 * n, 2 * n), np.inf)
    GAP[:n, :n] = gap
    wsum = np.zeros((2 * n, T))
    wsum[:n] = w[:, None] * d
    wtot = np.zeros(2 * n)
    wtot[:n] = w
    members = [[i] for i in range(n)] + [None] * n
    alive = np.zeros(2 * n, bool)
    alive[:n] = True
    steps = []
    for c in range(n, 2 * n):
        cand = np.where(GAP <= max_gap, R, -np.inf)
        a, b = np.unravel_index(np.argmax(cand), cand.shape)
        if cand[a, b] < CORR_THRESH:
            break
        steps.append({"a": members[a], "b": members[b], "r": float(cand[a, b]), "gap_px": float(GAP[a, b])})
        members[c] = members[a] + members[b]
        wsum[c], wtot[c] = wsum[a] + wsum[b], wtot[a] + wtot[b]
        vec[c] = residual_vectors((wsum[c] / wtot[c])[None], pcs)[0]
        alive[[a, b]] = False
        R[[a, b], :] = -np.inf
        R[:, [a, b]] = -np.inf
        GAP[c] = np.minimum(GAP[a], GAP[b])
        GAP[:, c] = GAP[c]
        r = np.where(alive, vec @ vec[c], -np.inf)
        R[c], R[:, c] = r, r
        R[c, c] = -np.inf
        alive[c] = True
    ids = [k for k in range(2 * n) if alive[k] and len(members[k]) > 1]
    return [sorted(members[k]) for k in ids], [wsum[k] / wtot[k] for k in ids], steps


# ----------------------------------------------------------------------------
# Merged ROIs
# ----------------------------------------------------------------------------
def union_roi(pix, group, leaves, Ly, Lx, diameter, stat):
    """stat entry for the union of the group's units (a shared pixel keeps its larger lam)."""
    y = np.concatenate([pix[u][0] for u in group])
    x = np.concatenate([pix[u][1] for u in group])
    lam = np.concatenate([pix[u][2] for u in group])
    key = y.astype(np.int64) * Lx + x
    order = np.lexsort((-lam, key))
    first = np.r_[True, key[order][1:] != key[order][:-1]]
    sel = order[first]
    new = {"ypix": y[sel], "xpix": x[sel], "lam": lam[sel],
           "imerge": np.array(sorted({k for u in group for k in leaves[u]})),
           "inmerge": -1, "chan2_prob": -1}
    new = roi_stats(np.array([new], dtype=object), Ly, Lx, diameter=diameter,
                    max_overlap=None, do_soma_crop=False)[0]
    # roi_stats normalises npix by the ROI itself; use the normaliser of the detected ROIs instead
    ref = next((s for s in stat if s.get("npix_norm", 0) > 0), None)
    if ref is not None:
        new["npix_norm"] = new["npix_soma"] * ref["npix_norm"] / ref["npix_soma"]
        new["npix_norm_no_crop"] = new["npix"] * ref["npix_norm_no_crop"] / ref["npix"]
    return new


def extract_rows(bin_path, stat_all, rows, ops, n_frames):
    """F, Fneu of stat_all[rows] from the registered binary, exactly as suite2p's extraction does it."""
    ext = ops["extraction"]
    Ly, Lx = int(ops["Ly"]), int(ops["Lx"])
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    with BinaryFile(Ly, Lx, str(bin_path)) as f:
        if f.shape[0] != n_frames:
            raise SystemExit(f"{bin_path} has {f.shape[0]} frames but F.npy has {n_frames}: wrong binary.")
        # allow_overlap=True, as in the run scripts
        cell_pix = create_cell_pix(stat_all, Ly, Lx, lam_percentile=ext["lam_percentile"])
        cell_masks = [create_cell_mask(stat_all[r], Ly, Lx, allow_overlap=True) for r in rows]
        neuropil_masks = create_neuropil_masks(
            ypixs=[stat_all[r]["ypix"] for r in rows], xpixs=[stat_all[r]["xpix"] for r in rows],
            cell_pix=cell_pix, inner_neuropil_radius=ext["inner_neuropil_radius"],
            min_neuropil_pixels=ext["min_neuropil_pixels"], circular=ext["circular_neuropil"])
        return extract_traces(f, cell_masks, neuropil_masks, batch_size=ext["batch_size"], device=device)


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------
def group_color(g):
    return SERIES[g] if g < len(SERIES) else INK2     # past 8 groups the number label carries the identity


def plot_pairs(path, r, gap, together):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(gap[~together], r[~together], s=10, color=MUTED, alpha=0.5, linewidths=0, label="other pairs")
    ax.scatter(gap[together], r[together], s=22, color=BLUE, linewidths=0, label="ended up in the same group")
    ax.axhline(CORR_THRESH, color=ORANGE, lw=1.2, ls="--")
    ax.text(0.99, CORR_THRESH, " CORR_THRESH", color=INK2, ha="right", va="bottom", transform=ax.get_yaxis_transform(), fontsize=9)
    if MAX_GAP_PX is not None:
        ax.axvline(MAX_GAP_PX, color=ORANGE, lw=1.2, ls="--")
        ax.text(MAX_GAP_PX, 1.0, " MAX_GAP_PX", color=INK2, ha="left", va="top", transform=ax.get_xaxis_transform(), fontsize=9)
    ax.set_xscale("symlog", linthresh=10)
    ax.set_xlim(left=0)
    ax.set_xlabel("gap between nearest pixels (px)")
    ax.set_ylabel(f"correlation of dF/F, {REMOVE_POP_PCS} population PC removed")
    ax.set_title(f"{len(r)} pairs of candidate ROIs  ·  single-ROI correlations; groups are decided on summed traces")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_map(path, ops, pix, groups):
    Ly, Lx = int(ops["Ly"]), int(ops["Lx"])
    img = np.asarray(ops.get("meanImg", np.zeros((Ly, Lx))))
    fig, ax = plt.subplots(figsize=(12, 12 * Ly / Lx + 0.8))
    ax.imshow(img, cmap="gray_r", vmin=np.percentile(img, 1), vmax=np.percentile(img, 99.5))

    def outline(u, **kw):
        m = np.zeros((Ly, Lx))
        m[pix[u][0], pix[u][1]] = 1
        ax.contour(m, levels=[0.5], **kw)

    in_group = {u for g in groups for u in g}
    for u in range(len(pix)):
        if u not in in_group:
            outline(u, colors=MUTED, linewidths=0.6)
    for k, g in enumerate(groups):
        for u in g:
            outline(u, colors=group_color(k), linewidths=1.6)
        y = np.concatenate([pix[u][0] for u in g])
        x = np.concatenate([pix[u][1] for u in g])
        ax.text(np.median(x), np.median(y), str(k + 1), ha="center", va="center", fontsize=9,
                color=INK, bbox=dict(boxstyle="round,pad=0.15", fc=SURFACE, ec="none", alpha=0.85))
    ax.set_xlim(0, Lx)
    ax.set_ylim(Ly, 0)
    ax.axis("off")
    ax.set_title(f"{len(groups)} merged groups (numbered, as in merge_traces.png) from "
                 f"{len(in_group)} of {len(pix)} candidate ROIs  ·  grey: ROIs left alone")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_traces(path, groups, best_dff, best_row, merged_dff, fs, estimated):
    n = min(N_EXAMPLES, len(groups))
    order = np.argsort([-len(g) for g in groups], kind="stable")[:n]
    t_min = (np.arange(merged_dff.shape[1]) + 0.5) * BIN_FRAMES / fs / 60
    fig, axes = plt.subplots(2 * n, 1, sharex=True, figsize=(12, 1.9 * n + 1.2), squeeze=False)
    axes = axes[:, 0]
    for i, k in enumerate(order):
        a, b = axes[2 * i], axes[2 * i + 1]
        both = np.r_[best_dff[k], merged_dff[k]]
        ylo, yhi = np.percentile(both, 0.05), np.percentile(both, 99.95)
        pad = 0.05 * (yhi - ylo)
        for ax, trace, color in ((a, best_dff[k], MUTED), (b, merged_dff[k], BLUE)):
            ax.plot(t_min, trace, color=color, lw=0.7)
            ax.set_ylim(ylo - pad, yhi + pad)
            ax.set_yticks([])
            for side in ("left", "bottom"):
                ax.spines[side].set_visible(False)
            ax.tick_params(bottom=False)
        snr, noise = event_snr(np.vstack([best_dff[k], merged_dff[k]]))
        a.set_title(f"{k + 1}: ROI {best_row[k]}, best of {len(groups[k])} segments   noise {noise[0]:.3f}   event SNR {snr[0]:.1f}",
                    fontsize=9, color=INK2)
        b.set_title(f"merged{' (estimate, not re-extracted)' if estimated else ''}   noise {noise[1]:.3f}   event SNR {snr[1]:.1f}",
                    fontsize=9, color=INK2)
    axes[-1].spines["bottom"].set_visible(True)
    axes[-1].tick_params(bottom=True)
    axes[-1].set_xlabel("Time (min)")
    axes[-1].set_xlim(t_min[0], t_min[-1])
    fig.subplots_adjust(left=0.03, right=0.99, top=0.97, bottom=0.07, hspace=0.55)
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------------------
def main():
    rec_dir = resolve_runtime_path(ARGS[0] if ARGS else REC_DIR)
    suite2p_dir = os.path.join(rec_dir, "suite2p")
    unmerged_dir = os.path.join(rec_dir, "suite2p_unmerged")
    plane = f"plane{PLANE}"
    live_dir = os.path.join(suite2p_dir, plane)
    if not os.path.isdir(live_dir):
        raise SystemExit(f"Not found: {live_dir}")

    print("=" * 70)
    print(f"MERGE ROIs  ·  {rec_dir}  ·  {plane}" + ("  ·  DRY RUN" if DRY_RUN else ""))
    print("=" * 70)
    if not DRY_RUN and not os.path.isdir(unmerged_dir):
        backup_unmerged(suite2p_dir, unmerged_dir)
    src_dir = os.path.join(unmerged_dir, plane)
    if not os.path.isdir(src_dir):
        src_dir = live_dir
    print(f"merging from: {src_dir}")
    report_dir = os.path.join(live_dir, "merge_report")
    os.makedirs(report_dir, exist_ok=True)

    stat, F, Fneu, spks, iscell, ops = load_plane(src_dir)
    fs = float(ops["fs"])
    Ly, Lx = int(ops["Ly"]), int(ops["Lx"])
    diameter = np.broadcast_to(np.atleast_1d(ops["diameter"]), (2,)).tolist()
    n_roi, n_frames = F.shape
    print(f"{n_roi} ROIs, {n_frames} frames @ {fs:g} Hz")

    # ---- decide the groups -------------------------------------------------
    rows, leaves = merge_units(stat, iscell)
    if len(rows) < 3:
        raise SystemExit(f"Only {len(rows)} candidate ROIs - nothing to merge.")
    pix = unit_pixels(stat, rows, leaves)
    dff, base = delta_f_over_f(F[rows], Fneu[rows], fs)
    w = np.array([sum(stat[k]["lam"].sum() for k in lv) for lv in leaves]) * np.maximum(np.median(base, axis=1), MIN_BASELINE_ADU)
    pcs, explained = population_pcs(dff, min(REMOVE_POP_PCS, len(rows) - 2))
    print(f"{len(rows)} candidate ROIs; population PCs removed: {len(pcs)}"
          + (f" (explain {', '.join(f'{100 * e:.0f}%' for e in explained)} of the variance)" if len(pcs) else ""))
    gap = pixel_gaps(pix)
    groups, est_dff, steps = agglomerate(dff, w, gap, pcs)
    print(f"\n{len(steps)} joins -> {len(groups)} merged groups "
          f"(CORR_THRESH {CORR_THRESH}, MAX_GAP_PX {MAX_GAP_PX}, REMOVE_POP_PCS {REMOVE_POP_PCS})")
    for s in steps:
        print(f"  ROIs {[int(rows[u]) for u in s['a']]} + {[int(rows[u]) for u in s['b']]}   r = {s['r']:.2f}   gap = {s['gap_px']:.0f} px")

    vec = residual_vectors(dff, pcs)
    iu = np.triu_indices(len(rows), 1)
    group_of = {u: k for k, g in enumerate(groups) for u in g}
    together = np.array([group_of.get(a, -1) == group_of.get(b, -2) for a, b in zip(*iu)])
    plot_pairs(os.path.join(report_dir, "merge_pairs.png"), (vec @ vec.T)[iu], gap[iu], together)
    plot_map(os.path.join(report_dir, "merge_map.png"), ops, pix, groups)
    best_u = [max(g, key=lambda u: event_snr(dff[[u]])[0][0]) for g in groups]
    best_dff = [dff[u] for u in best_u]

    if DRY_RUN or (not groups and src_dir == live_dir):
        if groups:
            plot_traces(os.path.join(report_dir, "merge_traces.png"), groups, best_dff,
                        [int(rows[u]) for u in best_u], np.array(est_dff), fs, estimated=True)
        print(f"\n{'Dry run' if DRY_RUN else 'No groups found'} - nothing changed; report in {report_dir}")
        return
    if not groups:
        print("\nNo groups found - suite2p/ is reset to the unmerged result.")

    # ---- build + extract the merged ROIs -------------------------------------
    new_stats = [union_roi(pix, g, leaves, Ly, Lx, diameter, stat) for g in groups]
    first_new = n_roi
    stat_all = np.concatenate([stat, np.array(new_stats, dtype=object)]) if groups else stat
    new_rows = list(range(first_new, first_new + len(groups)))
    grouped = {u for g in groups for u in g}
    # controls must be plain detected ROIs: the stored F of an earlier GUI merge is an average of its pieces, not an extraction
    free = [int(r) for u, r in enumerate(rows) if u not in grouped and leaves[u] == [int(r)]]
    check_rows = [free[i] for i in np.linspace(0, len(free) - 1, min(N_CHECK_ROIS, len(free))).astype(int)] if free else []

    if groups:
        bin_path = find_binary([live_dir, src_dir])
        print(f"\nExtracting {len(new_rows)} merged ROIs (+{len(check_rows)} check ROIs) from {bin_path}")
        t0 = time.time()
        F_all, Fneu_all = extract_rows(bin_path, stat_all, new_rows + check_rows, ops, n_frames)
        print(f"  done in {time.time() - t0:.0f} s")
        F_new, Fneu_new = F_all[:len(new_rows)], Fneu_all[:len(new_rows)]
        for r, Fc in zip(check_rows, F_all[len(new_rows):]):
            dev = float(np.abs(Fc - F[r]).mean() / np.abs(F[r]).mean())
            print(f"  check ROI {r}: mean |F_new - F_old| / F_old = {dev:.4f}")
            if dev > CHECK_TOL:
                raise SystemExit(f"Re-extraction does not reproduce F.npy for ROI {r} ({dev:.3f} > {CHECK_TOL}): "
                                 "the binary does not match this stat/F. Nothing was written.")
    else:
        F_new = Fneu_new = np.zeros((0, n_frames), np.float32)

    # ---- activity stats, spks, bookkeeping -----------------------------------
    neucoeff = ops.get("neucoeff", ops["extraction"]["neuropil_coefficient"])
    dF = F_new - neucoeff * Fneu_new
    for s, trace in zip(new_stats, dF):
        s["skew"] = float(skew(trace))
        s["std"] = float(trace.std())
        s["snr"] = float(1 - 0.5 * np.diff(trace).var() / (trace.var() + 1e-10))
    if len(groups) and ops.get("run", {}).get("do_deconvolution", True):
        spks_new = oasis(F=dF, batch_size=ops["batch_size"], tau=ops["tau"], fs=ops["fs"])
    else:
        spks_new = np.zeros((len(groups), n_frames), np.float32)

    iscell = iscell.copy()
    new_iscell = np.zeros((len(groups), 2))
    for k, g in enumerate(groups):
        hidden = sorted({int(rows[u]) for u in g} | {j for u in g for j in leaves[u]})
        for j in hidden:
            stat[j]["inmerge"] = first_new + k
            iscell[j, 0] = 0
        new_iscell[k] = [1, iscell[[rows[u] for u in g], 1].mean()]

    # ---- write ---------------------------------------------------------------
    outputs = {"stat.npy": stat_all, "F.npy": np.vstack([F, F_new]).astype(np.float32),
               "Fneu.npy": np.vstack([Fneu, Fneu_new]).astype(np.float32),
               "spks.npy": np.vstack([spks, spks_new]).astype(np.float32),
               "iscell.npy": np.vstack([iscell, new_iscell])}
    for fname, arr in outputs.items():
        np.save(os.path.join(live_dir, fname), arr)
    export_plane(live_dir)

    # ---- report: before / after ------------------------------------------------
    if groups:
        merged_dff, _ = delta_f_over_f(F_new, Fneu_new, fs)
        plot_traces(os.path.join(report_dir, "merge_traces.png"), groups, best_dff,
                    [int(rows[u]) for u in best_u], merged_dff, fs, estimated=False)
        snr_b, noise_b = event_snr(np.array(best_dff))
        snr_a, noise_a = event_snr(merged_dff)
        print("\ngroup  merged row  segments  npix   noise best -> merged   event SNR best -> merged")
        for k, g in enumerate(groups):
            print(f"{k + 1:5d}  {first_new + k:10d}  {len(g):8d}  {len(new_stats[k]['ypix']):5d}   "
                  f"{noise_b[k]:.3f} -> {noise_a[k]:.3f}        {snr_b[k]:5.1f} -> {snr_a[k]:5.1f}")
        print(f"median event SNR best segment {np.median(snr_b):.1f} -> merged {np.median(snr_a):.1f}")
    log = {"rec_dir": rec_dir, "plane": plane, "n_frames": n_frames,
           "params": {k: globals()[k] for k in PARAMS},
           "steps": [{"a": [int(rows[u]) for u in s["a"]], "b": [int(rows[u]) for u in s["b"]], "r": s["r"], "gap_px": s["gap_px"]} for s in steps],
           "groups": [{"merged_row": first_new + k, "imerge": [int(i) for i in new_stats[k]["imerge"]], "npix": int(len(new_stats[k]["ypix"]))}
                      for k in range(len(groups))]}
    with open(os.path.join(report_dir, "merge_log.json"), "w") as fh:
        json.dump(log, fh, indent=2)

    print(f"\nDone. {len(groups)} merged ROIs appended to {live_dir}  ({len(stat_all)} ROIs, {int(outputs['iscell.npy'][:, 0].sum())} cells)")
    print(f"Unmerged result: {unmerged_dir}")
    print(f"Report:          {report_dir}")
    if SHOW:
        plt.show()


if __name__ == "__main__":
    main()
