#!/usr/bin/env python3
"""
Export suite2p .npy outputs to a clean Fall.mat for MATLAB analysis.

Usage:
    python export_suite2p_to_mat.py <plane_dir>

Example:
    python export_suite2p_to_mat.py "Z:\...\AnimalRB6_250302_1711\suite2p\plane0"

The script loads stat.npy, F.npy, Fneu.npy, spks.npy, iscell.npy from the
given plane directory and writes Fall.mat in the same directory.

Why this exists: suite2p's GUI save_mat crashes when ops contains None values
(bin_size=None, classifier_path=None, etc.) and scipy.io.savemat cannot write
Python None to a MATLAB struct. This script bypasses ops entirely and writes
only the arrays needed for downstream MATLAB analysis.
"""

import sys
import os
import numpy as np
import scipy.io


STAT_FIELDS_FOR_MAT = [
    'ypix', 'xpix', 'lam',
    'med', 'medw',
    'npix', 'npix_soma', 'npix_norm',
    'compact', 'solidity', 'aspect_ratio',
    'skew', 'std', 'snr',
    'footprint', 'mrs', 'mrs0',
    'chan2_prob',
    'inmerge',
]


def sanitize_value(v):
    """Convert a stat field value to something scipy.io.savemat can write."""
    if v is None:
        return np.array([], dtype=np.float64)
    if isinstance(v, (list, tuple)):
        arr = np.array(v)
        return arr if arr.dtype != object else np.array([], dtype=np.float64)
    if isinstance(v, np.ndarray):
        return v
    if isinstance(v, (bool, np.bool_)):
        return int(v)
    return v


def stat_to_mat_cell(stat_arr):
    """
    Convert stat dicts to a (1, N) numpy object array.
    scipy.io.savemat writes object arrays as MATLAB cell arrays, giving
    Fall.stat{n}.med access — identical to suite2p's own save_mat output.
    """
    cell = np.empty((1, len(stat_arr)), dtype=object)
    for i, roi in enumerate(stat_arr):
        cell[0, i] = {f: sanitize_value(roi.get(f))
                      for f in STAT_FIELDS_FOR_MAT if f in roi}
    return cell


def load_ops_safe(plane_dir):
    """
    Load ops.npy and return a dict containing only MATLAB-serialisable fields.
    None values, Path objects, and mixed-type lists are silently dropped —
    these are the fields that crash scipy.io.savemat in the original suite2p GUI.
    Falls back to minimal stubs if ops.npy is absent.
    """
    ops_path = os.path.join(plane_dir, 'ops.npy')
    if not os.path.exists(ops_path):
        return {'max_proj': np.zeros((512, 512), dtype=np.float32),
                'meanImg':  np.zeros((512, 512), dtype=np.float32)}
    ops = np.load(ops_path, allow_pickle=True).item()
    safe = {}
    for k, v in ops.items():
        if v is None:
            continue
        if isinstance(v, np.ndarray):
            safe[k] = v
        elif isinstance(v, (int, float, bool, str,
                            np.integer, np.floating, np.bool_)):
            safe[k] = v
        # Path objects, dicts, lists-with-None, etc. are dropped
    return safe


def main():
    if len(sys.argv) < 2:
        print("Usage: python export_suite2p_to_mat.py <plane_dir>")
        print("  plane_dir: path to suite2p planeX folder containing .npy files")
        sys.exit(1)

    plane_dir = sys.argv[1]
    if not os.path.isdir(plane_dir):
        print(f"Error: directory not found: {plane_dir}")
        sys.exit(1)

    def load(fname):
        path = os.path.join(plane_dir, fname)
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing: {path}")
        return np.load(path, allow_pickle=True)

    print(f"Loading from: {plane_dir}")
    stat   = load("stat.npy")
    F      = load("F.npy")
    Fneu   = load("Fneu.npy")
    spks   = load("spks.npy")
    iscell = load("iscell.npy")   # shape (N, 2): [iscell_flag, probability]

    print(f"  {len(stat)} ROIs  |  {F.shape[1]} frames")
    print(f"  cells: {int(iscell[:, 0].sum())}  non-cells: {int((1 - iscell[:, 0]).sum())}")

    # Constituent ROIs (those absorbed into a merge) have inmerge > 0.
    # They are hidden visually in the GUI but their iscell flag is NOT cleared
    # by the merge operation. Zero it out here so MATLAB code that filters by
    # iscell(:,1)==1 only sees the actual merged ROI, not its hidden pieces.
    inmerge = np.array([s.get('inmerge', 0) for s in stat], dtype=float)
    constituent_mask = inmerge > 0
    n_constituents = int(constituent_mask.sum())
    if n_constituents > 0:
        iscell = iscell.copy()          # don't mutate the loaded array
        iscell[constituent_mask, 0] = 0
        print(f"  masked {n_constituents} hidden constituent ROI(s) (inmerge>0) → iscell set to 0 in Fall.mat")
        print(f"  real cells after masking: {int(iscell[:, 0].sum())}")

    print("Converting stat to MATLAB cell array...")
    stat_mat = stat_to_mat_cell(stat)

    print("Loading ops...")
    ops = load_ops_safe(plane_dir)

    redcell_path = os.path.join(plane_dir, 'redcell.npy')
    redcell = (np.load(redcell_path, allow_pickle=True)
               if os.path.exists(redcell_path)
               else np.zeros((len(stat), 2), dtype=np.float32))

    out_path = os.path.join(plane_dir, "Fall.mat")
    print(f"Writing {out_path} ...")
    scipy.io.savemat(out_path, {
        "stat":    stat_mat,
        "F":       F,
        "Fneu":    Fneu,
        "spks":    spks,
        "iscell":  iscell,
        "ops":     ops,
        "redcell": redcell,
    }, do_compression=True)

    size_mb = os.path.getsize(out_path) / 1e6
    print(f"Done — {size_mb:.1f} MB written to Fall.mat")


if __name__ == "__main__":
    main()
