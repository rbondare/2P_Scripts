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


def stat_to_mat_struct(stat_arr):
    """
    Convert a numpy array of stat dicts into a form scipy can write
    as a MATLAB struct array.
    """
    empty = np.array([], dtype=np.float64)

    # Collect fields present in ANY roi (union), filtered to our allow-list
    all_fields = []
    seen = set()
    for roi in stat_arr:
        for field in STAT_FIELDS_FOR_MAT:
            if field in roi and field not in seen:
                all_fields.append(field)
                seen.add(field)

    dtype = [(f, object) for f in all_fields]
    mat_struct = np.zeros(len(stat_arr), dtype=dtype)
    for i, roi in enumerate(stat_arr):
        for field in all_fields:
            mat_struct[i][field] = sanitize_value(roi.get(field, None)) if field in roi else empty
    return mat_struct


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

    print("Converting stat to MATLAB struct...")
    stat_mat = stat_to_mat_struct(stat)

    out_path = os.path.join(plane_dir, "Fall.mat")
    print(f"Writing {out_path} ...")
    scipy.io.savemat(out_path, {
        "stat":   stat_mat,
        "F":      F,
        "Fneu":   Fneu,
        "spks":   spks,
        "iscell": iscell,
    }, do_compression=True)

    size_mb = os.path.getsize(out_path) / 1e6
    print(f"Done — {size_mb:.1f} MB written to Fall.mat")


if __name__ == "__main__":
    main()
