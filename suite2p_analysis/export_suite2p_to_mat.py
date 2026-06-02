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

    Returns a numpy array of cleaned dicts, one per ROI.
    """
    cleaned = []
    for roi in stat_arr:
        d = {}
        for field in STAT_FIELDS_FOR_MAT:
            if field in roi:
                d[field] = sanitize_value(roi[field])
        cleaned.append(d)

    dtype = [(f, object) for f in cleaned[0].keys()]
    mat_struct = np.zeros(len(cleaned), dtype=dtype)
    for i, d in enumerate(cleaned):
        for f, v in d.items():
            mat_struct[i][f] = v
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
