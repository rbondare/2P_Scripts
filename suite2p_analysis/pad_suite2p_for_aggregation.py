#!/usr/bin/env python3
"""
Pad a PARTIAL suite2p run with NaN so that the aggregation accepts it.

When suite2p was run on only the first K TIFFs of a recording (LAST_TIF in run_suite2p_RB_axons_V2.py, e.g. because
the water in the objective ran out), F / Fneu / spks have fewer columns than the session has frames.
newAggregation/save_preprocessed.m cuts them with F(:, 1:<frames of ALL tifs>) and crashes on a shorter F.
This script rebuilds <recording>/suite2p/plane0/Fall.mat (the only file the aggregation reads) with F, Fneu and spks
padded with NaN at the end up to the full frame count. Column i of the padded matrices is frame i of the session.

What it does and does not touch
-------------------------------
* The suite2p .npy files (F, Fneu, spks, stat, iscell, ops, data.bin ...) are NOT modified, so the suite2p GUI keeps
  working on the unpadded data. Fall.mat is always rebuilt from them, so run this AFTER curating iscell in the GUI
  (the GUI writes iscell.npy) and again after any re-curation.
* An existing Fall.mat is kept once as Fall_unpadded_original.mat (a later run overwrites Fall.mat only).
* Writes plane0/padding_info.json (number of valid frames, total frames, last TIFF) for the analysis.

Checks (abort unless --force): the processed TIFFs (suite2p db.npy file_list) are exactly the first K TIFFs of the
recording, and their frame count equals the number of columns of F.

Effect downstream: the dF/F baseline (prctile) ignores NaN; the dF/F window median (movmedian) does not, so dF/F is
NaN for the ~15 s (window 15 s) before the cut, and from the cut on.

Usage (suite2p_venv python):
    python pad_suite2p_for_aggregation.py "Z:\\Group Members\\Rima\\DATA_2P\\AnimalRB24_261112_1408"
    options: --plane plane0   --iscell npy|mat (default npy: iscell.npy, the GUI's file)   --dry-run   --force
             --total-frames N (skip counting the TIFFs, tests only)
"""

import argparse
import datetime
import json
import os
import re
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import tifffile
from scipy.io import loadmat


def tif_number(name):
    m = re.search(r'_(\d+)\.tiff?$', name, re.IGNORECASE)
    return int(m.group(1)) if m else None


def count_frames(path):
    with tifffile.TiffFile(path) as tf:
        return len(tf.pages)


def clean_none(o):
    """suite2p's save_mat only replaces top-level None values; MATLAB cannot store nested None."""
    if isinstance(o, dict):
        return {k: clean_none(v) for k, v in o.items()}
    if o is None:
        return np.array([])
    if isinstance(o, (list, tuple)) and any(x is None for x in o):
        return [np.array([]) if x is None else clean_none(x) for x in o]
    return o


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('recording', help='recording folder (contains the TIFFs and suite2p/)')
    ap.add_argument('--plane', default='plane0')
    ap.add_argument('--iscell', choices=('npy', 'mat'), default='npy',
                    help="take the accept/reject flags from iscell.npy (default) or from the existing Fall.mat")
    ap.add_argument('--dry-run', action='store_true', help='do all checks, write nothing')
    ap.add_argument('--force', action='store_true', help='continue although a check fails')
    ap.add_argument('--total-frames', type=int, default=None, help='full frame count (skips counting the TIFFs; tests)')
    args = ap.parse_args()

    rec = args.recording
    s2p = os.path.join(rec, 'suite2p')
    plane = os.path.join(s2p, args.plane)
    if not os.path.isdir(plane):
        sys.exit(f"no suite2p plane folder: {plane}")

    def problem(msg):
        if args.force:
            print(f"  WARNING (continuing, --force): {msg}")
        else:
            sys.exit(f"ABORT: {msg}  (--force to ignore)")

    # ---- suite2p results (unpadded) -------------------------------------------------------------------------
    F = np.load(os.path.join(plane, 'F.npy'))
    Fneu = np.load(os.path.join(plane, 'Fneu.npy'))
    spks = np.load(os.path.join(plane, 'spks.npy'))
    stat = np.load(os.path.join(plane, 'stat.npy'), allow_pickle=True)
    ops = np.load(os.path.join(plane, 'ops.npy'), allow_pickle=True).item()
    n_roi, n_done = F.shape
    print(f"suite2p {args.plane}: {n_roi} ROIs x {n_done} frames (F), Fneu {Fneu.shape}, spks {spks.shape}, stat {len(stat)}")
    if not (Fneu.shape == F.shape == spks.shape and len(stat) == n_roi):
        problem("F / Fneu / spks / stat have different shapes")

    iscell_npy = np.load(os.path.join(plane, 'iscell.npy'))
    fall_path = os.path.join(plane, 'Fall.mat')
    iscell = iscell_npy
    if os.path.exists(fall_path):
        try:
            iscell_mat = loadmat(fall_path, variable_names=['iscell'])['iscell']
            if iscell_mat.shape == iscell_npy.shape:
                n_diff = int((iscell_mat[:, 0] != iscell_npy[:, 0]).sum())
                print(f"existing Fall.mat vs iscell.npy: accepted {int(iscell_mat[:, 0].sum())} vs "
                      f"{int(iscell_npy[:, 0].sum())}, {n_diff} ROIs differ -> using iscell.{args.iscell}")
                if n_diff and args.iscell == 'npy':
                    print("  WARNING: they disagree. iscell.npy is what the suite2p GUI saves on every click; if your "
                          "curation lives in the OLD Fall.mat instead, rerun with --iscell mat.")
                if args.iscell == 'mat':
                    iscell = iscell_mat
            else:
                print(f"existing Fall.mat has iscell {iscell_mat.shape}, iscell.npy {iscell_npy.shape} (different run?)")
                if args.iscell == 'mat':
                    problem("cannot take iscell from Fall.mat: shapes differ")
        except Exception as e:  # truncated / unreadable Fall.mat (seen before)
            print(f"existing Fall.mat unreadable ({e}); using iscell.npy")
            if args.iscell == 'mat':
                problem("cannot take iscell from an unreadable Fall.mat")
    elif args.iscell == 'mat':
        problem("no Fall.mat to take iscell from")
    print(f"accepted ROIs: {int(iscell[:, 0].sum())} of {n_roi}")

    # ---- which TIFFs were processed, how long is the session ------------------------------------------------
    all_tifs = sorted((f for f in os.listdir(rec) if f.lower().endswith(('.tif', '.tiff'))), key=tif_number)
    numbers = [tif_number(f) for f in all_tifs]
    if None in numbers:
        problem("TIFF names do not end in _<number>.tif")
    db = None
    for p in (os.path.join(plane, 'db.npy'), os.path.join(s2p, 'db.npy')):
        if os.path.exists(p):
            db = np.load(p, allow_pickle=True).item()
            break
    processed = [os.path.basename(f) for f in db['file_list']] if db is not None and db.get('file_list') else None
    last_tif = None
    if processed is None:
        print("  no db.npy / file_list found - cannot check which TIFFs were processed")
        problem("processed TIFF list unknown")
    else:
        if processed != all_tifs[:len(processed)]:
            problem(f"processed TIFFs are not the first {len(processed)} TIFFs of the recording")
        last_tif = tif_number(processed[-1])
        print(f"processed TIFFs: {len(processed)} of {len(all_tifs)} (1..{last_tif}), not processed: {len(all_tifs) - len(processed)}")

    if args.total_frames is not None:
        n_total = args.total_frames
        print(f"total frames given on the command line: {n_total}")
    else:
        print(f"counting frames of {len(all_tifs)} TIFFs ...")
        with ThreadPoolExecutor(4) as ex:
            counts = list(ex.map(count_frames, [os.path.join(rec, f) for f in all_tifs]))
        n_total = int(sum(counts))
        if processed is not None:
            n_proc = int(sum(counts[:len(processed)]))
            print(f"  frames in processed TIFFs: {n_proc}, columns of F: {n_done}")
            if n_proc != n_done:
                problem(f"frame count of the processed TIFFs ({n_proc}) != columns of F ({n_done})")
        print(f"  frames in the whole session: {n_total} ({n_total / 30.0 / 60:.1f} min at 30 Hz)")
    if n_total < n_done:
        problem(f"session has fewer frames ({n_total}) than F has columns ({n_done})")
    if n_total == n_done:
        print("F already covers the whole session - nothing to pad.")
        return
    print(f"padding {n_total - n_done} frames ({(n_total - n_done) / n_total * 100:.0f}% of the session) with NaN")

    if args.dry_run:
        print("dry run: nothing written")
        return

    # ---- build the padded Fall.mat --------------------------------------------------------------------------
    def pad(a):
        out = np.full((a.shape[0], n_total), np.nan, np.float32)
        out[:, :n_done] = a
        return out

    Fp, Fneup, spksp = pad(F), pad(Fneu), pad(spks)
    rc_path = os.path.join(plane, 'redcell.npy')
    redcell = np.load(rc_path) if os.path.exists(rc_path) else np.zeros((n_roi, 2), np.float32)

    ops = clean_none(ops)
    ops['pad_n_valid_frames'] = n_done
    ops['pad_n_total_frames'] = n_total
    tmp = os.path.join(plane, '_pad_tmp')
    os.makedirs(tmp, exist_ok=True)
    ops['save_path'] = tmp
    from suite2p.io.save import save_mat
    save_mat(ops, stat, Fp, Fneup, spksp, iscell, redcell)
    new_mat = os.path.join(tmp, 'Fall.mat')

    # verify what was written before it replaces anything
    chk = loadmat(new_mat, variable_names=['F', 'Fneu', 'spks', 'iscell'])
    for name, ref in (('F', F), ('Fneu', Fneu), ('spks', spks)):
        a = chk[name]
        ok = (a.shape == (n_roi, n_total) and np.array_equal(a[:, :n_done], ref.astype(a.dtype), equal_nan=True)
              and bool(np.isnan(a[:, n_done:]).all()))
        print(f"  verify {name}: shape {a.shape} {a.dtype}, first {n_done} columns identical, rest NaN -> {'OK' if ok else 'FAILED'}")
        if not ok:
            shutil.rmtree(tmp, ignore_errors=True)
            sys.exit(f"verification of {name} failed - Fall.mat NOT replaced")
    if not np.array_equal(chk['iscell'], iscell.astype(chk['iscell'].dtype)):
        shutil.rmtree(tmp, ignore_errors=True)
        sys.exit("verification of iscell failed - Fall.mat NOT replaced")

    # keep the previous Fall.mat once - but only if it is a genuine unpadded one (suite2p's own / the GUI's export),
    # not the padded file of an earlier run of this script
    backup = os.path.join(plane, 'Fall_unpadded_original.mat')
    if os.path.exists(fall_path) and not os.path.exists(backup):
        try:
            from scipy.io import whosmat
            f_shape = [s for (nm, s, _) in whosmat(fall_path) if nm == 'F'][0]
            is_unpadded = f_shape[1] != n_total
        except Exception:
            is_unpadded = True  # unreadable / truncated: keep it rather than lose it
        if is_unpadded:
            shutil.move(fall_path, backup)
            print(f"  kept the previous Fall.mat as {os.path.basename(backup)}")
    shutil.move(new_mat, fall_path)
    shutil.rmtree(tmp, ignore_errors=True)
    info = dict(n_valid_frames=int(n_done), n_total_frames=int(n_total), last_tif=last_tif,
                n_tifs_processed=len(processed) if processed else None, n_tifs_total=len(all_tifs),
                n_rois=int(n_roi), n_accepted=int(iscell[:, 0].sum()), iscell_source=args.iscell,
                valid_seconds_at_30hz=round(n_done / 30.0, 1),
                created=datetime.datetime.now().isoformat(timespec='seconds'))
    with open(os.path.join(plane, 'padding_info.json'), 'w') as fh:
        json.dump(info, fh, indent=2)
    print(f"wrote {fall_path} ({os.path.getsize(fall_path) / 1e6:.0f} MB) and padding_info.json")


if __name__ == '__main__':
    main()
