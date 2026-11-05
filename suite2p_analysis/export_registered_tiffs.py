#!/usr/bin/env python3
"""
Export the suite2p-registered movie as one motion-corrected TIFF per raw TIFF
(for manual ROI selection in Fiji instead of the suite2p ROIs).

suite2p registers all TIFFs of a recording against one reference and stores the result in
<recording>/suite2p/plane0/data.bin (int16, Ly x Lx, frames in sorted TIFF order). This script cuts that
binary back into the original TIFF boundaries, so
    AnimalRB22_261028_1631_00002.tif  ->  motion_corrected/AnimalRB22_261028_1631_00002_registered.tif
has the same frame count and frame size as the raw file, and ROI coordinates carry over 1:1. Because every
file was registered to the same reference, ROIs also line up across files.

Nothing in suite2p/ is touched and registration is not recomputed (the 41,805-frame RB22_261028_1631 registration
took hours and would be identical). The per-TIFF frame counts come from the raw TIFFs and must add up to the number of
frames in data.bin, otherwise the script stops.

Frames are int16 like the raw ScanImage files (Fiji shows them with a -32768 calibration function, same as the raw
TIFFs). Edges of each frame contain padding from the registration shifts (suite2p's valid region is in
ops.npy xrange/yrange; for RB22_261028_1631 x 16..1008, y 14..498). Open in Fiji with File > Import > Bio-Formats
(tick "Use virtual stack" for the big files).

Each file is written under a .part name and renamed when complete and verified, and finished files are skipped, so an
interrupted export can simply be run again.

Usage (suite2p_venv; only numpy + tifffile are needed):
    python export_registered_tiffs.py "Z:\\Group Members\\Rima\\DATA_2P\\AnimalRB22_261028_1631" [--out <folder>]
"""

import argparse
import glob
import os
import time

import numpy as np
import tifffile

CHUNK_FRAMES = 200      # frames loaded from data.bin per write (~100 MB per 100 frames at 512 x 1024)
N_CHECK = 3             # frames compared against data.bin after writing each file (first, middle, last)


def tiff_info(path):
    with tifffile.TiffFile(path) as tf:
        return len(tf.pages), tf.pages[0].shape, tf.pages[0].dtype


def verify(path, movie, first, n_frames):
    """Page count and a few frames of the written file must equal the data.bin slice."""
    with tifffile.TiffFile(path) as tf:
        if len(tf.pages) != n_frames:
            raise RuntimeError(f"{os.path.basename(path)}: {len(tf.pages)} pages written, expected {n_frames}")
        for k in sorted({0, n_frames // 2, n_frames - 1}):
            if not np.array_equal(tf.pages[k].asarray(), movie[first + k]):
                raise RuntimeError(f"{os.path.basename(path)}: frame {k} differs from data.bin frame {first + k}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("recording", help="recording folder with the raw TIFFs and suite2p/plane0/data.bin")
    ap.add_argument("--out", help="output folder (default: <recording>/motion_corrected)")
    args = ap.parse_args()

    rec = args.recording
    plane_dir = os.path.join(rec, "suite2p", "plane0")
    bin_path = os.path.join(plane_dir, "data.bin")
    out_dir = args.out or os.path.join(rec, "motion_corrected")

    # suite2p reads the TIFFs of the folder in sorted order (zero-padded numbering), not recursively
    raw = sorted(glob.glob(os.path.join(rec, "*.tif")))
    if not raw or not os.path.exists(bin_path):
        raise SystemExit(f"need raw *.tif in {rec} and {bin_path}")

    print(f"Counting frames of {len(raw)} raw TIFFs...")
    counts, shape, dtype = [], None, None
    for f in raw:
        n, shp, dt = tiff_info(f)
        if shape is not None and (shp, dt) != (shape, dtype):
            raise SystemExit(f"{os.path.basename(f)}: frame {shp} {dt} differs from the other files {shape} {dtype}")
        counts.append(n)
        shape, dtype = shp, dt
        print(f"  {os.path.basename(f)}: {n} frames")
    Ly, Lx = shape
    frame_bytes = 2 * Ly * Lx
    if dtype != np.int16:
        raise SystemExit(f"raw TIFFs are {dtype}, data.bin is int16 — not handled")

    n_bin = os.path.getsize(bin_path) / frame_bytes
    if n_bin != sum(counts):
        raise SystemExit(f"data.bin holds {n_bin:g} frames of {Ly}x{Lx}, the raw TIFFs add up to {sum(counts)} — "
                         f"the registration was not run on exactly these files, nothing written")
    movie = np.memmap(bin_path, dtype=np.int16, mode="r", shape=(sum(counts), Ly, Lx))
    edges = np.concatenate([[0], np.cumsum(counts)])

    os.makedirs(out_dir, exist_ok=True)
    print(f"\nWriting {len(raw)} registered TIFFs ({sum(counts) * frame_bytes / 1e9:.1f} GB) to {out_dir}")
    log = []
    for i, f in enumerate(raw):
        first, stop = int(edges[i]), int(edges[i + 1])
        name = os.path.splitext(os.path.basename(f))[0] + "_registered.tif"
        path = os.path.join(out_dir, name)
        log.append(f"{name}\t{os.path.basename(f)}\tdata.bin frames {first}-{stop - 1}\t{stop - first} frames")
        if os.path.exists(path):
            print(f"  {name} already there")
            continue
        t0 = time.time()
        with tifffile.TiffWriter(path + ".part", bigtiff=True) as tw:
            for s in range(first, stop, CHUNK_FRAMES):
                tw.write(np.asarray(movie[s:min(s + CHUNK_FRAMES, stop)]), contiguous=True)
        verify(path + ".part", movie, first, stop - first)
        os.replace(path + ".part", path)
        gb = (stop - first) * frame_bytes / 1e9
        print(f"  {name}: {stop - first} frames, {gb:.1f} GB, {gb * 1e3 / (time.time() - t0):.0f} MB/s", flush=True)

    with open(os.path.join(out_dir, "motion_corrected_info.txt"), "w") as fh:
        fh.write(f"Registered with suite2p ({plane_dir}), one file per raw TIFF, same frame count and size ({Ly}x{Lx}, int16).\n")
        fh.write("registered file\traw file\tsource in data.bin\tframes\n" + "\n".join(log) + "\n")
    print("\nDone.")


if __name__ == "__main__":
    main()
