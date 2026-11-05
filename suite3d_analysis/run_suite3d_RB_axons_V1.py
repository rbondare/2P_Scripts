#!/usr/bin/env python3
"""
Suite3D Axon Pipeline — V1
3D (z-aligned) alternative to run_suite2p_RB_axons_V2.py for multi-plane axon recordings

Why
---
In a single-plane recording an axon that moves in z changes brightness although nothing happened to the calcium.
Suite3D (Haydaroglu et al. 2025, https://pubmed.ncbi.nlm.nih.gov/40236252/) registers every VOLUME to a reference
volume in 3D (so the axial shift of each volume is measured) and detects/extracts ROIs in 3D, so a bouton that slides
between neighbouring planes is still captured as ONE ROI summed over z.

Run (conda env "suite3d", see environment_suite3d.yml and run_suite3d.bat):
    python run_suite3d_RB_axons_V1.py <recording_folder>                 # whole pipeline
    python run_suite3d_RB_axons_V1.py <recording_folder> --dry-run       # only read headers, print the plan
    python run_suite3d_RB_axons_V1.py <recording_folder> --from-stage corrmap --intensity-thresh 1.0
                                                                         # redo detection with another threshold
(on a second run of the same recording pass the same --job-id, and the same --out-dir if you changed it)

Everything about the acquisition (planes, flyback, volume rate, z step, pixel size) is read from the ScanImage
header of the first tif; override with --n-planes / --n-ch-tif / --fs / --z-step-um / --pixel-um if needed.

What is carried over from run_suite2p_RB_axons_V2.py (nothing about the axons changes)
--------------------------------------------------------------------------------------
  tau 1.0                 same (base ops of the working AnimalRB6_250302_1656 axon run)
  neuropil coefficient    0.7 -> npil_coeff 0.7
  nonrigid registration   128 x 128 px blocks, max 5 px -> block_size_3d (z, 128, 128), max_shift_nr (0, 5, 5)
                          (z range 0, see "z alignment" below: non-rigid z shifts are applied regardless of
                          --apply-z-shift and are noise on a 5-plane / 3 um stack)
  snr_thresh              0.8 (V2 value, Suite3D default is 1.2)
  max rigid shift         0.1 of the frame (V2 maxregshift) -> max_rigid_shift_pix
  fs                      V2 used the per-plane frame rate (30 Hz); Suite3D needs the VOLUME rate
                          (SI.hRoiManager.scanVolumeRate), read from the header
What has no counterpart in Suite3D and is therefore NOT carried over: diameter / sourcery / lam_percentile /
threshold_scaling / the suite2p classifier (suite2p-specific), min_neuropil_pixels (Suite3D builds its neuropil
masks with a fixed 1000-voxel minimum). Suite3D's own detection settings are used instead (see PARAMS below);
for boutons/axons its documentation says intensity_thresh must be lowered from 5 to ~1.5 (the default gives an
empty correlation map), so that is the one value to sweep (--intensity-thresh, --from-stage corrmap).

z alignment on a short stack (read this)
----------------------------------------
* By default the rigid z shift is MEASURED (3D phase correlation, +-2 planes = +-6 um at 3 um spacing; the hard limit
  is nz/2) but NOT APPLIED; --apply-z-shift applies it. Reasons: the shift is applied with trilinear interpolation
  and ZERO fill (shifting a 5-plane volume by one plane leaves an empty plane at one end, so ROIs in the outer planes
  dim), and Suite3D's own notes say the z estimate saturates / is unreliable on stacks of <= 4 planes.
* Tested on a synthetic 5-plane / 3 um bouton recording with known motion (axial SD 0.4 plane, lateral jitter +-4 px):
  x/y shifts were recovered (r 0.91-0.97 once the signal was high enough) but the z shift was NOT (r 0.12 at 4x and
  0.44 at 8x the baseline signal, estimate shrunk towards zero), and applying it reduced the number of detected
  boutons (synthetic_check/make_synth.py + eval_synth.py reproduce this). So on this kind of stack Suite3D's gain is the 3D ROI (a bouton sliding between neighbouring planes is
  one ROI summed over z), not an aligned volume. For a real z-aligned movie of a 5-plane axon stack see
  suite2p-in-depth (Baruchin et al. 2026, volume-based z-registration, designed for 5 planes 2-3 um apart).
  Check the same on your data: the measured z trace should look like real axial motion (steps at running onset),
  not like noise around 0.
* First real 5-plane recording (RB24_261104_1651, ScanImage: 5 planes + 1 flyback = LAST frame of each volume, 3 um,
  5 Hz volumes): the five planes are almost the same optical section (mean images correlate 0.94-0.97, bright
  objects change ~8-12% per plane step, monotonically), so there is no depth information to align on. Suite3D's rigid z
  was exactly 0 between lateral-motion bouts and flipped between +-2 planes during them (r = -0.05 with an independent
  estimate), and with max_shift_nr z = 1 its NON-RIGID z shifts (applied even with --apply-z-shift off) were non-zero in 99%
  of blocks (mean -1.3 planes) and turned planes 0-2 of the registered movie into near copies of each other. Hence z = 0
  in max_shift_nr here. Analysis scripts and details: motion_vs_signal/README.md.
* The measured z-shift trace is exported (Fall_suite3d.mat -> ops.shifts_zyx) and correlated with every ROI trace
  (motion_check.png / motion_check.csv). Axial motion correlates with running, so a correlation does not prove an
  artefact, but ROIs that follow the z trace with no other explanation deserve suspicion. The saturation fraction
  (volumes stuck at the +-2.5 plane window edge) is printed in registration_qc.png/.csv.
* Suite3D's documented fallback for short stacks, 3d_reg=False, crashes in this build (fuse_and_pad_gpu), so it is
  not offered.

Tifs: every tif is treated as its own acquisition (starts at plane 0). Tifs shorter than --min-volumes (10) are skipped and a
trailing partial volume is trimmed; nothing spills into the next tif (Suite3D's default would phase-shift all later
volumes after one stub tif). The skipped/trimmed files are listed in acquisition.json in the job folder.

Outputs (job folder <out-dir>\\s3d-<job-id>\\, results in s3d-results-<job-id>\\):
    rois\\F.npy Fneu.npy spks.npy stats.npy iscell.npy     raw Suite3D output (ROI x volume)
    s3d-results-<job-id>\\Fall_suite3d.mat                 F, Fneu, spks, iscell, stat{} (z/y/x pixels), ops (fs,
                                                           shifts_zyx, tif names + volumes per tif, mean images)
    s3d-results-<job-id>\\registration_qc.png/.csv, motion_check.png/.csv, viewer.html (ROI browser)
"""

import argparse
import ast
import csv
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

# Z: is mapped at different levels of the share depending on the session
Z_ROOTS = ("Z:\\group\\joeschgrp\\", "Z:\\joeschgrp\\", "Z:\\")

STAGES = ["init", "register", "corrmap", "segment", "extract", "export"]
MIN_VOLUMES = 10          # tifs with fewer volumes are rig stubs (RB24_261104_1651 has two 1-volume tifs)

# ============================================================================
# PARAMETERS
# ============================================================================
TAU = 1.0                  # s, GCaMP decay — same as the suite2p V2 run
NPIL_COEFF = 0.7           # neuropil coefficient — same as the suite2p V2 run
SNR_THRESH = 0.8           # nonrigid registration SNR gate — suite2p V2 value (Suite3D default 1.2)
MAXREGSHIFT_FRAC = 0.1     # suite2p V2 maxregshift, as a fraction of the frame
INTENSITY_THRESH = 1.5     # corrmap gate in noise-sigma; Suite3D default 5 is for bright somata and gives an empty
                           # map on boutons/axons (Suite3D debugging_tips.md: "high-zoom boutons sit around 1.5")


def resolve_runtime_path(path_str):
    """Resolve a Windows path whose Z: root may sit at a different level of the share."""
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


# ============================================================================
# ACQUISITION GEOMETRY (ScanImage header)
# ============================================================================

def find_tifs(folder):
    """Top-level tifs of the recording folder (not reg_tif/ or other subfolders), in acquisition order."""
    return sorted(str(p) for p in Path(folder).glob("*.tif") if "_header" not in p.name)


def read_acquisition(tif_path, args):
    """Read planes / flyback / rates / geometry from the ScanImage header of one tif; CLI overrides win."""
    import tifffile

    with tifffile.TiffFile(tif_path) as tf:
        md = tf.scanimage_metadata or {}
        fd = md.get("FrameData", {})
        ny, nx = tf.pages[0].shape[-2:]
        n_pages = len(tf.pages)

    def si(key, default=None):
        return fd.get("SI." + key, default)

    acq = {"tif": tif_path, "n_pages_first_tif": int(n_pages), "Ly": int(ny), "Lx": int(nx)}

    frame_rate = si("hRoiManager.scanFrameRate")
    vol_rate = si("hRoiManager.scanVolumeRate")
    n_slices = si("hStackManager.numSlices")
    n_flyback = si("hFastZ.numDiscardFlybackFrames")
    z_step = si("hStackManager.stackZStepSize")
    channels = np.atleast_1d(si("hChannels.channelSave", 1))
    fov = si("hRoiManager.imagingFovUm")
    acq["header"] = {
        "scanFrameRate": frame_rate, "scanVolumeRate": vol_rate, "numSlices": n_slices,
        "numDiscardFlybackFrames": n_flyback, "stackZStepSize": z_step,
        "numFramesPerVolumeWithFlyback": si("hStackManager.numFramesPerVolumeWithFlyback"),
        "channelSave": channels.tolist(), "scanZoomFactor": si("hRoiManager.scanZoomFactor"),
    }

    # frames per volume in the tif = frame rate / volume rate (planes + flyback frames)
    n_ch_tif = None
    if frame_rate and vol_rate:
        n_ch_tif = int(round(frame_rate / vol_rate))
    n_planes = int(n_slices) if n_slices else None
    pixel_um = None
    if fov is not None:
        fov = np.asarray(fov, float)
        pixel_um = (float(np.ptp(fov[:, 1]) / ny), float(np.ptp(fov[:, 0]) / nx))   # (y, x) um per px

    acq["n_ch_tif"] = args.n_ch_tif or n_ch_tif
    acq["n_planes"] = args.n_planes or n_planes
    acq["vol_rate"] = args.fs or (float(vol_rate) if vol_rate else None)
    acq["z_step_um"] = args.z_step_um or (abs(float(z_step)) if z_step else None)
    acq["pixel_um"] = (args.pixel_um, args.pixel_um) if args.pixel_um else pixel_um
    acq["n_colors"] = int(channels.size)

    missing = [k for k in ("n_ch_tif", "n_planes", "vol_rate", "z_step_um", "pixel_um") if not acq[k]]
    if missing:
        raise SystemExit(f"Could not read {missing} from the ScanImage header of {tif_path}; "
                         "give them with --n-ch-tif / --n-planes / --fs / --z-step-um / --pixel-um.")
    if acq["n_planes"] > acq["n_ch_tif"]:
        raise SystemExit(f"n_planes ({acq['n_planes']}) > frames per volume in the tif ({acq['n_ch_tif']})")
    return acq


def count_pages(path):
    """Exact number of frames in a ScanImage tif without walking every page over the network (~8 ms/page on Z:).

    All pages have the same size and sit back to back, so the count follows from the file size and the distance
    between the first two IFDs; it is confirmed by checking that the last IFD points nowhere and the one before it
    points at it. Falls back to tifffile's full walk if that check fails (checked equal on tifs of 3 to 9447 pages)."""
    import struct
    import tifffile

    size = os.path.getsize(path)
    with open(path, "rb") as f:
        head = f.read(16)
        bo = "<" if head[:2] == b"II" else ">"
        if struct.unpack(bo + "H", head[2:4])[0] == 43:        # BigTIFF
            ptr_fmt, ptr_sz, cnt_fmt, cnt_sz, tag_sz = "Q", 8, "Q", 8, 20
            ifd0 = struct.unpack(bo + "Q", head[8:16])[0]
        else:
            ptr_fmt, ptr_sz, cnt_fmt, cnt_sz, tag_sz = "I", 4, "H", 2, 12
            ifd0 = struct.unpack(bo + "I", head[4:8])[0]

        def next_ifd(off):
            f.seek(off)
            n_tags = struct.unpack(bo + cnt_fmt, f.read(cnt_sz))[0]
            f.seek(off + cnt_sz + n_tags * tag_sz)
            return struct.unpack(bo + ptr_fmt, f.read(ptr_sz))[0]

        try:
            ifd1 = next_ifd(ifd0)
            if ifd1 == 0:
                return 1
            stride = ifd1 - ifd0
            n_est = (size - ifd0) // stride
            for n in (n_est, n_est - 1, n_est + 1):
                if n >= 2 and next_ifd(ifd0 + (n - 1) * stride) == 0 \
                        and next_ifd(ifd0 + (n - 2) * stride) == ifd0 + (n - 1) * stride:
                    return int(n)
        except (struct.error, OSError):
            pass
    with tifffile.TiffFile(path) as tf:
        return len(tf.pages)


def count_volumes(tifs, n_ch_tif):
    """Exact page count of every tif -> whole volumes, trailing partial-volume frames, skipped stubs."""
    rows = []
    for i, p in enumerate(tifs):
        try:
            n, err = count_pages(p), None
        except Exception as e:      # unreadable / still being written
            n, err = 0, repr(e)
        rows.append({"tif": p, "frames": int(n), "volumes": int(n // n_ch_tif),
                     "trailing_frames": int(n % n_ch_tif), "error": err})
        if (i + 1) % 25 == 0:
            print(f"  counted {i + 1}/{len(tifs)} tifs", flush=True)
    return rows


# ============================================================================
# PARAMS
# ============================================================================

def build_params(acq, args):
    nz = acq["n_planes"]
    dz = acq["z_step_um"]
    py, px = acq["pixel_um"]
    if nz < 3:
        raise SystemExit(f"{nz} planes: Suite3D needs >= 3 planes for 3D registration; use run_suite2p_RB_axons_V2.py.")

    params = {
        # --- data: what the acquisition forces ---------------------------------------------------------------
        "fs": acq["vol_rate"],                       # VOLUME rate (not the per-plane frame rate)
        "tau": TAU,
        "voxel_size_um": (dz, py, px),               # (dz, dy, dx)
        "planes": np.arange(nz),                     # flyback frames come after the planes in each volume
        "n_ch_tif": acq["n_ch_tif"],                 # planes + flyback frames per volume in the tif
        "num_colors": acq["n_colors"],
        "functional_color_channel": args.functional_channel,
        "lbm": False, "faced": False,
        "subtract_crosstalk": False,                 # LBM only
        "fuse_strips": False,                        # mesoscope only
        # --- registration: 3D, GPU -----------------------------------------------------------------------------
        "3d_reg": True, "gpu_reg": True,
        "apply_z_shift": args.apply_z_shift,         # off by default, see the docstring
        "pc_size": np.asarray((min(2, nz // 2), 40, 40)),                  # z search is capped at nz/2 planes
        "max_rigid_shift_pix": int(round(MAXREGSHIFT_FRAC * max(acq["Ly"], acq["Lx"]))),
        "nonrigid": True,
        "block_size_3d": (min(5, nz), 128, 128),     # V2: 128 x 128 px blocks
        # V2: maxregshiftNR 5 px in x/y. z = 0: the NON-RIGID z shifts are applied to the registered movie even with
        # apply_z_shift=False, and on 5 planes @ 3 um (planes ~identical, no depth information) they are noise-driven:
        # on RB24_261104_1651 with z = 1 they were non-zero in 99% of blocks (mean -1.3 planes, range -2..+1) and made
        # planes 0/1/2 of the registered movie 0.9-0.99 correlated at single-pixel level (raw: 0.03)
        "max_shift_nr": (0, 5, 5),
        "nr_npad": (1, 3, 3),
        "snr_thresh": SNR_THRESH,
        # --- correlation map (detection) -------------------------------------------------------------------------
        "intensity_thresh": args.intensity_thresh,
        # z filters scaled to the plane spacing so they span the same number of planes as in the Suite3D
        # standard-2P demos (cell: 1 plane sigma, neuropil: 1.5 planes); the xy filters stay at the Suite3D
        # defaults (cell 1.5 um, neuropil 100 um), which is what Suite3D validated on axonal boutons
        "cell_filt_z_um": dz,
        "npil_filt_z_um": 1.5 * dz,
        # --- extraction ---------------------------------------------------------------------------------------------
        "npil_coeff": NPIL_COEFF,
    }
    return params


def apply_overrides(params, set_args):
    from suite3d.default_params import get_known_param_names

    known = get_known_param_names()
    for item in set_args or []:
        key, _, val = item.partition("=")
        if key not in known:
            raise SystemExit(f"--set {key}: not a Suite3D parameter")
        try:
            params[key] = ast.literal_eval(val)
        except (ValueError, SyntaxError):
            params[key] = val
    return params


# ============================================================================
# JOB
# ============================================================================

def create_job(args, tifs, acq, rows):
    """Create a new Suite3D job with exact per-tif frame bookkeeping (no spill-over between acquisitions)."""
    import suite3d.job as s3d_job
    from suite3d.job import Job

    n_ch = acq["n_ch_tif"]
    exact = {r["tif"]: r["frames"] for r in rows}
    # Job.preregister_tifs() estimates frame counts from file sizes; give it the exact ones instead
    s3d_job.get_frame_counts = lambda tif_list, safe_mode=False: {t: exact[t] for t in tif_list}

    params = apply_overrides(build_params(acq, args), args.set)
    job = Job(str(args.out_dir), args.job_id, tifs=tifs, params=params, create=True,
              overwrite=args.overwrite, verbosity=3)

    # every tif is its own acquisition: nothing carries over to the next tif, a partial last volume is trimmed
    job.params["previous_tif"] = {t: None for t in tifs}
    job.params["extra_frames"] = {t: exact[t] % n_ch for t in tifs}
    job.save_params()
    return job


def load_job(args):
    from suite3d.job import Job

    job = Job(str(args.out_dir), args.job_id, create=False)
    if args.set or args.intensity_thresh_given:
        overrides = apply_overrides({}, args.set)
        if args.intensity_thresh_given:
            overrides["intensity_thresh"] = args.intensity_thresh
        job.save_params(new_params=overrides)
        print(f"Updated params of the existing job: {overrides}")
    return job


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def stage_wanted(stage, from_stage):
    return STAGES.index(stage) >= STAGES.index(from_stage)


def extract_batch(job, target_gb=4.0):
    """Volumes held in RAM per extraction batch (whole multiples of the 100-volume on-disk chunk)."""
    mov = job.get_registered_movie()
    nz, nt, ny, nx = mov.shape
    chunk = mov.chunks[1][0]
    per_frame = nz * ny * nx * mov.dtype.itemsize * 2     # x2 for the shared-memory copy
    n_chunks = int(target_gb * 1024 ** 3 / (per_frame * chunk))
    return min(max(1, n_chunks) * chunk, nt)


# ============================================================================
# QC AND EXPORT
# ============================================================================

def load_shifts(job):
    """Per-volume rigid shifts (z, y, x in voxels) and volumes per registered batch (= per tif)."""
    d = job.dirs["registered_fused_data"]
    files = sorted(f for f in os.listdir(d) if f.startswith("offsets") and f.endswith(".npy"))
    offs = [np.load(os.path.join(d, f), allow_pickle=True).item() for f in files]
    shifts = np.concatenate([o["sub_pixel_shifts"] for o in offs])
    per_batch = [len(o["sub_pixel_shifts"]) for o in offs]
    return shifts, per_batch


def registration_qc(job, out_dir, shifts, per_batch):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fs = job.params["fs"]
    dz = job.params["voxel_size_um"][0]
    pc = np.asarray(job.params["pc_size"])
    nt = len(shifts)
    t = np.arange(nt) / fs
    edges = np.cumsum(per_batch)[:-1] / fs

    cap = pc + 0.5
    sat = [np.isclose(np.abs(shifts[:, ax]), cap[ax]) for ax in range(3)]
    sat_pct = [100.0 * s.mean() for s in sat]
    z = shifts[:, 0]

    with open(os.path.join(out_dir, "registration_qc.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["volume", "time_s", "z_planes", "z_um", "y_px", "x_px", "z_saturated"])
        for i in range(nt):
            w.writerow([i, f"{t[i]:.3f}", f"{z[i]:.3f}", f"{z[i] * dz:.3f}",
                        f"{shifts[i, 1]:.3f}", f"{shifts[i, 2]:.3f}", int(sat[0][i])])

    fig, ax = plt.subplots(3, 1, figsize=(12, 8), constrained_layout=True)
    ax[0].plot(t, z * dz, lw=0.6)
    for lim in (cap[0] * dz, -cap[0] * dz):
        ax[0].axhline(lim, color="r", ls="--", lw=0.8)
    ax[0].set_ylabel("z shift (um)")
    ax[0].set_title(f"rigid z shift per volume ({dz:g} um/plane); red = search-window edge (+-{cap[0]:.1f} planes), "
                    f"{sat_pct[0]:.2f}% of volumes saturated")
    ax[1].plot(t, shifts[:, 1], lw=0.5, label=f"y ({sat_pct[1]:.2f}% saturated)")
    ax[1].plot(t, shifts[:, 2], lw=0.5, label=f"x ({sat_pct[2]:.2f}% saturated)")
    ax[1].set_ylabel("xy shift (px)")
    ax[1].set_xlabel("time (s)")
    ax[1].legend(loc="upper right")
    for a in ax[:2]:
        for e in edges:
            a.axvline(e, color="k", lw=0.3, alpha=0.4)
    ax[2].hist(z, bins=np.arange(-cap[0] - 0.25, cap[0] + 0.5, 0.25))
    ax[2].set_xlabel("z shift (planes)")
    ax[2].set_ylabel("volumes")
    fig.savefig(os.path.join(out_dir, "registration_qc.png"), dpi=130)
    plt.close(fig)

    log(f"z shift: median {np.median(z) * dz:+.2f} um, 5-95% [{np.percentile(z, 5) * dz:+.2f}, "
        f"{np.percentile(z, 95) * dz:+.2f}] um, max |z| {np.abs(z).max() * dz:.2f} um")
    log(f"saturated at the search-window edge: z {sat_pct[0]:.2f}%  y {sat_pct[1]:.2f}%  x {sat_pct[2]:.2f}%")
    if sat_pct[0] >= 1.0:
        log("WARNING: >=1% of volumes hit the z search-window edge. That is the estimator failing (weak "
            "phase-correlation peak), not real motion; do not trust the z trace"
            + (" (and those volumes were shifted by the cap: re-run without --apply-z-shift)."
               if job.params["apply_z_shift"] else "."))
    return sat_pct


def delta_f_over_f(F, Fneu, fs):
    """Neuropil-subtracted dF/F with an 8th-percentile running baseline (30 s), denominator bounded at 1."""
    from scipy.ndimage import percentile_filter

    Fc = F.astype(np.float32) - NPIL_COEFF * Fneu.astype(np.float32)
    win = max(3, int(round(30 * fs)))
    base = percentile_filter(Fc, 8, size=(1, win), mode="nearest")
    return (Fc - base) / np.maximum(base, 1.0)


def motion_check(job, out_dir, shifts):
    """Correlate every ROI's dF/F with the measured z shift (and |xy| shift); circular-shift null."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rois = job.dirs["rois"]
    F, Fneu = np.load(os.path.join(rois, "F.npy")), np.load(os.path.join(rois, "Fneu.npy"))
    nroi, nt = F.shape
    if nt != len(shifts):
        log(f"motion check skipped: {nt} trace samples vs {len(shifts)} registration shifts")
        return
    fs = job.params["fs"]
    dff = delta_f_over_f(F, Fneu, fs)
    ok = np.isfinite(dff).all(axis=1)

    def corr_rows(X, y):
        y = (y - y.mean()) / (y.std() + 1e-12)
        X = (X - X.mean(axis=1, keepdims=True)) / (X.std(axis=1, keepdims=True) + 1e-12)
        return X @ y / X.shape[1]

    rng = np.random.default_rng(0)
    results = {}
    for name, trace in (("z", shifts[:, 0]), ("xy", np.hypot(shifts[:, 1], shifts[:, 2]))):
        if np.std(trace) < 1e-6:
            continue
        r = corr_rows(dff, trace)
        # null: same trace circularly shifted by >= 60 s (keeps its autocorrelation, breaks the alignment)
        min_shift = min(nt // 4, int(60 * fs))
        null = np.stack([corr_rows(dff, np.roll(trace, rng.integers(min_shift, nt - min_shift)))
                         for _ in range(20)])
        thr = np.percentile(np.abs(null[:, ok]), 99)
        results[name] = (r, thr)
        log(f"motion check [{name}]: median |r| {np.nanmedian(np.abs(r[ok])):.3f}, null 99th pct |r| {thr:.3f}, "
            f"{100 * np.mean(np.abs(r[ok]) > thr):.1f}% of {ok.sum()} ROIs beyond the null")

    if not results:
        log("motion check skipped: the registration shifts are constant")
        return
    with open(os.path.join(out_dir, "motion_check.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["roi", "r_dff_vs_z"] + (["r_dff_vs_absxy"] if "xy" in results else []))
        for i in range(nroi):
            w.writerow([i] + [f"{results[k][0][i]:.4f}" for k in ("z", "xy") if k in results])

    fig, axs = plt.subplots(1, len(results), figsize=(5.5 * len(results), 4), constrained_layout=True, squeeze=False)
    for a, (name, (r, thr)) in zip(axs[0], results.items()):
        a.hist(r[ok], bins=60)
        a.axvline(thr, color="r", ls="--")
        a.axvline(-thr, color="r", ls="--")
        a.set_xlabel(f"r(dF/F, {'z shift' if name == 'z' else '|xy shift|'})")
        a.set_ylabel("ROIs")
        a.set_title("red: 99th pct of the circular-shift null")
    fig.savefig(os.path.join(out_dir, "motion_check.png"), dpi=130)
    plt.close(fig)


def export_mat(job, out_dir, shifts, rows, acq):
    import scipy.io

    rois = job.dirs["rois"]
    F, Fneu = np.load(os.path.join(rois, "F.npy")), np.load(os.path.join(rois, "Fneu.npy"))
    spks = np.load(os.path.join(rois, "spks.npy"))
    stats = np.load(os.path.join(rois, "stats.npy"), allow_pickle=True)
    if len(stats) != F.shape[0]:
        raise RuntimeError(f"{len(stats)} stats vs {F.shape[0]} trace rows: ROI/trace pairing is ambiguous")

    stat = np.empty((1, len(stats)), dtype=object)
    for i, s in enumerate(stats):
        zz, yy, xx = (np.asarray(c, dtype=np.int32) for c in s["coords"])
        stat[0, i] = {"zpix": zz, "ypix": yy, "xpix": xx, "lam": np.asarray(s["lam"], dtype=np.float32),
                      "med": np.asarray(s["med"], dtype=np.float32), "npix": len(zz),
                      "n_planes": len(np.unique(zz))}

    ops = {"fs": float(job.params["fs"]), "tau": float(job.params["tau"]),
           "voxel_size_um": np.asarray(job.params["voxel_size_um"], dtype=float),
           "n_planes": int(len(job.params["planes"])), "n_ch_tif": int(job.params["n_ch_tif"]),
           "apply_z_shift": int(job.params["apply_z_shift"]), "shifts_zyx": shifts,
           "tif_names": np.array([os.path.basename(r["tif"]) for r in rows if r["volumes"] > 0], dtype=object),
           "tif_volumes": np.array([r["volumes"] for r in rows if r["volumes"] > 0]),
           "coords_note": "stat zpix/ypix/xpix and med are 0-based (z, y, x), like suite2p; F is ROI x volume"}
    try:
        for k, v in job.load_corr_map_results().items():
            ops[k] = np.asarray(v, dtype=np.float32)
        ops["ref_img_3d"] = np.asarray(job.load_summary()["ref_img_3d"], dtype=np.float32)
    except Exception as e:
        log(f"mean/vmap images not added to the .mat ({e!r})")

    path = os.path.join(out_dir, "Fall_suite3d.mat")
    scipy.io.savemat(path, {"stat": stat, "F": F, "Fneu": Fneu, "spks": spks,
                            "iscell": np.column_stack([np.ones(len(stats)), np.ones(len(stats))]),
                            "ops": ops}, do_compression=True)
    log(f"wrote {path} ({os.path.getsize(path) / 1e6:.1f} MB): {F.shape[0]} ROIs x {F.shape[1]} volumes")


# ============================================================================
# MAIN
# ============================================================================

def parse_args():
    p = argparse.ArgumentParser(description="Suite3D axon pipeline V1 (3D / z-aligned) — see the module docstring.",
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("recording", help="recording folder with the ScanImage tifs, e.g. Z:\\...\\AnimalRB25_261104_1200")
    p.add_argument("--out-dir", default=str(Path.home() / "Suite3D"),
                   help="where the job folder is written (default: %(default)s). Keep it on a LOCAL disk: the "
                        "registered movie is as big as the raw data")
    p.add_argument("--job-id", default=None, help="job name (default: the recording folder name)")
    p.add_argument("--from-stage", choices=STAGES, default=None,
                   help="re-run from this stage on, re-using the earlier stages of an existing job")
    p.add_argument("--overwrite", action="store_true", help="delete and recreate an existing job")
    p.add_argument("--dry-run", action="store_true", help="read headers and count frames, print the plan, stop")
    p.add_argument("--max-tifs", type=int, default=None, help="only use the first N tifs (quick test)")
    p.add_argument("--apply-z-shift", action="store_true",
                   help="also APPLY the measured rigid z shift (default: measure it, register x/y only)")
    p.add_argument("--intensity-thresh", type=float, default=None,
                   help=f"correlation-map gate in noise-sigma (default {INTENSITY_THRESH}); lower = more ROIs, noisier")
    p.add_argument("--min-volumes", type=int, default=MIN_VOLUMES,
                   help="skip tifs with fewer volumes than this (rig stubs; one broke the init pass on RB24_261104_1651)")
    p.add_argument("--functional-channel", type=int, default=0, help="0-based colour channel of the GCaMP signal")
    p.add_argument("--set", action="append", metavar="KEY=VALUE",
                   help="override any Suite3D parameter, e.g. --set peak_thresh=0.05 (repeatable)")
    for flag, typ, what in (("--n-planes", int, "imaged planes per volume"),
                            ("--n-ch-tif", int, "frames per volume in the tif (planes + flyback)"),
                            ("--fs", float, "volume rate (Hz)"), ("--z-step-um", float, "plane spacing (um)"),
                            ("--pixel-um", float, "um per pixel")):
        p.add_argument(flag, type=typ, default=None, help=f"override the header: {what}")
    args = p.parse_args()
    args.intensity_thresh_given = args.intensity_thresh is not None
    if args.intensity_thresh is None:
        args.intensity_thresh = INTENSITY_THRESH
    return args


def main():
    args = parse_args()
    recording = resolve_runtime_path(args.recording)
    if not os.path.isdir(recording):
        raise SystemExit(f"Recording folder not found: {recording}")
    args.job_id = args.job_id or Path(recording).name
    args.out_dir = Path(args.out_dir)
    job_dir = args.out_dir / f"s3d-{args.job_id}"
    existing = (job_dir / "params.npy").exists() and not args.overwrite

    print("\n" + "=" * 70)
    print("SUITE3D AXON PIPELINE V1")
    print("=" * 70)

    if existing and args.from_stage is None and not args.dry_run:
        raise SystemExit(f"Job {job_dir} already exists. Use --from-stage <stage> to continue it, or --overwrite.")
    from_stage = args.from_stage or "init"

    tifs = find_tifs(recording)
    if args.max_tifs:
        tifs = tifs[:args.max_tifs]
    if not tifs:
        raise SystemExit(f"No .tif files in {recording}")
    log(f"{len(tifs)} tifs in {recording}")

    acq = read_acquisition(tifs[0], args)
    print(f"\n  header of {os.path.basename(tifs[0])}:")
    for k, v in acq["header"].items():
        print(f"    {k:32s} {v}")
    py, px = acq["pixel_um"]
    print(f"  -> {acq['n_planes']} planes + {acq['n_ch_tif'] - acq['n_planes']} flyback frame(s) per volume "
          f"({acq['n_ch_tif']} frames), volume rate {acq['vol_rate']:.4f} Hz, dz {acq['z_step_um']:g} um, "
          f"pixel {py:.3f} x {px:.3f} um, frame {acq['Ly']} x {acq['Lx']}, {acq['n_colors']} colour channel(s)")
    n_flyback_header = acq["header"]["numDiscardFlybackFrames"]
    if n_flyback_header is not None and acq["n_ch_tif"] - acq["n_planes"] != n_flyback_header:
        print("  WARNING: frames per volume minus planes != numDiscardFlybackFrames in the header; "
              "check --n-ch-tif / --n-planes (the planes must be the FIRST frames of every volume)")

    rows = count_volumes(tifs, acq["n_ch_tif"])
    for r in rows:      # rig stubs (an acquisition that was started and stopped at once): Suite3D's init pass can pick one
        if r["volumes"] < args.min_volumes:
            r["volumes"], r["trailing_frames"] = 0, 0
    usable = [r for r in rows if r["volumes"] > 0]
    skipped = [r for r in rows if r["volumes"] == 0]
    trimmed = [r for r in usable if r["trailing_frames"]]
    tot_vol = sum(r["volumes"] for r in usable)
    log(f"{len(usable)} usable tifs, {tot_vol} volumes = {tot_vol / acq['vol_rate'] / 60:.1f} min; "
        f"{len(skipped)} tifs shorter than {args.min_volumes} volumes skipped, {len(trimmed)} with a trimmed partial last volume")
    for r in skipped:
        print(f"    skipped {os.path.basename(r['tif'])}: {r['frames']} frames {r['error'] or ''}")
    for r in trimmed:
        print(f"    trimmed {os.path.basename(r['tif'])}: {r['frames']} frames -> {r['trailing_frames']} dropped")
    if not usable:
        raise SystemExit("No tif holds a complete volume")
    usable_tifs = [r["tif"] for r in usable]

    preview = build_params(acq, args)
    print("\n  Suite3D parameters set by this script (everything else = Suite3D defaults):")
    for k, v in preview.items():
        print(f"    {k:26s} {v}")
    if args.set:
        print(f"    --set overrides: {args.set}")
    if args.dry_run:
        print("\nDry run: nothing was written.")
        return

    import cupy as cp
    if cp.cuda.runtime.getDeviceCount() < 1:
        raise SystemExit("No CUDA GPU: Suite3D 3D registration runs on the GPU only.")
    print(f"\n  GPU: {cp.cuda.runtime.getDeviceProperties(0)['name'].decode()}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    if existing:
        job = load_job(args)
    else:
        job = create_job(args, usable_tifs, acq, rows)
        with open(os.path.join(job.job_dir, "acquisition.json"), "w") as fh:
            json.dump({"acquisition": {k: v for k, v in acq.items() if k != "pixel_um"} | {"pixel_um": acq["pixel_um"]},
                       "tifs": rows}, fh, indent=1, default=str)

    t0 = time.time()
    if stage_wanted("init", from_stage):
        log("=== 1/6 initialization (reference volume, plane alignment) ===")
        job.run_init_pass()
    if stage_wanted("register", from_stage):
        log("=== 2/6 3D registration ===")
        job.register()
    if stage_wanted("corrmap", from_stage):
        log("=== 3/6 correlation map ===")
        job.calculate_corr_map()
    if stage_wanted("segment", from_stage):
        log("=== 4/6 segmentation ===")
        job.segment_rois()
    if stage_wanted("extract", from_stage):
        log("=== 5/6 neuropil masks + trace extraction + deconvolution ===")
        job.compute_npil_masks()
        job.extract_and_deconvolve(batchsize_frames=extract_batch(job))

    log("=== 6/6 export + QC ===")
    out_dir = job.export_results(job.job_dir, make_viewer=False)
    shifts, per_batch = load_shifts(job)
    registration_qc(job, out_dir, shifts, per_batch)
    motion_check(job, out_dir, shifts)
    export_mat(job, out_dir, shifts, rows, acq)
    try:
        # Suite3D's export copies the FULL segmentation (stats/iscell) next to the traces of the ROIs that were
        # extracted (a few are dropped in extraction): on RB24_261104_1651 9582 stats vs 9498 trace rows, so the browser
        # would pair ROI i with the wrong trace. job/rois (used for the .mat above) is paired correctly.
        n_stats = len(np.load(os.path.join(out_dir, "stats.npy"), allow_pickle=True))
        n_trace = np.load(os.path.join(out_dir, "F.npy"), mmap_mode="r").shape[0]
        if n_stats != n_trace:
            log(f"HTML viewer skipped: Suite3D's export has {n_stats} ROIs but {n_trace} traces (unpaired); "
                f"use Fall_suite3d.mat")
        else:
            job.make_html_viewer(export_path=out_dir)
    except Exception as e:     # the ROI browser is a convenience; the numbers above are what matters
        log(f"HTML viewer not written ({e!r})")

    print("\n" + "=" * 70)
    print(f"SUITE3D COMPLETE in {(time.time() - t0) / 60:.1f} min")
    print("=" * 70)
    print(f"  job folder : {job.job_dir}")
    print(f"  results    : {out_dir}")
    print("  Fall_suite3d.mat (MATLAB), registration_qc.png, motion_check.png, viewer.html")


if __name__ == "__main__":
    main()
