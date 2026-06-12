#!/usr/bin/env python3
"""
Suite2P Axon Detection Script — V2
Optimized for detecting faint, sparse axons as shown in fluorescence imaging

Changes from V1
---------------
1. lam_percentile: 0.0 → 50.0
   V1 set lam_percentile=0 to avoid an OOM that occurs when large axon ROIs
   trigger size=int(radius*5) inside suite2p's percentile_filter. Setting it
   to 0 fixes the OOM but includes all low-weight peripheral pixels in F,
   adding substantial noise. 50 matches the original working AnimalRB6_250302_1656
   ops and keeps only the top-50% brightest pixels per ROI mask. If an OOM
   reappears, try 25 as a middle ground.

2. spatial_taper: 3.45 → 50.0
   V1 set spatial_taper=3.45 (= 3 × smooth_sigma), which is almost no
   Fourier-edge tapering. The old working ops used the default 50 px. With
   3.45 px the phase-correlation registration has edge-ringing artefacts that
   degrade alignment and appear as jitter noise in F. Restored to 50 px.

3. Legacy flat keys from the 0.14.4 ops template are now explicitly dropped
   before the nested sub-dicts are set. Suite2p 1.0+ reads from the nested
   sub-dicts and ignores the flat equivalents, but retaining them caused
   confusion when reading saved ops (e.g. flat lam_percentile=50 coexisting
   with nested lam_percentile=0). Dropping them makes the nested values the
   single source of truth.

4. Combined folder creation removed. V1's build_combined_without_flyback() had
   issues saving the combined output reliably. V2 leaves outputs plane-separated
   (plane0/, plane1/) and does not attempt to merge them.

Everything else (paths, diameter, detection, classification, db) is unchanged.
"""

import numpy as np
import os
import torch
from pathlib import PureWindowsPath
from suite2p import run_s2p


# GPU check — abort early if CUDA is not available so you don't run a long job on CPU by accident
print("\n" + "="*70)
print("GPU CHECK")
print("="*70)
if torch.cuda.is_available():
    print(f"  GPU detected: {torch.cuda.get_device_name(0)}")
    print(f"  VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
    print("  Suite2p will use GPU for registration and detection.")
else:
    raise RuntimeError(
        "CUDA not available — aborting to avoid running on CPU.\n"
        "Fix: python -m pip install torch --index-url https://download.pytorch.org/whl/cu128"
    )
print("="*70 + "\n")


def resolve_runtime_path(path_str):
    """Resolve a Windows-style path on Windows or WSL/Linux."""
    if os.name == "nt":
        return path_str

    windows_path = PureWindowsPath(path_str)
    drive = windows_path.drive
    if drive:
        drive_letter = drive.rstrip(":").lower()
        return os.path.join("/mnt", drive_letter, *windows_path.parts[1:])

    return path_str.replace("\\", "/")


# Load ops from the previously SUCCESSFUL axon detection run (AnimalRB6_250302_1656)
# These settings are known to detect axons well — use as base, only override paths
settings_file = resolve_runtime_path(r"Z:\joeschgrp\Group Members\Rima\DATA_2P_OTHER\AnimalRB6\AnimalRB6_250302_1656\suite2p\combined\ops.npy")
data_path = resolve_runtime_path(r"Z:\joeschgrp\Rotation Students\Dow\DATA_2P\AnimalDG3_260320_1311")
classifier_path = resolve_runtime_path(r"Z:\joeschgrp\Group Members\Rima\DATA_2P_OTHER\AnimalRB6\AnimalRB6_250302_1656\suite2p\plane0\NAaxons.npy")

if os.path.exists(settings_file):
    ops = np.load(settings_file, allow_pickle=True).item()
    # Drop recording-specific keys — must be recomputed for the new recording
    for drop_key in (
        'version',
        'save_path0', 'save_path', 'fast_disk', 'data_path',
        'yrange', 'xrange',           # registration crop — shape differs per recording
        'meanImg', 'meanImgE',        # mean images — recording specific
        'refImg',                     # registration reference frame
        'max_proj',                   # max projection
        'Ly', 'Lx',                   # frame dimensions
        'tPC',                        # temporal PCs
        'ignore_flyback', 'nplanes',  # will be set via db dict
        # --- legacy flat keys now owned by nested sub-dicts ---
        # extraction
        'lam_percentile', 'neuropil_extract',
        'inner_neuropil_radius', 'min_neuropil_pixels', 'allow_overlap',
        # neucoeff is intentionally kept as a flat key — merge.py reads ops["neucoeff"] directly
        # detection
        'soma_crop', 'threshold_scaling', 'max_overlap', 'denoise',
        'connected', 'sparse_mode', 'nbinned', 'max_iterations',
        'roidetect', 'anatomical_only',
        # registration
        'spatial_taper', 'smooth_sigma', 'smooth_sigma_time',
        'nonrigid', 'maxregshiftNR', 'nimg_init', 'maxregshift',
        'norm_frames', 'th_badframes', 'subpixel', 'two_step_registration',
        'bidiphase', 'do_bidiphase', 'snr_thresh', 'reg_tif', 'reg_tif_chan2',
        'align_by_chan',
        # io
        'combined', 'save_mat', 'save_NWB', 'delete_bin', 'move_bin',
        # run
        'do_registration', 'multiplane_parallel', 'spikedetect',
        # classification
        'classifier_path', 'use_builtin_classifier', 'preclassify',
    ):
        ops.pop(drop_key, None)
    print("Loaded base settings from WORKING axon ops (AnimalRB6_250302_1656)")
else:
    print(f"Settings file not found: {settings_file}")
    exit(1)

# Validate custom classifier before doing anything else
print("\nValidating custom axon classifier...")
if not os.path.exists(classifier_path):
    print(f"  ERROR: classifier not found at {classifier_path}")
    exit(1)
_clf = np.load(classifier_path, allow_pickle=True).item()
_required_keys = {'stats', 'iscell', 'keys'}
if not _required_keys.issubset(_clf.keys()):
    print(f"  ERROR: classifier missing keys. Found: {list(_clf.keys())}, need: {_required_keys}")
    exit(1)
_n_rois = _clf['stats'].shape[0]
_n_cells = int(_clf['iscell'].sum())
print(f"  OK — {_n_rois} labeled ROIs ({_n_cells} axons, {_n_rois - _n_cells} non-axons)")
print(f"  Features: {_clf['keys']}")
del _clf

# ============================================================================
# AXON-SPECIFIC DETECTION PARAMETERS
# ============================================================================

print("\n" + "="*70)
print("APPLYING PARAMETERS — suite2p 1.0+ nested format")
print("="*70)

# suite2p 1.0+ uses NESTED settings. Flat keys like ops['combined'] or
# ops['threshold_scaling'] are silently ignored — they must go in the
# correct sub-dict. Legacy flat keys from the old ops template were dropped
# above so the nested values below are the single source of truth.

ops['save_path0'] = data_path
ops['diameter'] = 30              # flat top-level key — sigma=diameter/10=3px for sourcery (old ops used 20 on 2048-wide frame)

ops['io'] = {
    'combined': False,            # plane2 has no reg_outputs.npy — would crash if True
    'save_mat': True,             # backend save.py:save_mat() handles None values correctly
    'save_NWB': False,
    'save_ops_orig': True,
    'delete_bin': False,
    'move_bin': False,
}

ops['run'] = {
    'do_registration': 1,
    'do_regmetrics': False,
    'do_detection': True,
    'do_deconvolution': False,
    'multiplane_parallel': False,
}

# sourcery: grows ROI masks from SVD spatial components smoothed at sigma=diameter/10.
# No hardcoded 20% threshold — masks can extend to the full axon width.
ops['detection'] = {
    'algorithm': 'sourcery',
    'denoise': False,
    'block_size': (128, 128),
    'nbins': 5000,
    'bin_size': None,
    'highpass_time': 100,
    'threshold_scaling': 1.0,      # old working ops used 1.0; threshold_scaling only gates seeds, not mask size
    'npix_norm_min': 0.0,
    'npix_norm_max': 100,
    'max_overlap': 1.0,
    'soma_crop': False,
    'chan2_threshold': 0.25,
    'cellpose_chan2': False,
    'sparsery_settings': {
        'highpass_neuropil': 500,
        'max_ROIs': 5000,
        'spatial_scale': 2,
        'active_percentile': 0.0,
    },
    'sourcery_settings': {
        'connected': True,
        'max_iterations': 40,
        'smooth_masks': False,
    },
    'cellpose_settings': {
        'cellpose_model': 'cpsam',
        'img': 'max_proj / meanImg',
        'highpass_spatial': 0,
        'flow_threshold': 0.4,
        'cellprob_threshold': -2.0,
        'params': None,
        'params_chan2': None,
    }
}

ops['classification'] = {
    'classifier_path': None,
    'use_builtin_classifier': False,   # builtin classifier is trained on somas — it REJECTS axons
    'preclassify': 0.0,
}

ops['extraction'] = {
    'snr_threshold': 0.0,
    'batch_size': 500,
    'neuropil_extract': True,
    'neuropil_coefficient': 0.7,
    'inner_neuropil_radius': 2,
    'min_neuropil_pixels': 350,
    'lam_percentile': 50.0,        # V2: restored to match old working ops (was 0 in V1 to avoid OOM; try 25 if OOM recurs)
    'allow_overlap': True,
    'circular_neuropil': False,
}

ops['registration'] = {
    'align_by_chan2': False,
    'nimg_init': 300,
    'maxregshift': 0.1,
    'do_bidiphase': False,
    'bidiphase': 0.0,
    'batch_size': 100,
    'nonrigid': True,
    'maxregshiftNR': 5,
    'block_size': (128, 128),
    'smooth_sigma_time': 0,
    'smooth_sigma': 1.15,
    'spatial_taper': 50.0,         # V2: restored to default 50 px (was 3.45 in V1 — too small, caused registration artefacts)
    'th_badframes': 1.0,
    'norm_frames': True,
    'snr_thresh': 0.8,
    'subpixel': 10,
    'two_step_registration': False,
    'reg_tif': False,
    'reg_tif_chan2': False,
}

print(f"\n  algorithm: sourcery  |  diameter: {ops['diameter']} px (sigma={ops['diameter']/10:.0f}px)")
print(f"  threshold_scaling: {ops['detection']['threshold_scaling']}")
print(f"  combined: {ops['io']['combined']}  |  do_regmetrics: {ops['run']['do_regmetrics']}")


# ============================================================================
# DATABASE CONFIGURATION
# ============================================================================

db = {
    'data_path': [data_path],
    'save_path0': data_path,
    'nplanes': 4,                     # All 3 planes for correct alignment
    'nchannels': 1,
    'ignore_flyback': [3],         # Explicit plane index — plane2 is flyback (-1 does NOT work)
}

# Check for existing binary files — only check real planes, not flyback
suite2p_path = os.path.join(db['save_path0'], 'suite2p')
flyback_planes = set(db['ignore_flyback'])
real_planes = [p for p in range(db['nplanes']) if p not in flyback_planes]
binary_exists = False

print("\nChecking for existing binary files...")
if os.path.exists(suite2p_path):
    found = []
    for plane in real_planes:
        bin_file = os.path.join(suite2p_path, f'plane{plane}', 'data.bin')
        if os.path.exists(bin_file):
            found.append(plane)
            print(f"  Found: plane{plane}/data.bin")
        else:
            print(f"  Missing: plane{plane}/data.bin")

    # Always wipe stale ops/detection files — this runs whether or not binaries
    # exist, so cached parameters can never silently override the script's settings.
    stale = ['ops.npy', 'stat.npy', 'iscell.npy', 'F.npy', 'Fneu.npy', 'spks.npy', 'redcell.npy', 'Fall.mat']
    print("  Clearing stale detection outputs (keeping data.bin)...")
    for plane in range(db['nplanes']):
        plane_path = os.path.join(suite2p_path, f'plane{plane}')
        for fname in stale:
            fpath = os.path.join(plane_path, fname)
            if os.path.exists(fpath):
                os.remove(fpath)
                print(f"    Removed plane{plane}/{fname}")
    combined_path = os.path.join(suite2p_path, 'combined')
    if os.path.exists(combined_path):
        for fname in stale:
            fpath = os.path.join(combined_path, fname)
            if os.path.exists(fpath):
                os.remove(fpath)
                print(f"    Removed combined/{fname}")

    binary_exists = (len(found) == len(real_planes))
    if binary_exists:
        print("  All real-plane binaries found. Skipping TIFF conversion.")
        db['data_path'] = []
        db['look_one_level_down'] = False
        ops['fast_disk'] = suite2p_path
    else:
        print("  Some binaries missing — will convert from TIFFs.")

# ============================================================================
# RUN SUITE2P
# ============================================================================

print("\n" + "="*70)
print("FINAL OPS CHECK — values actually passed to run_s2p:")
print("="*70)
print(f"  algorithm:            {ops['detection']['algorithm']}")
print(f"  diameter:             {ops.get('diameter', 'NOT SET')}  (sigma={ops.get('diameter', 0)/10:.0f}px)")
print(f"  threshold_scaling:    {ops['detection']['threshold_scaling']}")
print(f"  npix_norm_max:        {ops['detection']['npix_norm_max']}")
print(f"  lam_percentile:       {ops['extraction']['lam_percentile']}")
print(f"  spatial_taper:        {ops['registration']['spatial_taper']}")
print(f"  combined:             {ops['io']['combined']}")
print(f"  do_regmetrics:        {ops['run']['do_regmetrics']}")
print(f"  use_builtin_clf:      {ops['classification']['use_builtin_classifier']}")
print(f"  nplanes (db):         {db.get('nplanes', 'NOT SET')}")
print(f"  ignore_flyback (db):  {db.get('ignore_flyback', 'NOT SET')}")
print(f"  total keys in ops:    {len(ops)}")
print("="*70 + "\n")

try:
    run_s2p(settings=ops, db=db)

    # Read back what suite2p actually saved — this is ground truth for what ran
    saved_ops_path = os.path.join(data_path, 'suite2p', 'plane0', 'ops.npy')
    if os.path.exists(saved_ops_path):
        saved = np.load(saved_ops_path, allow_pickle=True).item()
        saved_ext = saved.get('extraction', {})
        saved_reg = saved.get('registration', {})
        print("\n" + "="*70)
        print("SAVED ops.npy — what suite2p ACTUALLY used:")
        print(f"  diameter:          {saved.get('diameter', 'NOT SET')}")
        print(f"  spatial_hp_detect: {saved.get('spatial_hp_detect', 'NOT SET')}")
        print(f"  sparse_mode:       {saved.get('sparse_mode', 'NOT SET')}")
        print(f"  threshold_scaling: {saved.get('threshold_scaling', 'NOT SET')}")
        print(f"  spatial_scale:     {saved.get('spatial_scale', 'NOT SET')}")
        print(f"  lam_percentile (nested): {saved_ext.get('lam_percentile', 'NOT SET')}")
        print(f"  lam_percentile (flat):   {saved.get('lam_percentile', 'NOT SET')}")
        print(f"  spatial_taper (nested):  {saved_reg.get('spatial_taper', 'NOT SET')}")
        print(f"  spatial_taper (flat):    {saved.get('spatial_taper', 'NOT SET')}")
        print("="*70)

    print("\n" + "="*70)
    print("AXON DETECTION COMPLETE!")
    print("="*70)
    print("\nNext steps:")
    print("1. Inspect the detected axons in Suite2P GUI")
    print("2. Adjust threshold_scaling and cellprob_threshold if needed:")
    print("   - Lower values → more detections (including noise)")
    print("   - Higher values → fewer detections (missing faint axons)")
    print("3. Check iscell classifications")
    print()

except Exception as e:
    print(f"\n Error during processing: {e}")
    import traceback
    traceback.print_exc()
