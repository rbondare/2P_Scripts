#!/usr/bin/env python3
"""
Suite2P Drifting-Grating (DG) / Bouton Detection Script — V2
Sparsery-based bouton detection, tuned for AnimalDG3-type recordings.

Changes from V1
---------------
All detection/registration/extraction/classification parameters from V1 are
kept exactly as they were — these are the values known to work well for this
data (sparsery, threshold_scaling=0.2, lam_percentile=10, etc.). What's new in
V2 is the run-time checkpoint structure, ported from
run_suite2p_RB_axons_V2.py:

1. GPU check — aborts before running if CUDA is unavailable. V1 forced
   torch_device='cpu'; that explicit override is removed here so the default
   (CUDA) is used instead.
2. Classifier validation — confirms classifier_boutons_tony.npy exists and has
   the expected stats/iscell/keys structure before running.
3. Binary-file check — skips TIFF conversion if data.bin already exists for
   every real (non-flyback) plane.
4. Unconditional stale-output wipe — ops/stat/F/Fneu/spks/iscell/Fall.mat are
   removed before each run so cached parameters from a previous run can never
   silently persist.
5. Final-settings printout before run_s2p, and a saved-ops readback after, so
   what actually ran is always visible.
6. io['combined'] is set False for the run itself — V1's combined=True would
   make suite2p's native combined() crash on the flyback plane (it scans every
   planeX folder, including the unprocessed flyback one, for files that don't
   exist there). build_combined_without_flyback() (ported from
   run_suite2p_RB_axons_V1.py) is called after run_s2p to rebuild a combined
   Fall.mat from only the real planes, so the end result still matches V1's
   combined=True intent without the crash.

Paths/db are unchanged from V1: AnimalDG3_260320_1311, nplanes=4,
ignore_flyback=[3].
"""

import numpy as np
import os
import shutil
import tempfile
import torch
from pathlib import PureWindowsPath
from suite2p import run_s2p, default_settings
from suite2p import io as suite2p_io


def build_combined_without_flyback(source_save_folder, analysis_planes):
    """Rebuild combined folder from only the analysis planes (no flyback).

    Uses a temp folder so suite2p's combined() never sees the flyback plane
    directory. Always rebuilds — caller is responsible for calling this only
    after a fresh detection run.
    """
    combined_dir = os.path.join(source_save_folder, "combined")

    if os.path.exists(combined_dir):
        print("Removing existing combined folder to rebuild from fresh detection outputs...")
        shutil.rmtree(combined_dir)

    temp_root = tempfile.mkdtemp(prefix="suite2p_combined_")
    try:
        for plane in analysis_planes:
            src_plane = os.path.join(source_save_folder, f"plane{plane}")
            if not os.path.isdir(src_plane):
                continue
            dst_plane = os.path.join(temp_root, f"plane{plane}")
            os.makedirs(dst_plane, exist_ok=True)
            for fname in [
                "db.npy", "settings.npy",
                "reg_outputs.npy", "detect_outputs.npy",
                "stat.npy", "F.npy", "Fneu.npy",
                "spks.npy", "iscell.npy", "redcell.npy",
            ]:
                src = os.path.join(src_plane, fname)
                if os.path.exists(src):
                    shutil.copy2(src, os.path.join(dst_plane, fname))

        suite2p_io.combined(temp_root, save=True)

        temp_combined = os.path.join(temp_root, "combined")
        final_combined = os.path.join(source_save_folder, "combined")
        if os.path.exists(final_combined):
            shutil.rmtree(final_combined)
        shutil.copytree(temp_combined, final_combined)
        print(f"Created combined folder (planes {analysis_planes}): {final_combined}")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


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


# Old/alternate recording from V1 (kept for reference, not used):
# data_path = resolve_runtime_path(r"Y:\Rotation Students\Dow\DATA_2P\AnimalDG1_260302_1422")
data_path = resolve_runtime_path(r"Z:\joeschgrp\Rotation Students\Dow\DATA_2P\AnimalDG1_260224_1508")
classifier_path = resolve_runtime_path(r"Z:\joeschgrp\Group Members\Rima\2P_Scripts\suite2p\classifier_boutons_tony.npy")

# Validate bouton classifier before doing anything else
print("\nValidating bouton classifier...")
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
print(f"  OK — {_n_rois} labeled ROIs ({_n_cells} boutons, {_n_rois - _n_cells} non-boutons)")
print(f"  Features: {_clf['keys']}")
del _clf

# ============================================================================
# DG/BOUTON DETECTION PARAMETERS — unchanged from V1
# ============================================================================

print("\n" + "="*70)
print("APPLYING PARAMETERS (kept from DG_V1)")
print("="*70)

settings = default_settings()

# Top level
settings['tau']      = 0.7
settings['fs']       = 10.0
settings['diameter'] = [3.0, 3.0]
# torch_device left at default_settings()'s default (CUDA) — V1's explicit
# 'cpu' override is dropped now that the GPU check above requires CUDA.

# Run
settings['run']['do_registration']     = 1
settings['run']['do_regmetrics']       = True
settings['run']['do_detection']        = True
settings['run']['do_deconvolution']    = False
settings['run']['multiplane_parallel'] = False

# IO — combined left False for the run itself; see module docstring point 6.
# Rebuilt safely from real planes only via build_combined_without_flyback() below.
settings['io']['combined']      = False
settings['io']['save_mat']      = True
settings['io']['save_NWB']      = False
settings['io']['save_ops_orig'] = True
settings['io']['delete_bin']    = False
settings['io']['move_bin']      = False

# Registration
settings['registration']['nimg_init']           = 500
settings['registration']['batch_size']          = 250
settings['registration']['maxregshift']         = 0.12
settings['registration']['nonrigid']            = True
settings['registration']['block_size']          = [128, 128]
settings['registration']['maxregshiftNR']       = 5
settings['registration']['smooth_sigma']        = 1.5
settings['registration']['smooth_sigma_time']   = 0.5
settings['registration']['th_badframes']        = 0.8
settings['registration']['snr_thresh']          = 3.5
settings['registration']['norm_frames']         = True
settings['registration']['subpixel']            = 10
settings['registration']['align_by_chan2']      = False
settings['registration']['reg_tif']             = False
settings['registration']['reg_tif_chan2']       = False

# Detection
settings['detection']['algorithm']         = 'sparsery'
settings['detection']['denoise']           = True
settings['detection']['threshold_scaling'] = 0.2
settings['detection']['max_overlap']       = 0.75
settings['detection']['highpass_time']     = 50
settings['detection']['soma_crop']         = False
settings['detection']['sparsery_settings']['highpass_neuropil'] = 25
settings['detection']['sparsery_settings']['spatial_scale']     = 1
settings['detection']['sparsery_settings']['max_ROIs']          = 12000
settings['detection']['sparsery_settings']['active_percentile'] = 0.0

# Classification
settings['classification']['classifier_path']        = classifier_path
settings['classification']['use_builtin_classifier']  = False
settings['classification']['preclassify']             = 0.0

# Extraction
settings['extraction']['neuropil_extract']      = True
settings['extraction']['neuropil_coefficient']  = 0.7
settings['extraction']['inner_neuropil_radius'] = 1
settings['extraction']['min_neuropil_pixels']   = 50
settings['extraction']['lam_percentile']        = 10.0
settings['extraction']['allow_overlap']         = True

# Spike deconvolution preprocessing
settings['dcnv_preprocess']['baseline']         = 'maximin'
settings['dcnv_preprocess']['win_baseline']     = 60.0
settings['dcnv_preprocess']['sig_baseline']     = 10.0
settings['dcnv_preprocess']['prctile_baseline'] = 8.0

print(f"\n  algorithm: {settings['detection']['algorithm']}  |  diameter: {settings['diameter']}")
print(f"  threshold_scaling: {settings['detection']['threshold_scaling']}")
print(f"  combined (run-time): {settings['io']['combined']}  |  do_regmetrics: {settings['run']['do_regmetrics']}")


# ============================================================================
# DATABASE CONFIGURATION
# ============================================================================

db = {
    'data_path': [data_path],
    'save_path0': data_path,
    'nplanes':        4,
    'nchannels':      1,
    'ignore_flyback': [3],
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
        settings['fast_disk'] = suite2p_path
    else:
        print("  Some binaries missing — will convert from TIFFs.")

# ============================================================================
# RUN SUITE2P
# ============================================================================

print("\n" + "="*70)
print("FINAL SETTINGS CHECK — values actually passed to run_s2p:")
print("="*70)
print(f"  algorithm:              {settings['detection']['algorithm']}")
print(f"  diameter:               {settings.get('diameter', 'NOT SET')}")
print(f"  threshold_scaling:      {settings['detection']['threshold_scaling']}")
print(f"  highpass_neuropil:      {settings['detection']['sparsery_settings']['highpass_neuropil']}")
print(f"  spatial_scale:          {settings['detection']['sparsery_settings']['spatial_scale']}")
print(f"  max_ROIs:               {settings['detection']['sparsery_settings']['max_ROIs']}")
print(f"  lam_percentile:         {settings['extraction']['lam_percentile']}")
print(f"  neuropil_coefficient:   {settings['extraction']['neuropil_coefficient']}")
print(f"  classifier_path:        {settings['classification']['classifier_path']}")
print(f"  do_regmetrics:          {settings['run']['do_regmetrics']}")
print(f"  combined (run-time):    {settings['io']['combined']}  (rebuilt safely after run — see below)")
print(f"  nplanes (db):           {db.get('nplanes', 'NOT SET')}")
print(f"  ignore_flyback (db):    {db.get('ignore_flyback', 'NOT SET')}")
print(f"  total keys in settings: {len(settings)}")
print("="*70 + "\n")

try:
    run_s2p(settings=settings, db=db)

    # Rebuild combined Fall.mat from real planes only — run-time io['combined']
    # is False to avoid suite2p's native combined() crashing on the flyback plane.
    suite2p_dir = os.path.join(data_path, 'suite2p')
    print("\nBuilding combined folder from real planes only...")
    build_combined_without_flyback(suite2p_dir, real_planes)

    # Read back what suite2p actually saved — this is ground truth for what ran
    saved_settings_path = os.path.join(data_path, 'suite2p', 'plane0', 'ops.npy')
    if os.path.exists(saved_settings_path):
        saved = np.load(saved_settings_path, allow_pickle=True).item()
        saved_ext = saved.get('extraction', {})
        print("\n" + "="*70)
        print("SAVED ops.npy — what suite2p ACTUALLY used:")
        print(f"  diameter:          {saved.get('diameter', 'NOT SET')}")
        print(f"  spatial_hp_detect: {saved.get('spatial_hp_detect', 'NOT SET')}")
        print(f"  sparse_mode:       {saved.get('sparse_mode', 'NOT SET')}")
        print(f"  threshold_scaling: {saved.get('threshold_scaling', 'NOT SET')}")
        print(f"  spatial_scale:     {saved.get('spatial_scale', 'NOT SET')}")
        print(f"  lam_percentile (nested): {saved_ext.get('lam_percentile', 'NOT SET')}")
        print("="*70)

    print("\n" + "="*70)
    print("DG/BOUTON DETECTION COMPLETE!")
    print("="*70)
    print("\nNext steps:")
    print("1. Inspect the detected ROIs in Suite2P GUI")
    print("2. Check iscell classifications assigned via classifier_boutons_tony.npy")
    print("3. Adjust threshold_scaling / spatial_scale if needed:")
    print("   - Lower threshold_scaling → more detections (including noise)")
    print("   - Higher threshold_scaling → fewer detections (missing faint boutons)")
    print()

except Exception as e:
    print(f"\n Error during processing: {e}")
    import traceback
    traceback.print_exc()
