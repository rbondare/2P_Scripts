#!/usr/bin/env python3
"""
Suite2P Axon Detection Script
Optimized for detecting faint, sparse axons as shown in fluorescence imaging

Key differences from soma detection:
- Much smaller diameter (axons ~1-2 μm vs soma ~15-20 μm)
- Much lower detection thresholds
- Higher spatial_scale to detect fine structures
- Adjusted neuropil parameters
"""

import numpy as np
import os
import platform
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


# Load base settings from soma detection
settings_file = resolve_runtime_path(r"Z:\joeschgrp\Group Members\Rima\2P_Scripts\suite2p\settings.npy")
data_path = resolve_runtime_path(r"Z:\joeschgrp\Group Members\Rima\DATA_2P_OTHER\AnimalRB6\AnimalRB6_250302_1628")
classifier_path = resolve_runtime_path(r"Z:\joeschgrp\Group Members\Rima\DATA_2P_OTHER\AnimalRB6\AnimalRB6_250302_1656\suite2p\plane0\NAaxons.npy")

if os.path.exists(settings_file):
    ops = np.load(settings_file, allow_pickle=True).item()
    if 'version' in ops:
        ops.pop('version')
    print("Loaded base settings from file")
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
print("APPLYING AXON DETECTION PARAMETERS")
print("="*70)

# DETECTION SETTINGS (most critical for axons)
# ============================================================================
# sparse_mode = True is the MOST CRITICAL setting: activates the sparsery
# algorithm which grows masks along correlated-activity pixels rather than
# assuming round shapes. Without it, suite2p segments individual bouton
# hot-spots instead of the elongated process shown in the image.
ops['sparse_mode'] = True                 # CRITICAL: follow correlated activity along processes
ops['diameter'] = 20                      # Match axon TUBE WIDTH visible in image (~3-8 px), not soma size
ops['threshold_scaling'] = 0.3            # Lower = more sensitive; needed to catch dim segments along processes
ops['spatial_scale'] = 1                  # Detect at finest scale (1 = ~3 px); 0 = auto-estimate

# Cellpose parameters (not used in sparse_mode but kept for fallback)
ops['cellprob_threshold'] = 0.2
ops['flow_threshold'] = 1.2

# Spatial filtering for fine structures
ops['spatial_hp_detect'] = 12             # Lower than default 25: preserves thin tubular processes
ops['spatial_hp_reg'] = 42                # High-pass filter for registration

# Block and ROI parameters
ops['block_size'] = [64.0, 64.0]          # SMALLER blocks for finer detail (soma: 128)
ops['max_overlap'] = 0.75                 # Axons cross each other; allow high overlap
ops['soma_crop'] = False                  # Don't crop to soma region
ops['pre_smooth'] = 0                     # Pre-smoothing before registration

# ROI size constraints
# Elongated axon segments span many pixels — max must be large or long
# processes get clipped into fragments. Min is small because thin cross-sections
# have few pixels even at moderate length.
ops['npix_norm_min'] = 15.0               # Minimum ROI size
ops['npix_norm_max'] = 3000.0             # Large: winding axons accumulate many pixels

# Sparsery-specific detection
ops['connected'] = True                   # Keep masks as connected components along process
ops['max_iterations'] = 50                # More iterations to trace full axon length

# EXTRACTION SETTINGS
# ============================================================================
ops['neuropil_extract'] = True
ops['allow_overlap'] = True              # Match old settings
ops['inner_neuropil_radius'] = 1          # Minimal: axons have no meaningful neuropil halo
ops['min_neuropil_pixels'] = 50           # Much smaller surround needed for thin processes
ops['neucoeff'] = 0.5                     # Neuropil coefficient (old: neuropil_coefficient)
ops['circular_neuropil'] = False          # Allow irregular neuropil shape
ops['lam_percentile'] = 20.0              # Lower = mask grows further along the process (50 = compact blobs)

# Extraction batch size (different from registration batch)
ops['batch_size'] = 500                   # Extraction batch size

# Denoising
ops['denoise'] = True
ops['spikedetect'] = False                # No spike detection for axons
ops['roidetect'] = True                   # Enable ROI detection

# REGISTRATION SETTINGS
# ============================================================================
ops['do_registration'] = 1
ops['nonrigid'] = True                    # Use nonrigid registration
ops['norm_frames'] = True                 # 
ops['nimg_init'] = 400
ops['maxregshift'] = 0.12
ops['maxregshiftNR'] = 5
ops['smooth_sigma'] = 1.15                # Less blur → preserves thin tubular structures
ops['smooth_sigma_time'] = 0              # No temporal smoothing
ops['bidiphase'] = 0.0
ops['do_bidiphase'] = False
ops['align_by_chan'] = 0                  # Don't align by channel 2
ops['reg_tif'] = False
ops['reg_tif_chan2'] = False
ops['snr_thresh'] = 1.2
ops['spatial_taper'] = 3.45
ops['subpixel'] = 10
ops['th_badframes'] = 1.0
ops['two_step_registration'] = True
ops['1Preg'] = False                      # No 1P registration

# RUN SETTINGS
# ============================================================================
ops['do_detection'] = True                # Perform detection
ops['do_deconvolution'] = False
ops['do_regmetrics'] = False              # Skip metrics to avoid memory issues
ops['multiplane_parallel'] = False

# I/O SETTINGS
# ============================================================================
ops['combined'] = True
ops['delete_bin'] = False
ops['move_bin'] = False
ops['save_NWB'] = False
ops['save_mat'] = False                   # Disable MATLAB saving to avoid None conversion errors
ops['save_path0'] = data_path             # Set save path (data_path defined earlier)

# CLASSIFICATION SETTINGS
# ============================================================================
ops['preclassify'] = 0.0
ops['use_builtin_classifier'] = False     # Use custom axon classifier only
ops['classifier_path'] = classifier_path  # NAaxons.npy

# OTHER SETTINGS
# ============================================================================
ops['tau'] = 0.7                          # Time constant
ops['fs'] = 10.0                          # Sampling rate
ops['keep_movie_raw'] = False
ops['sparse_mode'] = True                 # Set above in detection section; repeated here for clarity
ops['pad_fft'] = True
ops['force_refImg'] = False

print("\nApplied axon-specific parameters:")
print(f"  - sparse_mode: {ops['sparse_mode']}  ← CRITICAL: sparsery algorithm for elongated structures")
print(f"  - Diameter: {ops['diameter']} px  (match tube width, not soma size)")
print(f"  - Threshold scaling: {ops['threshold_scaling']}")
print(f"  - Spatial scale: {ops['spatial_scale']}  (1 = finest ~3px scale)")
print(f"  - Spatial HP detect: {ops['spatial_hp_detect']}")
print(f"  - npix_norm_max: {ops['npix_norm_max']}  (large for elongated ROIs)")
print(f"  - max_overlap: {ops['max_overlap']}  (axons cross)")
print(f"  - Inner neuropil radius: {ops['inner_neuropil_radius']}")
print(f"  - norm_frames: {ops['norm_frames']}")
print(f"  - classifier: NAaxons.npy (custom, builtin disabled)")


# ============================================================================
# DATABASE CONFIGURATION
# ============================================================================

db = {
    'data_path': [data_path],
    'save_path0': data_path,
    'nplanes': 3,                      # All 3 planes for correct alignment
    'nchannels': 1,
    'ignore_flyback': [2],             # Explicit plane index — plane2 is flyback (-1 does NOT work)
}

# Check for existing binary files
suite2p_path = os.path.join(db['save_path0'], 'suite2p')
binary_exists = False

print("\nChecking for existing binary files...")
if os.path.exists(suite2p_path):
    for plane in range(db['nplanes']):
        plane_path = os.path.join(suite2p_path, f'plane{plane}')
        bin_file = os.path.join(plane_path, 'data.bin')
        if os.path.exists(bin_file):
            print(f"Found: plane{plane}/data.bin")
            binary_exists = True
        else:
            binary_exists = False
            break
    
    if binary_exists:
        print("\n  All binary files exist. Skipping conversion.")
        db['data_path'] = []
        db['look_one_level_down'] = False
        ops['fast_disk'] = suite2p_path

# ============================================================================
# RUN SUITE2P
# ============================================================================

print("\n" + "="*70)
print("RUNNING SUITE2P - AXON DETECTION MODE")
print("="*70)
print(f"  Data path: {db['data_path']}")
print(f"  Save path: {db['save_path0']}")
print(f"  Planes: {db['nplanes']}")
print(f"  Using existing binary: {binary_exists}")
print("\n  Detection settings:")
print(f"    - Algorithm: sparsery")
print(f"    - Diameter: {ops['diameter']} μm")
print(f"    - Cellprob threshold: {ops['cellprob_threshold']} (lower = more sensitive)")
print(f"    - Threshold scaling: {ops['threshold_scaling']} (lower = detects fainter)")
print("\n" + "="*70 + "\n")

try:
    run_s2p(settings=ops, db=db)
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
