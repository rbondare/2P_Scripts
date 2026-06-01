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
from pathlib import PureWindowsPath
from suite2p import run_s2p


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

if os.path.exists(settings_file):
    ops = np.load(settings_file, allow_pickle=True).item()
    if 'version' in ops:
        ops.pop('version')
    print("Loaded base settings from file")
else:
    print(f"Settings file not found: {settings_file}")
    exit(1)

# ============================================================================
# AXON-SPECIFIC DETECTION PARAMETERS
# ============================================================================

print("\n" + "="*70)
print("APPLYING AXON DETECTION PARAMETERS")
print("="*70)

# DETECTION SETTINGS (most critical for axons)
# ============================================================================
ops['diameter'] = 12                      # Axon diameter ~1-2 μm (soma: 15-20)
ops['threshold_scaling'] = 1.2            # Detection threshold (lower = more sensitive)
ops['cellprob_threshold'] = 0.0           # Cellpose probability threshold
ops['flow_threshold'] = 0.4               # Cellpose flow matching threshold

# Spatial filtering for fine structures
ops['spatial_hp_detect'] = 25             # High-pass filter for detection
ops['spatial_hp_reg'] = 42                # High-pass filter for registration

# Block and ROI parameters
ops['block_size'] = [64.0, 64.0]          # SMALLER blocks for finer detail (soma: 128)
ops['max_overlap'] = 0.7
ops['soma_crop'] = False                  # Don't crop to soma region
ops['pre_smooth'] = 0                     # Pre-smoothing before registration

# ROI size constraints
ops['npix_norm_min'] = 10.0               # Minimum ROI size
ops['npix_norm_max'] = 100.0              # Maximum ROI size (axons are small!)

# Sparsery-specific detection
ops['spatial_scale'] = 0                  # Spatial scale for detection
ops['connected'] = True                   # Use connected components
ops['max_iterations'] = 20                # Max iterations for detection

# EXTRACTION SETTINGS
# ============================================================================
ops['neuropil_extract'] = True
ops['allow_overlap'] = False              # Match old settings
ops['inner_neuropil_radius'] = 2          # Standard neuropil radius
ops['min_neuropil_pixels'] = 300          # Minimum neuropil pixels
ops['neucoeff'] = 0.5                     # Neuropil coefficient (old: neuropil_coefficient)
ops['circular_neuropil'] = False          # Allow irregular neuropil shape
ops['lam_percentile'] = 50.0

# Extraction batch size (different from registration batch)
ops['batch_size'] = 500                   # Extraction batch size

# Denoising
ops['denoise'] = False
ops['spikedetect'] = False                # No spike detection for axons
ops['roidetect'] = True                   # Enable ROI detection

# REGISTRATION SETTINGS
# ============================================================================
ops['do_registration'] = 1
ops['nonrigid'] = True                    # Use nonrigid registration
ops['norm_frames'] = True                 # ← CRITICAL: WAS MISSING (causes KeyError)
ops['nimg_init'] = 400
ops['maxregshift'] = 0.12
ops['maxregshiftNR'] = 5
ops['smooth_sigma'] = 1.5
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
ops['use_builtin_classifier'] = True

# OTHER SETTINGS
# ============================================================================
ops['tau'] = 0.7                          # Time constant
ops['fs'] = 10.0                          # Sampling rate
ops['keep_movie_raw'] = False
ops['sparse_mode'] = False
ops['pad_fft'] = True
ops['force_refImg'] = False

print("\nApplied axon-specific parameters:")
print(f"  - Diameter: {ops['diameter']}")
print(f"  - Threshold scaling: {ops['threshold_scaling']}")
print(f"  - Cellprob threshold: {ops['cellprob_threshold']}")
print(f"  - Spatial HP detect: {ops['spatial_hp_detect']}")
print(f"  - Block size: {ops['block_size']}")
print(f"  - Inner neuropil radius: {ops['inner_neuropil_radius']}")
print(f"  - norm_frames: {ops['norm_frames']} ✓ (required)")


# ============================================================================
# DATABASE CONFIGURATION
# ============================================================================

db = {
    'data_path': [data_path],
    'save_path0': data_path,
    'nplanes': 3,                      # All 3 planes for correct alignment
    'nchannels': 1,
    'ignore_flyback': [-1],             # Skip processing plane 2 (flyback)
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
