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
settings_file = resolve_runtime_path(r"Z:\group\joeschgrp\Group Members\Rima\2P_Scripts\suite2p\settings.npy")
data_path = resolve_runtime_path(r"Z:\group\joeschgrp\Group Members\Rima\DATA_2P_OTHER\AnimalRB6\AnimalRB6_250302_1628")

if os.path.exists(settings_file):
    ops = np.load(settings_file, allow_pickle=True).item()
    if 'version' in ops:
        ops.pop('version')
    print("✓ Loaded base settings from file")
else:
    print(f"✗ Settings file not found: {settings_file}")
    exit(1)

# ============================================================================
# AXON-SPECIFIC DETECTION PARAMETERS
# ============================================================================

print("\n" + "="*70)
print("APPLYING AXON DETECTION PARAMETERS")
print("="*70)

# DETECTION SETTINGS (most critical for axons)
# ============================================================================
ops['detection'] = {
    'algorithm': 'sparsery',
    'diameter': [3.0, 3.0],            # ← MUCH SMALLER: axons ~1-3 μm
    'threshold_scaling': 0.5,          # ← LOWER: detect fainter axons (default ~1.0)
    'cellprob_threshold': 0.1,         # ← MUCH LOWER: allow fainter detections (default ~0.5)
    'flow_threshold': 0.4,             # ← LOWER: less stringent flow matching
    
    # Spatial filtering for fine axon structures
    'spatial_hp_detect': 10.0,         # ← LOWER: preserve small structures (default ~25)
    'spatial_scale': 2,                # ← Higher spatial scale for finer detection (axons vs soma)
    
    # Cellpose settings (for axon shapes)
    'cellpose_chan2': False,
    
    # Block parameters
    'block_size': [64.0, 64.0],        # ← SMALLER blocks for finer detail (soma: 128)
    'max_overlap': 0.7,
    
    # ROI size constraints
    'max_ROIs': 10000,                 # Allow more ROIs for sparse axons
    'npix_norm_min': 0.0,              # ← LOWER: allow tiny axon segments
    'npix_norm_max': 100.0,            # ← LOWER: cap for axons (soma: 200+)
    
    # Neuropil and background
    'highpass_time': 100,              # Temporal high-pass
    'nbins': 5000,
    'denoise': False,
}

# EXTRACTION SETTINGS
# ============================================================================
ops['extraction'] = {
    'neuropil_extract': True,
    'allow_overlap': False,            # Keep False to avoid duplicate axons
    'inner_neuropil_radius': 1,        # ← SMALLER: axons don't have large neuropil
    'min_neuropil_pixels': 50,         # ← SMALLER: axons have less surrounding area (default 350)
    'circular_neuropil': False,        # Allow irregular neuropil shape
    'neuropil_coefficient': 0.5,       # Neuropil subtraction strength
    'lam_percentile': 50.0,
    'batch_size': 100,
    'snr_threshold': 0.0,
}

# REGISTRATION SETTINGS (keep most settings, but adjust for speed)
# ============================================================================
ops['registration'] = {
    'do_registration': 1,
    'nonrigid': True,                  # Use nonrigid for better local alignment
    'batch_size': 100,
    'block_size': [128.0, 128.0],
    'nimg_init': 400,
    'maxregshift': 0.12,
    'maxregshiftNR': 5,
    'smooth_sigma': 1.5,
    'bidiphase': 0.0,
    'do_bidiphase': False,
    'align_by_chan2': False,
    'reg_tif': False,
    'reg_tif_chan2': False,
    'snr_thresh': 1.2,
    'spatial_taper': 3.45,
    'subpixel': 10,
    'th_badframes': 1.0,
    'two_step_registration': True,
}

# DECONVOLUTION SETTINGS
# ============================================================================
ops['dcnv_preprocess'] = {
    'baseline': 'maximin',
    'win_baseline': 60.0,
    'sig_baseline': 10.0,
    'prctile_baseline': 8.0,
}

# RUN SETTINGS
# ============================================================================
ops['run'] = {
    'do_registration': 1,
    'do_detection': True,              # Perform detection
    'do_deconvolution': True,
    'do_regmetrics': False,            # Skip metrics to avoid memory issues
    'multiplane_parallel': False,
}

# I/O SETTINGS
# ============================================================================
ops['io'] = {
    'combined': True,
    'delete_bin': False,
    'move_bin': False,
    'save_NWB': False,
    'save_mat': True,
    'save_ops_orig': True,
}

# OTHER SETTINGS
# ============================================================================
ops['tau'] = 0.7
ops['fs'] = 10.0
ops['torch_device'] = 'cpu'
ops['version'] = '1.0.0.1'

print("\n✓ Applied axon-specific parameters:")
print(f"  - Diameter: {ops['detection']['diameter']}")
print(f"  - Threshold scaling: {ops['detection']['threshold_scaling']}")
print(f"  - Cellprob threshold: {ops['detection']['cellprob_threshold']}")
print(f"  - Spatial HP detect: {ops['detection']['spatial_hp_detect']}")
print(f"  - Block size: {ops['detection']['block_size']}")
print(f"  - Inner neuropil radius: {ops['extraction']['inner_neuropil_radius']}")

# ============================================================================
# DATABASE CONFIGURATION
# ============================================================================

db = {
    'data_path': [data_path],
    'save_path0': data_path,
    'nplanes': 4,
    'nchannels': 1,
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
            print(f"  ✓ Found: plane{plane}/data.bin")
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
print(f"    - Diameter: {ops['detection']['diameter']} μm")
print(f"    - Cellprob threshold: {ops['detection']['cellprob_threshold']} (lower = more sensitive)")
print(f"    - Threshold scaling: {ops['detection']['threshold_scaling']} (lower = detects fainter)")
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
