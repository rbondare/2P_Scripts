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

# Load ops directly from your settings file
settings_file = resolve_runtime_path(r"Z:\group\joeschgrp\Group Members\Rima\2P_Scripts\suite2p\settings.npy")
data_path = resolve_runtime_path(r"Z:\group\joeschgrp\Group Members\Rima\DATA_2P\AnimalRB19_260320_1609")

if os.path.exists(settings_file):
    ops = np.load(settings_file, allow_pickle=True).item()
    # Remove version if present
    if 'version' in ops:
        ops.pop('version')
    print("✓ Loaded settings from file")
else:
    print(f"✗ Settings file not found: {settings_file}")
    print(f"  Platform: {platform.system()} / {os.name}")
    print("  Please check the path!")
    exit(1)

# Database configuration
db = {
    'data_path': [data_path],
    'save_path0': data_path,
    'nplanes': 4,
    'nchannels': 1,
}

# Check if suite2p output already exists
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
        print("\n All binary files exist. Skipping conversion.")
        db['data_path'] = []
        db['look_one_level_down'] = False
        ops['fast_disk'] = suite2p_path

# Run Suite2P
print("\n" + "="*60)
print("Running Suite2P with:")
print(f"  Data path: {db['data_path']}")
print(f"  Save path: {db['save_path0']}")
print(f"  Planes: {db['nplanes']}")
print(f"  Using existing binary: {binary_exists}")
print("="*60 + "\n")

try:
    run_s2p(settings=ops, db=db)
    print("\n" + "="*60)
    print("Suite2P processing complete!")
    print("="*60)
except Exception as e:
    print(f"\nError during processing: {e}")
    import traceback
    traceback.print_exc()