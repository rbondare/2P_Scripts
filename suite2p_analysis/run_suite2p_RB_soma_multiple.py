import numpy as np
import os
import platform
from pathlib import PureWindowsPath
from suite2p import run_s2p
from suite2p import io as suite2p_io
from suite2p.run_s2p import run_plane
import time
import tempfile
import traceback
import glob
import subprocess
import shutil



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


def build_combined_without_flyback(source_save_folder, analysis_planes):
    """Create or complete the combined folder from only the analysis planes.
    
    If combined folder exists but Fall.mat is missing (e.g., from previous permission error),
    only regenerate the .mat file without rebuilding the entire combined folder.
    """
    combined_dir = os.path.join(source_save_folder, "combined")
    combined_mat = os.path.join(combined_dir, "Fall.mat")
    
    # If combined folder exists but Fall.mat is missing, just generate the .mat file
    if os.path.exists(combined_dir) and not os.path.exists(combined_mat):
        print(f"Combined folder exists but Fall.mat missing — regenerating .mat file only...")
        try:
            suite2p_io.save_mat(combined_dir)
            print(f"Fall.mat regenerated successfully in {combined_dir}")
            return
        except Exception as e:
            print(f"Error regenerating Fall.mat: {e}")
            raise
    
    # If combined folder and Fall.mat both exist, nothing to do
    if os.path.exists(combined_dir) and os.path.exists(combined_mat):
        print(f"Combined folder with Fall.mat already exists — skipping")
        return
    
    # Full rebuild: combined folder doesn't exist yet
    temp_root = tempfile.mkdtemp(prefix="suite2p_combined_")
    try:
        for plane in analysis_planes:
            src_plane = os.path.join(source_save_folder, f"plane{plane}")
            if not os.path.isdir(src_plane):
                continue
            dst_plane = os.path.join(temp_root, f"plane{plane}")
            os.makedirs(dst_plane, exist_ok=True)

            # Copy only the lightweight files needed by suite2p's combined()
            for fname in [
                "db.npy",
                "settings.npy",
                "reg_outputs.npy",
                "detect_outputs.npy",
                "stat.npy",
                "F.npy",
                "Fneu.npy",
                "spks.npy",
                "iscell.npy",
                "redcell.npy",
            ]:
                src = os.path.join(src_plane, fname)
                if os.path.exists(src):
                    shutil.copy2(src, os.path.join(dst_plane, fname))

        # Run Suite2P's native combine on the filtered temp folder
        suite2p_io.combined(temp_root, save=True)

        temp_combined = os.path.join(temp_root, "combined")
        final_combined = os.path.join(source_save_folder, "combined")
        if os.path.exists(final_combined):
            shutil.rmtree(final_combined)
        shutil.copytree(temp_combined, final_combined)
        
        print(f"Created combined folder without flyback plane: {final_combined}")
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)

# Load ops from your settings file
settings_file = resolve_runtime_path(r"C:\Users\rbondarenko\projects\settings.npy")

if os.path.exists(settings_file):
    ops = np.load(settings_file, allow_pickle=True).item()
    if 'version' in ops:
        ops.pop('version')
    print("Loaded settings from file\n")
else:
    print(f"Settings file not found: {settings_file}")
    print(f"Platform: {platform.system()} / {os.name}")
    exit(1)

# List of folders to process
data_folders = [
    resolve_runtime_path(r"Z:\group\joeschgrp\Group Members\Rima\DATA_2P\AnimalRB19_260303_1546"),
    resolve_runtime_path(r"Z:\group\joeschgrp\Group Members\Rima\DATA_2P\AnimalRB19_260303_1705"),
    resolve_runtime_path(r"Z:\group\joeschgrp\Group Members\Rima\DATA_2P\AnimalRB19_260303_1745")
    # Add more folders here
]

# Common settings
nplanes = 4
nchannels = 1
ignore_flyback = [nplanes - 1]
analysis_planes = [plane for plane in range(nplanes) if plane not in ignore_flyback]
use_local_copy = True #t to True to process locally

# Process each folder
for i, folder in enumerate(data_folders, 1):
    print("="*70)
    print(f"PROCESSING FOLDER {i}/{len(data_folders)}: {os.path.basename(folder)}")
    print("="*70)
    
    if not os.path.exists(folder):
        print(f"Folder not found, skipping: {folder}\n")
        continue

    if use_local_copy:
        # Copy remote folder to local workspace to avoid network locking issues
        LOCAL_COPY_ROOT = r"C:\Users\rbondarenko\projects\DATA_2P"
        os.makedirs(LOCAL_COPY_ROOT, exist_ok=True)
        folder_basename = os.path.basename(folder.rstrip('/\\'))
        local_dst = os.path.join(LOCAL_COPY_ROOT, folder_basename)
        if os.path.exists(local_dst):
            print(f"Local copy already exists: {local_dst} — skipping copy and using local folder")
            folder = local_dst
        else:
            print(f"Copying {folder} -> {local_dst} (this may take a while)")
            try:
                if os.name == 'nt':
                    # use robocopy for robust copying on Windows
                    cmd = [
                        'robocopy',
                        folder,
                        local_dst,
                        '/MIR',
                        '/Z',
                        '/J',
                        '/R:3',
                        '/W:5',
                        '/MT:8'
                    ]
                    # Ensure destination parent exists for robocopy
                    os.makedirs(os.path.dirname(local_dst), exist_ok=True)
                    res = subprocess.run(cmd, shell=False)
                    # robocopy exit codes < 8 indicate success or minor issues
                    if res.returncode >= 8:
                        print(f"robocopy failed with exit code {res.returncode}; falling back to shutil.copytree")
                        shutil.copytree(folder, local_dst)
                else:
                    shutil.copytree(folder, local_dst)
                print("Copy complete")
                folder = local_dst
            except Exception as e:
                print(f"Warning: failed to copy to local path: {e}")
                traceback.print_exc()
                print("Proceeding with original folder path")
    else:
        print("Processing directly from the group drive (no local copy)")
    
    # Database configuration for this folder
    db = {
        'data_path': [folder],
        'save_path0': folder,
        'nplanes': nplanes,
        'nchannels': nchannels,
        'ignore_flyback': ignore_flyback,
    }
    
    # Check if suite2p output already exists
    suite2p_path = os.path.join(folder, 'suite2p')
    binary_exists = False

    if os.path.exists(suite2p_path):
        for plane in analysis_planes:
            plane_path = os.path.join(suite2p_path, f'plane{plane}')
            bin_file = os.path.join(plane_path, 'data.bin')
            if os.path.exists(bin_file):
                binary_exists = True
            else:
                binary_exists = False
                break
        
        if binary_exists:
            # ensure outputs (.mat) exist for all planes before skipping
            all_mats = True
            for plane in analysis_planes:
                plane_dir = os.path.join(suite2p_path, f'plane{plane}')
                mat_files = glob.glob(os.path.join(plane_dir, '*.mat')) if os.path.exists(plane_dir) else []
                if not mat_files:
                    all_mats = False
                    break

            if all_mats:
                combined_dir = os.path.join(suite2p_path, 'combined')
                combined_mat = os.path.join(combined_dir, 'Fall.mat')
                if not os.path.exists(combined_dir) or not os.path.exists(combined_mat):
                    print("All analysis-plane .mat files present but combined folder missing — creating combined outputs now.")
                    try:
                        build_combined_without_flyback(suite2p_path, analysis_planes)
                        print("Combined outputs created.\n")
                    except Exception as e:
                        print(f"Warning: failed to create combined outputs: {e}")
                        traceback.print_exc()
                    # after attempting to create combined outputs, skip further processing of this folder
                    continue
                else:
                    print("Binary files, plane .mat outputs, and combined Fall.mat exist. Skipping this folder.\n")
                    continue
            else:
                print("Binary files exist for all planes but some .mat outputs are missing; continuing to (re)process missing planes.")
    
    # Run Suite2P
    print(f"Data path: {folder}")
    print(f"Planes: {nplanes}, Channels: {nchannels}\n")
    # quick writable test: try creating a temp file in the save folder
    def can_write_test(path):
        try:
            with tempfile.NamedTemporaryFile(dir=path, delete=True) as tf:
                tf.write(b"x")
            return True
        except Exception:
            return False

    if not os.path.exists(folder):
        print(f"Folder missing (unexpected): {folder}\n")
        continue

    if not can_write_test(folder):
        print(f"\nPermission problem: cannot create files in {folder}")
        print("Close MATLAB/other programs or check network share permissions and retry.\n")
        continue

    # Detect existing .mat files per plane and skip those planes
    suite2p_dir = os.path.join(folder, 'suite2p')
    planes_to_run = list(analysis_planes)
    if os.path.exists(suite2p_dir):
        planes_to_run = []
        for plane in analysis_planes:
            plane_dir = os.path.join(suite2p_dir, f'plane{plane}')
            mat_files = glob.glob(os.path.join(plane_dir, '*.mat')) if os.path.exists(plane_dir) else []
            if mat_files:
                print(f" Found existing .mat for plane{plane}, skipping this plane")
            else:
                planes_to_run.append(plane)

    if not planes_to_run:
        print(f"All planes already have .mat files — skipping folder {folder}\n")
        continue

    # If we need to run all analysis planes, call run_s2p without built-in combine.
    did_run_full = False
    if len(planes_to_run) == len(analysis_planes):
        max_attempts = 3
        for attempt in range(1, max_attempts+1):
            try:
                ops_no_combine = ops.copy()
                ops_no_combine.setdefault("io", {})
                ops_no_combine["io"] = ops_no_combine["io"].copy()
                ops_no_combine["io"]["combined"] = False
                run_s2p(settings=ops_no_combine, db=db)
                print(f"\nCompleted folder {i}/{len(data_folders)}\n")
                did_run_full = True
                break
            except Exception as e:
                err_str = str(e)
                print(f"\nError processing {folder} (attempt {attempt}/{max_attempts}): {e}")
                traceback.print_exc()
                if isinstance(e, PermissionError) or 'Permission denied' in err_str:
                    if attempt < max_attempts:
                        print("Permission denied detected — waiting 5s then retrying...")
                        time.sleep(5)
                        continue
                    else:
                        print("Permission error persisted after retries — skipping folder.\n")
                else:
                    print("Non-permission error — skipping folder.\n")
                break
    else:
        # Run only missing planes using run_plane
        for plane in planes_to_run:
            print(f"Running plane {plane} for folder {folder}")
            db_plane = db.copy()
            db_plane['iplane'] = plane

            # If a plane-specific db.npy exists from conversion, load it
            plane_db_path = os.path.join(suite2p_dir, f'plane{plane}', 'db.npy')
            if os.path.exists(plane_db_path):
                try:
                    loaded_db = np.load(plane_db_path, allow_pickle=True).item()
                    # prefer loaded values for frame counts and file_list
                    db_plane.update(loaded_db)
                    print(f"  Loaded plane db from {plane_db_path}")
                except Exception as e:
                    print(f"  Warning: failed to load {plane_db_path}: {e}")

            # ensure minimal required fields
            if 'nframes' not in db_plane:
                print(f"  Missing 'nframes' for plane {plane} — skipping plane")
                continue

            try:
                run_plane(db=db_plane, settings=ops)
                print(f"Completed plane {plane} for folder {folder}")
            except Exception as e:
                print(f"Error processing plane {plane} in {folder}: {e}")
                traceback.print_exc()
                print("Continuing to next plane...\n")

        # After per-plane processing, create combined outputs from the analysis planes only.
        if not did_run_full:
            try:
                print("Ensuring combined outputs from analysis planes only...")
                build_combined_without_flyback(suite2p_path, analysis_planes)
                print(f"\nCombined outputs created/updated for folder {folder}\n")
            except Exception as e:
                print(f"Warning: failed to create combined outputs: {e}")
                traceback.print_exc()

print("="*70)
print("ALL FOLDERS PROCESSED!")
print("="*70)