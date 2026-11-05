"""Compare a Suite3D job on the synthetic recording with the known truth (injected z motion, lateral jitter, activity).

    python make_synth.py <raw_dir> [signal_scale]      # 1 = dim; 4 and 8 are the numbers quoted in the V1 docstring
    run_suite3d.bat <raw_dir> --out-dir <jobs> --job-id synth --n-planes 5 --n-ch-tif 6 --fs 5 --z-step-um 3 --pixel-um 0.5
    python eval_synth.py <jobs>/s3d-synth <raw_dir>    # in the suite3d env
"""
import importlib.util
import os
import sys
from pathlib import Path

import numpy as np
import tifffile

job_dir, raw_dir = sys.argv[1], sys.argv[2]
spec = importlib.util.spec_from_file_location(
    "s3dscript", str(Path(__file__).resolve().parent.parent / "run_suite3d_RB_axons_V1.py"))
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)

truth = np.load(os.path.join(raw_dir, "truth.npz"))
fs, dz = float(truth["fs"]), float(truth["dz"])
z_true = truth["z_um"]
NT = len(z_true)

# --- registration: measured z shift vs truth --------------------------------------------------------------------
reg = os.path.join(job_dir, "registered_fused_data")
offs = sorted(f for f in os.listdir(reg) if f.startswith("offsets"))
shifts = np.concatenate([np.load(os.path.join(reg, f), allow_pickle=True).item()["sub_pixel_shifts"] for f in offs])
print("volumes registered:", len(shifts), "truth:", NT)
zs = shifts[:, 0]
r_z = np.corrcoef(zs, z_true)[0, 1]
slope = np.polyfit(z_true, zs * dz, 1)[0]
print(f"measured z shift vs true z: r = {r_z:+.3f}, slope = {slope:+.2f} (um per um; -1 or +1 = perfect), "
      f"saturated at the window edge: {100 * np.isclose(np.abs(zs), 2.5).mean():.1f}%")
r_y = np.corrcoef(shifts[:, 1], truth["sy"])[0, 1]
r_x = np.corrcoef(shifts[:, 2], truth["sx"])[0, 1]
print(f"measured y/x shift vs true lateral jitter: r_y = {r_y:+.3f}, r_x = {r_x:+.3f}")

# --- ROIs ---------------------------------------------------------------------------------------------------------------
rois = os.path.join(job_dir, "rois")
F = np.load(os.path.join(rois, "F.npy"))
Fneu = np.load(os.path.join(rois, "Fneu.npy"))
stats = np.load(os.path.join(rois, "stats.npy"), allow_pickle=True)
print("ROIs:", F.shape, "stats:", len(stats))
dff = S.delta_f_over_f(F, Fneu, fs)
meds = np.array([s["med"] for s in stats], dtype=float)    # z, y, x
npix = np.array([len(s["lam"]) for s in stats])
nz_span = np.array([len(np.unique(s["coords"][0])) for s in stats])
print(f"ROI size: median {np.median(npix):.0f} voxels, median planes spanned {np.median(nz_span):.1f}")


def corr(a, b):
    a = a - a.mean()
    b = b - b.mean()
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12))


bz_plane = truth["bz"] / dz + 2
matched, r_act, r_zc = [], [], []
for i in range(len(truth["by"])):
    d = np.hypot(meds[:, 1] - truth["by"][i], meds[:, 2] - truth["bx"][i])
    ok = np.where((d < 4) & (np.abs(meds[:, 0] - bz_plane[i]) <= 1.6))[0]
    if len(ok) == 0:
        continue
    best = max(ok, key=lambda j: corr(dff[j], truth["dff"][i]))
    matched.append(i)
    r_act.append(corr(dff[best], truth["dff"][i]))
    r_zc.append(corr(dff[best], z_true))
print(f"true boutons found: {len(matched)}/{len(truth['by'])}   (ROIs not matching a true bouton: "
      f"{len(stats) - len(matched)} of {len(stats)})")
if matched:
    print(f"  r(extracted dF/F, true activity): median {np.median(r_act):.2f}  (10th pct {np.percentile(r_act, 10):.2f})")
    print(f"  r(extracted dF/F, true axial motion): median |r| {np.median(np.abs(r_zc)):.2f}")
else:
    print("  no ROI matches a true bouton")
    sys.exit(0)

# --- control: naive single-plane trace from the raw tifs at the true bouton position, no registration ----------
tifs = sorted(os.path.join(raw_dir, f) for f in os.listdir(raw_dir) if f.endswith(".tif"))
frames = []
for t in tifs:
    a = tifffile.imread(t)
    if a.shape[0] >= 6:
        frames.append(a[: (a.shape[0] // 6) * 6].reshape(-1, 6, *a.shape[1:])[:, :5])
raw = np.concatenate(frames)            # (NT, 5, y, x)
print("raw volumes:", raw.shape[0])
naive_r_act, naive_r_z = [], []
for i in matched:
    p = int(np.clip(round(bz_plane[i]), 0, 4))
    y, x = int(round(truth["by"][i])), int(round(truth["bx"][i]))
    tr = raw[:, p, y - 1:y + 2, x - 1:x + 2].mean(axis=(1, 2))[None].astype(np.float32)
    nb = raw[:, p, max(0, y - 12):y + 13, max(0, x - 12):x + 13].mean(axis=(1, 2))[None].astype(np.float32)
    d = S.delta_f_over_f(tr, nb, fs)[0]
    naive_r_act.append(corr(d, truth["dff"][i]))
    naive_r_z.append(corr(d, z_true))
print(f"control (raw, nearest plane, no registration, 3x3 px): r(true activity) median {np.median(naive_r_act):.2f}, "
      f"median |r(true z)| {np.median(np.abs(naive_r_z)):.2f}")
