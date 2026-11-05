"""Sum the 5 x/y-registered planes of a Suite3D job into a single-plane movie (BigTIFF pieces) that plain suite2p can run on.

    python make_zprojection_tif.py <s3d_job_dir> <out_dir> [--volumes-per-tif 1000]

Why: on RB24_261104_1651 the five planes (nominally 3 um apart) are almost the same optical section (mean images correlate
0.94-0.97), so a plane-sum has ~sqrt(5) x the shot-noise SNR of one plane and is insensitive to small axial shifts, whereas
z-registration has no depth information to work with. The registered movie of Suite3D is registered in x/y only
(apply_z_shift is off). Output dtype int16 (sum of the planes, max pixel ~5000 counts, no saturation).
"""
import argparse
import os

import numpy as np
import tifffile

ap = argparse.ArgumentParser()
ap.add_argument("job")
ap.add_argument("out")
ap.add_argument("--volumes-per-tif", type=int, default=1000)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
d = os.path.join(a.job, "registered_fused_data")
files = sorted(f for f in os.listdir(d) if f.startswith("fused_reg_data") and f.endswith(".npy"))

buf, n_out, total = [], 0, 0


def flush():
    global buf, n_out
    if not buf:
        return
    mov = np.concatenate(buf)
    tifffile.imwrite(os.path.join(a.out, f"zsum_{n_out:03d}.tif"), mov, bigtiff=True)
    n_out += 1
    buf = []


have = 0
for i, f in enumerate(files):
    arr = np.load(os.path.join(d, f), mmap_mode="r")             # (nz, nt, ly, lx)
    s = np.zeros(arr.shape[1:], dtype=np.int32)
    for p in range(arr.shape[0]):
        s += np.asarray(arr[p], dtype=np.int32)
    buf.append(np.clip(s, -32768, 32767).astype(np.int16))
    have += s.shape[0]
    total += s.shape[0]
    if have >= a.volumes_per_tif:
        flush()
        have = 0
flush()
print(f"{total} volumes written as {n_out} tifs to {a.out}")
