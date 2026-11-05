"""2x2-bin Suite3D's x/y-registered movie (registered_fused_data/fused_reg_data*.npy, (planes, T, Ly, Lx)) into
<out>/reg_bin2.npy  int16 (T, planes, Ly/2, Lx/2) pixel SUMS.   python bin_registered.py <s3d_job_dir> <out_dir>"""
import os
import sys

import numpy as np

job, out = sys.argv[1:3]
d = os.path.join(job, "registered_fused_data")
files = sorted(f for f in os.listdir(d) if f.startswith("fused_reg_data") and f.endswith(".npy"))
parts = []
for i, f in enumerate(files):
    a = np.load(os.path.join(d, f), mmap_mode="r")             # (nz, nt, ly, lx)
    nz, nt, ly, lx = a.shape
    s = np.zeros((nt, nz, ly // 2, lx // 2), dtype=np.int32)
    for p in range(nz):
        blk = np.asarray(a[p], dtype=np.int32).reshape(nt, ly // 2, 2, lx // 2, 2)
        s[:, p] = blk.sum(axis=(2, 4))
    parts.append(np.clip(s, -32768, 32767).astype(np.int16))
    if (i + 1) % 10 == 0:
        print(f"{i + 1}/{len(files)}", flush=True)
mov = np.concatenate(parts)
np.save(os.path.join(out, "reg_bin2.npy"), mov)
print("registered binned movie", mov.shape)
