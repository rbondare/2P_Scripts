"""Raw 5-plane + flyback tifs -> x/y-registered (rigid only, NO z handling) outputs for the motion-vs-signal checks.

    python prep_rigid_registered.py <tif_folder> <s3d_job_dir> <out_dir> [--zsum-dir <dir>] [--bin 2]

Why not Suite3D's own registered movie: with 5 planes @ 3 um Suite3D's non-rigid z shifts are applied even with
apply_z_shift=False, and even with max_shift_nr z = 0 every block gets a fixed -1 plane shift with edge replication
(registered plane 0 == plane 1 in every volume, planes 1-4 = raw planes 0-3). Here only the rigid y/x shift per volume
(Suite3D registration offsets, sub_pixel_shifts[:, 1:]) is applied to all five planes, with bilinear interpolation.

Writes, in the same volume order as the Suite3D job (stub tifs < --min-volumes skipped, trailing partial volume dropped):
    <out>/regraw_bin2.npy   int16 (T, planes, Ly/b, Lx/b)  pixel SUMS of b x b pixels, rigid-registered
    <zsum-dir>/zsum_<k>.tif int16 (n, Ly, Lx)              sum of the five rigid-registered planes (for plain suite2p)
The sign of the shift is chosen automatically (the one that sharpens the mean image of 150 high-motion volumes).
"""
import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import tifffile
from scipy.ndimage import gaussian_filter, shift as ndshift


def load_shifts(job):
    d = os.path.join(job, "registered_fused_data")
    files = sorted(f for f in os.listdir(d) if f.startswith("offsets") and f.endswith(".npy"))
    return np.concatenate([np.load(os.path.join(d, f), allow_pickle=True).item()["sub_pixel_shifts"] for f in files])


def read_volumes(path, v0, nv, fpv, nz):
    a = tifffile.imread(path, key=range(v0 * fpv, (v0 + nv) * fpv))
    return a.reshape(nv, fpv, *a.shape[1:])[:, :nz].astype(np.float32)


def register(vols, sh, sign):
    out = np.empty_like(vols)
    for i in range(len(vols)):
        out[i] = ndshift(vols[i], (0, sign * sh[i, 1], sign * sh[i, 2]), order=1, mode="nearest")
    return out


def work(args):
    path, nv, sh, fpv, nz, b, sign = args
    vols = register(read_volumes(path, 0, nv, fpv, nz), sh, sign)
    zsum = np.clip(np.rint(vols.sum(1)), -32768, 32767).astype(np.int16)
    n, _, ly, lx = vols.shape
    s = vols.reshape(n, nz, ly // b, b, lx // b, b).sum(axis=(3, 5))
    return np.clip(np.rint(s), -32768, 32767).astype(np.int16), zsum


def sharpness(img):
    return float((img - gaussian_filter(img, 12)).std())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("job")
    ap.add_argument("out")
    ap.add_argument("--zsum-dir", default=None)
    ap.add_argument("--bin", type=int, default=2)
    ap.add_argument("--n-planes", type=int, default=5)
    ap.add_argument("--frames-per-volume", type=int, default=6)
    ap.add_argument("--min-volumes", type=int, default=10)
    ap.add_argument("--keep-stubs", action="store_true",
                    help="also output the rig-stub tifs (< --min-volumes volumes), with zero shift, so that the movie has "
                         "one frame per volume of EVERY tif, as the aggregation (tif-header timestamps) expects")
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    fpv, nz = a.frames_per_volume, a.n_planes
    S = load_shifts(a.job)

    tifs = sorted(str(f) for f in Path(a.folder).glob("*.tif"))
    jobs, names, vols, start, n_stub = [], [], [], 0, 0
    for t in tifs:
        with tifffile.TiffFile(t) as tf:
            nv = len(tf.pages) // fpv
        if nv < a.min_volumes:
            if not a.keep_stubs:
                continue
            n_stub += nv
            names.append(os.path.basename(t))
            vols.append(nv)
            jobs.append([t, nv, np.zeros((nv, 3)), fpv, nz, a.bin, None])
            continue
        names.append(os.path.basename(t))
        vols.append(nv)
        jobs.append([t, nv, S[start:start + nv], fpv, nz, a.bin, None])
        start += nv
    assert start == len(S), f"{start} volumes in the tifs vs {len(S)} Suite3D shifts"
    if n_stub:
        print(f"{n_stub} stub volumes kept with zero shift", flush=True)

    # sign test on the tif with the largest lateral excursion
    k = int(np.argmax([np.abs(j[2][:, 1:]).max() for j in jobs]))
    t, nv, sh = jobs[k][0], jobs[k][1], jobs[k][2]
    n_test = min(150, nv)
    raw = read_volumes(t, 0, n_test, fpv, nz)
    best = {}
    for sign in (1.0, -1.0):
        best[sign] = np.mean([sharpness(register(raw, sh[:n_test], sign)[:, p].mean(0)) for p in range(nz)])
    sign = max(best, key=best.get)
    print(f"shift sign {sign:+.0f} (sharpness {best[1.0]:.3f} vs {best[-1.0]:.3f}) on {os.path.basename(t)}, max |shift| there "
          f"{np.abs(sh[:, 1:]).max():.1f} px", flush=True)
    for j in jobs:
        j[6] = sign

    with ProcessPoolExecutor(a.workers) as ex:
        res = list(ex.map(work, jobs))
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    movie = np.concatenate([r[0] for r in res])
    np.save(out / f"regraw_bin{a.bin}.npy", movie)
    if a.zsum_dir:
        os.makedirs(a.zsum_dir, exist_ok=True)
        zs = np.concatenate([r[1] for r in res])
        for i, s0 in enumerate(range(0, len(zs), 1000)):
            tifffile.imwrite(os.path.join(a.zsum_dir, f"zsum_{i:03d}.tif"), zs[s0:s0 + 1000], bigtiff=True)
    json.dump({"tifs": names, "volumes": vols, "sign": sign}, open(out / "regraw_meta.json", "w"), indent=1)
    print("regraw movie", movie.shape, "; zsum", (len(zs), *zs.shape[1:]) if a.zsum_dir else "not written")


if __name__ == "__main__":
    main()
