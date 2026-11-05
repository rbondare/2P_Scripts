"""One pass over a 5-plane + flyback ScanImage recording -> a small movie for registration-free motion checks.

    python prep_binned_movie.py <tif_folder> <out_dir> [--bin 2] [--n-planes 5] [--frames-per-volume 6] [--min-volumes 10]

Per tif (= one acquisition; tifs with fewer than --min-volumes volumes are rig stubs and are skipped, a trailing partial
volume is dropped, the flyback frame(s) after the planes are ignored) it writes
    <out_dir>/movie_bin<b>.npy    int16 (T, planes, Ly/b, Lx/b), pixel SUMS over b x b pixels (not means)
    <out_dir>/mean_img.npy        float32 (planes, Ly, Lx), time mean at full resolution
    <out_dir>/meta.json           tif names, volumes per tif, volume index of each tif start
"""
import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import tifffile


def count_pages(path):
    with tifffile.TiffFile(path) as t:
        return len(t.pages)


def process(args):
    path, nv, fpv, nz, b = args
    a = tifffile.imread(path, key=range(0, nv * fpv))
    a = a.reshape(nv, fpv, *a.shape[1:])[:, :nz]                    # drop the flyback frame(s)
    mean = a.astype(np.float64).mean(axis=0).astype(np.float32)
    nvv, _, ly, lx = a.shape
    s = a.reshape(nvv, nz, ly // b, b, lx // b, b).astype(np.int32).sum(axis=(3, 5))
    return np.clip(s, -32768, 32767).astype(np.int16), mean


def main():
    p = argparse.ArgumentParser()
    p.add_argument("folder")
    p.add_argument("out")
    p.add_argument("--bin", type=int, default=2)
    p.add_argument("--n-planes", type=int, default=5)
    p.add_argument("--frames-per-volume", type=int, default=6)
    p.add_argument("--min-volumes", type=int, default=10)
    p.add_argument("--workers", type=int, default=8)
    a = p.parse_args()

    tifs = sorted(str(f) for f in Path(a.folder).glob("*.tif"))
    fpv = a.frames_per_volume
    jobs, names, vols = [], [], []
    for t in tifs:
        nv = count_pages(t) // fpv
        if nv < a.min_volumes:
            print("skip stub", os.path.basename(t), nv, "volumes")
            continue
        jobs.append((t, nv, fpv, a.n_planes, a.bin))
        names.append(os.path.basename(t))
        vols.append(nv)
    print(len(jobs), "tifs,", sum(vols), "volumes", flush=True)

    with ProcessPoolExecutor(a.workers) as ex:
        res = list(ex.map(process, jobs))
    movie = np.concatenate([r[0] for r in res])
    w = np.array(vols, dtype=np.float64)
    mean_img = sum(r[1] * wi for r, wi in zip(res, w)) / w.sum()

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    np.save(out / f"movie_bin{a.bin}.npy", movie)
    np.save(out / "mean_img.npy", mean_img.astype(np.float32))
    json.dump({"tifs": names, "volumes": vols, "start": np.cumsum([0] + vols[:-1]).tolist(), "bin": a.bin,
               "n_planes": a.n_planes, "frames_per_volume": fpv}, open(out / "meta.json", "w"), indent=1)
    print("movie", movie.shape, movie.dtype, "saturated sums:", int((np.abs(movie) >= 32767).sum()))


if __name__ == "__main__":
    main()
