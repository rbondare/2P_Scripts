#!/usr/bin/env python3
"""
fix_duplicate_merges.py  (v2 - corrected)

Restores from *_backup.npy files, then removes ONLY duplicate/subset
merged ROIs from stat/F/Fneu/spks/iscell.

Constituent ROIs are kept in place (so imerge indices stay valid).
inmerge on constituent ROIs is updated to point to the kept merged ROI.

Usage:
    conda activate suite2p_venv
    python fix_duplicate_merges.py "Z:/...path.../suite2p/plane0"
"""

import sys
import os
import shutil
import numpy as np


def main():
    if len(sys.argv) < 2:
        print("Usage: python fix_duplicate_merges.py <plane_dir>")
        sys.exit(1)

    plane_dir = sys.argv[1]
    fnames = ("stat.npy", "F.npy", "Fneu.npy", "spks.npy", "iscell.npy")

    # ── Step 1: restore from backups ──────────────────────────────────────
    print("Restoring from backups...")
    for fname in fnames:
        backup = os.path.join(plane_dir, fname.replace(".npy", "_backup.npy"))
        original = os.path.join(plane_dir, fname)
        if not os.path.exists(backup):
            print(f"  ERROR: backup not found: {backup}")
            sys.exit(1)
        shutil.copy2(backup, original)
        print(f"  restored {fname}")

    # ── Step 2: load restored files ────────────────────────────────────────
    def load(fname):
        return np.load(os.path.join(plane_dir, fname), allow_pickle=True)

    stat   = load("stat.npy")
    F      = load("F.npy")
    Fneu   = load("Fneu.npy")
    spks   = load("spks.npy")
    iscell = load("iscell.npy")
    N = len(stat)
    print(f"\nLoaded {N} ROIs, {F.shape[1]} frames.")

    # ── Step 3: identify merged ROIs ───────────────────────────────────────
    merged_idx    = []   # indices in stat
    merged_consts = []   # frozenset of constituent indices

    for n in range(N):
        imerge = stat[n].get("imerge", [])
        if imerge is not None and len(imerge) > 0:
            merged_idx.append(n)
            merged_consts.append(frozenset(int(k) for k in imerge))

    print(f"\nFound {len(merged_idx)} merged ROIs:")
    for idx, consts in zip(merged_idx, merged_consts):
        print(f"  stat[{idx}]  constituents={sorted(consts)}")

    if not merged_idx:
        print("No merged ROIs – nothing to do.")
        return

    # ── Step 4: mark duplicates/subsets for removal ────────────────────────
    # Keep a merged ROI if no other merged ROI has a superset of its constituents.
    # For exact duplicates, keep the one with the higher stat index (created later).
    to_remove_merged = set()
    for i in range(len(merged_idx)):
        for j in range(len(merged_idx)):
            if i == j:
                continue
            ci, cj = merged_consts[i], merged_consts[j]
            if ci < cj:                                          # strict subset
                to_remove_merged.add(merged_idx[i])
            elif ci == cj and merged_idx[i] < merged_idx[j]:   # exact duplicate, keep later
                to_remove_merged.add(merged_idx[i])

    to_keep_merged = [idx for idx in merged_idx if idx not in to_remove_merged]

    print(f"\nKeeping {len(to_keep_merged)} merged ROIs:")
    for idx in to_keep_merged:
        ci = merged_idx.index(idx)
        print(f"  stat[{idx}]  constituents={sorted(merged_consts[ci])}")
    print(f"Removing {len(to_remove_merged)} duplicate/subset merged ROIs: "
          f"{sorted(to_remove_merged)}")

    # ── Step 5: remove ONLY duplicate merged ROIs (not constituents) ───────
    to_keep = [n for n in range(N) if n not in to_remove_merged]

    stat_clean   = stat[to_keep]
    F_clean      = F[to_keep]
    Fneu_clean   = Fneu[to_keep]
    spks_clean   = spks[to_keep]
    iscell_clean = iscell[to_keep]

    # ── Step 6: rebuild old→new index map and update inmerge ───────────────
    old_to_new = {}
    for new_idx, old_idx in enumerate(to_keep):
        old_to_new[old_idx] = new_idx

    # For each constituent ROI, update inmerge to point to the surviving merged ROI
    # Build: old_constituent_idx → old_merged_ROI_idx (for kept merged ROIs only)
    constituent_to_merged = {}
    for kept_idx in to_keep_merged:
        ci = merged_idx.index(kept_idx)
        for c in merged_consts[ci]:
            constituent_to_merged[c] = kept_idx

    for new_n, old_n in enumerate(to_keep):
        s = stat_clean[new_n]
        imerge = s.get("imerge", [])
        if imerge is not None and len(imerge) > 0:
            # This is a kept merged ROI — no change to imerge (indices are constituent
            # ROI positions which are unchanged since we only removed merged ROIs)
            s["inmerge"] = 0   # merged ROIs themselves have inmerge=0
        else:
            # Regular or constituent ROI: update inmerge if it has one
            if old_n in constituent_to_merged:
                old_merged = constituent_to_merged[old_n]
                s["inmerge"] = old_to_new[old_merged]
            else:
                s["inmerge"] = 0

    # ── Step 7: save ───────────────────────────────────────────────────────
    np.save(os.path.join(plane_dir, "stat.npy"),   stat_clean)
    np.save(os.path.join(plane_dir, "F.npy"),      F_clean)
    np.save(os.path.join(plane_dir, "Fneu.npy"),   Fneu_clean)
    np.save(os.path.join(plane_dir, "spks.npy"),   spks_clean)
    np.save(os.path.join(plane_dir, "iscell.npy"), iscell_clean)

    print(f"\nDone. Cleaned stat has {len(stat_clean)} ROIs "
          f"({N - len(stat_clean)} duplicate merged ROIs removed).")
    print("Constituent ROIs kept in place — imerge indices are still valid.")


if __name__ == "__main__":
    main()
