"""Write the decisions of the review sheets into suite2p's iscell.npy (DRY RUN unless --apply).

usage:  python apply_curation.py [session ...]                 dry run: shows what would change, writes nothing
        python apply_curation.py --apply [session ...]         backs up iscell.npy, sets iscell = 1 for every ROI marked 'accept', writes a change log
        python apply_curation.py --apply --export-mat [...]    additionally backs up Fall.mat and rebuilds it from the npy files (suite2p_analysis/export_suite2p_to_mat.py)
        python apply_curation.py --revert [--apply] <session>  restores the most recent backup of iscell.npy (dry run without --apply)
Input: plots/summary_figures/rejected_roi_review/labels_review_<session>.csv (column 'decision': accept / reject / unsure; blank = untouched).
Only ROIs that are currently iscell 0 and are not hidden constituents of a merge (inmerge > 0) can be accepted; 'reject' and 'unsure' leave iscell as it is. The probability column of iscell.npy is not changed.
Backups are never overwritten: iscell_backup_original.npy (first run only) and iscell_backup_<timestamp>.npy in the same suite2p folder.

After --apply the downstream chain has to be rerun (nothing here does it for you):
  1. Fall.mat holds its own copy of iscell: rebuild it (--export-mat, or python suite2p_analysis/export_suite2p_to_mat.py <plane_dir>)
  2. delete the S2Presult intermediate of the session (it caches the old Fall) and re-aggregate the session in MATLAB
  3. re-export the session (cross_session/export_xs.m), rerun the cross-session analyses and paper_figures/make_all.py --clean
"""
import os, sys, csv, shutil, time, glob, numpy as np
from fig_common import *
import xs_common as xc

REV = os.path.join(OUT, 'rejected_roi_review')


def read_sheet(name):
    f = os.path.join(REV, 'labels_review_%s.csv' % name)
    if not os.path.exists(f): return None
    rows = []
    for r in csv.reader(open(f, newline='')):
        rows.append(r)
    head = rows[0]; k_idx = head.index('suite2p_roi_index'); k_dec = [i for i, h in enumerate(head) if h.startswith('decision')][0]; k_com = head.index('comment') if 'comment' in head else None
    return [(int(r[k_idx]), r[k_dec].strip().lower(), r[k_com] if k_com is not None and len(r) > k_com else '') for r in rows[1:] if len(r) > k_dec and r[k_idx].strip()]


def latest_backup(p):
    b = sorted(glob.glob(os.path.join(p, 'iscell_backup_2*.npy'))); return b[-1] if b else None


def main():
    args = sys.argv[1:]; apply = '--apply' in args; export = '--export-mat' in args; revert = '--revert' in args; names = [a for a in args if not a.startswith('--')]
    if not names: names = [os.path.basename(f)[len('labels_review_'):-4] for f in sorted(glob.glob(os.path.join(REV, 'labels_review_*.csv')))]
    ts = time.strftime('%Y%m%d_%H%M%S'); print('%s run (%s)' % ('APPLY' if apply else 'DRY', time.strftime('%Y-%m-%d %H:%M')))
    for name in names:
        p = os.path.join(xc.DATA_DIR, name, 'suite2p', 'plane0'); f_is = os.path.join(p, 'iscell.npy')
        if revert:
            b = latest_backup(p); print('%-26s revert -> %s' % (name, b or 'NO BACKUP FOUND'))
            if apply and b: shutil.copy2(f_is, os.path.join(p, 'iscell_before_revert_%s.npy' % ts)); shutil.copy2(b, f_is); print('   restored; rebuild Fall.mat and re-aggregate')
            continue
        sheet = read_sheet(name)
        if sheet is None: print('%-26s no decision sheet' % name); continue
        iscell = np.load(f_is); stat = np.load(os.path.join(p, 'stat.npy'), allow_pickle=True); new = iscell.copy(); acc, rej, uns, blank, bad, already = [], 0, 0, 0, [], []
        for idx, dec, com in sheet:
            if dec in ('', 'nan'): blank += 1; continue
            if dec == 'reject': rej += 1; continue
            if dec == 'unsure': uns += 1; continue
            if dec != 'accept': bad.append((idx, 'unknown decision "%s"' % dec)); continue
            if idx < 0 or idx >= len(iscell): bad.append((idx, 'index out of range')); continue
            if iscell[idx, 0] > 0: already.append(idx); continue
            if float(stat[idx].get('inmerge', 0)) > 0: bad.append((idx, 'hidden constituent of a merged ROI')); continue
            acc.append(idx); new[idx, 0] = 1
        print('%-26s accept %3d (iscell %d -> %d) | reject %3d, unsure %3d, undecided %3d | already accepted %d | refused %d%s' % (name, len(acc), int(iscell[:, 0].sum()), int(new[:, 0].sum()), rej, uns, blank, len(already), len(bad), ''.join('\n     roi %d: %s' % b for b in bad)))
        if not apply or not acc: continue
        if not os.path.exists(os.path.join(p, 'iscell_backup_original.npy')): shutil.copy2(f_is, os.path.join(p, 'iscell_backup_original.npy'))
        shutil.copy2(f_is, os.path.join(p, 'iscell_backup_%s.npy' % ts)); np.save(f_is, new)
        with open(os.path.join(p, 'iscell_review_log_%s.csv' % ts), 'w', newline='') as f:
            w = csv.writer(f); w.writerow(['suite2p_roi_index', 'old_iscell', 'new_iscell']); [w.writerow([i, int(iscell[i, 0]), 1]) for i in acc]
        print('   written %s (backup iscell_backup_%s.npy)' % (f_is, ts))
        if export:
            sys.path.insert(0, os.path.join(REPO, 'suite2p_analysis')); from export_suite2p_to_mat import export_plane
            fm = os.path.join(p, 'Fall.mat')
            if os.path.exists(fm): shutil.copy2(fm, os.path.join(p, 'Fall_backup_%s.mat' % ts))
            export_plane(p)
    if not apply and not revert: print('\n(dry run: nothing written; add --apply to write)')
    if apply and not revert: print('\nNext: rebuild Fall.mat if you did not use --export-mat, delete the S2Presult intermediate, re-aggregate in MATLAB, re-export (export_xs.m), rerun make_all.py --clean.')


if __name__ == '__main__':
    main()
