"""Roll the per-codec on/off levels of face_motion.ON_OFF (H.264: on above 1.45 / off below 1.40, 2026-11-10) out to the saved results.

usage:  python h264_threshold_rollout.py stage      re-classify ALL sessions of DLC\\work from their saved ME traces (no decoding) into a
                                                    staging folder and compare with what is saved now (regression test: MPEG-4 sessions
                                                    must come out identical, H.264 sessions must differ in whisk_on ONLY)
        python h264_threshold_rollout.py apply      copy the staged H.264 results to every place that holds a copy, after backing up
                                                    the version that is there (whisker_results_v2 next to it, once)
        python h264_threshold_rollout.py sync-stale copy the work version of the MPEG-4 sessions to stale copies (Camera folders that
                                                    still hold the pre-front_top results, backup whisker_results_v1 next to them)

Copies of the whisker results of a session (all must end up identical to DLC\\work\\<s>\\whisker_results):
    DLC\\work\\<s>\\whisker_results, DLC\\work\\quality_reports\\<s> (json + png), Z:\\Group Members\\Rima\\preCamera\\<animal>\\<s>\\whisker_results
    (after the session was aggregated) or Z:\\Group Members\\Rima\\Camera\\<s>\\whisker_results (before; the aggregation reads this one).
The aggregated files and DLC_filtered.mat are refreshed afterwards in MATLAB (newAggregation/refresh_whisker_behavior.m).
"""
import os, sys, glob, shutil, hashlib, argparse
import numpy as np, pandas as pd

sys.path.insert(0, r'C:\Users\rbondarenko\projects\2P_Scripts\DLC\face_model')
import face_motion as fm

WORK = r'C:\Users\rbondarenko\DLC\work'
PRE = r'Z:\Group Members\Rima\preCamera'
CAM = r'Z:\Group Members\Rima\Camera'
STAGE = r'C:\Users\rbondarenko\DLC\whisker_final_stage'
SUF = ('_whisker.csv', '_whisker.json', '_whisker_qc.png')


def sessions():
    return sorted(d for d in os.listdir(WORK) if d.startswith('Animal') and os.path.isfile(os.path.join(WORK, d, 'whisker_results', d + '_whisker.csv')))


def md5(p):
    return hashlib.md5(open(p, 'rb').read()).hexdigest() if os.path.isfile(p) else None


def copies(s):
    """folders (other than work) that hold the whisker_results of session s"""
    a = s.split('_')[0]
    out = []
    for p in (os.path.join(PRE, a, s, 'whisker_results'), os.path.join(CAM, s, 'whisker_results')):
        if os.path.isdir(p):
            out.append(p)
    return out


def stage():
    os.makedirs(STAGE, exist_ok=True)
    rows = []
    for s in sessions():
        old = os.path.join(WORK, s, 'whisker_results', s + '_whisker.csv')
        video = glob.glob(os.path.join(WORK, s, 'camera', '*.mp4'))
        if not video:
            print(s, 'NO VIDEO in work\\camera, skipped')
            continue
        res = fm.reclassify(video[0], old, s)
        fm.write_outputs(res, os.path.join(STAGE, s))
        O = pd.read_csv(old)
        N = pd.read_csv(os.path.join(STAGE, s, s + '_whisker.csv'))
        same_cols = [c for c in O.columns if c != 'whisk_on' and O[c].equals(N[c])]
        diff_cols = [c for c in O.columns if c != 'whisk_on' and not O[c].equals(N[c])]
        on_o, on_n = O['whisk_on'], N['whisk_on']
        art = O['artifact'].values.astype(bool)
        nan_same = bool((on_o.isna().values == on_n.isna().values).all()) and bool((on_o.isna().values == art).all())
        chg = int((on_o.fillna(-1).values != on_n.fillna(-1).values).sum())
        ok = ~art
        rows.append(dict(session=s, codec=res['codec'], frames=len(O), other_cols_identical=not diff_cols, nan_pattern_same=nan_same,
                         on_changed_frames=chg, whisking_old=round(float(on_o[ok].mean()), 3), whisking_new=round(float(on_n[ok].mean()), 3),
                         bouts_old=len(fm.runs(on_o.fillna(0).values.astype(bool))), bouts_new=len(fm.runs(on_n.fillna(0).values.astype(bool))),
                         diff_cols=','.join(diff_cols)))
        r = rows[-1]
        print('%-24s %-5s other cols identical %-5s  whisk_on changed in %7d frames  whisking %.3f -> %.3f  bouts %4d -> %4d' % (
            s, r['codec'], r['other_cols_identical'], chg, r['whisking_old'], r['whisking_new'], r['bouts_old'], r['bouts_new']), flush=True)
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(STAGE, 'stage_comparison.csv'), index=False)
    print('\nMPEG-4 sessions: %d, all identical incl. whisk_on: %s' % ((T.codec == 'mpeg4').sum(), bool(((T.codec == 'mpeg4') & (T.on_changed_frames == 0) & T.other_cols_identical).sum() == (T.codec == 'mpeg4').sum())))
    h = T[T.codec == 'h264']
    print('H.264 sessions: %d, all other columns identical: %s, NaN pattern identical: %s' % (len(h), bool(h.other_cols_identical.all()), bool(h.nan_pattern_same.all())))


def version_tag(folder):
    """'v2' = front_top version of 2026-11-10 (json has whisker_box), 'v1' = the original two-box version: names the backup"""
    import json
    js = glob.glob(os.path.join(folder, '*_whisker.json'))
    return 'v2' if js and 'whisker_box' in json.load(open(js[0])) else 'v1'


def backup_and_copy(dst_dir, src_dir, tag):
    """copy the staged / work files into dst_dir after backing the folder up once as whisker_results_<tag> next to it"""
    bak = os.path.join(os.path.dirname(dst_dir), 'whisker_results_' + tag)
    if not os.path.isdir(bak):
        shutil.copytree(dst_dir, bak)
    for suf in SUF:
        for f in glob.glob(os.path.join(src_dir, '*' + suf)):
            shutil.copy2(f, os.path.join(dst_dir, os.path.basename(f)))


def apply():
    T = pd.read_csv(os.path.join(STAGE, 'stage_comparison.csv'))
    h = T[T.codec == 'h264']
    assert h.other_cols_identical.all() and h.nan_pattern_same.all(), 'stage comparison failed'
    for s in h.session:
        st = os.path.join(STAGE, s)
        wr = os.path.join(WORK, s, 'whisker_results')
        backup_and_copy(wr, st, 'v2')
        qd = os.path.join(WORK, 'quality_reports', s)
        if os.path.isdir(qd):
            for suf in ('_whisker.json', '_whisker_qc.png'):
                shutil.copy2(os.path.join(st, s + suf), os.path.join(qd, s + suf))
        for c in copies(s):
            backup_and_copy(c, st, version_tag(c))
        print(s, 'updated work +', len(copies(s)), 'copies', flush=True)


def sync_stale():
    n = 0
    for s in sessions():
        wr = os.path.join(WORK, s, 'whisker_results')
        for c in copies(s):
            if md5(os.path.join(c, s + '_whisker.csv')) != md5(os.path.join(wr, s + '_whisker.csv')):
                backup_and_copy(c, wr, version_tag(c))   # normally the original pre-front_top version -> whisker_results_v1
                n += 1
                print(s, 'stale copy updated:', c, flush=True)
    print('synced', n)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['stage', 'apply', 'sync-stale'])
    a = ap.parse_args()
    {'stage': stage, 'apply': apply, 'sync-stale': sync_stale}[a.cmd]()
