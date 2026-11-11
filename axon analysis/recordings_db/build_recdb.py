"""Stage B of the recordings database: assemble the final tables (sessions / blocks / rois / axons / animals / notes) for RB22-24.

Run (after export_session_meta.m, see build_recdb.m):  miniforge3\\envs\\suite2p_venv\\python.exe build_recdb.py [--no-fov] [--no-class] [--snapshot]
Inputs (all read-only):
  <DB>\\_build\\m_sessions.csv, m_blocks.csv          from export_session_meta.m (stimulus blocks, aggregated-file QC)
  Z:\\...\\DATA_2P\\<session>\\*.tif (first header), suite2p\\plane0\\{iscell,stat,ops}.npy
  axon analysis\\plots\\cross_session\\session_qc.csv, plots\\fov_matching\\{fov_pose,master_axons,master_presence}.csv
  C:\\Users\\rbondarenko\\2P_cache\\xs\\<session>_xs.mat (only for the ROI signal/noise classes of the 19 exported sessions)
Outputs in <DB> = Z:\\Group Members\\Rima\\Recordings_DB:  sessions.csv blocks.csv rois.csv axons.csv animals.csv notes.csv fov\\<session>.png index.html README.md
Rules: auto columns are regenerated every run; notes.csv and animals.csv are hand-edited files that are created once and NEVER overwritten
(new sessions are appended to notes.csv with seeded values; existing cells are untouched).
"""
import os, re, sys, csv, ast, glob, json, shutil, argparse, datetime
import numpy as np

REPO = r'C:\Users\rbondarenko\projects\2P_Scripts'
AX = os.path.join(REPO, 'axon analysis')
ZROOT = r'Z:\joeschgrp\Group Members\Rima' if os.path.isdir(r'Z:\joeschgrp\Group Members\Rima') else r'Z:\Group Members\Rima'
DATA = os.path.join(ZROOT, 'DATA_2P')
AGG = os.path.join(ZROOT, 'Aggregated')
CAMERA, PRECAMERA, DLC = os.path.join(ZROOT, 'Camera'), os.path.join(ZROOT, 'preCamera'), os.path.join(ZROOT, 'DLCData')
DB = os.environ.get('RECDB_DIR') or os.path.join(ZROOT, 'Recordings_DB')        # RECDB_DIR / RECDB_BUILD: test overrides, normally unset
BUILD = os.environ.get('RECDB_BUILD') or os.path.join(DB, '_build')
SNAP = os.path.join(AX, 'recordings_db', 'snapshot')
XS_DIR = r'C:\Users\rbondarenko\2P_cache\xs'
QC_PATH = os.path.join(AX, 'plots', 'cross_session', 'session_qc.csv')
FOVM = os.path.join(AX, 'plots', 'fov_matching')
SESSION_RE = re.compile(r'^Animal(RB2[234])_(\d{6})_(\d{4})$')

# ---------------------------------------------------------------------------------------------------- small helpers
def read_csv(path):
    if not os.path.isfile(path): return []
    with open(path, newline='', encoding='utf-8-sig') as f: return list(csv.DictReader(f))


def write_csv(path, rows, cols):
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore'); w.writeheader()
        for r in rows: w.writerow({c: fmt(r.get(c)) for c in cols})


def fmt(v):
    if v is None: return ''
    if isinstance(v, (bool, np.bool_)): return 'yes' if v else 'no'
    if isinstance(v, (float, np.floating)):
        if not np.isfinite(v): return ''
        return ('%.6g' % v)
    return str(v)


def num(x):
    try:
        v = float(x); return v if np.isfinite(v) else None
    except (TypeError, ValueError): return None


def txt(x):
    return '' if x is None or str(x).strip().lower() == 'nan' else str(x)


def mtime_str(p):
    return datetime.datetime.fromtimestamp(os.path.getmtime(p)).strftime('%Y-%m-%d %H:%M:%S') if os.path.exists(p) else ''


def short(nm):
    m = SESSION_RE.match(nm); d, t = m.group(2), m.group(3)
    return '%s %s-%s %s' % (m.group(1), d[2:4], d[4:6], t)


# ---------------------------------------------------------------------------------------------------- hard-coded lists of the analysis code -> seeds (read with ast, nothing imported)
def literal_from(path, name):
    """value of the module-level assignment `name = <literal>` in a python file (None if absent / not a literal)"""
    try:
        tree = ast.parse(open(path, encoding='utf-8').read())
        for node in tree.body:
            if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == name for t in node.targets):
                return ast.literal_eval(node.value)
    except Exception as e:
        print('  seed %s from %s not read (%s)' % (name, os.path.basename(path), e))
    return None


def seeds():
    fc = os.path.join(AX, 'paper_figures', 'fig_common.py')
    flag = literal_from(fc, 'FLAG') or {}
    task = literal_from(fc, 'TASK_SESSIONS') or set()
    grey = literal_from(fc, 'GREY_SESSIONS') or set()
    excl = literal_from(os.path.join(AX, 'fov_matching', 'fov_common.py'), 'EXCLUDE_ROI_SESSIONS') or set()
    return dict(flag=dict(flag), task=set(task), grey=set(grey), exclude_roi=set(excl))


# ---------------------------------------------------------------------------------------------------- acquisition geometry from the first TIFF header (all sessions)
def read_header(nm, cache):
    files = sorted(glob.glob(os.path.join(DATA, nm, '*.tif')))
    if not files: return {'n_tifs': 0}
    key = '%s|%d|%.0f' % (os.path.basename(files[0]), len(files), os.path.getmtime(files[0]))
    if nm in cache and cache[nm].get('key') == key: return cache[nm]['val']
    import tifffile
    val = {'n_tifs': len(files), 'tif_first_mtime': mtime_str(files[0]), 'tif_last_mtime': mtime_str(files[-1])}
    fd = None
    for f in files[:3]:                         # the first header can be truncated (1540): fall back to the next acquisition
        try:
            with tifffile.TiffFile(f) as t: fd = t.scanimage_metadata['FrameData']
            break
        except Exception: fd = None
    if fd is None: val['acq_note'] = 'TIFF header unreadable'
    else:
        try:
            fov = np.array(fd['SI.hRoiManager.imagingFovUm'], float); fx, fy = 2 * abs(fov[0][0]), 2 * abs(fov[0][1])
            nx, ny = int(fd['SI.hRoiManager.pixelsPerLine']), int(fd['SI.hRoiManager.linesPerFrame'])
            ax = fd['SI.hMotors.axesPosition']; sp = fd['SI.hMotors.samplePosition']
            ns = int(fd.get('SI.hStackManager.numSlices', 1))
            val.update(zoom=float(fd['SI.hRoiManager.scanZoomFactor']), fov_x_um=round(fx, 1), fov_y_um=round(fy, 1), nx=nx, ny=ny,
                       um_per_px_x=round(fx / nx, 4), um_per_px_y=round(fy / ny, 4), axes_x=ax[0], axes_y=ax[1], axes_z=ax[2],
                       sample_dx=round(sp[0] - ax[0], 1), sample_dy=round(sp[1] - ax[1], 1), fastz_pos=fd['SI.hFastZ.position'],
                       power=fd['SI.hBeams.powers'][0], n_slices=ns, frame_rate=round(float(fd['SI.hRoiManager.scanFrameRate']), 3))
        except Exception as e:
            val['acq_note'] = 'header fields missing: %s' % e
    cache[nm] = {'key': key, 'val': val}
    return val


# ---------------------------------------------------------------------------------------------------- suite2p
def s2p_dir(nm): return os.path.join(DATA, nm, 'suite2p', 'plane0')


def load_s2p(nm):
    d = s2p_dir(nm); out = {}
    if os.path.isfile(os.path.join(d, 'iscell.npy')):
        out['iscell'] = np.load(os.path.join(d, 'iscell.npy'))
        out['iscell_mtime'] = os.path.getmtime(os.path.join(d, 'iscell.npy'))
    if os.path.isfile(os.path.join(d, 'stat.npy')): out['stat'] = np.load(os.path.join(d, 'stat.npy'), allow_pickle=True)
    if os.path.isfile(os.path.join(d, 'ops.npy')): out['ops'] = np.load(os.path.join(d, 'ops.npy'), allow_pickle=True).item()
    return out


def curation_files(nm):
    d = s2p_dir(nm); names = []
    if not os.path.isdir(d): return ''
    for f in sorted(os.listdir(d)):
        if f.startswith('iscell_') or f.startswith('user_curated') or f.startswith('Fall_backup') or f.startswith('iscell_review_log'): names.append(f)
    return '; '.join(names)


# ---------------------------------------------------------------------------------------------------- ROI classes (signal-like / noise-like / neuropil-limited) of the 19 exported sessions
def roi_classes(nm):
    """array of 'signal' / 'noise' / 'neuropil' per aggregated ROI row, or None when the session has no xs export. Same criterion and code as cross_session/xs_common.py"""
    p = os.path.join(XS_DIR, nm + '_xs.mat')
    if not os.path.isfile(p): return None
    sys.path.insert(0, os.path.join(AX, 'cross_session'))
    import xs_common as xc
    X = xc.load_xs(nm)
    tca = np.asarray(X['tca'], float); fs = 1 / np.median(np.diff(tca))
    F = np.atleast_2d(X['F']).astype(np.float32); Fn = np.atleast_2d(X['Fneu']).astype(np.float32); dp = np.atleast_2d(X['dFF_pipe']).astype(np.float32)
    sig, neu, _ = xc.signal_like_v3(dp, F - xc.NEUROPIL * Fn, fs)
    return np.where(neu, 'neuropil', np.where(sig, 'signal', 'noise'))


# ---------------------------------------------------------------------------------------------------- FOV picture
def fov_png(nm, s2, acc_cls, out_path, um_per_px_x, um_per_px_y, title):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy import ndimage as ndi
    ops, stat, ic = s2['ops'], s2['stat'], s2['iscell']
    img = np.asarray(ops['meanImg'], float)
    yr, xr = ops.get('yrange'), ops.get('xrange')
    sub = img[int(yr[0]):int(yr[1]), int(xr[0]):int(xr[1])] if yr is not None and xr is not None else img
    lo, hi = np.percentile(sub, 1), np.percentile(sub, 99.7)
    img = np.clip((img - lo) / max(hi - lo, 1e-6), 0, 1) ** 0.8
    h, w = img.shape
    rgb = np.dstack([img, img, img])
    col = {'signal': (0.15, 0.95, 0.95), 'noise': (1.0, 0.85, 0.2), 'neuropil': (1.0, 0.45, 0.45), None: (0.3, 1.0, 0.3)}
    k = 0
    for i, st in enumerate(stat):
        if ic[i, 0] <= 0: continue
        m = np.zeros((h, w), bool); m[st['ypix'], st['xpix']] = True
        edge = m & ~ndi.binary_erosion(m)
        c = col[acc_cls[k] if acc_cls is not None and k < len(acc_cls) else None]
        rgb[edge] = c; k += 1
    ar = (w * um_per_px_x) / (h * um_per_px_y) if um_per_px_x and um_per_px_y else 1.0
    fig, ax = plt.subplots(figsize=(6.2, 6.2 / max(ar, 0.3) * (1.0) + 0.5), dpi=110)
    ax.imshow(rgb, extent=[0, w * (um_per_px_x or 1), h * (um_per_px_y or 1), 0], interpolation='nearest')
    ax.set_xlabel('um' if um_per_px_x else 'px'); ax.set_ylabel('um' if um_per_px_y else 'px'); ax.set_title(title, fontsize=8)
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(labelsize=7)
    legend = 'cyan signal-like | yellow noise-like | red neuropil-limited' if acc_cls is not None else 'green = accepted ROIs, signal/noise class not available'
    ax.text(0.995, 0.005, legend, transform=ax.transAxes, ha='right', va='bottom', fontsize=5.5, color='w')
    fig.tight_layout(); fig.savefig(out_path); plt.close(fig)


# ---------------------------------------------------------------------------------------------------- notes.csv / animals.csv (hand-edited, never overwritten)
NOTE_COLS = ['session', 'pool', 'flag_reason', 'curation', 'tags', 'quality', 'visible_activity', 'comments',
             'notion_planes', 'notion_quality', 'notion_visible_activity', 'notion_comments', 'notion_suite2p', 'notion_stimuli_camera', 'seeded_from']
ANIMAL_COLS = ['animal', 'cohort', 'population', 'indicator', 'genotype', 'sex', 'date_birth', 'date_injection', 'date_implantation', 'virus', 'comments']


def update_notes(sessions, seed, in_qc):
    path = os.path.join(DB, 'notes.csv')
    old = {r['session']: r for r in read_csv(path)}
    rows = dict(old); added = []
    for nm in sessions:
        if nm in rows: continue
        r = {c: '' for c in NOTE_COLS}; r['session'] = nm; src = []
        if nm in seed['flag']: r['pool'] = 'no'; r['flag_reason'] = seed['flag'][nm]; src.append('fig_common.FLAG')
        elif nm in in_qc: r['pool'] = 'yes'; src.append('session_qc.csv (analysed, not flagged)')
        if nm in seed['task']: r['tags'] = (r['tags'] + ';task_analysis').strip(';'); src.append('fig_common.TASK_SESSIONS')
        if nm in seed['grey']: r['tags'] = (r['tags'] + ';grey_recorded').strip(';'); src.append('fig_common.GREY_SESSIONS')
        if nm in seed['exclude_roi']:
            r['tags'] = (r['tags'] + ';rois_unusable').strip(';'); src.append('fov_common.EXCLUDE_ROI_SESSIONS')
            if not r['flag_reason']: r['flag_reason'] = 'ROIs unusable (horizontal bands over the whole field); used for registration only'
        r['seeded_from'] = ' + '.join(src); rows[nm] = r; added.append(nm)
    ordered = sorted(rows)
    cols = list(NOTE_COLS) + [c for c in (next(iter(old.values())).keys() if old else []) if c not in NOTE_COLS]
    write_csv(path, [rows[k] for k in ordered], cols)
    return rows, added


def ensure_animals(animals):
    path = os.path.join(DB, 'animals.csv')
    if os.path.isfile(path): return {r['animal']: r for r in read_csv(path)}
    rows = []
    for a in animals:
        rows.append({'animal': a, 'cohort': '1st cohort of axo-GCaMP6 (RB22-RB24)', 'population': 'LC-NE axons (user statement)', 'indicator': '', 'genotype': '', 'sex': '',
                     'date_birth': '', 'date_injection': '', 'date_implantation': '', 'virus': '',
                     'comments': 'FILL IN: AnimalInfo in the aggregated files holds the dialog defaults (Ai148/Rorb, m, GCaMP6f, placeholder dates), not real values'})
    write_csv(path, rows, ANIMAL_COLS)
    return {r['animal']: r for r in rows}


# ---------------------------------------------------------------------------------------------------- main
SESSION_COLS = ['session', 'animal', 'date', 'time', 'short', 'date_actual', 'name_date_mismatch', 'stim_folder',
                'status', 'stale', 'stale_reason', 'has_stim', 'has_tif', 'has_suite2p', 'has_fall', 'has_camera', 'has_dlc', 'has_aggregated', 'in_session_qc', 'agg_mtime', 'iscell_mtime',
                'task_type', 'stimuli', 'block_sequence', 'n_blocks', 'n_blocks_interrupted', 'spout', 'n_lick_entries', 'spout_note', 'task_level', 'level_mode', 'level_counts',
                'reward_mode', 'reward_delay_s', 'square_contrast',
                'n_tifs', 'minutes', 'n_planes', 'multi_plane', 'frame_rate', 'n_frames', 'zoom', 'fov_x_um', 'fov_y_um', 'nx', 'ny', 'um_per_px_x', 'um_per_px_y', 'power',
                'axes_x', 'axes_y', 'axes_z', 'fastz_pos', 'sample_dx', 'sample_dy', 'pose_placed', 'pose_y_um', 'pose_x_um', 'pose_confidence',
                'n_roi_total', 'n_roi_accepted', 'n_agg_rois', 'accepted_fraction', 'n_signal', 'n_noise', 'n_neuropil', 'n_master_axons', 'curation_files',
                'timing_problems', 'red_sync', 'n_cam_triggers', 'video_codec', 'pupil_frac_bad', 'whisk_fraction', 'whisk_artifact_fraction',
                'ball_readings', 'ball_triggers', 'ball_sync_ok', 'ball_note', 'qc_face_whisk_r', 'qc_face_pupil_r', 'qc_whisk_pupil_r', 'qc_face_sd_log', 'qc_pc1_signal',
                'auto_tags', 'meta_error',
                'pool', 'flag_reason', 'curation', 'tags', 'all_tags', 'quality', 'visible_activity', 'comments',
                'notion_planes', 'notion_quality', 'notion_visible_activity', 'notion_comments', 'notion_suite2p', 'notion_stimuli_camera',
                'fov_png', 'built']
BLOCK_COLS = ['session', 'block', 'file', 'type', 'key', 'stim_desc', 'stim_short', 'loom_dir', 'addnote', 'trials_requested', 'n_trials', 'start_wallclock', 'duration_s', 'interrupted',
              't_start_2p', 't_end_2p', 'red_sync_ok', 'red_sync_method', 'projector', 'setup',
              # square (training) blocks only, from the raw stim file; blank for the other stimuli
              'stim_shape', 'stim_contrast_pct', 'stim_size', 'stim_location', 'gui_level', 'level_mode', 'n_L1', 'n_L2', 'n_L3', 'trials_logged',
              'reward_mode', 'reward_delay_s', 'lick_window_s', 'wait_range_s', 'spout', 'n_lick_entries', 'n_trials_with_licks', 'lick_log_trials', 'spout_note',
              'params', 'note']
ROI_COLS = ['session', 'roi', 'accepted', 'iscell_prob', 'agg_row', 'class', 'y_px', 'x_px', 'y_um', 'x_um', 'npix', 'skew', 'master_axon', 'master_cls']


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--no-fov', action='store_true', help='skip the FOV pictures'); ap.add_argument('--no-class', action='store_true', help='skip the ROI classes')
    ap.add_argument('--snapshot', action='store_true', help='also copy the CSVs to axon analysis/recordings_db/snapshot (for git)'); a = ap.parse_args()
    os.makedirs(os.path.join(DB, 'fov'), exist_ok=True)
    built = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
    ms = {r['session']: r for r in read_csv(os.path.join(BUILD, 'm_sessions.csv'))}
    mb = read_csv(os.path.join(BUILD, 'm_blocks.csv'))
    if not ms: sys.exit('m_sessions.csv missing: run export_session_meta.m first (build_recdb.m does both stages)')
    names = sorted(n for n in os.listdir(DATA) if SESSION_RE.match(n) and os.path.isdir(os.path.join(DATA, n)))
    qc = {r['session']: r for r in read_csv(QC_PATH)}
    pose = {r['session']: r for r in read_csv(os.path.join(FOVM, 'fov_pose.csv'))}
    axons = read_csv(os.path.join(FOVM, 'master_axons.csv'))
    pres = read_csv(os.path.join(FOVM, 'master_presence.csv'))
    master_of = {(r['session'], r['roi_id']): (r['master'], r['cls']) for r in pres if r.get('roi_id') not in (None, '')}
    seed = seeds()
    notes, added = update_notes(names, seed, set(qc))
    ensure_animals(sorted({SESSION_RE.match(n).group(1) for n in names}))
    cache_p = os.path.join(BUILD, 'acq_cache.json'); cache = json.load(open(cache_p)) if os.path.isfile(cache_p) else {}

    S_rows, R_rows = [], []
    for nm in names:
        print(nm, end=' ', flush=True)
        m = SESSION_RE.match(nm); mrow = ms.get(nm, {})
        r = {c: '' for c in SESSION_COLS}
        r.update(session=nm, animal=m.group(1), date='20%s-%s-%s' % (m.group(2)[:2], m.group(2)[2:4], m.group(2)[4:6]), time='%s:%s' % (m.group(3)[:2], m.group(3)[2:]), short=short(nm), built=built)
        hd = read_header(nm, cache)
        for k in ('n_tifs', 'zoom', 'fov_x_um', 'fov_y_um', 'nx', 'ny', 'um_per_px_x', 'um_per_px_y', 'power', 'axes_x', 'axes_y', 'axes_z', 'fastz_pos', 'sample_dx', 'sample_dy', 'frame_rate'):
            r[k] = hd.get(k, '')
        if hd.get('tif_first_mtime'): r['date_actual'] = hd['tif_first_mtime'][:10]; r['name_date_mismatch'] = (r['date_actual'] != r['date'])
        r['acq_note'] = hd.get('acq_note', '')
        # -- from the MATLAB stage
        for k in ('stim_folder', 'task_type', 'stimuli', 'block_sequence', 'n_blocks', 'n_blocks_interrupted', 'spout', 'n_lick_entries', 'spout_note', 'task_level', 'level_mode', 'level_counts',
                  'reward_mode', 'reward_delay_s', 'square_contrast',
                  'timing_problems', 'n_cam_triggers', 'video_codec', 'pupil_frac_bad', 'whisk_fraction', 'whisk_artifact_fraction', 'ball_readings', 'ball_triggers', 'ball_sync_ok', 'agg_mtime', 'meta_error', 'n_agg_rois'):
            v = mrow.get(k, ''); r[k] = '' if str(v).lower() == 'nan' else v
        if mrow.get('n_blocks_red_sync_checked') not in (None, '', 'NaN') and num(mrow.get('n_blocks_red_sync_checked')):
            r['red_sync'] = '%d/%d blocks' % (num(mrow['n_blocks_red_sync_ok']), num(mrow['n_blocks_red_sync_checked']))
        # -- processing status
        has_agg = os.path.isfile(os.path.join(AGG, nm + '_preprocessed.mat'))
        s2 = load_s2p(nm); has_s2 = 'iscell' in s2
        cam = os.path.isdir(os.path.join(CAMERA, nm)) or os.path.isdir(os.path.join(PRECAMERA, m.group(1), nm)) or os.path.isdir(os.path.join(PRECAMERA, 'Animal' + m.group(1), nm))
        r.update(has_stim=(mrow.get('stim_folder_found') == '1' and (num(mrow.get('n_stim_files')) or 0) > 0), has_tif=bool(hd.get('n_tifs')), has_suite2p=has_s2,
                 has_fall=os.path.isfile(os.path.join(DATA, nm, 'suite2p', 'plane0', 'Fall.mat')), has_camera=cam,
                 has_dlc=os.path.isfile(os.path.join(DLC, nm + '_DLC_filtered.mat')), has_aggregated=has_agg, in_session_qc=nm in qc)
        r['status'] = 'analysed' if nm in qc else ('aggregated' if has_agg else ('suite2p done' if has_s2 else ('suite2p partial' if os.path.isdir(os.path.join(DATA, nm, 'suite2p')) else 'raw only')))
        # -- duration / planes
        # ops['nframes'] is NOT the session length (8413 in every session): the real number of frames is the width of F.npy (header read only)
        nfr = None
        try: nfr = int(np.load(os.path.join(s2p_dir(nm), 'F.npy'), mmap_mode='r').shape[1])
        except Exception: pass
        if mrow.get('duration_min_2p') not in (None, '', 'NaN'): r['minutes'] = mrow['duration_min_2p']
        elif nfr and hd.get('frame_rate') and (num(hd.get('n_slices')) or 1) == 1: r['minutes'] = round(nfr / hd['frame_rate'] / 60, 2)
        r['n_frames'] = mrow.get('n_frames_2p') if mrow.get('n_frames_2p') not in (None, '', 'NaN') else (nfr or '')
        npl = num(mrow.get('n_planes')); r['n_planes'] = int(npl) if npl else (hd.get('n_slices', '') if hd.get('n_slices') else '')
        r['multi_plane'] = 'yes' if (num(r['n_planes']) or 1) > 1 else ('no' if r['n_planes'] != '' else '')
        # -- ROIs
        classes = None
        if has_s2:
            ic = s2['iscell']; n_tot = len(ic); acc = ic[:, 0] > 0; n_acc = int(acc.sum())
            r.update(n_roi_total=n_tot, n_roi_accepted=n_acc, accepted_fraction=round(n_acc / max(n_tot, 1), 3), iscell_mtime=mtime_str(os.path.join(s2p_dir(nm), 'iscell.npy')), curation_files=curation_files(nm))
            stale = []
            nagg = num(r['n_agg_rois'])
            if has_agg and nagg is not None and int(nagg) != n_acc: stale.append('aggregated file has %d ROIs, iscell.npy accepts %d' % (nagg, n_acc))
            if has_agg and s2['iscell_mtime'] > os.path.getmtime(os.path.join(AGG, nm + '_preprocessed.mat')): stale.append('iscell.npy newer than aggregated file')
            r['stale'] = bool(stale) if has_agg else ''; r['stale_reason'] = '; '.join(stale)
            if has_agg and not stale and not a.no_class:
                try:
                    classes = roi_classes(nm)
                    if classes is not None and len(classes) != n_acc: classes = None; r['stale_reason'] = (r['stale_reason'] + '; xs export has a different ROI count').strip('; ')
                except Exception as e:
                    print('[class failed: %s]' % e, end=' ')
            if classes is not None: r.update(n_signal=int((classes == 'signal').sum()), n_noise=int((classes == 'noise').sum()), n_neuropil=int((classes == 'neuropil').sum()))
            # ROI table
            pxy, pxx = num(r['um_per_px_y']), num(r['um_per_px_x']); k = 0; masters = set()
            for i, st in enumerate(s2['stat']):
                med = st.get('med', [np.nan, np.nan]); a_ = bool(acc[i])
                mm = master_of.get((nm, str(i)), ('', ''))
                if mm[0] != '' and a_ and mm[1] == 'accepted': masters.add(mm[0])
                R_rows.append(dict(session=nm, roi=i, accepted=a_, iscell_prob=float(ic[i, 1]), agg_row=(k + 1 if a_ and has_agg and not stale else ''),
                                   **{'class': (classes[k] if (a_ and classes is not None) else '')}, y_px=round(float(med[0]), 1), x_px=round(float(med[1]), 1),
                                   y_um=(round(float(med[0]) * pxy, 2) if pxy else ''), x_um=(round(float(med[1]) * pxx, 2) if pxx else ''),
                                   npix=int(st.get('npix', 0)), skew=round(float(st.get('skew', np.nan)), 2), master_axon=mm[0], master_cls=mm[1]))
                if a_: k += 1
            r['n_master_axons'] = len(masters) if masters else ''
            if not a.no_fov and 'ops' in s2:
                try:
                    out = os.path.join(DB, 'fov', nm + '.png')
                    title = '%s   %s   %s plane(s)   %d accepted / %d ROIs   zoom %s' % (nm, r['task_type'] or '', r['n_planes'] or '?', n_acc, n_tot, r['zoom'] or '?')
                    fov_png(nm, s2, classes, out, pxx, pxy, title); r['fov_png'] = 'fov/' + nm + '.png'
                except Exception as e: print('[fov failed: %s]' % e, end=' ')
        # -- other tables
        q = qc.get(nm)
        if q:
            r.update(ball_note=q.get('ball', ''), qc_face_whisk_r=q.get('beh_face_whisk_r', ''), qc_face_pupil_r=q.get('beh_face_pupil_r', ''), qc_whisk_pupil_r=q.get('beh_whisk_pupil_r', ''),
                     qc_face_sd_log=q.get('beh_face_sd_log', ''), qc_pc1_signal=q.get('pc1_signal', ''))
        pz = pose.get(nm)
        if pz: r.update(pose_placed=pz.get('placed', ''), pose_y_um=pz.get('y_um', ''), pose_x_um=pz.get('x_um', ''), pose_confidence=pz.get('confidence', ''))
        # -- tags
        blocks = [b for b in mb if b['session'] == nm]; keys = {b['key'] for b in blocks}
        tg = sorted(keys - {'spont'}) + (['spont'] if 'spont' in keys else [])
        if any(k.startswith('loom') for k in keys): tg.append('loom')
        if r['task_type'] == 'lick task': tg.append('task')
        if r['task_type'] == 'passive square': tg.append('passive')
        if r.get('spout') in ('yes', 'no'): tg.append('spout' if r['spout'] == 'yes' else 'no_spout')       # square sessions only: spout = a lick was logged at all (blank for the other sessions)
        if r.get('level_mode') == 'adaptive': tg.append('adaptive_level')
        if r['multi_plane'] == 'yes': tg.append('multi_plane')
        if r['video_codec']: tg.append(str(r['video_codec']))
        if any('pre_grey_t' in b.get('params', '') or 'post_grey_t' in b.get('params', '') for b in blocks if b['key'].startswith('loom')): tg.append('loom_grey_param')
        if r.get('n_roi_total') and r.get('n_roi_accepted') == r.get('n_roi_total'): tg.append('all_rois_accepted')
        if r['stale'] is True: tg.append('stale')
        if r['has_aggregated'] and mrow.get('timing_problems') in ('1',): tg.append('timing_problem')
        if r['n_blocks_interrupted'] not in ('', '0', 0) and num(r['n_blocks_interrupted']): tg.append('interrupted')
        r['auto_tags'] = ';'.join(dict.fromkeys(tg))
        nt = notes.get(nm, {})
        for k in ('pool', 'flag_reason', 'curation', 'tags', 'quality', 'visible_activity', 'comments', 'notion_planes', 'notion_quality', 'notion_visible_activity', 'notion_comments', 'notion_suite2p', 'notion_stimuli_camera'):
            r[k] = nt.get(k, '')
        r['all_tags'] = ';'.join(dict.fromkeys([t for t in (r['auto_tags'] + ';' + r['tags']).split(';') if t]))
        S_rows.append(r)
    print()
    json.dump(cache, open(cache_p, 'w'), default=lambda o: o.item() if hasattr(o, 'item') else str(o))

    # blocks
    B_rows = []
    for b in mb:
        if b['session'] not in names: continue
        d = dict(b); d['start_wallclock'] = b.get('start_wallclock', '') or b.get('start_wallclock_1', '')
        for k, v in list(d.items()):
            if isinstance(v, str) and v.lower() == 'nan': d[k] = ''
        B_rows.append(d)
    A_cols = ['master', 'n_sessions_accepted', 'n_sessions_in_view', 'y_um', 'x_um', 'n_cells', 'area_um2', 'sessions']
    write_csv(os.path.join(DB, 'sessions.csv'), S_rows, SESSION_COLS)
    write_csv(os.path.join(DB, 'blocks.csv'), B_rows, BLOCK_COLS)
    write_csv(os.path.join(DB, 'rois.csv'), R_rows, ROI_COLS)
    if axons: write_csv(os.path.join(DB, 'axons.csv'), axons, A_cols)
    print('wrote sessions %d, blocks %d, rois %d, axons %d to %s' % (len(S_rows), len(B_rows), len(R_rows), len(axons), DB))
    if added: print('notes.csv: %d new session rows seeded' % len(added))
    make_readme(); make_index(S_rows)
    if a.snapshot:
        os.makedirs(SNAP, exist_ok=True)
        for f in ('sessions.csv', 'blocks.csv', 'rois.csv', 'axons.csv', 'animals.csv', 'notes.csv'):
            if os.path.isfile(os.path.join(DB, f)): shutil.copy2(os.path.join(DB, f), os.path.join(SNAP, f))
        print('snapshot copied to', SNAP)


# ---------------------------------------------------------------------------------------------------- browse page + README
def make_index(rows):
    import html
    cols = [('session', 'session'), ('status', 'status'), ('task_type', 'type'), ('stimuli', 'stimuli'), ('spout', 'spout'), ('minutes', 'min'), ('n_planes', 'planes'), ('n_roi_accepted', 'acc ROIs'),
            ('n_signal', 'signal-like'), ('pool', 'pool'), ('stale', 'stale'), ('flag_reason', 'flag / reason'), ('all_tags', 'tags')]
    h = ['<!doctype html><html><head><meta charset="utf-8"><title>RB22-24 recordings</title><style>',
         'body{font:13px system-ui,sans-serif;margin:16px;background:#fff;color:#222}table{border-collapse:collapse}th,td{border:1px solid #ccc;padding:3px 6px;vertical-align:top}',
         'th{background:#f0f0f0;position:sticky;top:0;cursor:pointer}tr:hover{background:#f7fbff}img{height:110px;display:block}.no{color:#b00}.stale{background:#fff3cd}',
         '#q{padding:5px;width:380px;margin:0 0 8px}small{color:#666}</style></head><body>',
         '<h2>RB22-24 recordings (%d sessions)</h2><input id="q" placeholder="filter: e.g. RB24 loom_expand task 5plane"> <small>click a column header to sort; text filter = all words must appear in the row</small>' % len(rows),
         '<table id="t"><thead><tr><th>FOV</th>' + ''.join('<th>%s</th>' % c[1] for c in cols) + '</tr></thead><tbody>']
    for r in rows:
        img = '<img src="%s">' % r['fov_png'] if r.get('fov_png') else ''
        cells = ''.join('<td class="%s">%s</td>' % ('no' if (c[0] == 'pool' and r.get(c[0]) == 'no') else ('stale' if (c[0] == 'stale' and r.get(c[0]) is True) else ''),
                                                    html.escape(fmt(r.get(c[0])))) for c in cols)
        h.append('<tr><td>%s</td>%s</tr>' % (img, cells))
    h.append('</tbody></table><script>const q=document.getElementById("q"),rows=[...document.querySelectorAll("#t tbody tr")];'
             'q.oninput=()=>{const w=q.value.toLowerCase().split(/\\s+/).filter(Boolean);rows.forEach(r=>{const t=r.innerText.toLowerCase();r.style.display=w.every(x=>t.includes(x))?"":"none"})};'
             'document.querySelectorAll("#t th").forEach((th,i)=>th.onclick=()=>{const b=document.querySelector("#t tbody");const a=[...b.rows];const d=th.dataset.d=th.dataset.d==="1"?"-1":"1";'
             'a.sort((x,y)=>{const u=x.cells[i].innerText,v=y.cells[i].innerText,nu=parseFloat(u),nv=parseFloat(v);return (isNaN(nu)||isNaN(nv)?u.localeCompare(v):nu-nv)*d});a.forEach(r=>b.appendChild(r))});'
             '</script></body></html>')
    open(os.path.join(DB, 'index.html'), 'w', encoding='utf-8').write('\n'.join(h))


def make_readme():
    src = os.path.join(AX, 'recordings_db', 'README_db.md')
    if os.path.isfile(src): shutil.copy2(src, os.path.join(DB, 'README.md'))


if __name__ == '__main__':
    main()
