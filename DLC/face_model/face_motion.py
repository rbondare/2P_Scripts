"""Whisking on/off and amplitude, and face motion, from the raw camera video, in regions anchored on the nose tip.

    python face_motion.py <raw video.mp4> <face csv from analyze_face.py> --out-dir <folder> [--name <session name>]
    python face_motion.py --reclassify <raw video.mp4> <whisker csv> --out-dir <folder>    (new REST / rules, no decoding)

The face DLC model (analyze_face.py) tracks the nose tip. The head sits in a different place in every session (the
median nose position differs by up to ~120 px between sessions), so every region is placed relative to the nose. If the
face csv covers the whole video the anchor is the rolling median of the confident nose positions (ANCHOR_WINDOW_S);
if it covers only one clip the anchor is the median over that clip; if it holds several sampled clips
(analyze_face.py --sample) the anchor follows the clip medians over the session. Regions (dx0, dy0, dx1, dy1 from the
anchor, full-resolution px):
  front_top  whisker shafts in front of and above the snout, ABOVE the lick spout: the whisking index (WHISKER_BOX)
  fan_front  the whole fan in front of the snout; its lower part holds the spout tip in lick sessions (not used for
             the index, kept in the csv for comparison)
  fan_down   chin, lower jaw, tongue and the shafts hanging below the pad: JAW / LICK motion, not whisking
  pad        whisker pad skin right of the nose: face motion (moves with whisking, but also with chewing, breathing)
  ctrl       background without face or body, a noise reference
Why front_top only (RB24_261109_1506, 4,734 licks with the independent lick line): licking explains 45% of the variance of
the former index (mean of fan_front and fan_down) and 93% of its running-vs-still contrast, 87% of that through the
jaw box; front_top has 0.2% lick-explained variance and a lick effect of 4% of the running contrast (CI 0-9%) while keeping
the largest running-vs-still effect size (d 1.16 vs 1.02). See axon analysis/whisking_lick_check/.
Motion energy (ME) = mean |frame t - frame t-1| inside a region (both frames cropped at the same absolute position),
after a 1 px Gaussian blur. The raw video has a keyframe every 12 frames; the keyframe and the frame after it carry codec
noise, not motion (up to 4x background), so these two frames of every 12 are interpolated.

whisk_index = smoothed ME of front_top / REST, where REST is the resting ME of that region (median over the sessions
that have a clear rest period) for the codec of the video, so 1 = whiskers at rest. The camera files were MPEG-4
(~24 Mb/s) up to 2026-10-29 and are H.264 (~50 Mb/s) since 2026-10-30; the H.264 files keep much more sensor noise, so
the resting ME is about twice as high and each codec has its own REST (see calibrate_rest.py). A per-session floor
cannot be used: some mice whisk through the whole session and have no rest period.
whisk_on = hysteresis threshold on the index, ON_OFF per codec (short gaps closed, short bouts dropped): MPEG-4 on above
1.76 / off below 1.44 (1.6 x REST +-10 %), H.264 on above 1.45 / off below 1.40. whisk_amp = index - 1 (0 = rest).
Between 'held back' and 'clearly whisking' there is a continuum, so the threshold is a convention: whisk_amp carries the
graded information. ME follows the speed of the whiskers, not their angle, and cannot tell whisker motion from head motion.
The H.264 levels come from blind human labels (51 clips, 30 whisking / 21 still, RB24 H.264 sessions): the MPEG-4 rule
called 6 of them wrong (all whisking called rest, 0 false), 1.45 / 1.40 calls 3 wrong, and it makes the binary state
couple to the LC-axon dF/F as well as the continuous index (6 sessions). The H.264 noise floor needs about twice as much
whisker motion for the same index, so 1.76 was too insensitive there. Not changed for MPEG-4 (22 of 23 blind clips right).
See axon analysis/whisking_lick_check/fix2_rollout/h264_blind_check/.
Frames with an implausibly high index (a forepaw in front of the face when grooming) are flagged as artifact and
left empty in whisk_amp and whisk_on. The flag is raised when either the whisking index or the whole-face motion (the former
two-box index of fan_front and fan_down) exceeds ARTIFACT_X times its own 95th percentile, because grooming mostly shows in
the jaw / pad boxes and hardly in front_top.
"""
import argparse
import csv
import json
import os
import shutil
import subprocess
import sys

import cv2
import numpy as np
import pandas as pd
from scipy.ndimage import binary_dilation, gaussian_filter1d

FPS = 50.0
ANCHOR_WINDOW_S = 10.0
NOSE_MIN_LIKELIHOOD = 0.9
MIN_SAMPLE_FRAMES = 25  # confident nose frames a sampled clip needs to count
BLUR_SIGMA_PX = 1.0
BOXES = {
    'fan_front': (-260, -110, 60, 200),
    'front_top': (-260, -110, 60, 0),
    'fan_down': (60, 80, 300, 300),
    'pad': (60, -100, 260, 60),
    'ctrl': (-320, -300, -200, -180),
}
WHISKER_BOX = 'front_top'
# resting ME (gray levels/frame) of the whisker region (front_top) per video codec: median over the sessions of that codec
# that have a clear rest period (10th percentile within 15% of the 3rd) of the per-session 3rd percentile of the smoothed
# ME, calibrate_rest.py on the 43 RB22-24 sessions of 2026-10-17..11-09 (14 of 29 MPEG-4 sessions, range 0.30-0.58, CV 18%;
# 9 of 14 H.264 sessions, range 0.84-1.65, the tight cluster of 10-30..11-05 RB24 sessions is 0.84-0.98). MPEG-4 is the
# camera file format up to 2026-10-29, H.264 (50 Mb/s, ~2x the sensor noise) from 2026-10-30.
REST = {'mpeg4': 0.473, 'h264': 0.952}
SMOOTH_SIGMA_FRAMES = 3  # 60 ms
# whisk_on hysteresis on the index (1 = rest) per codec: (off below, on above). MPEG-4: 1.6 x REST +-10 % (unchanged since the
# first version); H.264: 1.40 / 1.45, set from the blind labels (docstring).
ON_OFF = {'mpeg4': (1.6 * (1 - 0.1), 1.6 * (1 + 0.1)), 'h264': (1.40, 1.45)}
MIN_GAP_S = 0.3  # off periods shorter than this inside a bout are closed
MIN_BOUT_S = 0.3  # on periods shorter than this are dropped
ARTIFACT_X = 2.0  # smoothed index above this multiple of its 95th percentile is not whisking
# A forepaw / grooming shows up mostly in the jaw and pad boxes and hardly at all in front_top (RB24_261105_1549: 3 grooming
# episodes, 8-12x median in fan_down, front_top <= 1.5x its 95th percentile), so the artifact flag also uses the former
# two-box motion index (mean of fan_front and fan_down over these fixed scales; only its ratio to its own 95th percentile
# matters). This keeps the flags of the earlier version, which were checked on frames.
ARTIFACT_BOXES = {'fan_front': 0.447, 'fan_down': 0.450}
ARTIFACT_PAD_S = 0.5  # flagged periods are widened by this much on both sides


def find_tool(name):
    exe = shutil.which(name)
    if exe:
        return exe
    for p in (os.path.join(os.path.dirname(sys.executable), 'Library', 'bin', name + '.exe'),
              os.path.join(os.path.dirname(sys.executable), 'Library', 'bin', name)):
        if os.path.isfile(p):
            return p
    raise FileNotFoundError(name + ' not found on PATH (activate the dlc-tf environment)')


def video_codec(video):
    """Codec name of the first video stream ('mpeg4', 'h264', ...)."""
    out = subprocess.run([find_tool('ffprobe'), '-v', 'error', '-select_streams', 'v:0', '-show_entries',
                          'stream=codec_name', '-of', 'csv=p=0', video], capture_output=True, text=True, check=True)
    return out.stdout.strip().split(',')[0]


def rest_for(codec):
    if codec not in REST:
        raise RuntimeError("no resting level for codec '%s': run calibrate_rest.py on sessions of this codec and add it to REST"
                           % codec)
    return REST[codec]


def on_off_for(codec):
    """(off below, on above) of whisk_on on the index for this codec"""
    if codec not in ON_OFF:
        raise RuntimeError("no on/off levels for codec '%s': add it to ON_OFF (and REST)" % codec)
    return ON_OFF[codec]


def keyframes(video):
    """Boolean array over frames, True where the video has a keyframe (read from the packet flags)."""
    out = subprocess.run([find_tool('ffprobe'), '-v', 'error', '-select_streams', 'v:0', '-show_entries',
                          'packet=flags', '-of', 'csv=p=0', video], capture_output=True, text=True, check=True)
    return np.array(['K' in line for line in out.stdout.split()], bool)


def drop_codec_frames(me, is_key):
    """Interpolate over the keyframe and the frame after it (codec noise); intra-only video is left alone."""
    n = len(me)
    key = is_key[:n]
    if key.mean() > 0.5:  # every frame is a keyframe: no inter-frame artefact to remove
        bad = np.zeros(n, bool)
    else:
        bad = key.copy()
        bad[1:] |= key[:-1]
    bad[0] = True  # frame 0 has no predecessor
    good = np.flatnonzero(~bad)
    return np.interp(np.arange(n), good, me[good])


def hysteresis(x, lo, hi):
    on = np.zeros(len(x), bool)
    state = False
    for i, v in enumerate(x):
        if state and v < lo:
            state = False
        elif not state and v > hi:
            state = True
        on[i] = state
    return on


def runs(mask):
    """(start, stop) of every run of True in mask (stop exclusive)."""
    d = np.diff(np.r_[False, mask, False].astype(int))
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)))


def clean_bouts(on, min_gap, min_len):
    on = on.copy()
    for a, b in runs(~on):
        if a > 0 and b < len(on) and b - a < min_gap:  # gap between two bouts
            on[a:b] = True
    for a, b in runs(on):
        if b - a < min_len:
            on[a:b] = False
    return on


def nose_anchor(face, n_frames):
    """Per-frame integer (x, y) anchor from the nose track. A track that covers the video gives a rolling median of
    the confident positions; a track made of several sampled clips ('sample' and 'source_frame' columns, see
    analyze_face.py --sample) gives each clip's median at the clip's centre, interpolated over the video (constant
    before the first and after the last clip); a single shorter clip gives one median for every frame."""
    ok = face['noseTip_likelihood'].values >= NOSE_MIN_LIKELIHOOD
    if ok.sum() < 100:
        raise RuntimeError('the nose is confident in only %d frames; cannot place the regions' % ok.sum())
    if 'sample' in face and face['sample'].nunique() > 1:
        confident = face[ok]
        groups = [g for _, g in confident.groupby('sample') if len(g) >= MIN_SAMPLE_FRAMES]
        if not groups:
            raise RuntimeError('no sampled clip has %d confident nose frames' % MIN_SAMPLE_FRAMES)
        centre = [g['source_frame'].mean() for g in groups]
        frames = np.arange(n_frames)
        return [np.round(np.interp(frames, centre, [g[c].median() for g in groups])).astype(int)
                for c in ('noseTip_x', 'noseTip_y')]
    if len(face) < 0.95 * n_frames:
        return [np.full(n_frames, int(round(float(np.median(face[c].values[ok]))))) for c in ('noseTip_x', 'noseTip_y')]
    out = []
    for c in ('noseTip_x', 'noseTip_y'):
        s = pd.Series(np.where(ok, face[c].values, np.nan))
        r = s.rolling(int(ANCHOR_WINDOW_S * FPS), center=True, min_periods=25).median()
        r = np.round(r.bfill().ffill().values).astype(int)
        out.append(np.resize(r, n_frames))
    return out


def video_frames(video):
    cap = cv2.VideoCapture(video)
    n, w, h = int(cap.get(7)), int(cap.get(3)), int(cap.get(4))
    cap.release()
    return n, w, h


def anchored_motion_energy(video, anchor, boxes=BOXES):
    """ME per frame for every box; frame 0 gives NaN. anchor = (ax, ay) integer arrays, one value per frame."""
    _, width, height = video_frames(video)
    ax, ay = anchor
    # one decoded window that holds every box at every anchor position
    x0 = max(int(ax.min() + min(b[0] for b in boxes.values())), 0)
    x1 = min(int(ax.max() + max(b[2] for b in boxes.values())), width)
    y0 = max(int(ay.min() + min(b[1] for b in boxes.values())), 0)
    y1 = min(int(ay.max() + max(b[3] for b in boxes.values())), height)
    w, h = x1 - x0, y1 - y0
    # exact=1: without it crop rounds the size down to even numbers on 4:2:0 video, the pipe then holds frames narrower
    # than w and every frame is read with a growing shear (huge fake motion) when the window width is odd
    cmd = [find_tool('ffmpeg'), '-v', 'error', '-i', video, '-vf', 'crop=%d:%d:%d:%d:exact=1,format=gray' % (w, h, x0, y0),
           '-f', 'rawvideo', 'pipe:1']
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=w * h * 32)
    me = {k: [] for k in boxes}
    prev = None
    i = 0
    while True:
        buf = proc.stdout.read(w * h)
        if len(buf) < w * h:
            if buf:
                raise RuntimeError('ffmpeg ended with a partial frame (%d of %d bytes): the crop size is not what was asked'
                                   % (len(buf), w * h))
            break
        g = cv2.GaussianBlur(np.frombuffer(buf, np.uint8).reshape(h, w), (0, 0), BLUR_SIGMA_PX).astype(np.float32)
        a = min(i, len(ax) - 1)
        for k, (bx0, by0, bx1, by1) in boxes.items():
            c0, c1 = max(ax[a] + bx0 - x0, 0), min(ax[a] + bx1 - x0, w)
            r0, r1 = max(ay[a] + by0 - y0, 0), min(ay[a] + by1 - y0, h)
            if prev is None or c1 - c0 < 20 or r1 - r0 < 20:
                me[k].append(np.nan)
            else:
                me[k].append(float(np.abs(g[r0:r1, c0:c1] - prev[r0:r1, c0:c1]).mean()))
        prev = g
        i += 1
    proc.wait()
    return {k: np.array(v) for k, v in me.items()}


def clean_me(me, video):
    """Interpolate over the keyframe + next frame (codec noise) in every trace."""
    n = len(next(iter(me.values())))
    is_key = keyframes(video)
    if len(is_key) != n:
        print('WARNING: %d decoded frames but %d packets; keyframe positions may be off' % (n, len(is_key)))
        is_key = np.resize(is_key, n)
    return {k: drop_codec_frames(v, is_key) for k, v in me.items()}


def classify(me, rest, codec):
    """From the cleaned ME traces to the whisking index, artifact flags and on/off state. rest = resting ME of WHISKER_BOX,
    codec ('mpeg4' / 'h264') selects the on/off levels (ON_OFF)."""
    sm = {k: gaussian_filter1d(v, SMOOTH_SIGMA_FRAMES) for k, v in me.items()}
    index = sm[WHISKER_BOX] / rest
    motion = np.mean([sm[k] / r for k, r in ARTIFACT_BOXES.items()], axis=0)  # whole-face motion (the former index)
    odd = (index > ARTIFACT_X * np.percentile(index, 95)) | (motion > ARTIFACT_X * np.percentile(motion, 95))
    artifact = binary_dilation(odd, iterations=int(ARTIFACT_PAD_S * FPS))
    off_below, on_above = on_off_for(codec)
    on = hysteresis(index, off_below, on_above)
    on = clean_bouts(on, int(MIN_GAP_S * FPS), int(MIN_BOUT_S * FPS))
    on[artifact] = False
    return dict(smooth=sm, index=index, artifact=artifact, on=on, amp=np.maximum(index - 1, 0), on_off=(off_below, on_above))


def process(video, face_csv, name=None):
    name = name or os.path.splitext(os.path.basename(video))[0]
    n_frames, width, height = video_frames(video)
    face = pd.read_csv(face_csv, index_col=0)
    anchor = nose_anchor(face, n_frames)
    codec = video_codec(video)
    rest = rest_for(codec)
    me = anchored_motion_energy(video, anchor)
    me = clean_me(me, video)
    res = classify(me, rest, codec)
    res.update(name=name, video=video, size=(width, height), anchor=anchor, me=me, face_frames=len(face), codec=codec, rest=rest)
    return res


def reclassify(video, whisker_csv, name=None):
    """Same result as process() but from the ME traces of an earlier run (the me_* columns of its csv and the anchor
    columns): a new REST or new rules without decoding the video again. The anchor is that run's per-frame anchor."""
    T = pd.read_csv(whisker_csv)
    me = {k: T['me_' + k].values.astype(float) for k in BOXES}
    anchor = [T['anchor_x'].values.astype(int), T['anchor_y'].values.astype(int)]
    info = json.load(open(whisker_csv.replace('_whisker.csv', '_whisker.json')))
    _, width, height = video_frames(video)
    codec = video_codec(video)
    rest = rest_for(codec)
    res = classify(me, rest, codec)
    res.update(name=name or info['session'], video=video, size=(width, height), anchor=anchor, me=me,
               face_frames=info['face_track_frames'], codec=codec, rest=rest)
    return res


def summary(res):
    on, art, n = res['on'], res['artifact'], len(res['on'])
    bouts = runs(on)
    idx = res['index'][~art]
    ax, ay = res['anchor']
    return dict(
        session=res['name'], frames=n, fps=FPS, face_track_frames=res['face_frames'],
        codec=res['codec'], whisker_box=WHISKER_BOX,
        anchor_median_xy_px=[int(np.median(ax)), int(np.median(ay))],
        anchor_range_px=[int(ax.max() - ax.min()), int(ay.max() - ay.min())],
        boxes_rel_to_nose_px={k: list(v) for k, v in BOXES.items()},
        fraction_whisking=round(float(on[~art].mean()), 3), fraction_artifact=round(float(art.mean()), 4),
        n_bouts=len(bouts), median_bout_s=round(float(np.median([b - a for a, b in bouts])) / FPS, 2) if bouts else None,
        index_p5=round(float(np.percentile(idx, 5)), 2), index_median=round(float(np.median(idx)), 2),
        index_p95=round(float(np.percentile(idx, 95)), 2),
        rest_level={WHISKER_BOX: res['rest']},
        whisk_off_below_index=round(float(res['on_off'][0]), 4), whisk_on_above_index=round(float(res['on_off'][1]), 4),
        ctrl_me_median=round(float(np.median(res['me']['ctrl'])), 3),
        params=dict(blur_sigma_px=BLUR_SIGMA_PX, smooth_sigma_frames=SMOOTH_SIGMA_FRAMES,
                    min_gap_s=MIN_GAP_S, min_bout_s=MIN_BOUT_S, artifact_x_p95=ARTIFACT_X, artifact_pad_s=ARTIFACT_PAD_S,
                    artifact_motion_boxes=dict(ARTIFACT_BOXES)))


def write_outputs(res, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    name, n = res['name'], len(res['on'])
    art = res['artifact']
    me = res['me']
    with open(os.path.join(out_dir, name + '_whisker.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['frame', 'time_s', 'whisk_index', 'whisk_amp', 'whisk_on', 'artifact'] + ['me_' + k for k in BOXES] +
                   ['anchor_x', 'anchor_y'])
        ax, ay = res['anchor']
        for i in range(n):
            a = min(i, len(ax) - 1)
            w.writerow([i, '%.2f' % (i / FPS), '%.3f' % res['index'][i],
                        '' if art[i] else '%.3f' % res['amp'][i], '' if art[i] else int(res['on'][i]), int(art[i])] +
                       ['%.4f' % me[k][i] for k in BOXES] + [int(ax[a]), int(ay[a])])
    info = summary(res)
    with open(os.path.join(out_dir, name + '_whisker.json'), 'w') as f:
        json.dump(info, f, indent=2)
    qc_figure(res, os.path.join(out_dir, name + '_whisker_qc.png'))
    return info


def qc_figure(res, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    idx, on, art, n = res['index'], res['on'], res['artifact'], len(res['index'])
    t = np.arange(n) / FPS
    fig = plt.figure(figsize=(16, 9.5))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 1.25])

    def shade(ax, lo, hi):
        for a, b in runs(on):
            if b / FPS > lo and a / FPS < hi:
                ax.axvspan(a / FPS, b / FPS, color='tab:green', alpha=.15, lw=0)
        for a, b in runs(art):
            if b / FPS > lo and a / FPS < hi:
                ax.axvspan(a / FPS, b / FPS, color='tab:red', alpha=.25, lw=0)
    ax = fig.add_subplot(gs[0, :])
    step = max(n // 6000, 1)
    ax.plot(t[::step], idx[::step], lw=.4)
    for lev in res['on_off']:
        ax.axhline(lev, c='tab:orange', ls='--', lw=1, label='on / off levels' if lev == res['on_off'][0] else None)
    ax.axhline(1, c='gray', ls=':', lw=1, label='rest (=1)')
    shade(ax, 0, t[-1])
    ax.set_ylim(0, np.percentile(idx, 99.8) * 1.1)
    ax.set_xlim(0, t[-1])
    ax.set_ylabel('whisking index (1 = rest)')
    ax.set_title('%s: green = whisking, red = artifact (grooming etc.)' % res['name'])
    ax.legend(loc='upper right', ncol=2)
    ax = fig.add_subplot(gs[1, :])
    lo = float(np.clip(t[-1] / 2 - 30, 0, max(t[-1] - 60, 0)))
    seg = slice(int(lo * FPS), int((lo + 60) * FPS))
    ax.plot(t[seg], idx[seg], lw=.8, c='tab:blue', label='whisking index')
    for k, col, lab in (('fan_down', 'tab:red', 'jaw / lick motion'), ('pad', 'tab:orange', 'pad (face) motion')):
        sm = gaussian_filter1d(res['me'][k], SMOOTH_SIGMA_FRAMES)
        ax.plot(t[seg], sm[seg] / np.percentile(sm, 5), lw=.6, c=col, alpha=.6, label=lab + ' (x its 5th pct)')
    for lev in res['on_off']:
        ax.axhline(lev, c='tab:orange', ls='--', lw=1)
    shade(ax, lo, lo + 60)
    ax.set_xlim(lo, lo + 60)
    ax.set_ylim(0, np.percentile(idx[seg], 99.5) * 1.1)
    ax.set_xlabel('time (s)')
    ax.set_title('60 s zoom')
    ax.legend(loc='upper right', ncol=2)
    ax = fig.add_subplot(gs[2, 0])
    ax.hist(np.log10(np.clip(idx[~art], 0.05, None)), 150, color='tab:blue')
    for lev in res['on_off']:
        ax.axvline(np.log10(lev), c='tab:orange', ls='--')
    ax.axvline(0, c='gray', ls=':')
    ax.set_xlabel('log10 whisking index')
    ax.set_title('distribution (rest dotted, threshold dashed)')
    ax = fig.add_subplot(gs[2, 1])
    cap = cv2.VideoCapture(res['video'])
    cap.set(cv2.CAP_PROP_POS_FRAMES, n // 2)
    ok, frame = cap.read()
    cap.release()
    if ok:
        ax.imshow(cv2.cvtColor(frame, cv2.COLOR_BGR2RGB))
        ax0, ay0 = res['anchor'][0][n // 2], res['anchor'][1][n // 2]
        for k, col in zip(BOXES, ('yellow', 'lime', 'magenta', 'orange', 'cyan')):
            x0, y0, x1, y1 = BOXES[k]
            ax.add_patch(plt.Rectangle((ax0 + x0, ay0 + y0), x1 - x0, y1 - y0, fill=False, ec=col, lw=3 if k == WHISKER_BOX else 1))
        ax.plot(ax0, ay0, 'r.', ms=8)
        ax.set_xlim(max(ax0 - 450, 0), min(ax0 + 450, res['size'][0]))
        ax.set_ylim(min(ay0 + 400, res['size'][1]), max(ay0 - 350, 0))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('regions on the middle frame (thick green = whisking index, red = nose anchor)', fontsize=9)
    ax = fig.add_subplot(gs[2, 2])
    ax.axis('off')
    info = summary(res)
    ax.text(0, 1, '\n'.join([
        'frames: %d (%.1f min)' % (info['frames'], info['frames'] / FPS / 60),
        'whisking: %.0f%% of valid frames' % (100 * info['fraction_whisking']),
        'bouts: %d, median %.1f s' % (info['n_bouts'], info['median_bout_s'] or 0),
        'index 5th / median / 95th pct:', '  %.1f / %.1f / %.1f' % (info['index_p5'], info['index_median'], info['index_p95']),
        'artifact (grooming etc.): %.1f%%' % (100 * info['fraction_artifact']),
        'nose anchor: (%d, %d), moved %s px' % (*info['anchor_median_xy_px'], info['anchor_range_px']),
        'background ME median: %.2f' % info['ctrl_me_median'],
        'codec %s, rest (index 1) = %.3f ME' % (info['codec'], res['rest'])]),
        va='top', family='monospace', fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('video')
    ap.add_argument('face_csv', help='face csv from analyze_face.py (with --reclassify: the whisker csv of an earlier run)')
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--name', help='name for the output files (default: the video file name)')
    ap.add_argument('--reclassify', action='store_true',
                    help='redo REST / index / artifact / on-off from the ME columns of an earlier whisker csv (no decoding)')
    args = ap.parse_args()
    res = reclassify(args.video, args.face_csv, args.name) if args.reclassify else process(args.video, args.face_csv, args.name)
    info = write_outputs(res, args.out_dir)
    print('whisking %.0f%%, %d bouts (median %.1f s), artifact %.1f%%, index 5/50/95 pct %.1f/%.1f/%.1f' % (
        100 * info['fraction_whisking'], info['n_bouts'], info['median_bout_s'] or 0, 100 * info['fraction_artifact'],
        info['index_p5'], info['index_median'], info['index_p95']))


if __name__ == '__main__':
    main()
