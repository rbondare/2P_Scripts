"""Checks after h264_threshold_rollout.py apply / sync-stale (all read-only):
 1. every copy of the whisker results (work, quality_reports, preCamera / Camera) is byte-identical to DLC\\work for all 43 sessions
 2. whisk_on of every H.264 session = an independent re-implementation of the rule (own hysteresis / gap / bout code, own gaussian
    smoothing of the saved ME trace, artifact mask from the csv) - exact equality; MPEG-4 sessions with the old rule 1.76 / 1.44 too
 3. json fields: codec, rest_level, whisk_off_below_index / whisk_on_above_index (H.264 1.40 / 1.45), whisker_box
 4. backups exist (whisker_results_v2 for the H.264 sessions, whisker_results_v1 for every stale Camera copy)
 5. the 51 human-labelled clips scored with the SAVED whisk_on (clip called whisking if on for > 50 % of its frames)
"""
import os, sys, json, glob, hashlib, csv
import numpy as np, pandas as pd
from scipy.ndimage import gaussian_filter1d

WORK = r'C:\Users\rbondarenko\DLC\work'
PRE = r'Z:\Group Members\Rima\preCamera'
CAM = r'Z:\Group Members\Rima\Camera'
FPS = 50.0
REST = {'mpeg4': 0.473, 'h264': 0.952}
LEVELS = {'mpeg4': (1.44, 1.76), 'h264': (1.40, 1.45)}
codec_of = {r[0]: r[1] for r in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'codec_table.json')))}
bad = []


def md5(p):
    return hashlib.md5(open(p, 'rb').read()).hexdigest() if os.path.isfile(p) else None


def my_rule(me_top, art, off_below, on_above, rest):
    """independent version of face_motion.classify's on/off: smooth, hysteresis, close gaps < 0.3 s, drop bouts < 0.3 s, clear artifact"""
    x = gaussian_filter1d(me_top, 3) / rest
    on = np.zeros(len(x), bool); st = False
    for i in range(len(x)):
        if st and x[i] < off_below: st = False
        elif (not st) and x[i] > on_above: st = True
        on[i] = st
    def runs(m):
        edges = np.flatnonzero(np.diff(np.r_[0, m.astype(np.int8), 0]))
        return list(zip(edges[::2], edges[1::2]))
    n = len(on); gap = int(0.3 * FPS)
    for a, b in runs(~on):
        if a > 0 and b < n and b - a < gap: on[a:b] = True
    for a, b in runs(on):
        if b - a < gap: on[a:b] = False
    on[art] = False
    return on


sessions = sorted(d for d in os.listdir(WORK) if d.startswith('Animal') and os.path.isfile(os.path.join(WORK, d, 'whisker_results', d + '_whisker.csv')))
print('%d sessions' % len(sessions))
print('\n1-4. copies, independent rule, json, backups')
n_copy = 0
for s in sessions:
    a = s.split('_')[0]
    wr = os.path.join(WORK, s, 'whisker_results')
    cp = [p for p in (os.path.join(PRE, a, s, 'whisker_results'), os.path.join(CAM, s, 'whisker_results')) if os.path.isdir(p)]
    for c in cp:
        for suf in ('_whisker.csv', '_whisker.json', '_whisker_qc.png'):
            n_copy += 1
            if md5(os.path.join(c, s + suf)) != md5(os.path.join(wr, s + suf)):
                bad.append('%s: %s differs from work' % (s, os.path.join(c, s + suf)))
    qd = os.path.join(WORK, 'quality_reports', s)
    if os.path.isdir(qd):
        for suf in ('_whisker.json', '_whisker_qc.png'):
            n_copy += 1
            if md5(os.path.join(qd, s + suf)) != md5(os.path.join(wr, s + suf)):
                bad.append('%s: quality_reports %s differs' % (s, suf))
    cod = codec_of[s]
    T = pd.read_csv(os.path.join(wr, s + '_whisker.csv'))
    J = json.load(open(os.path.join(wr, s + '_whisker.json')))
    art = T['artifact'].values.astype(bool)
    on_saved = T['whisk_on'].fillna(0).values.astype(bool)
    on_mine = my_rule(T['me_front_top'].values.astype(float), art, *LEVELS[cod], REST[cod])
    ndiff = int((on_saved != on_mine).sum())
    if ndiff: bad.append('%s: saved whisk_on differs from the independent rule in %d frames' % (s, ndiff))
    if not (T['whisk_on'].isna().values == art).all(): bad.append('%s: whisk_on NaN pattern is not the artifact mask' % s)
    if J['codec'] != cod: bad.append('%s: json codec %s, expected %s' % (s, J['codec'], cod))
    if J.get('whisker_box') != 'front_top': bad.append('%s: json whisker_box' % s)
    if abs(J['rest_level']['front_top'] - REST[cod]) > 1e-9: bad.append('%s: json rest' % s)
    if cod == 'h264':
        if (J.get('whisk_off_below_index'), J.get('whisk_on_above_index')) != (1.4, 1.45): bad.append('%s: json levels %s' % (s, (J.get('whisk_off_below_index'), J.get('whisk_on_above_index'))))
        for base in [wr] + cp:
            if not os.path.isdir(os.path.join(os.path.dirname(base), 'whisker_results_v2')) and not os.path.isdir(os.path.join(os.path.dirname(base), 'whisker_results_v1')):
                bad.append('%s: no backup folder next to %s' % (s, base))
        if not os.path.isdir(os.path.join(WORK, s, 'whisker_results_v2')): bad.append('%s: work backup v2 missing' % s)
    else:
        for c in cp:
            if os.path.basename(os.path.dirname(c)) == s and 'Camera' in c and not os.path.isdir(os.path.join(os.path.dirname(c), 'whisker_results_v1')) and md5(os.path.join(c, s + '_whisker.csv')) != md5(os.path.join(wr, s + '_whisker.csv')):
                bad.append('%s: stale Camera copy without backup' % s)
    print('  %-24s %-5s copies %d, independent rule differs in %d frames, json levels %s' % (s, cod, len(cp), ndiff, (J.get('whisk_off_below_index', 'old fields: threshold_x_rest %s' % J.get('threshold_x_rest')), J.get('whisk_on_above_index', ''))))
print('  files compared with work: %d' % n_copy)

print('\n5. human-labelled clips scored with the saved whisk_on')
rows = []
for tag, cl, clip_s in (('set1', r'C:\Users\rbondarenko\DLC\whisking_check_h264\clips', 5.0), ('strat', r'C:\Users\rbondarenko\DLC\whisking_check_h264_strat\clips', 3.0)):
    key = pd.read_csv(os.path.join(cl, 'answer_key.csv'))
    ans = {}
    for r in csv.reader(open(os.path.join(cl, 'your_answers.csv'))):
        cells = [c.strip(' ;\t"').lower() for c in r[1:]]
        cells = [c for c in cells if c]
        if r and r[0].startswith('clip') and cells: ans[r[0]] = cells[0][0]
    for _, r in key.iterrows():
        a = ans.get(r['clip'])
        if a not in ('y', 'n'): continue
        T = pd.read_csv(os.path.join(WORK, r['session'], 'whisker_results', r['session'] + '_whisker.csv'), usecols=['whisk_on'])
        f0 = int(float(r['start_s']) * FPS)
        on = T['whisk_on'].fillna(0).values[f0:f0 + int(clip_s * FPS)].mean() > 0.5
        rows.append((tag, r['clip'], r['session'], a == 'y', bool(on)))
D = pd.DataFrame(rows, columns=['set', 'clip', 'session', 'you', 'code'])
print('  judged clips: %d  (%d whisking, %d still)' % (len(D), D.you.sum(), (~D.you).sum()))
print('  wrong: %d  (missed whisking %d, false whisking %d)' % ((D.you != D.code).sum(), (D.you & ~D.code).sum(), (~D.you & D.code).sum()))
for _, r in D[D.you != D.code].iterrows():
    print('    %s %s %s  you=%s code=%s' % (r['set'], r['clip'], r['session'], 'whisking' if r['you'] else 'still', 'whisking' if r['code'] else 'rest'))

print('\n%s' % ('ALL CHECKS PASSED' if not bad else 'PROBLEMS:'))
for b in bad: print('  ', b)
