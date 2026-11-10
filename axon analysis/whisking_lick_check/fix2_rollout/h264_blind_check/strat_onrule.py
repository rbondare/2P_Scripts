"""Re-run the on/off rule of face_motion.classify with other THRESH_X / HYSTERESIS on the H.264 sessions (ME traces from the csv),
score against the 51 judged blind clips, and show what changes in whole sessions. Nothing is written back."""
import os, sys, json, csv, numpy as np, pandas as pd
sys.path.insert(0, r'C:\Users\rbondarenko\projects\2P_Scripts\DLC\face_model')
import face_motion as fm
W = r'C:\Users\rbondarenko\DLC\work'
FPS = 50.0
codec = json.load(open(r'C:\Users\rbondarenko\projects\2P_Scripts\axon analysis\whisking_lick_check\fix2_rollout\codec_table.json'))
h264 = [r[0] for r in codec if r[1] == 'h264']
print('H.264 sessions:', len(h264))
ME = {}
for s in h264:
    p = os.path.join(W, s, 'whisker_results', s + '_whisker.csv')
    T = pd.read_csv(p)
    ME[s] = {k: T['me_' + k].values.astype(float) for k in fm.BOXES}
def labels(cl):
    key = pd.read_csv(os.path.join(cl, 'answer_key.csv'))
    ans = {}
    for r in csv.reader(open(os.path.join(cl, 'your_answers.csv'))):
        cells = [c.strip(' ;\t"').lower() for c in r[1:]]
        cells = [c for c in cells if c]
        if r and r[0].startswith('clip') and cells:
            ans[r[0]] = cells[0][0]
    key['you'] = key['clip'].map(ans)
    return key
L = []
for cl, cs in ((r'C:\Users\rbondarenko\DLC\whisking_check_h264\clips', 5.0), (r'C:\Users\rbondarenko\DLC\whisking_check_h264_strat\clips', 3.0)):
    k = labels(cl); k['clip_s'] = cs; L.append(k)
L = pd.concat(L, ignore_index=True)
J = L[L.you.isin(['y', 'n'])].reset_index(drop=True)
y = (J.you == 'y').values
def run(thr, hyst, rest=0.952):
    fm.ON_OFF['h264'] = (thr * (1 - hyst), thr * (1 + hyst))   # was fm.THRESH_X / HYSTERESIS before 2026-11-10
    return {s: fm.classify(ME[s], rest, 'h264') for s in h264}
cur = run(1.6, 0.1)
print('\n%-28s %5s %6s %7s %7s | %s' % ('rule (on above / off below)', 'wrong', 'missed', 'false', 'unjudged-?', 'H.264 sessions: whisking fraction of non-artifact time (median, min-max), bouts/min'))
rows = []
for thr, h in ((1.6, 0.1), (1.5, 0.1), (1.45, 0.1), (1.4, 0.1), (1.35, 0.1), (1.5, 0.05), (1.45, 0.05), (1.4, 0.05), (1.35, 0.05), (1.3, 0.05)):
    R = run(thr, h)
    called = []
    for _, r in J.iterrows():
        a, n = int(r.start_s * FPS), int(r.clip_s * FPS)
        called.append(R[r.session]['on'][a:a + n].mean() > 0.5)
    called = np.array(called)
    q = L[L.you == '?']
    qf = [R[r.session]['on'][int(r.start_s*FPS):int(r.start_s*FPS)+int(r.clip_s*FPS)].mean() for _, r in q.iterrows()]
    fr = np.array([R[s]['on'][~R[s]['artifact']].mean() for s in h264])
    from face_motion import runs
    bpm = np.array([len(runs(R[s]['on'])) / (len(R[s]['on']) / FPS / 60) for s in h264])
    print('%.2f / %.2f (thr %.2f, h %.2f) %5d %6d %7d   ?-clips on-fraction %s | %.0f%% (%.0f-%.0f%%), %.1f bouts/min' % (
        thr * (1 + h), thr * (1 - h), thr, h, (called != y).sum(), (y & ~called).sum(), (~y & called).sum(), ' '.join('%.2f' % v for v in qf),
        100 * np.median(fr), 100 * fr.min(), 100 * fr.max(), np.median(bpm)))
    rows.append((thr, h, (called != y).sum(), np.median(fr)))
# per-session fractions: current vs two candidates
print('\nper session whisking fraction (current 1.6/0.1  ->  1.45/0.1  ->  1.4/0.05)')
a = run(1.45, 0.1); b = run(1.4, 0.05)
for s in h264:
    f = lambda R: 100 * R[s]['on'][~R[s]['artifact']].mean()
    print('  %-26s %5.1f -> %5.1f -> %5.1f' % (s.replace('Animal', ''), f(cur), f(a), f(b)))
