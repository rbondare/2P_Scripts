"""State from a slower-smoothed index (gaussian sigma in s) with hysteresis, scored on the 51 judged clips; whole-session effect."""
import os, sys, json, csv, numpy as np, pandas as pd
from scipy.ndimage import gaussian_filter1d
sys.path.insert(0, r'C:\Users\rbondarenko\projects\2P_Scripts\DLC\face_model')
import face_motion as fm
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'strat_onrule.py')).read().split("cur = run(1.6, 0.1)")[0])
FPS = 50.0
base = {s: fm.classify(ME[s], 0.952, 'h264') for s in h264}
def slow_on(s, sig, hi, lo):
    x = gaussian_filter1d(base[s]['index'], sig * FPS) if sig > 0 else base[s]['index']
    on = fm.hysteresis(x, lo, hi)
    on = fm.clean_bouts(on, int(0.3 * FPS), int(0.3 * FPS))
    on[base[s]['artifact']] = False
    return on
print('%-34s %5s %6s %6s  ?-clip on-fractions        | whisking fraction median (min-max)' % ('sigma s, on above / off below', 'wrong', 'missed', 'false'))
q = L[L.you == '?']
out = []
for sig in (0, 0.5, 1.0, 1.5):
    for hi, lo in ((1.76, 1.44), (1.5, 1.3), (1.45, 1.35), (1.45, 1.4), (1.4, 1.3), (1.4, 1.35), (1.35, 1.3), (1.5, 1.4)):
        cache = {s: slow_on(s, sig, hi, lo) for s in h264}
        called = np.array([cache[r.session][int(r.start_s*FPS):int(r.start_s*FPS)+int(r.clip_s*FPS)].mean() > 0.5 for _, r in J.iterrows()])
        qf = [cache[r.session][int(r.start_s*FPS):int(r.start_s*FPS)+int(r.clip_s*FPS)].mean() for _, r in q.iterrows()]
        fr = np.array([cache[s][~base[s]['artifact']].mean() for s in h264])
        wrong = (called != y)
        out.append((sig, hi, lo, wrong.sum(), (y & ~called).sum(), (~y & called).sum(), np.median(fr), fr.min(), fr.max(), qf, list(J['clip'][wrong] + '(' + J.session[wrong].str[-9:] + ')')))
out.sort(key=lambda r: (r[3], r[0]))
for r in out[:14]:
    print('sigma %.1f, %.2f / %.2f %12d %6d %6d  %s | %.0f%% (%.0f-%.0f%%)' % (r[0], r[1], r[2], r[3], r[4], r[5], ' '.join('%.2f' % v for v in r[9]), 100*r[6], 100*r[7], 100*r[8]))
print('\nwrong clips for the best rule:', out[0][10])
print('current rule (sigma 0, 1.76/1.44):', [r[3:6] for r in out if r[0] == 0 and r[1] == 1.76])
