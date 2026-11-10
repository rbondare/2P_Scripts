"""Lick independence of the front_top box exactly as face_motion.py computes it (rows -110..0, csv me_front_top), 1506."""
import os, numpy as np, pandas as pd
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lick_whisk2.py'), encoding='utf-8').read()
exec(src.split("print('== (a)")[0])
defs = src.split("WL = int(FPS)")[1].split("w2 = build(2, 1.0)")[0]
exec('WL = int(FPS)' + defs)
P1 = r'C:\Users\RBONDA~1\AppData\Local\Temp\claude\e64e36f16637d91c\c193d28f-132e-4dcd-81fe-9333a3f48817\scratchpad\fm2\pass1'
N = pd.read_csv(os.path.join(P1, sess, sess + '_whisker.csv'))
ME['top_csv'] = N['me_front_top'].values.astype(np.float32)[:n]
print('frames', len(N), n)
print('corr(cell-grid front_top, csv front_top) = %.4f' % np.corrcoef(ME['front_top'], ME['top_csv'])[0, 1])
print('\n%-10s %10s %12s %10s' % ('signal', 'state only', 'state+lick', 'unique lick'))
w2 = build(2, 1.0)
for k in ('top_csv', 'front_top', 'old_index'):
    a = r2cv(ME[k], Zbase)
    b = r2cv(ME[k], np.hstack([Zbase, Zl]))
    print('%-10s %10.3f %12.3f %10.3f   matched lick effect (%% of loco contrast): all %.1f, 1st half %.1f, 2nd half %.1f' % (
        k, a, b, b - a, eff(w2, ME[k]), eff(w2, ME[k], lambda d: d.half == 0), eff(w2, ME[k], lambda d: d.half == 1)))
cl = ~near1 & good
run = cl & (speed > np.percentile(speed, 80))
stl = cl & (speed < np.percentile(speed, 20))
for k in ('top_csv', 'front_top', 'old_index', 'jaw'):
    s = gaussian_filter1d(np.asarray(ME[k], np.float64), 3)
    a, b = s[run], s[stl]
    print('Cohen d run vs still (lick-free) %-10s %.2f' % (k, (a.mean() - b.mean()) / np.sqrt(0.5 * (a.var() + b.var()))))
