"""Mean ME (gray levels / frame) of candidate sub-regions by lick x running class, from the saved pixel maps.  usage: python subbox.py <session>"""
import sys, os, numpy as np

sess = sys.argv[1]
here = os.path.dirname(os.path.abspath(__file__))
z = np.load(os.path.join(here, 'maps_%s.npz' % sess))
X0, Y0, BS, ax, ay = [int(z[k]) for k in ('X0', 'Y0', 'BS', 'ax', 'ay')]
# regions relative to the nose anchor: (dx0, dy0, dx1, dy1)
reg = {
    'fan_front (all)': (-260, -110, 60, 200),
    'front_top (above spout)': (-260, -110, 60, 10),
    'front_top_near (x>-160)': (-160, -110, 60, 10),
    'front_spout_band': (-260, 10, 60, 200),
    'fan_down (jaw/tongue)': (60, 80, 300, 300),
    'pad': (60, -100, 260, 60),
    'ctrl (noise)': (-320, -300, -200, -180),
}
cls = ['still_nolick', 'still_lick', 'run_nolick', 'run_lick']
print('%-26s' % 'region', *['%13s' % c for c in cls], '   lick-effect still / run (abs, %)')
for k, (a, b, c, d) in reg.items():
    row = []
    for cl in cls:
        m = z[cl]
        x0 = max((ax + a - X0) // BS, 0)
        x1 = min((ax + c - X0) // BS, m.shape[1])
        y0 = max((ay + b - Y0) // BS, 0)
        y1 = min((ay + d - Y0) // BS, m.shape[0])
        row.append(m[y0:y1, x0:x1].mean() if (x1 > x0 and y1 > y0) else np.nan)
    e_s = row[1] - row[0]
    e_r = row[3] - row[2]
    print('%-26s' % k, *['%13.2f' % v for v in row], '   %+.2f (%+.0f%%) / %+.2f (%+.0f%%)' % (e_s, 100 * e_s / row[0], e_r, 100 * e_r / row[2]))
