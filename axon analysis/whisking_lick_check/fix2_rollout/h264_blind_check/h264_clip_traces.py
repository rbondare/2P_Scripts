"""Traces of the new index and the jaw / pad motion in the user-flagged H.264 clips (2, 7, 8, 13) with two correct clips for comparison."""
import os, csv, numpy as np, pandas as pd
from scipy.ndimage import gaussian_filter1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
CL = r'C:\Users\rbondarenko\DLC\whisking_check_h264\clips'
key = {r['clip']: r for r in csv.DictReader(open(os.path.join(CL, 'answer_key.csv')))}
ans = {}
for r in csv.reader(open(os.path.join(CL, 'your_answers.csv'))):
    cells = [c.strip(' ;\t"').lower() for c in r[1:]]
    cells = [c for c in cells if c]
    if r and r[0] in key and cells:
        ans[r[0]] = cells[0][0]
show = ['clip_02', 'clip_07', 'clip_08', 'clip_13', 'clip_05', 'clip_12']
fig, axs = plt.subplots(len(show), 1, figsize=(12, 2.6 * len(show)))
for ax, c in zip(axs, show):
    s, st = key[c]['session'], float(key[c]['start_s'])
    T = pd.read_csv(os.path.join(r'C:\Users\rbondarenko\DLC\work', s, 'whisker_results', s + '_whisker.csv'))
    rest = 0.952
    f0 = int(st * 50)
    a, b = f0 - 150, f0 + 400
    sl = slice(a, b)
    t = (np.arange(a, b) - f0) / 50
    idx = T['whisk_index'].values[sl]
    sm = lambda k: gaussian_filter1d(T['me_' + k].values, 3)
    jaw, pad = sm('fan_down'), sm('pad')
    allj = np.nanmedian(jaw), np.nanmedian(pad)
    ax.plot(t, idx, c='tab:blue', lw=1.2, label='whisking index (front_top)')
    ax.plot(t, jaw[sl] / allj[0], c='tab:red', lw=.8, alpha=.7, label='jaw ME / session median')
    ax.plot(t, pad[sl] / allj[1], c='tab:orange', lw=.8, alpha=.7, label='pad ME / session median')
    ax.axhline(1.6, c='gray', ls='--', lw=.8)
    ax.axvspan(0, 5, color='gold', alpha=.12)
    on = T['whisk_on'].values[sl]
    art = T['artifact'].values[sl]
    ax.set_title('%s %s at %.0f s: code %s, you %s, mean index %.2f, session floor (p3) %.2f' % (
        c, s.replace('Animal', ''), st, key[c]['code_says'], ans[c], float(key[c]['mean_index']),
        np.percentile(gaussian_filter1d(T['me_front_top'].values, 3), 3) / rest), fontsize=9)
    ax.set_ylim(0, 4)
    ax.set_xlim(t[0], t[-1])
    if c == show[0]:
        ax.legend(fontsize=7, ncol=3, loc='upper right')
    # print numbers
    w = slice(f0, f0 + 250)
    print('%s %s: index mean %.2f (1st half %.2f, 2nd half %.2f, max %.2f); jaw x median %.2f; pad x median %.2f; artifact frac in window %.2f' % (
        c, ans[c], T['whisk_index'].values[w].mean(), T['whisk_index'].values[f0:f0 + 125].mean(), T['whisk_index'].values[f0 + 125:f0 + 250].mean(),
        T['whisk_index'].values[w].max(), jaw[w].mean() / allj[0], pad[w].mean() / allj[1], T['artifact'].values[w].mean()))
fig.tight_layout()
fig.savefig(os.path.join(here, 'h264_clip_traces.png'), dpi=70)
