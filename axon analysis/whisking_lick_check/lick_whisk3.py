"""Summary figure LW2 + a test whether the background-noise region can be used to clean front_top.  Re-uses lick_whisk2 helpers by exec of its data section."""
import os, re, numpy as np
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lick_whisk2.py'), encoding='utf-8').read()
head = src.split("print('== (a)")[0]
exec(head)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

names = ['old_index', 'jaw', 'pad', 'front_all', 'front_low', 'front_top']
lab = ['old index\n(front_all+jaw)', 'jaw\n(old fan_down)', 'pad', 'front_all\n(old fan_front)', 'front_low\n(spout band)', 'front_top\n(above spout)']
uniq, dc, dcl, dch, dd = [], [], [], [], []
import pandas as pd
R = pd.read_csv(os.path.join(here, 'lw_table.csv'), index_col=0)
tabname = {'old_index': 'old index (front_all + jaw)', 'jaw': 'jaw', 'pad': 'pad', 'front_all': 'front_all', 'front_low': 'front_low', 'front_top': 'front_top'}
for k in names:
    a = r2cv(ME[k], Zbase)
    b = r2cv(ME[k], np.hstack([Zbase, Zl]))
    uniq.append(b - a)
    r = R.loc[tabname[k]]
    dc.append(100 * r.lick_dc / r.loco_contrast)
    dcl.append(100 * r.lick_dc_lo / r.loco_contrast)
    dch.append(100 * r.lick_dc_hi / r.loco_contrast)
cl = ~near1 & good
run = cl & (speed > np.percentile(speed, 80))
stl = cl & (speed < np.percentile(speed, 20))
for k in names:
    s = gaussian_filter1d(np.asarray(ME[k], np.float64), 3)
    a, b = s[run], s[stl]
    dd.append((a.mean() - b.mean()) / np.sqrt(0.5 * (a.var() + b.var())))
col = ['#c0504d', '#c0504d', '#e0a040', '#e0a040', '#e0a040', '#2e8b57']
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
x = np.arange(len(names))
ax[0].bar(x, 100 * np.array(uniq), color=col)
ax[0].set_ylabel('% of variance explained by licking\n(held-out, beyond ball speed + pupil)')
ax[0].set_title('how much of each signal is lick-driven')
ax[1].bar(x, dc, color=col, yerr=[np.array(dc) - np.array(dcl), np.array(dch) - np.array(dc)], capsize=3)
ax[1].set_ylabel('lick effect at matched speed + pupil\n(% of running-vs-still contrast)')
ax[1].set_title('lick-window minus lick-free window (95% CI)')
ax[2].bar(x, dd, color=col)
ax[2].set_ylabel("running vs still, lick-free frames (Cohen's d)")
ax[2].set_title('whisking-relevant dynamic range kept')
for a in ax:
    a.set_xticks(x)
    a.set_xticklabels(lab, fontsize=8)
    a.axhline(0, c='k', lw=.5)
fig.suptitle('RB24_261109_1506 (4,734 licks, aux3 lick line): lick contamination of the motion-energy boxes')
fig.tight_layout()
fig.savefig(os.path.join(here, 'LW2_AnimalRB24_261109_1506.png'), dpi=70)

# noise-region cleaning of front_top: regress out the background ME (common-mode codec / illumination noise) and see whether run-vs-still improves
print('\nnoise-region cleaning (common-mode):')
nz = np.asarray(ME['noise'], np.float64)
for k in ('front_top', 'front_all', 'pad'):
    y = np.asarray(ME[k], np.float64)
    m = good & ~near1
    beta = np.polyfit(nz[m], y[m], 1)
    yc = y - beta[0] * (nz - nz[m].mean())
    for nm, s in (('raw', y), ('noise-cleaned', yc)):
        s = gaussian_filter1d(s, 3)
        a, b = s[run], s[stl]
        print('%-10s %-14s slope on noise %.2f  d(run vs still) = %.2f  SD in still frames %.3f' % (k, nm, beta[0], (a.mean() - b.mean()) / np.sqrt(0.5 * (a.var() + b.var())), b.std()))
print('corr(noise, front_top) lick-free, 0.06 s smoothed: %.2f' % np.corrcoef(gaussian_filter1d(nz, 3)[m], gaussian_filter1d(np.asarray(ME['front_top'], np.float64), 3)[m])[0, 1])
