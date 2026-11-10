import numpy as np, pandas as pd, os, h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
sess = 'AnimalRB24_261109_1506'
A = np.load(os.path.join(here, 'align_%s.npz' % sess))
tcam, lick, tball = A['aux1'], A['aux3'], A['aux2']
d = pd.read_csv(r'C:\Users\rbondarenko\DLC\work\%s\whisker_results\%s_whisker.csv' % (sess, sess))
print('camera frames', len(tcam), 'csv frames', len(d), 'duration %.1f s' % (tcam[-1] - tcam[0]))
# lick events: drop chatter (< 30 ms after the previous event)
keep = np.r_[True, np.diff(lick) > 0.03]
print('lick events %d, after removing <30 ms chatter %d' % (len(lick), keep.sum()))
lk = lick[keep]
ili = np.diff(lk)
print('ILI percentiles 5/25/50/75/95 (s):', np.round(np.percentile(ili, [5, 25, 50, 75, 95]), 3))
bout_on = lk[np.r_[True, ili >= 1.0]]
bout_len = []
i0 = np.flatnonzero(np.r_[True, ili >= 1.0])
i1 = np.r_[i0[1:], len(lk)]
nl = i1 - i0
print('bouts (gap >= 1 s):', len(i0), ' licks per bout median %d, p90 %d, max %d' % (np.median(nl), np.percentile(nl, 90), nl.max()))
dur = lk[i1 - 1] - lk[i0]
print('bout duration median %.1f s, p90 %.1f s' % (np.median(dur), np.percentile(dur, 90)))
t_rel = lk - tcam[0]
# fraction of camera frames within +-0.5 s of a lick
fr = np.searchsorted(tcam, lk)
near = np.zeros(len(tcam), bool)
for f in fr:
    near[max(f - 25, 0):f + 26] = True
print('fraction of frames within 0.5 s of a lick: %.3f' % near.mean())
fig, ax = plt.subplots(3, 1, figsize=(15, 8), sharex=True)
ax[0].hist(ili[ili < 0.6], 120)
fig2, ax2 = plt.subplots(3, 1, figsize=(15, 7), sharex=True)
tt = np.arange(len(d)) / 50.0
sm = pd.Series(d.whisk_index.values).rolling(250, center=True, min_periods=1).median().values
ax2[0].plot(tt[::25], sm[::25], lw=.7)
ax2[0].set_ylabel('whisk index (5 s median)')
binw = 10.0
edges = np.arange(0, tt[-1] + binw, binw)
rate = np.histogram(t_rel, edges)[0] / binw
ax2[1].plot(edges[:-1] + binw / 2, rate, lw=.8)
ax2[1].set_ylabel('licks / s (10 s bins)')
ax2[2].plot(tt[::25], d.me_fan_down.values[::25], lw=.4, label='fan_down ME')
ax2[2].plot(tt[::25], d.me_fan_front.values[::25], lw=.4, label='fan_front ME')
ax2[2].legend()
ax2[2].set_xlabel('time (s)')
fig2.savefig(os.path.join(here, 'lick_overview_%s.png' % sess), dpi=70)
np.save(os.path.join(here, 'lick_clean_%s.npy' % sess), lk)
