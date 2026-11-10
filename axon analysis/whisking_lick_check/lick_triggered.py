"""Lick-triggered averages of the face-motion ME boxes + whisking index, for spout sessions and a no-spout control.
usage: python lick_triggered.py <session>   (needs xs/<session>_task.mat, whisker csv, aggregated TimeCamera)"""
import sys, os, numpy as np, pandas as pd, h5py, scipy.io as sio

sess = sys.argv[1]
W = r'C:\Users\rbondarenko\DLC\work\%s\whisker_results\%s_whisker.csv' % (sess, sess)
A = r'Z:\Group Members\Rima\Aggregated\%s_preprocessed.mat' % sess
T = sio.loadmat(r'C:\Users\rbondarenko\2P_cache\xs\%s_task.mat' % sess, squeeze_me=True)
lick = np.atleast_1d(T['lick_t']).astype(float)
on_t = np.atleast_1d(T['on_t']).astype(float)
d = pd.read_csv(W)
with h5py.File(A, 'r') as f:
    tc = f['Triggers/TimeCamera'][:].ravel()
n = min(len(d), len(tc))
tc = tc[:n]
fps = 50.0
print(sess, 'frames', n, 'licks', len(lick), 'first/last lick', lick[:3], lick[-3:])

# lick bouts: first lick after >= 1 s without licks
gap = np.diff(lick, prepend=-1e9)
bout_on = lick[gap >= 1.0]
print('lick bouts', len(bout_on))


def idx_of(t):
    i = np.searchsorted(tc, t)
    ok = (i > 0) & (i < n) & (np.abs(tc[np.clip(i, 0, n - 1)] - t) < 0.05)
    return i, ok


cols = ['me_fan_front', 'me_fan_down', 'me_pad', 'me_ctrl', 'whisk_index']
rel = np.arange(-int(2 * fps), int(3 * fps) + 1)
i0, ok = idx_of(bout_on)
i0 = i0[ok]
i0 = i0[(i0 + rel[0] >= 0) & (i0 + rel[-1] < n)]
print('bouts with camera', len(i0))
out = {}
for c in cols:
    x = d[c].values[:n].astype(float)
    M = np.stack([x[i + rel] for i in i0])
    base = np.nanmedian(M[:, : int(1.5 * fps)], axis=1, keepdims=True)
    out[c] = (np.nanmedian(M, axis=0), np.nanmedian(M - base, axis=0))
    pk = np.nanmedian(M[:, int(2 * fps): int(2 * fps) + int(1.5 * fps)])
    pre = np.nanmedian(M[:, : int(1.5 * fps)])
    print('%-13s median before bout (-2..-0.5 s) %.2f | median 0..1.5 s after first lick %.2f | ratio %.2f' % (c, pre, pk, pk / pre))

# the index is nonsensical to interpret without a lick-free control with the SAME running state: compare lick bouts to random times
rng = np.random.default_rng(0)
lick_all = np.zeros(n, bool)
il, okl = idx_of(lick)
lick_all[il[okl]] = True
near = np.convolve(lick_all.astype(float), np.ones(int(1.0 * fps)), 'same') > 0
print('fraction of frames within 0.5 s of a lick: %.3f' % near.mean())
for c in cols:
    x = d[c].values[:n].astype(float)
    print('%-13s median with lick %.2f  without lick %.2f' % (c, np.nanmedian(x[near]), np.nanmedian(x[~near])))
np.savez(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lt_%s.npz' % sess), rel=rel / fps, **{c + '_abs': out[c][0] for c in cols}, **{c + '_dlt': out[c][1] for c in cols}, bout_on=bout_on, i0=i0)
