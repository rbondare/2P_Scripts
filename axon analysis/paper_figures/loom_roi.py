"""ROI-level, trial-aligned loom data for Figure 10 (cached per session in plots/summary_figures/cache/loom_roi_<session>.pkl).
Per session and stimulus (loom_expand / loom_recede): signal-like ROIs x trials x time, best dF/F high-passed > 30 s (as in cross_session/stim_response.py), baseline (-2..0 s) subtracted, window -3..+10 s.
Trials of interrupted blocks are dropped exactly like in stim_response.py.   usage: python loom_roi.py [session ...]
"""
import os, sys, pickle, numpy as np
from fig_common import *
import xs_common as xc
from stim_response import hp, align, WIN, BWIN, AWIN

STIMS = ('loom_expand', 'loom_recede')


def build(name):
    f = os.path.join(CACHE, 'loom_roi_%s.pkl' % name)
    if os.path.exists(f): return pickle.load(open(f, 'rb'))
    S = xc.Session(name, with_rejected=False, verbose=False); X = S.X; out = dict(name=name, n_sig=int(S.sig.sum()), keys={})
    on_t = np.atleast_1d(np.asarray(X['on_t'], float)); on_key = np.atleast_1d(np.asarray(X['on_key'])).astype(str); on_blk = np.atleast_1d(np.asarray(X['on_blk'], int))
    binter = np.atleast_1d(np.asarray(X['blk_interrupted'], float)).astype(bool); keep = np.isfinite(on_t) & ~binter[on_blk - 1]
    on_t, on_key, on_blk = on_t[keep], on_key[keep], on_blk[keep]; e0 = S.edges[0]; D = np.stack([hp(r) for r in S.dff[S.sig]])
    for key in STIMS:
        m = on_key == key
        if m.sum() < 5: continue
        ons = on_t[m]; rel = None; A = []
        for r in D:
            rel, C = align(r, e0, ons, WIN[key]); A.append(C - np.nanmean(C[:, (rel >= BWIN[key][0]) & (rel < BWIN[key][1])], 1, keepdims=True))
        A = np.stack(A).astype(np.float32); out['keys'][key] = dict(rel=rel, A=A, blk=on_blk[m], on=ons)
    pickle.dump(out, open(f, 'wb')); return out


if __name__ == '__main__':
    R = load_pickle('stim_response_results.pkl')
    names = sys.argv[1:] or [n for n in ALL if any(k in R[n].get('keys', {}) for k in STIMS)]
    for n in names:
        d = build(n); print(n, {k: v['A'].shape for k, v in d['keys'].items()})
