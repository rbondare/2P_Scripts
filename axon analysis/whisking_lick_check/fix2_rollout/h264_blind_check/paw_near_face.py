"""Is 'forepaw close to the nose' (DLC body model) a usable grooming / paw-licking flag?  clip_13 (1423 at 2880.8 s) and the known
1549 grooming episodes (219.4-222.4, 465.2-468.1, 508.8-522.7 s of video time, old artifact flags)."""
import os, numpy as np, scipy.io as sio
from scipy.ndimage import gaussian_filter1d
here = os.path.dirname(os.path.abspath(__file__))
FPS = 50.0


def load(s):
    d = sio.loadmat(r'Z:\Group Members\Rima\DLCData\%s_DLC_filtered.mat' % s, variable_names=['Behav'], simplify_cells=True)['Behav']
    return d


def paw_dist(B):
    nose = B['noseTip'][:, :2].astype(float)
    out = {}
    for k in ('foreLimbRightTip', 'foreLimbRightWrist', 'foreLimbLeftTip', 'foreLimbLeftWrist'):
        if k in B:
            p = B[k][:, :2].astype(float)
            out[k] = np.hypot(*(p - nose).T)
    return out


for s, wins in (('AnimalRB24_261105_1423', {'clip_13 (paw licking, user)': (2880.8, 2885.8), 'clip_19 (whisking, user)': (2350.6, 2355.6), 'clip_16-like still 1744': None}),
                ('AnimalRB24_261105_1549', {'grooming 219.4-222.4': (219.4, 222.4), 'grooming 465.2-468.1': (465.2, 468.1), 'grooming 508.8-522.7': (508.8, 522.7), 'quiet 900-905': (900, 905)})):
    B = load(s)
    n = len(B['noseTip'])
    D = paw_dist(B)
    print('\n%s: %d frames; body-video px. nose NaN %.1f%%' % (s, n, 100 * np.isnan(B['noseTip'][:, 0]).mean()))
    for k, d in D.items():
        ok = ~np.isnan(d)
        print('  %-20s valid %.1f%%  distance to nose: p1 %.0f p5 %.0f p50 %.0f p95 %.0f' % (k, 100 * ok.mean(), *np.nanpercentile(d, [1, 5, 50, 95])))
    # nearest forepaw point
    stack = np.vstack(list(D.values()))
    near = np.nanmin(np.where(np.isnan(stack), np.inf, stack), axis=0)
    near[np.isinf(near)] = np.nan
    thr = np.nanpercentile(near, 3)
    print('  nearest forepaw point to nose: p3 %.0f px (3%% of frames closer than this)' % thr)
    for nm, w in wins.items():
        if w is None:
            continue
        a, b = int(w[0] * FPS), int(w[1] * FPS)
        seg = near[a:b]
        print('  %-28s nearest-paw distance: median %.0f, min %.0f, frames closer than session p3: %.0f%%, valid %.0f%%' % (
            nm, np.nanmedian(seg), np.nanmin(seg), 100 * np.nanmean(seg < thr), 100 * np.mean(~np.isnan(seg))))
    np.save(os.path.join(here, 'near_paw_%s.npy' % s[-9:]), near)
