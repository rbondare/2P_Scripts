"""Running-speed PROXY from the body video, for sessions whose ball log is truncated (open-loop blocks, see memory ballcam-truncation-rb22-24).

Features per time bin: frame-difference motion energy of the 480x264 DLC body video on a 6x8 grid (48 cells) + suite2p rigid brain-shift speed.
Model: ridge on log1p(features) -> log1p(ball speed in mm/s), fitted on REFERENCE session(s) with a complete ball log, then applied to the target.
The reference must be recorded with the same camera framing as the target (the video cells are position-specific): the same animal and day. A
reference from another day moved the activity map by several cells, and the prediction left the range of the training data (RB24_261023_1527 with the
RB24_261030_1744 reference: speeds to 5000 mm/s). Features are therefore clipped to the training range and the prediction to the 99.5th percentile of the
training speeds.

Validated (RB24_261023_1527 target, RB24_261023_1538 corridor block 3 as reference, 0.5 s bins): blocked 5-fold CV inside the reference r(log) 0.97,
running / still (> 20 mm/s) AUC 1.00; first ~9 s of true readings of the 8 open-loop blocks of 1538 (not used for fitting) r(log1p) 0.98 over 62 s, running /
still agreement 0.98; the 7 s of true readings at the start of 1527: median 3.0 mm/s true vs 3.4 predicted. It over-estimates slow speeds (head / face / limb
movement also moves the video), gives a scalar speed only, and its absolute scale beyond ~150 mm/s is not checked.

Workflow
  1. MATLAB:  export_proxy_inputs(<reference>, work_dir, <ball block>)  and  export_proxy_inputs(<target>, work_dir, 1)   (-> <session>_proxyin.mat)
  2. python (an environment with cv2, scipy; e.g. miniforge3\\envs\\dlc-tf):
       python locomotion_proxy_from_video.py --work work_dir --target AnimalRB24_261023_1527 --ref AnimalRB24_261023_1538:3
     computes the video features (cached as me_grid_<session>.npy in work_dir; ~25 s per 5 min of video), fits, applies, and writes
     <target>_locproxy.mat (t, speed_mm_s, bin_s) into this folder, where axon_analysis_spont_V1 looks for it.
"""
import argparse, os, time
import numpy as np
import scipy.io as sio

GR, GC = 6, 8
HERE = os.path.dirname(os.path.abspath(__file__))


def me_grid(session, work):
    """frame-difference motion energy on a GR x GC grid of the body video, [frames x 48]; cached"""
    f = os.path.join(work, f'me_grid_{session}.npy')
    if os.path.isfile(f):
        return np.load(f).astype(float)
    import cv2
    animal = session.split('_')[0]
    root = rf'Z:\Group Members\Rima\preCamera\{animal}\{session}\body_results\{session}_body_%05d.mp4'
    out, prev, t0 = [], None, time.time()
    for i in range(10):
        if not os.path.isfile(root % i):
            break
        cap = cv2.VideoCapture(root % i)
        while True:
            ok, fr = cap.read()
            if not ok:
                break
            g = fr[:, :, 0].astype(np.float32) if fr.ndim == 3 else fr.astype(np.float32)
            if prev is None:
                out.append(np.zeros(GR * GC, np.float32))
            else:
                d = np.abs(g - prev)
                H, W = d.shape
                out.append(d[:H // GR * GR, :W // GC * GC].reshape(GR, H // GR, GC, W // GC).mean(axis=(1, 3)).ravel().astype(np.float32))
            prev = g
        cap.release()
    me = np.stack(out)
    np.save(f, me)
    print(f'{session}: video features {me.shape} in {time.time() - t0:.0f} s')
    return me.astype(float)


def binned(t, v, edges):
    """mean of v inside each bin (NaN where empty)"""
    idx = np.digitize(t, edges) - 1
    ok = (idx >= 0) & (idx < len(edges) - 1) & np.isfinite(v)
    s = np.bincount(idx[ok], v[ok], len(edges) - 1)
    c = np.bincount(idx[ok], minlength=len(edges) - 1)
    with np.errstate(invalid='ignore', divide='ignore'):
        return s / c


def features(session, work, edges):
    """[bins x 49] log1p features: 48 video cells + rigid brain-shift speed (px/frame), on the 2P clock"""
    g = sio.loadmat(os.path.join(work, f'{session}_proxyin.mat'), squeeze_me=True)
    tcam, timeca = g['tcam'].astype(float), g['timeca'].astype(float)
    me = me_grid(session, work)
    n = min(len(me), len(tcam))
    bm = np.r_[0, np.hypot(np.diff(g['xoff'].astype(float)), np.diff(g['yoff'].astype(float)))]
    X = np.column_stack([binned(tcam[:n], me[:n, c], edges) for c in range(GR * GC)] + [binned(timeca[:len(bm)], bm, edges)])
    return np.log1p(np.clip(X, 0, None)), g


def ball_speed(g, edges):
    """ball speed (mm/s) per bin from the readings of the exported block (increments per trigger / dt)"""
    t = np.atleast_1d(g['t_ball']).astype(float)
    n = int(g['n_read']) if 'n_read' in g else len(t)
    t, fwd, side = t[:n], np.atleast_1d(g['fwd'])[:n], np.atleast_1d(g['side'])[:n]
    dt = np.diff(t, append=t[-1] + np.median(np.diff(t)))
    return binned(t, np.hypot(fwd, side) / dt, edges)


def ridge(X, y, lam):
    return np.linalg.solve(X.T @ X + lam * np.eye(X.shape[1]), X.T @ y)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--work', required=True, help='folder with <session>_proxyin.mat (and the cached video features)')
    ap.add_argument('--target', required=True, help='session whose speed is estimated, e.g. AnimalRB24_261023_1527')
    ap.add_argument('--ref', action='append', required=True, help='reference session:block with a complete ball log, e.g. AnimalRB24_261023_1538:3 (repeatable)')
    ap.add_argument('--bin', type=float, default=0.5, help='bin (s), default 0.5')
    ap.add_argument('--lam', type=float, default=100.0)
    ap.add_argument('--out', default=HERE, help='folder for <target>_locproxy.mat')
    a = ap.parse_args()
    B = a.bin

    Xs, ys = [], []
    for r in a.ref:
        sess = r.split(':')[0]
        g0 = sio.loadmat(os.path.join(a.work, f'{sess}_proxyin.mat'), squeeze_me=True)
        tb = np.atleast_1d(g0['t_ball']).astype(float)
        edges = np.arange(tb[0], tb[-1], B)
        X, g = features(sess, a.work, edges)
        y = ball_speed(g, edges)
        ok = np.isfinite(y) & np.all(np.isfinite(X), 1)
        Xs.append(X[ok]); ys.append(np.log1p(y[ok]))
        print(f'reference {r}: {ok.sum()} bins of {B} s, speed median {np.median(y[ok]):.1f} mm/s')
    X, y = np.vstack(Xs), np.concatenate(ys)
    lo, hi = np.percentile(X, 0.5, 0), np.percentile(X, 99.5, 0)
    Xc = np.clip(X, lo, hi)
    mu, sd = Xc.mean(0), Xc.std(0) + 1e-9
    w = ridge(np.c_[(Xc - mu) / sd, np.ones(len(Xc))], y, a.lam)
    ymax = np.percentile(y, 99.5)

    def predict(Xn):
        Z = (np.clip(Xn, lo, hi) - mu) / sd
        return np.expm1(np.minimum(np.c_[Z, np.ones(len(Z))] @ w, ymax))

    gt = sio.loadmat(os.path.join(a.work, f'{a.target}_proxyin.mat'), squeeze_me=True)
    edges = np.arange(0, gt['timeca'].astype(float)[-1], B)
    Xn, _ = features(a.target, a.work, edges)
    okn = np.all(np.isfinite(Xn), 1)
    speed = np.full(len(edges) - 1, np.nan)
    speed[okn] = predict(Xn[okn])
    tc = edges[:-1] + B / 2
    print(f'{a.target}: {okn.sum()} of {len(okn)} bins, predicted speed median {np.nanmedian(speed):.1f}, p90 {np.nanpercentile(speed, 90):.1f} mm/s, '
          f'share of bins > 20 mm/s {np.nanmean(speed[okn] > 20):.2f}')
    if 't_ball' in gt:                                         # the target's own true readings (a truncated log: only the first seconds)
        e2 = np.arange(np.atleast_1d(gt['t_ball'])[0], np.atleast_1d(gt['t_ball'])[min(int(gt['n_read']), len(np.atleast_1d(gt['t_ball']))) - 1], 1.0)
        yt = ball_speed(gt, e2)
        pp = np.array([np.nanmean(speed[(tc >= s) & (tc < s + 1)]) if np.any(np.isfinite(speed[(tc >= s) & (tc < s + 1)])) else np.nan for s in e2[:-1]])
        m = np.isfinite(yt) & np.isfinite(pp)
        if m.sum() > 2:
            print(f'check on the target\'s own {m.sum()} s of true readings: median true {np.median(yt[m]):.1f} vs predicted {np.median(pp[m]):.1f} mm/s')
    out = os.path.join(a.out, f'{a.target}_locproxy.mat')
    sio.savemat(out, {'t': tc, 'speed_mm_s': speed, 'bin_s': B, 'reference': np.array(a.ref, dtype=object)}, do_compression=True)
    print('written', out)


if __name__ == '__main__':
    main()
