"""What are the 'rejected' suite2p ROIs? Per-ROI feature table for every accepted and rejected ROI of a session (analysis behind Fig. 12).
usage: python rejected_rois.py [session ...]   (default: all unflagged sessions, 4 processes)  -> cache/rejroi_<session>.pkl

Per ROI (suite2p plane0 of the session, same F / Fneu / dF/F recipe for accepted and rejected ROIs; rejected: median Fc > 5 counts, at most CAP_REJ random ROIs):
  activity-likeness   skew of the 5 Hz best dF/F, V3 transient criterion on the grid (slow variance >= 4 x white-noise expectation, upward events >= 2 x downward events; approximate: accepted ROIs pass it 40-90%),
                      slow_ratio, up / down event counts
  neuropil-likeness   median F / median Fneu (puncta >> 1), r^2 of band-passed F vs Fneu
  shape / location    npix, compactness, aspect ratio (suite2p), local contrast (mean image inside the ROI / ring 2-6 px outside, accepted pixels excluded), overlap with and distance to accepted ROIs
  coupling            r (lag 0, band-passed 1-30 s) with whisking and with face motion and its circular-shift p (>= 120 s), for dF/F computed (a) as usual F - 0.5 Fneu, (b) from F alone, (c) from Fneu alone
"""
import os, sys, pickle, time, numpy as np
from fig_common import *
import xs_common as xc
from xs_common import DT
from scipy.ndimage import uniform_filter1d, binary_dilation
from scipy.stats import skew

SEP = int(120 / DT); CAP_REJ = 800; CAP_VAR = 220          # CAP_REJ >= the largest rejected set (621) so that the review queue covers every rejected ROI


def coupling(R, fy):
    r = np.full(len(R), np.nan); p = np.full(len(R), np.nan)
    for i, row in enumerate(R):
        if np.isfinite(row).mean() < 0.9: continue
        _, _, full = xc.xcorr_fft(row, fy, 1); n = len(full); r0 = full[0]; nul = full[SEP:n - SEP]; r[i] = r0; p[i] = (1 + (np.abs(nul) >= abs(r0)).sum()) / (1 + len(nul))
    return r, p


def v3_features(D, fs=1 / DT, thr=4.0, smooth_s=0.5):
    D = D.astype(np.float32); w = max(1, int(round(smooth_s * fs))); dsm = uniform_filter1d(D, w, axis=1, mode='nearest'); dev = dsm - np.median(dsm, axis=1, keepdims=True)
    neg = dev < 0; sig_neg = np.maximum(np.sqrt((dev ** 2 * neg).sum(1) / np.maximum(neg.sum(1), 1)), 1e-7); dd = np.diff(D, axis=1); hf = 1.4826 * np.median(np.abs(dd - np.median(dd, axis=1, keepdims=True)), axis=1) / np.sqrt(2)
    up = dev > 3 * sig_neg[:, None]; dn = dev < -3 * sig_neg[:, None]
    ev_pos = (np.diff(np.concatenate([np.zeros((len(up), 1), bool), up], 1).astype(np.int8), axis=1) == 1).sum(1); ev_neg = (np.diff(np.concatenate([np.zeros((len(dn), 1), bool), dn], 1).astype(np.int8), axis=1) == 1).sum(1)
    slow = dsm.var(axis=1) / (hf ** 2 / w + 1e-12)
    return dict(slow_ratio=slow, ev_pos=ev_pos, ev_neg=ev_neg, passes=(slow >= thr) & (ev_pos >= 2 * np.maximum(ev_neg, 1)))


def fill(D):
    med = np.nanmedian(D, axis=1, keepdims=True); return np.where(np.isfinite(D), D, med)


def analyse(name):
    t0 = time.time(); S = xc.Session(name, with_rejected=False, verbose=False); p = os.path.join(xc.DATA_DIR, name, 'suite2p', 'plane0')
    ic = np.load(os.path.join(p, 'iscell.npy'))[:, 0] > 0; acc = np.where(ic)[0]; rej_all = np.where(~ic)[0]
    if len(rej_all) == 0: return None
    rej = rej_all if len(rej_all) <= CAP_REJ else np.sort(np.random.default_rng(0).choice(rej_all, CAP_REJ, replace=False))
    F = np.load(os.path.join(p, 'F.npy'), mmap_mode='r'); Fn = np.load(os.path.join(p, 'Fneu.npy'), mmap_mode='r'); idx = np.concatenate([acc, rej]); grp = np.r_[np.zeros(len(acc), int), np.ones(len(rej), int)]
    Fs = np.asarray(F[idx]).astype(np.float32); Fns = np.asarray(Fn[idx]).astype(np.float32)
    # accepted ROIs of the xs file are the iscell ROIs in suite2p order (checked on the first ROI); signal-like flag from the pipeline criterion
    order_ok = bool(len(acc) == S.nroi and Fs.shape[1] == np.asarray(S.X['F']).shape[1] and np.corrcoef(np.asarray(S.X['F'])[0, :3000], Fs[0, :3000])[0, 1] > 0.99)
    sig_acc = np.zeros(len(acc), bool)
    if order_ok: sig_acc = S.sig.copy()
    Fg = xc.bin_rows(S.tca, Fs, S.edges); Fng = xc.bin_rows(S.tca, Fns, S.edges); dff, Fc, F0 = xc.best_dff_grid(Fg, Fng)
    keep = (grp == 0) | (np.nanmedian(Fc, 1) > 5); idx, grp, Fg, Fng, dff = idx[keep], grp[keep], Fg[keep], Fng[keep], dff[keep]
    sig_flag = np.zeros(len(idx), bool); sig_flag[:len(acc)] = sig_acc
    n = len(idx); out = dict(name=name, idx=idx, grp=grp, sig_acc=sig_flag, order_ok=order_ok, n_acc=int((grp == 0).sum()), n_rej=int((grp == 1).sum()))
    # ---- activity-likeness
    D = fill(dff); v3 = v3_features(D); out.update({k: v for k, v in v3.items()}); out['skew'] = skew(D, axis=1)
    # ---- neuropil-likeness
    mF = np.nanmedian(Fg, 1); mN = np.nanmedian(Fng, 1); out['F_over_Fneu'] = mF / np.maximum(mN, 1e-3)
    bF = xc.bandpass(Fg); bN = xc.bandpass(Fng); r2 = np.full(n, np.nan)
    for i in range(n):
        ok = np.isfinite(bF[i]) & np.isfinite(bN[i])
        if ok.sum() > 300: r2[i] = np.corrcoef(bF[i, ok], bN[i, ok])[0, 1] ** 2
    out['r2_F_Fneu'] = r2
    # ---- shape / location
    stat = np.load(os.path.join(p, 'stat.npy'), allow_pickle=True); ops = np.load(os.path.join(p, 'ops.npy'), allow_pickle=True).item(); mimg = np.asarray(ops['meanImg'], float); Ly, Lx = mimg.shape
    accmask = np.zeros((Ly, Lx), bool)
    for k in acc: accmask[stat[k]['ypix'], stat[k]['xpix']] = True
    accdil = binary_dilation(accmask, iterations=2); acc_med = np.array([stat[k]['med'] for k in acc], float)
    f = {k: np.full(n, np.nan) for k in ('npix', 'compact', 'aspect', 'overlap', 'adjacent', 'nn_dist', 'contrast', 'radius')}
    for j, k in enumerate(idx):
        s = stat[k]; y, x = s['ypix'], s['xpix']; f['npix'][j] = s['npix']; f['compact'][j] = s.get('compact', np.nan); f['aspect'][j] = s.get('aspect_ratio', np.nan); f['radius'][j] = s.get('radius', np.nan)
        if grp[j] == 1: f['overlap'][j] = accmask[y, x].mean(); f['adjacent'][j] = accdil[y, x].mean(); f['nn_dist'][j] = np.sqrt(((acc_med - np.asarray(s['med'], float)) ** 2).sum(1)).min()
        cy, cx = int(np.median(y)), int(np.median(x)); y0, y1, x0, x1 = max(cy - 25, 0), min(cy + 26, Ly), max(cx - 25, 0), min(cx + 26, Lx)
        m = np.zeros((y1 - y0, x1 - x0), bool); m[np.clip(y - y0, 0, y1 - y0 - 1), np.clip(x - x0, 0, x1 - x0 - 1)] = True; ring = binary_dilation(m, iterations=6) & ~binary_dilation(m, iterations=2) & ~accmask[y0:y1, x0:x1]
        if ring.sum() >= 10: f['contrast'][j] = mimg[y0:y1, x0:x1][m].mean() / max(mimg[y0:y1, x0:x1][ring].mean(), 1e-3)
    out.update(f); out['inmerge'] = np.array([float(stat[k].get('inmerge', 0)) for k in idx])       # > 0: hidden constituent of a merged ROI (still iscell 0): not an independent ROI
    # ---- coupling (and the neuropil test)
    behs = {b: xc.bandpass(S.beh[b]) for b in ('whisk', 'face') if b in S.beh}; out['beh'] = list(behs)
    sel = np.where((grp == 0) | (np.cumsum(grp == 1) <= CAP_VAR))[0]; out['sel_var'] = sel
    Z = np.zeros_like(Fng[sel]); dff_raw = xc.best_dff_grid(Fg[sel], Z)[0]; dff_neu = xc.best_dff_grid(Fng[sel], Z)[0]
    for b, fy in behs.items():
        out['r_' + b], out['p_' + b] = coupling(xc.bandpass(dff), fy)
        for tag, X_ in (('rawF', dff_raw), ('Fneu', dff_neu)):
            r, pv = coupling(xc.bandpass(X_), fy); rr = np.full(n, np.nan); pp = np.full(n, np.nan); rr[sel] = r; pp[sel] = pv; out['r_%s_%s' % (b, tag)] = rr; out['p_%s_%s' % (b, tag)] = pp
    out['seconds'] = time.time() - t0; return out


def run(name):
    f = os.path.join(CACHE, 'rejroi_%s.pkl' % name)
    if os.path.exists(f) and not os.environ.get('REFRESH'): return name, 'cached'
    try:
        r = analyse(name)
        if r is None: return name, 'no rejected ROIs'
        pickle.dump(r, open(f, 'wb')); return name, 'ok (%d accepted, %d rejected, order_ok %s, %.0f s)' % (r['n_acc'], r['n_rej'], r['order_ok'], r['seconds'])
    except Exception as e:
        import traceback; return name, 'FAILED: ' + traceback.format_exc()[-600:]


if __name__ == '__main__':
    names = sys.argv[1:] or UNFLAGGED
    if len(names) == 1: print(*run(names[0]))
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=4) as ex:
            for nm, msg in ex.map(run, names): print('%-26s %s' % (nm, msg), flush=True)
