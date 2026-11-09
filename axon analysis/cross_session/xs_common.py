"""Shared helpers of the cross-session axon analyses (state signal, stimulus responses).

Everything runs on ONE common 5 Hz grid (DT = 0.2 s) so that 30 Hz single-plane and 5 Hz 5-plane sessions are treated alike.
Input: <XS_DIR>\\<session>_xs.mat from export_xs.m (+ the suite2p folder of the session for the rejected ROIs and the rigid shifts).
dF/F: (Fc - F0) / max(F0, 5), Fc = F - 0.5 Fneu, F0 = running 10th percentile of Fc over +-60 s (on the grid), smoothed 10 s  (BestDFF recipe of 1317 / 1529)
      and the pipeline dF/F (V3 definition) binned to the grid for comparison.
ROI classes: signal-like (V3 criterion, from the native pipeline dF/F), noise-like, neuropil-limited (Fc <= 0 in > 1% of frames),
      rejected = suite2p ROIs the experimenter did not accept (control for common-mode artefacts: sign of the behaviour coupling is ~50/50 there).
"""
import os, numpy as np, scipy.io as sio
from scipy.ndimage import percentile_filter, gaussian_filter1d, uniform_filter1d

XS_DIR = r'C:\Users\rbondarenko\2P_cache\xs'
DATA_DIR = r'Z:\Group Members\Rima\DATA_2P'
DT = 0.2
NEUROPIL, FLOOR = 0.5, 5.0
BASE_P, BASE_WIN_S, BASE_SMOOTH_S = 10.0, 60.0, 10.0


def load_xs(session):
    d = sio.loadmat(os.path.join(XS_DIR, session + '_xs.mat'), squeeze_me=True)
    return {k: v for k, v in d.items() if not k.startswith('__')}


def make_grid(tca):
    t0 = float(tca[0]); n = int(np.floor((tca[-1] - t0) / DT))
    edges = t0 + DT * np.arange(n + 1)
    return edges, edges[:-1] + DT / 2


def bin_rows(t, X, edges):
    """mean of the finite samples of every row of X (rows x len(t), t sorted) in the bins between `edges`; NaN for empty bins"""
    X = np.atleast_2d(X); nb = len(edges) - 1
    idx = np.searchsorted(edges, t, side='right') - 1
    ok0 = (idx >= 0) & (idx < nb)
    out = np.full((X.shape[0], nb), np.nan, np.float32)
    for i, row in enumerate(X):
        ok = ok0 & np.isfinite(row)
        if not ok.any(): continue
        s = np.bincount(idx[ok], weights=row[ok].astype(np.float64), minlength=nb); c = np.bincount(idx[ok], minlength=nb)
        with np.errstate(invalid='ignore', divide='ignore'):
            out[i] = np.where(c > 0, s / c, np.nan)
    return out


def nsmooth(x, sigma_bins):
    """Gaussian smoothing that ignores NaN (normalised convolution); returns NaN where no finite sample is within reach"""
    m = np.isfinite(x).astype(float); y = np.where(np.isfinite(x), x, 0.0)
    num = gaussian_filter1d(y, sigma_bins, axis=-1, mode='nearest'); den = gaussian_filter1d(m, sigma_bins, axis=-1, mode='nearest')
    with np.errstate(invalid='ignore', divide='ignore'):
        r = num / den
    r[den < 0.05] = np.nan
    return r


def bandpass(x, lo_s=1.0, hi_s=30.0):
    """slow drifts (> ~hi_s) removed, fast noise (< ~lo_s) smoothed; sigma in seconds, NaN-safe"""
    return nsmooth(x, lo_s / DT) - nsmooth(x, hi_s / DT)


def best_dff_grid(F, Fn):
    """F, Fn: rows x bins on the grid (NaN allowed) -> dF/F (BestDFF recipe), Fc, F0"""
    Fc = F - NEUROPIL * Fn
    w = int(2 * BASE_WIN_S / DT) + 1
    F0 = np.empty_like(Fc)
    for i in range(Fc.shape[0]):
        r = Fc[i].copy(); bad = ~np.isfinite(r)
        if bad.all(): F0[i] = np.nan; continue
        if bad.any(): r[bad] = np.interp(np.where(bad)[0], np.where(~bad)[0], r[~bad])
        b = percentile_filter(r, BASE_P, size=w, mode='nearest'); F0[i] = gaussian_filter1d(b, BASE_SMOOTH_S / DT, mode='nearest')
    with np.errstate(invalid='ignore', divide='ignore'):
        dff = (Fc - F0) / np.maximum(F0, FLOOR)
    return dff.astype(np.float32), Fc, F0


def signal_like_v3(dFF, Fc_native, fs, slow_ratio_thr=4.0, smooth_s=0.5):
    """V3 section 5 criterion on the native pipeline dF/F (rows = ROIs): slow variance vs white-noise expectation, upward events, Fc > 0"""
    dFF = dFF.astype(np.float32); w = max(1, int(round(smooth_s * fs)))
    dsm = uniform_filter1d(dFF, w, axis=1, mode='nearest'); med = np.median(dsm, axis=1, keepdims=True); dev = dsm - med
    neg = dev < 0; sig_neg = np.maximum(np.sqrt((dev ** 2 * neg).sum(1) / np.maximum(neg.sum(1), 1)), 1e-7)
    dd = np.diff(dFF, axis=1); hf = 1.4826 * np.median(np.abs(dd - np.median(dd, axis=1, keepdims=True)), axis=1) / np.sqrt(2)
    up = dev > 3 * sig_neg[:, None]; dn = dev < -3 * sig_neg[:, None]
    ev_pos = (np.diff(np.concatenate([np.zeros((len(up), 1), bool), up], 1).astype(np.int8), axis=1) == 1).sum(1)
    ev_neg = (np.diff(np.concatenate([np.zeros((len(dn), 1), bool), dn], 1).astype(np.int8), axis=1) == 1).sum(1)
    slow_ratio = dsm.var(axis=1) / (hf ** 2 / w + 1e-12)
    frac_le0 = (Fc_native <= 0).mean(1)
    neuropil = frac_le0 > 0.01
    sig = (~neuropil) & (slow_ratio >= slow_ratio_thr) & (ev_pos >= 2 * np.maximum(ev_neg, 1))
    return sig, neuropil, slow_ratio


class Session:
    """one session on the common grid"""

    def __init__(self, name, with_rejected=True, verbose=True):
        self.name = name; X = load_xs(name); self.X = X
        self.tca = np.asarray(X['tca'], float); self.fs = 1 / np.median(np.diff(self.tca))
        self.edges, self.t = make_grid(self.tca); self.nb = len(self.t); self.dur = self.t[-1] - self.t[0]
        F = np.atleast_2d(X['F']).astype(np.float32); Fn = np.atleast_2d(X['Fneu']).astype(np.float32); dp = np.atleast_2d(X['dFF_pipe']).astype(np.float32)
        self.nroi = F.shape[0]
        self.sig, self.neuropil, self.slow_ratio = signal_like_v3(dp, F - NEUROPIL * Fn, self.fs)
        self.noise = (~self.sig) & (~self.neuropil)
        Fg = bin_rows(self.tca, F, self.edges); Fng = bin_rows(self.tca, Fn, self.edges)
        self.dff, self.Fc, self.F0 = best_dff_grid(Fg, Fng)
        self.dff_pipe = bin_rows(self.tca, dp, self.edges)
        self.Fg = Fg; self.Fng = Fng
        self.rej = None
        if with_rejected:
            try: self._load_rejected()
            except Exception as e:
                if verbose: print('  rejected ROIs not available:', e)
        self._behaviour()
        self._shift()
        if verbose: print('%s: %d ROIs (%d signal-like, %d noise-like, %d neuropil-limited), %.1f min, fs %.2f Hz%s' % (name, self.nroi, self.sig.sum(), self.noise.sum(), self.neuropil.sum(), self.dur / 60, self.fs,
                                                                                        '' if self.rej is None else ', %d rejected ROIs' % self.rej['dff'].shape[0]))

    # ---- rejected suite2p ROIs (control population)
    def _load_rejected(self):
        p = os.path.join(DATA_DIR, self.name, 'suite2p', 'plane0')
        ic = np.load(os.path.join(p, 'iscell.npy'))[:, 0] > 0
        F = np.load(os.path.join(p, 'F.npy'), mmap_mode='r'); Fn = np.load(os.path.join(p, 'Fneu.npy'), mmap_mode='r')
        if F.shape[1] != len(self.tca): raise ValueError('suite2p F has %d frames, session %d' % (F.shape[1], len(self.tca)))
        rid = np.where(~ic)[0]
        if len(rid) == 0: raise ValueError('no rejected ROIs')
        if len(rid) > 400: rid = np.random.default_rng(0).choice(rid, 400, replace=False); rid.sort()
        Fg = bin_rows(self.tca, np.asarray(F[rid]), self.edges); Fng = bin_rows(self.tca, np.asarray(Fn[rid]), self.edges)
        dff, Fc, F0 = best_dff_grid(Fg, Fng)
        ok = np.nanmedian(Fc, 1) > 3 * np.sqrt(np.maximum(np.nanmedian(F0, 1), 1)) * 0 + 5   # drop ROIs with ~no signal at all (Fc <= 5 counts): dF/F is meaningless there
        self.rej = dict(idx=rid[ok], dff=dff[ok], F0=np.nanmedian(F0[ok], 1))

    # ---- behaviour on the grid
    def _behaviour(self):
        X = self.X; self.beh = {}; self.beh_info = {}
        if int(np.asarray(X.get('has_behav', 0)).ravel()[0]) == 1 and np.size(X.get('tcam')) > 10:
            tc = np.asarray(X['tcam'], float)
            def g(name, log=False, mask=None):
                v = np.asarray(X.get(name, []), float)
                if v.size != tc.size: return None
                v = v.copy()
                if mask is not None: v[mask] = np.nan
                if log: v = np.log(np.maximum(v, 1e-6))
                return bin_rows(tc, v, self.edges)[0]
            bad = np.asarray(X.get('PupilBad', []), float)
            pm = (bad > 0) if bad.size == tc.size else None
            for nm, key, log, mk in (('pupil', 'PupilArea', False, pm), ('face', 'FacePadME', True, None), ('whisk', 'WhiskIndex', True, None), ('eyelid', 'EyeLidOpening', False, None)):
                v = g(key, log, mk)
                if v is not None and np.isfinite(v).mean() > 0.3: self.beh[nm] = v
            # saccade-free pupil dilation rate: derivative of the 1 s smoothed pupil
            if 'pupil' in self.beh:
                sp = nsmooth(self.beh['pupil'], 1.0 / DT); self.beh['dpupil'] = np.gradient(sp) / DT
        # ball: speed per bin from the readings, NaN where there are none
        lt = np.atleast_1d(np.asarray(X.get('loco_t', []), float)); self.loco_cov = 0.0; self.loco_ok = False; self.loco_note = 'no ball readings'
        if lt.size > 100:
            fw = np.asarray(X['loco_fwd'], float); sd = np.asarray(X['loco_side'], float)
            idx = np.searchsorted(self.edges, lt, side='right') - 1; ok = (idx >= 0) & (idx < self.nb) & np.isfinite(fw) & np.isfinite(sd)
            sf = np.bincount(idx[ok], fw[ok], minlength=self.nb); ss = np.bincount(idx[ok], sd[ok], minlength=self.nb); cn = np.bincount(idx[ok], minlength=self.nb)
            sp = np.where(cn > 0, np.hypot(sf, ss) / DT, np.nan)
            self.loco_cov = float(np.isfinite(sp).mean()); run_frac = float(np.nanmean(sp > 50)) if np.isfinite(sp).any() else np.nan
            self.loco_run_frac = run_frac
            # sanity: a mouse does not run 100% of 100 min; and the ball must relate to face motion (checked in the analysis, not here)
            self.loco_ok = bool(self.loco_cov >= 0.7 and run_frac < 0.7)     # 1317 / 1529 / RB22_261023_1427 read 79-82% 'running' for 20-100 min: implausible, ball not used
            self.loco_note = 'coverage %.0f%%, running %.0f%% of covered time' % (100 * self.loco_cov, 100 * run_frac)
            self.beh['loco_raw'] = np.log1p(sp)
            if self.loco_ok: self.beh['loco'] = np.log1p(sp)

    def _shift(self):
        self.shift = None
        try:
            p = os.path.join(DATA_DIR, self.name, 'suite2p', 'plane0', 'ops.npy')
            ops = np.load(p, allow_pickle=True).item(); xo, yo = ops['xoff'].astype(float), ops['yoff'].astype(float)
            if len(xo) == len(self.tca):
                w = int(100 * self.fs) | 1
                sh = np.hypot(xo - uniform_filter1d(xo, w, mode='nearest'), yo - uniform_filter1d(yo, w, mode='nearest'))
                self.shift = bin_rows(self.tca, sh, self.edges)[0]
        except Exception:
            pass
        if self.shift is not None: self.beh['shift'] = self.shift

    # ---- convenience
    def population(self, which='sig', dff='best'):
        D = self.dff if dff == 'best' else self.dff_pipe
        if which == 'sig': m = self.sig
        elif which == 'noise': m = self.noise
        elif which == 'all': m = np.ones(self.nroi, bool)
        else: m = which
        if m.sum() == 0: return None
        return np.nanmean(D[m], axis=0)


def xcorr_fft(x, y, maxlag):
    """circular cross-correlation r(lag) of two z-scored NaN-filled series: lag > 0 = y later than x (y[t+lag] ~ x[t]); returns lags (bins), r(lag) for |lag| <= maxlag and the full circular r"""
    ok = np.isfinite(x) & np.isfinite(y)
    def z(v):
        v = np.where(np.isfinite(v), v, np.nan); m = np.nanmean(v); s = np.nanstd(v)
        return np.where(np.isfinite(v), (v - m) / s, 0.0)
    xz, yz = z(np.where(ok, x, np.nan)), z(np.where(ok, y, np.nan))
    n = len(xz); fx = np.fft.rfft(xz); fy = np.fft.rfft(yz)
    full = np.fft.irfft(np.conj(fx) * fy, n) / max(ok.sum(), 1)      # full[k] = sum_t x[t] y[t+k] / N
    lags = np.arange(-maxlag, maxlag + 1)
    return lags, full[lags % n], full


def peak_stats(full, maxlag_bins, min_sep_bins):
    """observed peak |r| (and signed r, lag) within +-maxlag_bins, and the null distribution of the same statistic: windows of the same width at lags >= min_sep_bins away"""
    n = len(full); lags = np.arange(-maxlag_bins, maxlag_bins + 1); r = full[lags % n]
    k = int(np.argmax(np.abs(r))); obs = (lags[k], r[k])
    w = 2 * maxlag_bins + 1; starts = np.arange(min_sep_bins + maxlag_bins, n - min_sep_bins - maxlag_bins - w, max(1, w // 4))
    ext = np.concatenate([full, full[:w]])
    nul = np.array([np.abs(ext[s:s + w]).max() for s in starts]) if len(starts) else np.array([np.nan])
    return obs, nul
