"""Event-triggered averages, state tuning curves, spectral coherence and saccade-triggered averages of the axon population mean (Figures 2 and 3).
usage: python event_triggered.py [session ...]   -> cache/events_<session>.pkl (all 20 aggregated sessions by default)

Inputs (no ROI traces needed): cached state-analysis results (population mean P = mean best dF/F of the signal-like ROIs, behaviour on the 5 Hz grid) + xs export (stimulus onsets, eye position).
Events (hysteresis crossings, stimulus-free): bout of face motion / whisking, locomotion onset (ball, only where usable), pupil dilation, saccade.
Null for every event type: the SAME onset times circularly shifted by a random common offset (>= 60 s), 500 draws -> pointwise 95% band and a p value for the response metric.
Response metric R = mean(P, 0..+5 s) - mean(P, -6..-1 s) in dF/F (and in session z).
"""
import os, sys, pickle, numpy as np
from fig_common import *
import xs_common as xc
from xs_common import DT

PRE, POST = 10.0, 20.0
NPRE, NPOST = int(PRE / DT), int(POST / DT)
LAGS = np.arange(-NPRE, NPOST + 1) * DT
N_SHIFT = 500
BEH_LAG_S = {'face': 0.0, 'whisk': 0.0, 'loco': 0.0, 'pupil': 1.2}      # axon(t) vs behaviour(t + lag): the lags found in analysis 1


def zs(v): v = np.asarray(v, float); return (v - np.nanmean(v)) / np.nanstd(v)


def hyst_onsets(x, lo, mid, hi, quiet_bins, rise_bins, forbid=None):
    """upward crossings of `mid` that follow >= quiet_bins below `mid` (median < lo) and are followed by a rise above `hi` within rise_bins; x NaN-free where used"""
    n = len(x); out = []; i = quiet_bins
    while i < n - rise_bins - 1:
        if np.isfinite(x[i]) and x[i - 1] <= mid < x[i]:
            pre = x[i - quiet_bins:i]; post = x[i + 1:i + 1 + rise_bins]
            if np.isfinite(pre).all() and np.isfinite(post).all() and pre.max() < mid and np.median(pre) < lo and post.max() > hi and (forbid is None or not forbid[i]):
                out.append(i); i += quiet_bins; continue
        i += 1
    return np.array(out, int)


def stim_forbid(t, on_t, pre_s=2.0, post_s=6.0):
    """bins within [onset - pre, onset + post] of ANY stimulus / trial onset: not allowed as event onsets"""
    f = np.zeros(len(t), bool)
    if on_t is None or len(on_t) == 0: return f
    for o in on_t:
        a = np.searchsorted(t, o - pre_s); b = np.searchsorted(t, o + post_s); f[max(a, 0):min(b, len(t))] = True
    return f


def windows(x, idx):
    """event-aligned windows of x (n_events x (NPRE+NPOST+1)); events with incomplete windows dropped (NaN rows removed later)"""
    W = np.full((len(idx), NPRE + NPOST + 1), np.nan)
    for k, i in enumerate(idx):
        if i - NPRE >= 0 and i + NPOST < len(x): W[k] = x[i - NPRE:i + NPOST + 1]
    return W


def base(W):
    b = np.nanmean(W[:, int((PRE - 8) / DT):int((PRE - 3) / DT)], 1, keepdims=True); return W - b


def resp_metric(W, post=(0, 5), pre=(-6, -1)):
    """per event: mean(post window) - mean(pre window), windows in s relative to onset; W rows = events (NOT baseline-subtracted)"""
    a = np.nanmean(W[:, NPRE + int(post[0] / DT):NPRE + int(post[1] / DT) + 1], 1); b = np.nanmean(W[:, NPRE + int(pre[0] / DT):NPRE + int(pre[1] / DT) + 1], 1)
    return a - b


def trig(sig_dict, idx, P, nb, rng, n_shift=N_SHIFT, post=(0, 5), pre=(-6, -1)):
    """sig_dict: name -> series on the grid. returns triggered means/SEM of everything, observed metric and the circular-shift null for P"""
    out = dict(n=0)
    WP = windows(P, idx); ok = np.isfinite(WP).all(1)
    if ok.sum() < 3: return out
    idx = idx[ok]; out['n'] = int(len(idx)); out['onset_idx'] = idx
    Wp = WP[ok]; sd = np.nanstd(P)
    out['lags'] = LAGS
    out['P'] = (np.nanmean(base(Wp), 0), np.nanstd(base(Wp), 0) / np.sqrt(len(Wp)))
    out['Pz'] = (out['P'][0] / sd, out['P'][1] / sd)
    for k, s in sig_dict.items():
        W = windows(s, idx); out[k] = (np.nanmean(base(W), 0), np.nanstd(base(W), 0) / np.sqrt(max(1, np.isfinite(W).any(1).sum()))) if np.isfinite(W).any() else None
    R = float(np.nanmean(resp_metric(Wp, post, pre))); out['R'] = R; out['Rz'] = R / sd; out['R_events'] = resp_metric(Wp, post, pre); out['win'] = (post, pre); out['sd'] = float(sd)
    # circular-shift null (same onsets shifted by one random common offset)
    nul_tr = []; nul_R = []
    for _ in range(n_shift):
        sh = rng.integers(int(60 / DT), nb - int(60 / DT)); jj = (idx + sh) % nb; Wn = windows(P, jj); okn = np.isfinite(Wn).all(1)
        if okn.sum() < max(3, 0.6 * len(idx)): continue
        nul_tr.append(np.nanmean(base(Wn[okn]), 0)); nul_R.append(float(np.nanmean(resp_metric(Wn[okn], post, pre))))
    nul_tr = np.array(nul_tr); nul_R = np.array(nul_R)
    out['null_band'] = (np.percentile(nul_tr, 2.5, 0), np.percentile(nul_tr, 97.5, 0)); out['null_R'] = nul_R; out['null_traces_z'] = (nul_tr / sd).astype(np.float16)
    out['p'] = float((1 + (np.abs(nul_R) >= abs(R)).sum()) / (1 + len(nul_R)))
    return out


def tuning(P, b, lag_s, nbin=10):
    """mean P per decile of behaviour b at lag (axon(t) vs b(t + lag)); bin centres = percentile rank"""
    k = int(round(lag_s / DT)); bb = np.full_like(b, np.nan)
    if k >= 0: bb[:len(b) - k] = b[k:]
    else: bb[-k:] = b[:len(b) + k]
    ok = np.isfinite(P) & np.isfinite(bb)
    if ok.sum() < 500: return None
    edges = np.percentile(bb[ok], np.linspace(0, 100, nbin + 1)); which = np.clip(np.searchsorted(edges, bb, side='right') - 1, 0, nbin - 1)
    m = np.array([np.nanmean(P[ok & (which == j)]) for j in range(nbin)]); sd = np.nanstd(P[ok])
    return dict(x=(np.arange(nbin) + 0.5) / nbin * 100, y=m, yz=(m - np.nanmean(P[ok])) / sd, rng=float(m[-1] - m[0]), rngz=float((m[-1] - m[0]) / sd))


def coh_spectrum(P, b, seg_s=120.0, hop_s=60.0, rng=None, shift=False):
    """magnitude-squared coherence of two grid series (NaN-free Hann segments, linear detrend). returns f (Hz), coherence, n_segments"""
    n = int(seg_s / DT); hop = int(hop_s / DT); win = np.hanning(n); f = np.fft.rfftfreq(n, DT)
    if shift: b = np.roll(b, int(len(b) / 2 + (rng.integers(-200, 200) if rng is not None else 0)))
    Sxx = 0; Syy = 0; Sxy = 0; ns = 0
    for a in range(0, len(P) - n + 1, hop):
        x = P[a:a + n]; y = b[a:a + n]
        if np.isnan(x).mean() > 0.02 or np.isnan(y).mean() > 0.02: continue
        x = np.where(np.isfinite(x), x, np.nanmean(x)); y = np.where(np.isfinite(y), y, np.nanmean(y)); tt = np.arange(n)
        x = (x - np.polyval(np.polyfit(tt, x, 1), tt)) * win; y = (y - np.polyval(np.polyfit(tt, y, 1), tt)) * win
        X = np.fft.rfft(x); Y = np.fft.rfft(y); Sxx = Sxx + np.abs(X) ** 2; Syy = Syy + np.abs(Y) ** 2; Sxy = Sxy + X * np.conj(Y); ns += 1
    if ns < 3: return f, None, ns
    return f, np.abs(Sxy) ** 2 / (Sxx * Syy), ns


def detect_saccades(tc, az, el, bad):
    """saccade onsets from the 50 Hz eye position: displacement over 100 ms of the median-filtered (5 samples) position above max(4 deg, 8 x MAD of that displacement), blink / bad-pupil frames excluded; refractory 0.5 s"""
    from scipy.ndimage import median_filter
    fs = 1 / np.median(np.diff(tc)); a = median_filter(np.where(np.isfinite(az), az, np.nanmedian(az)), 5); e = median_filter(np.where(np.isfinite(el), el, np.nanmedian(el)), 5); k = max(1, int(round(0.1 * fs)))
    d = np.hypot(a[k:] - a[:-k], e[k:] - e[:-k]); d = np.concatenate([d, np.zeros(k)]); mad = 1.4826 * np.median(np.abs(d - np.median(d))); thr = max(4.0, np.median(d) + 8 * mad)
    badm = np.convolve(bad.astype(float), np.ones(int(0.5 * fs) | 1), 'same') > 0
    cand = np.where((d[1:] >= thr) & (d[:-1] < thr) & ~badm[1:])[0] + 1; keep = []; last = -1e9
    for c in cand:
        if tc[c] - last >= 0.5: keep.append(c); last = tc[c]
    return tc[np.array(keep, int)] if keep else np.array([]), thr


def analyse(name, st):
    r = st[name]; out = dict(name=name, usable=bool(r.get('usable', False)))
    if not out['usable']: return out
    t = np.asarray(r['t'], float); P = np.asarray(r['P'], float); bg = {k: np.asarray(v, float) for k, v in r['beh_grid'].items()}; nb = len(t)
    X = xc.load_xs(name); on_t = np.asarray(X['on_t'], float) if np.size(X.get('on_t', [])) else np.array([])
    forbid = stim_forbid(t, on_t); out['n_stim_forbid_frac'] = float(forbid.mean())
    import zlib; rng = np.random.default_rng(zlib.crc32(name.encode()))
    beh = {}
    if 'face' in bg: beh['face'] = zs(bg['face'])
    if 'whisk' in bg: beh['whisk'] = zs(bg['whisk'])
    if 'pupil' in bg: beh['pupil'] = zs(bg['pupil'])
    loco_ok = bool(r.get('loco_ok', False)) and 'loco' in bg
    if loco_ok: beh['loco'] = zs(bg['loco'])
    show = dict(beh)                                                  # traces shown in every triggered average
    ev = {}
    # --- face-motion bout and whisking bout onsets
    for key in ('face', 'whisk'):
        if key not in beh: continue
        s = xc.nsmooth(beh[key], 1.0 / DT); fin = s[np.isfinite(s)]
        if len(fin) < 500: continue
        lo, mid, hi = np.percentile(fin, [35, 50, 65]); idx = hyst_onsets(s, lo, mid, hi, int(5 / DT), int(3 / DT), forbid)
        ev[key + '_bout'] = trig(show, idx, P, nb, rng)
    # --- pupil dilation onsets (slow): crossing of the median after >= 8 s below the 35th percentile, rising above the 65th within 10 s
    if 'pupil' in beh:
        s = xc.nsmooth(beh['pupil'], 2.0 / DT); fin = s[np.isfinite(s)]
        if len(fin) > 500:
            lo, mid, hi = np.percentile(fin, [35, 50, 65]); idx = hyst_onsets(s, lo, mid, hi, int(8 / DT), int(10 / DT), forbid)
            ev['pupil_dilation'] = trig(show, idx, P, nb, rng)
    # --- locomotion onsets (ball): crossing 30 mm/s after >= 5 s below 30 (median < 10), above 50 within 3 s
    if loco_ok:
        sp = xc.nsmooth(np.expm1(bg['loco']), 1.0 / DT); idx = hyst_onsets(sp, 10, 30, 50, int(5 / DT), int(3 / DT), forbid)
        ev['loco_onset'] = trig(show, idx, P, nb, rng)
    # --- saccades
    try:
        tc = np.asarray(X['tcam'], float); ae = np.asarray(X['PupilAzElRaw'], float); bad = np.asarray(X['PupilBad'], float) > 0
        if ae.ndim == 2 and np.isfinite(ae).mean() > 0.5:
            ts, thr = detect_saccades(tc, ae[0], ae[1], bad); idx = np.searchsorted(t, ts); idx = idx[(idx > 0) & (idx < nb)]; idx = idx[~forbid[idx]]
            ev['saccade'] = trig(show, idx, P, nb, rng, post=(0, 2), pre=(-3, -1)); ev['saccade']['rate_per_min'] = float(len(ts) / ((tc[-1] - tc[0]) / 60)); ev['saccade']['thr_deg'] = float(thr)
    except Exception as e:
        out['saccade_error'] = str(e)
    out['events'] = ev
    # --- tuning curves
    out['tuning'] = {}
    for key in ('pupil', 'face', 'whisk', 'loco'):
        if key in beh:
            b = beh[key]; d = tuning(P, b, BEH_LAG_S[key])
            if d: out['tuning'][key] = d
    # --- spectral coherence (log bands)
    out['coh'] = {}; edges = np.geomspace(0.01, 2.0, 13)
    for key in ('face', 'whisk', 'pupil', 'loco'):
        if key not in beh: continue
        f, c, ns = coh_spectrum(P, beh[key])
        if c is None: continue
        _, cn, _ = coh_spectrum(P, beh[key], rng=rng, shift=True)
        band = lambda cc: np.array([np.nanmean(cc[(f >= edges[i]) & (f < edges[i + 1])]) if ((f >= edges[i]) & (f < edges[i + 1])).any() else np.nan for i in range(len(edges) - 1)])
        out['coh'][key] = dict(fc=np.sqrt(edges[:-1] * edges[1:]), c=band(c), c_null=band(cn) if cn is not None else None, n_seg=ns)
    return out


def load_events(name):
    p = os.path.join(CACHE, 'events_%s.pkl' % name)
    return pickle.load(open(p, 'rb')) if os.path.exists(p) else None


if __name__ == '__main__':
    st = load_state(); names = sys.argv[1:] or ALL
    for nm in names:
        r = analyse(nm, st); pickle.dump(r, open(os.path.join(CACHE, 'events_%s.pkl' % nm), 'wb'))
        ev = r.get('events', {}); print(short(nm), {k: (v['n'], round(v.get('R', np.nan), 3), round(v.get('p', np.nan), 3)) for k, v in ev.items() if v.get('n')}, 'tuning', list(r.get('tuning', {})), flush=True)
