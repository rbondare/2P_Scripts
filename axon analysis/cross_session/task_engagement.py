"""Cross-session analysis 4: LC-axon activity during the lick task as a question of CONTEXT / ENGAGEMENT (not of individual licks).

usage: python task_engagement.py [session ...]   (default TASK_SESSIONS; output task_results.pkl next to the figures)

Inputs per session: <XS_DIR>\\<session>_xs.mat (traces, behaviour) and <session>_task.mat (trial table + licks on the 2P clock, made by the task script of the repo:
trial start / square onset / square end / trial end, subsession spans, hit, latency, time-outs, licks before the square, pump kind).
Task: wait 8-13 s (licking in the wait = time-out, new wait) -> square 3 s (lick = hit) -> grey 8.5 s; 2-4 subsessions separated by ~100 s breaks without trials.
Population signal P = mean best dF/F of the signal-like ROIs (0.2 s grid); P_res = P minus a behaviour-only ridge model (face + whisking + locomotion if plausible, lags -1..+3 s) = what is left when the
concurrent movement (incl. the jaw movement of licking) is accounted for; behaviour = pupil, face motion, whisking, locomotion (ball only if plausible), licks.

A  Epochs of the trial cycle (wait / square / post-square, hit and miss trials separately) and the breaks: mean of P, P_res, behaviour, lick rate (z over the session) minus the mean over all epochs;
   null = the whole label set shifted circularly (>= 60 s; spacing kept).
B  Engagement over the session: per 2 min bin, lick rate (licks / min) against axon metrics (mean, SD, event rate, P_res) and behaviour; linear time trend removed from both; null = circular shift of the lick-rate series.
   Per subsession: hit rate, latency, time-outs, licks before the square, with the same axon metrics.
C  Pre-trial state and outcome: P_hp / P_res / pupil / face / whisking in [-6, -1] s before the square vs hit, log latency of the first lick, licks before the square (Spearman); null = outcome series shifted by >= 10 trials.
D  Hit vs miss: square response (0-3 s minus -3..-0.5 s) and the pre-square level of P_hp / P_res, difference hit - miss; null = hit labels permuted inside subsessions.
"""
import sys, os, pickle, time, numpy as np, scipy.io as sio
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import xs_common as xc
from xs_common import DT
from stim_response import fit_ridge, predict
from scipy.stats import spearmanr, rankdata
from scipy.signal import find_peaks

FIG = r'c:\Users\rbondarenko\projects\2P_Scripts\axon analysis\plots\cross_session'
TASK_SESSIONS = ['AnimalRB24_261029_1752', 'AnimalRB24_261030_1744', 'AnimalRB24_261105_1423', 'AnimalRB24_261105_1520', 'AnimalRB24_261106_1317', 'AnimalRB23_261029_1704']
FLAGGED = {'AnimalRB23_261029_1704': 'behaviour channels inconsistent (face-whisk r 0.13), lick clock from face motion (+-0.2 s), ball truncated: shown, not pooled'}
HP_S = 30.0; N_SH = 500; BIN_S = 60.0
rng = np.random.default_rng(11)
EPOCHS = ['wait', 'square', 'post', 'break']


def hp(x): return x - xc.nsmooth(x, HP_S / DT)
def zs(v): return (v - np.nanmean(v)) / np.nanstd(v)


def load_task(name):
    d = sio.loadmat(os.path.join(xc.XS_DIR, name + '_task.mat'), squeeze_me=True)
    T = {k: v for k, v in d.items() if not k.startswith('__')}
    for k, v in list(T.items()):
        if isinstance(v, np.ndarray) and v.ndim == 0: T[k] = v.item()
        elif isinstance(v, np.ndarray) and v.dtype != object and k not in ('session',): T[k] = np.atleast_1d(v)
    for k in ('sub_t0', 'sub_t1', 'lick_t', 'on_t', 'tr_t0', 'st_end', 'tr_t1', 'hit', 'lat1', 'n_to', 'n_pre', 'sub_of', 'rew_t', 'rew_kind', 'lvl_t', 'bout_on', 'wait_s'): T[k] = np.atleast_1d(np.asarray(T[k], float))
    return T


def window_mean(x, e0, t_a, t_b):
    i1 = np.floor((t_a - e0) / DT).astype(int); i2 = np.floor((t_b - e0) / DT).astype(int) + 1
    c = np.concatenate([[0], np.cumsum(np.where(np.isfinite(x), x, 0.0))]); n = np.concatenate([[0], np.cumsum(np.isfinite(x).astype(float))])
    ok = (i1 >= 0) & (i2 <= len(x)) & (i2 > i1); i1c = np.clip(i1, 0, len(x)); i2c = np.clip(i2, 0, len(x))
    with np.errstate(invalid='ignore', divide='ignore'): m = (c[i2c] - c[i1c]) / (n[i2c] - n[i1c])
    m[~ok] = np.nan; return m


def spear(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 8 or np.nanstd(a[ok]) == 0 or np.nanstd(b[ok]) == 0: return np.nan
    return float(spearmanr(a[ok], b[ok])[0])


def detrend(x):
    ok = np.isfinite(x); t = np.arange(len(x), dtype=float); y = x.copy()
    if ok.sum() > 3: p = np.polyfit(t[ok], x[ok], 1); y = x - np.polyval(p, t)
    return y


def analyse(name):
    S = xc.Session(name, with_rejected=False, verbose=True); T = load_task(name); e0 = S.edges[0]; nb = S.nb
    out = dict(name=name, flag=FLAGGED.get(name), n_sig=int(S.sig.sum()), dur_min=S.dur / 60, loco_ok=S.loco_ok, nTr=int(T['nTr']), rdel=float(T['rdel']), levels=sorted(set(np.unique(T['lvl_t']).tolist())))
    if S.sig.sum() < 3: out['why'] = 'fewer than 3 signal-like ROIs'; return out
    P = S.population('sig'); Php = hp(P)
    behs = {b: zs(S.beh[b]) for b in ('pupil', 'face', 'whisk', 'loco') if b in S.beh and (b != 'loco' or S.loco_ok) and np.isfinite(S.beh[b]).mean() > 0.3}
    out['behaviours'] = list(behs)
    # lick rate on the grid
    lick = np.histogram(T['lick_t'], bins=S.edges)[0].astype(float) / DT; lick_s = xc.nsmooth(lick, 1.0 / DT)
    # behaviour-only model (whole session, in-sample ridge): P_res and its high-passed version
    fam = [b for b in ('face', 'whisk', 'loco') if b in behs]; lags = (np.arange(-1.0, 3.01, 0.4) / DT).round().astype(int)
    Pres = Pres_hp = None; out['mv_train_r'] = np.nan
    if fam:
        Xd = np.concatenate([np.stack([np.roll(hp(behs[b]), -l) for l in lags], 1) for b in fam], 1)
        model = fit_ridge(Xd, Php)
        if model is not None:
            Phat = predict(model, Xd); Pres_hp = Php - Phat; ok = np.isfinite(Phat) & np.isfinite(Php); out['mv_train_r'] = float(np.corrcoef(Phat[ok], Php[ok])[0, 1])
            # level residual: regress P (not high-passed) on the unfiltered behaviours at the same lags (slow state included)
            Xl = np.concatenate([np.stack([np.roll(behs[b], -l) for l in lags], 1) for b in fam], 1); m2 = fit_ridge(Xl, P)
            Pres = P - predict(m2, Xl) if m2 is not None else None
    series = {'P': zs(P), 'Pres': zs(Pres) if Pres is not None else None, 'face': behs.get('face'), 'whisk': behs.get('whisk'), 'pupil': behs.get('pupil'), 'loco': behs.get('loco'), 'lick': zs(lick_s)}
    series = {k: v for k, v in series.items() if v is not None}
    # ------------------------------------------------------------------ A epochs
    lab_names = []; masks = []
    def add(name_, spans):
        m = np.zeros(nb, bool)
        for (a, b) in spans:
            if not (np.isfinite(a) and np.isfinite(b)) or b <= a: continue
            i1 = int(np.floor((a - e0) / DT)); i2 = int(np.floor((b - e0) / DT)) + 1; m[max(i1, 0):min(i2, nb)] = True
        lab_names.append(name_); masks.append(m)
    on, se, t0, t1, hit = T['on_t'], T['st_end'], T['tr_t0'], T['tr_t1'], T['hit'] > 0
    for nm, sel in (('hit', hit), ('miss', ~hit)):
        add('wait (%s)' % nm, [(t0[i] + 1.0, on[i] - 0.5) for i in np.where(sel)[0]]); add('square (%s)' % nm, [(on[i], se[i]) for i in np.where(sel)[0]]); add('post (%s)' % nm, [(se[i] + 0.5, t1[i]) for i in np.where(sel)[0]])
    brk = [(T['sub_t1'][s] + 5, T['sub_t0'][s + 1] - 5) for s in range(len(T['sub_t0']) - 1) if T['sub_t0'][s + 1] - T['sub_t1'][s] > 40]
    if brk: add('break', brk)
    M = np.array(masks).astype(float); allm = (M.sum(0) > 0)
    ser_keys = list(series); Z = np.stack([np.where(np.isfinite(series[k]), series[k], 0.0) for k in ser_keys]); Fin = np.stack([np.isfinite(series[k]).astype(float) for k in ser_keys])
    def eff(Zr, Fr, Mr):
        num = Mr @ (Zr * Fr).T; den = np.maximum(Mr @ Fr.T, 1)
        grand = (Mr.sum(0) > 0).astype(float) @ (Zr * Fr).T / np.maximum((Mr.sum(0) > 0).astype(float) @ Fr.T, 1)
        return num / den - grand            # [epochs x series]
    obs = eff(Z, Fin, M)
    nul = np.zeros((N_SH, M.shape[0], len(ser_keys))); lo, hi = int(60 / DT), nb - int(60 / DT)
    for i in range(N_SH):
        k = int(rng.integers(lo, hi)); nul[i] = eff(Z, Fin, np.roll(M, k, axis=1))
    out['epochs'] = dict(names=lab_names, series=ser_keys, obs=obs, null_mean=nul.mean(0), null_sd=nul.std(0), z=(obs - nul.mean(0)) / np.maximum(nul.std(0), 1e-9), n_bins=M.sum(1), n_trials=[int(hit.sum()), int((~hit).sum())])
    # ------------------------------------------------------------------ B engagement over the session (bins and subsessions)
    nbin = int(np.floor(S.dur / BIN_S)); bins = []
    thr = 2 * np.nanstd(xc.bandpass(P)); pk, _ = find_peaks(np.nan_to_num(xc.nsmooth(xc.bandpass(P), 1.0 / DT)), height=thr, distance=int(3 / DT)); ev = np.zeros(nb, bool); ev[pk] = True
    Pbp = xc.bandpass(P)
    for j in range(nbin):
        a = j * int(BIN_S / DT); b = (j + 1) * int(BIN_S / DT); sl = slice(a, b); ta, tb = e0 + a * DT, e0 + b * DT
        row = dict(t0=ta, t1=tb, lick_rate=float(((T['lick_t'] >= ta) & (T['lick_t'] < tb)).sum() / (BIN_S / 60)), P_mean=float(np.nanmean(P[sl])), P_sd=float(np.nanstd(Pbp[sl])), P_rate=float(ev[sl].sum() / (BIN_S / 60)))
        if Pres is not None: row['Pres_mean'] = float(np.nanmean(Pres[sl])); row['Pres_sd'] = float(np.nanstd(Pres_hp[sl]))
        for bn, v in behs.items(): row[bn] = float(np.nanmean(v[sl]))
        tt = (T['on_t'] >= ta) & (T['on_t'] < tb); row['n_trials'] = int(tt.sum()); row['hit_rate'] = float(np.mean(hit[tt])) if tt.sum() >= 2 else np.nan
        row['in_break'] = float(np.mean([any(a_ <= ta + x <= b_ for a_, b_ in brk) for x in (0, 15, 30, 45, 59)])) if brk else 0.0
        bins.append(row)
    out['bins'] = bins
    mets = [k for k in bins[0] if k not in ('t0', 't1', 'lick_rate', 'n_trials', 'in_break', 'hit_rate')]
    keep = np.array([r['in_break'] < 0.2 for r in bins])          # break bins are left out: the screen is dark / different there (pupil dilates) and nobody licks
    ll = np.array([r['lick_rate'] for r in bins])[keep]; res_b = {}
    for k in mets:
        v = np.array([r[k] for r in bins])[keep]; a_, b_ = detrend(ll), detrend(v); rho = spear(a_, b_)
        nl = np.array([spear(np.roll(a_, sh), b_) for sh in range(3, len(ll) - 2)])         # exhaustive circular shifts (the only distinct ones), p resolution 1/(n-4)
        res_b[k] = dict(rho=rho, null_sd=float(np.nanstd(nl)), null_mean=float(np.nanmean(nl)), p=float((1 + np.nansum(np.abs(nl - np.nanmean(nl)) >= abs(rho - np.nanmean(nl)))) / (1 + np.isfinite(nl).sum())), n=len(ll), n_shifts=int(np.isfinite(nl).sum()))
    out['bin_corr'] = res_b
    sub = []
    for s in range(int(T['subsessions'])):
        m = T['sub_of'] == s + 1
        if m.sum() == 0: continue
        a, b = T['sub_t0'][s], T['sub_t1'][s]; sl = slice(max(int((a - e0) / DT), 0), min(int((b - e0) / DT), nb))
        row = dict(sub=s + 1, trials=int(m.sum()), minutes=(b - a) / 60, hit=float(hit[m].mean()), lat=float(np.nanmedian(T['lat1'][m])), timeouts=float(T['n_to'][m].mean()), licks_pre=float(T['n_pre'][m].mean()),
                   lick_rate=float(((T['lick_t'] >= a) & (T['lick_t'] < b)).sum() / ((b - a) / 60)), P_sd=float(np.nanstd(Pbp[sl])), P_mean=float(np.nanmean(P[sl])), P_rate=float(ev[sl].sum() / ((b - a) / 60)))
        if Pres_hp is not None: row['Pres_sd'] = float(np.nanstd(Pres_hp[sl]))
        for bn, v in behs.items(): row[bn] = float(np.nanmean(v[sl]))
        sub.append(row)
    out['subs'] = sub
    # ------------------------------------------------------------------ C pre-trial state and outcome / D hit vs miss
    pre = (-6.0, -1.0)
    X_pre = {'P_hp': window_mean(Php, e0, on + pre[0], on + pre[1])}
    if Pres_hp is not None: X_pre['P_res'] = window_mean(Pres_hp, e0, on + pre[0], on + pre[1])
    for bn in ('pupil', 'face', 'whisk', 'loco'):
        if bn in behs: X_pre[bn] = window_mean(hp(behs[bn]), e0, on + pre[0], on + pre[1])
    y_hit = hit.astype(float); y_lat = np.where(np.isfinite(T['lat1']) & hit, np.log(np.maximum(T['lat1'], 0.05)), np.nan); y_pre = np.log1p(T['n_pre'])
    n = len(on); res_c = {}
    for ynm, y in (('hit', y_hit), ('log latency', y_lat), ('licks before square', y_pre)):
        res_c[ynm] = {}
        for xn, x in X_pre.items():
            if ynm == 'licks before square' and xn in ('face', 'whisk'): pass
            rho = spear(x, y); nl = np.array([spear(x, np.roll(y, int(rng.integers(10, max(n - 10, 11))))) for _ in range(N_SH)])
            res_c[ynm][xn] = dict(rho=rho, null_sd=float(np.nanstd(nl)), null_mean=float(np.nanmean(nl)), p=float((1 + np.nansum(np.abs(nl - np.nanmean(nl)) >= abs(rho - np.nanmean(nl)))) / (1 + np.isfinite(nl).sum())), n=int(np.isfinite(x * y).sum()))
    out['pre_corr'] = res_c
    # D
    resp = {}
    for nm_, ser in (('P_hp', Php), ('P_res', Pres_hp)):
        if ser is None: continue
        a = window_mean(ser, e0, on + 0.0, on + 3.0) - window_mean(ser, e0, on - 3.0, on - 0.5); lvl = X_pre['P_hp'] if nm_ == 'P_hp' else X_pre['P_res']
        sub_of = T['sub_of'].astype(int)
        for what, v in (('response', a), ('pre level', lvl)):
            ok = np.isfinite(v); d_obs = np.nanmean(v[hit & ok]) - np.nanmean(v[~hit & ok]) if (hit & ok).sum() > 2 and (~hit & ok).sum() > 2 else np.nan
            nul = []
            for _ in range(N_SH):
                h2 = hit.copy()
                for s_ in np.unique(sub_of): ii = np.where(sub_of == s_)[0]; h2[ii] = rng.permutation(hit[ii])
                if (h2 & ok).sum() > 2 and (~h2 & ok).sum() > 2: nul.append(np.nanmean(v[h2 & ok]) - np.nanmean(v[~h2 & ok]))
            nul = np.array(nul); resp['%s %s' % (nm_, what)] = dict(diff=float(d_obs), null_sd=float(nul.std()), z=float((d_obs - nul.mean()) / nul.std()), p=float((1 + (np.abs(nul - nul.mean()) >= abs(d_obs - nul.mean())).sum()) / (1 + len(nul))), n_hit=int((hit & ok).sum()), n_miss=int((~hit & ok).sum()))
    out['hit_miss'] = resp
    # event-related traces hit / miss (P_hp, baseline -3..-0.5)
    rel = np.arange(-8, 12.01, DT); idx = np.floor((on[:, None] + rel[None, :] - e0) / DT).astype(int); okk = (idx >= 0) & (idx < nb)
    def trace(ser):
        C = np.full(idx.shape, np.nan); C[okk] = ser[idx[okk]]; base = np.nanmean(C[:, (rel >= -3) & (rel < -0.5)], 1, keepdims=True); return C - base
    out['trace_rel'] = rel; out['trace_hit'] = trace(Php)[hit]; out['trace_miss'] = trace(Php)[~hit]
    if Pres_hp is not None: out['trace_res_hit'] = trace(Pres_hp)[hit]; out['trace_res_miss'] = trace(Pres_hp)[~hit]
    # ---- compact timeline for the overview figure (1 min bins, whole session)
    out['timeline'] = dict(t=np.array([0.5 * (r['t0'] + r['t1']) for r in bins]))
    nbm = int(np.floor(S.dur / 60)); tl = {}
    def bin1(x): return np.array([np.nanmean(x[j * int(60 / DT):(j + 1) * int(60 / DT)]) for j in range(nbm)])
    tl['t'] = e0 + 60 * (np.arange(nbm) + 0.5); tl['P'] = bin1(P); tl['P_sd'] = np.array([np.nanstd(Pbp[j * int(60 / DT):(j + 1) * int(60 / DT)]) for j in range(nbm)])
    tl['lick_rate'] = np.array([((T['lick_t'] >= e0 + 60 * j) & (T['lick_t'] < e0 + 60 * (j + 1))).sum() for j in range(nbm)], float)
    for bn, v in behs.items(): tl[bn] = bin1(v)
    out['timeline'] = tl
    out['trials'] = dict(on=on, hit=hit.astype(int), lat=T['lat1'], n_to=T['n_to'], n_pre=T['n_pre'], sub=T['sub_of'], rew_kind=T['rew_kind'], lvl=T['lvl_t'], sub_t0=T['sub_t0'], sub_t1=T['sub_t1'])
    return out


if __name__ == '__main__':
    names = sys.argv[1:] or TASK_SESSIONS; res = {}
    for nm in names:
        t0 = time.time()
        try:
            res[nm] = analyse(nm); r = res[nm]
            print('  -> %s | trials %d | %.0f s' % (r.get('why', 'ok'), r['nTr'], time.time() - t0), flush=True)
        except Exception as e:
            import traceback; traceback.print_exc(); res[nm] = dict(name=nm, why='error: %s' % e)
        pickle.dump(res, open(os.path.join(FIG, 'task_results.pkl'), 'wb'))
    print('TASK_DONE')
