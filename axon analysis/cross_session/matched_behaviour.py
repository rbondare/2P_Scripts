"""Behaviour-matched comparison of LC-axon activity: is there anything in the axons beyond what whisking / face motion / pupil / locomotion already say?

usage:  python matched_behaviour.py [session ...]        (default: every session in SESSIONS; one pickle per session in OUT/cache, REFRESH=1 recomputes)
        env  MB_VARIANT=<name>  runs one robustness variant (see VARIANTS) into OUT/cache_<name>

IDEA.  The earlier analyses REGRESS the behaviour out of the axon trace (linear, lagged ridge).  Here the comparison is NON-parametric: every "event" window (a stimulus,
a lick bout, a reward ...) is paired with 1-k control windows of the SAME session whose behaviour is the same - same mean pupil, face motion, whisking (and ball speed where
the ball is usable) - and the axon activity of the two windows is compared.  If the axons only follow behaviour, matched pairs have the same axon activity; a difference
is activity that the measured behaviour does not explain.  No functional form, no lag model: the match works for any monotone or non-monotone, additive or interacting relation.

WINDOWS.  response window [a, b] s after an event anchor (cue 0-3 s, post 3-10 s, loom 0-4 s, flash 0-2 s, lick bout 0-4 s ...).  Behaviour features of a window = the mean of
every behaviour (z-scored, 30 s high-passed like the axon trace) over the SAME window shifted by that behaviour's own lag to the axon (cross-correlation peak of the session, within +-4 s;
face / whisking ~0, pupil ~+1 s: the pupil at t+1 s is what the axon at t 'sees').  Simulation (matched_validate.py): this setting has the smallest bias on behaviour-only surrogates.
Axon outcome = mean of the 30 s high-passed best dF/F of the signal-like ROIs (and every ROI separately) over [a, b].

CONTRAST TYPES.
  vs_pool   event windows vs control windows of the same session on a 1 s grid, excluding every window that touches [-3, +12] s around ANY stimulus onset (or a lick bout for lick contrasts)
  vs_class  event class A vs class B (expanding vs receding loom, rewarded vs spontaneous lick bout, hit vs miss) - the control windows are the other class

MATCHING.  Optimal (Hungarian) 1:k matching without replacement on the Mahalanobis distance of the standardised features (pool covariance), with a caliper of CAL SD on EVERY
feature separately (so no single behaviour can differ by more than CAL SD inside a pair); focal windows without a partner inside the caliper are dropped and counted (the effect then
refers to the matched subset: coverage is reported).  Balance = standardised mean difference of every feature before / after matching.
ESTIMATES (all in dF/F of the session's own recipe, mean over focal windows):
  naive     mean(event) - mean(all control windows)                       (what the raw average shows)
  reg       mean(event - linear behaviour prediction trained on control windows)   (regression, no matching; the old approach, within-window features)
  matched   mean over pairs of (event - mean of its matched controls)      <- the quantity of interest
  adj       matched, bias-corrected for the residual feature imbalance with the control-window regression (Abadie & Imbens)
NULL.  vs_pool: all design events (focal and blockers) circularly shifted together by 60 s - T-60 s relative to the data and the WHOLE pipeline re-run (so autocorrelation, block
structure and residual matching bias are in the null); vs_class: the class labels are circularly shifted along the time-ordered event sequence (run structure kept).
Session = unit; pooled across sessions in matched_figures.py (animal-clustered bootstrap, unflagged sessions only, flagged drawn hollow, as in paper_figures).
"""
import os, sys, time, pickle, zlib, warnings, numpy as np
warnings.filterwarnings('ignore', category=RuntimeWarning)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import xs_common as xc
from xs_common import DT
from stim_response import hp as _hp30
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
if os.environ.get('MB_XS'): xc.XS_DIR = os.environ['MB_XS']      # e.g. C:\Users\rbondarenko\2P_cache\xs_v2: xs + task exports with the 2026-11-10 whisking index / JawME (default: ...\xs, 2026-11-08 export)

OUT = os.environ.get('MB_OUT', r'c:\Users\rbondarenko\projects\2P_Scripts\axon analysis\plots\matched_behaviour')
VARIANT = os.environ.get('MB_VARIANT', '')
CACHE = os.path.join(OUT, 'cache' + ('_' + VARIANT if VARIANT else ''))
os.makedirs(CACHE, exist_ok=True)

SESSIONS = ['AnimalRB24_261021_1630', 'AnimalRB24_261023_1538', 'AnimalRB24_261027_1456', 'AnimalRB24_261027_1534', 'AnimalRB24_261029_1752', 'AnimalRB24_261030_1744',
            'AnimalRB24_261104_1651', 'AnimalRB24_261105_1423', 'AnimalRB24_261105_1520', 'AnimalRB23_261023_1642', 'AnimalRB23_261027_1423', 'AnimalRB22_261027_1138',
            'AnimalRB22_261023_1427',
            'AnimalRB24_261023_1527', 'AnimalRB24_261104_1640', 'AnimalRB24_261106_1317', 'AnimalRB24_261106_1529', 'AnimalRB24_261106_1604', 'AnimalRB23_261027_1306', 'AnimalRB23_261029_1704', 'AnimalRB24_261109_1506']
TASK_SESSIONS = {'AnimalRB24_261029_1752', 'AnimalRB24_261030_1744', 'AnimalRB24_261105_1423', 'AnimalRB24_261105_1520', 'AnimalRB24_261106_1317', 'AnimalRB23_261029_1704', 'AnimalRB24_261109_1506'}

# ---- default parameters (VARIANTS change one at a time)
DEFAULT = dict(hp_s=30.0, lead=0.0, lag=0.0, lagfit=True, n_sur=200, pupil='pre', pupil_pre=(5.0, 1.0), hist=10.0, mincov=0.8, step=1.0, cal=0.5, k=2, fset='mean', beh=('face', 'whisk', 'loco'), n_null=300, lam=0.01,
               guard=(3.0, 12.0), lick_guard=(3.0, 8.0))
VARIANTS = {
    'cal025': dict(cal=0.25), 'cal100': dict(cal=1.0), 'k1': dict(k=1), 'k4': dict(k=4),
    'fsd': dict(fset='mean+sd'), 'fhist': dict(fset='mean+hist'), 'fdyn': dict(fset='mean+dyn'), 'fdyn_cal075': dict(fset='mean+dyn', cal=0.75),
    'pupilin': dict(pupil='in'), 'pupilnone': dict(pupil='none'), 'facewhisk': dict(beh=('face', 'whisk'), pupil='none'), 'noloco': dict(beh=('face', 'whisk')),
    'nohp': dict(hp_s=0.0), 'hp10': dict(hp_s=10.0), 'nolagfit': dict(lagfit=False), 'step5': dict(step=5.0), 'ext1': dict(lead=1.0, lag=2.0), 'ext05': dict(lead=0.5, lag=0.5),
}
BIG = 1e9
INCLUDE_LICK_POOL = bool(os.environ.get('MB_LICK_POOL'))
rng_g = np.random.default_rng(5)

# window per contrast (s relative to the anchor)
SPEC = {'cue': (0.0, 3.0), 'post': (3.0, 10.0), 'flash': (0.0, 2.0), 'bar': (0.5, 4.0), 'loom_expand': (0.0, 4.0), 'loom_recede': (0.0, 4.0), 'expand_vs_recede': (0.0, 4.0),
        'rew_vs_spont': (0.0, 4.0), 'rew_lick': (0.0, 4.0), 'spont_lick': (0.0, 4.0), 'hit_vs_miss': (0.0, 6.0)}
LABEL = {'cue': 'square cue (0-3 s)', 'post': 'square post (3-10 s)', 'flash': 'flash', 'bar': 'moving bar', 'loom_expand': 'expanding loom', 'loom_recede': 'receding loom',
         'expand_vs_recede': 'expanding - receding', 'rew_vs_spont': 'rewarded - spontaneous lick bout', 'rew_lick': 'rewarded lick bout', 'spont_lick': 'spontaneous lick bout',
         'hit_vs_miss': 'hit - miss trial'}


def zs(v): return (v - np.nanmean(v)) / np.nanstd(v)


# =====================================================================================================  prepared session
def _cum(x):
    ok = np.isfinite(x); y = np.where(ok, x, 0.0)
    return np.concatenate([[0], np.cumsum(y)]), np.concatenate([[0], np.cumsum(ok.astype(float))]), np.concatenate([[0], np.cumsum(y * y)])


def _wstat(cum, i1, i2, mincov, sd=False):
    cs, cn, cs2 = cum; n = len(cs) - 1; i1 = np.asarray(i1); i2 = np.asarray(i2)
    ok = (i1 >= 0) & (i2 <= n) & (i2 > i1); a = np.clip(i1, 0, n); b = np.clip(i2, 0, n)
    cnt = cn[b] - cn[a]; ok &= cnt >= mincov * np.maximum(i2 - i1, 1)
    with np.errstate(invalid='ignore', divide='ignore'):
        m = (cs[b] - cs[a]) / cnt
        if sd: m = np.sqrt(np.maximum((cs2[b] - cs2[a]) / cnt - m * m, 0))
    m = np.asarray(m, float); m[~ok] = np.nan; return m


class Rec:
    """one session prepared for matching (series on the 5 Hz grid, cumulative sums for fast window statistics)"""
    BEH_ALL = ('pupil', 'face', 'whisk', 'loco')      # behaviours taken from the session (task_deep_dive extends it by 'jaw' and 'lick' via a subclass)

    def __init__(self, name, par=None, S=None):
        self.par = dict(DEFAULT, **(par or {})); p = self.par; self.name = name
        S = S or xc.Session(name, with_rejected=False, verbose=False); self.S = S; self.nb = S.nb; self.e0 = S.edges[0]; self.n_sig = int(S.sig.sum()); self.ok = self.n_sig >= 3
        hpf = (lambda x: x - xc.nsmooth(x, p['hp_s'] / DT)) if p['hp_s'] else (lambda x: x - np.nanmean(x))
        self.P = hpf(S.population('sig')) if self.ok else None
        self.D = np.stack([hpf(r.astype(float)) for r in S.dff[S.sig]]) if self.ok else None
        self.beh = {b: hpf(zs(S.beh[b])) for b in self.BEH_ALL if b in S.beh and (b != 'loco' or S.loco_ok) and np.isfinite(S.beh[b]).mean() > 0.3}
        self.names = [b for b in (tuple(p['beh']) + (('pupil',) if p['pupil'] == 'in' else ())) if b in self.beh]       # behaviours matched on inside the window (+ the pupil only if pupil == 'in')
        core = [self.beh[b] / np.nanstd(self.beh[b]) for b in ('pupil', 'face', 'whisk') if b in self.beh]
        self.beh['arousal'] = np.nanmean(np.stack(core), 0) if core else None
        self.cum = {b: _cum(v) for b, v in self.beh.items() if v is not None}
        self.cumP = _cum(self.P) if self.ok else None
        self.blag = {b: 0 for b in self.beh}                      # lag (bins) of every behaviour: the behaviour at t + blag is what the axon at t "sees"
        if self.ok and p['lagfit']:
            for b in self.names:
                lags, r, _ = xc.xcorr_fft(self.P, self.beh[b], int(4 / DT)); self.blag[b] = int(lags[int(np.nanargmax(np.abs(r)))])
        self.cumD = None
        self._cache = {}

    # --- features and outcome at every bin as anchor, for one response window
    def table(self, spec):
        key = (spec, self.par['fset'], tuple(self.names))
        if key in self._cache: return self._cache[key]
        p = self.par; idx = np.arange(self.nb); a, b = int(round(spec[0] / DT)), int(round(spec[1] / DT)); lead, lag = int(round(p['lead'] / DT)), int(round(p['lag'] / DT))
        cols = [_wstat(self.cum[nm], idx + a + self.blag[nm] - lead, idx + b + self.blag[nm] + lag, p['mincov']) for nm in self.names]
        if 'sd' in p['fset']: cols += [_wstat(self.cum[nm], idx + a + self.blag[nm] - lead, idx + b + self.blag[nm] + lag, p['mincov'], sd=True) for nm in self.names]
        if p['pupil'] == 'pre' and 'pupil' in self.cum:      # pupil STATE just before the window (the light reflex of visual stimuli lives inside it)
            cols.append(_wstat(self.cum['pupil'], idx + a - int(round(p['pupil_pre'][0] / DT)), idx + a - int(round(p['pupil_pre'][1] / DT)), p['mincov']))
        if 'dyn' in p['fset']:      # within-window dynamics: second half minus first half of the window, per behaviour (rising vs falling pupil / face burst early vs late)
            m = (a + b) // 2
            for nm in self.names:
                L = self.blag[nm]; x1 = _wstat(self.cum[nm], idx + a + L - lead, idx + m + L, p['mincov']); x2 = _wstat(self.cum[nm], idx + m + L, idx + b + L + lag, p['mincov']); cols.append(x2 - x1)
        if 'hist' in p['fset'] and self.beh.get('arousal') is not None:
            h = int(round(p['hist'] / DT)); cols.append(_wstat(self.cum['arousal'], idx + a - lead - h, idx + a - lead, p['mincov']))
        F = np.stack(cols, 1); y = _wstat(self.cumP, idx + a, idx + b, p['mincov'])
        ok = np.isfinite(F).all(1) & np.isfinite(y); self._cache[key] = (F, y, ok); return self._cache[key]

    def surrogate_model(self):
        """behaviour-only description of the axon trace: lagged ridge (lags -1..+3 s) trained OUTSIDE every stimulus window [-3, +12] s (so event responses are not learned),
        prediction yhat everywhere, residual `res` (phase-randomised into surrogates).  Returns (yhat, res, train_fraction)"""
        if getattr(self, '_sur', None) is not None: return self._sur
        from stim_response import fit_ridge, predict
        lags = (np.arange(-1.0, 3.01, 0.4) / DT).round().astype(int)
        Xd = np.concatenate([np.stack([np.roll(self.beh[b], -l) for l in lags], 1) for b in self.names], 1)
        tr = np.isfinite(self.P).copy(); g0, g1 = [int(round(g / DT)) for g in self.par['guard']]
        if len(getattr(self, 'allon', [])): tr &= ~_blocked(self.nb, self.allon, g0, g1, 0, 1)
        model = fit_ridge(np.where(tr[:, None], Xd, np.nan), np.where(tr, self.P, np.nan))
        if model is None: self._sur = (None, None, float(tr.mean())); return self._sur
        yhat = predict(model, Xd); self._sur = (yhat, self.P - yhat, float(tr.mean())); return self._sur

    def feat_stats(self, spec):
        """mean / SD / covariance of the features over ALL valid windows of the session (anchor on the 1 s grid): the unit of the caliper and the whitening of the distance, identical for every
        contrast of the session (a small control class, e.g. 10 receding looms, has no usable SD of its own)"""
        F, y, ok = self.table(spec); X = F[ok][::max(1, int(round(self.par['step'] / DT)))]
        mu, sd = X.mean(0), X.std(0) + 1e-9; Z = (X - mu) / sd; C = np.atleast_2d(np.cov(Z.T)) + 1e-3 * np.eye(Z.shape[1]); return mu, sd, C

    def roi_outcome(self, idx, spec):
        """mean dF/F of every signal-like ROI in the response window at anchors idx -> nroi x len(idx)"""
        if self.cumD is None: self.cumD = [_cum(r) for r in self.D]
        a, b = int(round(spec[0] / DT)), int(round(spec[1] / DT)); return np.stack([_wstat(c, idx + a, idx + b, self.par['mincov']) for c in self.cumD])


# =====================================================================================================  matching and estimates
def do_match(Zf, Zr, cal, k, C):
    nf, nr = len(Zf), len(Zr); matches = [[] for _ in range(nf)]
    if nf == 0 or nr == 0: return matches
    w, V = np.linalg.eigh(C); Wm = V @ np.diag(w ** -0.5) @ V.T
    cost = cdist(Zf @ Wm, Zr @ Wm, 'sqeuclidean'); cost[cdist(Zf, Zr, 'chebyshev') > cal] = BIG
    used = np.zeros(nr, bool)
    for _ in range(k):
        c = cost.copy(); c[:, used] = BIG
        r, j = linear_sum_assignment(c); good = c[r, j] < BIG / 2
        if not good.any(): break
        for i, jj in zip(r[good], j[good]): matches[i].append(jj); used[jj] = True
    return matches


def ridge(Z, y, lam):
    n = len(y); ym = y.mean(); w = np.linalg.solve(Z.T @ Z + lam * n * np.eye(Z.shape[1]), Z.T @ (y - ym)); return w, ym


def estimate(Ff, yf, Fr, yr, par, keep_pairs=False, fs=None):
    """Ff/Fr raw feature rows (finite), yf/yr outcomes -> dict of estimates; features standardised by the control pool"""
    nf, nr = len(yf), len(yr); res = dict(n_focal=nf, n_ref=nr)
    if nf < 3 or nr < 10: return None
    mu, sd, C = fs if fs is not None else (Fr.mean(0), Fr.std(0) + 1e-9, np.atleast_2d(np.cov(((Fr - Fr.mean(0)) / (Fr.std(0) + 1e-9)).T)) + 1e-3 * np.eye(Fr.shape[1])); Zf, Zr = (Ff - mu) / sd, (Fr - mu) / sd
    beta, ym = ridge(Zr, yr, par['lam'])
    matches = do_match(Zf, Zr, par['cal'], par['k'], C); has = np.array([len(m) > 0 for m in matches]); nm = int(has.sum())
    res.update(n_matched=nm, naive=float(yf.mean() - yr.mean()), reg=float((yf - (Zf @ beta + ym)).mean()))
    if nm < 3:
        res.update(matched=np.nan, adj=np.nan, naive_m=np.nan, reg_m=np.nan, smd_before=np.full(Zf.shape[1], np.nan), smd_after=np.full(Zf.shape[1], np.nan)); return res
    idx = np.where(has)[0]; yref = np.array([yr[matches[i]].mean() for i in idx]); zref = np.stack([Zr[matches[i]].mean(0) for i in idx])
    dpair = yf[idx] - yref; res.update(matched=float(dpair.mean()), adj=float((dpair - (Zf[idx] - zref) @ beta).mean()), naive_m=float(yf[idx].mean() - yr.mean()),
                                       reg_m=float((yf[idx] - (Zf[idx] @ beta + ym)).mean()), smd_before=Zf[idx].mean(0), smd_after=(Zf[idx] - zref).mean(0),
                                       smd_before_all=Zf.mean(0), mean_dist=float(np.abs(Zf[idx] - zref).max(1).mean()))
    if keep_pairs: res['pairs'] = (idx, [np.array(matches[i]) for i in idx], dpair)
    return res


def _blocked(nb, events, g0, g1, a, b):
    """anchors whose response window [i+a, i+b) touches [e-g0, e+g1) of any event e (bins)"""
    d = np.zeros(nb + 2, np.int32)
    lo = np.clip(events - g0 - b + 1, 0, nb + 1); hi = np.clip(events + g1 - a, 0, nb + 1)
    np.add.at(d, lo, 1); np.add.at(d, hi, -1); return np.cumsum(d)[:nb] > 0


def null_stats(null, obs):
    null = np.asarray(null, float); null = null[np.isfinite(null)]
    if len(null) < 20 or not np.isfinite(obs): return dict(z=np.nan, p=np.nan, nmean=np.nan, nsd=np.nan)
    m, s = null.mean(), null.std()
    if s < 1e-7 or len(np.unique(np.round(null, 8))) < 8: return dict(z=np.nan, p=np.nan, nmean=float(m), nsd=float(s))
    s = s + 1e-12
    return dict(z=float((obs - m) / s), p=float((1 + (np.abs(null - m) >= abs(obs - m)).sum()) / (1 + len(null))), nmean=float(m), nsd=float(s))


# =====================================================================================================  contrasts
def to_bin(R, t):
    b = np.floor((np.asarray(t, float) - R.e0) / DT).astype(int); return b[(b >= 0) & (b < R.nb)]


def get_contrasts(R):
    """list of contrast dicts for this session: key, kind, focal (bins), ref (bins, vs_class) or blockers [(bins, g0, g1)]"""
    X = R.S.X; p = R.par; name = R.name; out = []
    on_t = np.atleast_1d(np.asarray(X['on_t'], float)); on_key = np.atleast_1d(np.asarray(X['on_key'])).astype(str); on_blk = np.atleast_1d(np.asarray(X['on_blk'], int))
    binter = np.atleast_1d(np.asarray(X['blk_interrupted'], float)).astype(bool)
    keep = np.isfinite(on_t) & ~binter[np.clip(on_blk, 1, len(binter)) - 1] if len(on_t) else np.zeros(0, bool)
    on_t, on_key = on_t[keep], on_key[keep]
    T = None
    if name in TASK_SESSIONS:
        from task_engagement import load_task
        T = load_task(name); sq = np.asarray(T['on_t'], float); sq = sq[np.isfinite(sq)]
        on_t = np.concatenate([on_t[on_key != 'square'], sq]); on_key = np.concatenate([on_key[on_key != 'square'], np.full(len(sq), 'square')])
    g0, g1 = [int(round(g / DT)) for g in p['guard']]; lg0, lg1 = [int(round(g / DT)) for g in p['lick_guard']]
    allon = to_bin(R, on_t); R.allon = allon
    for key in ('square', 'flash', 'bar', 'loom_expand', 'loom_recede'):
        m = on_key == key
        if m.sum() < 8: continue
        keys = ['cue', 'post'] if key == 'square' else [key]
        for ck in keys: out.append(dict(key=ck, kind='vs_pool', spec=SPEC[ck], focal=to_bin(R, on_t[m]), blockers=[(allon, g0, g1)]))
    me, mr = on_key == 'loom_expand', on_key == 'loom_recede'
    if me.sum() >= 6 and mr.sum() >= 6:
        out.append(dict(key='expand_vs_recede', kind='vs_class', spec=SPEC['expand_vs_recede'], focal=to_bin(R, on_t[me]), ref=to_bin(R, on_t[mr]), order_t=(on_t[me], on_t[mr])))
    if T is not None:
        co = np.asarray(T['cons_on'], float); co = co[np.isfinite(co)]; sp = np.asarray(T['sp_on'], float); sp = sp[np.isfinite(sp)]
        bo = np.asarray(T['bout_on'], float); bo = bo[np.isfinite(bo)]
        if INCLUDE_LICK_POOL:     # lick bouts vs non-lick windows: licking states have (almost) no non-licking match, kept switchable only to document that
            lickblock = [(allon, g0, g1), (to_bin(R, bo), lg0, lg1)]
            if len(co) >= 8: out.append(dict(key='rew_lick', kind='vs_pool', spec=SPEC['rew_lick'], focal=to_bin(R, co), blockers=lickblock))
            if len(sp) >= 8: out.append(dict(key='spont_lick', kind='vs_pool', spec=SPEC['spont_lick'], focal=to_bin(R, sp), blockers=lickblock))
        if len(co) >= 8 and len(sp) >= 8: out.append(dict(key='rew_vs_spont', kind='vs_class', spec=SPEC['rew_vs_spont'], focal=to_bin(R, co), ref=to_bin(R, sp), order_t=(co, sp)))
        hit = np.asarray(T['hit'], float); sq_t = np.asarray(T['on_t'], float)
        if np.nansum(hit == 1) >= 8 and np.nansum(hit == 0) >= 8:
            out.append(dict(key='hit_vs_miss', kind='vs_class', spec=SPEC['hit_vs_miss'], focal=to_bin(R, sq_t[hit == 1]), ref=to_bin(R, sq_t[hit == 0]), order_t=(sq_t[hit == 1], sq_t[hit == 0])))
    return out


def run_contrast(R, c, n_null, rng, with_extras=True):
    p = R.par; spec = c['spec']; F, y, ok = R.table(spec); fs = R.feat_stats(spec); nb = R.nb; a, b = int(round(spec[0] / DT)), int(round(spec[1] / DT))
    step = max(1, int(round(p['step'] / DT)))
    if c['kind'] == 'vs_pool':
        cand = np.arange(0, nb, step)
        def blocked_for(shift):
            m = np.zeros(nb, bool)
            for ev, g0, g1 in c['blockers']: m |= _blocked(nb, (ev + shift) % nb, g0, g1, a, b)
            return m
        def one(shift, keep_pairs=False, yv=None, okv=None):
            yv = y if yv is None else yv; okv = ok if okv is None else okv
            fa = (c['focal'] + shift) % nb; fa = fa[(fa >= 0) & okv[fa]]; pool = cand[okv[cand] & ~blocked_for(shift)[cand]]
            if len(fa) == 0 or len(pool) == 0: return None, fa, pool
            return estimate(F[fa], yv[fa], F[pool], yv[pool], p, keep_pairs, fs), fa, pool
    else:
        A, B = c['focal'], c['ref']; A = A[(A >= 0) & (A < nb)]; B = B[(B >= 0) & (B < nb)]
        allb = np.concatenate([A, B]); lab = np.concatenate([np.ones(len(A), bool), np.zeros(len(B), bool)]); o = np.argsort(allb, kind='stable'); allb, lab = allb[o], lab[o]
        def one(shift, keep_pairs=False, yv=None, okv=None):
            yv = y if yv is None else yv; okv = ok if okv is None else okv
            lb = np.roll(lab, shift) if np.isscalar(shift) else shift; fa, ra = allb[lb], allb[~lb]; fa = fa[okv[fa]]; ra = ra[okv[ra]]
            if len(fa) == 0 or len(ra) == 0: return None, fa, ra
            return estimate(F[fa], yv[fa], F[ra], yv[ra], p, keep_pairs, fs), fa, ra
    obs, fa, pool = one(0, True)
    if obs is None or not np.isfinite(obs.get('matched', np.nan)):
        return dict(key=c['key'], kind=c['kind'], spec=spec, n_events=len(c['focal']), n_focal=0 if obs is None else obs['n_focal'], n_matched=0, note='too few matched windows', obs=obs)
    res = dict(key=c['key'], kind=c['kind'], spec=spec, n_events=int(len(c['focal'])), obs=obs, n_focal=obs['n_focal'], n_matched=obs['n_matched'], n_ref=obs['n_ref'])
    # null
    nul = {k: [] for k in ('naive', 'reg', 'matched', 'adj', 'n_matched')}
    if c['kind'] == 'vs_pool':
        shifts = rng.integers(int(60 / DT), nb - int(60 / DT), n_null)
    else:
        Nn = len(lab); lo = max(3, min(len(A), len(B))); shifts = rng.integers(lo, max(Nn - lo, lo + 1), n_null)
        if c.get('strata_edges') is not None:      # optional (task_deep_dive): labels permuted INSIDE strata (e.g. subsessions) instead of rolled along the sequence (balanced classes leave a roll almost no distinct shifts)
            sid = np.searchsorted(np.asarray(c['strata_edges']), allb, side='right'); shifts = []
            for _ in range(n_null):
                lp = lab.copy()
                for s_ in np.unique(sid): ii = np.where(sid == s_)[0]; lp[ii] = rng.permutation(lab[ii])
                shifts.append(lp)
    for sft in shifts:
        e, _, _ = one(sft if not np.isscalar(sft) else int(sft))
        if e is None or not np.isfinite(e.get('matched', np.nan)): continue
        for k in nul: nul[k].append(e[k] if k != 'n_matched' else e['n_matched'])
    res['null'] = {k: np.array(v, float) for k, v in nul.items()}
    res['stats'] = {k: null_stats(res['null'][k], obs[k]) for k in ('naive', 'reg', 'matched', 'adj')}
    # behaviour-only benchmark: the SAME events and the SAME pipeline on surrogate axon traces = (lagged linear behaviour prediction trained outside the stimulus windows) + (phase-randomised residual)
    yhat, resid, _ = R.surrogate_model()
    if yhat is not None and p['n_sur'] > 0:
        r0 = np.where(np.isfinite(resid), resid, 0.0); f0 = np.fft.rfft(r0 - r0.mean()); sdr = np.std(r0); ben = {k: [] for k in ('naive', 'reg', 'matched', 'adj')}
        idx_all = np.arange(nb)
        for _ in range(p['n_sur']):
            ph = rng.uniform(0, 2 * np.pi, len(f0)); ph[0] = 0; sr = np.fft.irfft(np.abs(f0) * np.exp(1j * ph), nb); sr *= sdr / (sr.std() + 1e-12)
            ys = np.where(np.isfinite(yhat), yhat + sr, np.nan); ys = ys - np.nanmean(ys); yv = _wstat(_cum(ys), idx_all + a, idx_all + b, p['mincov'])
            e, _, _ = one(0, False, yv, ok & np.isfinite(yv))
            if e is None or not np.isfinite(e.get('matched', np.nan)): continue
            for k in ben: ben[k].append(e[k])
        res['bench'] = {k: np.array(v, float) for k, v in ben.items()}
        res['stats_sim'] = {}
        for k in ben:
            v = res['bench'][k]
            if len(v) >= 20 and np.isfinite(obs[k]):
                m, sd_ = v.mean(), v.std() + 1e-12; res['stats_sim'][k] = dict(z=float((obs[k] - m) / sd_), p=float((1 + (np.abs(v - m) >= abs(obs[k] - m)).sum()) / (1 + len(v))), nmean=float(m), nsd=float(sd_))
    if with_extras:
        idx, mlist, dpair = obs['pairs']; fidx = fa[idx]
        res['pair_d'] = dpair.astype(np.float32)
        # per-ROI matched effect with the same pairs
        yr_roi = R.roi_outcome(fidx, spec); rr = [pool[m] for m in mlist]
        ref_roi = np.stack([np.nanmean(R.roi_outcome(r, spec), 1) for r in rr], 1)
        d_roi = yr_roi - ref_roi; res['roi_effect'] = np.nanmean(d_roi, 1).astype(np.float32); res['roi_frac_pos'] = float(np.mean(res['roi_effect'] > 0))
        # aligned traces
        tau = np.arange(int(-5 / DT), int(15 / DT) + 1); res['tau'] = tau * DT
        def trace(x, anchors): ii = anchors[:, None] + tau[None, :]; o = np.full(ii.shape, np.nan); v = (ii >= 0) & (ii < nb); o[v] = x[ii[v]]; return o
        series = {'P': R.P}; series.update({b_: R.beh[b_] for b_ in R.names})
        tr = {}
        for nm, x in series.items():
            tf = trace(x, fidx); tref = np.stack([np.nanmean(trace(x, r), 0) for r in rr], 0)
            tr[nm] = dict(focal=np.nanmean(tf, 0).astype(np.float32), ref=np.nanmean(tref, 0).astype(np.float32), diff_mean=np.nanmean(tf - tref, 0).astype(np.float32), diff_sem=(np.nanstd(tf - tref, 0) / np.sqrt(np.isfinite(tf - tref).sum(0))).astype(np.float32))
        res['traces'] = tr; res['focal_bins'] = fidx; res['ref_bins'] = [np.asarray(r, np.int32) for r in rr]; res['pool_bins'] = pool.astype(np.int32)
        res['n_dropped_unmatched'] = int(obs['n_focal'] - obs['n_matched'])
    return res


def run_session(name, par=None, n_null=None, verbose=True, with_extras=True):
    t0 = time.time(); R = Rec(name, par); p = R.par; n_null = n_null or p['n_null']
    out = dict(name=name, xs_dir=xc.XS_DIR, n_sig=R.n_sig, n_roi=R.S.nroi, dur_min=R.S.dur / 60, behaviours=R.names, loco_ok=R.S.loco_ok, par=p, contrasts={})
    if not R.ok: out['why'] = 'fewer than 3 signal-like ROIs'; return out
    if len(R.names) < 2: out['why'] = 'fewer than 2 behaviours'; return out
    out['sd_P'] = float(np.nanstd(R.P))
    for c in get_contrasts(R):
        rng = np.random.default_rng(zlib.crc32((name + c['key']).encode())); rr = run_contrast(R, c, n_null, rng, with_extras)
        out['contrasts'][c['key']] = rr
        if verbose and rr.get('obs'):
            o = rr['obs']; s = rr.get('stats', {}); ss = rr.get('stats_sim', {})
            print('  %-18s n %3d valid %3d matched %3d | naive %+.3f reg %+.3f matched %+.3f (shift z %.1f p %.3f | behaviour-only bench %+.3f z %.1f p %.3f) | max|SMD| %.2f -> %.2f' % (
                c['key'], rr['n_events'], o['n_focal'], o['n_matched'], o['naive'], o['reg'], o['matched'], s.get('matched', {}).get('z', np.nan), s.get('matched', {}).get('p', np.nan),
                ss.get('matched', {}).get('nmean', np.nan), ss.get('matched', {}).get('z', np.nan), ss.get('matched', {}).get('p', np.nan),
                np.nanmax(np.abs(o['smd_before'])) if np.size(o.get('smd_before', 0)) else np.nan, np.nanmax(np.abs(o['smd_after'])) if np.size(o.get('smd_after', 0)) else np.nan), flush=True)
        elif verbose: print('  %-18s n %3d: %s' % (c['key'], rr['n_events'], rr.get('note', '')), flush=True)
    if verbose: print('  %.0f s' % (time.time() - t0), flush=True)
    return out


if __name__ == '__main__':
    par = VARIANTS.get(VARIANT) if VARIANT else None
    names = sys.argv[1:] or SESSIONS
    for nm in names:
        f = os.path.join(CACHE, nm + '.pkl')
        if os.path.exists(f) and not os.environ.get('REFRESH'): print('cached', nm); continue
        print(nm, VARIANT or 'default', flush=True)
        try:
            r = run_session(nm, par, n_null=int(os.environ.get('MB_NULL', 0)) or (100 if VARIANT else None), with_extras=not VARIANT)
            pickle.dump(r, open(f, 'wb'))
        except Exception as e:
            import traceback; traceback.print_exc()
    print('MB_DONE')
