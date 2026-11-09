"""Does the matching pipeline (matched_behaviour.py) do what it claims?  Synthetic axons on REAL behaviour and REAL event times.

usage: python matched_validate.py [session ...]      -> OUT/validation.pkl (+ printed table)

Axon surrogate = (what the real behaviour predicts) + (phase-randomised residual, i.e. same spectrum / SD but no event locking) [+ injected event response].
The behaviour-driven part is taken from a lagged ridge fit of the session's real signal-like population on its real behaviours (face, whisking, pupil, ball), so the surrogate has
the REAL confound: events of the real design coincide with real behaviour changes (startle / light reflex / licking), and the naive event-minus-control difference inherits them.

Scenarios (per replicate a new residual surrogate; n_rep replicates; the SAME shift null as in the real analysis gives the p-value):
  lin     axon = linear function of behaviour                                     truth: 0 effect.   Naive biased, regression / matching unbiased.
  nonlin  axon = expansive (|u|^1.8) function of the behaviour prediction         truth: 0 effect.   Regression (linear) biased, matching unbiased.
  noisy   as `lin`, but the matcher sees behaviour + measurement noise (SD 0.5 / 1.0 in z units)   truth: 0.  Residual confounding: how far does a noisy behavioural proxy break matching.
  resp    as `lin` + injected event response of known size (0.25 / 0.5 / 1.0 x SD of the session's window-mean noise)  truth: injected.  Recovery and power.
  resp_nl as `nonlin` + injected response.
"""
import os, sys, time, pickle, zlib, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matched_behaviour as mb
from matched_behaviour import Rec, get_contrasts, run_contrast, _cum, DT
from stim_response import fit_ridge, predict
from xs_common import nsmooth

SIM_SESSIONS = {'AnimalRB24_261030_1744': ['cue', 'post'], 'AnimalRB24_261023_1538': ['flash', 'loom_expand'], 'AnimalRB24_261021_1630': ['flash', 'loom_expand'], 'AnimalRB23_261027_1423': ['cue', 'post']}
LAGS = (np.arange(-1.0, 3.01, 0.4) / DT).round().astype(int)
N_REP, N_NULL = int(os.environ.get('MV_REP', 60)), int(os.environ.get('MV_NULL', 100))


def phase_surrogate(x, rng):
    x = np.where(np.isfinite(x), x, 0.0); n = len(x); f = np.fft.rfft(x - x.mean()); ph = rng.uniform(0, 2 * np.pi, len(f)); ph[0] = 0
    y = np.fft.irfft(np.abs(f) * np.exp(1j * ph), n); return y * (x.std() / (y.std() + 1e-12))


def behaviour_prediction(R):
    yhat, res, _ = R.surrogate_model()
    return yhat, res


def kernel(spec, tau0=1.5):
    """event response k(tau), tau >= 0 on the grid, scaled so that its mean over the response window `spec` is 1"""
    tt = np.arange(0, 20, DT); k = tt * np.exp(-tt / tau0); a, b = int(round(spec[0] / DT)), int(round(spec[1] / DT)); k = k / k[a:b].mean(); return k


def set_signal(R, y):
    R.P = y; R.cumP = _cum(y); R._cache = {}


def set_noise(R, sigma, rng, orig):
    """matcher sees behaviour + smooth measurement noise"""
    for b in R.names:
        e = nsmooth(rng.standard_normal(R.nb), 0.5 / DT); e = e / e.std() * sigma * np.nanstd(orig[b]); R.beh[b] = orig[b] + e
    R.cum = {b: _cum(v) for b, v in R.beh.items() if v is not None}; R._cache = {}


def simulate(name, keys, n_rep=N_REP, n_null=N_NULL, verbose=True, par=None, scen_only=None):
    R = Rec(name, dict(par or {}, n_sur=0)); orig = {b: v.copy() for b, v in R.beh.items() if v is not None}; P0 = R.P.copy()
    yhat, res = behaviour_prediction(R); sd_u = np.nanstd(yhat); out = dict(name=name, rows=[])
    contrasts = {c['key']: c for c in get_contrasts(R) if c['key'] in keys}
    u = (yhat - np.nanmean(yhat)) / sd_u; nl = np.sign(u) * np.abs(u) ** 1.8; nl = (nl - np.nanmean(nl)) / np.nanstd(nl) * sd_u; nl = nl + 0 * yhat
    # unit of the injected effect: SD of the window-mean of a purely behavioural trace across control-like windows (the noise level of a window mean)
    sd_win = {}
    for key, c in contrasts.items():
        F, y, ok = R.table(c['spec']); sd_win[key] = float(np.nanstd(y[ok]))
    scen = [('lin', 0.0, 0.0, 0.0), ('nonlin', 0.0, 0.0, 1.0), ('noisy05', 0.0, 0.5, 0.0), ('noisy10', 0.0, 1.0, 0.0), ('resp025', 0.25, 0.0, 0.0), ('resp050', 0.5, 0.0, 0.0), ('resp100', 1.0, 0.0, 0.0),
            ('resp050_nl', 0.5, 0.0, 1.0)]
    if scen_only: scen = [x for x in scen if x[0] in scen_only]
    for key, c in contrasts.items():
        spec = c['spec']; kern = kernel(spec); nb = R.nb
        evb = c['focal'] if c['kind'] == 'vs_pool' else c['focal']
        for scn, dfrac, noise, nonl in scen:
            t0 = time.time(); delta = dfrac * sd_win[key]; rows = []
            for rep in range(n_rep):
                rng = np.random.default_rng(zlib.crc32(('%s%s%s%d' % (name, key, scn, rep)).encode()))
                base = nl if nonl else yhat
                y = base + phase_surrogate(res, rng)
                if delta:
                    inj = np.zeros(nb)
                    for e in evb:
                        n = min(len(kern), nb - e)
                        if n > 0: inj[e:e + n] += kern[:n]
                    y = y + delta * inj
                y = np.where(np.isfinite(yhat), y, np.nan)
                for b_ in R.names: R.beh[b_] = orig[b_]
                R.cum = {b_: _cum(v) for b_, v in R.beh.items() if v is not None}
                set_signal(R, y - np.nanmean(y))
                if noise: set_noise(R, noise, rng, orig)
                r = run_contrast(R, c, n_null, rng, with_extras=False)
                if r.get('obs') is None or 'stats' not in r: continue
                o = r['obs']; s = r['stats']
                rows.append(dict(naive=o['naive'], reg=o['reg'], matched=o['matched'], adj=o['adj'], cover=o['n_matched'] / max(o['n_focal'], 1),
                                 p_naive=s['naive']['p'], p_reg=s['reg']['p'], p_matched=s['matched']['p'], p_adj=s['adj']['p'], nsd_matched=s['matched']['nsd'], nsd_naive=s['naive']['nsd']))
            R.beh = {b_: v for b_, v in orig.items()}; R.cum = {b_: _cum(v) for b_, v in R.beh.items() if v is not None}
            if not rows: continue
            M = {k: np.array([r[k] for r in rows]) for k in rows[0]}
            # expected window-mean effect of the injected response (events with a valid window; kernel mean over the window = 1 by construction, overlapping events add)
            truth = delta * float(np.mean(_expected_gain(c, kern, nb, spec)))
            row = dict(session=name, key=key, scen=scn, delta=delta, truth=truth, n=len(rows), sd_win=sd_win[key], **{k: M[k].mean() for k in ('naive', 'reg', 'matched', 'adj', 'cover')},
                       **{'sd_' + k: M[k].std() for k in ('naive', 'reg', 'matched', 'adj')}, **{'rej_' + k: float((M['p_' + k] < 0.05).mean()) for k in ('naive', 'reg', 'matched', 'adj')},
                       raw={k: M[k] for k in M})
            out['rows'].append(row)
            if verbose:
                print('%-22s %-12s %-10s truth %+.3f | mean naive %+.3f reg %+.3f matched %+.3f adj %+.3f (sd %.3f) | reject p<.05: naive %.2f reg %.2f matched %.2f adj %.2f | cover %.2f | %.0f s' % (
                    name[6:], key, scn, truth, row['naive'], row['reg'], row['matched'], row['adj'], row['sd_matched'], row['rej_naive'], row['rej_reg'], row['rej_matched'], row['rej_adj'], row['cover'], time.time() - t0), flush=True)
    return out


def _expected_gain(c, kern, nb, spec):
    """per focal anchor: the mean over the response window of the summed kernels of ALL events (own + neighbours' tails)"""
    a, b = int(round(spec[0] / DT)), int(round(spec[1] / DT)); inj = np.zeros(nb + len(kern))
    for e in c['focal']: inj[e:e + len(kern)] += kern
    cs = np.concatenate([[0], np.cumsum(inj)]); f = c['focal']; f = f[(f + b) < nb]
    return (cs[f + b] - cs[f + a]) / max(b - a, 1)


if __name__ == '__main__':
    names = sys.argv[1:] or list(SIM_SESSIONS); path = os.path.join(mb.OUT, 'validation.pkl')
    res = pickle.load(open(path, 'rb')) if os.path.exists(path) else {}
    for nm in names:
        print(nm, flush=True); res[nm] = simulate(nm, SIM_SESSIONS[nm]); pickle.dump(res, open(path, 'wb'))
    print('MV_DONE')
