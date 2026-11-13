"""Do the factor LOADINGS replicate across two sessions of the same field of view?  python fm_cross.py [domain] [sessionA sessionB]

ROIs are matched through the fov_matching master axons (Recordings_DB rois.csv, as ag_replicate.py). The factor TIME COURSES cannot be compared (different recordings); what can be compared is
which axons load on which factor: the loading matrix of the common ROIs (n_common x k) of A and B, rotation-invariantly (principal angles between the column spaces), against a null that
permutes the correspondence between the ROIs of B. Also the loading profile of factor 1 (Pearson / Spearman) and of the matched residual factors.
"""
import os, sys, pickle, numpy as np
from scipy.stats import spearmanr
from scipy.optimize import linear_sum_assignment
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'axon_groups'))
import ag_replicate as arep, fm_fit as fm                                       # noqa: E402


def load(name, domain):
    with open(os.path.join(fm.CACHE, '%s_fm_%s.pkl' % (name, domain)), 'rb') as fh: return pickle.load(fh)


def cross(A, B, domain='band', n_perm=2000, verbose=True):
    a, b = load(A, domain), load(B, domain); mA, mB = arep.master_of(A), arep.master_of(B); revB = {}
    for r, m in mB.items(): revB.setdefault(m, []).append(r)
    ia, ib = [], []
    for k, r in enumerate(a['acc_idx']):
        r = int(r)
        if r in mA:
            cand = [x for x in revB.get(mA[r], []) if x in set(int(q) for q in b['acc_idx'])]
            if cand: ia.append(k); ib.append(int(np.where(b['acc_idx'] == cand[0])[0][0]))
    ia, ib = np.array(ia), np.array(ib); k = min(a['kstar'], b['kstar']); out = dict(n_common=len(ia), k=k, A=A, B=B)
    if len(ia) < 8 or k < 1: return out
    LA, LB = a['L'][ia, :k], b['L'][ib, :k]; cos = fm.principal_cos(LA, LB); rng = np.random.default_rng(0)
    nul = np.array([fm.principal_cos(LA, LB[rng.permutation(len(ib))]) for _ in range(n_perm)]); p = [(1 + (nul[:, j] >= cos[j]).sum()) / (1 + n_perm) for j in range(k)]
    out.update(cos=cos, null_mean=nul.mean(0), null_p95=np.quantile(nul, 0.95, axis=0), p=p)
    out['f1_pearson'] = float(np.corrcoef(LA[:, 0], LB[:, 0])[0, 1]); out['f1_spearman'] = float(spearmanr(LA[:, 0], LB[:, 0])[0])
    if k >= 2:
        C = np.abs(np.corrcoef(LA[:, 1:].T, LB[:, 1:].T)[:k - 1, k - 1:]); r_, c_ = linear_sum_assignment(-C); out['resid_factor_match_absr'] = C[r_, c_]
    if verbose:
        print('%s -> %s [%s]: %d common ROIs, k = %d | principal-angle cosines %s vs permutation mean %s / p95 %s (p %s) | factor 1 loading profile r %.2f (Spearman %.2f)%s' % (
            A[-9:], B[-9:], domain, len(ia), k, np.round(cos, 2), np.round(out['null_mean'], 2), np.round(out['null_p95'], 2), np.round(p, 3), out['f1_pearson'], out['f1_spearman'],
            (' | best-matched residual factor |r| %s' % np.round(out['resid_factor_match_absr'], 2)) if 'resid_factor_match_absr' in out else ''))
    return out


if __name__ == '__main__':
    args = sys.argv[1:]; domain = 'band'
    if args and args[0] in ('band', 'dff'): domain = args[0]; args = args[1:]
    A, B = (args[:2] if len(args) >= 2 else ('AnimalRB24_261105_1423', 'AnimalRB24_261105_1520'))
    o1 = cross(A, B, domain); o2 = cross(B, A, domain)
    with open(os.path.join(fm.CACHE, 'cross_%s_%s_%s.pkl' % (A[-4:], B[-4:], domain)), 'wb') as fh: pickle.dump(dict(ab=o1, ba=o2), fh)
