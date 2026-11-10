import sys, numpy as np, h5py
sys.path.insert(0, r'C:\Users\rbondarenko\projects\2P_Scripts\axon analysis\cross_session')
import xs_common as xc
from xs_common import DT, nsmooth
def hp(x, s=30.0): return x - nsmooth(x, s / DT)
def corr(a, b):
    ok = np.isfinite(a) & np.isfinite(b); return np.corrcoef(a[ok], b[ok])[0, 1]
expected = {'AnimalRB24_261105_1423': 0.519, 'AnimalRB24_261105_1520': 0.510, 'AnimalRB24_261030_1744': 0.646, 'AnimalRB24_261104_1651': 0.577}
for s, exp in expected.items():
    ses = xc.Session(s, with_rejected=False, verbose=False)
    f = h5py.File(r'Z:\Group Members\Rima\Aggregated\%s_preprocessed.mat' % s, 'r')
    on = np.array(f['Behav']['WhiskOn']).ravel().astype(float)
    tcam = np.asarray(ses.X['tcam'], float)
    assert len(on) == len(tcam)
    g = xc.bin_rows(tcam, on, ses.edges)[0]
    y = hp(ses.population('sig'))
    x = np.full_like(g, np.nan); x[1:] = hp(g)[:-1]          # descriptor leads by 1 bin (0.2 s), as in coupling_check
    print('%s: r(aggregated WhiskOn, signal-like pop) = %.3f   (coupling_check candidate: %.3f)' % (s, corr(x, y), exp))
