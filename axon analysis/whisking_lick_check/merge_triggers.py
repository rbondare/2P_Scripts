"""Merge the per-TIFF trigger npz files of a session into camera frame times, lick events and ball readings on ONE clock (the ScanImage / 2P clock).
usage: python merge_triggers.py <session>     -> align_<session>.npz
Follows save_preprocessed.process_triggers: frames are de-duplicated by frame number (first occurrence kept); camera = aux1, ball = aux2, lick line = aux3."""
import sys, os, glob, numpy as np, h5py

sess = sys.argv[1]
here = os.path.dirname(os.path.abspath(__file__))
fs = sorted(glob.glob(os.path.join(here, 'trig_%s' % sess, '*.npz')))
print(len(fs), 'files')
FN, TS, E = [], [], {k: [] for k in range(4)}
page0 = 0
seen = set()
fn_all = []
page_of_event = {k: [] for k in range(4)}
offs = []
for fp in fs:
    z = np.load(fp)
    fn = z['fn']
    ts = z['ts']
    n = len(fn)
    FN.append(fn)
    TS.append(ts)
    for k in range(4):
        E[k].append(np.c_[z['aux%d_t' % k], z['aux%d_i' % k] + page0])
    offs.append(page0)
    page0 += n
FN = np.concatenate(FN)
TS = np.concatenate(TS)
print('pages', len(FN), 'bad (no frame number)', int(np.isnan(FN).sum()))
# cumulative frame numbers are continuous across files? check
good = np.isfinite(FN)
d = np.diff(FN[good])
print('frame number steps: min %g median %g max %g; negative steps %d; zero steps %d' % (d.min(), np.median(d), d.max(), (d < 0).sum(), (d == 0).sum()))
# keep first occurrence of every frame number
_, first = np.unique(FN, return_index=True)
keep = np.zeros(len(FN), bool)
keep[first] = True
keep &= good
print('unique frames kept', keep.sum(), 'of', len(FN))
out = {}
for k in range(4):
    a = np.concatenate(E[k]) if E[k] else np.zeros((0, 2))
    m = keep[a[:, 1].astype(int)]
    t = a[m, 0]
    out['aux%d' % k] = t
    dt = np.diff(t)
    print('aux%d: %d events%s' % (k, len(t), ('  monotonic %s, min dt %.4f median dt %.4f' % ((dt >= 0).all(), dt.min(), np.median(dt))) if len(t) > 1 else ''))
np.savez(os.path.join(here, 'align_%s.npz' % sess), fn=FN, ts=TS, **out)
