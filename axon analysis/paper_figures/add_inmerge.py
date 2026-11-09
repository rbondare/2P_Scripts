"""One-off patch: add the suite2p 'inmerge' flag (hidden constituent of a merged ROI) to the cached rejroi_<session>.pkl tables made before rejected_rois.py stored it.
usage: python add_inmerge.py [session ...]   (default: every cache that lacks it)"""
import os, sys, glob, pickle, numpy as np
from fig_common import *
import xs_common as xc

for f in sorted(glob.glob(os.path.join(CACHE, 'rejroi_*.pkl'))):
    d = pickle.load(open(f, 'rb')); name = d['name']
    if len(sys.argv) > 1 and name not in sys.argv[1:]: continue
    if 'inmerge' in d: continue
    stat = np.load(os.path.join(xc.DATA_DIR, name, 'suite2p', 'plane0', 'stat.npy'), allow_pickle=True)
    d['inmerge'] = np.array([float(stat[k].get('inmerge', 0)) for k in d['idx']]); pickle.dump(d, open(f, 'wb'))
    R = d['grp'] == 1
    print('%-26s rejected %3d, hidden merge constituents %3d (%.0f%%); of the activity-like ones %d' % (name, R.sum(), (d['inmerge'][R] > 0).sum(), 100 * (d['inmerge'][R] > 0).mean(), ((d['inmerge'][R] > 0) & d['passes'][R]).sum()))
