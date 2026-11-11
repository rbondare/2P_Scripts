"""Axon grouping layer, step 1: load ONE session as a "universe" of ROIs (accepted + rejected + hidden merge constituents) on the common 5 Hz grid.

Read-only. Nothing is written to Z:. All ROIs go through the SAME F / Fneu -> best dF/F recipe (xs_common.best_dff_grid), so an accepted and a rejected ROI
are directly comparable. Accepted ROIs are the suite2p iscell == 1 ROIs in suite2p order (the order of the xs file: checked, see order_ok).

  U = Universe('AnimalRB24_261021_1630')
  U.dff        N x T best dF/F on the 5 Hz grid (all ROIs)       U.src   0 accepted, 1 rejected, 2 hidden merge constituent (inmerge > 0, iscell 0)
  U.cls        '', 'sig', 'noise', 'neuropil' (accepted only; V3 criterion, xs_common.signal_like_v3)
  U.pix        list of (n_i x 2) pixel arrays [y, x];   U.um  um per pixel;   U.gap_um(...)  footprint gaps
  U.S          the xs_common.Session (behaviour on the grid: U.S.beh, rigid shift U.S.shift, ...)
"""
import os, sys, csv, numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'axon analysis', 'cross_session'))
import xs_common as xc                                          # noqa: E402
from xs_common import DT                                        # noqa: E402

DB_SESSIONS = r'Z:\Group Members\Rima\Recordings_DB\sessions.csv'
CACHE = r'C:\Users\rbondarenko\2P_cache\axon_groups'
OUT = os.path.join(REPO, 'axon analysis', 'plots', 'axon_groups')
BAD_SHIFT_UM = 0.93                     # rigid-shift deviation from the 100 s running mean above which the frames are not used: 6 px at 0.155 um/px (best_dff.py of the 1529 analysis used 6 px); in um so that the high-zoom sessions (0.10 um/px) are not masked more strictly
MIN_FC = 5.0                            # rejected ROIs with median Fc (F - 0.5 Fneu) <= this have no signal at all (as xs_common._load_rejected / rejected_rois)


def um_per_px(name):
    """um per pixel from the recordings database (header FOV / pixels); the fov_matching inventory only covers RB24"""
    with open(DB_SESSIONS, encoding='utf-8') as f:
        for r in csv.DictReader(f):
            if r['session'] == name:
                x, y = float(r['um_per_px_x']), float(r['um_per_px_y'])
                return 0.5 * (x + y)
    raise KeyError('%s not in %s' % (name, DB_SESSIONS))


def v3_features(D, fs=1 / DT, thr=4.0, smooth_s=0.5):
    """V3 transient criterion on grid dF/F (rows = ROIs, no NaN): copy of paper_figures/rejected_rois.v3_features (that module imports fig_common, which creates
    output folders and sets plot styles at import time, so it is not imported here). passes = slow variance >= thr x white-noise expectation and upward events >= 2 x downward."""
    D = D.astype(np.float32); w = max(1, int(round(smooth_s * fs))); dsm = uniform_filter1d(D, w, axis=1, mode='nearest'); dev = dsm - np.median(dsm, axis=1, keepdims=True)
    neg = dev < 0; sig_neg = np.maximum(np.sqrt((dev ** 2 * neg).sum(1) / np.maximum(neg.sum(1), 1)), 1e-7)
    dd = np.diff(D, axis=1); hf = 1.4826 * np.median(np.abs(dd - np.median(dd, axis=1, keepdims=True)), axis=1) / np.sqrt(2)
    up = dev > 3 * sig_neg[:, None]; dn = dev < -3 * sig_neg[:, None]
    cnt = lambda m: (np.diff(np.concatenate([np.zeros((len(m), 1), bool), m], 1).astype(np.int8), axis=1) == 1).sum(1)
    ev_pos, ev_neg = cnt(up), cnt(dn); slow = dsm.var(axis=1) / (hf ** 2 / w + 1e-12)
    return dict(slow_ratio=slow, ev_pos=ev_pos, ev_neg=ev_neg, passes=(slow >= thr) & (ev_pos >= 2 * np.maximum(ev_neg, 1)), hf=hf)


def fill_rows(D):
    med = np.nanmedian(D, axis=1, keepdims=True); return np.where(np.isfinite(D), D, med)


class Universe:
    def __init__(self, name, verbose=True):
        self.name = name
        S = xc.Session(name, with_rejected=False, verbose=False); self.S = S
        p = os.path.join(xc.DATA_DIR, name, 'suite2p', 'plane0'); self.plane_dir = p
        ic = np.load(os.path.join(p, 'iscell.npy'))[:, 0] > 0
        stat = np.load(os.path.join(p, 'stat.npy'), allow_pickle=True)
        acc = np.where(ic)[0]; rej = np.where(~ic)[0]
        F = np.load(os.path.join(p, 'F.npy'), mmap_mode='r'); Fn = np.load(os.path.join(p, 'Fneu.npy'), mmap_mode='r')
        if F.shape[1] != len(S.tca): raise ValueError('suite2p F has %d frames, session %d' % (F.shape[1], len(S.tca)))
        xsF = np.atleast_2d(S.X['F'])
        self.order_ok = bool(len(acc) == S.nroi and np.corrcoef(xsF[0, :3000], np.asarray(F[acc[0], :3000]))[0, 1] > 0.99)
        if not self.order_ok: raise ValueError('%s: xs ROI order does not match iscell==1 order of suite2p (%d xs rows, %d accepted)' % (name, S.nroi, len(acc)))
        idx = np.concatenate([acc, rej]); Fg = xc.bin_rows(S.tca, np.asarray(F[idx]).astype(np.float32), S.edges); Fng = xc.bin_rows(S.tca, np.asarray(Fn[idx]).astype(np.float32), S.edges)
        dff, Fc, F0 = xc.best_dff_grid(Fg, Fng)
        keep = np.r_[np.ones(len(acc), bool), np.nanmedian(Fc[len(acc):], 1) > MIN_FC]
        self.idx = idx[keep]; self.dff = dff[keep]; self.Fg = Fg[keep]; self.Fng = Fng[keep]; self.Fc = Fc[keep]; self.F0 = np.nanmedian(F0[keep], 1)
        self.N = len(self.idx); self.n_acc = len(acc); self.T = dff.shape[1]; self.t = S.t
        inm = np.array([float(stat[k].get('inmerge', 0)) for k in self.idx])
        self.inmerge = inm                                                                    # > 0: hidden constituent; value = suite2p index of the merged ROI that contains it
        self.src = np.where(np.arange(self.N) < self.n_acc, 0, np.where(inm > 0, 2, 1)).astype(int)
        self.cls = np.array(['sig' if S.sig[i] else ('neuropil' if S.neuropil[i] else 'noise') for i in range(self.n_acc)] + [''] * (self.N - self.n_acc), dtype=object)
        self.pix = [np.c_[np.asarray(stat[k]['ypix'], int), np.asarray(stat[k]['xpix'], int)] for k in self.idx]
        self.med = np.array([np.asarray(stat[k]['med'], float) for k in self.idx]); self.npix = np.array([len(q) for q in self.pix])
        self.skew = None
        self.um = um_per_px(name)
        self.F_over_Fneu = np.nanmedian(self.Fg, 1) / np.maximum(np.nanmedian(self.Fng, 1), 1e-3)
        self.mask_bad()
        bad = self.bad
        if verbose:
            print('%s: %d accepted (%d sig / %d noise / %d neuropil), %d rejected (%d hidden constituents), %.1f min, %.3f um/px, %.1f%% frames masked (shift)' % (
                name, self.n_acc, (self.cls == 'sig').sum(), (self.cls == 'noise').sum(), (self.cls == 'neuropil').sum(), (self.src > 0).sum(), (self.src == 2).sum(),
                S.dur / 60, self.um, 100 * bad.mean()))

    def mask_bad(self):
        """frames not used: large rigid-shift excursion (brain motion); recomputed after unpickling so that a changed rule needs no reload from Z:"""
        bad = np.zeros(self.T, bool); S = self.S
        if S.shift is not None:
            sh = np.where(np.isfinite(S.shift), S.shift, 0.0); b = sh * self.um > BAD_SHIFT_UM; bad = uniform_filter1d(b.astype(float), 11) > 0
        self.bad = bad

    def time_slice(self, a, b):
        """view of the grid bins a:b (footprint gaps are kept): halves for held-out tests. dff, bad, behaviour and time are cut; transient test / skew are recomputed on demand."""
        import copy
        V = copy.copy(self); V.dff = self.dff[:, a:b]; V.T = b - a; V.t = self.t[a:b]; V.bad = self.bad[a:b]
        S = copy.copy(self.S); S.beh = {k: (np.asarray(x)[a:b] if np.ndim(x) == 1 and len(x) == self.T else x) for k, x in self.S.beh.items()}; V.S = S
        V.skew = None
        if hasattr(V, '_v3'): del V._v3
        return V

    # ------------------------------------------------------------------ transient test / skew (grid dF/F, NaN median-filled)
    def v3(self):
        if not hasattr(self, '_v3'): self._v3 = v3_features(fill_rows(self.dff))
        return self._v3

    def get_skew(self):
        if self.skew is None:
            from scipy.stats import skew
            self.skew = skew(fill_rows(self.dff), axis=1)
        return self.skew

    # ------------------------------------------------------------------ footprint geometry
    def gap_um(self, dmax_um):
        """N x N matrix of the minimum pixel-to-pixel distance between footprints in um (inf beyond dmax_um; 0 = overlapping or touching). Window search per ROI + cKDTree,
        exact inside dmax. Cached per dmax."""
        key = float(dmax_um)
        if getattr(self, '_gap_key', None) == key: return self._gap
        r = dmax_um / self.um; N = self.N
        allp = np.concatenate(self.pix); lab = np.concatenate([np.full(len(q), i) for i, q in enumerate(self.pix)])
        order = np.argsort(allp[:, 0], kind='stable'); allp, lab = allp[order], lab[order]; ys = allp[:, 0]
        G = np.full((N, N), np.inf, np.float32)
        for i in range(N):
            P = self.pix[i]; y0, y1 = P[:, 0].min() - r, P[:, 0].max() + r; x0, x1 = P[:, 1].min() - r, P[:, 1].max() + r
            a, b = np.searchsorted(ys, y0), np.searchsorted(ys, y1, side='right'); W = allp[a:b]; L = lab[a:b]
            m = (W[:, 1] >= x0) & (W[:, 1] <= x1); W, L = W[m], L[m]
            if len(W) == 0: continue
            d, _ = cKDTree(P).query(W); d = d * self.um
            mn = np.full(N, np.inf); np.minimum.at(mn, L, d); G[i] = np.where(mn <= dmax_um, mn, np.inf)
        G = np.minimum(G, G.T); np.fill_diagonal(G, 0)
        self._gap, self._gap_key = G, key
        return G


GAP_SCAN_UM = 60.0


def load_universe(name, refresh=False, verbose=True):
    """Universe with the footprint-gap matrix (GAP_SCAN_UM), cached as a pickle in CACHE (loading from Z: takes 25-50 s, the gaps another 25 s)"""
    import pickle
    f = os.path.join(CACHE, name + '_U.pkl')
    if os.path.exists(f) and not refresh:
        with open(f, 'rb') as fh: U = pickle.load(fh)
        U.mask_bad(); return U
    U = Universe(name, verbose=verbose); U.gap_um(GAP_SCAN_UM); U.v3(); U.get_skew()
    os.makedirs(CACHE, exist_ok=True)
    with open(f, 'wb') as fh: pickle.dump(U, fh, protocol=4)
    return U
