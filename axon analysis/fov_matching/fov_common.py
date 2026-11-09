"""Shared helpers of the RB24 cross-session FOV / axon-location analysis (axon analysis/fov_matching).

Everything is read from the suite2p output of each session (ops.npy, stat.npy, iscell.npy) and the first TIFF header; no data.bin is touched.
Coordinates: images are resampled to ONE common grid of GRID um/px (separate x / y scale per session, because the 512x1024 frames have square pixels
but different zoom and the 512x512 frames a different pixel size). All positions below are in um unless a name ends in _px.
"""
import os, sys, csv, glob, pickle, numpy as np
from scipy import ndimage as ndi

REPO = r'c:\Users\rbondarenko\projects\2P_Scripts'
DATA = r'Z:\Group Members\Rima\DATA_2P'
HERE = os.path.join(REPO, 'axon analysis', 'fov_matching')
OUT = os.path.join(REPO, 'axon analysis', 'plots', 'fov_matching')
CACHE = os.path.join(OUT, 'cache')
for d in (OUT, CACHE): os.makedirs(d, exist_ok=True)
QC_PATH = os.path.join(REPO, 'axon analysis', 'plots', 'cross_session', 'session_qc.csv')

GRID = 0.25                       # um per pixel of the common grid (native 0.10-0.28 um/px; thin axons ~1 um)

# the 20 RB24 sessions with suite2p output that are in scope (261017_1620 = 3 planes 40 um apart, 261109_1506 = suite2p not finished: both excluded)
SESSIONS = ['AnimalRB24_261021_1630', 'AnimalRB24_261021_1826', 'AnimalRB24_261022_2026', 'AnimalRB24_261023_1453', 'AnimalRB24_261023_1527', 'AnimalRB24_261023_1538',
            'AnimalRB24_261027_1456', 'AnimalRB24_261027_1534', 'AnimalRB24_261029_1752', 'AnimalRB24_261030_1744', 'AnimalRB24_261104_1540', 'AnimalRB24_261104_1640',
            'AnimalRB24_261104_1651', 'AnimalRB24_261105_1423', 'AnimalRB24_261105_1520', 'AnimalRB24_261105_1549', 'AnimalRB24_261106_1317', 'AnimalRB24_261106_1529',
            'AnimalRB24_261106_1555', 'AnimalRB24_261106_1604']
# sessions recorded at the same stage position on the same day (identical axesPosition): same-day same-FOV groups = the ceiling of what matching can do
SAME_DAY = {'1023': ['AnimalRB24_261023_1453', 'AnimalRB24_261023_1527', 'AnimalRB24_261023_1538'],
            '1027': ['AnimalRB24_261027_1456', 'AnimalRB24_261027_1534'],
            '1104': ['AnimalRB24_261104_1640', 'AnimalRB24_261104_1651'],
            '1105': ['AnimalRB24_261105_1423', 'AnimalRB24_261105_1520', 'AnimalRB24_261105_1549'],
            '1106': ['AnimalRB24_261106_1529', 'AnimalRB24_261106_1555', 'AnimalRB24_261106_1604', 'AnimalRB24_261106_1317']}
# the stimulus-folder date of 1826 is 10-28, the DATA_2P folder says 10-21 (typo); calendar date used in the analysis comes from the TIFF time stamp (fov_inventory.csv)


def short(nm):
    d, t = nm.replace('AnimalRB24_', '').split('_')
    return '%s-%s %s' % (d[2:4], d[4:6], t)


def s2p_dir(nm): return os.path.join(DATA, nm, 'suite2p', 'plane0')


# ---------------------------------------------------------------- inventory (written by fov_inventory.py)
def read_inventory():
    p = os.path.join(OUT, 'fov_inventory.csv')
    rows = {}
    for r in csv.DictReader(open(p)): rows[r['session']] = r
    return rows


def pixel_size(nm, inv=None):
    """(um per pixel along y (lines), um per pixel along x) from the header FOV"""
    inv = inv or read_inventory(); r = inv[nm]
    return float(r['um_per_px_y']), float(r['um_per_px_x'])


# ---------------------------------------------------------------- one session's images and ROIs in full-frame pixel coordinates
class Sess:
    """meanImg (full frame), max_proj / Vcorr padded back into the full frame (suite2p crops them to ops.yrange / xrange), valid mask, ROI table."""

    def __init__(self, name, inv=None, rois=True):
        self.name = name
        p = s2p_dir(name)
        ops = np.load(os.path.join(p, 'ops.npy'), allow_pickle=True).item()
        self.ops = ops
        self.py, self.px = pixel_size(name, inv)
        m = np.asarray(ops['meanImg'], float)
        self.shape = m.shape
        yr, xr = ops['yrange'], ops['xrange']
        self.yr, self.xr = (int(yr[0]), int(yr[1])), (int(xr[0]), int(xr[1]))
        self.valid = np.zeros(m.shape, bool); self.valid[yr[0]:yr[1], xr[0]:xr[1]] = True
        self.mean = m
        self.max_proj = self._pad(ops.get('max_proj'))
        self.vcorr = self._pad(ops.get('Vcorr'))
        self.rois = None
        if rois: self._load_rois(p)

    def _pad(self, a):
        if a is None: return None
        a = np.asarray(a, float); out = np.full(self.shape, np.nan)
        y0, y1 = self.yr; x0, x1 = self.xr
        if a.shape != (y1 - y0, x1 - x0): raise ValueError('%s: cropped image %s does not fit yrange/xrange %s %s' % (self.name, a.shape, self.yr, self.xr))
        out[y0:y1, x0:x1] = a
        return out

    def _load_rois(self, p):
        stat = np.load(os.path.join(p, 'stat.npy'), allow_pickle=True)
        ic = np.load(os.path.join(p, 'iscell.npy'))
        if len(stat) != len(ic): raise ValueError('%s: stat %d rows, iscell %d rows' % (self.name, len(stat), len(ic)))
        self.n_total = len(stat); self.accepted = ic[:, 0] > 0; self.prob = ic[:, 1]
        self.stat = stat

    def roi_pix(self, i):
        s = self.stat[i]; return np.asarray(s['ypix'], int), np.asarray(s['xpix'], int), np.asarray(s['lam'], float)


# ---------------------------------------------------------------- resampling to the common um grid
def to_grid(a, py, px, order=1, cval=0.0, grid=GRID):
    """resample a full-frame pixel image (rows = y with py um/px, cols = x with px um/px) to a grid of `grid` um/px; pixel centres keep their um position
    (output pixel (i, j) is at y = i * grid, x = j * grid measured from the centre of input pixel (0, 0))"""
    ny = int(np.floor((a.shape[0] - 1) * py / grid)) + 1; nx = int(np.floor((a.shape[1] - 1) * px / grid)) + 1
    yy, xx = np.meshgrid(np.arange(ny) * grid / py, np.arange(nx) * grid / px, indexing='ij')
    return ndi.map_coordinates(a, [yy, xx], order=order, mode='constant', cval=cval, prefilter=False)


def dog(a, valid, s_small_um, s_large_um, py, px):
    """difference of Gaussians of a full-frame image in um units (sigmas in um), invalid border filled with the median first, normalised convolution at the border"""
    v = valid.astype(float); med = np.median(a[valid])
    b = np.where(valid, a, med)
    g1 = ndi.gaussian_filter(b, (s_small_um / py, s_small_um / px), mode='nearest')
    g2 = ndi.gaussian_filter(b, (s_large_um / py, s_large_um / px), mode='nearest')
    return g1 - g2


def robust_z(a, mask):
    x = a[mask]; med = np.median(x); mad = np.median(np.abs(x - med)) * 1.4826 + 1e-12
    return (a - med) / mad
