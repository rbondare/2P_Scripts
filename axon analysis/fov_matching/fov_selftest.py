"""Self-test of the registration: a real mean image (11-05 1520) is re-sampled at a DIFFERENT pixel size with a KNOWN offset, cropped, noise added; the pair registration
(fov_register.register_pair, coarse) + fine refinement (fov_refine.fine_edge) must recover the offset. Tests several offsets / pixel sizes / noise levels.
Run: miniforge3\\envs\\suite2p_venv\\python.exe fov_selftest.py
"""
import numpy as np
from types import SimpleNamespace
from scipy import ndimage as ndi
from fov_common import *
from fov_register import features, register_pair
from fov_refine import fine_edge, FG


def synth(A, t, px_b, shape_b, noise_frac, rng):
    """session B sampled from A's image at pixel size px_b, B's pixel (0,0) located at t (um) in A's frame"""
    ii, jj = np.meshgrid(np.arange(shape_b[0]), np.arange(shape_b[1]), indexing='ij')
    yy = (t[0] + ii * px_b) / A.py; xx = (t[1] + jj * px_b) / A.px
    # smooth to B's resolution first (box of one B pixel), then sample
    sm = ndi.gaussian_filter(A.mean, (0.5 * px_b / A.py, 0.5 * px_b / A.px))
    g = ndi.map_coordinates(sm, [yy, xx], order=1, mode='nearest')
    g = g + rng.normal(0, noise_frac * np.std(g), g.shape)
    valid = np.zeros(shape_b, bool); valid[5:-5, 5:-5] = True
    return SimpleNamespace(mean=g, valid=valid, py=px_b, px=px_b)


def main():
    inv = read_inventory(); A = Sess('AnimalRB24_261105_1520', inv, rois=False)
    rng = np.random.default_rng(0)
    FA = features(A)
    FAf = features(A, grid=FG)
    print('known offset (um)      pixel  noise | coarse peak (um)        err   NCC | fine (um)              err')
    errs = []
    for t, px_b, nf in [((5.0, 6.3), 0.138, 0.3), ((-4.4, 20.7), 0.211, 0.3), ((8.1, -3.9), 0.278, 0.3), ((2.2, 30.5), 0.138, 1.0), ((6.0, 12.0), 0.278, 1.0)]:
        shape_b = (int(40 / px_b), int(80 / px_b))
        B = synth(A, t, px_b, shape_b, nf, rng)
        FB = features(B)
        r = register_pair(FA, FB, rng, n_surr=2)
        p = r['peaks'][0]; ec = np.hypot(p['t'][0] - t[0], p['t'][1] - t[1])
        FBf = features(B, grid=FG); fe = fine_edge(FAf, FBf, np.array(p['t']))
        ef = np.hypot(fe[0][0] - t[0], fe[0][1] - t[1]) if fe else np.nan; errs.append(ef)
        print('(%5.1f, %5.1f)        %.3f  %.1f  | (%6.2f, %6.2f)   %.2f   %.2f | (%6.2f, %6.2f)   %.2f' % (t[0], t[1], px_b, nf, p['t'][0], p['t'][1], ec, p['r'], fe[0][0], fe[0][1], ef))
    print('median fine error %.2f um, max %.2f um' % (np.nanmedian(errs), np.nanmax(errs)))


if __name__ == '__main__':
    main()
