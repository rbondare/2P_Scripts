"""Resting motion energy (REST) of the whisking region, per video codec, from sessions already processed by face_motion.py.

    python calibrate_rest.py <whisker csv or folder of sessions> [...] [--max-p10-over-p3 1.15]

REST is the median over the sessions of one codec that have a clear rest period of the per-session 3rd percentile of the
smoothed ME of face_motion.WHISKER_BOX. A session has a clear rest period when its 10th percentile is close to its 3rd
(the lowest tenth of the session is a plateau at rest, not a continuous climb); a mouse that whisks or runs the whole
session has no rest and would bias REST upwards. The codec is read from the _whisker.json next to the csv (written by
face_motion.py). Prints one line per session and the REST line to paste into face_motion.REST.
Add --plot <png> for the distribution.
"""
import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import face_motion as fm


def session_stats(csv_path):
    T = pd.read_csv(csv_path)
    info = json.load(open(csv_path.replace('_whisker.csv', '_whisker.json')))
    sm = gaussian_filter1d(T['me_' + fm.WHISKER_BOX].values.astype(float), fm.SMOOTH_SIGMA_FRAMES)
    ok = T['artifact'].values == 0
    p3, p10, p50 = np.percentile(sm[ok], [3, 10, 50])
    return dict(session=info['session'], codec=info['codec'], frames=len(T), p3=p3, p10=p10, p50=p50, p10_over_p3=p10 / p3,
                ctrl=float(np.median(T['me_ctrl'])))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('paths', nargs='+')
    ap.add_argument('--max-p10-over-p3', type=float, default=1.15)
    ap.add_argument('--plot')
    args = ap.parse_args()
    files = []
    for p in args.paths:
        files += sorted(glob.glob(os.path.join(p, '**', '*_whisker.csv'), recursive=True)) if os.path.isdir(p) else [p]
    S = pd.DataFrame([session_stats(f) for f in files])
    S['clear_rest'] = S.p10_over_p3 <= args.max_p10_over_p3
    pd.set_option('display.width', 200)
    print(S.sort_values(['codec', 'session']).round(3).to_string(index=False))
    print()
    rest = {}
    for codec, g in S.groupby('codec'):
        r = g[g.clear_rest]
        print('%s: %d sessions, %d with a clear rest period; p3 of those: median %.3f, range %.3f-%.3f (CV %.0f%%); all sessions median p3 %.3f; ctrl ME median %.3f' % (
            codec, len(g), len(r), r.p3.median(), r.p3.min(), r.p3.max(), 100 * r.p3.std() / r.p3.mean(), g.p3.median(), g.ctrl.median()))
        rest[codec] = round(float(r.p3.median()), 3)
    print('\nREST = %s' % rest)
    if args.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 4.5))
        for codec, col in (('mpeg4', 'tab:blue'), ('h264', 'tab:red')):
            g = S[S.codec == codec]
            ax.scatter(g.p3, g.p10_over_p3, c=[col if c else 'none' for c in g.clear_rest], edgecolors=col, label=codec)
        ax.axhline(args.max_p10_over_p3, c='gray', ls=':')
        ax.set_xlabel('3rd percentile of the smoothed ME (front_top)')
        ax.set_ylabel('p10 / p3')
        ax.legend()
        fig.tight_layout()
        fig.savefig(args.plot, dpi=90)


if __name__ == '__main__':
    main()
