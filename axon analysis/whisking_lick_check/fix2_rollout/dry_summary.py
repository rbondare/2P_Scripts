"""Dry run: classify every session from the pass-1 ME with the calibrated REST (nothing written) and compare with the old results."""
import os, sys, glob, json
import numpy as np, pandas as pd
sys.path.insert(0, r'C:\Users\rbondarenko\projects\2P_Scripts\DLC\face_model')
import face_motion as fm
here = os.path.dirname(os.path.abspath(__file__))
rows = []
for jf in sorted(glob.glob(os.path.join(here, 'pass1', '*', '*_whisker.json'))):
    s = os.path.basename(jf).replace('_whisker.json', '')
    info = json.load(open(jf))
    T = pd.read_csv(jf.replace('.json', '.csv'))
    me = {k: T['me_' + k].values.astype(float) for k in fm.BOXES}
    res = fm.classify(me, fm.REST[info['codec']], info['codec'])
    idx, art, on = res['index'], res['artifact'], res['on']
    old = pd.read_csv(os.path.join(r'C:\Users\rbondarenko\DLC\work', s, 'whisker_results', s + '_whisker.csv'))
    oart = old.artifact.values.astype(bool)
    rows.append(dict(session=s.replace('Animal', ''), codec=info['codec'], p5=np.percentile(idx[~art], 5), p50=np.median(idx[~art]), p95=np.percentile(idx[~art], 95),
                     whisking=on[~art].mean(), n_bouts=len(fm.runs(on)), art_new=art.mean(), art_old=oart.mean(), old_not_in_new=(oart & ~art).mean(),
                     old_p5=old.whisk_index.values[~oart].min() if False else np.percentile(old.whisk_index.values[~oart], 5),
                     old_p50=np.median(old.whisk_index.values[~oart]), old_whisking=np.nanmean(old.whisk_on.values)))
D = pd.DataFrame(rows)
pd.set_option('display.width', 250)
print(D.round(3).to_string(index=False))
print()
for c, g in D.groupby('codec'):
    print('%s: new index p5 median %.2f (range %.2f-%.2f); p50 median %.2f; old index p5 median %.2f p50 median %.2f; artifact new %.2f%% old %.2f%%; old-flag frames not flagged now: max %.3f%%' % (
        c, g.p5.median(), g.p5.min(), g.p5.max(), g.p50.median(), g.old_p5.median(), g.old_p50.median(), 100 * g.art_new.mean(), 100 * g.art_old.mean(), 100 * g.old_not_in_new.max()))
D.to_csv(os.path.join(here, 'dry_summary.csv'), index=False)
