"""Per-pixel motion-energy maps (mean |frame t - frame t-1|) for 4 classes of 1 s windows: still/run x lick/no lick.  usage: python maps.py <session>"""
import sys, os, glob, json, numpy as np, pandas as pd, cv2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sess = sys.argv[1]
here = os.path.dirname(os.path.abspath(__file__))
vid = glob.glob(r'C:\Users\rbondarenko\DLC\work\%s\camera\*.mp4' % sess)[0]
J = json.load(open(r'C:\Users\rbondarenko\DLC\work\%s\whisker_results\%s_whisker.json' % (sess, sess)))
ax_, ay_ = J['anchor_median_xy_px']
boxes = J['boxes_rel_to_nose_px']
R = pd.read_pickle(os.path.join(here, 'win_%s.pkl' % sess))
rng = np.random.default_rng(5)
X0, X1, Y0, Y1 = max(ax_ - 340, 0), ax_ + 340, max(ay_ - 200, 0), ay_ + 330
BS = 8
classes = [('still', 'nolick'), ('still', 'lick'), ('run', 'nolick'), ('run', 'lick')]
cap = cv2.VideoCapture(vid)
maps = {}
for rc, lc in classes:
    sel = R[(R.runcls == rc) & (R.lickcls == lc)]
    if len(sel) > 70:
        sel = sel.iloc[rng.choice(len(sel), 70, replace=False)]
    acc = None
    cnt = 0
    for a in sel.a.values:
        cap.set(cv2.CAP_PROP_POS_FRAMES, int(a) - 1)
        prev = None
        for i in range(int(R.b.iloc[0] - R.a.iloc[0]) + 1):
            ok, f = cap.read()
            if not ok:
                break
            g = cv2.GaussianBlur(cv2.cvtColor(f[Y0:Y1, X0:X1], cv2.COLOR_BGR2GRAY), (0, 0), 1.0).astype(np.float32)
            if prev is not None:
                df = np.abs(g - prev)
                if acc is None:
                    acc = np.zeros_like(df)
                acc += df
                cnt += 1
            prev = g
    m = acc / max(cnt, 1)
    h, w = (m.shape[0] // BS) * BS, (m.shape[1] // BS) * BS
    maps[(rc, lc)] = m[:h, :w].reshape(h // BS, BS, w // BS, BS).mean(axis=(1, 3))
    print(rc, lc, 'windows', len(sel), 'frames', cnt)
np.savez(os.path.join(here, 'maps_%s.npz' % sess), X0=X0, Y0=Y0, BS=BS, ax=ax_, ay=ay_, **{'%s_%s' % k: v for k, v in maps.items()})
cap.set(cv2.CAP_PROP_POS_FRAMES, int(R.a.iloc[len(R) // 2]))
ok, frame = cap.read()
bg = cv2.cvtColor(frame[Y0:Y1, X0:X1], cv2.COLOR_BGR2GRAY)
cols = {'fan_front': 'lime', 'fan_down': 'magenta', 'pad': 'orange', 'ctrl': 'cyan'}
fig, axs = plt.subplots(2, 4, figsize=(18, 7.5))
vmax = max(np.percentile(m, 99.5) for m in maps.values())
for j, (rc, lc) in enumerate(classes):
    a = axs[0, j]
    im = a.imshow(maps[(rc, lc)], extent=(X0, X1, Y1, Y0), cmap='magma', vmin=0, vmax=vmax)
    a.set_title('%s, %s' % (rc, 'licking' if lc == 'lick' else 'no licks'))
for j, (c1, c2, ttl) in enumerate([(('still', 'lick'), ('still', 'nolick'), 'licking - no licking (still)'), (('run', 'lick'), ('run', 'nolick'), 'licking - no licking (running)'),
                                  (('run', 'nolick'), ('still', 'nolick'), 'running - still (no licks)')]):
    a = axs[1, j]
    D = maps[c1] - maps[c2]
    v = np.percentile(np.abs(D), 99.5)
    a.imshow(D, extent=(X0, X1, Y1, Y0), cmap='RdBu_r', vmin=-v, vmax=v)
    a.set_title(ttl)
axs[1, 3].imshow(bg, cmap='gray', extent=(X0, X1, Y1, Y0))
axs[1, 3].set_title('example frame')
for a in axs.ravel():
    for k, (x0, y0, x1, y1) in boxes.items():
        a.add_patch(plt.Rectangle((ax_ + x0, ay_ + y0), x1 - x0, y1 - y0, fill=False, ec=cols[k], lw=1.6))
    a.plot(ax_, ay_, 'r.')
    a.set_xticks([])
    a.set_yticks([])
    a.set_xlim(X0, X1)
    a.set_ylim(Y1, Y0)
fig.suptitle('%s: where the pixels move (green = fan_front, magenta = fan_down, orange = pad, cyan = ctrl)' % sess)
fig.tight_layout()
fig.savefig(os.path.join(here, 'maps_%s.png' % sess), dpi=70)
