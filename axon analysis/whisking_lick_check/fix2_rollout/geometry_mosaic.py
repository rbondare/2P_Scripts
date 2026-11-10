"""Contact sheet of every session: middle frame around the nose with the new whisking box (front_top) and the old boxes drawn,
to check that the spout never enters front_top. Also prints codec / bitrate per session."""
import os, sys, json, glob, subprocess
import numpy as np, cv2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
WORK = r'C:\Users\rbondarenko\DLC\work'
FF = r'C:\Users\rbondarenko\miniforge3\envs\dlc-tf\Library\bin\ffprobe.exe'
TOP = (-260, -110, 60, 0)
FRONT = (-260, -110, 60, 200)
DOWN = (60, 80, 300, 300)
PAD = (60, -100, 260, 60)

sessions = sorted(d for d in os.listdir(WORK) if d.startswith('Animal'))
rows = []
tiles = []
for s in sessions:
    vids = glob.glob(os.path.join(WORK, s, 'camera', '*.mp4'))
    js = glob.glob(os.path.join(WORK, s, 'whisker_results', '*_whisker.json'))
    if not vids:
        continue
    v = vids[0]
    out = subprocess.run([FF, '-v', 'error', '-select_streams', 'v:0', '-show_entries', 'stream=codec_name,bit_rate,nb_frames',
                          '-of', 'csv=p=0', v], capture_output=True, text=True).stdout.strip().split(',')
    codec, br = out[0], int(out[1]) / 1e6
    ax, ay = (None, None)
    if js:
        j = json.load(open(js[0]))
        ax, ay = j['anchor_median_xy_px']
    rows.append((s, codec, round(br, 1), ax, ay))
    print(s, codec, round(br, 1), ax, ay, flush=True)
    if ax is None:
        continue
    cap = cv2.VideoCapture(v)
    n = int(cap.get(7))
    cap.set(cv2.CAP_PROP_POS_FRAMES, n // 2)
    ok, fr = cap.read()
    cap.release()
    if not ok:
        continue
    x0, y0 = max(ax - 360, 0), max(ay - 200, 0)
    crop = fr[y0:ay + 340, x0:ax + 340]
    tiles.append((s, codec, crop, ax - x0, ay - y0))

nc = 6
nr = int(np.ceil(len(tiles) / nc))
fig, axs = plt.subplots(nr, nc, figsize=(nc * 4.2, nr * 3.4))
axs = np.atleast_2d(axs)
for a in axs.ravel():
    a.axis('off')
for a, (s, codec, crop, cx, cy) in zip(axs.ravel(), tiles):
    a.imshow(cv2.cvtColor(crop, cv2.COLOR_BGR2RGB))
    for box, col, lw in ((TOP, 'lime', 2), (FRONT, 'yellow', 1), (DOWN, 'magenta', 1), (PAD, 'orange', 1)):
        a.add_patch(plt.Rectangle((cx + box[0], cy + box[1]), box[2] - box[0], box[3] - box[1], fill=False, ec=col, lw=lw))
    a.plot(cx, cy, 'r.', ms=5)
    a.set_title('%s %s' % (s.replace('Animal', ''), codec), fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(here, 'geometry_mosaic.png'), dpi=60)
json.dump(rows, open(os.path.join(here, 'codec_table.json'), 'w'))
print('done')
