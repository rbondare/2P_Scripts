"""Read ONLY the per-frame ScanImage ImageDescription of every TIFF of a session (no pixel data) and keep the aux trigger lists.
usage: python read_triggers.py <session> [workers]
Output (one npz per TIFF in <here>/trig_<session>/): frame numbers, frame timestamps (s), aux0..aux3 event times (s) and the frame index each event came from.
aux0 = projector red-sync photodiode, aux1 = camera frame trigger, aux2 = ball camera, aux3 = lick line."""
import sys, os, re, time, glob
import numpy as np
import tifffile
from concurrent.futures import ThreadPoolExecutor

sess = sys.argv[1]
nw = int(sys.argv[2]) if len(sys.argv) > 2 else 12
here = os.path.dirname(os.path.abspath(__file__))
out_dir = os.path.join(here, 'trig_%s' % sess)
os.makedirs(out_dir, exist_ok=True)
files = sorted(glob.glob(r'Z:\Group Members\Rima\DATA_2P\%s\%s_*.tif' % (sess, sess)))
print(len(files), 'tifs', flush=True)
num = re.compile(r'^\s*(\w+)\s*=\s*(.*?)\s*$')


def parse(desc):
    r = {}
    for line in desc.replace('<output truncated>', '').split('\n'):
        m = num.match(line)
        if not m:
            continue
        k, v = m.group(1), m.group(2)
        if v.startswith('['):
            if not v.endswith(']'):
                break  # cut off: nothing after this can be trusted
            v = v.strip('[]').replace(',', ' ').split()
            r[k] = np.array(v, float)
        elif v.startswith('{'):
            continue
        else:
            try:
                r[k] = float(v)
            except ValueError:
                break
    return r


def do(fp):
    base = os.path.splitext(os.path.basename(fp))[0]
    op = os.path.join(out_dir, base + '.npz')
    if os.path.exists(op):
        return base, 0.0, 0
    t0 = time.time()
    fn, ts, ev = [], [], {0: ([], []), 1: ([], []), 2: ([], []), 3: ([], [])}
    nbad = 0
    with tifffile.TiffFile(fp) as tf:
        for i, pg in enumerate(tf.pages):
            d = pg.tags['ImageDescription'].value
            r = parse(d)
            if 'frameNumbers' not in r:
                nbad += 1
                fn.append(np.nan)
                ts.append(np.nan)
                continue
            fn.append(r['frameNumbers'])
            ts.append(r.get('frameTimestamps_sec', np.nan))
            for k in range(4):
                a = r.get('auxTrigger%d' % k)
                if a is not None and len(a):
                    ev[k][0].extend(a.tolist())
                    ev[k][1].extend([i] * len(a))
    np.savez(op, fn=np.array(fn), ts=np.array(ts), nbad=nbad,
             **{'aux%d_t' % k: np.array(ev[k][0]) for k in range(4)}, **{'aux%d_i' % k: np.array(ev[k][1]) for k in range(4)})
    return base, time.time() - t0, len(fn)


t00 = time.time()
done = 0
with ThreadPoolExecutor(nw) as ex:
    for base, dt, n in ex.map(do, files):
        done += 1
        print('%3d/%d %s  %d pages  %.0fs   (elapsed %.0f s)' % (done, len(files), base, n, dt, time.time() - t00), flush=True)
print('ALL DONE', flush=True)
