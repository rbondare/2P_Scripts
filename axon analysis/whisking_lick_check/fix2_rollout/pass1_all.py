"""Pass 1: motion energy of every box (incl. the new front_top) for every session -> scratch folder fm2/pass1/<session>/.
Nothing in the real whisker_results folders is touched. Skips sessions that are done. 3 parallel ffmpeg readers."""
import os, sys, glob, subprocess, time
from concurrent.futures import ThreadPoolExecutor

here = os.path.dirname(os.path.abspath(__file__))
WORK = r'C:\Users\rbondarenko\DLC\work'
PY = r'C:\Users\rbondarenko\miniforge3\envs\dlc-tf\python.exe'
FM = r'C:\Users\rbondarenko\projects\2P_Scripts\DLC\face_model\face_motion.py'
env = dict(os.environ)
env['PATH'] = r'C:\Users\rbondarenko\miniforge3\envs\dlc-tf\Library\bin;' + env['PATH']

priority = ['AnimalRB24_261030_1744', 'AnimalRB24_261105_1549', 'AnimalRB24_261109_1506', 'AnimalRB22_261029_1516', 'AnimalRB23_261029_1704']
sessions = sorted(d for d in os.listdir(WORK) if d.startswith('Animal') and glob.glob(os.path.join(WORK, d, 'face_samples', '*_face.csv')))
order = [s for s in priority if s in sessions] + [s for s in sessions if s not in priority]
log = open(os.path.join(here, 'pass1.log'), 'a', buffering=1)


def run(s):
    out = os.path.join(here, 'pass1', s)
    if os.path.isfile(os.path.join(out, s + '_whisker.json')):
        return s, 'skip', 0
    v = glob.glob(os.path.join(WORK, s, 'camera', '*.mp4'))[0]
    f = glob.glob(os.path.join(WORK, s, 'face_samples', '*_face.csv'))[0]
    t0 = time.time()
    r = subprocess.run([PY, '-u', FM, v, f, '--out-dir', out, '--name', s], capture_output=True, text=True, env=env)
    msg = (r.stdout.strip().splitlines() or [''])[-1] if r.returncode == 0 else 'FAILED %s' % r.stderr[-400:]
    log.write('%s  %.0f s  %s\n' % (s, time.time() - t0, msg))
    return s, msg, time.time() - t0


with ThreadPoolExecutor(3) as ex:
    for s, msg, dt in ex.map(run, order):
        pass
log.write('ALL DONE\n')
