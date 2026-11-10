"""Pass 2: classify with the calibrated REST (face_motion.py --reclassify, no decoding) and write the real whisker_results.
Local: DLC\\work\\<session>\\whisker_results (old version backed up once to whisker_results_v1) + work\\quality_reports\\<session>.
Archive (--archive): Z:\\Group Members\\Rima\\preCamera\\<animal>\\<session>\\whisker_results for sessions that have one (old version
backed up once to whisker_results_v1 next to it).
usage: python pass2_apply.py [--archive] [--only SESSION ...]"""
import os, sys, glob, shutil, subprocess, argparse, json

here = os.path.dirname(os.path.abspath(__file__))
WORK = r'C:\Users\rbondarenko\DLC\work'
PRE = r'Z:\Group Members\Rima\preCamera'
PY = r'C:\Users\rbondarenko\miniforge3\envs\dlc-tf\python.exe'
FM = r'C:\Users\rbondarenko\projects\2P_Scripts\DLC\face_model\face_motion.py'
env = dict(os.environ)
env['PATH'] = r'C:\Users\rbondarenko\miniforge3\envs\dlc-tf\Library\bin;' + env['PATH']

ap = argparse.ArgumentParser()
ap.add_argument('--archive', action='store_true')
ap.add_argument('--only', nargs='*')
args = ap.parse_args()

sessions = sorted(os.path.basename(os.path.dirname(p)) for p in glob.glob(os.path.join(here, 'pass1', '*', '*_whisker.json')))
if args.only:
    sessions = [s for s in sessions if s in args.only]
done = []
for s in sessions:
    animal = s.split('_')[0]
    wr = os.path.join(WORK, s, 'whisker_results')
    v1 = os.path.join(WORK, s, 'whisker_results_v1')
    if os.path.isdir(wr) and not os.path.isdir(v1):
        shutil.copytree(wr, v1)
    video = glob.glob(os.path.join(WORK, s, 'camera', '*.mp4'))[0]
    csv = os.path.join(here, 'pass1', s, s + '_whisker.csv')
    r = subprocess.run([PY, '-u', FM, '--reclassify', video, csv, '--out-dir', wr, '--name', s], capture_output=True, text=True, env=env)
    if r.returncode != 0:
        print(s, 'FAILED', r.stderr[-500:], flush=True)
        continue
    print(s, r.stdout.strip().splitlines()[-1], flush=True)
    qd = os.path.join(WORK, 'quality_reports', s)
    if os.path.isdir(qd):
        for suf in ('_whisker.json', '_whisker_qc.png'):
            shutil.copy2(os.path.join(wr, s + suf), os.path.join(qd, s + suf))
    if args.archive:
        ad = os.path.join(PRE, animal, s, 'whisker_results')
        if os.path.isdir(ad):
            av1 = os.path.join(PRE, animal, s, 'whisker_results_v1')
            if not os.path.isdir(av1):
                shutil.copytree(ad, av1)
            for suf in ('_whisker.csv', '_whisker.json', '_whisker_qc.png'):
                shutil.copy2(os.path.join(wr, s + suf), os.path.join(ad, s + suf))
            print('   archive updated', ad, flush=True)
    done.append(s)
print('done', len(done), 'of', len(sessions))
