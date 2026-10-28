"""Run the body + eye DeepLabCut pipeline on recordings that live on the drive.

    python run_dlc_to_drive.py <session_folder> [<session_folder> ...]
    python run_dlc_to_drive.py --scan "Z:\\Group Members\\Rima\\Camera"
    python run_dlc_to_drive.py --scan <root> --dry-run          (only list what would be processed)

(run_dlc.bat takes the same arguments and activates the dlc-tf environment first.)

A session folder is one recording (e.g. ...\\Camera\\AnimalRB22\\AnimalRB22_261023_1427) holding the raw
.mp4/.avi, either directly or in a camera\\ subfolder. For each one:
  1. the raw video is copied to the local work folder
  2. RB_DLC_analysis.py runs on the local copy (body model, eye model, eye ellipse fit)
  3. whisking on/off, amplitude and face motion from the raw video into whisker_results\\ (skip with
     --no-whisker): the face model (face_model\\analyze_face.py) tracks the nose in 5 short clips, then
     face_model\\face_motion.py measures motion energy in nose-anchored regions (see its docstring)
  4. the quality-metric report (summarize_eye_run.py) and the whisker QC figure are written to
     <work>\\quality_reports\\<session>; no labelled video is made
  5. body_results\\, eye_results\\, whisker_results\\ and success.log are copied into the session folder on
     the drive, laid out like the earlier sessions (e.g. AnimalRB18_260312_1803) plus whisker_results\\.
     Every file is size-checked, and success.log goes last so the MATLAB importer never picks up a
     half-copied session.

The raw video on the drive is not moved or changed, and nothing local is deleted: once the copy is
verified, delete <work>\\<session> yourself. Sessions that already have a success.log on the drive are
skipped, and a session whose local run already finished (success.log in the work folder) goes straight
to the copy, so an interrupted batch can simply be started again.
"""
import argparse
import datetime
import json
import os
import re
import shutil
import subprocess
import sys
import time

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_WORK_DIR = r'C:\Users\rbondarenko\DLC\work'
DEFAULT_CONFIG = os.path.join(HERE, 'experiment_config.yaml')
VIDEO_EXT = ('.mp4', '.avi')
TIFF_EXT = ('.tif', '.tiff')
SESSION_RE = re.compile(r'^.+_\d{6}_\d{4}$')  # AnimalRB22_261023_1427
RESULT_DIRS = ('body_results', 'eye_results')
WHISKER_DIR = 'whisker_results'
FACE_SAMPLES = 5  # clips the face model tracks per session to place the whisking regions on the nose
FACE_CLIP_S = 20  # seconds per clip
ANCHOR_RANGE_WARN_PX = 30  # warn when the nose anchor moves more than this over a session
COPY_BUFFER = 16 * 1024 * 1024

_log = None


def say(text='', end='\n'):
    sys.stdout.write(text + end)
    sys.stdout.flush()
    if _log:
        _log.write(text + end)
        _log.flush()


def stamp():
    return datetime.datetime.now().strftime('%H:%M:%S')


def list_files(folder, exts):
    if not os.path.isdir(folder):
        return []
    return sorted(f for f in os.listdir(folder)
                  if f.lower().endswith(exts) and os.path.isfile(os.path.join(folder, f)))


def find_raw_videos(session):
    """(folder, file names) of the raw videos in session\\camera, or in session itself."""
    for folder in (os.path.join(session, 'camera'), session):
        names = list_files(folder, VIDEO_EXT)
        if names:
            return folder, names
    return None, []


def has_tiffs(session):
    return any(list_files(f, TIFF_EXT) for f in (os.path.join(session, 'camera'), session))


def scan(root):
    """Recording folders under root, either directly in it or inside an animal folder (root\\Animal\\Recording)."""
    found = []
    for name in sorted(os.listdir(root)):
        path = os.path.join(root, name)
        if not os.path.isdir(path):
            continue
        if SESSION_RE.match(name):
            found.append(path)
            continue
        for sub in sorted(os.listdir(path)):
            subpath = os.path.join(path, sub)
            if SESSION_RE.match(sub) and os.path.isdir(subpath):
                found.append(subpath)
    return found


def plan(sessions):
    """Split sessions into (ready, skipped) with a reason for every skipped one."""
    ready, skipped, seen = [], [], {}
    for s in sessions:
        name = os.path.basename(os.path.normpath(s))
        if name in seen:
            skipped.append((s, 'same recording name as ' + seen[name]))
        elif os.path.isfile(os.path.join(s, 'success.log')):
            skipped.append((s, 'already has success.log on the drive'))
        elif not find_raw_videos(s)[1]:
            skipped.append((s, 'TIFF frames, not handled here' if has_tiffs(s) else 'no .mp4/.avi found'))
        else:
            seen[name] = s
            ready.append(s)
    return ready, skipped


def gb(n):
    return '%.2f GB' % (n / 1e9)


def copy_file(src, dst):
    """Copy via a .part file so an interrupted copy never looks complete; the size is checked afterwards."""
    size = os.path.getsize(src)
    if os.path.isfile(dst) and os.path.getsize(dst) == size:
        return size
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    part = dst + '.part'
    with open(src, 'rb') as fin, open(part, 'wb') as fout:
        shutil.copyfileobj(fin, fout, COPY_BUFFER)
    shutil.copystat(src, part)
    os.replace(part, dst)
    if os.path.getsize(dst) != size:
        raise IOError('size mismatch after copying %s -> %s' % (src, dst))
    return size


def copy_tree_files(src_dir, dst_dir):
    """Copy every file of src_dir into dst_dir; returns (n files, bytes)."""
    n, total = 0, 0
    for name in sorted(os.listdir(src_dir)):
        src = os.path.join(src_dir, name)
        if os.path.isfile(src):
            total += copy_file(src, os.path.join(dst_dir, name))
            n += 1
    return n, total


def run_streamed(cmd, keep=False):
    """Run cmd, echo its output live (progress bars included) and to the log; returns (exit code, text)."""
    env = dict(os.environ, PYTHONUTF8='1', PYTHONIOENCODING='utf-8')
    env.setdefault('CUDA_CACHE_MAXSIZE', '4294967296')  # keep TensorFlow's compiled GPU kernels
    proc = subprocess.Popen(cmd, cwd=HERE, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    text = []
    while True:
        chunk = proc.stdout.read1(4096)
        if not chunk:
            break
        piece = chunk.decode('utf-8', errors='replace').replace('\r\n', '\n')  # text-mode output adds \r again
        say(piece, end='')
        if keep:
            text.append(piece)
    return proc.wait(), ''.join(text)


def whisking(local, videos, session, whisker_dir):
    """Whisking on/off, amplitude and face motion for every raw video of the local session copy.
    The face model finds the nose in a few short clips (face_samples\\), face_motion.py places its regions on it and
    writes whisker_dir\\<out name>_whisker.csv/.json/_whisker_qc.png, where <out name> is the session name (the raw
    video is named after the camera, not the session; a session with several videos gets <session>_<video>).
    Returns the list of <out name>s, empty if anything failed."""
    face_dir = os.path.join(local, 'face_samples')
    names = []
    for v in videos:
        video, stem = os.path.join(local, 'camera', v), os.path.splitext(v)[0]
        out_name = session if len(videos) == 1 else '%s_%s' % (session, stem)
        face_csv = os.path.join(face_dir, stem + '_face.csv')
        if not os.path.isfile(face_csv):
            code, _ = run_streamed([sys.executable, '-u', os.path.join(HERE, 'face_model', 'analyze_face.py'), video,
                                    '--out-dir', face_dir, '--sample', str(FACE_SAMPLES), '--clip-s', str(FACE_CLIP_S)])
            if code != 0 or not os.path.isfile(face_csv):
                say('        face model failed on %s (exit code %s)' % (v, code))
                return []
        code, _ = run_streamed([sys.executable, '-u', os.path.join(HERE, 'face_model', 'face_motion.py'), video, face_csv,
                                '--out-dir', whisker_dir, '--name', out_name])
        if code != 0 or not os.path.isfile(os.path.join(whisker_dir, out_name + '_whisker.csv')):
            say('        face_motion.py failed on %s (exit code %s)' % (v, code))
            return []
        with open(os.path.join(whisker_dir, out_name + '_whisker.json')) as f:
            moved = max(json.load(f)['anchor_range_px'])
        if moved > ANCHOR_RANGE_WARN_PX:
            say('        WARNING: the nose position moved %d px over the session; check the regions in %s_whisker_qc.png'
                % (moved, out_name))
        names.append(out_name)
    return names


def process(session, args, options):
    name = os.path.basename(os.path.normpath(session))
    src_folder, videos = find_raw_videos(session)
    local = os.path.join(args.work_dir, name)
    qc_dir = os.path.join(args.qc_dir, name)
    t0 = time.time()

    say('[%s] 1/5 copy raw video to %s' % (stamp(), local))
    for v in videos:
        size = copy_file(os.path.join(src_folder, v), os.path.join(local, 'camera', v))
        say('        %s (%s)' % (v, gb(size)))

    if os.path.isfile(os.path.join(local, 'success.log')):
        say('[%s] 2/5 DLC already finished locally, skipping' % stamp())
    else:
        say('[%s] 2/5 run DLC' % stamp())
        # leftovers of an interrupted run in the work folder would be mixed into the new results
        for d in RESULT_DIRS:
            shutil.rmtree(os.path.join(local, d), ignore_errors=True)
        report = os.path.join(local, 'report.log')
        if os.path.isfile(report):
            os.remove(report)
        code, _ = run_streamed([sys.executable, '-u', os.path.join(HERE, 'RB_DLC_analysis.py'), args.config, local])
        if code != 0 or not os.path.isfile(os.path.join(local, 'success.log')):
            reason = open(report).read().strip() if os.path.isfile(report) else ''
            reason = reason or 'it stopped without a report, see the traceback above'
            raise RuntimeError('DLC pipeline did not finish (exit code %s): %s' % (code, reason))
    for d, suffix in (('body_results', '_body.csv'), ('eye_results', '_eye.csv')):
        if not os.path.isfile(os.path.join(local, d, name + suffix)):
            raise RuntimeError('missing %s\\%s%s after the DLC run' % (d, name, suffix))

    whisker_dir = os.path.join(local, WHISKER_DIR)
    whisker_names = []
    if args.no_whisker:
        say('[%s] 3/5 whisking skipped (--no-whisker)' % stamp())
    else:
        say('[%s] 3/5 whisking on/off, amplitude and face motion' % stamp())
        whisker_names = whisking(local, videos, name, whisker_dir)
        if not whisker_names:  # the DLC results are still worth copying
            say('        WARNING: whisking analysis failed; continuing without it')
    have_whisker = bool(whisker_names)

    say('[%s] 4/5 quality metrics -> %s' % (stamp(), qc_dir))
    code, text = run_streamed([sys.executable, '-u', os.path.join(HERE, 'summarize_eye_run.py'), local, qc_dir,
                               '--thr', str(options['eye_likelihood_threshold']),
                               '--fps', str(options['frames_per_second']),
                               '--batch', str(options['batch_size'])], keep=True)
    if code == 0:
        with open(os.path.join(qc_dir, 'summary.txt'), 'w', encoding='utf-8') as f:
            f.write(text)
    else:
        say('        WARNING: quality report failed (exit code %s); continuing with the copy' % code)
    if have_whisker:
        os.makedirs(qc_dir, exist_ok=True)
        for wname in whisker_names:
            for suffix in ('_whisker_qc.png', '_whisker.json'):
                shutil.copy2(os.path.join(whisker_dir, wname + suffix), os.path.join(qc_dir, wname + suffix))

    say('[%s] 5/5 copy results to %s' % (stamp(), session))
    n_files, n_bytes = 0, 0
    for d in RESULT_DIRS + ((WHISKER_DIR,) if have_whisker else ()):
        n, b = copy_tree_files(os.path.join(local, d), os.path.join(session, d))
        n_files += n
        n_bytes += b
        extra = sorted(set(os.listdir(os.path.join(session, d))) - set(os.listdir(os.path.join(local, d))))
        if extra:
            say('        NOTE: %s on the drive also holds files this run did not write: %s' % (d, ', '.join(extra)))
    copy_file(os.path.join(local, 'success.log'), os.path.join(session, 'success.log'))
    say('        %d files, %s copied and size-checked' % (n_files + 1, gb(n_bytes)))
    say('[%s] done: %s (%.1f min)' % (stamp(), name, (time.time() - t0) / 60))


def main():
    global _log
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('sessions', nargs='*', help='session folder(s) on the drive')
    ap.add_argument('--scan', metavar='ROOT', help='process every unprocessed recording under ROOT '
                    '(directly, or inside animal folders)')
    ap.add_argument('--work-dir', default=DEFAULT_WORK_DIR, help='local scratch folder for the raw video copies and '
                    'DLC output (default: %(default)s)')
    ap.add_argument('--qc-dir', help='where the quality reports go (default: <work-dir>\\quality_reports)')
    ap.add_argument('--config', default=DEFAULT_CONFIG, help='pipeline config (default: %(default)s)')
    ap.add_argument('--no-whisker', action='store_true', help='skip the whisking / face-motion step (face model on '
                    'a few clips + face_model\\face_motion.py); it is on by default')
    ap.add_argument('--dry-run', action='store_true', help='list what would be processed and stop')
    args = ap.parse_args()
    args.qc_dir = args.qc_dir or os.path.join(args.work_dir, 'quality_reports')
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')

    sessions = [os.path.abspath(s) for s in args.sessions]
    if args.scan:
        sessions += scan(args.scan)
    if not sessions:
        ap.error('give at least one session folder, or --scan ROOT')
    missing = [s for s in sessions if not os.path.isdir(s)]
    if missing:
        ap.error('not a folder: ' + ', '.join(missing))
    with open(args.config) as f:
        options = yaml.load(f, Loader=yaml.FullLoader)

    ready, skipped = plan(sessions)
    for s, why in skipped:
        say('skip   %s  (%s)' % (s, why))
    for s in ready:
        folder, videos = find_raw_videos(s)
        say('run    %s  (%s)' % (s, gb(sum(os.path.getsize(os.path.join(folder, v)) for v in videos))))
    say('%d to process, %d skipped' % (len(ready), len(skipped)))
    if args.dry_run or not ready:
        return 0

    os.makedirs(os.path.join(args.work_dir, 'logs'), exist_ok=True)
    log_path = os.path.join(args.work_dir, 'logs', 'dlc_%s.log' % datetime.datetime.now().strftime('%Y%m%d_%H%M%S'))
    _log = open(log_path, 'a', encoding='utf-8')
    say('log: ' + log_path)

    failed = []
    for k, s in enumerate(ready, 1):
        say('\n===== %d/%d  %s =====' % (k, len(ready), os.path.basename(s)))
        try:
            process(s, args, options)
        except Exception as e:  # keep going with the other sessions
            say('FAILED %s: %s' % (s, e))
            failed.append(s)

    say('\n%d of %d sessions done, %d failed' % (len(ready) - len(failed), len(ready), len(failed)))
    for s in failed:
        say('  failed: ' + s)
    say('quality reports: ' + args.qc_dir)
    say('local copies are still in %s; delete them once you have checked the results on the drive' % args.work_dir)
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
