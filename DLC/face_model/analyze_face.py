"""Run the face model on a raw camera video and write the nose-tip (and whisker-point) tracks in raw-video pixels.

Run via analyze_face.bat:
  analyze_face.bat "<raw video.mp4>" --out-dir "<folder>"                  (the whole video, ~1x real time)
  analyze_face.bat "<raw video.mp4>" --out-dir "<folder>" --sample 5       (5 clips of 20 s spread over the video)
The model sees the fixed 800x700 face window of the 1856x1020 frame (same window the labelling frames were cut from).
Writes <out-dir>\<video name>_face.csv with columns frame, sample, source_frame, <part>_x, <part>_y, <part>_likelihood
(x, y in the FULL frame; DLC itself reports them relative to the crop, the crop offset is added back here). For a whole
video sample is 0 and source_frame = frame; with --sample, frame counts within the clip and source_frame is the
(keyframe-aligned, so +-0.3 s) position of that clip frame in the video. DLC's own csv/h5 stay in <out-dir>\dlc_raw.
The sampled track is enough to place the whisking regions (face_motion.py), which only need where the nose sits.
"""
import argparse
import glob
import os
import subprocess
import sys
import time

TRAINING_DATA = r'C:\Users\rbondarenko\DLC\Training_Data'
FACE_WINDOW = (40, 840, 230, 930)  # x1, x2, y1, y2 of the face crop (must match make_face_project.py)

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('video')
ap.add_argument('--out-dir', required=True)
ap.add_argument('--config', help='project config.yaml (default: the face project in Training_Data)')
ap.add_argument('--sample', type=int, default=0, metavar='N',
                help='track N short clips spread evenly over the video instead of all of it')
ap.add_argument('--clip-s', type=float, default=20, help='length of each sampled clip in seconds (default: %(default)s)')
args = ap.parse_args()

if args.config:
    config = args.config
else:
    found = sorted(glob.glob(os.path.join(TRAINING_DATA, 'face-*', 'config.yaml')))
    if not found:
        sys.exit('ERROR: no face project in %s' % TRAINING_DATA)
    config = found[-1]
if not os.path.isfile(args.video):
    sys.exit('ERROR: video not found: %s' % args.video)

os.environ.setdefault('CUDA_CACHE_MAXSIZE', '4294967296')  # keep TensorFlow's compiled GPU kernels

import cv2
import pandas as pd
import deeplabcut as dlc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import face_motion as fm  # find_tool

raw_dir = os.path.join(args.out_dir, 'dlc_raw')
os.makedirs(raw_dir, exist_ok=True)
name, ext = os.path.splitext(os.path.basename(args.video))


def cut_clips(video, n, clip_s):
    """Stream-copy n clips of clip_s seconds, centred in n equal slices of the video; [(path, start frame), ...]."""
    cap = cv2.VideoCapture(video)
    frames, fps = int(cap.get(7)), cap.get(5) or fm.FPS
    cap.release()
    dur = frames / fps
    out = []
    for k in range(n):
        start = min(max(dur * (k + 0.5) / n - clip_s / 2, 0), max(dur - clip_s, 0))
        path = os.path.join(raw_dir, '%s_sample%d%s' % (name, k, ext))
        subprocess.run([fm.find_tool('ffmpeg'), '-v', 'error', '-y', '-ss', '%.3f' % start, '-i', video, '-t',
                        '%.3f' % clip_s, '-c', 'copy', path], check=True)
        out.append((path, int(round(start * fps))))
    return out


def read_track(video):
    csvs = sorted(glob.glob(os.path.join(raw_dir, os.path.splitext(os.path.basename(video))[0] + 'DLC_resnet50_face*.csv')))
    if not csvs:
        sys.exit('ERROR: no DLC csv found for %s in %s' % (video, raw_dir))
    return pd.read_csv(csvs[-1], header=[1, 2], index_col=0)


t0 = time.time()
if args.sample:
    clips = cut_clips(args.video, args.sample, args.clip_s)
    sources = [(p, s) for p, s in clips]
else:
    sources = [(args.video, 0)]
dlc.analyze_videos(config, [p for p, _ in sources], videotype=ext, shuffle=1, trainingsetindex=0, save_as_csv=True,
                   destfolder=raw_dir, cropping=list(FACE_WINDOW))
print('DLC analysis took %.1f min' % ((time.time() - t0) / 60))

parts = []
for k, (path, start) in enumerate(sources):
    df = read_track(path)
    out = pd.DataFrame(index=df.index)
    out['sample'] = k
    out['source_frame'] = df.index + start
    for part in dict.fromkeys(df.columns.get_level_values(0)):
        out[part + '_x'] = df[part]['x'] + FACE_WINDOW[0]
        out[part + '_y'] = df[part]['y'] + FACE_WINDOW[2]
        out[part + '_likelihood'] = df[part]['likelihood']
    parts.append(out)
out = pd.concat(parts)  # keeps each clip's own frame numbers as the index
out.index.name = 'frame'
dest = os.path.join(args.out_dir, name + '_face.csv')
out.to_csv(dest)
if args.sample:  # the cut clips are only needed for the analysis
    for path, _ in sources:
        os.remove(path)
print('wrote %s (%d frames%s)' % (dest, len(out), ', %d clips' % len(sources) if args.sample else ''))
