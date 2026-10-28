"""Create the DLC project for the face model (nose tip + one whisker) and pick the frames to label.

    python make_face_project.py [--camera-root "Z:\Group Members\Rima\Camera"] [--per-video 8]

Run in the dlc-tf environment (make_face_project.bat does that). Steps:
  1. from every new-setup video (AnimalRB22/23/24 folders), sample random frames inside the face window
     FACE_WINDOW (full resolution, fixed crop that holds the nose and the whisker fan in all sessions)
  2. per video keep --per-video frames: half from the still moments, half from the active moments (frame-to-frame
     change inside the window), each half spread over the appearance clusters (k-means), so the whiskers show up
     in many states and the head position varies
  3. write them as PNGs plus a stand-in video, create the DLC project <Training_Data>\face-RB-<date> with the
     bodyparts below, and put the PNGs in labeled-data\face_labelling_frames, ready for label_frames.bat
frames_provenance.csv in that folder says which session and frame each image came from.
"""
import argparse
import csv
import datetime
import glob
import os
import shutil
import subprocess
import sys

import cv2
import numpy as np

TRAINING_DATA = r'C:\Users\rbondarenko\DLC\Training_Data'
FACE_WINDOW = (40, 840, 230, 930)  # x0, x1, y0, y1 in full-resolution pixels (800 x 700)
BODYPARTS = ['noseTip', 'whiskerBase', 'whiskerMid', 'whiskerTip']
SKELETON = [['noseTip', 'whiskerBase'], ['whiskerBase', 'whiskerMid'], ['whiskerMid', 'whiskerTip']]
CANDIDATES_PER_VIDEO = 60
FOLDER = 'face_labelling_frames'


def find_ffmpeg():
    exe = shutil.which('ffmpeg')
    if exe:
        return exe
    raise FileNotFoundError('ffmpeg not on PATH (activate the dlc-tf environment)')


def pick_from_video(video, per_video, rng):
    """Returns [(frame index, activity, cropped BGR frame)] for the frames kept from this video."""
    from sklearn.cluster import KMeans
    x0, x1, y0, y1 = FACE_WINDOW
    cap = cv2.VideoCapture(video)
    n = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    # keyframe + the frame after it carry codec noise, so compare two frames that sit 3-9 frames after a keyframe
    idx = np.sort(rng.integers(int(0.01 * n), int(0.99 * n) - 12, CANDIDATES_PER_VIDEO))
    idx = idx - idx % 12 + rng.integers(3, 10, len(idx))
    cands = []
    for i in idx:
        cap.set(cv2.CAP_PROP_POS_FRAMES, int(i))
        ok0, f0 = cap.read()
        ok1, f1 = cap.read()
        if not (ok0 and ok1):
            continue
        a = cv2.GaussianBlur(cv2.cvtColor(f0[y0:y1, x0:x1], cv2.COLOR_BGR2GRAY), (0, 0), 1.0).astype(np.float32)
        b = cv2.GaussianBlur(cv2.cvtColor(f1[y0:y1, x0:x1], cv2.COLOR_BGR2GRAY), (0, 0), 1.0).astype(np.float32)
        feat = cv2.resize(a, (50, 44), interpolation=cv2.INTER_AREA).ravel()
        cands.append((int(i), float(np.abs(a - b).mean()), feat, f0[y0:y1, x0:x1].copy()))
    cap.release()
    cands.sort(key=lambda c: c[1])
    half = len(cands) // 2
    kept = []
    for group, k in ((cands[:half], per_video // 2), (cands[half:], per_video - per_video // 2)):
        feats = np.stack([c[2] for c in group])
        km = KMeans(n_clusters=k, n_init=5, random_state=0).fit(feats)
        for c in range(k):
            members = np.flatnonzero(km.labels_ == c)
            best = members[np.argmin(np.linalg.norm(feats[members] - km.cluster_centers_[c], axis=1))]
            kept.append((group[best][0], group[best][1], group[best][3]))
    return sorted(kept)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--camera-root', default=r'Z:\Group Members\Rima\Camera')
    ap.add_argument('--per-video', type=int, default=8)
    args = ap.parse_args()

    projects = glob.glob(os.path.join(TRAINING_DATA, 'face-*'))
    if projects:
        sys.exit('STOP: a face project already exists: %s\nDelete it first if you really want to start over.' % projects[0])
    videos = sorted(glob.glob(os.path.join(args.camera_root, 'AnimalRB2[234]', '*', '*.mp4')))
    if not videos:
        sys.exit('no videos found under ' + args.camera_root)

    stage = os.path.join(TRAINING_DATA, '_face_stage')
    os.makedirs(stage, exist_ok=True)
    rng = np.random.default_rng(0)
    rows, k = [], 0
    for v in videos:
        session = os.path.basename(os.path.dirname(v))
        cap = cv2.VideoCapture(v)
        n = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
        cap.release()
        if n < 3000:
            print('skip %s (%d frames)' % (session, n))
            continue
        print('%s: sampling %d candidate frames' % (session, CANDIDATES_PER_VIDEO), flush=True)
        for idx, act, img in pick_from_video(v, args.per_video, rng):
            k += 1
            cv2.imwrite(os.path.join(stage, 'img%04d.png' % k), img)
            rows.append(['img%04d.png' % k, session, idx, round(act, 3)])
    print('%d frames selected' % k)

    video = os.path.join(stage, FOLDER + '.mp4')
    subprocess.run([find_ffmpeg(), '-y', '-v', 'error', '-framerate', '5', '-i', os.path.join(stage, 'img%04d.png'),
                    '-c:v', 'libx264', '-crf', '4', '-pix_fmt', 'yuv420p', video], check=True)

    import deeplabcut as dlc
    from deeplabcut.utils import auxiliaryfunctions
    cfg = dlc.create_new_project('face', 'RB', [video], working_directory=TRAINING_DATA, copy_videos=True)
    auxiliaryfunctions.edit_config(cfg, {'bodyparts': BODYPARTS, 'skeleton': SKELETON, 'TrainingFraction': [0.9],
                                         'numframes2pick': k, 'cropping': False})
    folder = os.path.join(os.path.dirname(cfg), 'labeled-data', FOLDER)
    os.makedirs(folder, exist_ok=True)
    for name in sorted(os.listdir(stage)):
        if name.endswith('.png'):
            shutil.copy2(os.path.join(stage, name), os.path.join(folder, name))
    with open(os.path.join(folder, 'frames_provenance.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['image', 'session', 'frame', 'activity'])
        w.writerows(rows)
    shutil.rmtree(stage)
    print('\nDONE: project %s\n%d frames to label in %s' % (cfg, k, folder))


if __name__ == '__main__':
    main()
