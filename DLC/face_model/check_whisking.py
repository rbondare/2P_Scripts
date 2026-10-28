"""Blind check of the whisking on/off call: short unlabelled clips, you judge them, then compare with the code.

  1. make the clips (random 5 s stretches that the code calls whisking, and others it calls rest, from every session):
       python check_whisking.py make --out-dir <check folder> "<session1.mp4>|<session1_whisker.csv>" ["<session2.mp4>|<...csv>" ...]
     -> <check folder>\\clip_01.mp4 ... (half speed, the nose-anchored regions drawn on, no label), your_answers.csv
        (empty, to fill in) and answer_key.csv (what the code said; do not open it first).
  2. watch the clips and fill in your_answers.csv: whisking = y (whiskers moving), n (still), or ? (cannot tell)
  3. compare:   python check_whisking.py score --out-dir <check folder>
"""
import argparse
import csv
import json
import os
import subprocess
import sys

import cv2
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import face_motion as fm

CLIP_S = 5.0
PLAY_FPS = 25  # half speed
WIN = (-330, -140, 360, 340)  # dx0, dy0, dx1, dy1 around the nose anchor
TRACE_H = 90
TRACE_S = 4.0  # trace shows +-this many seconds
COLORS = {'fan_front': (0, 255, 0), 'fan_down': (255, 0, 255), 'pad': (0, 140, 255)}
SKIP_START_S = 30


def load(csv_path):
    idx, on, art = [], [], []
    with open(csv_path, newline='') as f:
        for row in csv.DictReader(f):
            idx.append(float(row['whisk_index']))
            on.append(row['whisk_on'] == '1')
            art.append(row['artifact'] == '1')
    return np.array(idx), np.array(on), np.array(art)


def pick_window(state, valid, rng, n):
    """Random start frame of an n-frame window where state is True throughout (and valid), or None."""
    ok = state & valid
    ok[:int(SKIP_START_S * fm.FPS)] = False
    cand = [(a, b) for a, b in fm.runs(ok) if b - a >= n]
    if not cand:
        return None
    a, b = cand[rng.integers(len(cand))]
    return int(rng.integers(a, b - n + 1))


def render(video, start, anchor, index, path, label, show_trace=False):
    ax, ay = anchor
    n = int(CLIP_S * fm.FPS)
    w, h = WIN[2] - WIN[0], WIN[3] - WIN[1]
    _, width, height = fm.video_frames(video)
    x0 = int(np.clip(ax + WIN[0], 0, width - w))  # the window stays inside the frame when the nose is near an edge
    y0 = int(np.clip(ay + WIN[1], 0, height - h))
    cmd = [fm.find_tool('ffmpeg'), '-v', 'error', '-ss', '%.3f' % (start / fm.FPS), '-i', video, '-frames:v', str(n),
           '-vf', 'crop=%d:%d:%d:%d' % (w, h, x0, y0), '-f', 'rawvideo', '-pix_fmt', 'bgr24', 'pipe:1']
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=w * h * 3 * 8)
    out = cv2.VideoWriter(path, cv2.VideoWriter_fourcc(*'mp4v'), PLAY_FPS, (w, h + (TRACE_H if show_trace else 0)))
    half = int(TRACE_S * fm.FPS)
    top = max(float(np.percentile(index[start:start + n], 99)), 3.0) * 1.1
    k = 0
    while k < n:
        buf = proc.stdout.read(w * h * 3)
        if len(buf) < w * h * 3:
            break
        img = np.frombuffer(buf, np.uint8).reshape(h, w, 3).copy()
        for name, col in COLORS.items():
            bx0, by0, bx1, by1 = fm.BOXES[name]
            cv2.rectangle(img, (ax + bx0 - x0, ay + by0 - y0), (ax + bx1 - x0, ay + by1 - y0), col, 1)
        cv2.putText(img, label, (8, 22), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (255, 255, 255), 2, cv2.LINE_AA)
        if show_trace:
            strip = np.zeros((TRACE_H, w, 3), np.uint8)
            t0 = start + k - half
            seg = index[max(t0, 0):start + k + half]
            xs = (np.arange(max(t0, 0), max(t0, 0) + len(seg)) - t0) / (2 * half) * (w - 1)
            ys = TRACE_H - 4 - np.clip(seg / top, 0, 1) * (TRACE_H - 10)
            cv2.polylines(strip, [np.stack([xs, ys], 1).astype(np.int32)], False, (255, 200, 120), 1, cv2.LINE_AA)
            yt = int(TRACE_H - 4 - min(fm.THRESH_X / top, 1) * (TRACE_H - 10))
            cv2.line(strip, (0, yt), (w - 1, yt), (0, 140, 255), 1)
            cv2.line(strip, (w // 2, 0), (w // 2, TRACE_H - 1), (255, 255, 255), 1)
            img = np.vstack([img, strip])
        out.write(img)
        k += 1
    out.release()
    proc.stdout.close()
    proc.wait()
    return k == n


def make(args):
    os.makedirs(args.out_dir, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    n = int(CLIP_S * fm.FPS)
    items = []
    for pair in args.sessions:
        video, csv_path = pair.split('|')
        idx, on, art = load(csv_path)
        info = json.load(open(csv_path.replace('_whisker.csv', '_whisker.json')))
        anchor = tuple(info['anchor_median_xy_px'])
        for state_name, state in (('whisking', on), ('rest', ~on)):
            for _ in range(args.per_state):
                s = pick_window(state, ~art, rng, n)
                if s is not None:
                    items.append((video, info['session'], anchor, idx, state_name, s))
    order = rng.permutation(len(items))
    key = []
    for rank, i in enumerate(order, 1):
        video, session, anchor, idx, state_name, s = items[i]
        name = 'clip_%02d' % rank
        ok = render(video, s, anchor, idx, os.path.join(args.out_dir, name + '.mp4'), name, args.show_trace)
        key.append(dict(clip=name, code_says=state_name, session=session, start_s=round(s / fm.FPS, 1),
                        mean_index=round(float(idx[s:s + n].mean()), 2), rendered=ok))
        print(name, 'done' if ok else 'SHORT READ', flush=True)
    with open(os.path.join(args.out_dir, 'answer_key.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(key[0]))
        w.writeheader()
        w.writerows(key)
    with open(os.path.join(args.out_dir, 'your_answers.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['clip', 'whisking (y/n/?)'])
        for k in key:
            w.writerow([k['clip'], ''])
    print('wrote %d clips to %s' % (len(key), args.out_dir))


def score(args):
    key = {r['clip']: r for r in csv.DictReader(open(os.path.join(args.out_dir, 'answer_key.csv')))}
    ans = {}
    for r in csv.reader(open(os.path.join(args.out_dir, 'your_answers.csv'))):
        # Excel with a ';' list separator saves 'clip_01,;y', so the answer is whatever is left once ';' and spaces go
        cells = [c.strip(' ;\t"').lower() for c in r[1:]]
        cells = [c for c in cells if c]
        if r and r[0] in key and cells:
            ans[r[0]] = cells[0][0]
    tab = {}
    for c, a in ans.items():
        tab[(key[c]['code_says'], a)] = tab.get((key[c]['code_says'], a), 0) + 1
    print('answered %d of %d clips\n' % (len(ans), len(key)))
    print('%-22s %8s %8s %8s' % ('code said / you said', 'whisking', 'still', 'cannot tell'))
    for s in ('whisking', 'rest'):
        print('%-22s %8d %8d %8d' % (s, tab.get((s, 'y'), 0), tab.get((s, 'n'), 0), tab.get((s, '?'), 0)))
    decided = [(key[c]['code_says'], a) for c, a in ans.items() if a in 'yn']
    agree = sum((s == 'whisking') == (a == 'y') for s, a in decided)
    print('\nagreement on the clips you could judge: %d of %d' % (agree, len(decided)))
    for c, a in sorted(ans.items()):
        if a in 'yn' and (key[c]['code_says'] == 'whisking') != (a == 'y'):
            print('  disagree: %s  code=%s you=%s  (%s at %.0f s, mean index %s)' % (
                c, key[c]['code_says'], a, key[c]['session'], float(key[c]['start_s']), key[c]['mean_index']))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    m = sub.add_parser('make')
    m.add_argument('--out-dir', required=True)
    m.add_argument('--per-state', type=int, default=1, help='clips per session for each of whisking / rest')
    m.add_argument('--seed', type=int, default=0)
    m.add_argument('--show-trace', action='store_true',
                   help='draw the whisking index under the clip (shows the code\'s answer, so not blind)')
    m.add_argument('sessions', nargs='+', help='"<raw video.mp4>|<session_whisker.csv>" (one argument per session)')
    m.set_defaults(func=make)
    s = sub.add_parser('score')
    s.add_argument('--out-dir', required=True)
    s.set_defaults(func=score)
    args = ap.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
