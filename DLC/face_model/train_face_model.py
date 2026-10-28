"""Train the face model (noseTip + one whisker) on the hand labels, then evaluate it on the held-out frames.

Run via train_face_model.bat. What it does:
  1. checks the labels (CollectedData_RB.h5) and holds out a seeded 10% of the labelled frames as TEST frames
  2. create_training_dataset (resnet_50 from the ImageNet weights shipped with DLC), iteration 0
  3. pose_cfg.yaml: learning-rate schedule scaled to --maxiters (DLC's default one is written for 1M iterations)
  4. train_network, then evaluate_network (pixel error on train and test frames)
The held-out frame names are written to <project>\dlc-models\iteration-0\heldout_test_frames.txt.
"""
import argparse
import glob
import os
import pickle
import sys
import time

import numpy as np

TRAINING_DATA = r'C:\Users\rbondarenko\DLC\Training_Data'
FOLDER = 'face_labelling_frames'

ap = argparse.ArgumentParser()
ap.add_argument('--config', help='project config.yaml (default: the face project in Training_Data)')
ap.add_argument('--maxiters', type=int, default=40000)
ap.add_argument('--saveiters', type=int, default=5000)
ap.add_argument('--displayiters', type=int, default=500)
ap.add_argument('--test-fraction', type=float, default=0.1)
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--skip-train', action='store_true', help='only evaluate the existing model')
args = ap.parse_args()

if args.config:
    config = args.config
else:
    found = sorted(glob.glob(os.path.join(TRAINING_DATA, 'face-*', 'config.yaml')))
    if not found:
        sys.exit('ERROR: no face project in %s; run make_face_project.bat first' % TRAINING_DATA)
    config = found[-1]
project = os.path.dirname(config)

os.environ.setdefault('CUDA_CACHE_MAXSIZE', '4294967296')  # keep TensorFlow's compiled GPU kernels

import pandas as pd
import yaml
import deeplabcut as dlc
from deeplabcut.utils import auxiliaryfunctions
from deeplabcut.generate_training_dataset.trainingsetmanipulation import merge_annotateddatasets

cfg = auxiliaryfunctions.read_config(config)
scorer = cfg['scorer']
parts = cfg['bodyparts']
frac = cfg['TrainingFraction'][0]
model_dir = os.path.join(project, 'dlc-models', 'iteration-%d' % cfg['iteration'])
train_dir = os.path.join(model_dir, '%s%s-trainset%dshuffle1' % (cfg['Task'], cfg['date'], round(frac * 100)), 'train')

if not args.skip_train:
    if os.path.isdir(model_dir):
        sys.exit('STOP: %s already exists. Delete it (and training-datasets\\iteration-%d) to train again.'
                 % (model_dir, cfg['iteration']))

    lab = pd.read_hdf(os.path.join(project, 'labeled-data', FOLDER, 'CollectedData_%s.h5' % scorer))
    print('labelled frames: %d' % len(lab))
    for p in parts:
        print('  %-12s %3d' % (p, int(lab[(scorer, p, 'x')].notna().sum())))

    tsf = os.path.join(project, auxiliaryfunctions.get_training_set_folder(cfg))
    os.makedirs(tsf, exist_ok=True)
    data = merge_annotateddatasets(cfg, tsf)[scorer]  # same (sorted) order create_training_dataset will use
    names = [i[-1] if isinstance(i, tuple) else os.path.basename(str(i).replace('\\', '/')) for i in data.index]
    n = len(data)
    n_test = int(round(args.test_fraction * n))
    if round((n - n_test) / n, 2) != frac:
        sys.exit('ERROR: %d test of %d frames gives train fraction %.3f, config TrainingFraction is %s'
                 % (n_test, n, (n - n_test) / n, frac))
    rng = np.random.default_rng(args.seed)
    test_idx = np.sort(rng.choice(n, n_test, replace=False))
    train_idx = np.setdiff1d(np.arange(n), test_idx)
    print('\n%d frames -> %d train / %d test' % (n, len(train_idx), len(test_idx)))

    dlc.create_training_dataset(config, num_shuffles=1, Shuffles=[1],
                                trainIndices=[train_idx.tolist()], testIndices=[test_idx.tolist()])

    doc = [f for f in os.listdir(tsf) if f.startswith('Documentation_data') and f.endswith('.pickle')][0]
    with open(os.path.join(tsf, doc), 'rb') as fh:
        _, tr, te, _ = pickle.load(fh)
    assert sorted(map(int, te)) == sorted(map(int, test_idx)), 'test indices in the dataset differ from the requested ones'
    with open(os.path.join(model_dir, 'heldout_test_frames.txt'), 'w') as fh:
        fh.write('\n'.join(sorted(names[i] for i in test_idx)) + '\n')

    pc = os.path.join(train_dir, 'pose_cfg.yaml')
    with open(pc, encoding='utf-8') as fh:
        p = yaml.safe_load(fh)
    m = args.maxiters
    p['multi_step'] = [[0.005, int(m * 0.125)], [0.02, int(m * 0.625)], [0.002, int(m * 0.875)], [0.001, m]]
    with open(pc, 'w', encoding='utf-8') as fh:
        yaml.safe_dump(p, fh, default_flow_style=False)
    print('pose_cfg: net_type=%s batch_size=%s global_scale=%s multi_step=%s'
          % (p['net_type'], p.get('batch_size'), p.get('global_scale'), p['multi_step']))

    import tensorflow as tf
    print('GPUs visible to TensorFlow:', tf.config.list_physical_devices('GPU'))
    print('training %d iterations (saving every %d); the first iterations are slow the very first time '
          '(one-off GPU kernel compile)\n' % (m, args.saveiters))
    t0 = time.time()
    dlc.train_network(config, shuffle=1, trainingsetindex=0, displayiters=args.displayiters,
                      saveiters=args.saveiters, maxiters=m, max_snapshots_to_keep=3, allow_growth=True)
    print('\nTRAINED: %d iterations in %.1f min' % (m, (time.time() - t0) / 60))

print('\nevaluating ...')
dlc.evaluate_network(config, Shuffles=[1], trainingsetindex=0, plotting=False, comparisonbodyparts='all')
res = glob.glob(os.path.join(project, 'evaluation-results', 'iteration-%d' % cfg['iteration'], '*', '*-results.csv'))
for f in res:
    print(f)
    print(pd.read_csv(f).to_string(index=False))
print('\nDONE')
