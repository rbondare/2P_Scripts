"""Open the face frames in the DLC labelling GUI (napari-deeplabcut).

Run via label_frames.bat. Labels are saved by the GUI (Ctrl+S) into
  <project>\labeled-data\face_labelling_frames\CollectedData_RB.h5   (+ .csv)
You can close it and come back any time; saved labels are loaded again.
"""
import argparse
import glob
import os
import sys

TRAINING_DATA = r'C:\Users\rbondarenko\DLC\Training_Data'
FOLDER = 'face_labelling_frames'

ap = argparse.ArgumentParser()
ap.add_argument('--config', help='project config.yaml (default: the face project in Training_Data)')
args = ap.parse_args()

if args.config:
    config = args.config
else:
    found = sorted(glob.glob(os.path.join(TRAINING_DATA, 'face-*', 'config.yaml')))
    if not found:
        sys.exit('ERROR: no face project in %s; run make_face_project.bat first' % TRAINING_DATA)
    config = found[-1]
if not os.path.isdir(os.path.join(os.path.dirname(config), 'labeled-data', FOLDER)):
    sys.exit('ERROR: no labeled-data\%s folder in the project' % FOLDER)

import deeplabcut as dlc
import napari

# image_folder must be given: with only the config, DLC opens the FIRST folder in labeled-data
dlc.label_frames(config, image_folder=FOLDER)
napari.run()  # label_frames only creates the window; this keeps it open until you close it
