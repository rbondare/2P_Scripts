"""Step 1: per-session inventory (header geometry, stage, z, power, ROI counts) -> plots/fov_matching/fov_inventory.csv

Header values come from the FIRST TIFF of each session (tifffile.scanimage_metadata); values are identical for all acquisitions of a session except in multi-plane sessions
(not relevant here). Accepted-ROI counts are cross-checked against plots/cross_session/session_qc.csv, which is what the aggregated files / analyses use.
Run: miniforge3\\envs\\suite2p_venv\\python.exe fov_inventory.py
"""
import os, glob, csv, datetime, numpy as np, tifffile
from fov_common import *


def header(nm):
    f = sorted(glob.glob(os.path.join(DATA, nm, '*.tif')))
    with tifffile.TiffFile(f[0]) as t: fd = t.scanimage_metadata['FrameData']
    return fd, f


def main():
    qc = {}
    for r in csv.DictReader(open(QC_PATH)): qc[r['session']] = r
    rows = []
    for nm in SESSIONS:
        fd, tifs = header(nm)
        fov = np.array(fd['SI.hRoiManager.imagingFovUm'], float)
        fx, fy = 2 * abs(fov[0][0]), 2 * abs(fov[0][1])
        nx, ny = int(fd['SI.hRoiManager.pixelsPerLine']), int(fd['SI.hRoiManager.linesPerFrame'])
        ax = fd['SI.hMotors.axesPosition']; sp = fd['SI.hMotors.samplePosition']
        ic = np.load(os.path.join(s2p_dir(nm), 'iscell.npy')); ops = np.load(os.path.join(s2p_dir(nm), 'ops.npy'), allow_pickle=True).item()
        n_acc = int((ic[:, 0] > 0).sum()); q = qc.get(nm, {})
        mt = datetime.datetime.fromtimestamp(os.path.getmtime(tifs[0]))
        rows.append(dict(session=nm, short=short(nm), tif_first_mtime=mt.strftime('%Y-%m-%d %H:%M'), n_tifs=len(tifs), zoom=fd['SI.hRoiManager.scanZoomFactor'],
                         fov_x_um=round(fx, 1), fov_y_um=round(fy, 1), nx=nx, ny=ny, um_per_px_x=round(fx / nx, 4), um_per_px_y=round(fy / ny, 4),
                         axes_x=ax[0], axes_y=ax[1], axes_z=ax[2], sample_dx=round(sp[0] - ax[0], 1), sample_dy=round(sp[1] - ax[1], 1), fastz_pos=fd['SI.hFastZ.position'],
                         power=fd['SI.hBeams.powers'][0], n_slices=fd.get('SI.hStackManager.numSlices', 1), frame_rate=round(float(fd['SI.hRoiManager.scanFrameRate']), 3),
                         n_frames=int(ops.get('nframes', 0)), minutes_qc=q.get('minutes', ''), n_roi_total=len(ic), n_roi_accepted=n_acc, n_accepted_qc=q.get('rois_accepted', ''),
                         roi_count_matches_qc=('' if not q else str(int(q.get('rois_accepted', -1)) == n_acc)),
                         beh_face_sd_log=q.get('beh_face_sd_log', ''), pc1_signal=q.get('pc1_signal', ''), plane_sum=('yes' if int(fd.get('SI.hStackManager.numSlices', 1)) > 1 else 'no')))
        print(rows[-1]['short'], rows[-1]['tif_first_mtime'], 'zoom', rows[-1]['zoom'], '%.3fx%.3f um/px' % (fy / ny, fx / nx), 'ROIs', len(ic), 'acc', n_acc, 'qc', q.get('rois_accepted'))
    p = os.path.join(OUT, 'fov_inventory.csv')
    with open(p, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print('written', p)
    bad = [r['short'] for r in rows if r['roi_count_matches_qc'] == 'False']
    print('accepted-count mismatches vs session_qc.csv:', bad if bad else 'none')


if __name__ == '__main__':
    main()
