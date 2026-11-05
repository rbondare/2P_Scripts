# Motion vs signal for 5-plane (+ flyback) axon recordings

Registration-free checks of "is this signal calcium or motion?" for ScanImage fastZ recordings (5 planes + 1 flyback frame
per volume), written for `RB24_261104_1651`. They run in the `suite3d` conda env (numpy / scipy / tifffile / matplotlib).
Everything works on a local copy of the tifs; nothing writes to Z:.

| step | script | what it does | output |
| --- | --- | --- | --- |
| 1 | `prep_binned_movie.py <tif_folder> <out>` | one pass over the tifs, drops the flyback frame and rig stubs, 2x2-bins to an int16 movie (unregistered; also the raw mean image) | `movie_bin2.npy`, `mean_img.npy` |
| 2 | Suite3D job (`../run_suite3d.bat`) | rigid x/y shift per volume (the z shift is measured, not applied), 3D ROIs | `~\Suite3D\s3d-<job>` |
| 3 | `prep_rigid_registered.py <tifs> <job> <out> --zsum-dir <dir>` | applies ONLY Suite3D's rigid y/x shifts to the raw planes: binned movie for the checks below + the plane-sum movie for suite2p | `regraw_bin2.npy`, `zsum_*.tif` |
| 4 | `axial_profile_check.py <out>` | brightness of bright objects across the 5 planes (raw mean image): are the planes different depths at all? | `axial_profile.png` |
| 5 | `z_motion_estimate.py <out> <job> --movie regraw_bin2.npy` | per-volume axial displacement from the plane profile (gain + shift model), split-half reliability, comparison with Suite3D's z | `z_motion.npz/.png` |
| 6 | `zreg_match_diagnostic.py <out> --movie regraw_bin2.npy` | steps 1-3 of suite2p-in-depth's volume z-registration (Baruchin et al. 2026) re-implemented: frame x reference-plane correlations, best-match plane, bounds | `zreg_match.png` |
| 7 | `roi_motion_audit.py <out> <job> [--s2p <suite2p root>] --out <dir>` | per ROI set: autocorrelation, distance-dependent correlation, fast-motion R2 vs circular-shift null, events in moving volumes, size of the artefact a one-plane stack shift would cause | `roi_audit/` |
| 8 | `event_triggered.py <out> <suite2p plane dir>` | triggered averages of population dF/F, lateral offset, contrast gain, z estimates at events / bout onsets (needs step 7 run with `--out <out>\audit_zsum`) | `event_triggered.png` |
| 9 | `cross_plane_consistency.py <out> <suite2p plane dir>` | fast-band correlation of the same ROI footprint between planes vs plane separation | stdout |
| 10 | `z_injection_test.py <out>` | injects known axial shifts / gain changes into the real movie and measures what the plane-profile estimator returns: linearity, leak of calcium-like gain into z, minimum detectable shift vs duration | stdout |

For aggregation, build the plane-sum with `prep_rigid_registered.py ... --keep-stubs` (stub tifs included with zero shift) so that
suite2p has one frame per volume of EVERY tif: the aggregator indexes the traces with the tif-header timestamps and fails or shifts
timing otherwise. The aggregation timing checks are `newAggregation/verify_aggregation_timing.m` and
`newAggregation/video_ball_lag_check.py`.

Do NOT build these from Suite3D's registered movie (`bin_registered.py`, `make_zprojection_tif.py` are kept only to document
the problem): on 5 planes @ 3 um Suite3D applies a non-rigid z shift to every block (mean -1.3 planes with `max_shift_nr`
z = 1, exactly -1.0 with z = 0) with edge replication, so registered plane 0 == plane 1 and the other planes are mixed or
shifted. Single-pixel fast-band correlation between planes: raw 0.02-0.04, Suite3D registered 0.9-1.0.

Plain suite2p on a 5-plane recording: `../../suite2p_analysis/run_suite2p_RB_axons_5plane_V1.py <tifs> --save-path <local>`
(V2 parameters, `nplanes` = 6 with the flyback frame ignored, per-plane fs = volume rate). For the plane-sum movie:
`--n-frames-per-volume 1 --flyback`.

Caveats
* Step 6 is NOT the suite2p-in-depth package, only its matching step on an already x/y-registered, binned movie.
* Lateral-motion / contrast-dip "moving" masks (step 7) are defined from Suite3D's shifts; if those fail the masks do.
* Axial displacement is expressed in planes; the planes' physical spacing was not verified (header: 3 um).
