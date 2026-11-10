# RB24: where are the axons in each recording, and what differs between days

Built 2026-11-09 by `axon analysis/fov_matching/` (Python, `miniforge3\envs\suite2p_venv\python.exe`; reads only `suite2p\plane0\{ops,stat,iscell}.npy` and the first TIFF header of each session on `Z:\Group Members\Rima\DATA_2P`; no data.bin, nothing modified, nothing committed).
20 RB24 sessions with suite2p output (261017_1620 = 3 planes 40 um apart and 261109_1506 = suite2p not finished are not included).

## What was done
1. `fov_inventory.py` -> `fov_inventory.csv`: zoom, FOV (um), um/px, stage (`axesPosition`), `hFastZ.position`, power, length, ROI counts (accepted counts equal `session_qc.csv` for every session that is in it).
2. `fov_register.py`, `fov_pose.py`, `fov_refine.py` -> `fov_pose_fine.csv`: mean images resampled to microns (pixel sizes 0.10-0.28 um/px, 512x512 or 512x1024 frames), masked normalised cross-correlation of axon-scale (DoG 0.8/6 um) and vessel-scale (3/15 um) features, 5 candidate peaks per pair, phase-randomised null, consensus placement (a session is only placed where candidates predicted by several already-placed sessions agree), sub-micron refinement.
3. `fov_roi_match.py` -> `roi_match_pairs.csv`: accepted ROIs rasterised on a 0.25 um grid (dilated 0.75 um), Hungarian matching at IoU >= 0.2, random-displacement control, same-day ceiling.
4. `fov_master.py` -> `master_axons.csv`, `master_presence.csv`: constrained clustering of matches into master axons (one ROI per session per axon); for every axon x session: accepted / rejected suite2p ROI present / structure but no ROI / nothing visible / out of FOV, using mean-image, Vcorr and max-projection contrast at the footprint versus 400 displaced background footprints of the same session.
5. `fov_analysis.py` -> `fov_session_summary.csv`, `fov_persistence.csv`, `fov_carryover.csv`; `fov_figures.py` -> F1-F6; `fov_selftest.py` (synthetic check).

## Registration: how far to trust it
- Self-test (real mean image re-sampled at 0.138 / 0.211 / 0.278 um/px with known offsets, noise 0.3-1.0 SD): recovered to <= 0.05 um (median 0.02).
- Local deformation between well-correlated pairs after the rigid shift: 0.31 um RMS, 0.56 um max (`fov_local_deformation.csv`) -> rigid translation is enough.
- Independent landmark (the bright blob with dark halo, ~(24, 97) um in the reference frame): found within 0.35 um of the prediction in 15 of the 19 sessions that contain it. Exceptions: 10-23 1538 (1.7 um), 10-30 1744 (0.8), 10-29 1752 (3.1; weakest placement, median NCC 0.38), 11-06 1604 (5.0; the blob has a dark centre in 1555/1604, so its maximum moved; 1604 vs 1555 are aligned to 0.3 um). 10-21 1826 does not contain the blob and has only 2 supporting edges (low confidence); 10-22 2026 also only 2 edges but its blob lands within 0.15 um.
- Reproduces the earlier ad-hoc 1423 -> 1520 shift: (0.5, 1.1) um vs (4, 10) px x 0.101 um.
- **Stage coordinates do not predict the image shift** (20-300 um differences, images line up within ~6 um on 11-05/06); `axesPosition` / `hFastZ` are descriptive columns only (units not verified).
- 11-04..11-06 sessions are one FOV (NCC 0.9-0.98); the 10-21..10-30 fields are larger (108-159 um) and contain the 104x52 um November window (F1, F2). The vessels look bright in 10-21/23/29/30 images and as dark shadows in 10-22/27 and in November (F2).

## Results
Sessions in time order; "confirmed" axon = accepted ROI in >= 2 OTHER sessions and fully inside this session's field (leave-one-out, so a session's own picks do not define the target set).

| | 10-21..10-30 (10 sessions, 2-53 accepted ROIs) | 11-04..11-06 (9 sessions, 63-160 accepted) |
|---|---|---|
| recall of confirmed axons (accepted again) | 0-0.18 | 0.30-0.58 |
| a suite2p ROI exists but was rejected | 0.03-0.22 | 0.17-0.45 (1604: no rejected class) |
| structure in the image, no ROI at all | 0.04-0.21 | 0.04-0.14 |
| nothing visible at that place | 0.56-0.86 | 0.11-0.29 (1604: 0.45) |

1. **Matching works and the same axons come back.** Matched fraction of the smaller evaluable ROI set (chance 0.07-0.15): same day / same stage position 0.71 (n 11 pairs, 0.50-0.82), 11-04..06 one day apart 0.64 (n 18), two days apart 0.56 (n 8), Oct vs Nov (>= 6 days) 0.15 (n 31, chance 0.12; only a few October ROIs lie inside the November window), different October days 0.28 (n 14, 0.00-0.78; 10-21 1630 / 10-23 1527 / 10-23 1538 share axons at 0.56-0.75) (F3). Carry-over P(accepted in t | accepted in s): 0.61 same day, 0.57 one day, 0.47 two days (0.66 / 0.59 without the accept-all 1604) (F5d).
2. **A large part of "different number of axons" is detection / curation, not visibility.** In November a session misses 42-70% of the axons the other sessions agree on, but 17-45% of them already have a suite2p ROI that was rejected, 4-14% show structure without an ROI, and only 11-29% show nothing at all (F5a, F4). The same-day ceiling is only ~0.7: even with the identical field, a third of the accepted axons are not re-accepted.
3. **Strongly axon-specific.** 194 axons are fully in view of all 9 November sessions and accepted at least once: 61 in a single session (random-subset model with the same per-session rates: 17 +- 8) and 21 in >= 7 sessions (model ~2); 68 of 465 master axons (212 confirmed) are accepted in >= 5 sessions (F5c, F4, F6). So there is a stable core and a sporadic set; the sporadic ones are the candidates for "not active at baseline" but also for noise / curation idiosyncrasy.
4. **Within 11-06 the field changes although stage and zoom are identical:** image similarity to 13:17 falls 1.0 -> 0.78 (15:30) -> 0.59 (15:56) -> 0.45 (16:04), recall of the confirmed axons 0.58 -> 0.53 -> 0.30 -> 0.41, "nothing visible" 0.14 -> 0.11 -> 0.29 -> 0.45 (F5f). Shifts are <= 3 um, so it is appearance (focus/z, brightness, the blob changing), not a lateral move. 1555 is only 4.2 min long (shorter recordings accept fewer axons: F5e), 1604 has no rejected class.
5. **October vs November:** the axons confirmed in November are mostly invisible at their location in the October images (contrast at background level, F5b, F6) and the October sessions accepted few ROIs. This could be real (expression/axon growth, different focal plane) or lower resolution (0.14-0.28 um/px) plus stricter curation; the data cannot separate these.

## What the data cannot say
- **Silent vs absent.** Only the green channel was saved (`channelSave=1`), so an axon that is dim at baseline and silent cannot be told from an axon that is absent or out of the focal plane. "Nothing visible" means exactly that, not "silent".
- z (focal plane) differs between days and is not corrected (`hFastZ.position` 30-370, motor z 4900-8500, units unverified; 5-plane tests earlier suggest ~+-6 um changes little).
- ROI sets come from different suite2p settings and curation passes (10-21..10-30 older/strict, 11-04+ merged / re-curated, 1604 accept-all of a curated pool). `iscell` is the user's curation, not a blind criterion.
- **261104_1640 is excluded from all ROI analyses**: its plane-sum ROIs are horizontal bands covering the field (accept-all of the default sourcery result); it is still used for registration. ROIs whose centre lies within 10 um of the blob landmark are dropped (it is not an axon).
- Sessions on one day are not independent; n is small; the within-day drift on 11-06 is described, not explained.

## Files
Figures: `F1_registration`, `F2_registered_gallery` (all 20 recordings, same window, red = accepted ROIs), `F3_roi_matching`, `F4_presence_matrix` (212 confirmed axons x 19 sessions), `F5_what_differs`, `F6_example_axons` (PNG + PDF); `Q1-Q3` are the first QC plots. Tables: `fov_inventory.csv`, `fov_pose_fine.csv`, `fov_landmark_check.csv`, `fov_local_deformation.csv`, `roi_match_pairs.csv`, `master_axons.csv` (positions in um in the 11-05 1520 frame), `master_presence.csv`, `fov_session_summary.csv`, `fov_persistence.csv`, `fov_carryover.csv`. `cache/` holds pickles (pairwise correlation candidates, poses, clusters).
Rebuild: `fov_inventory.py`, `fov_register.py` (~3 min), `fov_pose.py`, `fov_refine.py` (~3 min), `fov_roi_match.py` (~8 min), `fov_master.py`, `fov_analysis.py`, `fov_figures.py`.
