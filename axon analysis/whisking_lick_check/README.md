# Whisking vs licking check (RB24_261109_1506) and the new WhiskIndex (2026-11-10)

Why: the whisking index (`DLC/face_model/face_motion.py`, `Behav.WhiskIndex`) was the mean of the motion energy of two
boxes, one of which (fan_down) holds chin, jaw and tongue. With the independent lick line of RB24_261109_1506 (4,734 licks)
licking explained 45% of its variance and 93% of its running-vs-still contrast. The new index uses only the upper whisker
region (`front_top`, above the lick spout; lick effect 4% of the running contrast, CI 0-9%, 0.2% lick-explained variance,
Cohen d running vs still 1.15 vs 1.02 for the old index). `Behav.JawME` carries the old fan_down motion energy.

- `*.py` in this folder: the lick check itself (`read_triggers.py` reads the lick line from the TIFF headers, `cell_me.py`
  makes a 20 px cell grid of motion energy, `lick_whisk*.py`, `lick_triggered.py`, `classes.py`, `maps.py`, `subbox.py`,
  `frames_spout.py`); figures `LW1_*`, `LW2_*`, table `lw_table.csv`. The big caches (cell grid, ball, lick times) are not in
  the repo; the scripts expect them next to the script (they were made in a scratch folder).
- `fix2_rollout/`: how the new index was rolled out to the 43 RB22-24 sessions
  1. `pass1_all.py`  face_motion.py on every raw video into a scratch folder (motion energy of all boxes, ~1 h, 3 ffmpeg readers)
  2. `DLC/face_model/calibrate_rest.py` -> REST per codec (mpeg4 0.473, h264 0.952), set in `face_motion.REST`
  3. `validate_clips.py` scores the new index against the 24 blind clips (`DLC/whisking_check`); `dry_summary.py` classifies all
     sessions without writing; `geometry_mosaic.py` checks the spout geometry in every session (spout below `front_top`
     in all sessions where it is visible, RB22 / RB23 / RB24)
  4. `pass2_apply.py` `face_motion.py --reclassify` (no decoding) writes the real `whisker_results` (old version kept in
     `whisker_results_v1`; `--archive` does the same in `Z:\...\preCamera`)
  5. `run_refresh.m` + `newAggregation/refresh_whisker_behavior.m` replace the whisking fields in the 23 aggregated files
     (old fields backed up to `C:\Users\rbondarenko\2P_cache\backup_whisk_v1`), `verify_refresh.m` checks them
  6. (2026-11-10, second step) per-codec on/off levels: `face_motion.ON_OFF` = MPEG-4 on above 1.76 / off below 1.44 (unchanged),
     H.264 on above 1.45 / off below 1.40, set from 51 blind human clip labels (`h264_blind_check/`: the old rule missed 6
     whisking clips, the new one 3 wrong in 51) and checked against the LC-axon dF/F (coupling of the binary state, 6 sessions,
     `h264_blind_check/coupling_check.py`). Rolled out with `h264_threshold_rollout.py`: `stage` (all 43 sessions re-classified
     from the saved ME traces, regression test: the 29 MPEG-4 sessions identical in every column, the 14 H.264 sessions differ in
     `whisk_on` only, `h264_threshold_stage_comparison.csv`), `apply` (work, quality_reports, preCamera / Camera copies; previous
     version kept next to each as `whisker_results_v2`), `sync-stale` (18 Camera-folder copies of not-yet-aggregated MPEG-4 sessions
     still held the ORIGINAL two-box version and would have been imported by the next aggregation; backup `whisker_results_v1`),
     then MATLAB `run_refresh_h264.m` (12 aggregated H.264 files + DLC_filtered, backups `2P_cache\backup_whisk_v2\<s>_whisk_v2.mat`).
     Verification: `verify_h264_threshold.py` (215 copied files byte-identical to work, independent re-implementation of the
     rule equal in every frame of all 43 sessions, json levels, backups, the 51 clips), `verify_refresh_h264.m` (aggregated file
     == csv, only WhiskOn changed, new on-frames a superset of the old, DLC_filtered == csv, other Behav fields == DLC_filtered, size),
     `verify_h264_e2e_coupling.py` (axon coupling recomputed from the aggregated WhiskOn = the numbers of coupling_check.py).
     `h264_blind_check/strat_*.py` and `coupling_check.py` call `face_motion.classify(me, rest, codec)` (new argument) and override
     `fm.ON_OFF['h264']` for their candidate rules.
- Note: `whisk_index` in the csv files that `lick_triggered.py`, `classes.py` and `lick_whisk.py` read is the NEW index since
  2026-11-10; their results (figures, tables) were made with the old one (`whisker_results_v1`). `me_fan_front`, `me_fan_down`,
  `me_pad`, `me_ctrl` are unchanged.
