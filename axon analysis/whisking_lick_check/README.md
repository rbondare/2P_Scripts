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
- Note: `whisk_index` in the csv files that `lick_triggered.py`, `classes.py` and `lick_whisk.py` read is the NEW index since
  2026-11-10; their results (figures, tables) were made with the old one (`whisker_results_v1`). `me_fan_front`, `me_fan_down`,
  `me_pad`, `me_ctrl` are unchanged.
