# RB22-24 recordings database

One row per recording, linked tables, plain CSV (opens in Excel, `readtable` in MATLAB, `csv` in Python, diffable in git).
Master copy: `Z:\Group Members\Rima\Recordings_DB\`   Code: `2P_Scripts\axon analysis\recordings_db\`   Browse: open `index.html` (filter box, sortable, FOV pictures).

## Tables
| File | One row per | Key columns |
|---|---|---|
| `sessions.csv` | recording (every `DATA_2P` folder of RB22-24) | `session` (= `AnimalRB24_261105_1423`), `status`, `task_type`, `stimuli` (readable summary, see below), `spout`, `task_level`, `n_roi_accepted`, `pool`, `flag_reason`, `all_tags`, `stale`, acquisition + FOV geometry, behaviour / timing QC |
| `blocks.csv` | stimulus block (`stim*.mat`) | `session`, `block`, `key` (flash / fff / bar / loom_expand / loom_recede / square / spont / alleyCL), `stim_desc` (one readable line), `n_trials`, `interrupted`, `t_start_2p`, `red_sync_ok`, `params`; square blocks also the task columns below |
| `rois.csv` | suite2p ROI, accepted AND rejected | `session`, `roi` (suite2p index), `accepted`, `agg_row` (row in the aggregated `CaData`), `class` (signal / noise / neuropil), `y_um`, `x_um`, `master_axon` |
| `axons.csv` | master axon matched across RB24 sessions | `master`, position in the common um frame, `sessions` (tentative, see caveats) |
| `animals.csv` | animal | genotype, sex, dates, indicator, virus: **hand-filled** |
| `notes.csv` | session | **hand-edited**: `pool`, `flag_reason`, `curation`, `tags`, `quality`, `visible_activity`, `comments` + `notion_*` columns |
| `fov\<session>.png` | session | mean image with accepted ROI outlines (cyan = signal-like, yellow = noise-like, red = neuropil-limited, green = not classified) |

## Stimuli: what is saved per block
`blocks.csv` has `stim_desc` (full) and `stim_short` (what `sessions.stimuli` is built from) for every block, e.g.
`square | 100% contrast | centre | size 1 | level 2.5 adaptive (L2 59 / L3 61) | reward mode 1 (reward delay 1.5 s) | lick window 0.1-3 s | wait 8-13 s before the square | spout (4734 lick entries in 114 trials)`.
`sessions.stimuli` then reads e.g. `square 100% | centre | level 2.5 adaptive (L2 59 / L3 61) | spout x120` (x = completed trials) or `spont; flash x32; loom_expand (+grey) x10`. `(+grey)` = the loom block asks for 120 s of grey before / after (recorded outside the stimulus frames).
Square (training) blocks, read from the raw `stim*.mat` (so also for sessions that are not aggregated), have their own columns; blank for all other stimuli:
| Column | Meaning |
|---|---|
| `stim_shape`, `stim_contrast_pct`, `stim_size`, `stim_location` | stimulus type: `stimulusType` (square), `StimulusContrast` (30 / 100 %), `StimSize`, `locationMode` (centre) |
| `gui_level` / `level_mode` | the level set in the GUI: 1, 2, 3 or **2.5 = adaptive** (`level_mode = adaptive`: every trial runs level 2 or 3, drawn from the agent state). `specParams.level` in the file is only the level of the LAST trial, so it is not used |
| `n_L1`, `n_L2`, `n_L3` | trials per level: from `levelHistory` in adaptive blocks, all `trials_logged` at the GUI level in fixed blocks |
| `reward_mode`, `reward_delay_s`, `lick_window_s`, `wait_range_s` | `RewardMode` (0 = square shown without reward, 1 = spout + reward), `RewardDelay`, lick window min-max, wait before the square min-max |
| `spout` | **yes = at least one lick is in `lickData` of the block, no = none at all** (rig convention: no spout / reward mode 0 leaves `lickData` empty). Decided by the licks only, not by `reward_mode`; a disagreement is written to `spout_note` |
| `n_lick_entries`, `n_trials_with_licks`, `lick_log_trials` | raw lick-log entries (duplicates not removed, so slightly more than the lick count of the task analysis), trials with at least one, and trial columns in the log (can be one short: 1506) |
`sessions.csv` carries the summary (`spout` yes / no / mixed, `task_level`, `level_mode`, `level_counts`, `n_lick_entries`, `spout_note`); `spout` and `no_spout`, `adaptive_level` are also auto tags. Sessions without a square block have `spout` blank (unknown, not "no").
Before 2026-11-11 `sessions.reward_level` held `specParams.level`, which is wrong for adaptive sessions (3 or 2 instead of 2.5); it is now `task_level`.

## Select recordings
```matlab
names = recdb_select('animal', 'RB24', 'has_block', 'loom_expand', 'min_rois', 10, 'pool', true);   % cellstr of session names
names = recdb_select('task_type', 'lick task', 'spout', true, 'level', 2.5);                       % adaptive-level sessions with a spout
[names, T] = recdb_select('tags', {'task'}, 'exclude_tags', {'stale'}, 'date_from', '2026-10-27');   % T = matching sessions.csv rows
R = recdb_load('rois');   B = recdb_load('blocks');
```
```python
import sys; sys.path.insert(0, r'C:\Users\rbondarenko\projects\2P_Scripts\axon analysis\recordings_db'); import recdb
sel = recdb.select(animal='RB24', has_block='loom_expand', min_rois=10, pool=True); recdb.show(sel)
```
Filters: animal, status, task_type, pool, spout, level, tags / any_tags / exclude_tags, has_block (+ min_trials), min_rois, min_minutes, date_from / date_to, aggregated (Python also `where=lambda row: ...`).

## Rebuild
Run `build_recdb` (MATLAB; `.mlx` or `.m`). It runs `export_session_meta.m` (stimulus blocks and the QC that lives in the aggregated files; only sessions whose files changed) and then `build_recdb.py` (suite2p / TIFF header / FOV pictures / joins). Takes ~1 min when nothing changed.
Re-run after: aggregating a session, curating `iscell.npy`, a new recording, or editing `notes.csv` (so that `sessions.csv` / `index.html` pick it up).

## What is automatic and what is yours
* **Automatic (overwritten every build):** everything in `sessions.csv` / `blocks.csv` / `rois.csv` / `axons.csv` except the columns copied from `notes.csv`.
* **Yours (never overwritten):** `notes.csv`, `animals.csv`. New sessions get a seeded row (values taken once from the lists that used to be hard-coded in `paper_figures/fig_common.py`: `FLAG` -> `pool=no` + `flag_reason`, `TASK_SESSIONS` / `GREY_SESSIONS` -> `tags`). Edit freely; existing cells are never touched.
* `pool`: `yes` = goes into pooled statistics, `no` = drawn but not pooled, blank = not decided / not analysed.
* `tags` (yours) are merged with `auto_tags` (derived from the data) into `all_tags`, which is what the filters use.
* The Notion table: export it to CSV and run `import_notion.py <csv>`; its columns land in the `notion_*` columns of `notes.csv`.

## Caveats (read before trusting a column)
* Blank = unknown, never a guess. `class` (signal / noise / neuropil, drawn cyan / yellow / red in the FOV pictures) is computed for every aggregated, non-stale session from its cross-session export (`2P_cache\xs\<session>_xs.mat`, same criterion as `cross_session/xs_common.py`). The build makes a missing xs file itself (`export_session_meta.m`, option `xs`); an existing xs file is never re-made, so after re-aggregating a session with the SAME ROI count delete its xs file. Green outlines = no class available (not aggregated / stale / xs export failed). `master_axon` and `pose_*` only for RB24 (FOV matching 2026-11-09).
* `stale = yes`: the aggregated file was made from a different ROI selection than the current `iscell.npy` (different ROI count, or `iscell.npy` newer). `agg_row` and `class` are left blank then, re-aggregate first.
* `master_axon` ids come from `fov_matching`: no reference channel, so "silent" and "absent" axons cannot be told apart, and 1640's ROIs are junk bands. Treat as tentative.
* `AnimalInfo` inside the aggregated files holds the dialog defaults (Ai148 / Rorb, GCaMP6f, placeholder dates), so animal facts come only from `animals.csv`.
* `date_actual` (first TIFF time stamp) differs from the folder date for `AnimalRB24_261021_1826` (recorded 2026-10-28; its stimulus folder is `..._261028_1826`).
* `n_trials` = completed trials found in the data (same onset logic as `cross_session/export_xs.m`); raw-only sessions have it for loom / bar / fff blocks only. `trials_requested` is what the stimulus script was asked for.
* `ball_*` columns are what the aggregated file contains; whether the ball log is usable is in `ball_note` (from `session_qc.csv`, 20 sessions) and `ball_sync_ok` (training sessions only).
