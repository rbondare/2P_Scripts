# RB22-24 recordings database

One row per recording, linked tables, plain CSV (opens in Excel, `readtable` in MATLAB, `csv` in Python, diffable in git).
Master copy: `Z:\Group Members\Rima\Recordings_DB\`   Code: `2P_Scripts\axon analysis\recordings_db\`   Browse: open `index.html` (filter box, sortable, FOV pictures).

## Tables
| File | One row per | Key columns |
|---|---|---|
| `sessions.csv` | recording (every `DATA_2P` folder of RB22-24) | `session` (= `AnimalRB24_261105_1423`), `status`, `task_type`, `stimuli`, `n_roi_accepted`, `pool`, `flag_reason`, `all_tags`, `stale`, acquisition + FOV geometry, behaviour / timing QC |
| `blocks.csv` | stimulus block (`stim*.mat`) | `session`, `block`, `key` (flash / fff / bar / loom_expand / loom_recede / square / spont / alleyCL), `n_trials`, `interrupted`, `t_start_2p`, `red_sync_ok`, `params` |
| `rois.csv` | suite2p ROI, accepted AND rejected | `session`, `roi` (suite2p index), `accepted`, `agg_row` (row in the aggregated `CaData`), `class` (signal / noise / neuropil), `y_um`, `x_um`, `master_axon` |
| `axons.csv` | master axon matched across RB24 sessions | `master`, position in the common um frame, `sessions` (tentative, see caveats) |
| `animals.csv` | animal | genotype, sex, dates, indicator, virus: **hand-filled** |
| `notes.csv` | session | **hand-edited**: `pool`, `flag_reason`, `curation`, `tags`, `quality`, `visible_activity`, `comments` + `notion_*` columns |
| `fov\<session>.png` | session | mean image with accepted ROI outlines (cyan = signal-like, yellow = noise-like, red = neuropil-limited, green = not classified) |

## Select recordings
```matlab
names = recdb_select('animal', 'RB24', 'has_block', 'loom_expand', 'min_rois', 10, 'pool', true);   % cellstr of session names
[names, T] = recdb_select('tags', {'task'}, 'exclude_tags', {'stale'}, 'date_from', '2026-10-27');   % T = matching sessions.csv rows
R = recdb_load('rois');   B = recdb_load('blocks');
```
```python
import sys; sys.path.insert(0, r'C:\Users\rbondarenko\projects\2P_Scripts\axon analysis\recordings_db'); import recdb
sel = recdb.select(animal='RB24', has_block='loom_expand', min_rois=10, pool=True); recdb.show(sel)
```
Filters: animal, status, task_type, pool, tags / any_tags / exclude_tags, has_block (+ min_trials), min_rois, min_minutes, date_from / date_to, aggregated (Python also `where=lambda row: ...`).

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
* Blank = unknown, never a guess. `class` (signal / noise / neuropil) exists only for the sessions of the cross-session export (`2P_cache\xs`); `master_axon` and `pose_*` only for RB24 (FOV matching 2026-11-09).
* `stale = yes`: the aggregated file was made from a different ROI selection than the current `iscell.npy` (different ROI count, or `iscell.npy` newer). `agg_row` and `class` are left blank then, re-aggregate first.
* `master_axon` ids come from `fov_matching`: no reference channel, so "silent" and "absent" axons cannot be told apart, and 1640's ROIs are junk bands. Treat as tentative.
* `AnimalInfo` inside the aggregated files holds the dialog defaults (Ai148 / Rorb, GCaMP6f, placeholder dates), so animal facts come only from `animals.csv`.
* `date_actual` (first TIFF time stamp) differs from the folder date for `AnimalRB24_261021_1826` (recorded 2026-10-28; its stimulus folder is `..._261028_1826`).
* `n_trials` = completed trials found in the data (same onset logic as `cross_session/export_xs.m`); raw-only sessions have it for loom / bar / fff blocks only. `trials_requested` is what the stimulus script was asked for.
* `ball_*` columns are what the aggregated file contains; whether the ball log is usable is in `ball_note` (from `session_qc.csv`, 20 sessions) and `ball_sync_ok` (training sessions only).
