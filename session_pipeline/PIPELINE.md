# Session Pipeline — clean, per-session stimulus organization + QC

## Context / why this exists

The repo currently has three overlapping generations of analysis code:

1. **`newAggregation/`** — raw suite2p → `*_preprocessed.mat` (the *upstream* aggregation; not touched here).
2. **The 3-script pipeline** (`Preprocess_Stimulus_Data.m` → `analysis_functions/` → `Plot_Stimulus_Results.m`) — correct chunking/timing logic, but **baseline-vs-drug-pair-centric**: every extractor is called inside a pairwise `analyze_*_stimulus` function, so there is no reusable, cached "clean analysis of one session" the way `RB_DataGathering.m` intended.
3. **`Master_Stimulus_Analysis.m` + many one-offs** — the original monolith and poster/exploration scripts.

`RB_DataGathering.m` is the organizational pattern the user likes (`stim_responses.<stim_type>.<recording>` with parallel `param_values` / `responses{}` / `n_presentations_per_param`), **but it is a stub**: its per-presentation loop copies the whole block into every trial slot (`% In real implementation, would slice by presentation`) and never averages. It also has no red-frame/timing handling.

**This pipeline replaces `RB_DataGathering.m` done right**: each aggregated file is independently chunked, aligned, and averaged per stimulus into **one canonical struct**, saved to disk (`*_organized.mat`), with **hard-assert QC + an auto-generated QC verification figure**. Baseline-vs-drug (or any) comparison becomes a *thin* downstream layer that just loads two canonical structs. The heavy timing/RNG-reconstruction work is **reused verbatim** from the already-validated `analysis_functions/` extractors — nothing about timing is re-derived here.

### Design decisions (confirmed with user)
- **Unit of analysis:** per aggregated file → one canonical struct; comparison is a separate thin layer.
- **Intermediate:** saved to disk as a canonical `*_organized.mat` cache (with provenance metadata); plotting/exploration read the cache, never re-derive timing.
- **QC:** hard asserts that halt on data-integrity violations, plus recorded pass/fail scientific-validity flags and a per-session QC figure that must be reviewed before a session is trusted.

---

## Architecture

```
                 raw suite2p                (upstream, unchanged)
                     │  newAggregation/
                     ▼
        Z:\...\Aggregated\<session>_preprocessed.mat
                     │
   ┌─────────────────┴─────────────────┐
   │  SESSION PIPELINE (this folder)    │
   │                                    │
   │  run_organize.m   (config + batch entry point)
   │        │  for each session file:
   │        ▼
   │  organize_session.m  ── reuses ──►  analysis_functions/
   │        │   (extract_looming_trials, extract_direction_trials,
   │        │    get_flash_onsets + align_trials_to_onsets,
   │        │    extract_full_stimulus_responses, chunk_by_direction,
   │        │    reconstruct_direction_order, reorder_subtrials_by_direction)
   │        ▼
   │     S  (canonical per-session struct)
   │        │
   │        ├─►  qc_session.m         (asserts + scientific-validity flags)
   │        ├─►  qc_report_figure.m   (timing-alignment verification PNG)
   │        ▼
   │  organized_data\<session>_organized.mat   +   QC_<session>.png
   └────────────────────────────────────────────┘
                     │
        ┌────────────┴─────────────┐
        ▼                          ▼
  explore_session.m         compare_sessions.m
  (single-session plots)    (thin pairwise COMPUTE: matched-cell
                             raw/z, MI, signrank p  ->  R struct)
                                   │
                                   ▼
                            compare_plots.m
                            (violins / modindex / histograms /
                             scatter / summary -- PURE viz of R)
                                   │
                                   ▼
                            compare_notebook.m
                            (celled %% orchestration, run cell-by-cell)
                                   │
                                   ▼
                            compare_group.m   ← group_config.m (the cohort)
                            (N pairs -> per-pair R + ROI-match verification,
                             pooled AND per-recording stats, per condition)
                                   │
                         ┌─────────┴──────────┐
                         ▼                    ▼
                  group_plots.m       group_travelingwave.m
                  (traces/superplot/  (baseline|drug per direction,
                   mi/forest/summary)  cells pooled per condition)
                         │
                         ▼
                  group_notebook.m
```

**Separation of concerns** (the whole point):
- *Per-session* work (chunking, red-frame alignment, RNG direction reconstruction, trial averaging) happens **once per file** and is cached. This is inherently session-local — it does not need a drug pair or an ROI match.
- *Pairwise* work (matched-cell z-score, modulation index, baseline-vs-drug stats) is a thin layer over two cached structs + an ROI-match file. It re-derives nothing about timing.

---

## Canonical per-session struct `S`

`S = organize_session(preprocessed_file, cfg)`

### `S.meta`
| field | meaning |
|---|---|
| `session_id` | e.g. `AnimalRB19_260320_1444` (filename minus `_preprocessed`) |
| `source_file` | full path to the `*_preprocessed.mat` |
| `created` | timestamp string |
| `git_commit` | short hash of `2P_Scripts` at organize time (provenance) |
| `ca_type` | 1=`Ca_dFF`, 2=`Ca_deconvolved`, 3=`Ca_F` |
| `ca_field` | the resolved field name actually used |
| `selected_plane_idx`, `plane_z` | which imaging plane (1-based) and its centroid-Z value |
| `n_rois` | ROIs in the selected plane |
| `n_frames` | calcium frames |
| `fs_hz` | `1 / median(diff(TimeCa(plane,:)))` |
| `stim_types` | cellstr of stimulus types found |
| `qc_passed` | overall logical (all scientific checks passed) |
| `qc` | struct returned by `qc_session` (per-check detail) |

### `S.stimuli.<type>` — one field per stimulus type present
`<type>` is `matlab.lang.makeValidName(raw_type)`. Common fields on every entry:

| field | meaning |
|---|---|
| `type` | raw type string |
| `paradigm` | `'trial_onset'` \| `'direction'` \| `'block'` |
| `stim_index` | index into `B.Stimuli(...)` |
| `RedFrameSynchronized` | logical (from the stim block) |
| `n_rois` | == `S.meta.n_rois` |
| `qc` | per-stimulus QC sub-struct |

**`paradigm = 'trial_onset'`** (looming, sparse_local_global_flashes):
| field | shape | meaning |
|---|---|---|
| `t_com` | `1 × n_com` | trial-relative time grid (s), 0 = onset |
| `onsets` | `n_trials × 1` | onset times, camera clock |
| `trials` | `n_rois × n_com × n_trials` | per-trial aligned dF/F (raw, not averaged) |
| `avg` | `n_rois × n_com` | **median** across trials (per established convention) |
| `sem` | `n_rois × n_com` | s.e.m. across trials |
| `peak_window` | `1 × 2` | window used for amplitude metrics (s) |
| `n_trials` | scalar | |

**`paradigm = 'direction'`** (grating, moving_bar):
| field | shape | meaning |
|---|---|---|
| `t_com` | `1 × n_com` | subtrial-relative time grid (s) |
| `subtrials` | `n_rois × n_com × n_subtrials` | every subtrial aligned (raw) |
| `direction_of` | `1 × n_subtrials` | realized direction per subtrial (RNG-reconstructed) |
| `onsets` | `n_subtrials × 1` | subtrial onset times |
| `dir_values` | `1 × n_dir` | unique directions |
| `by_dir` | `n_rois × n_com × n_dir` | **median** of repeats per direction |
| `amp_by_dir` | `n_rois × n_dir` | peak-window mean per direction (tuning) |
| `pref_dir_idx` | `n_rois × 1` | index into `dir_values` of each cell's preferred direction |
| `pref_resp` | `n_rois × n_com` | each cell's own preferred-direction trace |
| `param_name` / `param_field` | e.g. `'orientation'` / `'bar_orientations'` | |
| `peak_window` | `1 × 2` | |
| `verify` | struct | `r_same`, `r_shifted`, `pass` from `verify_direction_reconstruction` |

**`paradigm = 'block'`** (spontaneous, checkers2, and generic fallback):
| field | shape | meaning |
|---|---|---|
| `block` | `n_rois × n_frames_block` | whole continuous block, native frame grid |
| `block_t` | `1 × n_frames_block` | real seconds within the block |
| `duration_sec` | scalar | |

---

## QC specification

`qc = qc_session(S, B, cfg)` — two tiers.

### Tier 1 — FATAL data-integrity (hard `assert`, halts this session)
- `length(TimeCa(plane,:)) == size(dff,2)` (time vector matches calcium frames).
- `fs_hz` within a plausible range (default `[2, 30]` Hz) — catches wrong-plane / corrupt-timebase.
- selected plane exists in centroid-Z.
- no `NaN`/`Inf` in the raw selected-plane dF/F beyond a tolerated fraction.
- for every trial/direction/block stimulus: `TimeStimulusFrame` monotonic non-decreasing, and its span lies within `[TimeCa(1), TimeCa(end)]`.

### Tier 2 — SCIENTIFIC-VALIDITY (recorded flags; set `qc_passed=false`, do **not** silently pass)
- `RedFrameSynchronized == true` for every stimulus block (the upstream fix must have applied).
- **onset count matches expectation**: looming → `Stimuli.trials`; direction → `n_dir × repeats`. Flashes have **no fixed count** (RNG-determined by cumulative random waits), so instead a **flash-timing** check asserts every inter-flash gap ∈ `[interval]` and the last onset ≤ `stimulus_trial_t` (± one display frame) — the falsifiable structural check for the reconstruction.
- all onsets fall inside `[TimeCa(1), TimeCa(end)]` with ≥ `pre+post` seconds of coverage.
- no all-NaN trial slots (every trial slice got ≥ 2 frames to interpolate from).
- direction reconstruction: `verify.r_same > verify.r_shifted + 0.05` (the existing PASS criterion).

Tier-2 failures don't halt (so a batch completes), but `S.meta.qc_passed=false`, the console prints a loud banner, and the QC figure title is stamped **QC FAILED** so it can't be missed.

### QC report figure (`QC_<session>.png`)
One row per stimulus. For each:
- **trial_onset / direction**: population-mean triggered average (mean over ROIs of `avg`/`pref_resp`) with `xline(0)` at onset + a raw heatmap (cells × time) with onset marked — a *correct* alignment shows the response band at the **same** relative latency for the population; a timing bug smears it.
- **block**: population-mean trace over the block + a raw heatmap.
- Title stamps `n_trials`/`n_dir`, `RedFrameSynchronized`, and PASS/FAIL per check.

This figure is the human gate: it makes a timing/alignment mistake *visible*, which is exactly what the "rigid, verified" requirement needs beyond numeric asserts.

---

## Files in this folder

| file | role |
|---|---|
| `run_organize.m` | **entry point.** Config block (list of session files, plane, ca_type, strict flag, output dir) → loops, organizes, QCs, saves `*_organized.mat` + `QC_*.png`. Wrapped per-file in try/catch so one bad file logs an ERROR and the batch continues. |
| `organize_session.m` | pure builder of `S` for one file; dispatches each stimulus type to the right reused extractor. |
| `qc_session.m` | Tier-1 asserts + Tier-2 flags → `qc` struct. |
| `qc_report_figure.m` | renders `QC_<session>.png`. |
| `explore_session.m` | load one `*_organized.mat`, plot per-stimulus averages / tuning / heatmaps for **single-session** exploration (no pair, no ROI match needed). |
| `compare_sessions.m` | thin pairwise **compute** layer: load two organized structs + ROI-match, emit matched-cell raw/z metrics + modulation index + paired `signrank` p per common stimulus as its own `R` struct. Verified (Step-7 walkthrough): matched pairs are the same neurons (centroid corr ≈ 1.0, dist ≪ shuffled), z-scores bit-exact vs by-hand recompute. |
| `compare_plots.m` | thin pairwise **plotting** library over `R`. `compare_plots(kind, R, 'metric','z'\|'raw', 'save',path)`; kinds: `violins`, `modindex`, `histograms`, `scatter`, `summary`. Pure visualization of values `compare_sessions` already computed (a plot can't disagree with the stats). Toolbox-free violin (manual Gaussian KDE) + robust [1,99]-pct axis clipping. |
| `compare_notebook.m` | celled (`%%`) orchestration — the exploratory comparison "notebook": CONFIG → compute `R` → per-session recap (`explore_session`) → violins / modindex / histograms / scatter / summary, each its own cell (run only what you want). |
| `group_config.m` | **the cohort**: one entry per baseline/drug pair (`baseline_file`, `drug_file`, `roi_match_file`, `condition`, optional `plane_base`/`plane_drug`). Single source of truth — the notebook and any batch driver both call it, so a figure and a table can't come from different cohorts. |
| `compare_group.m` | **group layer.** Runs `compare_sessions` on every pair, verifies each ROI match, then aggregates within condition two ways: **pooled** (all matched cells from all recordings, one signrank) and **per-recording** (each recording → one median, then summarised; `across_pair.underpowered` flags n<6). Also a between-condition `ranksum` on pooled MI. Re-derives nothing. |
| `group_plots.m` | figures over `G`; rows = condition, cols = stimulus. `traces` (grand mean across recordings ± s.e.m., thin line per recording), `superplot` (pooled violin + each recording's median), `mi`, `forest` (per-recording MI vs pooled median), `summary` (table). |
| `group_travelingwave.m` | paired traveling wave: per direction, the same matched cells drawn **baseline \| drug** in one cell order and one colour scale, cells pooled across the condition's recordings. Selection + latency ordering come from baseline repeat 1 only and the baseline panel shows the **remaining** repeats, so the diagonal is not an artefact of sorting the data being displayed. |
| `group_notebook.m` | celled orchestration for the group layer (compute → traces → superplot → mi → forest → summary → traveling wave → per-pair drill-down). |
| `organized_path.m` | the `<session>_plane<K-1>_organized.mat` naming convention, in one place. |
| `violin_kde.m`, `robust_lim.m` | shared plotting helpers (one copy, used by `compare_plots` and `group_plots`). |

## Reused from `analysis_functions/` (not reimplemented)
`extract_looming_trials`, `extract_direction_trials`, `get_flash_onsets`, `align_trials_to_onsets`, `extract_full_stimulus_responses`, `chunk_by_direction`, `reconstruct_direction_order`, `reorder_subtrials_by_direction`, `verify_direction_reconstruction`, `compute_baseline_zscore[_from_resp|_trace]`, `compute_whole_block_zscore`, `whole_block_zscore_resampled`, `compute_modulation_index`. Also `reconstructFlashes.m` (repo root).

## Keep / archive recommendations (not done automatically)
- **Keep & build on:** `newAggregation/`, `analysis_functions/`, this folder.
- **Deduplicate:** `RB_DataGathering.m` vs `organize_stimulus_by_parameters.m` are near-duplicate stubs — superseded by `organize_session.m`; archive both.
- **Archive (move to `legacy/`), don't delete:** `Master_Stimulus_Analysis.m`, `AggregationCode/`, `preprocess_and_organize_stimulus.m`, `stimulus_response_comparison_pipeline.m`, `Poster_*`, `Temp_*`, `Looming_Individual_Cells_0312.m`, `*.asv`.
- **Diagnostics:** fold `Diag_*` / `Test_*` intent into the QC figure where it belongs; keep the scripts as reference.

---

## How to run / verify (end-to-end)

```matlab
cd 2P_Scripts
% 1. organize one or more sessions (edit the config block first)
run_organize
%    -> session_pipeline/organized_data/<session>_organized.mat
%    -> session_pipeline/organized_data/QC_<session>.png   (INSPECT THIS)

% 2. explore a single session
explore_session('session_pipeline/organized_data/AnimalRB19_260320_1444_organized.mat')

% 3. compare a pair (thin compute layer)
R = compare_sessions( ...
      'session_pipeline/organized_data/AnimalRB19_260320_1444_organized.mat', ...
      'session_pipeline/organized_data/AnimalRB19_260320_1609_organized.mat', ...
      'roi_matching/roiMatch_plane0_0320_1444_1609_V1.mat');
compare_plots('violins', R);   % or 'modindex' / 'histograms' / 'scatter' / 'summary'

% 3b. OR the exploratory notebook: open compare_notebook.m, edit the CONFIG
%     cell (sess1/sess2/roi_match/metric), run cells with Ctrl+Enter.

% 4. GROUP: several pairs averaged per drug condition
%    edit group_config.m (the cohort), then open group_notebook.m and run cells
G = compare_group(group_config());
group_plots('traces', G);              % averaged across recordings
group_plots('superplot', G);           % pooled cells + per-recording medians
group_plots('forest', G);              % is the effect consistent across recordings?
T = group_plots('summary', G);
group_travelingwave(G, 'condition', 'yohimbine');
```

**Verification gates:**
1. `run_organize` prints a per-session QC table; every Tier-1 assert must pass (halts otherwise) and `qc_passed` must be `true`.
2. `QC_<session>.png` must show a consistent-latency response band per stimulus (timing sanity).
3. `compare_sessions` must reproduce the previously-validated p-values as a regression check that the cached path matches the live pipeline. On the **`1115_vs_1243`** pair (913 matched): moving_bar z `p≈2.0e-47`, grating z `p≈5.5e-31` (bit-exact). On **`1444_vs_1609`** (616 matched): moving_bar z `p≈6.1e-33`, grating z `p≈1.7e-12`, with matched-cell centroid corr ≈ 1.0 confirming the ROI pairing (Step-7 walkthrough).

---

## Session findings, conventions & gotchas — verification log (2026-07-29 → 08-02)

**Read this section to continue in a fresh chat.** Everything below was verified live in MATLAB R2024b against real Z: data (`AnimalRB19_260320_1444`/`1609`, `260312_1115`/`1243`) or ground-truthed against the presentation source `present_stimulus_addnew.m`. **Status: the pipeline is complete — walkthrough Steps 1–7 all confirmed, the cross-session comparison layer + notebook built and proven correct. The one open item is the "anticipatory-activity" research tests (below).**

### Timing conventions (load-bearing facts)
- **Sampling = 7.5 Hz per plane** for AnimalRB19 (3 functional planes): 30 Hz resonant scan ÷ 4 z-positions (3 planes + 1 flyback). `dt=0.1332 s`. A *2-plane* recording would be 10 Hz (30/3) — that's the source of the remembered "10 Hz".
- **Use each plane's OWN `TimeCa` row** (`TimeCa(selected_plane_idx,:)`), not row 1. Sub-volume offset ~33 ms/plane (66 ms for plane 3); only matters for planes 2/3. `organize_session` already does this.
- **`TimeStimulusFrame` and `Frameinfo.*` are 60 Hz DISPLAY-monitor frames, NOT calcium frames** (display:calcium ratio = 8.00 = 60/7.5). Every paradigm uses TSF correctly: blocks take its first/last as time *bounds*; direction reads TSF at the display indices where `subtrial_frame_count==1`. Nothing assumes TSF is 1:1 with calcium frames.

### Per-paradigm onset provenance (how each stimulus is aligned)
- **looming** (trial_onset): onsets = hardware red pulses `Triggers.TimeProjector` inside the stim window. 5 trials.
- **sparse_local_global_flashes** (trial_onset): onsets = **RNG replay** (`reconstructFlashes`→`get_flash_onsets`), anchored at `TimeStimulusFrame(1)`. **No hardware per-flash marker** (red frames are 5th-frame sync pulses). Count is **RNG-determined, varies by session** (17 on 1444, 15 on 1115) — not fixed. `s.randomseed` is a saved **twister state struct** (`Seed=0`), not a scalar (`rng(s.randomseed)` restores state; can't print with `%d`). Generator loops a 200-display-frame-padded window then keeps onsets ≤ `stimulus_trial_t`; on 1444, 18 wait-draws → 17 onsets → **0 dropped**. Gaps ∈ [5,15] s. Now guarded by `qc_session`'s `flash_timing` check.
- **grating / moving_bar** (direction): onsets = recorded subtrial starts (`Frameinfo.subtrial_frame_count==1` → `TimeStimulusFrame`). Direction order = `rng(randomseed)`+`randperm` — **ground-truthed** to `present_stimulus_addnew.m:2047` (bar) / `:2120` (grating). `reset_rng_each_trial=true` ⇒ same order every outer trial. **grating `trials=1`** (16 orientations shown once) ⇒ cross-repeat reconstruction is **N/A, not FAIL** (moving_bar's 3 repeats validate: r_same 0.18 > r_shifted −0.10). Subtrials: moving_bar 24 (8×3), grating 16 (16×1).
- **spontaneous / checkers2** (block): one continuous 300 s block, no onsets/RNG/chunking (`stimulus_trial_t=300`). spontaneous is internally 60×5 s display epochs, but the block is time-windowed on the calcium clock so those epochs correctly don't chunk it.

### dF/F & the "anticipatory activity" question (STILL OPEN)
- dF/F (`estimate_continuous_dff`, settings in `newAggregation/save_preprocessed.m`): **acausal** centered ±15 s window, **8th-percentile** baseline, windowed-median normalization, neuropil factor 0.5.
- Flashes show a small pre-onset rise. Ruled out: **NOT** the acausal baseline (causal recompute identical, +0.084 vs +0.086), **NOT** a plane offset, **NOT** misalignment. Best explanation: **7.5 Hz onset edge-smearing + sub-frame jitter** (benign; visibly fainter in 1115 than 1444, consistent with a quantization artifact). Only genuine-neural candidate: **hazard-rate temporal expectation** (bounded 5–15 s ISI).
- **Decisive tests OFFERED, NOT yet run** (the remaining next step): (1) onset-shuffle null, (2) jitter simulation, (3) hazard-rate / elapsed-time split.

### Cross-session comparison (verified trustworthy — Step 7)
- Matching via `roiMatchData.allSessionMapping`: **per-plane-local** suite2p indices used directly as dff row indices (col1=session1/baseline, col2=session2/drug). 1:1. **1444/1609=616 matched, 1115/1243=913 matched.**
- **Pairs proven to be the same physical neurons**: matched-cell centroid distance ≪ shuffled (14 px vs 264 px), centroid corr ≈ 1.0 (the ~14 px is a rigid FOV shift that roiMatchPub registered).
- **Z-scores** (all recompute bit-exact by hand): trial_onset = per-trial `(peak−base)/base_std` then median; direction = same on the preferred-direction trace (pref from baseline, carried to drug, common directions only); block = **population-relative** `(cell_mean−pop_mean)/pop_std` (self-z is ~0 by construction, so deliberately different — sensitive to *ranking* changes, blind to uniform shifts → use `raw` for those).
- **Modulation index** `(drug−base)/(|drug|+|base|)` ∈ [−1,1]; **drops NaN/zero-denominator cells**, so the MI vector is shorter than n_matched (treat as its own distribution).
- Result 1444(baseline)/1609(drug): drug **suppresses** grating (MI −0.60, p=1.7e-12) & moving_bar (−0.60, p=6e-33), **enhances** looming (+0.24, p=3.8e-8); spontaneous/checkers2 unchanged (internal control).

### Gotchas fixed this session (so you don't re-hit them)
- **MATLAB batch:** `run('$VAR/script.m')` mangles forward-slash paths → use `matlab -batch "cd('C:\abs\path'); scriptname"`.
- **grating perpetual QC-FAIL** → `verify.applicable` flag; NaN reconstruction = N/A not FAIL.
- **block truncation:** `time_end = max(stimulus_trial_t*trials, max(TimeStimulusFrame))` — the nominal product omits the per-direction multiplier for moving_bar/grating (would truncate to ~1/n_dir), so always also take the recorded TSF span.
- **roiMatch filename:** `roiMatch_plane0_0320_1444_1609_V1.mat` (convention `roiMatch_plane<P>_<date>_<t1>_<t2>_V1.mat`). The old `roiMatch_plane0_1444_V1.mat` still hard-coded in `Preprocess_Stimulus_Data.m:16` is **stale**.
- **violin outlier squashing** (grating/moving_bar z reach ~75) → `compare_plots` clips axes to robust [1,99] pct (display only; median/IQR/p unaffected).
- **`heatmap_sorted` duplication** — a local copy in `explore_session.m` shadowed the shared file; now one `session_pipeline/heatmap_sorted.m` used by both `explore_session` and `qc_report_figure`.

### Group layer — cohort, planes & stats (added 2026-08-02)
- **The cache is keyed by session AND plane**: `<session>_plane<K-1>_organized.mat` (`organized_path.m`; the number is the suite2p 0-based plane, so it reads like the roiMatch filenames). `run_organize`'s config is now `{file, plane_idx}` rows with `skip_existing`.
- **ROI-match files carry NO session provenance** — only `allSessionMapping` + mean frames. The filename is the only record of which sessions/plane, and it can be wrong or insufficient. So `compare_group` runs a **falsifiable match check per pair**: matched-cell centroid distance vs a shuffled re-pairing (`meta.roi_centroid`, stored by `organize_session`). Verified cohort: 4.5–14.4 px matched vs 282–300 px shuffled, corr(x,y) = 1.000. The check is **specific** — of the 9 plane combinations for 0310, exactly one passes.
- **Planes can differ between the two sessions of a pair.** 0310 was matched **1246 plane2 ↔ 1341 plane0** on purpose (z-drift); the filename says `plane1` and is simply not able to express that. Use `plane_base` / `plane_drug` in `group_config.m`. Applying the mapping at the wrong planes does **not** error — it silently pairs unrelated cells and still returns tiny p-values, which is exactly what the match check exists to catch.
- **Two aggregations, always reported together**: pooled over cells (power; cells nested in sessions) and per-recording medians (`across_pair`). With n = 2 recordings per drug the across-recording signrank **cannot** reach significance (min p = 0.5) — it is reported but flagged `underpowered`; use `consistent_sign` and the forest plot to judge reproducibility.
- **All 4 pairs are the same animal (AnimalRB19) on different days** — n is sessions, not animals.
- Stimulus sets differ by date: `0305_1327` (baseline) has **no moving_bar**, so guanfacine moving_bar rests on one recording; `0312` has no spontaneous. `compare_group` intersects per pair and records `n_pairs` per stimulus.

### Data & regeneration
- Aggregated inputs: `Z:\joeschgrp\Group Members\Rima\Aggregated\<session>_preprocessed.mat`.
- Cached organized structs: `session_pipeline/organized_data/<session>_plane<P>_organized.mat` (gitignored, ~100–140 MB; regenerate via `run_organize`). Currently cached: 1444/1609 (plane0 + plane2), 1115/1243 (plane0 + plane1), 1327/1409 (plane0), 1246 (plane2), 1341 (plane0).
- Group outputs: `session_pipeline/group_figures/` (`G_group.mat` + PNGs + `group_summary_z.csv`).
- ROI-match files: `roi_matching/roiMatch_plane0_0320_1444_1609_V1.mat`, `roiMatch_plane0_0312_1115_1243_V1.mat` (+ plane1/plane2 variants).
- **Ground-truth presentation scripts ARE readable** in `C:\Users\rbondarenko\Downloads\` (`present_stimulus_addnew.m`, `present_stimulus.m`, `DataGathering.m`, `save_preprocessed.m`) — use them to close any future timing question.

### To extend
- **New stimulus type:** add a case to `stim_plan()` in `organize_session.m`, then one branch each in `organize_session` / `qc_session` / `qc_report_figure` / `explore_session` / `compare_sessions` / `compare_plots`. The `switch E.paradigm` shape is identical in all of them.
- **New comparison view:** one local function in `compare_plots.m` + one `%%` cell in `compare_notebook.m`.
- **QC now has 31 checks** (Tier-1 fatal asserts + Tier-2 flags incl. `flash_timing`); fault-injection-tested (halts on integrity breaks, flags on validity breaks).
