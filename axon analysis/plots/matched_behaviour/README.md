# Behaviour-matched comparison of the LC-axon signal (2026-11-09)

Question: is the axon signal only a read-out of behaviour, or is there activity left when the behaviour is the same?
Method: every event window (square cue, flash, loom, lick bout, ...) is paired with control windows **of the same session** that have the same face motion, whisking, ball speed
(where the ball is usable) and the same pupil state, and the axon activity of the two windows is compared. No regression, no assumed functional form.

## Code (all in `axon analysis/cross_session/`, python = `miniforge3\envs\suite2p_venv`)
| file | what it does |
|---|---|
| `matched_behaviour.py` | the engine: windows, matching, estimates, shift null, behaviour-only benchmark. `python matched_behaviour.py [session ...]` -> `cache/<session>.pkl`. `MB_VARIANT=<name>` runs one robustness variant into `cache_<name>` (see `VARIANTS`). Needs `C:\Users\rbondarenko\2P_cache\xs\<session>_xs.mat` (from `export_xs.m`) and, for the task contrasts, `<session>_task.mat`. |
| `matched_validate.py` | simulation check on synthetic axons built from the REAL behaviour and event times (known truth) -> `validation.pkl` |
| `matched_pool.py` | pooling over sessions (session = unit, animal = cluster, unflagged + testable only), group difference with exact permutation |
| `matched_figures.py` | all figures + `matched_results.csv` (every session x contrast) + `matched_pooled.csv` (pooled, dF/F and in SD units) |

To run a new session: export it with `export_xs.m`, add its name to `SESSIONS` in `matched_behaviour.py`, run `matched_behaviour.py <session>`, then `matched_figures.py`.

## How a match is made
* Window = response window [a, b] after the event (cue 0-3 s, 3-10 s, flash 0-2 s, loom 0-4 s, bar 0.5-4 s, lick bout 0-4 s, hit/miss 0-6 s).
* Axon outcome = mean of the 30 s high-passed best dF/F of the signal-like ROIs in the window (every ROI separately for the ROI figure).
* Behaviour features = mean face motion, whisking, ball speed (only where the ball is plausible) over the same window shifted by that behaviour's own lag to the axon (cross-correlation peak of the session), plus the **pupil state 5-1 s BEFORE the window**. The pupil inside the window is deliberately not used for visual stimuli: the light reflex moves it out of the range of any control window (matching then fails or extrapolates - variant `pupilin` shows it: bias 0.02 dF/F, flash "excess" +0.22).
* Controls: windows of the same session on a 1 s grid that do not touch [-3, +12] s around any stimulus onset (class contrasts: the other class). 1:2 optimal matching without replacement, Mahalanobis distance, **caliper 0.5 SD on every feature**, SD = spread of the feature over all windows of the session. Event windows without a partner are dropped and counted (coverage = fraction matched; the effect refers to the matched subset). A session x contrast counts as *testable* only if >= 50% of its windows matched and >= 10 pairs exist.
* Estimates (dF/F, mean over event windows): `naive` = event - all controls; `reg` = event - linear behaviour prediction (no matching); `matched` = event - matched controls; `adj` = matched, bias-corrected for the residual imbalance.
* Two nulls: (1) **shift null** - all design events shifted together by >= 60 s against the data and the whole pipeline re-run (autocorrelation, block structure, matching bias inside); (2) **behaviour-only benchmark** - the same events and pipeline on 200 surrogate axon traces = linear lagged behaviour prediction (trained outside the stimulus windows) + phase-randomised residual, i.e. what the pipeline returns if the axon were a linear function of behaviour. A real effect must exceed both.

## Does the method work? (`M1_method_validation`, `validation.pkl`; 4 sessions, 7 contrasts, 60 surrogates each)
* Axon = behaviour only (linear or expansive non-linear): naive event-minus-control is biased by **+0.17 / +0.20 SD** of a window mean (up to +0.5) and rejects 23-28% at p < 0.05; matched bias **0.00 / +0.01 SD** (worst case 0.14), false positives 7-8%. With a noisy behaviour proxy (measurement noise 0.5 / 1 SD) matched stays within 0.04 SD on average (false positives 6 / 10%, worst session +0.03 dF/F): residual confounding exists but is small.
* A known event response is recovered with slope **1.00** (naive 1.22), power 81% at 0.5 SD of a window mean (42% at 0.25 SD, 99% at 1 SD) with 30-120 events; 10-15 events (looms) have about 35-60% power at 0.5 SD.
* Balance (44 testable session x contrasts): the largest feature imbalance is 0.05 SD after matching (median; 90th percentile 0.14, worst 0.28) versus 0.21 SD before (worst 0.72).
* Real limit: matching removes **measured** behaviour only. Axial (z) movement, neuropil and any state that face / whisker / ball / pupil do not carry stay in the "difference".

## Results (unflagged sessions, animal-clustered 95% CI, `M0_overview`, `matched_pooled.csv`)
| contrast | naive | matched [95% CI] | sessions / animals | p shift null / benchmark | reading |
|---|---|---|---|---|---|
| square cue 0-3 s, spout + reward | -0.029 | -0.016 [-0.040, +0.010] | 5 / 1 | 0.09 / 0.10 | small, 3 of 5 sessions negative (RB24 1744 -0.033 p 0.03, 1520 -0.060 p 0.02; 1506 +0.025); one animal |
| square 3-10 s, spout + reward | +0.003 | -0.001 [-0.025, +0.021] | 5 / 1 | 0.99 / 0.95 | nothing |
| square cue 0-3 s, no spout | -0.002 | -0.011 [-0.030, +0.005] | 3 / 2 | 0.57 / 0.71 | same small dip as with the spout, n.s. |
| square 3-10 s, no spout | | not testable | 0 / - | | too few control windows (trial period ~22 s; 27-41% of the windows matched in the unflagged sessions) |
| **task - passive, cue** | | **+0.000** (exact permutation p 0.96, smallest possible 0.018) | 5 vs 3 | | **no spout / no-spout difference** at matched behaviour |
| rewarded - spontaneous lick bout | +0.045 | **+0.011** [-0.003, +0.023] | 4 / 1 | 0.23 / 0.42 | most of the raw difference is lick movement (1744: +0.108 -> +0.024); no reward signal beyond behaviour |
| hit - miss trial | +0.07 (1506: -0.057, the others +0.06 to +0.13) | not testable | 0 | | hit trials contain licking, misses do not: only 8-58% of the windows have a match; no unflagged session reaches 10 pairs with a valid null |
| flash | +0.058 | **+0.030** [-0.002, +0.073] | 4 / 2 | 0.010 / 0.040 | small response survives in RB24 1538 (+0.084) and 1534 (+0.076), weak in 1630 (+0.025), none in RB22 1138 (-0.002); RB23 1642 not testable (15/38 matched) |
| expanding loom | +0.091 | +0.047 [-0.013, +0.119] | 3 / 2 | 0.056 / 0.040 | driven by RB24 1538 (+0.130, p 0.04 / 0.009) and 1630 (+0.084 n.s.); RB23 1642 -0.013 |
| receding loom | +0.045 | +0.042 [+0.001, +0.104] | 3 / 2 | 0.030 / 0.015 | same size as expanding: not direction (threat) specific; RB24 1538 +0.125, 1534 +0.038, RB23 1642 +0.001 |
| moving bar | +0.001 | +0.002 [-0.008, +0.011] | 5 / 3 | 0.72 / 0.83 | nothing |
| expanding - receding | | not testable | | | |

Short answer: **most of what looks like an event-related axon difference is behaviour** (rewarded-vs-spontaneous lick bout +0.045 -> +0.011; expanding loom +0.091 -> +0.047; task cue -0.029 -> -0.016).
What survives is small (about +0.03 to +0.05 dF/F, roughly 1 SD of a window mean for a single session): a visually evoked component for flashes and looms that is clear in only 2 sessions of one animal (RB24 1538 / 1534; weaker in 1630, absent in RB22 1138 and RB23 1642) and not shape-specific, plus a small cue dip (-0.016 task, -0.011 no spout) that is the same with and without the spout.
No difference between spout and no-spout sessions (task - passive cue +0.000, p 0.96), and no reward / outcome signal beyond licking movement.
These are 3-5 sessions in 1-2 animals each and the CIs of flash / expanding loom now touch zero: nothing here is established, the flash / loom component is the best candidate for a replication.

## Caveats
* Cross-session level comparisons (task vs passive "baseline" activity) are not possible with this dF/F definition (running 10th-percentile baseline +-60 s, 30 s high-pass absorbs level shifts; raw pupil / face units differ between sessions, e.g. pupil median 800-1480 px). Only within-session event contrasts are compared between conditions.
* Flash coverage is 65% (the rest of the flashes start with movement that no control window has): the effect applies to the matched flashes.
* Control windows are on a 1 s grid and may overlap each other (different pairs can share a control); this is inside the nulls because the same procedure is repeated, variant `step5` (controls >= 5 s apart) is the check (M4): similar estimates (task cue -0.020, flash +0.061, expanding loom +0.049, receding loom +0.045, rewarded - spontaneous +0.010) but fewer testable sessions (the control pool shrinks: flash 2). Dropping the pupil or the ball (`pupilnone`, `facewhisk`, `noloco`) or the lag fit (`nolagfit`) moves the pooled numbers by <= 0.02; the pupil inside the window (`pupilin`) breaks the flash / loom contrasts as expected (expanding loom +0.003).
* Flagged sessions (RB23 1306 / 1704, RB24 1527 / 1640 / 11-06 sessions) are drawn hollow and never pooled. RB23 1704 (flagged: inconsistent behaviour channels) shows rewarded - spontaneous +0.038 (p 0.003) - not interpreted.
* The ball is usable in about half of the sessions only (`loco_ok`), so locomotion is matched there only.

## Data versions (2026-11-11)
All 21 sessions are run from `C:/Users/rbondarenko/2P_cache/xs_v2` (set `MB_XS` to that folder): xs + task exports with the 2026-11-10 whisking index and JawME. The 7 task sessions were written by `task_deep_dive/export_task_events.py`, the other 14 were re-exported with `export_xs.m` on 2026-11-11 (new whisking differs from the old index: r 0.38-0.97 per session). Each pickle stores the folder it used (`xs_dir`).
The first 20-session results (old whisking index, `xs` folder) are kept in `backup_results_xs_v1_old_whisking/` (caches, csv, validation) and `backup_1317_old_whisking/`. Differences old -> new: task cue -0.021 -> -0.016 (the 20-session p 0.04 is gone), flash +0.036 -> +0.030, expanding loom +0.035 -> +0.047 (RB24 1534 dropped below 10 pairs), moving bar +0.014 -> +0.002, task - passive cue +0.002 -> +0.000.
RB24 1317 and the 11-06 sessions stay flagged (little behavioural range); RB23 sessions are drawn but their whisking is uninterpretable.

## Figures
`M0_overview` (naive vs regression vs matched, every contrast), `M1_method_validation` (example, feature space, balance, simulations), `M2_spout_vs_no_spout` (matched averages, per session, task vs passive), `M3_other_cases`, `M4_robustness` (17 settings, one changed at a time), `M5_single_roi`.
Tables: `matched_results.csv` (every session x contrast: counts, naive / reg / matched / adj, both nulls, balance), `matched_pooled.csv`.
