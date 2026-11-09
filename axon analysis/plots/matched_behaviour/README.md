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
| square cue 0-3 s, spout + reward | -0.045 | **-0.021** [-0.038, -0.003] | 4 / 1 | 0.04 / 0.04 | small dip, half of the naive dip is behaviour; one animal; 0/4 sessions above benchmark -> ambiguous |
| square 3-10 s, spout + reward | -0.002 | -0.009 [-0.025, +0.006] | 4 / 1 | 0.30 / 0.65 | nothing |
| square cue 0-3 s, no spout | -0.002 | -0.017 [-0.045, +0.002] | 3 / 2 | 0.45 / 0.64 | nothing |
| square 3-10 s, no spout | | not testable | 0 / - | | too few control windows (trial period ~22 s; 27-41% of the windows matched in the unflagged sessions) |
| **task - passive, cue** | | **+0.002** (exact permutation p 0.89, smallest possible 0.029) | 4 vs 3 | | **no spout / no-spout difference** at matched behaviour |
| rewarded - spontaneous lick bout | +0.052 | **+0.002** [-0.025, +0.016] | 3 / 1 | 0.88 / 0.65 | the whole raw difference is lick movement; no reward signal beyond behaviour |
| hit - miss trial | +0.10 | not testable | | | hit trials contain licking, misses do not: only 8-39% of the windows have a match in the unflagged sessions |
| flash | +0.058 | **+0.036** [+0.002, +0.082] | 4 / 2 | 0.003 / 0.03 | small response survives; 3 of 4 sessions positive (RB24 1630 / 1538 / 1534; RB22 1138 ~0) |
| expanding loom | +0.066 | +0.035 [-0.008, +0.109] | 4 / 2 | 0.04 / 0.035 | driven by RB24 1538 (+0.148) and 1534 (+0.053) |
| receding loom | +0.045 | +0.043 [-0.007, +0.112] | 3 / 2 | 0.05 / 0.025 | same size as expanding: not direction (threat) specific |
| moving bar | +0.001 | +0.014 [-0.004, +0.036] | 5 / 3 | 0.04 / 0.21 | nothing robust |
| expanding - receding | | not testable | | | |

Short answer: **most of what looks like an event-related axon difference is behaviour** (rewarded-vs-spontaneous lick bout +0.052 -> +0.002; loom +0.066 -> +0.035; task cue -0.045 -> -0.021).
What survives is small (about +0.035 dF/F, roughly 1 SD of a window mean for a single session): a visually evoked component for flashes and looms (not shape-specific, mostly RB24 1538 / 1534 / 1630) and a small cue dip in the lick-task sessions (one animal).
No difference between spout and no-spout sessions, and no reward / outcome signal beyond licking movement.
These are 3-4 sessions in 1-2 animals each: nothing here is established, the flash / loom component is the best candidate for a replication.

## Caveats
* Cross-session level comparisons (task vs passive "baseline" activity) are not possible with this dF/F definition (running 10th-percentile baseline +-60 s, 30 s high-pass absorbs level shifts; raw pupil / face units differ between sessions, e.g. pupil median 800-1480 px). Only within-session event contrasts are compared between conditions.
* Flash coverage is 65% (the rest of the flashes start with movement that no control window has): the effect applies to the matched flashes.
* Control windows are on a 1 s grid and may overlap each other (different pairs can share a control); this is inside the nulls because the same procedure is repeated, variant `step5` (controls >= 5 s apart) is the check (M4): same estimates (task cue -0.024, flash +0.035, expanding loom +0.069, receding loom +0.046, rewarded - spontaneous +0.004) but fewer testable sessions (the control pool shrinks: flash 2, passive cue 1).
* Flagged sessions (RB23 1306 / 1704, RB24 1527 / 1640 / 11-06 sessions) are drawn hollow and never pooled. RB23 1704 (flagged: inconsistent behaviour channels) shows rewarded - spontaneous +0.038 (p 0.003) - not interpreted.
* The ball is usable in about half of the sessions only (`loco_ok`), so locomotion is matched there only.

## Figures
`M0_overview` (naive vs regression vs matched, every contrast), `M1_method_validation` (example, feature space, balance, simulations), `M2_spout_vs_no_spout` (matched averages, per session, task vs passive), `M3_other_cases`, `M4_robustness` (17 settings, one changed at a time), `M5_single_roi`.
Tables: `matched_results.csv` (every session x contrast: counts, naive / reg / matched / adj, both nulls, balance), `matched_pooled.csv`.
