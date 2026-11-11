# Lick-task deep dive: RB24_261109_1506 (+ 6 other lick sessions), 2026-11-11

Code: `axon analysis/task_deep_dive/` (Python, `miniforge3\envs\suite2p_venv\python.exe`; nothing committed). Inputs: the aggregated files on Z: (read only) through a new export
(`export_task_events.py` -> `C:\Users\rbondarenko\2P_cache\xs_v2\<session>_task.mat` + `_xs.mat`, now with JawME and the new front_top WhiskIndex).
Outputs here: figures D1-D4 (PNG + PDF), `matched_contrasts.csv`, `fir_kernels.csv`, `hit_vs_miss_windows.csv`, `pooled_contrasts.csv` / `.txt`, `cache/` (pickles), run logs.
The task script `axon_analysis_task_analysis_V1.m/.mlx` was NOT changed. Shared engine files got small additive edits (`cross_session/export_xs.m` exports JawME,
`xs_common.py` reads it, `matched_behaviour.py`: configurable behaviour list + optional within-stratum label permutation null).

## Questions and short answers
1. **Do anticipatory licks (just before the square, sending the animal into a time-out), reward licks and post-reward licks differ in axon activity?**
   Not in a way the movement does not explain. In 1506 no class differs from another (FIR kernels and behaviour-matched bouts, every ROI group, noise-like control).
   Two weak patterns appear across sessions and are listed below as *hypotheses*, not findings.
2. **Hit vs miss with matched behaviour: anything left?** In 1506 no. The raw hit - miss difference (mid window 0.8-3 s: -0.071 [-0.121, -0.023]) is movement:
   after regressing the window's face / whisking / jaw / lick rate / ball it is 0 of 52 signal-like ROIs (raw 52%). The state before the square does not predict the outcome (AUC 0.48).
   The matched hit/miss contrasts are mostly *not testable* (only 33-61% of events find a partner: hits contain licking, misses do not).
3. **Other tests:** silent time-outs, free pump, trial-time level, per-ROI heterogeneity, stability over the session, first-lick vs square alignment. See section 4.

## 1. What the task does in 1506 (figure D1)
- **Time-out rule, read from `timed_Training.m` and verified on the rig's own counts:** at the END of every wait (8-13 s) the rig counts the licks of that wait; more than `maxErrorLicks x wait`
  (0.5 licks/s in 1506) = time-out (silent: `TimeOutSound` 0), a new wait is drawn, the count restarts. Verified on 361 wait segments (waitTimes + errorLicks): all 242 time-out waits exceed the limit,
  none of the 119 waits that ended in a square does; my lick count (2P clock) equals the rig's in 86% and is within +-1 in 99.7%. So a time-out is decided by the lick count of the WHOLE wait,
  not by "a lick just before the square", and nothing marks it for the animal.
- **Licking in the wait is flat:** ~0.85 Hz, no ramp towards the end of a wait (time-out waits 1.0-1.3 Hz, square-ending waits 0.2-0.4 Hz). "Anticipatory" bouts (last 3 s of a wait) are therefore the
  same kind of licking as earlier wait licks; in 1506 134 of them end in a time-out and only 9 before a square (the square comes after sparse licking).
- Bouts (>= 3 licks): early wait 351, anticipatory -> time-out 134, anticipatory -> square 9, reward window 73 (66 rewarded first bouts), post 3-6 s 45, post > 6 s 90. Hit rate 49% (level 2: 40%, level 3: 57%), 119 trials (the 120th has no lick log).
- **Free pump:** level-2 misses (pump 1.5 s after onset by the stated rule) show no extra licking after the pump compared with level-3 misses (no pump). Either the mouse is not responding to the pump or the rule is not what the script assumes (it is still unconfirmed, see memory note on the 1423 reward rule).

## 2. Lick classes (figures D2, D3)
**Classes** (by the lick's own time): early wait; anticipatory -> time-out / -> square (last 3 s before a wait end of that kind); reward window (0.1-3 s after onset; "reward bout" = the bout of the first lick); post 3-6 s; post > 6 s after the square onset.
**Encoding model (FIR, dd_fir.py):** one kernel (-1..+4 s) per lick class + square, pump, silent time-out + movement covariates (jaw, face, whisking, pupil, shift, ball). Held-out R2 of the population 0.27.
**Matched model (dd_matched.py):** every bout paired with bouts of another class whose face / whisking / jaw / lick rate / ball in the window and pupil state before it are the same (caliper 0.5 SD), shift-null and behaviour-only benchmark as in `matched_behaviour.py`.
Scale: SD of the population signal (30 s high-passed) is 0.163, a bout has a median of 4 licks.

1506 (population = 52 signal-like ROIs):
- Mean kernel 0-2 s per lick (x1e-3): early wait -2.0 [-4.5, +0.8], anticipatory -> time-out -4.3 [-7.3, -1.1], anticipatory -> square +2.0, reward window +0.1, post 3-6 s -0.6, post > 6 s +1.3. **No class difference** (all differences n.s.; ant_to - early wait -2.4 [-5.6, +1.0], reward - early wait +2.1 [-3.4, +7.2]).
- Matched contrasts, population: anticipatory vs early wait +0.014 (z +1.2), reward bout vs early wait +0.008, reward vs anticipatory +0.022, post 3-6 s vs early wait -0.010; the PRE-bout window (-2..0 s, before any lick) shows no anticipation signal (-0.014). The only nominal hit is *post > 6 s vs early wait* +0.039 (p 0.003; also in the noise-like control, +0.017, so partly common-mode).
- **Groups differ, not classes:** group 2 (17 ROIs) rises with licks of every class (post > 6 s +8.1e-3 per lick, z 4.7), group 3 (9 ROIs) falls with wait licks (-5.3e-3 early wait, -8.6e-3 anticipatory; z -4.8 / -4.7), groups 1 and 4 lean negative, the noise-like control is ~0. The sign pattern is stable in all four subsessions (D2e). But the per-ROI spread of lick kernels is the same in signal-like and noise-like ROIs (|z| > 3: 37% vs 44%, wait licks) and the sign structure shows only for wait licks (signal-like 15 negative / 4 positive, noise-like 9 / 9). Not seen in the same way in the other sessions (D2f; 1317 shows the same direction for wait licks but also in noise-like ROIs).
- **Pooled over the unflagged lick sessions (`pooled_contrasts.csv`, one animal, 2-5 sessions, bootstrap over sessions):**
  - trial time, 6-10 s after the square onset vs random wait times (behaviour matched, no licking condition): **+0.027 [+0.006, +0.044], 4 of 5 sessions positive** (1506 +0.030, 1423 +0.050, 1744 +0.046, 1752 +0.023, 1520 -0.015), noise-like control ~0 (1506 +0.011). 3-6 s: +0.014.
  - licks > 6 s after the square vs wait licks: +0.021 [-0.003, +0.039], 4 of 5 positive. Because the 6-10 s level is raised without any licking condition, this is most likely **a trial-time (post-square) state**, not a response to those licks (the window BEFORE those bouts is raised too: +0.013 [+0.004, +0.023]).
  - **licks just before a square (anticipatory -> square) vs earlier wait licks of the same kind of wait: -0.023 [-0.035, -0.006], 3 of 3 sessions negative** (1423 -0.027, 1744 -0.035, 1752 -0.006; noise-like control 0 / -0.008 / -0.006); in 1506 and 1520 there are too few such bouts (9). Anticipation with a time-out outcome vs earlier wait licks of the same kind of wait: +0.008 (3/3 positive, ~0). Only a hypothesis: sparse licking before the square and denser licking differ in more than the window's movement (surrounding lick state is not matched).
  - reward bout vs wait lick -0.003 [-0.019, +0.013] (2 of 4 positive); post 3-6 s vs wait lick -0.001; silent time-out 0..3 s after / -3..0 s before vs random wait time -0.000 / +0.003 (n 4, all ~0).

## 3. Hit vs miss (figure D4)
- Raw population (hit - miss, 95% CI): pre -0.016 [-0.068, +0.033], early lick-free 0..0.8 s -0.018, mid 0.8-3 s -0.071 [-0.121, -0.023], late 3-8 s -0.042 [-0.087, +0.001]. In the other sessions the raw mid window is positive in 1744, 1752, 1423, 1520 (+0.07 to +0.24) and slightly negative in 1317 (-0.02); the late window is positive in 1752, 1423, 1520 (+0.08 to +0.13) and negative in 1506 / 1317: the sign changes with the session / reward arrangement.
- Per ROI (hit - miss permuted inside subsessions): mid window 52% of signal-like and 59% of noise-like ROIs p < 0.05 raw, after regressing the window's movement 0% vs 2%. Late window: signal-like 25% -> 33% (all negative, mostly group 2: 76%), noise-like 15% -> 17%: the one place where the signal-like ROIs exceed the control in 1506; in 1317 (flagged) 56% vs 11%, but 1752 has the opposite sign and 1423 / 1520 are ~0, so no consistent effect.
- Behaviour-matched population contrasts (coverage 33-61%, i.e. not testable by the 50% rule): pre +0.026, early +0.047 (19 of 31 matched, z 0.9), late -0.014; none different from the shift null. Pooled: pre / late not testable, early -0.005 (2 sessions).
- Prediction of the outcome from the 5 s before the square (logistic, contiguous folds): movement + history AUC 0.65, + axon 0.62, axon alone 0.48 (dAUC -0.022, p 0.56); in all 7 sessions dAUC between -0.03 and +0.05 (p 0.08-0.69).

## 4. Other tests
- Silent time-out (no cue): no response and no ramp before it (population +0.002 / +0.006; 4 sessions all ~0).
- Free pump vs no pump in lick-free windows (1506 only, 12 of 34 matched): -0.037 (z -1.6), group 2 -0.069 (z -2.5, p 0.017); not significant after the number of tests; the pump rule is unconfirmed.
- Response to the square onset (FIR, population +0.005 [+0.001, +0.009] per event, noise-like +0.002) and the hit-specific part (+0.001): tiny, as in earlier analyses.
- First-lick vs square alignment of hit trials: peak of the average 0.033 (lick) vs 0.023 (square) in 1506, mixed in other sessions; response size vs first-lick latency |rho| <= 0.2 in 1506 (no preference).
- Lick coupling over the four subsessions: group-specific, stable (D2e).

## Caveats
- Exploratory: ~25 contrasts x outcomes were looked at; at p < 0.01 about 1 false hit per 100 is expected. Replication across days is the only arbiter, and all lick sessions are ONE animal (RB24); the only RB23 lick session is flagged. 1317 (11-06: no behaviour coupling, noise-like ROIs left in by design) is flagged and not pooled.
- Matching removes measured movement only (face, whisking, jaw, lick rate, ball, pupil state); z / brain motion, neuropil and unmeasured state stay. In 1506 the ball sync is not verified one-to-one (counts differ by 2) although validated against whisking.
- Whisking / jaw come from the 2026-11-10 front_top / JawME data. WhiskOn is not used (H.264 threshold caveat); the older xs caches of the other analyses are still pre-fix and were not touched (new exports are in `xs_v2`).
- Pump times are the script's rule (first lick + 1.5 s; level 2 without a lick: onset + 1.5 s), unconfirmed against the rig code.
- Groups of 1506 are those of the user's V1 run; groups of other sessions are a Ward tree on 1 - r (4 groups) and not comparable between sessions. Group results are only shown for 1506.
- 1744 / 1752 have only 8 signal-like and 3-4 noise-like ROIs (weak control); the lick log of 1506's last trial and the 242-row cap of one trial (trial 70) are excluded / noted by the task script.

## Rebuild
`python export_task_events.py <session...>` (headless MATLAB, 15-50 s per session) -> `dd_fir.py` (3 min per session), `dd_matched.py` (env `DD_OUTCOMES=P,g1,..`, `DD_NNULL`, `DD_NSUR`, `DD_TAG`), `dd_hitmiss.py`, `dd_extras.py` -> `dd_pool.py`, `dd_figures.py`.
