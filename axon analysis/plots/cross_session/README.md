# Cross-session LC-axon analysis (RB22-24), first pass: analyses 1-4 (start with SUMMARY_what_the_data_support.png)

Everything here is **exploratory**. 19 aggregated sessions of RB22 / RB23 / RB24 went through ONE standardised pipeline (code in `axon analysis/cross_session/`).
Common 5 Hz grid, dF/F = (Fc - F0)/F0 (Fc = F - 0.5 Fneu, F0 = running 10th percentile +-60 s), population signal = mean of the signal-like ROIs (V3 criterion) of a session.
Controls: the noise-like ROIs the experimenter left in on purpose (amber), and random sets of REJECTED suite2p ROIs of the same size (grey). Nulls are circular shifts (state analysis) or onsets shifted inside their own block (stimulus analysis).

## Analysis 1 - the shared behavioural-state signal  (S*.png)
| file | question | how to read it |
|---|---|---|
| S1_session_overview_coherence | is there one shared signal, more than in random non-accepted ROIs? | blue dot right of the grey range = more shared structure than random rejected sets |
| S2_coupling_heatmaps | how strongly does it follow pupil / face / whisking / eyelid / ball / brain shift? | peak r within +-10 s; * = p<0.01 vs lags >= 120 s away; compare the 3 panels: amber / grey as red as blue = common-mode |
| S3_lead_lag | who leads whom? | face, whisking, ball: lag ~ 0; pupil follows by ~1 s (slow pupil kinetics, not proof of a lead) |
| S4_sign_control | do all ROIs go the same way? | blue near +1 and grey near 0 = coherent signal vs mixed-sign artefact pattern |
| S5_variance_explained | how much does behaviour explain, and is it specific to axons? | blue vs amber / grey R2 |
| S6_examples | what does it look like? | 5 min excerpts |
| S7_states | size of the state effect | active - quiet, running - still, in dF/F units |
`session_qc.csv` has the numbers (and why a session is or is not informative: behaviour consistency, ball plausibility, ROI counts).

## Analysis 2 - stimulus responses  (L*.png, stim_summary.csv)
| file | question |
|---|---|
| L1_trial_averages | trial-averaged population dF/F and behaviour per stimulus (flash, expanding / receding loom, bar, square) |
| L2_amplitudes_raw_and_behaviour_adjusted | per-session z vs null: raw, after removing face + whisking + locomotion, after also removing pupil (secondary: pupil carries the light reflex) |
| L3_expanding_vs_receding | is the response specific to expanding discs? (within-session pairs) |
| L4_trial_relations_loom_expand | does a bigger response go with a bigger pupil / face change? |

## Cautions
- 11-06 sessions (1529, 1317): behaviour has little dynamic range (face-motion SD 0.17-0.21 vs 0.4-0.8) and the ball reads 'running' 80% of the time: weak coupling there is not necessarily absent coupling.
- RB23_261027_1306 and RB23_261029_1704: face / whisking / pupil disagree with each other (face-whisk r 0.09 / 0.13 vs 0.85-0.96): tracking problem, couplings not interpretable.
- 1527 (6 min) and 1640 (7 min): the circular null has few independent lags, p-values are weak.
- Ball locomotion is used only where it covers >= 70% of the time and reads 'running' < 70% (1744, 1640, 1651, 1423, 1520).
- Curation of accepted ROIs was not blind to the traces: the sign test (S4) supports a coherent signal but does not by itself exclude arousal-locked axial (z) motion, there is no reference channel.

## Analysis 3 - the recorded grey before / after the looming blocks  (G*.png, grey_windows.csv)  [code: cross_session/grey_compare.py, grey_figures.py]
6 sessions (RB24 1538, 1534, 1604; RB23 1642, 1306; RB22 1138), 21 loom blocks with a RECORDED 120 s grey before and after (block >= 121 s from its neighbours). Windows = 100 s each, 10 s guard.
| file | question |
|---|---|
| G1_grey_timecourse | the 120 s before / after every loom block: signal-like axons, noise-like ROIs, rejected ROIs, face, whisking, pupil, eyelid (z within session, mean over blocks) |
| G2_grey_pre_post_paired | block by block: before vs after (axon mean / SD / event rate / time active, movement-removed mean, behaviour, axon-face r) |
| G3_grey_contrasts | after - before across a loom block vs next-block's-before - this-after across PLAIN grey (elapsed-time control) vs their difference (the only contrast not confounded with time in session) |
| G4_grey_vs_other_epochs | same signal in spontaneous block / grey / loom blocks / other stimulus blocks (SD, coupling to face motion) |
| G5_grey_windows_timeline | every window in session time |
Cautions: no no-looming control block exists, so a loom effect cannot be separated from time in session except through the plain-grey control (only 5 sessions have consecutive loom blocks); sessions 1306 (4 signal-like ROIs) and 1534 (8) are thin; 1604 only has 2 bracketed blocks (the other 4 have no grey).

## Analysis 4 - the lick task as context / engagement  (T*.png, task_summary.csv)  [code: cross_session/task_engagement.py, task_figures.py; trial + lick export from the task script: 2P_cache/xs/*_task.mat]
5 RB24 sessions with a spout (1752, 1744, 1423, 1520, 1317) + RB23_1704 (shown, never pooled: behaviour channels inconsistent).
| file | question |
|---|---|
| T1_task_timeline | hits / misses, lick rate, axon SD / mean, pupil, face motion over every session (breaks = white gaps) |
| T2_task_subsessions | engagement (hit rate, latency, time-outs, licks / min) and axon signal per subsession |
| T3_task_epochs | wait / square / post-square (hit and miss trials) and breaks: axons, axons with movement removed, face motion, lick rate |
| T4_task_engagement_bins | lick rate per minute vs axon / behaviour per minute (time trend removed, circular-shift null) |
| T5_task_pre_trial_state | does the state 6-1 s before the square predict hit, first-lick latency, licks in the wait? |
| T6_task_hit_vs_miss | hit vs miss trials aligned to the square, with and without the concurrent movement removed |
Cautions: pupil dilates in the breaks (different screen), so it is not an arousal read-out there; 1752 has camera gaps inside 120 trials (the movement-removed version uses fewer trials); hit trials contain licking, miss trials do not, so every raw hit - miss difference is also a movement difference; 11-06 1317 has no behaviour coupling at all (movement-removed = raw); circular-shift p-values over minute bins have a floor of 1/(n bins - 4).
