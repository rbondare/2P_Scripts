# Factor model of the accepted ROIs (first version, 2026-11-12)

Instead of merging ROIs, fit a few shared factors across all accepted ROIs of a session: `z_i(t) = sum_k lambda_ik f_k(t) + e_i(t)` (factor analysis with a noise variance per ROI). Every ROI keeps its own
loading (amplitude); the number of shared signals is measured by cross-validation; nothing is decided by a correlation threshold. Uses the `axon_groups` loader (`ag_data.Universe`, `ag_pairs`).

```
miniforge3\envs\suite2p_venv\python.exe fm_run.py [--domain band|dff] [session ...]   # fit + printed summary, pickles in <2P_cache>/axon_factors/
miniforge3\envs\suite2p_venv\python.exe fm_cross.py [band|dff] [sessionA sessionB]     # do the loadings replicate across two sessions of the same field of view (default 1423 / 1520)
miniforge3\envs\suite2p_venv\python.exe fm_figures.py [--domain band] [session ...]     # plots/axon_factors/<session>_factors_<domain>.png
```
`band` = 1-30 s band-passed best dF/F (the grouping / coupling domain); `dff` = the unfiltered 5 Hz best dF/F (keeps everything the 2-min baseline leaves, i.e. up to ~2 min).

## What `analyse()` does (fm_fit.py)
1. accepted ROIs, valid frames (registration excursion > 0.93 um and missing bins dropped), every row z-scored (SD kept)
2. **rank**: leave-block-out (8 blocks, 60 s buffer) AND leave-ROI-out (the ROI and every ROI within 2 um are predicted from the factor scores of the others): held-out R2 per ROI for k = 0..8;
   k* = last k whose gain in the median R2 of the signal-like ROIs is >= 0.01
3. full fit at k*, scores, share of the standardised variance, varimax loadings for the maps
4. what each factor is: correlation with behaviours (face, whisking, pupil, ball), ridge R2 on lagged behaviours (held-out half), correlation with nuisance signals (rigid shift, mean neuropil of the accepted ROIs,
   PC1 of the rejected ROIs that do not touch an accepted ROI)
5. stability: principal angles between the loading subspaces of the two session halves vs random subspaces
6. residual: leave-one-ROI-out residual of every accepted ROI, skewness and the V3 asymmetry test (upward transients) vs the strict rejected noise pool
7. every rejected ROI regressed on the factor scores (held-out R2, and the part added by factors >= 2)

## First results (5 sessions; see the reply of 2026-11-12 for the reading)
| session | k* | variance share per factor (all ROIs) | split-half cosines | notes |
|---|---|---|---|---|
| RB24 1630 | 2 | 0.72, 0.03 | 0.99, 0.91 | factor 1 = state (face r 0.83, behaviour R2 0.72, also r 0.74 with rejected-ROI PC1); factor 2 uncorrelated with everything tested, but bipolar (minority sign 0.41) |
| RB22 1138 | 1 | 0.44 | 0.94 | factor 1 = PC1, r 0.96 with rejected-ROI PC1 (shared with non-axon ROIs); R2 keeps rising slowly to 0.77 at k = 8 |
| RB23 1642 | 6 | 0.22 0.13 0.10 0.09 0.05 0.03 | 0.96 0.94 0.92 0.88 0.38 0.01 | no factor related to behaviour or nuisance (RB23 behaviour is unreliable); last two unstable |
| RB24 1423 | 4 | 0.40 0.15 0.13 0.05 | 0.99 0.99 0.92 0.65 | f1 state; f2, f3 one-signed, stable, weakly / not behaviour related; f4 r -0.76 with rejected-ROI PC1, r -0.44 ball: nuisance-like |
| RB24 1520 | 3 | 0.29 0.15 0.09 | 0.96 0.94 0.60 | f1 state; f2 one-signed, stable, uncorrelated with anything tested; f3 r -0.86 with rejected-ROI PC1, r -0.65 ball: nuisance-like |

1423 <-> 1520 (40 common ROIs through the master axons): loading-subspace cosines 0.98, 0.94, 0.15 vs permutation p95 0.93, 0.45, 0.21; the nuisance-like third direction does not replicate.
Slow domain (`dff`) gives the same structure (k* 4 / 4 / 2) and a more specific residual test: upward transients remain in 20-43 % of the signal-like ROIs after the factors are removed vs 0 of 104 / 108 / 151 noise-pool blobs.

## Limits
Accepted (hand-curated) ROIs only; a factor is a shared time course, not a neuron; "unrelated to behaviour" only means unrelated to the measured, linear, lagged regressors; axial (z) motion is not measured
(rigid shift is lateral only) and gives bipolar "mirror" loadings, so bipolar factors are suspect; k* uses a conservative plateau rule; convergence warnings occur in the cross-validation folds at high k (the
final fits are converged: cosine 1.0 against fits run 20-200x longer); same-day replication shows stable structure, not axon identity.
