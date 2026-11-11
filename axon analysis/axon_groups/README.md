# Axon grouping layer (pilot)

Groups accepted suite2p ROIs into co-active, spatially close "units", attaches **rejected** ROIs to those units with a graded confidence, and builds a brightness-matched
**noise reference**. It is an **analysis layer only**: `iscell.npy`, `Fall.mat`, the aggregated files and every existing analysis are untouched, nothing is written to Z:.

The design driver: grouping by correlation must neither **erase real differences between ROIs** nor **manufacture homogeneity** (correlated ROIs selected, then averaged).
So the output is not only a grouping but a **measurement of what the grouping does** (report card, below). Group means are meant for later analyses only where those checks pass.

## Run

```
miniforge3\envs\suite2p_venv\python.exe ag_run.py <session> [--level pc1|raw|pc1beh] [--nrep 100] [--validate] [--refresh]
miniforge3\envs\suite2p_venv\python.exe ag_figures.py <session>      # map, unit traces, report card (needs --validate output)
```
Outputs: `C:\Users\rbondarenko\2P_cache\axon_groups\<session>_{roi_groups.csv, groups.pkl, units.npz, validation_<level>.pkl}`; figures in `axon analysis/plots/axon_groups/`.
Universe cache `<session>_U.pkl` (loading from Z: takes 25-60 s); delete it or pass `--refresh` after re-curating `iscell`.

## Modules

| file | step |
|---|---|
| `ag_data.py` | `Universe`: accepted (iscell == 1, xs order verified), rejected, hidden merge constituents (`inmerge > 0`) with footprints; same F/Fneu -> best dF/F recipe for all (`xs_common`); µm per pixel from the recordings DB; bad-registration frames masked (> 0.93 µm shift excursion) |
| `ag_pairs.py` | residual correlation of every ROI pair, three removal levels, Bartlett-sigma z + pooled circular-shift null, footprint gaps, BH-FDR, gap diagnostic |
| `ag_group.py` | greedy constrained agglomeration (units: `D_FRAG` 5 µm; clusters: `D_BUNDLE` 30 µm), homogeneity tests, tiers |
| `ag_attach.py` | rejected ROI -> unit, confidence high / medium / low / none, orphans |
| `ag_noise.py` | noise pool (strict control), brightness matching, top-20 + matched random draws |
| `ag_validate.py` | null runs, homogeneity calibration, ground-truth simulation of two different axons, parameter sweep |
| `ag_report.py` | held-out SNR / coherence, member coupling spread, in-sample vs held-out inflation, units vs matched noise |
| `ag_run.py`, `ag_figures.py` | runner, figures |

## What is computed

* **Residual r.** Band-passed (1-30 s) best dF/F with the population PC1 (of the signal-like accepted ROIs) projected out, because that shared signal is 30-80 % of the variance and makes every
  axon correlate with every other. Level `pc1beh` additionally removes the behaviour-predicted part (face / whisking / pupil / ball, lags -1..+3 s), so that grouping is blind to the readouts.
* **Significance.** z = r / sigma_ij (Bartlett, from both traces' autocorrelations) against a pooled null of circular shifts (>= 120 s); BH over candidate pairs. Edge: FDR < 0.01 and r >= 0.3.
* **Units (the only level that is ever averaged).** Footprints <= 5 µm apart, greedy join of the two groups whose mean residual traces correlate best, guard against bridging (every member of the
  smaller side must correlate significantly with the other side). **Clusters** (annotation only, never averaged): the same on unit traces within 30 µm.
* **Unit trace** = plain mean of the members' original dF/F (not the residual), so the shared signal stays in it.
* **Tiers.** `A` = >= 4 members, stable in both halves (mean within-unit residual r >= 0.15 each), the unit mean passes the V3 transient test, **and** the rank-1 test below passed -> homogeneity demonstrated.
  `B` = multi-member unit that is not A; `tier_why` says why (`untestable_n2`, `untestable_n3`, `heterogeneous`, `unstable`, `v3_fail`). `C` = single ROI passing V3, `D` = single ROI failing V3.
* **Rank-1 test.** Principal-axis one-factor fit of the members' residual correlation matrix; statistic = largest |R_ij - lambda_i lambda_j| / sigma_ij; threshold = 95th percentile of simulated
  homogeneous units of the same size that carry the **session's own noise** (rejected noise blobs). **Needs >= 4 members**: with 2 or 3 the model is exactly identified and nothing can be tested.
* **Attachment.** Rejected ROI vs the trace of each unit within 30 µm built **without members that overlap the ROI**; confidence `high` (q < 1e-3, r >= 0.4, r >= 0.25 in both halves, gap <= 5 µm),
  `medium` (q < 1e-2, r >= 0.3, both halves >= 0.15), `low` (q < 0.05, r >= 0.2). It never changes a unit trace. `orphan` = active on its own (V3 or skew > 0.5), attached to nothing, not touching an accepted ROI.
* **Noise reference.** Rejected, not a merge constituent, not within 2 px of an accepted ROI, fails V3, skew < 0.3 (Fig 12 strict control), not attached at any confidence; matched to the units' members on
  log F0, log npix, Fneu/F (nearest neighbour without replacement) -> top 20, plus 200 random draws from the matching caliper. A unit of n members is compared with the **mean of n** noise blobs.

## CSV columns (`<session>_roi_groups.csv`, one row per ROI)
`suite2p_idx, source (accepted|rejected|constituent), v3_class, unit_id, unit_size, cluster_id, tier, tier_why, attach_unit, attach_conf, attach_r, attach_q, attach_gap_um, orphan, noise_ref_top20, merged_parent_idx`

## Report card (what the pilot shows) - see the sections added after the pilot runs

RESULTS_PLACEHOLDER

## Caveats
* A unit is a co-active set, not proven to be one neuron: branches of one neuron correlate when apart, and two neurons driven by the same events can look alike.
* Removing PC1 (and behaviour) from the grouping makes groups reflect non-shared structure; the group count depends on the removal level (the sweep reports it).
* Selection is blind to behaviour coupling; the noise reference is selected on brightness, transient test and skew only. Accepted ROIs were curated by eye, so a units-vs-noise difference is a property of
  curated ROIs vs uncurated blobs, not a proof about the axons.
* Simulated fragments use independent real noise per fragment; real fragments share some noise (neuropil, motion), which would make the rank-1 test somewhat conservative, not liberal.
* The `shift` mask (rigid registration excursion) removes up to 18 % of the frames in the high-zoom November sessions; olive-green orphan outlines along vessel walls in October sessions are probably
  vessel pulsation, not axons: review before using them.
* Three animals: pool per session first (session = unit, animal = cluster), as in the paper figures.
