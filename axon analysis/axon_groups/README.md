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

## Pilot results (2026-11-11, level `pc1`, 5 sessions; `python ag_summary.py` reprints them)

Sessions: RB24 1630 (27 accepted ROIs), RB22 1138 (33), RB23 1642 (56), RB24 1423 (63), RB24 1520 (85; 1423 and 1520 share a field of view).

| session | units / multi-member | sizes of multi-member units | tier A | rejected ROIs attached high / med / low / none | orphans | noise pool (match F0, npix, Fneu/F) |
|---|---|---|---|---|---|---|
| 1630 | 14 / 7 | 5 4 3 2 2 2 2 | 0 | 10 / 46 / 79 / 183 | 65 | 66 (0.11, 0.07, -0.02) |
| 1138 | 27 / 3 | 4 3 2 | 0 | 13 / 116 / 90 / 121 | 16 | 55 (-0.49, 0.06, -0.43) |
| 1642 | 14 / 7 | 14 10 9 8 3 3 2 | 0 | 26 / 52 / 60 / 248 | 58 | 115 (-0.42, -0.29, -0.49) |
| 1423 | 20 / 9 | 13 11 6 6 6 3 3 2 2 | 1 | 39 / 188 / 62 / 104 | 41 | 21 (-1.35, 0.77, -1.66) |
| 1520 | 25 / 7 | 18 16 11 9 9 2 2 | 0 | 63 / 238 / 84 / 50 | 12 | 19 (-0.88, 0.23, -0.97) |

**1. Does the method create groups out of nothing?** No. Independent circular shifts or planted common mode: 0-1 edges among 160-941 candidate pairs (1520: 4-5 of 1568, 2 spurious pairs of units
in that session); real data 15-403 edges.

**2. Can it tell two different axons from fragments of one?** Not by correlation. Simulation (noisy fragments carrying the session's real noise, SNR 0.7-3): two axons whose true signal
correlation is >= 0.5 are joined in ~100 % of runs (only 0.3 separates them). A one-factor model is exactly identified for 2-3 members, so **pairs and triples can never be certified
homogeneous**. The rank-1 test (>= 4 members) flags mixed units at true correlation 0.95: 10-47 %, 0.9: 23-97 %, <= 0.8: 87-100 % (per session), at a 5 % false-flag rate. Only 1 of 33 multi-ROI units
(1423) passed it, so only one unit is tier A; 16 are `heterogeneous`, the rest `untestable`.

**3. Does grouping manufacture homogeneity?** Not measurably. PC1 share of the unit means is at or below the ROI-level share in the same half (ratios 0.73-1.08, in-sample and held-out alike), and the mean pair
correlation among units is at or below that among ROIs (within 0.02; e.g. 1630 0.48-0.53 vs 0.61) and below that of random merges of the same sizes (1630 0.69-0.72). Units stay coherent in the held-out half
(residual r among members median 0.27-0.64, above the near-ROI control's 95th percentile in 50-72 % of units).

**4. Does it improve the signal?** Little beyond averaging. Held-out SNR gain median 1.3-1.4 (ideal sqrt(n) 1.4-2.2) and above the control "mean of n near ungrouped ROIs" in only 1-5 units per
session (p95) because the shared state signal (PC1 share 0.3-0.8) is in every accepted ROI. The gain that grouping adds is in the residual (non-shared) signal, which is what held-out residual r shows.

**5. Do members of a unit differ?** Yes, in behaviour coupling, reproducibly across session halves (Spearman rho of member-minus-unit coupling half 1 vs half 2: face -0.5 (1642) to 0.8, whisking 0.3-0.9,
pupil 0.3-0.9, ball 0.4-0.7; SD of the member deviation 0.05-0.09 (1630, 1642), 0.07-0.16 (1423, 1520), 0.2-0.3 (1138)). Part of it follows member SNR (rho 0.1-0.55 in 1630 / 1642 / 1423), in 1138 and 1520 it does not (rho <= 0).
The unit mean hides this difference; per-ROI results stay the primary level.

**6. Attachments replicate only at high confidence** (1423 <-> 1520 through the fov_matching master axons; pairs of ROIs both accepted in both sessions within 30 µm): pairs sharing a unit in A share
one in B in 95 % (1423 -> 1520) and 60 % (reverse) vs 16 % / 1 % for pairs in different units. Rejected ROIs attached `high` to a unit end up with the same unit in B in 8 of 8 and 8 of 16 cases;
`medium` in 3 of 14 and 3 of 10 (about a quarter, so `medium` is weak evidence); n is small and the sessions are same-day, so this shows stable structure, not axon identity.

**7. Units vs the matched noise blobs** (mean of n blobs for a unit of n): unit |r| with behaviour above the 95th percentile of noise draws in 6 of 7 (1630), 0 of 3 (1138), 2 of 7 (1642), 4 of 9 (1423),
5 of 7 (1520) units; SNR and event rate above the noise almost everywhere. **The noise reference is only brightness-matched in 1630** (standardised difference <= 0.11): in the other sessions the strict pool
(rejected, inactive, not attached) holds only 19-115 blobs and they are dimmer / have a higher neuropil share than the axons (differences 0.4-1.7 SD), and in 1423 the pool is 21 blobs, so top-20 is essentially all of them.
Treat units-vs-noise differences in those sessions as confounded with brightness.

Parameter sweep (`validation_pc1.pkl`, `sweep`): without PC removal everything within 10 µm merges into one unit (1630: 24 of 27 ROIs); with PC1 removed the count of multi-member units is stable
(1630: 5-7 for d_frag 3-10 µm at r_min 0.3, 3-4 at 0.5); removing 2 PCs or behaviour changes it by 0-3 units.

**What this means for use.** Per-ROI analyses stay the default. Unit means are justified only for tier A (here 1 unit). Pair / triple units are co-activity annotations (`cluster_id`, `unit_id`) that
tell you which ROIs share non-state signal; averaging them is an assumption whose cost is quantified by the simulation (axons that share >= 50 % of their signal look the same to this test).
High-confidence attachments and orphans are candidates for the existing curation workflow (review_queue.py / apply_curation.py), not accepted ROIs.

## Caveats
* A unit is a co-active set, not proven to be one neuron: branches of one neuron correlate when apart, and two neurons driven by the same events can look alike.
* Removing PC1 (and behaviour) from the grouping makes groups reflect non-shared structure; the group count depends on the removal level (the sweep reports it).
* Selection is blind to behaviour coupling; the noise reference is selected on brightness, transient test and skew only. Accepted ROIs were curated by eye, so a units-vs-noise difference is a property of
  curated ROIs vs uncurated blobs, not a proof about the axons.
* Simulated fragments use independent real noise per fragment; real fragments share some noise (neuropil, motion), which would make the rank-1 test somewhat conservative, not liberal.
* The `shift` mask (rigid registration excursion) removes up to 18 % of the frames in the high-zoom November sessions; olive-green orphan outlines along vessel walls in October sessions are probably
  vessel pulsation, not axons: review before using them.
* Three animals: pool per session first (session = unit, animal = cluster), as in the paper figures.
