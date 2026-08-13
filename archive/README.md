# archive/

Superseded / one-off scripts moved here on 2026-07-29 to declutter the top
level while building the clean `session_pipeline/`. **Nothing here is used by
the active pipeline** — verified by dependency scan (only comment references
remained). Moving is fully reversible: `mv archive/<file> ..` to restore.

Kept OUT of here on purpose (still active):
- `reconstructFlashes.m` — real runtime dependency of `analysis_functions/get_flash_onsets.m`
- `Preprocess_Stimulus_Data.m`, `Plot_Stimulus_Results.m` — the current 3-script plotting pipeline
- `analysis_functions/` — extractors reused by `session_pipeline/`

Contents:
- **Legacy monolith / stubs:** `Master_Stimulus_Analysis.m`, `RB_DataGathering.m`,
  `organize_stimulus_by_parameters.m`, `preprocess_and_organize_stimulus.m`,
  `stimulus_response_comparison_pipeline.m`, `Stimulus_Response_*.m`
- **One-offs / posters / explorers:** `Poster_Looming_*.m`, `Looming_Individual_Cells_0312.m`,
  `Temp_Looming_0319_UnmatchedCheck.m`, `Explore_Heatmaps.*`, `Comprehensive_Data_Explorer.m`,
  `Plottingdff_from_preprocessed.m`, `MovingBar_AllDatasets_Comparison.m`, `Loom_Flashes_ZETA_Analysis.m`
- **Diagnostics / tests:** `Diag_*.m`, `Diag_outlog.txt`,
  `Test_Loom_Flashes_TriggeredAverage.m`, `Test_MovingBar_Chunking_Walkthrough.m`
- **Legacy aggregation:** `AggregationCode/` (predecessor of `newAggregation/`)
- `__pycache__/`
