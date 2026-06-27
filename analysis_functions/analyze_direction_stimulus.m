function result = analyze_direction_stimulus(stype, orientation_field, ...
        base_Stimuli, base_dff, base_TimeCa, drug_Stimuli, drug_dff, drug_TimeCa, ...
        base_match_idx_local, drug_match_idx_local, matched_rois_available)
    % Shared implementation for moving_bar and grating: reconstruct
    % realized direction per subtrial, average repeats of the same
    % direction, compare baseline vs drug at matched angles, then collapse
    % to each cell's own preferred direction (from baseline, for the
    % matched population; from its own session, for the "all cells"
    % population -- there's no cell-to-cell correspondence to carry a
    % direction across for the unpaired set anyway).
    %
    % No transient/triggered-average plot for this stimulus type (the
    % direction-resolved alignment is still being inspected) -- has_transient
    % is false. See analyze_looming_stimulus.m for the shared return-struct
    % contract.

    result = struct('t_com', [], 'resp_base', [], 'resp_drug', [], ...
        'Ztrace_base', [], 'Ztrace_drug', [], ...
        'raw_metric_base', [], 'raw_metric_drug', [], ...
        'z_metric_base', [], 'z_metric_drug', [], ...
        'raw_all_base', [], 'raw_all_drug', [], ...
        'z_all_base', [], 'z_all_drug', [], ...
        'has_transient', false, 'n_matched', 0, 'stim_label', stype, ...
        'common_dirs', [], 'amp_base', [], 'amp_drug', [], ...
        'dir_mean_base', [], 'dir_sem_base', [], 'dir_mean_drug', [], 'dir_sem_drug', [], ...
        't_com_block', [], 'resp_base_block', [], 'resp_drug_block', [], ...
        'Ztrace_base_block', [], 'Ztrace_drug_block', []);

    fprintf('    - Analyzing %s stimulus (direction-resolved, matched baseline vs drug)\n', stype);
    if ~matched_rois_available
        fprintf('      (No matched ROIs available -- skipping)\n');
        return;
    end

    pre_window_sec = 1; post_window_sec = 6; n_com = 100;
    peak_window = [0.5, min(4, post_window_sec)];

    [Ca_base, t_com, dir_base] = extract_direction_trials(base_Stimuli, base_dff, base_TimeCa, ...
        stype, orientation_field, pre_window_sec, post_window_sec, n_com);
    [Ca_drug, ~, dir_drug] = extract_direction_trials(drug_Stimuli, drug_dff, drug_TimeCa, ...
        stype, orientation_field, pre_window_sec, post_window_sec, n_com);

    if isempty(Ca_base) || isempty(Ca_drug)
        fprintf('      (%s missing in baseline or drug -- skipping)\n', stype);
        return;
    end
    result.t_com = t_com;
    peak_idx = t_com >= peak_window(1) & t_com <= peak_window(2);

    % ---- Matched cells: preferred direction from baseline, carried to drug ----
    [Ca_by_dir_base, dirs_base] = chunk_by_direction(Ca_base(base_match_idx_local, :, :), dir_base);
    [Ca_by_dir_drug, dirs_drug] = chunk_by_direction(Ca_drug(drug_match_idx_local, :, :), dir_drug);

    common_dirs = intersect(dirs_base, dirs_drug);
    if numel(common_dirs) < numel(union(dirs_base, dirs_drug))
        warning('%s: baseline/drug direction sets differ -- using %d common directions', stype, numel(common_dirs));
    end
    [~, ib]  = ismember(common_dirs, dirs_base);
    [~, idr] = ismember(common_dirs, dirs_drug);
    Ca_by_dir_base = Ca_by_dir_base(:, :, ib);
    Ca_by_dir_drug = Ca_by_dir_drug(:, :, idr);
    n_matched = size(Ca_by_dir_base, 1);
    result.n_matched = n_matched;

    % Per-orientation tuning curve + time-resolved population response
    % (population mean +/- SEM across matched cells, one column per
    % common direction) -- "how the response to each orientation is",
    % independent of any single cell's preferred direction.
    result.common_dirs = common_dirs;
    amp_per_cell_base = squeeze(mean(Ca_by_dir_base(:, peak_idx, :), 2, 'omitnan'));  % n_cells x n_dirs
    amp_per_cell_drug = squeeze(mean(Ca_by_dir_drug(:, peak_idx, :), 2, 'omitnan'));
    result.amp_base = mean(amp_per_cell_base, 1, 'omitnan');
    result.amp_drug = mean(amp_per_cell_drug, 1, 'omitnan');
    result.dir_mean_base = squeeze(mean(Ca_by_dir_base, 1, 'omitnan'))';   % n_dirs x n_com
    result.dir_sem_base  = squeeze(std(Ca_by_dir_base, 0, 1, 'omitnan'))' / sqrt(n_matched);
    result.dir_mean_drug = squeeze(mean(Ca_by_dir_drug, 1, 'omitnan'))';
    result.dir_sem_drug  = squeeze(std(Ca_by_dir_drug, 0, 1, 'omitnan'))' / sqrt(n_matched);

    [~, pref_dir_idx] = max(amp_per_cell_base, [], 2);

    resp_base_pref = nan(n_matched, numel(t_com));
    resp_drug_pref = nan(n_matched, numel(t_com));
    for c = 1:n_matched
        resp_base_pref(c, :) = Ca_by_dir_base(c, :, pref_dir_idx(c));
        resp_drug_pref(c, :) = Ca_by_dir_drug(c, :, pref_dir_idx(c));
    end
    result.resp_base = resp_base_pref;
    result.resp_drug = resp_drug_pref;

    result.raw_metric_base = mean(resp_base_pref(:, peak_idx), 2, 'omitnan');
    result.raw_metric_drug = mean(resp_drug_pref(:, peak_idx), 2, 'omitnan');
    result.z_metric_base = compute_baseline_zscore_from_resp(resp_base_pref, t_com, peak_window);
    result.z_metric_drug = compute_baseline_zscore_from_resp(resp_drug_pref, t_com, peak_window);

    % Full time-course z-score (per matched cell, at its own preferred
    % direction) for the z-scored heatmap.
    pre_idx = t_com < 0;
    bm = mean(resp_base_pref(:, pre_idx), 2, 'omitnan'); bs = std(resp_base_pref(:, pre_idx), 0, 2, 'omitnan');
    dm = mean(resp_drug_pref(:, pre_idx), 2, 'omitnan'); ds = std(resp_drug_pref(:, pre_idx), 0, 2, 'omitnan');
    result.Ztrace_base = (resp_base_pref - bm) ./ bs;
    result.Ztrace_drug = (resp_drug_pref - dm) ./ ds;

    % ---- All cells: each session's own preferred direction (independent, unpaired) ----
    [resp_all_base, ~] = preferred_direction_resp(Ca_base, dir_base, peak_idx, t_com);
    [resp_all_drug, ~] = preferred_direction_resp(Ca_drug, dir_drug, peak_idx, t_com);
    result.raw_all_base = mean(resp_all_base(:, peak_idx), 2, 'omitnan');
    result.raw_all_drug = mean(resp_all_drug(:, peak_idx), 2, 'omitnan');
    result.z_all_base = compute_baseline_zscore_from_resp(resp_all_base, t_com, peak_window);
    result.z_all_drug = compute_baseline_zscore_from_resp(resp_all_drug, t_com, peak_window);

    % ---- Continuous whole-block heatmap (matched cells) ----
    % resp_base/resp_drug above are each cell's own ~7s preferred-direction
    % snippet, not the continuous stimulus -- useful for the violin/MI
    % metrics, but it collapses away the natural repeat structure across
    % all directions/repeats that a population heatmap is meant to show.
    % This extracts the SAME entire block extract_full_stimulus_responses
    % gives spontaneous/checkers2's heatmap, resampled onto a shared
    % real-seconds axis (using baseline's own block duration), so the
    % heatmap shows the full ~tens-of-seconds block rather than one
    % direction's short window.
    resp_base_cell_block = extract_full_stimulus_responses(base_Stimuli, base_dff, stype, base_TimeCa);
    resp_drug_cell_block = extract_full_stimulus_responses(drug_Stimuli, drug_dff, stype, drug_TimeCa);
    if ~isempty(resp_base_cell_block) && ~isempty(resp_drug_cell_block)
        block_base = resp_base_cell_block{1};
        block_drug = resp_drug_cell_block{1};
        n_com_block = 300;
        dt_base = median(diff(base_TimeCa(1, :)), 'omitnan');
        dur_base_sec = (size(block_base, 2) - 1) * dt_base;
        t_com_block = linspace(0, dur_base_sec, n_com_block);
        result.t_com_block = t_com_block;

        [~, Ztrace_base_block, resp_base_block] = whole_block_zscore_resampled(...
            block_base, base_match_idx_local, t_com_block, dur_base_sec);

        dt_drug = median(diff(drug_TimeCa(1, :)), 'omitnan');
        dur_drug_sec = (size(block_drug, 2) - 1) * dt_drug;
        [~, Ztrace_drug_block, resp_drug_block] = whole_block_zscore_resampled(...
            block_drug, drug_match_idx_local, t_com_block, dur_drug_sec);

        result.resp_base_block = resp_base_block;
        result.resp_drug_block = resp_drug_block;
        result.Ztrace_base_block = Ztrace_base_block;
        result.Ztrace_drug_block = Ztrace_drug_block;
    end
end

function [resp_pref, pref_dir_idx] = preferred_direction_resp(Ca, direction_of, peak_idx, t_com)
    [Ca_by_dir, ~] = chunk_by_direction(Ca, direction_of);
    amp_per_cell = squeeze(mean(Ca_by_dir(:, peak_idx, :), 2, 'omitnan'));
    [~, pref_dir_idx] = max(amp_per_cell, [], 2);
    n_cells = size(Ca_by_dir, 1);
    resp_pref = nan(n_cells, numel(t_com));
    for c = 1:n_cells
        resp_pref(c, :) = Ca_by_dir(c, :, pref_dir_idx(c));
    end
end
