function result = analyze_blockwise_stimulus(stype, base_Stimuli, base_dff, base_TimeCa, ...
        drug_Stimuli, drug_dff, drug_TimeCa, base_match_idx_local, drug_match_idx_local, matched_rois_available)
    % Shared implementation for spontaneous and checkers2: no discrete
    % onset/trial structure (a single continuous block), so the
    % onset-relative baseline z-score used elsewhere doesn't apply.
    % "z-scored" here means each cell's own whole-block self z-score (see
    % compute_whole_block_zscore.m). No transient plot (nothing to anchor
    % one to). Baseline and drug blocks are independently resampled onto a
    % common normalized "fraction of block" axis (0-1) purely so the
    % returned matrices share a width for the heatmap/violin/histogram
    % helpers -- there's no shared onset to align sub-block timing to.

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

    fprintf('    - Analyzing %s (block-wise, matched baseline vs drug)\n', stype);
    if ~matched_rois_available
        fprintf('      (No matched ROIs available -- skipping)\n');
        return;
    end

    resp_base_cell = extract_full_stimulus_responses(base_Stimuli, base_dff, stype, base_TimeCa);
    resp_drug_cell = extract_full_stimulus_responses(drug_Stimuli, drug_dff, stype, drug_TimeCa);
    if isempty(resp_base_cell) || isempty(resp_drug_cell)
        fprintf('      (%s missing in baseline or drug -- skipping)\n', stype);
        return;
    end
    block_base = resp_base_cell{1};  % [n_rois x n_frames_in_block], first presentation
    block_drug = resp_drug_cell{1};

    n_com = 200;
    frac_axis = linspace(0, 1, n_com);
    result.t_com = frac_axis;

    [z_metric_base, Ztrace_base, resp_base_rs] = whole_block_zscore_resampled(block_base, base_match_idx_local, frac_axis, 1);
    [z_metric_drug, Ztrace_drug, resp_drug_rs] = whole_block_zscore_resampled(block_drug, drug_match_idx_local, frac_axis, 1);
    n_matched = numel(z_metric_base);
    result.n_matched = n_matched;
    result.resp_base = resp_base_rs;
    result.resp_drug = resp_drug_rs;
    result.Ztrace_base = Ztrace_base;
    result.Ztrace_drug = Ztrace_drug;
    result.z_metric_base = z_metric_base;
    result.z_metric_drug = z_metric_drug;
    result.raw_metric_base = mean(block_base(base_match_idx_local, :), 2, 'omitnan');
    result.raw_metric_drug = mean(block_drug(drug_match_idx_local, :), 2, 'omitnan');

    all_base_idx = 1:size(block_base, 1);
    all_drug_idx = 1:size(block_drug, 1);
    [z_all_base, ~, ~] = whole_block_zscore_resampled(block_base, all_base_idx, frac_axis, 1);
    [z_all_drug, ~, ~] = whole_block_zscore_resampled(block_drug, all_drug_idx, frac_axis, 1);
    result.z_all_base = z_all_base;
    result.z_all_drug = z_all_drug;
    result.raw_all_base = mean(block_base(all_base_idx, :), 2, 'omitnan');
    result.raw_all_drug = mean(block_drug(all_drug_idx, :), 2, 'omitnan');
end
