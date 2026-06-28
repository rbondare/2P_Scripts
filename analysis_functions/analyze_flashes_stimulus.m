function result = analyze_flashes_stimulus(base_Stimuli, base_dff, base_TimeCa, ...
        drug_Stimuli, drug_dff, drug_TimeCa, base_match_idx_local, drug_match_idx_local, matched_rois_available)
    % Flashes (sparse_local_global_flashes), baseline vs drug. Pure
    % data/metric function -- see analyze_looming_stimulus.m for the
    % shared return-struct contract.

    result = struct('t_com', [], 'resp_base', [], 'resp_drug', [], ...
        'Ztrace_base', [], 'Ztrace_drug', [], ...
        'raw_metric_base', [], 'raw_metric_drug', [], ...
        'z_metric_base', [], 'z_metric_drug', [], ...
        'raw_all_base', [], 'raw_all_drug', [], ...
        'z_all_base', [], 'z_all_drug', [], ...
        'has_transient', true, 'n_matched', 0, 'stim_label', 'Flashes', ...
        'common_dirs', [], 'amp_base', [], 'amp_drug', [], ...
        'dir_mean_base', [], 'dir_sem_base', [], 'dir_mean_drug', [], 'dir_sem_drug', [], ...
        't_com_block', [], 'resp_base_block', [], 'resp_drug_block', [], ...
        'Ztrace_base_block', [], 'Ztrace_drug_block', [], ...
        'resp_base_block_all', [], 'resp_drug_block_all', [], ...
        'resp_base_trialblocks', [], 'resp_drug_trialblocks', [], ...
        'trialblock_boundaries', [], 'trialblock_dir_labels', []);

    fprintf('    - Analyzing flashes stimulus (matched baseline vs drug)\n');
    if ~matched_rois_available
        fprintf('      (No matched ROIs available -- skipping)\n');
        return;
    end

    pre_window_sec = 1; post_window_sec = 4; n_com = 100;
    peak_window = [0 0.8];

    base_onsets = get_flash_onsets(base_Stimuli);
    drug_onsets = get_flash_onsets(drug_Stimuli);
    if isempty(base_onsets) || isempty(drug_onsets)
        fprintf('      (Flashes missing in baseline or drug -- skipping)\n');
        return;
    end

    Ca_base = align_trials_to_onsets(base_onsets, base_dff, base_TimeCa, pre_window_sec, post_window_sec, n_com);
    Ca_drug = align_trials_to_onsets(drug_onsets, drug_dff, drug_TimeCa, pre_window_sec, post_window_sec, n_com);
    t_com = linspace(-pre_window_sec, post_window_sec, n_com);
    result.t_com = t_com;
    peak_idx = t_com >= peak_window(1) & t_com <= peak_window(2);

    resp_base = median(Ca_base(base_match_idx_local, :, :), 3, 'omitnan');
    resp_drug = median(Ca_drug(drug_match_idx_local, :, :), 3, 'omitnan');
    n_matched = size(resp_base, 1);
    result.resp_base = resp_base;
    result.resp_drug = resp_drug;
    result.n_matched = n_matched;

    result.raw_metric_base = mean(resp_base(:, peak_idx), 2, 'omitnan');
    result.raw_metric_drug = mean(resp_drug(:, peak_idx), 2, 'omitnan');

    result.z_metric_base = compute_baseline_zscore(Ca_base, base_match_idx_local, t_com, peak_window);
    result.z_metric_drug = compute_baseline_zscore(Ca_drug, drug_match_idx_local, t_com, peak_window);

    result.Ztrace_base = compute_baseline_zscore_trace(Ca_base, base_match_idx_local, t_com);
    result.Ztrace_drug = compute_baseline_zscore_trace(Ca_drug, drug_match_idx_local, t_com);

    all_base_idx = 1:size(Ca_base, 1);
    all_drug_idx = 1:size(Ca_drug, 1);
    resp_all_base = median(Ca_base(all_base_idx, :, :), 3, 'omitnan');
    resp_all_drug = median(Ca_drug(all_drug_idx, :, :), 3, 'omitnan');
    result.raw_all_base = mean(resp_all_base(:, peak_idx), 2, 'omitnan');
    result.raw_all_drug = mean(resp_all_drug(:, peak_idx), 2, 'omitnan');
    result.z_all_base = compute_baseline_zscore(Ca_base, all_base_idx, t_com, peak_window);
    result.z_all_drug = compute_baseline_zscore(Ca_drug, all_drug_idx, t_com, peak_window);
end
