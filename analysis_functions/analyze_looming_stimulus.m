function result = analyze_looming_stimulus(base_Stimuli, base_dff, base_TimeCa, base_Triggers, ...
        drug_Stimuli, drug_dff, drug_TimeCa, drug_Triggers, ...
        base_match_idx_local, drug_match_idx_local, matched_rois_available)
    % Looming-triggered response, baseline vs drug. Looming fires exactly
    % one red sync frame per trial, at the true onset of that trial's
    % loom (unlike the other stimulus types, which have a continuous
    % periodic red-sync rhythm). Trials are NOT evenly spaced -- each has
    % a randomized wait period beforehand -- so onsets must come directly
    % from these red frames, not from stimulus_trial_t x trial_index.
    %
    % Pure data/metric function -- no plotting. Returns a struct with the
    % common fields Plot_Stimulus_Results.m expects from every
    % analyze_*_stimulus function (see analysis_functions/README context
    % in Plot_Stimulus_Results.m): t_com, resp_base/drug, Ztrace_base/drug,
    % raw_metric_base/drug, z_metric_base/drug (matched cells), and
    % raw_all_base/drug, z_all_base/drug (all cells, unpaired).

    result = struct('t_com', [], 'resp_base', [], 'resp_drug', [], ...
        'Ztrace_base', [], 'Ztrace_drug', [], ...
        'raw_metric_base', [], 'raw_metric_drug', [], ...
        'z_metric_base', [], 'z_metric_drug', [], ...
        'raw_all_base', [], 'raw_all_drug', [], ...
        'z_all_base', [], 'z_all_drug', [], ...
        'has_transient', true, 'n_matched', 0, 'stim_label', 'Looming', ...
        'common_dirs', [], 'amp_base', [], 'amp_drug', [], ...
        'dir_mean_base', [], 'dir_sem_base', [], 'dir_mean_drug', [], 'dir_sem_drug', [], ...
        't_com_block', [], 'resp_base_block', [], 'resp_drug_block', [], ...
        'Ztrace_base_block', [], 'Ztrace_drug_block', []);

    fprintf('    - Analyzing looming stimulus (matched baseline vs drug)\n');

    if ~matched_rois_available
        fprintf('      (No matched ROIs available -- skipping looming baseline/drug comparison)\n');
        return;
    end

    [Ca_base, t_com_base] = extract_looming_trials(base_Stimuli, base_dff, base_TimeCa, base_Triggers);
    [Ca_drug, t_com_drug] = extract_looming_trials(drug_Stimuli, drug_dff, drug_TimeCa, drug_Triggers);

    if isempty(Ca_base) || isempty(Ca_drug)
        fprintf('      (Looming stimulus missing in baseline or drug -- skipping)\n');
        return;
    end
    if ~isequal(size(t_com_base), size(t_com_drug)) || max(abs(t_com_base - t_com_drug)) > 1e-6
        warning('Baseline and drug looming analysis windows differ (stimulus_trial_t mismatch?) -- using baseline time axis');
    end
    t_com = t_com_base;
    result.t_com = t_com;

    loom_peak_window = [0.5 1.5];  % matches looming's sharp peak
    peak_idx = t_com >= loom_peak_window(1) & t_com <= loom_peak_window(2);

    % Restrict to matched cells. Use the MEDIAN across each condition's own
    % trials, not the mean -- with only a handful of trials per cell, a
    % single large-but-unreliable trial can otherwise dominate the mean
    % and flip a cell's responsive/non-responsive classification. The
    % median requires a majority of trials to be elevated instead.
    resp_base = median(Ca_base(base_match_idx_local, :, :), 3, 'omitnan');   % [n_matched x n_com]
    resp_drug = median(Ca_drug(drug_match_idx_local, :, :), 3, 'omitnan');   % [n_matched x n_com]
    n_matched = size(resp_base, 1);
    result.resp_base = resp_base;
    result.resp_drug = resp_drug;
    result.n_matched = n_matched;

    result.raw_metric_base = mean(resp_base(:, peak_idx), 2, 'omitnan');
    result.raw_metric_drug = mean(resp_drug(:, peak_idx), 2, 'omitnan');

    result.z_metric_base = compute_baseline_zscore(Ca_base, base_match_idx_local, t_com, loom_peak_window);
    result.z_metric_drug = compute_baseline_zscore(Ca_drug, drug_match_idx_local, t_com, loom_peak_window);

    result.Ztrace_base = compute_baseline_zscore_trace(Ca_base, base_match_idx_local, t_com);
    result.Ztrace_drug = compute_baseline_zscore_trace(Ca_drug, drug_match_idx_local, t_com);

    % All cells (unpaired, each session independently)
    all_base_idx = 1:size(Ca_base, 1);
    all_drug_idx = 1:size(Ca_drug, 1);
    resp_all_base = median(Ca_base(all_base_idx, :, :), 3, 'omitnan');
    resp_all_drug = median(Ca_drug(all_drug_idx, :, :), 3, 'omitnan');
    result.raw_all_base = mean(resp_all_base(:, peak_idx), 2, 'omitnan');
    result.raw_all_drug = mean(resp_all_drug(:, peak_idx), 2, 'omitnan');
    result.z_all_base = compute_baseline_zscore(Ca_base, all_base_idx, t_com, loom_peak_window);
    result.z_all_drug = compute_baseline_zscore(Ca_drug, all_drug_idx, t_com, loom_peak_window);

    % Responsive-cell diagnostic (console only, mirrors the prior
    % Master_Stimulus_Analysis.m console output; the bar/scatter figures
    % that used to accompany this were looming-specific one-offs that
    % don't fit the generic 6-stimulus plot suite, so they were dropped
    % rather than ported -- see report to user).
    post1s = t_com >= 1;
    pct_base = 100 * mean(mean(resp_base(:, post1s), 2) > 0.2);
    pct_drug = 100 * mean(mean(resp_drug(:, post1s), 2) > 0.2);
    fprintf('      Responsive cells (mean dF/F > 0.2, %.1fs to %.1fs post-onset): baseline %.1f%%, drug %.1f%%\n', ...
        t_com(find(post1s, 1)), t_com(end), pct_base, pct_drug);
end
