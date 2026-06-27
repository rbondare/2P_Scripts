function Z = compute_baseline_zscore_trace(Ca_trials, match_idx, t_com)
    % Same per-trial baseline normalization as compute_baseline_zscore,
    % but keeps the full time course instead of collapsing to one window
    % -- used for the z-scored heatmap (heatmap #2).
    pre_idx  = t_com < 0;
    Ca_sel   = Ca_trials(match_idx, :, :);
    n_trials = size(Ca_sel, 3);

    Z_trial = nan(size(Ca_sel));
    for tr = 1:n_trials
        base_mean = mean(Ca_sel(:, pre_idx, tr), 2, 'omitnan');
        base_std  = std(Ca_sel(:, pre_idx, tr), 0, 2, 'omitnan');
        Z_trial(:, :, tr) = (Ca_sel(:, :, tr) - base_mean) ./ base_std;
    end
    Z = median(Z_trial, 3, 'omitnan');
end
