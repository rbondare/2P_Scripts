function z = compute_baseline_zscore(Ca_trials, match_idx, t_com, peak_window)
    % Per-cell, per-trial: (peak-window mean - pre-onset baseline mean) /
    % pre-onset baseline std; then median across trials.
    pre_idx  = t_com < 0;
    peak_idx = t_com >= peak_window(1) & t_com <= peak_window(2);
    Ca_sel   = Ca_trials(match_idx, :, :);
    n_trials = size(Ca_sel, 3);
    z_trial = nan(size(Ca_sel, 1), n_trials);
    for tr = 1:n_trials
        base_mean = mean(Ca_sel(:, pre_idx, tr), 2, 'omitnan');
        base_std  = std(Ca_sel(:, pre_idx, tr), 0, 2, 'omitnan');
        peak_mean = mean(Ca_sel(:, peak_idx, tr), 2, 'omitnan');
        z_trial(:, tr) = (peak_mean - base_mean) ./ base_std;
    end
    z = median(z_trial, 2, 'omitnan');
end
