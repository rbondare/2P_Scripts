function z = compute_baseline_zscore_from_resp(resp, t_com, peak_window)
    % Same idea as compute_baseline_zscore, but for a trace that's already
    % been repeat-averaged upstream (e.g. chunk_by_direction's per-angle
    % median) -- z-scores the single merged trace against its own
    % pre-onset segment directly.
    pre_idx  = t_com < 0;
    peak_idx = t_com >= peak_window(1) & t_com <= peak_window(2);
    base_mean = mean(resp(:, pre_idx), 2, 'omitnan');
    base_std  = std(resp(:, pre_idx), 0, 2, 'omitnan');
    peak_mean = mean(resp(:, peak_idx), 2, 'omitnan');
    z = (peak_mean - base_mean) ./ base_std;
end
