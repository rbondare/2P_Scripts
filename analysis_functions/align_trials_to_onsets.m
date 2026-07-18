function Ca_trials = align_trials_to_onsets(onsets, dff, TimeCa, pre_window_sec, post_window_sec, n_com)
    t_com    = linspace(-pre_window_sec, post_window_sec, n_com);
    n_rois   = size(dff, 1);
    n_trials = numel(onsets);
    Ca_trials = nan(n_rois, n_com, n_trials);
    % Accept either a single plane's time vector (used directly) or the
    % full [n_planes x n_volumes] matrix (defaults to row 1, preserving
    % existing callers that haven't been updated to pass a specific plane).
    if isvector(TimeCa)
        time_vec = TimeCa(:)';
    else
        time_vec = TimeCa(1, :);
    end
    for k = 1:n_trials
        idx = find(time_vec > onsets(k) - pre_window_sec & time_vec < onsets(k) + post_window_sec);
        if numel(idx) < 2; continue; end
        t_rel = time_vec(idx) - onsets(k);
        Ca_trials(:, :, k) = interp1(t_rel(:), dff(:, idx)', t_com(:), 'linear', 'extrap')';
    end
end
