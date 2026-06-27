function [z_metric, Z_resampled, resp_resampled] = whole_block_zscore_resampled(block_dff, match_idx, query_axis, src_duration)
    % Whole-block self/population z-score (compute_whole_block_zscore),
    % then resample both the z-scored trace and the raw dF/F trace from
    % the block's own native frame grid onto a shared query_axis so
    % baseline and drug (which generally have slightly different frame
    % counts) end up with matching matrix widths for the heatmap helpers.
    %
    % src_duration sets the units of the block's native axis (and must
    % match query_axis's units): pass 1 for a normalized 0-1 "fraction of
    % block" axis (spontaneous/checkers2, which have no meaningful shared
    % real-time axis across conditions), or the block's actual duration in
    % seconds for a real elapsed-time axis (grating/moving_bar's
    % continuous whole-block heatmap, where 0-N seconds is meaningful).
    [z_metric, Z_trace] = compute_whole_block_zscore(block_dff, match_idx);
    n_frames = size(Z_trace, 2);
    src_axis = linspace(0, src_duration, n_frames);
    Z_resampled = interp1(src_axis, Z_trace', query_axis)';
    resp_resampled = interp1(src_axis, block_dff(match_idx, :)', query_axis)';
end
