function [Ca_by_dir, unique_dirs] = chunk_by_direction(Ca_subtrials, direction_of)
    % Averages (median) repeats of the SAME realized direction together --
    % the calcium data needs to correspond to the same direction in both
    % conditions: baseline's direction X gets compared against drug's
    % direction X, matched by actual angle, not subtrial position.
    unique_dirs = unique(direction_of);
    n_rois = size(Ca_subtrials, 1);
    n_com  = size(Ca_subtrials, 2);
    Ca_by_dir = nan(n_rois, n_com, numel(unique_dirs));
    for d = 1:numel(unique_dirs)
        cols = direction_of == unique_dirs(d);
        Ca_by_dir(:, :, d) = median(Ca_subtrials(:, :, cols), 3, 'omitnan');
    end
end
