function [Ca_concat, boundaries, dir_labels] = reorder_subtrials_by_direction(Ca_subtrials, direction_of, dir_order)
    % Reorders + concatenates subtrials so all repeats of the same
    % direction sit adjacent to each other, in dir_order's sequence,
    % instead of the natural chronological order (which interleaves every
    % direction once per repeat-cycle). Makes direction-tuning structure
    % visible as distinct vertical blocks in a heatmap, at the cost of no
    % longer reading left-to-right as real elapsed time.
    %
    % Ca_subtrials: [n_cells x n_com x n_subtrials] (one subtrial's
    % aligned snippet per 3rd-dim slice, e.g. from extract_direction_trials).
    % direction_of: [1 x n_subtrials], the realized direction per subtrial.
    % dir_order: the direction values to group by, in the order they
    % should appear left-to-right (pass baseline/drug's shared common_dirs
    % so both conditions end up with identical column structure).
    %
    % Returns:
    %   Ca_concat: [n_cells x total_columns], total_columns = sum over
    %     directions of (repeats_of_that_direction * n_com).
    %   boundaries: [1 x numel(dir_order)], the column index where each
    %     direction's block STARTS (1-indexed) -- for separator lines.
    %   dir_labels: [1 x numel(dir_order)], same order as boundaries.
    n_com = size(Ca_subtrials, 2);
    n_dir = numel(dir_order);

    Ca_concat = [];
    boundaries = zeros(1, n_dir);
    dir_labels = dir_order(:)';
    col = 0;
    for d = 1:n_dir
        cols = find(direction_of == dir_order(d));
        boundaries(d) = col + 1;
        for c = 1:numel(cols)
            Ca_concat = cat(2, Ca_concat, Ca_subtrials(:, :, cols(c)));
            col = col + n_com;
        end
    end
end
