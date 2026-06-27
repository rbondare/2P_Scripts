function [r_same, r_shifted] = verify_direction_reconstruction(Stimuli, dff, TimeCa, stype, orientation_field)
    % Step-by-step diagnostic for the RNG-replay direction reconstruction
    % (reconstruct_direction_order.m) used by grating/moving_bar. Re-confirms
    % correctness on every new dataset rather than trusting a one-time check.
    %
    % Logic: if the reconstructed per-subtrial direction labels are correct,
    % each cell's response-amplitude-by-direction profile should repeat
    % consistently across the stimulus's outer trial repeats (since the
    % SAME RNG order is replayed for every repeat). So: split subtrials into
    % repeat cycles by position, build each cell's [n_repeats x n_dir]
    % amplitude profile, and correlate cycle-to-cycle profiles at MATCHED
    % positions ("same-position") vs a circularly-shifted control. A real,
    % correct reconstruction should show same-position correlation clearly
    % above the shifted control (validated previously at ~0.20 vs ~0.00).

    r_same = NaN; r_shifted = NaN;
    pre_window_sec = 1; post_window_sec = 6; n_com = 100;
    peak_window = [0.5, min(4, post_window_sec)];

    [Ca, t_com, direction_of] = extract_direction_trials(Stimuli, dff, TimeCa, ...
        stype, orientation_field, pre_window_sec, post_window_sec, n_com);
    if isempty(Ca)
        fprintf('  [verify %s] stimulus not found -- skipping\n', stype);
        return;
    end

    unique_dirs = unique(direction_of);
    n_dir = numel(unique_dirs);
    n_subtrials = numel(direction_of);
    n_repeats = floor(n_subtrials / n_dir);
    if n_repeats < 2
        fprintf('  [verify %s] fewer than 2 full repeat cycles (%d subtrials, %d directions) -- cannot cross-validate\n', ...
            stype, n_subtrials, n_dir);
        return;
    end

    peak_idx = t_com >= peak_window(1) & t_com <= peak_window(2);
    amp_per_subtrial = squeeze(mean(Ca(:, peak_idx, :), 2, 'omitnan'));  % n_cells x n_subtrials
    pos_of = mod((0:n_subtrials-1), n_dir) + 1;

    % [n_cells x n_dir x n_repeats] amplitude profile, position-indexed
    n_cells = size(amp_per_subtrial, 1);
    profile = nan(n_cells, n_dir, n_repeats);
    for r = 1:n_repeats
        for p = 1:n_dir
            cols = find(pos_of == p, n_repeats);
            if r <= numel(cols)
                profile(:, p, r) = amp_per_subtrial(:, cols(r));
            end
        end
    end

    same_corrs = [];
    shifted_corrs = [];
    for r1 = 1:n_repeats
        for r2 = (r1+1):n_repeats
            for c = 1:n_cells
                v1 = profile(c, :, r1); v2 = profile(c, :, r2);
                if any(isnan(v1)) || any(isnan(v2)) || std(v1) == 0 || std(v2) == 0; continue; end
                same_corrs(end+1) = corr(v1', v2');              %#ok<AGROW>
                v2_shifted = circshift(v2, 1);
                shifted_corrs(end+1) = corr(v1', v2_shifted');   %#ok<AGROW>
            end
        end
    end

    r_same = mean(same_corrs, 'omitnan');
    r_shifted = mean(shifted_corrs, 'omitnan');

    fprintf('  [verify %s] same-position corr = %.3f, shifted-position control = %.3f (n_dir=%d, n_repeats=%d, n_cells=%d) -- %s\n', ...
        stype, r_same, r_shifted, n_dir, n_repeats, n_cells, ...
        verdict(r_same, r_shifted));
end

function v = verdict(r_same, r_shifted)
    if r_same > r_shifted + 0.05
        v = 'PASS (reconstruction looks correct)';
    else
        v = 'WARNING -- same-position correlation not clearly above shifted control, inspect manually';
    end
end
