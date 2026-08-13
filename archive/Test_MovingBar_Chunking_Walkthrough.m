%% Step-by-step walkthrough: how moving_bar trial chunking actually works
%
% Baseline recording only -- this is about understanding the MECHANICS of
% analyze_direction_stimulus.m, not a baseline-vs-drug comparison (that's
% what the main pipeline is for). Five sections, each one layer closer to
% the final per-cell-preferred-direction snippet used in the real
% pipeline:
%
%   1. Raw, untouched dF/F for the ENTIRE moving_bar block, ALL cells in
%      the plane (no matching, no trial structure imposed) -- what the
%      calcium signal actually looks like before anything is done to it.
%   2. Same raw traces, with reconstructed subtrial onsets marked --
%      sanity check that onsets line up with visible transients.
%   3. Per-subtrial extraction (Ca_subtrials): the same handful of cells,
%      all 24 individual subtrials, BEFORE any averaging -- trial-to-trial
%      variability is visible here.
%   4. Direction-chunked (chunk_by_direction): the 24 subtrials collapsed
%      to 8 direction-averaged (median) traces per cell.
%   5. Preferred-direction collapse: which ONE of those 8 traces becomes
%      that cell's row in the main pipeline's heatmap/violin/MI.
%
% Run section-by-section (Ctrl+Enter) rather than all at once -- each one
% builds on variables from the previous section.

addpath(fullfile(fileparts(mfilename('fullpath')), 'analysis_functions'));

%% ====================== CONFIGURATION ======================
baseline_file = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260320_1444_preprocessed.mat";
drug_file = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260320_1609_preprocessed.mat";
roi_match_file = "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane0_1444_V1.mat";
selected_plane_idx = 1;
n_example_cells = 6;   % how many cells to show in the small-multiples sections

test_outdir = fullfile(fileparts(mfilename('fullpath')), 'analysis_figures', '_test_movingbar_walkthrough');
if ~exist(test_outdir, 'dir'); mkdir(test_outdir); end

%% ====================== LOAD DATA (baseline + drug, for matching only) ======================
fprintf('Loading baseline: %s\n', baseline_file);
B = load(baseline_file);
base_dff_full = B.CaData(1).Ca_dFF;
centroid = B.CaData(1).Ca_centroid_voxel;
centroidZ = centroid(:, 3);
unique_planes = unique(centroidZ);
selected_roi_idx = find(centroidZ == unique_planes(selected_plane_idx));
base_dff = base_dff_full(selected_roi_idx, :);
fprintf('Selected plane %d: %d ROIs\n', selected_plane_idx, size(base_dff, 1));

% Drug + ROI matching loaded ONLY to get base_match_idx_local (which
% baseline ROIs have a drug-side match) -- the rest of this walkthrough is
% still baseline-only mechanics, matching is just a row-subset for Section 1b.
fprintf('Loading drug (for matching only): %s\n', drug_file);
D = load(drug_file);
drug_centroidZ = D.CaData(1).Ca_centroid_voxel(:, 3);
drug_unique_planes = unique(drug_centroidZ);
drug_selected_roi_idx = find(drug_centroidZ == drug_unique_planes(selected_plane_idx));

M = load(roi_match_file);
base_match_idx_global = M.roiMatchData.allSessionMapping(:, 1);
drug_match_idx_global = M.roiMatchData.allSessionMapping(:, 2);
[~, base_match_idx_local] = ismember(base_match_idx_global, selected_roi_idx);
[~, drug_match_idx_local] = ismember(drug_match_idx_global, drug_selected_roi_idx);
valid_matches = (base_match_idx_local > 0) & (drug_match_idx_local > 0);
base_match_idx_local = base_match_idx_local(valid_matches);
n_matched = numel(base_match_idx_local);
fprintf('Matched ROIs (baseline-local indices): %d\n', n_matched);

%% ====================== SECTION 1: raw, whole block, ALL cells ======================
% No trial structure, no matching -- literally what's in the calcium
% trace during the entire moving_bar stimulus, for every cell in the plane.
resp_cell = extract_full_stimulus_responses(B.Stimuli, base_dff, 'moving_bar', B.TimeCa);
block_base = resp_cell{1};   % [n_rois x n_frames_in_block]
n_rois_total = size(block_base, 1);
fprintf('\nSECTION 1: moving_bar block = %d frames, %d ROIs (all cells, unmatched)\n', size(block_base, 2), n_rois_total);

dt = median(diff(B.TimeCa(1, :)), 'omitnan');
block_t = (0:size(block_base,2)-1) * dt;

figure('Name', 'Section 1 -- raw block, all cells', 'NumberTitle', 'off', 'Position', [80 80 1000 600]);
% Raw dF/F is not sign-balanced -- grayscale (white=low, black=high)
% avoids a diverging map's white sitting at the middle of the range
% instead of at a meaningful zero.
imagesc(block_t, 1:n_rois_total, block_base, [0, prctile(block_base(:), 99)]);
colormap(flipud(gray));
colorbar; xlabel('Time (s)'); ylabel('ROI # (unsorted, all cells in plane)');
title('Section 1: raw dF/F, entire moving\_bar block, ALL cells (no matching, no trial structure)', 'Interpreter', 'none');
save_fig_local('Section1_raw_block_allcells_heatmap', test_outdir);

%% ====================== SECTION 1b: raw, whole block, MATCHED cells only ======================
% Same raw block, same grayscale convention, but restricted to the rows of
% base_dff that have a baseline<->drug match -- the population this
% walkthrough's later sections (and the main pipeline's violin/MI/heatmap)
% actually report on.
block_base_matched = block_base(base_match_idx_local, :);
[~, sort_order_matched] = sort(max(block_base_matched, [], 2), 'descend');

figure('Name', 'Section 1b -- raw block, matched cells only', 'NumberTitle', 'off', 'Position', [80 80 1000 600]);
imagesc(block_t, 1:n_matched, block_base_matched(sort_order_matched, :), [0, prctile(block_base_matched(:), 99)]);
colormap(flipud(gray));
colorbar; xlabel('Time (s)'); ylabel('Matched cell # (sorted by peak activity)');
title(sprintf('Section 1b: raw dF/F, entire moving\\_bar block, MATCHED cells only (n=%d)', n_matched), 'Interpreter', 'none');
save_fig_local('Section1b_raw_block_matchedcells_heatmap', test_outdir);

% Pick example cells by overall raw activity in this block (not matched --
% just the most active cells in the FULL population) so sections 3-5 have
% something visible to show.
[~, order_by_activity] = sort(max(block_base, [], 2), 'descend');
example_cells = order_by_activity(1:n_example_cells)';
fprintf('Example cells (top %d by raw peak activity in block): %s\n', n_example_cells, mat2str(example_cells));

figure('Name', 'Section 1 -- raw traces, example cells', 'NumberTitle', 'off', 'Position', [80 80 1000 700]);
for k = 1:n_example_cells
    subplot(n_example_cells, 1, k);
    plot(block_t, block_base(example_cells(k), :), 'k-', 'LineWidth', 1);
    ylabel(sprintf('ROI %d', example_cells(k)), 'FontSize', 9);
    set(gca, 'Box', 'off', 'FontSize', 8); xlim([block_t(1), block_t(end)]);
    if k < n_example_cells; set(gca, 'XTickLabel', []); end
end
xlabel('Time (s)');
sgtitle('Section 1: raw dF/F traces, example cells, no processing', 'FontWeight', 'bold');
save_fig_local('Section1_raw_traces_examplecells', test_outdir);

%% ====================== SECTION 2: mark reconstructed onsets ======================
% Same raw traces as Section 1, with vertical lines at each reconstructed
% subtrial onset -- do visible transients actually line up with these?
stim_idx = find(strcmp({B.Stimuli.type}, 'moving_bar'), 1);
s = B.Stimuli(stim_idx);
order = reconstruct_direction_order(s, 'bar_orientations');
n_dir = numel(order);
sfc = s.Frameinfo.subtrial_frame_count(~isnan(s.Frameinfo.frame_count));
subtrial_starts = find(sfc == 1);
n_subtrials = numel(subtrial_starts);
pos_of = mod((0:n_subtrials-1), n_dir) + 1;
direction_of = order(pos_of);
onsets_abs = s.TimeStimulusFrame(subtrial_starts);
onsets_rel = onsets_abs - s.TimeStimulusFrame(1);   % relative to block start, same axis as block_t

fprintf('\nSECTION 2: %d subtrials reconstructed, %d unique directions: %s\n', ...
    n_subtrials, n_dir, mat2str(unique(order)));

dir_colors = lines(n_dir);
[~, ~, dir_color_idx] = unique(direction_of);

figure('Name', 'Section 2 -- raw traces with onsets marked', 'NumberTitle', 'off', 'Position', [80 80 1000 700]);
for k = 1:n_example_cells
    subplot(n_example_cells, 1, k);
    hold on;
    plot(block_t, block_base(example_cells(k), :), 'k-', 'LineWidth', 1);
    for tr = 1:n_subtrials
        xline(onsets_rel(tr), '--', 'Color', dir_colors(dir_color_idx(tr), :), 'LineWidth', 1);
    end
    ylabel(sprintf('ROI %d', example_cells(k)), 'FontSize', 9);
    set(gca, 'Box', 'off', 'FontSize', 8); xlim([block_t(1), block_t(end)]);
    if k < n_example_cells; set(gca, 'XTickLabel', []); end
    hold off;
end
xlabel('Time (s)');
sgtitle('Section 2: raw traces with reconstructed subtrial onsets (color = direction)', 'FontWeight', 'bold');
save_fig_local('Section2_raw_traces_with_onsets', test_outdir);

%% ====================== SECTION 2b: our onset markers vs DataGathering.m's arithmetic formula ======================
% DataGathering.m (a different rig's pipeline) computes subtrial
% boundaries purely arithmetically: subtrial_start_time = (k-1) *
% (stimulus_trial_t + wait_intertrial) -- it trusts a perfectly regular
% clock and never looks at a real per-frame onset marker. Our method
% (above) reads the ACTUAL recorded onset frame for every subtrial via
% Frameinfo.subtrial_frame_count. This section checks whether that
% distinction actually matters for this data, by computing both and
% plotting them against each other.
wait_between_directions = s.specParams.wait_between_directions;
k = (1:n_subtrials)';
onsets_arithmetic = (k - 1) * (s.stimulus_trial_t + wait_between_directions);
drift = onsets_rel(:) - onsets_arithmetic;

fprintf('\nSECTION 2b: stimulus_trial_t=%.4f, wait_between_directions=%.4f\n', s.stimulus_trial_t, wait_between_directions);
fprintf('  Max drift (ours - arithmetic) across %d subtrials: %.3fs (at subtrial %d)\n', ...
    n_subtrials, max(abs(drift)), find(abs(drift) == max(abs(drift)), 1));

figure('Name', 'Section 2b -- onset marker vs arithmetic formula', 'NumberTitle', 'off', 'Position', [80 80 900 700]);
subplot(2,1,1);
hold on;
plot(k, onsets_rel, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'DisplayName', 'Ours: real onset marker (Frameinfo)');
plot(k, onsets_arithmetic, 'r--^', 'LineWidth', 1.5, 'DisplayName', 'DataGathering-style: (k-1)\times(trial\_t+wait)');
xlabel('Subtrial #'); ylabel('Onset time (s)');
legend('Location', 'best'); set(gca, 'Box', 'off');
title('Section 2b: subtrial onset time, two methods');

subplot(2,1,2);
bar(k, drift, 'FaceColor', [0.6 0.2 0.2]);
xlabel('Subtrial #'); ylabel('Drift: ours - arithmetic (s)');
set(gca, 'Box', 'off');
title(sprintf('Drift accumulates to %.3fs by subtrial %d', drift(end), n_subtrials));
save_fig_local('Section2b_onset_method_comparison', test_outdir);

%% ====================== SECTION 3: per-subtrial extraction (before averaging) ======================
% Same handful of cells, but now sliced into the 24 individual onset-aligned
% subtrials BEFORE any direction-averaging -- trial-to-trial variability
% is visible here, nothing has been smoothed over yet.
pre_window_sec = 1; post_window_sec = 6; n_com = 100;
[Ca_subtrials, t_com, direction_of_check] = extract_direction_trials(B.Stimuli, base_dff, B.TimeCa, ...
    'moving_bar', 'bar_orientations', pre_window_sec, post_window_sec, n_com);
fprintf('\nSECTION 3: Ca_subtrials = [%d cells x %d timepoints x %d subtrials]\n', size(Ca_subtrials));
assert(isequal(direction_of_check, direction_of), 'direction_of mismatch between Section 2''s manual reconstruction and extract_direction_trials -- investigate before trusting Section 4/5');

figure('Name', 'Section 3 -- individual subtrials, before averaging', 'NumberTitle', 'off', 'Position', [60 60 1500 850]);
for k = 1:n_example_cells
    subplot(n_example_cells, 1, k);
    hold on;
    cell_idx = example_cells(k);
    for tr = 1:n_subtrials
        plot(t_com, Ca_subtrials(cell_idx, :, tr), '-', 'Color', [dir_colors(dir_color_idx(tr), :), 0.5], 'LineWidth', 1);
    end
    xline(0, 'k--', 'LineWidth', 1);
    ylabel(sprintf('ROI %d', cell_idx), 'FontSize', 9);
    set(gca, 'Box', 'off', 'FontSize', 8); xlim([t_com(1), t_com(end)]);
    if k < n_example_cells; set(gca, 'XTickLabel', []); end
    hold off;
end
xlabel('Time from subtrial onset (s)');
sgtitle('Section 3: all 24 individual subtrials overlaid, color = direction (before averaging)', 'FontWeight', 'bold');
save_fig_local('Section3_subtrials_before_averaging', test_outdir);

%% ====================== SECTION 4: direction-chunked (median across repeats) ======================
[Ca_by_dir, unique_dirs] = chunk_by_direction(Ca_subtrials, direction_of);
fprintf('\nSECTION 4: Ca_by_dir = [%d cells x %d timepoints x %d directions], directions = %s\n', ...
    size(Ca_by_dir), mat2str(unique_dirs));

dir_colors2 = lines(numel(unique_dirs));
figure('Name', 'Section 4 -- direction-chunked (median of repeats)', 'NumberTitle', 'off', 'Position', [60 60 1500 850]);
for k = 1:n_example_cells
    subplot(n_example_cells, 1, k);
    hold on;
    cell_idx = example_cells(k);
    for d = 1:numel(unique_dirs)
        plot(t_com, Ca_by_dir(cell_idx, :, d), '-', 'Color', dir_colors2(d, :), 'LineWidth', 1.6, ...
            'DisplayName', sprintf('%g deg', unique_dirs(d)));
    end
    xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    ylabel(sprintf('ROI %d', cell_idx), 'FontSize', 9);
    set(gca, 'Box', 'off', 'FontSize', 8); xlim([t_com(1), t_com(end)]);
    if k < n_example_cells; set(gca, 'XTickLabel', []); end
    if k == 1; legend('Location', 'best', 'FontSize', 7); end
    hold off;
end
xlabel('Time from onset (s)');
sgtitle('Section 4: 8 direction-averaged (median of 3 repeats) traces per cell', 'FontWeight', 'bold');
save_fig_local('Section4_direction_chunked', test_outdir);

%% ====================== SECTION 5: preferred-direction collapse ======================
% Which ONE of the 8 direction-averaged traces becomes this cell's row in
% the main pipeline's resp_base / heatmap / violin / MI.
peak_window = [0.5, min(4, post_window_sec)];
peak_idx = t_com >= peak_window(1) & t_com <= peak_window(2);
amp_per_cell = squeeze(mean(Ca_by_dir(:, peak_idx, :), 2, 'omitnan'));  % n_cells x n_dir
[~, pref_dir_idx] = max(amp_per_cell, [], 2);

figure('Name', 'Section 5 -- preferred-direction collapse', 'NumberTitle', 'off', 'Position', [60 60 1500 850]);
for k = 1:n_example_cells
    subplot(n_example_cells, 1, k);
    hold on;
    cell_idx = example_cells(k);
    pdi = pref_dir_idx(cell_idx);
    for d = 1:numel(unique_dirs)
        if d == pdi
            plot(t_com, Ca_by_dir(cell_idx, :, d), '-', 'Color', dir_colors2(d, :), 'LineWidth', 3, ...
                'DisplayName', sprintf('%g deg (PREFERRED)', unique_dirs(d)));
        else
            plot(t_com, Ca_by_dir(cell_idx, :, d), '-', 'Color', [dir_colors2(d,:), 0.25], 'LineWidth', 1, ...
                'HandleVisibility', 'off');
        end
    end
    xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    ylabel(sprintf('ROI %d', cell_idx), 'FontSize', 9);
    set(gca, 'Box', 'off', 'FontSize', 8); xlim([t_com(1), t_com(end)]);
    if k < n_example_cells; set(gca, 'XTickLabel', []); end
    if k == 1; legend('Location', 'best', 'FontSize', 7); end
    hold off;
end
xlabel('Time from onset (s)');
sgtitle('Section 5: bold trace = this cell''s preferred direction (the one row that survives into the main pipeline)', 'FontWeight', 'bold');
save_fig_local('Section5_preferred_direction_collapse', test_outdir);

%% ====================== LOCAL FUNCTIONS ======================
function save_fig_local(fig_title, outdir)
    safe_name = regexprep(fig_title, '[^a-zA-Z0-9]+', '_');
    outfile = fullfile(outdir, [safe_name '.png']);
    saveas(gcf, outfile);
    fprintf('  Saved %s\n', outfile);
end
