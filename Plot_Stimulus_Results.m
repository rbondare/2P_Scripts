%% Plotting: violins, histograms, heatmaps, modulation index -- all stimuli
%
% Script 3 of 3 (Preprocess_Stimulus_Data -> analysis_functions/ -> Plot_Stimulus_Results).
% Calls Script 1 to load/match data, then analysis_functions/* (Script 2)
% per stimulus type, then produces the full figure suite:
%   - Transient (population mean +/- SEM): looming, flashes only
%   - Violins (6): {raw, z-score, raw with <0 filtered out} x {all cells, matched cells}
%   - Histograms (4): {raw, z-score} x {all cells, matched cells}
%   - Heatmaps (2): whole-duration raw (sorted by max), trial-averaged z-scored (sorted by max)
%   - Cross-stimulus: MI violin grid, MI spread histogram, MI CDF, MI linear regression
%
% grating/moving_bar deliberately have NO transient plot for now (the
% direction-resolved alignment is still being inspected) -- see
% analyze_direction_stimulus.m.

Preprocess_Stimulus_Data;
addpath(fullfile(fileparts(mfilename('fullpath')), 'analysis_functions'));

%% ====================== VERIFY DIRECTION RECONSTRUCTION ======================
fprintf('\n========== VERIFYING GRATING / MOVING_BAR DIRECTION RECONSTRUCTION ==========\n');
if any(strcmp(stimulus_types_to_analyze, 'grating')) && ...
        any(strcmp({B.Stimuli(:).type}, 'grating')) && any(strcmp({D.Stimuli(:).type}, 'grating'))
    fprintf('  Baseline:\n');
    verify_direction_reconstruction(B.Stimuli, base_dff, B.TimeCa, 'grating', 'grating_orientations');
    fprintf('  Drug:\n');
    verify_direction_reconstruction(D.Stimuli, drug_dff, D.TimeCa, 'grating', 'grating_orientations');
end
if any(strcmp(stimulus_types_to_analyze, 'moving_bar')) && ...
        any(strcmp({B.Stimuli(:).type}, 'moving_bar')) && any(strcmp({D.Stimuli(:).type}, 'moving_bar'))
    fprintf('  Baseline:\n');
    verify_direction_reconstruction(B.Stimuli, base_dff, B.TimeCa, 'moving_bar', 'bar_orientations');
    fprintf('  Drug:\n');
    verify_direction_reconstruction(D.Stimuli, drug_dff, D.TimeCa, 'moving_bar', 'bar_orientations');
end

%% ====================== PER-STIMULUS ANALYSIS + PLOTTING ======================
fprintf('\n========== PER-STIMULUS ANALYSIS ==========\n');

mi_metric_base = struct();
mi_metric_drug = struct();
raw_metric_base_all = struct();
raw_metric_drug_all = struct();

for i = 1:numel(stimulus_types_to_analyze)
    stype = stimulus_types_to_analyze{i};
    stim_exists_base = any(strcmp({B.Stimuli(:).type}, stype));
    stim_exists_drug = any(strcmp({D.Stimuli(:).type}, stype));
    if ~(stim_exists_base && stim_exists_drug)
        fprintf('Skipping "%s" (not in both baseline and drug)\n', stype);
        continue;
    end

    fprintf('Analyzing "%s"...\n', stype);
    switch stype
        case 'looming'
            result = analyze_looming_stimulus(B.Stimuli, base_dff, B.TimeCa, B.Triggers, ...
                D.Stimuli, drug_dff, D.TimeCa, D.Triggers, ...
                base_match_idx_local, drug_match_idx_local, matched_rois_available);
        case 'sparse_local_global_flashes'
            result = analyze_flashes_stimulus(B.Stimuli, base_dff, B.TimeCa, D.Stimuli, drug_dff, D.TimeCa, ...
                base_match_idx_local, drug_match_idx_local, matched_rois_available);
        case 'grating'
            result = analyze_grating_stimulus(B.Stimuli, base_dff, B.TimeCa, D.Stimuli, drug_dff, D.TimeCa, ...
                base_match_idx_local, drug_match_idx_local, matched_rois_available);
        case 'moving_bar'
            result = analyze_moving_bar_stimulus(B.Stimuli, base_dff, B.TimeCa, D.Stimuli, drug_dff, D.TimeCa, ...
                base_match_idx_local, drug_match_idx_local, matched_rois_available);
        case {'spontaneous', 'checkers2'}
            result = analyze_blockwise_stimulus(stype, B.Stimuli, base_dff, B.TimeCa, D.Stimuli, drug_dff, D.TimeCa, ...
                base_match_idx_local, drug_match_idx_local, matched_rois_available);
        otherwise
            fprintf('  (No analysis implemented for "%s" -- skipping)\n', stype);
            continue;
    end

    if isempty(result.resp_base)
        continue;  % analyze_* already printed why (no matched ROIs, stimulus missing, etc.)
    end

    stim_field = matlab.lang.makeValidName(stype);
    mi_metric_base.(stim_field) = result.z_metric_base;
    mi_metric_drug.(stim_field) = result.z_metric_drug;
    raw_metric_base_all.(stim_field) = result.raw_metric_base;
    raw_metric_drug_all.(stim_field) = result.raw_metric_drug;

    stim_outdir = fullfile(pair_outdir, stim_field);
    if ~exist(stim_outdir, 'dir'); mkdir(stim_outdir); end
    produce_stimulus_figures(result, stim_outdir);

    if ~isempty(result.common_dirs)
        produce_direction_figures(result, stim_outdir);
    end
end

%% ====================== CROSS-STIMULUS MODULATION INDEX ======================
fprintf('\n========== CROSS-STIMULUS MODULATION INDEX ==========\n');
cross_outdir = fullfile(pair_outdir, 'cross_stimulus');
if ~exist(cross_outdir, 'dir'); mkdir(cross_outdir); end

%% All-stimuli-in-one-figure violins: raw dF/F and z-score, matched cells
plot_all_stimuli_violin(raw_metric_base_all, raw_metric_drug_all, ...
    'All_stimuli_violin_raw', 'dF/F (per cell, matched)', cross_outdir);
plot_all_stimuli_violin(mi_metric_base, mi_metric_drug, ...
    'All_stimuli_violin_zscore', 'z-score (per cell, matched)', cross_outdir);

stim_names = fieldnames(mi_metric_base);
mi_by_stim = {};
stim_labels = {};
for i = 1:numel(stim_names)
    f = stim_names{i};
    zb = mi_metric_base.(f); zd = mi_metric_drug.(f);
    if numel(zb) < 2 || numel(zd) < 2; continue; end
    mi = compute_modulation_index(zb, zd);
    if numel(mi) < 2; continue; end
    mi_by_stim{end+1} = mi; %#ok<AGROW>
    stim_labels{end+1} = f; %#ok<AGROW>
    fprintf('  %-12s: n=%d, %.0f%% up, %.0f%% down, median MI=%.3f\n', f, numel(mi), ...
        100*mean(mi>0), 100*mean(mi<0), median(mi));
end

if numel(mi_by_stim) >= 1
    plot_mi_violin_grid(mi_by_stim, stim_labels, cross_outdir);
    plot_mi_spread_histogram(mi_by_stim, stim_labels, cross_outdir);
    plot_mi_cdf(mi_by_stim, stim_labels, cross_outdir);
end

%% Linear regression: baseline raw mean dF/F vs drug raw mean dF/F, matched cells, all stimuli together
fprintf('\n========== MODULATION REGRESSION (baseline vs drug raw mean dF/F) ==========\n');
plot_mi_regression_grid(stimulus_types_to_analyze, B, D, base_dff, drug_dff, ...
    base_match_idx_local, drug_match_idx_local, matched_rois_available, cross_outdir);

fprintf('\nAll figures saved under: %s\n', pair_outdir);


%% ====================== LOCAL FUNCTIONS ======================

function produce_stimulus_figures(result, outdir)
    label = result.stim_label;

    if result.has_transient
        plot_transient(result.t_com, result.resp_base, result.resp_drug, ...
            sprintf('%s_transient', label), outdir);
    end

    % ---- Violins (6): {raw, z-score, raw positive-filtered} x {matched, all cells} ----
    plot_violin(result.raw_metric_base, result.raw_metric_drug, ...
        sprintf('%s_violin_raw_matched', label), 'dF/F (per cell)', outdir, true);
    plot_violin(result.raw_all_base, result.raw_all_drug, ...
        sprintf('%s_violin_raw_allcells', label), 'dF/F (per cell)', outdir, false);

    plot_violin(result.z_metric_base, result.z_metric_drug, ...
        sprintf('%s_violin_zscore_matched', label), 'z-score (per cell)', outdir, true);
    plot_violin(result.z_all_base, result.z_all_drug, ...
        sprintf('%s_violin_zscore_allcells', label), 'z-score (per cell)', outdir, false);

    pos_matched = result.raw_metric_base >= 0 & result.raw_metric_drug >= 0;
    plot_violin(result.raw_metric_base(pos_matched), result.raw_metric_drug(pos_matched), ...
        sprintf('%s_violin_rawpositive_matched', label), 'dF/F (per cell, >=0 only)', outdir, true);
    plot_violin(result.raw_all_base(result.raw_all_base >= 0), result.raw_all_drug(result.raw_all_drug >= 0), ...
        sprintf('%s_violin_rawpositive_allcells', label), 'dF/F (per cell, >=0 only)', outdir, false);

    % ---- Histograms (4): {raw, z-score} x {matched, all cells} ----
    plot_histogram(result.raw_metric_base, result.raw_metric_drug, ...
        sprintf('%s_hist_raw_matched', label), 'dF/F (per cell)', outdir);
    plot_histogram(result.raw_all_base, result.raw_all_drug, ...
        sprintf('%s_hist_raw_allcells', label), 'dF/F (per cell)', outdir);
    plot_histogram(result.z_metric_base, result.z_metric_drug, ...
        sprintf('%s_hist_zscore_matched', label), 'z-score (per cell)', outdir);
    plot_histogram(result.z_all_base, result.z_all_drug, ...
        sprintf('%s_hist_zscore_allcells', label), 'z-score (per cell)', outdir);

    % ---- Heatmaps (4): {raw, z-score} x {sorted by baseline peak, unsorted} -- matched cells ----
    plot_heatmap(result.t_com, result.resp_base, result.resp_drug, ...
        sprintf('%s_heatmap_raw_sorted', label), outdir, 'dF/F', false, 'baseline');
    plot_heatmap(result.t_com, result.resp_base, result.resp_drug, ...
        sprintf('%s_heatmap_raw_unsorted', label), outdir, 'dF/F', false, 'none');
    plot_heatmap(result.t_com, result.Ztrace_base, result.Ztrace_drug, ...
        sprintf('%s_heatmap_zscore_sorted', label), outdir, 'z-score', true, 'baseline');
    plot_heatmap(result.t_com, result.Ztrace_base, result.Ztrace_drug, ...
        sprintf('%s_heatmap_zscore_unsorted', label), outdir, 'z-score', true, 'none');
end

function produce_direction_figures(result, outdir)
    % grating/moving_bar only: how the response looks at each orientation,
    % not just collapsed to each cell's preferred direction.
    label = result.stim_label;
    plot_tuning_curve(result.common_dirs, result.amp_base, result.amp_drug, ...
        sprintf('%s_tuning_curve', label), outdir, result.n_matched);
    plot_direction_grid(result.common_dirs, result.t_com, ...
        result.dir_mean_base, result.dir_sem_base, result.dir_mean_drug, result.dir_sem_drug, ...
        sprintf('%s_response_by_orientation', label), outdir, result.n_matched);

    % Continuous whole-block heatmap: the entire stimulus block (all
    % directions/repeats in their natural order), not collapsed to each
    % cell's preferred direction like resp_base/resp_drug above -- gives
    % the population-level repeat structure a per-direction collapse
    % erases.
    if ~isempty(result.resp_base_block)
        plot_heatmap(result.t_com_block, result.resp_base_block, result.resp_drug_block, ...
            sprintf('%s_heatmap_block_raw_sorted', label), outdir, 'dF/F', false, 'baseline');
        plot_heatmap(result.t_com_block, result.resp_base_block, result.resp_drug_block, ...
            sprintf('%s_heatmap_block_raw_unsorted', label), outdir, 'dF/F', false, 'none');
        plot_heatmap(result.t_com_block, result.resp_base_block, result.resp_drug_block, ...
            sprintf('%s_heatmap_block_raw_sorted_by_drug', label), outdir, 'dF/F', false, 'drug');
        plot_heatmap(result.t_com_block, result.Ztrace_base_block, result.Ztrace_drug_block, ...
            sprintf('%s_heatmap_block_zscore_sorted', label), outdir, 'z-score', true, 'baseline');
        plot_heatmap(result.t_com_block, result.Ztrace_base_block, result.Ztrace_drug_block, ...
            sprintf('%s_heatmap_block_zscore_unsorted', label), outdir, 'z-score', true, 'none');
    end

    % All cells (unmatched/full population), same continuous block, raw,
    % unsorted -- the un-curated population, not just the matched subset.
    if ~isempty(result.resp_base_block_all)
        plot_heatmap(result.t_com_block, result.resp_base_block_all, result.resp_drug_block_all, ...
            sprintf('%s_heatmap_block_raw_allcells_unsorted', label), outdir, 'dF/F', false, 'none');
    end

    % Trial-blocks-by-direction: subtrials reordered/concatenated so all
    % repeats of one direction sit adjacent, then the next direction --
    % trades real elapsed time for direction-tuning structure being
    % visible as distinct blocks.
    if ~isempty(result.resp_base_trialblocks)
        plot_trialblocks_heatmap(result.resp_base_trialblocks, result.resp_drug_trialblocks, ...
            result.trialblock_boundaries, result.trialblock_dir_labels, ...
            sprintf('%s_heatmap_trialblocks_by_direction', label), outdir);
    end
end

function plot_tuning_curve(common_dirs, amp_base, amp_drug, fig_title, outdir, n_matched)
    figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [100 100 700 500]);
    plot(common_dirs, amp_base, 'o-', 'Color', [0.2 0.6 1], 'LineWidth', 2, ...
        'MarkerFaceColor', [0.2 0.6 1], 'DisplayName', 'Baseline'); hold on;
    plot(common_dirs, amp_drug, 'o-', 'Color', [1 0.5 0.2], 'LineWidth', 2, ...
        'MarkerFaceColor', [1 0.5 0.2], 'DisplayName', 'Drug');
    xlabel('Direction (deg)', 'FontSize', 11); ylabel('Mean dF/F (peak window, population)', 'FontSize', 11);
    title(sprintf('%s (n=%d matched cells)', fig_title, n_matched), 'FontSize', 12, 'Interpreter', 'none');
    legend('Location', 'best'); set(gca, 'Box', 'off'); xticks(common_dirs);
    save_fig(fig_title, outdir);
end

function plot_direction_grid(common_dirs, t_com, dir_mean_base, dir_sem_base, dir_mean_drug, dir_sem_drug, fig_title, outdir, n_matched)
    % Small multiples: one subplot per orientation, population mean +/-
    % SEM, baseline vs drug overlaid -- the actual time-resolved response
    % shape at each direction, not just a collapsed amplitude.
    n_dir = numel(common_dirs);
    grid_cols = ceil(sqrt(n_dir));
    grid_rows = ceil(n_dir / grid_cols);
    figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [50 50 min(280*grid_cols, 1800) min(260*grid_rows, 1000)]);
    for d = 1:n_dir
        subplot(grid_rows, grid_cols, d);
        hold on;
        x_sh = [t_com, fliplr(t_com)];
        mb = dir_mean_base(d, :); sb = dir_sem_base(d, :);
        md = dir_mean_drug(d, :); sd = dir_sem_drug(d, :);
        fill(x_sh, [mb+sb, fliplr(mb-sb)], [0.2 0.6 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(t_com, mb, '-', 'Color', [0.2 0.6 1], 'LineWidth', 1.8, 'DisplayName', 'Baseline');
        fill(x_sh, [md+sd, fliplr(md-sd)], [1 0.5 0.2], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(t_com, md, '-', 'Color', [1 0.5 0.2], 'LineWidth', 1.8, 'DisplayName', 'Drug');
        xline(0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        title(sprintf('%g deg', common_dirs(d)), 'FontSize', 10);
        set(gca, 'Box', 'off', 'FontSize', 8); xlim([t_com(1), t_com(end)]);
        if d == 1; legend('Location', 'best', 'FontSize', 7); end
    end
    sgtitle(sprintf('%s (n=%d matched cells)', fig_title, n_matched), 'FontWeight', 'bold', 'Interpreter', 'none');
    save_fig(fig_title, outdir);
end

function save_fig(fig_title, outdir)
    safe_name = regexprep(fig_title, '[^a-zA-Z0-9]+', '_');
    outfile = fullfile(outdir, [safe_name '.png']);
    saveas(gcf, outfile);
    fprintf('      Saved %s\n', outfile);
end

function plot_transient(t_com, resp_base, resp_drug, fig_title, outdir)
    mean_base = mean(resp_base, 1, 'omitnan');
    sem_base  = std(resp_base, 0, 1, 'omitnan') / sqrt(size(resp_base, 1));
    mean_drug = mean(resp_drug, 1, 'omitnan');
    sem_drug  = std(resp_drug, 0, 1, 'omitnan') / sqrt(size(resp_drug, 1));

    figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [100 100 900 500]);
    hold on;
    x_sh = [t_com, fliplr(t_com)];
    y_sh_base = [(mean_base + sem_base), fliplr(mean_base - sem_base)];
    fill(x_sh, y_sh_base, [0.2 0.6 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t_com, mean_base, '-', 'Color', [0.2 0.6 1], 'LineWidth', 2.5, 'DisplayName', 'Baseline');
    y_sh_drug = [(mean_drug + sem_drug), fliplr(mean_drug - sem_drug)];
    fill(x_sh, y_sh_drug, [1 0.5 0.2], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t_com, mean_drug, '-', 'Color', [1 0.5 0.2], 'LineWidth', 2.5, 'DisplayName', 'Drug');
    xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlabel('Time from onset (s)', 'FontSize', 11);
    ylabel('dF/F (mean \pm SEM)', 'FontSize', 11);
    title(sprintf('%s (n=%d matched cells)', fig_title, size(resp_base,1)), 'FontSize', 12, 'Interpreter', 'none');
    legend('Location', 'best'); set(gca, 'Box', 'off'); xlim([t_com(1), t_com(end)]);
    save_fig(fig_title, outdir);
end

function plot_violin(vals_base, vals_drug, fig_title, ylab, outdir, paired)
    vals_base = vals_base(:); vals_drug = vals_drug(:);
    if paired
        valid = ~isnan(vals_base) & ~isnan(vals_drug);
        vals_base = vals_base(valid); vals_drug = vals_drug(valid);
    else
        vals_base = vals_base(~isnan(vals_base));
        vals_drug = vals_drug(~isnan(vals_drug));
    end
    if numel(vals_base) < 2 || numel(vals_drug) < 2
        fprintf('      (Skipping %s -- not enough data)\n', fig_title);
        return;
    end

    color_base = [0.6 0.6 0.6];      % grey
    color_drug = [0.40 0.65 0.40];   % green

    figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [100 100 600 650]);
    hold on;
    plot_one_violin(0, vals_base, color_base);
    plot_one_violin(1, vals_drug, color_drug);

    if paired && numel(vals_base) == numel(vals_drug)
        p = signrank(vals_base, vals_drug);
        test_name = 'signrank';
    else
        p = ranksum(vals_base, vals_drug);
        test_name = 'ranksum';
    end
    y_top  = max([vals_base; vals_drug]);
    y_span = range([vals_base; vals_drug]);
    y_bracket = y_top + 0.08 * y_span;
    plot([0 0 1 1], [y_bracket-0.02*y_span, y_bracket, y_bracket, y_bracket-0.02*y_span], ...
        'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
    text(0.5, y_bracket + 0.03*y_span, sig_stars(p), ...
        'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
    fprintf('      %s: %s p = %.4g (n_base=%d, n_drug=%d)\n', fig_title, test_name, p, numel(vals_base), numel(vals_drug));

    set(gca, 'XTick', [0 1], 'XTickLabel', {'Baseline', 'Drug'}, 'Box', 'off', 'XLim', [-0.6 1.6]);
    ylabel(ylab, 'FontSize', 11);
    title(fig_title, 'FontSize', 12, 'Interpreter', 'none');
    save_fig(fig_title, outdir);
end

function plot_one_violin(x0, vals, fill_color)
    if numel(vals) < 2; return; end
    [f, xi] = ksdensity(vals, 'NumPoints', 150);
    f = f / max(f) * 0.38;
    patch(x0 + [f, fliplr(-f)], [xi, fliplr(xi)], fill_color, 'FaceAlpha', 0.85, 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 1);
    plot(x0 + [-0.38 0.38], [median(vals) median(vals)], 'w-', 'LineWidth', 2.5);
    x_jitter = x0 + max(min(randn(numel(vals), 1)*0.10, 0.35), -0.35);
    scatter(x_jitter, vals, 14, [0.25 0.25 0.25], 'o', 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0);
end

function stars = sig_stars(p)
    if p < 0.001; stars = '***';
    elseif p < 0.01; stars = '**';
    elseif p < 0.05; stars = '*';
    else; stars = 'n.s.';
    end
end

function plot_histogram(vals_base, vals_drug, fig_title, xlab, outdir)
    vals_base = vals_base(~isnan(vals_base));
    vals_drug = vals_drug(~isnan(vals_drug));
    if isempty(vals_base) && isempty(vals_drug); return; end
    n_bins = 40;
    edges = linspace(min([vals_base; vals_drug]), max([vals_base; vals_drug]), n_bins+1);

    figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [100 100 700 450]);
    hold on;
    histogram(vals_base, edges, 'Normalization', 'probability', ...
        'FaceColor', [0.2 0.6 1], 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', sprintf('Baseline (n=%d)', numel(vals_base)));
    histogram(vals_drug, edges, 'Normalization', 'probability', ...
        'FaceColor', [1 0.5 0.2], 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', sprintf('Drug (n=%d)', numel(vals_drug)));
    xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlabel(xlab, 'FontSize', 11); ylabel('Probability', 'FontSize', 11);
    title(fig_title, 'FontSize', 12, 'Interpreter', 'none');
    legend('Location', 'best'); set(gca, 'Box', 'off');
    save_fig(fig_title, outdir);
end

function plot_heatmap(t_com, resp_base, resp_drug, fig_title, outdir, cbar_label, center_zero, sort_mode)
    % center_zero (z-score) data is sign-balanced around a meaningful 0, so
    % it gets the diverging blue-white-red colormap with limits symmetric
    % about 0. Raw dF/F is NOT sign-balanced (mostly-positive activity with
    % no meaningful "negative" side) -- a diverging map puts white at the
    % MIDDLE of the percentile range rather than at 0, which makes
    % quiescent baseline activity read as blue instead of neutral. Raw gets
    % a plain grayscale (white=low, black=high) sequential map instead --
    % matches Master_Stimulus_Analysis.m's original convention (flipud(gray)
    % for raw heatmaps, diverging only for z-score).
    % sort_mode = 'baseline': rows sorted by baseline's own peak value
    % (descending), that order carried to drug. 'drug': rows sorted by
    % drug's own peak instead (its own, independent ranking -- NOT the
    % same row-i-is-the-same-cell correspondence as 'baseline', since
    % baseline and drug panels would then show DIFFERENT cell orders; only
    % meant to be viewed side by side with a 'baseline'-sorted version of
    % the same data, not as a single internally-consistent figure).
    % 'none': natural matched-cell order, same row = same cell in both panels.
    if isempty(resp_base) || isempty(resp_drug); return; end
    switch sort_mode
        case 'baseline'
            base_max = max(resp_base, [], 2); base_max(isnan(base_max)) = -Inf;
            [~, order_base] = sort(base_max, 'descend');
            order_drug = order_base;
            ylab = 'Matched cell # (sorted by baseline peak)';
        case 'drug'
            drug_max = max(resp_drug, [], 2); drug_max(isnan(drug_max)) = -Inf;
            [~, order_drug] = sort(drug_max, 'descend');
            order_base = order_drug;
            ylab = 'Matched cell # (sorted by drug peak)';
        otherwise
            order_base = 1:size(resp_base, 1);
            order_drug = order_base;
            ylab = 'Matched cell # (unsorted)';
    end
    Zb = resp_base(order_base, :);
    Zd = resp_drug(order_drug, :);

    if center_zero
        clim = prctile([Zb(:); Zd(:)], [1 99]);
        clim = max(abs(clim)) * [-1 1];
        if any(~isfinite(clim)) || clim(1) == clim(2); clim = [-1 1]; end
    else
        clim = [0, prctile([Zb(:); Zd(:)], 99)];
        if any(~isfinite(clim)) || clim(1) == clim(2); clim = [0 1]; end
    end

    figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [100 100 1100 650]);
    if center_zero
        colormap(diverging_colormap());
        onset_color = 'k--';
    else
        colormap(flipud(gray));
        onset_color = 'r--';   % black would blend into high-activity (dark) cells on grayscale
    end

    ax1 = subplot(1,2,1);
    imagesc(t_com, 1:size(Zb,1), Zb, clim);
    hold on; xline(0, onset_color, 'LineWidth', 1); hold off;
    xlabel('Time / position'); ylabel(ylab); title('Baseline');
    pos1 = get(ax1, 'Position');

    ax2 = subplot(1,2,2);
    imagesc(t_com, 1:size(Zd,1), Zd, clim);
    hold on; xline(0, onset_color, 'LineWidth', 1); hold off;
    xlabel('Time / position'); title('Drug');
    pos2 = get(ax2, 'Position');   % capture BEFORE adding the colorbar
    cb = colorbar;
    cb.Label.String = cbar_label;
    set(ax2, 'Position', pos2);    % colorbar() shrinks the axes it's attached to -- restore
    set(ax1, 'Position', [pos1(1) pos2(2) pos2(3) pos2(4)]);  % match Baseline to Drug's (now-correct) size

    sgtitle(sprintf('%s (n=%d matched cells)', fig_title, size(Zb,1)), 'FontWeight', 'bold', 'Interpreter', 'none');
    save_fig(fig_title, outdir);
end

function plot_trialblocks_heatmap(resp_base, resp_drug, boundaries, dir_labels, fig_title, outdir)
    % X-axis is column index, NOT real elapsed time -- each direction's
    % repeats are concatenated end-to-end (reorder_subtrials_by_direction.m
    % built resp_base/resp_drug this way), with vertical separators +
    % direction labels marking where each direction's block starts.
    % Always sorted by baseline peak (matched cells only, raw dF/F).
    if isempty(resp_base) || isempty(resp_drug); return; end
    base_max = max(resp_base, [], 2); base_max(isnan(base_max)) = -Inf;
    [~, order] = sort(base_max, 'descend');
    Zb = resp_base(order, :);
    Zd = resp_drug(order, :);
    clim = [0, prctile([Zb(:); Zd(:)], 99)];
    if any(~isfinite(clim)) || clim(1) == clim(2); clim = [0 1]; end

    x_axis = 1:size(Zb, 2);
    xtick_labels = arrayfun(@(d) sprintf('%g deg', d), dir_labels, 'UniformOutput', false);

    figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [100 100 1200 650]);
    colormap(flipud(gray));

    ax1 = subplot(1,2,1);
    imagesc(x_axis, 1:size(Zb,1), Zb, clim);
    hold on; for b = boundaries(2:end); xline(b-0.5, 'r--', 'LineWidth', 1); end; hold off;
    set(gca, 'XTick', boundaries, 'XTickLabel', xtick_labels); xtickangle(45);
    ylabel('Matched cell # (sorted by baseline peak)'); title('Baseline');
    pos1 = get(ax1, 'Position');

    ax2 = subplot(1,2,2);
    imagesc(x_axis, 1:size(Zd,1), Zd, clim);
    hold on; for b = boundaries(2:end); xline(b-0.5, 'r--', 'LineWidth', 1); end; hold off;
    set(gca, 'XTick', boundaries, 'XTickLabel', xtick_labels); xtickangle(45);
    title('Drug');
    pos2 = get(ax2, 'Position');
    cb = colorbar; cb.Label.String = 'dF/F';
    set(ax2, 'Position', pos2);
    set(ax1, 'Position', [pos1(1) pos2(2) pos2(3) pos2(4)]);

    sgtitle(sprintf('%s (n=%d matched cells, each block = repeats of 1 direction)', fig_title, size(Zb,1)), ...
        'FontWeight', 'bold', 'Interpreter', 'none');
    save_fig(fig_title, outdir);
end

function cmap = diverging_colormap(n)
    % Blue (low) -> white (mid) -> red (high), matching the colormap
    % Loom_Flashes_ZETA_Analysis.m has always used.
    if nargin < 1; n = 256; end
    half = floor(n/2);
    neg = [linspace(0.05,1,half)', linspace(0.05,1,half)', ones(half,1)];      % blue -> white
    pos = [ones(n-half,1), linspace(1,0.05,n-half)', linspace(1,0.05,n-half)']; % white -> red
    cmap = [neg; pos];
end

function plot_all_stimuli_violin(metric_base_struct, metric_drug_struct, fig_title, ylab, outdir)
    % One figure, all stimuli together, each as a Baseline|Drug violin
    % pair -- same metric/aesthetic as the per-stimulus violins, just
    % combined so all stimuli can be compared side by side at a glance.
    stim_names = fieldnames(metric_base_struct);
    color_base = [0.6 0.6 0.6];      % grey
    color_drug = [0.40 0.65 0.40];   % green

    figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [100 100 max(320*numel(stim_names), 700) 650]);
    hold on;
    x0 = 0; x_tick_pos = []; x_tick_lab = {};
    for s = 1:numel(stim_names)
        f = stim_names{s};
        vb = metric_base_struct.(f)(:); vd = metric_drug_struct.(f)(:);
        valid = ~isnan(vb) & ~isnan(vd);
        vb = vb(valid); vd = vd(valid);
        if numel(vb) < 2 || numel(vd) < 2; continue; end

        plot_one_violin(x0, vb, color_base);
        plot_one_violin(x0 + 1, vd, color_drug);
        if numel(vb) == numel(vd)
            p = signrank(vb, vd);
            y_top = max([vb; vd]); y_span = range([vb; vd]);
            y_bracket = y_top + 0.08 * y_span;
            plot(x0 + [0 0 1 1], [y_bracket-0.02*y_span, y_bracket, y_bracket, y_bracket-0.02*y_span], ...
                'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
            text(x0 + 0.5, y_bracket + 0.04*y_span, sig_stars(p), ...
                'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
        end
        x_tick_pos(end+1) = x0 + 0.5; %#ok<AGROW>
        x_tick_lab{end+1} = f; %#ok<AGROW>
        x0 = x0 + 2.6;
    end
    yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    set(gca, 'XTick', x_tick_pos, 'XTickLabel', x_tick_lab, 'Box', 'off', 'TickLabelInterpreter', 'none');
    ylabel(ylab, 'FontSize', 11);
    title(fig_title, 'FontSize', 12, 'Interpreter', 'none');
    h_base = patch(NaN, NaN, color_base);
    h_drug = patch(NaN, NaN, color_drug);
    legend([h_base, h_drug], {'Baseline', 'Drug'}, 'Location', 'best');
    save_fig(fig_title, outdir);
end

function plot_mi_violin_grid(mi_by_stim, stim_labels, outdir)
    n_stim = numel(mi_by_stim);
    colors = [0.55 0.75 0.70; 0.70 0.60 0.78; 0.85 0.65 0.45; 0.55 0.70 0.85; 0.75 0.75 0.55; 0.65 0.80 0.55];

    figure('Name', 'Modulation index -- violin', 'NumberTitle', 'off', 'Position', [100 100 max(250*n_stim,600) 650]);
    hold on;
    for s = 1:n_stim
        mi = mi_by_stim{s};
        if numel(mi) < 2; continue; end
        plot_one_violin((s-1)*1.3, mi, colors(mod(s-1,size(colors,1))+1, :));
    end
    yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    set(gca, 'XTick', (0:n_stim-1)*1.3, 'XTickLabel', stim_labels, 'Box', 'off', 'TickLabelInterpreter', 'none');
    ylabel('Modulation index (drug vs baseline)', 'FontSize', 11);
    title('Modulation index per stimulus (bounded [-1, 1])', 'FontSize', 12);
    save_fig('Modulation_index_violin', outdir);
end

function plot_mi_spread_histogram(mi_by_stim, stim_labels, outdir)
    n_stim = numel(mi_by_stim);
    figure('Name', 'Modulation index -- spread histogram', 'NumberTitle', 'off', ...
        'Position', [100 100 min(420*n_stim, 1800) 480]);
    edges = linspace(-1, 1, 31);
    for s = 1:n_stim
        mi = mi_by_stim{s};
        subplot(1, n_stim, s);
        hold on;
        fill([-1 0 0 -1], [0 0 1 1], [0.20 0.45 0.80], 'FaceAlpha', 0.07, 'EdgeColor', 'none');
        fill([0 1 1 0],   [0 0 1 1], [0.85 0.25 0.15], 'FaceAlpha', 0.07, 'EdgeColor', 'none');
        if numel(mi) >= 2
            histogram(mi, edges, 'Normalization', 'probability', ...
                'FaceColor', [0.4 0.4 0.4], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
        end
        xline(0, 'k--', 'LineWidth', 1.5);
        pct_up = 100*mean(mi>0); pct_down = 100*mean(mi<0);
        xlim([-1 1]); xlabel('Modulation index', 'FontSize', 10); ylabel('Probability', 'FontSize', 10);
        title(sprintf('%s\n%.0f%% up | %.0f%% down', strrep(stim_labels{s},'_','\_'), pct_up, pct_down), ...
            'FontSize', 10, 'FontWeight', 'bold');
        grid on; box off;
        hold off;
    end
    sgtitle('Modulation index spread (bounded [-1, 1])', 'FontSize', 13, 'FontWeight', 'bold');
    save_fig('Modulation_index_spread_histogram', outdir);
end

function plot_mi_cdf(mi_by_stim, stim_labels, outdir)
    n_stim = numel(mi_by_stim);
    figure('Position', [100 50 min(420*n_stim, 1800) 540], ...
        'NumberTitle', 'off', 'Name', 'ModulationIndex_CDF', 'Color', 'w');
    for s = 1:n_stim
        mi = mi_by_stim{s};
        mi_sorted = sort(mi);
        cdf_y = (1:numel(mi_sorted))' / numel(mi_sorted);

        subplot(1, n_stim, s);
        hold on;
        fill([-1 0 0 -1], [0 0 1 1], [0.20 0.45 0.80], 'FaceAlpha', 0.07, 'EdgeColor', 'none');
        fill([0 1 1 0],   [0 0 1 1], [0.85 0.25 0.15], 'FaceAlpha', 0.07, 'EdgeColor', 'none');
        plot(mi_sorted, cdf_y, 'k-', 'LineWidth', 2.5);
        xline(0, 'r--', 'LineWidth', 1.8);
        yline(0.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
        hold off;
        xlabel('Modulation Index', 'FontSize', 10); ylabel('Cumulative fraction', 'FontSize', 10);
        xlim([-1 1]); ylim([0 1]);
        title(sprintf('%s\n%.0f%% up | %.0f%% down', strrep(stim_labels{s},'_','\_'), 100*mean(mi>0), 100*mean(mi<0)), ...
            'FontSize', 10, 'FontWeight', 'bold');
        grid on; box off; set(gca, 'LineWidth', 1.2, 'FontSize', 10);
    end
    sgtitle('Modulation Index CDF -- (drug - baseline) / (|drug| + |baseline|)', 'FontSize', 13, 'FontWeight', 'bold');
    save_fig('Modulation_index_CDF', outdir);
end

function plot_mi_regression_grid(stimulus_types_to_analyze, B, D, base_dff, drug_dff, ...
        base_match_idx_local, drug_match_idx_local, matched_rois_available, outdir)
    if ~matched_rois_available; return; end

    raw_base_by_stim = {}; raw_drug_by_stim = {}; stim_labels = {};
    for i = 1:numel(stimulus_types_to_analyze)
        stype = stimulus_types_to_analyze{i};
        if ~any(strcmp({B.Stimuli(:).type}, stype)) || ~any(strcmp({D.Stimuli(:).type}, stype))
            continue;
        end
        resp_base_cell = extract_full_stimulus_responses(B.Stimuli, base_dff, stype, B.TimeCa);
        resp_drug_cell = extract_full_stimulus_responses(D.Stimuli, drug_dff, stype, D.TimeCa);
        if isempty(resp_base_cell) || isempty(resp_drug_cell); continue; end

        base_per_neuron = [];
        for p = 1:numel(resp_base_cell)
            base_per_neuron = [base_per_neuron, mean(resp_base_cell{p}(base_match_idx_local, :), 2)]; %#ok<AGROW>
        end
        drug_per_neuron = [];
        for p = 1:numel(resp_drug_cell)
            drug_per_neuron = [drug_per_neuron, mean(resp_drug_cell{p}(drug_match_idx_local, :), 2)]; %#ok<AGROW>
        end
        base_per_neuron = mean(base_per_neuron, 2, 'omitnan');
        drug_per_neuron = mean(drug_per_neuron, 2, 'omitnan');

        ok = ~isnan(base_per_neuron) & ~isinf(base_per_neuron) & ~isnan(drug_per_neuron) & ~isinf(drug_per_neuron);
        if sum(ok) < 2; continue; end
        raw_base_by_stim{end+1} = base_per_neuron(ok); %#ok<AGROW>
        raw_drug_by_stim{end+1} = drug_per_neuron(ok); %#ok<AGROW>
        stim_labels{end+1} = stype; %#ok<AGROW>
    end

    n_grid = numel(stim_labels);
    if n_grid == 0; return; end
    grid_cols = ceil(sqrt(n_grid));
    grid_rows = ceil(n_grid / grid_cols);

    figure('Position', [100 100 1400 900], 'NumberTitle', 'off', 'Name', 'Baseline_vs_Drug_Modulation_Regression');
    for s = 1:n_grid
        bp = raw_base_by_stim{s}; dp = raw_drug_by_stim{s};
        subplot(grid_rows, grid_cols, s);
        scatter(bp, dp, 100, 'o', 'filled', 'MarkerFaceColor', [0.3 0.5 0.8], 'MarkerEdgeColor', 'black', ...
            'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.8, 'LineWidth', 1.5);
        hold on;
        coeffs = polyfit(bp, dp, 1);
        y_fit = polyval(coeffs, bp);
        ss_res = sum((dp - y_fit) .^ 2);
        ss_tot = sum((dp - mean(dp)) .^ 2);
        r_squared = 1 - (ss_res / ss_tot);
        x_range = [min(bp) - 0.5, max(bp) + 0.5];
        plot(x_range, polyval(coeffs, x_range), 'b-', 'LineWidth', 2.5, ...
            'DisplayName', sprintf('y = %.3f*x + %.3f (R^2 = %.3f)', coeffs(1), coeffs(2), r_squared));
        y_min = min([bp; dp]) - 0.5; y_max = max([bp; dp]) + 0.5;
        plot([y_min y_max], [y_min y_max], 'k--', 'LineWidth', 1.5, 'DisplayName', 'No change (y=x)');
        xlabel('Baseline dF/F', 'FontSize', 11); ylabel('Drug dF/F', 'FontSize', 11);
        title(sprintf('%s (n=%d)', strrep(stim_labels{s},'_','\_'), numel(bp)), 'FontSize', 12, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 8);
        grid on; set(gca, 'LineWidth', 1.2, 'FontSize', 9);
        axis equal; xlim([y_min y_max]); ylim([y_min y_max]);
        hold off;
        fprintf('  %-12s: slope=%.3f, R^2=%.3f\n', stim_labels{s}, coeffs(1), r_squared);
    end
    sgtitle('Baseline vs Drug -- Matched Neurons (whole-stimulus mean dF/F)', 'FontSize', 14, 'FontWeight', 'bold');
    save_fig('Modulation_regression_grid', outdir);
end
