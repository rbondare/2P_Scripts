%% Loom_Flashes_ZETA_Analysis.m
% Matched-cell baseline-vs-drug figures for looming and
% sparse_local_global_flashes. Saves everything into
% analysis_figures/<baseline>_vs_<drug>/ next to this script. Per
% stimulus:
%   1. Raw dF/F transient (time-resolved, population mean +/- SEM)
%   2. Histogram of per-cell peak-window response -- raw dF/F, and as a
%      baseline z-score
%   3. Violin of per-cell peak-window response -- raw dF/F, raw dF/F with
%      negative cells removed, and as a baseline z-score
%   4. Per-cell z-scored heatmap (matched cells, baseline | drug side by
%      side, same row order, shared color scale)
%
% Why z-score at all: raw windowed dF/F isn't zero-centered, so a violin/
% histogram of raw per-cell values piles up near (and below) zero -- see
% "Calcium Imaging and the Curse of Negativity" (Frontiers 2020).
% Z-scoring each cell against its own pre-onset baseline noise is the
% standard fix. Kept the raw-dF/F versions alongside it on request.
%
% No ZETA dependency. Kept separate from Master_Stimulus_Analysis.m until
% validated. See project memory "project-stimulus-timing-architecture" and
% "project-aggregation-pipeline-fixes" for the timing-fix background.

%% ====================== CONFIGURATION ======================
baseline_file  = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260320_1444_preprocessed.mat";
drug_file      = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260320_1609_preprocessed.mat";
roi_match_file = "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane0_1444_V1.mat";
selected_plane_idx = 1;

loom_pre  = 2; loom_post  = 5;    % trial-aligned window, looming (seconds around onset)
flash_pre = 1; flash_post = 4;    % trial-aligned window, flashes
n_com     = 100;                 % points in the common trial-aligned time grid
n_bins    = 40;                  % histogram bin count (finer than before, to see the shift)

% Per-cell peak window: matches the actual response transient, NOT the
% full post-onset window. A full-window mean cancels a peak-then-undershoot
% response toward zero (confirmed on flashes: 118/649 cells had a strong
% 0-0.8s peak >0.3 but a near-zero 0-4s mean) -- see
% project_aggregation_pipeline_fixes memory.
loom_peak_window  = [0.5 1.5];   % looming's sharp peak sits here (checked via population trace)
flash_peak_window = [0   0.8];   % flashes' peak before the undershoot begins

outdir = fullfile(fileparts(mfilename('fullpath')), 'analysis_figures', ...
    sprintf('%s_vs_%s', recording_name_from_path(baseline_file), recording_name_from_path(drug_file)));
if ~exist(outdir, 'dir'); mkdir(outdir); end
fprintf('Figures will be saved to: %s\n', outdir);

%% ====================== LOAD + PLANE SUBSET + ROI MATCHING ======================
fprintf('Loading baseline: %s\n', baseline_file);
B = load(baseline_file);
fprintf('Loading drug: %s\n', drug_file);
D = load(drug_file);

[base_dff, base_roi_idx] = subset_plane(B, selected_plane_idx);
[drug_dff, drug_roi_idx] = subset_plane(D, selected_plane_idx);

[base_match_idx_local, drug_match_idx_local, n_matched] = match_rois(roi_match_file, base_roi_idx, drug_roi_idx);
fprintf('Matched cells in selected plane: %d\n', n_matched);

%% ====================== ONSETS (red-frame ground truth) ======================
% TimeStimulusFrame is directly camera-clock-correct now (RedFrameSynchronized=true
% for every stimulus block, including looming) -- no manual anchor-correction
% against TimeProjector is needed for either stimulus.
loom_onsets_base  = get_looming_onsets(B.Stimuli, B.Triggers);
loom_onsets_drug  = get_looming_onsets(D.Stimuli, D.Triggers);
flash_onsets_base = get_flash_onsets(B.Stimuli);
flash_onsets_drug = get_flash_onsets(D.Stimuli);

fprintf('Looming onsets: baseline=%d, drug=%d\n', numel(loom_onsets_base), numel(loom_onsets_drug));
fprintf('Flash onsets:   baseline=%d, drug=%d\n', numel(flash_onsets_base), numel(flash_onsets_drug));

%% ====================== ONSET-VS-FORMULA CROSS-CHECK ======================
fprintf('\n--- Looming onset vs formula ---\n');
check_loom_onset_formula(B.Stimuli, loom_onsets_base, 'baseline');
check_loom_onset_formula(D.Stimuli, loom_onsets_drug, 'drug');

%% ====================== TRIAL-ALIGNED EXTRACTION ======================
[Ca_loom_base, t_com_loom]   = align_trials_to_onsets(loom_onsets_base, base_dff, B.TimeCa, loom_pre, loom_post, n_com);
[Ca_loom_drug, ~]            = align_trials_to_onsets(loom_onsets_drug, drug_dff, D.TimeCa, loom_pre, loom_post, n_com);
[Ca_flash_base, t_com_flash] = align_trials_to_onsets(flash_onsets_base, base_dff, B.TimeCa, flash_pre, flash_post, n_com);
[Ca_flash_drug, ~]           = align_trials_to_onsets(flash_onsets_drug, drug_dff, D.TimeCa, flash_pre, flash_post, n_com);

% Per-cell trial-aggregated response (median across trials), matched cells only.
resp_loom_base  = median(Ca_loom_base(base_match_idx_local, :, :), 3, 'omitnan');
resp_loom_drug  = median(Ca_loom_drug(drug_match_idx_local, :, :), 3, 'omitnan');
resp_flash_base = median(Ca_flash_base(base_match_idx_local, :, :), 3, 'omitnan');
resp_flash_drug = median(Ca_flash_drug(drug_match_idx_local, :, :), 3, 'omitnan');

%% ====================== 1. RAW dF/F TRANSIENT ======================
plot_matched_comparison(t_com_loom, resp_loom_base, resp_loom_drug, ...
    sprintf('Looming - raw dF-F transient (n = %d matched cells)', n_matched), outdir);
plot_matched_comparison(t_com_flash, resp_flash_base, resp_flash_drug, ...
    sprintf('Flashes - raw dF-F transient (n = %d matched cells)', n_matched), outdir);

%% ====================== PER-CELL METRICS ======================
% Raw dF/F: mean over the peak window, in dF/F units.
loom_raw_base  = mean(resp_loom_base(:,  t_com_loom  >= loom_peak_window(1)  & t_com_loom  <= loom_peak_window(2)),  2);
loom_raw_drug  = mean(resp_loom_drug(:,  t_com_loom  >= loom_peak_window(1)  & t_com_loom  <= loom_peak_window(2)),  2);
flash_raw_base = mean(resp_flash_base(:, t_com_flash >= flash_peak_window(1) & t_com_flash <= flash_peak_window(2)), 2);
flash_raw_drug = mean(resp_flash_drug(:, t_com_flash >= flash_peak_window(1) & t_com_flash <= flash_peak_window(2)), 2);

% Baseline z-score: (peak-window mean - pre-onset baseline mean) /
% pre-onset baseline std, per trial, then median across trials.
loom_z_base  = compute_baseline_zscore(Ca_loom_base,  base_match_idx_local, t_com_loom,  loom_peak_window);
loom_z_drug  = compute_baseline_zscore(Ca_loom_drug,  drug_match_idx_local, t_com_loom,  loom_peak_window);
flash_z_base = compute_baseline_zscore(Ca_flash_base, base_match_idx_local, t_com_flash, flash_peak_window);
flash_z_drug = compute_baseline_zscore(Ca_flash_drug, drug_match_idx_local, t_com_flash, flash_peak_window);

%% ====================== 2. HISTOGRAMS (raw dF/F and z-score) ======================
plot_response_histogram(loom_raw_base, loom_raw_drug, 'Looming - raw dF-F histogram', ...
    'Mean dF/F (peak window, per cell)', n_bins, outdir);
plot_response_histogram(flash_raw_base, flash_raw_drug, 'Flashes - raw dF-F histogram', ...
    'Mean dF/F (peak window, per cell)', n_bins, outdir);
plot_response_histogram(loom_z_base, loom_z_drug, 'Looming - z-score histogram', ...
    'Peak-window z-score (vs pre-onset baseline), per cell', n_bins, outdir);
plot_response_histogram(flash_z_base, flash_z_drug, 'Flashes - z-score histogram', ...
    'Peak-window z-score (vs pre-onset baseline), per cell', n_bins, outdir);

%% ====================== 3. VIOLINS (raw, raw positive-only, z-score) ======================
plot_violin_comparison(loom_raw_base, loom_raw_drug, 'Looming - raw dF-F violin', ...
    'Mean dF/F (peak window, per cell)', n_matched, outdir);
plot_violin_comparison(flash_raw_base, flash_raw_drug, 'Flashes - raw dF-F violin', ...
    'Mean dF/F (peak window, per cell)', n_matched, outdir);

loom_pos_mask  = loom_raw_base  >= 0 & loom_raw_drug  >= 0;
flash_pos_mask = flash_raw_base >= 0 & flash_raw_drug >= 0;
plot_violin_comparison(loom_raw_base(loom_pos_mask), loom_raw_drug(loom_pos_mask), ...
    'Looming - raw dF-F violin, negative cells removed', 'Mean dF/F (peak window, per cell)', sum(loom_pos_mask), outdir);
plot_violin_comparison(flash_raw_base(flash_pos_mask), flash_raw_drug(flash_pos_mask), ...
    'Flashes - raw dF-F violin, negative cells removed', 'Mean dF/F (peak window, per cell)', sum(flash_pos_mask), outdir);

plot_violin_comparison(loom_z_base, loom_z_drug, 'Looming - z-score violin', ...
    'Peak-window z-score (vs pre-onset baseline), per cell', n_matched, outdir);
plot_violin_comparison(flash_z_base, flash_z_drug, 'Flashes - z-score violin', ...
    'Peak-window z-score (vs pre-onset baseline), per cell', n_matched, outdir);

%% ====================== 4. Z-SCORED HEATMAPS ======================
loomZ_base_trace  = compute_baseline_zscore_trace(Ca_loom_base,  base_match_idx_local, t_com_loom);
loomZ_drug_trace  = compute_baseline_zscore_trace(Ca_loom_drug,  drug_match_idx_local, t_com_loom);
flashZ_base_trace = compute_baseline_zscore_trace(Ca_flash_base, base_match_idx_local, t_com_flash);
flashZ_drug_trace = compute_baseline_zscore_trace(Ca_flash_drug, drug_match_idx_local, t_com_flash);

plot_zscore_heatmap(t_com_loom, loomZ_base_trace, loomZ_drug_trace, loom_z_base, 'Looming - z-scored heatmap', outdir);
plot_zscore_heatmap(t_com_flash, flashZ_base_trace, flashZ_drug_trace, flash_z_base, 'Flashes - z-scored heatmap', outdir);

fprintf('\nAll figures saved to: %s\n', outdir);

%% ============= HELPER FUNCTIONS =============

function name = recording_name_from_path(filepath)
    [~, name, ~] = fileparts(filepath);
    name = regexprep(name, '_preprocessed$', '');
end

function [dff_plane, roi_idx] = subset_plane(S, plane_idx)
    centroid = S.CaData(1).Ca_centroid_voxel;
    centroidZ = centroid(:, 3);
    unique_planes = unique(centroidZ);
    roi_idx = find(centroidZ == unique_planes(plane_idx));
    dff_plane = S.CaData(1).Ca_dFF(roi_idx, :);
end

function [base_match_idx_local, drug_match_idx_local, n_matched] = match_rois(roi_match_file, base_roi_idx, drug_roi_idx)
    M = load(roi_match_file);
    base_match_idx_global = M.roiMatchData.allSessionMapping(:, 1);
    drug_match_idx_global = M.roiMatchData.allSessionMapping(:, 2);
    [~, base_match_idx_local] = ismember(base_match_idx_global, base_roi_idx);
    [~, drug_match_idx_local] = ismember(drug_match_idx_global, drug_roi_idx);
    valid_matches = (base_match_idx_local > 0) & (drug_match_idx_local > 0);
    base_match_idx_local = base_match_idx_local(valid_matches);
    drug_match_idx_local = drug_match_idx_local(valid_matches);
    n_matched = numel(base_match_idx_local);
end

function loom_onsets = get_looming_onsets(Stimuli, Triggers)
    % One red sync frame per trial, fired at the true onset of that
    % trial's loom (looming has no continuous periodic red-sync rhythm,
    % unlike the other stimulus types -- see project memory).
    loom_idx = find(strcmp({Stimuli.type}, 'looming'), 1);
    s = Stimuli(loom_idx);
    if ~s.RedFrameSynchronized
        warning('Looming: RedFrameSynchronized=false -- timing fell back to PC clock, expect an offset');
    end
    % Inclusive bounds: TimeStimulusFrame(1) lands exactly on the trial-1
    % pulse, so a strict > would silently drop it.
    PT = Triggers.TimeProjector;
    loom_onsets = PT(PT >= s.TimeStimulusFrame(1) & PT <= s.TimeStimulusFrame(end))';
end

function flash_onsets = get_flash_onsets(Stimuli)
    % Flash onset times are not recorded directly -- reconstructed by
    % replaying the same RNG draws used at presentation time
    % (reconstructFlashes.m), then anchored directly to TimeStimulusFrame(1),
    % which is now camera-clock-correct on its own (RedFrameSynchronized=true).
    flash_idx = find(strcmp({Stimuli.type}, 'sparse_local_global_flashes'), 1);
    s = Stimuli(flash_idx);
    if ~s.RedFrameSynchronized
        warning('Flashes: RedFrameSynchronized=false -- timing fell back to PC clock, expect an offset');
    end
    flash_onsets_rel = reconstructFlashes(s);
    flash_onsets = s.TimeStimulusFrame(1) + flash_onsets_rel;
end

function [Ca_trials, t_com] = align_trials_to_onsets(onsets, dff, TimeCa, pre_window_sec, post_window_sec, n_com)
    t_com     = linspace(-pre_window_sec, post_window_sec, n_com);
    n_rois    = size(dff, 1);
    n_trials  = numel(onsets);
    Ca_trials = nan(n_rois, n_com, n_trials);
    time_vec  = TimeCa(1, :);

    for k = 1:n_trials
        idx = find(time_vec > onsets(k) - pre_window_sec & time_vec < onsets(k) + post_window_sec);
        if numel(idx) < 2; continue; end
        t_rel = time_vec(idx) - onsets(k);
        Ca_trials(:, :, k) = interp1(t_rel(:), dff(:, idx)', t_com(:), 'linear', 'extrap')';
    end
end

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

function Z = compute_baseline_zscore_trace(Ca_trials, match_idx, t_com)
    % Same per-trial baseline normalization as compute_baseline_zscore,
    % but keeps the full time course instead of collapsing to one window
    % -- used for the heatmaps.
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

function plot_matched_comparison(t_com, resp_base, resp_drug, fig_title, outdir)
    % resp_base/resp_drug: [n_cells x n_com], already trial-aggregated per cell
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
    title(fig_title, 'FontSize', 12, 'Interpreter', 'none');
    legend('Location', 'best');
    set(gca, 'Box', 'off'); xlim([t_com(1), t_com(end)]);
    save_current_fig(fig_title, outdir);
end

function plot_violin_comparison(vals_base, vals_drug, fig_title, ylab, n_cells, outdir)
    % One violin per condition (baseline at x=0, drug at x=1). Muted fill
    % + white median bar + small low-alpha dots, plus a paired (signrank)
    % significance bracket.
    color_base = [0.55 0.75 0.70];   % muted teal
    color_drug = [0.70 0.60 0.78];   % muted purple

    valid = ~isnan(vals_base) & ~isnan(vals_drug);
    vals_base = vals_base(valid);
    vals_drug = vals_drug(valid);

    figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [100 100 600 650]);
    hold on;

    plot_one_violin(0, vals_base, color_base);
    plot_one_violin(1, vals_drug, color_drug);

    if numel(vals_base) >= 2
        p_signrank = signrank(vals_base, vals_drug);
        y_top = max([vals_base; vals_drug]);
        y_span = range([vals_base; vals_drug]);
        y_bracket = y_top + 0.08 * y_span;
        plot([0 0 1 1], [y_bracket-0.02*y_span, y_bracket, y_bracket, y_bracket-0.02*y_span], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
        text(0.5, y_bracket + 0.03*y_span, sig_stars(p_signrank), 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
        fprintf('  %s: paired signrank p = %.4g\n', fig_title, p_signrank);
    end

    set(gca, 'XTick', [0 1], 'XTickLabel', {'Baseline', 'Drug'}, 'Box', 'off', 'XLim', [-0.6 1.6]);
    ylabel(ylab, 'FontSize', 11);
    title(sprintf('%s (n = %d)', fig_title, n_cells), 'FontSize', 12, 'Interpreter', 'none');
    save_current_fig(fig_title, outdir);
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

function plot_response_histogram(vals_base, vals_drug, fig_title, xlab, n_bins, outdir)
    % Normalized-overlay histogram of the per-cell metric.
    vals_base = vals_base(~isnan(vals_base));
    vals_drug = vals_drug(~isnan(vals_drug));
    if isempty(vals_base) && isempty(vals_drug); return; end

    edges = linspace(min([vals_base; vals_drug]), max([vals_base; vals_drug]), n_bins+1);

    figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [100 100 700 450]);
    hold on;
    histogram(vals_base, edges, 'Normalization', 'probability', ...
        'FaceColor', [0.2 0.6 1], 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', sprintf('Baseline (n=%d)', numel(vals_base)));
    histogram(vals_drug, edges, 'Normalization', 'probability', ...
        'FaceColor', [1 0.5 0.2], 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', sprintf('Drug (n=%d)', numel(vals_drug)));
    xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

    xlabel(xlab, 'FontSize', 11);
    ylabel('Probability', 'FontSize', 11);
    title(fig_title, 'FontSize', 12, 'Interpreter', 'none');
    legend('Location', 'best');
    set(gca, 'Box', 'off');
    save_current_fig(fig_title, outdir);
end

function plot_zscore_heatmap(t_com, Z_base, Z_drug, sort_vals, fig_title, outdir)
    % Z_base/Z_drug: [n_cells x n_com], per-cell z-scored trial-aligned
    % traces, matched cells, same row order in both panels. Rows sorted
    % by sort_vals (baseline peak z-score), descending, so the strongest
    % baseline responders sit at the top of BOTH panels -- makes it easy
    % to see whether the same cells are still the strong responders under
    % drug, or whether the ranking shifts.
    [~, order] = sort(sort_vals, 'descend', 'MissingPlacement', 'last');
    Zb = Z_base(order, :);
    Zd = Z_drug(order, :);

    clim = prctile([Zb(:); Zd(:)], [1 99]);
    clim = max(abs(clim)) * [-1 1];   % symmetric around 0

    figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [100 100 1100 650]);
    colormap(diverging_colormap());

    subplot(1,2,1);
    imagesc(t_com, 1:size(Zb,1), Zb, clim);
    hold on; xline(0, 'k--', 'LineWidth', 1); hold off;
    xlabel('Time from onset (s)'); ylabel('Matched cell # (sorted by baseline peak z-score)');
    title('Baseline');

    subplot(1,2,2);
    imagesc(t_com, 1:size(Zd,1), Zd, clim);
    hold on; xline(0, 'k--', 'LineWidth', 1); hold off;
    xlabel('Time from onset (s)');
    title('Drug');
    cb = colorbar; cb.Label.String = 'z-score (vs pre-onset baseline)';

    sgtitle(sprintf('%s (n = %d matched cells)', fig_title, size(Zb,1)), 'FontWeight', 'bold');
    save_current_fig(fig_title, outdir);
end

function cmap = diverging_colormap(n)
    if nargin < 1; n = 256; end
    half = floor(n/2);
    neg = [linspace(0.05,1,half)', linspace(0.05,1,half)', ones(half,1)];   % blue -> white
    pos = [ones(n-half,1), linspace(1,0.05,n-half)', linspace(1,0.05,n-half)']; % white -> red
    cmap = [neg; pos];
end

function check_loom_onset_formula(Stimuli, onsets, label)
    % How far is each individual red-frame onset from the naive
    % stimulus_trial_t*(k-1) formula prediction? Expected to be large and
    % irregular for looming (that's why we use red frames at all here).
    loom_idx   = find(strcmp({Stimuli.type}, 'looming'), 1);
    trial_t    = Stimuli(loom_idx).stimulus_trial_t;
    time_start = Stimuli(loom_idx).TimeStimulusFrame(1);
    n_trials   = numel(onsets);

    formula_onsets = time_start + (0:n_trials-1)' * trial_t;
    onsets_col     = onsets(:);
    delta          = onsets_col - formula_onsets;

    fprintf('  %s:\n', label);
    disp(table((1:n_trials)', onsets_col, formula_onsets, delta, ...
        'VariableNames', {'trial', 'detected_onset_s', 'formula_onset_s', 'delta_s'}));
end

function save_current_fig(fig_title, outdir)
    safe_name = regexprep(fig_title, '[^a-zA-Z0-9]+', '_');
    outfile = fullfile(outdir, [safe_name '.png']);
    saveas(gcf, outfile);
    fprintf('Saved %s\n', outfile);
end
