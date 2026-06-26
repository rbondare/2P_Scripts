%% Loom_Flashes_ZETA_Analysis.m
% Standalone analysis for looming and sparse_local_global_flashes:
%   1. Matched-cell triggered average, baseline vs drug
%   2. Violin plots of each cell's post-onset response (all matched
%      cells, then restricted to ZETA-responsive cells only)
%   3. Onset-vs-formula cross-checks (does the red-frame/reconstructed
%      onset for each individual trial match what the naive
%      stimulus_trial_t formula would have predicted)
%   4. ZETA p-value distribution sanity check (should look broadly
%      uniform/flat, with an excess near 0 only for real responders --
%      a skewed or spiky distribution would flag a problem with the
%      test setup itself, not just "few responsive cells")
%   5. Responsive-vs-non-responsive triggered average, per condition --
%      the direct visual check that ZETA's call matches what the actual
%      traces look like
%
% Kept separate from Master_Stimulus_Analysis.m until validated -- will
% be merged in once confirmed. See project memory
% "project-stimulus-timing-architecture" for the full architecture
% writeup and open questions this script doesn't resolve.

addpath(genpath(fullfile(pwd, 'zetatest')));
addpath(genpath(fullfile(pwd, 'violin_plot_utils')));

%% ====================== CONFIGURATION ======================
baseline_file  = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1115_preprocessed.mat";
drug_file      = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1243_preprocessed.mat";
roi_match_file = "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane0_0312_V1.mat";
selected_plane_idx = 1;

zeta_alpha      = 0.05;
zeta_resamp_num = 100;     % toolbox default is 250; reduced for runtime, see note in script

loom_pre  = 2; loom_post  = 5;    % trial-aligned window, looming (seconds around onset)
flash_pre = 1; flash_post = 4;    % trial-aligned window, flashes
n_com     = 100;                 % points in the common trial-aligned time grid

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
loom_onsets_base  = get_looming_onsets(B.Stimuli, B.Triggers);
loom_onsets_drug  = get_looming_onsets(D.Stimuli, D.Triggers);
[flash_onsets_base, flash_anchor_base, flash_anchor_dt_base, flash_onsets_rel_base] = get_flash_onsets(B.Stimuli, B.Triggers);
[flash_onsets_drug, flash_anchor_drug, flash_anchor_dt_drug, flash_onsets_rel_drug] = get_flash_onsets(D.Stimuli, D.Triggers);

fprintf('Looming onsets: baseline=%d, drug=%d\n', numel(loom_onsets_base), numel(loom_onsets_drug));
fprintf('Flash onsets:   baseline=%d, drug=%d\n', numel(flash_onsets_base), numel(flash_onsets_drug));

%% ====================== ONSET-VS-FORMULA CROSS-CHECK ======================
% Looming: how far is each individual red-frame onset from the naive
% stimulus_trial_t*(k-1) formula? (Should be irregular/large -- that's
% the whole reason we use red frames instead of the formula for looming.)
fprintf('\n--- Looming onset vs formula ---\n');
check_loom_onset_formula(B.Stimuli, loom_onsets_base, 'baseline');
check_loom_onset_formula(D.Stimuli, loom_onsets_drug, 'drug');

% Flashes: how much did the red-frame anchor correction shift things,
% relative to just trusting TimeStimulusFrame(1) directly?
fprintf('\n--- Flashes anchor correction ---\n');
check_flash_anchor_correction(B.Stimuli, flash_anchor_base, flash_anchor_dt_base, flash_onsets_rel_base, 'baseline');
check_flash_anchor_correction(D.Stimuli, flash_anchor_drug, flash_anchor_dt_drug, flash_onsets_rel_drug, 'drug');

%% ====================== TRIAL-ALIGNED EXTRACTION ======================
[Ca_loom_base, t_com_loom]   = align_trials_to_onsets(loom_onsets_base, base_dff, B.TimeCa, loom_pre, loom_post, n_com);
[Ca_loom_drug, ~]            = align_trials_to_onsets(loom_onsets_drug, drug_dff, D.TimeCa, loom_pre, loom_post, n_com);
[Ca_flash_base, t_com_flash] = align_trials_to_onsets(flash_onsets_base, base_dff, B.TimeCa, flash_pre, flash_post, n_com);
[Ca_flash_drug, ~]           = align_trials_to_onsets(flash_onsets_drug, drug_dff, D.TimeCa, flash_pre, flash_post, n_com);

% Per-cell trial-aggregated response (median across trials -- robust to
% the kind of single-outlier-trial issue we found earlier)
resp_loom_base  = median(Ca_loom_base(base_match_idx_local, :, :), 3, 'omitnan');
resp_loom_drug  = median(Ca_loom_drug(drug_match_idx_local, :, :), 3, 'omitnan');
resp_flash_base = median(Ca_flash_base(base_match_idx_local, :, :), 3, 'omitnan');
resp_flash_drug = median(Ca_flash_drug(drug_match_idx_local, :, :), 3, 'omitnan');

%% ====================== PLOT: TRIGGERED AVERAGE, MATCHED CELLS, BASELINE VS DRUG ======================
plot_matched_comparison(t_com_loom, resp_loom_base, resp_loom_drug, ...
    sprintf('Looming: baseline vs drug (n = %d matched cells)', n_matched));
plot_matched_comparison(t_com_flash, resp_flash_base, resp_flash_drug, ...
    sprintf('Flashes: baseline vs drug (n = %d matched cells)', n_matched));

%% ====================== ZETA RESPONSIVENESS (matched cells only) ======================
fprintf('\nRunning ZETA: looming, baseline (%d cells)...\n', n_matched);
loom_p_base  = run_zeta_population(base_dff(base_match_idx_local, :), B.TimeCa(1,:)', loom_onsets_base,  loom_post, zeta_resamp_num);
fprintf('Running ZETA: looming, drug...\n');
loom_p_drug  = run_zeta_population(drug_dff(drug_match_idx_local, :), D.TimeCa(1,:)', loom_onsets_drug,  loom_post, zeta_resamp_num);
fprintf('Running ZETA: flashes, baseline...\n');
flash_p_base = run_zeta_population(base_dff(base_match_idx_local, :), B.TimeCa(1,:)', flash_onsets_base, flash_post, zeta_resamp_num);
fprintf('Running ZETA: flashes, drug...\n');
flash_p_drug = run_zeta_population(drug_dff(drug_match_idx_local, :), D.TimeCa(1,:)', flash_onsets_drug, flash_post, zeta_resamp_num);

% Responsive = significant in EITHER condition (so drug-induced silencing
% of a baseline-responsive cell doesn't drop it from the comparison) --
% used for the violin plots below.
loom_responsive  = (loom_p_base < zeta_alpha) | (loom_p_drug < zeta_alpha);
flash_responsive = (flash_p_base < zeta_alpha) | (flash_p_drug < zeta_alpha);

% Condition-specific masks -- used for the responsive-vs-non-responsive
% triggered-average check below, where we want to know "did THIS
% condition's own ZETA call match THIS condition's own traces" rather
% than mixing in the other condition's responsiveness.
loom_responsive_base  = loom_p_base  < zeta_alpha;
loom_responsive_drug  = loom_p_drug  < zeta_alpha;
flash_responsive_base = flash_p_base < zeta_alpha;
flash_responsive_drug = flash_p_drug < zeta_alpha;

fprintf('\nLooming ZETA-responsive:  %d / %d matched cells (either condition)\n', sum(loom_responsive), n_matched);
fprintf('  baseline only: %d, drug only: %d\n', sum(loom_responsive_base), sum(loom_responsive_drug));
fprintf('Flashes ZETA-responsive:  %d / %d matched cells (either condition)\n', sum(flash_responsive), n_matched);
fprintf('  baseline only: %d, drug only: %d\n', sum(flash_responsive_base), sum(flash_responsive_drug));

%% ====================== ZETA P-VALUE SANITY CHECK ======================
% Under the null (no real response), p-values should be roughly uniform
% on [0,1] -- a flat histogram. A real responsive subpopulation shows up
% as an excess near 0 on top of that flat background. A histogram that's
% skewed/spiky everywhere (not just near 0) would indicate a problem with
% the test setup itself (e.g. dblUseMaxDur or resampling count badly
% mismatched to the data), not just "few responsive cells".
plot_zeta_pvalue_check(loom_p_base, loom_p_drug, 'Looming: ZETA p-value distribution');
plot_zeta_pvalue_check(flash_p_base, flash_p_drug, 'Flashes: ZETA p-value distribution');

%% ====================== PER-CELL SCALAR METRIC FOR VIOLINS ======================
% Mean dF/F over the post-onset window, per cell, from the
% already-trial-aggregated (median) response
loom_val_base  = mean(resp_loom_base(:,  t_com_loom  >= 0), 2);
loom_val_drug  = mean(resp_loom_drug(:,  t_com_loom  >= 0), 2);
flash_val_base = mean(resp_flash_base(:, t_com_flash >= 0), 2);
flash_val_drug = mean(resp_flash_drug(:, t_com_flash >= 0), 2);

%% ====================== VIOLIN PLOTS ======================
plot_violin_comparison(loom_val_base, loom_val_drug, 'Looming -- all matched cells', n_matched);
plot_violin_comparison(loom_val_base(loom_responsive), loom_val_drug(loom_responsive), ...
    'Looming -- ZETA-responsive cells only', sum(loom_responsive));

plot_violin_comparison(flash_val_base, flash_val_drug, 'Flashes -- all matched cells', n_matched);
plot_violin_comparison(flash_val_base(flash_responsive), flash_val_drug(flash_responsive), ...
    'Flashes -- ZETA-responsive cells only', sum(flash_responsive));

%% ====================== RESPONSIVE VS NON-RESPONSIVE TRIGGERED AVERAGE ======================
% The direct validity check: do the cells ZETA calls significant in THIS
% condition actually show a visible bump in THIS condition's own
% triggered average, and do the non-responsive cells look flat by
% comparison? Uses each condition's own ZETA call (not the "either
% condition" mask used for the violins above).
plot_responsive_vs_nonresponsive(t_com_loom, resp_loom_base, loom_responsive_base, ...
    'Looming, baseline: ZETA-responsive vs non-responsive');
plot_responsive_vs_nonresponsive(t_com_loom, resp_loom_drug, loom_responsive_drug, ...
    'Looming, drug: ZETA-responsive vs non-responsive');
plot_responsive_vs_nonresponsive(t_com_flash, resp_flash_base, flash_responsive_base, ...
    'Flashes, baseline: ZETA-responsive vs non-responsive');
plot_responsive_vs_nonresponsive(t_com_flash, resp_flash_drug, flash_responsive_drug, ...
    'Flashes, drug: ZETA-responsive vs non-responsive');

%% ============= HELPER FUNCTIONS =============

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
    loom_idx   = find(strcmp({Stimuli.type}, 'looming'), 1);
    time_start = Stimuli(loom_idx).TimeStimulusFrame(1);
    time_end   = Stimuli(loom_idx).TimeStimulusFrame(end);
    PT = Triggers.TimeProjector;
    loom_onsets = PT(PT > time_start & PT < time_end)';
end

function [flash_onsets, anchor, anchor_dt, flash_onsets_rel] = get_flash_onsets(Stimuli, Triggers)
    % Flash onset times are not recorded directly -- reconstructed by
    % replaying the same RNG draws used at presentation time
    % (reconstructFlashes.m), then anchored to the block's real start via
    % the nearest red frame (this stimulus also uses PC-clock fallback
    % timing, same ~3.3s session-wide clock offset as the other stimuli --
    % confirmed root cause: stimulus PC and ScanImage PC are physically
    % separate, unsynced clocks; see project memory).
    flash_idx = find(strcmp({Stimuli.type}, 'sparse_local_global_flashes'), 1);
    s = Stimuli(flash_idx);
    block_start = s.TimeStimulusFrame(1);
    flash_onsets_rel = reconstructFlashes(s);

    PT = Triggers.TimeProjector;
    [anchor_dt, anchor_loc] = min(abs(PT - block_start));
    anchor = PT(anchor_loc);
    if anchor_dt > 1
        warning('Flashes: nearest red frame is %.2f s from TimeStimulusFrame(1) -- check block alignment', anchor_dt);
    end
    flash_onsets = anchor + flash_onsets_rel;
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

function vec_p = run_zeta_population(dff_matched, vecTime, onsets, dblUseMaxDur, intResampNum)
    % One zetatstest call per cell -- operates on the raw continuous
    % trace + event onset times directly (not the trial-aligned array).
    n_cells = size(dff_matched, 1);
    vec_p = nan(n_cells, 1);
    onsets = onsets(:);
    for c = 1:n_cells
        vec_p(c) = zetatstest(vecTime, dff_matched(c, :)', onsets, dblUseMaxDur, intResampNum, 0);
    end
end

function plot_matched_comparison(t_com, resp_base, resp_drug, fig_title)
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
end

function plot_violin_comparison(vals_base, vals_drug, fig_title, n_cells)
    % One violin per condition (baseline at x=0, drug at x=1), each dot
    % is one matched cell's post-onset mean dF/F.
    vals_base = vals_base(~isnan(vals_base));
    vals_drug = vals_drug(~isnan(vals_drug));

    figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [100 100 600 600]);
    hold on;

    if numel(vals_base) >= 2
        [f_bl, xi_bl] = ksdensity(vals_base, 'NumPoints', 150);
        f_bl = f_bl / max(f_bl) * 0.4;
        patch([f_bl, fliplr(-f_bl)], [xi_bl, fliplr(xi_bl)], [0.2 0.6 1], 'FaceAlpha', 0.8, 'EdgeColor', 'black', 'LineWidth', 1.5);
        plot([-0.4 0.4], [median(vals_base) median(vals_base)], 'k-', 'LineWidth', 2.5);
        x_jitter = max(min(randn(numel(vals_base), 1)*0.12, 0.4), -0.4);
        scatter(x_jitter, vals_base, 20, [0.3 0.3 0.3], 'o', 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.2);
    end

    if numel(vals_drug) >= 2
        [f_dr, xi_dr] = ksdensity(vals_drug, 'NumPoints', 150);
        f_dr = f_dr / max(f_dr) * 0.4;
        patch([f_dr + 1, fliplr(-f_dr + 1)], [xi_dr, fliplr(xi_dr)], [1 0.5 0.2], 'FaceAlpha', 0.8, 'EdgeColor', 'black', 'LineWidth', 1.5);
        plot([1-0.4 1+0.4], [median(vals_drug) median(vals_drug)], 'k-', 'LineWidth', 2.5);
        x_jitter = 1 + max(min(randn(numel(vals_drug), 1)*0.12, 0.4), -0.4);
        scatter(x_jitter, vals_drug, 20, [0.3 0.3 0.3], 'o', 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.2);
    end

    set(gca, 'XTick', [0 1], 'XTickLabel', {'Baseline', 'Drug'}, 'Box', 'off');
    ylabel('Mean dF/F (post-onset, per cell)', 'FontSize', 11);
    title(sprintf('%s (n = %d)', fig_title, n_cells), 'FontSize', 12, 'Interpreter', 'none');
    grid on;
end

function check_loom_onset_formula(Stimuli, onsets, label)
    % How far is each individual red-frame onset from the naive
    % stimulus_trial_t*(k-1) formula prediction? Expected to be large and
    % irregular for looming (that's why we use red frames at all here) --
    % this just makes that explicit and per-trial instead of only the
    % session-wide averages already written up in project memory.
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

function check_flash_anchor_correction(Stimuli, anchor, anchor_dt, flash_onsets_rel, label)
    % How much does the per-block red-frame anchor correction shift each
    % flash onset relative to just trusting TimeStimulusFrame(1) (the
    % uncorrected PC-clock-fallback timing) directly?
    flash_idx   = find(strcmp({Stimuli.type}, 'sparse_local_global_flashes'), 1);
    block_start = Stimuli(flash_idx).TimeStimulusFrame(1);

    uncorrected_onsets = block_start + flash_onsets_rel(:);
    corrected_onsets    = anchor + flash_onsets_rel(:);

    fprintf('  %s: block_start=%.3f s, nearest red frame=%.3f s, correction=+%.3f s applied uniformly to all %d onsets\n', ...
        label, block_start, anchor, anchor_dt, numel(flash_onsets_rel));
    disp(table((1:numel(flash_onsets_rel))', uncorrected_onsets, corrected_onsets, ...
        'VariableNames', {'flash', 'uncorrected_onset_s', 'corrected_onset_s'}));
end

function plot_zeta_pvalue_check(p_base, p_drug, fig_title)
    % Sanity check on the test itself, not the biology: under the null,
    % p-values should be roughly flat/uniform on [0,1]. A real responsive
    % subpopulation shows up as an excess near 0 on top of that flat
    % background. A histogram that's skewed/spiky everywhere (not just
    % near 0) would point to a problem with dblUseMaxDur/resampling, not
    % just "few responsive cells".
    n_bins = 20;

    figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [100 100 700 450]);
    hold on;
    histogram(p_base, n_bins, 'BinLimits', [0 1], 'Normalization', 'probability', ...
        'FaceColor', [0.2 0.6 1], 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Baseline');
    histogram(p_drug, n_bins, 'BinLimits', [0 1], 'Normalization', 'probability', ...
        'FaceColor', [1 0.5 0.2], 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Drug');
    yline(1/n_bins, 'k--', 'LineWidth', 1, 'DisplayName', 'Expected flat height under null');

    xlabel('ZETA p-value', 'FontSize', 11);
    ylabel('Probability', 'FontSize', 11);
    title(fig_title, 'FontSize', 12, 'Interpreter', 'none');
    legend('Location', 'best');
    set(gca, 'Box', 'off'); xlim([0 1]);
end

function plot_responsive_vs_nonresponsive(t_com, resp, responsive_mask, fig_title)
    % resp: [n_cells x n_com], already trial-aggregated per cell.
    % responsive_mask: logical, this condition's own ZETA call.
    resp_yes = resp(responsive_mask, :);
    resp_no  = resp(~responsive_mask, :);

    mean_yes = mean(resp_yes, 1, 'omitnan');
    sem_yes  = std(resp_yes, 0, 1, 'omitnan') / sqrt(max(size(resp_yes, 1), 1));
    mean_no  = mean(resp_no, 1, 'omitnan');
    sem_no   = std(resp_no, 0, 1, 'omitnan') / sqrt(max(size(resp_no, 1), 1));

    figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [100 100 900 500]);
    hold on;
    x_sh = [t_com, fliplr(t_com)];

    y_sh_no = [(mean_no + sem_no), fliplr(mean_no - sem_no)];
    fill(x_sh, y_sh_no, [0.6 0.6 0.6], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t_com, mean_no, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 2, ...
        'DisplayName', sprintf('Non-responsive (n=%d)', sum(~responsive_mask)));

    y_sh_yes = [(mean_yes + sem_yes), fliplr(mean_yes - sem_yes)];
    fill(x_sh, y_sh_yes, [0.85 0.1 0.1], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t_com, mean_yes, '-', 'Color', [0.85 0.1 0.1], 'LineWidth', 2.5, ...
        'DisplayName', sprintf('ZETA-responsive (n=%d)', sum(responsive_mask)));

    xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlabel('Time from onset (s)', 'FontSize', 11);
    ylabel('dF/F (mean \pm SEM)', 'FontSize', 11);
    title(fig_title, 'FontSize', 12, 'Interpreter', 'none');
    legend('Location', 'best');
    set(gca, 'Box', 'off'); xlim([t_com(1), t_com(end)]);
end
