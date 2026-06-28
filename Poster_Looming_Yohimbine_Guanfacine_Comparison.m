%% Poster figure: pooled baseline vs yohimbine vs guanfacine looming response
%
% One-off poster-prep script (not part of the 6-stimulus production
% pipeline in Plot_Stimulus_Results.m) -- same category as
% Test_MovingBar_Chunking_Walkthrough.m.
%
% Two recording pairs, same animal (AnimalRB19), opposite noradrenergic
% drugs:
%   Pair A: baseline 1115 vs drug 1243 = yohimbine (raises NA)  -> looming response DECREASES
%   Pair B: baseline 1327 vs drug 1409 = guanfacine (lowers NA) -> looming response INCREASES
% (confirmed against both pairs' own per-pair Looming_transient.png: pair A's
% drug peak ~0.17 sits below its own baseline ~0.33; pair B's drug peak
% ~0.52 sits above its own baseline ~0.33 -- consistent with the pooled
% comparison below, not an artifact of pooling.)
%
% Both pairs' looming stimulus_trial_t = 7.0s, so extract_looming_trials.m
% produces an IDENTICAL t_com = linspace(-2, 5, 100) for both -- their
% resp_base/resp_drug arrays share the same column structure and can be
% pooled directly with no interpolation. The two baseline sessions are
% pooled at the matched-cell level (not averaged as two pre-computed
% curves) into one larger reference population, since they're the same
% animal's untreated state on two different days.

clear; clc; close all;

addpath(fullfile(fileparts(mfilename('fullpath')), 'analysis_functions'));

%% ====================== CONFIGURATION ======================
baseline_file_A = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1115_preprocessed.mat";
drug_file_A      = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1243_preprocessed.mat";
roi_match_file_A = "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane0_0312_V1.mat";

baseline_file_B = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260305_1327_preprocessed.mat";
drug_file_B      = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260305_1409_preprocessed.mat";
roi_match_file_B = "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane0_V2.mat";

selected_plane_idx = 1;

poster_outdir = fullfile(fileparts(mfilename('fullpath')), 'analysis_figures', '_poster_figures');
if ~exist(poster_outdir, 'dir'); mkdir(poster_outdir); end

%% ====================== PAIR A: yohimbine ======================
fprintf('\n========== LOADING PAIR A (yohimbine) ==========\n');
pairA = load_pair(baseline_file_A, drug_file_A, roi_match_file_A, selected_plane_idx);
result_A = analyze_looming_stimulus(pairA.Stimuli_base, pairA.dff_base, pairA.TimeCa_base, pairA.Triggers_base, ...
    pairA.Stimuli_drug, pairA.dff_drug, pairA.TimeCa_drug, pairA.Triggers_drug, ...
    pairA.match_idx_local_base, pairA.match_idx_local_drug, pairA.matched_rois_available);
resp_base_A = result_A.resp_base;
resp_drug_A = result_A.resp_drug;   % yohimbine
t_com = result_A.t_com;
fprintf('Pair A: %d matched cells\n', size(resp_base_A, 1));

%% ====================== PAIR B: guanfacine ======================
fprintf('\n========== LOADING PAIR B (guanfacine) ==========\n');
pairB = load_pair(baseline_file_B, drug_file_B, roi_match_file_B, selected_plane_idx);
result_B = analyze_looming_stimulus(pairB.Stimuli_base, pairB.dff_base, pairB.TimeCa_base, pairB.Triggers_base, ...
    pairB.Stimuli_drug, pairB.dff_drug, pairB.TimeCa_drug, pairB.Triggers_drug, ...
    pairB.match_idx_local_base, pairB.match_idx_local_drug, pairB.matched_rois_available);
resp_base_B = result_B.resp_base;
resp_drug_B = result_B.resp_drug;   % guanfacine
fprintf('Pair B: %d matched cells\n', size(resp_base_B, 1));

assert(isequal(result_A.t_com, result_B.t_com), 'time axes differ between pairs -- pooling invalid');

%% ====================== POOL BASELINES ======================
resp_base_pooled = [resp_base_A; resp_base_B];
n_pooled = size(resp_base_pooled, 1);
n_yohimbine = size(resp_drug_A, 1);
n_guanfacine = size(resp_drug_B, 1);
fprintf('\nn_pooled=%d, n_yohimbine=%d, n_guanfacine=%d\n', n_pooled, n_yohimbine, n_guanfacine);

%% ====================== PLOT ======================
plot_transient_3cond(t_com, resp_base_pooled, resp_drug_A, resp_drug_B, ...
    'Looming_Yohimbine_Guanfacine_Comparison', poster_outdir, n_pooled, n_yohimbine, n_guanfacine);

%% ====================== LOCAL FUNCTIONS ======================
function pair = load_pair(baseline_file, drug_file, roi_match_file, selected_plane_idx)
    B = load(baseline_file);
    base_dff_full = B.CaData(1).Ca_dFF;
    centroidZ = B.CaData(1).Ca_centroid_voxel(:, 3);
    unique_planes = unique(centroidZ);
    selected_roi_idx = find(centroidZ == unique_planes(selected_plane_idx));
    base_dff = base_dff_full(selected_roi_idx, :);

    D = load(drug_file);
    drug_dff_full = D.CaData(1).Ca_dFF;
    drug_centroidZ = D.CaData(1).Ca_centroid_voxel(:, 3);
    unique_planes_drug = unique(drug_centroidZ);
    drug_selected_roi_idx = find(drug_centroidZ == unique_planes_drug(selected_plane_idx));
    drug_dff = drug_dff_full(drug_selected_roi_idx, :);

    M = load(roi_match_file);
    base_match_idx_global = M.roiMatchData.allSessionMapping(:, 1);
    drug_match_idx_global = M.roiMatchData.allSessionMapping(:, 2);
    [~, base_match_idx_local] = ismember(base_match_idx_global, selected_roi_idx);
    [~, drug_match_idx_local] = ismember(drug_match_idx_global, drug_selected_roi_idx);
    valid_matches = (base_match_idx_local > 0) & (drug_match_idx_local > 0);
    base_match_idx_local = base_match_idx_local(valid_matches);
    drug_match_idx_local = drug_match_idx_local(valid_matches);
    matched_rois_available = ~isempty(base_match_idx_local);

    pair.Stimuli_base = B.Stimuli;
    pair.Stimuli_drug = D.Stimuli;
    pair.dff_base = base_dff;
    pair.dff_drug = drug_dff;
    pair.TimeCa_base = B.TimeCa;
    pair.TimeCa_drug = D.TimeCa;
    pair.Triggers_base = B.Triggers;
    pair.Triggers_drug = D.Triggers;
    pair.match_idx_local_base = base_match_idx_local;
    pair.match_idx_local_drug = drug_match_idx_local;
    pair.matched_rois_available = matched_rois_available;
end

function plot_transient_3cond(t_com, resp_base, resp_yohimbine, resp_guanfacine, fig_title, outdir, n_pooled, n_yohimbine, n_guanfacine)
    % mean/std/sem computed ACROSS CELLS (dim=1: one row per matched cell,
    % resp_base is already the pooled [1562 x 100] array) at each of the
    % 100 timepoints -- standard population mean+-SEM convention, same
    % formula as Plot_Stimulus_Results.m's plot_transient. SEM = std /
    % sqrt(n_cells), so a larger pooled n (1562 vs 913/649) mechanically
    % produces a tighter band -- not a separate effect, just sqrt(n) at
    % work. Diagnostic print below shows the actual numbers at each
    % condition's own peak so this can be checked directly rather than
    % just trusted.
    mean_base = mean(resp_base, 1, 'omitnan');
    sem_base  = std(resp_base, 0, 1, 'omitnan') / sqrt(size(resp_base, 1));
    mean_yoh = mean(resp_yohimbine, 1, 'omitnan');
    sem_yoh  = std(resp_yohimbine, 0, 1, 'omitnan') / sqrt(size(resp_yohimbine, 1));
    mean_gua = mean(resp_guanfacine, 1, 'omitnan');
    sem_gua  = std(resp_guanfacine, 0, 1, 'omitnan') / sqrt(size(resp_guanfacine, 1));

    fprintf('\n  Peak diagnostics (mean +/- SEM at each condition''s own peak timepoint):\n');
    print_peak_stats('Baseline (pooled)', t_com, mean_base, sem_base, resp_base);
    print_peak_stats('Yohimbine', t_com, mean_yoh, sem_yoh, resp_yohimbine);
    print_peak_stats('Guanfacine', t_com, mean_gua, sem_gua, resp_guanfacine);

    color_base = [0.3 0.3 0.3];     % grey
    color_yoh  = [0.20 0.45 0.80];  % cool blue -- yohimbine DECREASES the loom response
    color_gua  = [0.85 0.25 0.15];  % warm red -- guanfacine INCREASES the loom response

    figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [100 100 900 500]);
    hold on;
    x_sh = [t_com, fliplr(t_com)];

    fill(x_sh, [(mean_base+sem_base), fliplr(mean_base-sem_base)], color_base, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t_com, mean_base, '-', 'Color', color_base, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Baseline (n=%d cells)', n_pooled));

    fill(x_sh, [(mean_yoh+sem_yoh), fliplr(mean_yoh-sem_yoh)], color_yoh, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t_com, mean_yoh, '-', 'Color', color_yoh, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Yohimbine (n=%d cells)', n_yohimbine));

    fill(x_sh, [(mean_gua+sem_gua), fliplr(mean_gua-sem_gua)], color_gua, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t_com, mean_gua, '-', 'Color', color_gua, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Guanfacine (n=%d cells)', n_guanfacine));

    xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlabel('Time from onset (s)', 'FontSize', 11);
    ylabel('dF/F (mean \pm SEM)', 'FontSize', 11);
    title(sprintf('Looming response: noradrenergic modulation (baseline n=%d, yohimbine n=%d, guanfacine n=%d)', ...
        n_pooled, n_yohimbine, n_guanfacine), 'FontSize', 12);
    legend('Location', 'best'); set(gca, 'Box', 'off'); xlim([t_com(1), t_com(end)]);
    save_fig_local(fig_title, outdir);
end

function print_peak_stats(label, t_com, mean_trace, sem_trace, resp)
    [pk_val, pk_idx] = max(mean_trace);
    n_cells = size(resp, 1);
    raw_std = std(resp(:, pk_idx), 0, 1, 'omitnan');
    fprintf('    %-20s n=%-5d  t=%.2fs  mean=%.4f  std(across cells)=%.4f  SEM=std/sqrt(n)=%.4f\n', ...
        label, n_cells, t_com(pk_idx), pk_val, raw_std, sem_trace(pk_idx));
end

function save_fig_local(fig_title, outdir)
    safe_name = regexprep(fig_title, '[^a-zA-Z0-9]+', '_');
    outfile = fullfile(outdir, [safe_name '.png']);
    saveas(gcf, outfile);
    fprintf('  Saved %s\n', outfile);
end
