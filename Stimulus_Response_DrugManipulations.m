%% STIMULUS RESPONSE COMPARISON PIPELINE
% Compares calcium responses between two aggregated recordings.
% Clear, section-by-section workflow for exploratory analysis.
% Plot aesthetics are kept similar to your original script.
%
% Requirements:
% - preprocessed files with fields: CaData, TimeCa, Stimuli
% - Statistics and Machine Learning Toolbox for ttest2/kstest2 (optional)

%% SECTION 1
clc;
close all;
clear;

% Recording labels for plots
label1 = 'Control';
label2 = 'Guanfacine';

% Input files
file1 = 'Z:\group\joeschgrp\Group Members\Rima\Aggregated\Animal25_240716_1453_preprocessed.mat';
file2 = 'Z:\group\joeschgrp\Group Members\Rima\Aggregated\Animal25_240717_1449_preprocessed.mat';

% Stimulus types to compare. Leave empty {} to auto-detect intersection.
stimulus_types = {'moving_bar', 'grating', 'spontaneous'};

% Plot/save options
save_figures = false;
out_dir = fullfile(pwd, 'comparison_figures');
num_bins_hist = 30;
num_bins_ratio = 200;
max_cells_heatmap = 1000;
max_frames_heatmap = 3000;

if save_figures && ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% SECTION 2 - LOAD DATA
data1 = load(file1);
data2 = load(file2);

calcium_activity1 = data1.CaData(1).Ca_dFF;
calcium_activity2 = data2.CaData(1).Ca_dFF;

time_ca_1 = get_time_vector(data1.TimeCa);
time_ca_2 = get_time_vector(data2.TimeCa);

fprintf('Loaded:\n');
fprintf('  %s: %d neurons x %d frames\n', label1, size(calcium_activity1, 1), size(calcium_activity1, 2));
fprintf('  %s: %d neurons x %d frames\n', label2, size(calcium_activity2, 1), size(calcium_activity2, 2));

%% SECTION 3 - WHOLE-RECORDING DISTRIBUTION
mean_activity1_all = mean(calcium_activity1, 2);
mean_activity2_all = mean(calcium_activity2, 2);

combined_data = [mean_activity1_all; mean_activity2_all];
bin_edges = linspace(min(combined_data), max(combined_data), num_bins_hist + 1);

figure('Name', 'All Neurons Mean dF/F Distribution', 'NumberTitle', 'off');
histogram(mean_activity1_all, 'Normalization', 'pdf', ...
    'FaceColor', 'r', 'FaceAlpha', 0.4, 'EdgeColor', 'r', 'BinEdges', bin_edges);
hold on;
histogram(mean_activity2_all, 'Normalization', 'pdf', ...
    'FaceColor', 'b', 'FaceAlpha', 0.4, 'EdgeColor', 'b', 'BinEdges', bin_edges);
xlabel('dF/F');
ylabel('Probability Density');
title('Distribution of Mean Activity (All Neurons)');
legend({label1, label2}, 'Location', 'best');
hold off;

if save_figures
    saveas(gcf, fullfile(out_dir, 'hist_mean_activity_all.png'));
end

%% SECTION 4 - PREPARE STIMULUS TYPES
if isempty(stimulus_types)
    types1 = get_stimulus_types(data1.Stimuli);
    types2 = get_stimulus_types(data2.Stimuli);
    stimulus_types = intersect(types1, types2, 'stable');
end

if isempty(stimulus_types)
    error('No common stimulus types found between datasets.');
end

fprintf('\nStimulus types to analyze:\n');
for i = 1:numel(stimulus_types)
    fprintf('  - %s\n', stimulus_types{i});
end

%% SECTION 5 - PER-STIMULUS ANALYSIS LOOP
results = struct();

for s = 1:numel(stimulus_types)
    stimulus_type = stimulus_types{s};
    fprintf('\n=== Stimulus: %s ===\n', stimulus_type);

    stim1 = filter_stimuli_by_type(data1.Stimuli, stimulus_type);
    stim2 = filter_stimuli_by_type(data2.Stimuli, stimulus_type);

    if isempty(stim1) || isempty(stim2)
        warning('Skipping %s: missing in one dataset.', stimulus_type);
        continue;
    end

    concat1 = concatenate_stimulus_segments(calcium_activity1, time_ca_1, stim1);
    concat2 = concatenate_stimulus_segments(calcium_activity2, time_ca_2, stim2);

    if isempty(concat1) || isempty(concat2)
        warning('Skipping %s: could not extract stimulus segments.', stimulus_type);
        continue;
    end

    mean1 = mean(concat1, 2);
    mean2 = mean(concat2, 2);

    std1 = std(concat1, 0, 2);
    std2 = std(concat2, 0, 2);
    ratio1 = std1 ./ max(mean1, eps);
    ratio2 = std2 ./ max(mean2, eps);

    % Store results
    results.(stimulus_type).mean1 = mean1;
    results.(stimulus_type).mean2 = mean2;
    results.(stimulus_type).ratio1 = ratio1;
    results.(stimulus_type).ratio2 = ratio2;
    results.(stimulus_type).concat1 = concat1;
    results.(stimulus_type).concat2 = concat2;

    %% SECTION 5A - Histogram of mean response for this stimulus
    combined_stim = [mean1; mean2];
    bin_edges_stim = linspace(min(combined_stim), max(combined_stim), num_bins_hist + 1);

    figure('Name', ['Mean dF/F - ' stimulus_type], 'NumberTitle', 'off');
    histogram(mean1, 'Normalization', 'pdf', ...
        'FaceColor', 'r', 'FaceAlpha', 0.4, 'EdgeColor', 'r', 'BinEdges', bin_edges_stim);
    hold on;
    histogram(mean2, 'Normalization', 'pdf', ...
        'FaceColor', 'b', 'FaceAlpha', 0.4, 'EdgeColor', 'b', 'BinEdges', bin_edges_stim);
    xlabel('dF/F');
    ylabel('Probability Density');
    title(['Mean dF/F Distribution - ' strrep(stimulus_type, '_', ' ')]);
    legend({label1, label2}, 'Location', 'best');
    hold off;

    if save_figures
        saveas(gcf, fullfile(out_dir, ['hist_mean_' stimulus_type '.png']));
    end

    %% SECTION 5B - Histogram of SD/Mean variability ratio
    combined_ratio = [ratio1; ratio2];
    bin_edges_ratio = linspace(min(combined_ratio), max(combined_ratio), num_bins_ratio + 1);

    figure('Name', ['SD-Mean Ratio - ' stimulus_type], 'NumberTitle', 'off');
    histogram(ratio1, 'Normalization', 'pdf', ...
        'FaceColor', 'r', 'FaceAlpha', 0.4, 'EdgeColor', 'r', 'BinEdges', bin_edges_ratio);
    hold on;
    histogram(ratio2, 'Normalization', 'pdf', ...
        'FaceColor', 'b', 'FaceAlpha', 0.4, 'EdgeColor', 'b', 'BinEdges', bin_edges_ratio);

    if strcmp(stimulus_type, 'grating') || strcmp(stimulus_type, 'full_field_flash')
        xlim([0 200]);
    else
        xlim([0 50]);
    end

    xlabel('SD/Mean Ratio');
    ylabel('Probability Density');
    title(['SD/Mean Ratio - ' strrep(stimulus_type, '_', ' ')]);
    legend({label1, label2}, 'Location', 'best');
    hold off;

    if save_figures
        saveas(gcf, fullfile(out_dir, ['hist_ratio_' stimulus_type '.png']));
    end

    %% SECTION 5C - Heatmap for each condition (same style as original)
    plot_stimulus_heatmap(concat1, stimulus_type, label1, max_cells_heatmap, max_frames_heatmap, save_figures, out_dir);
    plot_stimulus_heatmap(concat2, stimulus_type, label2, max_cells_heatmap, max_frames_heatmap, save_figures, out_dir);

    %% SECTION 5D - Stats for this stimulus
    [h_ttest, p_ttest] = safe_ttest2(mean1, mean2);
    [h_ks, p_ks, ks_stat] = safe_kstest2(mean1, mean2);

    results.(stimulus_type).stats = struct( ...
        'ttest_h', h_ttest, 'ttest_p', p_ttest, ...
        'ks_h', h_ks, 'ks_p', p_ks, 'ks_stat', ks_stat);

    fprintf('  n neurons: %s=%d, %s=%d\n', label1, numel(mean1), label2, numel(mean2));
    fprintf('  mean(dF/F): %s=%.4f, %s=%.4f\n', label1, mean(mean1), label2, mean(mean2));
    fprintf('  t-test p=%.3g, KS p=%.3g\n', p_ttest, p_ks);
end

%% SECTION 6 - CROSS-STIMULUS SUMMARY (SWARM STYLE)
valid_types = fieldnames(results);
if isempty(valid_types)
    warning('No valid stimulus results to summarize.');
else
    figure('Name', 'Cross-Stimulus Mean Activity Summary', 'NumberTitle', 'off');
    hold on;

    x_labels = {};
    x_idx = 1;

    for i = 1:numel(valid_types)
        st = valid_types{i};
        y1 = results.(st).mean1;
        y2 = results.(st).mean2;

        swarmchart(repmat(x_idx, numel(y1), 1), y1, 20, 'r', 'filled', 'MarkerFaceAlpha', 0.35);
        x_labels{end+1} = [strrep(st, '_', ' ') ' ' label1]; %#ok<AGROW>
        x_idx = x_idx + 1;

        swarmchart(repmat(x_idx, numel(y2), 1), y2, 20, 'b', 'filled', 'MarkerFaceAlpha', 0.35);
        x_labels{end+1} = [strrep(st, '_', ' ') ' ' label2]; %#ok<AGROW>
        x_idx = x_idx + 1;
    end

    title('Mean Activity by Stimulus and Condition');
    ylabel('Mean dF/F');
    set(gca, 'XTick', 1:numel(x_labels), 'XTickLabel', x_labels);
    xtickangle(35);
    xlim([0.5 numel(x_labels) + 0.5]);
    grid on;
    hold off;

    if save_figures
        saveas(gcf, fullfile(out_dir, 'summary_swarm_mean_activity.png'));
    end
end

%% SECTION 7 - PRINT SUMMARY TABLE
if ~isempty(valid_types)
    fprintf('\n===== SUMMARY TABLE =====\n');
    fprintf('%-22s  %-10s  %-10s  %-10s  %-10s\n', 'Stimulus', ['Mean ' label1], ['Mean ' label2], 'p(t-test)', 'p(KS)');
    for i = 1:numel(valid_types)
        st = valid_types{i};
        m1 = mean(results.(st).mean1);
        m2 = mean(results.(st).mean2);
        pt = results.(st).stats.ttest_p;
        pk = results.(st).stats.ks_p;
        fprintf('%-22s  %-10.4f  %-10.4f  %-10.3g  %-10.3g\n', st, m1, m2, pt, pk);
    end
end


%% ===================== LOCAL HELPER FUNCTIONS =====================
function t = get_time_vector(TimeCa)
% Handle either 1xN or 2xN TimeCa formats robustly.
    if isvector(TimeCa)
        t = TimeCa(:)';
    elseif size(TimeCa, 1) >= 2
        t = TimeCa(2, :);
    else
        t = TimeCa(1, :);
    end
end

function types = get_stimulus_types(Stimuli)
    if isempty(Stimuli)
        types = {};
        return;
    end
    types = {Stimuli.type};
    types = unique(types, 'stable');
end

function filtered = filter_stimuli_by_type(Stimuli, stimulus_type)
    if isempty(Stimuli)
        filtered = [];
        return;
    end
    idx = arrayfun(@(s) isfield(s, 'type') && strcmp(s.type, stimulus_type), Stimuli);
    filtered = Stimuli(idx);
end

function concat_data = concatenate_stimulus_segments(calcium_activity, time_ca, filtered_stimuli)
    concat_data = [];
    for i = 1:numel(filtered_stimuli)
        if ~isfield(filtered_stimuli(i), 'TimeStimulusFrame') || isempty(filtered_stimuli(i).TimeStimulusFrame)
            continue;
        end

        start_time = filtered_stimuli(i).TimeStimulusFrame(1);
        end_time = filtered_stimuli(i).TimeStimulusFrame(end);

        [~, start_idx] = min(abs(time_ca - start_time));
        [~, end_idx] = min(abs(time_ca - end_time));

        if end_idx < start_idx
            tmp = start_idx;
            start_idx = end_idx;
            end_idx = tmp;
        end

        segment = calcium_activity(:, start_idx:end_idx);
        concat_data = [concat_data, segment]; %#ok<AGROW>
    end
end

function plot_stimulus_heatmap(concat_data, stimulus_type, label_name, max_cells, max_frames, save_figures, out_dir)
    num_cells = size(concat_data, 1);
    num_frames = size(concat_data, 2);

    num_cells_plot = min(num_cells, max_cells);
    num_frames_plot = min(num_frames, max_frames);

    avg_activity = mean(concat_data, 2);
    [~, sort_order] = sort(avg_activity, 'descend');
    sorted_data = concat_data(sort_order, :);

    figure('Name', ['Heatmap - ' stimulus_type ' - ' label_name], 'NumberTitle', 'off');
    sampling_rate = 10;
    x_axis = (1:num_frames_plot) / sampling_rate;
    y_axis = 1:num_cells_plot;

    imagesc(x_axis, y_axis, sorted_data(1:num_cells_plot, 1:num_frames_plot));
    colormap(flipud(gray));
    set(gca, 'YDir', 'reverse');

    title(['Heatmap of Calcium Activity - ' strrep(stimulus_type, '_', ' ') ' - ' label_name]);
    xlabel('Time (s)');
    ylabel('Neuron Index');

    c = colorbar;
    caxis([0 5]);
    c.Label.String = 'dF/F';
    c.Label.FontSize = 12;

    if save_figures
        fname = sprintf('heatmap_%s_%s.png', stimulus_type, regexprep(lower(label_name), '\\s+', '_'));
        saveas(gcf, fullfile(out_dir, fname));
    end
end

function [h, p] = safe_ttest2(x, y)
    try
        [h, p] = ttest2(x, y);
    catch
        h = NaN;
        p = NaN;
    end
end

function [h, p, ks2stat] = safe_kstest2(x, y)
    try
        [h, p, ks2stat] = kstest2(x, y);
    catch
        h = NaN;
        p = NaN;
        ks2stat = NaN;
    end
end
