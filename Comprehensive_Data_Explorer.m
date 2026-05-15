clear; clc; close all;

%% ====================== USER CONFIGURATION ======================

% Data files
baseline_file = "Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260305_1327_preprocessed.mat";
drug_file     = "Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260305_1409_preprocessed.mat";

% ROI matching (optional, from ROIMatchPub)
roi_match_file = "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane0_V2.mat";

% Calcium data type (1=FVDff, 2=deconvolved, 3=F)
ca_type = 1;

% Stimulus configuration
selected_stimulus_type = 'sparse_local_global_flashes';  % e.g., 'grating', 'moving_bar', 'full_field_flash'
analysis_plane = 1;  % Plane index (1-based for internal indexing)

% Response window for SNR calculation (frames before and after stimulus onset)
response_window_frames = 50;  % frames to capture post-stimulus onset

% Visualization settings
response_threshold = 1.0;  % dF/F threshold to classify responding ROIs
n_top_rois_to_plot = 20;   % Number of top-SNR ROIs to visualize
figure_quality = 150;       % DPI for saved figures

%% ====================== DATA LOADING & VALIDATION ======================

fprintf('\n========== DATA LOADER ==========\n');

% Load baseline
fprintf('Loading baseline: %s\n', baseline_file);
B = load(baseline_file);
fprintf('  - Planes: %d\n', numel(B.CaData));
fprintf('  - Stimulus types: %s\n', strjoin(unique({B.Stimuli(:).type}), ', '));

% Load drug
fprintf('Loading drug: %s\n', drug_file);
D = load(drug_file);
fprintf('  - Planes: %d\n', numel(D.CaData));
fprintf('  - Stimulus types: %s\n', strjoin(unique({D.Stimuli(:).type}), ', '));

% Validate
if numel(B.CaData) < analysis_plane || numel(D.CaData) < analysis_plane
    error('Requested plane %d not available', analysis_plane);
end

% Extract data for selected plane and calcium type
base_ca = B.CaData(analysis_plane);
drug_ca = D.CaData(analysis_plane);

base_dff = get_calcium_data(base_ca, ca_type);
drug_dff = get_calcium_data(drug_ca, ca_type);

n_rois_base = size(base_dff, 1);
n_rois_drug = size(drug_dff, 1);

fprintf('\nData dimensions:\n');
fprintf('  Baseline: %d ROIs × %d frames\n', n_rois_base, size(base_dff, 2));
fprintf('  Drug:     %d ROIs × %d frames\n', n_rois_drug, size(drug_dff, 2));

%% ====================== STIMULUS ALIGNMENT ======================

fprintf('\n========== STIMULUS ALIGNMENT ==========\n');
fprintf('Searching for "%s" presentations...\n', selected_stimulus_type);

% Extract stimulus frame ranges
base_stim_info = extract_stimulus_presentations(B.Stimuli, selected_stimulus_type, response_window_frames);
drug_stim_info = extract_stimulus_presentations(D.Stimuli, selected_stimulus_type, response_window_frames);

n_stim_base = size(base_stim_info, 1);
n_stim_drug = size(drug_stim_info, 1);

fprintf('  Baseline: %d presentations\n', n_stim_base);
fprintf('  Drug:     %d presentations\n', n_stim_drug);

if n_stim_base == 0 || n_stim_drug == 0
    error('No "%s" stimulus presentations found', selected_stimulus_type);
end

%% ====================== POPULATION ANALYSIS ======================

fprintf('\n========== POPULATION ANALYSIS ==========\n');

% Calculate response metrics for all ROIs
[base_metrics] = calculate_roi_metrics(base_dff, base_stim_info, response_threshold);
[drug_metrics] = calculate_roi_metrics(drug_dff, drug_stim_info, response_threshold);

% Summary statistics
fprintf('\nBASELINE Population Statistics:\n');
print_population_summary(base_metrics, selected_stimulus_type);

fprintf('\nDRUG Population Statistics:\n');
print_population_summary(drug_metrics, selected_stimulus_type);

%% ====================== MATCHED ROI ANALYSIS (OPTIONAL) ======================

matched_rois_exist = ~isempty(roi_match_file) && isfile(roi_match_file);

if matched_rois_exist
    fprintf('\n========== MATCHED ROI ANALYSIS ==========\n');
    load(roi_match_file, 'allSessionMapping');
    
    % Extract matched ROI indices
    baseline_roi_idx = allSessionMapping(:, 1);
    drug_roi_idx = allSessionMapping(:, 2);
    
    n_matched = length(baseline_roi_idx);
    fprintf('Comparing %d matched ROIs\n', n_matched);
    
    % Get metrics for matched pairs
    matched_base_snr = base_metrics.snr(baseline_roi_idx);
    matched_drug_snr = drug_metrics.snr(drug_roi_idx);
    
    matched_base_response = base_metrics.mean_response(baseline_roi_idx);
    matched_drug_response = drug_metrics.mean_response(drug_roi_idx);
    
    % Calculate changes
    snr_change = matched_drug_snr - matched_base_snr;
    response_change = matched_drug_response - matched_base_response;
    
    fprintf('Matched ROI Statistics:\n');
    fprintf('  SNR change (drug - baseline):\n');
    fprintf('    Mean: %.3f\n', mean(snr_change));
    fprintf('    Median: %.3f\n', median(snr_change));
    fprintf('    Range: [%.3f, %.3f]\n', min(snr_change), max(snr_change));
    
    fprintf('  Response magnitude change (drug - baseline):\n');
    fprintf('    Mean: %.3f %%\n', 100 * mean(response_change));
    fprintf('    Median: %.3f %%\n', 100 * median(response_change));
end

%% ====================== VISUALIZATIONS ======================

fprintf('\n========== GENERATING VISUALIZATIONS ==========\n');

% Figure 1: Population Distributions
fig1 = figure('Position', [100 100 1400 600], 'NumberTitle', 'off', 'Name', 'Population Analysis');

% SNR distributions
subplot(2, 3, 1);
histogram(base_metrics.snr, 30, 'FaceColor', [1 0.4 0.4], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
hold on;
histogram(drug_metrics.snr, 30, 'FaceColor', [0.4 0.4 1], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
xlabel('SNR');
ylabel('# ROIs');
title('Signal-to-Noise Ratio Distribution');
legend('Baseline', 'Drug', 'Location', 'northeast');
set(gca, 'LineWidth', 1.5, 'FontSize', 10);

% Response magnitude distributions
subplot(2, 3, 2);
histogram(base_metrics.mean_response, 30, 'FaceColor', [1 0.4 0.4], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
hold on;
histogram(drug_metrics.mean_response, 30, 'FaceColor', [0.4 0.4 1], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
xlabel('Mean dF/F');
ylabel('# ROIs');
title('Response Magnitude Distribution');
legend('Baseline', 'Drug', 'Location', 'northeast');
set(gca, 'LineWidth', 1.5, 'FontSize', 10);

% Fraction responding
subplot(2, 3, 3);
data_fr = [base_metrics.fraction_responding, drug_metrics.fraction_responding];
bar(data_fr, 'BarWidth', 0.6);
set(gca, 'XTickLabel', {'Baseline', 'Drug'});
ylabel('Fraction of ROIs');
title(sprintf('Responsive ROIs (threshold: %.1f dF/F)', response_threshold));
ylim([0 1]);
set(gca, 'LineWidth', 1.5, 'FontSize', 10);

% Response latency distributions
subplot(2, 3, 4);
histogram(base_metrics.response_latency, 20, 'FaceColor', [1 0.4 0.4], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
hold on;
histogram(drug_metrics.response_latency, 20, 'FaceColor', [0.4 0.4 1], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
xlabel('Latency (frames)');
ylabel('# ROIs');
title('Response Latency Distribution');
legend('Baseline', 'Drug', 'Location', 'northeast');
set(gca, 'LineWidth', 1.5, 'FontSize', 10);

% Cumulative distribution of SNR
subplot(2, 3, 5);
[f_base, x_base] = ecdf(base_metrics.snr);
[f_drug, x_drug] = ecdf(drug_metrics.snr);
plot(x_base, f_base, 'Color', [1 0 0], 'LineWidth', 2.5);
hold on;
plot(x_drug, f_drug, 'Color', [0 0 1], 'LineWidth', 2.5);
xlabel('SNR');
ylabel('Cumulative Fraction');
title('SNR Cumulative Distribution');
legend('Baseline', 'Drug', 'Location', 'southeast');
grid on;
set(gca, 'LineWidth', 1.5, 'FontSize', 10);

% Stats box
subplot(2, 3, 6);
axis off;
stats_text = sprintf(['Population Summary (%s)\n\n' ...
    'Baseline:\n' ...
    '  Mean SNR: %.2f\n' ...
    '  Responding: %.1f%%\n' ...
    '\nDrug:\n' ...
    '  Mean SNR: %.2f\n' ...
    '  Responding: %.1f%%'], ...
    selected_stimulus_type, ...
    mean(base_metrics.snr), 100*base_metrics.fraction_responding, ...
    mean(drug_metrics.snr), 100*drug_metrics.fraction_responding);
text(0.1, 0.5, stats_text, 'FontSize', 11, 'VerticalAlignment', 'middle', ...
    'FontFamily', 'monospaced', 'BackgroundColor', [0.95 0.95 0.95], ...
    'EdgeColor', [0.5 0.5 0.5], 'Margin', 10);

sgtitle(sprintf('Population-Level Analysis: %s', selected_stimulus_type), ...
    'FontSize', 13, 'FontWeight', 'bold');

% Figure 2: Top SNR ROIs - Individual Responses
fig2 = figure('Position', [100 100 1400 800], 'NumberTitle', 'off', 'Name', 'Top ROI Responses');

[~, idx_base] = sort(base_metrics.snr, 'descend');
[~, idx_drug] = sort(drug_metrics.snr, 'descend');

top_idx_base = idx_base(1:min(n_top_rois_to_plot, n_rois_base));
top_idx_drug = idx_drug(1:min(n_top_rois_to_plot, n_rois_drug));

n_rows = ceil(n_top_rois_to_plot / 2);

for roi_num = 1:min(n_top_rois_to_plot, n_rois_base)
    % Baseline
    subplot(n_rows, 4, 2*roi_num - 1);
    roi_idx = top_idx_base(roi_num);
    plot_roi_responses(base_dff, base_stim_info, roi_idx, [1 0.4 0.4], 2.5);
    if roi_num == 1
        title(sprintf('Baseline Top %d ROIs', n_top_rois_to_plot), 'FontWeight', 'bold');
    end
    ylabel(sprintf('ROI %d\nSNR:%.2f', roi_idx, base_metrics.snr(roi_idx)));
end

for roi_num = 1:min(n_top_rois_to_plot, n_rois_drug)
    % Drug
    subplot(n_rows, 4, 2*roi_num);
    roi_idx = top_idx_drug(roi_num);
    plot_roi_responses(drug_dff, drug_stim_info, roi_idx, [0.4 0.4 1], 2.5);
    if roi_num == 1
        title(sprintf('Drug Top %d ROIs', n_top_rois_to_plot), 'FontWeight', 'bold');
    end
    ylabel(sprintf('ROI %d\nSNR:%.2f', roi_idx, drug_metrics.snr(roi_idx)));
end

sgtitle(sprintf('Top SNR Responses: %s', selected_stimulus_type), ...
    'FontSize', 13, 'FontWeight', 'bold');

% Figure 3: Matched ROI Analysis (if available)
if matched_rois_exist && n_matched > 0
    fig3 = figure('Position', [100 100 1200 800], 'NumberTitle', 'off', 'Name', 'Matched ROI Analysis');
    
    % SNR comparison scatter
    subplot(2, 3, 1);
    scatter(matched_base_snr, matched_drug_snr, 30, 'filled', ...
        'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k');
    hold on;
    lims = [min([matched_base_snr; matched_drug_snr]), max([matched_base_snr; matched_drug_snr])];
    plot(lims, lims, 'k--', 'LineWidth', 1.5);
    xlabel('Baseline SNR');
    ylabel('Drug SNR');
    title('SNR: Baseline vs Drug');
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    grid on;
    
    % Response magnitude comparison
    subplot(2, 3, 2);
    scatter(matched_base_response, matched_drug_response, 30, 'filled', ...
        'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k');
    hold on;
    lims = [min([matched_base_response; matched_drug_response]), max([matched_base_response; matched_drug_response])];
    plot(lims, lims, 'k--', 'LineWidth', 1.5);
    xlabel('Baseline Response (dF/F)');
    ylabel('Drug Response (dF/F)');
    title('Response Magnitude: Baseline vs Drug');
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    grid on;
    
    % SNR change distribution
    subplot(2, 3, 3);
    histogram(snr_change, 20, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'k', 'FaceAlpha', 0.7);
    xlabel('SNR Change (Drug - Baseline)');
    ylabel('# ROIs');
    title('SNR Changes');
    axline(0, 'k--', 'LineWidth', 1.5);
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    
    % Response change distribution
    subplot(2, 3, 4);
    histogram(response_change, 20, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'k', 'FaceAlpha', 0.7);
    xlabel('Response Change (Drug - Baseline, dF/F)');
    ylabel('# ROIs');
    title('Response Magnitude Changes');
    axline(0, 'k--', 'LineWidth', 1.5);
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    
    % Summary statistics
    subplot(2, 3, 5:6);
    axis off;
    summary_text = sprintf(['Matched ROI Summary (%d ROIs)\n\n' ...
        'SNR Change:\n' ...
        '  Mean: %.3f\n' ...
        '  Median: %.3f\n' ...
        '  % Increased: %.1f%%\n\n' ...
        'Response Change:\n' ...
        '  Mean: %.3f dF/F\n' ...
        '  Median: %.3f dF/F\n' ...
        '  % Enhanced: %.1f%%'], ...
        n_matched, ...
        mean(snr_change), median(snr_change), 100*sum(snr_change>0)/n_matched, ...
        mean(response_change), median(response_change), 100*sum(response_change>0)/n_matched);
    text(0.1, 0.5, summary_text, 'FontSize', 12, 'VerticalAlignment', 'middle', ...
        'FontFamily', 'monospaced', 'BackgroundColor', [0.95 0.95 0.95], ...
        'EdgeColor', [0.5 0.5 0.5], 'Margin', 15);
    
    sgtitle('Matched ROI Comparisons: Before vs After Drug', 'FontSize', 13, 'FontWeight', 'bold');
end

fprintf('Visualizations complete!\n');

%% ====================== EXPORT SUMMARY ======================

fprintf('\n========== ANALYSIS COMPLETE ==========\n');
fprintf('Figures generated:\n');
fprintf('  - Figure 1: Population distributions and statistics\n');
fprintf('  - Figure 2: Top SNR ROI responses\n');
if matched_rois_exist
    fprintf('  - Figure 3: Matched ROI comparisons\n');
end
fprintf('\nAnalysis summary saved.\n');

%% ====================== LOCAL FUNCTIONS ======================

function dff = get_calcium_data(ca_data, ca_type)
    % Extract calcium data by type
    % ca_type: 1=FVDff, 2=deconvolved, 3=F
    
    switch ca_type
        case 1
            dff = ca_data.Ca_dFF;
        case 2
            dff = ca_data.Ca_deconvolved;
        case 3
            dff = ca_data.Ca_F;
        otherwise
            error('Unknown calcium type: %d', ca_type);
    end
end

function stim_info = extract_stimulus_presentations(Stimuli, stim_type, response_window_frames)
    % Extract frame ranges for all presentations of a stimulus type
    % Returns: Nx2 array [start_frame, end_frame + response_window]
    
    stim_info = [];
    
    for i = 1:numel(Stimuli)
        if ~isfield(Stimuli(i), 'type') || ~strcmp(Stimuli(i).type, stim_type)
            continue;
        end
        
        if ~isfield(Stimuli(i), 'TimeStimulusFrame') || isempty(Stimuli(i).TimeStimulusFrame)
            continue;
        end
        
        stim_frames = Stimuli(i).TimeStimulusFrame;
        st_frame = min(stim_frames(:));
        en_frame = max(stim_frames(:));
        
        % Add response window
        en_frame_extended = en_frame + response_window_frames;
        
        stim_info = [stim_info; st_frame, en_frame_extended];
    end
end

function metrics = calculate_roi_metrics(dff, stim_info, response_threshold)
    % Calculate response metrics for all ROIs across stimulus presentations
    
    n_rois = size(dff, 1);
    n_stim = size(stim_info, 1);
    
    metrics.snr = zeros(n_rois, 1);
    metrics.mean_response = zeros(n_rois, 1);
    metrics.response_latency = zeros(n_rois, 1);
    metrics.fraction_responding = 0;
    
    % Calculate metrics for each ROI
    for roi = 1:n_rois
        snr_values = [];
        response_values = [];
        latencies = [];
        
        for stim = 1:n_stim
            st = stim_info(stim, 1);
            en = min(stim_info(stim, 2), size(dff, 2));
            
            if en > st
                response_window = dff(roi, st:en);
                
                % SNR: mean / std
                signal = mean(response_window);
                noise = std(response_window);
                snr_values = [snr_values; signal / (noise + eps)];
                
                response_values = [response_values; signal];
                
                % Latency: frame when response exceeds threshold
                exceeds = response_window >= response_threshold;
                if any(exceeds)
                    latencies = [latencies; find(exceeds, 1)];
                end
            end
        end
        
        % Average across presentations
        if ~isempty(snr_values)
            metrics.snr(roi) = mean(snr_values);
            metrics.mean_response(roi) = mean(response_values);
        end
        
        if ~isempty(latencies)
            metrics.response_latency(roi) = mean(latencies);
        else
            metrics.response_latency(roi) = NaN;
        end
    end
    
    % Fraction responding
    metrics.fraction_responding = sum(metrics.mean_response >= response_threshold) / n_rois;
end

function print_population_summary(metrics, stim_type)
    % Print population-level summary statistics
    
    fprintf('  Mean SNR: %.3f (± %.3f)\n', mean(metrics.snr), std(metrics.snr));
    fprintf('  Median SNR: %.3f\n', median(metrics.snr));
    fprintf('  SNR range: [%.3f, %.3f]\n', min(metrics.snr), max(metrics.snr));
    
    fprintf('  Mean response: %.3f dF/F (± %.3f)\n', mean(metrics.mean_response), std(metrics.mean_response));
    fprintf('  Fraction responding: %.1f%%\n', 100 * metrics.fraction_responding);
    
    latency_valid = metrics.response_latency(~isnan(metrics.response_latency));
    if ~isempty(latency_valid)
        fprintf('  Mean latency: %.1f frames\n', mean(latency_valid));
    end
end

function plot_roi_responses(dff, stim_info, roi_idx, color, linewidth)
    % Plot all responses for a single ROI across stimulus presentations
    
    hold on;
    for stim = 1:size(stim_info, 1)
        st = stim_info(stim, 1);
        en = min(stim_info(stim, 2), size(dff, 2));
        
        if en > st
            response = dff(roi_idx, st:en);
            plot(response, 'Color', color, 'LineWidth', linewidth, 'LineStyle', '-');
        end
    end
    
    % Overlay mean
    all_responses = [];
    for stim = 1:size(stim_info, 1)
        st = stim_info(stim, 1);
        en = min(stim_info(stim, 2), size(dff, 2));
        if en > st
            response = dff(roi_idx, st:en);
            % Pad to same length
            all_responses = [all_responses; response];
        end
    end
    
    if ~isempty(all_responses)
        mean_response = mean(all_responses, 1, 'omitnan');
        plot(mean_response, 'Color', color * 0.5, 'LineWidth', linewidth + 1, 'LineStyle', '-');
    end
    
    xlabel('Frame');
    ylabel('dF/F');
    set(gca, 'LineWidth', 1, 'FontSize', 9);
    grid on;
end
