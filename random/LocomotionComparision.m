%% Locomotion Comparison Across Conditions
function LocomotionComparision(Behav, LocomotionCal, params, conditions, Skip)
% LOCOMOTIONCOMPARISION - Analyzes locomotion speed patterns across different behavioral conditions
%
% INPUTS:
%   Behav        - Structure containing behavioral data for different conditions
%   LocomotionCal- Calibrated locomotion data structure
%   params       - Structure containing plotting parameters including colors for each condition
%   conditions   - Cell array of condition names to analyze (e.g., {'Naive', 'Beginner', 'Expert', 'NoSpout'})
%   Skip         - (optional) Structure containing recording indices to skip for each condition
%
% OUTPUTS:
%   • Figure 1 : Locomotion speed distribution violin plot across conditions
%   • Figure 2 : Alternative locomotion speed distribution plot
%
% EXAMPLE:
%   LocomotionComparision(Behav, LocomotionCal, params);
%   LocomotionComparision(Behav, LocomotionCal, params, {'Naive', 'Expert'});
%   LocomotionComparision(Behav, LocomotionCal, params, {'Naive', 'Expert'}, Skip);

%% FramRate Info
LocomotionFrameRate = 50; %Hz

%% Check inputs and set defaults
if nargin < 5 || isempty(Skip)
    Skip = struct(); 
end

if nargin < 4 || isempty(conditions)
    % Default to all available conditions if not specified
    conditions = {};
    if isfield(Behav, 'Naive'); conditions{end+1} = 'Naive'; end
    if isfield(Behav, 'Beginner'); conditions{end+1} = 'Beginner'; end
    if isfield(Behav, 'Expert'); conditions{end+1} = 'Expert'; end
    if isfield(Behav, 'NoSpout'); conditions{end+1} = 'NoSpout'; end
end

fprintf('Processing locomotion comparison for conditions: %s\n', strjoin(conditions, ', '));

% Print Skip information if provided
if ~isempty(fieldnames(Skip))
    fprintf('Skip information provided for conditions: %s\n', strjoin(fieldnames(Skip), ', '));
    for skip_field = fieldnames(Skip)'
        fprintf('  %s: skipping recordings %s\n', skip_field{1}, mat2str(Skip.(skip_field{1})));
    end
else
    fprintf('No recordings will be skipped (Skip structure is empty)\n');
end

%% Locomotion Analysis
fprintf('\n=== Processing Locomotion Analysis ===\n');

all_locomotion_data = [];
all_locomotion_labels = {};
all_locomotion_mean_data = [];

for c_idx = 1:length(conditions)
    condName = conditions{c_idx};
    fprintf('Processing locomotion data for condition: %s\n', condName);
    
    current_locomotion_data = [];
    current_mean_data = [];
    
    % Get recordings to skip for this condition
    recordingsToSkip = [];
    if isfield(Skip, condName); recordingsToSkip = Skip.(condName); end
    if ~isempty(recordingsToSkip)
        fprintf('  Skipping recordings %s for %s (locomotion analysis)\n', mat2str(recordingsToSkip), condName);
    end
    
    % Check if condition exists in LocomotionCal structure
    if isfield(LocomotionCal, condName) && ~isempty(LocomotionCal.(condName))
        condData = LocomotionCal.(condName);
        
        % Handle different data structures for different conditions
        if strcmp(condName, 'Beginner') || strcmp(condName, 'Expert')
            % These conditions have cell arrays with nested structures
            if iscell(condData)
                num_sessions = min(length(condData), 3); % Limit as in original code
                for i = 1:num_sessions
                    if ismember(i, recordingsToSkip); continue; end
                    if ~isempty(condData{1, i})
                        for a = 1:size(condData{1, i}, 2)
                            if isfield(condData{1, i}(a), 'Forward_mm') && ~isempty(condData{1, i}(a).Forward_mm)
                                current_mean_data = [current_mean_data, mean(condData{1, i}(a).Forward_mm * LocomotionFrameRate, 1)];
                                current_locomotion_data = [current_locomotion_data; condData{1, i}(a).Forward_mm * LocomotionFrameRate];
                            end
                        end
                    end
                end
            end
        elseif strcmp(condName, 'Naive') || strcmp(condName, 'NoSpout')
            % These conditions have direct structure arrays with trial segmentation
            for a = 1:size(condData, 2)
                if ismember(a, recordingsToSkip); continue; end
                if isfield(condData(a), 'Forward_mm') && ~isempty(condData(a).Forward_mm)
                    TrialLength = floor(length(condData(a).Forward_mm) / 150);
                    if TrialLength > 0
                        TrialMovement = 1:TrialLength:length(condData(a).Forward_mm);
                        current_locomotion_data = [current_locomotion_data; condData(a).Forward_mm * LocomotionFrameRate];
                        
                        % Extract means for each trial segment
                        for t = 1:150
                            if t < length(TrialMovement)
                                current_mean_data = [current_mean_data, mean(condData(a).Forward_mm(TrialMovement(t):TrialMovement(t+1))* LocomotionFrameRate, 1)];
                            end
                        end
                    end
                end
            end
        end
    end
    
    % Add data to combined arrays if we found any
    if ~isempty(current_locomotion_data)
 
        
        % Filter outliers
        current_locomotion_data(current_locomotion_data < 0) = NaN;
        % current_locomotion_data(current_locomotion_data > 120) = 1;
        current_mean_data(current_mean_data < 0) = NaN;
        % current_mean_data(current_mean_data > 25) = 1;
        
        all_locomotion_data = [all_locomotion_data; current_locomotion_data];
        all_locomotion_mean_data = [all_locomotion_mean_data, current_mean_data];
        
        % Create labels for this condition's data
        num_trials = size(current_locomotion_data, 1);
        all_locomotion_labels = [all_locomotion_labels; repmat({condName}, num_trials, 1)];
        
        fprintf('  Added %d locomotion trials for %s\n', num_trials, condName);
    else
        warning('No locomotion data found for condition: %s', condName);
    end
end

% Create locomotion plot if we have data
if ~isempty(all_locomotion_data)
    fprintf('Creating locomotion plot with %d total trials across %d conditions\n', size(all_locomotion_data, 1), length(conditions));
    
    % Get unique conditions and create color mapping
    unique_conds_locomotion = unique(all_locomotion_labels, 'stable');
    colors_locomotion = [];
    for i = 1:length(unique_conds_locomotion)
        if isfield(params, unique_conds_locomotion{i})
            colors_locomotion = [colors_locomotion; params.(unique_conds_locomotion{i})];
        else
            colors_locomotion = [colors_locomotion; 0.5 0.5 0.5]; % Default gray
        end
    end
    
    % Create group indices
    group_idx_locomotion = zeros(size(all_locomotion_data, 1), 1);
    for i = 1:length(all_locomotion_labels)
        for c_map = 1:length(unique_conds_locomotion)
            if strcmp(all_locomotion_labels{i}, unique_conds_locomotion{c_map})
                group_idx_locomotion(i) = c_map;
                break;
            end
        end
    end
    
    % Convert to proper format for plotting (flatten)
    locomotion_data_flat = all_locomotion_data(:);
    group_idx_flat = repmat(group_idx_locomotion, size(all_locomotion_data, 2), 1);
    
    fig = figure('Position', [134 160 1311 675], 'Color', 'w', 'Name', 'Locomotion Speed Comparison Across Conditions');
    h = daviolinplot(locomotion_data_flat, 'groups', group_idx_flat, 'colors', colors_locomotion, ...
                     'box', 2, 'boxcolors', 'same', 'xtlabels', unique_conds_locomotion, 'violin', 'full');
    ylabel('Speed [cm/sec]');
    title('Locomotion Speed Distribution Across Behavioral Conditions');
    xl = xlim; xlim([xl(1)-0.1, xl(2)+0.3]); % make more space for the legend
    set(gca, 'FontSize', 24);
    ylim([-1 15]);
    
    % Alternative violin plot
    figure('Color', 'w', 'Name', 'Alternative Locomotion Speed Comparison');
    violins = violinplot(group_idx_flat, locomotion_data_flat);
    ylim([-1 15]);
    xticklabels(unique_conds_locomotion);
    ylabel('Speed [cm/sec]');
    title('Locomotion Speed Distribution (Alternative View)');
    set(gca, 'FontSize', 14);
    
    %% New Enhanced Analysis Plots
    fprintf('\n=== Creating Enhanced Locomotion Analysis Plots ===\n');
    
    % Create new plots with detailed locomotion analysis
    createEnhancedLocomotionPlots(LocomotionCal, unique_conds_locomotion, colors_locomotion, conditions, Skip, LocomotionFrameRate);
    
else
    warning('No locomotion data available for any conditions');
end

fprintf('\nLocomotionComparision completed successfully!\n');
fprintf('Processed %d conditions: %s\n', length(conditions), strjoin(conditions, ', '));

% Summary of skip information
if ~isempty(fieldnames(Skip))
    fprintf('\nSummary of skipped recordings:\n');
    for skip_field = fieldnames(Skip)'
        fprintf('  %s: %d recording(s) skipped (%s)\n', skip_field{1}, length(Skip.(skip_field{1})), mat2str(Skip.(skip_field{1})));
    end
else
    fprintf('No recordings were skipped in this analysis.\n');
end

end

%% Helper function for enhanced locomotion analysis
function createEnhancedLocomotionPlots(LocomotionCal, unique_conditions, colors, conditions, Skip,LocomotionFrameRate)

% Initialize storage for new analyses
event_counts = struct('starts', [], 'stops', [], 'condition_names', {{}});
avg_speeds = [];
rest_moving_ratios = [];
condition_labels = {};

% Movement thresholds
movement_threshold = 1; % cm/s for defining movement vs rest
min_event_duration = 5; % minimum frames for an event to count

fprintf('Processing enhanced locomotion analysis...\n');

for c_idx = 1:length(unique_conditions)
    condName = unique_conditions{c_idx};
    fprintf('  Analyzing %s...\n', condName);
    
    % Get skip list
    recordingsToSkip = [];
    if isfield(Skip, condName); recordingsToSkip = Skip.(condName); end
    
    % Initialize condition-specific variables
    cond_starts = 0;
    cond_stops = 0;
    cond_speeds = [];
    cond_rest_time = 0;
    cond_moving_time = 0;
    
    if isfield(LocomotionCal, condName) && ~isempty(LocomotionCal.(condName))
        condData = LocomotionCal.(condName);
        
        % Process based on condition type (same logic as main function)
        if strcmp(condName, 'Beginner') || strcmp(condName, 'Expert')
            % Cell array structure
            if iscell(condData)
                num_sessions = min(length(condData), 3);
                for i = 1:num_sessions
                    if ismember(i, recordingsToSkip) || isempty(condData{1, i}); continue; end
                    
                    for a = 1:size(condData{1, i}, 2)
                        if isfield(condData{1, i}(a), 'Forward_mm') && ~isempty(condData{1, i}(a).Forward_mm)
                            locomotionTrace = condData{1, i}(a).Forward_mm * LocomotionFrameRate;
                            [starts, stops, avg_speed, rest_ratio] = analyzeLocomotionEvents(locomotionTrace, movement_threshold, min_event_duration);
                            
                            cond_starts = cond_starts + starts;
                            cond_stops = cond_stops + stops;
                            cond_speeds = [cond_speeds, avg_speed];
                            cond_rest_time = cond_rest_time + rest_ratio * length(locomotionTrace);
                            cond_moving_time = cond_moving_time + (1-rest_ratio) * length(locomotionTrace);
                        end
                    end
                end
            end
            
        elseif strcmp(condName, 'Naive') || strcmp(condName, 'NoSpout')
            % Structure array with trial segmentation
            for a = 1:size(condData, 2)
                if ismember(a, recordingsToSkip); continue; end
                
                if isfield(condData(a), 'Forward_mm') && ~isempty(condData(a).Forward_mm)
                    locomotionTrace = condData(a).Forward_mm ;
                    
                        % Apply scaling as in main function
                        scaled_trace = locomotionTrace * LocomotionFrameRate ;
                        [starts, stops, avg_speed, rest_ratio] = analyzeLocomotionEvents(scaled_trace, movement_threshold, min_event_duration);
                        
                        cond_starts = cond_starts + starts;
                        cond_stops = cond_stops + stops;
                        cond_speeds = [cond_speeds, avg_speed];
                        cond_rest_time = cond_rest_time + rest_ratio * length(scaled_trace);
                        cond_moving_time = cond_moving_time + (1-rest_ratio) * length(scaled_trace);
                    
                end
            end
        end
    end
    
    % Store results for this condition
    event_counts.starts(end+1) = cond_starts;
    event_counts.stops(end+1) = cond_stops;
    event_counts.condition_names{end+1} = condName;
    
    avg_speeds(end+1) = mean(cond_speeds, 'omitnan');
    
    total_time = cond_rest_time + cond_moving_time;
    if total_time > 0
        rest_moving_ratios(end+1) = cond_rest_time / total_time;
    else
        rest_moving_ratios(end+1) = NaN;
    end
    condition_labels{end+1} = condName;
end

%% Plot 1: Movement Events (Starts/Stops) per Condition
figure('Color', 'w', 'Name', 'Movement Events per Condition', 'Position', [100 100 800 600]);

subplot(2,1,1);
bar(1:length(unique_conditions), event_counts.starts, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'k');
set(gca, 'XTickLabel', unique_conditions, 'FontSize', 12);
ylabel('Number of Movement Starts');
title('Movement Start Events per Condition');
grid on;

subplot(2,1,2);
bar(1:length(unique_conditions), event_counts.stops, 'FaceColor', [0.8 0.4 0.2], 'EdgeColor', 'k');
set(gca, 'XTickLabel', unique_conditions, 'FontSize', 12);
ylabel('Number of Movement Stops');
title('Movement Stop Events per Condition');
grid on;

%% Plot 2: Average Speed per Condition
figure('Color', 'w', 'Name', 'Average Speed per Condition', 'Position', [200 100 600 500]);

bar_handle = bar(1:length(unique_conditions), avg_speeds, 'EdgeColor', 'k', 'LineWidth', 1.5);
for i = 1:length(unique_conditions)
    bar_handle.FaceColor = 'flat';
    bar_handle.CData(i,:) = colors(i,:);
end

set(gca, 'XTickLabel', unique_conditions, 'FontSize', 14);
ylabel('Average Speed [cm/s]', 'FontSize', 14);
title('Average Locomotion Speed per Condition', 'FontSize', 16);
grid on;

% Add value labels on bars
for i = 1:length(avg_speeds)
    text(i, avg_speeds(i) + 0.1, sprintf('%.2f', avg_speeds(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
end

%% Plot 3: Rest vs Moving Period Ratio per Condition
figure('Color', 'w', 'Name', 'Rest vs Moving Ratio per Condition', 'Position', [300 100 700 500]);

% Create stacked bar chart
moving_ratios = 1 - rest_moving_ratios;
stacked_data = [rest_moving_ratios(:), moving_ratios(:)]';

bar_handle = bar(1:length(unique_conditions), stacked_data', 'stacked', 'EdgeColor', 'k', 'LineWidth', 1.5);
bar_handle(1).FaceColor = [0.7 0.7 0.9]; % Rest periods - light blue
bar_handle(2).FaceColor = [0.9 0.7 0.7]; % Moving periods - light red

set(gca, 'XTickLabel', unique_conditions, 'FontSize', 14);
ylabel('Proportion of Time', 'FontSize', 14);
title('Rest vs Moving Time Proportions per Condition', 'FontSize', 16);
legend({'Rest', 'Moving'}, 'Location', 'best', 'FontSize', 12);
grid on;
ylim([0 1]);

% Add percentage labels
for i = 1:length(unique_conditions)
    rest_pct = rest_moving_ratios(i) * 100;
    moving_pct = moving_ratios(i) * 100;
    
    % Label for rest portion
    text(i, rest_moving_ratios(i)/2, sprintf('%.1f%%', rest_pct), ...
         'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
    
    % Label for moving portion
    text(i, rest_moving_ratios(i) + moving_ratios(i)/2, sprintf('%.1f%%', moving_pct), ...
         'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
end

%% Summary Statistics
fprintf('\n=== Enhanced Locomotion Analysis Summary ===\n');
for i = 1:length(unique_conditions)
    fprintf('%s:\n', unique_conditions{i});
    fprintf('  Movement starts: %d\n', event_counts.starts(i));
    fprintf('  Movement stops: %d\n', event_counts.stops(i));
    fprintf('  Average speed: %.2f cm/s\n', avg_speeds(i));
    fprintf('  Rest ratio: %.1f%% | Moving ratio: %.1f%%\n', ...
            rest_moving_ratios(i)*100, (1-rest_moving_ratios(i))*100);
    fprintf('\n');
end

end

%% Helper function to analyze locomotion events in a trace
function [starts, stops, avg_speed, rest_ratio] = analyzeLocomotionEvents(locomotionTrace, threshold, min_duration)

% Remove NaN values and outliers
clean_trace = locomotionTrace(~isnan(locomotionTrace));
clean_trace(clean_trace < 0) = 0;
clean_trace(clean_trace > 25) = 25; % Cap at reasonable maximum

if isempty(clean_trace)
    starts = 0; stops = 0; avg_speed = 0; rest_ratio = 1;
    return;
end

% Identify movement vs rest periods
is_moving = clean_trace > threshold;

% Find transitions
movement_starts = find(diff([0; is_moving]) == 1);
movement_stops = find(diff([is_moving; 0]) == -1);

% Filter out very short events
valid_starts = [];
valid_stops = [];

for i = 1:length(movement_starts)
    if i <= length(movement_stops)
        event_duration = movement_stops(i) - movement_starts(i) + 1;
        if event_duration >= min_duration
            valid_starts(end+1) = movement_starts(i);
            valid_stops(end+1) = movement_stops(i);
        end
    end
end

starts = length(valid_starts);
stops = length(valid_stops);

% Calculate average speed
avg_speed = mean(clean_trace, 'omitnan');

% Calculate rest ratio
rest_frames = sum(~is_moving);
total_frames = length(is_moving);
rest_ratio = rest_frames / total_frames;

end 