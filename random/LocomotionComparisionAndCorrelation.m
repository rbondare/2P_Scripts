function LocomotionComparisionAndCorrelation(ActivityData, Behav, WindowLength, Rec, Performance, LocomotionCal, params, conditions, Skip)
% LOCOMOTIONCOMPARISIONANDCORRELATION - Analyzes locomotion patterns across different behavioral conditions
%
% INPUTS:
%   ActivityData - Structure containing neural activity data
%   Behav        - Structure containing behavioral data for different conditions
%   WindowLength - Length of analysis window
%   Rec          - Recording metadata structure
%   Performance  - Performance metrics structure
%   LocomotionCal- Calibrated locomotion data structure
%   params       - Structure containing plotting parameters including colors for each condition
%   conditions   - Cell array of condition names to analyze (e.g., {'Naive', 'Beginner', 'Expert', 'NoSpout'})
%   Skip         - (optional) Structure containing recording indices to skip for each condition
%
% EXAMPLE:
%   LocomotionComparisionAndCorrelation(ActivityData, Behav, WindowLength, Rec, Performance, LocomotionCal, params, {'Naive', 'Expert'});
%   LocomotionComparisionAndCorrelation(ActivityData, Behav, WindowLength, Rec, Performance, LocomotionCal, params, {'Naive', 'Expert'}, Skip);

%% Check inputs and set defaults
if nargin < 9 || isempty(Skip)
    Skip = struct(); 
end

if nargin < 8 || isempty(conditions)
    % Default to all available conditions if not specified
    conditions = {};
    if isfield(Behav, 'Naive'); conditions{end+1} = 'Naive'; end
    if isfield(Behav, 'Beginner'); conditions{end+1} = 'Beginner'; end
    if isfield(Behav, 'Expert'); conditions{end+1} = 'Expert'; end
    if isfield(Behav, 'NoSpout'); conditions{end+1} = 'NoSpout'; end
end

fprintf('Processing conditions: %s\n', strjoin(conditions, ', '));

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
                                current_mean_data = [current_mean_data, mean(condData{1, i}(a).Forward_mm, 1)];
                                current_locomotion_data = [current_locomotion_data; condData{1, i}(a).Forward_mm];
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
                        current_locomotion_data = [current_locomotion_data; condData(a).Forward_mm];
                        
                        % Extract means for each trial segment
                        for t = 1:150
                            if t < length(TrialMovement)
                                current_mean_data = [current_mean_data, mean(condData(a).Forward_mm(TrialMovement(t):TrialMovement(t+1)), 1)];
                            end
                        end
                    end
                end
            end
        end
    end
    
    % Add data to combined arrays if we found any
    if ~isempty(current_locomotion_data)
        % Apply scaling and filtering as in original code
        if exist('TrialLength', 'var')
            current_locomotion_data = current_locomotion_data * TrialLength / 12;
            current_mean_data = current_mean_data * TrialLength / 12;
        end
        
        % Filter outliers
        current_locomotion_data(current_locomotion_data < 0) = NaN;
        current_locomotion_data(current_locomotion_data > 25) = 1;
        current_mean_data(current_mean_data < 0) = NaN;
        current_mean_data(current_mean_data > 25) = 1;
        
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
    ylabel('Speed [mm/sec]');
    title('Locomotion Speed Distribution Across Behavioral Conditions');
    xl = xlim; xlim([xl(1)-0.1, xl(2)+0.3]); % make more space for the legend
    set(gca, 'FontSize', 24);
    ylim([-1 15]);
    
    % Alternative violin plot
    figure('Color', 'w', 'Name', 'Alternative Locomotion Speed Comparison');
    violins = violinplot(group_idx_flat, locomotion_data_flat);
    ylim([-1 15]);
    xticklabels(unique_conds_locomotion);
    ylabel('Speed [mm/sec]');
    title('Locomotion Speed Distribution (Alternative View)');
    set(gca, 'FontSize', 14);
else
    warning('No locomotion data available for any conditions');
end

fprintf('\nLocomotionComparisionAndCorrelation completed successfully!\n');
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