function [Fig1_Data] = extract_percentile_based_data(conditions, Data, Behav, Skip, position, Performance, arousal_params, data_type)
% EXTRACT_PERCENTILE_BASED_DATA - Extract data for percentile-based arousal analysis
%
% This function replaces event-based pupil dilation detection with percentile-based
% arousal state analysis, comparing neural activity/entropy between high and low
% arousal states defined by pupil size percentiles.
%
% INPUTS:
%   conditions     - Cell array of condition names to analyze
%   Data   - Structure containing neural activity data (for activity analysis)
%   Behav          - Structure containing behavioral data including pupil size
%   Skip           - Structure containing recording indices to skip for each condition
%   position       - Position string ('All', 'P1', 'P3')
%   Performance    - Performance metrics structure for trial selection
%   arousal_params - Structure containing percentile analysis parameters:
%                   .use_percentile_binning: true to enable percentile mode
%                   .high_arousal_percentile: percentile for high arousal (default 90)
%                   .low_arousal_percentile: percentile for low arousal (default 10)
%                   .min_frames_per_bin: minimum frames required per bin (default 50)
%                   .pooling_level: 'recording' (default), 'condition', or 'trial'
%   data_type      - 'activity' or 'entropy' to specify analysis type
%
% OUTPUTS:
%   Fig1_Data - Structure containing percentile-based arousal data:
%              For activity analysis:
%                .HighArousalActivity: neural activity data for high arousal frames
%                .LowArousalActivity: neural activity data for low arousal frames
%              For entropy analysis:
%                .HighArousalEntropy: entropy data for high arousal frames
%                .LowArousalEntropy: entropy data for low arousal frames
%              Common fields:
%                .HighArousalConditionLabels: condition labels for high arousal data
%                .LowArousalConditionLabels: condition labels for low arousal data
%                .pupil_percentiles: percentile thresholds used
%                .analysis_params: copy of arousal_params used

fprintf('=== Extracting Percentile-Based Arousal Data for %s (%s) ===\n', position, upper(data_type));

% Note: This function now only performs percentile-based analysis

% Set default parameters
if ~isfield(arousal_params, 'min_frames_per_bin')
    arousal_params.min_frames_per_bin = 50;
end

% Calculate percentile thresholds
pupil_percentiles = calculate_pupil_percentiles(conditions, Behav, Skip, position, Performance, arousal_params);

% Initialize output structure
Fig1_Data = struct();
Fig1_Data.pupil_percentiles = pupil_percentiles;
Fig1_Data.analysis_params = arousal_params;

% Initialize data arrays based on analysis type
if strcmp(data_type, 'activity')
    Fig1_Data.HighArousalActivity = [];
    Fig1_Data.LowArousalActivity = [];
elseif strcmp(data_type, 'entropy')
    Fig1_Data.HighArousalEntropy = [];
    Fig1_Data.LowArousalEntropy = [];
else
    error('data_type must be ''activity'' or ''entropy''');
end

Fig1_Data.HighArousalConditionLabels = {};
Fig1_Data.LowArousalConditionLabels = {};

% Extraction counters
total_high_frames = 0;
total_low_frames = 0;

% Process each condition
for c_idx = 1:length(conditions)
    condName = conditions{c_idx};
    fprintf('Processing percentile-based %s data for condition: %s\n', data_type, condName);
    
    % Get recordings to skip for this condition
    recordingsToSkip = [];
    if isfield(Skip, condName)
        recordingsToSkip = Skip.(condName);
    end
    if ~isempty(recordingsToSkip)
        fprintf('  Skipping recordings %s for %s\n', mat2str(recordingsToSkip), condName);
    end
    
    % Check if condition exists in required structures
    if ~isfield(Behav, condName)
        warning('Condition %s not found in Behav structure', condName);
        continue;
    end
    
    % Get condition data based on analysis type
    if strcmp(data_type, 'activity')
        if ~isfield(Data, condName)
            warning('Condition %s not found in ActivityData structure', condName);
            continue;
        end
        condAnalysisData = Data.(condName);
    else % entropy
        if ~isfield(Data, condName) % Data parameter actually contains Entropy for entropy analysis
            warning('Condition %s not found in Entropy structure', condName);
            continue;
        end
        condAnalysisData = Data.(condName); % This is actually Entropy structure
    end
    
    condBehavData = Behav.(condName);
    
    condition_high_frames = 0;
    condition_low_frames = 0;
    
    % Process each recording
    for rec = 1:length(condBehavData)
        if rec > length(condBehavData) 
            continue;
        end
        if ismember(rec, recordingsToSkip)
            continue;
        end
        
        % Check required fields
        if ~isfield(condBehavData(rec), 'allPupilSize') || isempty(condBehavData(rec).allPupilSize)
            continue;
        end
        
        % Get data based on analysis type
        if strcmp(data_type, 'activity')
            if ~isfield(condAnalysisData(rec), 'all') || isempty(condAnalysisData(rec).all)
                continue;
            end
            analysisData = condAnalysisData(rec).all; % [neurons x frames x trials]
        else % entropy
            if ~isfield(condAnalysisData, 'AllNeurons') || ~isfield(condAnalysisData.AllNeurons, 'RecordingEntropy_Raw')
                continue;
            end
            entropyData = condAnalysisData.AllNeurons.RecordingEntropy_Raw;
    
            analysisData = entropyData{rec,1}; % [trials x frames]
        end
        
        pupilData = condBehavData(rec).allPupilSize / 1000; % [trials x frames]
        
        % Get position-specific trials
        validTrials = getValidTrials(position, Performance, condName, rec, size(pupilData, 1));
        
        fprintf('  Recording %d: Processing %d valid trials...\n', rec, length(validTrials));
        
        % Process each trial
        for trial = validTrials
            if trial > size(pupilData, 1)
                continue;
            end
            
            trialPupil = pupilData(trial, :);
            
            % Get trial analysis data based on type
            if strcmp(data_type, 'activity')
                if trial > size(analysisData, 3)
                    continue;
                end
                trialAnalysisData = squeeze(analysisData(:, :, trial)); % [neurons x frames]
            else % entropy
                if trial > size(analysisData, 1)
                    continue;
                end
                trialAnalysisData = analysisData(trial, :); % [1 x frames]
            end
            
            % Align data lengths
            lenFrames = min(length(trialPupil), size(trialAnalysisData, 2));
            if lenFrames < 25
                continue; % Skip very short trials
            end
            
            trialPupil = trialPupil(1:lenFrames);
            if strcmp(data_type, 'activity')
                trialAnalysisData = trialAnalysisData(:, 1:lenFrames);
            else
                trialAnalysisData = trialAnalysisData(1:lenFrames);
            end
            
            % Identify high and low arousal frames
            high_arousal_frames = find(trialPupil > pupil_percentiles.high_threshold);
            low_arousal_frames = find(trialPupil < pupil_percentiles.low_threshold);
            
            % Extract corresponding analysis data
            if ~isempty(high_arousal_frames)
                if strcmp(data_type, 'activity')
                    high_arousal_data = mean(trialAnalysisData(:, high_arousal_frames), 1); % Average across neurons
                    Fig1_Data.HighArousalActivity = [Fig1_Data.HighArousalActivity, high_arousal_data];
                else
                    high_arousal_data = trialAnalysisData(high_arousal_frames);
                    Fig1_Data.HighArousalEntropy = [Fig1_Data.HighArousalEntropy, high_arousal_data];
                end
                
                % Add condition labels
                for i = 1:length(high_arousal_frames)
                    Fig1_Data.HighArousalConditionLabels{end+1} = condName;
                end
                
                condition_high_frames = condition_high_frames + length(high_arousal_frames);
            end
            
            if ~isempty(low_arousal_frames)
                if strcmp(data_type, 'activity')
                    low_arousal_data = mean(trialAnalysisData(:, low_arousal_frames), 1); % Average across neurons
                    Fig1_Data.LowArousalActivity = [Fig1_Data.LowArousalActivity, low_arousal_data];
                else
                    low_arousal_data = trialAnalysisData(low_arousal_frames);
                    Fig1_Data.LowArousalEntropy = [Fig1_Data.LowArousalEntropy, low_arousal_data];
                end
                
                % Add condition labels
                for i = 1:length(low_arousal_frames)
                    Fig1_Data.LowArousalConditionLabels{end+1} = condName;
                end
                
                condition_low_frames = condition_low_frames + length(low_arousal_frames);
            end
        end
    end
    
    fprintf('  Added %d high arousal and %d low arousal frames for %s\n', ...
            condition_high_frames, condition_low_frames, condName);
    total_high_frames = total_high_frames + condition_high_frames;
    total_low_frames = total_low_frames + condition_low_frames;
end

% Convert to column vectors for consistency
if strcmp(data_type, 'activity')
    Fig1_Data.HighArousalActivity = Fig1_Data.HighArousalActivity(:);
    Fig1_Data.LowArousalActivity = Fig1_Data.LowArousalActivity(:);
else
    Fig1_Data.HighArousalEntropy = Fig1_Data.HighArousalEntropy(:);
    Fig1_Data.LowArousalEntropy = Fig1_Data.LowArousalEntropy(:);
end

fprintf('  Total percentile-based arousal data extracted:\n');
fprintf('    High arousal frames: %d\n', total_high_frames);
fprintf('    Low arousal frames: %d\n', total_low_frames);
fprintf('    Percentile thresholds: High=%.3f mm, Low=%.3f mm\n', ...
        pupil_percentiles.high_threshold, pupil_percentiles.low_threshold);

% Validate minimum data requirements
if total_high_frames < arousal_params.min_frames_per_bin
    warning('High arousal bin has insufficient data (%d < %d required)', ...
            total_high_frames, arousal_params.min_frames_per_bin);
end

if total_low_frames < arousal_params.min_frames_per_bin
    warning('Low arousal bin has insufficient data (%d < %d required)', ...
            total_low_frames, arousal_params.min_frames_per_bin);
end

end

%% ======================== HELPER FUNCTIONS ========================

function [validTrials] = getValidTrials(position, Performance, condName, rec, nTrials)
    % Get position-specific trials if Performance data is available
    
    validTrials = [];
    
    if exist('Performance', 'var') && ~isempty(Performance) && isfield(Performance, condName) && ...
       length(Performance.(condName)) >= rec && ~strcmp(position, 'All')
        
        % Get position-specific trials
        if strcmp(position, 'P1')
            hitField = 'hitP1'; missField = 'missP1';
        elseif strcmp(position, 'P3')
            hitField = 'hitP3'; missField = 'missP3';
        else
            hitField = 'hit'; missField = 'miss';
        end
        
        if isfield(Performance.(condName)(rec), hitField) && isfield(Performance.(condName)(rec), missField)
            hitTrials = Performance.(condName)(rec).(hitField);
            missTrials = Performance.(condName)(rec).(missField);
            validTrials = [hitTrials, missTrials];
        else
            validTrials = 1:nTrials;
        end
    else
        % Use all trials for 'All' position or when Performance data not available
        validTrials = 1:nTrials;
    end
end