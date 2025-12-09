function [pupil_percentiles] = calculate_pupil_percentiles(conditions, Behav, Skip, position, Performance, arousal_params)
% CALCULATE_PUPIL_PERCENTILES - Calculate percentile thresholds for pupil size binning
%
% This function calculates pupil size percentiles for arousal state analysis,
% supporting different pooling levels (recording, condition, trial) to account
% for baseline variations across sessions.
%
% INPUTS:
%   conditions      - Cell array of condition names to analyze
%   Behav          - Structure containing behavioral data including pupil size
%   Skip           - Structure containing recording indices to skip for each condition
%   position       - Position string ('All', 'P1', 'P3')
%   Performance    - Performance metrics structure for trial selection
%   arousal_params - Structure containing percentile calculation parameters:
%                   .pooling_level: 'recording' (default), 'condition', or 'trial'
%                   .high_arousal_percentile: percentile for high arousal (default 90)
%                   .low_arousal_percentile: percentile for low arousal (default 10)
%
% OUTPUTS:
%   pupil_percentiles - Structure containing percentile thresholds:
%                      .high_threshold: threshold for high arousal state
%                      .low_threshold: threshold for low arousal state
%                      .all_pupil_sizes: collected pupil data used for calculation
%                      .pooling_level: pooling method used
%                      .n_recordings: number of recordings processed
%                      .n_frames: total number of frames processed

% Set default parameters
if ~isfield(arousal_params, 'pooling_level')
    arousal_params.pooling_level = 'recording'; % Default: recording-level percentiles
end

if ~isfield(arousal_params, 'high_arousal_percentile')
    arousal_params.high_arousal_percentile = 90; % Top 10%
end

if ~isfield(arousal_params, 'low_arousal_percentile')
    arousal_params.low_arousal_percentile = 10; % Bottom 10%
end

fprintf('=== Calculating Pupil Percentiles for %s ===\n', position);
fprintf('  Pooling level: %s\n', arousal_params.pooling_level);
fprintf('  High arousal percentile: %.1f%%\n', arousal_params.high_arousal_percentile);
fprintf('  Low arousal percentile: %.1f%%\n', arousal_params.low_arousal_percentile);

% Initialize output structure
pupil_percentiles = struct();
pupil_percentiles.pooling_level = arousal_params.pooling_level;
pupil_percentiles.n_recordings = 0;
pupil_percentiles.n_frames = 0;

switch arousal_params.pooling_level
    case 'recording'
        % Calculate percentiles per recording (DEFAULT - accounts for session-specific baselines)
        pupil_percentiles = calculateRecordingLevelPercentiles(conditions, Behav, Skip, position, Performance, arousal_params, pupil_percentiles);
        
    case 'condition'
        % Pool across all recordings/trials within each condition
        pupil_percentiles = calculateConditionLevelPercentiles(conditions, Behav, Skip, position, Performance, arousal_params, pupil_percentiles);
        
    case 'trial'
        % Calculate percentiles per trial (most conservative approach)
        pupil_percentiles = calculateTrialLevelPercentiles(conditions, Behav, Skip, position, Performance, arousal_params, pupil_percentiles);
        
    otherwise
        error('Invalid pooling_level: %s. Must be ''recording'', ''condition'', or ''trial''', arousal_params.pooling_level);
end

% Validate percentile separation
if pupil_percentiles.high_threshold <= pupil_percentiles.low_threshold
    warning('High arousal threshold (%.3f) is not greater than low arousal threshold (%.3f)', ...
            pupil_percentiles.high_threshold, pupil_percentiles.low_threshold);
end

fprintf('  Percentile calculation complete:\n');
fprintf('    High arousal threshold (%.1f%%): %.3f mm\n', arousal_params.high_arousal_percentile, pupil_percentiles.high_threshold);
fprintf('    Low arousal threshold (%.1f%%): %.3f mm\n', arousal_params.low_arousal_percentile, pupil_percentiles.low_threshold);
fprintf('    Threshold separation: %.3f mm\n', pupil_percentiles.high_threshold - pupil_percentiles.low_threshold);
fprintf('    Total recordings processed: %d\n', pupil_percentiles.n_recordings);
fprintf('    Total frames processed: %d\n', pupil_percentiles.n_frames);

end

%% ======================== HELPER FUNCTIONS ========================

function [pupil_percentiles] = calculateRecordingLevelPercentiles(conditions, Behav, Skip, position, Performance, arousal_params, pupil_percentiles)
    % Calculate percentiles per recording (accounts for session-specific baselines)
    
    fprintf('  → Using recording-level percentiles (accounts for session-specific baselines)\n');
    
    all_high_thresholds = [];
    all_low_thresholds = [];
    
    % Process each condition
    for c_idx = 1:length(conditions)
        condName = conditions{c_idx};
        
        % Get recordings to skip for this condition
        recordingsToSkip = [];
        if isfield(Skip, condName)
            recordingsToSkip = Skip.(condName);
        end
        
        % Check if condition exists
        if ~isfield(Behav, condName)
            warning('Condition %s not found in Behav structure', condName);
            continue;
        end
        
        condBehavData = Behav.(condName);
        
        % Process each recording in this condition
        for rec = 1:length(condBehavData)
            if ismember(rec, recordingsToSkip)
                continue;
            end
            
            % Check required fields
            if ~isfield(condBehavData(rec), 'allPupilSize') || isempty(condBehavData(rec).allPupilSize)
                continue;
            end
            
            pupilData = condBehavData(rec).allPupilSize / 1000; % Convert to mm [trials x frames]
            
            % Get position-specific trials
            validTrials = getValidTrials(position, Performance, condName, rec, size(pupilData, 1));
            
            % Collect pupil data for this recording
            recording_pupil_data = [];
            for trial = validTrials
                if trial > size(pupilData, 1)
                    continue;
                end
                
                trialPupil = pupilData(trial, :);
                if length(trialPupil) < 10
                    continue; % Skip very short trials
                end
                
                recording_pupil_data = [recording_pupil_data; trialPupil(:)];
            end
            
            % Calculate percentiles for this recording
            if length(recording_pupil_data) >= 50 % Minimum data requirement
                high_threshold = prctile(recording_pupil_data, arousal_params.high_arousal_percentile);
                low_threshold = prctile(recording_pupil_data, arousal_params.low_arousal_percentile);
                
                all_high_thresholds = [all_high_thresholds; high_threshold];
                all_low_thresholds = [all_low_thresholds; low_threshold];
                
                pupil_percentiles.n_recordings = pupil_percentiles.n_recordings + 1;
                pupil_percentiles.n_frames = pupil_percentiles.n_frames + length(recording_pupil_data);
            end
        end
    end
    
    % Average thresholds across recordings
    if ~isempty(all_high_thresholds)
        pupil_percentiles.high_threshold = mean(all_high_thresholds);
        pupil_percentiles.low_threshold = mean(all_low_thresholds);
        pupil_percentiles.all_pupil_sizes = []; % Not stored for recording-level percentiles
    else
        error('No valid recordings found for percentile calculation');
    end
end

function [pupil_percentiles] = calculateConditionLevelPercentiles(conditions, Behav, Skip, position, Performance, arousal_params, pupil_percentiles)
    % Pool across all recordings within each condition
    
    fprintf('  → Using condition-level percentiles (may confound baseline differences)\n');
    
    all_pupil_sizes = [];
    
    % Process each condition
    for c_idx = 1:length(conditions)
        condName = conditions{c_idx};
        
        % Get recordings to skip for this condition
        recordingsToSkip = [];
        if isfield(Skip, condName)
            recordingsToSkip = Skip.(condName);
        end
        
        % Check if condition exists
        if ~isfield(Behav, condName)
            warning('Condition %s not found in Behav structure', condName);
            continue;
        end
        
        condBehavData = Behav.(condName);
        
        % Process each recording in this condition
        for rec = 1:length(condBehavData)
            if ismember(rec, recordingsToSkip)
                continue;
            end
            
            % Check required fields
            if ~isfield(condBehavData(rec), 'allPupilSize') || isempty(condBehavData(rec).allPupilSize)
                continue;
            end
            
            pupilData = condBehavData(rec).allPupilSize / 1000; % Convert to mm [trials x frames]
            pupil_percentiles.n_recordings = pupil_percentiles.n_recordings + 1;
            
            % Get position-specific trials
            validTrials = getValidTrials(position, Performance, condName, rec, size(pupilData, 1));
            
            % Collect pupil data for this recording
            for trial = validTrials
                if trial > size(pupilData, 1)
                    continue;
                end
                
                trialPupil = pupilData(trial, :);
                if length(trialPupil) < 10
                    continue; % Skip very short trials
                end
                
                all_pupil_sizes = [all_pupil_sizes; trialPupil(:)];
                pupil_percentiles.n_frames = pupil_percentiles.n_frames + length(trialPupil);
            end
        end
    end
    
    % Calculate percentiles across all conditions
    if length(all_pupil_sizes) >= 100 % Minimum data requirement
        pupil_percentiles.high_threshold = prctile(all_pupil_sizes, arousal_params.high_arousal_percentile);
        pupil_percentiles.low_threshold = prctile(all_pupil_sizes, arousal_params.low_arousal_percentile);
        pupil_percentiles.all_pupil_sizes = all_pupil_sizes;
    else
        error('Insufficient data for percentile calculation (need at least 100 frames, got %d)', length(all_pupil_sizes));
    end
end

function [pupil_percentiles] = calculateTrialLevelPercentiles(conditions, Behav, Skip, position, Performance, arousal_params, pupil_percentiles)
    % Calculate percentiles per trial (most conservative approach)
    
    fprintf('  → Using trial-level percentiles (most conservative approach)\n');
    
    % This approach would calculate percentiles for each trial individually
    % For now, implement as average of recording-level percentiles
    % (can be expanded to store per-trial thresholds if needed)
    
    pupil_percentiles = calculateRecordingLevelPercentiles(conditions, Behav, Skip, position, Performance, arousal_params, pupil_percentiles);
    
    fprintf('  → Note: Trial-level percentiles currently implemented as recording-level average\n');
end

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