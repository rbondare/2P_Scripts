function LocomotionCorrelation_Entropy(Entropy, WindowLength, params, Rec, Performance, LocomotionCal, conditions, Skip)
% LOCOMOTIONCORRELATION_ENTROPY - Analyze entropy around locomotion periods
% across different behavioral conditions. Uses trial-based locomotion analysis.
%
% INPUTS:
%   Entropy      - Structure containing entropy data for different conditions
%   WindowLength - Structure containing window length information
%   params       - Structure containing plotting parameters including colors
%   Rec          - Recording metadata structure
%   Performance  - Performance metrics structure containing trial data
%   LocomotionCal- Structure containing calibrated locomotion data for different conditions
%   conditions   - Cell array of condition names to analyze (e.g., {'Expert', 'Beginner'})
%   Skip         - (optional) Structure containing recording indices to skip for each condition
%
% OUTPUTS:
%   For each position (All, P1, P3):
%   • Figure 1: Pre vs Post entropy scatter around locomotion events – Wilcoxon signed‑rank (paired)
%   • Figure 2: Pre vs Post entropy violin plot
%   • Figure 3: Histogram of Post/Pre entropy ratios – sign‑rank against median of 1
%   • Figure 4: Pre value vs Change scatter (entropy) – Spearman ρ with p‑value
%   • Figure 5: Locomotion triggered plots (Population Entropy time series)
%   • Figure 6: Locomotion vs entropy frame-by-frame scatter plot with correlation analysis

%% FramRate Info
LocomotionFrameRate = 50; %Hz

threshold = 2.0; % cm/s

%% Input validation
if nargin < 8 || isempty(Skip)
    Skip = struct(); 
end

if nargin < 7 || isempty(conditions)
    conditions = {};
    if isfield(Entropy, 'Expert'); conditions{end+1} = 'Expert'; end
    if isfield(Entropy, 'Beginner'); conditions{end+1} = 'Beginner'; end
    if isfield(Entropy, 'Naive'); conditions{end+1} = 'Naive'; end
    if isfield(Entropy, 'NoSpout'); conditions{end+1} = 'NoSpout'; end
end

% positions = {'All', 'P1', 'P3'};
positions = {'All'};


fprintf('Processing locomotion correlations for ENTROPY in conditions: %s\n', strjoin(conditions, ', '));
fprintf('Positions: %s\n', strjoin(positions, ', '));

% Print Skip information if provided
if ~isempty(fieldnames(Skip))
    fprintf('Skip information provided for conditions: %s\n', strjoin(fieldnames(Skip), ', '));
    for skip_field = fieldnames(Skip)'
        fprintf('  %s: skipping recordings %s\n', skip_field{1}, mat2str(Skip.(skip_field{1})));
    end
else
    fprintf('No recordings will be skipped (Skip structure is empty)\n');
end

%% Loop through each position with "Data First, Plots Second" approach
for pos_idx = 1:length(positions)
    position = positions{pos_idx};
    fprintf('\n\n======== PROCESSING POSITION: %s ========\n', position);
    
    %% ==================== PHASE 1: DATA EXTRACTION ====================
    fprintf('\n*** PHASE 1: DATA EXTRACTION ***\n');
    
    % Step 1: Extract locomotion event data
    fprintf('\n--- Step 1: Locomotion Event Entropy Data Extraction ---\n');
    Fig1_Data = extractLocomotionEventEntropyData(conditions, Entropy, LocomotionCal, Skip, position, Performance, LocomotionFrameRate, threshold);
    
    % Step 2: Calculate entropy ratios
    fprintf('\n--- Step 2: Entropy Ratio Calculation ---\n');
    Fig3_Data = calculateLocomotionEntropyRatios(Fig1_Data);
    
    % Step 3: Calculate entropy changes
    fprintf('\n--- Step 3: Entropy Change Calculation ---\n');
    Fig4_Data = calculateLocomotionEntropyChanges(Fig1_Data);
    
    % Step 4: Extract triggered data for peri-event plots
    fprintf('\n--- Step 4: Locomotion Triggered Data Extraction ---\n');
    Fig5_Data = extractLocomotionTriggeredEntropyData(conditions, Entropy, LocomotionCal, Skip, position, Performance,LocomotionFrameRate, threshold);
    
    % Step 5: Extract frame-by-frame locomotion vs entropy data
    fprintf('\n--- Step 5: Frame-by-Frame Locomotion vs Entropy Data Extraction ---\n');
    Fig6_Data = extractFrameByFrameEntropyData(conditions, Entropy, LocomotionCal, Skip, position, Performance,LocomotionFrameRate, threshold);
    
    % Step 6: Validate data consistency
    fprintf('\n--- Step 6: Data Validation ---\n');
    validateLocomotionEntropyDataConsistency(Fig1_Data, Fig3_Data, Fig4_Data);
    
    %% ==================== PHASE 2: PLOT GENERATION ====================
    fprintf('\n*** PHASE 2: PLOT GENERATION ***\n');
    
    % Check if we have any data
    if isempty(Fig1_Data.PreLocomotionEntropy)
        warning('No locomotion event data found for position %s', position);
        continue;
    end
    
    fprintf('Data extraction complete for %s:\n', position);
    fprintf('  Total locomotion events: %d\n', length(Fig1_Data.PreLocomotionEntropy));
    fprintf('  Conditions: %s\n', strjoin(unique(Fig1_Data.AllConditionLabels), ', '));
    
    %% Print condition breakdown
    unique_conditions = unique(Fig1_Data.AllConditionLabels);
    for i = 1:length(unique_conditions)
        count = sum(strcmp(Fig1_Data.AllConditionLabels, unique_conditions{i}));
        fprintf('  %s: %d events (%.1f%%)\n', unique_conditions{i}, count, 100*count/length(Fig1_Data.AllConditionLabels));
    end

    %% ------------------------- Figure 1: Entropy scatter -------------------- %%
    fprintf('\n--- Generating Figure 1: Pre vs Post Entropy Scatter ---\n');
    
    scatterColormap = gray;
    scatterMarkerSize = 30;
    
    figure('Color', 'w', 'Name', sprintf('Pre vs Post Locomotion Entropy Scatter - %s', position)); 
    hold on;
    dscatter(Fig1_Data.PreLocomotionEntropy, Fig1_Data.PostLocomotionEntropy, 'MSIZE', scatterMarkerSize);
    colormap(scatterColormap);
    axis square;
    % Set equal axis limits
    allData = [Fig1_Data.PreLocomotionEntropy; Fig1_Data.PostLocomotionEntropy];
    dataLim = [min(allData), max(allData)];
    xlim(dataLim);
    % Set nice axis ticks at half and full numbers
    setNiceAxisTicks(gca, dataLim);
    ylim(dataLim);
    ref = refline(1); 
    ref.LineWidth = 5; 
    ref.Color = 'r';
    
    % Add fitted trend line
    fit = lsline(gca);
    fit.LineStyle = '--';
    fit.LineWidth = 2;
    fit.Color = [0.5 0.5 0.5];
    
    [pEnt, ~, statsEnt] = signrank(Fig1_Data.PreLocomotionEntropy, Fig1_Data.PostLocomotionEntropy);
    ratioEnt = mean(Fig1_Data.PreLocomotionEntropy ./ Fig1_Data.PostLocomotionEntropy, 'omitnan');
    [rhoEnt, pRhoEnt] = corr(Fig1_Data.PreLocomotionEntropy, Fig1_Data.PostLocomotionEntropy, 'Type', 'Spearman', 'Rows', 'complete');
    
    ax = gca;
    ax.FontSize = 20;
    xlabel('Pre-Locomotion Entropy');
    ylabel('Post-Locomotion Entropy');
    
    title({sprintf('Pre vs Post Locomotion Entropy - %s', position), sprintf('ratio = %.3f,  p_{signrank} = %.3g,  rho = %.3f,  p_{rho} = %.3g', ratioEnt, pEnt, rhoEnt, pRhoEnt), ...
           sprintf('n = %d events from %d conditions', length(Fig1_Data.PreLocomotionEntropy), length(conditions))});
    
    %% ------------------------- Figure 2: Violin plot -------------------- %%
    fprintf('\n--- Generating Figure 2: Pre vs Post Entropy Violin Plot ---\n');
    
    % Calculate effect size (Hedges' g)
    pre_data = Fig1_Data.PreLocomotionEntropy;
    post_data = Fig1_Data.PostLocomotionEntropy;
    
    mean_pre = mean(pre_data);
    mean_post = mean(post_data);
    std_pre = std(pre_data);
    std_post = std(post_data);
    n_pre = length(pre_data);
    n_post = length(post_data);
    
    % Pooled standard deviation
    pooled_std = sqrt(((n_pre - 1) * std_pre^2 + (n_post - 1) * std_post^2) / (n_pre + n_post - 2));
    
    % Cohen's d
    cohens_d = (mean_post - mean_pre) / pooled_std;
    
    % Hedges' g (bias correction for Cohen's d)
    correction_factor = 1 - (3 / (4 * (n_pre + n_post - 2) - 1));
    hedges_g = cohens_d * correction_factor;
    
    % Effect size interpretation
    if abs(hedges_g) < 0.2
        effect_interp = 'negligible';
    elseif abs(hedges_g) < 0.5
        effect_interp = 'small';
    elseif abs(hedges_g) < 0.8
        effect_interp = 'medium';
    else
        effect_interp = 'large';
    end
    
    figure('Color', 'w', 'Name', sprintf('Pre vs Post Locomotion Entropy Violin Plot - %s', position));
    axesFont = 24;
    daviolinplot({[Fig1_Data.PreLocomotionEntropy, Fig1_Data.PostLocomotionEntropy]}, ...
                 'violin', 'full', 'groups', ones(numel(Fig1_Data.PreLocomotionEntropy), 1), ...
                 'colors', params.colour(13, :), 'box', 3, 'boxcolor', 'w', ...
                 'linkline', 0, 'withinlines', 0, ...
                 'xtlabels', {"Before Locomotion Event", "After Locomotion Event"});
    set(gca, 'FontSize', axesFont);
    ylabel('Entropy');
    title(gca, {sprintf('Locomotion Entropy %s: p = %.3g, g = %.3f (%s)', position, pEnt, hedges_g, effect_interp), ...
                sprintf('n = %d', length(Fig1_Data.PreLocomotionEntropy))});
    
    fprintf('Effect size: Hedges'' g = %.3f (%s effect)\n', hedges_g, effect_interp);
    
    %% ------------------------- Figure 3: Ratio Histogram --------------- %%
    fprintf('\n--- Generating Figure 3: Post/Pre Entropy Ratio Histogram ---\n');
    
    figure('Color', 'w', 'Name', sprintf('Post/Pre Locomotion Entropy Ratio Histogram - %s', position));
    
    histogram(Fig3_Data.EntropyRatio, 'Normalization', 'probability', ...
              'BinMethod', 'sturges', 'FaceColor', params.colour(9, :), ...
              'EdgeColor', 'none', 'FaceAlpha', 0.7);
    [pRatioEnt] = signrank(Fig3_Data.EntropyRatio - 1);
    ax = gca;
    ax.FontSize = axesFont;
    ax.Box = 'on';
    ax.XLabel.String = 'Post / Pre Entropy';
    ax.YLabel.String = 'Probability';
    
    title({sprintf('Locomotion Entropy Ratio - %s', position), sprintf('p = %.3g', pRatioEnt), sprintf('n = %d', length(Fig3_Data.EntropyRatio))});
    
    %% ------------------------- Figure 4: Change vs Pre ------------------ %%
    fprintf('\n--- Generating Figure 4: Entropy Change vs Pre-Event Values ---\n');
    
    figure('Color', 'w', 'Name', sprintf('Locomotion Entropy Change vs Pre-Event Values - %s', position));
    
    dscatter(Fig1_Data.PreLocomotionEntropy, Fig4_Data.EntropyChange, 'MSIZE', scatterMarkerSize);
    colormap(scatterColormap);
    % Add fitted trend line
    fit = lsline(gca);
    fit.LineStyle = '--';
    fit.LineWidth = 2;
    fit.Color = [0.5 0.5 0.5];
    axis square; 
    grid on;
    [rhoEnt, pRhoEnt] = corr(Fig1_Data.PreLocomotionEntropy, Fig4_Data.EntropyChange, ...
                             'Type', 'Spearman', 'Rows', 'complete');
    ax = gca;
    ax.FontSize = axesFont;
    ax.XLabel.String = 'Pre Locomotion Entropy';
    ax.YLabel.String = '\Delta Entropy (Post - Pre)';
    
    title({sprintf('Locomotion Entropy Change vs Pre - %s', position), sprintf('rho = %.2f,  p = %.3g', rhoEnt, pRhoEnt), ...
           sprintf('n = %d', length(Fig1_Data.PreLocomotionEntropy))});
    
    %% ------------------------- Figure 5: Triggered Plots --------------- %%
    fprintf('\n--- Generating Figure 5: Locomotion Triggered Plot ---\n');
    
    if ~isempty(Fig5_Data.TriggeredLocomotion) && ~isempty(Fig5_Data.TriggeredEntropy)
        
        % Create subplot figure for triggered data
        figure('Color', 'w', 'Name', sprintf('Locomotion-Triggered Plot - %s', position));

        timeVector = Fig5_Data.TimeVector;
      
        % Triggered Population Entropy
        meanEntropy = nanmean(Fig5_Data.TriggeredEntropy, 1);
        stdEntropy = nanstd(Fig5_Data.TriggeredEntropy, [], 1);
        
        plot(timeVector, meanEntropy, 'LineWidth', 3, 'Color', params.colour(9,:));
        hold on;
        fill([timeVector, fliplr(timeVector)], [meanEntropy + stdEntropy, fliplr(meanEntropy - stdEntropy)], ...
             params.colour(9,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        
        xline(0, '--k', 'LineWidth', 2);
        xlabel('Time from Locomotion Event (s)');
        ylabel('Population Entropy');
        title(sprintf('Population Entropy Triggered by Locomotion Events - %s (n=%d events)', position, size(Fig5_Data.TriggeredLocomotion,1)));
        grid on;
        ax = gca; ax.FontSize = 16;
        
        fprintf('Generated locomotion triggered plots with %d events\n', size(Fig5_Data.TriggeredLocomotion,1));
    else
        fprintf('No triggered data available for plotting\n');
    end
    
    %% ------------------------- Figure 6: Locomotion vs Entropy Scatter ------ %%
    fprintf('\n--- Generating Figure 6: Locomotion vs Entropy Scatter Plot ---\n');
    
    if ~isempty(Fig6_Data.AllLocomotionMagnitudes) && ~isempty(Fig6_Data.AllEntropyValues)
        
        % Validate data before plotting
        validLocomotionIdx = ~isnan(Fig6_Data.AllLocomotionMagnitudes) & ~isinf(Fig6_Data.AllLocomotionMagnitudes) & (Fig6_Data.AllLocomotionMagnitudes >= 0);
        validEntropyIdx = ~isnan(Fig6_Data.AllEntropyValues) & ~isinf(Fig6_Data.AllEntropyValues);
        validIdx = validLocomotionIdx & validEntropyIdx;
        
        if sum(validIdx) > 10
            figure('Color', 'w', 'Name', sprintf('Locomotion vs Entropy Scatter - %s', position));
            
            dscatter(Fig6_Data.AllLocomotionMagnitudes(validIdx), Fig6_Data.AllEntropyValues(validIdx));
            colormap(scatterColormap);
            
            % Calculate correlation
            [rho, pval] = corr(Fig6_Data.AllLocomotionMagnitudes(validIdx), Fig6_Data.AllEntropyValues(validIdx), 'Type', 'Spearman');
            
            % Add trend line
            hold on;
            % Add fitted trend line
    fit = lsline(gca);
    fit.LineStyle = '--';
    fit.LineWidth = 2;
    fit.Color = [0.5 0.5 0.5];
            
            ax = gca;
            ax.FontSize = 16;
            xlabel('Locomotion Speed (cm/s)');
            ylabel('Population Entropy');
            title({sprintf('Locomotion vs Population Entropy - %s', position), ...
                   sprintf('Spearman rho = %.3f, p = %.3g', rho, pval), ...
                   sprintf('n = %d frames from %d conditions', sum(validIdx), length(conditions))});
            grid on;
            
            fprintf('Generated locomotion vs entropy scatter with %d valid data points\n', sum(validIdx));
            fprintf('Spearman correlation: rho = %.3f, p = %.3g\n', rho, pval);
        else
            fprintf('Insufficient valid data points for correlation analysis (%d valid points)\n', sum(validIdx));
            if sum(~validLocomotionIdx) > 0
                fprintf('  Invalid locomotion data: %d NaN/inf/negative values\n', sum(~validLocomotionIdx));
            end
            if sum(~validEntropyIdx) > 0
                fprintf('  Invalid entropy data: %d NaN/inf values\n', sum(~validEntropyIdx));
            end
        end
    else
        fprintf('No frame-by-frame data available for plotting\n');
    end
    
    %% Summary output for this position
    fprintf('\nLocomotionCorrelation_Entropy completed for position %s!\n', position);
    fprintf('Processed %d conditions: %s\n', length(conditions), strjoin(conditions, ', '));
    fprintf('Total locomotion events analyzed: %d\n', length(Fig1_Data.PreLocomotionEntropy));

end

fprintf('\n\n======== ALL POSITIONS COMPLETED ========\n');
fprintf('Analyzed positions: %s\n', strjoin(positions, ', '));

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

%% ======================== HELPER FUNCTIONS ========================

function [Fig1_Data] = extractLocomotionEventEntropyData(conditions, Entropy, LocomotionCal, Skip, position, Performance, LocomotionFrameRate, threshold)
    % EXTRACTLOCOMOTIONEVENTNENTROPYDATA - Extract locomotion event data for entropy analysis
    
    fprintf('=== Extracting Locomotion Event Entropy Data for %s ===\n', position);
    
    % Initialize output structure
    Fig1_Data = struct();
    Fig1_Data.PreLocomotionEntropy = [];
    Fig1_Data.PostLocomotionEntropy = [];
    Fig1_Data.AllConditionLabels = {};
    
    % Process each condition
    for c_idx = 1:length(conditions)
        condName = conditions{c_idx};
        fprintf('Processing locomotion event entropy data for condition: %s\n', condName);
        
        % Get recordings to skip for this condition
        recordingsToSkip = [];
        if isfield(Skip, condName); recordingsToSkip = Skip.(condName); end
        if ~isempty(recordingsToSkip)
            fprintf('  Skipping recordings %s for %s\n', mat2str(recordingsToSkip), condName);
        end
        
        % Check if condition exists in all required structures
        if ~isfield(Entropy, condName) || ~isfield(LocomotionCal, condName)
            warning('Condition %s not found in one or more required data structures', condName);
            continue;
        end
        
        % Get condition data
        condEntropy = Entropy.(condName);
        condLocomotionData = LocomotionCal.(condName);
        
        locomotion_event_count_this_condition = 0;
        
        % Handle different data structures for different conditions
        if strcmp(condName, 'Expert') || strcmp(condName, 'Beginner')
            % Expert/Beginner: use the same structure as the working activity function
            if isfield(condEntropy, 'AllNeurons') && isfield(condEntropy.AllNeurons, 'RecordingEntropy_Raw')
                entropyData = condEntropy.AllNeurons.RecordingEntropy_Raw;
                
                for i = 1: length(entropyData)
                    if ismember(i, recordingsToSkip) || isempty(condLocomotionData(i).RunForward); continue; end
                    
                    if iscell(entropyData) && i <= length(entropyData) && ~isempty(entropyData{i, 1})
                        recEntropyData = entropyData{i, 1}; % [trials x frames]
                        
                        fprintf('  Processing recording %d...\n', i);
                        
                        % Get position-specific trials if Performance data is available
                        validTrials = [];
                        if exist('Performance', 'var') && ~isempty(Performance) && isfield(Performance, condName) && ...
                           length(Performance.(condName)) >= i && ~strcmp(position, 'All')
                            
                            % Get position-specific trials
                            if strcmp(position, 'P1')
                                hitField = 'hitP1'; missField = 'missP1';
                            elseif strcmp(position, 'P3')
                                hitField = 'hitP3'; missField = 'missP3';
                            else
                                hitField = 'hit'; missField = 'miss';
                            end
                            
                            if isfield(Performance.(condName)(i), hitField) && isfield(Performance.(condName)(i), missField)
                                hitTrials = Performance.(condName)(i).(hitField);
                                missTrials = Performance.(condName)(i).(missField);
                                validTrials = [hitTrials, missTrials];
                            else
                                validTrials = 1:min([size(recEntropyData, 1), 150]);
                            end
                        else
                            % Use all trials for 'All' position or when Performance data not available
                            validTrials = 1:min([size(recEntropyData, 1), 150]);
                        end
                        
                        fprintf('    Valid trials for this recording: %d\n', length(validTrials));
                        
                        % Process trials with locomotion data
                        for t = validTrials
                                locomotionTrace = condLocomotionData(i).RunForward(t,:)' * LocomotionFrameRate;
                                
                                high_movement = locomotionTrace > threshold;
                                onsets = find(diff([0; high_movement]) == 1);
                                
                                % Process locomotion events
                                for onset_idx = 1:length(onsets)
                                    onset_frame = onsets(onset_idx);
                                    neural_frame = round(onset_frame * size(recEntropyData, 2) / length(locomotionTrace));
                                    
                                    if neural_frame > 5 && neural_frame < (size(recEntropyData, 2) - 5)
                                        preEntropy = mean(recEntropyData(t, neural_frame-2:neural_frame-1), 2);
                                        postEntropy = mean(recEntropyData(t, neural_frame+1:neural_frame+2), 2);
                                        
                                        Fig1_Data.PreLocomotionEntropy = [Fig1_Data.PreLocomotionEntropy; preEntropy];
                                        Fig1_Data.PostLocomotionEntropy = [Fig1_Data.PostLocomotionEntropy; postEntropy];
                                        Fig1_Data.AllConditionLabels{end+1} = condName;
                                        locomotion_event_count_this_condition = locomotion_event_count_this_condition + 1;
                                    end
                                end
                            
                        end
                    end
                end
            end
            
        else
            % Naive/NoSpout: different structure 
            if isfield(condEntropy, 'AllNeurons') && isfield(condEntropy.AllNeurons, 'RecordingEntropy_Raw')
                entropyData = condEntropy.AllNeurons.RecordingEntropy_Raw;
                
                for a = 1:size(condLocomotionData, 2)
                    if ismember(a, recordingsToSkip); continue; end
                    
                    if isfield(condLocomotionData(a), 'Forward_mm') && ...
                       iscell(entropyData) && a <= length(entropyData) && ~isempty(entropyData{a, 1})
                        locomotionTrace = condLocomotionData(a).Forward_mm * LocomotionFrameRate ;
                        recEntropyData = entropyData{a, 1}; % [trials x frames]
                        
                        fprintf('  Processing recording %d...\n', a);
                        
                        % Trial segmentation
                        TrialLength = floor(length(locomotionTrace) / 150);
                        if TrialLength > 0
                            TrialMovement = 1:TrialLength:length(locomotionTrace);
                            
                            % Get position-specific trials
                            validTrialSegments = [];
                            if exist('Performance', 'var') && ~isempty(Performance) && isfield(Performance, condName) && ...
                               length(Performance.(condName)) >= a && ~strcmp(position, 'All')
                                
                                if strcmp(position, 'P1')
                                    hitField = 'hitP1'; missField = 'missP1';
                                elseif strcmp(position, 'P3')
                                    hitField = 'hitP3'; missField = 'missP3';
                                else
                                    hitField = 'hit'; missField = 'miss';
                                end
                                
                                if isfield(Performance.(condName)(a), hitField) && isfield(Performance.(condName)(a), missField)
                                    hitTrials = Performance.(condName)(a).(hitField);
                                    missTrials = Performance.(condName)(a).(missField);
                                    validTrialSegments = [hitTrials, missTrials];
                                else
                                    validTrialSegments = 1:min([150, length(TrialMovement)-1, size(recEntropyData, 1)]);
                                end
                            else
                                % Use all trial segments for 'All' position or when Performance data not available
                                validTrialSegments = 1:min([150, length(TrialMovement)-1, size(recEntropyData, 1)]);
                            end
                            
                            % Process trial segments
                            for t = validTrialSegments
                                
                                trial_start = TrialMovement(t);
                                trial_end = TrialMovement(t+1) - 1;
                                trial_locomotion = locomotionTrace(trial_start:trial_end) ;
                                
                                high_movement = trial_locomotion > threshold;
                                onsets = find(diff([0; high_movement]) == 1);
                                
                                for onset_idx = 1:length(onsets)
                                    onset_frame = onsets(onset_idx);
                                    neural_frame = round(onset_frame * size(recEntropyData, 2) / length(trial_locomotion));
                                    
                                    if neural_frame > 5 && neural_frame < (size(recEntropyData, 2) - 5) && t <= size(recEntropyData, 1)
                                        preEntropy = mean(recEntropyData(t, neural_frame-2:neural_frame-1), 2);
                                        postEntropy = mean(recEntropyData(t, neural_frame+1:neural_frame+2), 2);
                                        
                                        Fig1_Data.PreLocomotionEntropy = [Fig1_Data.PreLocomotionEntropy; preEntropy];
                                        Fig1_Data.PostLocomotionEntropy = [Fig1_Data.PostLocomotionEntropy; postEntropy];
                                        Fig1_Data.AllConditionLabels{end+1} = condName;
                                        locomotion_event_count_this_condition = locomotion_event_count_this_condition + 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        fprintf('  Added %d locomotion events for %s\n', locomotion_event_count_this_condition, condName);
    end
    
    fprintf('  Total locomotion events extracted: %d\n', length(Fig1_Data.PreLocomotionEntropy));
end

function [Fig3_Data] = calculateLocomotionEntropyRatios(Fig1_Data)
    % CALCULATELOCOMOTIONENTROPYRATIOS - Calculate Post/Pre ratios for Figure 3
    
    fprintf('=== Calculating Entropy Ratios ===\n');
    
    Fig3_Data = struct();
    
    if ~isempty(Fig1_Data.PreLocomotionEntropy)
        Fig3_Data.EntropyRatio = Fig1_Data.PostLocomotionEntropy ./ Fig1_Data.PreLocomotionEntropy;
        fprintf('  Entropy ratios calculated: %d\n', length(Fig3_Data.EntropyRatio));
    else
        Fig3_Data.EntropyRatio = [];
        fprintf('  No data available for ratio calculation\n');
    end
end

function [Fig4_Data] = calculateLocomotionEntropyChanges(Fig1_Data)
    % CALCULATELOCOMOTIONENTROPYCHANGES - Calculate Post-Pre changes for Figure 4
    
    fprintf('=== Calculating Entropy Changes ===\n');
    
    Fig4_Data = struct();
    
    if ~isempty(Fig1_Data.PreLocomotionEntropy)
        Fig4_Data.EntropyChange = Fig1_Data.PostLocomotionEntropy - Fig1_Data.PreLocomotionEntropy;
        fprintf('  Entropy changes calculated: %d\n', length(Fig4_Data.EntropyChange));
    else
        Fig4_Data.EntropyChange = [];
        fprintf('  No data available for change calculation\n');
    end
end

function [Fig5_Data] = extractLocomotionTriggeredEntropyData(conditions, Entropy, LocomotionCal, Skip, position, Performance, LocomotionFrameRate, threshold)
    % EXTRACTLOCOMOTIONTRIGGEREDENTROPYDATA - Extract locomotion-triggered time series data for peri-event plots
    
    fprintf('=== Extracting Locomotion Triggered Time Series Entropy Data for %s ===\n', position);
    
    % Calculate sampling rates for both data types
    locomotionSamplingRate =  size(LocomotionCal.Expert(2).RunForward,2)/18.5; % Hz 
    entropySamplingRate = size(Entropy.Expert(1).AllNeurons.Entropy_Raw,2)/18.5; % Hz
    
    % Parameters for triggered analysis 
    preWindowSec = 2;   % 2 seconds before event
    postWindowSec = 4;  % 4 seconds after event
    totalWindowSec = preWindowSec + postWindowSec;
    
    locomotionPreWindowFrames = round(preWindowSec * locomotionSamplingRate);
    locomotionPostWindowFrames = round(postWindowSec * locomotionSamplingRate);
    entropyPreWindowFrames = round(preWindowSec * entropySamplingRate);
    entropyPostWindowFrames = round(postWindowSec * entropySamplingRate);
    
    % Use entropy sampling rate for the time vector (since that's what we'll plot)
    totalEntropyFrames = entropyPreWindowFrames + entropyPostWindowFrames + 1;
    
    % Initialize output structure
    Fig5_Data = struct();
    Fig5_Data.TriggeredLocomotion = [];
    Fig5_Data.TriggeredEntropy = [];
    Fig5_Data.TimeVector = linspace(-preWindowSec, postWindowSec, totalEntropyFrames);
    
    fprintf('Locomotion sampling rate: %.1f Hz (pre: %d, post: %d frames)\n', locomotionSamplingRate, locomotionPreWindowFrames, locomotionPostWindowFrames);
    fprintf('Entropy sampling rate: %.1f Hz (pre: %d, post: %d frames)\n', entropySamplingRate, entropyPreWindowFrames, entropyPostWindowFrames);
    fprintf('Window: -%.1f to +%.1f seconds\n', preWindowSec, postWindowSec);
    fprintf('Total entropy frames per event: %d\n', totalEntropyFrames);
    
    % Process each condition
    for c_idx = 1:length(conditions)
        condName = conditions{c_idx};
        fprintf('Processing triggered entropy data for condition: %s\n', condName);
        
        % Get recordings to skip for this condition
        recordingsToSkip = [];
        if isfield(Skip, condName); recordingsToSkip = Skip.(condName); end
        if ~isempty(recordingsToSkip)
            fprintf('  Skipping recordings %s for %s\n', mat2str(recordingsToSkip), condName);
        end
        
        % Check if condition exists in all required structures
        if ~isfield(Entropy, condName) || ~isfield(LocomotionCal, condName)
            warning('Condition %s not found in one or more required data structures', condName);
            continue;
        end
        
        % Get condition data
        condEntropy = Entropy.(condName);
        condLocomotionData = LocomotionCal.(condName);
        
        triggered_event_count = 0;
        
        % Handle different data structures for different conditions
        if strcmp(condName, 'Expert') || strcmp(condName, 'Beginner')
            % Expert/Beginner: use the same structure as the working functions
            if isfield(condEntropy, 'AllNeurons') && isfield(condEntropy.AllNeurons, 'RecordingEntropy_Raw')
                entropyData = condEntropy.AllNeurons.RecordingEntropy_Raw;
                
                for i = 1:min(length(condLocomotionData), length(entropyData))
                    if ismember(i, recordingsToSkip) || isempty(condLocomotionData(i).RunForward); continue; end
                    
                    if iscell(entropyData) && i <= length(entropyData) && ~isempty(entropyData{i, 1})
                        recEntropyData = entropyData{i, 1}; % [trials x frames]
                        
                        % Get position-specific trials
                        validTrials = [];
                        if exist('Performance', 'var') && ~isempty(Performance) && isfield(Performance, condName) && ...
                           length(Performance.(condName)) >= i && ~strcmp(position, 'All')
                            
                            if strcmp(position, 'P1')
                                hitField = 'hitP1'; missField = 'missP1';
                            elseif strcmp(position, 'P3')
                                hitField = 'hitP3'; missField = 'missP3';
                            else
                                hitField = 'hit'; missField = 'miss';
                            end
                            
                            if isfield(Performance.(condName)(i), hitField) && isfield(Performance.(condName)(i), missField)
                                hitTrials = Performance.(condName)(i).(hitField);
                                missTrials = Performance.(condName)(i).(missField);
                                validTrials = [hitTrials, missTrials];
                            else
                                validTrials = 1:min([size(recEntropyData, 1), 150]);
                            end
                        else
                            validTrials = 1:min([size(recEntropyData, 1), 150]);
                        end
                        
                        % Process trials with locomotion data
                        for t = validTrials
                                locomotionTrace = condLocomotionData(i).RunForward(t,:)' * LocomotionFrameRate;
             
                                high_movement = locomotionTrace > threshold;
                                onsets = find(diff([0; high_movement]) == 1);
                                
                                % Process locomotion events
                                for onset_idx = 1:min(10, length(onsets))
                                    onset_frame = onsets(onset_idx);
                                    neural_frame = round(onset_frame * size(recEntropyData, 2) / length(locomotionTrace));
                                    
                                    if neural_frame > entropyPreWindowFrames && neural_frame <= (size(recEntropyData, 2) - entropyPostWindowFrames)
                                        % Extract locomotion data with padding if necessary
                                        locStart = onset_frame - locomotionPreWindowFrames;
                                        locEnd = onset_frame + locomotionPostWindowFrames;
                                        
                                        if locStart < 1 || locEnd > length(locomotionTrace)
                                            triggeredLocomotion = nan(1, locomotionPreWindowFrames + locomotionPostWindowFrames + 1);
                                            validStart = max(1, locStart);
                                            validEnd = min(length(locomotionTrace), locEnd);
                                            offsetStart = validStart - locStart + 1;
                                            offsetEnd = offsetStart + (validEnd - validStart);
                                            triggeredLocomotion(offsetStart:offsetEnd) = locomotionTrace(validStart:validEnd);
                                        else
                                            triggeredLocomotion = locomotionTrace(locStart:locEnd);
                                        end
                                        
                                        % Extract entropy data
                                        entStart = neural_frame - entropyPreWindowFrames;
                                        entEnd = neural_frame + entropyPostWindowFrames;
                                        triggeredEntropy = recEntropyData(t, entStart:entEnd);
                                        
                                        % Resample locomotion data to match entropy sampling rate
                                        if length(triggeredLocomotion) ~= totalEntropyFrames
                                            % Check if we have enough points for interpolation
                                            if length(triggeredLocomotion) < 2
                                                fprintf('    Warning: Insufficient locomotion data points (%d) for interpolation, skipping event\n', length(triggeredLocomotion));
                                                continue;
                                            end
                                            
                                            locomotionTimeVector = linspace(-preWindowSec, postWindowSec, length(triggeredLocomotion));
                                            entropyTimeVector = Fig5_Data.TimeVector;
                                            triggeredLocomotion = interp1(locomotionTimeVector, triggeredLocomotion, entropyTimeVector, 'linear', 'extrap');
                                        end
                                        
                                        % Add to collection
                                        Fig5_Data.TriggeredLocomotion = [Fig5_Data.TriggeredLocomotion; triggeredLocomotion];
                                        Fig5_Data.TriggeredEntropy = [Fig5_Data.TriggeredEntropy; triggeredEntropy];
                                        triggered_event_count = triggered_event_count + 1;
                                    end
                                end
                            
                        end
                    end
                end
            end
            
        else
            % Naive/NoSpout: different structure 
            if isfield(condEntropy, 'AllNeurons') && isfield(condEntropy.AllNeurons, 'RecordingEntropy_Raw')
                entropyData = condEntropy.AllNeurons.RecordingEntropy_Raw;
                
                for a = 1:size(condLocomotionData, 2)
                    if ismember(a, recordingsToSkip); continue; end
                    
                    if isfield(condLocomotionData(a), 'Forward_mm') && ...
                       iscell(entropyData) && a <= length(entropyData) && ~isempty(entropyData{a, 1})
                        locomotionTrace = condLocomotionData(a).Forward_mm * LocomotionFrameRate;
                        recEntropyData = entropyData{a, 1}; % [trials x frames]
                        
                        % Trial segmentation
                        TrialLength = floor(length(locomotionTrace) / 150);
                        if TrialLength > 0
                            TrialMovement = 1:TrialLength:length(locomotionTrace);
                            
                            % Get position-specific trials
                            validTrialSegments = [];
                            if exist('Performance', 'var') && ~isempty(Performance) && isfield(Performance, condName) && ...
                               length(Performance.(condName)) >= a && ~strcmp(position, 'All')
                                
                                if strcmp(position, 'P1')
                                    hitField = 'hitP1'; missField = 'missP1';
                                elseif strcmp(position, 'P3')
                                    hitField = 'hitP3'; missField = 'missP3';
                                else
                                    hitField = 'hit'; missField = 'miss';
                                end
                                
                                if isfield(Performance.(condName)(a), hitField) && isfield(Performance.(condName)(a), missField)
                                    hitTrials = Performance.(condName)(a).(hitField);
                                    missTrials = Performance.(condName)(a).(missField);
                                    validTrialSegments = [hitTrials, missTrials];
                                else
                                    validTrialSegments = 1:min([150, length(TrialMovement)-1, size(recEntropyData, 1)]);
                                end
                            else
                                validTrialSegments = 1:min([150, length(TrialMovement)-1, size(recEntropyData, 1)]);
                            end
                            
                            % Process trial segments
                            for t = validTrialSegments
                                
                                trial_start = TrialMovement(t);
                                trial_end = TrialMovement(t+1) - 1;
                                trial_locomotion = locomotionTrace(trial_start:trial_end) ;
                                
                                high_movement = trial_locomotion > threshold;
                                onsets = find(diff([0; high_movement]) == 1);
                                
                                for onset_idx = 1:min(5, length(onsets))
                                    onset_frame = onsets(onset_idx);
                                    neural_frame = round(onset_frame * size(recEntropyData, 2) / length(trial_locomotion));
                                    
                                    if neural_frame > entropyPreWindowFrames && neural_frame <= (size(recEntropyData, 2) - entropyPostWindowFrames) && t <= size(recEntropyData, 1)
                                        % Extract locomotion data with padding if necessary
                                        locStart = onset_frame - locomotionPreWindowFrames;
                                        locEnd = onset_frame + locomotionPostWindowFrames;
                                        
                                        if locStart < 1 || locEnd > length(trial_locomotion)
                                            triggeredLocomotion = nan(1, locomotionPreWindowFrames + locomotionPostWindowFrames + 1);
                                            validStart = max(1, locStart);
                                            validEnd = min(length(trial_locomotion), locEnd);
                                            offsetStart = validStart - locStart + 1;
                                            offsetEnd = offsetStart + (validEnd - validStart);
                                            triggeredLocomotion(offsetStart:offsetEnd) = trial_locomotion(validStart:validEnd);
                                        else
                                            triggeredLocomotion = trial_locomotion(locStart:locEnd);
                                        end
                                        
                                        % Extract entropy data
                                        entStart = neural_frame - entropyPreWindowFrames;
                                        entEnd = neural_frame + entropyPostWindowFrames;
                                        triggeredEntropy = recEntropyData(t, entStart:entEnd);
                                        
                                        % Resample locomotion data to match entropy sampling rate
                                        if length(triggeredLocomotion) ~= totalEntropyFrames
                                            % Check if we have enough points for interpolation
                                            if length(triggeredLocomotion) < 2
                                                fprintf('    Warning: Insufficient locomotion data points (%d) for interpolation, skipping event\n', length(triggeredLocomotion));
                                                continue;
                                            end
                                            
                                            locomotionTimeVector = linspace(-preWindowSec, postWindowSec, length(triggeredLocomotion));
                                            entropyTimeVector = Fig5_Data.TimeVector;
                                            triggeredLocomotion = interp1(locomotionTimeVector, triggeredLocomotion, entropyTimeVector, 'linear', 'extrap');
                                        end
                                        
                                        % Add to collection
                                        Fig5_Data.TriggeredLocomotion = [Fig5_Data.TriggeredLocomotion; triggeredLocomotion];
                                        Fig5_Data.TriggeredEntropy = [Fig5_Data.TriggeredEntropy; triggeredEntropy];
                                        triggered_event_count = triggered_event_count + 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        fprintf('  Added %d triggered events for %s\n', triggered_event_count, condName);
    end
    
    fprintf('  Total triggered events extracted: %d\n', size(Fig5_Data.TriggeredLocomotion, 1));
    fprintf('  Time series length: %d frames (%.1f seconds)\n', totalEntropyFrames, totalEntropyFrames/entropySamplingRate);
end

function [Fig6_Data] = extractFrameByFrameEntropyData(conditions, Entropy, LocomotionCal, Skip, position, Performance, LocomotionFrameRate, threshold)
    % EXTRACTFRAMEBYFRAMEENTROPYDATA - Extract all locomotion and entropy data for frame-by-frame correlation
    
    fprintf('=== Extracting Frame-by-Frame Entropy Data for %s ===\n', position);
    
    % Initialize output structure
    Fig6_Data = struct();
    Fig6_Data.AllLocomotionMagnitudes = [];
    Fig6_Data.AllEntropyValues = [];
    
    total_frames = 0;
    
    % Process each condition
    for c_idx = 1:length(conditions)
        condName = conditions{c_idx};
        fprintf('Processing frame-by-frame entropy data for condition: %s\n', condName);
        
        % Get recordings to skip for this condition
        recordingsToSkip = [];
        if isfield(Skip, condName); recordingsToSkip = Skip.(condName); end
        if ~isempty(recordingsToSkip)
            fprintf('  Skipping recordings %s for %s\n', mat2str(recordingsToSkip), condName);
        end
        
        % Check if condition exists in all required structures
        if ~isfield(Entropy, condName) || ~isfield(LocomotionCal, condName)
            warning('Condition %s not found in one or more required data structures', condName);
            continue;
        end
        
        % Get condition data
        condEntropy = Entropy.(condName);
        condLocomotionData = LocomotionCal.(condName);
        
        frames_this_condition = 0;
        
        % Handle different data structures for different conditions
        if strcmp(condName, 'Expert') || strcmp(condName, 'Beginner')
            % Expert/Beginner: use the same structure as the working functions
            if isfield(condEntropy, 'AllNeurons') && isfield(condEntropy.AllNeurons, 'RecordingEntropy_Raw')
                entropyData = condEntropy.AllNeurons.RecordingEntropy_Raw;
                
                for i = 1:min(length(condLocomotionData), length(entropyData))
                    if ismember(i, recordingsToSkip) || isempty(condLocomotionData(i).RunForward); continue; end
                    
                    if iscell(entropyData) && i <= length(entropyData) && ~isempty(entropyData{i, 1})
                        recEntropyData = entropyData{i, 1}; % [trials x frames]
                        
                        % Get position-specific trials
                        validTrials = [];
                        if exist('Performance', 'var') && ~isempty(Performance) && isfield(Performance, condName) && ...
                           length(Performance.(condName)) >= i && ~strcmp(position, 'All')
                            
                            if strcmp(position, 'P1')
                                hitField = 'hitP1'; missField = 'missP1';
                            elseif strcmp(position, 'P3')
                                hitField = 'hitP3'; missField = 'missP3';
                            else
                                hitField = 'hit'; missField = 'miss';
                            end
                            
                            if isfield(Performance.(condName)(i), hitField) && isfield(Performance.(condName)(i), missField)
                                hitTrials = Performance.(condName)(i).(hitField);
                                missTrials = Performance.(condName)(i).(missField);
                                validTrials = [hitTrials, missTrials];
                            else
                                validTrials = 1:min([size(recEntropyData, 1), 150]);
                            end
                        else
                            validTrials = 1:min([size(recEntropyData, 1), 150]);
                        end
                        
                        % Process position-specific trials
                        for trial = validTrials(1:min(20, length(validTrials)))
                            if trial > size(recEntropyData, 1) || trial > size(condLocomotionData(i).RunForward, 1); continue; end
                            
                            % Get trial locomotion and entropy data
                            locomotionTrace = condLocomotionData(i).RunForward(trial,:)' * LocomotionFrameRate;
                            trialEntropy = recEntropyData(trial, :); % [frames]
                            
                            % Interpolate locomotion to match entropy sampling rate
                            locomotionFrames = length(locomotionTrace);
                            entropyFrames = length(trialEntropy);
                            
                            if locomotionFrames ~= entropyFrames
                                % Check if we have enough points for interpolation
                                if locomotionFrames < 2
                                    fprintf('    Warning: Insufficient locomotion frames (%d) for interpolation, skipping trial\n', locomotionFrames);
                                    continue;
                                end
                                
                                locomotionTimeVector = linspace(0, 1, locomotionFrames);
                                entropyTimeVector = linspace(0, 1, entropyFrames);
                                trialLocomotion = interp1(locomotionTimeVector, locomotionTrace, entropyTimeVector, 'linear', 'extrap');
                            else
                                trialLocomotion = locomotionTrace;
                            end
                            
                            lenFrames = min(length(trialEntropy), length(trialLocomotion));
                            if lenFrames < 10; continue; end % Skip very short trials
                            
                            trialEntropy = trialEntropy(1:lenFrames);
                            trialLocomotion = trialLocomotion(1:lenFrames);
                            
                            % Add to collection
                            Fig6_Data.AllLocomotionMagnitudes = [Fig6_Data.AllLocomotionMagnitudes; trialLocomotion(:)];
                            Fig6_Data.AllEntropyValues = [Fig6_Data.AllEntropyValues; trialEntropy(:)];
                            frames_this_condition = frames_this_condition + lenFrames;
                        end
                    end
                end
            end
            
        else
            % Naive/NoSpout: different structure 
            if isfield(condEntropy, 'AllNeurons') && isfield(condEntropy.AllNeurons, 'RecordingEntropy_Raw')
                entropyData = condEntropy.AllNeurons.RecordingEntropy_Raw;
                
                for a = 1:size(condLocomotionData, 2)
                    if ismember(a, recordingsToSkip); continue; end
                    
                    if isfield(condLocomotionData(a), 'Forward_mm') && ...
                       iscell(entropyData) && a <= length(entropyData) && ~isempty(entropyData{a, 1})
                        locomotionTrace = condLocomotionData(a).Forward_mm * LocomotionFrameRate;
                        recEntropyData = entropyData{a, 1}; % [trials x frames]
                        
                        % Trial segmentation
                        TrialLength = floor(length(locomotionTrace) / 150);
                        if TrialLength > 0
                            TrialMovement = 1:TrialLength:length(locomotionTrace);
                            
                            % Get position-specific trials
                            validTrialSegments = [];
                            if exist('Performance', 'var') && ~isempty(Performance) && isfield(Performance, condName) && ...
                               length(Performance.(condName)) >= a && ~strcmp(position, 'All')
                                
                                if strcmp(position, 'P1')
                                    hitField = 'hitP1'; missField = 'missP1';
                                elseif strcmp(position, 'P3')
                                    hitField = 'hitP3'; missField = 'missP3';
                                else
                                    hitField = 'hit'; missField = 'miss';
                                end
                                
                                if isfield(Performance.(condName)(a), hitField) && isfield(Performance.(condName)(a), missField)
                                    hitTrials = Performance.(condName)(a).(hitField);
                                    missTrials = Performance.(condName)(a).(missField);
                                    validTrialSegments = [hitTrials, missTrials];
                                else
                                    validTrialSegments = 1:min([150, length(TrialMovement)-1, size(recEntropyData, 1)]);
                                end
                            else
                                validTrialSegments = 1:min([150, length(TrialMovement)-1, size(recEntropyData, 1)]);
                            end
                            
                            % Process trial segments
                            for t = validTrialSegments
                                
                                trial_start = TrialMovement(t);
                                trial_end = TrialMovement(t+1) - 1;
                                trial_locomotion = locomotionTrace(trial_start:trial_end) * TrialLength / 12;
                                
                                % Get trial entropy data
                                trialEntropy = recEntropyData(t, :); % [frames]
                                
                                % Interpolate to match sampling rates
                                locomotionFrames = length(trial_locomotion);
                                entropyFrames = length(trialEntropy);
                                
                                if locomotionFrames ~= entropyFrames
                                    % Check if we have enough points for interpolation
                                    if locomotionFrames < 2
                                        fprintf('    Warning: Insufficient locomotion frames (%d) for interpolation, skipping trial\n', locomotionFrames);
                                        continue;
                                    end
                                    
                                    locomotionTimeVector = linspace(0, 1, locomotionFrames);
                                    entropyTimeVector = linspace(0, 1, entropyFrames);
                                    trial_locomotion = interp1(locomotionTimeVector, trial_locomotion, entropyTimeVector, 'linear', 'extrap');
                                end
                                
                                lenFrames = min(length(trialEntropy), length(trial_locomotion));
                                if lenFrames < 10; continue; end % Skip very short trials
                                
                                trialEntropy = trialEntropy(1:lenFrames);
                                trial_locomotion = trial_locomotion(1:lenFrames);
                                
                                % Add to collection
                                Fig6_Data.AllLocomotionMagnitudes = [Fig6_Data.AllLocomotionMagnitudes; trial_locomotion(:)];
                                Fig6_Data.AllEntropyValues = [Fig6_Data.AllEntropyValues; trialEntropy(:)];
                                frames_this_condition = frames_this_condition + lenFrames;
                            end
                        end
                    end
                end
            end
        end
        
        fprintf('  Added %d frames for %s\n', frames_this_condition, condName);
        total_frames = total_frames + frames_this_condition;
    end
    
    fprintf('  Total frames extracted: %d\n', total_frames);
    fprintf('  Valid data points: %d\n', sum(~isnan(Fig6_Data.AllLocomotionMagnitudes) & ~isnan(Fig6_Data.AllEntropyValues)));
end

function validateLocomotionEntropyDataConsistency(Fig1_Data, Fig3_Data, Fig4_Data)
    % VALIDATELOCOMOTIONENTROPYDATACONSISTENCY - Validate consistency across extractions
    
    fprintf('\n=== LOCOMOTION ENTROPY DATA CONSISTENCY VALIDATION ===\n');
    
    % Check consistency of array lengths
    nPrePost = length(Fig1_Data.PreLocomotionEntropy);
    nLabels = length(Fig1_Data.AllConditionLabels);
    nRatios = length(Fig3_Data.EntropyRatio);
    nChanges = length(Fig4_Data.EntropyChange);
    
    fprintf('Data array lengths:\n');
    fprintf('  Pre/Post entropy: %d\n', nPrePost);
    fprintf('  Condition labels: %d\n', nLabels);
    fprintf('  Entropy ratios: %d\n', nRatios);
    fprintf('  Entropy changes: %d\n', nChanges);
    
    % Check consistency
    consistency_check = true;
    
    if nPrePost ~= nLabels
        fprintf('ERROR: Pre/Post entropy and condition labels length mismatch!\n');
        consistency_check = false;
    end
    
    if nPrePost ~= nRatios
        fprintf('ERROR: Pre/Post entropy and ratios length mismatch!\n');
        consistency_check = false;
    end
    
    if nPrePost ~= nChanges
        fprintf('ERROR: Pre/Post entropy and changes length mismatch!\n');
        consistency_check = false;
    end
    
    if consistency_check
        fprintf('✓ All locomotion entropy data extractions are CONSISTENT!\n');
    else
        fprintf('✗ Locomotion entropy data INCONSISTENCIES detected - please review\n');
    end
    
    fprintf('=== END VALIDATION ===\n\n');
end

% Helper function to set nice axis ticks at half and full numbers
function setNiceAxisTicks(ax, dataLim)
    % Generate nice tick values at half and full numbers (e.g., 1.5, 2, 2.5, etc.)
    range = dataLim(2) - dataLim(1);
    
    if range <= 0
        return; % Invalid range
    end
    
    % Determine appropriate tick spacing
    if range <= 2
        tickSpacing = 0.5;
    elseif range <= 5
        tickSpacing = 0.5;
    elseif range <= 10
        tickSpacing = 1;
    else
        tickSpacing = 2;
    end
    
    % Find start and end ticks
    startTick = ceil(dataLim(1) / tickSpacing) * tickSpacing;
    endTick = floor(dataLim(2) / tickSpacing) * tickSpacing;
    
    % Generate tick values
    ticks = startTick:tickSpacing:endTick;
    
    % Ensure we have at least a few ticks
    if length(ticks) < 3
        % Fall back to automatic ticking if our algorithm produces too few ticks
        return;
    end
    
    % Apply ticks to both X and Y axes
    set(ax, 'XTick', ticks, 'YTick', ticks);
end 