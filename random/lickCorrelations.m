%% Lick correlations
function lickCorrelations(ActivityData, Entropy, WindowLength, params, Rec, Performance, conditions, Skip)
% LICKCORRELATIONS - Analyse pre- and post-lick activity/entropy and produce
% correlation and distribution plots across different behavioral conditions.
%
% INPUTS:
%   ActivityData - Structure containing neural activity data for different conditions
%   Entropy      - Structure containing entropy data for different conditions
%   WindowLength - Structure containing window length information
%   params       - Structure containing plotting parameters including colors
%   Rec          - Recording metadata structure
%   Performance  - Performance metrics structure containing lick data
%   conditions   - Cell array of condition names to analyze (e.g., {'Expert', 'Beginner'})
%   Skip         - (optional) Structure containing recording indices to skip for each condition
%
% OUTPUTS:
%   • Figure 1 : Pre vs Post dF/F scatter   – Wilcoxon signed‑rank (paired)
%   • Figure 2 : Pre vs Post Entropy scatter – Wilcoxon signed‑rank (paired)
%   • Figure 3 : Histograms of Post/Pre ratios (separate tiles) – sign‑rank
%                against a median of 1
%   • Figure 4 : Pre value vs Change scatter (dF/F & Entropy) – Spearman ρ
%                with p‑value; also signed‑rank of change ≠ 0
%
% EXAMPLE:
%   lickCorrelations(ActivityData, Entropy, WindowLength, params, Rec, Performance);
%   lickCorrelations(ActivityData, Entropy, WindowLength, params, Rec, Performance, {'Expert'});
%   lickCorrelations(ActivityData, Entropy, WindowLength, params, Rec, Performance, {'Expert', 'Beginner'}, Skip);

%% Check inputs and set defaults
if nargin < 8 || isempty(Skip)
    Skip = struct(); 
end

if nargin < 7 || isempty(conditions)
    % Default to all available conditions if not specified
    conditions = {};
    if isfield(ActivityData, 'Expert'); conditions{end+1} = 'Expert'; end
    if isfield(ActivityData, 'Beginner'); conditions{end+1} = 'Beginner'; end
    if isfield(ActivityData, 'Naive'); conditions{end+1} = 'Naive'; end
    if isfield(ActivityData, 'NoSpout'); conditions{end+1} = 'NoSpout'; end
end

fprintf('Processing lick correlations for conditions: %s\n', strjoin(conditions, ', '));

% Print Skip information if provided
if ~isempty(fieldnames(Skip))
    fprintf('Skip information provided for conditions: %s\n', strjoin(fieldnames(Skip), ', '));
    for skip_field = fieldnames(Skip)'
        fprintf('  %s: skipping recordings %s\n', skip_field{1}, mat2str(Skip.(skip_field{1})));
    end
else
    fprintf('No recordings will be skipped (Skip structure is empty)\n');
end

%% ------------------------- Data Accumulation ------------------------ %%
% Time axis (–8 s to +10.5 s around lick)
timeReference = linspace(-8, 10.5, WindowLength.All);
ifi = 18.5/WindowLength.All; %#ok<NASGU>

% Initialize data arrays
PreExtraneousLickActivity = [];
PostExtraneousLickActivity = [];
PreExtraneousLickEntropy = [];
PostExtraneousLickEntropy = [];
AllConditionLabels = {};

fprintf('\n=== Processing Lick Data Extraction ===\n');

% Process each condition
for c_idx = 1:length(conditions)
    condName = conditions{c_idx};
    fprintf('Processing lick data for condition: %s\n', condName);
    
    % Get recordings to skip for this condition
    recordingsToSkip = [];
    if isfield(Skip, condName); recordingsToSkip = Skip.(condName); end
    if ~isempty(recordingsToSkip)
        fprintf('  Skipping recordings %s for %s\n', mat2str(recordingsToSkip), condName);
    end
    
    % Check if condition exists in all required structures
    if ~isfield(ActivityData, condName) || ~isfield(Entropy, condName) || ~isfield(Performance, condName) || ~isfield(Rec, condName)
        warning('Condition %s not found in one or more required data structures', condName);
        continue;
    end
    
    % Get condition data
    condActivityData = ActivityData.(condName);
    condEntropy = Entropy.(condName);
    condPerformance = Performance.(condName);
    condRec = Rec.(condName);
    
    % Determine number of recordings for this condition
    if isfield(condRec, 'AnimalID')
        numRecs = length(condRec.AnimalID);
    else
        numRecs = size(condRec, 1);
    end
    
    lick_count_this_condition = 0;
    
    % Process each recording
    for r = 1:numRecs
        if ismember(r, recordingsToSkip); continue; end
        
        % Check if this recording has the required data
        if r > length(condActivityData) || r > length(condPerformance)
            continue;
        end
        
        if ~isfield(condActivityData(r), 'all') || ~isfield(condPerformance(r), 'LicksInFrame')
            continue;
        end
        
        % Check if entropy data exists for this recording
        if ~isfield(condEntropy, 'AllNeurons') || ~isfield(condEntropy.AllNeurons, 'RecordingEntropy_Raw')
            warning('Entropy data not found for condition %s', condName);
            continue;
        end
        
        entropyData = condEntropy.AllNeurons.RecordingEntropy_Raw;
        if iscell(entropyData)
            if r > length(entropyData) || isempty(entropyData{r, 1})
                continue;
            end
            recEntropyData = entropyData{r, 1};
        else
            continue;
        end
        
        % Process each trial for this recording
        for t = 1:size(condPerformance(r).LicksInFrame, 2)
            % Find extraneous licks (licks between -7.9 and -0.2 seconds)
            extraneousLicks = condPerformance(r).LicksInFrame( ...
                condPerformance(r).LicksInFrame(:, t) < -0.2 & ...
                condPerformance(r).LicksInFrame(:, t) > -7.9, t);
            
            % Process each extraneous lick
            for e = 1:numel(extraneousLicks)
                StartFrame = find(extraneousLicks(e) < timeReference, 1);
                
                % Check if we have enough frames before and after
                if StartFrame > 2 && StartFrame < (size(condActivityData(r).all, 2) - 2)
                    % Extract pre and post activity
                    preActivity = mean(condActivityData(r).all(:, StartFrame-2, t), 1);
                    postActivity = mean(condActivityData(r).all(:, StartFrame+2, t), 1);
                    
                    % Extract pre and post entropy
                    if t <= size(recEntropyData, 1) && StartFrame-2 > 0 && StartFrame+2 <= size(recEntropyData, 2)
                        preEntropy = recEntropyData(t, StartFrame-2);
                        postEntropy = recEntropyData(t, StartFrame+2);
                        
                        % Add to arrays
                        PreExtraneousLickActivity = [PreExtraneousLickActivity; preActivity]; %#ok<AGROW>
                        PostExtraneousLickActivity = [PostExtraneousLickActivity; postActivity]; %#ok<AGROW>
                        PreExtraneousLickEntropy = [PreExtraneousLickEntropy; preEntropy]; %#ok<AGROW>
                        PostExtraneousLickEntropy = [PostExtraneousLickEntropy; postEntropy]; %#ok<AGROW>
                        AllConditionLabels{end+1} = condName; %#ok<AGROW>
                        lick_count_this_condition = lick_count_this_condition + 1;
                    end
                end
            end
        end
    end
    
    fprintf('  Added %d extraneous licks for %s\n', lick_count_this_condition, condName);
end

%% Check if we have any data
if isempty(PreExtraneousLickActivity)
    warning('No extraneous lick data found for any conditions');
    return;
end

fprintf('Total extraneous licks found: %d across %d conditions\n', length(PreExtraneousLickActivity), length(conditions));

%% Convenience variables
ActivityRatio = PostExtraneousLickActivity ./ PreExtraneousLickActivity;
EntropyRatio = PostExtraneousLickEntropy ./ PreExtraneousLickEntropy;
ActivityChange = PostExtraneousLickActivity - PreExtraneousLickActivity;
EntropyChange = PostExtraneousLickEntropy - PreExtraneousLickEntropy;

%% ------------------------- Figure 1  dF/F scatter -------------------- %%
figure('Color', 'w', 'Name', 'Pre vs Post dF/F Scatter'); 
hold on;
scatter(PreExtraneousLickActivity, PostExtraneousLickActivity, 30, [0 0 0], 'filled');
axis square;
ref = refline(1); 
ref.LineWidth = 5; 
ref.Color = 'r';

[pAct, ~, statsAct] = signrank(PreExtraneousLickActivity, PostExtraneousLickActivity);
ratioAct = mean(PreExtraneousLickActivity ./ PostExtraneousLickActivity, 'omitnan');

ax = gca;
ax.FontSize = 20;
xlabel('Pre-Lick dF/F');
ylabel('Post-Lick dF/F');

title({'Pre vs Post dF/F', sprintf('ratio = %.3f,  p = %.3g', ratioAct, pAct), ...
       sprintf('n = %d licks from %d conditions', length(PreExtraneousLickActivity), length(conditions))});

%% Half‑violin (dF/F)
figure('Color', 'w', 'Name', 'Pre vs Post dF/F Violin Plot');
axesFont = 24;
daviolinplot({[PreExtraneousLickActivity, PostExtraneousLickActivity]}, ...
             'groups', ones(numel(PreExtraneousLickActivity), 1), ...
             'colors', params.colour(13, :), 'box', 3, 'boxcolor', 'w', ...
             'scatter', 2, 'jitter', 1, 'scattercolor', 'same', 'scattersize', 20, ...
             'scatteralpha', 0.7, 'linkline', 0, 'withinlines', 0, ...
             'xtlabels', {"Before Lick", "After Lick"});
set(gca, 'FontSize', axesFont);
ylabel('dF/F');
title(gca, sprintf('dF/F: p = %.3g (n = %d)', pAct, length(PreExtraneousLickActivity)));

%% ------------------------- Figure 2  Entropy scatter ----------------- %%
figure('Color', 'w', 'Name', 'Pre vs Post Entropy Scatter'); 
hold on;
scatter(PreExtraneousLickEntropy, PostExtraneousLickEntropy, 30, [0 0 0], 'filled');
axis square;
ref = refline(1); 
ref.LineWidth = 5; 
ref.Color = 'r';

[pEnt, ~, statsEnt] = signrank(PreExtraneousLickEntropy, PostExtraneousLickEntropy);
ratioEnt = mean(PreExtraneousLickEntropy ./ PostExtraneousLickEntropy, 'omitnan');

ax = gca;
ax.FontSize = 24;
xlabel('Pre-Lick Entropy');
ylabel('Post-Lick Entropy');

title(ax, {'Pre vs Post Entropy', sprintf('ratio = %.3f,  p = %.3g', ratioEnt, pEnt), ...
          sprintf('n = %d licks from %d conditions', length(PreExtraneousLickEntropy), length(conditions))});

%% Half‑violin (Entropy)
figure('Color', 'w', 'Name', 'Pre vs Post Entropy Violin Plot');
daviolinplot({[PreExtraneousLickEntropy, PostExtraneousLickEntropy]}, ...
             'groups', ones(numel(PreExtraneousLickEntropy), 1), ...
             'colors', params.colour(13, :), 'box', 3, 'boxcolor', 'w', ...
             'scatter', 2, 'jitter', 1, 'scattercolor', 'same', 'scattersize', 20, ...
             'scatteralpha', 0.7, 'linkline', 1, 'withinlines', 0, ...
             'xtlabels', {"Before Lick", "After Lick"});
set(gca, 'FontSize', axesFont);
ylabel('Entropy');
title(gca, sprintf('Entropy: p = %.3g (n = %d)', pEnt, length(PreExtraneousLickEntropy)));

%% ------------------------- Figure 3  Ratio Histograms --------------- %%
figure('Color', 'w', 'Name', 'Post/Pre Ratio Histograms');
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% dF/F ratio tile
ax1 = nexttile;
histogram(ax1, ActivityRatio, 'Normalization', 'probability', ...
          'BinMethod', 'sturges', 'FaceColor', params.colour(7, :), ...
          'EdgeColor', 'none', 'FaceAlpha', 0.7);
[pRatioAct] = signrank(ActivityRatio - 1);
ax1.FontSize = axesFont;
ax1.Box = 'on';
ax1.XLabel.String = 'Post / Pre dF/F';
ax1.YLabel.String = 'Probability';

title(ax1, {'dF/F Ratio', sprintf('p = %.3g', pRatioAct), sprintf('n = %d', length(ActivityRatio))});

% Entropy ratio tile
ax2 = nexttile;
histogram(ax2, EntropyRatio, 'Normalization', 'probability', ...
          'BinMethod', 'sturges', 'FaceColor', params.colour(9, :), ...
          'EdgeColor', 'none', 'FaceAlpha', 0.7);
[pRatioEnt] = signrank(EntropyRatio - 1);
ax2.FontSize = axesFont;
ax2.Box = 'on';
ax2.XLabel.String = 'Post / Pre Entropy';
ax2.YLabel.String = 'Probability';

title(ax2, {'Entropy Ratio', sprintf('p = %.3g', pRatioEnt), sprintf('n = %d', length(EntropyRatio))});

%% ------------------------- Figure 4  Change vs Pre ------------------ %%
figure('Color', 'w', 'Name', 'Change vs Pre-Lick Values');

% dF/F change
ax3 = subplot(1, 2, 1);
scatter(ax3, PreExtraneousLickActivity, ActivityChange, 30, params.colour(7, :), 'filled');
lsline(ax3);
axis(ax3, 'square'); 
grid(ax3, 'on');
[rhoAct, pRhoAct] = corr(PreExtraneousLickActivity, ActivityChange, ...
                         'Type', 'Spearman', 'Rows', 'complete');
ax3.FontSize = axesFont;
ax3.XLabel.String = 'Pre dF/F';
ax3.YLabel.String = '\Delta dF/F (Post - Pre)';

title(ax3, {'Activity Change vs Pre', sprintf('rho = %.2f,  p = %.3g', rhoAct, pRhoAct), ...
           sprintf('n = %d', length(PreExtraneousLickActivity))});

% Entropy change
ax4 = subplot(1, 2, 2);
scatter(ax4, PreExtraneousLickEntropy, EntropyChange, 30, params.colour(9, :), 'filled');
lsline(ax4);
axis(ax4, 'square'); 
grid(ax4, 'on');
[rhoEnt, pRhoEnt] = corr(PreExtraneousLickEntropy, EntropyChange, ...
                         'Type', 'Spearman', 'Rows', 'complete');
ax4.FontSize = axesFont;
ax4.XLabel.String = 'Pre Entropy';
ax4.YLabel.String = '\Delta Entropy (Post - Pre)';

title(ax4, {'Entropy Change vs Pre', sprintf('rho = %.2f,  p = %.3g', rhoEnt, pRhoEnt), ...
           sprintf('n = %d', length(PreExtraneousLickEntropy))});

%% Summary output
fprintf('\nlickCorrelations completed successfully!\n');
fprintf('Processed %d conditions: %s\n', length(conditions), strjoin(conditions, ', '));
fprintf('Total extraneous licks analyzed: %d\n', length(PreExtraneousLickActivity));

% Print condition breakdown
unique_conditions = unique(AllConditionLabels);
for i = 1:length(unique_conditions)
    count = sum(strcmp(AllConditionLabels, unique_conditions{i}));
    fprintf('  %s: %d licks (%.1f%%)\n', unique_conditions{i}, count, 100*count/length(AllConditionLabels));
end

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
