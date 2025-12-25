%% Plot Number of Timeouts Across Trials

TimeOuts = data.options.timeoutHistory  ;

% Separate by opto condition
TimeOuts_control = TimeOuts(optoTrials == 0);
trials_control     = find(optoTrials == 0);

TimeOuts_opto    = TimeOuts(optoTrials == 1);
trials_opto        = find(optoTrials == 1);

% Keep only trials with >0 errors
keep_control = TimeOuts_control > 0;
TimeOuts_control = TimeOuts_control(keep_control);
trials_control = trials_control(keep_control);

keep_opto = TimeOuts_opto > 0;
TimeOuts_opto = TimeOuts_opto(keep_opto);
trials_opto = trials_opto(keep_opto);


%Plot 1: Scatter plot showing progression across all trials
figure('Color', 'w', 'Position', [100 100 800 500]);
hold on;

% Shade opto blocks in background
if any(optoTrials == 1)
    yLims = [0 max(nErrorLicks) + 1];
    for t = 1:nTrials
        if optoTrials(t) == 1
            fill([t-0.5 t+0.5 t+0.5 t-0.5], [yLims(1) yLims(1) yLims(2) yLims(2)], ...
                [0.4 0.6 1], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        end
    end
end

% Plot data
scatter(trials_control, TimeOuts_control, 50, [0.5 0.5 0.5], 'filled', ...
    'MarkerFaceAlpha', 0.6, 'DisplayName', 'Control');
scatter(trials_opto, TimeOuts_opto, 50, [0.4 0.6 1], 'filled', ...
    'MarkerFaceAlpha', 0.6, 'DisplayName', 'Opto');

% Add smooth trend line
if length(trials_control) > 3
    p_control = polyfit(trials_control, TimeOuts_control, 1);
    plot(trials_control, polyval(p_control, trials_control), '-', ...
        'Color', [0.3 0.3 0.3], 'LineWidth', 2.5);
end
if length(trials_opto) > 3
    p_opto = polyfit(trials_opto, TimeOuts_opto, 1);
    plot(trials_opto, polyval(p_opto, trials_opto), '-', ...
        'Color', [0.2 0.4 0.8], 'LineWidth', 2.5);
end

xlabel('Trial Number', 'FontSize', 14);
ylabel('Number of Timeouts', 'FontSize', 14);
title('Timeouts', 'FontSize', 15, 'FontWeight', 'bold');
set(gca, 'FontSize', 13, 'LineWidth', 1.5, 'Box', 'off');
grid on;
grid minor;

%% NORMALISED number of trials with timeouts

% Count trials with timeouts (>0)
nTimeouts_control = sum(TimeOuts_control > 0);
nTimeouts_opto = sum(TimeOuts_opto > 0);

% Total trials per condition
nTrials_control = length(TimeOuts_control);
nTrials_opto = length(TimeOuts_opto);

% Calculate normalized rates (%)
timeoutRate_control = (nTimeouts_control / nTrials_control) * 100;
timeoutRate_opto = (nTimeouts_opto / nTrials_opto) * 100;

% Plot
figure('Color', 'w', 'Position', [100 100 450 500]);

% Data
barData = [timeoutRate_control, timeoutRate_opto];

% Create bar plot
b = bar(barData, 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 1.5);
b.CData(1,:) = [0.3 0.3 0.3];  % Control = dark gray
b.CData(2,: ) = [0.4 0.6 1];     % Opto = blue

hold on;

% Add 95% confidence intervals
% Wilson score interval for proportions
z = 1.96; % 95% CI
p_control = timeoutRate_control / 100;
p_opto = timeoutRate_opto / 100;

ci_control = z * sqrt(p_control * (1-p_control) / nTrials_control) * 100;
ci_opto = z * sqrt(p_opto * (1-p_opto) / nTrials_opto) * 100;

errorbar(1:2, barData, [ci_control, ci_opto], 'k.', ... 
    'LineWidth', 2, 'CapSize', 12);

% Add actual numbers on bars
text(1, timeoutRate_control + ci_control + 3, ... 
    sprintf('%d/%d', nTimeouts_control, nTrials_control), ...
    'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
text(2, timeoutRate_opto + ci_opto + 3, ... 
    sprintf('%d/%d', nTimeouts_opto, nTrials_opto), ...
    'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');

% Statistics (Fisher's exact test)
contingencyTable = [nTimeouts_control, nTrials_control - nTimeouts_control; ... 
                    nTimeouts_opto, nTrials_opto - nTimeouts_opto];
[~, p_fisher] = fishertest(contingencyTable);

% Add significance
if p_fisher < 0.001
    sigText = '***';
elseif p_fisher < 0.01
    sigText = '**';
elseif p_fisher < 0.05
    sigText = '*';
else
    sigText = 'n.s.';
end

% Significance line
yMax = max(barData + [ci_control, ci_opto]);
yPos = yMax * 1.2;
plot([1 2], [yPos yPos], 'k-', 'LineWidth', 2);
text(1.5, yPos*1.05, sigText, 'FontSize', 16, ... 
    'HorizontalAlignment', 'center', 'FontWeight', 'bold');

% Labels
set(gca, 'XTickLabel', {'Control', 'Opto'}, 'FontSize', 14);
ylabel('Timeout Rate (%)', 'FontSize', 15);
title('Normalized Timeout Rate', 'FontSize', 16, 'FontWeight', 'bold');
set(gca, 'FontSize', 13, 'LineWidth', 1.5, 'Box', 'off');
ylim([0 yMax*1.25]);
grid on; grid minor;
%% Timeouts across opto block positions

% Parameters
OPTO_BLOCK_OFF = 10;
OPTO_BLOCK_ON = 5;
BLOCK_SIZE = OPTO_BLOCK_OFF + OPTO_BLOCK_ON;

% Get timeout data
TimeOuts = data.options.timeoutHistory;
nTrials = length(TimeOuts);

% Initialize arrays for OPTO ON blocks
opto_position = [];
timeout_count = [];
block_number = [];

% Initialize arrays for OPTO OFF blocks (control)
control_position = [];
control_timeout_count = [];
control_block_number = [];

% Find all blocks and extract timeout data
blockNum = 1;
for startBlock = 1 : BLOCK_SIZE : nTrials
    % OPTO OFF trials (first 10 trials of block)
    control_start = startBlock;
    control_end = min(startBlock + OPTO_BLOCK_OFF - 1, nTrials);
    control_trials = control_start: control_end;
    control_timeouts = TimeOuts(control_trials);
    
    positions = 1:length(control_timeouts);
    control_position = [control_position; positions(: )];
    control_timeout_count = [control_timeout_count; control_timeouts(:)];
    control_block_number = [control_block_number; repmat(blockNum, length(control_timeouts), 1)];
    
    % OPTO ON trials (next 5 trials of block)
    opto_start = startBlock + OPTO_BLOCK_OFF;
    if opto_start <= nTrials
        opto_end = min(opto_start + OPTO_BLOCK_ON - 1, nTrials);
        opto_trials = opto_start:opto_end;
        opto_timeouts = TimeOuts(opto_trials);
        
        positions = 1:length(opto_timeouts);
        opto_position = [opto_position; positions(:)];
        timeout_count = [timeout_count; opto_timeouts(:)];
        block_number = [block_number; repmat(blockNum, length(opto_timeouts), 1)];
    end
    
    blockNum = blockNum + 1;
end

%% Plot 1: Individual blocks overlaid (Control vs Opto)
figure('Color', 'w', 'Position', [100 100 900 550]);
hold on;

nBlocks = max(block_number);
cmap_control = winter(nBlocks);
cmap_opto = autumn(nBlocks);

% Plot control blocks (blue tones)
for b = 1:nBlocks
    idx = control_block_number == b;
    pos = control_position(idx);
    counts = control_timeout_count(idx);
    plot(pos, counts, '--o', 'Color', cmap_control(b,:), 'LineWidth', 1, ...
        'MarkerSize', 6, 'MarkerFaceColor', cmap_control(b,:));
end

% Plot opto blocks (red tones)
for b = 1:nBlocks
    idx = block_number == b;
    pos = opto_position(idx);
    counts = timeout_count(idx);
    plot(pos + OPTO_BLOCK_OFF, counts, '-o', 'Color', cmap_opto(b,:), 'LineWidth', 1.5, ...
        'MarkerSize', 8, 'MarkerFaceColor', cmap_opto(b,:));
end

% Add mean trend lines
mean_control = zeros(OPTO_BLOCK_OFF, 1);
for p = 1:OPTO_BLOCK_OFF
    mean_control(p) = mean(control_timeout_count(control_position == p));
end
plot(1:OPTO_BLOCK_OFF, mean_control, 'b-', 'LineWidth', 4, 'DisplayName', 'Control Mean');

mean_opto = zeros(OPTO_BLOCK_ON, 1);
for p = 1:OPTO_BLOCK_ON
    mean_opto(p) = mean(timeout_count(opto_position == p));
end
plot((1:OPTO_BLOCK_ON) + OPTO_BLOCK_OFF, mean_opto, 'r-', 'LineWidth', 4, 'DisplayName', 'Opto Mean');

% Add vertical line to separate control and opto
yl = ylim;
plot([OPTO_BLOCK_OFF + 0.5, OPTO_BLOCK_OFF + 0.5], yl, 'k--', 'LineWidth', 2, 'DisplayName', 'Opto Onset');

xlabel('Trial Position in Block', 'FontSize', 14);
ylabel('Number of Timeouts', 'FontSize', 14);
title('Timeouts Across Block (Control vs Opto)', 'FontSize', 15, 'FontWeight', 'bold');
set(gca, 'FontSize', 13, 'LineWidth', 1.5, 'Box', 'off');
xlim([0.5 BLOCK_SIZE + 0.5]);
xticks([1 5 10 11 12 13 14 15]);
xticklabels({'1', '5', '10', 'O1', 'O2', 'O3', 'O4', 'O5'});
grid on;

% Plot 2: Summary comparison with error bars
figure('Color', 'w', 'Position', [100 100 900 550]);
hold on;

% Calculate means and SEMs
mean_control = zeros(OPTO_BLOCK_OFF, 1);
sem_control = zeros(OPTO_BLOCK_OFF, 1);
for p = 1:OPTO_BLOCK_OFF
    data_points = control_timeout_count(control_position == p);
    mean_control(p) = mean(data_points);
    sem_control(p) = std(data_points) / sqrt(length(data_points));
end

mean_opto = zeros(OPTO_BLOCK_ON, 1);
sem_opto = zeros(OPTO_BLOCK_ON, 1);
for p = 1:OPTO_BLOCK_ON
    data_points = timeout_count(opto_position == p);
    mean_opto(p) = mean(data_points);
    sem_opto(p) = std(data_points) / sqrt(length(data_points));
end

% Plot control with error bars
errorbar(1:OPTO_BLOCK_OFF, mean_control, sem_control, 'o-', ...
    'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', [0.5 0.7 0.9], ...
    'Color', [0.2 0.4 0.7], 'CapSize', 8, 'DisplayName', 'Control');

% Plot opto with error bars
errorbar((1:OPTO_BLOCK_ON) + OPTO_BLOCK_OFF, mean_opto, sem_opto, 'o-', ...
    'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', [1 0.6 0.6], ...
    'Color', [0.8 0.2 0.2], 'CapSize', 8, 'DisplayName', 'Opto');

% Add vertical line
yl = ylim;
plot([OPTO_BLOCK_OFF + 0.5, OPTO_BLOCK_OFF + 0.5], yl, 'k--', 'LineWidth', 2);

xlabel('Trial Position in Block', 'FontSize', 14);
ylabel('Number of Timeouts', 'FontSize', 14);
title('Mean Timeouts: Control vs Opto Blocks', 'FontSize', 15, 'FontWeight', 'bold');
set(gca, 'FontSize', 13, 'LineWidth', 1.5, 'Box', 'off');
xlim([0.5 BLOCK_SIZE + 0.5]);
xticks([1 5 10 11 12 13 14 15]);
xticklabels({'1', '5', '10', 'O1', 'O2', 'O3', 'O4', 'O5'});
grid on;

%% False Alarm Rate *premature licks* calculation (as defined by Arnau's -2s lick period relative to stim)

% This script analyzes false alarm rates from behavioral data
% Run RB_OptoAnalysis_hitmiss.m first to load and process data

FA_WINDOW_START = -2;  % Start of false alarm window (seconds before stimulus)
FA_WINDOW_END = 0;      % End of false alarm window (stimulus onset)

%CALCULATE FALSE ALARMS
[FalseAlarms, FalseAlarmRate] = calculateFalseAlarms(LicksInFrame, discardedTrials, ... 
    FA_WINDOW_START, FA_WINDOW_END);

% DISPLAY RESULTS
fprintf('\n===== FALSE ALARM ANALYSIS =====\n');
fprintf('False Alarm Window: %.1f to %.1f seconds\n', FA_WINDOW_START, FA_WINDOW_END);
fprintf('Overall False Alarm Rate: %.2f%% (%d/%d trials)\n\n', ...
    FalseAlarmRate * 100, sum(FalseAlarms), sum(~isDiscarded));

% False alarms by condition (Control vs Opto)
if any(optoTrials == 1)
    FA_Control = sum(FalseAlarms & optoTrials == 0 & ~isDiscarded);
    FA_Opto = sum(FalseAlarms & optoTrials == 1 & ~isDiscarded);
    nControl = sum(optoTrials == 0 & ~isDiscarded);
    nOpto = sum(optoTrials == 1 & ~isDiscarded);
    
    if nControl > 0
        FA_Rate_Control = FA_Control / nControl;
        fprintf('Control False Alarm Rate:   %.2f%% (%d/%d trials)\n', ...
            FA_Rate_Control*100, FA_Control, nControl);
    end
    
    if nOpto > 0
        FA_Rate_Opto = FA_Opto / nOpto;
        fprintf('Opto False Alarm Rate:  %.2f%% (%d/%d trials)\n', ...
            FA_Rate_Opto*100, FA_Opto, nOpto);
    end
    
    % Statistical comparison
    if nControl > 0 && nOpto > 0
        fprintf('\nDifference:   %.2f%%\n', (FA_Rate_Opto - FA_Rate_Control)*100);
    end
end

% False alarms by location (if applicable)
if exist('stimLocations', 'var') && ~isempty(stimLocations)
    fprintf('\n--- By Location ---\n');
    
    FA_Right = sum(FalseAlarms & stimLocations == 6 & ~isDiscarded);
    FA_Left = sum(FalseAlarms & stimLocations == 7 & ~isDiscarded);
    nRight = sum(stimLocations == 6 & ~isDiscarded);
    nLeft = sum(stimLocations == 7 & ~isDiscarded);
    
    if nRight > 0
        fprintf('Right Stimulus FA Rate: %.2f%% (%d/%d trials)\n', ...
            (FA_Right/nRight)*100, FA_Right, nRight);
    end
    
    if nLeft > 0
        fprintf('Left Stimulus FA Rate:   %.2f%% (%d/%d trials)\n', ...
            (FA_Left/nLeft)*100, FA_Left, nLeft);
    end
end

fprintf('================================\n\n');

% PLOT FALSE ALARM RASTER
plotFalseAlarmRaster(LicksInFrame, FalseAlarms, nTrials, optoTrials, ... 
    FA_WINDOW_START, FA_WINDOW_END);

% PLOT FALSE ALARM RATES BY CONDITION
if any(optoTrials == 1)
    plotFalseAlarmsByCondition(FalseAlarms, optoTrials, isDiscarded);
end

%PLOT FALSE ALARM RATES BY CONDITION AND LOCATION
if exist('stimLocations', 'var') && ~isempty(stimLocations) && any(optoTrials == 1)
    plotFalseAlarmsByConditionAndLocation(FalseAlarms, optoTrials, stimLocations, isDiscarded);
end

% FUNCTIONS

function [FalseAlarms, FalseAlarmRate] = calculateFalseAlarms(LicksInFrame, ... 
    discardedTrials, faWindowStart, faWindowEnd)
    % Calculate false alarm trials (premature licks between faWindowStart and faWindowEnd)
    % Returns a logical array indicating FA trials and the overall FA rate
    
    nTrials = size(LicksInFrame, 2);
    FalseAlarms = false(1, nTrials);
    
    for t = 1:nTrials
        % Skip discarded trials
        if ~isempty(discardedTrials) && ismember(t, discardedTrials)
            continue
        end
        
        % Get all licks for this trial
        allLicks = LicksInFrame(LicksInFrame(: ,t) ~= 0, t);
        
        % Check if any licks fall in the false alarm window
        prematureLicks = allLicks(allLicks >= faWindowStart & allLicks < faWindowEnd);
        
        if ~isempty(prematureLicks)
            FalseAlarms(t) = true;
        end
    end
    
    % Calculate false alarm rate (excluding discarded trials)
    validTrials = true(1, nTrials);
    if ~isempty(discardedTrials)
        validTrials(discardedTrials) = false;
    end
    
    nValidTrials = sum(validTrials);
    if nValidTrials > 0
        FalseAlarmRate = sum(FalseAlarms) / nValidTrials;
    else
        FalseAlarmRate = 0;
    end
end


function plotFalseAlarmRaster(LicksInFrame, FalseAlarms, nTrials, optoTrials, ... 
    faWindowStart, faWindowEnd)
    % Plot lick raster highlighting false alarm trials
    
    [x, y] = convertLicksToXY(LicksInFrame);
    
    figure('Color', 'w', 'Position', [100 100 800 600]);
    set(0, 'defaultaxesfontname', 'arial');
    hold on;
    
    % Color code licks:  FA trials in red, others in black
    c = repmat([0.1 0.1 0.1], length(x), 1);
    
    scatter(x, y, 20, c, 'filled');
    ylim([0 nTrials]);
    xlim([-5, 5]);
    
    % Highlight false alarm window
    rectangle('Position', [faWindowStart, 0, faWindowEnd-faWindowStart, nTrials], ...
        'FaceColor', [1 0.2 0.2 0.2], 'EdgeColor', [0.8 0.2 0.2], 'LineWidth', 2);
    
    % Highlight opto trials
    if ~isempty(optoTrials)
        xLims = get(gca, 'XLim');
        for t = 1:nTrials
            if optoTrials(t) == 1
                fill([xLims(1) xLims(2) xLims(2) xLims(1)], [t-0.5 t-0.5 t+0.5 t+0.5], ... 
                    [0.6 0.8 1], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            end
        end
    end
    
    % Add vertical line at stimulus onset
    plot([0 0], [0 nTrials], 'k--', 'LineWidth', 2);
    
    xlabel('Time relative to stimulus (s)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Trials', 'FontSize', 14, 'FontWeight', 'bold');
    set(gca, 'YAxisLocation', 'right', 'FontSize', 13);
    title('False Alarm Raster (Red = FA trials)', 'FontSize', 15, 'FontWeight', 'bold'); 
    hold off;
end


function plotFalseAlarmsByCondition(FalseAlarms, optoTrials, isDiscarded)
    % Plot false alarm rates for Control vs Opto conditions
    
    if nargin < 3
        isDiscarded = false(size(optoTrials));
    end
    
    idxOpto = (optoTrials == 1) & ~isDiscarded;
    idxNoOpto = (optoTrials == 0) & ~isDiscarded;
    
    denomNoOpto = sum(idxNoOpto);
    denomOpto = sum(idxOpto);
    
    if denomNoOpto == 0 && denomOpto == 0
        warning('No valid trials for either condition.');
        return;
    end
    
    FA_NoOpto = sum(FalseAlarms(idxNoOpto)) / max(denomNoOpto, 1);
    FA_Opto = sum(FalseAlarms(idxOpto)) / max(denomOpto, 1);
    
    figure('Color', 'w', 'Position', [100 100 500 400]);
    hold on;
    
    b = bar([FA_NoOpto, FA_Opto], 'FaceColor', [0.8 0.3 0.3], ... 
        'EdgeColor', 'k', 'LineWidth', 1.2, 'BarWidth', 0.6);
    
    set(gca, 'XTick', 1:2, 'XTickLabel', {'Control', 'Opto'}, ... 
        'FontSize', 12, 'LineWidth', 1.2, 'Box', 'off');
    ylabel('False Alarm Rate', 'FontSize', 13, 'FontWeight', 'bold');
    
    yMax = min(1, max([FA_NoOpto, FA_Opto]) * 1.3);
    if yMax == 0
        yMax = 0.1;
    end
    ylim([0 yMax]);
    
    grid on;
    set(gca, 'GridAlpha', 0.2);
    title('False Alarm Rate by Condition', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Add text labels on bars
    if FA_NoOpto > 0
        text(1, FA_NoOpto, sprintf('%.1f%%', FA_NoOpto*100), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 11, 'FontWeight', 'bold');
    end
    
    if FA_Opto > 0
        text(2, FA_Opto, sprintf('%.1f%%', FA_Opto*100), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 11, 'FontWeight', 'bold');
    end
    
    % Add sample sizes below bars
    text(1, -0.05*yMax, sprintf('n=%d', denomNoOpto), ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
    text(2, -0.05*yMax, sprintf('n=%d', denomOpto), ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
    
    hold off;
end


function plotFalseAlarmsByConditionAndLocation(FalseAlarms, optoTrials, stimLocations, isDiscarded)
    % Plot FA rates for Control Right/Left and Opto Right/Left
    
    if nargin < 4
        isDiscarded = false(size(optoTrials));
    end
    
    n = numel(optoTrials);
    if numel(FalseAlarms) ~= n || numel(stimLocations) ~= n || numel(isDiscarded) ~= n
        error('All inputs must be same length.');
    end
    
    % Create masks excluding discarded trials
    ctrl_right = (optoTrials == 0) & (stimLocations == 6) & ~isDiscarded;
    ctrl_left  = (optoTrials == 0) & (stimLocations == 7) & ~isDiscarded;
    opto_right = (optoTrials == 1) & (stimLocations == 6) & ~isDiscarded;
    opto_left  = (optoTrials == 1) & (stimLocations == 7) & ~isDiscarded;
    
    masks = {ctrl_right, ctrl_left, opto_right, opto_left};
    FA_rates = zeros(1, 4);
    n_trials = zeros(1, 4);
    
    for k = 1:4
        n_trials(k) = sum(masks{k});
        if n_trials(k) > 0
            FA_rates(k) = sum(FalseAlarms(masks{k})) / n_trials(k);
        end
    end
    
    groups = {'Control Right', 'Control Left', 'Opto Right', 'Opto Left'};
    
    figure('Color', 'w', 'Position', [100 100 700 400]);
    hold on;
    
    colors = [0.7 0.3 0.3; 0.9 0.5 0.5; 0.5 0.2 0.2; 0.8 0.3 0.3];
    b = bar(FA_rates, 'EdgeColor', 'k', 'LineWidth', 1.2, 'BarWidth', 0.6);
    
    % Color each bar differently
    b.FaceColor = 'flat';
    for k = 1:4
        b.CData(k,: ) = colors(k,:);
    end
    
    set(gca, 'XTick', 1:4, 'XTickLabel', groups, 'FontSize', 11, ... 
        'LineWidth', 1.2, 'Box', 'off');
    ylabel('False Alarm Rate', 'FontSize', 13, 'FontWeight', 'bold');
    
    yMax = min(1, max(FA_rates) * 1.3);
    if yMax == 0
        yMax = 0.1;
    end
    ylim([0 yMax]);
    
    grid on;
    set(gca, 'GridAlpha', 0.2);
    title('False Alarm Rate by Condition and Location', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Add text labels
    for k = 1:4
        if FA_rates(k) > 0
            text(k, FA_rates(k), sprintf('%.1f%%', FA_rates(k)*100), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ... 
                'FontSize', 10, 'FontWeight', 'bold');
        end
        text(k, -0.05*yMax, sprintf('n=%d', n_trials(k)), ...
            'HorizontalAlignment', 'center', 'FontSize', 9);
    end
    
    hold off;
end


function [x, y] = convertLicksToXY(LicksInFrame)
    % Convert a LicksInFrame matrix to x (times) and y (trial indices) vectors
    x = [];
    y = [];
    for t = 1:size(LicksInFrame, 2)
        licks = LicksInFrame(LicksInFrame(:,t) ~= 0, t);
        if ~isempty(licks)
            x = [x; licks];
            y = [y; repmat(t, length(licks), 1)];
        end
    end
end
