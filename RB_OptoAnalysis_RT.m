%% Reaction Time Analysis - Opto and Location

% Assuming there are:
% FirstLicks - vector of first lick times per trial
% optoTrials - vector where 1=opto ON, 0=opto OFF
% stimLocations - vector where 6=Right, 7=Left
%% Detect session type
hasOpto = any(optoTrials == 1);

% PLOT 1: Reaction Time by Condition
if hasOpto
    plotReactionTimeByOpto(LicksHit, optoTrials);
end

% PLOT 2: Reaction Time by Location
plotReactionTimeByLocation(LicksHit, stimLocations);

% PLOT 3: Reaction Time by Condition and Location
if hasOpto
    plotReactionTimeByOptoAndLocation(LicksHit, optoTrials, stimLocations);
end

%% VIOLIN PLOT OF ALL LICKS 
LICK_WINDOW_PRE_OPTO = -3;
LICK_WINDOW_POST_OPTO = 5;
relLicks = [];
groupTag = [];

for t = 1:nTrials
    lickTimes = LicksInFrame(:, t);
    lickTimes = lickTimes(~isnan(lickTimes) & lickTimes ~= 0);
    
    % Filter to keep only licks between -3s and 5s
    lickTimes = lickTimes(lickTimes >= LICK_WINDOW_PRE_OPTO & lickTimes <= LICK_WINDOW_POST_OPTO);
    
    % Each lick gets tagged with this trial's opto condition
    relLicks = [relLicks; lickTimes];
    groupTag = [groupTag; repmat(optoTrials(t), numel(lickTimes), 1)];
end


% Create group labels for plotting
groupLabels = cell(size(groupTag));
groupLabels(groupTag == 0) = {'Control'};
groupLabels(groupTag == 1) = {'Opto'};

% Plot violin
figure('Color', 'w', 'Position', [100 100 500 600]);
hold on;

% Reference lines
plot([0.5 2.5], [-2 -2], 'Color', [0.8 0.2 0.2], 'LineWidth', 2, 'LineStyle', '--');
plot([0.5 2.5], [0 0], 'Color', [0.2 0.2 0.2], 'LineWidth', 2, 'LineStyle', '-');

% Violin plot with custom colors
colors = [0.7 0.7 0.7; 0.4 0.6 1];
vp = violinplot(relLicks, groupLabels, 'ShowMean', true, 'ShowData', true, ...
    'ViolinColor', colors);

% Enhance mean lines
for i = 1:length(vp)
    if isfield(vp(i), 'MeanPlot') && ~isempty(vp(i).MeanPlot)
        set(vp(i).MeanPlot, 'LineWidth', 3, 'Color', [0 0 0]);
    end
    % Make scatter points smaller and more transparent
    if isfield(vp(i), 'ScatterPlot') && ~isempty(vp(i).ScatterPlot)
        set(vp(i).ScatterPlot, 'SizeData', 15, 'MarkerFaceAlpha', 0.4);
    end
end

% Labels and aesthetics
ylabel('Time from stimulus onset (s)', 'FontSize', 14);
set(gca, 'FontSize', 13, 'LineWidth', 1.5, 'Box', 'off');
ylim([LICK_WINDOW_PRE_OPTO LICK_WINDOW_POST_OPTO]);

% Add text labels for reference lines
text(2.6, -2, 'Opto', 'FontSize', 9, 'Color', [0.8 0.2 0.2]);
text(2.6, 0, 'Stim', 'FontSize', 9, 'Color', [0.2 0.2 0.2]);

% Find median for Opto trials between -2 and 0 seconds
idx_opto = strcmp(groupLabels, 'Opto');  % Find Opto trials
opto_licks = relLicks(idx_opto);         % Get licks from Opto trials

% Filter to window between -2 and 0 seconds
opto_licks_window = opto_licks(opto_licks >= -2 & opto_licks <= 0);

% Calculate median
median_opto_window = median(opto_licks_window);

fprintf('Median lick time for Opto: %.3f s\n', 2 + median_opto_window);

% % for Control
% idx_control = strcmp(groupLabels, 'Control');
% control_licks = relLicks(idx_control);
% control_licks_window = control_licks(control_licks >= -2 & control_licks <= 0);
% median_control_window = median(control_licks_window);
% fprintf('Median lick time for Control (between -2 and 0s): %.3f s\n', median_control_window);

%% VIOLIN PLOT WITH INTERLICK PERIOD (REACTION TIMES ONLY)

LICK_WINDOW_PRE_OPTO  = -3;
LICK_WINDOW_POST_OPTO = 2;
INTERLICK_THRESH = 0.3;  % 300 ms

relLicks = [];
groupTag = [];

for t = 1:nTrials

    lickTimes = LicksInFrame(:, t);
    lickTimes = lickTimes(~isnan(lickTimes) & lickTimes ~= 0);
    lickTimes = lickTimes(lickTimes >= LICK_WINDOW_PRE_OPTO & ...
                          lickTimes <= LICK_WINDOW_POST_OPTO);

    if isempty(lickTimes)
        continue
    end

    lickTimes = sort(lickTimes);

    % Identify bout starts
    isBoutStart = [true; diff(lickTimes) > INTERLICK_THRESH];
    boutStartTimes = lickTimes(isBoutStart);

    relLicks = [relLicks; boutStartTimes];
    groupTag = [groupTag; repmat(optoTrials(t), numel(boutStartTimes), 1)];

end


groupLabels = cell(size(groupTag));
groupLabels(groupTag == 0) = {'Control'};
groupLabels(groupTag == 1) = {'Opto'};

figure('Color', 'w', 'Position', [100 100 500 600]); hold on;

plot([0.5 2.5], [-2 -2], 'Color', [0.8 0.2 0.2], 'LineWidth', 2, 'LineStyle', '--');
plot([0.5 2.5], [0 0], 'Color', [0.2 0.2 0.2], 'LineWidth', 2);

colors = [0.7 0.7 0.7; 0.4 0.6 1];
vp = violinplot(relLicks, groupLabels, ...
    'ShowMean', true, 'ShowData', true, ...
    'ViolinColor', colors);

ylabel('Reaction time from stimulus (s)', 'FontSize', 14);
ylim([LICK_WINDOW_PRE_OPTO LICK_WINDOW_POST_OPTO]);

%% plot reaction times split by OPTO 
x = [];
y = [];
for t = 1 : size(LicksInFrame,2)
    for o = 1 : size(LicksInFrame,1)
        if LicksInFrame(o,t) ~= 0
            x = [x, LicksInFrame(o,t)];
            y = [y, t];
        end
    end
end

% --- Build matrix of lick times per trial ---
% FIX: Use nTrials instead of max(y)
well = zeros(nTrials, 100);  % preallocate with max possible licks
rownum = min(y);
colnum = 1;

for i = 1:length(x)
    if rownum ~= y(i)
        colnum = 1;
    end
    if x(i) > 0.1
        well(y(i), colnum) = x(i);
        colnum = colnum + 1;
    end
    rownum = y(i);
end

% --- Get first valid lick per trial ---
FirstLicks = NaN(nTrials, 1);  % FIX: Use nTrials

for t = 1:nTrials
    if ~ismember(t, discardedTrials) && well(t,1) ~= 0 && well(t,1) < lickWindowMax
        FirstLicks(t) = well(t,1);
    end
end

% --- Split by opto condition ---
rt_opto    = FirstLicks(optoTrials == 1);
rt_nonopto = FirstLicks(optoTrials == 0);

% --- Remove NaNs and zeros ---
rt_opto    = rt_opto(~isnan(rt_opto) & rt_opto ~= 0);
rt_nonopto = rt_nonopto(~isnan(rt_nonopto) & rt_nonopto ~= 0);

allRT = [rt_nonopto; rt_opto];
allGroup = [repmat("Non-opto", numel(rt_nonopto), 1);
            repmat("Opto",     numel(rt_opto),     1)];
allGroup = categorical(allGroup, ["Non-opto","Opto"]);

% --- Plot violin plot ---
figure('Color','w','Position',[200 200 400 500]); hold on;
colors = [0.7 0.7 0.7; 0.4 0.6 1];  % gray = Non-opto, blue = Opto
vp = violinplot(allRT, allGroup, ...
    'ShowMean', true, ...
    'ViolinColor', colors);

set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.2, ...
         'FontSize', 11, 'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2]);

ylim([min(allRT)*-0.2, max(allRT)*1.3]);
xlabel('Condition');
ylabel('Reaction time (s)');
box on;

% --- Prepare figure ---
figure('Color','w','Position',[200 200 500 400]); hold on;

% --- Define bin edges ---
nbins=10;  % adjust bin width as needed

% --- Plot histograms ---
h1 = histogram(rt_nonopto, nbins, 'Normalization', 'probability', ...
    'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h2 = histogram(rt_opto, nbins, 'Normalization', 'probability', ...
    'FaceColor', [0.4 0.6 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% --- Format axes ---
xlabel('Reaction time (s)');
ylabel('Probability');
xlim([0 lickWindowMax]);
ylim([0 max([h1.Values h2.Values])*1.2]);

set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.2, ...
         'FontSize', 11, 'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2]);

legend({'Non-opto','Opto'}, 'Location', 'northeast');
title('Reaction time distribution by condition');

%% FUNCTIONS

function plotReactionTimeByOpto(FirstLicks, optoTrials)
    
    allRT = [];
    allGroup = [];
    
    for opto = 0:1
        idx = optoTrials == opto;
        rt = FirstLicks(idx);
        rt = rt(~isnan(rt));
        
        if ~isempty(rt)
            optoLabel = {'Opto OFF', 'Opto ON'};
            allRT = [allRT; rt(:)];
            allGroup = [allGroup; repmat(optoLabel(opto+1), numel(rt), 1)];
        end
    end
    
    figure('Color', 'w');
    hold on;
    
    violinplot(allRT, allGroup, 'ShowMean', true, 'ShowData', true);
    
    ylabel('Reaction Time (s)', 'FontSize', 13);
    title('Reaction Time by Optogenetic Condition', 'FontSize', 14, 'FontWeight', 'bold');
    set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'Box', 'on');
end

function plotReactionTimeByLocation(FirstLicks, stimLocations)
    
    locCodes = [6, 7];
    locLabels = {'Right', 'Left'};
    
    allRT = [];
    allGroup = [];
    
    for i = 1:length(locCodes)
        idx = stimLocations == locCodes(i);
        rt = FirstLicks(idx);
        rt = rt(~isnan(rt));
        
        if ~isempty(rt)
            allRT = [allRT; rt(:)];
            allGroup = [allGroup; repmat(locLabels(i), numel(rt), 1)];
        end
    end
    
    figure('Color', 'w');
    hold on;
    
    violinplot(allRT, allGroup, 'ShowMean', true, 'ShowData', true);
    
    ylabel('Reaction Time (s)', 'FontSize', 13);
    title('Reaction Time by Stimulus Location', 'FontSize', 14, 'FontWeight', 'bold');
    set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'Box', 'on');
end

function plotReactionTimeByOptoAndLocation(FirstLicks, optoTrials, stimLocations)
    
    locCodes = [6, 7];
    locLabels = {'Right', 'Left'};
    
    allRT = [];
    allGroup = [];
    
    for opto = 0:1
        for i = 1:length(locCodes)
            idx = (optoTrials == opto) & (stimLocations == locCodes(i));
            rt = FirstLicks(idx);
            rt = rt(~isnan(rt));
            
            if ~isempty(rt)
                optoLabel = {'OFF', 'ON'};
                groupLabel = sprintf('%s Opto %s', locLabels{i}, optoLabel{opto+1});
                
                allRT = [allRT; rt(:)];
                allGroup = [allGroup; repmat({groupLabel}, numel(rt), 1)];
            end
        end
    end
    
    figure('Color', 'w');
    hold on;
    
    violinplot(allRT, allGroup, 'ShowMean', true, 'ShowData', true);
    
    ylabel('Reaction Time (s)', 'FontSize', 13);
    title('Reaction Time by Location and Optogenetic Condition', 'FontSize', 14, 'FontWeight', 'bold');
    set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'Box', 'on');
    xtickangle(45);
end
