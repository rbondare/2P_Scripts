%% First version of the script to plot and analyse behaviour data

%can create a loop to process several files here
% Load your file
data = load('Z:\Group Members\Rima\Stimulus\AnimalRB15\Training\stim_001__training_20251209_161704');

%{ 
red frames is a windows time, so we substract the trial time start from the stimulus time onset (arduino time) ms * 1000 
to get time in seconds from beggining of trial until the stimulus is shown the first red frame appears at the beginning 
of the trial, so we add that physical time to the red frame to get a reference time (in windows clock) of when the
stimulus starts to align licks to that windows timing 
%} 

%% Lick Plotting (as from Florina's lickPlotting Script)
% physical time of when the stimulus appears

StimTime = data.options.Trialtsecbegin;
ntrials = data.options.trials * data.options.subsessions;
lickData = data.lickData;
redFrames = data.redFrames;
lickWindowMax = data.lickWindowMax;

% Convert lick timestamps from string to numeric, removing empty cells
LickTimes = str2double(lickData(~cellfun('isempty',lickData)));


% If multiple rows of frame data exist, take only first row
if size(redFrames,1) > 1 
    redFrames = redFrames(1,:);
end

RedFramesAsNumber = str2double(redFrames);
RedFramesAsNumber = RedFramesAsNumber + ((StimTime - data.options.tsecbegin) * 1000)';

% Initialize arrays for analysis
Limits = NaN(2,ntrials);
maxLicks = 100; % random upper bound for licks
LicksInFrame = zeros(maxLicks, ntrials);

%%Process licks within time window around each stimulus *(level 3)*
% For each frame, find licks within -16s to +10s window
for j = 1 : length(RedFramesAsNumber)
    % Define time window boundaries (in milliseconds)
    Limits(1,j) = RedFramesAsNumber(j) - 16000; % 16s before stimulus
    Limits(2,j) = RedFramesAsNumber(j) + 10000; % 10s after stimulus
    
    % Find licks within the time window
    [col] = find(LickTimes>Limits(1,j) & LickTimes<Limits(2,j));
    temp = LickTimes(col);
    
    % Convert lick times to seconds relative to stimulus onset
    if ~isempty(temp)
        LicksInFrame(1:length(temp),j) = (temp - RedFramesAsNumber(j))/1000;
    end
end

%%Prepare data for plotting
% Convert lick matrix to x,y coordinates for scatter plot
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

%%Visualisation 
% Set up figure with specific styling
fire = figure;
set(0,'defaultaxesfontname','arial');
grid(gca,'minor');
hAx=gca;
set(hAx,'xminorgrid','off','yminorgrid','off');
hold on

% Plot licks as scatter points
scatter(x,y,20,[0.1 0.1 0.1],'filled');
ylim([0 ntrials]);
xlim([-10, 5]); 

% Add shaded response window
rectangle('Position',[0,0,lickWindowMax,ntrials],'FaceColor',[0.5 .5 .5 .3],'EdgeColor',[0.5 .5 .5 0]);

offBlock = 10;   % number of non-opto trials
onBlock  = 5;    % number of opto trials
blockSize = offBlock + onBlock;

% --- add shaded opto window ---
for startTrial = offBlock + 1 : blockSize : ntrials
    endTrial = min(startTrial + onBlock - 1, ntrials);
    fill([xlim fliplr(xlim)], [startTrial startTrial endTrial endTrial], ...
         [0.6 0.8 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end

% Label axes and set ticks
xlabel({'Time in seconds'},'FontSize',16,'FontWeight','bold');
set(gca, 'YAxisLocation', 'right');
ylabel('Trials', 'FontSize',16,'FontWeight','bold');
yticks([10 20 30 40 50 60 70 80 90 100 110 120]);
xticks([-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3]);
ax=gca;
ax.FontSize = 15;
%scatter(x1,y1, 20 , [0.3 0.3 0.3], 'filled');


%% SEPARATED BY LEVEL (FOR LEVEL 2.5)

levelHistory = options.levelHistory;
idx_validtrials = find(levelHistory == 3);
idx_nonvalidtrials = find(levelHistory == 2);
validTrials.frames = RedFramesAsNumber(idx_validtrials);
nonvalidTrials.frames = RedFramesAsNumber(idx_nonvalidtrials);

%%Prepare data for plotting
x1 = [];
y1 = [];
c = []; % color matrix (will be N x 3)

for t = 1:size(LicksInFrame,2)
    for o = 1:size(LicksInFrame,1)
        if LicksInFrame(o,t) ~= 0
            x1 = [x1, LicksInFrame(o,t)];
            y1 = [y1, t];
            
            if ismember(t, idx_validtrials)
                % valid trial (level 3)
                c = [c; 0.1 0.1 0.1];  % dark gray
            elseif ismember(t, idx_nonvalidtrials)
                % nonvalid trial (level 2)
                c = [c; 1 0.1 0.1];    % red
            else
                % other trials
                c = [c; 0.5 0.5 0.5];  % medium gray
            end
        end
    end
end

fire = figure;
set(0, 'defaultaxesfontname', 'arial');
hAx = gca;
set(hAx, 'xminorgrid', 'off', 'yminorgrid', 'off');
hold on;

% Plot licks as scatter points with per-point colors
scatter(x1, y1, 20, c, 'filled');

% Define limits
numTrials = length(RedFramesAsNumber);
ylim([0, numTrials]);
xlim([-10, 5]);

% Add shaded response window
rectangle('Position', [0, 0, lickWindowMax, numTrials], ...
    'FaceColor', [0.5, 0.5, 0.5, 0.3], 'EdgeColor', 'none');

% Label axes and set ticks
xlabel({'Time in seconds'}, 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Trials', 'FontSize', 16, 'FontWeight', 'bold');
set(gca, 'YAxisLocation', 'right');
yticks(0:10:numTrials);
xticks(-10:1:3);
ax = gca;
ax.FontSize = 15;


%% groups into hit, miss and unattended 

time_before_stim = -0.2; %disregard trials where there was a lick before this time (0 is the stimulus)
time_after_stim =0.1; %disregard trials with reaction time faster than 0.1s (not physiological)

discardedTrials = [];
for trials = 1:length(x)
    if x(trials) > time_before_stim & x(trials) < time_after_stim;
        discardedTrials = [discardedTrials, y(trials)];
    end
end
discardedTrials = unique(discardedTrials);
discardedTrials

ntrials = size(LicksInFrame, 2);     % number of trials
maxLicks = size(LicksInFrame, 1);    % max number of licks per trial

% Preallocate matrix to hold valid licks
AllValidLicks = NaN(maxLicks, ntrials);  

for t = 1:ntrials
    if ismember(t, discardedTrials)
        continue
    end

    % Get all licks for this trial
    trial_licks = LicksInFrame(:, t);

    % Keep only positive, non-zero licks
    valid_licks = trial_licks(trial_licks > 0);

    % Fill the column with valid licks (top-aligned)
    n = length(valid_licks);
    if n > 0
        AllValidLicks(1:n, t) = valid_licks;
    end
end

% Preallocate arrays
ntrials = size(LicksInFrame, 2);  % number of columns = trials
LicksHit        = NaN(1, ntrials);
LicksMiss       = NaN(1, ntrials);
LicksUnattended = NaN(1, ntrials);

for t = 1:ntrials
    % Skip discarded trials
    if ismember(t, discardedTrials)
        continue
    end

    % Get all licks for this trial
    trial_licks = LicksInFrame(:, t);

    % Keep only positive, non-zero licks
    valid_licks = trial_licks(trial_licks > 0);

    % Skip if no valid licks
    if isempty(valid_licks)
        continue
    end

    % Take the first valid lick
    rt = valid_licks(1);

    % Classify trial based on first lick
    if rt < 2
        LicksHit(t) = rt;                % Hit
    elseif rt < 5
        LicksMiss(t) = rt;               % Miss
    else
        LicksUnattended(t) = rt;         % Unattended
    end
end

% Count trials in each category
corrected_trials =  ntrials -length(discardedTrials);
numHits        = sum(~isnan(LicksHit));
numMisses      = sum(~isnan(LicksMiss));
numUnattended  = corrected_trials - (numMisses + numHits);

fprintf('Hits: %d, Misses: %d, Unattended: %d\n', numHits, numMisses, numUnattended);
HitRate = numHits / ntrials *100

%% sepapare hit/miss/unattended based on OPTO

opto_trials = data.options.opto_trial;

% --- Select trials by condition and remove NaNs ---
idx_opto    = find(opto_trials == 1);   % indices where opto is 1
idx_nonopto = find(opto_trials == 0);   % indices where opto is 0

isUnattended = isnan(LicksHit) & isnan(LicksMiss);

% --- Count for each condition ---
counts_nonopto = [ ...
    sum(~isnan(LicksHit(idx_nonopto))), ...
    sum(~isnan(LicksMiss(idx_nonopto))), ...
    sum(isUnattended(idx_nonopto))];

counts_opto = [ ...
    sum(~isnan(LicksHit(idx_opto))), ...
    sum(~isnan(LicksMiss(idx_opto))), ...
    sum(isUnattended(idx_opto))];

% --- Normalize to total number of valid trials in each condition ---
n_nonopto = numel(idx_nonopto);
n_opto = numel(idx_opto);

prop_nonopto = counts_nonopto / n_nonopto;
prop_opto    = counts_opto / n_opto;

% --- Combine for plotting ---
combined = [prop_nonopto; prop_opto];   % rows = conditions, cols = Hit/Miss/Unattended

% --- Plot stacked (filled) bar chart ---
figure('Color', 'w'); hold on;

% Nice, soft colors (Hit = green, Miss = amber, Unattended = light gray)
colors = [0.3 0.75 0.4;    % Hit
          1.0 0.6 0.2;     % Miss
          0.75 0.75 0.75]; % Unattended

b = bar(combined, 'stacked', 'BarWidth', 0.6, 'EdgeColor', 'none');
for i = 1:numel(b)
    b(i).FaceColor = colors(i,:);
end

% Aesthetics
set(gca, 'XTick', 1:2, ...
         'XTickLabel', {'No-Stim', 'Stim'}, ...
         'FontSize', 12, ...
         'LineWidth', 1.2, ...
         'Box', 'off');

ylabel('Proportion of Trials', 'FontSize', 13);
title('Trial outcome distribution by condition', 'FontSize', 14, 'FontWeight', 'bold');

ylim([0 1]);
yticks(0:0.2:1);
ax = gca;
ax.YGrid = 'on';
ax.GridAlpha = 0.2;

legend({'Hit', 'Miss', 'Unattended'}, 'Location', 'northeastoutside', ...
    'FontSize', 11, 'Box', 'off');


%% SPLIT HIT/MISS/UNATTENDED BY OPTO AND |LOCATION|

% --- Identify conditions ---
opto_trials = data.options.opto_trial;            % 1=Opto, 0=Non-opto
stim_locations = data.options.stimLocationHistory; % e.g., 6=Right, 7=Left

% Identify outcomes
isHit = ~isnan(LicksHit);
isMiss = ~isnan(LicksMiss);
isUnattended = isnan(LicksHit) & isnan(LicksMiss);

% Define groups
groups = {'No-opto Right', 'No-opto Left', 'Opto Right', 'Opto Left'};

data = zeros(4,3); % rows = groups, cols = Hit/Miss/Unattended

% Loop over each group and count outcomes
for i = 1:4
    switch i
        case 1
            idx = (opto_trials==0 & stim_locations==6); % Non-opto Right
        case 2
            idx = (opto_trials==0 & stim_locations==7); % Non-opto Left
        case 3
            idx = (opto_trials==1 & stim_locations==6); % Opto Right
        case 4
            idx = (opto_trials==1 & stim_locations==7); % Opto Left
    end
    data(i,1) = sum(isHit(idx));
    data(i,2) = sum(isMiss(idx));
    data(i,3) = sum(isUnattended(idx));
    
    % Normalize to total trials in that group
    n = sum(idx);
    if n>0
        data(i,:) = data(i,:) / n;
    end
end

% --- Plot stacked bar chart ---
figure('Color','w'); hold on;
colors = [0.3 0.75 0.4; 1.0 0.6 0.2; 0.75 0.75 0.75]; % Hit, Miss, Unattended

b = bar(data, 'stacked', 'BarWidth', 0.6, 'EdgeColor', 'none');
for i = 1:numel(b)
    b(i).FaceColor = colors(i,:);
end

% --- Aesthetics ---
set(gca, 'XTick', 1:4, 'XTickLabel', groups, 'FontSize', 12, 'LineWidth', 1.2, 'Box', 'off');
ylabel('Proportion of Trials', 'FontSize', 13);
ylim([0 1]);
yticks(0:0.2:1);
ax = gca;
ax.YGrid = 'on';
ax.GridAlpha = 0.2;

legend({'Hit','Miss','Unattended'}, 'Location','northeastoutside','FontSize',11,'Box','off');

%% plot reaction time for opto and split by location
% FirstLicks      : vector of first lick times per trial (NaN for invalid trials)
% Stim_locations  : vector of stimulus location per trial (e.g., 6=Left,7=Right)
% opto_trials     : vector, 0 = non-opto, 1 = opto

Stim_locations = options.stimLocationHistory(:);   % ensure column vector
FirstLicks_column = FirstLicks(:);                % column vector of first licks
opto_trials = options.opto_trial(:);             % 1 = opto, 0 = non-opto

locCodes  = [6, 1, 7];                           % Left, Center, Right
locLabels = ["Left", "Center", "Right"];

condCodes  = [0, 1];                             % Non-opto = 0, Opto = 1
condLabels = ["Non-opto", "Opto"];

allRT    = [];
allGroup = [];
allCond  = [];

% Map numeric codes to labels
locMap = containers.Map([6,7], {'Left','Right'});
locLabels = cellfun(@(x) locMap(x), num2cell(Stim_locations), 'UniformOutput', false);


%% Plot reaction times for each subsession

figure('Color','w'); hold on;

allRT = [];
allGroup = [];

for i = 1:ntrials
    varName = sprintf('FirstLicks_subsession%d', i);
    if evalin('base', sprintf('exist("%s", ''var'')', varName))
        data = evalin('base', varName);
        data = data(~isnan(data));
        
        allRT = [allRT; data(:)];                 % Reaction times
        allGroup = [allGroup; repmat(i, numel(data), 1)];  % Numeric subsession index
    end
end

% --- Boxplot ---
%boxplot(allRT, allGroup, 'Colors', [0 0 1], 'Whisker', 1.5, 'Symbol', ''); % 'Symbol' empty to hide default points

% --- Scatter individual data points on top ---
uniqueGroups = unique(allGroup);
for i = 1:numel(uniqueGroups)
    x = allGroup(allGroup == uniqueGroups(i));
    y = allRT(allGroup == uniqueGroups(i));
    jitter = (rand(size(y)) - 0.5) * 0.2;  % small horizontal jitter
    scatter(x + jitter, y, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
end

vp = violinplot(allRT, allGroup, 'ShowMean', true, 'ViolinColor', repmat([0.6 0.6 1], nSubsessions, 1));

% --- Aesthetics ---
xlim([0.5, nSubsessions + 0.5]);
xticks(1:nSubsessions);
xticklabels(arrayfun(@(x)sprintf('Subsession %d', x), 1:nSubsessions, 'UniformOutput', false));
ylabel('Reaction time (s)');
title('Reaction time distribution per subsession');
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'Box', 'on');

%% calculates hit & miss rate per subsession (for when opto is subsession based)

nSubsessions = options.subsessions;
ntrials = options.trials;

% Preallocate
HitRate = nan(1, nSubsessions);
MissRate = nan(1, nSubsessions);


for i = 1 : nSubsessions
    start_trial = (i-1)*nTrials + 1;
    end_trial = i* nTrials; 
    FirstLicks_subsession = FirstLicks(1, start_trial:end_trial);

    lick_var_name = sprintf('FirstLicks_subsession%d', i); 
    assignin('base', lick_var_name, FirstLicks_subsession);

    nHits = sum(~isnan(FirstLicks_subsession));
    nMiss = sum(isnan(FirstLicks_subsession));
    total = nHits + nMiss;
    HitRate(i) = nHits / total;
    MissRate(i) = nMiss / total;
    fprintf('Subsession %d: Hits = %d, Misses = %d, HitRate = %.2f\n', ...
    i, nHits, nMiss, HitRate(i));

end 

%% plot reaction time based on location (all trials, location separated)
Stim_locations = options.stimLocationHistory;  % stimulus location per trial
FirstLicks_column = FirstLicks(:);                    % ensure column vector

center_idx = find(Stim_locations == 1);
left_idx   = find(Stim_locations == 6);
right_idx  = find(Stim_locations == 7);

RT_center = FirstLicks_column(center_idx);
RT_left   = FirstLicks_column(left_idx);
RT_right  = FirstLicks_column(right_idx);

allRT = [RT_left; RT_center; RT_right];
allGroup = [repmat({'Left'}, numel(RT_left), 1); ...
            repmat({'Center'}, numel(RT_center), 1); ...
            repmat({'Right'}, numel(RT_right), 1)];

% --- Boxplot ---
%figure('Color','w'); hold on;
%boxplot(allRT, allGroup, 'Colors', [0 0 1], 'Whisker', 1.5, 'Symbol', '');  % hide default points

figure; hold on;


vp = violinplot(allRT, allGroup, 'ShowMean', true, 'ViolinColor', repmat([0.6 0.6 1], 3, 1), 'ShowData', true);

ylabel('Reaction time (s)');
title('Reaction time by stimulus location');
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'Box', 'on');

%% plot reaction times by subsession *and* location
figure('Color','w'); hold on;

Stim_locations = options.stimLocationHistory;  % stimulus location per trial
FirstLicks_column = FirstLicks(:);                    % ensure column vector

locCodes = [6, 1, 7];
locLabels = {'Left', 'Center', 'Right'};

allRT = [];
allGroup = [];
allSub = [];

for li = 1:numel(locCodes)
    % indices for this stimulus location
    idx_loc = find(Stim_locations == locCodes(li));
    RT_loc = FirstLicks_column(idx_loc);
    RT_loc = RT_loc(~isnan(RT_loc));  % remove NaNs

    % split by subsession
    for s = 1:nSubsessions
        startIdx = (s-1)*nTrials + 1;
        endIdx = min(s*nTrials, numel(Stim_locations));
        subMask = idx_loc >= startIdx & idx_loc <= endIdx;

        RT_sub = FirstLicks_column(idx_loc(subMask));
        RT_sub = RT_sub(~isnan(RT_sub));

        if ~isempty(RT_sub)
            allRT = [allRT; RT_sub];
            allGroup = [allGroup; repmat(locLabels(li), numel(RT_sub), 1)];
            allSub = [allSub; repmat(s, numel(RT_sub), 1)];
        end
    end
end

% --- Grouped boxplot ---
boxplot(allRT, {allGroup, allSub}, 'Colors', [0 0 1], ...
        'Whisker', 1.5, 'Symbol', '', 'FactorSeparator', 1);

% --- Scatter individual data points ---
uniqueLoc = unique(allGroup, 'stable');
uniqueSub = unique(allSub);

positions = [];
count = 1;
for li = 1:numel(uniqueLoc)
    for si = 1:numel(uniqueSub)
        idx = strcmp(allGroup, uniqueLoc{li}) & (allSub == uniqueSub(si));
        y = allRT(idx);
        if ~isempty(y)
            xpos = count;
            jitter = (rand(size(y)) - 0.5) * 0.2;
            scatter(xpos + jitter, y, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
        end
        positions = [positions; count];
        count = count + 1;
    end
end

% --- Adjust x-axis labels ---
xticks(mean(reshape(positions, numel(uniqueSub), []).', 2));
xticklabels(uniqueLoc);
xlabel('Stimulus location');
ylabel('Reaction time (s)');
title('Reaction time by location and subsession');
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'Box', 'on');

% --- Violinplot (commented out for when available) ---
% vp = violinplot(allRT, strcat(allGroup, "_S", string(allSub)), ...
%     'ShowMean', true, 'ViolinColor', repmat([0.6 0.6 1], numel(uniqueLoc)*numel(uniqueSub), 1), 'ShowData', true);



%% plot timeouts (erranous licks) seconds correspond to number of timeouts? 

%% old scripts for reaction time calculations

%% Lick Plotting (as from Florina's lickPlotting Script)
% physical time of when the stimulus appears

StimTime = data.options.Trialtsecbegin;
ntrials = data.options.trials * data.options.subsessions;
lickData = data.lickData;
redFrames = data.redFrames;
lickWindowMax = data.lickWindowMax;
opto_trials = data.options.opto_trial;

% Convert lick timestamps from string to numeric, removing empty cells
LickTimes = str2double(lickData(~cellfun('isempty',lickData)));


% If multiple rows of frame data exist, take only first row
if size(redFrames,1) > 1 
    redFrames = redFrames(1,:);
end

RedFramesAsNumber = str2double(redFrames);
RedFramesAsNumber = RedFramesAsNumber + ((StimTime - data.options.tsecbegin) * 1000)';

% Initialize arrays for analysis
Limits = NaN(2,ntrials);
maxLicks = 100; % random upper bound for licks
LicksInFrame = zeros(maxLicks, ntrials);

%%Process licks within time window around each stimulus *(level 3)*
% For each frame, find licks within -16s to +10s window
for j = 1 : length(RedFramesAsNumber)
    % Define time window boundaries (in milliseconds)
    Limits(1,j) = RedFramesAsNumber(j) - 16000; % 16s before stimulus
    Limits(2,j) = RedFramesAsNumber(j) + 10000; % 10s after stimulus
    
    % Find licks within the time window
    [col] = find(LickTimes>Limits(1,j) & LickTimes<Limits(2,j));
    temp = LickTimes(col);
    
    % Convert lick times to seconds relative to stimulus onset
    if ~isempty(temp)
        LicksInFrame(1:length(temp),j) = (temp - RedFramesAsNumber(j))/1000;
    end
end

%%Prepare data for plotting
% Convert lick matrix to x,y coordinates for scatter plot
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

%%Visualisation 
% Set up figure with specific styling
fire = figure;
set(0,'defaultaxesfontname','arial');
grid(gca,'minor');
hAx=gca;
set(hAx,'xminorgrid','off','yminorgrid','off');
hold on

% Plot licks as scatter points
scatter(x,y,20,[0.1 0.1 0.1],'filled');
ylim([0 ntrials]);
xlim([-10, 5]); 

% Add shaded response window
rectangle('Position',[0,0,lickWindowMax,ntrials],'FaceColor',[0.5 .5 .5 .3],'EdgeColor',[0.5 .5 .5 0]);

offBlock = 10;   % number of non-opto trials
onBlock  = 5;    % number of opto trials
blockSize = offBlock + onBlock;

% --- add shaded opto window ---
for startTrial = offBlock + 1 : blockSize : ntrials
    endTrial = min(startTrial + onBlock - 1, ntrials);
    fill([xlim fliplr(xlim)], [startTrial startTrial endTrial endTrial], ...
         [0.6 0.8 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end

% Label axes and set ticks
xlabel({'Time in seconds'},'FontSize',16,'FontWeight','bold');
set(gca, 'YAxisLocation', 'right');
ylabel('Trials', 'FontSize',16,'FontWeight','bold');
yticks([10 20 30 40 50 60 70 80 90 100 110 120]);
xticks([-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3]);
ax=gca;
ax.FontSize = 15;
%scatter(x1,y1, 20 , [0.3 0.3 0.3], 'filled');

%% plot ALL reaction times
% --- Parameters ---
preTime  = -2;    % start of window (s)
postTime = 3;     % end of window (s)

% --- Collect all lick times relative to stimulus ---
relLicks = [];
groupTag = [];

for t = 1:ntrials
    lickTimes = LicksInFrame(:, t);
    lickTimes = lickTimes(~isnan(lickTimes) & lickTimes ~= 0);  % remove NaN and exact 0

    % keep only within -2 to 5 s
    lickTimes = lickTimes(lickTimes >= preTime & lickTimes <= postTime);

    % store data
    relLicks = [relLicks; lickTimes];
    groupTag = [groupTag; repmat(opto_trials(t), numel(lickTimes), 1)];
end

% --- Split by condition ---
licks_nonopto = relLicks(groupTag == 0);
licks_opto    = relLicks(groupTag == 1);

% --- Plot violin plot ---
figure('Color','w','Position',[200 200 400 500]); hold on;
colors = [0.6 0.6 0.6; 0.2 0.4 1];  % gray = non-opto, blue = opto

vp = violinplot({licks_nonopto, licks_opto}, {'Non-opto','Opto'}, ...
    'ViolinColor', colors, 'ShowMean', true, 'ShowData', true, 'MedianColor', [0 0 0]);

yline(0, '--k', 'Stimulus');  % stimulus onset
xlabel('Condition');
ylabel('Lick time relative to stimulus (s)');
ylim([preTime postTime]);

set(gca, 'Box','off', 'TickDir','out', 'FontSize',11, 'LineWidth',1.2);
title('Reaction times by condition');



%% plot reaction times split by OPTO 

opto_trials = data.options.opto_trial;


% time_before_stim = -0.2; %disregard trials where there was a lick before this time (0 is the stimulus)
% time_after_stim =0.1; %disregard trials with reaction time faster than 0.1s (not physiological)
% 
% discardedTrials = [];
% for trials = 1:length(x)
%     if x(trials) > time_before_stim & x(trials) < time_after_stim;
%         discardedTrials = [discardedTrials, y(trials)];
%     end
% end
% discardedTrials = unique(discardedTrials);
% discardedTrials


% --- Build matrix of lick times per trial ---
well = zeros(max(y), 1);
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
FirstLicks = NaN(size(well,1),1);
for t = 1:size(well,1)
    if ~ismember(t, discardedTrials) && well(t,1) ~= 0 && well(t,1) < lickWindowMax
        FirstLicks(t) = well(t,1);
    end
end

% --- Split by opto condition ---
rt_opto    = FirstLicks(opto_trials == 1);
rt_nonopto = FirstLicks(opto_trials == 0);

% --- Remove NaNs and zeros ---
rt_opto    = rt_opto(~isnan(rt_opto) & rt_opto ~= 0);
rt_nonopto = rt_nonopto(~isnan(rt_nonopto) & rt_nonopto ~= 0);

% --- Combine for plotting ---
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

