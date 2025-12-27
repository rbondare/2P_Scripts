%% ANIMAL SUMMARY ANALYSIS
% Variable naming: dataRB##_MMDD_opto or dataRB##_MMDD_baseline
% Gray/Blue = Control/Opto in opto sessions


% PARAMETERS
ANIMAL_ID = 'RB15';  % Change this for different animals
ILI_THRESHOLD = 0.3;
LICK_WINDOW_PRE = -3;
LICK_WINDOW_POST = 2;
HIT_THRESHOLD = 2;
MISS_THRESHOLD = 5;
FA_WINDOW_START = -2;
FA_WINDOW_END = 0;
TIME_BEFORE_STIM = -0.2;
TIME_AFTER_STIM = 0.1;
LICK_WINDOW_PRE_FULL = -16;
LICK_WINDOW_POST_FULL = 10;

% FIND ALL SESSIONS
allVars = who(['data' ANIMAL_ID '*']);
allVars = sort(allVars);
nSessions = length(allVars);

fprintf('Found %d sessions:\n', nSessions);
for i = 1:length(allVars)
    fprintf('  %d. %s\n', i, allVars{i});
end
fprintf('\n');

% COLLECT DATA FROM ALL SESSIONS
allLicks = [];
allLabels = [];
sessionData = struct();

for s = 1:nSessions
    sessionName = allVars{s};
    fprintf('Processing %s...\n', sessionName);
    
    try
        data = eval(sessionName);
        
        % Determine session type
        isBaselineSession = contains(sessionName, 'baseline');
        
        % Extract basic info
        nTrials = data.i;
        StimTime = data.options.Trialtsecbegin;
        LickTimes = str2double(data.lickData(~cellfun('isempty', data.lickData)));
        redFrames = data.redFrames(1, :);
        RedFramesNum = str2double(redFrames) + ((StimTime - data.options.tsecbegin) * 1000)';
        optoTrials = data.options.opto_trial(1:nTrials);
        stimLocations = data.options.stimLocationHistory(1:nTrials);
        
        % Extract licks in full window for trial classification
        [LicksInFrame, discardedTrials] = extractLicksInWindow(LickTimes, RedFramesNum, ...
            LICK_WINDOW_PRE_FULL, LICK_WINDOW_POST_FULL, TIME_BEFORE_STIM, TIME_AFTER_STIM, nTrials);
        
        % Classify trials as Hit/Miss/Unattended
        [isHit, isMiss, isUnattended, isDiscarded] = classifyTrials(LicksInFrame, ...
            discardedTrials, HIT_THRESHOLD, MISS_THRESHOLD, nTrials);
        
        % Calculate false alarms
        [FalseAlarms, ~] = calculateFalseAlarms(LicksInFrame, discardedTrials, ...
            FA_WINDOW_START, FA_WINDOW_END, nTrials);
        
        % Extract date label
        dateParts = regexp(sessionName, 'dataRB\d+_(\d{2})(\d{2})', 'tokens');
        if ~isempty(dateParts)
            month = dateParts{1}{1};
            day = dateParts{1}{2};
            sessionLabel = sprintf('%s/%s', month, day);
        else
            sessionLabel = sprintf('S%d', s);
        end
        
        % Get bout start times for violin plot
        relLicks_control = [];
        relLicks_opto = [];
        
        for t = 1:min(nTrials, length(RedFramesNum))
            limMin = RedFramesNum(t) + LICK_WINDOW_PRE * 1000;
            limMax = RedFramesNum(t) + LICK_WINDOW_POST * 1000;
            
            licksInWindow = LickTimes(LickTimes > limMin & LickTimes < limMax);
            if ~isempty(licksInWindow)
                licksRel = (licksInWindow - RedFramesNum(t)) / 1000;
                licksRel = sort(licksRel);
                boutIdx = [true; diff(licksRel) > ILI_THRESHOLD];
                boutStartTimes = licksRel(boutIdx);
                
                if optoTrials(t) == 0
                    relLicks_control = [relLicks_control; boutStartTimes];
                else
                    relLicks_opto = [relLicks_opto; boutStartTimes];
                end
            end
        end
        
        % Store session data
        sessionData(s).name = sessionName;
        sessionData(s).label = sessionLabel;
        sessionData(s).isBaseline = isBaselineSession;
        sessionData(s).nTrials = nTrials;
        
        if isBaselineSession
            % Baseline session - all trials together
            idx_valid = ~isDiscarded;
            sessionData(s).hitRate = sum(isHit(idx_valid)) / max(sum(idx_valid), 1);
            sessionData(s).missRate = sum(isMiss(idx_valid)) / max(sum(idx_valid), 1);
            sessionData(s).unattendedRate = sum(isUnattended(idx_valid)) / max(sum(idx_valid), 1);
            sessionData(s).FARate = sum(FalseAlarms(idx_valid)) / max(sum(idx_valid), 1);
            
            % By location (baseline)
            idx_right = (stimLocations == 6) & ~isDiscarded;
            idx_left = (stimLocations == 7) & ~isDiscarded;
            sessionData(s).hitRate_right = sum(isHit(idx_right)) / max(sum(idx_right), 1);
            sessionData(s).missRate_right = sum(isMiss(idx_right)) / max(sum(idx_right), 1);
            sessionData(s).unattendedRate_right = sum(isUnattended(idx_right)) / max(sum(idx_right), 1);
            sessionData(s).hitRate_left = sum(isHit(idx_left)) / max(sum(idx_left), 1);
            sessionData(s).missRate_left = sum(isMiss(idx_left)) / max(sum(idx_left), 1);
            sessionData(s).unattendedRate_left = sum(isUnattended(idx_left)) / max(sum(idx_left), 1);
            
            % Add to violin plot
            allLicks = [allLicks; relLicks_control];
            allLabels = [allLabels; repmat({[sessionLabel ' Base']}, length(relLicks_control), 1)];
            
        else
            % Opto session - separate control and opto
            idx_ctrl = (optoTrials == 0) & ~isDiscarded;
            idx_opto = (optoTrials == 1) & ~isDiscarded;
            
            sessionData(s).control.hitRate = sum(isHit(idx_ctrl)) / max(sum(idx_ctrl), 1);
            sessionData(s).control.missRate = sum(isMiss(idx_ctrl)) / max(sum(idx_ctrl), 1);
            sessionData(s).control.unattendedRate = sum(isUnattended(idx_ctrl)) / max(sum(idx_ctrl), 1);
            sessionData(s).control.FARate = sum(FalseAlarms(idx_ctrl)) / max(sum(idx_ctrl), 1);
            
            sessionData(s).opto.hitRate = sum(isHit(idx_opto)) / max(sum(idx_opto), 1);
            sessionData(s).opto.missRate = sum(isMiss(idx_opto)) / max(sum(idx_opto), 1);
            sessionData(s).opto.unattendedRate = sum(isUnattended(idx_opto)) / max(sum(idx_opto), 1);
            sessionData(s).opto.FARate = sum(FalseAlarms(idx_opto)) / max(sum(idx_opto), 1);
            
            % By location - Control
            idx_ctrl_right = (optoTrials == 0) & (stimLocations == 6) & ~isDiscarded;
            idx_ctrl_left = (optoTrials == 0) & (stimLocations == 7) & ~isDiscarded;
            sessionData(s).control.hitRate_right = sum(isHit(idx_ctrl_right)) / max(sum(idx_ctrl_right), 1);
            sessionData(s).control.missRate_right = sum(isMiss(idx_ctrl_right)) / max(sum(idx_ctrl_right), 1);
            sessionData(s).control.unattendedRate_right = sum(isUnattended(idx_ctrl_right)) / max(sum(idx_ctrl_right), 1);
            sessionData(s).control.hitRate_left = sum(isHit(idx_ctrl_left)) / max(sum(idx_ctrl_left), 1);
            sessionData(s).control.missRate_left = sum(isMiss(idx_ctrl_left)) / max(sum(idx_ctrl_left), 1);
            sessionData(s).control.unattendedRate_left = sum(isUnattended(idx_ctrl_left)) / max(sum(idx_ctrl_left), 1);
            
            % By location - Opto
            idx_opto_right = (optoTrials == 1) & (stimLocations == 6) & ~isDiscarded;
            idx_opto_left = (optoTrials == 1) & (stimLocations == 7) & ~isDiscarded;
            sessionData(s).opto.hitRate_right = sum(isHit(idx_opto_right)) / max(sum(idx_opto_right), 1);
            sessionData(s).opto.missRate_right = sum(isMiss(idx_opto_right)) / max(sum(idx_opto_right), 1);
            sessionData(s).opto.unattendedRate_right = sum(isUnattended(idx_opto_right)) / max(sum(idx_opto_right), 1);
            sessionData(s).opto.hitRate_left = sum(isHit(idx_opto_left)) / max(sum(idx_opto_left), 1);
            sessionData(s).opto.missRate_left = sum(isMiss(idx_opto_left)) / max(sum(idx_opto_left), 1);
            sessionData(s).opto.unattendedRate_left = sum(isUnattended(idx_opto_left)) / max(sum(idx_opto_left), 1);
            
            % Add to violin plot
            if ~isempty(relLicks_control)
                allLicks = [allLicks; relLicks_control];
                allLabels = [allLabels; repmat({[sessionLabel ' C']}, length(relLicks_control), 1)];
            end
            if ~isempty(relLicks_opto)
                allLicks = [allLicks; relLicks_opto];
                allLabels = [allLabels; repmat({[sessionLabel ' O']}, length(relLicks_opto), 1)];
            end
        end
        
        fprintf('  Done.\n');
        
    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        sessionData(s).name = sessionName;
        sessionData(s).error = true;
    end
end


%% PLOT 1: VIOLIN PLOT WITH ILI (Bout Start Times)
plotViolinSummary(allLicks, allLabels, ANIMAL_ID, LICK_WINDOW_PRE, LICK_WINDOW_POST);

%% PLOT 2: HIT/MISS/UNATTENDED - OVERALL (Control vs Opto)
plotOutcomesOverall(sessionData, ANIMAL_ID);

%% PLOT 3: AGGREGATED HIT/MISS/UNATTENDED - AVERAGE ACROSS ALL SESSIONS
plotOutcomesAggregated(sessionData, ANIMAL_ID);

%% PLOT 4: FALSE ALARM RATES
plotFalseAlarms(sessionData, ANIMAL_ID);

%% ========== HELPER FUNCTIONS ==========

function [LicksInFrame, discardedTrials] = extractLicksInWindow(LickTimes, RedFramesNum, ...
    windowPre, windowPost, minRT, maxRT, nTrials)
    maxLicks = 100;
    LicksInFrame = zeros(maxLicks, nTrials);
    
    for j = 1:length(RedFramesNum)
        limMin = RedFramesNum(j) + windowPre * 1000;
        limMax = RedFramesNum(j) + windowPost * 1000;
        licksInWindow = LickTimes(LickTimes > limMin & LickTimes < limMax);
        if ~isempty(licksInWindow)
            licksRel = (licksInWindow - RedFramesNum(j)) / 1000;
            LicksInFrame(1:length(licksRel), j) = licksRel;
        end
    end
    
    [x, y] = convertLicksToXY(LicksInFrame);
    discardedTrials = unique(y(x > minRT & x < maxRT));
end

function [x, y] = convertLicksToXY(LicksInFrame)
    x = []; y = [];
    for t = 1:size(LicksInFrame, 2)
        licks = LicksInFrame(LicksInFrame(:,t) ~= 0, t);
        if ~isempty(licks)
            x = [x; licks];
            y = [y; repmat(t, length(licks), 1)];
        end
    end
end

function [isHit, isMiss, isUnattended, isDiscarded] = classifyTrials(LicksInFrame, ...
    discardedTrials, hitThresh, missThresh, nTrials)
    
    isDiscarded = false(1, nTrials);
    if ~isempty(discardedTrials)
        discardedTrials = discardedTrials(discardedTrials >= 1 & discardedTrials <= nTrials);
        isDiscarded(discardedTrials) = true;
    end
    
    LicksHit = NaN(1, nTrials);
    LicksMiss = NaN(1, nTrials);
    
    for t = 1:nTrials
        if isDiscarded(t), continue; end
        
        validLicks = LicksInFrame(LicksInFrame(:,t) > 0, t);
        if isempty(validLicks), continue; end
        
        rt = validLicks(1);
        if rt < hitThresh
            LicksHit(t) = rt;
        elseif rt < missThresh
            LicksMiss(t) = rt;
        end
    end
    
    isHit = ~isnan(LicksHit) & ~isDiscarded;
    isMiss = ~isnan(LicksMiss) & ~isDiscarded;
    isUnattended = ~isHit & ~isMiss & ~isDiscarded;
end

function [FalseAlarms, FalseAlarmRate] = calculateFalseAlarms(LicksInFrame, ...
    discardedTrials, faWindowStart, faWindowEnd, nTrials)
    
    FalseAlarms = false(1, nTrials);
    
    for t = 1:nTrials
        if ~isempty(discardedTrials) && ismember(t, discardedTrials)
            continue
        end
        allLicks = LicksInFrame(LicksInFrame(:,t) ~= 0, t);
        prematureLicks = allLicks(allLicks >= faWindowStart & allLicks < faWindowEnd);
        if ~isempty(prematureLicks)
            FalseAlarms(t) = true;
        end
    end
    
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

%% ========== PLOTTING FUNCTIONS ==========

function plotViolinSummary(allLicks, allLabels, animalID, windowPre, windowPost)
    if isempty(allLicks), return; end
    
    uniqueLabels = unique(allLabels, 'stable');
    nGroups = length(uniqueLabels);
    
    figure('Color', 'w', 'Position', [50 50 max(900, 80*nGroups) 600]);
    hold on;
    
    % Reference lines
    plot([0.5 nGroups+0.5], [-2 -2], 'Color', [0.8 0.2 0.2], 'LineWidth', 2, 'LineStyle', '--');
    plot([0.5 nGroups+0.5], [0 0], 'Color', [0.2 0.2 0.2], 'LineWidth', 2, 'LineStyle', '-');
    
    % Assign colors
    colors = zeros(nGroups, 3);
    for i = 1:nGroups
        if contains(uniqueLabels{i}, ' C')
            colors(i, :) = [0.7 0.7 0.7];  % Gray - Control
        elseif contains(uniqueLabels{i}, ' O')
            colors(i, :) = [0.4 0.6 1];     % Blue - Opto
        elseif contains(uniqueLabels{i}, ' Base')
            colors(i, :) = [0.5 0.8 0.5];   % Green - Baseline
        end
    end
    
    % Violin plot
    vp = violinplot(allLicks, allLabels, 'ShowMean', true, 'ShowData', true, ...
        'ViolinColor', colors);
    
    % Enhance aesthetics
    for i = 1:length(vp)
        if isfield(vp(i), 'MeanPlot') && ~isempty(vp(i).MeanPlot)
            set(vp(i).MeanPlot, 'LineWidth', 3, 'Color', [0 0 0]);
        end
        if isfield(vp(i), 'ScatterPlot') && ~isempty(vp(i).ScatterPlot)
            set(vp(i).ScatterPlot, 'SizeData', 15, 'MarkerFaceAlpha', 0.4);
        end
    end
    
    ylabel('Time from stimulus onset (s)', 'FontSize', 14);
    xlabel('Session (Month/Day)', 'FontSize', 14);
    title(sprintf('%s - Bout Start Times (ILI=0.3s) Across Sessions', animalID), ...
        'FontSize', 16, 'FontWeight', 'bold');
    ylim([windowPre windowPost]);
    set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'Box', 'off');
    xtickangle(45);
    grid on;
    set(gca, 'GridAlpha', 0.2);
    
    % Legend
    h_gray = patch(NaN, NaN, [0.7 0.7 0.7]);
    h_blue = patch(NaN, NaN, [0.4 0.6 1]);
    h_green = patch(NaN, NaN, [0.5 0.8 0.5]);
    legend([h_gray h_blue h_green], {'Control', 'Opto', 'Baseline'}, ...
        'Location', 'northeast', 'FontSize', 10);
    
    text(nGroups+0.8, -2, 'Opto', 'FontSize', 9, 'Color', [0.8 0.2 0.2]);
    text(nGroups+0.8, 0, 'Stim', 'FontSize', 9, 'Color', [0.2 0.2 0.2]);
    
    hold off;
end

function plotOutcomesOverall(sessionData, animalID)
    % Extract opto sessions only
    optoSessions = sessionData(~[sessionData.isBaseline]);
    if isempty(optoSessions), return; end
    
    nSessions = length(optoSessions);
    ctrl_data = zeros(nSessions, 3);
    opto_data = zeros(nSessions, 3);
    
    for i = 1:nSessions
        ctrl_data(i,:) = [optoSessions(i).control.hitRate, ...
                          optoSessions(i).control.missRate, ...
                          optoSessions(i).control.unattendedRate];
        opto_data(i,:) = [optoSessions(i).opto.hitRate, ...
                          optoSessions(i).opto.missRate, ...
                          optoSessions(i).opto.unattendedRate];
    end
    
    colors = [0.3 0.75 0.4; 1.0 0.6 0.2; 0.75 0.75 0.75]; % Hit, Miss, Unattended
    
    figure('Color', 'w', 'Position', [100 100 900 500]);
    
    % Control
    subplot(1, 2, 1);
    b = bar(ctrl_data, 'stacked', 'EdgeColor', 'none');
    for i = 1:3
        b(i).FaceColor = colors(i,:);
    end
    ylabel('Proportion', 'FontSize', 13);
    title('Control Trials', 'FontSize', 14, 'FontWeight', 'bold');
    ylim([0 1]);
    xlabel('Session', 'FontSize', 12);
    set(gca, 'FontSize', 11, 'Box', 'off');
    legend({'Hit', 'Miss', 'Unattended'}, 'Location', 'best');
    grid on;
    
    % Opto
    subplot(1, 2, 2);
    b = bar(opto_data, 'stacked', 'EdgeColor', 'none');
    for i = 1:3
        b(i).FaceColor = colors(i,:);
    end
    ylabel('Proportion', 'FontSize', 13);
    title('Opto Trials', 'FontSize', 14, 'FontWeight', 'bold');
    ylim([0 1]);
    xlabel('Session', 'FontSize', 12);
    set(gca, 'FontSize', 11, 'Box', 'off');
    grid on;
    
    sgtitle(sprintf('%s - Trial Outcomes Across Sessions', animalID), ...
        'FontSize', 15, 'FontWeight', 'bold');
end

function plotOutcomesAggregated(sessionData, animalID)
    % Aggregate and plot average hit/miss/unattended rates across all sessions
    % Shows overall averages and by location (left/right)
    
    % Extract opto sessions only
    optoSessions = sessionData(~[sessionData.isBaseline]);
    if isempty(optoSessions), return; end
    
    nSessions = length(optoSessions);
    
    % Collect rates for averaging
    ctrl_hit = zeros(nSessions, 1);
    ctrl_miss = zeros(nSessions, 1);
    ctrl_unatt = zeros(nSessions, 1);
    
    opto_hit = zeros(nSessions, 1);
    opto_miss = zeros(nSessions, 1);
    opto_unatt = zeros(nSessions, 1);
    
    % By location
    ctrl_hit_right = zeros(nSessions, 1);
    ctrl_miss_right = zeros(nSessions, 1);
    ctrl_unatt_right = zeros(nSessions, 1);
    ctrl_hit_left = zeros(nSessions, 1);
    ctrl_miss_left = zeros(nSessions, 1);
    ctrl_unatt_left = zeros(nSessions, 1);
    
    opto_hit_right = zeros(nSessions, 1);
    opto_miss_right = zeros(nSessions, 1);
    opto_unatt_right = zeros(nSessions, 1);
    opto_hit_left = zeros(nSessions, 1);
    opto_miss_left = zeros(nSessions, 1);
    opto_unatt_left = zeros(nSessions, 1);
    
    for i = 1:nSessions
        % Overall
        ctrl_hit(i) = optoSessions(i).control.hitRate;
        ctrl_miss(i) = optoSessions(i).control.missRate;
        ctrl_unatt(i) = optoSessions(i).control.unattendedRate;
        
        opto_hit(i) = optoSessions(i).opto.hitRate;
        opto_miss(i) = optoSessions(i).opto.missRate;
        opto_unatt(i) = optoSessions(i).opto.unattendedRate;
        
        % By location
        ctrl_hit_right(i) = optoSessions(i).control.hitRate_right;
        ctrl_miss_right(i) = optoSessions(i).control.missRate_right;
        ctrl_unatt_right(i) = optoSessions(i).control.unattendedRate_right;
        ctrl_hit_left(i) = optoSessions(i).control.hitRate_left;
        ctrl_miss_left(i) = optoSessions(i).control.missRate_left;
        ctrl_unatt_left(i) = optoSessions(i).control.unattendedRate_left;
        
        opto_hit_right(i) = optoSessions(i).opto.hitRate_right;
        opto_miss_right(i) = optoSessions(i).opto.missRate_right;
        opto_unatt_right(i) = optoSessions(i).opto.unattendedRate_right;
        opto_hit_left(i) = optoSessions(i).opto.hitRate_left;
        opto_miss_left(i) = optoSessions(i).opto.missRate_left;
        opto_unatt_left(i) = optoSessions(i).opto.unattendedRate_left;
    end
    
    % Calculate means
    ctrl_avg = [mean(ctrl_hit), mean(ctrl_miss), mean(ctrl_unatt)];
    opto_avg = [mean(opto_hit), mean(opto_miss), mean(opto_unatt)];
    
    ctrl_right_avg = [mean(ctrl_hit_right), mean(ctrl_miss_right), mean(ctrl_unatt_right)];
    ctrl_left_avg = [mean(ctrl_hit_left), mean(ctrl_miss_left), mean(ctrl_unatt_left)];
    
    opto_right_avg = [mean(opto_hit_right), mean(opto_miss_right), mean(opto_unatt_right)];
    opto_left_avg = [mean(opto_hit_left), mean(opto_miss_left), mean(opto_unatt_left)];
    
    % Calculate SEMs
    ctrl_sem = [std(ctrl_hit)/sqrt(nSessions), std(ctrl_miss)/sqrt(nSessions), std(ctrl_unatt)/sqrt(nSessions)];
    opto_sem = [std(opto_hit)/sqrt(nSessions), std(opto_miss)/sqrt(nSessions), std(opto_unatt)/sqrt(nSessions)];
    
    colors = [0.3 0.75 0.4; 1.0 0.6 0.2; 0.75 0.75 0.75]; % Hit, Miss, Unattended
    
    figure('Color', 'w', 'Position', [100 100 1200 500]);
    
    % SUBPLOT 1: Overall (Control vs Opto)
    subplot(1, 3, 1);
    data_overall = [ctrl_avg; opto_avg];
    b = bar(data_overall, 'stacked', 'EdgeColor', 'k', 'LineWidth', 1);
    for i = 1:3
        b(i).FaceColor = colors(i,:);
    end
    set(gca, 'XTickLabel', {'Control', 'Opto'}, 'FontSize', 12);
    ylabel('Proportion', 'FontSize', 13);
    title('Overall Average', 'FontSize', 14);
    ylim([0 1]);
    legend({'Hit', 'Miss', 'Unattended'}, 'Location', 'best', 'FontSize', 10);
    grid on;
    set(gca, 'Box', 'off');
    
    % SUBPLOT 2: Right Stimulus (Control vs Opto)
    subplot(1, 3, 2);
    data_right = [ctrl_right_avg; opto_right_avg];
    b = bar(data_right, 'stacked', 'EdgeColor', 'k', 'LineWidth', 1);
    for i = 1:3
        b(i).FaceColor = colors(i,:);
    end
    set(gca, 'XTickLabel', {'Control', 'Opto'}, 'FontSize', 12);
    ylabel('Proportion', 'FontSize', 13);
    title('Right Stimulus', 'FontSize', 14);
    ylim([0 1]);
    grid on;
    set(gca, 'Box', 'off');
    
    % SUBPLOT 3: Left Stimulus (Control vs Opto)
    subplot(1, 3, 3);
    data_left = [ctrl_left_avg; opto_left_avg];
    b = bar(data_left, 'stacked', 'EdgeColor', 'k', 'LineWidth', 1);
    for i = 1:3
        b(i).FaceColor = colors(i,:);
    end
    set(gca, 'XTickLabel', {'Control', 'Opto'}, 'FontSize', 12);
    ylabel('Proportion', 'FontSize', 13);
    title('Left Stimulus', 'FontSize', 14);
    ylim([0 1]);
    grid on;
    set(gca, 'Box', 'off');
    
    sgtitle(sprintf('%s - Average Trial Outcomes (n=%d sessions)', animalID, nSessions), ...
        'FontSize', 15);
end

function plotFalseAlarms(sessionData, animalID)
    % Extract opto sessions only
    optoSessions = sessionData(~[sessionData.isBaseline]);
    if isempty(optoSessions), return; end
    
    nSessions = length(optoSessions);
    ctrl_FA = zeros(nSessions, 1);
    opto_FA = zeros(nSessions, 1);
    
    for i = 1:nSessions
        ctrl_FA(i) = optoSessions(i).control.FARate * 100;
        opto_FA(i) = optoSessions(i).opto.FARate * 100;
    end
    
    figure('Color', 'w', 'Position', [100 100 700 500]);
    
    x = 1:nSessions;
    width = 0.35;
    
    bar(x - width/2, ctrl_FA, width, 'FaceColor', [0.5 0.5 0.5], ...
        'EdgeColor', 'k', 'LineWidth', 1);
    hold on;
    bar(x + width/2, opto_FA, width, 'FaceColor', [0.4 0.6 1], ...
        'EdgeColor', 'k', 'LineWidth', 1);
    
    ylabel('False Alarm Rate (%)', 'FontSize', 14);
    xlabel('Session', 'FontSize', 14);
    title(sprintf('%s - False Alarm Rates Across Sessions', animalID), ...
        'FontSize', 15, 'FontWeight', 'bold');
    legend({'Control', 'Opto'}, 'Location', 'best');
    set(gca, 'FontSize', 12, 'Box', 'off');
    grid on;
    
    hold off;
end