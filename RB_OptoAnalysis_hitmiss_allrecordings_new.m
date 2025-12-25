%% PLOTTING ALL VIOLIN PLOTS OF REACTION TIME ACROSS SESSIONS
% Loops through all dataRB* variables and creates one violin plot
% NOW WITH BOUT DETECTION (ILI = 0.3s)

% Variable naming: dataRB##_MMDD_opto or dataRB##_MMDD_baseline
% Gray/Blue = Control/Opto in opto sessions
% Green = Baseline sessions

% PARAMETERS
ILI_THRESHOLD = 0.3;
LICK_WINDOW_PRE = -3;
LICK_WINDOW_POST = 2;

% FIND ALL SESSIONS
allVars = who('dataRB12*');
allVars = sort(allVars);

nSessions = length(allVars);
fprintf('\n\n=== PLOTTING ALL SESSIONS ===\n');
fprintf('Found %d sessions (sorted):\n', nSessions);
for i = 1:length(allVars)
    fprintf('  %d. %s\n', i, allVars{i});
end

% COLLECT DATA
allLicks = [];
allLabels = [];

for s = 1:nSessions
    fprintf('\nProcessing %s...\n', allVars{s});
    
    try
        data = eval(allVars{s});
        
        % Determine session type from variable name
        isBaselineSession = contains(allVars{s}, 'baseline');
        
        % Extract data
        nTrials = data.i;
        StimTime = data.options.Trialtsecbegin;
        LickTimes = str2double(data.lickData(~cellfun('isempty', data.lickData)));
        redFrames = data.redFrames(1, :);
        RedFramesNum = str2double(redFrames) + ((StimTime - data.options.tsecbegin) * 1000)';
        optoTrials = data.options.opto_trial(1:nTrials);
        
        % Get licks relative to stimulus
        relLicks_control = [];
        relLicks_opto = [];
        
        for t = 1:min(nTrials, length(RedFramesNum))
            limMin = RedFramesNum(t) + LICK_WINDOW_PRE * 1000;
            limMax = RedFramesNum(t) + LICK_WINDOW_POST * 1000;
            
            licksInWindow = LickTimes(LickTimes > limMin & LickTimes < limMax);
            if ~isempty(licksInWindow)
                licksRel = (licksInWindow - RedFramesNum(t)) / 1000;
                
                % BOUT DETECTION - Only keep first lick of each bout
                licksRel = sort(licksRel);
                boutIdx = [true; diff(licksRel) > ILI_THRESHOLD];
                boutStartTimes = licksRel(boutIdx);
                
                % Separate by opto condition
                if optoTrials(t) == 0
                    relLicks_control = [relLicks_control; boutStartTimes];
                else
                    relLicks_opto = [relLicks_opto; boutStartTimes];
                end
            end
        end
        
        % Extract date (MMDD format -> display as MM/DD)
        dateParts = regexp(allVars{s}, 'dataRB\d+_(\d{2})(\d{2})', 'tokens');
        if ~isempty(dateParts)
            month = dateParts{1}{1};
            day = dateParts{1}{2};
            sessionLabel = sprintf('%s/%s', month, day);
        else
            sessionLabel = sprintf('S%d', s);
        end
        
        % Create labels based on session type
        if isBaselineSession
            % BASELINE SESSION → 1 GREEN VIOLIN (all trials)
            allLicks = [allLicks; relLicks_control];
            allLabels = [allLabels; repmat({[sessionLabel ' Base']}, length(relLicks_control), 1)];
            fprintf('  -> Baseline: %d violin (green)\n', length(relLicks_control));
        else
            % OPTO SESSION → 2 VIOLINS (gray for opto=0, blue for opto=1)
            if ~isempty(relLicks_control)
                allLicks = [allLicks; relLicks_control];
                allLabels = [allLabels; repmat({[sessionLabel ' C']}, length(relLicks_control), 1)];
            end
            if ~isempty(relLicks_opto)
                allLicks = [allLicks; relLicks_opto];
                allLabels = [allLabels; repmat({[sessionLabel ' O']}, length(relLicks_opto), 1)];
            end
            fprintf('  -> Opto: %d control (gray) + %d opto bout starts (blue)\n', ...
                length(relLicks_control), length(relLicks_opto));
        end
        
    catch ME
        fprintf('  ERROR: %s\n', ME.message);
    end
end

% PLOT
if ~isempty(allLicks)
    uniqueLabels = unique(allLabels, 'stable');
    nGroups = length(uniqueLabels);
    
    fprintf('\nPlotting %d violins:\n', nGroups);
    for i = 1:nGroups
        fprintf('  %d. %s\n', i, uniqueLabels{i});
    end
    
    figure('Color', 'w', 'Position', [50 50 max(900, 80*nGroups) 600]);
    hold on;
    
    % Reference lines
    plot([0.5 nGroups+0.5], [-2 -2], 'Color', [0.8 0.2 0.2], 'LineWidth', 2, 'LineStyle', '--');
    plot([0.5 nGroups+0.5], [0 0], 'Color', [0.2 0.2 0.2], 'LineWidth', 2, 'LineStyle', '-');
    
    % Assign colors based on label
    colors = zeros(nGroups, 3);
    for i = 1:nGroups
        if contains(uniqueLabels{i}, ' C')
            colors(i, :) = [0.7 0.7 0.7];  % Gray - Control trials
        elseif contains(uniqueLabels{i}, ' O')
            colors(i, :) = [0.4 0.6 1];     % Blue - Opto trials
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
    title('Lick Time Progression Across Sessions', 'FontSize', 16, 'FontWeight', 'bold');
    ylim([LICK_WINDOW_PRE LICK_WINDOW_POST]);
    set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'Box', 'off');
    xtickangle(45);
    grid on;
    set(gca, 'GridAlpha', 0.2);
    
    % Legend
    h_gray = patch(NaN, NaN, [0.7 0.7 0.7]);
    h_blue = patch(NaN, NaN, [0.4 0.6 1]);
    h_green = patch(NaN, NaN, [0.5 0.8 0.5]);
    legend([h_gray h_blue h_green], {'Control', 'Opto', 'Baseline session'}, ...
        'Location', 'northeast', 'FontSize', 10);
    
    % Reference line labels
    text(nGroups+0.8, -1, 'Opto', 'FontSize', 9, 'Color', [0.8 0.2 0.2]);
    text(nGroups+0.8, 0, 'Stim', 'FontSize', 9, 'Color', [0.2 0.2 0.2]);
    
    hold off;
end