%% LICK PROBABILITY PSTH - MULTIPLE SESSIONS
% Creates PSTH-style plot of lick probability over time
% Aligned to stimulus onset
% Plots all selected sessions as separate lines

% PARAMETERS
ANIMAL_ID = 'RB15';  % For title only

% INCLUSION LIST
INCLUDED_SESSIONS = {
    'AnimalRB15_1030_1558_opto';
    'AnimalRB15_1103_1413_opto';
    'AnimalRB15_1105_1307_opto';
    'AnimalRB15_1106_1329_opto';
};

% Time window
TIME_WINDOW_PRE = -3;   % seconds before stimulus
TIME_WINDOW_POST = 3;   % seconds after stimulus

% Binning parameters
BIN_SIZE = 0.1;  % 100ms bins
SLIDE_STEP = 0.05;  % 50ms slide (50% overlap)

% CHECK DATA
if ~exist('recordings', 'var')
    error('recordings struct not found in workspace');
end

nSessions = length(INCLUDED_SESSIONS);
fprintf('\n=== PROCESSING %d SESSIONS ===\n', nSessions);

% CREATE TIME BINS
binEdges = TIME_WINDOW_PRE: SLIDE_STEP:(TIME_WINDOW_POST - BIN_SIZE);
nBins = length(binEdges);
binCenters = binEdges + BIN_SIZE/2;

% STORAGE FOR ALL SESSIONS
allData = struct();

% LOOP THROUGH SESSIONS
for s = 1:nSessions
    sessionName = INCLUDED_SESSIONS{s};
    fprintf('\n[%d/%d] Processing:  %s\n', s, nSessions, sessionName);
    
    if ~isfield(recordings, sessionName)
        warning('Session not found, skipping:  %s', sessionName);
        continue;
    end
    
    data = recordings.(sessionName);
    
    % EXTRACT DATA
    nTrials = data.i;
    StimTime = data. options. Trialtsecbegin;
    LickTimes = str2double(data.lickData(~cellfun('isempty', data.lickData)));
    redFrames = data. redFrames(1, :);
    RedFramesNum = str2double(redFrames) + ((StimTime - data.options.tsecbegin) * 1000)';
    optoTrials = data.options.opto_trial(1:nTrials);
    
    isBaselineSession = contains(sessionName, 'baseline');
    
    fprintf('  Total trials: %d | Control:  %d | Opto: %d\n', ... 
        nTrials, sum(optoTrials == 0), sum(optoTrials == 1));
    
    % BUILD LICK MATRIX
    maxLicks = 100;
    LicksInFrame = zeros(maxLicks, nTrials);
    
    for t = 1:min(nTrials, length(RedFramesNum))
        limMin = RedFramesNum(t) + TIME_WINDOW_PRE * 1000;
        limMax = RedFramesNum(t) + TIME_WINDOW_POST * 1000;
        
        licksInWindow = LickTimes(LickTimes > limMin & LickTimes < limMax);
        
        if ~isempty(licksInWindow)
            licksRel = (licksInWindow - RedFramesNum(t)) / 1000;
            LicksInFrame(1: length(licksRel), t) = licksRel;
        end
    end
    

    
    
    % SEPARATE CONTROL VS OPTO
    idx_ctrl = (optoTrials == 0);
    idx_opto = (optoTrials == 1);
    
    LicksCtrl = LicksInFrame(: , idx_ctrl);
    LicksOpto = LicksInFrame(:, idx_opto);
    
    nCtrl = sum(idx_ctrl);
    nOpto = sum(idx_opto);
    
    fprintf('  Valid - Control: %d | Opto: %d\n', nCtrl, nOpto);
    
    % CALCULATE LICK PROBABILITY
    probCtrl = zeros(1, nBins);
    probOpto = zeros(1, nBins);
    
    for b = 1:nBins
        binStart = binEdges(b);
        binEnd = binEdges(b) + BIN_SIZE;
        
        % Control
        if nCtrl > 0
            hasLickCtrl = false(1, nCtrl);
            for t = 1:nCtrl
                licks = LicksCtrl(LicksCtrl(:,t) ~= 0, t);
                if any(licks >= binStart & licks < binEnd)
                    hasLickCtrl(t) = true;
                end
            end
            probCtrl(b) = sum(hasLickCtrl) / nCtrl;
        end
        
        % Opto
        if nOpto > 0
            hasLickOpto = false(1, nOpto);
            for t = 1:nOpto
                licks = LicksOpto(LicksOpto(:,t) ~= 0, t);
                if any(licks >= binStart & licks < binEnd)
                    hasLickOpto(t) = true;
                end
            end
            probOpto(b) = sum(hasLickOpto) / nOpto;
        end
    end
    
    % STORE DATA
    % Extract date label
    dateParts = regexp(sessionName, 'RB\d*_(\d{2})(\d{2})_', 'tokens');
    if ~isempty(dateParts)
        sessionLabel = sprintf('%s/%s', dateParts{1}{1}, dateParts{1}{2});
    else
        sessionLabel = sprintf('S%d', s);
    end
    
    allData(s).name = sessionName;
    allData(s).label = sessionLabel;
    allData(s).isBaseline = isBaselineSession;
    allData(s).probCtrl = probCtrl;
    allData(s).probOpto = probOpto;
    allData(s).nCtrl = nCtrl;
    allData(s).nOpto = nOpto;
end

% Colorblind-friendly palette
colors_ctrl = [
    0.00, 0.45, 0.74;   % Blue
    0.85, 0.33, 0.10;   % Red-orange
    0.93, 0.69, 0.13;   % Yellow
    0.49, 0.18, 0.56;   % Purple
    0.47, 0.67, 0.19;   % Green
    0.30, 0.75, 0.93;   % Cyan
    0.64, 0.08, 0.18;   % Dark red
];

colors_opto = colors_ctrl;  % Use same palette for consistency
% FIGURE 1: SIDE-BY-SIDE SUBPLOTS 
fig1 = figure('Color', 'w', 'Position', [100 100 1000 450]);
set(fig1, 'PaperPositionMode', 'auto');
set(fig1, 'Renderer', 'painters');  % Vector graphics

% Set default font
set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');

% SUBPLOT A: CONTROL TRIALS
subplot(1, 2, 1);
hold on;

maxProb_ctrl = 0;
validSessions_ctrl = 0;

for s = 1:nSessions
    if isempty(allData(s).name) || allData(s).nCtrl == 0
        continue;
    end
    validSessions_ctrl = validSessions_ctrl + 1;
    
    colorIdx = mod(s-1, size(colors_ctrl, 1)) + 1;
    
    plot(binCenters, allData(s).probCtrl * 100, ...
        'Color', colors_ctrl(colorIdx,:), ... 
        'LineWidth', 2.5, ... 
        'DisplayName', allData(s).label);
    
    maxProb_ctrl = max(maxProb_ctrl, max(allData(s).probCtrl * 100));
end

% Stimulus onset line
ymax = ceil(maxProb_ctrl / 10) * 10 + 10;
plot([0 0], [0 ymax], 'k--', 'LineWidth', 1.5);

% Formatting
xlabel('Time from stimulus (s)', 'FontSize', 14, 'FontWeight', 'normal');
ylabel('Lick probability (%)', 'FontSize', 14, 'FontWeight', 'normal');
title('Control', 'FontSize', 16, 'FontWeight', 'bold');

xlim([TIME_WINDOW_PRE TIME_WINDOW_POST]);
ylim([0 ymax]);

% Clean legend
leg = legend({'10/30', '11/03','11/05','11/06'},...
    'Location', 'eastoutside', 'FontSize', 10, 'Box', 'off');
leg.ItemTokenSize = [20, 10];

% Axes
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
ax.TickDir = 'out';
ax.TickLength = [0.015 0.015];
ax.Box = 'off';
grid on;
ax.GridAlpha = 0.15;
ax.GridLineStyle = '-';

hold off;

% SUBPLOT B: OPTO TRIALS
subplot(1, 2, 2);
hold on;

maxProb_opto = 0;
validSessions_opto = 0;

for s = 1:nSessions
    if isempty(allData(s).name) || allData(s).isBaseline || allData(s).nOpto == 0
        continue;
    end
    validSessions_opto = validSessions_opto + 1;
    
    colorIdx = mod(s-1, size(colors_opto, 1)) + 1;
    
    plot(binCenters, allData(s).probOpto * 100, ... 
        'Color', colors_opto(colorIdx,:), ...
        'LineWidth', 2.5, ...
        'DisplayName', allData(s).label);
    
    maxProb_opto = max(maxProb_opto, max(allData(s).probOpto * 100));
end

% Stimulus onset line
ymax = ceil(maxProb_opto / 10) * 10 + 10;
plot([0 0], [0 ymax], 'k--', 'LineWidth', 1.5);

% Formatting
xlabel('Time from stimulus (s)', 'FontSize', 14, 'FontWeight', 'normal');
ylabel('Lick probability (%)', 'FontSize', 14, 'FontWeight', 'normal');
title('Opto Trials', 'FontSize', 16, 'FontWeight', 'bold');

xlim([TIME_WINDOW_PRE TIME_WINDOW_POST]);
ylim([0 ymax]);

% Legend
leg = legend({'10/30', '11/03','11/05','11/06'},...
    'Location', 'eastoutside', 'FontSize', 10, 'Box', 'off');
leg.ItemTokenSize = [20, 10];

% Axes
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.5;
ax.TickDir = 'out';
ax.TickLength = [0.015 0.015];
ax.Box = 'off';
grid on;
ax.GridAlpha = 0.15;
ax.GridLineStyle = '-';

hold off;

% Overall title
sgtitle(sprintf('%s - Licking across sessions', ANIMAL_ID), ...
    'FontSize', 18, 'FontWeight', 'bold');

% FIGURE 2: OVERLAY PLOT (SUPPLEMENTARY)
fig2 = figure('Color', 'w', 'Position', [150 150 800 600]);
set(fig2, 'PaperPositionMode', 'auto');
set(fig2, 'Renderer', 'painters');

hold on;

maxProb_all = 0;

% Plot all control (solid lines)
for s = 1:nSessions
    if isempty(allData(s).name) || allData(s).nCtrl == 0
        continue;
    end
    
    colorIdx = mod(s-1, size(colors_ctrl, 1)) + 1;
    
    plot(binCenters, allData(s).probCtrl * 100, ...
        'Color', colors_ctrl(colorIdx,:), ...
        'LineWidth', 2, ...
        'LineStyle', '-', ...
        'DisplayName', sprintf('%s Control', allData(s).label));
    
    maxProb_all = max(maxProb_all, max(allData(s).probCtrl * 100));
end

% Plot all opto (dashed lines)
for s = 1:nSessions
    if isempty(allData(s).name) || allData(s).isBaseline || allData(s).nOpto == 0
        continue;
    end
    
    colorIdx = mod(s-1, size(colors_opto, 1)) + 1;
    
    plot(binCenters, allData(s).probOpto * 100, ... 
        'Color', colors_opto(colorIdx,:), ...
        'LineWidth', 2, ...
        'LineStyle', '--', ...
        'DisplayName', sprintf('%s Opto', allData(s).label));
    
    maxProb_all = max(maxProb_all, max(allData(s).probOpto * 100));
end

% Stimulus onset
ymax = ceil(maxProb_all / 10) * 10 + 10;
plot([0 0], [0 ymax], 'k-', 'LineWidth', 2);

% Formatting
xlabel('Time from stimulus (s)', 'FontSize', 16, 'FontWeight', 'normal');
ylabel('Lick probability (%)', 'FontSize', 16, 'FontWeight', 'normal');
title(sprintf('%s - All sessions (solid:  control, dashed: opto)', ANIMAL_ID), ...
    'FontSize', 18, 'FontWeight', 'bold');

xlim([TIME_WINDOW_PRE TIME_WINDOW_POST]);
ylim([0 ymax]);

% Legend
leg = legend({'10/30', '11/03','11/05','11/06'},...
    'Location', 'eastoutside', 'FontSize', 10, 'Box', 'off');
leg.ItemTokenSize = [30, 10];



% Axes
ax = gca;
ax.FontSize = 14;
ax.LineWidth = 1.5;
ax.TickDir = 'out';
ax. TickLength = [0.015 0.015];
ax. Box = 'off';
grid on;
ax.GridAlpha = 0.15;
ax.GridLineStyle = '-';

hold off;
