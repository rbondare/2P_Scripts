%% this script is for random timing of optogenetics (binning it to early and late opto) 
% and for different contrasts (need to adjust it different contrasts) 

stimContrasts = data.options.contrastHistory;
optoTrials = data.options.opto_trial;

% Calculate opto timing relative to stimulus
TotalTrialWaitTime = zeros(max(size(data.waitTimes,1)), 1);
TotalTrialWaitTime = sum(data.waitTimes, 2);
TotalTrialWaitTime = -TotalTrialWaitTime;

opto_trials = data.options.opto_trial(:);
parallel_opto_start = data.options.parallel_opto_start(:);
TotalTrialWaitTime = TotalTrialWaitTime(:);

% Find indices of opto trials
opto_idx = find(opto_trials == 1);

% Preallocate output vector
opto_time_relative_to_stim = nan(size(TotalTrialWaitTime));

% Subtract opto start times from trial times (only for opto trials)
opto_time_relative_to_stim(opto_idx) = TotalTrialWaitTime(opto_idx) + parallel_opto_start;


%% ========== OPTO TIMING BINNING ==========
% Define time bins (relative to stimulus onset)
BIN_EARLY_START = -2;   % Early opto: -2 to 0 seconds (close to stimulus)
BIN_EARLY_END = 0;
BIN_LATE_START = -6;    % Late opto: -5 to -2 seconds (farther from stimulus)
BIN_LATE_END = 0;

% Create binary masks for each bin
opto_early = false(size(opto_trials));  % -2 to 0 sec
opto_late = false(size(opto_trials));   % -5 to -2 sec
opto_other = false(size(opto_trials));  % Outside defined bins

for t = 1:length(opto_trials)
    if opto_trials(t) == 1 && ~isnan(opto_time_relative_to_stim(t))
        opto_time = opto_time_relative_to_stim(t);
        
        if opto_time >= BIN_EARLY_START && opto_time < BIN_EARLY_END
            opto_early(t) = true;
        elseif opto_time >= BIN_LATE_START && opto_time < BIN_LATE_END
            opto_late(t) = true;
        else
            opto_other(t) = true;
        end
    end
end


%%
optoTrials_early = optoTrials;  % Copy original
for t = 1:length(optoTrials)
    if optoTrials(t) == 1 && ~opto_early(t)
        optoTrials_early(t) = 0;  % Turn non-early opto to control
    end
end

optoTrials_late = optoTrials;  % Copy original
for t = 1:length(optoTrials)
    if optoTrials(t) == 1 && ~opto_late(t)
        optoTrials_late(t) = 0;  % Turn non-late opto to control
    end
end

% Print summary
fprintf('Vector size: %d\n', length(optoTrials));
fprintf('Original - Control: %d, Opto: %d\n', sum(optoTrials == 0), sum(optoTrials == 1));
fprintf('Early filtered - Control: %d, Early opto: %d\n', sum(optoTrials_early == 0), sum(optoTrials_early == 1));
fprintf('Late filtered - Control: %d, Late opto: %d\n', sum(optoTrials_late == 0), sum(optoTrials_late == 1));

% plot hit&miss trials 
if any(optoTrials == 1)
    plotOutcomesByConditionAndContrast(isHit, isMiss, isUnattended, optoTrials_late, stimContrasts, isDiscarded);
else
    fprintf('Baseline session - no contrast analysis\n');
end


%%
function plotOutcomesByConditionAndContrast(isHit, isMiss, isUnattended, optoTrials, stimContrasts, isDiscarded)
% Plot stacked bars for Control/Opto at Low/High contrast, excluding discarded trials.
    if nargin < 6
        isDiscarded = false(size(optoTrials));
    end
    
    n = numel(optoTrials);
    if numel(isHit) ~= n || numel(isMiss) ~= n || numel(isUnattended) ~= n || numel(stimContrasts) ~= n || numel(isDiscarded) ~= n
        error('All inputs must be same length.');
    end
    
    % Create masks excluding discarded trials
    low_contrast_ctrl = (optoTrials == 0) & (stimContrasts == 12) & ~isDiscarded;
    low_contrast_opto = (optoTrials == 1) & (stimContrasts == 12) & ~isDiscarded;
    high_contrast_ctrl = (optoTrials == 0) & (stimContrasts == 60) & ~isDiscarded;
    high_contrast_opto = (optoTrials == 1) & (stimContrasts == 60) & ~isDiscarded;
    
    masks = {low_contrast_ctrl, low_contrast_opto, high_contrast_ctrl, high_contrast_opto};
    propData = zeros(4,3);
    
    for k = 1:4
        d = sum(masks{k});
        if d > 0
            propData(k,:) = [sum(isHit(masks{k})), sum(isMiss(masks{k})), sum(isUnattended(masks{k}))] / d;
        else
            propData(k,:) = [0 0 0];
        end
    end
    
    if any(isnan(propData(:))) || any(isinf(propData(:)))
        warning('NaN/Inf detected in propData — skipping plotOutcomesByConditionAndContrast.');
        return;
    end
    
    groups = {'Control Low', 'Opto Low', 'Control High', 'Opto High'};
    figure('Color', 'w');
    hold on;
    
    colors = [0.3 0.75 0.4; 1.0 0.6 0.2; 0.75 0.75 0.75];
    b = bar(propData, 'stacked', 'BarWidth', 0.6, 'EdgeColor', 'none');
    
    try
        for i = 1:min(numel(b), size(propData,2))
            if isgraphics(b(i))
                b(i).FaceColor = colors(i,:);
            end
        end
    catch
    end
    
    set(gca, 'XTick', 1:4, 'XTickLabel', groups, 'FontSize', 11, 'LineWidth', 1.2, 'Box', 'off');
    ylabel('Proportion of Trials', 'FontSize', 13);
    ylim([0 1]);
    yticks(0:0.2:1);
    grid on;
    set(gca, 'GridAlpha', 0.2);
    legend({'Hit', 'Miss', 'Unattended'}, 'Location', 'northeastoutside', 'FontSize', 11, 'Box', 'off');
    title('Trial Outcomes by Condition and Contrast', 'FontSize', 14, 'FontWeight', 'bold');
    hold off;
end
