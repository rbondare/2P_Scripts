%% Lick correlations - Entropy Analysis
function lickCorrelations_Entropy(Entropy, WindowLength, params, Rec, Performance, conditions, Skip, positions)
% LICKCORRELATIONS_ENTROPY - Analyse pre- and post-lick entropy and produce
% correlation and distribution plots across different behavioral conditions and stimulus positions.
%
% INPUTS:
%   Entropy      - Structure containing entropy data for different conditions
%   WindowLength - Structure containing window length information
%   params       - Structure containing plotting parameters including colors
%   Rec          - Recording metadata structure
%   Performance  - Performance metrics structure containing lick data
%   conditions   - Cell array of condition names to analyze (e.g., {'Expert', 'Beginner'})
%   Skip         - (optional) Structure containing recording indices to skip for each condition
%   positions    - (optional) Cell array of positions to analyze ({'All', 'P1', 'P3'})
%
% OUTPUTS:
%   • Figure 1 : Pre vs Post Entropy scatter (3 lick types: Extraneous | Stimulus | Reward) (Raw top | Z-scored bottom)
%   • Figure 2 : Pre vs Post Entropy violin plot (3 lick types: Extraneous | Stimulus | Reward) (Raw left | Z-scored right)
%   • Figure 3 : Histogram of Post/Pre entropy ratios (3 lick types) (Raw top | Z-scored bottom) – sign‑rank against median of 1
%   • Figure 4 : Pre value vs Change scatter (3 lick types) (Raw top | Z-scored bottom) – Spearman ρ with p‑value
%   • Figure 5 : Extraneous Licks - Lick-Triggered Average (Raw | Z-scored) (time-locked to extraneous licks, smart NaN padding)
%   • Figure 6 : Stimulus Licks - Lick-Triggered Average (Raw | Z-scored) (time-locked to first licks, sorted by reaction time)*
%   • Figure 7 : Stimulus Licks - Lick-Triggered Modulation (Raw | Z-scored) (baseline-corrected by miss trials)*
%   • Figure 8 : Stimulus Licks - Reward-Triggered Average (Raw | Z-scored) (time-locked to reward delivery @ 350ms after first lick)*
%   • Figure 9 : Stimulus Licks - Reward Modulation (Raw | Z-scored) (pre-reward @ 300ms, post-reward @ 600ms after first lick)*
%   • Figure 10: Stimulus Licks - Previous trial history analysis (Raw | Z-scored) (Hit vs Miss trial effects)*

%   
%   Note: ALL figures require lick timing data and spout presence. This function automatically
%         skips 'Naive' and 'NoSpout' conditions as they have no spout present during recording.
%         Figure 8 specifically analyzes reward delivery at 350ms after first lick.
%
% EXAMPLE:
%   lickCorrelations_Entropy(Entropy, WindowLength, params, Rec, Performance);
%   lickCorrelations_Entropy(Entropy, WindowLength, params, Rec, Performance, {'Expert'});
%   lickCorrelations_Entropy(Entropy, WindowLength, params, Rec, Performance, {'Expert', 'Beginner'}, Skip, {'All', 'P1'});

%% Check inputs and set defaults
if nargin < 8 || isempty(positions)
    positions = {'All', 'P1', 'P3'};
end

if nargin < 7 || isempty(Skip)
    Skip = struct(); 
end

if nargin < 6 || isempty(conditions)
    % Default to all available conditions if not specified
    conditions = {};
    if isfield(Entropy, 'Expert'); conditions{end+1} = 'Expert'; end
    if isfield(Entropy, 'Beginner'); conditions{end+1} = 'Beginner'; end
end

fprintf('Processing lick correlations for ENTROPY\n');
fprintf('Conditions: %s\n', strjoin(conditions, ', '));
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

%% Initialize summary data collection for cross-position plot
summary_ExtraneousLTA = [];
summary_P3_StimulusLTA = [];
summary_P3_RewardRTA = [];
summary_timeAxis_LTA = [];
summary_timeAxis_StimulusLTA = [];
summary_timeAxis_RTA = [];

%% Loop through each position with "Data First, Plots Second" approach
for pos_idx = 1:length(positions)
    position = positions{pos_idx};
    fprintf('\n\n======== PROCESSING POSITION: %s ========\n', position);
    
    %% ==================== PHASE 1: DATA EXTRACTION ====================
    fprintf('\n*** PHASE 1: DATA EXTRACTION ***\n');
    
    % Time axis (–8 s to +10.5 s around lick)
    timeReference = linspace(-8, 10.5, WindowLength.All);
    ifi = 18.5/WindowLength.All; %#ok<NASGU>
    
    % Analysis window parameters
    LTAwindow = 40; % Frames around lick event for stimulus analysis
    RTAwindow = 30; % Frames around reward event
    
    % Check if any conditions have lick data (exclude Naive/NoSpout which have no spout)
    conditionsWithLicks = {};
    for c_idx = 1:length(conditions)
        condName = conditions{c_idx};
        if ~strcmp(condName, 'Naive') && ~strcmp(condName, 'NoSpout')
            conditionsWithLicks{end+1} = condName; %#ok<AGROW>
        end
    end
    
    if isempty(conditionsWithLicks)
        fprintf('Skipping analyses for %s - no conditions with spout data available\n', position);
        fprintf('(Naive and NoSpout conditions have no spout present, so licks cannot be recorded)\n');
        continue;
    else
        fprintf('Performing analyses for %s using conditions with spout: %s\n', position, strjoin(conditionsWithLicks, ', '));
        fprintf('(Excluding Naive/NoSpout: no spout present during recording)\n');
    end
    
    % Step 1: Extract unified valid lick trials
    fprintf('\n--- Step 1: Unified Trial Selection ---\n');
    validLickTrials = extractValidLickTrials_Entropy(conditionsWithLicks, position, Entropy, Performance, Rec, Skip, timeReference, LTAwindow, RTAwindow);
    
    % Step 2: Extract data for Figures 1-4 (Pre/Post entropy)
    fprintf('\n--- Step 2: Pre/Post Entropy Data ---\n');
    Fig1_Data = extractPrePostEntropy(validLickTrials, Entropy, Performance, timeReference);
    Fig3_Data = calculateEntropyRatios(Fig1_Data);
    Fig4_Data = calculateEntropyChanges(Fig1_Data);
    
    % Step 3: Extract data for Figure 5 (Extraneous LTA)  
    fprintf('\n--- Step 3: Extraneous LTA Data ---\n');
    Fig5_Data = extractExtraneousLTA_Entropy(conditionsWithLicks, Entropy, Performance, timeReference, 30);
    
    % Step 4: Extract data for Figure 6 (Stimulus LTA)
    fprintf('\n--- Step 4: Stimulus LTA Data ---\n');
    Fig6_Data = extractStimulusLTA_Entropy(validLickTrials, timeReference, LTAwindow);
    
    % Step 5: Extract data for Figure 7 (Stimulus LTM)
    fprintf('\n--- Step 5: Stimulus LTM Data ---\n');
    Fig7_Data = extractStimulusLTM_Entropy(validLickTrials, Entropy, Performance, timeReference, LTAwindow, position);
    
    % Step 6: Extract data for Figure 8 (Reward RTA)
    fprintf('\n--- Step 6: Reward RTA Data ---\n');
    Fig8_Data = extractRewardRTA_Entropy(validLickTrials, timeReference, RTAwindow);
    
    % Step 7: Extract data for Figure 9 (Reward Modulation)
    fprintf('\n--- Step 7: Reward Modulation Data ---\n');
    Fig9_Data = extractRewardModulation_Entropy(validLickTrials, timeReference);
    
    % Step 8: Extract data for Figure 10 (Previous Trial History)
    fprintf('\n--- Step 8: Previous Trial History Data ---\n');
    Fig10_Data = extractPreviousTrialHistory_Entropy(conditionsWithLicks, Entropy, Performance, Rec, Skip, position);
    
    % Step 9: Validate trial count consistency
    fprintf('\n--- Step 9: Validation ---\n');
    validateTrialCounts_Entropy(validLickTrials, Fig1_Data, Fig6_Data, Fig8_Data, Fig9_Data, position);
    
    %% ==================== PHASE 2: PLOT GENERATION ====================
    fprintf('\n*** PHASE 2: PLOT GENERATION ***\n');
    
    % Now all data has been extracted - proceed with plotting using pre-extracted data
    
    %% Check if we have any data for this position
    if isempty(Fig1_Data.Extraneous.Pre)
        fprintf('No extraneous lick entropy data found for position %s\n', position);
        fprintf('Skipping all plots for this position\n');
        continue;
    end
    
    fprintf('Data extraction complete for %s:\n', position);
    fprintf('  Extraneous licks: %d\n', size(Fig1_Data.Extraneous.Pre, 1));
    fprintf('  Stimulus licks: %d\n', size(Fig1_Data.Stimulus.Pre, 1));
    fprintf('  Reward timepoints: %d\n', size(Fig1_Data.Reward.Pre, 1));


    
    %% ------------------------- Figure 1: All Lick Types Entropy scatter -------------------- %%
    fprintf('\n--- Generating Figure 1: Pre vs Post Entropy Scatter ---\n');
    
    % Create the expanded figure with 2x3 layout
    figure('Color', 'w', 'Name', sprintf('All Lick Types: Pre vs Post Entropy Scatter - %s', position)); 
    tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Row 1: Raw data
    % Extraneous licks (raw)
    ax1 = nexttile;
    dscatter(Fig1_Data.Extraneous.Pre, Fig1_Data.Extraneous.Post, 'MSIZE', 30);
    colormap(ax1, gray);
    axis(ax1, 'square');
    % Set equal axis limits
    allData1 = [Fig1_Data.Extraneous.Pre; Fig1_Data.Extraneous.Post];
    dataLim1 = [min(allData1), max(allData1)];
    xlim(ax1, dataLim1);
    ylim(ax1, dataLim1);
    % Set nice axis ticks at half and full numbers
    setNiceAxisTicks(ax1, dataLim1);
    ref1 = refline(ax1, 1); 
    ref1.LineWidth = 3; 
    ref1.Color = 'r';
    % Add fitted trend line
    fit1 = lsline(ax1);
    fit1.LineStyle = '--';
    fit1.LineWidth = 2;
    fit1.Color = [0.5 0.5 0.5];
    
    [pEnt, ~, ~] = signrank(Fig1_Data.Extraneous.Pre, Fig1_Data.Extraneous.Post);
    [rhoEnt, pRhoEnt] = corr(Fig1_Data.Extraneous.Pre, Fig1_Data.Extraneous.Post, 'Type', 'Spearman', 'Rows', 'complete');
    
    % Calculate Hedges' g effect size
    mean_pre = mean(Fig1_Data.Extraneous.Pre);
    mean_post = mean(Fig1_Data.Extraneous.Post);
    std_pre = std(Fig1_Data.Extraneous.Pre);
    std_post = std(Fig1_Data.Extraneous.Post);
    n_pre = length(Fig1_Data.Extraneous.Pre);
    n_post = length(Fig1_Data.Extraneous.Post);
    pooled_std = sqrt(((n_pre - 1) * std_pre^2 + (n_post - 1) * std_post^2) / (n_pre + n_post - 2));
    cohens_d = (mean_post - mean_pre) / pooled_std;
    correction_factor = 1 - (3 / (4 * (n_pre + n_post - 2) - 1));
    hedges_g = cohens_d * correction_factor;
    if abs(hedges_g) < 0.2; effect_interp = 'negligible'; elseif abs(hedges_g) < 0.5; effect_interp = 'small'; elseif abs(hedges_g) < 0.8; effect_interp = 'medium'; else; effect_interp = 'large'; end
    
    ax1.FontSize = 16;
    xlabel(ax1, 'Pre-Extraneous-Lick Entropy');
    ylabel(ax1, 'Post-Extraneous-Lick Entropy');
    title(ax1, {'Extraneous Licks (Raw)', sprintf('n=%d, p_{signrank}=%.3g, rho=%.3f, p_{rho}=%.3g', length(Fig1_Data.Extraneous.Pre), pEnt, rhoEnt, pRhoEnt), sprintf('Hedges'' g=%.3f (%s)', hedges_g, effect_interp)});
    
    % Stimulus licks (raw)
    ax2 = nexttile;
    if ~isempty(Fig1_Data.Stimulus.Pre)
        dscatter(Fig1_Data.Stimulus.Pre, Fig1_Data.Stimulus.Post, 'MSIZE', 30);
        colormap(ax2, gray);
        axis(ax2, 'square');
        % Set equal axis limits
        allData2 = [Fig1_Data.Stimulus.Pre; Fig1_Data.Stimulus.Post];
        dataLim2 = [min(allData2), max(allData2)];
        xlim(ax2, dataLim2);
        % Set nice axis ticks at half and full numbers
        setNiceAxisTicks(ax2, dataLim2);
        ylim(ax2, dataLim2);
        ref2 = refline(ax2, 1); 
        ref2.LineWidth = 3; 
        ref2.Color = 'r';
        % Add fitted trend line
        fit2 = lsline(ax2);
        fit2.LineStyle = '--';
        fit2.LineWidth = 2;
        fit2.Color = [0.5 0.5 0.5];
        
        [pStim, ~, ~] = signrank(Fig1_Data.Stimulus.Pre, Fig1_Data.Stimulus.Post);
        [rhoStim, pRhoStim] = corr(Fig1_Data.Stimulus.Pre, Fig1_Data.Stimulus.Post, 'Type', 'Spearman', 'Rows', 'complete');
        
        % Calculate Hedges' g effect size
        mean_pre = mean(Fig1_Data.Stimulus.Pre);
        mean_post = mean(Fig1_Data.Stimulus.Post);
        std_pre = std(Fig1_Data.Stimulus.Pre);
        std_post = std(Fig1_Data.Stimulus.Post);
        n_pre = length(Fig1_Data.Stimulus.Pre);
        n_post = length(Fig1_Data.Stimulus.Post);
        pooled_std = sqrt(((n_pre - 1) * std_pre^2 + (n_post - 1) * std_post^2) / (n_pre + n_post - 2));
        cohens_d = (mean_post - mean_pre) / pooled_std;
        correction_factor = 1 - (3 / (4 * (n_pre + n_post - 2) - 1));
        hedges_g = cohens_d * correction_factor;
        if abs(hedges_g) < 0.2; effect_interp = 'negligible'; elseif abs(hedges_g) < 0.5; effect_interp = 'small'; elseif abs(hedges_g) < 0.8; effect_interp = 'medium'; else; effect_interp = 'large'; end
        
        ax2.FontSize = 16;
        xlabel(ax2, 'Pre-Stimulus-Lick Entropy');
        ylabel(ax2, 'Post-Stimulus-Lick Entropy');
        title(ax2, {'Stimulus Licks (Raw)', sprintf('n=%d, p_{signrank}=%.3g, rho=%.3f, p_{rho}=%.3g', length(Fig1_Data.Stimulus.Pre), pStim, rhoStim, pRhoStim), sprintf('Hedges'' g=%.3f (%s)', hedges_g, effect_interp)});
    else
        text(ax2, 0.5, 0.5, 'No stimulus lick data', 'HorizontalAlignment', 'center');
        title(ax2, 'Stimulus Licks (Raw)');
        ax2.FontSize = 16;
    end
    
    % Reward timepoints (raw)
    ax3 = nexttile;
    if ~isempty(Fig1_Data.Reward.Pre)
        dscatter(Fig1_Data.Reward.Pre, Fig1_Data.Reward.Post, 'MSIZE', 30);
        colormap(ax3, gray);
        axis(ax3, 'square');
        % Set equal axis limits
        allData3 = [Fig1_Data.Reward.Pre; Fig1_Data.Reward.Post];
        dataLim3 = [min(allData3), max(allData3)];
        xlim(ax3, dataLim3);
        % Set nice axis ticks at half and full numbers
        setNiceAxisTicks(ax3, dataLim3);
        ylim(ax3, dataLim3);
        ref3 = refline(ax3, 1); 
        ref3.LineWidth = 3; 
        ref3.Color = 'r';
        % Add fitted trend line
        fit3 = lsline(ax3);
        fit3.LineStyle = '--';
        fit3.LineWidth = 2;
        fit3.Color = [0.5 0.5 0.5];
        
        [pRewardTP, ~, ~] = signrank(Fig1_Data.Reward.Pre, Fig1_Data.Reward.Post);
        [rhoRewardTP, pRhoRewardTP] = corr(Fig1_Data.Reward.Pre, Fig1_Data.Reward.Post, 'Type', 'Spearman', 'Rows', 'complete');
        
        % Calculate Hedges' g effect size
        mean_pre = mean(Fig1_Data.Reward.Pre);
        mean_post = mean(Fig1_Data.Reward.Post);
        std_pre = std(Fig1_Data.Reward.Pre);
        std_post = std(Fig1_Data.Reward.Post);
        n_pre = length(Fig1_Data.Reward.Pre);
        n_post = length(Fig1_Data.Reward.Post);
        pooled_std = sqrt(((n_pre - 1) * std_pre^2 + (n_post - 1) * std_post^2) / (n_pre + n_post - 2));
        cohens_d = (mean_post - mean_pre) / pooled_std;
        correction_factor = 1 - (3 / (4 * (n_pre + n_post - 2) - 1));
        hedges_g = cohens_d * correction_factor;
        if abs(hedges_g) < 0.2; effect_interp = 'negligible'; elseif abs(hedges_g) < 0.5; effect_interp = 'small'; elseif abs(hedges_g) < 0.8; effect_interp = 'medium'; else; effect_interp = 'large'; end
        
        ax3.FontSize = 16;
        xlabel(ax3, 'Pre-Reward Entropy');
        ylabel(ax3, 'Post-Reward Entropy');
        title(ax3, {'Reward Timepoints (Raw)', sprintf('n=%d, p_{signrank}=%.3g, rho=%.3f, p_{rho}=%.3g', length(Fig1_Data.Reward.Pre), pRewardTP, rhoRewardTP, pRhoRewardTP), sprintf('Hedges'' g=%.3f (%s)', hedges_g, effect_interp)});
    else
        text(ax3, 0.5, 0.5, 'No reward timepoint data', 'HorizontalAlignment', 'center');
        title(ax3, 'Reward Timepoints (Raw)');
        ax3.FontSize = 16;
    end
    
    % Row 2: Z-scored data
    % Extraneous licks (z-scored)
    ax4 = nexttile;
    dscatter(Fig1_Data.Extraneous.Pre_zscore, Fig1_Data.Extraneous.Post_zscore, 'MSIZE', 30);
    colormap(ax4, gray);
    axis(ax4, 'square');
    % Set equal axis limits
    allData4 = [Fig1_Data.Extraneous.Pre_zscore; Fig1_Data.Extraneous.Post_zscore];
    dataLim4 = [min(allData4), max(allData4)];
    xlim(ax4, dataLim4);
    % Set nice axis ticks at half and full numbers
    setNiceAxisTicks(ax4, dataLim4);
    ylim(ax4, dataLim4);
    ref4 = refline(ax4, 1); 
    ref4.LineWidth = 3; 
    ref4.Color = 'r';
    % Add fitted trend line
    fit4 = lsline(ax4);
    fit4.LineStyle = '--';
    fit4.LineWidth = 2;
    fit4.Color = [0.5 0.5 0.5];
    
    [pEnt_zscore, ~, ~] = signrank(Fig1_Data.Extraneous.Pre_zscore, Fig1_Data.Extraneous.Post_zscore);
    [rhoEnt_zscore, pRhoEnt_zscore] = corr(Fig1_Data.Extraneous.Pre_zscore, Fig1_Data.Extraneous.Post_zscore, 'Type', 'Spearman', 'Rows', 'complete');
    
    % Calculate Hedges' g effect size
    mean_pre = mean(Fig1_Data.Extraneous.Pre_zscore);
    mean_post = mean(Fig1_Data.Extraneous.Post_zscore);
    std_pre = std(Fig1_Data.Extraneous.Pre_zscore);
    std_post = std(Fig1_Data.Extraneous.Post_zscore);
    n_pre = length(Fig1_Data.Extraneous.Pre_zscore);
    n_post = length(Fig1_Data.Extraneous.Post_zscore);
    pooled_std = sqrt(((n_pre - 1) * std_pre^2 + (n_post - 1) * std_post^2) / (n_pre + n_post - 2));
    cohens_d = (mean_post - mean_pre) / pooled_std;
    correction_factor = 1 - (3 / (4 * (n_pre + n_post - 2) - 1));
    hedges_g = cohens_d * correction_factor;
    if abs(hedges_g) < 0.2; effect_interp = 'negligible'; elseif abs(hedges_g) < 0.5; effect_interp = 'small'; elseif abs(hedges_g) < 0.8; effect_interp = 'medium'; else; effect_interp = 'large'; end
    
    ax4.FontSize = 16;
    xlabel(ax4, 'Pre-Extraneous-Lick Entropy (Z)');
    ylabel(ax4, 'Post-Extraneous-Lick Entropy (Z)');
    title(ax4, {'Extraneous Licks (Z-scored)', sprintf('n=%d, p_{signrank}=%.3g, rho=%.3f, p_{rho}=%.3g', length(Fig1_Data.Extraneous.Pre_zscore), pEnt_zscore, rhoEnt_zscore, pRhoEnt_zscore), sprintf('Hedges'' g=%.3f (%s)', hedges_g, effect_interp)});
    
    % Stimulus licks (z-scored)
    ax5 = nexttile;
    if ~isempty(Fig1_Data.Stimulus.Pre) && isfield(Fig1_Data.Stimulus, 'Pre_zscore')
        dscatter(Fig1_Data.Stimulus.Pre_zscore, Fig1_Data.Stimulus.Post_zscore, 'MSIZE', 30);
        colormap(ax5, gray);
        axis(ax5, 'square');
        % Set equal axis limits
        allData5 = [Fig1_Data.Stimulus.Pre_zscore; Fig1_Data.Stimulus.Post_zscore];
        dataLim5 = [min(allData5), max(allData5)];
        xlim(ax5, dataLim5);
        % Set nice axis ticks at half and full numbers
        setNiceAxisTicks(ax5, dataLim5);
        ylim(ax5, dataLim5);
        ref5 = refline(ax5, 1); 
        ref5.LineWidth = 3; 
        ref5.Color = 'r';
        % Add fitted trend line
        fit5 = lsline(ax5);
        fit5.LineStyle = '--';
        fit5.LineWidth = 2;
        fit5.Color = [0.5 0.5 0.5];
        
        [pStim_zscore, ~, ~] = signrank(Fig1_Data.Stimulus.Pre_zscore, Fig1_Data.Stimulus.Post_zscore);
        [rhoStim_zscore, pRhoStim_zscore] = corr(Fig1_Data.Stimulus.Pre_zscore, Fig1_Data.Stimulus.Post_zscore, 'Type', 'Spearman', 'Rows', 'complete');
        
        % Calculate Hedges' g effect size
        mean_pre = mean(Fig1_Data.Stimulus.Pre_zscore);
        mean_post = mean(Fig1_Data.Stimulus.Post_zscore);
        std_pre = std(Fig1_Data.Stimulus.Pre_zscore);
        std_post = std(Fig1_Data.Stimulus.Post_zscore);
        n_pre = length(Fig1_Data.Stimulus.Pre_zscore);
        n_post = length(Fig1_Data.Stimulus.Post_zscore);
        pooled_std = sqrt(((n_pre - 1) * std_pre^2 + (n_post - 1) * std_post^2) / (n_pre + n_post - 2));
        cohens_d = (mean_post - mean_pre) / pooled_std;
        correction_factor = 1 - (3 / (4 * (n_pre + n_post - 2) - 1));
        hedges_g = cohens_d * correction_factor;
        if abs(hedges_g) < 0.2; effect_interp = 'negligible'; elseif abs(hedges_g) < 0.5; effect_interp = 'small'; elseif abs(hedges_g) < 0.8; effect_interp = 'medium'; else; effect_interp = 'large'; end
        
        ax5.FontSize = 16;
        xlabel(ax5, 'Pre-Stimulus-Lick Entropy (Z)');
        ylabel(ax5, 'Post-Stimulus-Lick Entropy (Z)');
        title(ax5, {'Stimulus Licks (Z-scored)', sprintf('n=%d, p_{signrank}=%.3g, rho=%.3f, p_{rho}=%.3g', length(Fig1_Data.Stimulus.Pre_zscore), pStim_zscore, rhoStim_zscore, pRhoStim_zscore), sprintf('Hedges'' g=%.3f (%s)', hedges_g, effect_interp)});
    else
        text(ax5, 0.5, 0.5, 'No stimulus lick data', 'HorizontalAlignment', 'center');
        title(ax5, 'Stimulus Licks (Z-scored)');
        ax5.FontSize = 16;
    end
    
    % Reward timepoints (z-scored)
    ax6 = nexttile;
    if ~isempty(Fig1_Data.Reward.Pre) && isfield(Fig1_Data.Reward, 'Pre_zscore')
        dscatter(Fig1_Data.Reward.Pre_zscore, Fig1_Data.Reward.Post_zscore, 'MSIZE', 30);
        colormap(ax6, gray);
        axis(ax6, 'square');
        % Set equal axis limits
        allData6 = [Fig1_Data.Reward.Pre_zscore; Fig1_Data.Reward.Post_zscore];
        dataLim6 = [min(allData6), max(allData6)];
        xlim(ax6, dataLim6);
        % Set nice axis ticks at half and full numbers
        setNiceAxisTicks(ax6, dataLim6);
        ylim(ax6, dataLim6);
        ref6 = refline(ax6, 1); 
        ref6.LineWidth = 3; 
        ref6.Color = 'r';
        % Add fitted trend line
        fit6 = lsline(ax6);
        fit6.LineStyle = '--';
        fit6.LineWidth = 2;
        fit6.Color = [0.5 0.5 0.5];
        
        [pRewardTP_zscore, ~, ~] = signrank(Fig1_Data.Reward.Pre_zscore, Fig1_Data.Reward.Post_zscore);
        [rhoRewardTP_zscore, pRhoRewardTP_zscore] = corr(Fig1_Data.Reward.Pre_zscore, Fig1_Data.Reward.Post_zscore, 'Type', 'Spearman', 'Rows', 'complete');
        
        % Calculate Hedges' g effect size
        mean_pre = mean(Fig1_Data.Reward.Pre_zscore);
        mean_post = mean(Fig1_Data.Reward.Post_zscore);
        std_pre = std(Fig1_Data.Reward.Pre_zscore);
        std_post = std(Fig1_Data.Reward.Post_zscore);
        n_pre = length(Fig1_Data.Reward.Pre_zscore);
        n_post = length(Fig1_Data.Reward.Post_zscore);
        pooled_std = sqrt(((n_pre - 1) * std_pre^2 + (n_post - 1) * std_post^2) / (n_pre + n_post - 2));
        cohens_d = (mean_post - mean_pre) / pooled_std;
        correction_factor = 1 - (3 / (4 * (n_pre + n_post - 2) - 1));
        hedges_g = cohens_d * correction_factor;
        if abs(hedges_g) < 0.2; effect_interp = 'negligible'; elseif abs(hedges_g) < 0.5; effect_interp = 'small'; elseif abs(hedges_g) < 0.8; effect_interp = 'medium'; else; effect_interp = 'large'; end
        
        ax6.FontSize = 16;
        xlabel(ax6, 'Pre-Reward Entropy (Z)');
        ylabel(ax6, 'Post-Reward Entropy (Z)');
        title(ax6, {'Reward Timepoints (Z-scored)', sprintf('n=%d, p_{signrank}=%.3g, rho=%.3f, p_{rho}=%.3g', length(Fig1_Data.Reward.Pre_zscore), pRewardTP_zscore, rhoRewardTP_zscore, pRhoRewardTP_zscore), sprintf('Hedges'' g=%.3f (%s)', hedges_g, effect_interp)});
    else
        text(ax6, 0.5, 0.5, 'No reward timepoint data', 'HorizontalAlignment', 'center');
        title(ax6, 'Reward Timepoints (Z-scored)');
        ax6.FontSize = 16;
    end
    
    sgtitle(sprintf('All Lick Types: Pre vs Post Entropy - %s', position));
    
    %% ------------------------- Figure 2: All Lick Types Entropy violin plot --------------- %%
    
    fprintf('\n--- Generating Figure 2: Pre vs Post Entropy Violin Plot ---\n');
    
    % Calculate effect sizes (Hedges' g) for all lick types
    % Extraneous licks - raw data
    mean_pre_extr = mean(Fig1_Data.Extraneous.Pre);
    mean_post_extr = mean(Fig1_Data.Extraneous.Post);
    std_pre_extr = std(Fig1_Data.Extraneous.Pre);
    std_post_extr = std(Fig1_Data.Extraneous.Post);
    n_pre_extr = length(Fig1_Data.Extraneous.Pre);
    n_post_extr = length(Fig1_Data.Extraneous.Post);
    pooled_std_extr = sqrt(((n_pre_extr - 1) * std_pre_extr^2 + (n_post_extr - 1) * std_post_extr^2) / (n_pre_extr + n_post_extr - 2));
    cohens_d_extr = (mean_post_extr - mean_pre_extr) / pooled_std_extr;
    correction_factor_extr = 1 - (3 / (4 * (n_pre_extr + n_post_extr - 2) - 1));
    hedges_g_extr_raw = cohens_d_extr * correction_factor_extr;
    
    % Extraneous licks - z-scored data
    mean_pre_extr_z = mean(Fig1_Data.Extraneous.Pre_zscore);
    mean_post_extr_z = mean(Fig1_Data.Extraneous.Post_zscore);
    std_pre_extr_z = std(Fig1_Data.Extraneous.Pre_zscore);
    std_post_extr_z = std(Fig1_Data.Extraneous.Post_zscore);
    n_pre_extr_z = length(Fig1_Data.Extraneous.Pre_zscore);
    n_post_extr_z = length(Fig1_Data.Extraneous.Post_zscore);
    pooled_std_extr_z = sqrt(((n_pre_extr_z - 1) * std_pre_extr_z^2 + (n_post_extr_z - 1) * std_post_extr_z^2) / (n_pre_extr_z + n_post_extr_z - 2));
    cohens_d_extr_z = (mean_post_extr_z - mean_pre_extr_z) / pooled_std_extr_z;
    correction_factor_extr_z = 1 - (3 / (4 * (n_pre_extr_z + n_post_extr_z - 2) - 1));
    hedges_g_extr_zscore = cohens_d_extr_z * correction_factor_extr_z;
    
    % Stimulus licks effect sizes
    if ~isempty(Fig1_Data.Stimulus.Pre)
        mean_pre_stim = mean(Fig1_Data.Stimulus.Pre);
        mean_post_stim = mean(Fig1_Data.Stimulus.Post);
        std_pre_stim = std(Fig1_Data.Stimulus.Pre);
        std_post_stim = std(Fig1_Data.Stimulus.Post);
        n_pre_stim = length(Fig1_Data.Stimulus.Pre);
        n_post_stim = length(Fig1_Data.Stimulus.Post);
        pooled_std_stim = sqrt(((n_pre_stim - 1) * std_pre_stim^2 + (n_post_stim - 1) * std_post_stim^2) / (n_pre_stim + n_post_stim - 2));
        cohens_d_stim = (mean_post_stim - mean_pre_stim) / pooled_std_stim;
        correction_factor_stim = 1 - (3 / (4 * (n_pre_stim + n_post_stim - 2) - 1));
        hedges_g_stim_raw = cohens_d_stim * correction_factor_stim;
        
        mean_pre_stim_z = mean(Fig1_Data.Stimulus.Pre_zscore);
        mean_post_stim_z = mean(Fig1_Data.Stimulus.Post_zscore);
        std_pre_stim_z = std(Fig1_Data.Stimulus.Pre_zscore);
        std_post_stim_z = std(Fig1_Data.Stimulus.Post_zscore);
        n_pre_stim_z = length(Fig1_Data.Stimulus.Pre_zscore);
        n_post_stim_z = length(Fig1_Data.Stimulus.Post_zscore);
        pooled_std_stim_z = sqrt(((n_pre_stim_z - 1) * std_pre_stim_z^2 + (n_post_stim_z - 1) * std_post_stim_z^2) / (n_pre_stim_z + n_post_stim_z - 2));
        cohens_d_stim_z = (mean_post_stim_z - mean_pre_stim_z) / pooled_std_stim_z;
        correction_factor_stim_z = 1 - (3 / (4 * (n_pre_stim_z + n_post_stim_z - 2) - 1));
        hedges_g_stim_zscore = cohens_d_stim_z * correction_factor_stim_z;
    else
        hedges_g_stim_raw = NaN;
        hedges_g_stim_zscore = NaN;
    end
    
    % Reward licks effect sizes
    if ~isempty(Fig1_Data.Reward.Pre)
        mean_pre_reward = mean(Fig1_Data.Reward.Pre);
        mean_post_reward = mean(Fig1_Data.Reward.Post);
        std_pre_reward = std(Fig1_Data.Reward.Pre);
        std_post_reward = std(Fig1_Data.Reward.Post);
        n_pre_reward = length(Fig1_Data.Reward.Pre);
        n_post_reward = length(Fig1_Data.Reward.Post);
        pooled_std_reward = sqrt(((n_pre_reward - 1) * std_pre_reward^2 + (n_post_reward - 1) * std_post_reward^2) / (n_pre_reward + n_post_reward - 2));
        cohens_d_reward = (mean_post_reward - mean_pre_reward) / pooled_std_reward;
        correction_factor_reward = 1 - (3 / (4 * (n_pre_reward + n_post_reward - 2) - 1));
        hedges_g_reward_raw = cohens_d_reward * correction_factor_reward;
        
        mean_pre_reward_z = mean(Fig1_Data.Reward.Pre_zscore);
        mean_post_reward_z = mean(Fig1_Data.Reward.Post_zscore);
        std_pre_reward_z = std(Fig1_Data.Reward.Pre_zscore);
        std_post_reward_z = std(Fig1_Data.Reward.Post_zscore);
        n_pre_reward_z = length(Fig1_Data.Reward.Pre_zscore);
        n_post_reward_z = length(Fig1_Data.Reward.Post_zscore);
        pooled_std_reward_z = sqrt(((n_pre_reward_z - 1) * std_pre_reward_z^2 + (n_post_reward_z - 1) * std_post_reward_z^2) / (n_pre_reward_z + n_post_reward_z - 2));
        cohens_d_reward_z = (mean_post_reward_z - mean_pre_reward_z) / pooled_std_reward_z;
        correction_factor_reward_z = 1 - (3 / (4 * (n_pre_reward_z + n_post_reward_z - 2) - 1));
        hedges_g_reward_zscore = cohens_d_reward_z * correction_factor_reward_z;
    else
        hedges_g_reward_raw = NaN;
        hedges_g_reward_zscore = NaN;
    end
    
    fprintf('Effect sizes calculated:\n');
    fprintf('  Extraneous (Raw): g = %.3f\n', hedges_g_extr_raw);
    fprintf('  Extraneous (Z-scored): g = %.3f\n', hedges_g_extr_zscore);
    if ~isempty(Fig1_Data.Stimulus.Pre)
        fprintf('  Stimulus (Raw): g = %.3f\n', hedges_g_stim_raw);
        fprintf('  Stimulus (Z-scored): g = %.3f\n', hedges_g_stim_zscore);
    end
    if ~isempty(Fig1_Data.Reward.Pre)
        fprintf('  Reward (Raw): g = %.3f\n', hedges_g_reward_raw);
        fprintf('  Reward (Z-scored): g = %.3f\n', hedges_g_reward_zscore);
    end
    
    figure('Color', 'w', 'Name', sprintf('All Lick Types: Pre vs Post Entropy Violin Plot - %s', position));
    axesFont = 18;
    tl = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Calculate p-values for violin plot titles
    [pAct, ~, ~] = signrank(Fig1_Data.Extraneous.Pre, Fig1_Data.Extraneous.Post);
    [pAct_zscore, ~, ~] = signrank(Fig1_Data.Extraneous.Pre_zscore, Fig1_Data.Extraneous.Post_zscore);
    
    if ~isempty(Fig1_Data.Stimulus.Pre)
        [pStim, ~, ~] = signrank(Fig1_Data.Stimulus.Pre, Fig1_Data.Stimulus.Post);
        [pStim_zscore, ~, ~] = signrank(Fig1_Data.Stimulus.Pre_zscore, Fig1_Data.Stimulus.Post_zscore);
    end
    
    if ~isempty(Fig1_Data.Reward.Pre)
        [pRewardTP, ~, ~] = signrank(Fig1_Data.Reward.Pre, Fig1_Data.Reward.Post);
        [pRewardTP_zscore, ~, ~] = signrank(Fig1_Data.Reward.Pre_zscore, Fig1_Data.Reward.Post_zscore);
    end
    
    % Row 1: Extraneous Licks
    % Raw data (left)
    ax1 = nexttile;
    daviolinplot({[Fig1_Data.Extraneous.Pre, Fig1_Data.Extraneous.Post]}, ...
                 'violin', 'full', 'groups', ones(numel(Fig1_Data.Extraneous.Pre), 1), ...
                 'colors', [0 0 0], 'box', 3, 'boxcolor', 'w', ...
                 'scatter', 0, 'linkline', 1, ...
                 'xtlabels', {"Before", "After"});
    set(ax1, 'FontSize', axesFont);
    ylabel(ax1, 'Entropy (Raw)');
    title(ax1, sprintf('Extraneous Licks (Raw): p = %.3g, g = %.3f, n = %d', pAct, hedges_g_extr_raw, length(Fig1_Data.Extraneous.Pre)));
    
    % Z-scored data (right)
    ax2 = nexttile;
    daviolinplot({[Fig1_Data.Extraneous.Pre_zscore, Fig1_Data.Extraneous.Post_zscore]}, ...
                 'violin', 'full', 'groups', ones(numel(Fig1_Data.Extraneous.Pre_zscore), 1), ...
                 'colors', [0 0 0], 'box', 3, 'boxcolor', 'w', ...
                 'scatter', 0, 'linkline', 1, ...
                 'xtlabels', {"Before", "After"});
    set(ax2, 'FontSize', axesFont);
    ylabel(ax2, 'Entropy (Z-scored)');
    title(ax2, sprintf('Extraneous Licks (Z-scored): p = %.3g, g = %.3f, n = %d', pAct_zscore, hedges_g_extr_zscore, length(Fig1_Data.Extraneous.Pre_zscore)));
    
    % Row 2: Stimulus Licks
    % Raw data (left)
    ax3 = nexttile;
    if ~isempty(Fig1_Data.Stimulus.Pre)
        daviolinplot({[Fig1_Data.Stimulus.Pre, Fig1_Data.Stimulus.Post]}, ...
                     'violin', 'full', 'groups', ones(numel(Fig1_Data.Stimulus.Pre), 1), ...
                     'colors', [0.2 0.6 0.8], 'box', 3, 'boxcolor', 'w', ...
                     'scatter', 0, 'linkline', 1, ...
                     'xtlabels', {"Before", "After"});
        set(ax3, 'FontSize', axesFont);
        ylabel(ax3, 'Entropy (Raw)');
        title(ax3, sprintf('Stimulus Licks (Raw): p = %.3g, g = %.3f, n = %d', pStim, hedges_g_stim_raw, length(Fig1_Data.Stimulus.Pre)));
    else
        text(ax3, 0.5, 0.5, 'No stimulus lick data', 'HorizontalAlignment', 'center', 'FontSize', 14);
        set(ax3, 'XTick', [], 'YTick', []);
        title(ax3, 'Stimulus Licks (Raw): No data');
    end
    
    % Z-scored data (right)
    ax4 = nexttile;
    if ~isempty(Fig1_Data.Stimulus.Pre) && isfield(Fig1_Data.Stimulus, 'Pre_zscore')
        daviolinplot({[Fig1_Data.Stimulus.Pre_zscore, Fig1_Data.Stimulus.Post_zscore]}, ...
                     'violin', 'full', 'groups', ones(numel(Fig1_Data.Stimulus.Pre_zscore), 1), ...
                     'colors', [0.2 0.6 0.8], 'box', 3, 'boxcolor', 'w', ...
                     'scatter', 0, 'linkline', 1, ...
                     'xtlabels', {"Before", "After"});
        set(ax4, 'FontSize', axesFont);
        ylabel(ax4, 'Entropy (Z-scored)');
        title(ax4, sprintf('Stimulus Licks (Z-scored): p = %.3g, g = %.3f, n = %d', pStim_zscore, hedges_g_stim_zscore, length(Fig1_Data.Stimulus.Pre_zscore)));
    else
        text(ax4, 0.5, 0.5, 'No stimulus lick data', 'HorizontalAlignment', 'center', 'FontSize', 14);
        set(ax4, 'XTick', [], 'YTick', []);
        title(ax4, 'Stimulus Licks (Z-scored): No data');
    end
    
    % Row 3: Reward Timepoints
    % Raw data (left)
    ax5 = nexttile;
    if ~isempty(Fig1_Data.Reward.Pre)
        daviolinplot({[Fig1_Data.Reward.Pre, Fig1_Data.Reward.Post]}, ...
                     'violin', 'full', 'groups', ones(numel(Fig1_Data.Reward.Pre), 1), ...
                     'colors', [0.8 0.4 0.2], 'box', 3, 'boxcolor', 'w', ...
                     'scatter', 0, 'linkline', 1, ...
                     'xtlabels', {"Before", "After"});
        set(ax5, 'FontSize', axesFont);
        ylabel(ax5, 'Entropy (Raw)');
        title(ax5, sprintf('Reward Timepoints (Raw): p = %.3g, g = %.3f, n = %d', pRewardTP, hedges_g_reward_raw, length(Fig1_Data.Reward.Pre)));
    else
        text(ax5, 0.5, 0.5, 'No reward timepoint data', 'HorizontalAlignment', 'center', 'FontSize', 14);
        set(ax5, 'XTick', [], 'YTick', []);
        title(ax5, 'Reward Timepoints (Raw): No data');
    end
    
    % Z-scored data (right)
    ax6 = nexttile;
    if ~isempty(Fig1_Data.Reward.Pre) && isfield(Fig1_Data.Reward, 'Pre_zscore')
        daviolinplot({[Fig1_Data.Reward.Pre_zscore, Fig1_Data.Reward.Post_zscore]}, ...
                     'violin', 'full', 'groups', ones(numel(Fig1_Data.Reward.Pre_zscore), 1), ...
                     'colors', [0.8 0.4 0.2], 'box', 3, 'boxcolor', 'w', ...
                     'scatter', 0, 'linkline', 1, ...
                     'xtlabels', {"Before", "After"});
        set(ax6, 'FontSize', axesFont);
        ylabel(ax6, 'Entropy (Z-scored)');
        title(ax6, sprintf('Reward Timepoints (Z-scored): p = %.3g, g = %.3f, n = %d', pRewardTP_zscore, hedges_g_reward_zscore, length(Fig1_Data.Reward.Pre_zscore)));
    else
        text(ax6, 0.5, 0.5, 'No reward timepoint data', 'HorizontalAlignment', 'center', 'FontSize', 14);
        set(ax6, 'XTick', [], 'YTick', []);
        title(ax6, 'Reward Timepoints (Z-scored): No data');
    end
    
    sgtitle(sprintf('All Lick Types: Pre vs Post Entropy - %s', position));
    
    %% ------------------------- Figure 3: All Lick Types Entropy Ratio Histogram ------- %%
    
    fprintf('\n--- Generating Figure 3: Post/Pre Entropy Ratio Histogram ---\n');
    
    figure('Color', 'w', 'Name', sprintf('All Lick Types: Post/Pre Entropy Ratio Histogram - %s', position));
    tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Row 1: Raw data
    % Extraneous licks (raw)
    ax1 = nexttile;
    if isfield(Fig3_Data, 'Extraneous') && ~isempty(Fig3_Data.Extraneous.Ratio)
        histogram(ax1, Fig3_Data.Extraneous.Ratio, 'Normalization', 'probability', ...
                  'BinMethod', 'sturges', 'FaceColor', [0 0 0], ...
                  'EdgeColor', 'none', 'FaceAlpha', 0.7);
        [pRatioAct] = signrank(Fig3_Data.Extraneous.Ratio - 1);
        set(ax1, 'FontSize', 16);
        box(ax1, 'on');
        xlabel(ax1, 'Post / Pre Entropy');
        ylabel(ax1, 'Probability');
        title(ax1, {'Extraneous Licks (Raw)', sprintf('p = %.3g, n = %d', pRatioAct, length(Fig3_Data.Extraneous.Ratio))});
    else
        text(ax1, 0.5, 0.5, 'No extraneous lick data', 'HorizontalAlignment', 'center', 'FontSize', 14);
        set(ax1, 'XTick', [], 'YTick', []);
        title(ax1, 'Extraneous Licks (Raw): No data');
    end
    
    % Stimulus licks (raw)
    ax2 = nexttile;
    if isfield(Fig3_Data, 'Stimulus') && ~isempty(Fig3_Data.Stimulus.Ratio)
        histogram(ax2, Fig3_Data.Stimulus.Ratio, 'Normalization', 'probability', ...
                  'BinMethod', 'sturges', 'FaceColor', [0.2 0.6 0.8], ...
                  'EdgeColor', 'none', 'FaceAlpha', 0.7);
        [pRatioStim] = signrank(Fig3_Data.Stimulus.Ratio - 1);
        set(ax2, 'FontSize', 16);
        box(ax2, 'on');
        xlabel(ax2, 'Post / Pre Entropy');
        ylabel(ax2, 'Probability');
        title(ax2, {'Stimulus Licks (Raw)', sprintf('p = %.3g, n = %d', pRatioStim, length(Fig3_Data.Stimulus.Ratio))});
    else
        text(ax2, 0.5, 0.5, 'No stimulus lick data', 'HorizontalAlignment', 'center', 'FontSize', 14);
        set(ax2, 'XTick', [], 'YTick', []);
        title(ax2, 'Stimulus Licks (Raw): No data');
    end
    
    % Reward timepoints (raw)
    ax3 = nexttile;
    if isfield(Fig3_Data, 'Reward') && ~isempty(Fig3_Data.Reward.Ratio)
        histogram(ax3, Fig3_Data.Reward.Ratio, 'Normalization', 'probability', ...
                  'BinMethod', 'sturges', 'FaceColor', [0.8 0.4 0.2], ...
                  'EdgeColor', 'none', 'FaceAlpha', 0.7);
        [pRatioRewardTP] = signrank(Fig3_Data.Reward.Ratio - 1);
        set(ax3, 'FontSize', 16);
        box(ax3, 'on');
        xlabel(ax3, 'Post / Pre Entropy');
        ylabel(ax3, 'Probability');
        title(ax3, {'Reward Timepoints (Raw)', sprintf('p = %.3g, n = %d', pRatioRewardTP, length(Fig3_Data.Reward.Ratio))});
    else
        text(ax3, 0.5, 0.5, 'No reward timepoint data', 'HorizontalAlignment', 'center', 'FontSize', 14);
        set(ax3, 'XTick', [], 'YTick', []);
        title(ax3, 'Reward Timepoints (Raw): No data');
    end
    
    % Row 2: Z-scored data
    % Extraneous licks (z-scored)
    ax4 = nexttile;
    if isfield(Fig3_Data, 'Extraneous') && ~isempty(Fig3_Data.Extraneous.Ratio_zscore)
        histogram(ax4, Fig3_Data.Extraneous.Ratio_zscore, 'Normalization', 'probability', ...
                  'BinMethod', 'sturges', 'FaceColor', [0 0 0], ...
                  'EdgeColor', 'none', 'FaceAlpha', 0.7);
        [pRatioAct_zscore] = signrank(Fig3_Data.Extraneous.Ratio_zscore - 1);
        set(ax4, 'FontSize', 16);
        box(ax4, 'on');
        xlabel(ax4, 'Post / Pre Entropy (Z)');
        ylabel(ax4, 'Probability');
        title(ax4, {'Extraneous Licks (Z-scored)', sprintf('p = %.3g, n = %d', pRatioAct_zscore, length(Fig3_Data.Extraneous.Ratio_zscore))});
    else
        text(ax4, 0.5, 0.5, 'No extraneous lick data', 'HorizontalAlignment', 'center', 'FontSize', 14);
        set(ax4, 'XTick', [], 'YTick', []);
        title(ax4, 'Extraneous Licks (Z-scored): No data');
    end
    
    % Stimulus licks (z-scored)
    ax5 = nexttile;
    if isfield(Fig3_Data, 'Stimulus') && ~isempty(Fig3_Data.Stimulus.Ratio_zscore)
        histogram(ax5, Fig3_Data.Stimulus.Ratio_zscore, 'Normalization', 'probability', ...
                  'BinMethod', 'sturges', 'FaceColor', [0.2 0.6 0.8], ...
                  'EdgeColor', 'none', 'FaceAlpha', 0.7);
        [pRatioStim_zscore] = signrank(Fig3_Data.Stimulus.Ratio_zscore - 1);
        set(ax5, 'FontSize', 16);
        box(ax5, 'on');
        xlabel(ax5, 'Post / Pre Entropy (Z)');
        ylabel(ax5, 'Probability');
        title(ax5, {'Stimulus Licks (Z-scored)', sprintf('p = %.3g, n = %d', pRatioStim_zscore, length(Fig3_Data.Stimulus.Ratio_zscore))});
    else
        text(ax5, 0.5, 0.5, 'No stimulus lick data', 'HorizontalAlignment', 'center', 'FontSize', 14);
        set(ax5, 'XTick', [], 'YTick', []);
        title(ax5, 'Stimulus Licks (Z-scored): No data');
    end
    
    % Reward timepoints (z-scored)
    ax6 = nexttile;
    if isfield(Fig3_Data, 'Reward') && ~isempty(Fig3_Data.Reward.Ratio_zscore)
        histogram(ax6, Fig3_Data.Reward.Ratio_zscore, 'Normalization', 'probability', ...
                  'BinMethod', 'sturges', 'FaceColor', [0.8 0.4 0.2], ...
                  'EdgeColor', 'none', 'FaceAlpha', 0.7);
        [pRatioRewardTP_zscore] = signrank(Fig3_Data.Reward.Ratio_zscore - 1);
        set(ax6, 'FontSize', 16);
        box(ax6, 'on');
        xlabel(ax6, 'Post / Pre Entropy (Z)');
        ylabel(ax6, 'Probability');
        title(ax6, {'Reward Timepoints (Z-scored)', sprintf('p = %.3g, n = %d', pRatioRewardTP_zscore, length(Fig3_Data.Reward.Ratio_zscore))});
    else
        text(ax6, 0.5, 0.5, 'No reward timepoint data', 'HorizontalAlignment', 'center', 'FontSize', 14);
        set(ax6, 'XTick', [], 'YTick', []);
        title(ax6, 'Reward Timepoints (Z-scored): No data');
    end
    
    sgtitle(sprintf('All Lick Types: Post/Pre Entropy Ratios - %s', position));
    
    %% ------------------------- Figure 4: All Lick Types Entropy Change vs Pre --------- %%
    
    fprintf('\n--- Generating Figure 4: Entropy Change vs Pre-Lick Values ---\n');
    
    figure('Color', 'w', 'Name', sprintf('All Lick Types: Entropy Change vs Pre-Lick Values - %s', position));
    tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Row 1: Raw data
    % Extraneous licks (raw)
    ax1 = nexttile;
    dscatter(Fig1_Data.Extraneous.Pre, Fig4_Data.Extraneous.Change, 'MSIZE', 30);
    colormap(ax1, gray);
    lsline(ax1);
    axis(ax1, 'square'); 
    grid(ax1, 'on');
    [rhoAct, pRhoAct] = corr(Fig1_Data.Extraneous.Pre, Fig4_Data.Extraneous.Change, ...
                             'Type', 'Spearman', 'Rows', 'complete');
    set(ax1, 'FontSize', 16);
    xlabel(ax1, 'Pre Entropy');
    ylabel(ax1, '\Delta Entropy');
    title(ax1, {'Extraneous Licks (Raw)', sprintf('rho = %.2f, p = %.3g', rhoAct, pRhoAct), sprintf('n = %d', length(Fig1_Data.Extraneous.Pre))});
    
    % Stimulus licks (raw)
    ax2 = nexttile;
    if ~isempty(Fig1_Data.Stimulus.Pre) && isfield(Fig4_Data, 'Stimulus') && ~isempty(Fig4_Data.Stimulus.Change)
        dscatter(Fig1_Data.Stimulus.Pre, Fig4_Data.Stimulus.Change, 'MSIZE', 30);
        colormap(ax2, gray);
        lsline(ax2);
        axis(ax2, 'square'); 
        grid(ax2, 'on');
        [rhoStim, pRhoStim] = corr(Fig1_Data.Stimulus.Pre, Fig4_Data.Stimulus.Change, ...
                                   'Type', 'Spearman', 'Rows', 'complete');
        set(ax2, 'FontSize', 16);
        xlabel(ax2, 'Pre Entropy');
        ylabel(ax2, '\Delta Entropy');
        title(ax2, {'Stimulus Licks (Raw)', sprintf('rho = %.2f, p = %.3g', rhoStim, pRhoStim), sprintf('n = %d', length(Fig1_Data.Stimulus.Pre))});
    else
        text(ax2, 0.5, 0.5, 'No stimulus lick data', 'HorizontalAlignment', 'center', 'FontSize', 14);
        set(ax2, 'XTick', [], 'YTick', []);
        title(ax2, 'Stimulus Licks (Raw): No data');
    end
    
    % Reward timepoints (raw)
    ax3 = nexttile;
    if ~isempty(Fig1_Data.Reward.Pre) && isfield(Fig4_Data, 'Reward') && ~isempty(Fig4_Data.Reward.Change)
        dscatter(Fig1_Data.Reward.Pre, Fig4_Data.Reward.Change, 'MSIZE', 30);
        colormap(ax3, gray);
        lsline(ax3);
        axis(ax3, 'square'); 
        grid(ax3, 'on');
        [rhoRewardTP, pRhoRewardTP] = corr(Fig1_Data.Reward.Pre, Fig4_Data.Reward.Change, ...
                                           'Type', 'Spearman', 'Rows', 'complete');
        set(ax3, 'FontSize', 16);
        xlabel(ax3, 'Pre Entropy');
        ylabel(ax3, '\Delta Entropy');
        title(ax3, {'Reward Timepoints (Raw)', sprintf('rho = %.2f, p = %.3g', rhoRewardTP, pRhoRewardTP), sprintf('n = %d', length(Fig1_Data.Reward.Pre))});
    else
        text(ax3, 0.5, 0.5, 'No reward timepoint data', 'HorizontalAlignment', 'center', 'FontSize', 14);
        set(ax3, 'XTick', [], 'YTick', []);
        title(ax3, 'Reward Timepoints (Raw): No data');
    end
    
    % Row 2: Z-scored data
    % Extraneous licks (z-scored)
    ax4 = nexttile;
    dscatter(Fig1_Data.Extraneous.Pre_zscore, Fig4_Data.Extraneous.Change_zscore, 'MSIZE', 30);
    colormap(ax4, gray);
    lsline(ax4);
    axis(ax4, 'square'); 
    grid(ax4, 'on');
    [rhoAct_zscore, pRhoAct_zscore] = corr(Fig1_Data.Extraneous.Pre_zscore, Fig4_Data.Extraneous.Change_zscore, ...
                                           'Type', 'Spearman', 'Rows', 'complete');
    set(ax4, 'FontSize', 16);
    xlabel(ax4, 'Pre Entropy (Z)');
    ylabel(ax4, '\Delta Entropy (Z)');
    title(ax4, {'Extraneous Licks (Z-scored)', sprintf('rho = %.2f, p = %.3g', rhoAct_zscore, pRhoAct_zscore), sprintf('n = %d', length(Fig1_Data.Extraneous.Pre_zscore))});
    
    % Stimulus licks (z-scored)
    ax5 = nexttile;
    if ~isempty(Fig1_Data.Stimulus.Pre) && isfield(Fig1_Data.Stimulus, 'Pre_zscore') && isfield(Fig4_Data, 'Stimulus') && ~isempty(Fig4_Data.Stimulus.Change_zscore)
        dscatter(Fig1_Data.Stimulus.Pre_zscore, Fig4_Data.Stimulus.Change_zscore, 'MSIZE', 30);
        colormap(ax5, gray);
        lsline(ax5);
        axis(ax5, 'square'); 
        grid(ax5, 'on');
        [rhoStim_zscore, pRhoStim_zscore] = corr(Fig1_Data.Stimulus.Pre_zscore, Fig4_Data.Stimulus.Change_zscore, ...
                                                 'Type', 'Spearman', 'Rows', 'complete');
        set(ax5, 'FontSize', 16);
        xlabel(ax5, 'Pre Entropy (Z)');
        ylabel(ax5, '\Delta Entropy (Z)');
        title(ax5, {'Stimulus Licks (Z-scored)', sprintf('rho = %.2f, p = %.3g', rhoStim_zscore, pRhoStim_zscore), sprintf('n = %d', length(Fig1_Data.Stimulus.Pre_zscore))});
    else
        text(ax5, 0.5, 0.5, 'No stimulus lick data', 'HorizontalAlignment', 'center', 'FontSize', 14);
        set(ax5, 'XTick', [], 'YTick', []);
        title(ax5, 'Stimulus Licks (Z-scored): No data');
    end
    
    % Reward timepoints (z-scored)
    ax6 = nexttile;
    if ~isempty(Fig1_Data.Reward.Pre) && isfield(Fig1_Data.Reward, 'Pre_zscore') && isfield(Fig4_Data, 'Reward') && ~isempty(Fig4_Data.Reward.Change_zscore)
        dscatter(Fig1_Data.Reward.Pre_zscore, Fig4_Data.Reward.Change_zscore, 'MSIZE', 30);
        colormap(ax6, gray);
        lsline(ax6);
        axis(ax6, 'square'); 
        grid(ax6, 'on');
        [rhoRewardTP_zscore, pRhoRewardTP_zscore] = corr(Fig1_Data.Reward.Pre_zscore, Fig4_Data.Reward.Change_zscore, ...
                                                         'Type', 'Spearman', 'Rows', 'complete');
        set(ax6, 'FontSize', 16);
        xlabel(ax6, 'Pre Entropy (Z)');
        ylabel(ax6, '\Delta Entropy (Z)');
        title(ax6, {'Reward Timepoints (Z-scored)', sprintf('rho = %.2f, p = %.3g', rhoRewardTP_zscore, pRhoRewardTP_zscore), sprintf('n = %d', length(Fig1_Data.Reward.Pre_zscore))});
    else
        text(ax6, 0.5, 0.5, 'No reward timepoint data', 'HorizontalAlignment', 'center', 'FontSize', 14);
        set(ax6, 'XTick', [], 'YTick', []);
        title(ax6, 'Reward Timepoints (Z-scored): No data');
    end
    
    sgtitle(sprintf('All Lick Types: Entropy Change vs Pre-Lick Values - %s', position));
    
    %% --------------- Figure 5: Extraneous Licks - Lick-Triggered Average (LTA) --------------- %%
    fprintf('\n--- Generating Figure 5: Extraneous Lick-Triggered Average ---\n');
    
    % Create Extraneous LTA plot
    if ~isempty(Fig5_Data.allTrials)
        [~, LickTimeSortedIndex] = sort(Fig5_Data.lickTimes);
        
        figure('Color', 'w', 'Name', sprintf('Extraneous Licks: Lick-Triggered Average Analysis - %s', position));
        f = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        % Raw data heatmap (left top)
        ax1 = nexttile([2 1]);
        imagesc(ax1, Fig5_Data.timeAxis, 1:size(Fig5_Data.allTrials, 1), ...
                Fig5_Data.allTrials(LickTimeSortedIndex, :));
        colormap(ax1, params.cmap);
        colorbar(ax1);
        ylabel(ax1, 'Trials');
        xlabel(ax1, 'Time to extraneous lick (s)');
        xlim(ax1, [-3 3]);
        title(ax1, 'Raw Data (sorted by lick time)');
        
        % Z-scored data heatmap (right top)
        ax2 = nexttile([2 1]);
        imagesc(ax2, Fig5_Data.timeAxis, 1:size(Fig5_Data.allTrials_zscore, 1), ...
                Fig5_Data.allTrials_zscore(LickTimeSortedIndex, :));
        colormap(ax2, params.cmap);
        colorbar(ax2);
        ylabel(ax2, 'Trials');
        xlabel(ax2, 'Time to extraneous lick (s)');
        xlim(ax2, [-3 3]);
        title(ax2, 'Z-scored Data (sorted by lick time)');
        
        % Raw data average trace (left bottom)
        ax3 = nexttile;
        LTA_mean = mean(Fig5_Data.allTrials, 1, 'omitnan');
        plot(ax3, Fig5_Data.timeAxis, LTA_mean, 'LineWidth', 4, 'Color', 'k');
        xlim(ax3, [-3 3]);
        xlabel(ax3, 'Time to extraneous lick (s)');
        ylabel(ax3, 'Entropy (Raw)');
        title(ax3, 'Average Raw Response');
        
        % Z-scored average trace (right bottom)
        ax4 = nexttile;
        LTA_mean_zscore = mean(Fig5_Data.allTrials_zscore, 1, 'omitnan');
        plot(ax4, Fig5_Data.timeAxis, LTA_mean_zscore, 'LineWidth', 4, 'Color', 'k');
        xlim(ax4, [-3 3]);
        xlabel(ax4, 'Time to extraneous lick (s)');
        ylabel(ax4, 'Entropy (Z-scored)');
        title(ax4, 'Average Z-scored Response');
        
        fontsize(f, 20, 'points');
        sgtitle(sprintf('Extraneous Licks: Lick-triggered averages - %s (%d licks)', position, size(Fig5_Data.allTrials, 1)));
        
        fprintf('Extraneous lick-triggered average for %s: %d licks analyzed\n', position, size(Fig5_Data.allTrials, 1));
    else
        fprintf('Warning: No data available for extraneous lick-triggered average analysis for %s\n', position);
    end
    
    %% --------------- Figure 6: Stimulus Licks - Lick-Triggered Average (LTA) --------------- %%
    fprintf('\n--- Generating Figure 6: Stimulus Lick-Triggered Average ---\n');
    
    % Create LTA plot
    if ~isempty(Fig6_Data.allTrials)
        [~, ReactionTimeSortedIndex] = sort(Fig6_Data.reactionTimes);
        
        figure('Color', 'w', 'Name', sprintf('Stimulus Licks: Lick-Triggered Average Analysis - %s', position));
        f = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        % Raw data heatmap (left top)
        ax1 = nexttile([2 1]);
        imagesc(ax1, Fig6_Data.timeAxis, ...
                1:size(Fig6_Data.allTrials, 1), Fig6_Data.allTrials(ReactionTimeSortedIndex, :));
        colormap(ax1, params.cmap);
        colorbar(ax1);
        ylabel(ax1, 'Trials');
        xlabel(ax1, 'Time to first lick (s)');
        xlim(ax1, [-3 3]);
        title(ax1, 'Raw Data (sorted by reaction time)');
        
        % Z-scored data heatmap (right top)
        ax2 = nexttile([2 1]);
        imagesc(ax2, Fig6_Data.timeAxis, ...
                1:size(Fig6_Data.allTrials_zscore, 1), Fig6_Data.allTrials_zscore(ReactionTimeSortedIndex, :));
        colormap(ax2, params.cmap);
        colorbar(ax2);
        ylabel(ax2, 'Trials');
        xlabel(ax2, 'Time to first lick (s)');
        xlim(ax2, [-3 3]);
        title(ax2, 'Z-scored Data (sorted by reaction time)');
        
        % Raw data average trace (left bottom)
        ax3 = nexttile;
        LTA_mean = mean(Fig6_Data.allTrials, 1);
        plot(ax3, Fig6_Data.timeAxis, LTA_mean, ...
             'LineWidth', 4, 'Color', 'k');
        xlim(ax3, [-3 3]);
        xlabel(ax3, 'Time to first lick (s)');
        ylabel(ax3, 'Entropy (Raw)');
        title(ax3, 'Average Raw Response');
        
        % Z-scored average trace (right bottom)
        ax4 = nexttile;
        LTA_mean_zscore = mean(Fig6_Data.allTrials_zscore, 1);
        plot(ax4, Fig6_Data.timeAxis, LTA_mean_zscore, ...
             'LineWidth', 4, 'Color', 'k');
        xlim(ax4, [-3 3]);
        xlabel(ax4, 'Time to first lick (s)');
        ylabel(ax4, 'Entropy (Z-scored)');
        title(ax4, 'Average Z-scored Response');
        
        fontsize(f, 24, 'points');
        sgtitle(sprintf('Stimulus Licks: Lick-triggered averages - %s (%d trials, all valid licks)', position, size(Fig6_Data.allTrials, 1)));
        
        fprintf('Lick-triggered average for %s: %d trials analyzed (all valid licks)\n', position, size(Fig6_Data.allTrials, 1));
    else
        fprintf('Warning: No data available for lick-triggered average analysis for %s\n', position);
    end
    
    %% ----------- Figure 7: Stimulus Licks - Lick-Triggered Modulation (LTM) --------------- %%
    fprintf('\n--- Generating Figure 7: Stimulus Lick-Triggered Modulation ---\n');
    
    % Create LTM plot
    if ~isempty(Fig7_Data.allTrials)
        [~, ReactionTimeSortedIndex] = sort(Fig7_Data.reactionTimes);
        
        figure('Color', 'w', 'Name', sprintf('Stimulus Licks: Lick-Triggered Modulation Analysis - %s', position));
        f = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        % Raw data heatmap (left top)
        ax1 = nexttile([2 1]);
        imagesc(ax1, Fig7_Data.timeAxis, ...
                1:size(Fig7_Data.allTrials, 1), Fig7_Data.allTrials(ReactionTimeSortedIndex, :));
        colormap(ax1, params.cmap);
        colorbar(ax1);
        clim(ax1, [-2 2]);
        ylabel(ax1, 'Trials');
        xlabel(ax1, 'Time to first lick (s)');
        xlim(ax1, [-4 4]);
        title(ax1, 'Raw Data (baseline-corrected)');
        
        % Z-scored data heatmap (right top)
        ax2 = nexttile([2 1]);
        imagesc(ax2, Fig7_Data.timeAxis, ...
                1:size(Fig7_Data.allTrials_zscore, 1), Fig7_Data.allTrials_zscore(ReactionTimeSortedIndex, :));
        colormap(ax2, params.cmap);
        colorbar(ax2);
        clim(ax2, [-2 2]);
        ylabel(ax2, 'Trials');
        xlabel(ax2, 'Time to first lick (s)');
        xlim(ax2, [-4 4]);
        title(ax2, 'Z-scored Data (baseline-corrected)');
        
        % Raw data average trace (left bottom)
        ax3 = nexttile;
        LTM_mean = mean(Fig7_Data.allTrials, 1);
        plot(ax3, Fig7_Data.timeAxis, LTM_mean, ...
             'LineWidth', 4, 'Color', 'k');
        xlim(ax3, [-4 4]);
        xlabel(ax3, 'Time to first lick (s)');
        ylabel(ax3, 'Entropy (Raw)');
        title(ax3, 'Average Raw Modulation');
        
        % Z-scored average trace (right bottom)
        ax4 = nexttile;
        LTM_mean_zscore = mean(Fig7_Data.allTrials_zscore, 1);
        plot(ax4, Fig7_Data.timeAxis, LTM_mean_zscore, ...
             'LineWidth', 4, 'Color', 'k');
        xlim(ax4, [-4 4]);
        xlabel(ax4, 'Time to first lick (s)');
        ylabel(ax4, 'Entropy (Z-scored)');
        title(ax4, 'Average Z-scored Modulation');
        
        fontsize(f, 24, 'points');
        sgtitle(sprintf('Stimulus Licks: Lick-triggered modulation - %s (%d trials)', position, size(Fig7_Data.allTrials, 1)));
        
        fprintf('Lick-triggered modulation for %s: %d trials analyzed\n', position, size(Fig7_Data.allTrials, 1));
    else
        fprintf('Warning: No data available for lick-triggered modulation analysis for %s\n', position);
    end
    
    %% ------------- Figure 8: Stimulus Licks - Reward-Triggered Average (RTA) -------------- %%
    fprintf('\n--- Generating Figure 8: Reward-Triggered Average ---\n');
    
    % Create RTA plot
    if ~isempty(Fig8_Data.allTrials)
        [~, ReactionTimeSortedIndex] = sort(Fig8_Data.reactionTimes);
        
        figure('Color', 'w', 'Name', sprintf('Stimulus Licks: Reward-Triggered Average Analysis - %s', position));
        f = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        % Raw data heatmap (left top)
        ax1 = nexttile([2 1]);
        imagesc(ax1, Fig8_Data.timeAxis, 1:size(Fig8_Data.allTrials, 1), Fig8_Data.allTrials(ReactionTimeSortedIndex, :));
        colormap(ax1, params.cmap);
        colorbar(ax1);
        ylabel(ax1, 'Trials');
        xlabel(ax1, 'Time to reward delivery (s)');
        xlim(ax1, [-3 3]);
        title(ax1, 'Raw Data (sorted by reaction time)');
        
        % Z-scored data heatmap (right top)
        ax2 = nexttile([2 1]);
        imagesc(ax2, Fig8_Data.timeAxis, 1:size(Fig8_Data.allTrials_zscore, 1), Fig8_Data.allTrials_zscore(ReactionTimeSortedIndex, :));
        colormap(ax2, params.cmap);
        colorbar(ax2);
        ylabel(ax2, 'Trials');
        xlabel(ax2, 'Time to reward delivery (s)');
        xlim(ax2, [-3 3]);
        title(ax2, 'Z-scored Data (sorted by reaction time)');
        
        % Raw data average trace (left bottom)
        ax3 = nexttile;
        RTA_mean = mean(Fig8_Data.allTrials, 1);
        plot(ax3, Fig8_Data.timeAxis, RTA_mean, 'LineWidth', 4, 'Color', [0.8 0.4 0.2]);
        xlim(ax3, [-3 3]);
        xlabel(ax3, 'Time to reward delivery (s)');
        ylabel(ax3, 'Entropy (Raw)');
        title(ax3, 'Average Raw Response');
        % Add vertical line at reward delivery
        xline(ax3, 0, '--', 'Reward', 'LineWidth', 2, 'Color', 'r');
        
        % Z-scored average trace (right bottom)
        ax4 = nexttile;
        RTA_mean_zscore = mean(Fig8_Data.allTrials_zscore, 1);
        plot(ax4, Fig8_Data.timeAxis, RTA_mean_zscore, 'LineWidth', 4, 'Color', [0.8 0.4 0.2]);
        xlim(ax4, [-3 3]);
        xlabel(ax4, 'Time to reward delivery (s)');
        ylabel(ax4, 'Entropy (Z-scored)');
        title(ax4, 'Average Z-scored Response');
        % Add vertical line at reward delivery
        xline(ax4, 0, '--', 'Reward', 'LineWidth', 2, 'Color', 'r');
        
        fontsize(f, 20, 'points');
        sgtitle(sprintf('Stimulus Licks: Reward-triggered averages - %s (%d hit trials)', position, size(Fig8_Data.allTrials, 1)));
        
        fprintf('Reward-triggered average for %s: %d hit trials analyzed\n', position, size(Fig8_Data.allTrials, 1));
    else
        fprintf('Warning: No data available for reward-triggered average analysis for %s\n', position);
    end
    
    %% ------------- Figure 9: Stimulus Licks - Reward Modulation Analysis -------------- %%
    fprintf('\n--- Generating Figure 9: Reward Modulation Analysis ---\n');
    
    % Create reward modulation plot
    if ~isempty(Fig9_Data.PreReward)
        % Calculate statistics
        [pReward, ~, ~] = signrank(Fig9_Data.PreReward, Fig9_Data.PostReward);
        [pReward_zscore, ~, ~] = signrank(Fig9_Data.PreReward_zscore, Fig9_Data.PostReward_zscore);
        
        ratioReward = mean(Fig9_Data.PreReward ./ Fig9_Data.PostReward, 'omitnan');
        ratioReward_zscore = mean(Fig9_Data.PreReward_zscore ./ Fig9_Data.PostReward_zscore, 'omitnan');
        
        % Calculate Hedges' g effect sizes
        % Raw data
        mean_pre_raw = mean(Fig9_Data.PreReward);
        mean_post_raw = mean(Fig9_Data.PostReward);
        std_pre_raw = std(Fig9_Data.PreReward);
        std_post_raw = std(Fig9_Data.PostReward);
        n_pre_raw = length(Fig9_Data.PreReward);
        n_post_raw = length(Fig9_Data.PostReward);
        pooled_std_raw = sqrt(((n_pre_raw - 1) * std_pre_raw^2 + (n_post_raw - 1) * std_post_raw^2) / (n_pre_raw + n_post_raw - 2));
        cohens_d_raw = (mean_post_raw - mean_pre_raw) / pooled_std_raw;
        correction_factor_raw = 1 - (3 / (4 * (n_pre_raw + n_post_raw - 2) - 1));
        hedges_g_raw = cohens_d_raw * correction_factor_raw;
        if abs(hedges_g_raw) < 0.2; effect_interp_raw = 'negligible'; elseif abs(hedges_g_raw) < 0.5; effect_interp_raw = 'small'; elseif abs(hedges_g_raw) < 0.8; effect_interp_raw = 'medium'; else; effect_interp_raw = 'large'; end
        
        % Z-scored data
        mean_pre_z = mean(Fig9_Data.PreReward_zscore);
        mean_post_z = mean(Fig9_Data.PostReward_zscore);
        std_pre_z = std(Fig9_Data.PreReward_zscore);
        std_post_z = std(Fig9_Data.PostReward_zscore);
        n_pre_z = length(Fig9_Data.PreReward_zscore);
        n_post_z = length(Fig9_Data.PostReward_zscore);
        pooled_std_z = sqrt(((n_pre_z - 1) * std_pre_z^2 + (n_post_z - 1) * std_post_z^2) / (n_pre_z + n_post_z - 2));
        cohens_d_z = (mean_post_z - mean_pre_z) / pooled_std_z;
        correction_factor_z = 1 - (3 / (4 * (n_pre_z + n_post_z - 2) - 1));
        hedges_g_z = cohens_d_z * correction_factor_z;
        if abs(hedges_g_z) < 0.2; effect_interp_z = 'negligible'; elseif abs(hedges_g_z) < 0.5; effect_interp_z = 'small'; elseif abs(hedges_g_z) < 0.8; effect_interp_z = 'medium'; else; effect_interp_z = 'large'; end
        
        figure('Color', 'w', 'Name', sprintf('Stimulus Licks: Reward Modulation Analysis - %s', position));
        tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        % Scatter plot: Pre vs Post Reward (Raw vs Z-scored)
        ax1 = nexttile;
        dscatter(Fig9_Data.PreReward, Fig9_Data.PostReward, 'MSIZE', 30);
        colormap(ax1, gray);
        axis(ax1, 'square');
        % Set equal axis limits
        allRewardData1 = [Fig9_Data.PreReward; Fig9_Data.PostReward];
        rewardLim1 = [min(allRewardData1), max(allRewardData1)];
        xlim(ax1, rewardLim1);
        ylim(ax1, rewardLim1);
        ref1 = refline(ax1, 1); 
        ref1.LineWidth = 5; 
        ref1.Color = 'r';
        set(ax1, 'FontSize', 20);
        xlabel(ax1, 'Pre-Reward Entropy (Raw)');
        ylabel(ax1, 'Post-Reward Entropy (Raw)');
        title(ax1, {'Raw Data', sprintf('ratio = %.3f, p = %.3g', ratioReward, pReward), sprintf('Hedges'' g=%.3f (%s)', hedges_g_raw, effect_interp_raw)});
        
        ax2 = nexttile;
        dscatter(Fig9_Data.PreReward_zscore, Fig9_Data.PostReward_zscore, 'MSIZE', 30);
        colormap(ax2, gray);
        axis(ax2, 'square');
        % Set equal axis limits
        allRewardData2 = [Fig9_Data.PreReward_zscore; Fig9_Data.PostReward_zscore];
        rewardLim2 = [min(allRewardData2), max(allRewardData2)];
        xlim(ax2, rewardLim2);
        ylim(ax2, rewardLim2);
        ref2 = refline(ax2, 1); 
        ref2.LineWidth = 5; 
        ref2.Color = 'r';
        set(ax2, 'FontSize', 20);
        xlabel(ax2, 'Pre-Reward Entropy (Z-scored)');
        ylabel(ax2, 'Post-Reward Entropy (Z-scored)');
        title(ax2, {'Z-scored Data', sprintf('ratio = %.3f, p = %.3g', ratioReward_zscore, pReward_zscore), sprintf('Hedges'' g=%.3f (%s)', hedges_g_z, effect_interp_z)});
        
        % Violin plot: Pre vs Post Reward (Raw vs Z-scored)
        ax3 = nexttile;
        daviolinplot({[Fig9_Data.PreReward, Fig9_Data.PostReward]}, ...
                     'violin', 'full', 'groups', ones(numel(Fig9_Data.PreReward), 1), ...
                     'colors', [0.8 0.4 0.2], 'box', 3, 'boxcolor', 'w', ...
                     'scatter', 0, 'linkline', 1, ...
                     'xtlabels', {"Pre-Reward", "Post-Reward"});
        set(ax3, 'FontSize', 20);
        ylabel(ax3, 'Entropy (Raw)');
        title(ax3, sprintf('Raw: p = %.3g', pReward));
        
        ax4 = nexttile;
        daviolinplot({[Fig9_Data.PreReward_zscore, Fig9_Data.PostReward_zscore]}, ...
                     'violin', 'full', 'groups', ones(numel(Fig9_Data.PreReward_zscore), 1), ...
                     'colors', [0.8 0.4 0.2], 'box', 3, 'boxcolor', 'w', ...
                     'scatter', 0, 'linkline', 1, ...
                     'xtlabels', {"Pre-Reward", "Post-Reward"});
        set(ax4, 'FontSize', 20);
        ylabel(ax4, 'Entropy (Z-scored)');
        title(ax4, sprintf('Z-scored: p = %.3g', pReward_zscore));
        
        sgtitle(sprintf('Stimulus Licks: Reward Modulation - %s (n = %d hit trials, reward @ 350ms)', position, length(Fig9_Data.PreReward)));
        
        fprintf('Reward modulation analysis for %s: %d hit trials analyzed\n', position, length(Fig9_Data.PreReward));
        fprintf('  Reward timing: 350ms after first lick\n');
        fprintf('  Pre-reward sampling: 300ms after first lick\n');
        fprintf('  Post-reward sampling: 600ms after first lick\n');
    else
        fprintf('Warning: No reward modulation data available for %s\n', position);
        fprintf('(This analysis requires hit trials with valid first lick timing)\n');
    end
    
    %% ------------- Figure 10: Stimulus Licks - Previous Trial History Analysis -------------- %%
    fprintf('\n--- Generating Figure 10: Previous Trial History Analysis ---\n');
    
    % Create previous trial history plot
    if ~isempty(Fig10_Data.PreviousHit) && ~isempty(Fig10_Data.PreviousMiss)
        figure('Color', 'w', 'Name', sprintf('Stimulus Licks: Previous Trial History Analysis - %s', position));
        tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        % Raw data (left)
        ax1 = nexttile;
        hold(ax1, 'on');
        aline_raw(2) = stdshade(Fig10_Data.PreviousMiss, 0.15, [0.8 0.2 0.2], Fig10_Data.timeAxis);
        aline_raw(1) = stdshade(Fig10_Data.PreviousHit, 0.15, [0.2 0.6 0.8], Fig10_Data.timeAxis);
        aline_raw(2).LineWidth = 3;
        aline_raw(1).LineWidth = 3;
        legend(ax1, aline_raw, 'Previous Hit', 'Previous Miss');
        axis(ax1, 'square');
        xlabel(ax1, 'Time in seconds');
        ylabel(ax1, 'Entropy (Raw)');
        title(ax1, 'Raw Data');
        xlim(ax1, [-2 4]);
        set(ax1, 'FontSize', 24);
        
        % Z-scored data (right)
        ax2 = nexttile;
        hold(ax2, 'on');
        aline_zscore(2) = stdshade(Fig10_Data.PreviousMiss_zscore, 0.15, [0.8 0.2 0.2], Fig10_Data.timeAxis);
        aline_zscore(1) = stdshade(Fig10_Data.PreviousHit_zscore, 0.15, [0.2 0.6 0.8], Fig10_Data.timeAxis);
        aline_zscore(2).LineWidth = 3;
        aline_zscore(1).LineWidth = 3;
        legend(ax2, aline_zscore, 'Previous Hit', 'Previous Miss');
        axis(ax2, 'square');
        xlabel(ax2, 'Time in seconds');
        ylabel(ax2, 'Entropy (Z-scored)');
        title(ax2, 'Z-scored Data');
        xlim(ax2, [-2 4]);
        set(ax2, 'FontSize', 24);
        
                 sgtitle(sprintf('Stimulus Licks: Entropy Based on Previous Trial Outcome - %s', position));
         
         fprintf('Previous trial history for %s: %d hit-following, %d miss-following trials\n', ...
                 position, size(Fig10_Data.PreviousHit, 1), size(Fig10_Data.PreviousMiss, 1));
     else
         fprintf('Warning: Insufficient data for previous trial history analysis for %s\n', position);
     end
     
     %% Collect data for summary plot
    % Collect extraneous LTA data from all positions
    if ~isempty(Fig5_Data.allTrials)
        summary_ExtraneousLTA = [summary_ExtraneousLTA; Fig5_Data.allTrials];
        if isempty(summary_timeAxis_LTA)
            summary_timeAxis_LTA = Fig5_Data.timeAxis;
        end
    end
    
    % Collect P3-specific stimulus and reward data
    if strcmp(position, 'P3')
        if ~isempty(Fig6_Data.allTrials)
            % Accumulate trials 
            if isempty(summary_P3_StimulusLTA)
                summary_P3_StimulusLTA = Fig6_Data.allTrials;
                summary_timeAxis_StimulusLTA = Fig6_Data.timeAxis;
            else
                summary_P3_StimulusLTA = [summary_P3_StimulusLTA; Fig6_Data.allTrials];
            end
            fprintf('Entropy function - P3 stimulus trials collected: %d (total: %d)\n', size(Fig6_Data.allTrials, 1), size(summary_P3_StimulusLTA, 1));
        end
        if ~isempty(Fig8_Data.allTrials)
            summary_P3_RewardRTA = Fig8_Data.allTrials;
            if isempty(summary_timeAxis_RTA)
                summary_timeAxis_RTA = Fig8_Data.timeAxis;
            end
        end
    end
     
     %% Summary output for this position
    fprintf('\nlickCorrelations_Entropy completed for position %s!\n', position);
    fprintf('Processed %d conditions: %s\n', length(conditionsWithLicks), strjoin(conditionsWithLicks, ', '));
    fprintf('Total extraneous licks analyzed: %d\n', length(Fig1_Data.Extraneous.Pre));
    fprintf('Total stimulus licks analyzed: %d\n', length(Fig1_Data.Stimulus.Pre));
    fprintf('Total reward timepoints analyzed: %d\n', length(Fig1_Data.Reward.Pre));
end

fprintf('\n\n======== ALL POSITIONS COMPLETED ========\n');
fprintf('Analyzed positions: %s\n', strjoin(positions, ', '));

%% ======================== SUMMARY PLOT ========================
fprintf('\n--- Generating Summary Plot: Cross-Position Lick-Triggered Entropy ---\n');

% Create summary plot only if we have data
if ~isempty(summary_ExtraneousLTA) || ~isempty(summary_P3_StimulusLTA) || ~isempty(summary_P3_RewardRTA)
    
    figure('Color', 'w', 'Name', 'Summary: Cross-Position Lick-Triggered Entropy');
    tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Tile 1: Extraneous licks (all positions combined)
    ax1 = nexttile;
    if ~isempty(summary_ExtraneousLTA) && ~isempty(summary_timeAxis_LTA)
        % Calculate mean and SEM
        mean_trace = mean(summary_ExtraneousLTA, 1, 'omitnan');
        std_trace = std(summary_ExtraneousLTA, 0, 1, 'omitnan');
        
        % Ensure all vectors are row vectors for consistent dimensions
        timeAxis = summary_timeAxis_LTA(:)'; % Force row vector
        meanTrace = mean_trace(:)'; % Force row vector
        stdTrace = std_trace(:)'; % Force row vector
        
        % Plot with shaded error bars
        hold(ax1, 'on');
        fill(ax1, [timeAxis, fliplr(timeAxis)], ...
             [meanTrace + stdTrace, fliplr(meanTrace - stdTrace)], ...
             [0 0 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        plot(ax1, timeAxis, meanTrace, 'LineWidth', 3, 'Color', [0 0 0]);
        
        xlim(ax1, [-3 3]);
        ylim(ax1,[0.75 0.95]);
        xlabel(ax1, 'Time to extraneous lick (s)');
        ylabel(ax1, 'Entropy (Raw)');
        title(ax1, sprintf('Extraneous Licks\n(All Positions, n=%d)', size(summary_ExtraneousLTA, 1)));
        grid(ax1, 'on');
        set(ax1, 'FontSize', 16);
        
        % Add vertical line at lick time
        xline(ax1, 0, '--', 'Lick', 'LineWidth', 2, 'Color', 'r');
    else
        text(ax1, 0.5, 0.5, 'No extraneous lick data', 'HorizontalAlignment', 'center', 'FontSize', 14);
        title(ax1, 'Extraneous Licks (All Positions)');
        set(ax1, 'XTick', [], 'YTick', []);
    end
    
    % Tile 2: Combined Stimulus + Reward licks (P3 only)
    ax2 = nexttile;
    if ~isempty(summary_P3_StimulusLTA) && ~isempty(summary_timeAxis_StimulusLTA)
        % Calculate mean and std for stimulus data
        mean_trace = mean(summary_P3_StimulusLTA, 1, 'omitnan');
        std_trace = std(summary_P3_StimulusLTA, 0, 1, 'omitnan');
        
        % Ensure all vectors are row vectors for consistent dimensions
        timeAxis = summary_timeAxis_StimulusLTA(:)'; % Force row vector - use correct time axis
        meanTrace = mean_trace(:)'; % Force row vector
        stdTrace = std_trace(:)'; % Force row vector
        
        % Check if vectors have same length
        if length(timeAxis) ~= length(meanTrace) || length(timeAxis) ~= length(stdTrace)
            % Skip fill and just plot the line if possible
            if length(timeAxis) == length(meanTrace)
                plot(ax2, timeAxis, meanTrace, 'LineWidth', 3, 'Color', [0.2 0.6 0.8]);
            else
                text(ax2, 0, 0, 'Dimension mismatch in stimulus data', 'HorizontalAlignment', 'center');
            end
        else
            % Plot with shaded error bars
            hold(ax2, 'on');
            fill(ax2, [timeAxis, fliplr(timeAxis)], ...
                 [meanTrace + stdTrace, fliplr(meanTrace - stdTrace)], ...
                 [0.2 0.6 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            plot(ax2, timeAxis, meanTrace, 'LineWidth', 3, 'Color', [0.2 0.6 0.8]);
        end
        
        xlim(ax2, [-4 4]); % Extended window by 1 second
        ylim(ax1,[0.75 0.95]);
        xlabel(ax2, 'Time to first lick (s)');
        ylabel(ax2, 'Entropy (Raw)');
        title(ax2, sprintf('Stimulus Licks + Reward\n(P3 Only, n=%d)', size(summary_P3_StimulusLTA, 1)));
        grid(ax2, 'on');
        set(ax2, 'FontSize', 16);
        
        % Add vertical line at lick time
        xline(ax2, 0, '--', 'Lick', 'LineWidth', 2, 'Color', 'r');
        % Add vertical line at reward time (350ms = 0.35s after lick)
        xline(ax2, 0.35, '--', 'Reward', 'LineWidth', 2, 'Color', [0.8 0.4 0.2]);
    else
        text(ax2, 0.5, 0.5, 'No P3 stimulus lick data', 'HorizontalAlignment', 'center', 'FontSize', 14);
        title(ax2, 'Stimulus Licks + Reward (P3 Only)');
        set(ax2, 'XTick', [], 'YTick', []);
    end
    
    sgtitle('Summary: Lick-Triggered Entropy Across Conditions');
    
    fprintf('Summary plot generated successfully!\n');
    fprintf('  Extraneous licks (all positions): %d trials\n', size(summary_ExtraneousLTA, 1));
    fprintf('  P3 stimulus licks + reward (combined): %d trials\n', size(summary_P3_StimulusLTA, 1));
else
    fprintf('No data available for summary plot generation\n');
end

% Note about experimental conditions
fprintf('\nExperimental setup note:\n');
fprintf('- Conditions with spout: Expert, Beginner (full lick-triggered analyses performed)\n');
fprintf('- Conditions without spout: Naive, NoSpout (basic analyses only - no lick data available)\n');

% Summary of skip information
if ~isempty(fieldnames(Skip))
    fprintf('\nSummary of skipped recordings:\n');
    for skip_field = fieldnames(Skip)'
        fprintf('  %s: %d recording(s) skipped (%s)\n', skip_field{1}, length(Skip.(skip_field{1})), mat2str(Skip.(skip_field{1})));
    end
else
    fprintf('No recordings were skipped in this analysis.\n');
end

%% ======================== CONTROL FIGURE ========================
fprintf('\n--- Generating Control Figure: Pre-stimulus vs Stimulus Onset (Frame 80 vs 83) ---\n');

% Extract entropy data at specific frames for P1 and P3 positions
controlFig_Data = struct();
controlFig_Data.P1 = struct('PreStim', [], 'StimOnset', []);
controlFig_Data.P3 = struct('PreStim', [], 'StimOnset', []);

% Process each condition to extract frame 80 and 83 data
for c_idx = 1:length(conditionsWithLicks)
    condName = conditionsWithLicks{c_idx};
    fprintf('Extracting control data for condition: %s\n', condName);
    
    % Get recordings to skip for this condition
    recordingsToSkip = [];
    if isfield(Skip, condName); recordingsToSkip = Skip.(condName); end
    
    % Check if condition exists
    if ~isfield(Entropy, condName) || ~isfield(Performance, condName) || ~isfield(Rec, condName)
        continue;
    end
    
    condEntropyData = Entropy.(condName);
    condPerformanceData = Performance.(condName);
    
    % Determine number of recordings
    if isfield(Rec.(condName), 'AnimalID')
        numRecs = length(Rec.(condName).AnimalID);
    else
        numRecs = size(Rec.(condName), 1);
    end
    fprintf('  Total recordings for %s: %d\n', condName, numRecs);
    
    % Extract data for each recording
    for rec = 1:numRecs
        if rec > length(condPerformanceData) || ismember(rec, recordingsToSkip)
            fprintf('  Skipping recording %d (performance data or skip list)\n', rec);
            continue;
        end
        
        % Check if entropy data exists for this recording
        if ~isfield(condEntropyData, 'AllNeurons') || ~isfield(condEntropyData.AllNeurons, 'RecordingEntropy_Raw')
            fprintf('  Recording %d: Missing AllNeurons.RecordingEntropy_Raw field\n', rec);
            continue;
        end
        
        % Extract entropy data from cell array
        if iscell(condEntropyData.AllNeurons.RecordingEntropy_Raw) && rec <= size(condEntropyData.AllNeurons.RecordingEntropy_Raw, 1)
            if ~isempty(condEntropyData.AllNeurons.RecordingEntropy_Raw{rec, 1})
                entropyData = condEntropyData.AllNeurons.RecordingEntropy_Raw{rec, 1}; % [trials x frames]
                fprintf('  Recording %d: Found entropy data [%d trials x %d frames]\n', rec, size(entropyData, 1), size(entropyData, 2));
            else
                fprintf('  Recording %d: Empty entropy data cell\n', rec);
                continue;
            end
        else
            fprintf('  Recording %d: Invalid cell array or index out of bounds\n', rec);
            continue;
        end
        
        % Get P1 and P3 trials
        for pos_name = {'P1', 'P3'}
            position = pos_name{1};
            
            % Get hit and miss trials for this position
            hitField = sprintf('hit%s', position);
            missField = sprintf('miss%s', position);
            
            if isfield(condPerformanceData(rec), hitField) && isfield(condPerformanceData(rec), missField)
                hitTrials = condPerformanceData(rec).(hitField);
                missTrials = condPerformanceData(rec).(missField);
                validTrials = [hitTrials, missTrials];
                fprintf('    %s - Hit: %d, Miss: %d, Total: %d trials\n', position, length(hitTrials), length(missTrials), length(validTrials));
                
                % Extract frame 80 and 83 data for valid trials
                validCount = 0;
                for trial = validTrials
                    if trial <= size(entropyData, 1) && size(entropyData, 2) >= 83
                        % Frame 80 (pre-stimulus) and 83 (stimulus onset)
                        preStimEntropy = entropyData(trial, 80);
                        stimOnsetEntropy = entropyData(trial, 83);
                        
                        % Store data
                        controlFig_Data.(position).PreStim(end+1) = preStimEntropy;
                        controlFig_Data.(position).StimOnset(end+1) = stimOnsetEntropy;
                        validCount = validCount + 1;
                    end
                end
                fprintf('    %s - Valid data points extracted: %d\n', position, validCount);
            else
                fprintf('    %s - Missing hit/miss fields (%s, %s)\n', position, hitField, missField);
            end
        end
    end
end

% Generate control figure with tiledlayout
if ~isempty(controlFig_Data.P1.PreStim) || ~isempty(controlFig_Data.P3.PreStim)
    figure('Color', 'w', 'Name', 'Control: Pre-stimulus vs Stimulus Onset Entropy (Frame 80 vs 83)');
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Define colormap and marker size for scatter plots
    scatterColormap = gray;
    scatterMarkerSize = 30;
    
    % Left tile: P1 data
    nexttile;
    if ~isempty(controlFig_Data.P1.PreStim)
        % Create scatter plot - ensure column vectors
        preStimP1 = controlFig_Data.P1.PreStim(:);
        stimOnsetP1 = controlFig_Data.P1.StimOnset(:);
        dscatter(preStimP1, stimOnsetP1, 'MSIZE', scatterMarkerSize);
        colormap(scatterColormap);
        axis square;
        
        % Set equal axis limits
        allData = [preStimP1; stimOnsetP1];
        dataLim = [min(allData), max(allData)];
        xlim(dataLim);
        ylim(dataLim);
        setNiceAxisTicks(gca, dataLim);
        
        % Add unity line
        ref = refline(1);
        ref.LineWidth = 5;
        ref.Color = 'r';
        
        % Add fitted trend line
        fit = lsline(gca);
        fit.LineStyle = '--';
        fit.LineWidth = 2;
        fit.Color = [0.5 0.5 0.5];
        
        % Statistical tests
        [pP1, ~, ~] = signrank(preStimP1, stimOnsetP1);
        [rhoP1, pRhoP1] = corr(preStimP1, stimOnsetP1, 'Type', 'Spearman', 'Rows', 'complete');
        
        % Labels and title
        ax = gca;
        ax.FontSize = 16;
        xlabel('Frame 80 (Pre-stimulus) Entropy');
        ylabel('Frame 83 (Stimulus Onset) Entropy');
        title({sprintf('P1 Position Control'), ...
               sprintf('p_{signrank} = %.3g,  ρ = %.3f,  p_{rho} = %.3g', pP1, rhoP1, pRhoP1), ...
               sprintf('n = %d trials', length(preStimP1))});
        
        fprintf('P1 Control: n=%d, p_signrank=%.3g, rho=%.3f, p_rho=%.3g\n', ...
                length(preStimP1), pP1, rhoP1, pRhoP1);
    else
        % No data message
        text(0.5, 0.5, 'No P1 data available', 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', 'FontSize', 14);
        xlim([0 1]); ylim([0 1]);
        title('P1 Position Control (No Data)');
    end
    
    % Right tile: P3 data
    nexttile;
    if ~isempty(controlFig_Data.P3.PreStim)
        % Create scatter plot - ensure column vectors
        preStimP3 = controlFig_Data.P3.PreStim(:);
        stimOnsetP3 = controlFig_Data.P3.StimOnset(:);
        dscatter(preStimP3, stimOnsetP3, 'MSIZE', scatterMarkerSize);
        colormap(scatterColormap);
        axis square;
        
        % Set equal axis limits
        allData = [preStimP3; stimOnsetP3];
        dataLim = [min(allData), max(allData)];
        xlim(dataLim);
        ylim(dataLim);
        setNiceAxisTicks(gca, dataLim);
        
        % Add unity line
        ref = refline(1);
        ref.LineWidth = 5;
        ref.Color = 'r';
        
        % Add fitted trend line
        fit = lsline(gca);
        fit.LineStyle = '--';
        fit.LineWidth = 2;
        fit.Color = [0.5 0.5 0.5];
        
        % Statistical tests
        [pP3, ~, ~] = signrank(preStimP3, stimOnsetP3);
        [rhoP3, pRhoP3] = corr(preStimP3, stimOnsetP3, 'Type', 'Spearman', 'Rows', 'complete');
        
        % Labels and title
        ax = gca;
        ax.FontSize = 16;
        xlabel('Frame 80 (Pre-stimulus) Entropy');
        ylabel('Frame 83 (Stimulus Onset) Entropy');
        title({sprintf('P3 Position Control'), ...
               sprintf('p_{signrank} = %.3g,  ρ = %.3f,  p_{rho} = %.3g', pP3, rhoP3, pRhoP3), ...
               sprintf('n = %d trials', length(preStimP3))});
        
        fprintf('P3 Control: n=%d, p_signrank=%.3g, rho=%.3f, p_rho=%.3g\n', ...
                length(preStimP3), pP3, rhoP3, pRhoP3);
    else
        % No data message
        text(0.5, 0.5, 'No P3 data available', 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', 'FontSize', 14);
        xlim([0 1]); ylim([0 1]);
        title('P3 Position Control (No Data)');
    end
    
    fprintf('Control figure generated: Pre-stimulus (Frame 80) vs Stimulus Onset (Frame 83) Entropy\n');
else
    fprintf('No control data available for control figure generation\n');
end

end 

%% ======================== HELPER FUNCTIONS ========================

function [validLickTrials] = extractValidLickTrials_Entropy(conditionsWithLicks, position, Entropy, Performance, Rec, Skip, timeReference, LTAwindow, RTAwindow)
    % EXTRACTVALIDLICKTRIALS_ENTROPY - Unified trial selection logic for all entropy analyses
    % This function creates master lists of valid trials ensuring consistency across all figures
    
    validLickTrials = struct();
    
    fprintf('=== Unified Trial Selection for %s ===\n', position);
    
    for c_idx = 1:length(conditionsWithLicks)
        condName = conditionsWithLicks{c_idx};
        
        % Get recordings to skip for this condition
        recordingsToSkip = [];
        if isfield(Skip, condName); recordingsToSkip = Skip.(condName); end
        
        if ~isfield(Entropy, condName) || ~isfield(Performance, condName)
            continue;
        end
        
        condEntropy = Entropy.(condName);
        condPerformance = Performance.(condName);
        
        % Determine number of recordings
        if isfield(Rec.(condName), 'AnimalID')
            numRecs = length(Rec.(condName).AnimalID);
        else
            numRecs = size(Rec.(condName), 1);
        end
        
        % Initialize condition structure
        validLickTrials.(condName) = struct();
        
        % Process each recording
        for r = 1:numRecs
            if ismember(r, recordingsToSkip); continue; end
            
            if r > length(condPerformance)
                continue;
            end
            
            % Skip if no LicksInFrame data
            if ~isfield(condPerformance(r), 'LicksInFrame')
                continue;
            end
            
            % Check if entropy data exists for this recording
            if ~isfield(condEntropy, 'AllNeurons') || ~isfield(condEntropy.AllNeurons, 'RecordingEntropy_Raw')
                continue;
            end
            
            if iscell(condEntropy.AllNeurons.RecordingEntropy_Raw) && r <= size(condEntropy.AllNeurons.RecordingEntropy_Raw, 1)
                if ~isempty(condEntropy.AllNeurons.RecordingEntropy_Raw{r, 1})
                    recEntropyData = condEntropy.AllNeurons.RecordingEntropy_Raw{r, 1};
                else
                    continue;
                end
            else
                continue;
            end
            
            % Get position-specific trials
            hitTrials = [];
            missTrials = [];
            
            hitFieldName = sprintf('hit%s', position);
            missFieldName = sprintf('miss%s', position);
            
            if isfield(condPerformance(r), hitFieldName)
                hitTrials = condPerformance(r).(hitFieldName);
            elseif isfield(condPerformance(r), 'hit')
                hitTrials = condPerformance(r).hit;
            end
            
            if isfield(condPerformance(r), missFieldName)
                missTrials = condPerformance(r).(missFieldName);
            elseif isfield(condPerformance(r), 'miss')
                missTrials = condPerformance(r).miss;
            end
            
            % Combine hit and miss trials for this position
            allPositionTrials = [hitTrials, missTrials];
            
            if isempty(allPositionTrials)
                continue; % Skip if no trials for this position
            end
            
            % Extract reaction times and apply unified filtering
            validAllTrials = [];
            validHitTrials = [];
            allReactionTimes = [];
            hitReactionTimes = [];
            
            % Process all position trials
            for i = 1:length(allPositionTrials)
                t = allPositionTrials(i);
                if t <= size(condPerformance(r).LicksInFrame, 2) && t <= size(recEntropyData, 1)
                    firstLick = find(condPerformance(r).LicksInFrame(:, t) > 0.2, 1);
                    if ~isempty(firstLick)
                        reactionTime = condPerformance(r).LicksInFrame(firstLick, t);
                        
                        % Apply reaction time filtering (< 2 seconds)
                        if reactionTime < 2
                            % Apply unified bounds checking
                            lickFrame = find(reactionTime < timeReference, 1);
                            rewardFrame = lickFrame + 4; % 350ms after lick
                            
                            % Use same criteria as activity data function - only check stimulus window for trial inclusion
                            if lickFrame > LTAwindow && (lickFrame + LTAwindow) <= size(recEntropyData, 2)
                                
                                % This trial passes stimulus criteria
                                validAllTrials(end+1) = t; %#ok<AGROW>
                                allReactionTimes(end+1) = reactionTime; %#ok<AGROW>
                                
                                % Check if this is also a hit trial 
                                if ismember(t, hitTrials)
                                    validHitTrials(end+1) = t; %#ok<AGROW>
                                    hitReactionTimes(end+1) = reactionTime; %#ok<AGROW>
                                end
                            end
                        end
                    end
                end
            end
            
            % Store results for this recording
            if ~isempty(validAllTrials)
                validLickTrials.(condName).(sprintf('rec%d', r)).allTrials = validAllTrials;
                validLickTrials.(condName).(sprintf('rec%d', r)).allReactionTimes = allReactionTimes;
                validLickTrials.(condName).(sprintf('rec%d', r)).hitTrials = validHitTrials;
                validLickTrials.(condName).(sprintf('rec%d', r)).hitReactionTimes = hitReactionTimes;
                validLickTrials.(condName).(sprintf('rec%d', r)).recEntropyData = recEntropyData;
                
                fprintf('  %s rec %d: %d total valid licks, %d hit licks\n', ...
                        condName, r, length(validAllTrials), length(validHitTrials));
            end
        end
    end
end

function [Fig1_Data] = extractPrePostEntropy(validLickTrials, Entropy, Performance, timeReference)
    % EXTRACTPREPOSTENTROPY - Extract pre/post entropy data for Figures 1-4
    % Extracts entropy before and after lick events for different lick types
    
    fprintf('=== Extracting Pre/Post Entropy Data ===\n');
    
    % Initialize output structure
    Fig1_Data = struct();
    Fig1_Data.Extraneous = struct('Pre', [], 'Post', []);
    Fig1_Data.Stimulus = struct('Pre', [], 'Post', []);
    Fig1_Data.Reward = struct('Pre', [], 'Post', []);
    
    condNames = fieldnames(validLickTrials);
    
    for c_idx = 1:length(condNames)
        condName = condNames{c_idx};
        condEntropy = Entropy.(condName);
        condPerformance = Performance.(condName);
        
        recNames = fieldnames(validLickTrials.(condName));
        
        for r_idx = 1:length(recNames)
            recName = recNames{r_idx};
            recNum = str2double(recName(4:end)); % Extract number from 'rec#'
            
            if recNum > length(condPerformance)
                continue;
            end
            
            allTrials = validLickTrials.(condName).(recName).allTrials;
            allReactionTimes = validLickTrials.(condName).(recName).allReactionTimes;
            hitTrials = validLickTrials.(condName).(recName).hitTrials;
            hitReactionTimes = validLickTrials.(condName).(recName).hitReactionTimes;
            recEntropyData = validLickTrials.(condName).(recName).recEntropyData;
            
            % Extract stimulus lick data (all valid trials)
            for i = 1:length(allTrials)
                t = allTrials(i);
                reactionTime = allReactionTimes(i);
                
                lickFrame = find(reactionTime < timeReference, 1);
                if ~isempty(lickFrame) && lickFrame > 2 && lickFrame < (size(recEntropyData, 2) - 2)
                    preEntropy = recEntropyData(t, lickFrame-2);
                    postEntropy = recEntropyData(t, lickFrame+2);
                    
                    Fig1_Data.Stimulus.Pre = [Fig1_Data.Stimulus.Pre; preEntropy]; %#ok<AGROW>
                    Fig1_Data.Stimulus.Post = [Fig1_Data.Stimulus.Post; postEntropy]; %#ok<AGROW>
                end
            end
            
            % Extract reward timepoint data (all valid trials)
            for i = 1:length(allTrials)
                t = allTrials(i);
                reactionTime = allReactionTimes(i);
                
                lickFrame = find(reactionTime < timeReference, 1);
                if ~isempty(lickFrame) && lickFrame + 6 <= size(recEntropyData, 2)
                    preRewardEntropy = recEntropyData(t, lickFrame+1); % Just before reward
                    postRewardEntropy = recEntropyData(t, lickFrame+5); % Just after reward
                    
                    Fig1_Data.Reward.Pre = [Fig1_Data.Reward.Pre; preRewardEntropy]; %#ok<AGROW>
                    Fig1_Data.Reward.Post = [Fig1_Data.Reward.Post; postRewardEntropy]; %#ok<AGROW>
                end
            end
        end
    end
    
    % Extract extraneous lick data (separate logic for pre-stimulus licks)
    % This uses all trials and looks for licks between -7.9 and -0.2 seconds
    for c_idx = 1:length(condNames)
        condName = condNames{c_idx};
        condEntropy = Entropy.(condName);
        condPerformance = Performance.(condName);
        
        recNames = fieldnames(validLickTrials.(condName));
        
        for r_idx = 1:length(recNames)
            recName = recNames{r_idx};
            recNum = str2double(recName(4:end));
            
            if recNum > length(condPerformance)
                continue;
            end
            
            recEntropyData = validLickTrials.(condName).(recName).recEntropyData;
            
            % Process each trial for extraneous licks
            for t = 1:size(condPerformance(recNum).LicksInFrame, 2)
                if t > size(recEntropyData, 1)
                    continue;
                end
                
                % Find extraneous licks (licks between -7.9 and -0.2 seconds)
                extraneousLicks = condPerformance(recNum).LicksInFrame( ...
                    condPerformance(recNum).LicksInFrame(:, t) < -0.2 & ...
                    condPerformance(recNum).LicksInFrame(:, t) > -7.9, t);
                
                % Process each extraneous lick
                for e = 1:numel(extraneousLicks)
                    lickTime = extraneousLicks(e);
                    lickFrame = find(lickTime < timeReference, 1);
                    
                    % Check if we have enough frames before and after
                    if ~isempty(lickFrame) && lickFrame > 2 && lickFrame < (size(recEntropyData, 2) - 2)
                        preEntropy = recEntropyData(t, lickFrame-2);
                        postEntropy = recEntropyData(t, lickFrame+2);
                        
                        Fig1_Data.Extraneous.Pre = [Fig1_Data.Extraneous.Pre; preEntropy]; %#ok<AGROW>
                        Fig1_Data.Extraneous.Post = [Fig1_Data.Extraneous.Post; postEntropy]; %#ok<AGROW>
                    end
                end
            end
        end
    end
    
    % Compute z-scored versions
    % Extraneous licks
    if ~isempty(Fig1_Data.Extraneous.Pre)
        allExtraneousData = [Fig1_Data.Extraneous.Pre; Fig1_Data.Extraneous.Post];
        meanExtraneous = mean(allExtraneousData);
        stdExtraneous = std(allExtraneousData);
        stdExtraneous(stdExtraneous == 0) = 1;
        
        Fig1_Data.Extraneous.Pre_zscore = (Fig1_Data.Extraneous.Pre - meanExtraneous) ./ stdExtraneous;
        Fig1_Data.Extraneous.Post_zscore = (Fig1_Data.Extraneous.Post - meanExtraneous) ./ stdExtraneous;
    end
    
    % Stimulus licks
    if ~isempty(Fig1_Data.Stimulus.Pre)
        allStimulusData = [Fig1_Data.Stimulus.Pre; Fig1_Data.Stimulus.Post];
        meanStimulus = mean(allStimulusData);
        stdStimulus = std(allStimulusData);
        stdStimulus(stdStimulus == 0) = 1;
        
        Fig1_Data.Stimulus.Pre_zscore = (Fig1_Data.Stimulus.Pre - meanStimulus) ./ stdStimulus;
        Fig1_Data.Stimulus.Post_zscore = (Fig1_Data.Stimulus.Post - meanStimulus) ./ stdStimulus;
    end
    
    % Reward timepoints
    if ~isempty(Fig1_Data.Reward.Pre)
        allRewardData = [Fig1_Data.Reward.Pre; Fig1_Data.Reward.Post];
        meanReward = mean(allRewardData);
        stdReward = std(allRewardData);
        stdReward(stdReward == 0) = 1;
        
        Fig1_Data.Reward.Pre_zscore = (Fig1_Data.Reward.Pre - meanReward) ./ stdReward;
        Fig1_Data.Reward.Post_zscore = (Fig1_Data.Reward.Post - meanReward) ./ stdReward;
    end
    
    fprintf('  Extraneous licks: %d\n', size(Fig1_Data.Extraneous.Pre, 1));
    fprintf('  Stimulus licks: %d\n', size(Fig1_Data.Stimulus.Pre, 1));
    fprintf('  Reward timepoints: %d\n', size(Fig1_Data.Reward.Pre, 1));
end 

function [Fig3_Data] = calculateEntropyRatios(Fig1_Data)
    % CALCULATEENTROPYRATOS - Calculate Post/Pre ratios for Figure 3
    
    fprintf('=== Calculating Entropy Ratios ===\n');
    
    Fig3_Data = struct();
    
    % Extraneous licks
    if ~isempty(Fig1_Data.Extraneous.Pre)
        Fig3_Data.Extraneous.Ratio = Fig1_Data.Extraneous.Post ./ Fig1_Data.Extraneous.Pre;
        Fig3_Data.Extraneous.Ratio_zscore = Fig1_Data.Extraneous.Post_zscore ./ Fig1_Data.Extraneous.Pre_zscore;
        Fig3_Data.Extraneous.Ratio_zscore = Fig3_Data.Extraneous.Ratio_zscore(isfinite(Fig3_Data.Extraneous.Ratio_zscore));
    end
    
    % Stimulus licks
    if ~isempty(Fig1_Data.Stimulus.Pre)
        Fig3_Data.Stimulus.Ratio = Fig1_Data.Stimulus.Post ./ Fig1_Data.Stimulus.Pre;
        Fig3_Data.Stimulus.Ratio_zscore = Fig1_Data.Stimulus.Post_zscore ./ Fig1_Data.Stimulus.Pre_zscore;
        Fig3_Data.Stimulus.Ratio_zscore = Fig3_Data.Stimulus.Ratio_zscore(isfinite(Fig3_Data.Stimulus.Ratio_zscore));
    end
    
    % Reward timepoints
    if ~isempty(Fig1_Data.Reward.Pre)
        Fig3_Data.Reward.Ratio = Fig1_Data.Reward.Post ./ Fig1_Data.Reward.Pre;
        Fig3_Data.Reward.Ratio_zscore = Fig1_Data.Reward.Post_zscore ./ Fig1_Data.Reward.Pre_zscore;
        Fig3_Data.Reward.Ratio_zscore = Fig3_Data.Reward.Ratio_zscore(isfinite(Fig3_Data.Reward.Ratio_zscore));
    end
end

function [Fig4_Data] = calculateEntropyChanges(Fig1_Data)
    % CALCULATEENTROPYCHANGES - Calculate Post-Pre changes for Figure 4
    
    fprintf('=== Calculating Entropy Changes ===\n');
    
    Fig4_Data = struct();
    
    % Extraneous licks
    if ~isempty(Fig1_Data.Extraneous.Pre)
        Fig4_Data.Extraneous.Change = Fig1_Data.Extraneous.Post - Fig1_Data.Extraneous.Pre;
        Fig4_Data.Extraneous.Change_zscore = Fig1_Data.Extraneous.Post_zscore - Fig1_Data.Extraneous.Pre_zscore;
    end
    
    % Stimulus licks
    if ~isempty(Fig1_Data.Stimulus.Pre)
        Fig4_Data.Stimulus.Change = Fig1_Data.Stimulus.Post - Fig1_Data.Stimulus.Pre;
        Fig4_Data.Stimulus.Change_zscore = Fig1_Data.Stimulus.Post_zscore - Fig1_Data.Stimulus.Pre_zscore;
    end
    
    % Reward timepoints
    if ~isempty(Fig1_Data.Reward.Pre)
        Fig4_Data.Reward.Change = Fig1_Data.Reward.Post - Fig1_Data.Reward.Pre;
        Fig4_Data.Reward.Change_zscore = Fig1_Data.Reward.Post_zscore - Fig1_Data.Reward.Pre_zscore;
    end
end

function [Fig5_Data] = extractExtraneousLTA_Entropy(conditionsWithLicks, Entropy, Performance, timeReference, LTAwindow)
    % EXTRACTEXTRANEOULA_ENTROPY - Extract extraneous lick-triggered averages for Figure 5
    
    fprintf('=== Extracting Extraneous LTA Data ===\n');
    
    Fig5_Data = struct();
    Fig5_Data.allTrials = [];
    Fig5_Data.lickTimes = [];
    
    for c_idx = 1:length(conditionsWithLicks)
        condName = conditionsWithLicks{c_idx};
        condEntropy = Entropy.(condName);
        condPerformance = Performance.(condName);
        
        % Process each recording
        for r = 1:length(condPerformance)
            if ~isfield(condPerformance(r), 'LicksInFrame')
                continue;
            end
            
            % Check if entropy data exists for this recording
            if ~isfield(condEntropy, 'AllNeurons') || ~isfield(condEntropy.AllNeurons, 'RecordingEntropy_Raw')
                continue;
            end
            
            if iscell(condEntropy.AllNeurons.RecordingEntropy_Raw) && r <= size(condEntropy.AllNeurons.RecordingEntropy_Raw, 1)
                if ~isempty(condEntropy.AllNeurons.RecordingEntropy_Raw{r, 1})
                    recEntropyData = condEntropy.AllNeurons.RecordingEntropy_Raw{r, 1};
                else
                    continue;
                end
            else
                continue;
            end
            
            % Process each trial for extraneous licks
            for t = 1:size(condPerformance(r).LicksInFrame, 2)
                if t > size(recEntropyData, 1)
                    continue;
                end
                
                % Find extraneous licks (licks between -7.9 and -0.2 seconds)
                extraneousLicks = condPerformance(r).LicksInFrame( condPerformance(r).LicksInFrame(:, t) < -3 );
                
                % Process each extraneous lick
                for e = 1:numel(extraneousLicks)
                    lickTime = extraneousLicks(e);
                    lickFrame = find(lickTime < timeReference, 1);
                    
                    if isempty(lickFrame)
                        continue;
                    end
                    
                    % Calculate window boundaries with smart padding
                    leftBound = max(1, lickFrame - LTAwindow);
                    rightBound = min(size(recEntropyData, 2), lickFrame + LTAwindow);
                    
                    % Don't let the window extend past stimulus onset
                    stimulusFrame = find(0 < timeReference, 1);
                    if ~isempty(stimulusFrame)
                        rightBound = min(rightBound, stimulusFrame - 1);
                    end
                    
                    % Check if we have enough data
                    if rightBound <= leftBound
                        continue;
                    end
                    
                    % Extract the available entropy data
                    availableEntropy = recEntropyData(t, leftBound:rightBound);
                    
                    % Create full window with NaN padding
                    fullWindow = NaN(1, 2*LTAwindow + 1);
                    
                    % Calculate where to place the available data
                    windowCenter = LTAwindow + 1;
                    dataStartInWindow = windowCenter - (lickFrame - leftBound);
                    dataEndInWindow = dataStartInWindow + length(availableEntropy) - 1;
                    
                    % Place available data in the appropriate position
                    if dataStartInWindow >= 1 && dataEndInWindow <= length(fullWindow)
                        fullWindow(dataStartInWindow:dataEndInWindow) = availableEntropy;
                        Fig5_Data.allTrials = [Fig5_Data.allTrials; fullWindow]; %#ok<AGROW>
                        Fig5_Data.lickTimes = [Fig5_Data.lickTimes, lickTime]; %#ok<AGROW>
                    end
                end
            end
        end
    end
    
    % Compute z-scored data
    if ~isempty(Fig5_Data.allTrials)
        validData = Fig5_Data.allTrials(~isnan(Fig5_Data.allTrials));
        if ~isempty(validData)
            meanExtraneous = mean(validData);
            stdExtraneous = std(validData);
            if stdExtraneous == 0; stdExtraneous = 1; end
            Fig5_Data.allTrials_zscore = (Fig5_Data.allTrials - meanExtraneous) ./ stdExtraneous;
        else
            Fig5_Data.allTrials_zscore = Fig5_Data.allTrials;
        end
        
        % Create time axis
        ifi = 18.5/(2*LTAwindow + 1);
        Fig5_Data.timeAxis = linspace(-LTAwindow*ifi, LTAwindow*ifi, size(Fig5_Data.allTrials, 2));
    end
    
    fprintf('  Extraneous LTA trials: %d\n', size(Fig5_Data.allTrials, 1));
end

function [Fig6_Data] = extractStimulusLTA_Entropy(validLickTrials, timeReference, LTAwindow)
    % EXTRACTSTIMULUSLTA_ENTROPY - Extract stimulus lick-triggered averages for Figure 6
    
    fprintf('=== Extracting Stimulus LTA Data ===\n');
    
    Fig6_Data = struct();
    Fig6_Data.allTrials = [];
    Fig6_Data.reactionTimes = [];
    
    condNames = fieldnames(validLickTrials);
    
    for c_idx = 1:length(condNames)
        condName = condNames{c_idx};
        
        recNames = fieldnames(validLickTrials.(condName));
        
        for r_idx = 1:length(recNames)
            recName = recNames{r_idx};
            
            allTrials = validLickTrials.(condName).(recName).allTrials;
            allReactionTimes = validLickTrials.(condName).(recName).allReactionTimes;
            recEntropyData = validLickTrials.(condName).(recName).recEntropyData;
            
            % Extract LTA for each valid trial
            for i = 1:length(allTrials)
                t = allTrials(i);
                reactionTime = allReactionTimes(i);
                
                lickFrame = find(reactionTime < timeReference, 1);
                if ~isempty(lickFrame)
                    % Extract entropy around lick
                    lickTriggeredEntropy = recEntropyData(t, lickFrame-LTAwindow:lickFrame+LTAwindow);
                    Fig6_Data.allTrials = [Fig6_Data.allTrials; lickTriggeredEntropy]; %#ok<AGROW>
                    Fig6_Data.reactionTimes = [Fig6_Data.reactionTimes, reactionTime]; %#ok<AGROW>
                end
            end
        end
    end
    
    % Compute z-scored data
    if ~isempty(Fig6_Data.allTrials)
        meanLTA = mean(Fig6_Data.allTrials(:));
        stdLTA = std(Fig6_Data.allTrials(:));
        if stdLTA == 0; stdLTA = 1; end
        Fig6_Data.allTrials_zscore = (Fig6_Data.allTrials - meanLTA) ./ stdLTA;
        
        % Create time axis
        ifi = 18.5/(2*LTAwindow + 1);
        Fig6_Data.timeAxis = linspace(-LTAwindow*ifi, LTAwindow*ifi, size(Fig6_Data.allTrials, 2));
    end
    
    fprintf('  Stimulus LTA trials: %d\n', size(Fig6_Data.allTrials, 1));
end

function [Fig7_Data] = extractStimulusLTM_Entropy(validLickTrials, Entropy, Performance, timeReference, LTAwindow, position)
    % EXTRACTSTIMULUSLTA_ENTROPY - Extract stimulus lick-triggered modulation for Figure 7
    % This requires miss trial averages for baseline correction
    
    fprintf('=== Extracting Stimulus LTM Data ===\n');
    
    Fig7_Data = struct();
    Fig7_Data.allTrials = [];
    Fig7_Data.reactionTimes = [];
    
    condNames = fieldnames(validLickTrials);
    
    for c_idx = 1:length(condNames)
        condName = condNames{c_idx};
        condEntropy = Entropy.(condName);
        condPerformance = Performance.(condName);
        
        recNames = fieldnames(validLickTrials.(condName));
        
        for r_idx = 1:length(recNames)
            recName = recNames{r_idx};
            recNum = str2double(recName(4:end));
            
            if recNum > length(condPerformance)
                continue;
            end
            
            recEntropyData = validLickTrials.(condName).(recName).recEntropyData;
            
            % Get miss trials for baseline using position-specific fields
            missTrials = [];
            missFieldName = sprintf('miss%s', position);
            if isfield(condPerformance(recNum), missFieldName)
                missTrials = condPerformance(recNum).(missFieldName);
            elseif isfield(condPerformance(recNum), 'miss')
                missTrials = condPerformance(recNum).miss;
            end
            
            if isempty(missTrials)
                continue; % Skip if no miss trials available
            end
            
            % Calculate miss trial average
            if max(missTrials) <= size(recEntropyData, 1)
                missAverage = mean(recEntropyData(missTrials, :), 1);
            else
                continue;
            end
            
            allTrials = validLickTrials.(condName).(recName).allTrials;
            allReactionTimes = validLickTrials.(condName).(recName).allReactionTimes;
            
            % Extract LTM for each valid trial
            for i = 1:length(allTrials)
                t = allTrials(i);
                reactionTime = allReactionTimes(i);
                
                lickFrame = find(reactionTime < timeReference, 1);
                if ~isempty(lickFrame)
                    % Extract entropy around lick and subtract miss baseline
                    lickTriggeredEntropy = recEntropyData(t, lickFrame-LTAwindow:lickFrame+LTAwindow);
                    baselineEntropy = missAverage(lickFrame-LTAwindow:lickFrame+LTAwindow);
                    modulatedEntropy = lickTriggeredEntropy - baselineEntropy;
                    
                    Fig7_Data.allTrials = [Fig7_Data.allTrials; modulatedEntropy]; %#ok<AGROW>
                    Fig7_Data.reactionTimes = [Fig7_Data.reactionTimes, reactionTime]; %#ok<AGROW>
                end
            end
        end
    end
    
    % Compute z-scored data
    if ~isempty(Fig7_Data.allTrials)
        meanLTM = mean(Fig7_Data.allTrials(:));
        stdLTM = std(Fig7_Data.allTrials(:));
        if stdLTM == 0; stdLTM = 1; end
        Fig7_Data.allTrials_zscore = (Fig7_Data.allTrials - meanLTM) ./ stdLTM;
        
        % Create time axis
        ifi = 18.5/(2*LTAwindow + 1);
        Fig7_Data.timeAxis = linspace(-LTAwindow*ifi, LTAwindow*ifi, size(Fig7_Data.allTrials, 2));
    end
    
    fprintf('  Stimulus LTM trials: %d\n', size(Fig7_Data.allTrials, 1));
end

function [Fig8_Data] = extractRewardRTA_Entropy(validLickTrials, timeReference, RTAwindow)
    % EXTRACTREWARDRTA_ENTROPY - Extract reward-triggered averages for Figure 8
    
    fprintf('=== Extracting Reward RTA Data ===\n');
    
    Fig8_Data = struct();
    Fig8_Data.allTrials = [];
    Fig8_Data.reactionTimes = [];
    
    condNames = fieldnames(validLickTrials);
    
    for c_idx = 1:length(condNames)
        condName = condNames{c_idx};
        
        recNames = fieldnames(validLickTrials.(condName));
        
        for r_idx = 1:length(recNames)
            recName = recNames{r_idx};
            
                 hitTrials = validLickTrials.(condName).(recName).allTrials;
            hitReactionTimes = validLickTrials.(condName).(recName).allReactionTimes;
            recEntropyData = validLickTrials.(condName).(recName).recEntropyData;
            
            % Extract RTA for each valid hit trial
            for i = 1:length(hitTrials)
                t = hitTrials(i);
                reactionTime = hitReactionTimes(i);
                
                firstLickFrame = find(reactionTime < timeReference, 1);
                if ~isempty(firstLickFrame)
                    rewardFrame = firstLickFrame + 4; % 350ms after lick
                    
                    % Extract entropy around reward delivery
                    rewardTriggeredEntropy = recEntropyData(t, rewardFrame-RTAwindow:rewardFrame+RTAwindow);
                    Fig8_Data.allTrials = [Fig8_Data.allTrials; rewardTriggeredEntropy]; %#ok<AGROW>
                    Fig8_Data.reactionTimes = [Fig8_Data.reactionTimes, reactionTime]; %#ok<AGROW>
                end
            end
        end
    end
    
    % Compute z-scored data
    if ~isempty(Fig8_Data.allTrials)
        meanRTA = mean(Fig8_Data.allTrials(:));
        stdRTA = std(Fig8_Data.allTrials(:));
        if stdRTA == 0; stdRTA = 1; end
        Fig8_Data.allTrials_zscore = (Fig8_Data.allTrials - meanRTA) ./ stdRTA;
        
        % Create time axis
        ifi = 18.5/(2*RTAwindow + 1);
        Fig8_Data.timeAxis = linspace(-RTAwindow*ifi, RTAwindow*ifi, size(Fig8_Data.allTrials, 2));
    end
    
    fprintf('  Reward RTA trials: %d\n', size(Fig8_Data.allTrials, 1));
end

function [Fig9_Data] = extractRewardModulation_Entropy(validLickTrials, timeReference)
    % EXTRACTREWARDMODULATION_ENTROPY - Extract reward modulation data for Figure 9
    
    fprintf('=== Extracting Reward Modulation Data ===\n');
    
    Fig9_Data = struct();
    Fig9_Data.PreReward = [];
    Fig9_Data.PostReward = [];
    
    condNames = fieldnames(validLickTrials);
    
    for c_idx = 1:length(condNames)
        condName = condNames{c_idx};
        
        recNames = fieldnames(validLickTrials.(condName));
        
        for r_idx = 1:length(recNames)
            recName = recNames{r_idx};
            
            hitTrials = validLickTrials.(condName).(recName).hitTrials;
            hitReactionTimes = validLickTrials.(condName).(recName).hitReactionTimes;
            recEntropyData = validLickTrials.(condName).(recName).recEntropyData;
            
            % Extract reward modulation for each hit trial
            for i = 1:length(hitTrials)
                t = hitTrials(i);
                reactionTime = hitReactionTimes(i);
                
                lickFrame = find(reactionTime < timeReference, 1);
                if ~isempty(lickFrame) && lickFrame + 6 <= size(recEntropyData, 2)
                    % Pre-reward: 3 frames after lick (300ms), Post-reward: 6 frames after lick (600ms)
                    preRewardEntropy = recEntropyData(t, lickFrame+3);
                    postRewardEntropy = recEntropyData(t, lickFrame+6);
                    
                    Fig9_Data.PreReward = [Fig9_Data.PreReward; preRewardEntropy]; %#ok<AGROW>
                    Fig9_Data.PostReward = [Fig9_Data.PostReward; postRewardEntropy]; %#ok<AGROW>
                end
            end
        end
    end
    
    % Compute z-scored data
    if ~isempty(Fig9_Data.PreReward)
        allRewardData = [Fig9_Data.PreReward; Fig9_Data.PostReward];
        meanReward = mean(allRewardData(:));
        stdReward = std(allRewardData(:));
        if stdReward == 0; stdReward = 1; end
        
        Fig9_Data.PreReward_zscore = (Fig9_Data.PreReward - meanReward) ./ stdReward;
        Fig9_Data.PostReward_zscore = (Fig9_Data.PostReward - meanReward) ./ stdReward;
    end
    
    fprintf('  Reward modulation trials: %d\n', size(Fig9_Data.PreReward, 1));
end

function [Fig10_Data] = extractPreviousTrialHistory_Entropy(conditionsWithLicks, Entropy, Performance, Rec, Skip, position)
    % EXTRACTPREVIOUSTRIALHISTORY_ENTROPY - Extract previous trial history data for Figure 10
    
    fprintf('=== Extracting Previous Trial History Data ===\n');
    
    Fig10_Data = struct();
    Fig10_Data.PreviousHit = [];
    Fig10_Data.PreviousMiss = [];
    
    % Find minimum timepoints across all recordings
    minTimepoints = inf;
    for c_idx = 1:length(conditionsWithLicks)
        condName = conditionsWithLicks{c_idx};
        if isfield(Entropy, condName)
            condEntropy = Entropy.(condName);
            if isfield(condEntropy, 'AllNeurons') && isfield(condEntropy.AllNeurons, 'RecordingEntropy_Raw')
                if iscell(condEntropy.AllNeurons.RecordingEntropy_Raw)
                    for r = 1:size(condEntropy.AllNeurons.RecordingEntropy_Raw, 1)
                        if ~isempty(condEntropy.AllNeurons.RecordingEntropy_Raw{r, 1})
                            minTimepoints = min(minTimepoints, size(condEntropy.AllNeurons.RecordingEntropy_Raw{r, 1}, 2));
                        end
                    end
                end
            end
        end
    end
    
    for c_idx = 1:length(conditionsWithLicks)
        condName = conditionsWithLicks{c_idx};
        
        % Get recordings to skip
        recordingsToSkip = [];
        if isfield(Skip, condName); recordingsToSkip = Skip.(condName); end
        
        if ~isfield(Entropy, condName) || ~isfield(Performance, condName)
            continue;
        end
        
        condEntropy = Entropy.(condName);
        condPerformance = Performance.(condName);
        
        % Determine number of recordings
        if isfield(Rec.(condName), 'AnimalID')
            numRecs = length(Rec.(condName).AnimalID);
        else
            numRecs = size(Rec.(condName), 1);
        end
        
        % Process each recording
        for r = 1:numRecs
            if ismember(r, recordingsToSkip); continue; end
            
            if r > length(condPerformance)
                continue;
            end
            
            % Get entropy data
            if iscell(condEntropy.AllNeurons.RecordingEntropy_Raw) && r <= size(condEntropy.AllNeurons.RecordingEntropy_Raw, 1)
                if ~isempty(condEntropy.AllNeurons.RecordingEntropy_Raw{r, 1})
                    recEntropyData = condEntropy.AllNeurons.RecordingEntropy_Raw{r, 1};
                else
                    continue;
                end
            else
                continue;
            end
            
            % Get hit/miss trials using position-specific fields
            hitTrials = [];
            missTrials = [];
            
            hitFieldName = sprintf('hit%s', position);
            missFieldName = sprintf('miss%s', position);
            
            if isfield(condPerformance(r), hitFieldName)
                hitTrials = condPerformance(r).(hitFieldName);
            elseif isfield(condPerformance(r), 'hit')
                hitTrials = condPerformance(r).hit;
            end
            
            if isfield(condPerformance(r), missFieldName)
                missTrials = condPerformance(r).(missFieldName);
            elseif isfield(condPerformance(r), 'miss')
                missTrials = condPerformance(r).miss;
            end
            
            if isempty(hitTrials) %|| length(hitTrials) < 2
                continue;
            end
            
            % Find trials that come after hits (exclude first trial)
            hitTrials = hitTrials(hitTrials > 1);
            previousHitTrials = hitTrials - 1;
            
            % Find trials that come after misses
            if ~isempty(missTrials)
                missTrials = missTrials(missTrials > 1);
                previousMissTrials = missTrials - 1;
            else
                previousMissTrials = [];
            end
            
            % Extract entropy for trials following hits vs misses
            if ~isempty(previousHitTrials) && max(previousHitTrials) <= size(recEntropyData, 1)
                % Average across trials, truncate to minTimepoints
                temp = mean(recEntropyData(previousHitTrials, 1:minTimepoints), 1); % [1 x minTimepoints]
                
                % Ensure it's a row vector [1 x minTimepoints] - force reshape
                prevHitEntropy = reshape(temp, 1, []);
                
                Fig10_Data.PreviousHit = [Fig10_Data.PreviousHit; prevHitEntropy]; %#ok<AGROW>
            end
            
            if ~isempty(previousMissTrials) && max(previousMissTrials) <= size(recEntropyData, 1)
                % Average across trials, truncate to minTimepoints
                temp = mean(recEntropyData(previousMissTrials, 1:minTimepoints), 1); % [1 x minTimepoints]
                
                % Ensure it's a row vector [1 x minTimepoints] - force reshape
                prevMissEntropy = reshape(temp, 1, []);
                
                Fig10_Data.PreviousMiss = [Fig10_Data.PreviousMiss; prevMissEntropy]; %#ok<AGROW>
            end
        end
    end
    
    % Compute z-scored data
    if ~isempty(Fig10_Data.PreviousHit) && ~isempty(Fig10_Data.PreviousMiss)
        allPrevData = [Fig10_Data.PreviousHit; Fig10_Data.PreviousMiss];
        meanPrev = mean(allPrevData(:));
        stdPrev = std(allPrevData(:));
        if stdPrev == 0; stdPrev = 1; end
        
        Fig10_Data.PreviousHit_zscore = (Fig10_Data.PreviousHit - meanPrev) ./ stdPrev;
        Fig10_Data.PreviousMiss_zscore = (Fig10_Data.PreviousMiss - meanPrev) ./ stdPrev;
        
        % Create time axis
        Fig10_Data.timeAxis = linspace(-8, 10.5, size(Fig10_Data.PreviousHit, 2));
    end
    
    fprintf('  Previous trial history: %d hit-following, %d miss-following\n', ...
            size(Fig10_Data.PreviousHit, 1), size(Fig10_Data.PreviousMiss, 1));
end

function validateTrialCounts_Entropy(validLickTrials, Fig1_Data, Fig6_Data, Fig8_Data, Fig9_Data, position)
    % VALIDATETRIALCOUNTS_ENTROPY - Validate consistency of trial counts across different analyses
    
    fprintf('\n=== TRIAL COUNT VALIDATION for %s ===\n', position);
    
    % Count total valid trials across all conditions
    totalValidTrials = 0;
    totalHitTrials = 0;
    
    condNames = fieldnames(validLickTrials);
    for c_idx = 1:length(condNames)
        condName = condNames{c_idx};
        recNames = fieldnames(validLickTrials.(condName));
        
        for r_idx = 1:length(recNames)
            recName = recNames{r_idx};
            totalValidTrials = totalValidTrials + length(validLickTrials.(condName).(recName).allTrials);
            totalHitTrials = totalHitTrials + length(validLickTrials.(condName).(recName).hitTrials);
        end
    end
    
    % Compare across analyses
    fprintf('Expected from unified trial selection:\n');
    fprintf('  Total valid trials (all licks): %d\n', totalValidTrials);
    fprintf('  Total hit trials: %d\n', totalHitTrials);
    
    fprintf('\nActual from data extraction:\n');
    fprintf('  Fig 1 Stimulus licks: %d\n', size(Fig1_Data.Stimulus.Pre, 1));
    fprintf('  Fig 1 Reward timepoints: %d\n', size(Fig1_Data.Reward.Pre, 1));
    fprintf('  Fig 6 Stimulus LTA: %d\n', size(Fig6_Data.allTrials, 1));
    fprintf('  Fig 8 Reward RTA: %d\n', size(Fig8_Data.allTrials, 1));
    fprintf('  Fig 9 Reward Modulation: %d\n', size(Fig9_Data.PreReward, 1));
    
    % Check consistency
    consistency_check = true;
    
    if size(Fig1_Data.Stimulus.Pre, 1) ~= totalValidTrials
        fprintf('WARNING: Fig 1 Stimulus count mismatch!\n');
        consistency_check = false;
    end
    
    if size(Fig1_Data.Reward.Pre, 1) ~= totalHitTrials
        fprintf('WARNING: Fig 1 Reward count mismatch!\n');
        consistency_check = false;
    end
    
    if size(Fig6_Data.allTrials, 1) ~= totalValidTrials
        fprintf('WARNING: Fig 6 LTA count mismatch!\n');
        consistency_check = false;
    end
    
    if size(Fig8_Data.allTrials, 1) ~= totalHitTrials
        fprintf('WARNING: Fig 8 RTA count mismatch!\n');
        consistency_check = false;
    end
    
    if size(Fig9_Data.PreReward, 1) ~= totalHitTrials
        fprintf('WARNING: Fig 9 Reward Modulation count mismatch!\n');
        consistency_check = false;
    end
    
    if consistency_check
        fprintf('✓ All trial counts are CONSISTENT!\n');
    else
        fprintf('✗ Trial count INCONSISTENCIES detected - please review\n');
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