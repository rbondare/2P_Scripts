data= load('\\fs.ista.ac.at\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB6_250302_1656_preprocessed.mat');
%% plotting Positioning of ROIs

%extract ROI position (both planes here)
centroid = data.CaData(1).Ca_centroid_voxel;
centroidX = centroid(:, 1);
centroidY = centroid(:, 2);
centroidZ = centroid(:, 3);

figure;
scatter(centroidY, centroidX)
set(gca, 'YDir', 'reverse'); % This inverts the y-axis
%% Normalising Calcium Activity

dFF = data.CaData(1).Ca_dFF;

%different normalisations 

%normalize by min&max - worst 
%dFF_norm = (dFF - min(dFF,[],2)) ./ (max(dFF,[],2) - min(dFF,[],2)); 
%-zsocre
dFF_z = (dFF - mean(dFF, 2)) ./ std(dFF, 0, 2); 
% 0 to max normalisation 
dFF_norm = (dFF - min(dFF,[],2)) ./ (max(dFF,[],2) - min(dFF,[],2)) .* max(dFF,[],2);
% mean normlisation 
dFF_norm_mean = dFF ./ mean(dFF, 2);
%% Plotting Heatmap of Calcium Activity


roi_indices = (1:11); %ALL ROIs  
%roi_indices = setdiff(1:11, 3); % All ROIs from 1 to 11 except 3

dFF_selected = dFF(roi_indices, 3000:4000);
%dFF_norm_selected = dFF_norm(roi_indices, 3000:4000); 
dFF_z_selected = dFF_z(roi_indices, 3000:4000);
dFF_norm_selected = dFF_norm(roi_indices, 3000:4000);
dFF_norm_mean_selected = dFF_norm_mean(roi_indices, 3000:4000);

%sort ROIs by maximum
[~, sort_idx] = sort(max(dFF_selected, [], 2), 'descend'); % Change to 'ascend' if you want lowest to highest

dFF_selected = dFF_selected(sort_idx, :);
%dFF_norm_selected = dFF_norm_selected(sort_idx, :);
dFF_z_selected = dFF_z_selected(sort_idx, :);
dFF_norm_selected = dFF_norm_selected(sort_idx, :);
dFF_norm_mean_selected = dFF_norm_mean_selected(sort_idx, :);

original_colormap = gray;  % Get the default gray colormap
reversed_colormap = flipud(original_colormap);

figure;
imagesc(dFF_selected)
title("Raw")
colormap('turbo'); % Use a colormap like 'hot' or 'jet'
colorbar;  % Add a color scale

figure;
imagesc(dFF_z_selected)
title("z-scored")
colormap('turbo'); 
colorbar;  % Add a color scale

figure;
imagesc(dFF_norm_selected)
title("normalized by max")
colormap(turbo); 
colorbar;
c = colorbar;  % Add a color scale
caxis([0 12]); % Adjust as needed based on the range of your data
%% Aligning Calcium Activity with Behaviour
% Plot Calcium Activity traces with pupil dynamics and locomotion 

% calcium activity is aligned to the stimulus using the same time vector 
% Get stimulus timing
stim_total_time = data.Stimuli(2).stimulus_trial_t * data.Stimuli(2).trials
stim_start = data.Stimuli(2).TimeStimulusFrame(1)
stim_end = stim_start + stim_total_time

% Find frames of calcium for defined stimulus timing
Ca_index = find(data.TimeCa(1, :) > stim_start & data.TimeCa(1, :) < stim_end);
Ca_selected = dFF_z(:,Ca_index);

% Find behavior data for defined stimulus timing
Beh_index = find(data.Triggers.TimeCamera > stim_start & data.Triggers.TimeCamera < stim_end);
Beh_selected = data.Behav.PupilArea(Beh_index);

% Find locomotion data for defined stimulus timing
Loc_index = find(data.Triggers.TimeBall{4, 1} > stim_start & data.Triggers.TimeBall{4, 1} < stim_end);
Loc_time = data.Triggers.TimeBall{4, 1}(Loc_index); % Extract timestamps
Loc_position = data.LocomotionCal(1).Forward_mm(Loc_index); % Extract position

% Compute locomotion speed (mm/s)
Loc_velocity = diff(Loc_position) ./ diff(Loc_time); 
Loc_time_velocity = Loc_time(1:end-1); % Adjust time vector for velocity

% Create time vector
time_vector = linspace(stim_start, stim_end, size(Ca_selected,2));

figure;

% Plot calcium data (first ROI only)
subplot(3,1,1);
plot(time_vector, Ca_selected(4,:), 'LineWidth', 1);
title('Calcium Activity (ROI 1)');
ylabel('dF/F');
ylim([-3 15]);
ax1 = gca; % Save axis handle for linking

% Plot behavioral data (Pupil Area)
subplot(3,1,2);
plot(linspace(stim_start, stim_end, length(Beh_selected)), Beh_selected, 'r', 'LineWidth', 1);
title('Pupil Area');
ylabel('Pupil Area');
ax2 = gca; % Save axis handle for linking

% Interpolate locomotion velocity to match calcium time vector
Loc_velocity_resampled = interp1(Loc_time_velocity, Loc_velocity, time_vector, 'linear', 'extrap');

% Plot locomotion speed
subplot(3,1,3);
plot(time_vector, Loc_velocity_resampled, 'k', 'LineWidth', 1);
title('Locomotion Speed');
xlabel('Time (s)');
ylabel('Speed (mm/s)');
ax3 = gca; % Save axis handle for linking

% Link the x-axes so they zoom/pan together
linkaxes([ax1, ax2, ax3], 'x');

% Adjust spacing and figure size
set(gcf, 'Position', [100, 100, 900,700]); % Make figure larger
%% Generating fff_stimulus 

ifi = data.Stimuli(1).ifi(1) 
stim_start = data.Stimuli(2).TimeStimulusFrame(1);
options = data.Stimuli(2);%.specParams;


[fff_brightness, single_pulse] = mk_fff_curve(options, ifi, stim_start);

%% 
% Aligning Calcium Data with locomotion and Pupil 

% Get stimulus timing for a particular stimulus
stim_idx = 2;
roi_id=4;
roi_plane=data.CaData(1).Ca_centroid_voxel(roi_id,3);
stim_total_time = data.Stimuli(stim_idx).stimulus_trial_t * data.Stimuli(stim_idx).trials
stim_start = data.Stimuli(stim_idx).TimeStimulusFrame(1)
stim_end = stim_start + stim_total_time

% Find frames of calcium for defined stimulus timing
Ca_index = find(data.TimeCa(roi_plane, :) > stim_start & data.TimeCa(roi_plane, :) < stim_end);
Ca_selected = data.CaData(roi_plane).Ca_dFF(:, Ca_index);

% Find behavior data for defined stimulus timing
Beh_index = find(data.Triggers.TimeCamera > stim_start & data.Triggers.TimeCamera < stim_end);
Beh_selected = data.Behav.PupilArea(Beh_index);

% Find locomotion data for defined stimulus timing
 Loc_index = find(data.Triggers.TimeBall{stim_idx *2} > stim_start & data.Triggers.TimeBall{stim_idx *2} < stim_end);
 Loc_time = data.Triggers.TimeBall{stim_idx *2}(Loc_index); % Extract timestamps
 Loc_position = data.LocomotionCal(1).Forward_mm(Loc_index); % Extract position

% Compute locomotion speed (mm/s), handling edge case
if length(Loc_time) > 1
    Loc_velocity = diff(Loc_position) ./ diff(Loc_time); 
    Loc_time_velocity = Loc_time(1:end-1); % Adjust time vector for velocity
else
    Loc_velocity = zeros(1, length(Loc_time)); % Prevent errors
    Loc_time_velocity = Loc_time;
end

% Create time vector for consistent x-axis across all plots
time_vector = linspace(stim_start, stim_end, size(Ca_selected,2));

% Create figure for Stimulus 2
figure('Name', 'Stimulus 2');

%plot the stimulus
subplot(3,1,1);
num_repeats = 5;  % Number of stimulus repetitions
fff_brightness_repeated = repmat(fff_brightness, 1, num_repeats); 
plot(linspace(stim_start, stim_end, length(fff_brightness_repeated)), fff_brightness_repeated, 'k', 'LineWidth', 1.5);
title ("Stimulus")
axis off; 
set(gca, 'Box', 'off');          
xlim([stim_start, stim_end])
        

% Plot calcium data (1 ROI only)
subplot(3,1,2);
plot(time_vector, Ca_selected(roi_id,:), 'LineWidth', 1);
title('Calcium Activity (ROI 1)');
ylabel('dF/F');
ylim([-3 15]);
set(gca, 'Box', 'off');  
xlim([stim_start, stim_end])

% Plot behavioral data (Pupil Area)
subplot(3,1,3);
plot(linspace(stim_start, stim_end, length(Beh_selected)), Beh_selected, 'r', 'LineWidth', 1);
title('Pupil Area');
ylabel('Pupil Area');
set(gca, 'Box', 'off');  
xlim([stim_start, stim_end])

% Interpolate locomotion velocity to match calcium time vector
if length(Loc_time_velocity) > 1
    Loc_velocity_resampled = interp1(Loc_time_velocity, Loc_velocity, time_vector, 'linear', 'extrap');
else
    Loc_velocity_resampled = zeros(size(time_vector)); % No movement case
end

% % Plot locomotion speed
% subplot(4,1,4);
% plot(time_vector, Loc_velocity_resampled, 'k', 'LineWidth', 1);
% title('Locomotion Speed');
% xlabel('Time (s)');
% ylabel('Speed (mm/s)');
% ax3 = gca; % Save axis handle for linking

% Link the x-axes so they zoom/pan together

% Adjust spacing and figure size
set(gcf, 'Position', [100, 100, 800, 600]); % Make figure larger
%% 
% Chunking data by trial 

num_trials = data.Stimuli.trials
no_stimuli = 4;
chunk_size = length(data.Stimuli(1).TimeStimulusFrame)/num_trials;

for idx = 1:no_stimuli 
    for no_trial = 1:num_trials
      desired_vector = data.Stimuli(idx).TimeStimulusFrame;
      desired_vector = desired_vector(chunk_size * (no_trial -1) + 1 : chunk_size * no_trial);
      data.Stimuli(idx).TrialTimes{no_trial} = desired_vector;
    end
end

%%
num_stimuli = length(data.Stimuli); % Number of stimuli
trial_number = 5;  %number of trials 

for stim_idx = 1:num_stimuli
    for trial_times =1:trial_number
        % Get stimulus timing
        trial_time = data.Stimuli(stim_idx).stimulus_trial_t;
        stim_start = data.Stimuli(stim_idx).TrialTimes{1, trial_times}(1);
        stim_end = stim_start + trial_time;
    
        % Find frames of calcium for defined stimulus timing
        Ca_index = find(data.TimeCa(1, :) > stim_start & data.TimeCa(1, :) < stim_end);
        Ca_selected = data.CaData(1).Ca_dFF(4, Ca_index);
        
    end
end 
%%
% Generate two plots for each stimulus:
% 1) All trials in grey with average in black
% 2) Average with standard deviation cloud

num_stimuli = length(data.Stimuli); % Number of stimuli
trial_number = 5;  % number of trials 
roi_number = 4;    % ROI number to plot

% Create figures for each stimulus
for stim_idx = 1:num_stimuli
    % Storage for all trials to calculate average later
    all_trials_ca = cell(1, trial_number);
    all_trials_time = cell(1, trial_number);
    
    % First collect all trials for this stimulus
    for trial_times = 1:trial_number
        % Get stimulus timing
        trial_time = data.Stimuli(stim_idx).stimulus_trial_t;
        stim_start = data.Stimuli(stim_idx).TrialTimes{1, trial_times}(1);
        stim_end = stim_start + trial_time;
    
        % Find frames of calcium for defined stimulus timing
        Ca_index = find(data.TimeCa(1, :) > stim_start & data.TimeCa(1, :) < stim_end);
        Ca_selected = data.CaData(1).Ca_dFF(roi_number, Ca_index);
        
        % Create time vector for this trial
        time_vector = linspace(0, trial_time, length(Ca_selected));
        
        % Store data for averaging later
        all_trials_ca{trial_times} = Ca_selected;
        all_trials_time{trial_times} = time_vector;
    end
    
    % Create a common time vector for resampling
    common_time = linspace(0, trial_time, 1000); % 1000 points for smooth curve
    resampled_ca = zeros(trial_number, length(common_time));
    
    % Resample all trials to common time base
    for trial_times = 1:trial_number
        resampled_ca(trial_times, :) = interp1(all_trials_time{trial_times}, all_trials_ca{trial_times}, common_time, 'linear', 'extrap');
    end
    
    % Calculate average and standard deviation across trials
    avg_ca = mean(resampled_ca, 1);
    std_ca = std(resampled_ca, 0, 1);
    
    % Calculate upper and lower bounds for shaded area
    upper_bound = avg_ca + std_ca;
    lower_bound = avg_ca - std_ca;
    
    %% PLOT 1: All trials in grey with average in black
    figure('Name', sprintf('Calcium Activity - All Trials - Stimulus %d', stim_idx));
    hold on;
    
    % Plot each trial in grey
    for trial_times = 1:trial_number
        plot(common_time, resampled_ca(trial_times, :), 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    end
    
    % Plot average in thick black line
    plot(common_time, avg_ca, 'k', 'LineWidth', 2);
    
    % Add labels and title
    xlabel('Time (s)');
    ylabel('dF/F');
    title(sprintf('Calcium Activity (ROI %d) - Stimulus %d', roi_number, stim_idx));
    
    % Add legend
    legend('Individual Trials', 'Average', 'Location', 'best');
    
    % Set axes limits and appearance
    ylim([-3 15]);  % Based on your previous code
    set(gca, 'Box', 'off');
    grid on;
    
    % Adjust spacing and figure size
    set(gcf, 'Position', [100, 100, 800, 400]);
    
%     %% PLOT 2: Average with standard deviation cloud
%     figure('Name', sprintf('Average Calcium Activity - Stimulus %d', stim_idx));
%     hold on;
%     
%     % Create shaded area for standard deviation
%     x_shade = [common_time, fliplr(common_time)];
%     y_shade = [upper_bound, fliplr(lower_bound)];
%     
%     % Plot the shaded area
%     fill(x_shade, y_shade, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%     
%     % Plot average in thick black line
%     plot(common_time, avg_ca, 'k', 'LineWidth', 2);
%     
%     % Add labels and title
%     xlabel('Time (s)');
%     ylabel('dF/F');
%     title(sprintf('Average Calcium Activity (ROI %d) - Stimulus %d', roi_number, stim_idx));
%     
%     % Add legend
%     legend('Standard Deviation', 'Average', 'Location', 'best');
%     
%     % Set axes limits and appearance
%     ylim([-3 15]);  % Based on your previous code
%     set(gca, 'Box', 'off');
%     grid on;
%     
%     % Adjust spacing and figure size
%     set(gcf, 'Position', [100, 100, 800, 400]);
end

%%

% Generate two plots for each stimulus:
% 1) Average of all ROIs for each trial (in grey) with overall average in black
% 2) Overall average with standard deviation cloud

num_stimuli = length(data.Stimuli); % Number of stimuli
trial_number = 5;  % number of trials 

% Get total number of ROIs from the data structure
num_rois = size(data.CaData(1).Ca_dFF, 1);

% Create figures for each stimulus
for stim_idx = 1:num_stimuli
    % Storage for trial averages (each trial being an average across all ROIs)
    trial_averages = cell(1, trial_number);
    all_trials_time = cell(1, trial_number);
    
    % First calculate the average across all ROIs for each trial
    for trial_times = 1:trial_number
        % Get stimulus timing
        trial_time = data.Stimuli(stim_idx).stimulus_trial_t;
        stim_start = data.Stimuli(stim_idx).TrialTimes{1, trial_times}(1);
        stim_end = stim_start + trial_time;
    
        % Find frames of calcium for defined stimulus timing
        Ca_index = find(data.TimeCa(1, :) > stim_start & data.TimeCa(1, :) < stim_end);
        
        % Get calcium data for all ROIs for this trial
        Ca_all_rois = data.CaData(1).Ca_dFF(:, Ca_index);
        
        % Calculate average across all ROIs for this trial
        Ca_trial_avg = mean(Ca_all_rois, 1);
        
        % Create time vector for this trial
        time_vector = linspace(0, trial_time, length(Ca_trial_avg));
        
        % Store the ROI-averaged data for this trial
        trial_averages{trial_times} = Ca_trial_avg;
        all_trials_time{trial_times} = time_vector;
    end
    
    % Create a common time vector for resampling
    common_time = linspace(0, trial_time, 1000); % 1000 points for smooth curve
    
    % Resample all trial averages to common time base
    resampled_trial_avgs = zeros(trial_number, length(common_time));
    for trial_times = 1:trial_number
        resampled_trial_avgs(trial_times, :) = interp1(all_trials_time{trial_times}, ...
            trial_averages{trial_times}, common_time, 'linear', 'extrap');
    end
    
    % Calculate overall average and standard deviation across trials
    overall_avg = mean(resampled_trial_avgs, 1);
    overall_std = std(resampled_trial_avgs, 0, 1);
    
    % Calculate upper and lower bounds for shaded area
    upper_bound = overall_avg + overall_std;
    lower_bound = overall_avg - overall_std;
    
    %% PLOT 1: Average across all ROIs for each trial (in grey) with overall average in black
    figure('Name', sprintf('Calcium Activity by Trial - Stimulus %d', stim_idx));
    hold on;
    
         
    
    % Plot each trial's average (across all ROIs) in grey
    for trial_times = 1:trial_number
        plot(common_time, resampled_trial_avgs(trial_times, :), 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    end
    
    % Plot overall average in thick black line
    plot(common_time, overall_avg, 'k', 'LineWidth', 2);
    
    % Add labels and title
    xlabel('Time (s)');
    ylabel('dF/F');
    title(sprintf('Average Calcium Activity Across All ROIs - Stimulus %d', stim_idx));
    
    % Add legend
    legend('Trial Averages (All ROIs)', 'Overall Average', 'Location', 'best');
    
    % Set axes limits and appearance
    ylim([-3 15]);  % Based on your previous code
    set(gca, 'Box', 'off');
    grid on;
    
    % Adjust spacing and figure size
    set(gcf, 'Position', [100, 100, 800, 400]);
    
    %% PLOT 2: Overall average with standard deviation cloud
    figure('Name', sprintf('Overall Average Calcium Activity - Stimulus %d', stim_idx));
    hold on;
    
    % Create shaded area for standard deviation
    x_shade = [common_time, fliplr(common_time)];
    y_shade = [upper_bound, fliplr(lower_bound)];
    
    % Plot the shaded area
    fill(x_shade, y_shade, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    % Plot overall average in thick black line
    plot(common_time, overall_avg, 'k', 'LineWidth', 2);
    
    % Add labels and title
    xlabel('Time (s)');
    ylabel('dF/F');
    title(sprintf('Overall Average Calcium Activity Across All ROIs - Stimulus %d', stim_idx));
    
    % Add legend
    legend('Standard Deviation', 'Overall Average', 'Location', 'best');
    
    % Set axes limits and appearance
    ylim([-3 15]);  % Based on your previous code
    set(gca, 'Box', 'off');
    grid on;
    
    % Adjust spacing and figure size
    set(gcf, 'Position', [100, 100, 800, 400]);
end
%% 
% Loop through all the stimuli and trial 
% 
% Plot each trial for each stimulus 



num_stimuli = length(data.Stimuli); % Number of stimuli
trial_number = 5; 

for stim_idx = 1:num_stimuli
    for trial_times =1:trial_number
        % Get stimulus timing
        trial_time = data.Stimuli(stim_idx).stimulus_trial_t;
        stim_start = data.Stimuli(stim_idx).TrialTimes{1, trial_times}(1);
        stim_end = stim_start + trial_time;
    
        % Find frames of calcium for defined stimulus timing
        Ca_index = find(data.TimeCa(1, :) > stim_start & data.TimeCa(1, :) < stim_end);
        Ca_selected = data.CaData(1).Ca_dFF(4, Ca_index);
        
        %z-score and apply moving mean
        %Ca_selected = (Ca_selected - mean(Ca_selected)) / std(Ca_selected);
    
        % Find behavior data for defined stimulus timing
        Beh_index = find(data.Triggers.TimeCamera > stim_start & data.Triggers.TimeCamera < stim_end);
        Beh_selected = data.Behav.PupilArea(Beh_index);
        Whisker_selected = data.Behav.whiskerFollicleA1(Beh_index, 1);

        %z-score and apply moving mean
       % Beh_selected = (Beh_selected - mean(Beh_selected)) / std(Beh_selected);
       % Whisker_selected = (Whisker_selected - mean(Whisker_selected)) / std(Whisker_selected);

        
        % Find locomotion data for defined stimulus timing
        Loc_index = find(data.Triggers.TimeBall{stim_idx *2} > stim_start & data.Triggers.TimeBall{stim_idx *2} < stim_end);
        Loc_time = data.Triggers.TimeBall{stim_idx *2}(Loc_index); % Extract timestamps
        Loc_position = data.LocomotionCal(1).Forward_mm(Loc_index); % Extract position
    
        % Compute locomotion speed (mm/s)
        if length(Loc_time) > 1
            Loc_velocity = diff(Loc_position) ./ diff(Loc_time); 
            Loc_time_velocity = Loc_time(1:end-1); % Adjust time vector for velocity
        else
            Loc_velocity = zeros(1, length(Loc_time)); % Avoid errors if no movement
            Loc_time_velocity = Loc_time;
        end
        
        % Interpolate locomotion velocity to match calcium time vector
        if length(Loc_time_velocity) > 1
            Loc_velocity_resampled = interp1(Loc_time_velocity, Loc_velocity, time_vector, 'linear', 'extrap');
        else
            Loc_velocity_resampled = zeros(size(time_vector)); % No movement case
        end
  
        % Create time vector for consistent x-axis across all plots
        time_vector = linspace(stim_start, stim_end, size(Ca_selected,2));
    
        % Create a new figure for each stimulus
        figure('Name', sprintf('Stimulus %d', stim_idx));
    
        %plot the stimulus
        subplot(4,1,1);
        plot(linspace(stim_start, stim_end, length(fff_brightness)), fff_brightness, 'k', 'LineWidth', 1.5);
        title ("Stimulus")
        axis off; 
        set(gca, 'Box', 'off');          
        
        % Plot calcium data 
        subplot(4,1,2);
        plot(time_vector, Ca_selected(:,:), 'LineWidth', 1.5);
        title(sprintf('Calcium Activity (ROI 1) - Stimulus %d', stim_idx));
        ylabel('dF/F');
        ylim([-3 15]);
        ax1 = gca; % Save axis handle for linking
    
   
        % Plot behavioral data (Pupil Area)
        subplot(4,1,3);
        plot(linspace(stim_start, stim_end, length(Beh_selected)), Beh_selected, 'r', 'LineWidth', 1.5);
        title(sprintf('Pupil Area - Stimulus %d', stim_idx));
        ylabel('Pupil Area');
        ax2 = gca; % Save axis handle for linking
        
        
        % Plot behavioral data (Whisker Movement)
        subplot(4,1,4);
        plot(linspace(stim_start, stim_end, length(Whisker_selected)), Whisker_selected, 'k', 'LineWidth', 1.5);
        title(sprintf('Whisker - Stimulus %d', stim_idx));
        ylabel('Whisker Movement');
        ax3 = gca; % Save axis handle for linking
        
        
%         % Plot locomotion speed
%         subplot(4,1,4);
%         plot(time_vector, Loc_velocity_resampled, 'k', 'LineWidth', 1.5);
%         title(sprintf('Locomotion Speed - Stimulus %d', stim_idx));
%         xlabel('Time (s)');
%         ylabel('Speed (mm/s)');
%         ax3 = gca; % Save axis handle for linking
        

    
        % Link the x-axes so they zoom/pan together
        linkaxes([ax1, ax2, ax3], 'x');
        
        % Adjust spacing and figure size
        set(gcf, 'Position', [100, 100, 800, 600]); % Make figure larger
    end 
end 
    
%%
%%Loop through all the stimulus

num_stimuli = length(data.Stimuli); % Number of stimuli

for stim_idx = 1:num_stimuli
    % Get stimulus timing
    stim_total_time = data.Stimuli(stim_idx).stimulus_trial_t * data.Stimuli(stim_idx).trials;
    stim_start = data.Stimuli(stim_idx).TimeStimulusFrame(1);
    stim_end = stim_start + stim_total_time;

    % Find frames of calcium for defined stimulus timing
    Ca_index = find(data.TimeCa(1, :) > stim_start & data.TimeCa(1, :) < stim_end);
    Ca_selected = data.CaData(1).Ca_dFF(:, Ca_index);

    % Find behavior data for defined stimulus timing
    Beh_index = find(data.Triggers.TimeCamera > stim_start & data.Triggers.TimeCamera < stim_end);
    Beh_selected = data.Behav.PupilArea(Beh_index);

    % Find locomotion data for defined stimulus timing
    Loc_index = find(data.Triggers.TimeBall{stim_idx *2} > stim_start & data.Triggers.TimeBall{stim_idx *2} < stim_end);
    Loc_time = data.Triggers.TimeBall{stim_idx *2}(Loc_index); % Extract timestamps
    Loc_position = data.LocomotionCal(1).Forward_mm(Loc_index); % Extract position

    % Compute locomotion speed (mm/s)
    if length(Loc_time) > 1
        Loc_velocity = diff(Loc_position) ./ diff(Loc_time); 
        Loc_time_velocity = Loc_time(1:end-1); % Adjust time vector for velocity
    else
        Loc_velocity = zeros(1, length(Loc_time)); % Avoid errors if no movement
        Loc_time_velocity = Loc_time;
    end

    % Create time vector for consistent x-axis across all plots
    time_vector = linspace(stim_start, stim_end, size(Ca_selected,2));

    % Create a new figure for each stimulus
    figure('Name', sprintf('Stimulus %d', stim_idx));

    % Plot calcium data (first ROI only)
    subplot(3,1,1);
    plot(time_vector, Ca_selected(4,:), 'LineWidth', 1.5);
    title(sprintf('Calcium Activity (ROI 1) - Stimulus %d', stim_idx));
    ylabel('dF/F');
    ylim([-3 15]);
    ax1 = gca; % Save axis handle for linking

    % Plot behavioral data (Pupil Area)
    subplot(3,1,2);
    plot(linspace(stim_start, stim_end, length(Beh_selected)), Beh_selected, 'r', 'LineWidth', 1.5);
    title(sprintf('Pupil Area - Stimulus %d', stim_idx));
    ylabel('Pupil Area');
    ax2 = gca; % Save axis handle for linking

    % Interpolate locomotion velocity to match calcium time vector
    if length(Loc_time_velocity) > 1
        Loc_velocity_resampled = interp1(Loc_time_velocity, Loc_velocity, time_vector, 'linear', 'extrap');
    else
        Loc_velocity_resampled = zeros(size(time_vector)); % No movement case
    end

    % Plot locomotion speed
    subplot(3,1,3);
    plot(time_vector, Loc_velocity_resampled, 'k', 'LineWidth', 1.5);
    title(sprintf('Locomotion Speed - Stimulus %d', stim_idx));
    xlabel('Time (s)');
    ylabel('Speed (mm/s)');
    ax3 = gca; % Save axis handle for linking

    % Link the x-axes so they zoom/pan together
    linkaxes([ax1, ax2, ax3], 'x');

    % Adjust spacing and figure size
    set(gcf, 'Position', [100, 100, 800, 600]); % Make figure larger
end

%%
dFF_norm = (dFF - min(dFF,[],2)) ./ (max(dFF,[],2) - min(dFF,[],2)) .* max(dFF,[],2);


% Get stimulus timing
trial_time = data.Stimuli(2).stimulus_trial_t;
stim_start = data.Stimuli(2).TrialTimes{1, 2}(1);
stim_end = stim_start + trial_time;
    
% Find frames of calcium for defined stimulus timing
Ca_index = find(data.TimeCa(1, :) > stim_start & data.TimeCa(1, :) < stim_end)
Ca_selected = dFF_norm(:, Ca_index);

[~, sort_idx] = sort(max(Ca_selected, [], 2), 'descend'); % Change to 'ascend' if you want lowest to highest
Ca_selected = Ca_selected(sort_idx, :);
    
% Find behavior data for defined stimulus timing
Beh_index = find(data.Triggers.TimeCamera > stim_start & data.Triggers.TimeCamera < stim_end);
Pupil_selected = data.Behav.PupilArea(Beh_index);
Whisker_selected = data.Behav.whiskerFollicleA1(Beh_index, 1);

Pupil_selected_norm = (Pupil_selected - min(Pupil_selected)) / (max(Pupil_selected) - min(Pupil_selected));
Whisker_selected_norm = (Whisker_selected - min(Whisker_selected)) / (max(Whisker_selected) - min(Whisker_selected));


% Create time vector for consistent x-axis across all plots
time_vector = linspace(stim_start, stim_end, size(Ca_selected,2));
    
% Create a figure with two subplots
figure;

% === Calcium Activity Heatmap ===

%plot the stimulus
subplot(3,1,1);
plot(linspace(stim_start, stim_end, length(fff_brightness)), fff_brightness, 'k', 'LineWidth', 1.5);
title ("Stimulus")
axis off; 
set(gca, 'Box', 'off');  

subplot(3,1,2); % 
imagesc(time_vector, [1 size(Ca_selected,1)], Ca_selected); % Use the time_vector
colormap turbo;
xlabel('Time (s)');
ylabel('Neuron #');
title('Calcium Activity Heatmap');


% === Subplot 2: Pupil and Whisker Traces ===
subplot(3,1,3); 
axis off; 
set(gca, 'Box', 'off');  
hold on;
plot(linspace(stim_start, stim_end, length(Pupil_selected)), Pupil_selected_norm, 'k', 'LineWidth', 1); 
plot(linspace(stim_start, stim_end, length(Whisker_selected)), Whisker_selected_norm, 'r', 'LineWidth', 1);
hold off;

%% 
% Computing correlations

%computing correlations

trial_time = data.Stimuli(2).stimulus_trial_t;
stim_start = data.Stimuli(2).TrialTimes{1, 2}(1);
stim_end = stim_start + trial_time;

%for the whole stimulus 

% size(Beh_selected);
% Ca_sel = Ca_selected(4,:);
% 
% [p, q] = rat(10/50);  % p/q = 1/5
% pupil_data_resampled = resample(Beh_selected, p, q);
% 
% min_length = min(length(Ca_sel), length(pupil_data_resampled));
% calcium_data_trimmed = Ca_sel(1:min_length);
% pupil_data_trimmed = pupil_data_resampled(1:min_length);
% 
% 
% [r, lags] = xcorr(calcium_data_trimmed, pupil_data_resampled, 'coeff');

%for individual trials 

%num_stimuli = length(data.Stimuli); % Number of stimuli

% Get stimulus timing

trial_number = 5; 
Ca_traces = struct(); 

%for stim_idx = 1:num_stimuli
 for trial_times =1:trial_number
    % Get stimulus timing
    trial_time = data.Stimuli(2).stimulus_trial_t;
    stim_start = data.Stimuli(2).TrialTimes{1, trial_times}(1);
    stim_end = stim_start + trial_time;
    
%     stim_start = stim_start + 2; % Add 2 seconds to start
%     stim_end = stim_end - 2;     % Subtract 2 seconds from end

    % Find frames of calcium for defined stimulus timing
    Ca_index = find(data.TimeCa(1, :) > stim_start & data.TimeCa(1, :) < stim_end);
    Ca_selected = data.CaData(1).Ca_dFF(:, Ca_index);

    % Find behavior data for defined stimulus timing
    Beh_index = find(data.Triggers.TimeCamera > stim_start & data.Triggers.TimeCamera < stim_end);
    Beh_selected = data.Behav.PupilArea(Beh_index);
    
    % Store calcium trace for this trial
    Ca_traces.(['trial_', num2str(trial_times)]) = Ca_selected;
    Beh_traces.(['trial_', num2str(trial_times)]) = Beh_selected;  
    
end 
%end 

Beh_selected = Beh_traces.trial_2;    
Ca_sel = Ca_traces.trial_2(4,:);

% stimulus intensity vector 

[p, q] = rat(10/50);  % p/q = 1/5

pupil_data_resampled = resample(Beh_selected, p, q);

min_length = min(length(Ca_sel), length(pupil_data_resampled));
calcium_data_trimmed = Ca_sel(1:min_length);
pupil_data_trimmed = pupil_data_resampled(1:min_length);

% z-score signals
calcium_norm = (calcium_data_trimmed - mean(calcium_data_trimmed)) / std(calcium_data_trimmed);
pupil_norm = (pupil_data_trimmed - mean(pupil_data_trimmed)) / std(pupil_data_trimmed);

%smooth using mean
calcium_norm = movmean(calcium_norm, 10);
pupil_norm = movmean(pupil_norm, 10);

% Calculate cross-correlation
[r, lags] = xcorr(calcium_data_trimmed, pupil_data_trimmed, 'coeff');

% Plot results
figure;
plot(lags/10, r);  % Convert lags to seconds
xlabel('Lag (seconds)');
ylabel('Correlation coefficient');
title('Cross-correlation between calcium activity and pupil dynamics');

% Create figure with subplots
figure('Position', [100, 100, 1000, 800]);

% Plot 1: Original signals (time-aligned and normalized)
subplot(3,1,1);
t_calcium = (0:length(calcium_norm)-1)/10; % Time in seconds
plot(t_calcium, calcium_norm, 'b', 'LineWidth', 1.5);
hold on;
plot(t_calcium, pupil_norm, 'r', 'LineWidth', 1.5);
legend('Calcium Signal (normalized)', 'Pupil Signal (normalized)');
xlabel('Time (s)');
% ylim([0 1]);
set(gca, 'Box', 'off');  

% Plot 2: Cross-correlation function
subplot(3,1,2);
lag_time = lags/10; % Convert lags to seconds
plot(lag_time, r, 'k', 'LineWidth', 1.5);
title('Cross-correlation between Calcium and Pupil Signals');
xlabel('Lag (s)');
ylabel('Correlation Coefficient');
grid on;

% Find and mark the peak correlation
[max_corr, max_idx] = max(abs(r));
lag_max = lag_time(max_idx);
hold on;
plot(lag_max, r(max_idx), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
%text(lag_max, r(max_idx), ['  Peak: ' num2str(r(max_idx), '%.3f') ' at lag = ' num2str(lag_max, '%.2f') 's'], 'VerticalAlignment', 'bottom');

% Plot 3: Zoomed view of cross-correlation around zero lag
subplot(3,1,3);
% Define range to zoom (e.g., ±10 seconds around zero)
zoom_range = 10; % seconds
zoom_indices = find(lag_time >= -zoom_range & lag_time <= zoom_range);
plot(lag_time(zoom_indices), r(zoom_indices), 'k', 'LineWidth', 1.5);
title('Zoomed Cross-correlation (±10s around zero lag)');
xlabel('Lag (s)');
ylabel('Correlation Coefficient');
grid on;

% Add zero lag line for reference
hold on;
line([0 0], ylim, 'Color', 'r', 'LineStyle', '--');
text(0.2, mean(ylim), 'Zero lag', 'Color', 'r');



%%
% 1. Calculate the actual correlation
calcium_norm = (calcium_data_trimmed - mean(calcium_data_trimmed)) / std(calcium_data_trimmed);
pupil_norm = (pupil_data_trimmed - mean(pupil_data_trimmed)) / std(pupil_data_trimmed);

calcium_norm = movmean(calcium_norm, 10);
pupil_norm = movmean(pupil_norm, 10);

[r_actual, lags] = xcorr(calcium_norm, pupil_norm, 'coeff');
max_actual_corr = max(abs(r_actual));
lag_at_max = lags(find(abs(r_actual) == max_actual_corr, 1));

% 2. Perform permutation test
num_permutations = 10;  
max_shuffled_corrs = zeros(num_permutations, 1);

%shuffle pupil 
for i = 1:num_permutations
    % Shuffle the pupil data randomly
    pupil_shuffled = pupil_norm(randperm(length(pupil_norm)));
    
    % Calculate cross-correlation with shuffled data
    [r_shuffled, ~] = xcorr(calcium_norm, pupil_shuffled, 'coeff');
    
    % Store maximum correlation for this permutation
    max_shuffled_corrs(i) = max(abs(r_shuffled));
end

% 3. Calculate p-value
% Count how many permutations produced a correlation as high or higher than the actual
p_value = sum(max_shuffled_corrs >= max_actual_corr) / num_permutations;

% Find the maximum correlation and its lag
[max_actual_corr, max_idx] = max(abs(r_actual));
lag_at_max = lags(max_idx);
max_r_value = r_actual(max_idx);  % Get the actual correlation value (not just abs)

% 4. Visualize results
figure('Position', [100, 100, 1000, 800]);

subplot(2,1,1);
t_calcium = (0:length(calcium_norm)-1)/10; % Time in seconds
plot(t_calcium, calcium_norm, 'b', 'LineWidth', 1.5);
hold on;
plot(t_calcium, pupil_norm, 'r', 'LineWidth', 1.5);
legend('Calcium Signal (normalized)', 'Pupil Signal (normalized)');
xlabel('Time (s)');
% ylim([0 1]);
set(gca, 'Box', 'off');  


% Original cross-correlation plot
subplot(2,1,2);
lag_time = lags/10; % Convert lags to seconds
plot(lag_time, r_actual, 'k', 'LineWidth', 1.5);
title('Cross-correlation between Calcium and Pupil Signals');
xlabel('Lag (s)');
ylabel('Correlation Coefficient');
set(gca, 'Box', 'off');  
hold on;

% Create confidence bands at each lag point
shuffled_corrs_at_lags = zeros(num_permutations, length(lags));
for i = 1:num_permutations
    pupil_shuffled = pupil_norm(randperm(length(pupil_norm))); % Changed from randsample to randperm for consistency
    [shuffled_corrs_at_lags(i,:), ~] = xcorr(calcium_norm, pupil_shuffled, 'coeff');
end

plot(lag_time, shuffled_corrs_at_lags, 'Color', [.7 .7 .7], 'LineWidth', 1.5);


% Plot percentile bounds of shuffled correlations (e.g., 95% confidence interval)
upper_bound = prctile(shuffled_corrs_at_lags, 97.5, 1);
lower_bound = prctile(shuffled_corrs_at_lags, 2.5, 1);
% fill([lag_time, fliplr(lag_time)], [upper_bound, fliplr(lower_bound)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot the peak correlation point
plot(lag_at_max/10, max_r_value, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlim([-50, 50]);
legend('Original correlation', 'Shuffled signal', 'Peak correlation');

%%
% Calculate R-squared between calcium and pupil signals
calcium_norm = (calcium_data_trimmed - mean(calcium_data_trimmed)) / std(calcium_data_trimmed);
pupil_norm = (pupil_data_trimmed - mean(pupil_data_trimmed)) / std(pupil_data_trimmed);
pupil_norm = pupil_norm';

% Calculate actual R-squared values at different lags
max_lag = 50; % Adjust based on your data
[r_squared_actual, lags] = xcorr_r_squared(calcium_norm, pupil_norm, max_lag);
max_actual_r_squared = max(r_squared_actual);
lag_at_max = lags(find(r_squared_actual == max_actual_r_squared, 1));

% Perform permutation test
num_permutations = 5000;
max_shuffled_r_squared = zeros(num_permutations, 1);

for i = 1:num_permutations
    % Shuffle the pupil data randomly
    pupil_shuffled = pupil_norm(randperm(length(pupil_norm)));
    
    % Calculate R-squared with shuffled data
    [r_squared_shuffled, ~] = xcorr_r_squared(calcium_norm, pupil_shuffled, max_lag);
    
    % Store maximum R-squared for this permutation
    max_shuffled_r_squared(i) = max(r_squared_shuffled);
end

% Calculate p-value using the Laplace correction
p_value = (sum(max_shuffled_r_squared >= max_actual_r_squared) + 1) / (num_permutations + 1);

% Calculate Z-score
z_score = (max_actual_r_squared - mean(max_shuffled_r_squared)) / std(max_shuffled_r_squared);

% Visualize results
figure('Position', [100, 100, 1000, 800]);

% Plot R-squared values at different lags
subplot(2,1,1);
lag_time = lags/10; % Convert lags to seconds (assuming 10Hz sampling)
plot(lag_time, r_squared_actual, 'k', 'LineWidth', 1.5);
title(['R-squared between Calcium and Pupil Signals (p = ' num2str(p_value, '%.4e') ', Z = ' num2str(z_score, '%.2f') ')']);
xlabel('Lag (s)');
ylabel('R-squared');
grid on;
hold on;

% Plot the peak R-squared point
plot(lag_at_max/10, max_actual_r_squared, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
text(lag_at_max/10, max_actual_r_squared, ['  Max R² = ' num2str(max_actual_r_squared, '%.3f') ' at lag = ' num2str(lag_at_max/10, '%.2f') 's'], 'VerticalAlignment', 'bottom');

% Plot histogram of permutation results
subplot(2,1,2);
histogram(max_shuffled_r_squared, 30, 'Normalization', 'probability');
hold on;
xline(max_actual_r_squared, 'r', 'LineWidth', 2);
title(['Distribution of Maximum R-squared Under Null Hypothesis']);
xlabel('Maximum R-squared');
ylabel('Probability');

% Try fitting extreme value distribution for better p-value estimation
    pd = fitdist(max_shuffled_r_squared, 'GeneralizedExtremeValue');
    p_value_fitted = 1 - cdf(pd, max_actual_r_squared);
    fprintf('P-value from fitted distribution: %e\n', p_value_fitted);
    
    % Plot the fitted distribution
    x = linspace(min(max_shuffled_r_squared), max(max_actual_r_squared*1.1, max(max_shuffled_r_squared)), 1000);
    y = pdf(pd, x);
    plot(x, y, 'g-', 'LineWidth', 2);
    legend('Permutation results', 'Observed R-squared', 'GEV fit');
%%
% Load pupil area from your data
pupil_area = data.Behav.PupilArea;

% Create appropriate time vectors for both original and downsampled data
t_original = (1:length(pupil_area)) / 50;  
t_downsampled = (1:size(dFF,2)) / 10;      

% Downsample pupil size using interpolation
pupil_downsampled = interp1(t_original, pupil_area, t_downsampled, 'linear');

% Z-score the downsampled pupil data
pupil_z = (pupil_downsampled - mean(pupil_downsampled)) / std(pupil_downsampled);

% Check if z-scoring worked correctly
disp(['Mean of z-scored pupil: ' num2str(mean(pupil_z))]);
disp(['Std of z-scored pupil: ' num2str(std(pupil_z))]);

% Z-score the calcium data (dFF)
dFF_z = (dFF - mean(dFF, 2)) ./ std(dFF, 0, 2); 

% Plot 1: Original vs Downsampled (raw values)
figure;
subplot(2,1,1);
plot(t_original, pupil_area, 'r', 'LineWidth', 1); hold on;
plot(t_downsampled, pupil_downsampled, 'b', 'LineWidth', 2);
xlabel('Time (seconds)'); 
ylabel('Pupil Area (raw)');
legend('Original Pupil Area (60 Hz)', 'Downsampled Pupil Area (10 Hz)');
title('Raw Pupil Area: Original vs Downsampled');

% Plot 2: Z-scored downsampled pupil data (in separate figure to ensure it works)
figure;
plot(t_downsampled, pupil_z, 'b', 'LineWidth', 2);
xlabel('Time (seconds)'); 
ylabel('Pupil Area (z-scored)');
title('Z-scored Downsampled Pupil Area');
%%
% Load pupil area from your data
pupil_area = data.Behav.PupilArea;

% Create appropriate time vectors for both original and downsampled data
t_original = (1:length(pupil_area)) / 60;  % Time in seconds for 60 Hz data
t_downsampled = (1:size(dFF,2)) / 10;      % Time in seconds for 10 Hz data

% Downsample pupil size using interpolation
pupil_downsampled = interp1(t_original, pupil_area, t_downsampled, 'linear');

% Normalize (Z-score) if needed
pupil_z = (pupil_downsampled - mean(pupil_downsampled)) / std(pupil_downsampled); 
dFF_z = (dFF - mean(dFF, 2)) ./ std(dFF, 0, 2); 

% Plot both original and downsampled data
figure;
plot(t_original, pupil_area, 'k', 'LineWidth', 1); hold on;
plot(t_downsampled, pupil_downsampled, 'b', 'LineWidth', 1);
xlabel('Time (seconds)'); 
ylabel('Pupil Area');
legend('Original Pupil Area (60 Hz)', 'Downsampled Pupil Area (10 Hz)');
title('Comparison of Original and Downsampled Pupil Area');
hold off;
%% Extracting Timing

% Define stimulus time window
fs = 60;
num_trials = data.Stimuli.trials;
trial_duration = 30; 
stim_total_time = trial_duration * num_trials; 

start_time = data.Stimuli(1).TimeStimulusFrame(1);
end_time = start_time + stim_total_time;
[~, start_indx] = min(abs(data.TimeCa(2,:) - start_time));
[~, end_indx]   = min(abs(data.TimeCa(2,:) - end_time));
%%
dFF = data.CaData(1).Ca_dFF;
time_Ca = data.TimeCa(2,:);
dFF_norm = (dFF - min(dFF, [], 2)) ./ (max(dFF, [], 2) - min(dFF, [], 2));
roi_indices = 1:11;
num_trials = data.Stimuli.trials;
start_time = data.Stimuli(1).TimeStimulusFrame(1);
end_time = start_time + stim_total_time;
[~, start_indx] = min(abs(data.TimeCa(2,:) - start_time));
[~, end_indx]   = min(abs(data.TimeCa(2,:) - end_time));
    
dFF_norm_selected = dFF_norm(roi_indices,  start_indx:end_indx);
%% 
% *Plotting Heatmap with traces*


stim_total_time = data.Stimuli(2).stimulus_trial_t * data.Stimuli(2).trials;
stim_start = data.Stimuli(2).TimeStimulusFrame(1);
stim_end = stim_start + stim_total_time;

% Find frames of calcium for defined stimulus timing
Ca_index = find(data.TimeCa(1, :) > stim_start & data.TimeCa(1, :) < stim_end);
Ca_selected = data.CaData(1).Ca_dFF(:, Ca_index);

% Find behavior data for defined stimulus timing
Beh_index = find(data.Triggers.TimeCamera > stim_start & data.Triggers.TimeCamera < stim_end);
Beh_selected = data.Behav.PupilArea(Beh_index);

[~, sort_idx] = sort(max(Ca_selected, [], 2), 'descend'); % Change to 'ascend' if you want lowest to highest

Ca_selected = Ca_selected(sort_idx, :);

%dFF_norm_selected = (Ca_selected - mean(Ca_selected, 2)) ./ std(Ca_selected, 0, 2);
dFF_norm_selected = (Ca_selected - mean(Ca_selected, 2)) ./ std(Ca_selected, 0, 2)  .* max(Ca_selected,[],2);



% Create a figure with two vertically stacked subplots
figure;
subplot(4,1,1)
plot(linspace(stim_start, stim_end, length(Beh_selected)), Beh_selected, 'r');
axis off; 
xlim([stim_start, stim_end])
set(gca, 'Box', 'off');
set(gca, 'Position', [0.13 0.80 0.775 0.15]);

subplot(4,1,2)

Whisker_selected = data.Behav.whiskerFollicleA1(Beh_index, 1);
Whisker_selected_norm = (Whisker_selected - mean(Whisker_selected)) / std(Whisker_selected);
plot(linspace(stim_start, stim_end, length(Whisker_selected_norm)), Whisker_selected_norm, 'k')
axis off; 
xlim([stim_start, stim_end])
set(gca, 'Box', 'off');  
set(gca, 'Position', [0.13 0.62 0.775 0.15]);

subplot(4,1,3)

num_repeats = 5;  % Number of stimulus repetitions
fff_brightness_repeated = repmat(fff_brightness, 1, num_repeats); 
plot(linspace(stim_start, stim_end, length(fff_brightness_repeated)), fff_brightness_repeated, 'k', 'LineWidth', 1.5);
title ("Stimulus")
axis off; 
set(gca, 'Box', 'off');          
xlim([stim_start, stim_end])
        
subplot(4,1,4)

original_colormap = gray;  % Get the default gray colormap
reversed_colormap = flipud(original_colormap);
imagesc('XData',linspace(stim_start, stim_end, size(dFF_norm_selected, 2)), 'CData', dFF_norm_selected)
colormap(reversed_colormap);
caxis([0 3])
set(gca, 'YDir', 'reverse');
xlim([stim_start, stim_end])
xlabel('Time (seconds)'); 
ylabel('ROIs');



%% PUPIL-CALCIUM SETUP  (run this cell first)
% =========================================================================
%  Extracts calcium (10 Hz) and pupil (~50 Hz) for Stimulus 2 over all
%  trials, then resamples pupil onto the calcium time grid using real
%  timestamps so the two signals are properly aligned.
% =========================================================================

% Stimulus 2 time window
stim_total_time = data.Stimuli(2).stimulus_trial_t * data.Stimuli(2).trials;
stim_start_pc   = data.Stimuli(2).TimeStimulusFrame(1);
stim_end_pc     = stim_start_pc + stim_total_time;

% --- Calcium (10 Hz) ---
Ca_idx_pc  = find(data.TimeCa(1,:) > stim_start_pc & data.TimeCa(1,:) < stim_end_pc);
Ca_all_pc  = data.CaData(1).Ca_dFF(:, Ca_idx_pc);   % [nROIs × nCa]
Ca_avg_pc  = mean(Ca_all_pc, 1);                     % mean across ROIs: [1 × nCa]
Ca_time_pc = data.TimeCa(1, Ca_idx_pc);              % real timestamps at 10 Hz

% --- Pupil (behavioral camera, ~50 Hz) ---
Beh_idx_pc   = find(data.Triggers.TimeCamera > stim_start_pc & ...
                    data.Triggers.TimeCamera < stim_end_pc);
pupil_raw_pc = double(data.Behav.PupilArea(Beh_idx_pc));
pupil_time_pc = data.Triggers.TimeCamera(Beh_idx_pc);

% Resample pupil onto calcium timepoints using shared real timestamps
pupil_at_ca = interp1(pupil_time_pc(:), pupil_raw_pc(:), Ca_time_pc(:), ...
                      'linear', 'extrap');

% Remove NaNs — force column vectors so histogram2 / boxplot don't complain
valid_pc    = ~(isnan(pupil_at_ca) | isnan(Ca_avg_pc(:)));
ca_clean    = Ca_avg_pc(valid_pc);    ca_clean    = ca_clean(:);
pupil_clean = pupil_at_ca(valid_pc);  pupil_clean = pupil_clean(:);
t_clean     = Ca_time_pc(valid_pc);   t_clean     = t_clean(:);

% Correlations
[R_pc, P_pc]  = corrcoef(pupil_clean(:), ca_clean(:));
r_pearson_pc  = R_pc(1,2);
p_pearson_pc  = P_pc(1,2);
rho_spear_pc  = corr(pupil_clean(:), ca_clean(:), 'type', 'Spearman');

fprintf('Pearson r=%.3f (p=%.2e)  |  Spearman rho=%.3f\n', ...
        r_pearson_pc, p_pearson_pc, rho_spear_pc);

% --- Pupil as % of pre-stimulus baseline ---
% Uses 2 s before stimulus 2 onset as the reference.  Falls back to the
% within-window mean if no pre-stimulus pupil data is available.
pct_bsl_dur = 2;
bsl_beh_idx = find(data.Triggers.TimeCamera > (stim_start_pc - pct_bsl_dur) & ...
                   data.Triggers.TimeCamera < stim_start_pc);
if numel(bsl_beh_idx) > 1
    pupil_baseline = mean(double(data.Behav.PupilArea(bsl_beh_idx)), 'omitnan');
else
    pupil_baseline = mean(pupil_raw_pc, 'omitnan');
    fprintf('WARNING: no pre-stimulus pupil data — using window mean as baseline.\n');
end
fprintf('Pupil baseline = %.1f px²\n', pupil_baseline);

pupil_pct_pc    = (pupil_at_ca / pupil_baseline) * 100;   % % of baseline, [nCa × 1]
pupil_pct_clean = pupil_pct_pc(valid_pc);
pupil_pct_clean = pupil_pct_clean(:);

% Pupil bins (reused by plots below)
n_bins_pc    = 5;
bin_edges_pc = linspace(min(pupil_clean), max(pupil_clean), n_bins_pc + 1);
bin_ctr_pc   = (bin_edges_pc(1:end-1) + bin_edges_pc(2:end)) / 2;

ca_by_bin_pc = cell(n_bins_pc, 1);
bin_n_pc     = zeros(n_bins_pc, 1);
bin_mean_pc  = zeros(n_bins_pc, 1);
bin_sem_pc   = zeros(n_bins_pc, 1);

for ib = 1:n_bins_pc
    if ib < n_bins_pc
        mask = pupil_clean >= bin_edges_pc(ib) & pupil_clean < bin_edges_pc(ib+1);
    else
        mask = pupil_clean >= bin_edges_pc(ib) & pupil_clean <= bin_edges_pc(ib+1);
    end
    ca_by_bin_pc{ib} = ca_clean(mask);
    bin_n_pc(ib)     = sum(mask);
    if bin_n_pc(ib) > 1
        bin_mean_pc(ib) = mean(ca_by_bin_pc{ib});
        bin_sem_pc(ib)  = std(ca_by_bin_pc{ib}) / sqrt(bin_n_pc(ib));
    end
end

%% PLOT 1: Scatter — Pupil Area vs Mean Calcium
% Each point is one calcium timepoint. Regression line shows overall trend.

figure('Name', 'P1: Scatter Pupil vs Calcium', 'Position', [100 100 650 560]);

scatter(pupil_clean, ca_clean, 8, [0.2 0.4 0.8], 'filled', 'MarkerFaceAlpha', 0.25);
hold on;

p1_coeff = polyfit(pupil_clean, ca_clean, 1);
x1_fit   = linspace(min(pupil_clean), max(pupil_clean), 200);
plot(x1_fit, polyval(p1_coeff, x1_fit), 'r-', 'LineWidth', 2);

xlabel('Pupil Area (px²)',    'FontSize', 12);
ylabel('Mean Calcium (ΔF/F)', 'FontSize', 12);
title('Scatter: Pupil Area vs Calcium Activity', 'FontSize', 13);
text(0.05, 0.93, sprintf('r = %.3f  (p = %.2e)', r_pearson_pc, p_pearson_pc), ...
    'Units', 'normalized', 'FontSize', 11, 'BackgroundColor', 'w', 'EdgeColor', 'k');
legend({'Data', 'Linear fit'}, 'Location', 'best');
set(gca, 'Box', 'off');

%% PLOT 2: Mean ± SEM Calcium by Pupil Size Bins
% Bins pupil into 5 equal-width ranges; shows mean calcium in each bin.

figure('Name', 'P2: Mean Calcium by Pupil Bin', 'Position', [100 100 580 480]);

errorbar(bin_ctr_pc, bin_mean_pc, bin_sem_pc, 'o-', ...
    'LineWidth', 2, 'MarkerSize', 8, ...
    'Color', [0.2 0.6 0.9], 'MarkerFaceColor', [0.2 0.6 0.9]);
hold on;

for ib = 1:n_bins_pc
    if bin_n_pc(ib) > 0
        text(bin_ctr_pc(ib), bin_mean_pc(ib) + bin_sem_pc(ib) + 0.05, ...
            sprintf('n=%d', bin_n_pc(ib)), ...
            'HorizontalAlignment', 'center', 'FontSize', 9);
    end
end

xlabel('Pupil Area (px²) — bin centre', 'FontSize', 12);
ylabel('Mean Calcium (ΔF/F)',            'FontSize', 12);
title('Mean ± SEM Calcium by Pupil Size', 'FontSize', 13);
set(gca, 'Box', 'off');  grid on;

%% PLOT 3: 2D Density Heatmap — Joint distribution
% Colour = number of timepoints falling in each (pupil, calcium) cell.

figure('Name', 'P3: 2D Density Heatmap', 'Position', [100 100 680 540]);

histogram2(pupil_clean, ca_clean, [20 20], 'DisplayStyle', 'tile', 'EdgeColor', 'none');
colormap('hot');  colorbar;
xlabel('Pupil Area (px²)',  'FontSize', 12);
ylabel('Calcium (ΔF/F)',    'FontSize', 12);
title('2D Density: Joint Distribution', 'FontSize', 13);
set(gca, 'Box', 'off');

%% PLOT 4: Time Series — Calcium (black) and Pupil (blue) on dual axes
% Lets you see whether calcium and pupil co-vary over the recording.

figure('Name', 'P4: Time Series Calcium + Pupil', 'Position', [100 100 1100 380]);

yyaxis left
plot(t_clean, ca_clean, 'k-', 'LineWidth', 1.2);
ylabel('Mean Calcium (ΔF/F)', 'FontSize', 12);

yyaxis right
plot(t_clean, pupil_clean, '-', 'Color', [0.1 0.45 0.85], 'LineWidth', 1.2);
ylabel('Pupil Area (px²)', 'FontSize', 12);

xlabel('Time (s)', 'FontSize', 12);
title('Calcium (black) and Pupil Size (blue) Over Time', 'FontSize', 13);
set(gca, 'Box', 'off');

%% PLOT 5: Box Plots — Calcium distribution per pupil size bin
% Each box shows median, IQR, and outliers for one pupil size range.

figure('Name', 'P5: Box Plots by Pupil Bin', 'Position', [100 100 680 520]);

grp_data   = [];
grp_labels = [];
for ib = 1:n_bins_pc
    if bin_n_pc(ib) > 0
        grp_data   = [grp_data;   ca_by_bin_pc{ib}(:)];          %#ok<AGROW>
        grp_labels = [grp_labels; repmat(ib, bin_n_pc(ib), 1)];  %#ok<AGROW>
    end
end

tick_labels = arrayfun(@(x) sprintf('%.0f', x), bin_ctr_pc, 'UniformOutput', false);
boxplot(grp_data, grp_labels, 'Labels', tick_labels, 'Colors', 'b', 'Symbol', '+');

xlabel('Pupil Size Bin Centre (px²)', 'FontSize', 12);
ylabel('Calcium (ΔF/F)',              'FontSize', 12);
title('Calcium Distribution by Pupil Size', 'FontSize', 13);
set(gca, 'Box', 'off');  grid on;

%% PLOT 6: Scatter — Pupil % Baseline vs Mean Calcium
% Same as Plot 1 but x-axis shows pupil relative to its pre-stimulus baseline.
% 100% = baseline level; >100% = dilated; <100% = constricted.

[R_pct, P_pct] = corrcoef(pupil_pct_clean, ca_clean);
r_pct = R_pct(1,2);  p_pct = P_pct(1,2);

figure('Name', 'P6: Scatter Pupil % Baseline vs Calcium', 'Position', [100 100 650 560]);
scatter(pupil_pct_clean, ca_clean, 8, [0.2 0.4 0.8], 'filled', 'MarkerFaceAlpha', 0.25);
hold on;
pf6 = polyfit(pupil_pct_clean, ca_clean, 1);
xf6 = linspace(min(pupil_pct_clean), max(pupil_pct_clean), 200);
plot(xf6, polyval(pf6, xf6), 'r-', 'LineWidth', 2);
xline(100, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);   % baseline level
xlabel('Pupil Size (% of baseline)', 'FontSize', 12);
ylabel('Mean Calcium (ΔF/F)',         'FontSize', 12);
title(sprintf('Scatter: Pupil %% Baseline vs Calcium  (baseline = %.0f px²)', pupil_baseline), ...
      'FontSize', 13);
text(0.05, 0.93, sprintf('r = %.3f  (p = %.2e)', r_pct, p_pct), ...
    'Units', 'normalized', 'FontSize', 11, 'BackgroundColor', 'w', 'EdgeColor', 'k');
legend({'Data', 'Linear fit', 'Baseline (100%)'}, 'Location', 'best');
set(gca, 'Box', 'off');

%% PLOT 7: Per-Axon Correlation — Pupil % Baseline vs Individual ROI Calcium
% Computes Pearson r between pupil % and each axon bouton (ROI) individually.
% Left panel: ROIs in anatomical order.  Right panel: sorted by r-value.
% Orange bars = positive correlation; blue = negative.  * p<0.05, ** p<0.01, *** p<0.001.

nROIs_pc  = size(Ca_all_pc, 1);
r_roi_pct = NaN(nROIs_pc, 1);
p_roi_pct = NaN(nROIs_pc, 1);

for roi = 1:nROIs_pc
    ca_roi    = Ca_all_pc(roi, :)';
    valid_roi = valid_pc(:) & ~isnan(ca_roi);
    if sum(valid_roi) > 5
        [R_r, P_r]    = corrcoef(pupil_pct_pc(valid_roi), ca_roi(valid_roi));
        r_roi_pct(roi) = R_r(1,2);
        p_roi_pct(roi) = P_r(1,2);
    end
end

n_sig_pct = sum(p_roi_pct < 0.05 & ~isnan(p_roi_pct));
fprintf('Pupil %% vs Ca:  %d / %d ROIs significant (p<0.05)\n', n_sig_pct, nROIs_pc);

figure('Name', 'P7: Per-ROI Correlation Pupil % vs Ca', 'Position', [100 100 960 480]);

for panel = 1:2
    subplot(1, 2, panel);  hold on;

    if panel == 1
        r_plot = r_roi_pct;
        p_plot = p_roi_pct;
        x_vals = 1:nROIs_pc;
        xlabel('ROI (axon bouton index)', 'FontSize', 11);
        title('Per-ROI r  —  anatomical order', 'FontSize', 12);
    else
        [r_plot, sort_idx] = sort(r_roi_pct, 'descend');
        p_plot = p_roi_pct(sort_idx);
        x_vals = 1:nROIs_pc;
        xlabel('Rank (sorted by r)', 'FontSize', 11);
        title(sprintf('Sorted  —  %d significant (p<0.05)', n_sig_pct), 'FontSize', 12);
    end

    for k = 1:nROIs_pc
        if ~isnan(r_plot(k))
            if r_plot(k) >= 0
                clr = [0.85 0.35 0.1];
            else
                clr = [0.1 0.45 0.85];
            end
            bar(x_vals(k), r_plot(k), 0.7, 'FaceColor', clr, 'EdgeColor', 'none');
            if p_plot(k) < 0.05
                if p_plot(k) < 0.001,    mk = '***';
                elseif p_plot(k) < 0.01, mk = '**';
                else,                    mk = '*';
                end
                text(x_vals(k), r_plot(k) + sign(r_plot(k)) * 0.04, mk, ...
                     'HorizontalAlignment', 'center', 'FontSize', 7, 'FontWeight', 'bold');
            end
        end
    end

    yline(0, 'k-', 'LineWidth', 0.8);
    ylabel('Pearson r  (pupil % vs Ca dF/F)', 'FontSize', 11);
    ylim([-1, 1]);  xlim([0, nROIs_pc + 1]);
    set(gca, 'Box', 'off');
end

%% =========================================================================
%  PER-STIMULUS PUPIL ANALYSIS  (Plots 4x–7x)
%  Repeats Plots 4–7 across all stimuli.
%  Requires: trial-chunking cell (~line 228) run first.
%  Run PER-STIM PUPIL SETUP, then any plot cell independently.
%  =========================================================================

%% PER-STIM PUPIL SETUP — extract full block for every stimulus

pct_bsl_dur_ps = 2;   % seconds of pre-stimulus baseline for pupil %
n_bins_ps      = 5;   % pupil bins used in Plot 5x
num_stim_ps    = length(data.Stimuli);
pup_stim       = cell(num_stim_ps, 1);

for si = 1:num_stim_ps
    blk_dur_si = data.Stimuli(si).stimulus_trial_t * data.Stimuli(si).trials;
    stim_blk_s = data.Stimuli(si).TimeStimulusFrame(1);
    stim_blk_e = stim_blk_s + blk_dur_si;
    lbl_si     = sprintf('%s | %s', data.Stimuli(si).type, data.Stimuli(si).addnote);

    % Calcium full block
    Ca_idx_ps = find(data.TimeCa(1,:) > stim_blk_s & data.TimeCa(1,:) < stim_blk_e);
    t_ca_ps   = data.TimeCa(1, Ca_idx_ps) - stim_blk_s;     % row, block-relative
    Ca_all_ps = data.CaData(1).Ca_dFF(:, Ca_idx_ps);         % [nROIs × nCa]
    Ca_avg_ps = mean(Ca_all_ps, 1, 'omitnan');                % [1 × nCa]

    % Pupil full block, resampled to Ca grid
    Beh_idx_ps  = find(data.Triggers.TimeCamera > stim_blk_s & data.Triggers.TimeCamera < stim_blk_e);
    t_beh_ps    = data.Triggers.TimeCamera(Beh_idx_ps) - stim_blk_s;
    pup_raw_ps  = double(data.Behav.PupilArea(Beh_idx_ps));
    pup_atca_ps = interp1(t_beh_ps(:), pup_raw_ps(:), t_ca_ps(:), 'linear', 'extrap');  % col

    % Pupil % of pre-stimulus baseline
    bsl_beh_ps = find(data.Triggers.TimeCamera > (stim_blk_s - pct_bsl_dur_ps) & ...
                      data.Triggers.TimeCamera < stim_blk_s);
    if numel(bsl_beh_ps) > 1
        pup_bsl_ps = mean(double(data.Behav.PupilArea(bsl_beh_ps)), 'omitnan');
    else
        pup_bsl_ps = mean(pup_raw_ps, 'omitnan');
    end
    pup_pct_ps = (pup_atca_ps / pup_bsl_ps) * 100;   % col, [nCa × 1]

    % NaN mask — force column vectors
    valid_ps  = ~(isnan(pup_atca_ps) | isnan(Ca_avg_ps(:)));
    ca_v_ps   = Ca_avg_ps(valid_ps)';   ca_v_ps   = ca_v_ps(:);
    pup_v_ps  = pup_atca_ps(valid_ps);  pup_v_ps  = pup_v_ps(:);
    pct_v_ps  = pup_pct_ps(valid_ps);   pct_v_ps  = pct_v_ps(:);
    t_v_ps    = t_ca_ps(valid_ps)';     t_v_ps    = t_v_ps(:);

    % Pupil bins (Plot 5x)
    bin_edg_ps  = linspace(min(pup_v_ps), max(pup_v_ps), n_bins_ps + 1);
    bin_ctr_ps  = (bin_edg_ps(1:end-1) + bin_edg_ps(2:end)) / 2;
    bin_mean_ps = zeros(n_bins_ps, 1);
    bin_sem_ps  = zeros(n_bins_ps, 1);
    bin_n_ps    = zeros(n_bins_ps, 1);
    bin_dat_ps  = cell(n_bins_ps, 1);
    for ib = 1:n_bins_ps
        if ib < n_bins_ps
            msk_ps = pup_v_ps >= bin_edg_ps(ib) & pup_v_ps < bin_edg_ps(ib+1);
        else
            msk_ps = pup_v_ps >= bin_edg_ps(ib) & pup_v_ps <= bin_edg_ps(ib+1);
        end
        bin_dat_ps{ib}  = ca_v_ps(msk_ps);
        bin_n_ps(ib)    = sum(msk_ps);
        if bin_n_ps(ib) > 1
            bin_mean_ps(ib) = mean(bin_dat_ps{ib});
            bin_sem_ps(ib)  = std(bin_dat_ps{ib}) / sqrt(bin_n_ps(ib));
        end
    end

    % Per-ROI Pearson r with pupil %
    nROIs_ps = size(Ca_all_ps, 1);
    r_roi_ps = NaN(nROIs_ps, 1);
    p_roi_ps = NaN(nROIs_ps, 1);
    for roi = 1:nROIs_ps
        ca_roi_ps = Ca_all_ps(roi, :)';
        v_roi_ps  = valid_ps(:) & ~isnan(ca_roi_ps);
        if sum(v_roi_ps) > 5
            [Rr_ps, Pr_ps]  = corrcoef(pup_pct_ps(v_roi_ps), ca_roi_ps(v_roi_ps));
            r_roi_ps(roi)   = Rr_ps(1,2);
            p_roi_ps(roi)   = Pr_ps(1,2);
        end
    end

    pup_stim{si}.label       = lbl_si;
    pup_stim{si}.t_v         = t_v_ps;
    pup_stim{si}.ca_v        = ca_v_ps;
    pup_stim{si}.pup_v       = pup_v_ps;
    pup_stim{si}.pct_v       = pct_v_ps;
    pup_stim{si}.blk_dur     = blk_dur_si;
    pup_stim{si}.pup_bsl     = pup_bsl_ps;
    pup_stim{si}.bin_ctr     = bin_ctr_ps;
    pup_stim{si}.bin_mean    = bin_mean_ps;
    pup_stim{si}.bin_sem     = bin_sem_ps;
    pup_stim{si}.bin_n       = bin_n_ps;
    pup_stim{si}.bin_dat     = bin_dat_ps;
    pup_stim{si}.r_roi       = r_roi_ps;
    pup_stim{si}.p_roi       = p_roi_ps;
    pup_stim{si}.nROIs       = nROIs_ps;
end
fprintf('Per-stim pupil setup done: %d stimuli.\n', num_stim_ps);

%% PLOT 4x: Per-Stimulus — Time Series Calcium + Pupil (dual axis)

for si = 1:num_stim_ps
    t_p   = pup_stim{si}.t_v;
    ca_p  = pup_stim{si}.ca_v;
    pup_p = pup_stim{si}.pup_v;
    blk_d = pup_stim{si}.blk_dur;

    figure('Name', sprintf('P4x: Stim %d — Time Series  [%s]', si, pup_stim{si}.label), ...
           'Position', [100 100 1100 380]);
    yyaxis left
    plot(t_p, ca_p, 'k-', 'LineWidth', 1.2);
    ylabel('Mean Calcium (ΔF/F)', 'FontSize', 12);
    yyaxis right
    plot(t_p, pup_p, '-', 'Color', [0.1 0.45 0.85], 'LineWidth', 1.2);
    ylabel('Pupil Area (px²)', 'FontSize', 12);
    xlabel('Time within block (s)', 'FontSize', 12);
    title(sprintf('Stim %d — Calcium and Pupil Over Time  [%s]', si, pup_stim{si}.label), ...
          'FontSize', 12, 'Interpreter', 'none');
    xlim([0, blk_d]);  set(gca, 'Box', 'off');
end

%% PLOT 5x: Per-Stimulus — Calcium Distribution by Pupil Bin (boxplot)

for si = 1:num_stim_ps
    b_dat = pup_stim{si}.bin_dat;
    b_ctr = pup_stim{si}.bin_ctr;
    b_n   = pup_stim{si}.bin_n;

    grp_d = [];  grp_l = [];
    for ib = 1:n_bins_ps
        if b_n(ib) > 0
            grp_d = [grp_d; b_dat{ib}(:)];               %#ok<AGROW>
            grp_l = [grp_l; repmat(ib, b_n(ib), 1)];     %#ok<AGROW>
        end
    end

    tick_lbl = arrayfun(@(x) sprintf('%.0f', x), b_ctr, 'UniformOutput', false);

    figure('Name', sprintf('P5x: Stim %d — Ca by Pupil Bin  [%s]', si, pup_stim{si}.label), ...
           'Position', [100 100 680 520]);
    boxplot(grp_d, grp_l, 'Labels', tick_lbl, 'Colors', 'b', 'Symbol', '+');
    xlabel('Pupil Size Bin Centre (px²)', 'FontSize', 12);
    ylabel('Calcium (ΔF/F)',              'FontSize', 12);
    title(sprintf('Stim %d — Ca Distribution by Pupil Size  [%s]', si, pup_stim{si}.label), ...
          'FontSize', 12, 'Interpreter', 'none');
    set(gca, 'Box', 'off');  grid on;
end

%% PLOT 6x: Per-Stimulus — Scatter Pupil % Baseline vs Mean Calcium

for si = 1:num_stim_ps
    pct_p   = pup_stim{si}.pct_v;
    ca_p    = pup_stim{si}.ca_v;
    bsl_p   = pup_stim{si}.pup_bsl;

    [R6x, P6x] = corrcoef(pct_p, ca_p);
    r6x = R6x(1,2);  p6x = P6x(1,2);

    figure('Name', sprintf('P6x: Stim %d — Scatter Pupil %% vs Ca  [%s]', si, pup_stim{si}.label), ...
           'Position', [100 100 650 560]);
    scatter(pct_p, ca_p, 8, [0.2 0.4 0.8], 'filled', 'MarkerFaceAlpha', 0.25);  hold on;
    pf6x = polyfit(pct_p, ca_p, 1);
    xf6x = linspace(min(pct_p), max(pct_p), 200);
    plot(xf6x, polyval(pf6x, xf6x), 'r-', 'LineWidth', 2);
    xline(100, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
    xlabel('Pupil Size (% of baseline)', 'FontSize', 12);
    ylabel('Mean Calcium (ΔF/F)',         'FontSize', 12);
    title(sprintf('Stim %d — Pupil %% Baseline vs Ca  (bsl = %.0f px²)', si, bsl_p), ...
          'FontSize', 12, 'Interpreter', 'none');
    text(0.05, 0.93, sprintf('r = %.3f  (p = %.2e)', r6x, p6x), ...
        'Units', 'normalized', 'FontSize', 11, 'BackgroundColor', 'w', 'EdgeColor', 'k');
    legend({'Data', 'Linear fit', 'Baseline (100%)'}, 'Location', 'best');
    set(gca, 'Box', 'off');
end

%% PLOT 7x: Per-Stimulus — Per-ROI Correlation Pupil % vs Axon Ca
% Left: ROIs in anatomical order.  Right: sorted by r-value.
% Orange = positive correlation, blue = negative.

for si = 1:num_stim_ps
    r_p  = pup_stim{si}.r_roi;
    p_p  = pup_stim{si}.p_roi;
    nR   = pup_stim{si}.nROIs;
    nsig = sum(p_p < 0.05 & ~isnan(p_p));

    figure('Name', sprintf('P7x: Stim %d — Per-ROI Correlation  [%s]', si, pup_stim{si}.label), ...
           'Position', [100 100 960 480]);

    for panel = 1:2
        subplot(1, 2, panel);  hold on;

        if panel == 1
            r_pl = r_p;  p_pl = p_p;
            xlabel('ROI (axon bouton index)', 'FontSize', 11);
            title(sprintf('Stim %d — anatomical order', si), 'FontSize', 12);
        else
            [r_pl, srt] = sort(r_p, 'descend');
            p_pl = p_p(srt);
            xlabel('Rank (sorted by r)', 'FontSize', 11);
            title(sprintf('Sorted  —  %d / %d significant (p<0.05)', nsig, nR), 'FontSize', 12);
        end

        for k = 1:nR
            if ~isnan(r_pl(k))
                clr = [0.85 0.35 0.1] * (r_pl(k) >= 0) + [0.1 0.45 0.85] * (r_pl(k) < 0);
                bar(k, r_pl(k), 0.7, 'FaceColor', clr, 'EdgeColor', 'none');
                if ~isnan(p_pl(k)) && p_pl(k) < 0.05
                    if p_pl(k) < 0.001,    mk = '***';
                    elseif p_pl(k) < 0.01, mk = '**';
                    else,                  mk = '*';
                    end
                    text(k, r_pl(k) + sign(r_pl(k)) * 0.04, mk, ...
                         'HorizontalAlignment', 'center', 'FontSize', 7, 'FontWeight', 'bold');
                end
            end
        end

        yline(0, 'k-', 'LineWidth', 0.8);
        ylabel('Pearson r  (pupil % vs Ca dF/F)', 'FontSize', 11);
        ylim([-1, 1]);  xlim([0, nR + 1]);
        set(gca, 'Box', 'off');
    end
end

%% =========================================================================
%  PUPIL-CALCIUM EXTENDED ANALYSIS — ALL STIMULI
%  Groups stimuli by data.Stimuli.type + data.Stimuli.addnote so only
%  truly comparable stimuli are averaged together in Cell D.
%  REQUIRES: trial-chunking cell (~line 228) run first (populates TrialTimes).
%  Run EXTENDED SETUP first, then any plot cell independently.
%  =========================================================================

%% EXTENDED SETUP — extract all stimuli, group by type + addnote

% ---- ROI SELECTION FLAG — edit here, then re-run any plot cell ----
roi_mode = 'all';   % 'all' = mean across all ROIs  |  'single' = one ROI only
roi_id   = 4;       % only used when roi_mode = 'single'
% Function handle lets plot cells check the flag without triggering
% MATLAB's "unreachable code" static-analysis warning on the else branches.
use_single = @() strcmp(roi_mode, 'single');
% -------------------------------------------------------------------

if ~isfield(data.Stimuli(1), 'TrialTimes') || isempty(data.Stimuli(1).TrialTimes)
    error('TrialTimes missing — run the trial-chunking cell (~line 228) first.');
end

num_stim_x   = length(data.Stimuli);
num_trials_x = data.Stimuli(1).trials;
n_com_x      = 500;   % resampling points per trial (increase for smoother curves)
nROIs_x      = size(data.CaData(1).Ca_dFF, 1);

if use_single() && (roi_id < 1 || roi_id > nROIs_x)
    error('roi_id=%d is out of range — available ROIs: 1 to %d.', roi_id, nROIs_x);
end

% Group stimuli by type + addnote — only stimuli with identical BOTH fields are averaged
stim_types_x    = {data.Stimuli.type};
stim_addnotes_x = {data.Stimuli.addnote};
stim_labels_x   = cellfun(@(t,a) sprintf('%s | %s', t, a), ...
                           stim_types_x, stim_addnotes_x, 'UniformOutput', false);
[uniq_labels_x, ~, label_grp_x] = unique(stim_labels_x);

fprintf('\nStimulus groups by type + addnote:\n');
for g = 1:length(uniq_labels_x)
    fprintf('  Group %d  [%s]  →  Stim %s\n', g, uniq_labels_x{g}, ...
            num2str(find(label_grp_x == g)'));
end
fprintf('ROI mode: %s', roi_mode);
if use_single(),  fprintf(' (ROI %d)', roi_id);  end
fprintf('\n');

% Extract per-trial Ca and pupil for every stimulus onto a trial-relative grid
stim_data = cell(num_stim_x, 1);

for si = 1:num_stim_x
    trial_dur    = data.Stimuli(si).stimulus_trial_t;
    t_com        = linspace(0, trial_dur, n_com_x);  % [1 × n_com_x], row, trial-relative

    Ca_trials_x  = NaN(nROIs_x, n_com_x, num_trials_x);
    pup_trials_x = NaN(n_com_x, num_trials_x);

    for tr = 1:num_trials_x
        stim_s = data.Stimuli(si).TrialTimes{1, tr}(1);
        stim_e = stim_s + trial_dur;

        % Calcium (10 Hz) — shift to trial-relative time
        Ca_idx = find(data.TimeCa(1,:) > stim_s & data.TimeCa(1,:) < stim_e);
        if length(Ca_idx) > 1
            t_ca   = data.TimeCa(1, Ca_idx) - stim_s;
            Ca_raw = data.CaData(1).Ca_dFF(:, Ca_idx);   % [nROIs × nCa]
            for roi = 1:nROIs_x
                Ca_trials_x(roi, :, tr) = interp1(t_ca, Ca_raw(roi,:), t_com, 'linear', 'extrap');
            end
        end

        % Pupil (~50 Hz) — shift to trial-relative time, resample onto calcium grid
        Beh_idx = find(data.Triggers.TimeCamera > stim_s & data.Triggers.TimeCamera < stim_e);
        if length(Beh_idx) > 1
            t_pup   = data.Triggers.TimeCamera(Beh_idx) - stim_s;
            pup_raw = double(data.Behav.PupilArea(Beh_idx));
            pup_trials_x(:, tr) = interp1(t_pup(:), pup_raw(:), t_com(:), 'linear', 'extrap');
        end
    end

    stim_data{si}.Ca_trials  = Ca_trials_x;        % [nROIs × n_com × nTrials]
    stim_data{si}.pup_trials = pup_trials_x;        % [n_com × nTrials]
    stim_data{si}.t_com      = t_com;               % [1 × n_com], trial-relative
    stim_data{si}.trial_dur  = trial_dur;
    stim_data{si}.label      = stim_labels_x{si};   % 'type | addnote' string
    stim_data{si}.label_grp  = label_grp_x(si);
end

fprintf('Done: %d stimuli × %d trials × %d ROIs × %d timepoints.\n', ...
        num_stim_x, num_trials_x, nROIs_x, n_com_x);

%% CELL A: Per-Stimulus — Individual Trial Traces + Mean
% roi_mode='all':    grey = per-ROI trial mean (nROIs lines), black = grand mean
% roi_mode='single': grey = individual trials for roi_id, black = trial mean
% Bottom panel: individual pupil trials (pink) + mean ± SEM (red).

for si = 1:num_stim_x
    t_com   = stim_data{si}.t_com;
    Ca_all  = stim_data{si}.Ca_trials;
    pup_all = stim_data{si}.pup_trials;

    pup_mean_A = mean(pup_all, 2, 'omitnan')';
    pup_sem_A  = (std(pup_all, 0, 2, 'omitnan') / sqrt(num_trials_x))';

    fig_name_A = sprintf('A: Stim %d — Traces  [%s]', si, stim_data{si}.label);
    figure('Name', fig_name_A, 'Position', [100 100 900 500]);

    subplot(2, 1, 1);  hold on;
    if ~use_single()
        Ca_roi_mean = mean(Ca_all, 3, 'omitnan');   % [nROIs × n_com]
        Ca_grand    = mean(Ca_roi_mean, 1);
        for roi = 1:nROIs_x
            plot(t_com, Ca_roi_mean(roi,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
        end
        plot(t_com, Ca_grand, 'k-', 'LineWidth', 2.5);
        title(sprintf('Stim %d — All ROI means (grey) + grand mean (black)', si), ...
              'FontSize', 11, 'Interpreter', 'none');
    else
        Ca_sel   = squeeze(Ca_all(roi_id, :, :));   % [n_com × nTrials]
        Ca_mean_A = mean(Ca_sel, 2)';
        for tr = 1:num_trials_x
            plot(t_com, Ca_sel(:,tr)', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
        end
        plot(t_com, Ca_mean_A, 'k-', 'LineWidth', 2.5);
        title(sprintf('Stim %d — ROI %d: individual trials (grey) + mean (black)', si, roi_id), ...
              'FontSize', 11, 'Interpreter', 'none');
    end
    ylabel('Calcium (ΔF/F)', 'FontSize', 11);
    xlim([0, stim_data{si}.trial_dur]);  set(gca, 'Box', 'off');

    subplot(2, 1, 2);  hold on;
    for tr = 1:num_trials_x
        plot(t_com, pup_all(:,tr)', 'Color', [1 0.7 0.7], 'LineWidth', 0.8);
    end
    x_sh_A = [t_com, fliplr(t_com)];
    y_sh_A = [(pup_mean_A + pup_sem_A), fliplr(pup_mean_A - pup_sem_A)];
    fill(x_sh_A, y_sh_A, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(t_com, pup_mean_A, 'r-', 'LineWidth', 2.5);
    xlabel('Time within trial (s)', 'FontSize', 11);
    ylabel('Pupil Area (px²)', 'FontSize', 11);
    xlim([0, stim_data{si}.trial_dur]);  set(gca, 'Box', 'off');
end

%% CELL B: Per-Stimulus — Scatter (Pupil vs Ca) + XCorr
% All trials concatenated per stimulus.  XCorr limited to ±10 s lag.

for si = 1:num_stim_x
    Ca_all  = stim_data{si}.Ca_trials;
    pup_all = stim_data{si}.pup_trials;
    fs_si   = n_com_x / stim_data{si}.trial_dur;

    % Apply ROI flag → Ca_sel: [n_com × nTrials]
    if use_single()
        Ca_sel = squeeze(Ca_all(roi_id, :, :));
    else
        Ca_sel = squeeze(mean(Ca_all, 1, 'omitnan'));
    end

    ca_B  = Ca_sel(:);
    pup_B = pup_all(:);
    valid_B = ~(isnan(ca_B) | isnan(pup_B));
    ca_B  = ca_B(valid_B);
    pup_B = pup_B(valid_B);

    [R_B, P_B] = corrcoef(pup_B, ca_B);
    r_B = R_B(1,2);  p_B = P_B(1,2);

    ca_z_B  = (ca_B  - mean(ca_B))  / std(ca_B);
    pup_z_B = (pup_B - mean(pup_B)) / std(pup_B);
    [xc_B, lags_B] = xcorr(ca_z_B, pup_z_B, round(fs_si * 10), 'coeff');
    lag_s_B = lags_B / fs_si;
    [~, pk_B] = max(abs(xc_B));

    figure('Name', sprintf('B: Stim %d — Scatter + XCorr  [%s]', si, stim_data{si}.label), ...
           'Position', [100 100 950 430]);

    subplot(1, 2, 1);
    scatter(pup_B, ca_B, 5, [0.2 0.4 0.8], 'filled', 'MarkerFaceAlpha', 0.15);
    hold on;
    pf_B = polyfit(pup_B, ca_B, 1);
    xf_B = linspace(min(pup_B), max(pup_B), 200);
    plot(xf_B, polyval(pf_B, xf_B), 'r-', 'LineWidth', 2);
    text(0.05, 0.93, sprintf('r = %.3f   p = %.2e', r_B, p_B), ...
        'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k');
    xlabel('Pupil Area (px²)',  'FontSize', 11);
    ylabel('Calcium (ΔF/F)',    'FontSize', 11);
    title(sprintf('Stim %d — Scatter', si), 'FontSize', 12);
    legend({'Data', 'Linear fit'}, 'Location', 'best');
    set(gca, 'Box', 'off');

    subplot(1, 2, 2);
    plot(lag_s_B, xc_B, 'k-', 'LineWidth', 1.5);  hold on;
    xline(0, 'r--', 'LineWidth', 1);
    plot(lag_s_B(pk_B), xc_B(pk_B), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    xlabel('Lag (s)',                 'FontSize', 11);
    ylabel('Correlation Coefficient', 'FontSize', 11);
    title(sprintf('Stim %d — XCorr  (peak at %.2f s)', si, lag_s_B(pk_B)), 'FontSize', 12);
    set(gca, 'Box', 'off');  grid on;
end

%% CELL C1: Per-Stimulus — Dual-Axis Time Trace (Ca + Pupil)
% Trial-averaged Ca (black, left axis) and pupil (blue, right axis).

for si = 1:num_stim_x
    t_com   = stim_data{si}.t_com;
    Ca_all  = stim_data{si}.Ca_trials;
    pup_all = stim_data{si}.pup_trials;

    if use_single()
        Ca_sel = squeeze(Ca_all(roi_id, :, :));
    else
        Ca_sel = squeeze(mean(Ca_all, 1, 'omitnan'));
    end
    Ca_avg_C1  = mean(Ca_sel, 2)';
    pup_avg_C1 = mean(pup_all, 2, 'omitnan')';

    figure('Name', sprintf('C1: Stim %d — Dual-Axis Trace  [%s]', si, stim_data{si}.label), ...
           'Position', [100 100 1000 400]);
    yyaxis left
    plot(t_com, Ca_avg_C1, 'k-', 'LineWidth', 1.8);
    ylabel('Calcium (ΔF/F)', 'FontSize', 12);
    yyaxis right
    plot(t_com, pup_avg_C1, '-', 'Color', [0.1 0.45 0.85], 'LineWidth', 1.8);
    ylabel('Pupil Area (px²)', 'FontSize', 12);
    xlabel('Time within trial (s)', 'FontSize', 12);
    title(sprintf('Stim %d — Trial-Averaged Ca (black) and Pupil (blue)', si), ...
          'FontSize', 13, 'Interpreter', 'none');
    xlim([0, stim_data{si}.trial_dur]);  set(gca, 'Box', 'off');
end

%% CELL C2: Per-Stimulus — Histogram Overlay (Ca vs Pupil distributions)
% Z-scored histograms overlaid to compare distribution shapes.

for si = 1:num_stim_x
    Ca_all  = stim_data{si}.Ca_trials;
    pup_all = stim_data{si}.pup_trials;

    if use_single()
        Ca_sel = squeeze(Ca_all(roi_id, :, :));
    else
        Ca_sel = squeeze(mean(Ca_all, 1, 'omitnan'));
    end
    Ca_avg_C2  = mean(Ca_sel, 2)';
    pup_avg_C2 = mean(pup_all, 2, 'omitnan')';

    ca_z_C2  = (Ca_avg_C2  - mean(Ca_avg_C2(:)))  / std(Ca_avg_C2(:));
    pup_z_C2 = (pup_avg_C2 - mean(pup_avg_C2(:))) / std(pup_avg_C2(:));

    figure('Name', sprintf('C2: Stim %d — Histograms  [%s]', si, stim_data{si}.label), ...
           'Position', [100 100 700 500]);
    hold on;
    histogram(ca_z_C2,  25, 'Normalization', 'probability', ...
              'FaceColor', 'k', 'FaceAlpha', 0.45, 'EdgeColor', 'none');
    histogram(pup_z_C2, 25, 'Normalization', 'probability', ...
              'FaceColor', [0.1 0.45 0.85], 'FaceAlpha', 0.45, 'EdgeColor', 'none');
    xlabel('Z-scored value', 'FontSize', 12);
    ylabel('Probability',    'FontSize', 12);
    title(sprintf('Stim %d — Value Distributions (z-scored)', si), 'FontSize', 13);
    legend({'Calcium', 'Pupil'}, 'Location', 'best');
    set(gca, 'Box', 'off');
end

%% CELL D: Grouped Average — Stimuli with Same type + addnote
% Only stimuli sharing BOTH type AND addnote are averaged together.
% Grey lines = per-stimulus grand mean.  Black/red = group mean ± SEM.
% Groups with mismatched trial durations use normalised time [0,1].

for g = 1:length(uniq_labels_x)
    members = find(label_grp_x == g);
    if length(members) < 2
        fprintf('Group %d [%s]: single member (Stim %d) — skipping.\n', ...
                g, uniq_labels_x{g}, members);
        continue
    end

    n_mem = length(members);

    % Use real time if all durations match, else normalise to [0,1]
    trial_durs_D = arrayfun(@(k) stim_data{k}.trial_dur, members);
    if range(trial_durs_D) < 0.01
        t_com_D   = stim_data{members(1)}.t_com;
        x_lbl_D   = 'Time within trial (s)';
        xlim_D    = stim_data{members(1)}.trial_dur;
    else
        t_com_D   = linspace(0, 1, n_com_x);
        x_lbl_D   = 'Normalised trial time (0–1)';
        xlim_D    = 1;
        fprintf('WARNING: Group %d durations differ (%s s) — using normalised time.\n', ...
                g, num2str(trial_durs_D'));
    end

    Ca_per_stim_D  = zeros(n_mem, n_com_x);
    pup_per_stim_D = zeros(n_mem, n_com_x);

    for k = 1:n_mem
        si_k     = members(k);
        Ca_all_k = stim_data{si_k}.Ca_trials;

        if use_single()
            Ca_sel_k = squeeze(Ca_all_k(roi_id, :, :));
        else
            Ca_sel_k = squeeze(mean(Ca_all_k, 1, 'omitnan'));
        end
        Ca_per_stim_D(k,:)  = mean(Ca_sel_k, 2)';
        pup_per_stim_D(k,:) = mean(stim_data{si_k}.pup_trials, 2, 'omitnan')';
    end

    Ca_grp_mean_D  = mean(Ca_per_stim_D,  1);
    Ca_grp_sem_D   = std(Ca_per_stim_D,  0, 1) / sqrt(n_mem);
    pup_grp_mean_D = mean(pup_per_stim_D, 1);
    pup_grp_sem_D  = std(pup_per_stim_D, 0, 1) / sqrt(n_mem);
    x_sh_D = [t_com_D, fliplr(t_com_D)];

    fig_lbl_D = sprintf('D: Group %d — Stim [%s]  |  %s', ...
                        g, num2str(members'), uniq_labels_x{g});
    figure('Name', fig_lbl_D, 'Position', [100 100 950 520]);

    % --- Calcium ---
    subplot(2, 1, 1);  hold on;
    for k = 1:n_mem
        plot(t_com_D, Ca_per_stim_D(k,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);
    end
    hf1 = fill(x_sh_D, [(Ca_grp_mean_D + Ca_grp_sem_D), fliplr(Ca_grp_mean_D - Ca_grp_sem_D)], ...
               'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    set(hf1, 'HandleVisibility', 'off');
    plot(t_com_D, Ca_grp_mean_D, 'k-', 'LineWidth', 2.5);
    ylabel('Calcium (ΔF/F)', 'FontSize', 11);
    title(fig_lbl_D, 'FontSize', 11, 'Interpreter', 'none');
    xlim([0, xlim_D]);  set(gca, 'Box', 'off');
    leg_lbl_D = [arrayfun(@(x) sprintf('Stim %d', x), members, 'UniformOutput', false); ...
                 {'Group mean ± SEM'}];
    legend(leg_lbl_D, 'Location', 'best', 'FontSize', 9);

    % --- Pupil ---
    subplot(2, 1, 2);  hold on;
    for k = 1:n_mem
        plot(t_com_D, pup_per_stim_D(k,:), 'Color', [1 0.7 0.7], 'LineWidth', 1.2);
    end
    hf2 = fill(x_sh_D, [(pup_grp_mean_D + pup_grp_sem_D), fliplr(pup_grp_mean_D - pup_grp_sem_D)], ...
               'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    set(hf2, 'HandleVisibility', 'off');
    plot(t_com_D, pup_grp_mean_D, 'r-', 'LineWidth', 2.5);
    xlabel(x_lbl_D, 'FontSize', 11);
    ylabel('Pupil Area (px²)', 'FontSize', 11);
    xlim([0, xlim_D]);  set(gca, 'Box', 'off');
end

%% =========================================================================
%  LOCOMOTION & WHISKING ANALYSIS  (Cells E1–E4)
%  Requires EXTENDED SETUP (stim_data, num_stim_x, num_trials_x, nROIs_x,
%  use_single, roi_id) to have been run first.
%  Run LOCO-WHISK SETUP, then any plot cell independently.
%  =========================================================================

%% LOCO-WHISK SETUP — extend stim_data with full-block behavioural signals
% Extracts the full stimulus block (all repeats as one continuous trace).
% Whisker and locomotion are resampled onto the calcium time grid via interp1.

loco_smooth_win = 10;   % movmean window (samples) applied to raw diff velocity
loco_run_thresh = 0.5;    % mm/s — locomotion above this = "running" (used in Cell E4)

for si = 1:num_stim_x
    blk_dur_si = data.Stimuli(si).stimulus_trial_t * data.Stimuli(si).trials;
    stim_blk_s = data.Stimuli(si).TimeStimulusFrame(1);
    stim_blk_e = stim_blk_s + blk_dur_si;

    % Calcium: full block, all ROIs
    Ca_idx_full  = find(data.TimeCa(1,:) > stim_blk_s & data.TimeCa(1,:) < stim_blk_e);
    t_full       = data.TimeCa(1, Ca_idx_full) - stim_blk_s;   % row, block-relative
    Ca_raw_full  = data.CaData(1).Ca_dFF(:, Ca_idx_full);      % [nROIs × nCa]
    ca_full_mean = mean(Ca_raw_full, 1, 'omitnan');             % [1 × nCa]
    ca_full_sem  = std(Ca_raw_full, 0, 1) / sqrt(nROIs_x);     % [1 × nCa]

    % Whisker + Pupil: column 1 only / PupilArea, ~50 Hz — same timestamps
    Beh_idx_full = find(data.Triggers.TimeCamera > stim_blk_s & data.Triggers.TimeCamera < stim_blk_e);
    t_beh        = data.Triggers.TimeCamera(Beh_idx_full) - stim_blk_s;
    whisk_raw    = double(data.Behav.whiskerFollicleA1(Beh_idx_full, 1));
    valid_wh     = ~isnan(whisk_raw);
    if sum(valid_wh) > 1
        whisk_full = interp1(t_beh(valid_wh), whisk_raw(valid_wh), t_full(:), 'linear', 'extrap');
    else
        whisk_full = zeros(numel(t_full), 1);
    end
    pup_raw_full = double(data.Behav.PupilArea(Beh_idx_full));
    pup_full_e   = interp1(t_beh(:), pup_raw_full(:), t_full(:), 'linear', 'extrap');

    % Locomotion: own timestamps, velocity from diff, smoothed
    Loc_idx_full = find(data.Triggers.TimeBall{si*2} > stim_blk_s & data.Triggers.TimeBall{si*2} < stim_blk_e);
    if length(Loc_idx_full) > 1
        t_loc     = data.Triggers.TimeBall{si*2}(Loc_idx_full) - stim_blk_s;
        p_loc     = data.LocomotionCal(1).Forward_mm(Loc_idx_full);
        v_raw     = diff(p_loc(:)) ./ diff(t_loc(:));
        t_vel     = t_loc(1:end-1);
        v_smooth  = movmean(v_raw, loco_smooth_win);
        loco_full = interp1(t_vel(:), v_smooth(:), t_full(:), 'linear', 'extrap');
    else
        loco_full = zeros(numel(t_full), 1);
    end

    % Trial boundary times (relative to block start) for vertical markers
    trial_bdry = zeros(1, data.Stimuli(si).trials);
    for tr = 1:data.Stimuli(si).trials
        trial_bdry(tr) = data.Stimuli(si).TrialTimes{1, tr}(1) - stim_blk_s;
    end

    stim_data{si}.t_full       = t_full(:)';
    stim_data{si}.ca_raw_full  = Ca_raw_full;
    stim_data{si}.ca_full_mean = ca_full_mean(:)';
    stim_data{si}.ca_full_sem  = ca_full_sem(:)';
    stim_data{si}.whisk_full   = whisk_full(:)';
    stim_data{si}.pup_full     = pup_full_e(:)';
    stim_data{si}.loco_full    = loco_full(:)';
    stim_data{si}.trial_bdry   = trial_bdry;
    stim_data{si}.blk_dur      = blk_dur_si;
end
fprintf('Loco-whisk setup done: %d stimuli.\n', num_stim_x);

%% CELL E1: Per-Stimulus — Continuous Full-Block Time Series
% 4 stacked panels: Ca mean ± SEM, pupil (blue), whisker (green), locomotion (orange).
% Dashed grey lines mark trial boundaries.  All panels linked for panning.

for si = 1:num_stim_x
    t     = stim_data{si}.t_full;
    ca_m  = stim_data{si}.ca_full_mean;
    ca_s  = stim_data{si}.ca_full_sem;
    pup   = stim_data{si}.pup_full;
    whisk = stim_data{si}.whisk_full;
    loco  = stim_data{si}.loco_full;
    bdry  = stim_data{si}.trial_bdry;
    blk_d = stim_data{si}.blk_dur;

    figure('Name', sprintf('E1: Stim %d — Full Block  [%s]', si, stim_data{si}.label), ...
           'Position', [100 100 1100 780]);

    ax1_E1 = subplot(4, 1, 1);  hold on;
    h_sh_E1 = fill([t, fliplr(t)], [(ca_m + ca_s), fliplr(ca_m - ca_s)], ...
                   [0.7 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    set(h_sh_E1, 'HandleVisibility', 'off');
    plot(t, ca_m, 'k-', 'LineWidth', 1.8);
    for k = 2:length(bdry)
        xline(bdry(k), '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
    ylabel('Calcium (ΔF/F)', 'FontSize', 11);
    title(sprintf('Stim %d — Full block (%d repeats)  [%s]', si, data.Stimuli(si).trials, stim_data{si}.label), ...
          'FontSize', 12, 'Interpreter', 'none');
    xlim([0, blk_d]);  set(gca, 'Box', 'off');

    ax2_E1 = subplot(4, 1, 2);  hold on;
    plot(t, pup, '-', 'Color', [0.1 0.45 0.85], 'LineWidth', 1.2);
    for k = 2:length(bdry)
        xline(bdry(k), '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
    ylabel('Pupil Area (px²)', 'FontSize', 11);
    xlim([0, blk_d]);  set(gca, 'Box', 'off');

    ax3_E1 = subplot(4, 1, 3);  hold on;
    plot(t, whisk, '-', 'Color', [0.2 0.6 0.2], 'LineWidth', 1.2);
    for k = 2:length(bdry)
        xline(bdry(k), '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
    ylabel('Whisker A1 (px)', 'FontSize', 11);
    xlim([0, blk_d]);  set(gca, 'Box', 'off');

    ax4_E1 = subplot(4, 1, 4);  hold on;
    plot(t, loco, '-', 'Color', [0.85 0.45 0.1], 'LineWidth', 1.2);
    yline(loco_run_thresh, ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 1);
    for k = 2:length(bdry)
        xline(bdry(k), '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
    xlabel('Time within block (s)', 'FontSize', 11);
    ylabel('Velocity (mm/s)', 'FontSize', 11);
    xlim([0, blk_d]);  set(gca, 'Box', 'off');

    linkaxes([ax1_E1, ax2_E1, ax3_E1, ax4_E1], 'x');
end

%% CELL E2: Per-Stimulus — Ca vs Whisker: Scatter + XCorr
% Respects roi_mode / use_single flag.  XCorr limited to ±10 s lag.

for si = 1:num_stim_x
    t     = stim_data{si}.t_full;
    whisk = stim_data{si}.whisk_full;
    Ca_rf = stim_data{si}.ca_raw_full;
    fs_si = numel(t) / stim_data{si}.blk_dur;

    if use_single()
        ca_vec = Ca_rf(roi_id, :);
    else
        ca_vec = mean(Ca_rf, 1, 'omitnan');
    end

    valid_E2 = ~(isnan(ca_vec) | isnan(whisk));
    ca_E2    = ca_vec(valid_E2)';
    wh_E2    = whisk(valid_E2)';

    [R_E2, P_E2] = corrcoef(wh_E2, ca_E2);
    r_E2 = R_E2(1,2);  p_E2 = P_E2(1,2);

    ca_z_E2  = (ca_E2 - mean(ca_E2)) / std(ca_E2);
    wh_z_E2  = (wh_E2 - mean(wh_E2)) / std(wh_E2);
    [xc_E2, lags_E2] = xcorr(ca_z_E2, wh_z_E2, round(fs_si * 10), 'coeff');
    lag_s_E2 = lags_E2 / fs_si;
    [~, pk_E2] = max(abs(xc_E2));

    figure('Name', sprintf('E2: Stim %d — Ca vs Whisker  [%s]', si, stim_data{si}.label), ...
           'Position', [100 100 950 430]);

    subplot(1, 2, 1);
    scatter(wh_E2, ca_E2, 5, [0.2 0.6 0.2], 'filled', 'MarkerFaceAlpha', 0.15);  hold on;
    pf_E2 = polyfit(wh_E2, ca_E2, 1);
    xf_E2 = linspace(min(wh_E2), max(wh_E2), 200);
    plot(xf_E2, polyval(pf_E2, xf_E2), 'r-', 'LineWidth', 2);
    text(0.05, 0.93, sprintf('r = %.3f   p = %.2e', r_E2, p_E2), ...
        'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k');
    xlabel('Whisker A1 (px)', 'FontSize', 11);
    ylabel('Calcium (ΔF/F)',  'FontSize', 11);
    title(sprintf('Stim %d — Scatter', si), 'FontSize', 12);
    set(gca, 'Box', 'off');

    subplot(1, 2, 2);
    plot(lag_s_E2, xc_E2, '-', 'Color', [0.2 0.6 0.2], 'LineWidth', 1.5);  hold on;
    xline(0, 'r--', 'LineWidth', 1);
    plot(lag_s_E2(pk_E2), xc_E2(pk_E2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    xlabel('Lag (s)',                 'FontSize', 11);
    ylabel('Correlation Coefficient', 'FontSize', 11);
    title(sprintf('Stim %d — XCorr Ca–Whisker  (peak %.2f s)', si, lag_s_E2(pk_E2)), 'FontSize', 12);
    set(gca, 'Box', 'off');  grid on;
end

%% CELL E3: Per-Stimulus — Ca vs Locomotion: Scatter + XCorr
% Respects roi_mode / use_single flag.  XCorr limited to ±10 s lag.

for si = 1:num_stim_x
    t    = stim_data{si}.t_full;
    loco = stim_data{si}.loco_full;
    Ca_rf = stim_data{si}.ca_raw_full;
    fs_si = numel(t) / stim_data{si}.blk_dur;

    if use_single()
        ca_vec = Ca_rf(roi_id, :);
    else
        ca_vec = mean(Ca_rf, 1, 'omitnan');
    end

    valid_E3 = ~(isnan(ca_vec) | isnan(loco));
    ca_E3    = ca_vec(valid_E3)';
    lo_E3    = loco(valid_E3)';

    [R_E3, P_E3] = corrcoef(lo_E3, ca_E3);
    r_E3 = R_E3(1,2);  p_E3 = P_E3(1,2);

    ca_z_E3  = (ca_E3 - mean(ca_E3)) / std(ca_E3);
    lo_z_E3  = (lo_E3 - mean(lo_E3)) / std(lo_E3);
    [xc_E3, lags_E3] = xcorr(ca_z_E3, lo_z_E3, round(fs_si * 10), 'coeff');
    lag_s_E3 = lags_E3 / fs_si;
    [~, pk_E3] = max(abs(xc_E3));

    figure('Name', sprintf('E3: Stim %d — Ca vs Locomotion  [%s]', si, stim_data{si}.label), ...
           'Position', [100 100 950 430]);

    subplot(1, 2, 1);
    scatter(lo_E3, ca_E3, 5, [0.85 0.45 0.1], 'filled', 'MarkerFaceAlpha', 0.15);  hold on;
    pf_E3 = polyfit(lo_E3, ca_E3, 1);
    xf_E3 = linspace(min(lo_E3), max(lo_E3), 200);
    plot(xf_E3, polyval(pf_E3, xf_E3), 'r-', 'LineWidth', 2);
    text(0.05, 0.93, sprintf('r = %.3f   p = %.2e', r_E3, p_E3), ...
        'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k');
    xlabel('Velocity (mm/s)',  'FontSize', 11);
    ylabel('Calcium (ΔF/F)',   'FontSize', 11);
    title(sprintf('Stim %d — Scatter', si), 'FontSize', 12);
    set(gca, 'Box', 'off');

    subplot(1, 2, 2);
    plot(lag_s_E3, xc_E3, '-', 'Color', [0.85 0.45 0.1], 'LineWidth', 1.5);  hold on;
    xline(0, 'r--', 'LineWidth', 1);
    plot(lag_s_E3(pk_E3), xc_E3(pk_E3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    xlabel('Lag (s)',                 'FontSize', 11);
    ylabel('Correlation Coefficient', 'FontSize', 11);
    title(sprintf('Stim %d — XCorr Ca–Loco  (peak %.2f s)', si, lag_s_E3(pk_E3)), 'FontSize', 12);
    set(gca, 'Box', 'off');  grid on;
end

%% CELL E4: Per-Stimulus — Running State Analysis
% Panel 1: Ca trace with running epochs shaded orange.
% Panel 2: Histogram overlay — Ca values during running (orange) vs stationary (grey).
% Panel 3: Locomotion-onset-triggered Ca average ± SEM.
%          Only bouts longer than min_bout_dur are counted as valid onsets.
% Adjust loco_run_thresh and min_bout_dur in LOCO-WHISK SETUP / here as needed.

min_bout_dur = 1;   % seconds — ignore brief locomotion blips
ota_pre_s    = 3;   % seconds before onset to include
ota_post_s   = 7;   % seconds after onset to include

for si = 1:num_stim_x
    t     = stim_data{si}.t_full;
    ca_m  = stim_data{si}.ca_full_mean;
    loco  = stim_data{si}.loco_full;
    blk_d = stim_data{si}.blk_dur;
    fs_si = numel(t) / blk_d;

    % Identify running bouts, filter by minimum duration
    running_bin = loco > loco_run_thresh;
    run_s_all = find(diff([false, running_bin]) ==  1);
    run_e_all = find(diff([running_bin, false])  == -1);
    n_all     = min(numel(run_s_all), numel(run_e_all));
    run_s_all = run_s_all(1:n_all);
    run_e_all = run_e_all(1:n_all);
    long_bout = (run_e_all - run_s_all + 1) >= round(min_bout_dur * fs_si);
    run_starts = run_s_all(long_bout);
    run_ends   = run_e_all(long_bout);
    n_bouts    = numel(run_starts);

    % Running mask from long bouts only
    running_mask = false(size(loco));
    for k = 1:n_bouts
        running_mask(run_starts(k):run_ends(k)) = true;
    end

    if sum(running_mask) < 5
        fprintf('Stim %d: no sustained running above %.1f mm/s — skipping E4.\n', si, loco_run_thresh);
        continue
    end

    ca_run  = ca_m(running_mask);
    ca_stat = ca_m(~running_mask);

    ca_pad = 0.5 * std(ca_m);
    y_lo   = min(ca_m) - ca_pad;
    y_hi   = max(ca_m) + ca_pad;

    figure('Name', sprintf('E4: Stim %d — Running State  [%s]', si, stim_data{si}.label), ...
           'Position', [100 100 1100 750]);

    % --- Panel 1: Ca trace with running epochs ---
    subplot(3, 1, 1);  hold on;
    for k = 1:n_bouts
        fill([t(run_starts(k)), t(run_ends(k)), t(run_ends(k)), t(run_starts(k))], ...
             [y_lo, y_lo, y_hi, y_hi], ...
             [1 0.7 0.3], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    plot(t, ca_m, 'k-', 'LineWidth', 1.8);
    ylim([y_lo, y_hi]);  xlim([0, blk_d]);
    ylabel('Calcium (ΔF/F)', 'FontSize', 11);
    title(sprintf('Stim %d — Ca with running epochs (orange, >%.1f mm/s, >%.0f s)  [%s]', ...
          si, loco_run_thresh, min_bout_dur, stim_data{si}.label), 'FontSize', 11, 'Interpreter', 'none');
    set(gca, 'Box', 'off');

    % --- Panel 2: Histogram overlay running vs stationary ---
    subplot(3, 1, 2);  hold on;
    histogram(ca_stat, 30, 'Normalization', 'probability', ...
              'FaceColor', [0.6 0.6 0.6], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    histogram(ca_run,  30, 'Normalization', 'probability', ...
              'FaceColor', [1 0.6 0.1],     'FaceAlpha', 0.6, 'EdgeColor', 'none');
    xline(median(ca_stat), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xline(median(ca_run),  '--', 'Color', [0.85 0.35 0],  'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlabel('Calcium (ΔF/F)', 'FontSize', 11);
    ylabel('Probability',    'FontSize', 11);
    title(sprintf('Distribution  —  stationary (grey) vs running (orange)  |  median: stat=%.3f  run=%.3f', ...
          median(ca_stat), median(ca_run)), 'FontSize', 10);
    legend({'Stationary', 'Running'}, 'Location', 'best', 'FontSize', 9);
    set(gca, 'Box', 'off');

    % --- Panel 3: Onset-triggered Ca average ---
    win_pre  = round(ota_pre_s  * fs_si);
    win_post = round(ota_post_s * fs_si);
    t_ota    = linspace(-ota_pre_s, ota_post_s, win_pre + win_post + 1);
    valid_ons = run_starts(run_starts > win_pre & run_starts <= numel(t) - win_post);

    subplot(3, 1, 3);  hold on;
    if numel(valid_ons) >= 2
        ota_mat = zeros(numel(valid_ons), win_pre + win_post + 1);
        for oi = 1:numel(valid_ons)
            ota_mat(oi, :) = ca_m((valid_ons(oi) - win_pre) : (valid_ons(oi) + win_post));
        end
        ota_mean = mean(ota_mat, 1);
        ota_sem  = std(ota_mat, 0, 1) / sqrt(numel(valid_ons));
        x_sh_ota = [t_ota, fliplr(t_ota)];
        y_sh_ota = [(ota_mean + ota_sem), fliplr(ota_mean - ota_sem)];
        h_ota = fill(x_sh_ota, y_sh_ota, [0.7 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
        set(h_ota, 'HandleVisibility', 'off');
        plot(t_ota, ota_mean, 'k-', 'LineWidth', 1.8);
        xline(0, '--', 'Color', [0.85 0.45 0.1], 'LineWidth', 1.5, 'Label', 'Onset', ...
              'LabelVerticalAlignment', 'bottom');
        xlabel('Time from locomotion onset (s)', 'FontSize', 11);
        ylabel('Calcium (ΔF/F)', 'FontSize', 11);
        title(sprintf('Stim %d — Onset-triggered Ca mean ± SEM  (n = %d onsets)', si, numel(valid_ons)), ...
              'FontSize', 12);
    else
        text(0.5, 0.5, sprintf('Fewer than 2 valid onsets (n=%d)\nLower loco_run_thresh or min_bout_dur', ...
             numel(valid_ons)), 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 11);
        title(sprintf('Stim %d — Onset-triggered Ca  (insufficient onsets)', si), 'FontSize', 12);
    end
    xlim([-ota_pre_s, ota_post_s]);  set(gca, 'Box', 'off');
end