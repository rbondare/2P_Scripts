
data= load('\\fs.ista.ac.at\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB6_250302_1656_preprocessed.mat'); 



%%
dFF = data.CaData(1).Ca_dFF;
centroid = data.CaData(1).Ca_centroid_voxel;

dFF_norm = (dFF - min(dFF,[],2)) ./ (max(dFF,[],2) - min(dFF,[],2)); % Normalize each row
roi_indices = [1,2,3, 4, 5, 6,7,8,9, 10,11]; % Replace with your desired ROI indices
dFF_selected = dFF_norm(roi_indices, 1000:2000); % Extract only these ROIs

figure;
imagesc(dFF_selected)
colormap('turbo'); % Use a colormap like 'hot' or 'jet'
colorbar;  % Add a color scale


%%

% Define stimulus time window
fs = 60;
num_trials = data.Stimuli.trials
trial_duration = 30; 
stim_total_time = trial_duration * num_trials; 

start_time = data.Stimuli(1).TimeStimulusFrame(1);
end_time = start_time + stim_total_time;
[~, start_indx] = min(abs(data.TimeCa(2,:) - start_time));
[~, end_indx]   = min(abs(data.TimeCa(2,:) - end_time));









% Select ROIs and extract calcium data for the stimulus window.
roi_indices = [1,2,3,4,5,6,7,8,9,10,11];  % example ROI indices
dFF_selected = dFF_norm(roi_indices, start_indx:end_indx); 

% Create a figure with two subplots
figure;
subplot(2,1,1);
imagesc(dFF_selected);
colormap('turbo');
colorbar;
xlabel('Time (frames)');
ylabel('ROI');
title('Calcium Activity Heatmap');

% Extract durations from your data structure (in seconds)
durations = data.Stimuli(1).specParams.fff_durations;  
stim_signal_all = [];
for trial = 1:num_trials
    stim_signal_trial = [];
    for i = 1:length(durations)
        % Create time vector for current segment
        t_seg = linspace(0, durations(i), durations(i)*fs);
        % Generate the segment based on the stimulus type
        switch i
            case 1  % Segment 1: black
                seg = zeros(size(t_seg)); 
            case 2  % Segment 2: white
                seg = ones(size(t_seg));
            case 3  % Segment 3: grey
                seg = 0.5 * ones(size(t_seg));
            case 4  % Segment 4: black
                seg = zeros(size(t_seg));
            case 5  % Segment 5: grey
                seg = 0.5 * ones(size(t_seg));
            case 6  % Segment 6: sine amplitude (one cycle sine wave modulated around grey)
                seg = 0.5 + 0.5 * sin(2 * pi * (1/durations(i)) * t_seg);
            case 7  % Segment 7: grey
                seg = 0.5 * ones(size(t_seg));
            case 8  % Segment 8: sine frequency (sine wave with frequency factor 2)
                seg = 0.5 + 0.5 * sin(2 * pi * 2 * t_seg);
            case 9  % Segment 9: grey
                seg = 0.5 * ones(size(t_seg));
        end
        % Append the segment to the trial's stimulus signal
        stim_signal_trial = [stim_signal_trial, seg];
    end
    % Concatenate the trial's stimulus signal to the overall signal
    stim_signal_all = [stim_signal_all, stim_signal_trial];
end

% Create a time vector for the concatenated stimulus trace
total_time = trial_duration * num_trials;  % total seconds across all trials
t_stim_all = linspace(0, total_time, length(stim_signal_all));

subplot(2,1,2);
plot(t_stim_all, stim_signal_all, 'k', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Luminance');
title('Visual Stimulus (Full Field Flashes) across 5 Trials');
ylim([-0.1 1.1]); % set y-limits to clearly show black, grey, and white
grid on;


%%







%%





% Load pupil area from your data
pupil_area = data.Behav.PupilArea;

% Downsample pupil size using interpolation
pupil_downsampled = interp1(1:length(pupil_area), pupil_area, linspace(1, length(pupil_area), size(dFF,2)), 'linear');

% Normalize (Z-score)
pupil_z = (pupil_downsampled - mean(pupil_downsampled)) / std(pupil_downsampled); 
dFF_z = (dFF - mean(dFF, 2)) ./ std(dFF, 0, 2); 

% Compute correlation between pupil size and calcium activity for each ROI
num_ROIs = size(dFF, 1);
correlations = zeros(1, num_ROIs);  % Fixed typo in `zeros`

for i = 1:num_ROIs 
    correlations(i) = corr(pupil_z', dFF_z(i, :)');  % Ensure both inputs are column vectors
end

% Display correlation coefficients
disp('Correlation coefficients for each ROI:')
disp(correlations)

% Plot time-series of pupil size and calcium activity (ROI 1)
figure;
plot(pupil_z, 'b', 'LineWidth', 1); hold on;
plot(dFF_z(1,:), 'r', 'LineWidth', 1);
xlabel('Time (frames)'); 
ylabel('Z-scored activity');
title('Pupil Size vs Calcium Activity (ROI 1)');
legend('Pupil Size', 'Calcium Activity');
hold off;


%% Plot Calcium Activity with Pupil Diameter Below
figure;
t = 1:size(dFF_z,2); % Time axis (frames)

% Top subplot: Calcium activity (example ROI 1)
subplot(2,1,1);
plot(t, dFF_z(1,:), 'k', 'LineWidth', 1);
xlabel('Time (frames)');
ylabel('Ca^{2+} Activity (Z-scored)');
title('Calcium Activity (ROI 1)');
legend('Calcium Activity');
ylim([min(dFF_z(1,:)) max(dFF_z(1,:))]); % Adjust based on data

% Bottom subplot: Pupil diameter (Z-scored)
subplot(2,1,2);
plot(t, pupil_z, 'b', 'LineWidth', 1);
xlabel('Time (frames)');
ylabel('Pupil Diameter (Z-scored)');
title('Pupil Diameter Over Time');
legend('Pupil Diameter');
ylim([min(pupil_z) max(pupil_z)]); % Adjust based on data

%% Plot Absolute Pupil Area Over Time
figure;
plot(1:length(pupil_area), pupil_area, 'k', 'LineWidth', 1);
xlabel('Time (frames)');
ylabel('Pupil Area (absolute units)');
title('Absolute Pupil Area Over Time');
legend('Pupil Area');
grid on;


% Scatter plot for ROI 1
figure;
scatter(pupil_z, dFF_z(1,:), 10, 'filled');
xlabel('Pupil Size (z-scored)'); ylabel('Calcium Activity (z-scored)');
title(sprintf('Scatter Plot (ROI 1): r = %.2f', correlations(1)));
grid on;

% Bar plot of correlations across ROIs
figure;
bar(correlations);
xlabel('ROI'); ylabel('Correlation (r)');
title('Correlation Between Pupil Size and Calcium Activity Across ROIs');
yline(0, '--k'); % Add reference line at r=0
grid on;

% Load pupil area from your data
pupil_area = data.Behav.PupilArea;

% Downsample pupil size using interpolation to match calcium imaging frames
pupil_downsampled = interp1(1:length(pupil_area), pupil_area, linspace(1, length(pupil_area), size(dFF,2)), 'linear');

% Normalize both signals (Z-score)
pupil_z = (pupil_downsampled - mean(pupil_downsampled)) / std(pupil_downsampled); 
dFF_z = (dFF - mean(dFF, 2)) ./ std(dFF, 0, 2); 

% Time axis (frames) for calcium imaging
t = (1:size(dFF_z,2)) / 10; 

% Select an example ROI to visualize
example_ROI = 4;


figure;
% Top subplot: Calcium activity
subplot(2,1,1);
plot(t, dFF_z(example_ROI, :), 'k', 'LineWidth', 1);
ylabel('Ca^{2+} Activity (Z-scored)');
title(sprintf('Calcium Activity (ROI %d)', example_ROI));
xlim([t(1) t(end)]);
ylim([min(dFF_z(example_ROI,:)) max(dFF_z(example_ROI,:))]); % Adjust based on data
hold on;
grid on;

% Bottom subplot: Pupil area (absolute or z-scored)
subplot(2,1,2);
plot(t, pupil_z, 'b', 'LineWidth', 1);
xlabel('Time (seconds)');
ylabel('Pupil Diameter (Z-scored)');
title('Pupil Diameter Over Time');
xlim([t(1) t(end)]);
ylim([min(pupil_z) max(pupil_z)]); % Adjust based on data
grid on;

% Link axes for synchronized zooming
linkaxes(findall(gcf,'Type','axes'),'x');


%%

figure;
hold on;
plot(dFF(1,1000:2000), 'b', 'LineWidth', 1); 
plot(dFF(4,1000:2000), 'r', 'LineWidth', 1); 
plot(dFF(11,1000:2000), 'k', 'LineWidth', 1); 
axis tight

% Formatting
xlabel('Frame');
ylabel('\DeltaF/F');
legend({'ROI 1', 'ROI 2'}, 'Location', 'Best');
grid on;
hold off;