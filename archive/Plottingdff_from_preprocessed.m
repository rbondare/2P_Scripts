%% This script plots the calcium data as traces, heatmap and xy coordinates for initial data visualisation
%%load preprocessed data file

filename=load('Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260305_1327_preprocessed.mat');

%%

% Extract calcium data and time vector
CalciumData = filename.CaData(1).Ca_dFF;
TimeCa = filename.TimeCa;

% Normalize for visualization
CalciumData_norm = CalciumData - min(CalciumData, [], 'all');
CalciumData_norm = CalciumData_norm / max(CalciumData_norm, [], 'all');

% Extract neuron metadata
Centroid = filename.CaData(1).Ca_centroid_voxel;
numCells = size(CalciumData, 1);
numTimepoints = size(CalciumData, 2);

fprintf('Data loaded successfully:\n');
fprintf('  Number of neurons: %d\n', numCells);
fprintf('  Number of timepoints: %d\n', numTimepoints);
%% Plotting individual dff traces
figure;

% Plot the first 5 cells as an example
for i = 1:min(numCells, 5)
    subplot(5, 1, i);
    plot(CalciumData(i, :));
    title(['Cell ', num2str(i)]);
    xlabel('Time');
    ylabel('dF/F');
end

% Adjust subplot spacing
sgtitle('Calcium Imaging Data (dF/F)');

%%

% Set the number of cells to plot
numCellsToPlot = min(numCells, 5);

% Determine the number of rows and columns for subplots
numRows = ceil(sqrt(numCellsToPlot)); % Number of rows
numCols = ceil(numCellsToPlot / numRows); % Number of columns

% Create a figure
figure;

% Plot each cell's activity in a subplot
for i = 1:numCellsToPlot
    subplot(numRows, numCols, i);
    plot(TimeCa(1, :), CalciumData(i, :)); % Use the time vector for x-axis
    title(['Cell ', num2str(i)]);
    xlabel('Time (s)');
    ylabel('dF/F');
end
sgtitle('Calcium  Data (dF/F)');

%%
figure('Name', 'Population Heatmap - Time', 'NumberTitle', 'off');

numCellsToPlot = 200;

% Calculate average activity per neuron
avg_activity = mean(CalciumData, 2);
[~, sort_order] = sort(avg_activity, 'descend');
Ca_sorted = CalciumData(sort_order, :);

% Reversed grayscale colormap for better contrast
original_colormap = gray;
reversed_colormap = flipud(original_colormap);

imagesc(Ca_sorted(1:numCellsToPlot, 5500:7000));  % Plot data for the first numCellsToPlot cells
colormap(reversed_colormap);
set(gca, 'YDir', 'reverse');

c = colorbar;
c.Label.String = 'dF/F';
c.Label.FontSize = 11;
%caxis([0 max(CalciumData_norm, [], 'all')]);
caxis([0 8])

xlabel('Time (s)', 'FontSize', 11);
ylabel('Neuron Index (sorted by avg activity)', 'FontSize', 11);
title('Population Heatmap: All Neurons', 'FontSize', 12, 'FontWeight', 'bold');
grid off;

%%
% Plot the heatmap

% Z-score normalization 
%mean_dFF = mean(calcium_activity1, 2); % Calculate mean across time for each cell
%std_dFF = std(calcium_activity1, 0, 2); % Calculate standard deviation across time for each cell
%normalized_data = (calcium_activity1 - mean_dFF) ./ std_dFF;

%CalciumDataZ = zscore(calcium_activity1')';

numCells = size(CalciumData, 1);

numCellsToPlot = min(numCells, 400);

avg_activity = mean(CalciumData, 2);

% Sort neurons by peak activity
[~, sort_order] = sort(avg_activity, 'descend'); % Sort in descending order

Ca_sorted = CalciumData(sort_order, :);

original_colormap = gray;  % Get the default gray colormap
reversed_colormap = flipud(original_colormap);

figure;
imagesc(CalciumData(1:numCellsToPlot, 5500:6500));  % Plot data for the first numCellsToPlot cells
colorbar;  % Add colorbar for reference
title('Heatmap of Calcium Activity');
xlabel('Frame number');
ylabel('Neuron Index');

% Adjust colormap (optional)
%colormap(gray);  % Use 'jet' colormap for better visualization
colormap(reversed_colormap);

% Reverse the y-axis so the first cell is at the top
set(gca, 'YDir', 'reverse');

% Optionally, adjust axis limits for better visualization (e.g., xlim, ylim)

% Optionally, add labels to the colorbar and adjust its appearance
c = colorbar;
caxis([0 5]); % Adjust as needed based on the range of your data
c.Label.String = 'dF/F';
c.Label.FontSize = 12;

% filename = 'Heatmap_Time_.png';  % Specify your filename and format
% saveas(gcf, filename);

%% plotting location of the cells (x,y grid) 
figure('Name', 'Neuron Spatial Locations', 'NumberTitle', 'off');

% Filter neurons for first plane (z == 1)
centroidX = Centroid(:, 1);
centroidY = Centroid(:, 2);
centroidZ = Centroid(:, 3);
plane_indices = find(centroidZ == 2);

x_coords = centroidX(plane_indices);
y_coords = centroidY(plane_indices);

scatter(y_coords, x_coords, 50, 'filled', 'MarkerFaceColor', '#4DBEEE', ...
        'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.3);

xlabel('Y Position (µm)', 'FontSize', 11);
ylabel('X Position (µm)', 'FontSize', 11);
title(sprintf('Neuron Locations: %d neurons', length(plane_indices)), ...
      'FontSize', 12, 'FontWeight', 'bold');
axis equal;
grid on;
grid minor;

