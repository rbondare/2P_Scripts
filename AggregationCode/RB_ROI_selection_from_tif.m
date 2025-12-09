%% nLight ROI detection
clear; clc; close all;

%% --- USER INPUTS ---
tiffFolder = uigetdir(pwd, 'Select folder with TIFF stack');
roi_size = 30;          % ROI side length (pixels)
step = 20;              % step between ROIs
topN = 10;              % number of top ROIs
min_baseline_frac = 0.1; % discard ROIs whose mean F0 < 10% of global max

%% --- LOAD ALL TIFFS (multi-frame corrected) ---
fprintf('Loading TIFFs from %s...\n', tiffFolder);
files = dir(fullfile(tiffFolder, '*.tif'));
[~, order] = sort({files.name});
files = files(order);

% Get info from the first file
info = imfinfo(fullfile(tiffFolder, files(1).name));
nx = info(1).Height;
ny = info(1).Width;
framesPerFile = numel(info);

nt = framesPerFile * numel(files);  % total frames
movie = zeros(nx, ny, nt, 'single');

frameIdx = 1;
for f = 1:numel(files)
    infoF = imfinfo(fullfile(tiffFolder, files(f).name));
    for k = 1:numel(infoF)
        movie(:,:,frameIdx) = imread(fullfile(tiffFolder, files(f).name), k);
        frameIdx = frameIdx + 1;
    end
end
fprintf('Loaded %d frames (%dx%d)\n', nt, nx, ny);

%% --- NORMALIZE AND DFF ---
F0 = prctile(movie, 10, 3);
movie_dff = (movie - F0) ./ max(F0, eps);
movie_dff = movie_dff - mean(movie_dff,3); % remove global mean

%% --- ROI SCAN ---
x_centers = 1+roi_size/2 : step : nx - roi_size/2;
y_centers = 1+roi_size/2 : step : ny - roi_size/2;

fprintf('Scanning ROIs...\n');
activity = zeros(length(x_centers), length(y_centers));
baseline = zeros(length(x_centers), length(y_centers));

for ix = 1:length(x_centers)
    for iy = 1:length(y_centers)
        x0 = round(x_centers(ix) - roi_size/2);
        y0 = round(y_centers(iy) - roi_size/2);
        roi = squeeze(mean(mean(movie(x0:x0+roi_size-1, y0:y0+roi_size-1, :),1),2));

        % baseline fluorescence (dark areas will have low values)
        baseline(ix,iy) = prctile(roi, 10);

        % transient activity = energy of positive fluorescence increases
        dtrace = diff(roi);
        posEnergy = mean(dtrace(dtrace>0).^2);
        activity(ix,iy) = posEnergy;
    end
end

%% --- FILTER ROIS ---
baseline_thresh = min_baseline_frac * max(baseline(:));
validMask = baseline > baseline_thresh;

% Only keep ROIs with decent baseline and strong transients
activity_filtered = activity .* validMask;

%% --- SELECT NON-OVERLAPPING TOP ROIS ---
roi_mask = false(size(activity_filtered));
selected_idx = [];
[sortedActivity, sortedIdx] = sort(activity_filtered(:), 'descend');
for k = 1:length(sortedIdx)
    [ixk, iyk] = ind2sub(size(activity_filtered), sortedIdx(k));
    
    % Skip if already masked
    if roi_mask(ixk, iyk)
        continue
    end
    
    % Select this ROI
    selected_idx(end+1) = sortedIdx(k); %#ok<SAGROW>
    
    % Mask nearby ROIs to avoid overlap
    [X, Y] = meshgrid(1:size(activity_filtered,2), 1:size(activity_filtered,1));
    dist2 = (X - iyk).^2 + (Y - ixk).^2;
    roi_mask(dist2 <= (roi_size/step)^2) = true;
    
    % Stop if reached desired number
    if numel(selected_idx) >= topN
        break
    end
end

[ix, iy] = ind2sub(size(activity_filtered), selected_idx);
centers = [x_centers(ix)', y_centers(iy)'];

%% --- OPTIONAL: SCROLLABLE MOVIE WITH ROIS ---
fig = figure('Name','Scrollable Movie with ROIs','NumberTitle','off');
ax = axes('Parent',fig);
hIm = imshow(movie(:,:,1), [], 'Parent', ax); hold on;

% Draw ROI rectangles
roiRects = gobjects(numel(ix),1);
for n = 1:numel(ix)
    roiRects(n) = rectangle('Position',[centers(n,2)-roi_size/2, centers(n,1)-roi_size/2, roi_size, roi_size], ...
                            'EdgeColor','r','LineWidth',1.5);
    text(centers(n,2), centers(n,1), num2str(n), 'Color','y','FontSize',9);
end

% Slider to scroll frames
uicontrol('Style', 'slider', ...
    'Min',1, 'Max', nt, 'Value',1, 'SliderStep',[1/(nt-1) 10/(nt-1)], ...
    'Position', [150 10 300 20], ...
    'Callback', @(src,evt) set(hIm,'CData',movie(:,:,round(src.Value))));

%% --- EXTRACT TRACES ---
F = zeros(numel(ix), nt);
for n = 1:numel(ix)
    x0 = round(centers(n,1) - roi_size/2);
    y0 = round(centers(n,2) - roi_size/2);
    roi = movie(x0:x0+roi_size-1, y0:y0+roi_size-1, :);
    F(n,:) = squeeze(mean(roi, [1 2]));
end
F0 = prctile(F,10,2);
dFF = (F - F0) ./ F0;

%% --- VISUALIZE ROI TRACES ---
figure('Position',[100 100 1000 400]);
subplot(1,2,1);
imagesc(mean(movie,3)); axis image off; colormap(gray);
title('Top non-overlapping ROIs (filtered for true fluorescence)');
hold on;
for n = 1:numel(ix)
    rectangle('Position',[centers(n,2)-roi_size/2, centers(n,1)-roi_size/2, roi_size, roi_size], ...
        'EdgeColor','r','LineWidth',1.5);
    text(centers(n,2), centers(n,1), num2str(n), 'Color','y','FontSize',9);
end
hold off;

subplot(1,2,2);
plot(dFF');
xlabel('Frame'); ylabel('\DeltaF/F');
title('ROI fluorescence traces');
sgtitle('nLight ROI detection (transient + baseline filtered, non-overlapping)');

fprintf('\n✅ Done. Displaying %d filtered non-overlapping ROIs.\n', numel(ix));
