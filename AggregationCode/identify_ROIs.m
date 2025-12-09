clear; clc; close all;
%% Step 1: Load all TIFF frames into a 3D stack
tiffFolder = uigetdir(pwd, 'Select folder with TIFFs');
files = dir(fullfile(tiffFolder, '*.tif'));
[~,order] = sort({files.name});
files = files(order);

% Get image size from first TIFF
info1 = imfinfo(fullfile(tiffFolder, files(1).name));
nx = info1(1).Height;
ny = info1(1).Width;

% Count total frames
totalFrames = 0;
for f = 1:numel(files)
    infoF = imfinfo(fullfile(tiffFolder, files(f).name));
    totalFrames = totalFrames + numel(infoF);
end

% Preallocate stack
stack = zeros(nx, ny, totalFrames, 'single');

% Load frames
frameIdx = 1;
fprintf('Loading TIFF frames...\n');
for f = 1:numel(files)
    fname = fullfile(tiffFolder, files(f).name);
    t = Tiff(fname, 'r');
    while true
        stack(:,:,frameIdx) = single(t.read());
        frameIdx = frameIdx + 1;
        if t.lastDirectory()
            break;
        else
            t.nextDirectory();
        end
    end
    t.close();
end
fprintf('Done! Stack size: %d x %d x %d\n', size(stack));

size(stack)  % should be 512 x 512 x ~29,500
imshow(max(stack,[],3), []); title('Max projection')


%% Step 2: Define square tiles
roi_size = 10;     % pixels
step = roi_size;   % non-overlapping

x_centers = round((roi_size/2) : step : (nx - roi_size/2));
y_centers = round((roi_size/2) : step : (ny - roi_size/2));
nX = numel(x_centers);
nY = numel(y_centers);

fprintf('Tile grid: %d x %d = %d tiles\n', nX, nY, nX*nY);


%% Step 3: Sum pixels within each tile to get temporal trace
tileTraces = zeros(nX, nY, totalFrames, 'single');

for ix = 1:nX
    x0 = max(1, x_centers(ix) - floor(roi_size/2));
    x1 = min(nx, x0 + roi_size - 1);
    for iy = 1:nY
        y0 = max(1, y_centers(iy) - floor(roi_size/2));
        y1 = min(ny, y0 + roi_size - 1);
        roiStack = stack(x0:x1, y0:y1, :);
        tileTraces(ix, iy, :) = mean(reshape(roiStack, [], totalFrames), 1); % mean over pixels
    end
end

% Quick check: plot one tile trace
figure; plot(squeeze(tileTraces(10,10,:))); xlabel('Frame'); ylabel('Fluorescence'); title('Example tile trace');


%% Step 4: Compute max delta from baseline
baseline = prctile(tileTraces, 10, 3); % baseline per tile
maxDF = max(tileTraces - baseline, [], 3); % max ΔF per tile

% Visualize max ΔF across the grid
figure; imagesc(maxDF'); axis image; colorbar; title('Max ΔF per tile');


%% Step 5: Threshold to remove low-activity tiles
thresholdFactor = 1; % times the median ΔF across all tiles
thresh = thresholdFactor * median(maxDF(:));

activeMask = maxDF > thresh;

figure; imagesc(activeMask'); axis image; title('Active tiles after thresholding');


%% Step 6: Select top N non-overlapping ROIs 
topN = 8;
roi_mask = false(size(maxDF));
selected_idx = [];

[~, sortedIdx] = sort(maxDF(:), 'descend');

for k = 1:length(sortedIdx)
    idx = sortedIdx(k);
    [ix, iy] = ind2sub(size(maxDF), idx);

    % Skip if tile is already masked or below threshold
    if roi_mask(ix, iy) || ~activeMask(ix, iy)
        continue;
    end

    selected_idx(end+1) = sub2ind(size(maxDF), ix, iy); %#ok<SAGROW>

    % Block immediate neighbors (3x3 square around selected tile)
    for dx = -1:1
        for dy = -1:1
            nx_ = ix + dx;
            ny_ = iy + dy;
            if nx_ >=1 && nx_ <= size(maxDF,1) && ny_ >=1 && ny_ <= size(maxDF,2)
                roi_mask(nx_, ny_) = true;
            end
        end
    end

    if numel(selected_idx) >= topN
        break;
    end
end

[ix_sel, iy_sel] = ind2sub(size(maxDF), selected_idx);
centers = [x_centers(ix_sel)', y_centers(iy_sel)'];
fprintf('Selected %d ROIs\n', size(centers,1));


%% Step 7: Overlay ROIs on max projection
maxProj = max(stack, [], 3);
figure; imshow(maxProj, []); hold on;
title('Selected ROIs');

half = floor(roi_size/2);
for r = 1:size(centers,1)
    cx = centers(r,1); cy = centers(r,2);
    rectangle('Position',[cy-half, cx-half, roi_size, roi_size], 'EdgeColor','r','LineWidth',1.2);
    text(cy, cx, num2str(r), 'Color','y','FontSize',10,'HorizontalAlignment','center');
end
hold off;

%% Step 8: Interactive viewer
fig = figure('Name','ROI Viewer','NumberTitle','off','Position',[100 100 800 600]);
ax = axes('Parent', fig);
hIm = imshow(stack(:,:,1), [], 'Parent', ax); colormap(gray);
title(ax,'Frame 1');

hold(ax,'on');
for r = 1:size(centers,1)
    cx = centers(r,1); cy = centers(r,2);
    rectangle(ax, 'Position',[cy-half, cx-half, roi_size, roi_size], 'EdgeColor','r','LineWidth',1.2);
    text(cy, cx, num2str(r), 'Color','y','FontSize',10,'HorizontalAlignment','center','Parent',ax);
end
hold(ax,'off');

uicontrol('Style','slider','Min',1,'Max',totalFrames,'Value',1,...
    'Units','normalized','Position',[0.2 0.02 0.6 0.04],...
    'SliderStep',[1/(totalFrames-1),50/(totalFrames-1)],...
    'Callback', @(src,~) set(hIm,'CData',stack(:,:,round(src.Value))));

