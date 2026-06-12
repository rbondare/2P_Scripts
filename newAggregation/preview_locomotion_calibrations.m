function preview_locomotion_calibrations(input_path)
%PREVIEW_LOCOMOTION_CALIBRATIONS Compare all calibration models before aggregation.
%
%   preview_locomotion_calibrations(input_path)
%
%   input_path can be EITHER:
%     (a) The folder containing stim*.mat files directly, e.g.:
%           'Z:\joeschgrp\Group Members\Rima\Stimulus_Axons\AnimalRB6\AnimalRB6_250302_1711'
%     (b) The 2P recording folder (same as save_preprocessed input), in which
%         case the stim path is derived automatically via get_user_config, e.g.:
%           'Z:\...\DATA_2P\AnimalRB6_250302_1711\'
%
%   Plots locomotion under every CalibrationModel.mat entry so you can pick
%   the right one before running aggregation. Auto-selected model is highlighted.

if nargin < 1 || isempty(input_path)
    input_path = uigetdir('', 'Select stim folder or recording folder');
    if isequal(input_path, 0); return; end
end

basedir = input_path;
if basedir(end) ~= filesep; basedir(end+1) = filesep; end

parts = strsplit(basedir, filesep);
parts = parts(~cellfun(@isempty, parts));
RecordingName = parts{end};
AnimalName = strsplit(RecordingName, '_'); AnimalName = AnimalName{1};

fprintf('Recording : %s\n', RecordingName);
fprintf('Animal    : %s\n', AnimalName);

%% Locate stim files — check given folder first, then derive via config
StimFiles = dir(fullfile(basedir, 'stim*.mat'));
if ~isempty(StimFiles)
    StimPath = basedir;
    fprintf('Stim path : %s (given folder)\n', StimPath);
else
    % Derive from get_user_config (works for non-axon / standard pipeline)
    config   = get_user_config(AnimalName);
    StimPath = fullfile(config.StimBasePath, AnimalName, RecordingName, filesep);
    StimFiles = dir(fullfile(StimPath, 'stim*.mat'));
    if isempty(StimFiles)
        [parent, recname] = fileparts(StimPath(1:end-1));
        StimFiles = dir(fullfile(parent, [recname '*'], 'stim*.mat'));
    end
    if ~isempty(StimFiles)
        fprintf('Stim path : %s (from config)\n', StimPath);
    end
end

if isempty(StimFiles)
    error('No stim*.mat files found.\nTried: %s\nPass the stim folder directly if the path is non-standard.', basedir);
end
fprintf('Found %d stim file(s)\n', numel(StimFiles));

%% Extract raw BallCamData values from stim files
BallValues = {};   % one cell per stimulus epoch with an N×4 matrix
epoch_labels = {};

for F = 1:numel(StimFiles)
    stim_data = load(fullfile(StimFiles(F).folder, StimFiles(F).name));
    if ~isfield(stim_data, 'BallCamData')
        fprintf('  Stim %d: no BallCamData field, skipping\n', F);
        continue;
    end

    raw = stim_data.BallCamData;   % cell array of readings (each is 1×5)
    if iscell(raw)
        good = cellfun(@numel, raw) == 5;
        if sum(good) == 0; continue; end
        V = cell2mat(raw(good)');  % N×5
    elseif isnumeric(raw)
        V = raw;
    else
        continue;
    end

    BallValues{end+1} = V(:, 1:4);  % only first 4 channels used by calibration
    epoch_labels{end+1} = sprintf('Stim%d', F);
    fprintf('  Stim %d: %d ball samples\n', F, size(V,1));
end

if isempty(BallValues)
    error('No valid BallCamData found in stim files.');
end

%% Load all calibration models
cal_file = fullfile(fileparts(mfilename('fullpath')), 'CalibrationModel.mat');
if ~exist(cal_file, 'file')
    [f, p] = uigetfile('*CalibrationModel.mat', 'Select CalibrationModel.mat');
    if isequal(f,0); error('No calibration file selected'); end
    cal_file = fullfile(p, f);
end
tmp = load(cal_file, 'ModelParams');
ModelParams = tmp.ModelParams;
nModels = numel(ModelParams);
fprintf('\nFound %d calibration model(s) in CalibrationModel.mat\n', nModels);

%% Determine auto-selected calibration based on recording date
% Parse date from folder name: AnimalXX_YYMMDD_HHMM or AnimalXX_YYYYMMDD_HHMM
rec_date = [];
tokens = regexp(RecordingName, '_(\d{6,8})_', 'tokens');
if ~isempty(tokens)
    datestr_raw = tokens{1}{1};
    try
        if numel(datestr_raw) == 6        % YYMMDD
            rec_date = datetime(['20' datestr_raw], 'InputFormat', 'yyyyMMdd');
        else                               % YYYYMMDD
            rec_date = datetime(datestr_raw, 'InputFormat', 'yyyyMMdd');
        end
    catch
        rec_date = [];
    end
end

auto_idx = nModels;  % default: last entry
if ~isempty(rec_date)
    model_dates = cat(1, ModelParams(:).Date);
    mask = rec_date >= model_dates;
    if any(mask)
        auto_idx = find(mask, 1, 'last');
    else
        auto_idx = 1;
    end
end
fprintf('Auto-selected calibration for this recording: #%d (%s)\n', ...
    auto_idx, datestr(ModelParams(auto_idx).Date));

%% Apply each calibration model to all epochs
% Results: cal_results{model_idx}{epoch_idx} = struct with Forward_mm etc.
cal_results = cell(nModels, 1);
for m = 1:nModels
    C = ModelParams(m).Calibration;
    if size(C,2) == 4; C = [zeros(3,1), C]; end  % add zero intercept if missing

    cal_results{m} = cell(numel(BallValues), 1);
    for ep = 1:numel(BallValues)
        X = BallValues{ep};   % N×4
        A = [ones(size(X,1),1), X];  % N×5
        cal_results{m}{ep}.Forward_mm   = A * C(1,:)';
        cal_results{m}{ep}.SideRight_mm = A * C(2,:)';
        cal_results{m}{ep}.RotRight_deg = A * C(3,:)';
    end
end

%% Build legend entries
colors = lines(nModels);
legend_strs = cell(nModels, 1);
for m = 1:nModels
    marker = '';
    if m == auto_idx; marker = '  <<< auto-selected'; end
    legend_strs{m} = sprintf('#%d  %s%s', m, datestr(ModelParams(m).Date, 'yyyy-mmm-dd'), marker);
end

%% Plot — one figure per epoch, 3 rows (Forward / SideRight / RotRight)
ball_fs    = 66;   % ball camera acquisition rate (Hz)
dim_names  = {'Forward\_mm', 'SideRight\_mm', 'RotRight\_deg'};
dim_fields = {'Forward_mm', 'SideRight_mm', 'RotRight_deg'};
dim_units  = {'mm', 'mm', 'deg'};

for ep = 1:numel(BallValues)
    N = size(BallValues{ep}, 1);
    t = (0:N-1)' / ball_fs;   % seconds

    fig = figure('Name', sprintf('Locomotion calibration preview — %s (%s)', ...
        RecordingName, epoch_labels{ep}), ...
        'NumberTitle', 'off', 'Color', 'w', ...
        'Position', [100, 100, 1100, 700]);

    for d = 1:3
        ax = subplot(3, 1, d);
        hold on;

        for m = 1:nModels
            lw = 1.2;
            ls = '-';
            if m == auto_idx; lw = 2.5; ls = '-'; end

            plot(t, cal_results{m}{ep}.(dim_fields{d}), ls, ...
                'Color', colors(m,:), 'LineWidth', lw, ...
                'DisplayName', legend_strs{m});
        end

        hold off;
        ylabel([dim_names{d} ' (' dim_units{d} ')'], 'FontSize', 10);
        if d == 1
            title(sprintf('%s — %s', strrep(RecordingName,'_','\_'), epoch_labels{ep}), ...
                'FontSize', 11, 'FontWeight', 'bold');
        end
        if d == 3
            xlabel('Time (s)', 'FontSize', 10);
        end
        grid on; box off;
        if d == 1
            legend('show', 'Location', 'best', 'FontSize', 8);
        end
        xlim([0, (N-1)/ball_fs]);
    end

    linkaxes(findobj(fig, 'Type', 'axes'), 'x');
end

%% Plot — raw sensors (separate subplots) vs last 5 calibrated models (one figure per epoch)
sensor_colors = lines(4);

last5_idx = max(1, nModels - 4) : nModels;   % last 5 (or all if fewer)
n5        = numel(last5_idx);
muted_pal = [
    0.25  0.50  0.70;   % steel blue
    0.40  0.65  0.42;   % sage green
    0.65  0.42  0.60;   % dusty mauve
    0.78  0.57  0.28;   % warm tan
    0.50  0.50  0.50;   % mid grey
];
colors5   = muted_pal(mod((0:n5-1)', size(muted_pal, 1)) + 1, :);
legend5   = legend_strs(last5_idx);
cal_alpha = 0.55;   % base transparency for calibration lines

for ep = 1:numel(BallValues)
    N = size(BallValues{ep}, 1);
    t = (0:N-1)' / ball_fs;   % seconds

    figure('Name', sprintf('Raw sensors vs Calibrated — %s (%s)', ...
        RecordingName, epoch_labels{ep}), ...
        'NumberTitle', 'off', 'Color', 'w', ...
        'Position', [150, 30, 1300, 820]);

    sgtitle(sprintf('%s  –  %s', strrep(RecordingName, '_', '\_'), epoch_labels{ep}), ...
        'FontSize', 12, 'FontWeight', 'bold');

    % ---- Left column (col 1): one subplot per raw sensor channel ----
    ax_raw = gobjects(4, 1);
    for ch = 1:4
        ax_raw(ch) = subplot(4, 2, (ch - 1) * 2 + 1);
        plot(t, BallValues{ep}(:, ch), '-', ...
            'Color', sensor_colors(ch, :), 'LineWidth', 1.0);
        ylabel(sprintf('Ch%d (ADC)', ch), 'FontSize', 9);
        if ch == 1
            title('Raw sensor channels', 'FontSize', 10, 'FontWeight', 'bold');
        end
        if ch == 4
            xlabel('Time (s)', 'FontSize', 9);
        end
        grid on; box off;
        xlim([0, (N-1)/ball_fs]);
    end

    % ---- Right column (col 2, rows 1-3): calibrated outputs, last 5 models ----
    ax_cal = gobjects(3, 1);
    for d = 1:3
        ax_cal(d) = subplot(4, 2, (d - 1) * 2 + 2);
        hold on;
        % non-auto models first so auto-selected renders on top
        for mi = 1:numel(last5_idx)
            m = last5_idx(mi);
            if m == auto_idx; continue; end
            plot(t, cal_results{m}{ep}.(dim_fields{d}), '-', ...
                'Color', [colors5(mi, :), cal_alpha], 'LineWidth', 1.2, ...
                'DisplayName', legend5{mi});
        end
        % auto-selected last → sits on top
        auto_mi = find(last5_idx == auto_idx, 1);
        if ~isempty(auto_mi)
            plot(t, cal_results{auto_idx}{ep}.(dim_fields{d}), '-', ...
                'Color', [colors5(auto_mi, :), 0.90], 'LineWidth', 2.5, ...
                'DisplayName', legend5{auto_mi});
        end
        hold off;
        ylabel([dim_names{d} ' (' dim_units{d} ')'], 'FontSize', 9);
        if d == 1
            title('Calibrated (last 5 models)', 'FontSize', 10, 'FontWeight', 'bold');
            legend('show', 'Location', 'best', 'FontSize', 7);
        end
        if d == 3
            xlabel('Time (s)', 'FontSize', 9);
        end
        grid on; box off;
        xlim([0, (N-1)/ball_fs]);
    end

    linkaxes([ax_raw; ax_cal], 'x');
end

%% Plot — per-second sums (area under curve per 1-s bin, one figure per epoch)
% Sums all frame values within each 1-second window → displacement / ADC per second

for ep = 1:numel(BallValues)
    N      = size(BallValues{ep}, 1);
    n_sec  = floor(N / ball_fs);
    if n_sec < 1; continue; end
    t_sec  = (0:n_sec-1)';

    % Helper: sum-bin a column vector into n_sec × 1
    bin1s = @(x) sum(reshape(x(1:n_sec*ball_fs), ball_fs, n_sec), 1)';

    % Bin raw channels
    raw_sec = zeros(n_sec, 4);
    for ch = 1:4
        raw_sec(:, ch) = bin1s(BallValues{ep}(:, ch));
    end

    % Bin calibrated outputs for last 5 models
    cal_sec = cell(numel(last5_idx), numel(dim_fields));
    for mi = 1:numel(last5_idx)
        m = last5_idx(mi);
        for d = 1:3
            cal_sec{mi, d} = bin1s(cal_results{m}{ep}.(dim_fields{d}));
        end
    end

    figure('Name', sprintf('Per-second sums — %s (%s)', RecordingName, epoch_labels{ep}), ...
        'NumberTitle', 'off', 'Color', 'w', 'Position', [200, 30, 1300, 820]);

    sgtitle(sprintf('%s  –  %s  |  per-second sums', ...
        strrep(RecordingName, '_', '\_'), epoch_labels{ep}), ...
        'FontSize', 12, 'FontWeight', 'bold');

    % Left column: raw sensor channels (summed per second)
    ax_raw2 = gobjects(4, 1);
    for ch = 1:4
        ax_raw2(ch) = subplot(4, 2, (ch - 1) * 2 + 1);
        bar(t_sec, raw_sec(:, ch), 1, 'FaceColor', sensor_colors(ch, :), 'EdgeColor', 'none');
        ylabel(sprintf('Ch%d (ADC·s)', ch), 'FontSize', 9);
        if ch == 1
            title('Raw sensor channels', 'FontSize', 10, 'FontWeight', 'bold');
        end
        if ch == 4
            xlabel('Time (s)', 'FontSize', 9);
        end
        grid on; box off;
        xlim([-0.5, n_sec - 0.5]);
    end

    % Right column: calibrated outputs (summed per second), last 5 models
    ax_cal2 = gobjects(3, 1);
    for d = 1:3
        ax_cal2(d) = subplot(4, 2, (d - 1) * 2 + 2);
        hold on;
        for mi = 1:numel(last5_idx)
            m  = last5_idx(mi);
            if m == auto_idx; continue; end
            plot(t_sec, cal_sec{mi, d}, '-', ...
                'Color', [colors5(mi, :), cal_alpha], 'LineWidth', 1.2, ...
                'DisplayName', legend5{mi});
        end
        auto_mi = find(last5_idx == auto_idx, 1);
        if ~isempty(auto_mi)
            plot(t_sec, cal_sec{auto_mi, d}, '-', ...
                'Color', [colors5(auto_mi, :), 0.90], 'LineWidth', 2.5, ...
                'DisplayName', legend5{auto_mi});
        end
        hold off;
        ylabel([dim_names{d} ' (' dim_units{d} '·s)'], 'FontSize', 9);
        if d == 1
            title('Calibrated (last 5 models)', 'FontSize', 10, 'FontWeight', 'bold');
            legend('show', 'Location', 'best', 'FontSize', 7);
        end
        if d == 3
            xlabel('Time (s)', 'FontSize', 9);
        end
        grid on; box off;
        xlim([-0.5, n_sec - 0.5]);
    end

    linkaxes([ax_raw2; ax_cal2], 'x');
end

%% Print calibration coefficient table
fprintf('\n=== Calibration coefficients ===\n');
fprintf('%-5s  %-14s  %-10s  %-10s  %-10s  %-10s  %-10s\n', ...
    'Entry', 'Date', 'Intercept', 'Sensor1', 'Sensor2', 'Sensor3', 'Sensor4');
dim_short = {'Fwd','Side','Rot'};
for m = 1:nModels
    C = ModelParams(m).Calibration;
    if size(C,2) == 4; C = [zeros(3,1), C]; end
    marker = ''; if m == auto_idx; marker = ' <<<'; end
    for d = 1:3
        if d == 1
            entry_label = sprintf('#%d%s', m, marker);
            date_label  = datestr(ModelParams(m).Date, 'yyyy-mm-dd');
        else
            entry_label = '';
            date_label  = dim_short{d};
        end
        fprintf('%-5s  %-14s  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f\n', ...
            entry_label, date_label, C(d,1), C(d,2), C(d,3), C(d,4), C(d,5));
    end
    fprintf('\n');
end
