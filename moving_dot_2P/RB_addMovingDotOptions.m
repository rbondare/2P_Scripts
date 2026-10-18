% RB_addMovingDotOptions  Append the 'moving_dot_4dir' entries to default_options_mouse.mat
%
% Run ONCE on the rig PC with MouseTrainingScripts as current folder (or on the
% path). It adds two entries to the end of the options array, so no existing index
% changes:
%     'dot 4dir, grey background'      specParams.checkers = 0
%     'dot 4dir, checkers background'  specParams.checkers = 1
% and prints their indices - those are the numbers to enter in the stimulus
% selection of the GUI. A time-stamped backup of the mat file is written first.
% Edit the specParams / timing values below to change the defaults, or edit them
% on the options struct after loading it in FS_potentScript.
%
% specParams
%   dot_speed_deg_s, dot_radius_deg  dot speed and radius
%   dot_start_t, dot_stop_t          s of rest before / after each sweep
%   visual_area_10_20                [bottom top left right] on a 10 x 20 grid of the canvas
%   dot_color, dot_background_col    RGB 0..1, multiplied by maxintensity*relativeintensity
%   checkers                         0 = grey background, 1 = moving checkers behind the dot
%   checker_ncol                     checker columns across the canvas
%   checker_move_ratio               random shifts of the pattern are multiples of 1/move_ratio checker
%   checker_nframes                  new pattern every n frames (1 = every frame)
%   checker_fraction                 probability that a checker has checker_max_col
%   checker_min_col, checker_max_col RGB 0..1 of the two checker values

optfile = which('default_options_mouse.mat');
if isempty(optfile)
    error('default_options_mouse.mat is not on the path - cd to MouseTrainingScripts first');
end
S = load(optfile, 'options');
options = S.options;

done = find(strcmp({options.type}, 'moving_dot_4dir'));
if ~isempty(done)
    error('moving_dot_4dir entries already exist (indices %s) - nothing changed', mat2str(done));
end

backup = strrep(optfile, '.mat', ['_backup_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) '.mat']);
copyfile(optfile, backup);
fprintf('backup written: %s\n', backup);

% a moving_bar entry already has every field the pipeline expects
new = options(find(strcmp({options.type}, 'moving_bar'), 1));
new.type = 'moving_dot_4dir';
new.stimulus_trial_t = 10;          % placeholder, present_stimulus_addnew sets it from the trajectories
new.trials = 3;                     % repetitions of the whole area
new.update_each_nth_frame = 1;      % must stay 1
new.reset_rng_each_trial = 0;       % 0 = fresh checker noise in every repetition
new.waitrange_between_trials = [0 0];
new.wait_init = 5;
new.wait_exit = 10;
new.randomseed = [];

sp = struct( ...
    'dot_speed_deg_s',   50, ...
    'dot_radius_deg',    2.5, ...
    'dot_start_t',       0.5, ...
    'dot_stop_t',        2, ...
    'visual_area_10_20', [0 5 10 20], ...
    'dot_color',         [0 0 0], ...
    'dot_background_col',[0 0.5 0.5], ...
    'checkers',          0, ...
    'checker_ncol',      30, ...
    'checker_move_ratio',10, ...
    'checker_nframes',   1, ...
    'checker_fraction',  0.5, ...
    'checker_min_col',   [0 0 0], ...
    'checker_max_col',   [0 0.5 0.5]);

first = numel(options)+1;

new.specParams = sp;
new.addnote = 'dot 4dir, grey background';
options(end+1) = new;

new.specParams.checkers = 1;
new.addnote = 'dot 4dir, checkers background';
options(end+1) = new;

save(optfile, 'options');
fprintf('added options(%d) "%s" and options(%d) "%s" to %s\n', ...
    first, options(first).addnote, first+1, options(first+1).addnote, optfile);
