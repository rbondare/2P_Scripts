function C = reconstruct_checkers2(s)
% RECONSTRUCT_CHECKERS2  Regenerate the checkers2 white-noise movie from its seed.
%
%   C = reconstruct_checkers2(s)   where s = B.Stimuli(i) with type 'checkers2'
%
% checkers2 never saves its frames: present_stimulus_addnew.m only writes
% stimdata.pattern when project==0 (simulation), and a real experiment runs
% project>0. The pattern is therefore recovered the same way flash onsets are in
% reconstructFlashes.m -- by replaying the saved RNG stream.
%
% GROUND TRUTH (present_stimulus_addnew.m, verified for this stimulus):
%   setup, case 'checkers2' (l.238-308)
%     checkers_num >= size(available_colors,1)  ->  shuffle_pattern_distribution
%     colI       = 1 + interp1(1:nc, 1:nc, linspace(1,nc,checkers_num), 'nearest')
%     patternI   = colI                (index 1 = background, k>=2 = available_colors(k-1))
%     patternLin = colors(patternI,:)  -- a FIXED colour multiset
%     change_pattern_every_n_frames = round((1/ifi)/pattern_change_rate)
%   frame loop, case 'checkers2' (l.1682-1697), reached every frame because
%   update_each_nth_frame == 1:
%     if mod(frame_count, change_pattern_every_n_frames) == 0
%         pattern       = reshape(patternLin(randperm(numel)), ...)   % draw 1
%         current_shift = rand(1,2) .* maxshift                       % draw 2
%     end
%   frame_count is 0-based, so the FIRST frame also triggers a change and the
%   setup pattern is never actually shown.
%
% Why the replay is exact: present_stimulus calls rng(rseed) at the top of every
% trial when reset_rng_each_trial is true (l.1459-1461). checkers2 has trials==1,
% so the generator state at frame 0 is exactly s.randomseed regardless of
% anything the setup phase drew. Within the loop the switch on options.type means
% only the two draws above are consumed -- every other rand/randperm in that loop
% body belongs to a different stimulus case.
%
% Returns C with
%   .contrast        [n_rows x n_cols x n_patterns] single, background-relative
%                    contrast (-1 = black, +1 = bright, 0 = background)
%   .onset_t         [1 x n_patterns] pattern onset, same clock as TimeStimulusFrame
%   .shift           [n_patterns x 2] sub-element jitter (row, col), in ELEMENTS
%   .vis_rows/.vis_cols  indices of the on-screen sub-rectangle
%   .az_deg/.el_deg  centre of each visible column/row, degrees, 0 = screen centre
%   .meta            parameters echoed for QC

    assert(strcmp(s.type, 'checkers2'), 'not a checkers2 stimulus: %s', s.type);
    sp = s.specParams;

    % ---- geometry ---------------------------------------------------------
    ne = double(sp.checkers_elements);          % [n_rows n_cols], border included
    n_el = prod(ne);
    assert(~isfield(sp, 'restrict_area') && ~isfield(sp, 'define_area'), ...
        ['restrict_area/define_area set: the visible sub-rectangle and the ' ...
         'degree axes below assume the full-screen path.']);

    % base_source_rect = [0 0 fliplr(ne)-2], then +1 in the draw call, so the
    % visible window is one element in from each edge; the spare border is what
    % the sub-element jitter slides into.
    vis_rows = 2:(ne(1)-1);
    vis_cols = 2:(ne(2)-1);

    % ---- colour -> background-relative contrast ---------------------------
    % relativeintensity kills the red channel on this rig ([0 1 1]), so weight
    % the channels by it rather than taking a plain mean.
    relint = double(s.relativeintensity(:));
    if all(relint == 0); relint = ones(size(relint)); end
    lum = @(c) (double(c) * relint) ./ sum(relint);

    avail   = sp.available_colors;               % [n_col x 3]
    bg_lum  = lum(sp.background_col);
    col_lum = [bg_lum; lum(avail)];              % index 1 = background
    if bg_lum > 0
        contrast_by_index = (col_lum - bg_lum) ./ bg_lum;
    else
        contrast_by_index = col_lum;             % no background to modulate around
    end
    m = max(abs(contrast_by_index));
    if m > 0; contrast_by_index = contrast_by_index ./ m; end

    % ---- the fixed colour multiset (setup phase) --------------------------
    n_col = size(avail, 1);
    checkers_num = min(double(sp.checkers_num), n_el);
    assert(checkers_num >= n_col, ...
        ['checkers_num (%d) < n available colours (%d): the stimulus took the ' ...
         'non-shuffle branch, which draws a different RNG sequence and is not ' ...
         'implemented here.'], checkers_num, n_col);

    patternI = ones(n_el, 1);
    if n_col > 1
        colI = 1 + interp1(1:n_col, 1:n_col, linspace(1, n_col, checkers_num), 'nearest');
    else
        colI = 1 + ones(1, checkers_num);
    end
    patternI(1:numel(colI)) = colI;
    contrast_lin = contrast_by_index(patternI);        % [n_el x 1]

    % ---- when the pattern changes -----------------------------------------
    ifi = double(s.ifi);
    change_every = round((1/ifi) / double(sp.pattern_change_rate));
    assert(s.update_each_nth_frame == 1, ...
        ['update_each_nth_frame = %d. The pattern-change test sits inside ' ...
         '"if mod(frame_count, update_each_nth_frame)==0", so the effective ' ...
         'period is lcm(%d,%d), not %d -- not handled.'], ...
        s.update_each_nth_frame, s.update_each_nth_frame, change_every, change_every);

    n_frames  = numel(s.TimeStimulusFrame);
    frame_ids = 0:(n_frames - 1);                      % frame_count is 0-based
    pat_frame = frame_ids(mod(frame_ids, change_every) == 0);
    n_pat     = numel(pat_frame);

    if strcmp(sp.subshift, 'pixel'); maxshift = [1 1]; else; maxshift = [0 0]; end
    maxshift(fliplr(ne) == 1) = 0;

    % drift would make the pattern circshift within a presentation; this rig
    % records drift_speed = [0 0], and the reconstruction below assumes it.
    has_drift = isfield(sp, 'drift_speed') && any(abs(sp.drift_speed) > 0);
    assert(~has_drift, 'drift_speed is non-zero; within-presentation circshift not implemented.');

    % ---- replay -----------------------------------------------------------
    rng(s.randomseed);
    contrast = zeros([ne, n_pat], 'single');
    shift    = zeros(n_pat, 2);
    for p = 1:n_pat
        idx = randperm(n_el);                          % draw 1
        contrast(:, :, p) = single(reshape(contrast_lin(idx), ne));
        shift(p, :) = rand(1, 2) .* maxshift;          % draw 2
    end

    C = struct();
    C.contrast = contrast;
    C.onset_t  = s.TimeStimulusFrame(pat_frame + 1);   % ->1-based
    C.shift    = shift;
    C.vis_rows = vis_rows;
    C.vis_cols = vis_cols;

    % Elements are stretched to fill the screen exactly, so degrees per element
    % come from the screen extent, not from the nominal checker_size.
    deg_col = double(s.full_screen_az) / numel(vis_cols);
    deg_row = double(s.full_screen_el) / numel(vis_rows);
    C.az_deg = ((1:numel(vis_cols)) - 0.5) * deg_col - double(s.full_screen_az)/2;
    C.el_deg = ((1:numel(vis_rows)) - 0.5) * deg_row - double(s.full_screen_el)/2;

    C.meta = struct('n_patterns', n_pat, 'change_every_frames', change_every, ...
        'pattern_change_rate', double(sp.pattern_change_rate), 'ifi', ifi, ...
        'checkers_elements', ne, 'checkers_num', checkers_num, ...
        'deg_per_col', deg_col, 'deg_per_row', deg_row, ...
        'contrast_levels', unique(contrast_lin(:))', ...
        'maxshift', maxshift, 'n_frames', n_frames, ...
        'block_t0', s.TimeStimulusFrame(1), 'block_t1', s.TimeStimulusFrame(end));
end
