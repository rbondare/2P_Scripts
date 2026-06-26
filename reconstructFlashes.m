function [flash_onsets, Object_in_frame, Frametime] = reconstructFlashes(options_flashes, ~)
% Function used to reconstruct sparse_local_global_flashes onset times
% Used for visualization purposes and to recover trial-onset timing
%
% This is a faithful replay of the 'sparse_local_global_flashes' case in
% present_stimulus_addnew.m (the actual PsychToolbox stimulus-presentation
% script), not a guess -- confirmed against that source directly. It
% reproduces Object_in_frame/Frametime exactly as that function computes
% them, by reseeding to the same recorded randomseed and drawing in the
% same order: one rand() for the inter-flash gap, then one randi() for
% which flash "type" to show (always 1 here since specParams.types has
% only one element, but the draw still advances the RNG stream, so it
% must be replayed even though its result doesn't vary).
%
% INPUT:
% options_flashes: the Stimuli() struct entry for sparse_local_global_flashes
%                   (needs .specParams.interval, .specParams.duration,
%                   .specParams.types, .stimulus_trial_t, .ifi, .randomseed)
% fs: unused (kept for interface parity with reconstructChirp.m); pass [] or omit
%
% OUTPUT:
% flash_onsets: vector of flash onset times (seconds, relative to
%               stimulus start, i.e. relative to TimeStimulusFrame(1))
% Object_in_frame, Frametime: same outputs present_stimulus_addnew.m
%               stores in stimdata.Object_in_frame for this stimulus
%

debug = false; % if true does plotting

ifi = options_flashes.ifi(1);

% Replay the exact RNG state recorded at presentation time -- same call
% present_stimulus_addnew.m makes at its top
rng(options_flashes.randomseed);

Object_in_frame = zeros(200 + ceil(options_flashes.stimulus_trial_t/ifi), 1);
Frametime       = (1:ceil(options_flashes.stimulus_trial_t/ifi)+200) * ifi;
t = 0;
while t <= Frametime(end)
    t = t + rand(1)*diff(options_flashes.specParams.interval) + options_flashes.specParams.interval(1);
    Object_in_frame(Frametime >= t & Frametime < t + options_flashes.specParams.duration) = ...
        randi(size(options_flashes.specParams.types, 2), 1);
end

% Derive onset times from Object_in_frame: first frame of each contiguous
% nonzero run
is_on = Object_in_frame > 0;
onset_mask = diff([0; is_on]) == 1;
flash_onsets = Frametime(onset_mask)';
flash_onsets = flash_onsets(flash_onsets <= options_flashes.stimulus_trial_t);

if debug == true
   figure;
   plot(Frametime, Object_in_frame > 0)
   xlim([0 options_flashes.stimulus_trial_t])
   ylim([-0.5 1.5])
   xlabel('Time (s)'); ylabel('Flash on/off');
   title(sprintf('Reconstructed flashes (n = %d)', numel(flash_onsets)));
end

end
