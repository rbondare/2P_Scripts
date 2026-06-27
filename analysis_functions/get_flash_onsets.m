function get_onsets = get_flash_onsets(Stimuli)
    % Flash onsets are reconstructed by replaying the RNG draws used at
    % presentation time (reconstructFlashes.m), anchored directly to
    % TimeStimulusFrame(1), which is camera-clock-correct on its own
    % (RedFrameSynchronized=true) -- no manual anchor-correction needed.
    flash_idx = find(strcmp({Stimuli.type}, 'sparse_local_global_flashes'), 1);
    if isempty(flash_idx); get_onsets = []; return; end
    s = Stimuli(flash_idx);
    if isfield(s, 'RedFrameSynchronized') && ~s.RedFrameSynchronized
        warning('Flashes: RedFrameSynchronized=false -- timing fell back to PC clock, expect an offset');
    end
    flash_onsets_rel = reconstructFlashes(s);
    get_onsets = s.TimeStimulusFrame(1) + flash_onsets_rel;
end
