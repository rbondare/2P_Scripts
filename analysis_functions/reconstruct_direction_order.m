function order = reconstruct_direction_order(Stimopts, orientation_field)
    % Direction order is NOT saved anywhere in the stim file (only the
    % available orientation list, e.g. specParams.bar_orientations, is
    % saved -- never the realized per-trial order). It is only recoverable
    % by replaying the RNG from the saved seed, matching
    % present_stimulus_addnew.m's project>0 branch exactly (confirmed via
    % this lab's data: SetUp==1 -> project==1):
    %   for R=1:trials
    %       new_trial=true; if reset_rng_each_trial; rng(rseed); end
    %       ... case 'moving_bar'/'grating': if new_trial
    %             orientations = specParams.X_orientations(randperm(numel(...)));
    %   end
    % Zero other random draws happen between the rng() reset and this
    % randperm() call, and reset_rng_each_trial=true for both stimulus
    % types in this data, so EVERY outer trial repeat gets the IDENTICAL
    % order -- one randperm() call covers all repeats. Validated
    % empirically (not just assumed): same-subtrial-position responses
    % correlate ~0.20 on average across moving_bar's 3 repeats, vs ~0.00
    % for a shifted-position control -- see verify_direction_reconstruction.m.
    rng(Stimopts.randomseed);
    orientations = Stimopts.specParams.(orientation_field);
    order = orientations(randperm(numel(orientations)));
end
