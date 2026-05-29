%% TEST: Run organize_stimulus_by_parameters with debug output

clear all; close all;

rec1_file = 'Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1115_preprocessed.mat';
rec2_file = 'Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1243_preprocessed.mat';

fprintf('Test 1: Check if files exist\n');
fprintf('  Rec1 exists: %d\n', isfile(rec1_file));
fprintf('  Rec2 exists: %d\n', isfile(rec2_file));

fprintf('\nTest 2: Try to run organize_stimulus_by_parameters\n');
try
    % Run with plane 1 (first/top plane)
    stim_responses = organize_stimulus_by_parameters(rec1_file, rec2_file, 1);
    fprintf('SUCCESS: Function completed\n');
catch ME
    fprintf('ERROR: %s\n', ME.message);
    fprintf('\nStack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('  Line %d in %s\n', ME.stack(i).line, ME.stack(i).name);
    end
end

fprintf('\nTest complete.\n');
