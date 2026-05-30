%% QUICK TEST: organize_stimulus_by_parameters with plane 1

clear all; close all;

fprintf('\n========== QUICK TEST: organize_stimulus_by_parameters ==========\n\n');

rec1_file = 'Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1115_preprocessed.mat';
rec2_file = 'Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1243_preprocessed.mat';

try
    fprintf('Starting organization...\n');
    
    % Call function with plane 1 (not 0!)
    stim_responses = organize_stimulus_by_parameters(rec1_file, rec2_file, 1);
    
    fprintf('\n✓ SUCCESS! Function completed.\n\n');
    
    % Show what was organized
    fprintf('Results:\n');
    resp_fields = fieldnames(stim_responses);
    for i = 1:length(resp_fields)
        field = resp_fields{i};
        fprintf('  - %s\n', field);
    end
    
catch ME
    fprintf('\n✗ ERROR:\n');
    fprintf('%s\n', ME.message);
    fprintf('\nFull stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('  [%d] Line %d in %s\n', i, ME.stack(i).line, ME.stack(i).name);
    end
end

fprintf('\nTest complete.\n');
