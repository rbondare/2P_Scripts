%% TEST: Run preprocess_and_organize_stimulus with corrected plane indexing

clear all; close all;

fprintf('\n========== TEST: PREPROCESS AND ORGANIZE STIMULUS ==========\n\n');

rec1_file = 'Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1115_preprocessed.mat';
rec2_file = 'Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1243_preprocessed.mat';

% Run the preprocessing script
try
    % Note: Using plane 1 (not 0) since planes are 1-indexed
    preprocess_and_organize_stimulus;
    fprintf('\n✓ PREPROCESSING COMPLETE\n');
catch ME
    fprintf('\n✗ ERROR during preprocessing:\n');
    fprintf('%s\n', ME.message);
    fprintf('\nStack trace:\n');
    for i = 1:min(5, length(ME.stack))
        fprintf('  Line %d in %s\n', ME.stack(i).line, ME.stack(i).name);
    end
end
