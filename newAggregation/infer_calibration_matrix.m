% Infer the calibration matrix from an already-loaded preprocessed file.
% Assumes the variable 'data' is already in your workspace (matfile object or struct).
%
% Run as a script: just press Run or type 'infer_calibration_matrix' in the command window.

NonCal = data.LocomotionNonCal;
Cal    = data.LocomotionCal;

% Stack all epochs
all_X   = [];
all_fwd = [];
all_sid = [];
all_rot = [];

for s = 1:numel(NonCal)
    X   = NonCal(s).Values(:, 1:4);
    fwd = Cal(s).Forward_mm(:);
    sid = Cal(s).SideRight_mm(:);
    rot = Cal(s).RotRight_deg(:);

    N = min([size(X,1), numel(fwd), numel(sid), numel(rot)]);
    X=X(1:N,:); fwd=fwd(1:N); sid=sid(1:N); rot=rot(1:N);

    valid = ~any(isnan(X),2) & ~isnan(fwd) & ~isnan(sid) & ~isnan(rot);
    all_X   = [all_X;   X(valid,:)];
    all_fwd = [all_fwd; fwd(valid)];
    all_sid = [all_sid; sid(valid)];
    all_rot = [all_rot; rot(valid)];
end

fprintf('Using %d data points from %d epochs\n', size(all_X,1), numel(NonCal));

% Design matrix: [intercept | sensor1 | sensor2 | sensor3 | sensor4]
A = [ones(size(all_X,1),1), all_X];

% Solve each output row (exact inversion — no noise in this transform)
row1 = A \ all_fwd;
row2 = A \ all_sid;
row3 = A \ all_rot;

Calibration = [row1'; row2'; row3'];  % 3x5

% Residuals — should be ~1e-10 if this is the exact calibration
res_fwd = max(abs(A*row1 - all_fwd));
res_sid = max(abs(A*row2 - all_sid));
res_rot = max(abs(A*row3 - all_rot));

fprintf('\n=== Inferred Calibration (3x5) ===\n');
fprintf('           Intercept     Sensor1      Sensor2      Sensor3      Sensor4\n');
labels = {'Forward_mm  ','SideRight_mm','RotRight_deg'};
for d = 1:3
    fprintf('%s  %11.6f  %11.6f  %11.6f  %11.6f  %11.6f\n', ...
        labels{d}, Calibration(d,1), Calibration(d,2), ...
        Calibration(d,3), Calibration(d,4), Calibration(d,5));
end

fprintf('\nMax residuals (should be ~1e-10):\n');
fprintf('  Forward_mm:    %.2e\n', res_fwd);
fprintf('  SideRight_mm:  %.2e\n', res_sid);
fprintf('  RotRight_deg:  %.2e\n', res_rot);

% Compare against CalibrationModel.mat if available
cal_file = fullfile(fileparts(which('transform_locomotion')), 'CalibrationModel.mat');
if exist(cal_file, 'file')
    tmp = load(cal_file, 'ModelParams');
    fprintf('\n=== Match against CalibrationModel.mat entries ===\n');
    for k = 1:numel(tmp.ModelParams)
        C_ref = tmp.ModelParams(k).Calibration;
        if size(C_ref,2) == 4; C_ref = [zeros(3,1), C_ref]; end
        d = norm(Calibration - C_ref, 'fro');
        match = '';
        if d < 1e-4; match = '  <-- MATCH'; end
        fprintf('  Entry %d (%s): distance = %.2e%s\n', k, datestr(tmp.ModelParams(k).Date), d, match);
    end
end
