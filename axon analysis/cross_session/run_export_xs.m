% usage: set part (1..nparts) and nparts before running, e.g. part = 1; nparts = 3; run_export_xs
sessions = {'AnimalRB24_261021_1630', 'AnimalRB24_261023_1527', 'AnimalRB24_261023_1538', 'AnimalRB24_261027_1456', 'AnimalRB24_261027_1534', 'AnimalRB24_261029_1752', ...
            'AnimalRB24_261030_1744', 'AnimalRB24_261104_1640', 'AnimalRB24_261104_1651', 'AnimalRB24_261105_1423', 'AnimalRB24_261105_1520', 'AnimalRB24_261106_1317', 'AnimalRB24_261106_1529', ...
            'AnimalRB23_261023_1642', 'AnimalRB23_261027_1306', 'AnimalRB23_261027_1423', 'AnimalRB23_261029_1704', 'AnimalRB22_261023_1427', 'AnimalRB22_261027_1138'};
addpath('C:\Users\rbondarenko\projects\2P_Scripts\axon analysis\cross_session');
for k = part:nparts:numel(sessions)
    t0 = tic;
    try
        export_xs(sessions{k});
    catch ME
        fprintf(2, 'FAILED %s: %s\n', sessions{k}, getReport(ME));
    end
    fprintf('  (%d/%d, %.1f min)\n', k, numel(sessions), toc(t0) / 60);
end
fprintf('EXPORT_DONE part %d\n', part);
