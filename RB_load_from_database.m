load('RB_Opto_Database_latest.mat');  

recordings = struct();

for i = 1:height(RB_Opto_Database.DB)
    animal = RB_Opto_Database.DB.MouseID{i};
    date = RB_Opto_Database.DB.Date{i};
    filepath = RB_Opto_Database.DB.FilePath{i};
    condition = RB_Opto_Database.DB.condition{i};
   
    if isnumeric(date)
        date = num2str(date);
    end
    month_day = date(5:8);
    
    [~, filename, ~] = fileparts(filepath);
    time_pattern = regexp(filename, '_(\d{6})$', 'tokens');
    if ~isempty(time_pattern)
        hour_min = time_pattern{1}{1}(1:4);
    else
        hour_min = '0000';
    end
    
    new_name = sprintf('%s_%s_%s_%s', animal, month_day, hour_min, condition);
    
    try
        recordings.(new_name) = load(filepath);
        fprintf('Loaded: %s\n', new_name);
    catch ME
        fprintf('Error:  %s\n', ME.message);
    end
end