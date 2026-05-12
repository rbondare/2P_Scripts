function [LocomotionCal,ModelParams]=FS_transform_locomotion(LocomotionNonCal,SI_Info,varargin)


if isempty(varargin) || isempty(varargin{1})
    % Try config-based path first
    try
        cfg = get_user_config('AnimalFS0');
        cal_path = cfg.CalibrationFile;
    catch
        cal_path = '';
    end

    if ~isempty(cal_path) && exist(cal_path, 'file')
        load(cal_path, 'ModelParams');
    elseif exist(fullfile(fileparts(mfilename('fullpath')), 'CalibrationModel.mat'), 'file')
        load(fullfile(fileparts(mfilename('fullpath')), 'CalibrationModel.mat'), 'ModelParams');
    else
        [file, path] = uigetfile('*CalibrationModel.mat', 'select calibration file');
        if isequal(file, 0); error('No calibration file selected'); end
        load(fullfile(path, file), 'ModelParams');
    end
else
    ModelParams=varargin{1};
end
if ~isempty(SI_Info)
CalibrationNumber=find(SI_Info.Clockstart>=cat(1,ModelParams(:).Date),1,'last');
else
    CalibrationNumber=numel(ModelParams);
end
if isempty(CalibrationNumber)
    CalibrationNumber=1;%numel(ModelParams);
end
fprintf('Using calibration #%d (date: %s)\n', CalibrationNumber, datestr(ModelParams(CalibrationNumber).Date));
Calibration=ModelParams(CalibrationNumber).Calibration;
if size(Calibration,2)==4; Calibration=[zeros(3,1) Calibration];end
if isstruct(LocomotionNonCal)
    LocomotionCal=LocomotionNonCal;
    for s=1:numel(LocomotionCal)
        d=1;LocomotionCal(s).Forward_mm=Calibration(d,1)+sum(LocomotionNonCal(s).Values(:,1:4).*Calibration(d,2:end),2);
        d=2;LocomotionCal(s).SideRight_mm=Calibration(d,1)+sum(LocomotionNonCal(s).Values(:,1:4).*Calibration(d,2:end),2);
        d=3;LocomotionCal(s).RotRight_deg=Calibration(d,1)+sum(LocomotionNonCal(s).Values(:,1:4).*Calibration(d,2:end),2);
    end
else
    d=1;LocomotionCal.Forward_mm=Calibration(d,1)+sum(LocomotionNonCal(:,1:4).*Calibration(d,2:end),2);
    d=2;LocomotionCal.SideRight_mm=Calibration(d,1)+sum(LocomotionNonCal(:,1:4).*Calibration(d,2:end),2);
    d=3;LocomotionCal.RotRight_deg=Calibration(d,1)+sum(LocomotionNonCal(:,1:4).*Calibration(d,2:end),2);
end










