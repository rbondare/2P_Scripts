function [LocomotionCal,ModelParams]=FS_transform_locomotion(LocomotionNonCal,SI_Info,varargin)


if isempty(varargin) || isempty(varargin{1})
  
    if exist('A:\Florian\BallCalibration\CalibrationModel.mat','file')
        load('A:\Florian\BallCalibration\CalibrationModel.mat','ModelParams');
    elseif exist('S:\fs3-joeschgrp\Toni\BallCalibration\CalibrationModel.mat','file')
        load('S:\fs3-joeschgrp\Toni\BallCalibration\CalibrationModel.mat','ModelParams');
    elseif exist('K:\Toni\BallCalibration\CalibrationModel.mat','file')
        load('K:\Toni\BallCalibration\CalibrationModel.mat','ModelParams');
    elseif exist('K:\fs3-joeschgrp\Toni\BallCalibration\CalibrationModel.mat','file')
        load('K:\fs3-joeschgrp\Toni\BallCalibration\CalibrationModel.mat','ModelParams');
    elseif exist('B:\fs3-joeschgrp\Toni\BallCalibration\CalibrationModel.mat','file')
        load('B:\fs3-joeschgrp\Toni\BallCalibration\CalibrationModel.mat','ModelParams');
    else
        [file,path]=uigetfile('*CalibrationModel.mat','select calibration file');
        load([path '\' file],'ModelParams');
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










