function RB_import_VideoBehav
params.movmed_filt=3;
params.pupil_min_p=.95;
params.body_min_p=.8;
params.eye_r_mm=1.5;
params.eye_px_per_mm=23.4;
params.body_vid=480;
params.overwrite=true;
params.delete_camera_tifs=false;
%parentdir='K:\Rimma\Camera\';
parentdir = 'Z:\Group Members\Rima\Camera\'; 
archivedir = 'Z:\Group Members\Rima\preCamera';
DLCdir='Z:\Group Members\Rima\DLCData\';
preprocesseddir='Z:\Group Members\Rima\ToAnalyse\';
%archivedir='K:\Florian\preCamera\';

data_exists=true;
new_data=true;

%%
while data_exists
    
    
    animallist=dir(parentdir);animallist(1:3)=[];

    if isempty(animallist)
        data_exists=false;
        fprintf('no data available anymore\n')
        break;
    end
    for a=1:numel(animallist)
        explist=dir([animallist(a).folder '\' animallist(a).name]);explist(1:2)=[];

        if isempty(explist)
            rmdir([animallist(a).folder '\' animallist(a).name]);
        else
            for e=1:numel(explist)
                clearvars -except params new_data data_exists parentdir preprocesseddir DLCdir animallist a explist e archivedir
                if exist([explist(e).folder '\' explist(e).name '\success.log'],'file')
                    %exist([explist(e).folder '\' explist(e).name '\eye_results'],'dir')
                    new_data=true;
                    finalFile=[preprocesseddir explist(e).name '_preprocessed.mat'];
                    preprocessed_exists=exist(finalFile,'file');
                    if preprocessed_exists
                        DLCnotimported=isempty(who('-file',finalFile,'VideoParams'));
                    else
                        DLCnotimported=true;
                    end
                    if ~preprocessed_exists || ((preprocessed_exists && DLCnotimported) || params.overwrite)
                        if ~exist([DLCdir explist(e).name '_DLC_filtered.mat'],'file') || params.overwrite
                            eyename=[explist(e).folder '\' explist(e).name '\eye_results\' explist(e).name '_eye.csv'];
                            if ~exist(eyename,'file')
                                tempf=dir([explist(e).folder '\' explist(e).name '\eye_results\*_eye.csv']);
                                 eyename=fullfile(tempf.folder,tempf.name);
                            end
                            bodyname=[explist(e).folder '\' explist(e).name '\body_results\' explist(e).name '_body.csv'];
                             if ~exist(bodyname,'file')
                                tempf=dir([explist(e).folder '\' explist(e).name '\body_results\*_body.csv']);
                                 bodyname=fullfile(tempf.folder,tempf.name);
                            end
                            fprintf('importing %s DLC file\n',explist(e).name)
                            [DLC,VideoFrameTimes] = import_DLC_eye(eyename,params);
                            [DLCbody,DLCraw] = import_DLC_body(bodyname,params);
                            DLCraw.PupilPoints=DLC.Pupil.PointsRaw;
                            VideoParams.Eye.Radius_mm=DLC.eye_r_mm;
                            VideoParams.Eye.Center=DLC.Eye.center;
                            VideoParams.Eye.EllipseParams=DLC.Eye.A;
                            VideoParams.Eye.Semiaxes=DLC.Eye.semiaxes;
                            VideoParams.Eye.Rotation=DLC.Eye.phi;
                            VideoParams.Eye.medfilt=DLC.Eye.denoise.movmed_points;
                            VideoParams.Eye.CornealReflection=DLC.cornealReflection;
                            VideoParams.Eye.px_per_mm=DLC.eye_px_per_mm;
                            VideoParams.Pupil.p_thresh=DLC.Pupil.denoise.p_thresh;
                            VideoParams.Body.HeadplatRBC=DLCbody.headplateRBC;
                            VideoParams.Body.BallLoc_RAC=DLCbody.BallRAC;
                            VideoParams.Body.EyeLocFull=DLCbody.EyeLocFull;
                            VideoParams.Body.p_trhesh=DLCbody.prob_cutoff;
                            VideoParams.Body.medfilt=DLCbody.medfilt;
                            
                            DLCSupporting.Probability.EyeLid=DLC.Eye.p;
                            DLCSupporting.Probability.Pupil=DLC.Pupil.p;
                            DLCSupporting.Probability.EllipseRMSE=DLC.Pupil.RMSE;
                            DLCSupporting.Probability.BallApex=DLCbody.ballApex(:,3);
                            DLCSupporting.Probability.BallRostral=DLCbody.ballRostral(:,3);
                            DLCSupporting.Probability.BallCaudal=DLCbody.ballCaudal(:,3);
                            
                            DLCSupporting.Pupil.Eccentricity=DLC.Pupil.eccentricity;
                            DLCSupporting.Pupil.EllipsePhi=DLC.Pupil.Rotation;
                            DLCSupporting.Pupil.EllipseSemiAxes=DLC.Pupil.semiaxes;
                            DLCSupporting.Pupil.EllipseCenter=DLC.Pupil.center;
                            DLCSupporting.Pupil.Points=DLC.Pupil.PointsRaw;
                            DLCSupporting.CropParams=DLC.CropParams;
                            DLCSupporting.CropParams.Body_vid_width=params.body_vid;
                            DLCSupporting.CropParams.Body_scale=DLCbody.bodyvid_scale;
                            DLCSupporting.Ball.Apex=DLCbody.ballApex(:,1:2);
                            DLCSupporting.Ball.Rostral=DLCbody.ballRostral(:,1:2);
                            DLCSupporting.Ball.Caudal=DLCbody.ballCaudal(:,1:2);
                            
                            Behav.VideoNum=DLC.Vid(:,1);
                            Behav.FrameNum=DLC.Vid(:,2);
                            Behav.PupilAzEl=DLC.Pupil.az_el;
                            Behav.PupilArea=DLC.Pupil.area;
                            FN=fieldnames(DLCbody);
                            FN=FN(1:19);
                            for fn=1:numel(FN)
                                Behav.(FN{fn})=DLCbody.(FN{fn})(:,1:2);
                                DLCSupporting.Probability.(FN{fn})=DLCbody.(FN{fn})(:,3);
                            end
                            Behav.EyeLidArea=DLC.Eye.area;
                            Behav.EyeLidCenter=DLC.Eye.position;
                            
                            fprintf('saving %s DLC file\n',explist(e).name)
                            
                            save([DLCdir explist(e).name '_DLC_raw.mat'],'DLCraw')
                            save([DLCdir explist(e).name '_DLC_filtered.mat'],'DLCSupporting','Behav','VideoParams','VideoFrameTimes')
                            clear DLC DLCraw DLCbody FN
                        else
                            fprintf('imported dlc exists... loading')
                            load([DLCdir explist(e).name '_DLC_filtered.mat'],'DLCSupporting','Behav','VideoParams','VideoFrameTimes')
                        end
                        if preprocessed_exists
                            fprintf('correct timing %s\n',explist(e).name)
                            
                            Data=matfile(finalFile,'Writable',true);
                            Triggers=Data.Triggers;
                            Stimuli=Data.Stimuli;
                            if isfield(VideoParams,'selected_frameIdx')
                                %load(['K:\archive-joeschgrp\Toni\Run eye correction\' explist(e).name '\frame_to_include.mat'])
                                selected_frameIdx=VideoParams.selected_frameIdx;
                            else
                                [selected_frameIdx,Triggers]=correct_videoFrames_OS(Triggers,VideoFrameTimes,Behav);
                                Data.Triggers=Triggers;
                                VideoParams.selected_frameIdx=selected_frameIdx;
                                save([DLCdir explist(e).name '_DLC_filtered.mat'],'VideoParams','-append')
                            end
                            FN=fieldnames(DLCSupporting.Probability);
                            for fn=1:numel(FN)
                                DLCSupporting.Probability.(FN{fn})=DLCSupporting.Probability.(FN{fn})(selected_frameIdx,:,:);
                            end
                            FN=fieldnames(DLCSupporting.Pupil);
                            for fn=1:numel(FN)
                                DLCSupporting.Pupil.(FN{fn})=DLCSupporting.Pupil.(FN{fn})(selected_frameIdx,:,:);
                            end
                            FN=fieldnames(DLCSupporting.Ball);
                            for fn=1:numel(FN)
                                DLCSupporting.Ball.(FN{fn})=DLCSupporting.Ball.(FN{fn})(selected_frameIdx,:,:);
                            end
                            FN=fieldnames(Behav);
                            for fn=1:numel(FN)
                                Behav.(FN{fn})=Behav.(FN{fn})(selected_frameIdx,:,:);
                            end
                            
                            Data.DLCSupporting=DLCSupporting;
                            Data.Behav=Behav;
                            Data.VideoParams=VideoParams;
                            
                        else
                            warning([explist(e).name ' no preprocessed found... skipping'])
                        end
                    else
                        warning([explist(e).name ' already imported... skipping'])
                    end
                    if params.delete_camera_tifs && exist([explist(e).folder '\' explist(e).name '\camera'],'dir')
                        rmdir([explist(e).folder '\' explist(e).name '\camera'],'s')
                    end
                    
                    mkdir([archivedir animallist(a).name '\' explist(e).name]);
                    stat=movefile([explist(e).folder '\' explist(e).name '\*'],[archivedir animallist(a).name '\' explist(e).name]);
                    if stat
                        rmdir([explist(e).folder '\' explist(e).name])
                    end
                else
                   % warning([explist(e).name ' not yet dlcd... skipping'])
                end
            end
        end
    end
    if new_data
        fprintf('waiting for dlc\n')
        new_data=false;
    end
end

%%
function [DLC,VideoFrameTimes] = import_DLC_eye(filename,params)

delimiter = ',';

%% get general variables
fileID = fopen(filename,'r');
 startRow = 2;
    endRow = 2;
    formatSpec = '%C%f%f%f%f%f%f%f%f%f%f%f%f%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

    dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end
    output_general = table2struct(table(dataArray{1:end-1}, 'VariableNames', {'Experiment_Name','eye_center_x','eye_center_y','eye_width','eye_height','eye_theta','pupil_jackknife_hold','eye_R','eye_area','crop_size','original_width','original_height','original_eye_center_x','original_eye_center_y'}));
fclose(fileID);
%% get frame variables
fileID = fopen(filename,'r');
 startRow = 4;
    endRow = inf;

formatSpec = '%s%f%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%C%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end


fclose(fileID);

output_frame = table2struct(table(dataArray{1:end-1}, 'VariableNames', {'videoref','frame','tiff_name','time_stamp','eye3D_x','eye3D_y','eye3D_z','eye_width_2','eye_heigth','eye_theta','eye_area','pupil3d_x','pupil3d_y','pupil3d_z','pupil_width','pupil_height','pupil_theta','pupil_el','pupil_az','pupil_std_x','pupil_std_y','jackknife_inds','pupilN','pupilN1','pupilN2','pupilNE','pupilNE1','pupilNE2','pupilE','pupilE1','pupilE2','pupilSE','pupilSE1','pupilSE2','pupilS','pupilS1','pupilS2','pupilSW','pupilSW1','pupilSW2','pupilW','pupilW1','pupilW2','pupilNW','pupilNW1','pupilNW2','cornealReflection','cornealReflection1','cornealReflection2','eyeN','eyeN1','eyeN2','eyeNE','eyeNE1','eyeNE2','eyeE','eyeE1','eyeE2','eyeSE','eyeSE1','eyeSE2','eyeS','eyeS1','eyeS2','eyeSW','eyeSW1','eyeSW2','eyeW','eyeW1','eyeW2','eyeNW','eyeNW1','eyeNW2'}));

%%
temp=struct2cell(output_frame);
frame_temp=temp(1,:)';
DLC.Vid=cat(2,str2double(cellfun(@(x) x{1}(end-4:end),frame_temp,'UniformOutput',0)),...
    cat(1,output_frame(:).frame));

DLC.Eye.denoise.movmed_points=params.movmed_filt;
DLC.Eye.center=[output_general.eye_center_x output_general.eye_center_y];
DLC.Eye.semiaxes=[output_general.eye_width output_general.eye_height]/2;
DLC.Eye.phi=output_general.eye_theta;
DLC.Eye.area=cat(1,output_frame(:).eye_area);
DLC.Eye.position=cat(2,cat(1,output_frame(:).eye3D_x),cat(1,output_frame(:).eye3D_y));
DLC.Eye.p=mean(cell2mat(temp(end-21:3:end,:)'),2);
DLC.cornealReflection=median(cell2mat(temp(end-26:end-25,:)'),1);

DLC.Pupil.center=cat(2,cat(1,output_frame(:).pupil3d_x),cat(1,output_frame(:).pupil3d_y));
DLC.Pupil.semiaxes=cat(2,cat(1,output_frame(:).pupil_width),cat(1,output_frame(:).pupil_height))/2;
DLC.Pupil.Rotation=cat(1,output_frame(:).pupil_theta);
DLC.Pupil.RMSE=mean(cat(2,cat(1,output_frame(:).pupil_std_x),cat(1,output_frame(:).pupil_std_y)),2);
DLC.Pupil.p=cell2mat(temp(25:3:46,:)');
DLC.Pupil.denoise.p_thresh=params.pupil_min_p;
DLC.Pupil.denoise.movmed_points=params.movmed_filt;
DLC.Pupil.PointsRaw=reshape(cell2mat(temp(23:46,:)'),size(temp,2),3,8);
p=repmat(reshape(DLC.Pupil.p,[],1,8)<DLC.Pupil.denoise.p_thresh,1,2,1);
% temp_pupil=DLC.Pupil.PointsRaw(:,1:2,:);
% temp_pupil(p)=nan;
% temp_pupil=fillmissing(temp_pupil,'linear');
%DLC.Pupil.Points=movmedian(temp_pupil,DLC.Pupil.denoise.movmed_points);
DLC.Pupil.eccentricity=sqrt(1-(DLC.Pupil.semiaxes(:,1).^2)./(DLC.Pupil.semiaxes(:,2).^2));
DLC.Pupil.area=prod(DLC.Pupil.semiaxes,2)*pi;

DLC.eye_r_mm=params.eye_r_mm;
DLC.eye_px_per_mm=params.eye_px_per_mm;

DLC.Pupil.az_el=cat(2,cat(1,output_frame(:).pupil_az),cat(1,output_frame(:).pupil_el));

DLC.CropParams.Size=output_general.crop_size;
DLC.CropParams.original_width=output_general.original_width;
DLC.CropParams.original_height=output_general.original_height;
DLC.CropParams.original_eye_center_x=output_general.original_eye_center_x;
DLC.CropParams.original_eye_center_y=output_general.original_eye_center_y;
DLC.Eye.A=nan(1,8);
videoframetimestemp=cellfun(@(x) x{1},temp(4,:),'UniformOutput',0)';
if str2double(videoframetimestemp{1}(1:4))==0
    VideoFrameTimes=nan;
else
    VideoFrameTimes=datetime(videoframetimestemp,'InputFormat','yyyy_MM_dd HH:mm:ss:SSS');
end
end

%%
function [DLC,DLCraw] = import_DLC_body(filename,params)
%%
delimiter = ',';
startRow = 4;
endRow = inf;

formatSpec = '%C%f%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

fileID = fopen(filename,'r');


dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);


dataArray=dataArray(1:end-1);
NameOut = import_body_part_name(filename)';
[NameUn,~,ib]=unique(NameOut,'stable');
%Output.frameN=dataArray{1};


for n=1:numel(NameUn)
    
DLC.(NameUn{n})=cat(2,dataArray{ib==n});
end 

DLC.headplateRBC=[median(DLC.headplateRostralCorner(:,1:2),1);median(DLC.headplateBendCorner(:,1:2),1);nan(1,2)];
DLC.BallRAC=[median(DLC.ballRostral(:,1:2),1);median(DLC.ballApex(:,1:2),1);median(DLC.ballCaudal(:,1:2),1)];
%DLC.EyeLocFull=mean([median(DLC.eyeRostralCorner(:,1:2));median(DLC.eyeBottomCorner(:,1:2));median(DLC.eyeTopCorner(:,1:2));median(DLC.eyeCaudalCorner(:,1:2))],1);
DLC.EyeLocFull=median(DLC.eye(:,1:2));
FN=fieldnames(DLC);
DLC=rmfield(DLC,FN([1:4 8:9]));
FN=fieldnames(DLC);
DLCraw=DLC;
DLC.prob_cutoff=params.body_min_p;
DLC.medfilt=params.movmed_filt;
DLC.bodyvid_scale = import_DLC_eyescale(filename);
fprintf('filtering\n')
for fn=1:19
            DLC.(FN{fn})(DLC.(FN{fn})(:,end)<DLC.prob_cutoff,1:end-1)=nan;
            DLC.(FN{fn})(:,1:end-1)=fillmissing(DLC.(FN{fn})(:,1:end-1),'linear');
            DLC.(FN{fn})(:,1:end-1)=movmedian(DLC.(FN{fn})(:,1:end-1),DLC.medfilt);
            DLC.(FN{fn})(DLC.(FN{fn})(:,end)<DLC.prob_cutoff,1:end-1)=nan;
end

end

%%
function NameOut = import_body_part_name(filename)

delimiter = ',';
startRow = 3;
endRow = 3;
formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';
fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

fclose(fileID);
nonemptystrings=cellfun(@(x) strlength(x),dataArray)>0;
NameOut = [dataArray{nonemptystrings}];

end

%%
function scale = import_DLC_eyescale(filename)
    startRow = 2;
    endRow = 2;

delimiter = ',';
formatSpec = '%*s%*s%*s%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';


fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    dataArray{1} = [dataArray{1};dataArrayBlock{1}];
end

fclose(fileID);

scale = [dataArray{1:end-1}];
end



end




