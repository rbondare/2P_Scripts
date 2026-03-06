function [Saccade,direction,speed,selBehavTimes, params]=get_saccades(StimulusFrameTimes,VideoFrameTimes,PupilAz,PupilEl,params, stim_onoff_sr)
%% function to get saccade times and parameters
%% inputs:
% StimulusFrameTimes: time indices (in seconds) for stimulus frames to include
    %(includeTimes fs needs to be >= VideoFrameTimes, useful if only
    %considering parts of recording. empty if all recording
% VideoFrameTimes: time indices in seconds for behavior video frames
% PupilAz: position of eye in azimuth (degrees), useful to have median subtracted
% PupilEl: position of eye in elevation (degrees), useful to have median subtracted
% params.x.x=default;
% params.Behavior.min_saccade_interval=.25          : min interval between saccades
% params.Behavior.min_saccade_speed=45              : min speed of saccades (filtered deg/sec)
% params.Behavior.min_saccade_dur=2/50              : min duration of saccade (seconds or Nframes/fs)
% params.Behavior.max_saccade_dur=.25               : max duration of saccade (seconds)
% params.Behavior.min_saccade_amp=2                 : min amplitude of saccade (degrees)
% params.Behavior.saccade_median_filter_length=.7   : filter length in seconds
%% outputs:
%Saccade: struct with Saccade parameters for each Saccade
%optional:
    %direction: direction of eye movements per sample in rad
    %speed: speed of eye movements per sample in rad/sec
    %selBehavTimes: selected timepoints of VideoFrameTimes
debug_plot=true;

if ~exist("params","var") || isempty(params)
%     params.x.x=default;
    params.Behavior.min_saccade_interval=.25;          %: min interval between saccades
    params.Behavior.min_saccade_speed=45;              %: min speed of saccades (filtered deg/sec)
    params.Behavior.min_saccade_dur=2/30;%50;              %: min duration of saccade (seconds or Nframes/fs)
    params.Behavior.max_saccade_dur=.25;               %: max duration of saccade (seconds)
    params.Behavior.min_saccade_amp=2;                 %: min amplitude of saccade (degrees)
    params.Behavior.saccade_median_filter_length=.7;   %: filter length in seconds
end

%% subselect timepoints
if ~exist("stim_onoff_sr","var") || isempty(stim_onoff_sr)
    if isempty(StimulusFrameTimes)
        selBehavTimes=VideoFrameTimes(:); %if empty use all VideoTimePoints
        TestIdx=1:length(selBehavTimes);
    elseif length(StimulusFrameTimes)==2 % if start & end time of stimulus is given
        TestIdx=VideoFrameTimes>=StimulusFrameTimes(1) & VideoFrameTimes<=StimulusFrameTimes(2);
        selBehavTimes=VideoFrameTimes(TestIdx);
        selBehavTimes=selBehavTimes(:);
    else
        TestIdx=unique(interp1(VideoFrameTimes,1:numel(VideoFrameTimes),StimulusFrameTimes,'nearest')); %find indices of VideoFrameTimes close to StimulusFrameTimes
        TestIdx(isnan(TestIdx))=[]; % purge indices out of range
        selBehavTimes=VideoFrameTimes(TestIdx); %subselect VideoFrameTimes by indices
        selBehavTimes=selBehavTimes(:);
    end
else %select videoframes that happend during stim on
    VideoFrameTimes_s = uint32(VideoFrameTimes*stim_onoff_sr);
    TestIdx = find(StimulusFrameTimes(VideoFrameTimes_s)); %VideoFrameTimes during the stimulus
    selBehavTimes=VideoFrameTimes(TestIdx); %subselect VideoFrameTimes by indices
    selBehavTimes=selBehavTimes(:);
end
%% prepare variables
relativeSamples=50; %samples before/after saccade to check for changes
min_interval_samples=floor(params.Behavior.min_saccade_interval/median(diff(selBehavTimes)));
if ~isempty(PupilEl)
Pupil=cat(2,PupilAz(TestIdx),PupilEl(TestIdx)); %rearrange eye position in samples x 2 
else
Pupil=cat(2,PupilAz(TestIdx),zeros(size(TestIdx),'single')); %rearrange eye position in samples x 2 
end
if mean(any(isnan(Pupil),2))>.1
    error('pupil trace contains >10% missing values')
else
    Pupil=fillmissing(Pupil,'linear',1); % fill missing (NaN) values
end
PupilF=movmedian(Pupil,params.Behavior.saccade_median_filter_length,1,'SamplePoints',selBehavTimes);% filter pupil trace
%PupilF=cat(2,smooth(PupilF1(:,1),15,'loess'),smooth(PupilF1(:,2),15,'loess')); %alternate filtering option
Vel=cat(1,[0 0],diff(PupilF,1,1)./diff(selBehavTimes)); %compute velocity from position
[direction,speed]=cart2pol(Vel(:,1),Vel(:,2)); %transform velocity to polar coordinates

%% find peaks in eye movements
[saccade_peak_h,saccade_peak_loc]=findpeaks(speed,...
    'MinPeakHeight',params.Behavior.min_saccade_speed,... %select only peaks of sufficient magnitude
    'MinPeakDistance',min_interval_samples,... %select only peaks of sufficient interval
    'MinPeakProminence',params.Behavior.min_saccade_speed/2); %select only peaks of sufficient prominence (to avoid 2-step peaks)

% if debug_plot
%     fig=figure;
%     plot(Pupil(:,1),'Color',[.7 .7 .7]);hold on;plot(PupilF(:,1),'k');
%     plot(Pupil(:,2),'Color',[.6 .6 .95]);hold on;plot(PupilF(:,2),'b');
%     plot(saccade_peak_loc,PupilF(saccade_peak_loc,1),'r*');
%     yyaxis right
%     plot(speed,'r');hold on
%     yyaxis left
% %     pause
% %     close(fig) 
% end

% get eye speed trace relative to preliminary saccades
sacc_around_idx=saccade_peak_loc+(-relativeSamples:relativeSamples); % relative indices
oori=sacc_around_idx<1 | sacc_around_idx>numel(speed);
sacc_around_idx(sacc_around_idx<1)=1; % purge indices out of range
sacc_around_idx(sacc_around_idx>numel(speed))=numel(speed); % purge indices out of range
speed=reshape(speed,1,[]);
sacc_around_speed=speed(sacc_around_idx);% get relative speed
sacc_around_speed(oori)=0;
speed=reshape(speed,[],1);

%% compute saccade parameters
%preallocate parameters
saccade_start_end=zeros(numel(saccade_peak_loc),2);
saccade_start_pos=zeros(numel(saccade_peak_loc),2);
saccade_end_pos=zeros(numel(saccade_peak_loc),2);
saccade_amp=zeros(numel(saccade_peak_loc),1);
saccade_amp_azel=zeros(numel(saccade_peak_loc),2);
saccade_dir=zeros(numel(saccade_peak_loc),1);
for s=1:numel(saccade_peak_loc)
    saccade_start_end(s,1)=find(sacc_around_speed(s,1:relativeSamples+1)<saccade_peak_h(s)/2,1,'last')-relativeSamples-1;%saccade onset relative (crossing half magnitude)
    saccade_start_end(s,2)=find(sacc_around_speed(s,relativeSamples+1:end)<saccade_peak_h(s)/2,1,'first');%saccade offset relative (crossing half magnitude)
    saccade_start_end(s,:)=saccade_start_end(s,:)+saccade_peak_loc(s,:); %saccade on/offset in samples
    saccade_start_end(saccade_start_end-5>size(Pupil,1))=size(Pupil,1)-5; %crop offset if outside of samples
    saccade_end_pos(s,:)=median(Pupil(saccade_start_end(s,2):min(saccade_start_end(s,2)+5,size(Pupil,1)),:),1);%find stable end position
    saccade_start_pos(s,:)=median(PupilF(max(1,saccade_start_end(s,1)-5):saccade_start_end(s,1),:),1);%find stable start position
    temp=saccade_end_pos(s,:)-saccade_start_pos(s,:);
    [saccade_dir(s),saccade_amp(s)]=cart2pol(temp(1),temp(2)); % polar coordinate of saccade
    saccade_amp_azel(s,:)=temp(:)'; %saccade magnitude in each dim
end

%% purge/correct preliminary saccades
saccade_start_end(saccade_start_end>numel(selBehavTimes))=numel(selBehavTimes);%crop if out of range
saccade_dur=diff(selBehavTimes(saccade_start_end),1,2);%duration of saccade
sel=saccade_amp>=params.Behavior.min_saccade_amp & ... %select saccades of minimal amplitude
    saccade_dur>=params.Behavior.min_saccade_dur & ... %select saccades of minimal duration
    saccade_dur<=params.Behavior.max_saccade_dur;       %select saccades of maximal duration

%% assign saccades to output structure
Saccade.PeakTime=selBehavTimes(saccade_peak_loc(sel));
Saccade.PeakSpeed=speed(saccade_peak_loc(sel));
Saccade.StartTime=selBehavTimes(saccade_start_end(sel,1));
Saccade.EndTime=selBehavTimes(saccade_start_end(sel,2));
Saccade.Duration=saccade_dur(sel);
Saccade.Direction=saccade_dir(sel);
Saccade.Amplitude=saccade_amp(sel);
Saccade.AmpAzEl=saccade_amp_azel(sel,:);
Saccade.StartPos=saccade_start_pos(sel,:);
Saccade.EndPos= saccade_end_pos(sel,:);




















