%%Function to align video frames with trigger events and behaviour data 

%.. It corrects discrepancies between the number of video frames and triggers, 
%.. adjusts timing mismatches, and identifies a set of valid frame indices 
%..(selected_frameIdx) to be used in subsequent analyses.

%INPUTS
% Triggers: Struct containing timing data for events, including TimeCamera.
% VideoFrameTimes: Timestamps of the video frames.
% Behav: Struct containing behavioral data, including FrameNum (total number of behavioral frames).
% prec: Optional input (default is true) to specify the tolerance for timing precision.

% OUTPUTS:
% selected_frameIdx: Indices of valid video frames aligned with triggers.
% Triggers: Updated Triggers struct with corrected TimeCamera if necessary.
%%
function [selected_frameIdx,Triggers]=correct_videoFrames_OS(Triggers,VideoFrameTimes,Behav,prec)

if ~exist('prec','var')
    prec=true;
end
% 
% tS=fillmissing(Stimuli(end).TimeStimulusFrame(end-20:end),'spline');
% tS=tS(end);
frames2rem=find(diff([-inf Triggers.TimeCamera])>.5);
if numel(Behav.FrameNum)-numel(Triggers.TimeCamera)==1
    selected_frameIdx=2:numel(Behav.FrameNum);
elseif numel(Behav.FrameNum)-numel(Triggers.TimeCamera)==0
    selected_frameIdx=1:numel(Behav.FrameNum);
    %selected_frameIdx(frames2rem)=[];
elseif numel(Behav.FrameNum)-numel(Triggers.TimeCamera(Triggers.TimeCamera>-1))==0
    Triggers.TimeCamera=Triggers.TimeCamera(Triggers.TimeCamera>-1);
     selected_frameIdx=1:numel(Behav.FrameNum);
elseif numel(VideoFrameTimes)==1
     Triggers.TimeCamera=Triggers.TimeCamera(1:min(numel(Triggers.TimeCamera),numel(Behav.FrameNum)));
     selected_frameIdx=1:numel(Triggers.TimeCamera);
else
    dt=VideoFrameTimes;
    if prec
    tol=0.33;
    else
        tol=0.66;
    end
    if isempty(Triggers.TimeCamera)
        dt=seconds(dt-dt(2))+Triggers.TimeBall{1}(1);
    else
        dt=seconds(dt-dt(1));
        if abs(diff(dt(1:2))-diff(Triggers.TimeCamera(1:2)))>1
            dt=dt-dt(2)+Triggers.TimeCamera(1);
        end
    end
    dt=reshape(dt,1,[]);
    
    
    if numel(dt)>numel(Triggers.TimeCamera)+1 && ~isempty(Triggers.TimeCamera)
        dfC=diff([-inf dt]);
        frames2rem2=find(dfC>1.1 & dfC([2:end end])>.5);
        if  ~(numel(frames2rem2)==numel(Behav.FrameNum)-numel(Triggers.TimeCamera)) %&& ~use_write_times
            dc=median(Triggers.TimeCamera(1:500)-dt(2:501));
            selected_frameIdx=1:numel(dt);
            tm=1;tc=1;
            
            while 1
                temp=-dt(selected_frameIdx(1:min(numel(Triggers.TimeCamera),numel(selected_frameIdx))))...
                    +Triggers.TimeCamera(1:min(numel(Triggers.TimeCamera),numel(selected_frameIdx)))-dc;
                wr=find(temp>tol,tm);
                
                if ~isempty(wr)
                    wr=wr(end);
                    poff=find(temp(wr+1:end)<-5*tol);
                    tempsf=selected_frameIdx;
                    tempsf(wr)=[];
                    temp2=-dt(tempsf(1:min(numel(Triggers.TimeCamera),numel(tempsf))))+Triggers.TimeCamera(1:min(numel(Triggers.TimeCamera),numel(tempsf)))-dc;
                    pnoff=find(temp2(wr:end)<-5*tol);
%                                                         plot(max(1,wr-5000):min(wr+5000,numel(temp)),temp(max(1,wr-5000):min(wr+5000,numel(temp))));hold on;
%                                                         plot(max(1,wr-5000):min(wr+5000,numel(temp)),temp2(max(1,wr-5000):min(wr+5000,numel(temp))));hold off
%                                                                                 plot(temp);hold on;plot(temp2);hold off
%                                                           plot(cumsum(temp));hold on;plot(cumsum(temp2));hold off
%                                                                                 pause
                    if mean(temp2(wr+(0:tc-1)))<mean(temp(wr+(0:tc-1)))% && numel(pnoff)<=numel(poff)
                        selected_frameIdx(wr)=[];
                        
                    else
                        tm=tm+1;
                    end
                    
                end
                
                if numel(find(temp>tol))<=tm
                    break;
                end
                
            end
            fig= figure;
            plot(-dt(selected_frameIdx(1:min(numel(Triggers.TimeCamera),numel(selected_frameIdx))))+Triggers.TimeCamera(1:min(numel(Triggers.TimeCamera),numel(selected_frameIdx)))-dc)
            [~,~,b]=ginput(1);
            close(fig)
            if b==3
                check;
            end
            if numel(selected_frameIdx)-numel(Triggers.TimeCamera)==1
                selected_frameIdx(end)=[];
            end
            selected_frameIdx(selected_frameIdx>numel(Behav.FrameNum))=[];
            if numel(selected_frameIdx)<numel(Triggers.TimeCamera)
                if dt(selected_frameIdx(end))-Triggers.TimeCamera(end)<-5
                    Triggers.TimeCamera=Triggers.TimeCamera(1:numel(selected_frameIdx));
                    
                else
                    error('correction mismatch: more triggers than frames')
                end
            end
            if numel(selected_frameIdx)>numel(Triggers.TimeCamera)
                if dt(selected_frameIdx(end))-Triggers.TimeCamera(end)>5
                    selected_frameIdx=selected_frameIdx(1:numel(Triggers.TimeCamera));
                else
                    error('correction mismatch: more frames than triggers')
                end
            end
            
            Triggers.TimeCamera= Triggers.TimeCamera(1:numel(selected_frameIdx));
            %                     Behav.PupilCenter=Behav.PupilCenter(selected_frameIdx,:);
            %                     Behav.PupilRadius=Behav.PupilRadius(selected_frameIdx,:);
        else
            %if use_write_times;frames2rem2=[];end
            selected_frameIdx=1:numel(Behav.FrameNum);
            selected_frameIdx(frames2rem2)=[];
            %                     Behav.PupilCenter(frames2rem2,:)=[];
            %                     Behav.PupilRadius(frames2rem2,:)=[];
        end
    elseif numel(dt)==numel(Triggers.TimeCamera)+1
        selected_frameIdx=2:numel(Behav.FrameNum);
    elseif numel(dt)==numel(Triggers.TimeCamera)
        selected_frameIdx=1:numel(Behav.FrameNum);
    elseif ~isempty(Triggers.TimeCamera)
        warning('frame number lower that trigger number, using approx timing')
        dc=mean(Triggers.TimeCamera(1:500)-dt(2:501));
        dt=dt+dc;
        ddt=diff([dt(1) dt]);
        sw=find(ddt);
        ds=dt(sw);
        dt2=interp1(sw,ds,1:numel(dt),'linear','extrap');
        Triggers.TimeCamera=dt2;
        selected_frameIdx=1:numel(Behav.FrameNum);
        
        Triggers.ApproximateTiming=true;
    else
        selected_frameIdx=1:numel(Behav.FrameNum);
        Triggers.TimeCamera=dt;
        Triggers.ApproximateTiming=true;
    end
% else
%     if numel(frames2rem)>numel(Behav.FrameNum)-numel(Triggers.TimeCamera)
%         frames2rem=frames2rem(1:(numel(Behav.FrameNum)-numel(Triggers.TimeCamera)));
%         selected_frameIdx=1:numel(Behav.FrameNum);
%         selected_frameIdx(frames2rem)=[];
%         %                     Behav.PupilCenter(frames2rem,:)=[];
%         %                     Behav.PupilRadius(frames2rem,:)=[];
%         
%     else
%         selected_frameIdx=1:numel(Behav.FrameNum);
%         selected_frameIdx(frames2rem)=[];
%         selected_frameIdx=selected_frameIdx(1:numel(Triggers.TimeCamera));
%         %                     Behav.PupilCenter(frames2rem,:)=[];
%         %                     Behav.PupilRadius(frames2rem,:)=[];
%         %                     Behav.PupilCenter=Behav.PupilCenter(1:numel(Triggers.TimeCamera),:);
%         %                     Behav.PupilRadius=Behav.PupilRadius(1:numel(Triggers.TimeCamera),:);
%     end
%     if body_tracking;BehavVid.ApproximateTiming=true;end
%     Triggers.ApproximateTiming=true;
%     
    % end
    
    %             Behav.PupilCenter(frames2rem,:)=[];
    %             Behav.PupilRadius(frames2rem,:)=[];
end
% mkdir(['K:\archive-joeschgrp\Toni\Run eye correction\' name '\'])
% save(['K:\archive-joeschgrp\Toni\Run eye correction\' name '\frame_to_include.mat'],'selected_frameIdx')


