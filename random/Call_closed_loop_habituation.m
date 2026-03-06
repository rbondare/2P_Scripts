%%%%%% Main Script to run CL habituation
%%%%%% A.Sumser 2018/12/18
%% startup
load([fileparts(which('callstim')) '\default_options_mouse.mat']);
global teensy ballcam includeserial abort_now distortion debug_mode
includeserial=true; %false: do not use serial port (for debugging)
distortion=true;
SkipSync=0;
screenid=[];
startupLCr(274,200,758)
%% set up mouse
mouse_number=input('vlgn mouse number: ');
Animal_Name = strcat('Gad2-axoGC-VLGN-',num2str(mouse_number));
 Recording_Name=[Animal_Name '_' datestr(now,'yymmdd_HHMM')]
 targetlogfolder=['D:\Data\Habituation\' Animal_Name '\' Recording_Name '\'];mkdir(targetlogfolder);
%
 Stimnumber_times=[  
     332 10; %12.5deg/cy corridor closed loop w sometimes patch no rot gain 1
%       332 7; %12.5deg/cy corridor closed loop w sometimes patch no rot gain 1

    333 3; %12.5deg/cy corridor open loop w sometimes patch 30deg/sec
    ] ;
predefine_rng={[];[]};%same rng in cell
force_change=false;
[order,options]=mk_order(Stimnumber_times,options,predefine_rng,force_change);

%% RUN
sca
%

try % if not catching errors, we get problems with the serial ports being undefined
    copyfile(which('stimulus'),targetlogfolder)
    quitkey=initialize_ptb(options(order(1)),SkipSync,screenid); % initialize PTB, open serial ports etc
    for S=1:numel(order) % run through stimuli
        fprintf('stimulus %u of %u \n',S,numel(order))
        stimulus(options(order(S)),sprintf('%sstim_%03u_',targetlogfolder,S),1);%run stimulus
        if abort_now
            if includeserial
                BallCamData=get_BallCamSerialInfo; %get ball movements from serial input buffer
                if~strcmp(teensy.Status,'open');fopen(teensy);end
                fprintf(teensy,'3');%send stop signal
                fclose(teensy);fclose(ballcam); 
            end
            break ;
        end
    end
    
    if includeserial
        if ~strcmp(teensy.Status,'open');fopen(teensy);end
        fprintf(teensy,'3');%send stop signal
        fclose(teensy);fclose(ballcam);
    end
     
catch ME %error info in ME
    
    if includeserial
        BallCamData=get_BallCamSerialInfo; %get ball movements from serial input buffer
        if ~strcmp(teensy.Status,'open');fopen(teensy);end
        fprintf(teensy,'3');%send stop signal
        fclose(teensy);fclose(ballcam); 
    end
    sca;%close all windows and return
    ME
    ME.stack(1)
end

sca;%close all windows and return




    