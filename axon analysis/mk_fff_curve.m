function [fff_brightness,single_pulse]=mk_fff_curve(options,ifi,timestart)
single_pulse=[];
if ~iscell(options.specParams.fff_type) && strcmp(options.specParams.fff_type,'sin_pulses_with_interval')
    %
    tb=1000;
    tl=linspace(0,options.stimulus_trial_t,options.stimulus_trial_t*tb);%timeline in ms (or other)
    options.specParams.fff_pulse_duration=options.specParams.fff_pulse_duration*tb;
   
    if ~isfield(options.specParams,'phase_off');options.specParams.phase_off=zeros(size(options.specParams.fff_avail_freqs));end
    options.specParams.phase_off=options.specParams.phase_off*tb;

    if options.specParams.expand_conditions
        single_pulse=cell(size(options.specParams.fff_pulse_duration,1),size(options.specParams.fff_avail_amps,1),size(options.specParams.fff_avail_freqs,1));
        if numel(options.specParams.phase_off)==1
            options.specParams.phase_off=repmat(options.specParams.phase_off,numel(options.specParams.fff_avail_freqs),1);
        end
        for d=1:size(options.specParams.fff_pulse_duration,1)
            for a=1:size(options.specParams.fff_avail_amps,1)
                for f=1:size(options.specParams.fff_avail_freqs,1)
                    if options.specParams.fff_pulse_duration(d)>0
                        tl_indiv=linspace(0,((options.specParams.fff_pulse_duration(d))*min(options.specParams.fff_avail_freqs(f,:))+options.specParams.phase_off(f))/tb,options.specParams.fff_pulse_duration(d)+round(options.specParams.phase_off(f)));
                    else
                        tl_indiv=linspace(0,1,tb./min(options.specParams.fff_avail_freqs(f,:)));
                    end
                    if sum(~isnan(options.specParams.fff_avail_amps(a,:)))>1
                        amp_indiv=  interp1(0:sum(~isnan(options.specParams.fff_avail_amps(a,:)))-1,    options.specParams.fff_avail_amps(a,~isnan(options.specParams.fff_avail_amps(a,:))),    linspace(0,sum(~isnan(options.specParams.fff_avail_amps(a,:)))-1,numel(tl_indiv)));
                    else
                        amp_indiv=repmat(options.specParams.fff_avail_amps(f,~isnan(options.specParams.fff_avail_amps(a,:))),1,numel(tl_indiv));
                    end
                    if sum(~isnan(options.specParams.fff_avail_freqs(f,:)))>1
                        freq_indiv= interp1(0:sum(~isnan(options.specParams.fff_avail_freqs(f,:)))-1,   options.specParams.fff_avail_freqs(f,~isnan(options.specParams.fff_avail_freqs(f,:))),  linspace(0,sum(~isnan(options.specParams.fff_avail_freqs(f,:)))-1,numel(tl_indiv)));
                    else
                        freq_indiv=repmat(options.specParams.fff_avail_freqs(f,~isnan(options.specParams.fff_avail_freqs(f,:))),1,numel(tl_indiv));
                    end
                    single_pulse{d,a,f}=(-amp_indiv/2).*sin(-pi*freq_indiv.*tl_indiv) + options.specParams.fff_baseline_lum;
                    
                end
            end
        end
    else
        conds=max([size(options.specParams.fff_pulse_duration,1) size(options.specParams.fff_avail_amps,1) size(options.specParams.fff_avail_freqs,1)]);
        if size(options.specParams.fff_pulse_duration,1)<conds
            options.specParams.fff_pulse_duration=repmat(options.specParams.fff_pulse_duration(1,:),conds,1);
        end
        if size(options.specParams.fff_avail_amps,1)<conds
            options.specParams.fff_avail_amps=repmat(options.specParams.fff_avail_amps(1,:),conds,1);
        end
        if size(options.specParams.fff_avail_freqs,1)<conds
            options.specParams.fff_avail_freqs=repmat(options.specParams.fff_avail_freqs(1,:),conds,1);
        end
        
        single_pulse=cell(size(options.specParams.fff_pulse_duration,1),1);
        for d=1:size(options.specParams.fff_pulse_duration,1)
            a=d;f=d;
            if options.specParams.fff_pulse_duration(d)>0
               if timestart<datetime(2021,09,23)
                tl_indiv=linspace(0,((options.specParams.fff_pulse_duration(d))*min(options.specParams.fff_avail_freqs(f,:))+options.specParams.phase_off(f))/tb,options.specParams.fff_pulse_duration(d)+round(options.specParams.phase_off(f)));
               else
                tl_indiv=linspace(0,((options.specParams.fff_pulse_duration(d))+options.specParams.phase_off(f))/tb,options.specParams.fff_pulse_duration(d)+round(options.specParams.phase_off(f)));
               end
            else
                tl_indiv=linspace(0,1,tb./min(options.specParams.fff_avail_freqs(f,:)));
            end
            if sum(~isnan(options.specParams.fff_avail_amps(a,:)))>1
                amp_indiv=  interp1(0:sum(~isnan(options.specParams.fff_avail_amps(a,:)))-1,    options.specParams.fff_avail_amps(a,~isnan(options.specParams.fff_avail_amps(a,:))),    linspace(0,sum(~isnan(options.specParams.fff_avail_amps(a,:)))-1,numel(tl_indiv)));
            else
                amp_indiv=repmat(options.specParams.fff_avail_amps(f,~isnan(options.specParams.fff_avail_amps(a,:))),1,numel(tl_indiv));
            end
            if sum(~isnan(options.specParams.fff_avail_freqs(f,:)))>1
                freq_indiv= interp1(0:sum(~isnan(options.specParams.fff_avail_freqs(f,:)))-1,   options.specParams.fff_avail_freqs(f,~isnan(options.specParams.fff_avail_freqs(f,:))),  linspace(0,sum(~isnan(options.specParams.fff_avail_freqs(f,:)))-1,numel(tl_indiv)));
            else
                freq_indiv=repmat(options.specParams.fff_avail_freqs(f,~isnan(options.specParams.fff_avail_freqs(f,:))),1,numel(tl_indiv));
            end
            single_pulse{d}=(-amp_indiv/2).*sin(-2*pi*freq_indiv.*tl_indiv) + options.specParams.fff_baseline_lum;
            %single_pulse{d}=(-amp_indiv/2).*sin(-pi*freq_indiv.*tl_indiv) + options.specParams.fff_baseline_lum;
        end
    end
    %
    fff_brightnessT=zeros(1,numel(tl))+options.specParams.fff_baseline_lum;
    current_time=1;
    if options.specParams.randomize_order;sporder=randperm(numel(single_pulse));else;sporder=1:numel(single_pulse);end
    current_pulse=1;
    while current_time<=numel(tl)
        current_time=current_time+options.specParams.pre_interval*tb;
        fff_brightnessT(current_time:current_time+numel(single_pulse{sporder(current_pulse)})-1)=single_pulse{sporder(current_pulse)};
        current_time=current_time+numel(single_pulse{sporder(current_pulse)});
        fff_brightnessT(current_time:current_time+options.specParams.post_interval*1000-1)=single_pulse{sporder(current_pulse)}(end);
        current_time=current_time+options.specParams.post_interval*tb;
        current_pulse=current_pulse+1;
        if current_pulse>numel(sporder);if options.specParams.randomize_order;sporder=randperm(numel(single_pulse));end;current_pulse=1;end
    end
%     figure, 
%     plot(tl,fff_brightnessT(1:numel(tl)))
%     close(gcf)
    
elseif ~iscell(options.specParams.fff_type) && strcmp(options.specParams.fff_type,'linearcolorlevels')
    levels=linspace(0,1,options.specParams.linearcolorlevels);
    [intsA,intsB] = meshgrid(levels,levels);
    start_and_end=[reshape(intsA,1,[]);reshape(intsB,1,[])]';
    start_and_end(diff(start_and_end,1,2)==0,:)=[];
    start_and_end=unique(start_and_end,'rows');
    %  start_and_end=[combnk(levels,2); fliplr(combnk(levels,2))]';
    
    
    incremental=nan(size(start_and_end,1),options.specParams.fff_durations*1000);
    for c=1:size(start_and_end,1)
        incremental(c,:)=linspace(start_and_end(c,1),start_and_end(c,2),options.specParams.fff_durations*1000);
    end
    incremental=incremental';
    if isfield(options.specParams,'rand_order') && options.specParams.rand_order
        k=1;
        while k==1
            temp=incremental(:,randperm(size(incremental,2)));
            fff_brightnessT=temp(:);
            if all(diff(fff_brightnessT)~=0)
                k=0;
            end
        end
    else
        fff_brightnessT=incremental(:);
    end
    fff_brightnessT=reshape(fff_brightnessT,1,[]);
    tl=linspace(0,numel(fff_brightnessT)/1000,numel(fff_brightnessT));%timeline in ms
    
else
    % options.fff_params.durations=[2 1 3 2 8 2 8 2 2]; %seconds
    % options.fff_params.type={'black','white','black','grey','sineamp','grey','sinefreq','grey','black'};
    tl=linspace(0,sum(options.specParams.fff_durations),sum(options.specParams.fff_durations)*1000);%timeline in ms
    fff_brightnessT=zeros(size(tl));%initialize curve (over time) with black
    sections=[1 1+cumsum(options.specParams.fff_durations(1:end-1))*1000;cumsum(options.specParams.fff_durations)*1000];%define sections of curve
    for s=1:numel(options.specParams.fff_durations)
        switch options.specParams.fff_type{s}
            case 'white'
                fff_brightnessT(sections(1,s):sections(2,s))=1;%options.maxintensity(options.color_ch); %set section to white value
            case 'bright'
                fff_brightnessT(sections(1,s):sections(2,s))=0.75;
            case {'grey','gray'}
                fff_brightnessT(sections(1,s):sections(2,s))=0.5;%options.greylevel(options.color_ch);%set section to grey value
            case 'dim'
                fff_brightnessT(sections(1,s):sections(2,s))=0.25;
            case 'sinefreq'
                t1 = linspace(0,options.specParams.fff_durations(s),1000*options.specParams.fff_durations(s));%time for sine
                f1 = linspace(0,(1/options.specParams.fff_freq_factor)*options.specParams.fff_durations(s),numel(t1));%freq for sine
                fff_brightnessT(sections(1,s):sections(2,s)) =-0.5.*sin(-2*pi*f1.*t1) + 0.5;% -options.greylevel(options.color_ch).*sin(-2*pi*f1.*t1) + options.greylevel(options.color_ch);%set section to frequency modulated sine
            case 'stepfreq'
                t1 = linspace(0,options.specParams.fff_durations(s),1000*options.specParams.fff_durations(s));%time for sine
                f1 = linspace(0,(1/options.specParams.fff_freq_factor)*options.specParams.fff_durations(s),numel(t1));%freq for sine
                fff_brightnessT(sections(1,s):sections(2,s)) =-0.5.*sign(sin(-2*pi*f1.*t1)) + 0.5;% -options.greylevel(options.color_ch).*sin(-2*pi*f1.*t1) + options.greylevel(options.color_ch);%set section to frequency modulated sine
            case 'sineamp'
                t2 = linspace(0,options.specParams.fff_durations(s),1000*options.specParams.fff_durations(s));%time for sine
                A2 = linspace(0,0.5,numel(t2));% linspace(0,options.greylevel(options.color_ch),numel(t2));%amp for sine
                fff_brightnessT(sections(1,s):sections(2,s)) = -A2.*sin(-2*pi*options.specParams.fff_sineamp.*t2) + 0.5;%options.greylevel(options.color_ch);       %set section to amplitude modulated sine
            case 'black'
                fff_brightnessT(sections(1,s):sections(2,s))=0;
                %             case 'green'
                %                 fff_brightnessT(sections(1,s):sections(2,s))=0.5;%options.greylevel(options.color_ch);%set section to grey value
                %
                
                
        end
    end
end
%[~,framenums] = unique(ceil(tl(2:end)/ifi)); %find frame numbers in time vector
% framenums=interp1(tl,1:numel(tl),0:ifi:tl(end),'nearest');
% %fff_brightness=uint8([fff_brightnessT(framenums) ones(1,100)*options.greylevel(options.color_ch)]); %remap time vector to frame vector
% fff_brightness=[fff_brightnessT(framenums) ones(1,100)*0.5]; %remap time vector to frame vector
fff_brightness=[interp1(tl,fff_brightnessT(1:numel(tl)),0:ifi:tl(end),'linear') ones(1,100)*0.5];
% figure,plot(0:ifi:tl(end),fff_brightness(1:end-100))

