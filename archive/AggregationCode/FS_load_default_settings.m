function [params,g]=FS_load_default_settings(params)


params.Stimulus.screen_px=[0 342;0 608]; %float: screen size in px
params.Stimulus.screen_deg=[80 -40;170 -100];%float: screen extent in deg from mouse
params.Stimulus.force_pattern=true;%bool: make full pattern even if sparse stim
params.Stimulus.RF_deg_per_px=1;%float: downsample stimulus to resolution in vis deg

params.Neurons.use_neurons=true;
params.Neurons.activity_type='dFF';%dFF deconvolved both
params.Neurons.SNR_cutoff=5;%float: min SNR for neurons (signal=params.high_prc of neuron, noise= std of neg. part of dFF)
params.Neurons.high_prc=99;%float(1:100) high percentile of activity considered signal strength

params.Calc.interpfun='linear';%interpolation function for neuronal act
params.Calc.Act_interp_func='linear';%interpolation function for neuronal act
params.Calc.relativeTime=[-2:.1:2];%vector: relative time of stimulus to neurons in sec
params.Calc.useGPU=gpuDeviceCount>0;
if params.Calc.useGPU
    g = gpuDevice(1);
    reset(g);
else
    g=[];
end
params.Calc.ActivityMean='mean';%'prctile95'
params.Calc.filtbarT=21;
params.Calc.bar_latency=15;
params.Calc.spatial_sigma=30;
params.Calc.fitfunc='poly22';
params.Calc.loess_span=.5;
params.Calc.use_compact=true;

params.Behavior.use_behav=true;
params.Behavior.include='default';
params.Behavior.median_filt_length=[0.05 0.05 0.05];%seconds
params.Behavior.pupil_relative=true;
params.Behavior.min_saccade_speed=45; %deg/sec
params.Behavior.min_saccade_amp=2; %deg
params.Behavior.min_saccade_dur=2/50; %sec
params.Behavior.max_saccade_dur=.25; %sec
params.Behavior.min_az_saccade=3;
params.Behavior.min_saccade_interval=.25;
params.Behavior.saccade_median_filter_length=.7;
params.colour = [215/255, 65/255,155/255;... % beginner
    255/255,166/255,230/255;... % beginner no reward
    0/255, 158/255,115/255;... % expert
    154/255, 224/255, 116/255;... % expert no reward
33/255,220/255,216/255;... % expert random
213/255,94/255,0;... %Hit
    0,114/255,178/255;...Miss
    204/255,121/255,167/255;...
    0,158/255,115/255;...
    213/255,94/255,0;...
    0.5 0.5 0.5;...
    1 0 0 ;...
    0 0 0; ...
    0 0 1;...
    0.2 0.2 0.2];
params.cmap=get_blue_red_cmap;
%params.vir =viridis;
params.Calc.do_plot=true;
params.Naive = [0.3373, 0.7059, 0.9137];
params.NaiveIndividual = [0.3373, 0.7059, 0.9137];

params.Beginner = [0.8431, 0.2549, 0.6078];
params.BeginnerIndividual = [0.8431, 0.2549, 0.6078];

params.Expert = [0, 0.6196, 0.4510];
params.ExpertIndividual = [0, 0.6196, 0.4510];

params.ExpertRandom = [ .2, .2, .2];
params.ExpertRandomIndividual = [ .2, .2, .2];

params.ExpertAll = [ 0, 0, 0];
params.ExpertAllIndividual = [ 0, 0, 0];

params.NoSpout = [0.8353, 0.3686, 0];
params.NoSpoutIndividual = [0.8353, 0.3686, 0];

params.Scrambled = [0.5020, 0.5020, 0.5020]; 
params.Hit = [0.8500, 0.3250, 0.0980];
params.HitIndividual = [0.8500, 0.3250, 0.0980];

params.Miss = [0, 0.4470, 0.7410];
params.MissIndividual = [0, 0.4470, 0.7410];