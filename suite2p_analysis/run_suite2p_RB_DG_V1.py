from suite2p import run_s2p, default_settings

settings = default_settings()

# Top level
settings['torch_device'] = 'cpu'
settings['tau']           = 0.7
settings['fs']            = 10.0
settings['diameter']      = [3.0, 3.0]

# Run
settings['run']['do_registration']     = 1
settings['run']['do_regmetrics']       = True
settings['run']['do_detection']        = True
settings['run']['do_deconvolution']    = False
settings['run']['multiplane_parallel'] = False

# IO
settings['io']['combined']      = True
settings['io']['save_mat']      = True
settings['io']['save_NWB']      = False
settings['io']['delete_bin']    = False

# Registration
settings['registration']['nimg_init']           = 500
settings['registration']['batch_size']          = 250
settings['registration']['maxregshift']         = 0.12
settings['registration']['nonrigid']            = True
settings['registration']['block_size']          = [128, 128]
settings['registration']['maxregshiftNR']       = 5
settings['registration']['smooth_sigma']        = 1.5
settings['registration']['smooth_sigma_time']   = 0.5
settings['registration']['th_badframes']        = 0.8
settings['registration']['snr_thresh']          = 3.5
settings['registration']['norm_frames']         = True
settings['registration']['subpixel']            = 10
settings['registration']['align_by_chan2']      = False
settings['registration']['reg_tif']             = False
settings['registration']['reg_tif_chan2']       = False

# Detection
settings['detection']['algorithm']         = 'sparsery'
settings['detection']['denoise']           = True
settings['detection']['threshold_scaling'] = 0.2
settings['detection']['max_overlap']       = 0.75
settings['detection']['highpass_time']     = 50
settings['detection']['soma_crop']         = False
settings['detection']['sparsery_settings']['highpass_neuropil'] = 25
settings['detection']['sparsery_settings']['spatial_scale']     = 1
settings['detection']['sparsery_settings']['max_ROIs']          = 12000
settings['detection']['sparsery_settings']['active_percentile'] = 0.0

# Classification
settings['classification']['classifier_path']        = r"Y:\Group Members\Rima\2P_Scripts\suite2p\classifier_boutons_tony.npy"
settings['classification']['use_builtin_classifier'] = False
settings['classification']['preclassify']            = 0.0

# Extraction
settings['extraction']['neuropil_extract']      = True
settings['extraction']['neuropil_coefficient']  = 0.7
settings['extraction']['inner_neuropil_radius'] = 1
settings['extraction']['min_neuropil_pixels']   = 50
settings['extraction']['lam_percentile']        = 10.0
settings['extraction']['allow_overlap']         = True

# Spike deconvolution preprocessing
settings['dcnv_preprocess']['baseline']         = 'maximin'
settings['dcnv_preprocess']['win_baseline']     = 60.0
settings['dcnv_preprocess']['sig_baseline']     = 10.0
settings['dcnv_preprocess']['prctile_baseline'] = 8.0

db = {
#    'data_path': [r"Y:\Rotation Students\Dow\DATA_2P\AnimalDG1_260302_1422"],
#    'save_path0': r"Y:\Rotation Students\Dow\DATA_2P\AnimalDG1_260302_1422",
     'data_path': [r"Y:\Rotation Students\Dow\DATA_2P\AnimalDG3_260318_1235"],
     'save_path0': r"Y:\Rotation Students\Dow\DATA_2P\AnimalDG3_260318_1235",
     'nplanes':        4,
     'nchannels':      1,
     'ignore_flyback': [-1],
}

run_s2p(db=db, settings=settings)