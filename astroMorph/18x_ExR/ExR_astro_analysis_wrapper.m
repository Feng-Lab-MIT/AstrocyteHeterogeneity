%%  WRAPPER for analyzing multiExR astrocyte data from October - December 2024 %%

% Last modified by Margaret Schroeder on 3/18/25

% Image processing up to this point:
% - background subtraction in Fiji

% NOTE: It's best to run this wrapper script section by section and be
% mindful of the outputs

%% Set up and parameter setting

clear all

%the directory containing all of your image files, named according to the
%following convention: [FIELD OF VIEW NAME]_[ROUND]_[CHANNEL].tif, as
%output from the ExSeqProcessing registration pipeline
params.parentfolder = 'A:/Margaret/Astrocytes/Exp194_ExR_astro_rdegs/preprocessed/';

%set parameters
params.xystep = 0.1625/18; %physical pixel size divided by expansoin factor, um/voxel in x and y
params.zstep = 0.25/18; %physical z-step size divided by expansion factor, um/voxel in z
params.dofilt = 1; %1 to do median filter to eliminate high spatial frequency noise
params.filt_size = [5 5 3]; %filter size [x y z] for 3D median filter
params.filt = 'med'; %method for filtering image
params.vol_converter = ((2/3)*params.xystep + (1/3)*params.zstep)^3; %weighted average to convert voxels to um^3
params.nchannels=3;
params.gfilt=1; %whether or not to blur the synaptic reference channel for synapse segmentation
params.sizefilt=1; %whether or not to do size filtration

params.savetargetmasks=1;
params.maskfolder='A:/Margaret/Astrocytes/Exp194_ExR_astro_rdegs/gfp_masks/'; %for astrocyte process segmentations
params.outfolder ='A:/Margaret/Astrocytes/Exp194_ExR_astro_rdegs/syn_masks/'; %for synapse segmentations
params.imgbit = 2e16; %maximum pixel value to convert to grayscale


%% Make and save the GFP channel masks
params.thresh_method = 'zscore';
params.gfp_thresh = 0.75;

params.lowerlim = 200; %for GFP channel
params.gfp_sigma=10;
params.gfp_folder = params.parentfolder;
params.gfp_channel = 3;
params.savegfp = 1;

make_gfp_mask_exr_auto(params);

%% Make and save the synapse masks

params.cavsynfolder = params.parentfolder;
params.cav_channel = 1;
params.sigma=3;
params.lowerlim=100; %for synapse identification
params.savemasks=1;
params.use_adaptive_thresh=0;

params.thresh_method = 'zscore';
params.thresh_multiplier = 3;

params.doplot=1;
params.upperlim = 5000; %objects larger than this will get filtered out
params.nbuffer = 30;%number of pixels by which to buffer each bounding box (total)

% data_syn = make_synaptic_mask_exr(params);

%% Run quantification in the target channel using previously generated masks - needs updating

params.target_channel = 2;

data_target = run_target_quant_exr(params);
save('A:/Margaret/Astrocytes/Exp194_ExR_astro_rdegs/analyzed_target_data_20250413.mat','data_target')
