%%  WRAPPER for analyzing multiExR astrocyte data from October - December 2024 %%

% Last modified by Margaret Schroeder on 12/30/24

% Image processing up to this point:
% - background subtraction in Fiji and registration

% NOTE: It's best to run this wrapper script section by section and be
% mindful of the outputs

%% Set up and parameter setting - Cortex

clear all

%the directory containing all of your image files, named according to the
%following convention: [FIELD OF VIEW NAME]_[ROUND]_[CHANNEL].tif, as
%output from the ExSeqProcessing registration pipeline
params.parentfolder = 'A:/Margaret/Astrocytes/Exp184_multiExR_astros/registered/';

%set parameters
params.xystep = 0.1625/17; %physical pixel size divided by expansoin factor, um/voxel in x and y
params.zstep = 0.25/17; %physical z-step size divided by expansion factor, um/voxel in z
params.dofilt = 1; %1 to do median filter to eliminate high spatial frequency noise
params.filt_size = [5 5 3]; %filter size [x y z] for 3D median filter
params.filt = 'med'; %method for filtering image
params.vol_converter = ((2/3)*params.xystep + (1/3)*params.zstep)^3; %weighted average to convert voxels to um^3
params.nchannels=3;
params.gfilt=1; %whether or not to blur the synaptic reference channel for synapse segmentation
params.sizefilt=1; %whether or not to do size filtration

params.savetargetmasks=1;
params.maskfolder='A:/Margaret/Astrocytes/Exp184_multiExR_astros/gfp_masks/'; %for astrocyte process segmentations
params.outfolder ='A:/Margaret/Astrocytes/Exp184_multiExR_astros/syn_masks/'; %for synapse segmentations
params.imgbit = 2e16; %maximum pixel value to convert to grayscale


%% Make and save the GFP channel masks
params.thresh_method = 'stdev';
params.gfp_thresh = 5;
params.lower_percentile = 75;

params.lowerlim = 200; %for GFP channel
params.gfp_sigma=10;
params.gfp_folder = [params.parentfolder 'round1.1'];
params.gfp_channel = 'ch03';
params.savegfp = 1;

make_gfp_mask_auto(params);

%% Make and save the synapse masks

params.cavsynfolder = params.gfp_folder;
params.cav_channel = 'ch01';
params.sigma=3;
params.syn_channels = {'2-3','5-3'}; %either one of these needs to be colocalized with Cav2.1 to call it a synapse
params.lowerlim=100; %for synapse identification
params.savemasks=1;
params.use_adaptive_thresh=0;

params.thresh_method = 'stdev';
params.cav_thresh_multiplier = 5;
params.thresh_multiplier = 7;
params.lower_percentile = 75;

params.doplot=1;
params.upperlim = 5000; %objects larger than this will get filtered out
params.nbuffer = 30;%number of pixels by which to buffer each bounding box (total)

data_syn = make_synaptic_mask(params);

%% Run quantification in the target channel using previously generated masks - needs updating

params.gfp_folder = [params.parentfolder 'round1.1'];
params.sigma=1;
params.cavsynfolder = params.gfp_folder;
params.cav_channel = 'ch01';
params.target_sizefilt = 1;
params.target_names = {'Cav2.1','GLAST','SPARC','PSD95','GluA1','NR2B','mGluR3','Trmp3','Gat3'};
params.target_channels = {'2-1';%this is Cav2.1, should be round 1.1 so direct to alternate folder
    '1-2';
    '2-2';'2-3';
    '3-2';'3-3';
    '4-2';
    '5-2';
    '6-2';
    };
params.savemasks=1;
params.lowerlim=100; %for target objects
params.upperlim = 5000; %objects larger than this will get filtered out
params.use_adaptive_thresh = 0;
params.in_over_out_thresh = 2;
params.thresh_multiplier=1;

params.thresh_multiplier = 4; 
params.thresh_method = 'zscore';

params.targetoutfolder ='A:/Margaret/Astrocytes/Exp184_multiExR_astros/target_masks/'; %for target rDEG segmentations
data_target = run_target_quant(params);

save('A:/Margaret/Astrocytes/Exp184_multiExR_astros/analyzed_target_data_20250203.mat','data_target')
