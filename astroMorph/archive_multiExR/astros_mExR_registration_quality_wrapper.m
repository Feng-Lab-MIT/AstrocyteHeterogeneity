%% SCRIPT TO EVALUATE QUALITY OF REGISTRATION 

% Uses Dan Goodwin's method from ExSeq paper (Alon et al., Science 2021)
% Modified from code used in the manuscript manuscript by Kang*, Schroeder* et al., 2024
% Last modified by Margaret Schroeder on 8/18/24

% Image processing up to this point:
% - background subtraction in Fiji and
% - registration using modified ExSeq pipeline

%% Set up and set parameters
clear all

%folder where the field of view image volumes are stored (either full or
%cropped z-stacks to mutually overlapping areas)
parentfolder = 'A:/Margaret/Astrocytes/Exp174_Aldh1l1-Cre_CAG-FLEX-GFP_18x-exp_multiExR/registered/round1.1/';

fovs = {
    'mouse1-A-fovroi1';
    'mouse1-A-fovroi3';
    'mouse1-A-fovroi4';
    'mouse1-A-fovroi5';
    'mouse1-A-fovroi6';

    'mouse1-B-fovroi1';
    'mouse1-B-fovroi2';
    % 'mouse1-B-fovroi3';
    'mouse1-B-fovroi4';
    'mouse1-B-fovroi5';
    'mouse1-B-fovroi6';

    'mouse1-C-fovroi1';
    'mouse1-C-fovroi2';
    'mouse1-C-fovroi4';
    'mouse1-C-fovroi5';
    'mouse1-C-fovroi6';

    'mouse3-A-fovroi1';
    'mouse3-A-fovroi2';
    'mouse3-A-fovroi4';
    'mouse3-A-fovroi6';

    'mouse3-B-fovroi1';
    'mouse3-B-fovroi3';
    'mouse3-B-fovroi4';
    'mouse3-B-fovroi5';
    'mouse3-B-fovroi6';

    'mouse3-C-fovroi1';
    'mouse3-C-fovroi2';
    'mouse3-C-fovroi3';
    'mouse3-C-fovroi4';
    'mouse3-C-fovroi6';

    };

params.xystep = 0.1625/16; %physical pixel size divided by expansoin factor, um/voxel in x and y
params.zstep = 0.25/16; %physical z-step size divided by expansion factor, um/voxel in z
params.nchannels=3; %number of channels in each stack
params.parentfolder = parentfolder;
params.error_channel = {'3','3'};%,'4','4','2','2'}; %channel index, channel on which to calculate error
params.subtract_morph = 0; %subtract the morphology channel? no, because we are using it for reg quality analysis
params.morph_channel = nan; %morphology channel used for registration, for subtracting out
params.rounds = {'01','02'};%,'12','13','14','15','16','17','18','19'};%,'6','7','8','9','10','11','12'};%rounds to analyze for registration error

%Params for Dan Goodwin's registration quality evaluation
params.subvol_dim=100;%length of side of sub-volume, in pixels to use in Dan's registration quality measure
params.xrange = 1:2048;%in pixels, extent of volume to analyze
params.yrange = 1:2048;%in pixels, extent of volume to analyze
params.N = 1000; %number of sub-volumes for DG method
params.pct_thresh = 99;%percentage of pixel intensity values for thresholding
params.doplot = 1; %plotting flag (1 if output plots desired)
params.nonzero_thresh=.2*2048*2048*60;%number of nonzero pixels required to run registration - any
params.doplotrgb=1;

%% Running section

for fovidx = 1:length(fovs)
    fov = fovs{fovidx};
    disp(fov);
    error_DG{fovidx} = measure_round_alignment_mExR(params,fov);
end

%% Organize data for easy copy/paste
% modify as needed for convenient copy/paste into Excel or Prism for plotting
fovid = 2;
cp = [error_DG{1,fovid}{1,2}
    ];