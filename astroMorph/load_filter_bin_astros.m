%% Load, filter, andd binarize astrocytes, and calculate their surface area and volume
%Last modified by MES on 8/9/24

%Change filenames here!!
parentdir = 'A:/Margaret/Astrocytes/Exp177_Aldh1l1-Cre_CAG-FLEX-GFP_4x-exp_astro_morph/striatum/preprocessed/cropped/';
savedir = 'A:/Margaret/Astrocytes/Exp177_Aldh1l1-Cre_CAG-FLEX-GFP_4x-exp_astro_morph/striatum/segmentations_1std/';

files = dir([parentdir '*tif*']);
params.thresh_method = 'zscore';
params.thresh_multiplier=1; %standard deviations above the mean
params.sigma=2; %for gaussian filter BEFORE thresholding
params.lowerlim=1e7; %in voxels, for size filter

vols = zeros(length(files),1);
SA = zeros(length(files),1);

for fidx = 1:length(files)
    img = mat2gray(loadtiff([parentdir files(fidx).name]));
    %blur the image
    imgblur = imgaussfilt3(img,params.sigma);
    %binarize the image
    imgbin = binarize_intensity_threshold(imgblur,params);
    %filter the image
    imgfilt = bwareaopen(imgbin,params.lowerlim);
    %find connected components
    cc = bwconncomp(imgfilt,26);
    disp(cc.NumObjects)
    if cc.NumObjects > 1 %only use the largest object
        volume=regionprops3(cc,'Volume');
        lowerlim = max(volume.Volume);
        imgfilt = bwareaopen(imgbin,lowerlim-10); %assume the objects are more than 10 voxels apart in size
        cc = bwconncomp(imgfilt,26);
        disp(cc.NumObjects)
    end
    %get volume
    volume=regionprops3(cc,'Volume');
    vols(fidx,1)=volume.Volume;
    %get surface area
    satemp=regionprops3(cc,'SurfaceArea');
    SA(fidx,1)=satemp.SurfaceArea;
    options.overwrite=true;
    %save down the segmentation
    saveastiff(uint8(imgfilt*255),[savedir files(fidx).name(1:end-4) '_filt.tif'],options);
end

%Change filename here!
save([savedir 'striatum_vol_SA_workspace.mat'])