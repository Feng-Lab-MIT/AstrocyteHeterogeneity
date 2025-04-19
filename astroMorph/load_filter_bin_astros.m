%% Load, filter, andd binarize astrocytes, and calculate their surface area and volume
%Last modified by MES on 3/3/25

%Change filenames here!!
clear all
parentdir = 'A:/Margaret/Astrocytes/Exp189_4xexp_astromorph_202502/in DIW 2025.3.18/preprocessed_gfp/';
savedir = 'A:/Margaret/Astrocytes/Exp189_4xexp_astromorph_202502/in DIW 2025.3.18/segmentations_1std/';

files = dir([parentdir '*tif*']);
params.thresh_method = 'zscore';
params.thresh_multiplier=1; %standard deviations above the mean
params.thresh_multiplier_reduction = 0.3333; %how much to lower the thresh multiplier for dimmer images
params.sigma=2; %for gaussian filter BEFORE thresholding
params.lowerlim=1e6; %in voxels, for size filter

vols = zeros(length(files),1);
SA = zeros(length(files),1);
eqdis = zeros(length(files),1);
ars = zeros(length(files),1);

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
    if cc.NumObjects == 0 %if we have filtered out our largest object, the intensity threshold is likely too high
        %binarize the image using a different thresh multiplier
        thresh_multiplier_old = params.thresh_multiplier;
        %adjust the thresh multiplier
        params.thresh_multiplier = params.thresh_multiplier_reduction * params.thresh_multiplier;
        %binarize the image
        imgbin = binarize_intensity_threshold(imgblur,params);
        %filter the image
        imgfilt = bwareaopen(imgbin,params.lowerlim);
        cc = bwconncomp(imgfilt,26);
        disp(cc.NumObjects)
        %reset the thresh multiplier
        params.thresh_multiplier = thresh_multiplier_old;
    end
    if cc.NumObjects > 1 %if we have more than one object, we only want the largest objects
        volume=regionprops3(cc,'Volume');
        lowerlim = max(volume.Volume);
        imgfilt = bwareaopen(imgbin,lowerlim-10); %assume the objects are more than 10 voxels apart in size
        cc = bwconncomp(imgfilt,26);
        disp(cc.NumObjects)
    end
    if cc.NumObjects == 0 %if we still didn't detect any astrocytes
        vols(fidx,1) = NaN;
        SA(fidx,1) = NaN;
        eqdis(fidx,1) = NaN;
        ars(fidx,1) = NaN;
        disp([files(fidx).name ' has no detectable astrocyte - skipping ...'])
    else
        %get volume
        volume=regionprops3(cc,'Volume');
        vols(fidx,1)=volume.Volume;
        %get surface area
        satemp=regionprops3(cc,'SurfaceArea');
        SA(fidx,1)=satemp.SurfaceArea;
        options.overwrite=true;
        %get equivalent diamater
        eqdi=regionprops3(cc,'EquivDiameter');
        eqdis(fidx,1)=eqdi.EquivDiameter;
        %get aspect ratio
        paxs=regionprops3(cc,'PrincipalAxisLength');
        ars(fidx,1) = paxs.PrincipalAxisLength(1)/paxs.PrincipalAxisLength(2);
        %save down the segmentation
        saveastiff(uint8(imgfilt*255),[savedir files(fidx).name(1:end-4) '_filt.tif'],options);
    end
end

%Change filename here!
save([savedir 'Exp189_DIW_SoRa_vol_SA_feret_lwr_workspace_1std.mat'])