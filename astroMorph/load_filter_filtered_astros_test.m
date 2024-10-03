%% Load and filter binarized astrocytes
%Last modified by MES on 7/24/24

parentdir = 'A:/Margaret/Astrocytes/Exp173_Aldh1l1-Cre_CAG-FLEX-GFP_4x-exp_astrocyte_morphology/mouse8/bin/';

files = dir([parentdir '*tif*']);

for fidx = 1:length(files)
    imgbin = mat2gray(loadtiff([parentdir files(fidx).name]));
    imgbin = ~imgbin;
    imgfilt = bwareaopen(imgbin,1e6);
    cc = bwconncomp(imgfilt,26);
    disp(cc.NumObjects)
    volume=regionprops3(cc,'Volume');
    vol=volume.Volume;
    options.overwrite=true;
    saveastiff(double(imgfilt),[parentdir files(fidx).name(1:end-4) '_filt.tif'],options);
end