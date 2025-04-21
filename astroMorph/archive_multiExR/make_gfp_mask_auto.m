function make_gfp_mask_auto(params)

%Extract parameters
parentfolder = params.gfp_folder;
fnames = dir([parentfolder filesep '*' params.gfp_channel '*.tif']);
nfovs = length(fnames);

counter = 0;

%loop through each field of view in the folder
for fidx = 1:nfovs
    img = loadtiff([parentfolder filesep fnames(fidx).name]);
    img = double(img);
   
    %% Segment GFP
    if params.gfilt
        ch3img_blur = imgaussfilt3(img,params.gfp_sigma);
    else
        ch3img_blur=img;
    end

    params.thresh_multiplier=params.gfp_thresh;
    imbingfp = binarize_intensity_threshold(ch3img_blur,params);

    if strcmp(params.filt,'med')
        imbin_filt_gfp = medfilt3(imbingfp,params.filt_size);
    else
        imbin_filt_gfp = imbingfp;
    end

    if params.sizefilt
        mask_gfp = bwareaopen(imbin_filt_gfp, params.lowerlim); %filtration, removes binary objects LESS than this size
    else
        mask_gfp = imbin_filt_gfp;
    end

    if params.savegfp
        figure();
        imagesc(max(mask_gfp,[],3));
        options.overwrite=true;
        saveastiff(uint16(mask_gfp),[params.maskfolder fnames(fidx).name(1:end-4) '_gfp_seg.tif'],options);
    end


end
