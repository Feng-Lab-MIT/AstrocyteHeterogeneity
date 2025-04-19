function  data = make_synaptic_mask_exr(params)

%Extract parameters
parentfolder = params.parentfolder;
cavfolder = params.cavsynfolder;

fnames_cav = dir([cavfolder '*.tif']);

nfovs = length(fnames_cav);
nsyn = zeros(nfovs,1);

fovnames = {};

%loop through each field of view in the folder
for fidx = 1:nfovs
    cavimg = loadtiff([cavfolder filesep fnames_cav(fidx).name]);
    cavimg = cavimg(:,:,params.cav_channel:3:(end-(3-params.cav_channel)));
    cavimg = double(cavimg);

    %% Create Cav2.1 segmentation
    if params.gfilt
        cavimg_blur = imgaussfilt3(cavimg,params.sigma);
    else
        cavimg_blur= cavimg;
    end

    %Use Cav2.1 channel to segment synapses
    imbin = binarize_intensity_threshold(cavimg_blur,params);

    if strcmp(params.filt,'med')
        imbin_filt = medfilt3(imbin,params.filt_size);
    else
        imbin_filt = imbin;
    end

    if params.sizefilt
        %with help from chatGPT on 12/31/24
        cc = bwconncomp(imbin_filt);
        vols = regionprops3(cc,'Volume');
        valididx = find([vols.Volume] >= params.lowerlim & [vols.Volume] <= params.upperlim);
        mask = ismember(labelmatrix(cc), valididx);
    else
        mask = imbin_filt;
    end

    if params.doplot
        figure;
        subplot(1,2,1)
        imagesc(max(cavimg,[],3))
        title('Cav2.1 Original MIP')
        subplot(1,2,2)
        imagesc(max(mask,[],3))
        title('Cav2.1 Segmentation MIP')
    end

    cc = bwconncomp(mask,26);
    nsyn(fidx,1) = cc.NumObjects;

    splits = strsplit(fnames_cav(fidx).name,"_");
    fovname = strcat(splits(1),'-',splits(2),'-',splits(3),'-',splits(4));
    disp(fovname);
    fovnames{fidx,1} = fovname{1};

    if (cc.NumObjects > 0)
        bboxes_pos = regionprops3(cc,'BoundingBox');
        
       %From MATLAB documentation: Smallest cuboid containing the region,
        % returned as a 1-by-6 vector of the form [ulf_x ulf_y ulf_z width_x width_y width_z].
        % ulf_x, ulf_y, and ulf_z specify the upper-left front corner of the cuboid.
        % width_x, width_y, and width_z specify the width of the cuboid along each dimension.
    
        % Below with help from ChatGPT 4.0 by MES on 7/7/24
        segimg_combined = false(size(mask));
    
        buffer_half = params.nbuffer/2;
        zbuffer_half = 2;%floor(buffer_half/2);
        for k = 1:cc.NumObjects
            bbox = bboxes_pos(k,1).BoundingBox;
            ulf_x = bbox(1);
            ulf_y = bbox(2);
            ulf_z = bbox(3);
            width_x = bbox(4);
            width_y = bbox(5);
            width_z = bbox(6);
    
            % Calculate the start and end indices for each dimension
            xStart = max(1, floor(ulf_x) - buffer_half);
            xEnd = min(size(mask,2), floor(ulf_x + width_x) + buffer_half);
            yStart = max(1, floor(ulf_y) - buffer_half);
            yEnd = min(size(mask,1), floor(ulf_y + width_y) + buffer_half);
            zStart = max(1, floor(ulf_z) - zbuffer_half);
            zEnd = min(size(mask,3), floor(ulf_z + width_z) + zbuffer_half);
            
            % Set the bounding box region to 1 (white)
            segimg_combined(yStart:yEnd, xStart:xEnd, zStart:zEnd) = true;
        end
    
        if params.savemasks
            options.overwrite=true;
            outfile1 = strcat(params.outfolder,fovname,'_syn_seg.tif');
            outfile2 = strcat(params.outfolder,fovname,'_syn_bbox_seg.tif');
            saveastiff(uint8(mask),outfile1{1},options);
            saveastiff(uint8(segimg_combined),outfile2{1},options);
        end
    else
        disp(['No objects detected for ' fovname{1}]);
        nsyn(fidx,1) = 0;
    end

end

data.nsyn = nsyn;
data.fovnames = fovnames;

