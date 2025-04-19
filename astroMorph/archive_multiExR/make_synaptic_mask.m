function  data = make_synaptic_mask(params)

%Extract parameters
parentfolder = params.parentfolder;
cavfolder = params.cavsynfolder;

fnames_cav = dir([cavfolder filesep '*' params.cav_channel '*.tif']);
psd95_splits = strsplit(params.syn_channels{1},'-');
geph_splits = strsplit(params.syn_channels{2},'-');

fnames_psd95 = dir([parentfolder filesep '*round00' psd95_splits{1} '*ch0' psd95_splits{2} '*.tif']);
fnames_geph = dir([parentfolder filesep '*round00' geph_splits{1} '*ch0' geph_splits{2} '*.tif']);

nfovs = length(fnames_cav);

nsyn_excit = zeros(nfovs,1);
nsyn_inh = zeros(nfovs,1);
nsyn_filtered= zeros(nfovs,1);
fovnames = {};

%loop through each field of view in the folder
for fidx = 1:nfovs
    cavimg = loadtiff([cavfolder filesep fnames_cav(fidx).name]);
    cavimg = double(cavimg);

    %% Create Cav2.1 segmentation
    if params.gfilt
        cavimg_blur = imgaussfilt3(cavimg,params.sigma);
    else
        cavimg_blur= cavimg;
    end

    %Use Cav2.1 channel to segment synapses
    thresh_old = params.thresh_multiplier;
    params.thresh_multiplier = params.cav_thresh_multiplier;
    imbin = binarize_intensity_threshold(cavimg_blur,params);
    params.thresh_multiplier = thresh_old;

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

    %% Create PSD95 and gephyrin segmentations
    splits = strsplit(fnames_cav(fidx).name,'_');
    fovname = splits{1};

    fname_psd95_idx = find(contains({fnames_psd95.name}, fovname)==1);
    fname_geph_idx = find(contains({fnames_geph.name}, fovname)==1);

    if ~isempty(fname_psd95_idx) && ~isempty(fname_geph_idx)
        fname_psd95 = fnames_psd95(fname_psd95_idx).name;
        fname_geph = fnames_geph(fname_geph_idx).name;
    
        psd95img = loadtiff([parentfolder filesep fname_psd95]);
        psd95img = mat2gray(psd95img);
    
        gephimg = loadtiff([parentfolder filesep fname_geph]);
        gephimg = mat2gray(gephimg);
    
        if params.gfilt
            psd95img_blur = imgaussfilt3(psd95img,params.sigma);
            gephimg_blur = imgaussfilt3(gephimg,params.sigma);
        else
            psd95img_blur = psd95img;
            gephimg_blur = gephimg;
        end
    
        %binarize psd95 and gephyrin channels
        if params.use_adaptive_thresh %use adaptive threshold?
            if mean(psd95img_blur(:)) >= 3e-3 %if we have a higher-intensity image, SNR usually higher, use lower threshold
                params.thresh_multiplier = 3;
                psd95_imbin = binarize_intensity_threshold(psd95img_blur,params);
            else
                params.thresh_multiplier = 6;
                psd95_imbin = binarize_intensity_threshold(psd95img_blur,params);
                disp(['PSD95 in low threshold'])
            end

            if mean(gephimg_blur(:)) >= 3e-3%if we have a higher-intensity image, SNR usually higher, use lower threshold
                params.thresh_multiplier = 3;
                geph_imbin = binarize_intensity_threshold(gephimg_blur,params);
                disp(['Gephyrin in low threshold'])
            else
                params.thresh_multiplier = 6;
                geph_imbin = binarize_intensity_threshold(gephimg_blur,params);
            end

        else
            params.thresh_multiplier=3;
            psd95_imbin = binarize_intensity_threshold(psd95img_blur,params);
            params.thresh_multiplier=5; %Gephyrin channel is super noisy, need to increase this
            geph_imbin = binarize_intensity_threshold(gephimg_blur,params);
        end

        if strcmp(params.filt,'med')
            psd95_imbin_filt = medfilt3(psd95_imbin,params.filt_size);
            geph_imbin_filt = medfilt3(geph_imbin,params.filt_size);
        else
            psd95_imbin_filt = psd95_imbin;
            geph_imbin_filt = geph_imbin;
        end
    
        %% Merge the synaptic masks
        combined_mask1 = mask & psd95_imbin_filt; %excitatory synapses
        combined_mask2 = mask & geph_imbin_filt;
    
        cc1 = bwconncomp(combined_mask1);
        nsyn_excit(fidx,1) = cc1.NumObjects;
        cc2 = bwconncomp(combined_mask2);
        nsyn_inh(fidx,1) = cc2.NumObjects;
    
        combined_mask = combined_mask1 | combined_mask2;
    
        if params.doplot
            figure;
            subplot(1,2,1)
            imagesc(max(mask,[],3))
            title('All Cav2.1')
            subplot(1,2,2)
            imagesc(max(combined_mask,[],3))
            title('Cav2.1 + PSD95 or Gephyrin')
        end
    
        %% Create bounding boxes
       
        %get rid of any small puncta remaining
        combined_mask = bwareaopen(combined_mask, params.lowerlim);
        
        CC_combined = bwconncomp(combined_mask,26);%find all connected components in the combined image
        nsyn_filtered_temp = CC_combined.NumObjects;
        bboxes_pos = regionprops3(CC_combined,'BoundingBox');
        nsyn_filtered(fidx,1) = nsyn_filtered_temp;
        
        fovnames{fidx,1} = fovname;
        
       %From MATLAB documentation: Smallest cuboid containing the region,
        % returned as a 1-by-6 vector of the form [ulf_x ulf_y ulf_z width_x width_y width_z].
        % ulf_x, ulf_y, and ulf_z specify the upper-left front corner of the cuboid.
        % width_x, width_y, and width_z specify the width of the cuboid along each dimension.
    
        % Below with help from ChatGPT 4.0 by MES on 7/7/24
        segimg_combined = false(size(combined_mask));
    
        buffer_half = params.nbuffer/2;
        zbuffer_half = 2;%floor(buffer_half/2);
        for k = 1:nsyn_filtered_temp
            bbox = bboxes_pos(k,1).BoundingBox;
            ulf_x = bbox(1);
            ulf_y = bbox(2);
            ulf_z = bbox(3);
            width_x = bbox(4);
            width_y = bbox(5);
            width_z = bbox(6);
    
            % Calculate the start and end indices for each dimension
            xStart = max(1, floor(ulf_x) - buffer_half);
            xEnd = min(size(combined_mask,2), floor(ulf_x + width_x) + buffer_half);
            yStart = max(1, floor(ulf_y) - buffer_half);
            yEnd = min(size(combined_mask,1), floor(ulf_y + width_y) + buffer_half);
            zStart = max(1, floor(ulf_z) - zbuffer_half);
            zEnd = min(size(combined_mask,3), floor(ulf_z + width_z) + zbuffer_half);
            
            % Set the bounding box region to 1 (white)
            segimg_combined(yStart:yEnd, xStart:xEnd, zStart:zEnd) = true;
        end
    
        if params.savemasks
            options.overwrite=true;
            saveastiff(uint8(combined_mask),[params.outfolder fovname '_syn_seg.tif'],options);
            saveastiff(uint8(segimg_combined),[params.outfolder fovname '_syn_bbox_seg.tif'],options);
        end

    else
        disp([fovname ' does not have the requisite rounds - skipping']);
        nsyn_excit(fidx,1) = NaN;
        nsyn_inh(fidx,1) = NaN;
        nsyn_filtered(fidx,1) = NaN;
        fovnames{fidx,1} = fovname;
    end

end

data.nsyn_inh = nsyn_inh;
data.nsyn_excit = nsyn_excit;
data.nsyn_filtered = nsyn_filtered;
data.fovnames = fovnames;

