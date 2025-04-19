function  data = run_target_quant(params)

%Extract parameters
parentfolder = params.parentfolder;
cavfolder = params.cavsynfolder;

fnames_cav = dir([cavfolder filesep '*' params.cav_channel '*.tif']);
nfovs = length(fnames_cav);
ntargets = length(params.target_channels);

ncoloc = zeros(nfovs,ntargets);
totalvol = zeros(nfovs,ntargets);
meanint_object = zeros(nfovs,ntargets);
meanint = zeros(nfovs,ntargets);
meanvol = zeros(nfovs,ntargets);
enrich_ratio = zeros(nfovs,ntargets);
fovnames = {};
inoverout = zeros(nfovs,ntargets);
enrich_ratio_int = zeros(nfovs,ntargets);

%loop through each field of view in the folder
for fidx = 1:nfovs
    cavimg = loadtiff([cavfolder filesep fnames_cav(fidx).name]);
    cavimg = double(cavimg)/params.imgbit;
    %cavimg = mat2gray(cavimg);

    %% Create Cav2.1 segmentation
    if params.gfilt
        cavimg_blur = imgaussfilt3(cavimg,params.sigma);
    else
        cavimg_blur= cavimg;
    end

    %load the appropriate astrocyte morphology mask and only count objects
    %in the intersection
    splits = strsplit(fnames_cav(fidx).name,'_');
    fovname = splits{1};
    morph_fnames = dir([params.maskfolder filesep '*' fovname '*.tif']);
    morph_fname = morph_fnames(1).name;
    morphimg = logical(loadtiff([params.maskfolder filesep morph_fname]));

    if params.use_adaptive_thresh
        middleslice = cavimg_blur(:,:,ceil(size(cavimg_blur,3)/2));
        middleslice_morph_mask = morphimg(:,:,ceil(size(cavimg_blur,3)/2));
        
        masked_target_middleZ = middleslice_morph_mask .* middleslice;
        invmasked_target_middleZ = ~middleslice_morph_mask .* middleslice;
        
        nonzero_std_in = std(masked_target_middleZ(masked_target_middleZ(:)>0));
        nonzero_std_out = std(invmasked_target_middleZ(invmasked_target_middleZ(:)>0));

        in_over_out = nonzero_std_in/nonzero_std_out;
        
        %std_over_mean = std(middleslice(:))/mean(middleslice(:));
        %sigmaovermu(fidx,1) = std_over_mean;
        inoverout(fidx,1) = in_over_out;
        % 
        % if in_over_out < params.in_over_out_thresh
        %     params.lower_percentile = 98;
        %     params.thresh_multiplier = 10;
        % else
        %     params.lower_percentile = 75;
        %     params.thresh_multiplier = 7;
        % end
        if strcmp(params.thresh_mult_values{1},'hi')
            params.thresh_multiplier = 3;
        else
            params.thresh_multiplier = 1.5;
        end
    end
    
    imbin = binarize_intensity_threshold(cavimg_blur,params);

    if strcmp(params.filt,'med')
        imbin_filt = medfilt3(imbin,params.filt_size);
    else
        imbin_filt = imbin;
    end

    if params.target_sizefilt
        %with help from chatGPT on 12/31/24
        cc = bwconncomp(imbin_filt);
        vols = regionprops3(cc,'Volume');
        valididx = find([vols.Volume] >= params.lowerlim & [vols.Volume] <= params.upperlim);
        mask = ismember(labelmatrix(cc), valididx);
    else
        mask = imbin_filt;
    end

    mask_inastro = morphimg & mask; %take intersection of Cav2.1 and GFP masks

    %quantify number, intensity and volume of target objects
    %colocalized with GFP
    %number
    cc_target = bwconncomp(mask_inastro,26);
    ncoloc(fidx,1) = cc_target.NumObjects;

    %total volume
    totalvol(fidx,1) = nnz(mask_inastro)*params.vol_converter;

    %mean volume
    voltemp = regionprops3(cc_target,'Volume');
    voltemp = voltemp.Volume;
    meanvol(fidx,1) = mean(voltemp)*params.vol_converter;

    %mean intensity = average grayscale value in GFP-masked region for that
    %channel
    masked_target_raw = cavimg .* morphimg;
    nonzero_mean = sum(masked_target_raw(:))/nnz(morphimg);
    meanint(fidx,1) = nonzero_mean;

    %mean intensity object = mean intensity in masked region of target channel
    %(within objects)
    masked_target = cavimg .* mask_inastro;
    nonzero_mean = sum(masked_target(:))/nnz(masked_target);
    meanint_object(fidx,1) = nonzero_mean;

    %enrichment ratio - what is the ratio of volume of target within vs.
    %outside of the astrocyte?
    invmask_target = ~morphimg & mask;
    enrich_ratio(fidx,1) = nnz(mask_inastro)/nnz(invmask_target);

    %intensity enrichment ratio - what is the ratio of raw target
    %signal within vs. outside of the astrocyte?
    invmask_target_int = ~morphimg .* cavimg;
    enrich_ratio_int(fidx,1) = sum(masked_target_raw(:))/sum(invmask_target_int(:));

    %save Cav2.1 mask
    if params.savemasks
        options.overwrite=true;
        saveastiff(uint8(mask),[params.targetoutfolder fovname '_Cav2.1_seg.tif'],options);
        saveastiff(uint8(mask_inastro),[params.targetoutfolder fovname '_Cav2.1_GFP_seg.tif'],options);
    end

    %% Create segmentations for other target channels
    for chidx = 2:ntargets
        chid = params.target_channels{chidx};
        target_splits = strsplit(chid,'-');
        fnames_target = dir([parentfolder filesep '*round00' target_splits{1} '*ch0' target_splits{2} '*.tif']);
        fname_target_idx = find(contains({fnames_target.name}, fovname)==1);

    if ~isempty(fname_target_idx)
        fname_target = fnames_target(fname_target_idx).name;
    
        targetimg = loadtiff([parentfolder filesep fname_target]);
        targetimg = double(targetimg)/params.imgbit;
        %targetimg = mat2gray(targetimg);

        if params.gfilt
            targetimg_blur = imgaussfilt3(targetimg,params.sigma);
        else
            targetimg_blur = targetimg;
        end
    
        if params.use_adaptive_thresh
            middleslice = targetimg_blur(:,:,ceil(size(targetimg_blur,3)/2));
            
            masked_target_middleZ = middleslice_morph_mask .* middleslice;
            invmasked_target_middleZ = ~middleslice_morph_mask .* middleslice;
            
            nonzero_std_in = std(masked_target_middleZ(masked_target_middleZ(:)>0));
            nonzero_std_out = std(invmasked_target_middleZ(invmasked_target_middleZ(:)>0));
    
            in_over_out = nonzero_std_in/nonzero_std_out;
            
            %std_over_mean = std(middleslice(:))/mean(middleslice(:));
            %sigmaovermu(fidx,1) = std_over_mean;
            inoverout(fidx,chidx) = in_over_out;

            % if in_over_out < params.in_over_out_thresh
            %     params.lower_percentile = 98;
            %     params.thresh_multiplier = 10;
            % else
            %     params.lower_percentile = 75;
            %     params.thresh_multiplier = 7;
            % end
            if strcmp(params.thresh_mult_values{chidx},'hi')
                params.thresh_multiplier = 3;
            else
                params.thresh_multiplier = 1.5;
            end
        end
    
        %binarize target channel
        target_imbin = binarize_intensity_threshold(targetimg_blur,params);

        if strcmp(params.filt,'med')
            target_imbin_filt = medfilt3(target_imbin,params.filt_size);
        else
            target_imbin_filt = target_imbin;
        end

        if params.target_sizefilt
            %with help from chatGPT on 12/31/24
            cc = bwconncomp(target_imbin_filt);
            vols = regionprops3(cc,'Volume');
            valididx = find([vols.Volume] >= params.lowerlim & [vols.Volume] <= params.upperlim);
            targetmask = ismember(labelmatrix(cc), valididx);
        else
            targetmask = target_imbin_filt;
        end
        
        %quantify number, intensity and volume of target objects
        %colocalized with GFP
        %number
        mask_inastro = morphimg & targetmask;
        cc_target = bwconncomp(mask_inastro,26);
        ncoloc(fidx,chidx) = cc_target.NumObjects;

        %total volume
        totalvol(fidx,chidx) = nnz(mask_inastro)*params.vol_converter;

        %mean volume
        voltemp = regionprops3(cc_target,'Volume');
        voltemp = voltemp.Volume;
        meanvol(fidx,chidx) = mean(voltemp)*params.vol_converter;

        %mean intensity = average grayscale value in GFP-masked region for that
        %channel
        masked_target_raw = targetimg .* morphimg;
        nonzero_mean = sum(masked_target_raw(:))/nnz(morphimg);
        meanint(fidx,chidx) = nonzero_mean;
    
        %mean intensity object = mean intensity in masked region of target channel
        %(within objects)
        masked_target = targetimg .* mask_inastro;
        nonzero_mean = sum(masked_target(:))/nnz(masked_target);
        meanint_object(fidx,chidx) = nonzero_mean;

        %enrichment ratio - what is the ratio of volume of target within vs.
        %outside of the astrocyte?
        invmask_target = ~morphimg & targetmask;
        enrich_ratio(fidx,chidx) =  nnz(mask_inastro)/nnz(invmask_target);

        %intensity enrichment ratio - what is the ratio of raw target
        %signal within vs. outside of the astrocyte?
        invmask_target_int = ~morphimg .* targetimg;
        enrich_ratio_int(fidx,chidx) = sum(masked_target_raw(:))/sum(invmask_target_int(:));

        fovnames{fidx,1} = fovname;

        %save target mask
        if params.savemasks
            options.overwrite=true;
            targetname = params.target_names{chidx};
            if params.use_adaptive_thresh
                disp([targetname ' in astro std: out astro std = ' num2str(in_over_out)])
            end
            saveastiff(uint8(targetmask),[params.targetoutfolder fovname '_' targetname '_seg.tif'],options);
            saveastiff(uint8(mask_inastro),[params.targetoutfolder fovname '_' targetname '_GFP_seg.tif'],options);
        end

    else
        disp([fovname ' does not have the requisite rounds - skipping']);
        
        ncoloc(fidx,chidx) = NaN;
        totalvol(fidx,chidx) = NaN;
        meanvol(fidx,chidx) = NaN;
        meanint_object(fidx,chidx) = NaN;
        meanint(fidx,chidx)=NaN;
        fovnames{fidx,1} = fovname;
        inoverout(fidx,chidx) = NaN;
        enrich_ratio(fidx,chidx) = NaN;
        enrich_ratio_int(fidx,chidx) = NaN;

    end
    end
end

data.ncoloc = ncoloc;
data.totalvol = totalvol;
data.meanvol = meanvol;
data.meanint_object = meanint_object;
data.meanint = meanint;
data.fovnames = fovnames;
data.inoverout = inoverout;
data.enrich_ratio = enrich_ratio;
data.enrich_ratio_int = enrich_ratio_int;

