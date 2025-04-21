function  data = run_target_quant_exr(params)

%Extract parameters
parentfolder = params.parentfolder;
fnames = dir([parentfolder '*.tif']);
nfovs = length(fnames);

fovnames = {};
meanint = zeros(nfovs,1);
enrich_ratio_int = zeros(nfovs,1);
meanint_cav = zeros(nfovs,1);
enrich_ratio_int_cav = zeros(nfovs,1);

%loop through each field of view in the folder
for fidx = 1:nfovs
    img = loadtiff([parentfolder filesep fnames(fidx).name]);
    targetimg = img(:,:,params.target_channel:3:(end-(3-params.target_channel)));
    targetimg = double(targetimg)/params.imgbit;
    
    cavimg = img(:,:,params.cav_channel:3:(end-(3-params.cav_channel)));
    cavimg = double(cavimg)/params.imgbit;

    img = targetimg;

    %load the appropriate astrocyte morphology mask and only count objects
    %in the intersection
    splits = strsplit(fnames(fidx).name,'_');
    fovname = [splits{1} '_' splits{2} '_' splits{3} '_' splits{4}];
    disp(fovname)
    morph_fnames = dir([params.maskfolder '*' fovname '*.tif']);
    morph_fname = morph_fnames(1).name;
    morphimg = logical(loadtiff([params.maskfolder filesep morph_fname]));

    %% For target protein
    %mean intensity = average grayscale value in GFP-masked region for that
    %channel
    masked_target_raw = img .* morphimg;
    nonzero_mean = sum(masked_target_raw(:))/nnz(morphimg);
    meanint(fidx,1) = nonzero_mean;

    %intensity enrichment ratio - what is the ratio of raw target
    %signal within vs. outside of the astrocyte?
    invmask_target_int = ~morphimg .* img;
    enrich_ratio_int(fidx,1) = mean(masked_target_raw(masked_target_raw>0))/mean(invmask_target_int(invmask_target_int>0));

    %% For synaptic reference (Cav2.1)
    %mean intensity = average grayscale value in GFP-masked region for that
    %channel
    masked_target_raw = cavimg .* morphimg;
    nonzero_mean = sum(masked_target_raw(:))/nnz(morphimg);
    meanint_cav(fidx,1) = nonzero_mean;

    %intensity enrichment ratio - what is the ratio of raw target
    %signal within vs. outside of the astrocyte?
    invmask_target_int = ~morphimg .* cavimg;
    enrich_ratio_int_cav(fidx,1) = mean(masked_target_raw(masked_target_raw>0))/mean(invmask_target_int(invmask_target_int>0));

    fovnames{fidx,1} = fovname;

end


data.meanint_target = meanint;
data.meanint_cav = meanint_cav;
data.fovnames = fovnames;
data.enrich_ratio_int_target = enrich_ratio_int;
data.enrich_ratio_int_cav = enrich_ratio_int_cav;

