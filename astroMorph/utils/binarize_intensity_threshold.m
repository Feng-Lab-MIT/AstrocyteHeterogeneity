function bin = binarize_intensity_threshold(img,params)
%Binarizes an image using an intensity-based threshold

if strcmp(params.thresh_method,'pct')
    threshold = prctile(img(:),params.thresh_pctl);
    bin = imbinarize(img,threshold);
elseif strcmp(params.thresh_method,'zscore')
    meanint = mean(img(:));
    stdint = std(img(:));
    threshold = meanint + params.thresh_multiplier*stdint;
    bin = imbinarize(img,threshold);
elseif strcmp(params.thresh_method,'stdev')
    pct_low = prctile(img(:),params.lower_percentile);
    imgflat = img(:);
    stdint_low = std(imgflat(imgflat<pct_low));
    threshold = params.thresh_multiplier*stdint_low;
    bin = imbinarize(img,threshold);
elseif strcmp(params.thresh_method,'absolute')
    bin = imbinarize(img,params.threshold);
elseif strcmp(params.thresh_method,'otsu_mult')
    grayimg = mat2gray(img);
    otsuthresh = graythresh(grayimg(500:1500,500:1500,10:70)); %take the middle of the image only to avoid edge artifacts (blood vessels, black spaces after registration)
    threshold = otsuthresh*max(img(:)) * params.thresh_multiplier;
    bin = imbinarize(img,threshold);
end
end

