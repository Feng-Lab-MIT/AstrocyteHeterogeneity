%% Run the box-counting algorithm on astrocyte segmentations (binary)
%Last modified by MES in August 2024

%change directories here!
parentdir = 'A:/Margaret/Astrocytes/Exp177_Aldh1l1-Cre_CAG-FLEX-GFP_4x-exp_astro_morph/striatum/segmentations_1std/';
savedir = 'A:/Margaret/Astrocytes/Exp177_Aldh1l1-Cre_CAG-FLEX-GFP_4x-exp_astro_morph/striatum/boxcount_out/';

files = dir([parentdir '*tif*']);
n = zeros((length(files)),12);
r = zeros((length(files)),12);
mean_df = zeros(length(files),1);

for fidx = 1:length(files)
    img = mat2gray(loadtiff([parentdir files(fidx).name]));
    [ntemp,rtemp] = boxcount(img,'plot');
    df = -diff(log(ntemp))./diff(log(rtemp));
    n(fidx,1:12) = ntemp;
    r(fidx,1:12) = rtemp;
    mean_df(fidx,1)=mean(df);
    savefig([savedir files(fidx).name(1:end-4) '_boxcount.fig'])
end

%change directories here!
save([savedir 'boxcount_workspace.mat'])