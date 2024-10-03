// ************* FIJI SCRIPT FOR RUNNING SHOLL ANALYSIS IN BATCH ************* \\
// Last modified by MES on 8/6/24
run("Close All");
setBatchMode(true);

//CHANGE FILENAMES HERE
//Note: manually add a point and save as an overlay in the approximate center of the Soma before this (see the Sholl plugin instructions on Fiji website)
parent = "A:/Margaret/Astrocytes/Exp177_Aldh1l1-Cre_CAG-FLEX-GFP_4x-exp_astro_morph/thalamus/forsholl/"; //parent folder where raw microscope images are stored (e.g. in .nd2 format)
parentout = "A:/Margaret/Astrocytes/Exp177_Aldh1l1-Cre_CAG-FLEX-GFP_4x-exp_astro_morph/thalamus/sholl_results/"; //where to save the images
	    	
folder = parent;
//print(folder);
imgList = getFileList(folder);
//calculate pixel width/height/depth based on expansion factor
//.1625/4.22 

for (k = 0; k < imgList.length; k++){
	if (endsWith(imgList[k], ".tif")){ //CHANGE HERE for your filetype
		open(parent + imgList[k]);
		//CHANGE to biological (post-expansion units here) - should be different for every batch based on expansion factor
		run("Properties...", "unit=um pixel_width=0.039062500 pixel_height=0.039062500 voxel_depth=0.120192308");
		run("Legacy: Sholl Analysis (From Image)...", "starting=5 ending=70 radius_step=2.5 enclosing=1 #_primary=2 infer linear polynomial=[Best fitting degree] most normalizer=Volume save directory=" + parentout);
		close();
		run("Close All");
		setBatchMode(false);
		run("Collect Garbage");
		setBatchMode(true);
		
	}
}
	
