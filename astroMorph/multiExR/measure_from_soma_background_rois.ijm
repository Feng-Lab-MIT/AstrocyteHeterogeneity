// ************* FIJI SCRIPT FOR AUTOMATICALLY MEASURING PROPERTIES IN SOMA/PROCESS ROIs ************* \\
// Last modified by MES on 1/27/25
// Please read through this script and take note of anything marged "CHANGE HERE". Change accordingly
// Lots of this script is hard-coded and will need to be modified based on your file saving/naming conventions
// Some familiarity with Fiji / ijm is best
// I have tried to compensate for this with good commenting

//CHANGE ALL OF THE BELOW AS NEEDED
parent = "/run/user/1000/gvfs/smb-share:server=exr.mit.edu,share=exr/Margaret/Astrocytes/Exp184_multiExR_astros/registered/singleZ/renamed_files/"; //where whole field of view images are stored
roidir = "/run/user/1000/gvfs/smb-share:server=exr.mit.edu,share=exr/Margaret/Astrocytes/Exp184_multiExR_astros/registered/singleZ/renamed_files/soma_process_rois/"; //where the .zip files for Fiji ROIs are saved
				
image_folder = parent;
roilist = getFileList(roidir);	//all files in roi directory					
imagelist = getFileList(image_folder); //all files in image directory

run("Set Measurements...", "area mean standard display redirect=None decimal=3"); //set measurements

for (j=0; j<roilist.length; j++) {		//loop through all ROIs
	
	//parse the ROI name
	roi_fname = roilist[j];
	fname_splits = split(roi_fname, ".");
	roiname = fname_splits[0]; 
	
	for (i = 0; i<imagelist.length; i++) { //loop through all images (not super efficient but will work)
		if (indexOf(imagelist[i], roiname) >= 0){ //only open images that match the field of view name
					
		    roiManager("Open", roidir + roiname + ".zip"); //open the manually-identified ROIs. NOTE: they must be named the same as the prefix for the .tif stacks
		    
		    open(image_folder + imagelist[i]);
		   
			// Iterate all ROIs in ROI Manager
		    for (n=0; n<roiManager("count"); ++n) {
		    	roiManager("Select", n);
		    	run("Measure");//measure
		    }
			
			selectWindow(imagelist[i]);
		    run("Close");
		    roiManager("Delete");
		    selectWindow("ROI Manager");
		    run("Close");
		}
	}
}
