// ************* FIJI SCRIPT FOR PRE-PROCESSING EXR IMAGES ************* \\
// Last modified by MES on 3/14/25
// Please read through this script and take note of anything marged "CHANGE HERE". Change accordingly
// Lots of this script is hard-coded and will need to be modified based on your file saving/naming conventions
// Some familiarity with Fiji / ijm is best
// I have tried to compensate for this with good commenting

run("Close All");

parent = "/run/user/1000/gvfs/smb-share:server=exr.mit.edu,share=exr/Jinyoung/2025.3.5 Astrocyte/"; //parent folder where raw microscope images are stored (e.g. in .nd2 format)
parentout = "/run/user/1000/gvfs/smb-share:server=exr.mit.edu,share=exr/Margaret/Astrocytes/Exp194_ExR_astro_rdegs/"; //where to save the images

regList = newArray("PFC/",
					"STR/",
					"THAL/"
					); 

chanList = newArray(
	"GFP+Gat3_546+Ca++_633/",
	"GFP+GLAST_546+Ca++_633/",
	"GFP+mGluR3_546+Ca++_633/"
	);
	
mouseList = newArray("mouse1 40x/","mouse2 40x/");

setBatchMode(true); //true for speed

for (m = 0; m < regList.length; m++){
	for (j = 0; j < chanList.length; j++){
		for (i=0; i < mouseList.length; i++){
	
			folder = parent+regList[m] + chanList[j] + mouseList[i];
			//print(folder);
			imgList = getFileList(folder);
			
			//regsplits = split(condsplits[1],"/"); // uncomment this if you have region information in the directory structure, which will allow you to forgo hard-coding
		
			//CHANGE HERE per the naming convention of your files. This assumes "ROI[X] 40x.tif"
			//roundind = indexOf(regList[j],"R");
			//slashind = indexOf(regList[j],"/");
			//roundno = substring(regList[j], roundind+1, slashind);
			roundno = toString(m + 1);
			
			for (k = 0; k < imgList.length; k++){
		    	if (endsWith(imgList[k], ".nd2")){ //CHANGE HERE for your filetype
		    		if (indexOf(imgList[k],"Max") < 0){ //CHANGE HERE: provide a string that is ONLY in the images you're interested in
		    			
	    			regsplits = split(regList[m],"/");
	    			regname = regsplits[0];
	    			
	    			chansplits = split(chanList[j],"+");
	    			targetsplits = split(chansplits[1],"_");
	    			targetname = targetsplits[0];
	    			
	    			mousename = "mouse" + toString(i + 1);
	    			
	    			dotIndex = indexOf(imgList[k], ".");	
					title = substring(imgList[k], 0, dotIndex); 
					titlesplit = split(title,"_");

					newname = mousename + "_" + regname + "_" + targetname + "_fov" + titlesplit[0];
					print(newname); //output file 

	    			// you can also skip this extra for-loop
	    			preprocessImage(imgList[k],folder,newname);

					}
				}
			}
		}
	}
}


function preprocessImage(imageFile,folder,newname)
{
	run("Bio-Formats Importer", "open=[" + folder + imageFile + "] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	//filename = getTitle();
	getDimensions(width, height, channels, slices, frames);

	//subtract off the background
	run("Subtract Background...", "rolling=50");


	saveAs("Tiff", parentout + newname + "_pp.tif");;

	close("*");
	run("Collect Garbage");
}
