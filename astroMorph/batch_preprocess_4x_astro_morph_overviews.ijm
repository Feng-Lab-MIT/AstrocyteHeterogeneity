// ************* FIJI SCRIPT FOR PRE-PROCESSING 4x expanded EXR ASTROCYTE IMAGES ************* \\
// Last modified by MES 2/27/25
// Please read through this script and take note of anything marged "CHANGE HERE". Change accordingly
// Lots of this script is hard-coded and will need to be modified based on your file saving/naming conventions
// Some familiarity with Fiji / ijm is best
// I have tried to compensate for this with good commenting

run("Close All");
setBatchMode(true);

parent = "A:/Margaret/Astrocytes/Exp189_4xexp_astromorph_202502/4x/"; //parent folder where raw microscope images are stored (e.g. in .nd2 format)
parentout = "A:/Margaret/Astrocytes/Exp189_4xexp_astromorph_202502/4x/GFP_mip/"; //where to save the images
//reg = "ctx"; //specify the brain region here (you can delete this from the final filename in the last section if you want)

	
setBatchMode(true); //true for speed

	    	
folder = parent;
//print(folder);
imgList = getFileList(folder);


for (k = 0; k < imgList.length; k++){
	if (endsWith(imgList[k], ".nd2")){ //CHANGE HERE for your filetype
		if (indexOf(imgList[k],"Max") < 0){ //CHANGE HERE: provide a string that is ONLY in the images you're interested in
			if (indexOf(imgList[k],"4x") >= 0){
			//we have "40" because all of our high-mag images are named with "40x"
			// you can also skip this extra for-loop
				preprocessImage(imgList[k],folder);
			}
		}
	}
}


function preprocessImage(imageFile,folder)
{
	run("Bio-Formats Importer", "open=[" + folder + imageFile + "] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	//filename = getTitle();
	getDimensions(width, height, channels, slices, frames);

	dotIndex = indexOf(imageFile, ".");	
	title = substring(imageFile, 0, dotIndex); 								
	titlesplit = split(title,"x");
	roino = titlesplit[1];
	
	newname = title;
	print(newname); //output file 
	
	//we only want to process the first channel (GFP)
	run("Split Channels");
	
	//Max project and save
	run("Z Project...", "projection=[Max Intensity]");
	saveAs("Tiff", parentout + newname + "_GFP_MIP.tif");
	close();
	
	//close the other channel
	close();
	
	run("Collect Garbage");
}
