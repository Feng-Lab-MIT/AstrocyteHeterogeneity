// ************* FIJI SCRIPT FOR MAKING MIP JPEGS FOR ASTRO EXR IMAGES ************* \\
// Last modified by MES on 8/8/24
// Please read through this script and take note of anything marged "CHANGE HERE". Change accordingly
// Lots of this script is hard-coded and will need to be modified based on your file saving/naming conventions
// Some familiarity with Fiji / ijm is best
// I have tried to compensate for this with good commenting

run("Close All");
setBatchMode(true);

parent = "A:/Margaret/Astrocytes/Exp177_Aldh1l1-Cre_CAG-FLEX-GFP_4x-exp_astro_morph/thalamus/40x/"; //parent folder where raw microscope images are stored (e.g. in .nd2 format)
parentout = "A:/Margaret/Astrocytes/Exp177_Aldh1l1-Cre_CAG-FLEX-GFP_4x-exp_astro_morph/thalamus/MIP_jpeg/"; //where to save the images

setBatchMode(true); //true for speed

	    	
folder = parent;
//print(folder);
imgList = getFileList(folder);


for (k = 0; k < imgList.length; k++){
	if (endsWith(imgList[k], ".nd2")){ //CHANGE HERE for your filetype
		if (indexOf(imgList[k],"Max") < 0){ //CHANGE HERE: provide a string that is ONLY in the images you're interested in
			if (indexOf(imgList[k],"40x") >= 0){
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
	run("Subtract Background...", "rolling=50 stack");
	run("Z Project...", "projection=[Max Intensity]");
	Stack.setDisplayMode("color");
	Stack.setChannel(1);
	run("Cyan");
	run("Enhance Contrast", "saturated=0.35");
	Stack.setChannel(2);
	run("Magenta");
	run("Enhance Contrast", "saturated=0.35");
	Stack.setChannel(3);
	run("Yellow");
	run("Enhance Contrast", "saturated=0.35");
	Stack.setDisplayMode("composite");
	saveAs("Jpeg", parentout + newname + "_pp_MIP.jpb");
	close();
	close();
	run("Collect Garbage");
}
