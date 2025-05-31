// ************* FIJI SCRIPT FOR PRE-PROCESSING MULTI-EXR IMAGES ************* \\
// Last modified by MES on 8/9/24
// Please read through this script and take note of anything marged "CHANGE HERE". Change accordingly
// Lots of this script is hard-coded and will need to be modified based on your file saving/naming conventions
// Some familiarity with Fiji / ijm is best
// I have tried to compensate for this with good commenting

run("Close All");

parent = "A:/Margaret/Astrocytes/Exp174_Aldh1l1-Cre_CAG-FLEX-GFP_18x-exp_multiExR/"; //parent folder where raw microscope images are stored (e.g. in .nd2 format)
parentout = "A:/Margaret/Astrocytes/Exp174_Aldh1l1-Cre_CAG-FLEX-GFP_18x-exp_multiExR/preprocessed/"; //where to save the images

roundList = newArray(//"R1.1/"
					"R1.2/",
					"R2/"
					); //round folders - again, will need to be modified depending on how you saved your data

setBatchMode(true); //true for speed

for (m = 0; m < roundList.length; m++){
	
	folder = parent+roundList[m] + "40x/";
	//print(folder);
	imgList = getFileList(folder);
	
	//regsplits = split(condsplits[1],"/"); // uncomment this if you have region information in the directory structure, which will allow you to forgo hard-coding

	//CHANGE HERE per the naming convention of your files. This assumes "ROI[X] 40x.tif"
	//roundind = indexOf(roundList[j],"R");
	//slashind = indexOf(roundList[j],"/");
	//roundno = substring(roundList[j], roundind+1, slashind);
	roundno = toString(m + 1);
	
	for (k = 0; k < imgList.length; k++){
    	if (endsWith(imgList[k], ".nd2")){ //CHANGE HERE for your filetype
    		if (indexOf(imgList[k],"Max") < 0){ //CHANGE HERE: provide a string that is ONLY in the images you're interested in
    			if (indexOf(imgList[k],"40x") >= 0){
    			//we have "40" because all of our high-mag images are named with "40x"
    			// you can also skip this extra for-loop
    				preprocessImage(imgList[k],folder,roundno);
    			}
			}
		}
	}
}


function preprocessImage(imageFile,folder,roundno)
{
	run("Bio-Formats Importer", "open=[" + folder + imageFile + "] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	//filename = getTitle();
	getDimensions(width, height, channels, slices, frames);

	//subtract off the background
	run("Subtract Background...", "rolling=50");

	dotIndex = indexOf(imageFile, ".");	
	title = substring(imageFile, 0, dotIndex); 
	titlesplit = split(title,"_");
	
	mouseno = titlesplit[0];
	regid = titlesplit[1];
	roino = titlesplit[3];
	
	if (parseInt(roundno) < 10) {
		roundid = "round00" + roundno;
	} else if (parseInt(roundno) >= 10) {
		roundid = "round0" + roundno;
	}
	//roundid="round002"; hard-coded for the round without GFAP and SMI (R1.1)
	//newname = title;
	newname = mouseno + "-" + regid + "-" + "fov"+ roino + "_" + roundid;
	print(newname); //output file 

	saveAs("Tiff", parentout + newname + "_pp.tif");;

	close("*");
	run("Collect Garbage");
}
