//FIJI SCRIPT FOR TAKING MAX PROJECTIONS AND SEPARATING CHANNELS FOR HIPLEX IMAGES
//Last updated by MES August 2024

parent = "/media/mschro/MS_HD3/Astrocytes/Exp179_hiplex_cjadult_tail-order-swap_20240829/renamed_tif/"; // CHANGE HERE
target = "/ssd/Dropbox/Glia_SingleCell/Data/FISH+IHC/Exp179_hiplex_cjadult_roundflip_20240831/MIP_sep_channels/"; //CHANGE HERE

///// DON'T CHANGE BELOW HERE (usually) ////////

run("Close All");
image_folder = parent;
setBatchMode(true);

list = getFileList(image_folder);	 //find all files in current data folder					
extension = "tif";		 //set file types to look at files ending in .tif

for (i=0; i<list.length; i++) {								 //go through all files in the current data folder
	if (endsWith(toLowerCase(list[i]), "." + extension)){ 
		name =  list[i];												//record the file name
		dotIndex = indexOf(name, ".tif");		//find the position of the extension
		title = substring(name, 0, dotIndex); 	//record the abbreviated file name, excluding the extension
		//print(title);
		
		open(image_folder + list[i]);
		run("Z Project...", "projection=[Max Intensity]"); //take maximum intensity projection in Z
		run("Split Channels"); // split channels
	
		saveAs("Tiff", target+ title +  "_647.tif"); //save down channel 1 as 647 (always first)
		close(); // close it
	
		saveAs("Tiff", target + title +  "_546.tif"); //do the same for all the other channels in order
		close();
	
		saveAs("Tiff", target + title +  "_488.tif");
		close();
	
		saveAs("Tiff", target + title +  "_DAPI.tif"); //last channel is alwasy 405 (DAPI)
		close();
		
		run("Close All");
	}
}
run("Close All");

run("Collect Garbage");