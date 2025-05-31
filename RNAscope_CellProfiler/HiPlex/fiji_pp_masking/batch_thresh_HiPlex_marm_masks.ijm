/// SCRIPT for creating lipofuscin autofluorescence masks based on minimum intensity projections ////

parent = "/media/mschro/hd3/Glia_SingleCell/Data/HiPlex/marm/grouped/adult_flipped/masks/min_proj/"; //change here
target = "/media/mschro/hd3/Glia_SingleCell/Data/HiPlex/marm/grouped/adult_flipped/masks/bin/"; //change here

run("Close All");
image_folder = parent;
setBatchMode(true);

list = getFileList(image_folder);									 //find all files in current data folder					
extension = "tif";							 //set file types to look at 'tif' stacks

for (i=0; i<list.length; i++) {								 //go through all files in the current data folder
	if (endsWith(toLowerCase(list[i]), "." + extension)){ 
		name =  list[i];				//record the full file name 
		dotIndex = indexOf(name, ".tif");	//get location of file extension
		title = substring(name, 0, dotIndex);  //get abbreviated filename
		
		open(image_folder + list[i]);
		getStatistics(area, mean, min, max, std, histogram); //get image statistics, including mean and std, which we will use for thresholding
		run("Gaussian Blur...", "sigma=2"); //run gaussian blurring with sigma = 2
		
		lowerlim=mean + 1.5*std; //threshold at X (change the number) standard deviations above the mean
		upperlim=max; //upper limit is maximum intensity value
		setThreshold(lowerlim,upperlim);
		
		run("Convert to Mask");
		run("Median...", "radius=2"); //get rid of really small puncta
		run("8-bit");
		
		saveAs("Tiff", target + title +  "_bin.tif");

		run("Close All");
	}
}
run("Close All");

run("Collect Garbage");