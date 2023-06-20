#@ File (label = "Input directory", style = "directory") input
#@ String (label = "File suffix", value = ".tif") suffix
#@ String (label = "Pixels per distance", value = "25") pix_len
#@ String (label = "Pixel distance unit", value = "mm") res_unit
setBatchMode(true);
/* Macro to take a parent (input) directory containing directories with .tif-images
 *  
 *  The images will be stacked together based on the center of the images and have values of 0 added to pad the image and make all images the same size
 *  They will then be scaled linearly by the image with the highest intensity in the stack and recieve a calibration bar (LUT scale bar) and a distance scale bar
 *  You must supply: 
 *  	- An input directory containing directories of .tif-files
 *  	- A file suffix (the macro has been written for .tif-files which is the default)
 *  	- A pixel length number (see metaSPACE metadata)
 *  	- A pixel length unit (see metaSPACE metadata)
 *  The output will be an .svg file for each input image, where the LUT of the images is "viridis" and both the calibration and scale bars are vector elements
 *  that can be edited in Inkscape or Illustrator
 *  
 *  
 *  You will need the BioVoxxel figure toolbox https://zenodo.org/badge/latestdoi/542074524
 *  Go to "Help" --> "Update..." --> "Manage update sites" and tick the "BioVoxxel" update site. Click "close", "apply changes" and restart ImageJ.
 *  
 *  
 *  Run the macro with ctrl+r
 *  
 *  Created by Mark Lyng, May 2023
*/

processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		image_list = getFileList(input + File.separator + list[i]);
		processFile(input, list[i], image_list);
		}
}

function processFile(input, parent_dir, file) {
	// Do the processing here by adding your own code.
	// Leave the print statements until things work, then remove them.
	dir = input + File.separator + list[i];
	
	for (i = 0; i < file.length; i++) {
		open(dir + File.separator + file[i]);
	}
	
	run("Set Scale...", "distance=" + pix_len + " known=1 pixel=1 unit=" + res_unit + " global");
	run("Images to Stack", "name=stack fill=Black use");
	run("Duplicate...", "title=blur duplicate");
	run("Mean...", "radius=1 stack");
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	contrast = max * 0.8;
	close("blur");
	
	selectWindow("stack");
	run("mpl-viridis");
	//run("Brightness/Contrast...");
	setMinAndMax(0.00, contrast);
	
	run("Stack to Images");
	for (i = 1; i <= file.length; i++) {
		selectWindow("stack-000" + i);
		rename(file[i-1]);
		run("Calibration Bar...", "location=[Upper Right] fill=None label=White number=5 decimal=0 font=12 zoom=1 bold overlay");
		run("Scale Bar...", "width=3 height=3 thickness=10 font=20 color=White background=None location=[Lower Right] horizontal bold hide overlay");
	}

	output = input + "_viridis";
	output_dir = File.getName(dir);
	
			
	if (!File.exists(output)) {
		File.makeDirectory(output);
		}
		
	if (!File.exists(output + File.separator + output_dir)) {
		File.makeDirectory(output + File.separator + output_dir);
		}	

	
// Save as .svg
    while(nImages > 0){
	title = getTitleStripExtension();
	run("Export SVG", "filename=" + title + " folder=[" + output + File.separator + output_dir + "] exportchannelsseparately=None interpolationrange=0.0 locksensitiverois=false");
	close();
	}
		
	print("Saving to: " + output + output_dir);
}

//close("*");
function getTitleStripExtension() {
  t = getTitle();
  t = replace(t, ".tif", "");        
  t = replace(t, ".tiff", "");      
  t = replace(t, ".lif", "");      
  t = replace(t, ".lsm", "");    
  t = replace(t, ".czi", "");      
  t = replace(t, ".nd2", "");    
  return t;
}