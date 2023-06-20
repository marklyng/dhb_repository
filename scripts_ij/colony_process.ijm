/*
 * Macro template to process multiple images in a folder
 */

#@ File (label = "Input directory", style = "directory") input
//#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".czi") suffix

// See also Process_Folder.py for a version of this code
// in the Python scripting language.
setBatchMode(true);

processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, list[i]);
	}
}

function processFile(input, file) {
	// Do the processing here by adding your own code.
	// Leave the print statements until things work, then remove them.

	run("Bio-Formats Importer", 
		"open=["+input + File.separator + file+"] autoscale color_mode=Composite open_files rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT"); //Import files
	
	output_parent = File.getParent(input);
	Table.open(output_parent + "/contrast.txt");
	
	mag_min = Table.get("value", 0);
	mag_max = Table.get("value", 1);
	gfp_min = Table.get("value", 2);
	gfp_max = Table.get("value", 3);
	grey_min = Table.get("value", 4);
	grey_max = Table.get("value", 5);
	
	
/*   	Table.open(output_parent + "/crop.txt");
   	dim0 = Table.get("dim", 0);
	dim1 = Table.get("dim", 1);	
	dim2 = Table.get("dim", 2);
	dim3 = Table.get("dim", 3);
*/
		
	title = getTitle();
	title_noext = getTitleStripExtension();
	getPixelSize(unit, pixelWidth, pixelHeight);
	scale = 1/pixelWidth*1000;
	run("Set Scale...", "distance="+scale+" known=1 unit=mm");
	getDimensions(width, height, channels, slices, frames);

	if(channels <= 1) {
		setMinAndMax(grey_min, grey_max);
		run("Grays");
		} else if(channels == 2) {
			if(matches(title_noext, ".*_bs.*")) {
		
				setMinAndMax(mag_min, mag_max);
				run("Magenta");
				run("Next Slice [>]");
				setMinAndMax(grey_min, grey_max);
				run("Grays");
		} else if(matches(title_noext, ".*_ps.*")) {
				
				setMinAndMax(gfp_min, gfp_max);
				run("Green");
				run("Next Slice [>]");
				setMinAndMax(grey_min, grey_max);
				run("Grays");			
		} else if(matches(title_noext, ".*_co.*")) {
				setMinAndMax(gfp_min, gfp_max);
				run("Green");
				run("Next Slice [>]");
				setMinAndMax(grey_min, grey_max);
				run("Grays");
		}
		} else if(channels == 3) {
			setMinAndMax(mag_min, mag_max);
			run("Magenta");
			run("Next Slice [>]");
			setMinAndMax(gfp_min, gfp_max);
			run("Green");
			run("Next Slice [>]");
			setMinAndMax(grey_min, grey_max);
			run("Grays");
			}
	
			input_dir = File.getName(input);
			output = output_parent + File.separator + input_dir + "_processed";
			if (!File.exists(output)) {
				File.makeDirectory(output);
			}
		

// Save as .jpeg - can be changed to "png" or "tiff"
    while(nImages > 0){
//    	run("RGB Color");
//    	makeRectangle(dim0, dim1, dim2, dim3);
//		run("Crop");
	title = getTitleStripExtension();
	saveAs("tiff", output + File.separator + title);
	close("*");
	}

}

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