#@ File (label = "Input directory", style = "directory") input
#@ String (label = "File suffix", value = ".czi") suffix
setBatchMode(true);
run("Set Measurements...", "area display redirect=None decimal=3");

processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i])) {
			processFolder(input + File.separator + list[i]);
	}
		if(endsWith(list[i], suffix)) {
			dir = File.getName(input);
			parent = File.getParent(input);
			processFile(input, list[i]);			
		}
	}
	saveAs("Results", parent + File.separator + dir + "_area.csv");
	run("Clear Results");
}



function processFile(input, file) {
	run("Bio-Formats Importer", 
		"open=["+input + File.separator + file+"] autoscale color_mode=Default open_all_series rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT"); //Import .czi-files


	title = getTitleStripExtension();
	title_ext = getTitle();
	getDimensions(width, height, channels, slices, frames);
		
	if(channels > 2) {
		run("Duplicate...", "title=" + title + "_bs.czi duplicate channels=1");
		selectWindow(title_ext);
		run("Duplicate...", "title=" + title + "_ps.czi duplicate channels=2");
		segment_and_area(title + "_bs.czi");
		segment_and_area(title + "_ps.czi");
		close("*");
		}
		
	else {
		run("Duplicate...", "title=" + title + "_1.czi duplicate channels=1");
		segment_and_area(title + "_1.czi");
		close("*");
	}
	
}


function segment_and_area(img) {
	selectWindow(img);
	run("Gaussian Blur...", "sigma=2");
	setAutoThreshold("Otsu dark");
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("Close-");
	run("Fill Holes");
	run("Analyze Particles...", "size=1000000-Infinity display exclude"); //Remove particles below colony size (if you have some big straggler blobs that "Close" can't remove)
}

// Use this function to strip any number of extensions
// off images.
// Returns the title without the extension.
//====================================
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