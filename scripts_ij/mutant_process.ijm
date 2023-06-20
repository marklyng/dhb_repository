#@ File (label = "Input directory", style = "directory") input
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
		if(endsWith(list[i], suffix)) {
			processFile(input, list[i]);
		}
	}
	
}


function processFile(input, file) {
	run("Bio-Formats", "open=[" + input + File.separator + file + "] autoscale color_mode=Composite open_files rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");

	image = getTitleStripExtension();

	run("Set Scale...", "distance=68.4 known=1 unit=mm global");

	makeRectangle(326, 241, 962, 962);
	run("Crop");
	
	setMinAndMax(400, 11000);
	run("Green");
	run("Next Slice [>]");
	setMinAndMax(800, 10000);
	run("Grays");
	
	while(nImages > 0){
    	name = getTitleStripExtension();
    	run("Scale Bar...", "width=3 height=14 font=20 color=White background=None location=[Lower Right] bold hide");
    	saveAs("png", input + File.separator + name);
    	close();
    }
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