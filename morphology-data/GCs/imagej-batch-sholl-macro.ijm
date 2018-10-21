/*
 * Macro template to process multiple images in a folder
 */

input = "C:/Users/Justas/Downloads/morphology/all-gcs";
output = "C:/Users/Justas/Downloads/sholl/gcs";
suffix = ".swc";

processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
}

function processFile(input, output, file) {
	// Do the processing here by adding your own code.
	// Leave the print statements until things work, then remove them.
	print("Processing: " + input + File.separator + file);
	dir=input;
	run("Sholl Analysis (Tracings)...", "traces/(e)swc=["+dir+"/"+file+"] image=[] load center=[Start of main path] radius=10 enclosing=1 #_primary=0 infer linear polynomial=[Best fitting degree] normalizer=Area/Volume directory=[]");
	selectWindow(file+" no-filtering (0.000,0.000,0.000)_Sholl-Profiles");
	saveAs("Results",output+"/"+file+".csv");
	close("*");
	selectWindow(file+".csv"); 
	run("Close"); 
	print("Saving to: " + output);
}
