title = getTitleStripExtension();
rename(title);

run("Deinterleave", "how=2 keep");

selectWindow(title+ " #2");

run("Duplicate...", "duplicate");
selectWindow(title + " #2-1");

run("Mean 3D...", "x=0 y=0 z=50");

run("Merge Channels...", "c1=["+title+" #1] c2=["+title+" #2] c3=["+title+" #2-1] create");

selectWindow("Composite");
run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
run("Correct 3D drift", "channel=3 multi_time_scale only=0 lowest=1 highest=1");

selectWindow("registered time points");
run("Split Channels");
selectWindow("C1-registered time points");
run("Grays");
selectWindow("C2-registered time points");
run("Grays");
selectWindow("C3-registered time points");
run("Grays");


selectWindow("C2-registered time points");
removeLightFrames();

run("Mean 3D...", "x=0 y=0 z=100");
getDimensions(width, height, channels, slices, frames)

run("Size...", "width="+ width +" height=" + height +" depth=1400 constrain average interpolation=Bilinear");
imageCalculator("Subtract create 32-bit stack", "C2-registered time points","C2-registered time points-1");
imageCalculator("Divide create 32-bit stack", "Result of C2-registered time points","C2-registered time points-1");
run("Size...", "width="+ width +" height=" + height +" depth=700 constrain average interpolation=Bilinear");
run("Enhance Contrast", "saturated=0.35");

selectWindow("C2-registered time points-1");
run("Z Project...", "projection=[Average Intensity]");
run("Enhance Contrast", "saturated=0.35");

selectWindow("C1-registered time points");
removeLightFrames();
run("Z Project...", "projection=[Average Intensity]");
run("Enhance Contrast", "saturated=0.35");

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



function removeLightFrames() { 

	run("Duplicate...", "duplicate");
	
	run("Make Substack...", "delete slices=1338-1370");
	close();
	run("Make Substack...", "delete slices=1282-1314");
	close();
	run("Make Substack...", "delete slices=1226-1259");
	close();
	run("Make Substack...", "delete slices=1171-1203");
	close();
	run("Make Substack...", "delete slices=1115-1147");
	close();
	run("Make Substack...", "delete slices=1059-1091");
	close();
	run("Make Substack...", "delete slices=1003-1036");
	close();
	run("Make Substack...", "delete slices=948-980");
	close();
	run("Make Substack...", "delete slices=892-924");
	close();
	run("Make Substack...", "delete slices=836-868");
	close();
	run("Make Substack...", "delete slices=780-813");
	close();
	run("Make Substack...", "delete slices=725-757");
	close();
	run("Make Substack...", "delete slices=669-701");
	close();
	run("Make Substack...", "delete slices=613-645");
	close();
	run("Make Substack...", "delete slices=557-590");
	close();
	run("Make Substack...", "delete slices=502-534");
	close();
	run("Make Substack...", "delete slices=446-478");
	close();
	run("Make Substack...", "delete slices=390-422");
	close();
	run("Make Substack...", "delete slices=334-367");
	close();
	run("Make Substack...", "delete slices=279-311");
	close();
	run("Make Substack...", "delete slices=223-255");
	close();
	run("Make Substack...", "delete slices=167-199");
	close();
	run("Make Substack...", "delete slices=111-144");
	close();
	run("Make Substack...", "delete slices=56-88");
	close();

}
