# DS code git folder
 
Current files (kind of in the order you would run the analysis from raw data -> results):

1) RegisterAndCalcDFOF_withGFP_v2.ijm   (FIJI macro)
This is the first code I run in my analysis pipeline. It takes a two-channel tif stack and uses one of the channels (usually the cell-GFP channel though this works on calcium channel) for motion correction. It then uses a predetermined set of trial times (these don't change in my entire dataset) to obtain Fbackground and calculates dF/F using this. 

2) FindDS.m  (Matlab code)
This is the second code I run in my analysis pipeline. This code lets you find DS ROIs within a field of view. The inputs are the dF/F movie generated by using the code in 1) and the text file that has the trial direction information (generated during acquisition). Briefly, the code moves around every pixels in a frame, takes an average of its neighbors in 2D, then compute the DS of that region. If it is DS, it uses the vector sum to figure out its preferred direction. Pixels over a DS cells will generally all be DS and have similar preferred directions.

3) lookAtCells.m and lookAtCells_redraw.m  (Matlab code)
This is the third code I run in my analysis pipeline. Once I draw ROIs in FIJI over DS regions using the output of code 2), I run this code which has a GUI that allows the user to cycle through all ROIs. It outputs things like the entire dF/F trace, the mean response to each directions, a tuning curve, and DSI and VS values. I use this code to manually classify cells as ON-OFF and ON DSGCs. Anything else gets classified as BAD (common examples: no responses to any trials, clear noise caused by edge of the FOV effect, very variable responses to repeated trials).

3b) CalculateDS.m  (Matlab function)
This is a standalone function that can be used to calculate DS and VS.

4) mainDS.m (Matlab code)
Using the cell IDs obtained using code 3), this code does many things. For example it will plot the local DS map for this FOV for ON-OFF and ON DSGCs and it will also plot the sorted mean responses to trial directions. More importantly, this code is used to export a matrix called roiInt where every row is a cell and the columns contain the dF/F trace of that cell. This code also exports the text file with direction information. This information is exported to a main folder that will be used to generate a mega table where every row in the table is one cell and the columns contain metadata of that cell. The table ultimately contains every single cell from the entire dataset.






A) shuffleTrialDirections.m
This code will block shuffle trial directions. It used during a permutation test to figure out if cells are statistically significantly DS.




