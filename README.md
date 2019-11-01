# DS code git folder
 
Current files:
1) FindDS.m
This code lets you find DS ROIs within a field of view. It moves around every pixels in a frame, takes an average of its neighbors in 2D, then compute the DS of that region. If it is DS, it uses the vector sum to figure out its preferred direction. Pixels over a DS cells will generally all be DS and have similar preferred directions.

2) shuffleTrialDirections.m
