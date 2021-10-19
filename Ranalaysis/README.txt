The Folder Ranalaysis contains the main program main.r

Here you can specify the numParts for the simulations by editing  numparts = c(3,5,1) etc.

there are different scripts for different purposes. 
To run the code the user needs to change the path in read_files.r and also the seeds and the tree density.
You can make all this as an input in Main (needs to be done).
To run them individually use, However make sure that the argument is given to these as numparts are given as arguemtns as
e.g.  commandArgs <- function() 1
./deviationMain.r
./radiation_analysis.r
./radiationAnalysisAllSeeds.r
./deviationForAllSeeds.r
./CVPlot.r

else use ./main.r

All the functions are written in functionsInR.r and this script is sourced whenever needed.








