The R analysis that was written to analysie the data for forest radiation are as follows.

1). read_files.r : First make the changes here. You need to change the path for the files (where they are in the computer for the code to read them) then change the tree density  with numOfTrees and then the seed values depending on the values that is used while running the lignum model.
2).functionsInR.r : This file is a collection of various functions used in the analysis.
3).deviationMainR.r: This file plots all the deviation for all the seeds that is given for the same tree density. Please change the path here. This includes the paths for the scripts and also the output.
4).deviationForAllSeeds.r : This file performs analysis of all the seeds put together.Please change the path here. This includes the paths for the scripts and also the output.
5).radiation_analysis.r : This script plots the radiation i.e. the Qacc vs Qvox for all the seeds mentioned.Please change the path here. This includes the paths for the scripts and also the output.
6). radiationAnalysisAllSeeds.r : Change the path here also as required. Radiation analysis for all the seeds.


