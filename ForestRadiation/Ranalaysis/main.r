#!/usr/bin/env Rscript

numparts = c(3)

#Please change the paths in read_files.r suitable to your machine along with the seeds and the tree density.
#That will be the main code where all the paths are changed.
for(counter in numparts)
{
   commandArgs <- function() counter 
   source("deviationMain.r") 
   source("radiation_analysis.r")  
   source("radiationAnalysisAllSeeds.r") 
   source("deviationForAllSeeds.r")
   source("CVPlot.r")

}






