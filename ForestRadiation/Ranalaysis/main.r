#!/usr/bin/env Rscript

numparts = c(3,5)



for(counter in numparts)
{
   commandArgs <- function() counter
 #  source("read_files_general.r")
   source("deviationMain.r") 
   source("radiation_analysis.r")  
   source("radiationAnalysisAllSeeds.r") 
   source("deviationForAllSeeds.r")

}






