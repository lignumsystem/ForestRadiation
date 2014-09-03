#!/usr/bin/env Rscript

numparts = c(3,5) # Please give the required nparts here.



for(counter in numparts) # loop for numparts
{
   commandArgs <- function() counter
 #  source("read_files_general.r")
   source("deviationMain.r")                # deviation for each seeds
   source("radiation_analysis.r")           # radiation analysis for each seeds
   source("radiationAnalysisAllSeeds.r")    # Radiation analysis for all seeds
   source("deviationForAllSeeds.r")         # Deviation for all seeds. 

}






