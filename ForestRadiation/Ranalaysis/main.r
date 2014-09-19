#!/usr/bin/env Rscript

numparts = c(3)

rootpath = "/home/likewise-open/IN/gopalkri/Developer/core-model/ForestRadiation/Resultapp"

for(counter in numparts)
{
   commandArgs <- function() counter 
   source("deviationMain.r") 
   source("radiation_analysis.r")  
   source("radiationAnalysisAllSeeds.r") 
   source("deviationForAllSeeds.r")
   source("CVPlot.r")

}






