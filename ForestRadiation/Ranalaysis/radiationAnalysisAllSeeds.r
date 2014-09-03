#!/usr/bin/env Rscript
args <- commandArgs()
commandArgs <- function() args[1]
source("read_files.r")
rootpath = "/home/likewise-open/IN/gopalkri/Developer/core-model/ForestRadiation/Resultapp"
path = paste(rootpath,"/Radiation",as.character(numOfTrees),"numParts_",as.character(args[1]),"/","forAllSeeds",sep ='')
dir.create(path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
cnt = 1
for(seedInt in seeds){
   
dirData      = cbind(get(paste("dirDataVoxBox",as.character(cnt),sep='')))
meanData     = cbind(get(paste("meanDataVoxBox",as.character(cnt),sep='')))
accurateData = cbind(get(paste("accurateDataVoxBox",as.character(cnt),sep='')))
meanBDYesData= cbind(get(paste("meanBDYes",as.character(cnt),sep ='')))



cnt = cnt +1
}


custom_plotting<-function(accurateDataVoxBox1,meanDataVoxBox1,dirDataVoxBox1,titleString){

matplot(accurateDataVoxBox1[7],meanDataVoxBox1$V7[meanDataVoxBox1$V1== 0.1],type="p",pch = 1,col = 2:6,xlab = "Q accurate (MJ/m2)", ylab = " Q mean(MJ/m2)" , xlim= c(min(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.1],accurateDataVoxBox1$V7), max(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.1],accurateDataVoxBox1$V7)),ylim=c(min(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.1],accurateDataVoxBox1$V7), max(meanDataVoxBox1$V7[meanDataVoxBox1$V1== 0.1],accurateDataVoxBox1$V7)),main = "Radiation mean vs accurate Vox = 0.1 m", )
abline(0,1)

matplot(accurateDataVoxBox1[7],dirDataVoxBox1$V7[dirDataVoxBox1$V1== 0.1],type="p",pch = 1,col = 2:6,xlab = "Q accurate (MJ/m2)", ylab = " Q dir(MJ/m2)" ,  xlim= c(min(dirDataVoxBox1$V7[dirDataVoxBox1$V1== 0.1],accurateDataVoxBox1$V7), max(dirDataVoxBox1$V7[dirDataVoxBox1$V1== 0.1],accurateDataVoxBox1$V7)),ylim=c(min(dirDataVoxBox1$V7[dirDataVoxBox1$V1== 0.1],accurateDataVoxBox1$V7), max(dirDataVoxBox1$V7[dirDataVoxBox1$V1== 0.1],accurateDataVoxBox1$V7)),main = "Radiation directional vs accurate Vox = 0.1 m", )
abline(0,1)

matplot(accurateDataVoxBox1[7],meanDataVoxBox1$V7[meanDataVoxBox1$V1== 0.2],type="p",pch = 1,col = 2:6,xlab = "Q accurate (MJ/m2)", ylab = " Q mean(MJ/m2)" , xlim= c(min(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.2],accurateDataVoxBox1$V7), max(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.2],accurateDataVoxBox1$V7)),ylim=c(min(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.2],accurateDataVoxBox1$V7), max(meanDataVoxBox1$V7[meanDataVoxBox1$V1== 0.2],accurateDataVoxBox1$V7)),main = "Radiation mean vs accurate Vox = 0.2 m", )
abline(0,1)

matplot(accurateDataVoxBox1[7],dirDataVoxBox1$V7[dirDataVoxBox1$V1== 0.2],type="p",pch = 1,col = 2:6,xlab = "Q accurate (MJ/m2)", ylab = " Q dir(MJ/m2)" ,  xlim= c(min(dirDataVoxBox1$V7[dirDataVoxBox1$V1== 0.2],accurateDataVoxBox1$V7), max(dirDataVoxBox1$V7[dirDataVoxBox1$V1== 0.2],accurateDataVoxBox1$V7)),ylim=c(min(dirDataVoxBox1$V7[dirDataVoxBox1$V1== 0.2],accurateDataVoxBox1$V7), max(dirDataVoxBox1$V7[dirDataVoxBox1$V1== 0.2],accurateDataVoxBox1$V7)),main = "Radiation directional vs accurate Vox = 0.2 m", )
abline(0,1)

matplot(accurateDataVoxBox1[7],meanDataVoxBox1$V7[meanDataVoxBox1$V1== 0.3],type="p",pch = 1,col = 2:6,xlab = "Q accurate (MJ/m2)", ylab = " Q mean(MJ/m2)" , xlim= c(min(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.3],accurateDataVoxBox1$V7), max(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.3],accurateDataVoxBox1$V7)),ylim=c(min(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.3],accurateDataVoxBox1$V7), max(meanDataVoxBox1$V7[meanDataVoxBox1$V1== 0.3],accurateDataVoxBox1$V7)),main = "Radiation mean vs accurate Vox = 0.3 m", )
abline(0,1)

matplot(accurateDataVoxBox1[7],dirDataVoxBox1$V7[dirDataVoxBox1$V1== 0.3],type="p",pch = 1,col = 2:6,xlab = "Q accurate (MJ/m2)", ylab = " Q dir(MJ/m2)" ,  xlim= c(min(dirDataVoxBox1$V7[dirDataVoxBox1$V1== 0.3],accurateDataVoxBox1$V7), max(dirDataVoxBox1$V7[dirDataVoxBox1$V1== 0.3],accurateDataVoxBox1$V7)),ylim=c(min(dirDataVoxBox1$V7[dirDataVoxBox1$V1== 0.3],accurateDataVoxBox1$V7), max(dirDataVoxBox1$V7[dirDataVoxBox1$V1== 0.3],accurateDataVoxBox1$V7)),main = "Radiation directional vs accurate Vox = 0.3 m", )
abline(0,1)

matplot(accurateDataVoxBox1[7],meanDataVoxBox1$V7[meanDataVoxBox1$V1== 0.4],type="p",pch = 1,col = 2:6,xlab = "Q accurate (MJ/m2)", ylab = " Q mean(MJ/m2)" , xlim= c(min(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.4],accurateDataVoxBox1$V7), max(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.4],accurateDataVoxBox1$V7)),ylim=c(min(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.4],accurateDataVoxBox1$V7), max(meanDataVoxBox1$V7[meanDataVoxBox1$V1== 0.4],accurateDataVoxBox1$V7)),main = "Radiation mean vs accurate Vox = 0.4 m", )
abline(0,1)

matplot(accurateDataVoxBox1[7],dirDataVoxBox1$V7[dirDataVoxBox1$V1== 0.4],type="p",pch = 1,col = 2:6,xlab = "Q accurate (MJ/m2)", ylab = " Q dir(MJ/m2)" ,  xlim= c(min(dirDataVoxBox1$V7[dirDataVoxBox1$V1== 0.4],accurateDataVoxBox1$V7), max(dirDataVoxBox1$V7[dirDataVoxBox1$V1== 0.4],accurateDataVoxBox1$V7)),ylim=c(min(dirDataVoxBox1$V7[dirDataVoxBox1$V1== 0.4],accurateDataVoxBox1$V7), max(dirDataVoxBox1$V7[dirDataVoxBox1$V1== 0.4],accurateDataVoxBox1$V7)),main = "Radiation directional vs accurate Vox = 0.4 m", )
abline(0,1)
title(main= titleString,outer=T)
dev.off()
}
#***********************************************************************************************************************************

custom_plotting2<-function(accurateDataVoxBox1,meanDataVoxBox1,titleString){
#Function to plot the boxdir = yes part here no directional calculations are performed.

matplot(accurateDataVoxBox1[7],meanDataVoxBox1$V7[meanDataVoxBox1$V1== 0.1],type="p",pch = 1,col = 2:6,xlab = "Q accurate (MJ/m2)", ylab = " Q mean(MJ/m2)" , xlim= c(min(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.1],accurateDataVoxBox1$V7), max(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.1],accurateDataVoxBox1$V7)),ylim=c(min(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.1],accurateDataVoxBox1$V7), max(meanDataVoxBox1$V7[meanDataVoxBox1$V1== 0.1],accurateDataVoxBox1$V7)),main = "BoxDir Yes mean vs accurate Vox = 0.1 m", )
abline(0,1)

matplot(accurateDataVoxBox1[7],meanDataVoxBox1$V7[meanDataVoxBox1$V1== 0.2],type="p",pch = 1,col = 2:6,xlab = "Q accurate (MJ/m2)", ylab = " Q mean(MJ/m2)" , xlim= c(min(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.2],accurateDataVoxBox1$V7), max(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.2],accurateDataVoxBox1$V7)),ylim=c(min(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.2],accurateDataVoxBox1$V7), max(meanDataVoxBox1$V7[meanDataVoxBox1$V1== 0.2],accurateDataVoxBox1$V7)),main = "BoxDir Yes mean vs accurate Vox = 0.2 m", )
abline(0,1)

matplot(accurateDataVoxBox1[7],meanDataVoxBox1$V7[meanDataVoxBox1$V1== 0.3],type="p",pch = 1,col = 2:6,xlab = "Q accurate (MJ/m2)", ylab = " Q mean(MJ/m2)" , xlim= c(min(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.3],accurateDataVoxBox1$V7), max(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.3],accurateDataVoxBox1$V7)),ylim=c(min(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.3],accurateDataVoxBox1$V7), max(meanDataVoxBox1$V7[meanDataVoxBox1$V1== 0.3],accurateDataVoxBox1$V7)),main = "BoxDir Yes mean vs accurate Vox = 0.3 m", )
abline(0,1)

matplot(accurateDataVoxBox1[7],meanDataVoxBox1$V7[meanDataVoxBox1$V1== 0.4],type="p",pch = 1,col = 2:6,xlab = "Q accurate (MJ/m2)", ylab = " Q mean(MJ/m2)" , xlim= c(min(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.4],accurateDataVoxBox1$V7), max(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.4],accurateDataVoxBox1$V7)),ylim=c(min(meanDataVoxBox1$V7[meanDataVoxBox1$V1==0.4],accurateDataVoxBox1$V7), max(meanDataVoxBox1$V7[meanDataVoxBox1$V1== 0.4],accurateDataVoxBox1$V7)),main = "BoxDir Yes mean vs accurate Vox = 0.4 m", )
abline(0,1)

title(main= titleString,outer=T)
dev.off()
}
#**********************************************************************************************************************************

pdf(paste(rootpath,"/Radiation",as.character(numOfTrees),"numParts_",as.character(args[1]),"/","forAllSeeds","/","ForAllSeeds",".","pdf",sep = ''))
attach(mtcars)
par(mfrow=c(4,2))
custom_plotting(accurateData,meanData,dirData,"For all Seeds")

path = paste(rootpath,"/Radiation",as.character(numOfTrees),"numParts_",as.character(args[1]),"/","forAllSeeds","/","bdYes",sep='')
dir.create(path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
pdf(paste(rootpath,"/Radiation",as.character(numOfTrees),"numParts_",as.character(args[1]),"/","forAllSeeds","/","bdYes/","bdYes","AllSeedYMean",".pdf",sep = ''))
attach(mtcars)
par(mfrow=c(2,2))
custom_plotting2(accurateData,meanBDYesData,"AllSeedYMean")



