#!/usr/bin/env Rscript
# Plotting Script for all the seeds together
source("read_files.r")
seeds =c(787237, 536199,  1676779,  546327,  235663)
path = paste("/home/likewise-open/IN/gopalkri/Developer/core-model/ForestRadiation/Resultapp/Rad",as.character(numOfTrees),"ForAllSeeds", sep ='')
dir.create(path, showWarnings = TRUE, recursive = FALSE, mode = "0777")

dirData      = cbind(dirDataVoxBox1, dirDataVoxBox2,dirDataVoxBox3,dirDataVoxBox4,dirDataVoxBox5)
meanData     = cbind(meanDataVoxBox1, meanDataVoxBox2,meanDataVoxBox3,meanDataVoxBox4,meanDataVoxBox5)
accurateData = cbind(accurateDataVoxBox1, accurateDataVoxBox2,accurateDataVoxBox3,accurateDataVoxBox4,accurateDataVoxBox5)
meanBDYesData= cbind(meanBDYes1, meanBDYes2,meanBDYes3,meanBDYes4,meanBDYes5)

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
pdf(paste("Rad",as.character(numOfTrees),"ForAllSeeds/","ForAllSeeds",".pdf",sep = ''))
attach(mtcars)
par(mfrow=c(4,2))
custom_plotting(accurateData,meanData,dirData,"For all Seeds")

path = paste("Rad",as.character(numOfTrees),"ForAllSeeds/bdYes",sep='')
dir.create(path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
pdf(paste("Rad",as.character(numOfTrees),"ForAllSeeds/bdYes/","AllSeedYMean",".pdf",sep = ''))
attach(mtcars)
par(mfrow=c(2,2))
custom_plotting2(accurateData,meanBDYesData,"AllSeedYMean")


