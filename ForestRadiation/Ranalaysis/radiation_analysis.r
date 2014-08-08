#!/usr/bin/env Rscript
# Plotting script for each seed and for each voxel size.
source("read_files.r")
seeds = c(787237, 536199, 1676779, 546327, 235663)
path = paste("/home/likewise-open/IN/gopalkri/Developer/core-model/ForestRadiation/Resultapp/Rad",as.character(numOfTrees),"New",sep ='') #Change the path here.
dir.create(path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
# Function to plot and save the files in the desired directory. Could be more concised by using a loop but can be improved later.

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


#***********************************************************************************************************************************
pdf(paste("Rad",as.character(numOfTrees),"New/","Seed",as.character(seeds[1]),".pdf",sep = ''))
attach(mtcars)
par(mfrow=c(4,2))
custom_plotting(accurateDataVoxBox1,meanDataVoxBox1,dirDataVoxBox1,paste("Seed", as.character(seeds[1])))

pdf(paste("Rad",as.character(numOfTrees),"New/","Seed",as.character(seeds[2]),".pdf",sep = ''))
attach(mtcars)
par(mfrow=c(4,2))
custom_plotting(accurateDataVoxBox2,meanDataVoxBox2,dirDataVoxBox2,paste("Seed", as.character(seeds[2])))


pdf(paste("Rad",as.character(numOfTrees),"New/","Seed",as.character(seeds[3]),".pdf",sep = ''))
attach(mtcars)
par(mfrow=c(4,2))
custom_plotting(accurateDataVoxBox3,meanDataVoxBox3,dirDataVoxBox3,paste("Seed", as.character(seeds[3])))


pdf(paste("Rad",as.character(numOfTrees),"New/","Seed",as.character(seeds[4]),".pdf",sep = ''))
attach(mtcars)
par(mfrow=c(4,2))
custom_plotting(accurateDataVoxBox4,meanDataVoxBox4,dirDataVoxBox4,paste("Seed", as.character(seeds[4])))


pdf(paste("Rad",as.character(numOfTrees),"New/","Seed",as.character(seeds[5]),".pdf",sep = ''))
attach(mtcars)
par(mfrow=c(4,2))
custom_plotting(accurateDataVoxBox5,meanDataVoxBox5,dirDataVoxBox5,paste("Seed", as.character(seeds[5])))


#**********************************************************************************************************************************
# code to plot the boxdir yes part
path = paste("Rad",as.character(numOfTrees),"New/bdYes",sep='')
dir.create(path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
pdf(paste("Rad",as.character(numOfTrees),"New/bdYes/","SeedYMean",as.character(seeds[1]),".pdf",sep = ''))
attach(mtcars)
par(mfrow=c(2,2))
custom_plotting2(accurateDataVoxBox1,meanBDYes1,paste("Seed", as.character(seeds[1])))

pdf(paste("Rad",as.character(numOfTrees),"New/bdYes/","SeedYMean",as.character(seeds[2]),".pdf",sep = ''))
attach(mtcars)
par(mfrow=c(2,2))
custom_plotting2(accurateDataVoxBox2,meanBDYes2,paste("Seed", as.character(seeds[2])))


pdf(paste("Rad",as.character(numOfTrees),"New/bdYes/","SeedYMean",as.character(seeds[3]),".pdf",sep = ''))
attach(mtcars)
par(mfrow=c(2,2))
custom_plotting2(accurateDataVoxBox3,meanBDYes3,paste("Seed", as.character(seeds[3])))


pdf(paste("Rad",as.character(numOfTrees),"New/bdYes/","SeedYMean",as.character(seeds[4]),".pdf",sep = ''))
attach(mtcars)
par(mfrow=c(2,2))
custom_plotting2(accurateDataVoxBox4,meanBDYes4,paste("Seed", as.character(seeds[4])))


pdf(paste("Rad",as.character(numOfTrees),"New/bdYes/","SeedYMean",as.character(seeds[5]),".pdf",sep = ''))
attach(mtcars)
par(mfrow=c(2,2))
custom_plotting2(accurateDataVoxBox5,meanBDYes5,paste("Seed", as.character(seeds[5])))


