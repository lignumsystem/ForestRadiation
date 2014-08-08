#!/usr/bin/env Rscript
seeds = c(787237, 536199,  1676779,  546327,  235663) # These are the seeds to be changed
numOfTrees =  100 # this is the no of trees you want to plot the analysis on.
# Add all the files that needs to be read in to this script. This will enable the user to read all the files needed in one script.
# Change the paths according to your computer. read files that are in this folder you get them running the .sh files. Please change the paths where this data is and usr them here
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

pathDir      = paste("/home/likewise-open/IN/gopalkri/Developer/core-model/ForestRadiation/Resultapp/",as.character(numOfTrees),"/Directional/",sep = '')
pathBoxDirNo = paste("/home/likewise-open/IN/gopalkri/Developer/core-model/ForestRadiation/Resultapp/",as.character(numOfTrees),"/MeanStar/BoxDirNo/",sep ='')
pathAccurate = paste("/home/likewise-open/IN/gopalkri/Developer/core-model/ForestRadiation/Resultapp/accurate",as.character(numOfTrees),"/",sep='')
pathBoxDirYes= paste("/home/likewise-open/IN/gopalkri/Developer/core-model/ForestRadiation/Resultapp/",as.character(numOfTrees),"/MeanStar/BoxDirYes/",sep='')
pathdataSeed = paste("/home/likewise-open/IN/gopalkri/Developer/core-model/ForestRadiation/Resultapp/",as.character(numOfTrees),"/",sep='')

dirDataVoxBox1<- read.table(file.path(pathDir,paste("boxdirNdir","_",as.character(seeds[1]),"_","den",as.character(numOfTrees),".dat",sep = '')))

dirDataVoxBox2<- read.table(file.path(pathDir,paste("boxdirNdir","_",as.character(seeds[2]),"_","den",as.character(numOfTrees),".dat",sep = '')))

dirDataVoxBox3<- read.table(file.path(pathDir,paste("boxdirNdir","_",as.character(seeds[3]),"_","den",as.character(numOfTrees),".dat",sep = '')))

dirDataVoxBox4<- read.table(file.path(pathDir,paste("boxdirNdir","_",as.character(seeds[4]),"_","den",as.character(numOfTrees),".dat",sep = '')))

dirDataVoxBox5<- read.table(file.path(pathDir,paste("boxdirNdir","_",as.character(seeds[5]),"_","den",as.character(numOfTrees),".dat",sep = '')))

#*********************************************************************************************************************************
meanDataVoxBox1<- read.table(file.path(pathBoxDirNo,paste("boxdirNmean","_",as.character(seeds[1]),"_","den",as.character(numOfTrees),".dat",sep = '')))

meanDataVoxBox2<- read.table(file.path(pathBoxDirNo,paste("boxdirNmean","_",as.character(seeds[2]),"_","den",as.character(numOfTrees),".dat",sep = '')))

meanDataVoxBox3<- read.table(file.path(pathBoxDirNo,paste("boxdirNmean","_",as.character(seeds[3]),"_","den",as.character(numOfTrees),".dat",sep = '')))

meanDataVoxBox4<- read.table(file.path(pathBoxDirNo,paste("boxdirNmean","_",as.character(seeds[4]),"_","den",as.character(numOfTrees),".dat",sep = '')))

meanDataVoxBox5<- read.table(file.path(pathBoxDirNo,paste("boxdirNmean","_",as.character(seeds[5]),"_","den",as.character(numOfTrees),".dat",sep = '')))


#*************************************************************************************************************************************

accurateDataVoxBox1<- read.table(file.path(pathAccurate,paste("accurate","_",as.character(seeds[1]),"_","den",as.character(numOfTrees),".dat",sep = '')))

accurateDataVoxBox2<- read.table(file.path(pathAccurate,paste("accurate","_",as.character(seeds[2]),"_","den",as.character(numOfTrees),".dat",sep = '')))

accurateDataVoxBox3<- read.table(file.path(pathAccurate,paste("accurate","_",as.character(seeds[3]),"_","den",as.character(numOfTrees),".dat",sep = '')))

accurateDataVoxBox4<- read.table(file.path(pathAccurate,paste("accurate","_",as.character(seeds[4]),"_","den",as.character(numOfTrees),".dat",sep = '')))

accurateDataVoxBox5<- read.table(file.path(pathAccurate,paste("accurate","_",as.character(seeds[5]),"_","den",as.character(numOfTrees),".dat",sep = '')))
#************************************************************************************************************************************


meanBDYes1<- read.table(file.path(pathBoxDirYes,paste("boxdirYmean","_",as.character(seeds[1]),"_","den",as.character(numOfTrees),".dat",sep = '')))

meanBDYes2<- read.table(file.path(pathBoxDirYes,paste("boxdirYmean","_",as.character(seeds[2]),"_","den",as.character(numOfTrees),".dat",sep = '')))

meanBDYes3<- read.table(file.path(pathBoxDirYes,paste("boxdirYmean","_",as.character(seeds[3]),"_","den",as.character(numOfTrees),".dat",sep = '')))

meanBDYes4<- read.table(file.path(pathBoxDirYes,paste("boxdirYmean","_",as.character(seeds[4]),"_","den",as.character(numOfTrees),".dat",sep = '')))

meanBDYes5<- read.table(file.path(pathBoxDirYes,paste("boxdirYmean","_",as.character(seeds[5]),"_","den",as.character(numOfTrees),".dat",sep = '')))
#***********************************************************************************************************************************
dataSeed1 = read.table(file.path(pathdataSeed,paste("runfile","_",as.character(seeds[1]),"_",as.character(numOfTrees),"_","BDN","_","0.1",".","dat",sep = '')),header = TRUE)

dataSeed2 = read.table(file.path(pathdataSeed,paste("runfile","_",as.character(seeds[2]),"_",as.character(numOfTrees),"_","BDN","_","0.1",".","dat",sep = '')),header = TRUE)


dataSeed3 = read.table(file.path(pathdataSeed,paste("runfile","_",as.character(seeds[3]),"_",as.character(numOfTrees),"_","BDN","_","0.1",".","dat",sep = '')),header = TRUE)

dataSeed4 = read.table(file.path(pathdataSeed,paste("runfile","_",as.character(seeds[4]),"_",as.character(numOfTrees),"_","BDN","_","0.1",".","dat",sep = '')),header = TRUE)


dataSeed5 = read.table(file.path(pathdataSeed,paste("runfile","_",as.character(seeds[5]),"_",as.character(numOfTrees),"_","BDN","_","0.1",".","dat",sep = '')),header = TRUE)



