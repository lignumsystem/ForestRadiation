#!/usr/bin/env Rscript
args <- commandArgs()
batch_args <- args[1]
print(args[1])
numOfTrees = 22
seeds = c(787237, 536199 , 1676779,  546327,  235663 ,787231, 536197, 1576775, 446327, 135663)

#Paths for data
rootpath="/home/likewise-open/IN/gopalkri/Developer/core-model/ForestRadiation/Resultapp"
pathDir = paste(rootpath,"/",as.character(numOfTrees),"/numParts_",as.character(args[1]), "/Directional",sep = '')

pathBoxDirNo = paste(rootpath,"/",as.character(numOfTrees),"/numParts_",as.character(args[1]),"/MeanStar/BoxDirNo/",sep ='')

pathAccurate = paste(rootpath,"/","accurate",as.character(numOfTrees),"/numparts_",as.character(args[1]),sep='')

pathBoxDirYes= paste(rootpath,"/",as.character(numOfTrees),"/numParts_",as.character(args[1]), "/MeanStar/BoxDirYes/",sep='')

pathdataSeed = paste(rootpath,"/",as.character(numOfTrees),"/",sep='')

print(pathBoxDirYes)
cnt = 1

for(seedInt in seeds)
{
 dirDataVB<- read.table(file.path(pathDir,paste("boxdirNdir","_",as.character(seedInt),"_","den",as.character(numOfTrees),".dat",sep = '')))
 assign(paste("dirDataVoxBox", cnt, sep = ''),dirDataVB)

#print(dirDataVoxBox1)
 
meanDataVB <- read.table(file.path(pathBoxDirNo,paste("boxdirNmean","_",as.character(seedInt),"_","den",as.character(numOfTrees),".dat",sep = '')))
 assign(paste("meanDataVoxBox", cnt, sep = ""),meanDataVB)

accurateDataVB<-  read.table(file.path(pathAccurate,paste("accurate","_",as.character(seedInt),"_","den",as.character(numOfTrees),".dat",sep = '')))
 assign(paste("accurateDataVoxBox", cnt, sep = ""),accurateDataVB)

meanBDVB<- read.table(file.path(pathBoxDirYes,paste("boxdirYmean","_",as.character(seedInt),"_","den",as.character(numOfTrees),".dat",sep = '')))
 assign(paste("meanBDYes", cnt, sep = ""),meanBDVB)
cnt = cnt +1
   
}
