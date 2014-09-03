#!/usr/bin/env Rscript
args <- commandArgs()
batch_args <- args[1]
print(args[1])
numOfTrees = 150
seeds = c(787237, 536199 , 1676779,  546327,  235663 ,787231, 536197, 1576775, 446327, 135663)

#Paths for data
pathDir = paste("/home/likewise-open/IN/gopalkri/Developer/core-model/ForestRadiation/Resultapp/",as.character(numOfTrees),"/numParts_",as.character(args[1]), "/Directional",sep = '')

pathBoxDirNo = paste("/home/likewise-open/IN/gopalkri/Developer/core-model/ForestRadiation/Resultapp/",as.character(numOfTrees),"/numParts_",as.character(args[1]),"/MeanStar/BoxDirNo/",sep ='')

pathAccurate = paste("/home/likewise-open/IN/gopalkri/Developer/core-model/ForestRadiation/Resultapp/accurate",as.character(numOfTrees),"/numparts_",as.character(args[1]),sep='')

pathBoxDirYes= paste("/home/likewise-open/IN/gopalkri/Developer/core-model/ForestRadiation/Resultapp/",as.character(numOfTrees),"/numParts_",as.character(args[1]), "/MeanStar/BoxDirYes/",sep='')

pathdataSeed = paste("/home/likewise-open/IN/gopalkri/Developer/core-model/ForestRadiation/Resultapp/",as.character(numOfTrees),"/",sep='')

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
