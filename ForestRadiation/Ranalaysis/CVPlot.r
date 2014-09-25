#!/usr/bin/env Rscript
args <- commandArgs()
commandArgs <- function() args[1]
source("read_files.r")
source("functionsInR.r") #Functions are obtained here
voxBoxNo = c(0.1,0.2,0.3,0.4)
path1 = paste(rootpath,"/Deviation/","CV_",as.character(args[1]),sep='')
dir.create(path1, showWarnings = TRUE, recursive = FALSE, mode = "0777")


pathmain = paste(rootpath,"/Deviation/","CV_",as.character(args[1]),"/",as.character(numOfTrees),sep='')
dir.create(pathmain, showWarnings = TRUE, recursive = FALSE, mode = "0777")


cnt = 1
for(seedInt in seeds){
   meanData     = cbind (get(paste("meanDataVoxBox",as.character(cnt),sep='')))
   accurateData = cbind(get(paste("accurateDataVoxBox",as.character(cnt),sep='')))
   cnt = cnt +1
}

plotCVCoeff<-function(data,vox,htlev,titleString){

   ty = list()
   for (i in 1:length(htlev)){
    realData = data$V7[data$V4==htlev[i]&data$V1==vox]
    ty[i]  = calculateCV(realData)
  }

  relHgtM1 = relativeHeight(data$V4[data$V1 == vox])
  newRelHt = levels(as.factor(relHgtM1)) 

  plot(newRelHt,ty,type="l",xlab = "Relative Height", ylab = "Relative Coefficient of Variation")
  title(main= titleString,outer=T)
 }

# To get the data for the levels and the classification/enumerator type.

for (box in voxBoxNo){
pdf(paste(rootpath,"/Deviation/CV_",as.character(args[1]),"/",as.character(numOfTrees),"/","CVversesRelativeHt","Vox",as.character(box),".","pdf",sep = '' ))
attach(mtcars)
htlev = levels(as.factor(meanData$V4))
plotCVCoeff(meanData,box,htlev,paste("All Seeds Vox =",as.character(box),sep =''))
}
  

pdf(paste(rootpath,"/Deviation/CV_",as.character(args[1]),"/",as.character(numOfTrees),"/","CVversesRelativeHtAccurate","Vox",as.character(0.2),".","pdf",sep = '' ))
attach(mtcars)
htlev = levels(as.factor(accurateData$V4))
 plotCVCoeff(accurateData,0.2,htlev,paste("All Seeds Vox =",as.character(0.2),sep =''))



