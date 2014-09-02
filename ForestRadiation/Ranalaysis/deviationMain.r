#!/usr/bin/env Rscript
#commandArgs <- function() 5
args <- commandArgs()
commandArgs <- function() args[1]
source("read_files.r")
path = "/home/likewise-open/IN/gopalkri/Developer/core-model/ForestRadiation/Resultapp/Deviation"  # Change the path where you want the result

dir.create(path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
pathmain = paste("Deviation/",as.character(numOfTrees),sep='')
dir.create(pathmain, showWarnings = TRUE, recursive = FALSE, mode = "0777")
source("functionsInR.r")

voxBoxNo = c(0.1,0.2,0.3,0.4)

CalculateRelHtandDev<- function(meanDataVoxBox,accurateDataVoxBox,dirDataVoxBox,meanBDYes,titleString,vox){ 


relDevM1 = relativeDeviation(meanDataVoxBox$V7[meanDataVoxBox$V1 == vox],accurateDataVoxBox[7])
print(nrow(relDevM1))             
print(nrow(meanDataVoxBox$V4[meanDataVoxBox$V1 == vox]))
relHgtM1 = relativeHeight(meanDataVoxBox$V4[meanDataVoxBox$V1 == vox])
total_num1= index_summary(meanDataVoxBox$V7[meanDataVoxBox$V1 == vox],accurateDataVoxBox[7]) 
total_num1=specify_decimal(total_num1,3)

#For Directional STAR values
relDevD1 = relativeDeviation(dirDataVoxBox$V7[dirDataVoxBox$V1 == vox],accurateDataVoxBox[7])
relHgtD1 = relativeHeight(dirDataVoxBox$V4[dirDataVoxBox$V1 == vox])
total_num2= index_summary(dirDataVoxBox$V7[dirDataVoxBox$V1 == vox],accurateDataVoxBox[7]) 
total_num2=specify_decimal(total_num2,3)

relDevBDM1 = relativeDeviation(meanBDYes$V7[meanBDYes$V1 == vox],accurateDataVoxBox[7])
relHgtBDM1 = relativeHeight(meanBDYes$V4[meanBDYes$V1 == vox])
total_num3= index_summary(meanBDYes$V7[meanBDYes$V1 == vox],accurateDataVoxBox[7]) 
total_num3=specify_decimal(total_num3,3)



matplot(relHgtM1,relDevM1,type="p",pch = 1,col = 2:6,xlab = "Relative Height", ylab = "Relative Deviation for Mean" ,ylim=c(-0.5,0.5),main= paste("Vox =",as.character(vox),"RMSE =",as.character(total_num1))) 
y<- yvalue(accurateDataVoxBox,meanDataVoxBox,vox)
par(new=T)
plot(seq(from=0.0,to=1.0,by=1/9),y,type="l",ylim=c(-0.5,0.5),xlim=c(0,1),xlab = "Relative Height", ylab = "Relative Deviation for Mean")
abline(0,0) 

matplot(relHgtD1,relDevD1,type="p",pch = 1,col = 2:6,xlab = "Relative Height", ylab = "Relative Deviation Directional" ,ylim=c(-0.5,0.5),main= paste("Vox =",as.character(vox),"RMSE =",as.character(total_num2))) 
 yd<- yvalue(accurateDataVoxBox,dirDataVoxBox,vox)
par(new=T)
plot(seq(from=0.0,to=1.0,by=1/9),yd,type="l",ylim=c(-0.5,0.5),xlim=c(0,1),xlab = "Relative Height", ylab = "Relative Deviation Directional")
abline(0,0)

matplot(relHgtBDM1,relDevBDM1,type="p",pch = 1,col = 2:6,xlab = "Relative Height", ylab = "Relative Deviation for Mean BDY" ,ylim=c(-0.5,0.5),main= paste("Vox =",as.character(vox),"RMSE =",as.character(total_num3)))  
yy<- yvalue(accurateDataVoxBox,meanDataVoxBox,vox)
par(new=T)
plot(seq(from=0.0,to=1.0,by=1/9),yy,type="l",ylim=c(-0.5,0.5),xlim=c(0,1),xlab = "Relative Height", ylab = "Relative Deviation for Mean BDY")
abline(0,0) 

title(main= titleString,outer=T)
dev.off()
}



#*********************************************************************
                                   

path = paste("Deviation/",as.character(numOfTrees),"/","numParts_",as.character(args[1]),sep='')
dir.create(path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
cnt = 1

for(seedInt in seeds){
path = paste("Deviation/",as.character(numOfTrees),"/","numParts_",as.character(args[1]),"/",as.character(seedInt),sep='')
dir.create(path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
for (box in voxBoxNo){

pdf(paste("Deviation/",as.character(numOfTrees),"/","numParts_",as.character(args[1]),"/",as.character(seedInt),"/","deviationSeed",as.character(seedInt),"Vox",as.character(box),".","pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(get(paste("meanDataVoxBox",as.character(cnt),sep='')),get(paste("accurateDataVoxBox",as.character(cnt),sep='')),get(paste("dirDataVoxBox",as.character(cnt),sep='')),get(paste("meanBDYes",as.character(cnt),sep='')),paste("SEED",as.character(seedInt)),box)

}
cnt = cnt+1
}
