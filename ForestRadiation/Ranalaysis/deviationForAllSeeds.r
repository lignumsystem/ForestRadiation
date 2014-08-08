#!/usr/bin/env Rscript
source("read_files.r") # Read in the files of data
seeds =c(787237, 536199,  1676779,  546327,  235663) # Specify the seeds here
#numOfTrees =  100
source("functionsInR.r") #Functions are obtained here
pathmain = paste("Deviation/","AllSeeds",as.character(numOfTrees),sep='')
dir.create(pathmain, showWarnings = TRUE, recursive = FALSE, mode = "0777")

dirData      = cbind(dirDataVoxBox1, dirDataVoxBox2,dirDataVoxBox3,dirDataVoxBox4,dirDataVoxBox5)
meanData     = cbind(meanDataVoxBox1, meanDataVoxBox2,meanDataVoxBox3,meanDataVoxBox4,meanDataVoxBox5)
accurateData = cbind(accurateDataVoxBox1, accurateDataVoxBox2,accurateDataVoxBox3,accurateDataVoxBox4,accurateDataVoxBox5)
meanBDYesData= cbind(meanBDYes1, meanBDYes2,meanBDYes3,meanBDYes4,meanBDYes5)


CalculateRelHtandDev<- function(meanDataVoxBox,accurateDataVoxBox,dirDataVoxBox,meanBDYes,titleString,vox){ 
# calculates the relative deviation for mean STAR values
relDevM1 = relativeDeviation(meanDataVoxBox$V7[meanDataVoxBox$V1 == vox],accurateDataVoxBox[7])
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
#*******************************************************************************************************************************

pdf(paste("Deviation/AllSeeds",as.character(numOfTrees),"/","deviationAllSeed","Vox0.1.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanData,accurateData,dirData,meanBDYesData,"All Seeds Vox= 0.1",0.1)

pdf(paste("Deviation/AllSeeds", as.character(numOfTrees),"/","deviationAllSeed","Vox0.2.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanData,accurateData,dirData,meanBDYesData,"All Seeds Vox= 0.2",0.2)


pdf(paste("Deviation/AllSeeds",as.character(numOfTrees),"/","deviationAllSeed","Vox0.3.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanData,accurateData,dirData,meanBDYesData,"All Seeds Vox= 0.3",0.3)


pdf(paste("Deviation/AllSeeds",as.character(numOfTrees),"/","deviationAllSeed","Vox0.4.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanData,accurateData,dirData,meanBDYesData,"All Seeds Vox= 0.4",0.4)




