#!/usr/bin/env Rscript
seeds =c(787237, 536199,  1676779,  546327,  235663)
source("read_files.r") # source the file to be read
path = "/home/likewise-open/IN/gopalkri/Developer/core-model/ForestRadiation/Resultapp/Deviation"  # Change the path where you want the results here.
dir.create(path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
pathmain = paste("Deviation/",as.character(numOfTrees),sep='')
dir.create(pathmain, showWarnings = TRUE, recursive = FALSE, mode = "0777")
source("functionsInR.r")
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

# FOR SEED1
path = paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[1]),sep='')
dir.create(path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[1]),"/","deviationSeed",as.character(seeds[1]),"Vox0.1.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanDataVoxBox1,accurateDataVoxBox1,dirDataVoxBox1,meanBDYes1,paste("SEED",as.character(seeds[1])),0.1)

pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[1]),"/","deviationSeed",as.character(seeds[1]),"Vox0.2.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanDataVoxBox1,accurateDataVoxBox1,dirDataVoxBox1,meanBDYes1,paste("SEED",as.character(seeds[1])),0.2)


pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[1]),"/","deviationSeed",as.character(seeds[1]),"Vox0.3.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanDataVoxBox1,accurateDataVoxBox1,dirDataVoxBox1,meanBDYes1,paste("SEED",as.character(seeds[1])),0.3)


pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[1]),"/","deviationSeed",as.character(seeds[1]),"Vox0.4.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanDataVoxBox1,accurateDataVoxBox1,dirDataVoxBox1,meanBDYes1,paste("SEED",as.character(seeds[1])),0.4)




#FOR SEED 2
path = paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[2]),sep='')
dir.create(path, showWarnings = TRUE, recursive = FALSE, mode = "0777")

pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[2]),"/","deviationSeed",as.character(seeds[2]),"Vox0.1.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanDataVoxBox2,accurateDataVoxBox2,dirDataVoxBox2,meanBDYes2,paste("SEED",as.character(seeds[2])),0.1)

pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[2]),"/","deviationSeed",as.character(seeds[2]),"Vox0.2.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanDataVoxBox2,accurateDataVoxBox2,dirDataVoxBox2,meanBDYes2,paste("SEED",as.character(seeds[2])),0.2)


pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[2]),"/","deviationSeed",as.character(seeds[2]),"Vox0.3.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanDataVoxBox2,accurateDataVoxBox2,dirDataVoxBox2,meanBDYes2,paste("SEED",as.character(seeds[2])),0.3)


pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[2]),"/","deviationSeed",as.character(seeds[2]),"Vox0.4.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanDataVoxBox2,accurateDataVoxBox2,dirDataVoxBox2,meanBDYes2,paste("SEED",as.character(seeds[2])),0.4)


#SEED 3
path = paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[3]),sep='')
dir.create(path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[3]),"/","deviationSeed",as.character(seeds[3]),"Vox0.1.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanDataVoxBox3,accurateDataVoxBox3,dirDataVoxBox3,meanBDYes3,paste("SEED",as.character(seeds[3])),0.1)

pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[3]),"/","deviationSeed",as.character(seeds[3]),"Vox0.2.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanDataVoxBox3,accurateDataVoxBox3,dirDataVoxBox3,meanBDYes3,paste("SEED",as.character(seeds[3])),0.2)


pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[3]),"/","deviationSeed",as.character(seeds[3]),"Vox0.3.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanDataVoxBox3,accurateDataVoxBox3,dirDataVoxBox3,meanBDYes3,paste("SEED",as.character(seeds[3])),0.3)


pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[3]),"/","deviationSeed",as.character(seeds[3]),"Vox0.4.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanDataVoxBox3,accurateDataVoxBox3,dirDataVoxBox3,meanBDYes3,paste("SEED",as.character(seeds[3])),0.4)

#FOR SEED 4
path = paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[4]),sep='')
dir.create(path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[4]),"/","deviationSeed",as.character(seeds[4]),"Vox0.1.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanDataVoxBox4,accurateDataVoxBox4,dirDataVoxBox4,meanBDYes4,paste("SEED",as.character(seeds[4])),0.1)

pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[4]),"/","deviationSeed",as.character(seeds[4]),"Vox0.2.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanDataVoxBox4,accurateDataVoxBox4,dirDataVoxBox4,meanBDYes4,paste("SEED",as.character(seeds[4])),0.2)


pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[4]),"/","deviationSeed",as.character(seeds[4]),"Vox0.3.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanDataVoxBox4,accurateDataVoxBox4,dirDataVoxBox4,meanBDYes4,paste("SEED",as.character(seeds[4])),0.3)


pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[4]),"/","deviationSeed",as.character(seeds[4]),"Vox0.4.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanDataVoxBox4,accurateDataVoxBox4,dirDataVoxBox4,meanBDYes4,paste("SEED",as.character(seeds[4])),0.4)

#FOR SEED 5
path =  paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[5]),sep='')
dir.create(path,showWarnings = TRUE, recursive = FALSE, mode = "0777")
pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[5]),"/","deviationSeed",as.character(seeds[5]),"Vox0.1.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanDataVoxBox5,accurateDataVoxBox5,dirDataVoxBox5,meanBDYes5,paste("SEED",as.character(seeds[5])),0.1)


pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[5]),"/","deviationSeed",as.character(seeds[5]),"Vox0.2.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanDataVoxBox5,accurateDataVoxBox5,dirDataVoxBox5,meanBDYes5,paste("SEED",as.character(seeds[5])),0.2)



pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[5]),"/","deviationSeed",as.character(seeds[5]),"Vox0.3.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanDataVoxBox5,accurateDataVoxBox5,dirDataVoxBox5,meanBDYes5,paste("SEED",as.character(seeds[5])),0.3)


pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[5]),"/","deviationSeed",as.character(seeds[5]),"Vox0.4.pdf",sep = '' ))
attach(mtcars)
par(mfrow=c(2,2))
CalculateRelHtandDev(meanDataVoxBox5,accurateDataVoxBox5,dirDataVoxBox5,meanBDYes5,paste("SEED",as.character(seeds[5])),0.4)
#*********************************************************************************************************************************

pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[1]),"/","TreeArrangement",as.character(seeds[1]),".","pdf",sep = '' ))
plotTreeArrangement(dataSeed1,paste("SEED",as.character(seeds[1])))

pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[2]),"/","TreeArrangement",as.character(seeds[2]),".","pdf",sep = '' ))
plotTreeArrangement(dataSeed2,paste("SEED",as.character(seeds[2])))

pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[3]),"/","TreeArrangement",as.character(seeds[3]),".","pdf",sep = '' ))
plotTreeArrangement(dataSeed3,paste("SEED",as.character(seeds[3])))

pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[4]),"/","TreeArrangement",as.character(seeds[4]),".","pdf",sep = '' ))
plotTreeArrangement(dataSeed4,paste("SEED",as.character(seeds[4])))

pdf(paste("Deviation/",as.character(numOfTrees),"/",as.character(seeds[5]),"/","TreeArrangement",as.character(seeds[5]),".","pdf",sep = '' ))
plotTreeArrangement(dataSeed5,paste("SEED",as.character(seeds[5])))
#*********************************************************************************************************************************
