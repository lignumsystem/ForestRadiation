#!/usr/bin/env Rscript

source("read_files.r") # source the file to be read

relativeDeviation<-function(Qvox,Qacc){
# Function calculated the relative error between the radiations used by the 2 methods
         (Qvox-Qacc)/Qacc
        }

relativeHeight<- function(h){
# This funcion calculates  the relative height of the trees
  minHeight = min(h)
  maxHeight = max(h)
  (h - minHeight)/(maxHeight - minHeight)
}
 
yvalue<- function(a,b,box){
# a clever way to group heights and then to select correspond values of the height sets to plot the mean
    df<-factor(cut(b$V4[b$V1==0.1],10))
    y= 1
    a$df=df
    b$df=df
    for (z in levels(df)){ y<- c(y,mean((b$V7[b$df==z&b$V1==box]-a$V7[a$df==z]) /a$V7[a$df ==z]))}
    y=y[2:11] 

}

# calculatin the index in this code

index_summary<- function(Qvox,Qacc)
{
   N = 130
   sqrt(sum((Qvox-Qacc)^2)/N)
}

CalculateRelHtandDev<- function(meanDataVoxBox1,accurateDataVoxBox1,dirDataVoxBox1,titleString){ 

# calculates the relative deviation for mean STAR values
relDevM1 = relativeDeviation(meanDataVoxBox1$V7[meanDataVoxBox1$V1 == 0.1],accurateDataVoxBox1[7])
relHgtM1 = relativeHeight(meanDataVoxBox1$V4[meanDataVoxBox1$V1 == 0.1])

relDevM2 = relativeDeviation(meanDataVoxBox1$V7[meanDataVoxBox2$V1 == 0.2],accurateDataVoxBox1[7])
relHgtM2 = relativeHeight(meanDataVoxBox1$V4[meanDataVoxBox1$V1 == 0.2])

relDevM3 = relativeDeviation(meanDataVoxBox1$V7[meanDataVoxBox1$V1 == 0.3],accurateDataVoxBox1[7])
relHgtM3 = relativeHeight(meanDataVoxBox1$V4[meanDataVoxBox1$V1 == 0.3])

relDevM4 = relativeDeviation(meanDataVoxBox1$V7[meanDataVoxBox1$V1 == 0.4],accurateDataVoxBox1[7])
relHgtM4 = relativeHeight(meanDataVoxBox1$V4[meanDataVoxBox1$V1 == 0.4])


#For Directional STAR values
relDevD1 = relativeDeviation(dirDataVoxBox1$V7[dirDataVoxBox1$V1 == 0.1],accurateDataVoxBox1[7])
relHgtD1 = relativeHeight(dirDataVoxBox1$V4[dirDataVoxBox1$V1 == 0.1])

relDevD2 = relativeDeviation(dirDataVoxBox1$V7[dirDataVoxBox1$V1 == 0.2],accurateDataVoxBox1[7])
relHgtD2 = relativeHeight(dirDataVoxBox1$V4[dirDataVoxBox1$V1 == 0.2])

relDevD3 = relativeDeviation(dirDataVoxBox1$V7[dirDataVoxBox1$V1 == 0.3],accurateDataVoxBox1[7])
relHgtD3 = relativeHeight(dirDataVoxBox1$V4[dirDataVoxBox1$V1==0.3])

relDevD4 = relativeDeviation(dirDataVoxBox1$V7[dirDataVoxBox1$V1 == 0.4],accurateDataVoxBox1[7])
relHgtD4 = relativeHeight(dirDataVoxBox1$V4[dirDataVoxBox1$V1==0.4])




matplot(relHgtM1,relDevM1,type="p",pch = 1,col = 2:6,xlab = "Relative Height", ylab = "Relative Deviation for Mean" ,ylim=c(-0.2,0.2),main= "Vox =0.1")  
y<- yvalue(accurateDataVoxBox1,meanDataVoxBox1,0.1)
par(new=T)
plot(seq(from=0.0,to=1.0,by=1/9),y,type="l",ylim=c(-0.2,0.2),xlim=c(0,1),xlab = "Relative Height", ylab = "Relative Deviation for Mean")
abline(0,0) 

matplot(relHgtM1,relDevD1,type="p",pch = 1,col = 2:6,xlab = "Relative Height", ylab = "Relative Deviation Directional" ,ylim=c(-0.2,0.2),main= "Vox =0.1")  
 y<- yvalue(accurateDataVoxBox1,dirDataVoxBox1,0.1)
par(new=T)
plot(seq(from=0.0,to=1.0,by=1/9),y,type="l",ylim=c(-0.2,0.2),xlim=c(0,1),xlab = "Relative Height", ylab = "Relative Deviation Directional")

abline(0,0)

matplot(relHgtM2,relDevM2,type="p",pch = 1,col = 2:6,xlab = "Relative Height", ylab = "Relative Deviation for Mean" ,ylim=c(-0.2,0.2),main= "Vox =0.2") 
 y<- yvalue(accurateDataVoxBox1,meanDataVoxBox1,0.2)
par(new=T)
plot(seq(from=0.0,to=1.0,by=1/9),y,type="l",ylim=c(-0.2,0.2),xlim=c(0,1),xlab = "Relative Height", ylab = "Relative Deviation for Mean")

abline(0,0)

#directional for 0.2
matplot(relHgtM2,relDevD2,type="p",pch = 1,col = 2:6,xlab = "Relative Height", ylab = "Relative Deviation Directional" ,ylim=c(-0.2,0.2),main= "Vox =0.2") 

y<- yvalue(accurateDataVoxBox1,dirDataVoxBox1,0.2)
par(new=T)
plot(seq(from=0.0,to=1.0,by=1/9),y,type="l",ylim=c(-0.2,0.2),xlim=c(0,1),xlab = "Relative Height", ylab = "Relative Deviation Directional") 
abline(0,0)

#mean for 0.3

matplot(relHgtM3,relDevM3,type="p",pch = 1,col = 2:6,xlab = "Relative Height", ylab = "Relative Deviation for Mean" ,ylim=c(-0.2,0.2),main= "Vox =0.3") 
 y<- yvalue(accurateDataVoxBox1,meanDataVoxBox1,0.3)
par(new=T)
plot(seq(from=0.0,to=1.0,by=1/9),y,type="l",ylim=c(-0.2,0.2),xlim=c(0,1),xlab = "Relative Height", ylab = "Relative Deviation for Mean")
abline(0,0)

#directinoal for 0.3
matplot(relHgtM3,relDevD3,type="p",pch = 1,col = 2:6,xlab = "Relative Height", ylab = "Relative Deviation Directional" ,ylim=c(-0.2,0.2),main= "Vox =0.3") 
y<- yvalue(accurateDataVoxBox1,dirDataVoxBox1,0.3)
par(new=T)
plot(seq(from=0.0,to=1.0,by=1/9),y,type="l",ylim=c(-0.2,0.2),xlim=c(0,1),xlab = "Relative Height", ylab = "Relative Deviation Directional")  
abline(0,0)



#mean for 0.4

matplot(relHgtM4,relDevM4,type="p",pch = 1,col = 2:6,xlab = "Relative Height", ylab = "Relative Deviation for Mean" ,ylim=c(-0.2,0.2),main= "Vox =0.4") 
y<- yvalue(accurateDataVoxBox1,meanDataVoxBox1,0.4)
par(new=T)
plot(seq(from=0.0,to=1.0,by=1/9),y,type="l",ylim=c(-0.2,0.2),xlim=c(0,1),xlab = "Relative Height", ylab = "Relative Deviation for Mean")
abline(0,0)

#directional for 0.4
matplot(relHgtM4,relDevD4,type="p",pch = 1,col = 2:6,xlab = "Relative Height", ylab = "Relative Deviation Directional" ,ylim=c(-0.2,0.2),main= "Vox =0.4")
y<- yvalue(accurateDataVoxBox1,dirDataVoxBox1,0.4)
par(new=T)
plot(seq(from=0.0,to=1.0,by=1/9),y,type="l",ylim=c(-0.2,0.2),xlim=c(0,1),xlab = "Relative Height", ylab = "Relative Deviation Directional")    
abline(0,0)
#title(main="Relative height  vs Relative Deviation",outer=T)
title(main= titleString,outer=T)
dev.off()

}


CalculateRelHtandDev2<- function(meanDataVoxBox1,accurateDataVoxBox1,titleString){

relDevBDM1 = relativeDeviation(meanDataVoxBox1$V7[meanDataVoxBox1$V1 == 0.1],accurateDataVoxBox1[7])
relHgtBDM1 = relativeHeight(meanDataVoxBox1$V4[meanDataVoxBox1$V1 == 0.1])

relDevBDM2 = relativeDeviation(meanDataVoxBox1$V7[meanDataVoxBox1$V1 == 0.2],accurateDataVoxBox1[7])
relHgtBDM2 = relativeHeight(meanDataVoxBox1$V4[meanDataVoxBox1$V1 == 0.2])

relDevBDM3 = relativeDeviation(meanDataVoxBox1$V7[meanDataVoxBox1$V1 == 0.3],accurateDataVoxBox1[7])
relHgtBDM3 = relativeHeight(meanDataVoxBox1$V4[meanDataVoxBox1$V1 == 0.3])

relDevBDM4 = relativeDeviation(meanDataVoxBox1$V7[meanDataVoxBox1$V1 == 0.4],accurateDataVoxBox1[7])
relHgtBDM4 = relativeHeight(meanDataVoxBox1$V4[meanDataVoxBox1$V1 == 0.4])


matplot(relHgtBDM1,relDevBDM1,type="p",pch = 1,col = 2:6,xlab = "Relative Height", ylab = "Relative Deviation for Mean BDY" ,ylim=c(-0.3,0.3),main= "Vox =0.1")  
y<- yvalue(accurateDataVoxBox1,meanDataVoxBox1,0.1)
par(new=T)
plot(seq(from=0.0,to=1.0,by=1/9),y,type="l",ylim=c(-0.3,0.3),xlim=c(0,1),xlab = "Relative Height", ylab = "Relative Deviation for Mean BDY")
abline(0,0) 

matplot(relHgtBDM2,relDevBDM2,type="p",pch = 1,col = 2:6,xlab = "Relative Height", ylab = "Relative Deviation for Mean BDY" ,ylim=c(-0.3,0.3),main= "Vox =0.2")  
y<- yvalue(accurateDataVoxBox1,meanDataVoxBox1,0.2)
par(new=T)
plot(seq(from=0.0,to=1.0,by=1/9),y,type="l",ylim=c(-0.3,0.3),xlim=c(0,1),xlab = "Relative Height", ylab = "Relative Deviation for Mean BDY")
abline(0,0) 

matplot(relHgtBDM3,relDevBDM3,type="p",pch = 1,col = 2:6,xlab = "Relative Height", ylab = "Relative Deviation for Mean BDY" ,ylim=c(-0.3,0.3),main= "Vox =0.3")  
y<- yvalue(accurateDataVoxBox1,meanDataVoxBox1,0.3)
par(new=T)
plot(seq(from=0.0,to=1.0,by=1/9),y,type="l",ylim=c(-0.3,0.3),xlim=c(0,1),xlab = "Relative Height", ylab = "Relative Deviation for Mean BDY")
abline(0,0) 

matplot(relHgtBDM4,relDevBDM4,type="p",pch = 1,col = 2:6,xlab = "Relative Height", ylab = "Relative Deviation for Mean BDY" ,ylim=c(-0.3,0.3),main= "Vox =0.4")  
y<- yvalue(accurateDataVoxBox1,meanDataVoxBox1,0.4)
par(new=T)
plot(seq(from=0.0,to=1.0,by=1/9),y,type="l",ylim=c(-0.3,0.3),xlim=c(0,1),xlab = "Relative Height", ylab = "Relative Deviation for Mean BDY")
abline(0,0) 
title(main= titleString,outer=T)
dev.off()
}

pdf("Deviation/deviationSeed26497.pdf")
attach(mtcars)
par(mfrow=c(4,2))
CalculateRelHtandDev(meanDataVoxBox1,accurateDataVoxBox1,dirDataVoxBox1,"SEED 26497")

pdf("Deviation/meanyes/deviationSeed26497.pdf")
attach(mtcars)
par(mfrow=c(4,2))
CalculateRelHtandDev2(meanBDYes1,accurateDataVoxBox1,"SEED 26497")


pdf("Deviation/deviationSeed26557.pdf")
attach(mtcars)
par(mfrow=c(4,2))
CalculateRelHtandDev(meanDataVoxBox2,accurateDataVoxBox2,dirDataVoxBox2,"Seed 26557")

pdf("Deviation/meanyes/deviationSeed26557.pdf")
attach(mtcars)
par(mfrow=c(4,2))
CalculateRelHtandDev2(meanBDYes2,accurateDataVoxBox2,"Seed 26557")

pdf("Deviation/deviationSeed35577.pdf")
attach(mtcars)
par(mfrow=c(4,2))
CalculateRelHtandDev(meanDataVoxBox3,accurateDataVoxBox3,dirDataVoxBox3,"Seed 35577")

pdf("Deviation/meanyes/deviationSeed35577.pdf")
attach(mtcars)
par(mfrow=c(4,2))
CalculateRelHtandDev2(meanBDYes3,accurateDataVoxBox3,"Seed 35577")


pdf("Deviation/deviationSeed36577.pdf")
attach(mtcars)
par(mfrow=c(4,2))
CalculateRelHtandDev(meanDataVoxBox4,accurateDataVoxBox4,dirDataVoxBox4,"Seed 36577")

pdf("Deviation/meanyes/deviationSeed36577.pdf")
attach(mtcars)
par(mfrow=c(4,2))
CalculateRelHtandDev2(meanBDYes4,accurateDataVoxBox4,"Seed 36577")

pdf("Deviation/deviationSeed39577.pdf")
attach(mtcars)
par(mfrow=c(4,2))
CalculateRelHtandDev(meanDataVoxBox5,accurateDataVoxBox5,dirDataVoxBox5,"Seed 39577")

pdf("Deviation/meanyes/deviationSeed39577.pdf")
attach(mtcars)
par(mfrow=c(4,2))
CalculateRelHtandDev2(meanBDYes5,accurateDataVoxBox5,"Seed 39577")

