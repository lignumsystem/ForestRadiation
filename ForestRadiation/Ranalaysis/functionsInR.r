#!/usr/bin/env Rscript
# All the functions required in analysis are defined here.

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
    df<-factor(cut(b$V4[b$V1==box],10))
    y= 1
    a$df=df
    b$df=df
    for (z in levels(df)){y<- c(y,mean((b$V7[b$df==z&b$V1==box]-a$V7[a$df==z]) /a$V7[a$df ==z]))}
    y=y[2:11] 

}

# calculating the index in this code
index_summary<- function(Qvox,Qacc)
{
   N = 130
   sqrt(sum((Qvox-Qacc)^2)/N)
}

#specify the decimal number
specify_decimal <- function(x, k) format(round(x, k), nsmall=k)


plotTreeArrangement<- function(data,titleString)
{
plot(data$x,data$y, xlim= c(0,20),ylim=c(0,20),col = "red",xlab = "x", ylab = "y" )
par(new=T)
plot(10,10,type ="p",pch = 17,col = "black",xlab = "x", ylab = "y",xlim= c(0,20),ylim=c(0,20))
title(main= titleString,outer=T)
}
