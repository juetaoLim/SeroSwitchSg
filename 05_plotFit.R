####################model fit#####################
####################model fit#####################
####################model fit#####################
rm(list=ls())
library(quantmod)
setwd("~/Google Drive/work_nea/sero_switch_sg/")
load("out/01outVAR.RData")
#create lag matrix, given number of lags.
xMat<- apply(DF,MARGIN=1,as.numeric)
xMat <- t(xMat)
xMatL_1 <- apply(xMat,MARGIN=2,FUN=Lag,k=1)
xMatL_2 <- apply(xMat,MARGIN=2,FUN=Lag,k=2)
xMatL_3 <- apply(xMat,MARGIN=2,FUN=Lag,k=3)
xMatL_4 <- apply(xMat,MARGIN=2,FUN=Lag,k=4)
xMat <- cbind(xMatL_1,xMatL_2,xMatL_3,xMatL_4); rm(xMatL_1,xMatL_2,xMatL_3,xMatL_4)
#not much difference between samples, just use var1 for model fit
beta <- apply(var1$BDraws,MARGIN=c(2,3),quantile,probs=0.5,na.rm=T)
yPred1 <- xMat%*%beta[,1]
yPred2 <- xMat%*%beta[,2]
yPred3 <- xMat%*%beta[,3]
yPred4 <- xMat%*%beta[,4]

#plot fit
plotterTS <- function(yT,yFit,fitCol,pointCol,ylabel,yaxt=F){
  limiter <- cbind(yT,yFit)
  plot(x=c(0,length(yT)),y=c(min(limiter,na.rm=T),max(limiter,na.rm=T))*1.2,col="white",xlab="",ylab=ylabel,xaxt='n')
  points(yT,pch=16,cex=0.5,col=pointCol)
  lines(yFit,lwd=2,col=fitCol)
  legend("topleft",legend=c("Total Cases","Fitted"),lty=c(NA,1),pch=c(16,NA),lwd=2,bty="n",col=c(pointCol,fitCol),cex=1.2)
  if(yaxt==T){axis(side=1,
                   at=seq(from=26,
                          to=length(yT)+26,
                          length.out=length(seq(from=2006,to=2021,by=1))),
                          labels=c(seq(from=2006,to=2020,by=1),""),
                          tick=F)}
  if(yaxt==T){axis(side=1,
                   at=seq(from=1,
                          to=length(yT),
                          length.out=length(seq(from=2006,to=2020,by=1))),
                   labels=rep("",length.out=length(seq(from=2006,to=2020,by=1))),
                   tick=T)}
  box()
}
pdf("plot/plotFit.pdf",width=8.27,height=7)
par(mfrow=c(4,1),las=1,mar=c(0,4,0,2),oma=c(4,2,4,2),cex.lab=1.2)
plotterTS(yT=DF[,1],yFit=yPred1,fitCol="darkorchid3",pointCol="tomato1",ylabel="Serotype 1")
plotterTS(yT=DF[,2],yFit=yPred2,fitCol="cadetblue4",pointCol="tomato1",ylabel="Serotype 2")
plotterTS(yT=DF[,3],yFit=yPred3,fitCol="thistle4",pointCol="tomato1",ylabel="Serotype 3")
plotterTS(yT=DF[,4],yFit=yPred4,fitCol="peru",pointCol="tomato1",ylabel="Serotype 4",yaxt=T)
dev.off()
