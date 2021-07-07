rm(list=ls())
library(VARsignR)
library(mvtnorm)
library(xtable)
setwd("~/Google Drive/work_nea/sero_switch_sg/")
data <- read.csv("data/data_sg_noms2.csv")
#get sero%

S1 <- data$D1
S2 <- data$D2
S3 <- data$D3
S4 <- data$D4
# # #first difference log
# S1 <- diff(data$S1)
# S2 <- diff(data$S2)
# S3 <- diff(data$S3)
# S4 <- diff(data$S4)
# 
# #first difference data
# diffDF <- cbind(S1,S2,S3,S4)
# diffDF <- as.ts(diffDF)

DF     <- cbind(S1,S2,S3,S4)
DF     <- as.ts(DF)

#setConstraints
constr1 <- c(+1,-2,-3,-4)
constr2 <- c(+2,-1,-3,-4)
constr3 <- c(+3,-2,-1,-4)
constr4 <- c(+4,-2,-3,-1)

#runMCMC using uhlig.reject
var_lag1 <- uhlig.reject(Y=DF,nlags=1,constrained = constr1,steps=300,draws=10000)
var_lag2 <- uhlig.reject(Y=DF,nlags=2,constrained = constr1,steps=300,draws=10000)
var_lag3 <- uhlig.reject(Y=DF,nlags=3,constrained = constr1,steps=300,draws=10000)
var_lag4 <- uhlig.reject(Y=DF,nlags=4,constrained = constr1,steps=300,draws=10000)


#create dataframe
require(quantmod)
xMat<- apply(DF,MARGIN=1,as.numeric)
xMat <- t(xMat)
xMatL_1 <- apply(xMat,MARGIN=2,FUN=Lag,k=1)
xMatL_2 <- apply(xMat,MARGIN=2,FUN=Lag,k=2)
xMatL_3 <- apply(xMat,MARGIN=2,FUN=Lag,k=3)
xMatL_4 <- apply(xMat,MARGIN=2,FUN=Lag,k=4)
xMat1 <- cbind(xMatL_1,1); 
xMat2 <- cbind(xMatL_1,xMatL_2,1); 
xMat3 <- cbind(xMatL_1,xMatL_2,xMatL_3,1); 
xMat4 <- cbind(xMatL_1,xMatL_2,xMatL_3,xMatL_4,1); 
rm(xMatL_1,xMatL_2,xMatL_3,xMatL_4)
#not much difference between samples, just use var1 for model fit
beta1 <- apply(var_lag1$BDraws,MARGIN=c(2,3),quantile,probs=0.5,na.rm=T)
beta2 <- apply(var_lag2$BDraws,MARGIN=c(2,3),quantile,probs=0.5,na.rm=T)
beta3 <- apply(var_lag3$BDraws,MARGIN=c(2,3),quantile,probs=0.5,na.rm=T)
beta4 <- apply(var_lag4$BDraws,MARGIN=c(2,3),quantile,probs=0.5,na.rm=T)

varcov1 <- apply(var_lag1$SDraws,MARGIN=c(2,3),quantile,probs=0.5,na.rm=T)
varcov2 <- apply(var_lag2$SDraws,MARGIN=c(2,3),quantile,probs=0.5,na.rm=T)
varcov3 <- apply(var_lag3$SDraws,MARGIN=c(2,3),quantile,probs=0.5,na.rm=T)
varcov4 <- apply(var_lag4$SDraws,MARGIN=c(2,3),quantile,probs=0.5,na.rm=T)

assessFunc <- function(xMat,beta,varcov){
yPred1 <- xMat%*%beta[,1]
yPred2 <- xMat%*%beta[,2]
yPred3 <- xMat%*%beta[,3]
yPred4 <- xMat%*%beta[,4]
#compute percentage errors, R2 for all
MAE1 <- mean(abs(yPred1 - DF[,1]),na.rm=T)
MAE2 <- mean(abs(yPred2 - DF[,2]),na.rm=T)
MAE3 <- mean(abs(yPred3 - DF[,3]),na.rm=T)
MAE4 <- mean(abs(yPred4 - DF[,4]),na.rm=T)

rs1 <- cor(yPred1,DF[,1],use="pairwise.complete.obs")
rs2 <- cor(yPred2,DF[,2],use="pairwise.complete.obs")
rs3 <- cor(yPred3,DF[,3],use="pairwise.complete.obs")
rs4 <- cor(yPred4,DF[,4],use="pairwise.complete.obs")

#logBF, DIC computation 
#compute acfs for each time series and store
acf1 <- acf(DF[,1]-yPred1,na.action=na.pass,plot=F)
acf2 <- acf(DF[,2]-yPred2,na.action=na.pass,plot=F)
acf3 <- acf(DF[,3]-yPred3,na.action=na.pass,plot=F)
acf4 <- acf(DF[,4]-yPred4,na.action=na.pass,plot=F)

assessList <- list(MAE1,MAE2,MAE3,MAE4,rs1,rs2,rs3,rs4,acf1,acf2,acf3,acf4,yPred1,yPred2,yPred3,yPred4,xMat)
names(assessList) <- c("MAE1","MAE2","MAE3","MAE4","rs1","rs2","rs3","rs4",
                       "acf1","acf2","acf3","acf4","yPred1","yPred2","yPred3","yPred4","xMat")
return(assessList)
}



dic <- function(y,varOutput){
  
  beta_p <- varOutput$BDraws
  s_p    <- varOutput$SDraws
  nlags <- (dim(beta_p)[2]-1)/4
  
  outX <- list()
  for (i in 1:nlags){
    outX[[i]] <- apply(y,MARGIN=2,FUN=Lag,k=i)
  }
  
  outX <- cbind(do.call(cbind,outX),1)
  llh <- list()
  for (i in 1:dim(beta_p)[1]){
    ##get predictions
    beta_temp <- beta_p[i,,]
    yt1 <- outX%*%beta_temp[,1]
    yt2 <- outX%*%beta_temp[,2]
    yt3 <- outX%*%beta_temp[,3]
    yt4 <- outX%*%beta_temp[,4]
    
    #compute llh value
    temp <- list()
    for (j in 1:nrow(DF)){
    temp[[j]] <- dmvnorm(DF[j,],
                        mean=cbind(yt1,yt2,yt3,yt4)[j,], 
                        sigma=s_p[i,,], log=T)
    
    }
    llh[[i]] <- sum(unlist(temp),na.rm=T)  ##llh for ith posterior draw
  }
  
  ##compute posterior mean llh 
  
  beta_m <- apply(varOutput$BDraws,MARGIN=c(2,3),mean)
  s_m    <- apply(varOutput$SDraws,MARGIN=c(2,3),mean)
  nlags <- (dim(beta_m)[1]-1)/4
  outX <- list()
  for (i in 1:nlags){
    outX[[i]] <- apply(y,MARGIN=2,FUN=Lag,k=i)
  }
  outX <- cbind(do.call(cbind,outX),1)
    ##get predictions
    yt1 <- outX%*%beta_m[,1]
    yt2 <- outX%*%beta_m[,2]
    yt3 <- outX%*%beta_m[,3]
    yt4 <- outX%*%beta_m[,4]
    
    #compute llh value
    temp <- list()
    for (j in 1:nrow(DF)){
      temp[[j]] <- dmvnorm(DF[j,],
                           mean=cbind(yt1,yt2,yt3,yt4)[j,], 
                           sigma=s_m, log=T)
      
    }
    llh_m <- sum(unlist(temp),na.rm=T)  ##llh for ith posterior draw
  
    #compute dic
    deviance_all <- -2* unlist(llh)
    deviance_exp <- -2*llh_m
    p_d <- mean(deviance) - deviance_exp
    dic <- p_d + mean(deviance) 
    return(dic)
}

assess1 <- assessFunc(xMat = xMat1,beta = beta1,varcov=varcov1)
assess2 <- assessFunc(xMat = xMat2,beta = beta2,varcov=varcov2)
assess3 <- assessFunc(xMat = xMat3,beta = beta3,varcov=varcov3)
assess4 <- assessFunc(xMat = xMat4,beta = beta4,varcov=varcov4)

dic1 <- dic(y=DF,varOutput=var_lag1)
dic2 <- dic(y=DF,varOutput=var_lag2) - dic1
dic3 <- dic(y=DF,varOutput=var_lag3) - dic1
dic4 <- dic(y=DF,varOutput=var_lag4) - dic1
dic1 <- 1
require(stargazer)
#create Xtable
lag1 <- unlist(assess1[1:8])
lag2 <- unlist(assess2[1:8])
lag3 <- unlist(assess3[1:8])
lag4 <- unlist(assess4[1:8])
df <- cbind(lag1,lag2,lag3,lag4)
df <- rbind(df,c(dic1,dic2,dic3,dic4))
df <- round(df,digits=3)
stargazer(df)

pdf("plot_appendix/qq.pdf",width=8.27,height=8.27)
##########################plot all in appendix#############################
pdf("~/nBox/sero_switch_th/plot_appendix/qq.pdf",width=8.27,height=8.27)
##########qqplots for appendix
par(mfrow=c(2,2),pty="s")
#lag1 model
qqplot(assess1$yPred1,assess1$xMat[,1],xlab="Prediction DENV-1",ylab="Actual DENV-1",main="1 Lag Model")
abline(a=0,b=1,lty=2)
qqplot(assess1$yPred2,assess1$xMat[,2],xlab="Prediction DENV-2",ylab="Actual DENV-2",main="1 Lag Model")
abline(a=0,b=1,lty=2)
qqplot(assess1$yPred3,assess1$xMat[,3],xlab="Prediction DENV-3",ylab="Actual DENV-3",main="1 Lag Model")
abline(a=0,b=1,lty=2)
qqplot(assess1$yPred4,assess1$xMat[,4],xlab="Prediction DENV-4",ylab="Actual DENV-4",main="1 Lag Model")
abline(a=0,b=1,lty=2)

#lag2 model
qqplot(assess2$yPred1,assess2$xMat[,1],xlab="Prediction DENV-1",ylab="Actual DENV-1",main="2 Lag Model")
abline(a=0,b=1,lty=2)
qqplot(assess2$yPred2,assess2$xMat[,2],xlab="Prediction DENV-2",ylab="Actual DENV-2",main="2 Lag Model")
abline(a=0,b=1,lty=2)
qqplot(assess2$yPred3,assess2$xMat[,3],xlab="Prediction DENV-3",ylab="Actual DENV-3",main="2 Lag Model")
abline(a=0,b=1,lty=2)
qqplot(assess2$yPred4,assess2$xMat[,4],xlab="Prediction DENV-4",ylab="Actual DENV-4",main="2 Lag Model")
abline(a=0,b=1,lty=2)

#lag3 model
qqplot(assess3$yPred1,assess3$xMat[,1],xlab="Prediction DENV-1",ylab="Actual DENV-1",main="3 Lag Model")
abline(a=0,b=1,lty=2)
qqplot(assess3$yPred2,assess3$xMat[,2],xlab="Prediction DENV-2",ylab="Actual DENV-2",main="3 Lag Model")
abline(a=0,b=1,lty=2)
qqplot(assess3$yPred3,assess3$xMat[,3],xlab="Prediction DENV-3",ylab="Actual DENV-3",main="3 Lag Model")
abline(a=0,b=1,lty=2)
qqplot(assess3$yPred4,assess3$xMat[,4],xlab="Prediction DENV-4",ylab="Actual DENV-4",main="3 Lag Model")
abline(a=0,b=1,lty=2)

#lag4 model
qqplot(assess4$yPred1,assess4$xMat[,1],xlab="Prediction DENV-1",ylab="Actual DENV-1",main="4 Lag Model")
abline(a=0,b=1,lty=2)
qqplot(assess4$yPred2,assess4$xMat[,2],xlab="Prediction DENV-2",ylab="Actual DENV-2",main="4 Lag Model")
abline(a=0,b=1,lty=2)
qqplot(assess4$yPred3,assess4$xMat[,3],xlab="Prediction DENV-3",ylab="Actual DENV-3",main="4 Lag Model")
abline(a=0,b=1,lty=2)
qqplot(assess4$yPred4,assess4$xMat[,4],xlab="Prediction DENV-4",ylab="Actual DENV-4",main="4 Lag Model")
abline(a=0,b=1,lty=2)

#acf
plot(assess1$acf1,main="DENV-1, 1 Lag Model")
plot(assess1$acf2,main="DENV-2, 1 Lag Model")
plot(assess1$acf3,main="DENV-3, 1 Lag Model")
plot(assess1$acf4,main="DENV-4, 1 Lag Model")

plot(assess2$acf1,main="DENV-1, 2 Lag Model")
plot(assess2$acf2,main="DENV-2, 2 Lag Model")
plot(assess2$acf3,main="DENV-3, 2 Lag Model")
plot(assess2$acf4,main="DENV-4, 2 Lag Model")

plot(assess3$acf1,main="DENV-1, 3 Lag Model")
plot(assess3$acf2,main="DENV-2, 3 Lag Model")
plot(assess3$acf3,main="DENV-3, 3 Lag Model")
plot(assess3$acf4,main="DENV-4, 3 Lag Model")

plot(assess4$acf1,main="DENV-1, 4 Lag Model")
plot(assess4$acf2,main="DENV-2, 4 Lag Model")
plot(assess4$acf3,main="DENV-3, 4 Lag Model")
plot(assess4$acf4,main="DENV-4, 4 Lag Model")

dev.off()