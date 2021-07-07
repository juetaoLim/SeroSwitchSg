####################model fit#####################
####################model fit#####################
####################model fit#####################
rm(list=ls())
setwd("~/Google Drive/work_nea/sero_switch_sg/")
load("out/01outVAR.RData")
require(quantmod)
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
acf1 <- acf(DF[,1]-yPred1,na.action=na.pass)
acf2 <- acf(DF[,2]-yPred2,na.action=na.pass)
acf3 <- acf(DF[,3]-yPred3,na.action=na.pass)
acf4 <- acf(DF[,4]-yPred4,na.action=na.pass)

# require(mvtnorm)
#predictive simulation
# BETA <- apply(var1$BDraws,MARGIN=c(2,3),quantile,probs=0.5,na.rm=T)
# VARCOV <- apply(var1$SDraws,MARGIN=c(2,3),quantile,probs=0.5,na.rm=T)

# pps <- function(beta=BETA,varcov=VARCOV,iter=2000,subiter=600){
#   #initialize
#   store <- list()
#   for (j in 1:iter){
#   x <- matrix(data=c(rep(0,nrow(beta)-1),1),
#                 ncol=nrow(beta),
#                 nrow=1)
#   mu <- x %*% beta
#   temp <- list()
#   
#   #iterate
#   
#   for (i in 1:subiter){
#   temp[[i]] <- rmvnorm(n=1,mean=mu,sigma=varcov)
#   
#   #old lag 1-3 is new lag 2-4
#   
#   oldx1       <- x[1:3]
#   oldx2       <- x[5:7]
#   oldx3       <- x[9:11]
#   oldx4       <- x[13:15]
#   
#   #construct new design matrix
#   
#   newX        <- temp[[i]]
#   x           <- c(newX[1],
#                    oldx1,
#                    newX[2],
#                    oldx2,
#                    newX[3],
#                    oldx3,
#                    newX[4],
#                    oldx4,
#                    1)
#    mu <- x %*% beta
#   }
#   scaleTemp  <- do.call(rbind,temp)
#   scaleTemp  <- apply(scaleTemp,MARGIN=2,scale)
#   store[[j]] <- scaleTemp
#   }
#   return(store)
# }
# 
# ppsStore<- pps()

save.image("~/nBox/sero_switch_th/out/02assesment.RData")
