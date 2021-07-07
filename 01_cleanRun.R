rm(list=ls())
library(VARsignR)
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
var1 <- uhlig.reject(Y=DF,nlags=4,constrained = constr1,steps=300,draws=10000,constant=F, nkeep=10000)
var2 <- uhlig.reject(Y=DF,nlags=4,constrained = constr2,steps=300,draws=10000,constant=F, nkeep=10000)
var3 <- uhlig.reject(Y=DF,nlags=4,constrained = constr3,steps=300,draws=10000,constant=F, nkeep=10000)
var4 <- uhlig.reject(Y=DF,nlags=4,constrained = constr4,steps=300,draws=10000,constant=F, nkeep=10000)


#get IRFs
irfs1 <- var1$IRFS[,1:150,]
irfs2 <- var2$IRFS[,1:150,]
irfs3 <- var3$IRFS[,1:150,]
irfs4 <- var4$IRFS[,1:150,]
irfplot(irfs1)
irfplot(irfs2)
irfplot(irfs3)
irfplot(irfs4)

save.image("out/01outVAR.RData")

#with zika analysis
# 
# rm(list=ls())
# library(VARsignR)
# data <- read.csv("nBox/sero_switch_th/data/data_sg_wZika.csv")
# data <- data[1:220,]
# #get sero%
# 
# S1 <- data$S1/data$TotalSero*data$Dengue.cases
# S2 <- data$S2/data$TotalSero*data$Dengue.cases
# S3 <- data$S3/data$TotalSero*data$Dengue.cases
# S4 <- data$S4/data$TotalSero*data$Dengue.cases
# S5 <- data$Zika
# # # #first difference log
# # S1 <- diff(data$S1)
# # S2 <- diff(data$S2)
# # S3 <- diff(data$S3)
# # S4 <- diff(data$S4)
# # 
# # #first difference data
# # diffDF <- cbind(S1,S2,S3,S4)
# # diffDF <- as.ts(diffDF)
# 
# DF     <- cbind(S1,S2,S3,S4,S5)
# DF     <- as.ts(DF)
# 
# #setConstraints
# constr1 <- c(+1,-2,-3,-4,-5)
# constr2 <- c(+2,-1,-3,-4,-5)
# constr3 <- c(+3,-2,-1,-4,-5)
# constr4 <- c(+4,-2,-3,-1,-5) #DENV protective against ZIKA in cohort, but can enhance ZIKA infeciton, so no constrain 
# constr5 <- c(+5,+1,+2,+3,+4) #zika increases DENV risk
# #runMCMC using uhlig.reject
# var1 <- uhlig.reject(Y=DF,nlags=4,constrained = constr1,steps=300,draws=1000)
# var2 <- uhlig.reject(Y=DF,nlags=4,constrained = constr2,steps=300,draws=1000)
# var3 <- uhlig.reject(Y=DF,nlags=4,constrained = constr3,steps=300,draws=1000)
# var4 <- uhlig.reject(Y=DF,nlags=4,constrained = constr4,steps=300,draws=1000)
# var5 <- uhlig.reject(Y=DF,nlags=4,constrained = constr5,steps=300,draws=1000)
# 
# #get IRFs
# irfs1 <- var1$IRFS[,1:150,]
# irfs2 <- var2$IRFS[,1:150,]
# irfs3 <- var3$IRFS[,1:150,]
# irfs4 <- var4$IRFS[,1:150,]
# irfs5 <- var5$IRFS[,1:150,]
# irfplot(irfs1)
# irfplot(irfs2)
# irfplot(irfs3)
# irfplot(irfs4)
# irfplot(irfs5)
# 
# save.image("~/nBox/sero_switch_th/out/01outVAR_z.RData")
