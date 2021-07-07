####################model fit#####################
####################model fit#####################
####################model fit#####################
rm(list=ls())
setwd("~/Google Drive/work_nea/sero_switch_sg/")
load("out/01outVAR.RData")
source("functions/hisDecomp.R")
require(quantmod)


#generate historical decompostion
hd <- VARhd(Estimation=var1,DF=as.matrix(DF))
#clean data
s1_hd <- -hd[,,1]
s2_hd <- -hd[,,2]
s3_hd <- -hd[,,3]
s4_hd <- -hd[,,4]

save.image("out/04_aHistDecomps.RData")
#add a b c d, add big bottom top labels
####### time series decomposition
####### time series decomposition
####### time series decomposition
####### time series decomposition

plotterHD <- function(df,   #time series of serotype of interest
                      HD,   #own serotypes
                      HDex, #all other serotypes
                      whatifCol,
                      whatifCol2,
                      labTS,
                      ylabTS,
                      labBP){
  whatif <- df-rowSums(HDex)
  plot(x=c(0,length(df)),
       y=c(min(c(0,whatif,rowSums(HDex),df),na.rm = T),
           max(df,na.rm=T)),col="white",xlab="",ylab=ylabTS,xaxt='n')
  polygon(y=c(0,df,0,0),
              x=c(0,seq(1,length(df)),length(df),0),col="black")
  lines(HD,col=whatifCol)
  lines(rowSums(HDex),col=whatifCol2)
  axis(side=1,at=seq(0,length(df),by=52),labels=F)
  axis(side=1,at=seq(0,length(df),by=52)+26,labels=c(paste("0",seq(6,9),sep=""),
                                                  seq(10,21)),tick=F)
  legend('topleft',
         legend=c("Observed","Decomposition (Own Serotype)","Decomposition (Other Serotypes)"),bty='n',
         col = c("black",whatifCol,whatifCol2),pch=15
         )
  mtext(labTS,side=3,adj=0)
  boxplot(cbind(df,HD,rowSums(HDex)),outline=F,col=c("black",whatifCol,whatifCol2),border="grey50",names=c("All","Own","Others"))
  mtext(labBP,side=3,adj=0)
}


pdf("plot/histDecomp3.pdf",width=8.27,height=10)
tester<-layout(matrix(c(1,1,2,
                        3,3,4,
                        5,5,6,
                        7,7,8), 4, 3, byrow = TRUE))
par(las=1,mar=c(2.5, 4.1, 2.5, 1.1),oma=c(2,2,2,2))

plotterHD(df=DF[,1],HD=s1_hd[,1],HDex=s1_hd[,-1],whatifCol = "lightblue3",whatifCol2="red",labTS = "A: Serotype 1",ylabTS="Dengue Case Counts",labBP="B")
plotterHD(df=DF[,2],HD=s2_hd[,2],HDex=s2_hd[,-2],whatifCol = "lightblue3",whatifCol2="red",labTS = "C: Serotype 2",ylabTS="Dengue Case Counts",labBP="D")
plotterHD(df=DF[,3],HD=s3_hd[,3],HDex=s3_hd[,-3],whatifCol = "lightblue3",whatifCol2="red",labTS = "E: Serotype 3",ylabTS="Dengue Case Counts",labBP="F")
plotterHD(df=DF[,4],HD=s4_hd[,4],HDex=s4_hd[,-4],whatifCol = "lightblue3",whatifCol2="red",labTS = "G: Serotype 4",ylabTS="Dengue Case Counts",labBP="H")
dev.off()

#BOXPLOT OF all historical decompositions#
#BOXPLOT OF all historical decompositions#
#BOXPLOT OF all historical decompositions#
pdf("plot/histDecomp4.pdf",width=8.27,height=4)
par(las=1,cex.axis=0.7)
boxplot(cbind(s1_hd,s2_hd,s3_hd,s4_hd),whisklty=1,outline=F,names=c("D1/1","D1/2","D1/3","D1/4",
                                                         "D2/1","D2/2","D2/3","D2/4",
                                                         "D3/1","D3/2","D3/3","D3/4",
                                                         "D4/1","D4/2","D4/3","D4/4"),
        col=c(rep("lightblue4",4),
              rep("lightblue3",4),
              rep("lightblue2",4),
              rep("lightblue",4)),ylab='Attributable Case Counts')
abline(v=4.5,lty=2,col="grey")
abline(v=8.5,lty=2,col="grey")
abline(v=12.5,lty=2,col="grey")
legend(x="topright",pch=15,legend=c("Historical Decomposition: DENV1",
                                    "Historical Decomposition: DENV2",
                                    "Historical Decomposition: DENV3",
                                    "Historical Decomposition: DENV4"), cex=0.70,
       col=c("lightblue4","lightblue3","lightblue2","lightblue",4),bty='n')
dev.off()


# 
# ###ridgeline plot very ugly
# ridgeDF <- cbind(s1_hd,s2_hd,s3_hd,s4_hd)
# ridgeDF <- apply(ridgeDF,MARGIN=2,function(x)density(x,na.rm=T))
# #normalize all density to 0-1 
# ridgeDF_y <- lapply(ridgeDF,function(x)(x$y - min(x$y))/(max(x$y)-min(x$y)))
# ridgeDF_x<- lapply(ridgeDF,function(x)x$x)
# ridgeDF_x <- do.call(cbind,ridgeDF_x)
# 
# xlim <- 200
# plot(x=c(-xlim,xlim),
#      y=c(0,length(ridgeDF_y)),col='white',ylab='Value',xlab='')
# 
# 
# ind <- seq(1,length(ridgeDF_y),by=2)
# polyLim <- 1000
# for (i in ind){
#   #draw squares
#   polygon(x=c(-polyLim,polyLim,polyLim,-polyLim,-polyLim),
#           y=c(i-1,i-1,i,i,i-1),col="grey89",border=F)}
# 
# for (i in rev(1:length(ridgeDF_y))){
#   #draw lines
#   lines(y=ridgeDF_y[[i]]+i-1,
#       x=ridgeDF_x[,i])
#   #draw shades
#   # polygon()  
# }
# 

# plotterTS <- function(yT,yDecomp,yDecompL,yDecompH,MainTitle,decompCol,decompCol2,leg=F,ylab=F,ylabel){
#   limiter <- cbind(yT,yDecomp,yDecompL,yDecompH)
#   plot(x=c(0,length(yT)),y=c(min(limiter),max(limiter))*1.2,col="white",xlab="",ylab="",xaxt='n')
#   
#   lines(yDecompL,lwd=6,col=decompCol2)
#   lines(yDecompH,lwd=6,col=decompCol2)
#   lines(yDecomp,lwd=1,col=decompCol)
#   lines(yT,lwd=1,col="black")
#   legend("topleft",legend=c("Cases","Decomposition"),lty=1,lwd=2,bty="n",col=c("black",decompCol))
#   mtext(text=MainTitle,side=3,cex=0.8)
#   mtext(text="weeks",side=1,padj=2.8,cex=0.8)
#   if (ylab==T) {par(las=0);mtext(ylabel,side=2,cex=0.8,padj=-4);par(las=1)}
#   box()
# }
# 
# pdf("~/nBox/sero_switch_th/plot/histDecomp2.pdf",width=8.27,height=10.69)
# par(mfrow=c(4,4),las=1,mar=c(1,2,1,1),pty='s',oma=c(3,2,3,1)) 
# 
# indDecomp <- seq(nrow(decomp1_m)/2,nrow(decomp1_m))
# indDF     <- seq(nrow(DF)-length(indDecomp)+1,nrow(DF))  
# plotterTS(yT=DF[indDF,1],yDecomp = decomp1_m[indDecomp,1],yDecompL = decomp1_t016[indDecomp,1],yDecompH = decomp1_t084[indDecomp,1],MainTitle = "S1 Decomposition",decompCol = "red2",decompCol2 = "pink",ylab=T,ylabel="Serotype 1")
# plotterTS(yT=DF[indDF,1],yDecomp = decomp1_m[indDecomp,2],yDecompL = decomp1_t016[indDecomp,2],yDecompH = decomp1_t084[indDecomp,2],MainTitle = "S2 Decomposition",decompCol = "cadetblue4",decompCol2 = "cadetblue2")
# plotterTS(yT=DF[indDF,1],yDecomp = decomp1_m[indDecomp,3],yDecompL = decomp1_t016[indDecomp,3],yDecompH = decomp1_t084[indDecomp,3],MainTitle = "S3 Decomposition",decompCol = "thistle4",decompCol2 = "thistle3")
# plotterTS(yT=DF[indDF,1],yDecomp = decomp1_m[indDecomp,4],yDecompL = decomp1_t016[indDecomp,4],yDecompH = decomp1_t084[indDecomp,4],MainTitle = "S4 Decomposition",decompCol = "peachpuff2",decompCol2 = "peru")
# 
# plotterTS(yT=DF[indDF,2],yDecomp = decomp2_m[indDecomp,1],yDecompL = decomp2_t016[indDecomp,1],yDecompH = decomp2_t084[indDecomp,1],MainTitle = "S1 Decomposition",decompCol = "red2",decompCol2 = "pink",ylab=T,ylabel="Serotype 2")
# plotterTS(yT=DF[indDF,2],yDecomp = decomp2_m[indDecomp,2],yDecompL = decomp2_t016[indDecomp,2],yDecompH = decomp2_t084[indDecomp,2],MainTitle = "S2 Decomposition",decompCol = "cadetblue4",decompCol2 = "cadetblue2")
# plotterTS(yT=DF[indDF,2],yDecomp = decomp2_m[indDecomp,3],yDecompL = decomp2_t016[indDecomp,3],yDecompH = decomp2_t084[indDecomp,3],MainTitle = "S3 Decomposition",decompCol = "thistle4",decompCol2 = "thistle3")
# plotterTS(yT=DF[indDF,2],yDecomp = decomp2_m[indDecomp,4],yDecompL = decomp2_t016[indDecomp,4],yDecompH = decomp2_t084[indDecomp,4],MainTitle = "S4 Decomposition",decompCol = "peachpuff2",decompCol2 = "peru")
# 
# plotterTS(yT=DF[indDF,3],yDecomp = decomp3_m[indDecomp,1],yDecompL = decomp3_t016[indDecomp,1],yDecompH = decomp3_t084[indDecomp,1],MainTitle = "S1 Decomposition",decompCol = "red2",decompCol2 = "pink",ylab=T,ylabel="Serotype 3")
# plotterTS(yT=DF[indDF,3],yDecomp = decomp3_m[indDecomp,2],yDecompL = decomp3_t016[indDecomp,2],yDecompH = decomp3_t084[indDecomp,2],MainTitle = "S2 Decomposition",decompCol = "cadetblue4",decompCol2 = "cadetblue2")
# plotterTS(yT=DF[indDF,3],yDecomp = decomp3_m[indDecomp,3],yDecompL = decomp3_t016[indDecomp,3],yDecompH = decomp3_t084[indDecomp,3],MainTitle = "S3 Decomposition",decompCol = "thistle4",decompCol2 = "thistle3")
# plotterTS(yT=DF[indDF,3],yDecomp = decomp3_m[indDecomp,4],yDecompL = decomp3_t016[indDecomp,4],yDecompH = decomp3_t084[indDecomp,4],MainTitle = "S4 Decomposition",decompCol = "peachpuff2",decompCol2 = "peru")
# 
# plotterTS(yT=DF[indDF,4],yDecomp = decomp4_m[indDecomp,1],yDecompL = decomp4_t016[indDecomp,1],yDecompH = decomp4_t084[indDecomp,1],MainTitle = "S1 Decomposition",decompCol = "red2",decompCol2 = "pink",ylab=T,ylabel="Serotype 4")
# plotterTS(yT=DF[indDF,4],yDecomp = decomp4_m[indDecomp,2],yDecompL = decomp4_t016[indDecomp,2],yDecompH = decomp4_t084[indDecomp,2],MainTitle = "S2 Decomposition",decompCol = "cadetblue4",decompCol2 = "cadetblue2")
# plotterTS(yT=DF[indDF,4],yDecomp = decomp4_m[indDecomp,3],yDecompL = decomp4_t016[indDecomp,3],yDecompH = decomp4_t084[indDecomp,3],MainTitle = "S3 Decomposition",decompCol = "thistle4",decompCol2 = "thistle3")
# plotterTS(yT=DF[indDF,4],yDecomp = decomp4_m[indDecomp,4],yDecompL = decomp4_t016[indDecomp,4],yDecompH = decomp4_t084[indDecomp,4],MainTitle = "S4 Decomposition",decompCol = "peachpuff2",decompCol2 = "peru")
# dev.off()
