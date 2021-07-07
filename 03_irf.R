rm(list=ls())
setwd("~/Google Drive/work_nea/sero_switch_sg/")
load("out/01outVAR.RData")
#get median IRF for all
irf1_m <- apply(var1$IRFS,MARGIN=c(2,3),quantile,probs=c(0.5))
irf2_m <- apply(var2$IRFS,MARGIN=c(2,3),quantile,probs=c(0.5))
irf3_m <- apply(var3$IRFS,MARGIN=c(2,3),quantile,probs=c(0.5))
irf4_m <- apply(var4$IRFS,MARGIN=c(2,3),quantile,probs=c(0.5))
#q1 IRF for all
irf1_q1 <- apply(var1$IRFS,MARGIN=c(2,3),quantile,probs=c(0.16))
irf2_q1 <- apply(var2$IRFS,MARGIN=c(2,3),quantile,probs=c(0.16))
irf3_q1 <- apply(var3$IRFS,MARGIN=c(2,3),quantile,probs=c(0.16))
irf4_q1 <- apply(var4$IRFS,MARGIN=c(2,3),quantile,probs=c(0.16))
#q3 IRF for all
irf1_q3 <- apply(var1$IRFS,MARGIN=c(2,3),quantile,probs=c(0.84))
irf2_q3 <- apply(var2$IRFS,MARGIN=c(2,3),quantile,probs=c(0.84))
irf3_q3 <- apply(var3$IRFS,MARGIN=c(2,3),quantile,probs=c(0.84))
irf4_q3 <- apply(var4$IRFS,MARGIN=c(2,3),quantile,probs=c(0.84))

plotter <- function(q1,q3,m,samps,MainTitle,sampCol,qCol,leg=F){
plot(x=c(0,50),y=c(min(q1),max(q3))*1.2,col="white",xlab="",ylab="",xaxt='n')
for (i in 1:50){
  ind <- round(runif(1,1,nrow(samps)))
  lines(samps[ind,],col=sampCol) 
}
  lines(m,lwd=2)
  lines(q1,col=qCol,lwd=2)
  lines(q3,col=qCol,lwd=2)
  abline(h=0,lty=2,col="black")
  abline(v=0,lty=2,col="black")
  if (leg==T){
  legend("bottomright",legend=c("Median","16, 84 Tilde"),lty=1,lwd=2,bty="n",col=c("black",qCol))}
  mtext(text=MainTitle,side=3,cex=0.8)
  mtext(text="weeks",side=1,padj=2.8,cex=0.8)
  axis(side=1,at=c(0,10,20,30,40,50),labels=c(0,10,20,30,40,50)*2)
  box()
}

pdf("plot/irf.pdf",width=8.27,height=10.69)
par(mfrow=c(4,4),las=1,mar=c(1,2,1,1),pty='s',oma=c(3,2,3,0))
plotter(q1=irf1_q1[,1],q3=irf1_q3[,1],m=irf1_m[,1],samps=var1$IRFS[,,1],MainTitle = "Shock Serotype 1",sampCol="pink",qCol="red2")
mtext("Case Counts (Serotype 1)",side=2,las=0,cex=0.8,padj=-3.1)
plotter(q1=irf2_q1[,1],q3=irf2_q3[,1],m=irf2_m[,1],samps=var2$IRFS[,,1],MainTitle = "Shock Serotype 2",sampCol="cadetblue2",qCol="cadetblue4")
plotter(q1=irf3_q1[,1],q3=irf3_q3[,1],m=irf3_m[,1],samps=var3$IRFS[,,1],MainTitle = "Shock Serotype 3",sampCol="thistle1",qCol="thistle4")
plotter(q1=irf4_q1[,1],q3=irf4_q3[,1],m=irf4_m[,1],samps=var4$IRFS[,,1],MainTitle = "Shock Serotype 4",sampCol="navajowhite2",qCol="peru",leg=T)

plotter(q1=irf1_q1[,2],q3=irf1_q3[,2],m=irf1_m[,2],samps=var1$IRFS[,,2],MainTitle = "Shock Serotype 1",sampCol="pink",qCol="red2")
mtext("Case Counts (Serotype 2)",side=2,las=0,cex=0.8,padj=-3.1)
plotter(q1=irf2_q1[,2],q3=irf2_q3[,2],m=irf2_m[,2],samps=var2$IRFS[,,2],MainTitle = "Shock Serotype 2",sampCol="cadetblue2",qCol="cadetblue4")
plotter(q1=irf3_q1[,2],q3=irf3_q3[,2],m=irf3_m[,2],samps=var3$IRFS[,,2],MainTitle = "Shock Serotype 3",sampCol="thistle1",qCol="thistle4")
plotter(q1=irf4_q1[,2],q3=irf4_q3[,2],m=irf4_m[,2],samps=var4$IRFS[,,2],MainTitle = "Shock Serotype 4",sampCol="navajowhite2",qCol="peru",leg=T)

plotter(q1=irf1_q1[,3],q3=irf1_q3[,3],m=irf1_m[,3],samps=var1$IRFS[,,3],MainTitle = "Shock Serotype 1",sampCol="pink",qCol="red2")
mtext("Case Counts (Serotype 3)",side=2,las=0,cex=0.8,padj=-3.1)
plotter(q1=irf2_q1[,3],q3=irf2_q3[,3],m=irf2_m[,3],samps=var2$IRFS[,,3],MainTitle = "Shock Serotype 2",sampCol="cadetblue2",qCol="cadetblue4")
plotter(q1=irf3_q1[,3],q3=irf3_q3[,3],m=irf3_m[,3],samps=var3$IRFS[,,3],MainTitle = "Shock Serotype 3",sampCol="thistle1",qCol="thistle4")
plotter(q1=irf4_q1[,3],q3=irf4_q3[,3],m=irf4_m[,3],samps=var4$IRFS[,,3],MainTitle = "Shock Serotype 4",sampCol="navajowhite2",qCol="peru",leg=T)


plotter(q1=irf1_q1[,4],q3=irf1_q3[,4],m=irf1_m[,4],samps=var1$IRFS[,,4],MainTitle = "Shock Serotype 1",sampCol="pink",qCol="red2")
mtext("Case Counts (Serotype 4)",side=2,las=0,cex=0.8,padj=-3.1)
plotter(q1=irf2_q1[,4],q3=irf2_q3[,4],m=irf2_m[,4],samps=var2$IRFS[,,4],MainTitle = "Shock Serotype 2",sampCol="cadetblue2",qCol="cadetblue4")
plotter(q1=irf3_q1[,4],q3=irf3_q3[,4],m=irf3_m[,4],samps=var3$IRFS[,,4],MainTitle = "Shock Serotype 3",sampCol="thistle1",qCol="thistle4")
plotter(q1=irf4_q1[,4],q3=irf4_q3[,4],m=irf4_m[,4],samps=var4$IRFS[,,4],MainTitle = "Shock Serotype 4",sampCol="navajowhite2",qCol="peru",leg=T)
dev.off()
