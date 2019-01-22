#########################
### LINEAR REGRESSION ###
#########################
library(scales)
library(viridis)
path.to.git.repository <- "~/Documents/DATA/Codes/fsc-size-calibration"
setwd(path.to.git.repository)



### BEADS
beads740 <- read.csv("740-summary.csv")
beads740 <- beads740[order(beads740$size),]
beads740$fsc <- 10^((beads740$fsc.corr.high/2^16)*3.5)
id.740 <- which(beads740$size == 1) # find 1 micron beads

beads989<- read.csv("989-summary.csv")
beads989 <- beads989[order(beads989$size),]
beads989$fsc <- 10^((beads989$fsc.corr.high/2^16)*3.5)
id.989 <- which(beads989$size == 1) # find 1 micron beads

beads751 <- read.csv("751-summary.csv")
beads751 <- beads751[order(beads751$size),]
beads751$fsc <- 10^((beads751$fsc.corr.high/2^16)*3.5)
id.751<- which(beads751$size == 1) # find 1 micron beads

mie.b <- t(read.csv("~/Documents/DATA/Codes/mie-theory/meidata-beads.csv" ,header=F)) # beads
id <- findInterval(c(0.3, 0.5, 0.75, 1, 1.83, 3.1, 5.7) , mie.b[,1]) # find beads
id1 <- findInterval(1 , mie.b[,1]) # find 1 micron beads

png("Mie-beads-scatter.png",width=6, height=12, unit='in', res=100)

par(mfrow=c(3,1),pty='s')

inst <- 740; c <- 450*0.9; b <- 1.18
plot((mie.b[,2]/c)^b, mie.b[,1], col='red3', type='l', lwd=2, ylim=c(0.3,7), xlim=c(0.005,10), xaxt='n', yaxt='n', log='xy',xlab="Normalized scatter (dimensionless)", main=paste("#",inst),ylab=substitute(paste("Cell diameter (",mu,"m)")))
points((beads740$fsc/mean(beads740[id.1,'fsc'])), beads740$size, bg=alpha(viridis(nrow(beads740)),0.5),cex=2, pch=21)
axis(2, at=c(0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,1000), labels=c(0.005,0.01, 0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,1000), las=1)
axis(1, at=c(0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10),labels=c(0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10))
legend("topleft",legend=c(paste(unique(beads740$size), 'µm-beads'), "Mie-based model (index refraction = 1.6)"), bty='n',
        pch=c(rep(21,nrow(beads740)/2), NA), lwd=c(rep(NA,nrow(beads740)/2), 2),col=c(rep(1,nrow(beads740)/2),'red3'),
        pt.bg=alpha(c(viridis(nrow(beads740)/2), 'red3'),0.5))

inst <- 751; c <- 450; b <- 1
plot((mie.b[,2]/c)^b, mie.b[,1], col='red3', type='l', lwd=2, ylim=c(0.3,7), xlim=c(0.005,10), xaxt='n', yaxt='n', log='xy',xlab="Normalized scatter (dimensionless)", main=paste("#",inst),ylab=substitute(paste("Cell diameter (",mu,"m)")))
points((beads751$fsc/mean(beads751[id.1,'fsc'])), beads751$size, bg=alpha(viridis(nrow(beads740)),0.5),cex=2, pch=21)
axis(2, at=c(0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,1000), labels=c(0.005,0.01, 0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,1000), las=1)
axis(1, at=c(0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10),labels=c(0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10))
legend("topleft",legend=c(paste(unique(beads740$size), 'µm-beads'), "Mie-based model (index refraction = 1.6)"), bty='n',
        pch=c(rep(21,nrow(beads751)/2), NA), lwd=c(rep(NA,nrow(beads751)/2), 2),col=c(rep(1,nrow(beads751)/2),'red3'),
        pt.bg=alpha(c(viridis(nrow(beads751)/2), 'red3'),0.5))

inst <- 989; c <- 450; b <- 1.05
plot((mie.b[,2]/c)^b, mie.b[,1], col='red3', type='l', lwd=2, ylim=c(0.3,7), xlim=c(0.005,10), xaxt='n', yaxt='n', log='xy', xlab="Normalized scatter (dimensionless)", main=paste("#",inst),ylab=substitute(paste("Cell diameter (",mu,"m)")))
points((beads989$fsc/mean(beads989[id.1,'fsc'])), beads989$size, bg=alpha(viridis(nrow(beads740)),0.5),cex=2, pch=21)
axis(2, at=c(0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,1000), labels=c(0.005,0.01, 0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,1000), las=1)
axis(1, at=c(0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10),labels=c(0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10))
legend("topleft",legend=c(paste(unique(beads740$size), 'µm-beads'), "Mie-based model (index refraction = 1.6)"), bty='n',
        pch=c(rep(21,nrow(beads989)/2), NA), lwd=c(rep(NA,nrow(beads989)/2), 2),col=c(rep(1,nrow(beads989)/2),'red3'),
        pt.bg=alpha(c(viridis(nrow(beads989)/2), 'red3'),0.5))

dev.off()



### CULTURES


culture <- read.csv("scatter_calibration.csv")
culture$norm.fsc <- culture$fsc / culture$fsc.beads
culture$norm.fsc.sd <- culture$fsc.sd / culture$fsc.beads
culture$volume <- 4/3 * pi * (culture$diameter/2)^3
culture$volume.sd <- culture$volume * culture$diameter.sd/culture$diameter

culture2 <- aggregate(culture, by=list(culture$species), FUN=mean)
culture2 <- culture2[order(culture2$norm.fsc),]

mie <- read.csv("calibrated-mie.csv")
inst <- 740


png(paste0(inst,"-Size-scatter.png"),width=12, height=12, unit='in', res=100)

  par(mfrow=c(1,1), pty='s',cex=1.4)
  plot(culture2$norm.fsc, culture2$diameter, log='xy', pch=NA,ylab=substitute(paste("Cell diameter (",mu,"m)")), xlab="Normalized scatter (dimensionless)",cex=2, xaxt='n', yaxt='n', xlim=c(0.002,10), ylim=c(0.5,20), main=paste("#",inst))
  with(culture2, arrows(norm.fsc, diameter-culture2$diameter.sd, norm.fsc, diameter + culture2$diameter.sd,  code = 3, length=0, col='grey',lwd=2))
  with(culture2, arrows(norm.fsc-norm.fsc.sd, diameter, norm.fsc+norm.fsc.sd, diameter,  code = 3, length=0,col='grey',lwd=2))
  lines(mie$scatter, mie[,paste0("diam_",inst,"_mid")], col='red3', lwd=2)
  lines(mie$scatter, mie[,paste0("diam_",inst,"_upr")], col='grey', lwd=2)
  lines(mie$scatter, mie[,paste0("diam_",inst,"_lwr")], col='grey', lwd=2)
  axis(1, at=c(0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5))
  axis(2, at=c(0.1,0.2,0.5,1,2,5,10,20),las=1)
  points(culture2$norm.fsc, culture2$diameter, bg=alpha(viridis(nrow(culture2)),0.5), pch=21,cex=2)
  legend("topleft",legend=c(as.vector(culture2$Group.1),"Mie-based model (index refraction = 1.031 +/- 0.014)","beads"), pch=c(rep(21,nrow(culture2)),NA, 13), lwd=c(rep(NA,nrow(culture2)),2, NA), bty='n',
            pt.bg=alpha(viridis(nrow(culture2)),0.5), col=c(rep(1,nrow(culture2)),'red3','red3'), text.font=c(rep(3,nrow(culture2)),1,1))

dev.off()
