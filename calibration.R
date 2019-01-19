#########################
### LINEAR REGRESSION ###
#########################
library(scales)
library(viridis)
path.to.git.repository <- "~/Documents/DATA/Codes/fsc-size-calibration"
setwd(path.to.git.repository)

culture <- read.csv("scatter_calibration.csv")
culture$norm.fsc <- culture$fsc / culture$fsc.beads
culture$norm.fsc.sd <- culture$fsc.sd / culture$fsc.beads
culture$volume <- 4/3 * pi * (culture$diameter/2)^3
culture$volume.sd <- culture$volume * culture$diameter.sd/culture$diameter

culture2 <- aggregate(culture, by=list(culture$species), FUN=mean)
culture2 <- culture2[order(culture2$norm.fsc),]

beads <-




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

             points(1,1,pch=13, col='red3',lwd=2, cex=1.5) # beads location based on  corrected Mie theory

dev.off()
