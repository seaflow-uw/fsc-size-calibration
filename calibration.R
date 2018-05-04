#########################
### LINEAR REGRESSION ###
#########################
library(lmodel2)
path.to.git.repository <- "~/Documents/DATA/Codes/fsc-volume-calibration"
setwd(path.to.git.repository)
culture <- read.csv("volume_calibration.csv")
### WARNING !!! ###
# Calibration is only perform for Cyanobacteria in cultures (Pro & Syn), picoeuks have been excluded from the analysis
# linear regression type II
reg <- lmodel2(volume ~ norm.fsc, data=log(culture[15:22,c("volume","norm.fsc")],10))


png("Seaflow-volume-scatter.png",width=12, height=12, unit='in', res=100)

par(mfrow=c(1,1), pty='s',cex=1.4)
plot(culture$norm.fsc, culture$volume, log='xy', bg='grey', pch=21,ylab=substitute(paste("Cell Volume (",mu,"m"^{3},")")), xlab="Normalized scatter (dimensionless)",cex=2, xaxt='n', yaxt='n')
points(culture$norm.fsc[15:22], culture$volume[15:22], bg='red3', pch=21,cex=2)
with(culture, arrows(norm.fsc, volume - merge$volume.sd, norm.fsc, volume + merge$volume.sd,  code = 3, length=0))
with(culture, arrows(norm.fsc-norm.fsc.sd, volume, norm.fsc+norm.fsc.sd, volume,  code = 3, length=0))
axis(1, at=c(0.05,0.1,0.5,1,5), labels=c(0.05,0.1,0.5,1,5))
axis(2, at=c(1,5,10,50,100,500,1000), labels=c(1,5,10,50,100,500,1000))
par(new=T)
plot(log10(culture$norm.fsc), log10(culture$volume), yaxt='n',xaxt='n',xlab=NA, ylab=NA,pch=NA, bty='n')
abline(b=reg$regression.results$Slope[1], a=reg$regression.results$Intercept[1], col=2,lwd=2)
abline(b=reg$confidence.intervals[1,4], a=reg$confidence.intervals[1,2], col='grey',lwd=2)
abline(b=reg$confidence.intervals[1,5], a=reg$confidence.intervals[1,3], col='grey',lwd=2)
legend("topleft", legend=bquote(paste("Volume=",.(round(10^reg$regression.results$Intercept[1],3)),"(scatter"^{.(round(reg$regression.results$Slope[1],3))},")")), bty='n',cex=2)

dev.off()
