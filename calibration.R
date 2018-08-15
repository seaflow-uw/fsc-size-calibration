#########################
### LINEAR REGRESSION ###
#########################
library(scales)
path.to.git.repository <- "~/Documents/DATA/Codes/fsc-volume-calibration"
setwd(path.to.git.repository)
culture <- read.csv("scatter_calibration.csv")
culture <- culture[order(culture$norm.fsc),]

# convert diameter to volume, assuming spherical shape
culture$volume <- 4/3 * pi * (culture$diameter/2)^3
culture$volume.sd <- culture$volume * culture$diameter.sd/culture$diameter


df <- t(read.table("meidata.txt"))
df <- data.frame(cbind(diam=df[1:100,], scatter=df[101:200,]))
write.csv(df, "meidata.csv", row.names=F, quote=F)
mie <- read.csv("meidata.csv")

par(mfrow=c(1,1))
plot(mie$scatter/17000, mie$diam, log='xy',type='l', ylim=c(0.1,15), xlab="FSC",ylab="Diameter (µm)")
f <- mie[which(mie$diam == 1),'scatter']
with(culture, arrows(norm.fsc, diameter-culture$diameter.sd, norm.fsc, diameter + culture$diameter.sd,  code = 3, length=0, col='grey',lwd=2))
with(culture, arrows(norm.fsc-norm.fsc.sd, diameter, norm.fsc+norm.fsc.sd, diameter,  code = 3, length=0,col='grey',lwd=2))
points(culture$norm.fsc,culture$diameter,pch=16, col=2,cex=2)



reg <- lm(diameter ~ norm.fsc, data=log(culture[1:8,c("diameter","norm.fsc")],10))
reg2 <- lm(diameter ~ norm.fsc, data=log(culture[9:25,c("diameter","norm.fsc")],10))

png("Seaflow-diameter-scatter.png",width=12, height=12, unit='in', res=100)

par(mfrow=c(1,1), pty='s',cex=1.4)
plot(culture$norm.fsc, culture$diameter, log='xy', bg=alpha('grey',0.3), pch=21,ylab=substitute(paste("Cell diameter (",mu,"m)")), xlab="Normalized scatter (dimensionless)",cex=2, xaxt='n', yaxt='n')
with(culture, arrows(norm.fsc, diameter- culture$diameter.sd, norm.fsc, diameter + culture$diameter.sd,  code = 3, length=0))
with(culture, arrows(norm.fsc-norm.fsc.sd, diameter, norm.fsc+norm.fsc.sd, diameter,  code = 3, length=0))
axis(1, at=c(0.05,0.1,0.5,1,5), labels=c(0.05,0.1,0.5,1,5))
axis(2, at=c(0.1,1,2,5,10,20), labels=c(0.1,1,2,5,10,20))
par(new=T)
plot(log10(culture$norm.fsc), log10(culture$diameter), yaxt='n',xaxt='n',xlab=NA, ylab=NA,pch=NA, bty='n')
lines(x=log10(culture$norm.fsc),predict(reg, newdata=data.frame(norm.fsc=log10(culture$norm.fsc)),interval='predict')[,"fit"], col='red3',lwd=2 )
lines(x=log10(culture$norm.fsc),predict(reg, newdata=data.frame(norm.fsc=log10(culture$norm.fsc)),interval='predict')[,"upr"], col='grey',lwd=2 )
lines(x=log10(culture$norm.fsc),predict(reg, newdata=data.frame(norm.fsc=log10(culture$norm.fsc)),interval='predict')[,"lwr"], col='grey',lwd=2 )
legend("topleft", legend=bquote(paste("Diam_cyano=",.(round(10^reg$coefficients[1],3)),"(scatter"^{.(round(reg$coefficients[2],3))},")")), bty='n',cex=1.5, text.col='red3')

lines(x=log10(culture$norm.fsc),predict(reg2, newdata=data.frame(norm.fsc=log10(culture$norm.fsc)),interval='predict')[,"fit"], col='seagreen3',lwd=2 )
lines(x=log10(culture$norm.fsc),predict(reg2, newdata=data.frame(norm.fsc=log10(culture$norm.fsc)),interval='predict')[,"upr"], col='grey',lwd=2 )
lines(x=log10(culture$norm.fsc),predict(reg2, newdata=data.frame(norm.fsc=log10(culture$norm.fsc)),interval='predict')[,"lwr"], col='grey',lwd=2 )
legend(y=log10(10), x=log10(0.022), legend=bquote(paste("Diam_pico=",.(round(10^reg2$coefficients[1],3)),"(scatter"^{.(round(reg2$coefficients[2],3))},")")), bty='n',cex=1.5, text.col='seagreen3')

dev.off()





reg <- lm(volume ~ norm.fsc, data=log(culture[1:8,c("volume","norm.fsc")],10))
reg2 <- lm(volume ~ norm.fsc, data=log(culture[9:25,c("volume","norm.fsc")],10))

png("Seaflow-volume-scatter.png",width=12, height=12, unit='in', res=100)

par(mfrow=c(1,1), pty='s',cex=1.4)
plot(culture$norm.fsc, culture$volume, log='xy', bg=alpha('grey',0.3), pch=21,ylab=substitute(paste("Cell Volume (",mu,"m"^{3},")")), xlab="Normalized scatter (dimensionless)",cex=2, xaxt='n', yaxt='n')
with(culture, arrows(norm.fsc, volume - culture$volume.sd, norm.fsc, volume + culture$volume.sd,  code = 3, length=0))
with(culture, arrows(norm.fsc-norm.fsc.sd, volume, norm.fsc+norm.fsc.sd, volume,  code = 3, length=0))
axis(1, at=c(0.05,0.1,0.5,1,5), labels=c(0.05,0.1,0.5,1,5))
axis(2, at=c(1,5,10,50,100,500,1000), labels=c(1,5,10,50,100,500,1000))
par(new=T)
plot(log10(culture$norm.fsc), log10(culture$volume), yaxt='n',xaxt='n',xlab=NA, ylab=NA,pch=NA, bty='n')
lines(x=log10(culture$norm.fsc),predict(reg, newdata=data.frame(norm.fsc=log10(culture$norm.fsc)),interval='predict')[,"fit"], col='red3',lwd=2 )
lines(x=log10(culture$norm.fsc),predict(reg, newdata=data.frame(norm.fsc=log10(culture$norm.fsc)),interval='predict')[,"upr"], col='grey',lwd=2 )
lines(x=log10(culture$norm.fsc),predict(reg, newdata=data.frame(norm.fsc=log10(culture$norm.fsc)),interval='predict')[,"lwr"], col='grey',lwd=2 )
legend("topleft", legend=bquote(paste("Vol_cyano=",.(round(10^reg$coefficients[1],3)),"(scatter"^{.(round(reg$coefficients[2],3))},")")), bty='n',cex=1.5, text.col='red3')

lines(x=log10(culture$norm.fsc),predict(reg2, newdata=data.frame(norm.fsc=log10(culture$norm.fsc)),interval='predict')[,"fit"], col='seagreen3',lwd=2 )
lines(x=log10(culture$norm.fsc),predict(reg2, newdata=data.frame(norm.fsc=log10(culture$norm.fsc)),interval='predict')[,"upr"], col='grey',lwd=2 )
lines(x=log10(culture$norm.fsc),predict(reg2, newdata=data.frame(norm.fsc=log10(culture$norm.fsc)),interval='predict')[,"lwr"], col='grey',lwd=2 )
legend(y=log10(500), x=log10(0.022), legend=bquote(paste("Vol_pico=",.(round(10^reg2$coefficients[1],3)),"(scatter"^{.(round(reg2$coefficients[2],3))},")")), bty='n',cex=1.5, text.col='seagreen3')


dev.off()
