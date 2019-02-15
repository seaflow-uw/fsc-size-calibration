library(scales)
library(viridis)
path.to.git.repository <- "~/Documents/DATA/Codes/fsc-size-calibration"
setwd(path.to.git.repository)

#############################
### scatter normalization ###
#############################

beads740 <- read.csv("740-summary.csv")
beads740$fsc <- 10^((beads740$fsc.corr.high/2^16)*3.5)
id.740 <- which(beads740$size == 1) # find 1 micron beads
beads740$normalized.fsc <- beads740$fsc/mean(beads740[id.740,'fsc'])
beads740 <- beads740[order(beads740$size),]

beads989<- read.csv("989-summary.csv")
beads989$fsc <- 10^((beads989$fsc.corr.high/2^16)*3.5)
id.989 <- which(beads989$size == 1) # find 1 micron beads
beads989$normalized.fsc <- beads989$fsc/mean(beads989[id.989,'fsc'])
beads989 <- beads989[order(beads989$size),]

beads751 <- read.csv("751-summary.csv")
beads751$fsc <- 10^((beads751$fsc.corr.high/2^16)*3.5)
id.751<- which(beads751$size == 1) # find 1 micron beads
beads751$normalized.fsc <- beads751 $fsc/mean(beads751 [id.751,'fsc'])
beads751 <- beads751[order(beads751$size),]


####################
### OPTIMIZATION ###
####################
library(DEoptim)

# Mie theory fitting
# n <- c(1.35/1.3371, 1.41/1.3371) # range given in Lehmuskero et al. Progr Oceanogr 2018
mie2 <- t(read.csv("meidata-1017.csv" ,header=F)) # low
mie1 <- t(read.csv("meidata-1031.csv" ,header=F)) # fit
mie3 <- t(read.csv("meidata-1045.csv" ,header=F)) # high
mie4 <- t(read.csv("meidata-beads.csv" ,header=F)) # particle of 1 µm, index of refraction of 1.6033

# optimization routine
sigma.lsq <- function(mie, beads, params){

     c <- params[1]
     b <- params[2]
     id <- findInterval(beads$size, mie[,1])
     scatt <- (mie[id,2]/c)^b

           df <- data.frame(obs=beads$normalized.fsc, pred=scatt)
      sigma <- mean(abs(df$obs - df$pred)/df$obs,na.rm=T)
      return(sigma)

  }




##########################
#### CALIBRATION BEADS ### Mie theory optimization
##########################
png("Mie-beads-scatter.png",width=12, height=6, unit='in', res=400)

par(mfrow=c(1,3), pty='s', cex=1.2)

for(inst in c(740,751,989)){

# inst <- 740


### Optimization

if(inst == 740) beads <- beads740
if(inst == 751) beads <- beads751
if(inst == 989) beads <- beads989

f <- function(params) sigma.lsq(mie=mie4, beads=beads, params)
res <- DEoptim(f, lower=c(0,0), upper=c(10000,2), control=DEoptim.control(itermax=1000, reltol=1e-8,trace=100, steptol=100, strategy=2, parallelType=0))
print(res$optim$bestval)
params <- res$optim$bestmem # optimized 'c' and 'b' values
print(params)


### CREATE LOOKUP TABLE
d <- 0.220 #Shalapyonok et al. 2001; 0.220 (Li et al. 1992, Veldhuis et al. 2004, and more studies agreed with 220 fg C um-3)
max.scatter <- 20
min.scatter <- 0.0001

b <- params[2]
c <- params[1]
scatter <- 10^(seq(log10(min(mie1[,2]/c)), log10(max(mie3[,2]/c)),length.out=10000))

s1 <- approx((mie1[,2]/c)^b, mie1[,1], xout=scatter)
s2 <- approx((mie2[,2]/c)^b, mie2[,1], xout=scatter)
s3 <- approx((mie3[,2]/c)^b, mie3[,1], xout=scatter)
s4 <- approx((mie1[,2]/c)^b, d*(4/3*pi*(0.5*mie1[,1])^3), xout=scatter)
s5 <- approx((mie2[,2]/c)^b, d*(4/3*pi*(0.5*mie2[,1])^3), xout=scatter)
s6 <- approx((mie3[,2]/c)^b, d*(4/3*pi*(0.5*mie3[,1])^3), xout=scatter)

		if(inst == 740){mie_740 <- data.frame(cbind(scatter=s1$x,
                                                  diam_740_mid=s1$y,diam_740_upr=s2$y,diam_740_lwr = s3$y,
                                                  Qc_740_mid=s4$y, Qc_740_upr=s5$y, Qc_740_lwr=s6$y))
                      	mie_740 <- subset(mie_740, scatter >= min.scatter & scatter <= max.scatter)
                     	}

        if(inst == 751){mie_751 <- data.frame(cbind(scatter=s1$x,
                                                  diam_751_mid=s1$y,diam_751_upr=s2$y,diam_751_lwr = s3$y,
                                                  Qc_751_mid=s4$y, Qc_751_upr=s5$y, Qc_751_lwr=s6$y))
                      	mie_751 <- subset(mie_751, scatter >= min.scatter & scatter <= max.scatter)
                      	}
	    if(inst == 989){mie_989 <- data.frame(cbind(scatter=s1$x,
                                                  diam_989_mid=s1$y,diam_989_upr=s2$y,diam_989_lwr = s3$y,
                                                  Qc_989_mid=s4$y, Qc_989_upr=s5$y, Qc_989_lwr=s6$y))
                      	mie_989 <- subset(mie_989, scatter >= min.scatter & scatter <= max.scatter)
						}



### PLOT Observartions vs Predicted cell size
plot(beads$normalized.fsc, beads$size,log='xy', xaxt='n',main=paste(inst), xlim=c(0.005,10), ylim=c(0.3,7), bg=alpha(viridis(nrow(beads)),0.5),cex=2, pch=21, xlab="scatter", ylab="size (µm)")
axis(1, at=c(0.01, 0.1, 1,10))
lines((mie4[,2]/params[1])^params[2], mie4[,1], col='red3')
legend("topleft",cex=0.5, legend=c(paste(unique(beads$size), 'µm-beads'), "Mie-based model (n = 1.6)"), bty='n', pch=c(rep(21,nrow(beads)/2), NA), lwd=c(rep(NA,nrow(beads)/2), 2),col=c(rep(1,nrow(beads)/2),'red3'), pt.bg=alpha(c(viridis(nrow(beads)/2), 'red3'),0.5))

}
dev.off()



mie <- data.frame(cbind(mie_740, mie_751[,-1], mie_989[,-1]))

write.csv(mie, "calibrated-mie.csv", row.names=F, quote=F)





#####################
### PHYTOPLANKTON ### validation
#####################


### SIZE
culture <- read.csv("scatter_calibration.csv")
culture$norm.fsc <- culture$fsc / culture$fsc.beads
culture$norm.fsc.sd <- culture$fsc.sd / culture$fsc.beads
culture$volume <- 4/3 * pi * (culture$diameter/2)^3
culture$volume.sd <- culture$volume * culture$diameter.sd/culture$diameter

culture2 <- aggregate(culture, by=list(culture$species), FUN=mean)
culture2 <- subset(culture2, Group.1 !="Phaeodactylum tricornutum") # remove non-spherical cells
culture2 <- culture2[order(culture2$norm.fsc),]

mie <- read.csv("calibrated-mie.csv")
inst <- 740


png(paste0(inst,"-Size-scatter.png"),width=6, height=6, unit='in', res=200)

  par(mfrow=c(1,1), pty='s',cex=1.2)
  plot(culture2$norm.fsc, culture2$diameter, log='xy', pch=NA,ylab=substitute(paste("Cell diameter (",mu,"m)")), xlab="Normalized scatter (dimensionless)",cex=2, xaxt='n', yaxt='n', xlim=c(0.002,10), ylim=c(0.5,20), main=paste(inst))
  with(culture2, arrows(norm.fsc, diameter-culture2$diameter.sd, norm.fsc, diameter + culture2$diameter.sd,  code = 3, length=0, col='grey',lwd=2))
  with(culture2, arrows(norm.fsc-norm.fsc.sd, diameter, norm.fsc+norm.fsc.sd, diameter,  code = 3, length=0,col='grey',lwd=2))
  lines(mie$scatter, mie[,paste0("diam_",inst,"_mid")], col='red3', lwd=2)
  lines(mie$scatter, mie[,paste0("diam_",inst,"_upr")], col='grey', lwd=2)
  lines(mie$scatter, mie[,paste0("diam_",inst,"_lwr")], col='grey', lwd=2)
  axis(1, at=c(0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5))
  axis(2, at=c(0.1,0.2,0.5,1,2,5,10,20),las=1)
  points(culture2$norm.fsc, culture2$diameter, bg=alpha(viridis(nrow(culture2)),0.5), pch=21,cex=2)
  legend("topleft",legend=c(as.vector(culture2$Group.1), "Mie-based model (n = 1.031 +/- 0.014)"), cex=0.5, pch=c(rep(21,nrow(culture2)),NA), lwd=c(rep(NA,nrow(culture2)),2), bty='n',
            pt.bg=alpha(viridis(nrow(culture2)),0.5), col=c(rep(1,nrow(culture2)),'red3'), text.font=c(rep(3,nrow(culture2)),1))

dev.off()
