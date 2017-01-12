ppi <- 300
tiff("figure2.tiff", width=6*ppi, height=6*ppi, res=ppi)

par(mfrow=c(2,2))

labels = c('Int.', 'RUG', 'DTHM', 'CVFR', 'CVBS', 'LGRT', 'AREA') 

# initial occupancy
mles = c(-1.27,-0.42,-1.52,0.12,0.24,0.63,NA)
names(mles) = labels
conf_int = matrix(c(-2.13, -0.41,
-1.91, 1.06,
-2.86, -0.17,
-0.63, 0.88,
-0.57, 1.04,
-0.64, 1.89,
0,0),byrow=T,ncol=2)
row.names(conf_int) = labels
plot(mles, 1:7, ylab='',yaxt="n",xlim=range(c(mles,conf_int[,1],conf_int[,2]),na.rm=T),main='Initial occupancy',xlab='')
axis(2,at=1:7,labels=labels,las=1)
segments(conf_int[,1], 1:7, conf_int[,2], 1:7)
abline(v=0, lty=2, col='red')

# colonization
mles = c(-7.26,3.38,-3.07,0.6,NA,NA,NA)
names(mles) = labels
conf_int = matrix(c(-16.19, 1.66,
-1, 7.76,
-8.49, 2.34,
-1.14, 2.35,
0,0,
0,0,
0,0),byrow=T,ncol=2)
row.names(conf_int) = labels
plot(mles, 1:7, ylab='',yaxt="n",xlim=range(c(mles,conf_int[,1],conf_int[,2]),na.rm=T),main='Colonization',xlab='')
axis(2,at=1:7,labels=labels,las=1)
segments(conf_int[,1], 1:7, conf_int[,2], 1:7)
abline(v=0, lty=2, col='red')

# extinction
mles = c(-7.6,NA,-5.82,NA,NA,-1.48,NA)
names(mles) = labels
conf_int = matrix(c(-13.91, -1.29,
0,0,
-11.02, -0.62,
0,0,
0,0,
-4.65, 1.68,
0,0),byrow=T,ncol=2)
row.names(conf_int) = labels
plot(mles, 1:7, ylab='',yaxt="n",xlim=range(c(mles,conf_int[,1],conf_int[,2]),na.rm=T),main='Extinction',xlab='')
axis(2,at=1:7,labels=labels,las=1)
segments(conf_int[,1], 1:7, conf_int[,2], 1:7)
abline(v=0, lty=2, col='red')

# detection
mles = c(-2.36, 0.62,0.58,0.22,NA,NA,-1.4)
names(mles) = labels
conf_int = matrix(c(-2.8, -1.92,
0.1, 1.14,
-0.31, 1.47,
-0.13, 0.58,
0,0,
0,0,
-2.18, -0.62),byrow=T,ncol=2)
row.names(conf_int) = labels
plot(mles, 1:7, ylab='',yaxt="n",xlim=range(c(mles,conf_int[,1],conf_int[,2]),na.rm=T),main='Detection',xlab='')
axis(2,at=1:7,labels=labels,las=1)
segments(conf_int[,1], 1:7, conf_int[,2], 1:7)
abline(v=0, lty=2, col='red')


dev.off()