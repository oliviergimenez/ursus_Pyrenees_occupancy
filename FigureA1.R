# packages
library(rgdal)
library(classInt)
library(RColorBrewer)

# read in det/non-det data
data_occ <- read.csv('Bear_OccSM_0814_mod.csv', header=FALSE, sep=",")

# sites ids
site_list <- data_occ[,1]

# get subsections and associated covariates 
sousmassif.rg <- readOGR(".", "sousmassif_og")

# get covariates
data_cov <- sousmassif.rg@data

# filter to keep only subsections with monitoring
cov <- subset(data_cov,select=c('Numero','alt_moy','tri_moy','dens_my','prc_frt','prc_shr','prc_rds','cnnct_f','diff_hm','Area'))
cov2 = NULL
test_suivi <- rep(0, dim(cov)[1])
for (i in 1:nrow(cov)){
	if (sum(cov[i,'Numero'] == site_list) != 0){
		cov2 <- rbind(cov2,cov[i,])
		test_suivi[i] <- 1
	}
}

# standardize covariates
RUG <- cov2$tri_moy	# roughness
RUG <- (RUG-mean(RUG))/sd(RUG)
DTHM <- cov2$dens_my # human density
DTHM <- (DTHM-mean(DTHM))/sd(DTHM)
CVFR <- cov2$prc_frt # forest cover
CVFR <- (CVFR-mean(CVFR))/sd(CVFR)
CVBS <- cov2$prc_shr # shrub cover
CVBS <- (CVBS-mean(CVBS))/sd(CVBS)
LGRT <- cov2$prc_rds # road length
LGRT <- (LGRT-mean(LGRT))/sd(LGRT)
AREA <- (cov2$Area/1000000 - mean(cov2$Area/1000000))/sd(cov2$Area/1000000) # subsection size

ppi <- 300
tiff("figureA1.tiff", width=6*ppi, height=6*ppi, res=ppi)

par(mar=c(0.2, 0.2, 0.2, 0.2), mfrow=c(3,2),
     oma = c(0.2, 0.2, 0.2, 0.2))

cols <- brewer.pal(n=9,name="Greys")

# roughness
lcols <- cut(data_cov$tri_moy,breaks=quantile(data_cov$tri_moy,probs=seq(0, 1, 0.11)),labels=cols)
plot(sousmassif.rg, axes=TRUE,border="gray",xaxt="n", yaxt="n")
plot(sousmassif.rg, col=as.character(lcols), add=TRUE) 
title('Roughness', line=-1.5, cex.main=2)

# forest cover
lcols <- cut(data_cov$prc_frt,breaks=quantile(data_cov$prc_frt,probs=seq(0, 1, 0.11)),labels=cols)
plot(sousmassif.rg, axes=TRUE,border="gray",xaxt="n", yaxt="n")
plot(sousmassif.rg, col=as.character(lcols), add=TRUE)
title('Forest cover', line=-1.5, cex.main=2)

# shrub cover
lcols <- cut(data_cov$prc_shr,breaks=quantile(data_cov$prc_shr,probs=seq(0, 1, 0.11)),labels=cols)
plot(sousmassif.rg, axes=TRUE,border="gray",xaxt="n", yaxt="n")
plot(sousmassif.rg, col=as.character(lcols), add=TRUE)
title('Shrub cover', line=-1.5, cex.main=2)

# road length
lcols <- cut(data_cov$prc_rds,breaks=quantile(data_cov$prc_rds,probs=seq(0, 1, 0.11)),labels=cols)
plot(sousmassif.rg, axes=TRUE,border="gray",xaxt="n", yaxt="n")
plot(sousmassif.rg, col=as.character(lcols), add=TRUE)
title('Road length', line=-1.5, cex.main=2)

# human density
lcols <- cut(data_cov$dens_my,breaks=quantile(data_cov$dens_my,probs=seq(0, 1, 0.11)),labels=cols)
plot(sousmassif.rg, axes=TRUE,border="gray",xaxt="n", yaxt="n")
plot(sousmassif.rg, col=as.character(lcols), add=TRUE)
title('Human density', line=-1.5, cex.main=2)

dev.off()
