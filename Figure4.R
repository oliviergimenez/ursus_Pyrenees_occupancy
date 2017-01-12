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

# psi/Occupancy 
logit_psi <- -1.27 -1.52 * DTHM + 0.63*LGRT + 0.24*CVBS -0.42*RUG + 0.12*CVFR
data_psi <- 1/(1+exp(-logit_psi))
data_psi <- as.vector(data_psi)

# gam/Colonization 
logit_gamma <- -7.26 + 3.38*RUG -3.07*DTHM +0.6*CVFR
data_gamma <- 1/(1+exp(-logit_gamma))
data_gamma <- as.vector(data_gamma)

# eps/Extinction 
logit_eps <- -7.6 - 5.82*DTHM - 1.48*LGRT
data_eps <- 1/(1+exp(-logit_eps))
data_eps <- as.vector(data_eps)

# p/Détection 
logit_p <- -2.36 + 0.62*RUG + 0.58*DTHM +0.22*CVFR - 1.4*AREA
data_p <- 1/(1+exp(-logit_p))
data_p <- as.vector(data_p)

psi <- gamma <- eps <- p <- rep(NA, dim(data_cov)[1])
a <- 1
for (i in 1:dim(data_cov)[1]){
	if (test_suivi[i] == 1){
		psi[i] <- data_psi[a]
		gamma[i] <- data_gamma[a]
		eps[i] <- data_eps[a]
		p[i] <- data_p[a]
		a <- a+1
	}	
}

data_cov$psi <- psi
data_cov$gamma <- gamma
data_cov$eps <- eps
data_cov$p <- p

ppi <- 300
tiff("figure4.tiff", width=6*ppi, height=6*ppi, res=ppi)

par(mar=c(0.2, 0.2, 0.2, 0.2), mfrow=c(2,2),
     oma = c(0.2, 0.2, 0.2, 0.2))

# Detection

min <- 0 
cex_test <- 0.9
yi_test <- 0.65

nclr <- 9
plotclr <- brewer.pal(nclr,"Greys")
class <- classIntervals(data_cov$p, nclr, style="quantile")
brks <- round(class$brks, 2)
class <- classIntervals(data_cov$p, nclr, style = "fixed",
   fixedBreaks = brks)
colcode <- findColours(class, plotclr)

plot(sousmassif.rg, axes=TRUE,border="gray", xaxt="n", yaxt="n")
plot(sousmassif.rg, col=colcode, add=TRUE)
title("A - Detection", line=-1.5, cex.main=2)
legend("bottomleft", legend=names(attr(colcode, "table")), title="Probability", fill=attr(colcode, "palette"), cex=cex_test, bty="n",y.intersp=yi_test)


# Initial Occupancy

min <- 0 
cex_test <- 0.9
yi_test <- 0.65

nclr <- 9
plotclr <- brewer.pal(nclr,"Greys")
data_cov$psi[which.max(data_cov$psi)] = signif(data_cov$psi[which.max(data_cov$psi)],2)
class <- classIntervals(data_cov$psi, nclr, style="quantile")
brks <- round(class$brks, 2)
class <- classIntervals(data_cov$psi, nclr, style = "fixed",
   fixedBreaks = brks)
colcode <- findColours(class, plotclr)

plot(sousmassif.rg, axes=TRUE,border="gray", xaxt="n", yaxt="n")
plot(sousmassif.rg, col=colcode, add=TRUE)
title("B - Initial occupancy", line=-1.5, cex.main=1.8)
legend("bottomleft", legend=names(attr(colcode, "table")), title="Probability", fill=attr(colcode, "palette"), cex=cex_test, bty="n",y.intersp=yi_test)

# Colonisation

min <- 0 
cex_test <- 0.9
yi_test <- 0.65

nclr <- 9
plotclr <- brewer.pal(nclr,"Greys")
class <- classIntervals(data_cov$gamma, nclr, style="quantile")
brks <- round(class$brks, 6)
class <- classIntervals(data_cov$gamma, nclr, style = "fixed",
   fixedBreaks = brks)
colcode <- findColours(class, plotclr)

plot(sousmassif.rg, axes=TRUE,border="gray", xaxt="n", yaxt="n")
plot(sousmassif.rg, col=colcode, add=TRUE)
title("C - Colonization", line=-1.5, cex.main=2)
legend("bottomleft", legend=names(attr(colcode, "table")), title="Probability", fill=attr(colcode, "palette"), cex=cex_test, bty="n",y.intersp=yi_test)

# Extinction

min <- 0 
cex_test <- 0.9
yi_test <- 0.65

nclr <- 9
plotclr <- brewer.pal(nclr,"Greys")
data_cov$eps[which.max(data_cov$eps)] = signif(data_cov$eps[which.max(data_cov$eps)],2)
class <- classIntervals(data_cov$eps, nclr, style="quantile")
brks <- round(class$brks, 8)
class <- classIntervals(data_cov$eps, nclr, style = "fixed",
   fixedBreaks = brks)
colcode <- findColours(class, plotclr)

plot(sousmassif.rg, axes=TRUE,border="gray",xaxt="n", yaxt="n")
plot(sousmassif.rg, col=colcode, add=TRUE)
title("D - Extinction", line=-1.5, cex.main=2)
legend("bottomleft", legend=names(attr(colcode, "table")), title="Probability", fill=attr(colcode, "palette"), cex=cex_test, bty="n",y.intersp=yi_test)


dev.off()


