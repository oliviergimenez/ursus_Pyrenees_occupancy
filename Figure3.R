# packages
library(rgdal)

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
AREA <- (cov2$Area/1000000 - mean(cov2$Area/1000000))/sd(cov2$Area/1000000) # subsection size

# function that backtransforms and plot param-cov relationship with conf int
plot_bear <- function(cov,cov_name,param_name,estmean,estvar,panel_name){
	grid = seq(min(cov),max(cov),length=100)
	logit_psi <- estmean[1] + estmean[2] * grid
	logit_psi_se <- rep(NA,length(grid))
	for (i in 1:length(grid)){
		grad = c(1,grid[i])
		logit_psi_se[i] = sqrt(t(grad)%*%estvar%*%grad) # delta-method
	}
	psi <- 1/(1+exp(-logit_psi))
	psi_ci_lower <- 1/(1+exp(-(logit_psi-2*logit_psi_se)))
	psi_ci_upper <- 1/(1+exp(-(logit_psi+2*logit_psi_se)))
	plot(grid, psi,type='l',ylim=c(0,1),xlab="",ylab="")
	polygon(c(rev(grid), grid), c(rev(psi_ci_upper), psi_ci_lower), col = 'grey80', border = NA)
	lines(grid, psi, col="black",ylim=c(0,1),type='l',lwd=2,pch=21)
	points(cov, rep(0,84), pch='|')
	title(panel_name, line=-1.5, cex.main=2)
	title(xlab=cov_name,ylab=param_name, line=2)
}

ppi <- 300
tiff("figure3.tiff", width=6*ppi, height=6*ppi, res=ppi)

par(mar=c(3, 3, 0.2, 0.2), mfrow=c(2,2),
     oma = c(0.2, 1, 0.2, 0.2))

## Occupancy x DTHM
cov = DTHM
cov_name = 'human density'
param_name = 'Initial occupancy'
estmean <- c(-1.27,-1.52)
estvar <- diag(c(0.44,0.69)^2)
panel_name <- ''
plot_bear(cov=cov,cov_name=cov_name,param_name=param_name,estmean=estmean,estvar=estvar,panel_name=panel_name)

## extinction x DTHM
cov = DTHM
cov_name = 'human density'
param_name = 'Extinction'
estmean <- c(-7.6,-5.82)
estvar <- diag(c(3.22,2.65)^2)
panel_name <- ''
plot_bear(cov,cov_name,param_name,estmean,estvar,panel_name)

## detection x RUG
cov = RUG
cov_name = 'roughness'
param_name = 'Detection'
estmean <- c(-2.36,0.62)
estvar <- diag(c(0.22,0.27)^2)
panel_name <- ''
plot_bear(cov,cov_name,param_name,estmean,estvar,panel_name)

## detection x AREA
cov = AREA
cov_name = 'subsection area'
param_name = 'Detection'
estmean <- c(-2.36,-1.4)
estvar <- diag(c(0.22,0.4)^2)
panel_name <- ''
plot_bear(cov,cov_name,param_name,estmean,estvar,panel_name)

dev.off()

