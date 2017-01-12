####### packages
library(rgdal)  
library(unmarked)
library(AICcmodavg)


# read in data 
data_occ <- read.csv('Bear_OccSM_0814_mod.csv', header=FALSE, sep=",")

# ids of sites
site_list <- data_occ[,1]

# detections/non-detections only
data_occ <- data_occ[,2:36]

# number of sites sites
nb_sites <- nrow(data_occ)

# nb of seasons
nb_sec <- 5

# nb of surveyrs
nb_prim <- 7

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

# put all covariates together
sites.covs <- data.frame(RUG=RUG,DTHM=DTHM,CVFR=CVFR,CVBS=CVBS,LGRT=LGRT,AREA = AREA)

# explore correlations
cor(sites.covs)

# format det/non-det and covariates for occupancy analysis with unmarked
simUMF <- unmarkedMultFrame(y= data_occ, siteCovs=sites.covs, numPrimary=nb_prim)
summary(simUMF)

# covariates considered on each parameter (see Table 1):
# detection: roughness, human density, forest cover and area (the latter being always included) 
# colonization: roughness, human density, forest cover
# extinction: human density, road length
# occupancy: roughness, human density, forest cover, shrub cover, road length
# 2^13 = 8192 models in total

# specify combinations to be tested 
ind_det = list('1+AREA','RUG+ AREA','DTHM+ AREA','CVFR+ AREA','RUG+DTHM+ AREA','RUG+CVFR+ AREA','DTHM+CVFR+ AREA','RUG+DTHM+CVFR+ AREA') # + AREA
ind_col = list(1,'RUG','DTHM','CVFR','RUG+DTHM','RUG+CVFR','DTHM+CVFR','RUG+DTHM+CVFR')
ind_ext = list(1,'DTHM','LGRT','DTHM+LGRT')
ind_psi = list(1,'RUG','DTHM','CVFR','CVBS','LGRT',
"RUG+DTHM",
"RUG+CVFR",
"RUG+CVBS",
"RUG+LGRT",
"DTHM+CVFR",
"DTHM+CVBS",
"DTHM+LGRT",
"CVFR+CVBS",
"CVFR+LGRT",
"CVBS+LGRT",
"RUG+DTHM+CVFR",
"RUG+DTHM+CVBS",
"RUG+DTHM+LGRT",
"RUG+CVFR+CVBS",
"RUG+CVFR+LGRT",
"RUG+CVBS+LGRT",
"DTHM+CVFR+CVBS",
"DTHM+CVFR+LGRT",
"DTHM+CVBS+LGRT",
"CVFR+CVBS+LGRT",
"RUG+DTHM+CVFR+CVBS",
"RUG+DTHM+CVFR+LGRT",
"RUG+DTHM+CVBS+LGRT",
"RUG+CVFR+CVBS+LGRT",
"DTHM+CVFR+CVBS+LGRT",
"RUG+DTHM+CVFR+CVBS+LGRT")

# fit all 8192 models - DO NOT RUN! Takes between 8 and 10 hours.
all_aic = list()
all_formula = list()
index = 1
for (i in 1:length(ind_det)){
	for (j in 1:length(ind_col)){
		for (k in 1:length(ind_ext)){
			for (t in 1:length(ind_psi)){
				f_psi <- as.formula(paste('', ind_psi[[t]], sep="~"))
				f_p <- as.formula(paste('', ind_det[[i]], sep="~"))
				f_col <- as.formula(paste('', ind_col[[j]], sep="~"))
				f_ext <- as.formula(paste('', ind_ext[[k]], sep="~"))
				fit <- do.call("colext", list(psiformula=f_psi, gammaformula =f_col, epsilonformula = f_ext, pformula = f_p, data=simUMF, method='BFGS',se=F))
				#all_mod[[index]] <- fit
				all_aic[[index]] <- fit@AIC
				all_formula[[index]] <- fit@formula
				index = index + 1
				print(index)
			}
		}
	}
}
save(all_aic,all_formula,file='modselursus.RData')

# get AICs
load('modselursus.RData')
ord=order(unlist(all_aic))
head(unlist(all_aic)[ord],100)
tail(unlist(all_aic)[ord],100)

# much model uncertainty
# select models with delta-AIC <= 2

# look for models with delta-AIC <= 2
all_aic_ord = unlist(all_aic)[ord]
all_formula_ord = unlist(all_formula)[ord]
mask = (all_aic_ord - min(all_aic_ord)) <= 2

# get Table 3
formulas = all_formula_ord[mask]
as.character(formulas)

all_mod = list()
index = 1
for (i in 1:length(formulas)){
# for (i in 1:30) print(length(as.character(formulas[[i]])))
	temp_f = as.character(formulas[[i]])
	split_tilde = unlist(strsplit(temp_f[2], "~"))
	psi_formula = as.formula(paste('~', split_tilde[2],sep=''))
	gam_formula = as.formula(paste('~', split_tilde[3],sep=''))
	eps_formula = as.formula(paste('~', split_tilde[4],sep=''))
	p_formula = as.formula(paste('~',temp_f[3],sep=''))
	fit = colext(psiformula = psi_formula, gammaformula = gam_formula, 
	epsilonformula = eps_formula, pformula = p_formula, data = simUMF, method = "BFGS", 
    se = T)
    	all_mod[[index]] <- fit
    	index = index + 1
	print(index)
}
save(all_mod,file='allmodursus.RData')

# perform model-averaging (estimates needed for Figure 2)

# psi
modavg(all_mod,modnames=as.character(formulas),parm="(Intercept)",parm.type="psi")
modavg(all_mod,modnames=as.character(formulas),parm="DTHM",parm.type="psi")
modavg(all_mod,modnames=as.character(formulas),parm="LGRT",parm.type="psi")
modavg(all_mod,modnames=as.character(formulas),parm="CVBS",parm.type="psi")
modavg(all_mod,modnames=as.character(formulas),parm="RUG",parm.type="psi")
modavg(all_mod,modnames=as.character(formulas),parm="CVFR",parm.type="psi")

# gamma
modavg(all_mod,modnames=as.character(formulas),parm="(Intercept)",parm.type="gamma")
modavg(all_mod,modnames=as.character(formulas),parm="RUG",parm.type="gamma")
modavg(all_mod,modnames=as.character(formulas),parm="DTHM",parm.type="gamma")
modavg(all_mod,modnames=as.character(formulas),parm="CVFR",parm.type="gamma")

# epsilon
modavg(all_mod,modnames=as.character(formulas),parm="(Intercept)",parm.type="epsilon")
modavg(all_mod,modnames=as.character(formulas),parm="DTHM",parm.type="epsilon")
modavg(all_mod,modnames=as.character(formulas),parm="LGRT",parm.type="epsilon")

# detection
modavg(all_mod,modnames=as.character(formulas),parm="(Intercept)",parm.type="detect")
modavg(all_mod,modnames=as.character(formulas),parm="RUG",parm.type="detect")
modavg(all_mod,modnames=as.character(formulas),parm="DTHM",parm.type="detect")
modavg(all_mod,modnames=as.character(formulas),parm="CVFR",parm.type="detect")
modavg(all_mod,modnames=as.character(formulas),parm="AREA",parm.type="detect")
