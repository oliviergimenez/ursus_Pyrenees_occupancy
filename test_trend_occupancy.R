# packages
library(rgdal)  ## pour occupancy
#library(RColorBrewer)  ## pour colorer maps
library(unmarked)
#library(classInt)
library(spdep)


# get 30 models with delta-AIC < 2 (see Table3_covariate_selection.R)
load('allmodursus.RData')
all_mod

# get covariates
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

# test trend using each of the 30 models
result = NULL
for (i in 1:30){ # loop on models
	
# psi/Occupancy 
labels <- names(all_mod[[i]]['psi']@estimates)[-1] # drop intercept
nb_cov <- length(labels)
logit_intercept <- all_mod[[i]]['psi']@estimates[1]
logit_slopes <- all_mod[[i]]['psi']@estimates[2:(2+nb_cov-1)] * t(sapply(labels, function(x) eval(parse(text=x))))
logit_par = logit_intercept + apply(logit_slopes,2,sum)
data_psi <- 1/(1+exp(-logit_par))

# gam/Colonization 
labels <- names(all_mod[[i]]['col']@estimates)[-1] # drop intercept
nb_cov <- length(labels)
logit_intercept <- all_mod[[i]]['col']@estimates[1]
logit_slopes <- all_mod[[i]]['col']@estimates[2:(2+nb_cov-1)] * t(sapply(labels, function(x) eval(parse(text=x))))
logit_par = logit_intercept + apply(logit_slopes,2,sum)
data_gamma <- 1/(1+exp(-logit_par))

# eps/Extinction 
labels <- names(all_mod[[i]]['ext']@estimates)[-1] # drop intercept
nb_cov <- length(labels)
logit_intercept <- all_mod[[i]]['ext']@estimates[1]
logit_slopes <- all_mod[[i]]['ext']@estimates[2:(2+nb_cov-1)] * t(sapply(labels, function(x) eval(parse(text=x))))
logit_par = logit_intercept + apply(logit_slopes,2,sum)
data_eps <- 1/(1+exp(-logit_par))

psi <- gamma <- eps <- rep(NA, dim(data_cov)[1])
a <- 1
for (j in 1:dim(data_cov)[1]){
	if (test_suivi[j] == 1){
		psi[j] <- data_psi[a]
		gamma[j] <- data_gamma[a]
		eps[j] <- data_eps[a]
		a <- a+1
	}	
}

# get parameters
data_cov$psi <- psi
data_cov$psi2 <- psi * (1-eps) + (1-psi) * gamma
data_cov$psi3 <- data_cov$psi2 * (1-eps) + (1-data_cov$psi2) * gamma
data_cov$psi4 <- data_cov$psi3 * (1-eps) + (1-data_cov$psi3) * gamma
data_cov$psi5 <- data_cov$psi4 * (1-eps) + (1-data_cov$psi4) * gamma
data_cov$psi6 <- data_cov$psi5 * (1-eps) + (1-data_cov$psi5) * gamma
data_cov$psi7 <- data_cov$psi6 * (1-eps) + (1-data_cov$psi6) * gamma
data_cov$gamma <- gamma
data_cov$eps <- eps

# dataset with occ prob estimates by subsections/years
data_trend = data.frame(subsection=rep(cov$Numero,7),occupancy=c(data_cov$psi,data_cov$psi2,data_cov$psi3,data_cov$psi4,data_cov$psi5,data_cov$psi6,data_cov$psi7),year=c(rep(2008,138),rep(2009,138),rep(2010,138),rep(2011,138),rep(2012,138),rep(2013,138),rep(2014,138)))

# posterior distributions oflatent occurrence) using empirical Bayes methods 
re <- ranef(all_mod[[i]]) 

# stores the estimated posterior distributions of the latent occurrence
# "safer to use the posterior mean even though this will not be an integer in general", see ?bup
res = bup(re, stat="mean") # posterior mean

# build adjancency matrix
nb.r = poly2nb(sousmassif.rg, queen=F) # two areas are neighbors if share common edges with length >0
mat = nb2mat(nb.r, style="B") # mat is the 0/1 adjacency matrix
mat2 = matrix(as.numeric(mat),nrow=138,byrow=T) 
mat3 = mat2[!is.na(data_cov$psi),!is.na(data_cov$psi)] 

# test time effect
library(spaMM) # help(spaMM)
stat = fixedLRT(null.formula=occupancy~1+adjacency(1|subsection),
       formula=occupancy~year+adjacency(1|subsection), adjMatrix=mat3,family=gaussian(),
       HLmethod='ML',data= data_trend)
detach(package:spaMM) # spaMM overrides unmarked, and generates conflicts in using the raneff function
result = rbind(result,stat$basicLRT)
}
result

        LR2 df       pvalue
1  25.65404  1 4.084385e-07
2  23.56641  1 1.206760e-06
3  26.44434  1 2.712417e-07
4  27.87521  1 1.293973e-07
5  35.46650  1 2.594754e-09
6  76.49419  1 0.000000e+00
7  23.52359  1 1.233915e-06
8  27.85504  1 1.307531e-07
9  26.48704  1 2.653127e-07
10 35.35988  1 2.740745e-09
11 21.27237  1 3.984325e-06
12 18.10081  1 2.095128e-05
13 77.09445  1 0.000000e+00
14 36.11472  1 1.860360e-09
15 25.91195  1 3.573507e-07
16 98.22300  1 0.000000e+00
17 87.23120  1 0.000000e+00
18 34.38891  1 4.512878e-09
19 31.33306  1 2.173469e-08
20 23.73614  1 1.104885e-06
21 16.54959  1 4.739426e-05
22 72.11309  1 0.000000e+00
23 21.98990  1 2.740893e-06
24 21.14715  1 4.253332e-06
25 18.20122  1 1.987510e-05
26 22.10090  1 2.586884e-06
27 24.05342  1 9.369987e-07
28 37.54229  1 8.945223e-10
29 97.99203  1 0.000000e+00
30 37.92658  1 7.345774e-10
