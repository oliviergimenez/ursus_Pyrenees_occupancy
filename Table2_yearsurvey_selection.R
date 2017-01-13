# packages
library(unmarked)

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

# nb of surveys
nb_prim <- 7

# year effect
year <- matrix(c('08','09','10','11','12','13','14'), 
	nrow(data_occ), ncol=nb_prim, byrow=T)
year

# survey effect
survey <- matrix(rep(c('1','2','3','4','5'), each=nb_sites), nb_sites, nb_sec*nb_prim)
survey

# format det/non-det and covariates for occupancy analysis with unmarked
simUMF <- unmarkedMultFrame(y= data_occ, yearlySiteCovs=list(year=year),obsCovs=list(survey=survey), numPrimary=nb_prim)
summary(simUMF)

# fit 8 models
mod1 = colext(psiformula = ~ 1, gammaformula =  ~ 1, epsilonformula =  ~ 1, pformula =  ~ 1, data = simUMF, method = "BFGS", se = F)
mod2 = colext(psiformula = ~ 1, gammaformula =  ~ 1, epsilonformula =  ~ year, pformula =  ~ 1, data = simUMF, method = "BFGS", se = F)
mod3 = colext(psiformula = ~ 1, gammaformula =  ~ 1, epsilonformula =  ~ 1, pformula =  ~ survey, data = simUMF, method = "BFGS", se = F)
mod4 = colext(psiformula = ~ 1, gammaformula =  ~ 1, epsilonformula =  ~ year, pformula =  ~ survey, data = simUMF, method = "BFGS", se = F)
mod5 = colext(psiformula = ~ 1, gammaformula =  ~ year, epsilonformula =  ~ 1, pformula =  ~ 1, data = simUMF, method = "BFGS", se = F)
mod6 = colext(psiformula = ~ 1, gammaformula =  ~ year, epsilonformula =  ~ 1, pformula =  ~ survey, data = simUMF, method = "BFGS", se = F)
mod7 = colext(psiformula = ~ 1, gammaformula =  ~ year, epsilonformula =  ~ year, pformula =  ~ 1, data = simUMF, method = "BFGS", se = F)
mod8 = colext(psiformula = ~ 1, gammaformula =  ~ year, epsilonformula =  ~ year, pformula =  ~ survey, data = simUMF, method = "BFGS", se = F)

# get Table 2
models <- fitList('psi(.)gam(.)eps(.)p(.)' = mod1,'psi(.)gam(.)eps(year)p(.)' = mod2,'psi(.)gam(.)eps(.)p(survey)' = mod3,'psi(.)gam(.)eps(year)p(survey)' = mod4,'psi(.)gam(year)eps(.)p(.)' = mod5,'psi(.)gam(year)eps(.)p(survey)' = mod6,'psi(.)gam(year)eps(year)p(.)' = mod7,
'psi(.)gam(year)eps(year)p(survey)' = mod8)ms <- modSel(models) # this is Table 2ms
