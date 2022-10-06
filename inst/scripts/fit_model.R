library(publicWW)
library(INLA)
library(tibble)
library(dplyr)

#####################
###  user input
#####################
#  specify the directory in which the model fit will be saved
save.dir <- '/Volumes/WorkSpace/publicWW_results/'

#  region is included as fixed ('fixed') or random effect ('re') or not included at all ('none')
region_effect <- 're'

#  LOD threshold
L <- 160

#  include the two genomic covariates as in the paper?
#  this is by default set to FALSE (excluding the genomic covariates) as we do not have permission to make the data public
include_genomic_covars <- FALSE

#  covariates to be included
covars <- tibble(covariates =  c('IMD_score','old_prop','young_prop','bame_proportion','population_density','f_industrial','pct_genome_coverage','num_snp'), #  covariates to be included
                 aggregation= c('wmean','wmean','wmean','wmean','wmean','wmean','none','none') # how to aggregate from LSOA to catchment, taking weighted mean or sum, none means that the covariate is spatio temporal and should note be aggregated.
                )
if (!include_genomic_covars) {
	covars <- tibble(covariates =  c('IMD_score','old_prop','young_prop','bame_proportion','population_density','f_industrial'), #  covariates to be included
                 aggregation= c('wmean','wmean','wmean','wmean','wmean','wmean') # how to aggregate from LSOA to catchment, taking weighted mean or sum, none means that the covariate is spatio temporal and should note be aggregated.
                )
}
analysis_unit <- list(time='week',space='uwwCode')  #  space and time units of analysis (they should be columns in the dataset)

# Defaults to syears = NULL, smonths = NULL. To select specific time regions change these values.
syears  <- NULL #2021  #  which year of data to be modelled
smonths <- NULL #6:9

# construct catchment level covariates from LSOA values?
construct_catchment_covariates <- FALSE 
#####################
###  end user input
#####################

#  read wastewater data in
dat <- read_public_data('30Mar2022')

#  gather and format data for fitting
ww.weekly <- make_weekly_data(dat = dat, analysis_unit = analysis_unit, covars = covars, 
                              LOD = L, syears = syears, smonths = smonths,
                              construct_catchment_covariates=construct_catchment_covariates)

# set prior
priors <- set_priors()
hyper.prec.time <- priors$hyper.prec.time
hyper.prec.space <- priors$hyper.prec.space
if (region_effect=='re') {
  hyper.prec.region <- priors$hyper.prec.region
}

# fit inla
outputs <- fit_inla(ww.weekly=ww.weekly, covars=covars, y='log_e_gc', priors=priors, mesh=NULL, region_effect=region_effect)

res <- outputs[[1]]
stk <- outputs[[2]]
dm <- outputs[[3]]
A.indexs.spde <- outputs[[4]]

summary(res)

########################
###  save results
########################
dir.create(save.dir, showWarnings=FALSE, recursive=TRUE) #creates directory if doesn't exist
filename <- format(Sys.time(), "ww_model_fit-%Y%m%d-%H%M%S.RData")
datafile <- sub('model_fit','data',filename)
# save fit
save(file=paste0(save.dir,filename),res,stk,dm, A.indexs.spde, ww.weekly, covars, region_effect)
# save data
save(file=paste0(save.dir,datafile),ww.weekly)
