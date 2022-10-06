##########################
### running CV
##########################
library(publicWW)
library(INLA)
library(tibble)
library(dplyr)

#####################
###  user input
#####################
start_fold <- 1
end_fold <- 10

#  result directory
save.fit <- TRUE  #  saving INLA fit for each CV run
save.dir  <- '/Volumes/WorkSpace/publicWW_results/CVrun/'  #  specify the directory in which the model fit will be saved

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
syears <- NULL  # 2021  #  which year of data to be modelled
smonths <- NULL #6:9

# construct catchment level covariates from LSOA values?
construct_catchment_covariates <- FALSE

nsims <- 1000  # number of samples drawn from marginal posterior

#####################
###  end user input
#####################

#  read wastewater data in
dat <- read_public_data('30Mar2022')

#  gather and format data for fitting
ww.weekly <- make_weekly_data(dat = dat, analysis_unit = analysis_unit, covars = covars, 
                              LOD = L, syears = syears, smonths = smonths,
                              construct_catchment_covariates=construct_catchment_covariates)

ntimes <- length(unique(ww.weekly$time_index))

# prior specification
priors <- set_priors()
hyper.prec.time <- priors$hyper.prec.time
hyper.prec.space <- priors$hyper.prec.space
if (region_effect=='re') {
  hyper.prec.region <- priors$hyper.prec.region
}

######################################################
# looping through the CV folds (10 in total)
######################################################
for (fold_i in start_fold:end_fold){
  #################
  ####   fit   ####
  #################
  #  forming the cv and the training datasets for this CV run
  dat <- setup_cv_training_datasets(ww.weekly, covars, fold_i, readinCV_sites=TRUE, data_used='30Mar2022')
  training <- dat$training
  cv <- dat$cv

  print_with_time('Fitting cv model...')
  outputs <- fit_inla(ww.weekly=training, covars=covars, y='log_e_gc' , priors=priors, mesh=NULL, region_effect=region_effect)
  #  save model fit if required
  if (save.fit) {
    print_with_time(paste0('Saving cv model fit fold',fold_i,' ...'))
    save.file <- paste0(save.dir,'fit_fold',fold_i,'.RData')
    save(file=save.file,outputs,priors,cv,training)
  }
  #########################
  ####   predictions   ####
  #########################
  print_with_time('Making predictions...')
  #  form the design matrix for the CV sites
  pred.dm  <- create_covariate_model_matrix(cv, covars, region_effect)
  #  locations of prediction
  cvsites <- as.data.frame(unique(cv[, c('eastings', 'northings','uwwCode')]))
  pred.geo <- as.matrix(cvsites[,c('eastings', 'northings')])
  ck <- summary(pred.geo[,'eastings']-cv$eastings[which(cv$time_index==1)])
  if (any(ck>0)) {
    stop('Error: order of prediction locations does not match with the data')
  }
  #  predict for all CV sites
  sims.pred <- predict_sites(outputs, nsims, pred.dm, pred.geo, ntimes, cv=TRUE, region_effect, region_index=cv$region_index)

  #  summarise CV predictions
  summ.pred <- summarise_CV_sims(sims.pred,cv,cvsites,ntimes)

  #  save CV predictions
  print_with_time('Saving predictions...')
  save.file <- paste0(save.dir,'pred_fold',fold_i,'.RData')
  save(file=save.file,summ.pred,sims.pred,region_effect)
}
