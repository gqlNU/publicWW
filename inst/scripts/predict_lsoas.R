library(publicWW)

# directory where the fit result is saved
save.dir <- '/Volumes/WorkSpace/publicWW_results/'

#model_name <- 'ww_model_fit-20220711-191020_allcovars'
model_name <- 'ww_model_fit-20220726-003507'  #  region as fixed effects
model_name <- 'ww_model_fit-20220801-195150'  #  region as random effects

nsims <- 1000  # number of samples to obtain from the joint posterior distribution

####################
#  end user inputs
####################
#  set working directory
setwd(save.dir)

print_with_time(' loading INLA fits ...')
################################################################################
# Model fits are saved from fit_model.R. They include the following objects:
# A.indexs.spde, covars, dm, res, stk, ww.weekly.
fit.file <- paste0(model_name,'.RData')
load(fit.file)
outputs <- list(res,stk,dm,A.indexs.spde) #  gather all objects together for the predict_sites function

print_with_time(' preparing LSOA data for prediction ...')
pred.data <- prepare_lsoa_data(ww.weekly, covars, region_effect)
ntimes <- length(unique(ww.weekly$time_index))

print_with_time(' obtaining weekly LSOA prediction ...')
sims.pred <- predict_sites(outputs, nsims, pred.data$pred.dm, pred.data$pred.geo, ntimes, cv=FALSE, region_effect, region_index=pred.data$region_index)

print_with_time(' saving predictions: iterations ')
#  saving LSOA prediction over iterations
save.file <- paste0(model_name,'_LSOApred_sims.RData')
save(file=save.file,sims.pred,pred.data,region_effect)
