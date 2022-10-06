library(publicWW)

##############################
#  user input
##############################
model <- 'ww_model_fit-20220726-003507'   # regional fixed effects
model <- 'ww_model_fit-20220801-195150'   # regional random effects

aggregate.method <- 'simple' # simple sum
aggregate.method <- 'pw'     # population weighted

geography <- 'LAD'  #  local authority district (i.e. LTLA)
geography <- 'RGN'  #  regions
geography <- 'CCG'  #  Clinical commissioning groups
#geography <- 'CTRY'  #  country (i.e. England)

#  working directory
workdir <- '/Volumes/WorkSpace/publicWW_results/'

##############################
#  end user input
##############################

#  set working directory
setwd(workdir)

################################################
#  load lookup tables
################################################
lookup <- get.geo.lookups()
nlsoas <- length(lookup$lookup$LSOA11CD)

################################################
#  extract predictions
#  Note: this does not work with batch run
#        of LSOA prediction yet
################################################
print_with_time(' loading LSOA prediction ... ')
file <- paste0(model,'_LSOApred_sims.RData')
load(file)

nsims <-  nrow(sims.pred)
ntimes <- ncol(sims.pred)/nlsoas
if (ntimes!=floor(ntimes)) {
	stop('Error: ntimes is not an integer')
}

################################################
#  get LSOA to LAD/region/England lookup
#  this is where other LSOA-to-other geography
#  come in (e.g. LSOA-to-CCG, LSOA-to-LTLA)
################################################
#  get covariates and geo indicators (e.g. LAD/CCG) for all LSOAs
all.lsoa <- get.LSOA.covariates()
names(all.lsoa)[grep('lsoa11cd',names(all.lsoa),ignore.case=TRUE)] <- 'LSOA11CD'
names(all.lsoa)[grep('lsoa11nm',names(all.lsoa),ignore.case=TRUE)] <- 'LSOA11NM'

#  identify the columns of area indicators over which LSOA predictions are aggregated to
#   + population weights
agg <- get_lsoa_agg_popw(geography,all.lsoa,aggregate.method)
save.file <- paste0(model,'_',agg$save.file)
print_with_time(save.file)

agg.sp <- agg$agg
code <- paste0(geography,'21CD')
print_with_time(paste0(' start constructing predictions over time for ',geography))
uareas <- unique(agg.sp[[code]])
nareas <- length(uareas)

################################################
#  aggregating the LSOA predictions
################################################
refresh <- c(50,10,1,1)  # number of areas to print update
names(refresh) <- c('LAD','CCG','RGN','CTRY')

#  population weighted or not
popw <- agg$agg$popw
if (aggregate.method=='simple') {
  popw <- NULL
}

#  creating storage array
agg.sims <- array(0,c(nsims,nareas,ntimes))

time_index <- rep(1:ntimes,each=nlsoas)

#  aggregation starts
for (tt in 1:ntimes) {
	time.ids <- which(time_index==tt)
	sm <- sims.pred[,time.ids]
	for (iarea in 1:nareas) {
		agg.sims[,iarea,tt] <- aggregate_lsoa(area=uareas[iarea],m=sm,agg.index=agg.sp[[code]],pw=popw)
		if (iarea%%refresh[geography]==0) {
		  msg <- paste0(' t = ',tt, ' (out of ',ntimes,' weeks), ',iarea,' out of ',nareas,' ',geography,'s done ...')
		  print_with_time(msg)
		}
	}
}
dimnames(agg.sims)[[2]] <- uareas

#  save aggregation
save(file=save.file,agg.sims)
