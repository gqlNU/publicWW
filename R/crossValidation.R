###########################################################################
#'  A wrapper function to create the training and the CV datasets for a CV run
#'
#' @return
#' @export
setup_cv_training_datasets <- function(ww.weekly, covars, fold_i, readinCV_sites=TRUE, data_used='30Mar2022') {
	#  getting the CV sites for the 10 CV folds
	if (readinCV_sites) {
		#  pre-generated
		data.file <- system.file("extdata", paste0("cv_sites_10folds_",data_used,".csv"), package = "publicWW")
		sites <- utils::read.csv(data.file,header=TRUE)
	} else {
		#  generate a set (the setting given below gives exactly the same setting as in the package file cv_sites_10folds.csv)
		sites <- setup_CV_sites(seed=1,prop=0.1,file=NULL,data_used,plot.sites=paste0('cv_sites_10folds_',data_used,'.pdf'))
	}
	#  select the CV sites in a specific fold
	fold_sites <- dplyr::filter(sites, fold == fold_i)

	#  labelling observations of each site as cv or train
	cv_df <- get_data_cross_validation(ww.weekly, sitecodes = fold_sites$uwwCode, fold=fold_i)

	if (!is.null(covars)) {
		#  covariates are included
		covar_names <- covars$covariates
		#  forming the training data
		training <- scale_covariates(dplyr::filter(cv_df, cv=='train'), covar_names)
		#  forming the cv data
		cv <- scale_covariates_by_training_set(dplyr::filter(cv_df, cv=='cv'), training, covar_names)
	} else {
		#  no covariates in this model fit
		training <- dplyr::filter(cv_df, cv=='train')
		cv <- dplyr::filter(cv_df, cv=='cv')
	}
	return(list(training=training,cv=cv))
}

###########################################################################
#' set up cross validation with 10% sites removed (use seed=1 to ensure everyone is using the same set of CV sites)
#'
#' @export
setup_CV_sites <- function(seed=1,prop=0.1,file=NULL,data_used='30Mar2022',
                           plot.sites=NULL,analysis_unit=list(time='week',space='uwwCode')) {
  #  format everything based on the 07Mar2022 data so that the CV setting for
	#  the 30Mar2022 runs is the same as that for 07Mar2022 with the extra site goes
  #  in the last CV fold
  dat <- read_public_data(upto='07Mar2022')

	#  transform the wide format data to long
	ww.long <- ww_wide_to_long(dat,analysis_unit$space,show.warning=FALSE)

	#  read in the locations of the sewage wastewater plants
	site_geo <- get_sites_geo(ww.long$site_uniq)

	u <- unique(site_geo[c('uwwCode','eastings','northings')])

	site.codes <- u$uwwCode
	nsites <- length(site.codes)

	set.seed(seed)
	#  number of sites to hide
	ncv.sites <- as.integer(nsites*prop)
	#  number of CV folds
	nfolds <- 1/prop

	cv.ids <- fold <- NULL
	available.ids <- 1:nsites
	for (ifold in 1:nfolds) {
	  if (ifold<nfolds) {
			#  sample ids without replacement from the remaining available ids
			tmp <- sample(available.ids,ncv.sites,replace=FALSE)
			fold <- c(fold,rep(ifold,ncv.sites))
			#  remove ids that have been chosen
			available.ids <- setdiff(available.ids,tmp)
		} else {
			#  put all the remaining sites into the last fold
			tmp <- available.ids
			fold <- c(fold,rep(ifold,length(tmp)))
		}
		cv.ids <- c(cv.ids,tmp)
	}
  if (data_used=='30Mar2022') {
		#  add the one extra site in the 30Mar2022 data compared with the 07Mar2022 dataset
		#  read 30Mar2022 wastewater data in
		dat <- read_public_data('30Mar2022')
		#  transform the wide format data to long
		ww.long <- ww_wide_to_long(dat,analysis_unit$space,show.warning=FALSE)
		ids <- sapply(ww.long$site_uniq$uwwCode,function(x){which(u$uwwCode==x)})
		ll <- unlist(lapply(ids,length))
  	extra.site <- names(ll)[which(ll==0)]
		site.codes <- c(site.codes,extra.site)
		cv.ids <- c(cv.ids,length(site.codes))
		fold <- c(fold,nfolds)
  }
	# summary(as.numeric(table(cv.ids))) # run this to check if each ID appears only once
	out <- data.frame(uwwCode=site.codes[cv.ids],fold=fold)
	if (!is.null(plot.sites)) {
		if (used_data=='30Mar2022') {
 	  	stop('Error: the plot functionality does not work for the 30Mar2022 data yet')
		}
		pdf(file=plot.sites,width=7,height=7)
		for (ifold in 1:nfolds) {
			ii <- sapply(subset(out,out$fold==ifold)$uwwCode,function(x){which(u$uwwCode==x)})
			plot(u$eastings[-c(ii)],u$northings[-c(ii)],pch=19,col='grey',cex=0.5,xlim=c(0,700),ylim=c(0,700))
			points(u$eastings[ii],u$northings[ii],pch=4,col='red',lwd=2)
		}
		dev.off()
	}
	if (!is.null(file)) {
		write.csv(out,file,row.names=FALSE,col.names=TRUE,quote=FALSE)
		return(paste0(getwd(),file))
	} else {
		return(out)
	}
}
#################
#  example:
#################
#setup_CV_sites(seed=1,prop=0.1,file='cv_sites_10folds_07Mar2022.csv',data_used='07Mar2022',plot.sites='cv_sites_10folds_07Mar2022.pdf')
#setup_CV_sites(seed=1,prop=0.1,file='cv_sites_10folds_30Mar2022.csv',data_used='30Mar2022',plot.sites=NULL)

###########################################################################
#'  select sites to be cross validated (not used)
#'
#' @return
#' @export
select_cv_sites <- function(ww.weekly, stratify_column=NULL, cv_size = .2, folds = 1) {

  sites_df <- ww.weekly %>%
    select(region, uwwCode) %>%
    unique()

  df <- tibble()
  for (fold in 1:folds){
    if (is.null(stratify_column)){
      # add random by stratify_column = NULL
      dfnew <- dplyr::slice_sample(sites_df, prop = cv_size)
    }
    else{
      dfnew <- stratified(sites_df, stratify_column, size = cv_size)

    }
    dfnew$fold <- fold
    df <- bind_rows(df, dfnew)
  }

  return(select(df, -region))
}

###########################################################################
#'  Plot CV sites (not used)
#'
#' @return
#' @export
plot_cv_site <- function(site_name, pred.summary, cv, fold, dirname){

  df <- dplyr::filter(pred.summary,uwwCode==site_name)
  data <- dplyr::filter(cv,uwwCode==site_name)

  jpeg(paste0(paste0(dirname,'cv_fold_',fold,'_',site_name,'.png')))

  plot(df$time, df$mean, ylim = c(0,14), type = "l",main=paste0('site ',site_name),xlab='week',ylab='log(flow norm gc/L)')
  plot(data$time_index, data$log_e_gc, ylim = c(0,14),pch=19)

  #make polygon where coordinates start with lower limit and
  # then upper limit in reverse order
  lines(df$time, df$mean, lwd = 1)
  #add red lines on borders of polygon
  lines(df$time, df$q975, col="grey",lty=3, lwd=2)
  lines(df$time, df$q025, col="grey",lty=3, lwd=2)

  title(site_name,
        cex.main = 1,   font.main= 1, col.main= "blue")

  dev.off()
}


###########################################################################
#'  drop sites to be used in cross validation
#'
#' @return
#' @export
get_data_cross_validation <- function(ww.weekly, sitecodes, fold){
  out <- ww.weekly %>%
    mutate(fold = fold,
           cv = case_when(
             !uwwCode %in% sitecodes ~ "train",
             uwwCode %in% sitecodes ~ "cv"
           ))
  return (out)
}

###########################################################################
#'  drop sites to be used in cross validation
#'
#' @return
#' @export
summarise_CV_sims <- function(sims.pred,cv,cvsites,ntimes) {
	pred_time_index <- rep(1:ntimes,each=nrow(cvsites))
	m <- apply(sims.pred,2,mean)
	md <- apply(sims.pred,2,median)
	s <- apply(sims.pred,2,sd)
	q1 <- apply(sims.pred,2,function(x){quantile(x,0.025)})
	q2 <- apply(sims.pred,2,function(x){quantile(x,0.975)})
	# make sure site names being the same order as the one in the arrays
	pred.summary <- data.frame(uwwCode=cv$uwwCode,
														 time=pred_time_index,
														 week_start=unique(cv[c('time_index','day1week')])$day1week[pred_time_index],
														 mean=m,sd=s,median=md,
														 q025=q1,q975=q2)
	#  prediction error
	pred.errors <- calculate_prediction_error(pred.summary, data = cv)
	pred.errors$fold <- fold_i
	return(pred.summary)
}


###########################################################################
#' Given a sims list and raw data observation and prediction, return a dataframe of error metrics
#'
#'
#' @return
#' @export
calculate_prediction_error <- function(pred.summary, data){
	errors <- lapply(unique(data$uwwCode), site_prediction_error, data = data, preds = pred.summary) %>% bind_rows()
	return(errors)
}


###########################################################################
#' Given a site observation and prediction, return a dataframe of error metrics
#'
#' @return
#' @export
site_prediction_error <- function(site_name, data, preds){

	site_preds <- preds %>% dplyr::filter(uwwCode==site_name) %>% rename(mn = mean)
	site_preds$log_e_gc <- data %>% dplyr::filter(uwwCode==site_name) %>% dplyr::select(log_e_gc)


	errors <- site_preds %>% tidyr::drop_na() %>% summarise(uwwCode = first(uwwCode),
	                                                 correlation = cor(mn, log_e_gc, method = 'pearson'),
	                                                 coverage = mean(log_e_gc > q025 & log_e_gc < q975),
	                                                 rmse = sqrt(mean(as.matrix(mn-log_e_gc)^2)))


	return (errors)
}
