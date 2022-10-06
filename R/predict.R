###########################################################################
#' A wrapper function to predict (applied to both CV and LSOA)
#'
#'
#' @return
#' @export
predict_sites <- function(outputs, nsims=1000, pred.dm, pred.geo, ntimes, cv=FALSE, region_effect, region_index) {
  #  extract components from INLA outputs
  res <- outputs[[1]]
  stk <- outputs[[2]]
  dm <- outputs[[3]]
  A.indexs.spde <- outputs[[4]]

  #  extract samples for relevant parameters
  sims.list <- extract.posterior.samples.from.joint(nsims,res,dm,A.indexs.spde,fit.object=TRUE, region_effect)

  #  project the SPDE mesh to the prediction locations
  pred.mesh <- INLA::inla.mesh.projector(A.indexs.spde$spde$mesh,loc=pred.geo)
  print_with_time(' calculating space-time interactions over the predicted locations ... ')

	#  space-time interaction
  eta <- get_eta_interaction(nsims,sims.list,pred.mesh)

	#  unstructured spatial variability
  vsp <- get_between_site_variability(nsims, nrow(pred.geo), ntimes, sims.list)

	#  sampling variability (added for CV only)
	epsilon <- 0
	if (cv) {
		epsilon <- sim_epsilon(nsims, nrow(pred.geo), ntimes, sims.list)
	}

  pred_time_index <- rep(1:ntimes,each=nrow(pred.geo))
  #  putting everything together
  print_with_time(' putting everything together to form site predictions ... ')
  sims.pred <- sims.list$intercept +
               sims.list$fixed%*% t(pred.dm) +
               sims.list$time[,pred_time_index] +
               vsp +
               eta +
							 epsilon

  #  region random effects
  if (region_effect=='re') {
    rm(eta)
    rm(vsp)
    rm(res);gc()
    rgn <- sims.list$region[,region_index]
    sims.pred <- sims.pred + rgn
  }
  return (sims.pred)
}


###########################################################################
#' Obtain posterior samples from the INLA approximate joint posterior
#'
#' Note: The smoothing parameter, nu, in the Matérn covariance function is fixed at 1. (this note is no longer relevant)
#'       This comes from the setting alpha=2 in INLA::inla.spde2.pcmatern (called in spde.make.A in the file 03_run_INLA.R).
#'       See the note there for explanation.
#'
#'
#' @return
#' @export
extract.posterior.samples.from.joint <- function(nsims, fits, dm, A.indexs.spde, fit.object=TRUE, region_effect) {
	if (fit.object) {
		#  obtain random samples from the approximate joint posterior
		samples <- INLA::inla.posterior.sample(n=nsims,fits)
	} else {
		#  the object fits is already the results from inla.posterior.sample
		samples <- fits
	}
	sims.list.joint <- list()
	#  extract samples for intercept
	sims.list.joint$intercept <- INLA::inla.posterior.sample.eval('b0',samples)[1,]
	#  extract samples for the common time trend
	sims.list.joint$time <- t(INLA::inla.posterior.sample.eval('time.index',samples))
	#  extract samples for the common spatial pattern
	sims.list.joint$site <- t(INLA::inla.posterior.sample.eval('site.index',samples))

	#  extract samples for between-site variability (SD)
	sims.list.joint$sigma.site <- sqrt(1/extract.joint.hyperpar(samples,'Precision for site.index'))
  #  extract samples for the variance of the temporal random effects  (SD)
  sims.list.joint$sigma.time <- sqrt(1/extract.joint.hyperpar(samples,'Precision for time.index'))

  if (region_effect=='re') {
    #  extract samples for the regional random effects
	  sims.list.joint$region <- t(INLA::inla.posterior.sample.eval('region.index',samples))
    #  extract samples for the variance of the regional random effects  (SD)
    sims.list.joint$sigma.region <- sqrt(1/extract.joint.hyperpar(samples,'Precision for region.index'))
  }

	#  extract samples for the fixed effects
	sims.list.joint$fixed <- array(0,c(nsims,ncol(dm)))
	for (i in 1:ncol(dm)) {
		sims.list.joint$fixed[,i] <- INLA::inla.posterior.sample.eval(colnames(dm)[i],samples)[1,]
	}

	#  extract samples for error precision
	sims.list.joint$sigma.residual <- sqrt(1/extract.joint.hyperpar(samples,'Precision for the Gaussian observations'))

	#  extract samples for AR1 coefficient in the space-time interactions
	sims.list.joint$rho <- extract.joint.hyperpar(samples,'GroupRho for s')

	###  Matern related
	#  extract samples for the SD in the Matern function
	sims.list.joint$sigma.matern <- extract.joint.hyperpar(samples,'Stdev for s')
	#  extract samples for the range in the Matern function
	sims.list.joint$range <- extract.joint.hyperpar(samples,'Range for s')

	#  extract space-time interactions defined on the vertices of the SPDE mesh
	#  Note that the name of the random effects is called s,
	#  i.e. in the fitting code: f(s, model = A.indexs.spde$spde ...
	cont <- attr(samples, which = ".contents", exact = TRUE)
	id <- which(cont$tag=='s')  #  's' corresponds to the name associated with the space-time interaction in the fitted model: f(s, model = A.indexs.spde$spde ...
	start <- cont$start[id]
	end <- cont$length[id] + start - 1
	ss <- lapply(samples,function(x){x$latent[start:end,1]})
	ss.array <- list2array(ss)
	ntimes <- ncol(sims.list.joint$time)
	s.by.time <- lapply(1:ntimes,function(it) {
		                    extract.1d.sptm.interactions(ss.array,itm=it,isp=NULL,
		                          site.index=A.indexs.spde$indices$s,
		                          time.index=A.indexs.spde$indices$s.group)
		                  })
	sims.list.joint$s.by.time <- list2array(s.by.time)
	return(sims.list.joint)
}


###########################################################################
#' Extract hyperparameters from the posterior draws of the INLA approximate joint posterior
#'
#'
#' @return
#' @export
extract.joint.hyperpar <- function(samples,para) {
	out <- unlist(lapply(samples,function(x){x$hyperpar[para]}))
	return(out)
}



###########################################################################
#' Extract space-time interaction samples for a given site or a given time point
#'
#'
#' @return
#' @export
extract.1d.sptm.interactions <- function(interactions,itm=NULL,isp=NULL,site.index,time.index) {
	if (!is.null(itm) & !is.null(isp)) {
		stop('ERROR in extract.1d.sptm.interactions: This function can only extract one dimension, space or time but not both')
	}
	if (!is.null(itm)) {
		# extracting all terms at time itm
		ids <- which(time.index==itm)
		if (is_consecutive(site.index[ids],unique=FALSE)) {
			out <- interactions[,ids]
		} else {
			stop('ERROR in extract.1d.sptm.interactions: order of site indices is wrong')
		}
	}
	if (!is.null(isp)) {
		# extracting all terms associated with site isp
		ids <- which(site.index==isp)
		if (is_consecutive(time.index[ids],unique=FALSE)) {
			out <- interactions[,ids]
		} else {
			stop('ERROR in extract.1d.sptm.interactions: order of time indices is wrong')
		}
	}
	return(out)
}


###########################################################################
#' Get eta interaction
#'
#'
#' @return
#' @export
get_eta_interaction <-function(nsims,sims.list,pred.mesh) {
 	eta.lsoa.sims <- lapply(1:nsims,function(isim) {
    x <- sims.list$s.by.time[,isim,]
    #  this projects onto a (nlsoas x ntimes) matrix of eta values over the LSOA centroids
    m <- apply(x,1,function(vs){INLA::inla.mesh.project(pred.mesh,vs)})
    #  format the matrix to long in the form of eta_11,..., eta_N1,eta_12,....,eta_N2,...,eta_1T,..., eta_NT
    #  eg c(matrix(1:12,nrow=3,byrow=FALSE))
    c(m)
  })
  eta <- list2array(eta.lsoa.sims)

  return (eta)
}


###########################################################################
#' Get between-site variability (the exchangeable site-level random effects)
#'
#'
#' @return
#' @export
get_between_site_variability <- function(nsims, npred.locs, ntimes, sims.list) {
  vsp <- sapply(1:nsims,function(x){rnorm(npred.locs,0,sims.list$sigma.site[x])}) # per predicted location
  vsp <- do.call(rbind,replicate(ntimes,vsp,simplify=FALSE)) # repeat over ntimes
  vsp <- t(vsp)
  return (vsp)
}

###########################################################################
#' Simulate epsilon, the sampling variability
#'
#'
#' @return
#' @export
sim_epsilon <- function(nsims, npred.locs, ntimes, sims.list) {
  eps <- sapply(1:nsims,function(x){rnorm(npred.locs * ntimes,0,sims.list$sigma.residual[x])})  # per space-time cell
  return (t(eps))
}

###########################################################################
#' Prepare LSOA data for prediction
#'
#'
#'
#' @return
#' @export
prepare_lsoa_data <- function(ww.weekly, covars, region_effect) {
	ntimes <- length(unique(ww.weekly$time_index))
	#  read population-weighted LSOA centroids
  pc <- read_LSOA_centroids(england.only=TRUE)
	#  get spatial only covariates for all LSOAs
  lsoa.X <- get.LSOA.covariates()

  lsoa.covas <- add.region.to.lsoas(lsoa.X)	 #  add region indicator
  nl <- names(lsoa.covas)
  names(lsoa.covas)[grep('LSOA11CD',nl)] <- 'lsoa11cd'
  names(lsoa.covas)[grep('LSOA11NM',nl)] <- 'lsoa11nm'
  uniq.lsoa.covas <- dplyr::distinct(lsoa.covas,lsoa11nm,.keep_all=TRUE)
  all.lsoa <- merge(pc,uniq.lsoa.covas,by=c('lsoa11cd','lsoa11nm'),all.x=TRUE,all.y=TRUE)
  all.lsoa$region <- sapply(all.lsoa$RGN11NM,function(x){gsub(' ','',x)})

  pred.geo <- as.matrix(all.lsoa[c('eastings','northings')])
	nlsoas <- nrow(pred.geo)

	#  replicate the design matrix ntimes times then rbind them so that it can be added
	#  to the space-time interactions that will be created later on
	all.lsoa.sptm <- do.call(rbind,replicate(ntimes,all.lsoa,simplify=FALSE)) #  require some RAM to store this

  all.lsoa.sptm$time_index <- rep(1:ntimes,each=nlsoas)
	if (any(covars$aggregation=='none')) {
		#  add genomic data (temporal only)
		gmX.nms <- covars$covariates[grep('none',covars$aggregation)] # names of time-only covariates
		gmX <- unique(ww.weekly[c('time_index',gmX.nms)])
		gmX.sptm <- gmX[all.lsoa.sptm$time_index,]
		all.lsoa.sptm <- cbind(all.lsoa.sptm,gmX.sptm[gmX.nms])
	}

	#  scale each continuous covariate based on the mean and sd in the fitting data
  all.lsoa.sptm <- scale_covariates_by_training_set(all.lsoa.sptm, ww.weekly, covars$covariates)

	#  create the design matrix for all covariates
	pred.dm  <- create_covariate_model_matrix(all.lsoa.sptm, covars, region_effect)

  #  sptm region index
  region_index <- all.lsoa.sptm$region_index
  names(region_index) <- all.lsoa.sptm$RGN21NM

  out <- list(pred.dm=pred.dm,pred.geo=pred.geo,region_index=region_index)
  return (out)
}
