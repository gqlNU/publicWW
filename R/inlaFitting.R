###########################################################################
#' fit inla model from data, covariates, and predictor.
#'
#'
#'
#' @return
#' @export
fit_inla <- function(ww.weekly, covars, y, priors, mesh = NULL, region_effect){
  if (!any(region_effect==c('fixed','none','re'))) {
    stop('region_effect must be a string as one of the following: fixed, re or none')
  }

  #  create the design matrix for all covariates
  dm <- create_covariate_model_matrix(ww.weekly, covars, region_effect)

  if (is.null(mesh)){
    mesh <- make_mesh(ww.weekly)
  }

  #  preparation for INLA SPDE
  A.indexs.spde <- spde.make.A(fitdata=ww.weekly,mesh=mesh,geo_columns=c('eastings','northings'),priors=priors)

  #  stacking data
  stk <- form.fitdata.stk(fitdata=ww.weekly, A.indexs=A.indexs.spde, X=dm, y_column='log_e_gc', region_effect)

  #  model specification
  f <- specify_space_time_model(colnames(dm),region_effect)
  formula <- as.formula(f)

  #print(paste0('Fitting INLA formula: ', formula))
  print('Fitting the model below in INLA: ')
  print(formula)

  #  INLA fitting
  t_start <- Sys.time()
  res <- inla(formula,
              data = inla.stack.data(stk),
              family='gaussian',
              # priors on fixed effects and intercept
              control.fixed = priors$fixed,
              control.predictor = list(compute = TRUE,A = inla.stack.A(stk)),
              # need this to sample from posterior later
              control.compute = list(config = TRUE,
                                     # need this to get marginals for the fitted values
                                     return.marginals.predictor=TRUE)
  )
  t_end <- Sys.time()
  print(t_end - t_start)
  return(list(res, stk, dm, A.indexs.spde))
}


###########################################################################
#' Makes model matrix with all covariates. On input, all covariates in dat are already space-time.
#'
#'
#'
#' @return
#' @export
create_covariate_model_matrix <- function(dat, covars, region_effect) {
	if (length(grep('time_index',names(dat)))==0) {
		stop('Error in create_covariate_model_matrix: replicate the covariate matrix over time before entering this function and add the column time_index')
	}
	if (!is.null(covars)) {
    bm <- '~'
    if (region_effect=='fixed') { #  regional fixed effect
      bm <- '~ region + '
    }
		#  with covariates
		cov_model_formula <- as.formula(paste(bm, paste0("scale_",covars$covariates,collapse=" + ")))
		dm <- stats::model.matrix(cov_model_formula,data=dat)[,-1]
		colnames(dm) <- sapply(colnames(dm),function(x){gsub(' ','',x)})
		colnames(dm) <- sapply(colnames(dm),function(x){gsub('\\(','_',x)})
		colnames(dm) <- sapply(colnames(dm),function(x){gsub(')','',x)})
	} else {
		#  without covariates
    dm <- NULL
    if (region_effect=='fixed') {                     #  region as fixed
      cov_model_formula <- as.formula("~ region ")
      dm <- stats::model.matrix(cov_model_formula,data=dat)[,-1]
      colnames(dm) <- sapply(colnames(dm),function(x){gsub(' ','',x)})
    }
	}
  if (FALSE) {
  	if (!is.null(ntimes)){
    	dm <- dm[lsoa.ids,]
  	}
  	if (!is.null(ntimes)){
    	dm <- do.call(rbind,replicate(ntimes,dm,simplify=FALSE))
  	}
	}
  return (dm)
}

###########################################################################
#' make inla mesh
#'
#'
#'
#' @return
#' @export
make_mesh <- function(ww.weekly){
  coords <- get.england.boundary()
  site.geo <- unique(ww.weekly[c('eastings','northings')])

  mesh <- INLA::inla.mesh.2d(
    loc.domain = coords,
    loc = site.geo,
    offset = c(50, 100),
    cutoff = 1, max.edge = c(30, 60)
  )
  return(mesh)
}

###########################################################################
#'
#'
#'
#' @return TRUE/FALSE
#' @export
#'
spde.make.A <- function(fitdata,mesh,geo_columns=c('eastings','northings'),priors) {
	#  get coordinates of sites
	site_geo <- as.matrix(fitdata[,geo_columns])
	#  time index
	if (!any(colnames(fitdata)=='time_index')) stop(' Error: time_index not in the dataset, run the function add_time_site_indices first')
	time_index <- fitdata$time_index
	times <- length(unique(time_index))

	spde <- INLA::inla.spde2.pcmatern(
		mesh = mesh,
		alpha = 2, # alpha is related to the smoothing parameter nu where alpha=nu + D/2 with D=2 (dimensions) for a spatial case
		           # thus alpha=2 => nu is fixed at 1, a special case of Matern correlation function
		           # (see Whittle 1954 Eq. 68 p. 448 and Diggle and Ribeiro 2010, p. 52)
		constr = TRUE,
  	prior.range = priors$matern$range,  # P(range < range[1]) = range[2]
  	prior.sigma = priors$matern$sigma   # P(sigma > sigma[1]) = sigma[2]
	)

	indexs <- INLA::inla.spde.make.index("s",
  		n.spde = spde$n.spde,
  		n.group = times
	)

	A <- INLA::inla.spde.make.A(mesh = mesh, loc = site_geo, group = time_index)
	out <- list(A=A,indices=indexs,spde=spde)
	return(out)
}


###########################################################################
#'  y_column can be 'log_gc' (for the original data)
#'  or 'log_gcX' where X=1,...,nsims, one of the simulated datasets
#'  with LOD entries simulated from U(0,L)
#'
#' @return
#' @export
form.fitdata.stk <- function(fitdata, A.indexs, X=NULL, y_column='log_e_gc', region_effect) {
	###  extracting the outcome column
	if (!is.null(y_column)) {
    outcome <- fitdata[[y_column]]
    tag <- 'est'
  } else {
    outcome <- NA
    tag <- 'pred'
  }

	###  creating a dataframe containing intercept, fixed effects (if any) and indices of random effects
	X_and_REindex <- data.frame(b0 = rep(1, nrow(fitdata)),
                              time.index=fitdata$time_index,
                              site.index=fitdata$site_index)
  #  adding fixed effects/covariates
  if (!is.null(X)) {
  	X_and_REindex <- cbind(X_and_REindex,X)
  }

  #  adding a numeric regional index
  if (region_effect=='re') {
    rgn <- data.frame(region.index=fitdata$region_index)
    X_and_REindex <- cbind(X_and_REindex,rgn)
  }
  ###  stacking everything
	stk <- INLA::inla.stack(tag = tag,
                          data = list(y = outcome),
                          A = list(1, A.indexs$A),
                          effects = list(X_and_REindex,s = A.indexs$indices)
                         )
	return(stk)
}



###########################################################################
#' make hard-coded model.
#'
#'
#'
#' @return
#' @export
specify_space_time_model <- function(fixed_effects = NULL, region_effect){
  #
  rgn_effect <- ''
  if (region_effect=='re') {   # region random effect
    rgn_effect <- ' + f(region.index,model=\'iid\',hyper = hyper.prec.region)'
  }
  fx <- ''
  if (!is.null(fixed_effects)) {
    fx <- paste0(paste(fixed_effects,collapse='+'),'+')
  }
  # fixed effects are a list
  f <- paste('y ~ 0 + b0 + ',                                           # intercept
             fx,                                                        # covariates/fixed effects
             '   f(time.index,model=\'rw1\',hyper = hyper.prec.time)',  # overall time profile
             ' + f(site.index,model=\'iid\',hyper = hyper.prec.space)', # overall site-effects
             rgn_effect,                                                # region random effect
             # spatially-temporally correlated space-time interactions
             ' + f(s, model = A.indexs.spde$spde,
                   group = s.group,
                   control.group = list(model = "ar1", hyper = pcprior.rho()))'
            )
  return(f)
}
