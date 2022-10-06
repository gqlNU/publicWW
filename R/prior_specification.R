###########################################################################
#'  Setting PC prior on the AR1 coefficient that Pr(|rho|>a)=pb
#'  By default, a=0.1 and pb=0.9 such that the prior assumes Pr(|rho|>0.1)=0.9
#'
#' @return
#' @export
pcprior.rho <- function(prior='pc.cor0',a=0.1,pb=0.9) {
	pr <- list(theta = list(prior = prior, param = c(a, pb)))
	return(pr)
}

###########################################################################
#'  A function to set prior specification (apart from the AR1 coefficient; for that, see pcprior.rho)
#'
#' @description
#' See https://github.com/gqlNU/publicWW/issues/42 for the discussion of prior specification.
#' The default prior is used for the residual variance, a Gamma prior with 1 and 0.00005 on the Precision
#' (see inla.doc("gaussian") and also https://groups.google.com/g/r-inla-discussion-group/c/MB53BG4WP68).
#'
#' @return
#' @export
set_priors <- function(fixed=list(mean.intercept=0,prec.intercept=0.0001,
                                  mean.fixed=0,prec.fixed=0.0001),
                       prec.overall.time=c(10,0.05),
                       prec.overall.space=c(10,0.05),
											 prec.overall.region=c(10,0.05),
                       matern=list(range=c(10,0.05),
                                   sigma=c(10,0.05))
                      ) {
    #  prior for intercept
    intercept <- list(mean.intercept=fixed$mean.intercept,
                      prec.intercept=fixed$prec.intercept)

    #  prior for fixed effects
    fixed <- list(mean=fixed$mean.fixed,prec=fixed$prec.fixed)

    #  PC prior for precision of overall time assuming Pr(sd>param[1])=param[2]
    time <- list(prec = list(prior="pc.prec", param=prec.overall.time))

    #  PC prior for precision of overall space assuming Pr(sd>param[1])=param[2]
    space <- list(prec = list(prior="pc.prec", param=prec.overall.space))

		#  PC prior for precision of regional random effect assuming Pr(sd>param[1])=param[2]
		region <- list(prec = list(prior="pc.prec", param=prec.overall.region))

    #  PC priors for the hyperparameters in the space-time interactions
    #  for detail, see Thm 2.6 in https://www.tandfonline.com/doi/epub/10.1080/01621459.2017.1415907?needAccess=true
    matern <- list(range=matern$range,  #  prior assumes Pr(range < range[1])=range[2]
                   sigma=matern$sigma)  #  prior assumes Pr(sigma > sigma[1])=sigma[2]

    out <- list(fixed=c(intercept,fixed),
                hyper.prec.time=time,
                hyper.prec.space=space,
								hyper.prec.region=region,
                matern=matern)
    return(out)
}
