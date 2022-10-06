###########################################################################
#' Read the publicly available data into R
#t This function checks if the unique values in an array are n consecutive numbers - adopted from https://stackoverflow.com/questions/16118050/how-to-check-if-a-vector-contains-n-consecutive-numbers
#'
#'
#' @return TRUE/FALSE
#' @export
is_consecutive <- function(x,unique=TRUE) {
	if (unique) {
		ux <- unique(x)
		res <- rle(diff(sort(ux),lag=1))
	} else {
		ux <- x
		res <- rle(diff(ux,lag=1))
	}
	#out <- any(res$lengths>=2 & res$values==1)
	out <- !any(res$values!=1)
	return(out)
}
# is_consecutive(1:10)  # return TRUE
# is_consecutive(c(1:10,10,2),unique=TRUE)    # return TRUE
# is_consecutive(c(1:10,10,2),unique=FALSE)   # return FALSE
# is_consecutive(c(1:10,13,2))    # return FALSE


###########################################################################
#' Scale each covariate (x-mean(x))/sd(x) in a list
#'
#'
#'
#' @return
#' @export
scale_covariates <- function(data, covariates){
  # could write our own scale function here
  df <- data %>% mutate(across(.col = covariates, .fns = scale, .names = "scale_{.col}"))
  return(df)
}

###########################################################################
#' Scale the covariates of the prediction sites (CV or LSOA) based on the training data
#'
#' @return
#' @export
scale_covariates_by_training_set <- function(data2scale, training.data ,covariates){
  scale.values <- lapply(covariates,function(x){
    ax <- training.data[[x]]
    out <- c(mean(ax),sd(ax))
    names(out) <- c('mean','sd')
    return(out)
  })
  names(scale.values) <- covariates

  for (i in covariates) {
    sc <- scale.values[[i]]
    data2scale[[paste0('scale_',i)]] <- (data2scale[[i]] - sc['mean']) / sc['sd']
  }
  return(data2scale)
}

###########################################################################
#'  To check if a point (x, y) is within a circle centred at (c1, c2) with radius of r
#'
#' @return
#' @export
inCircle <- function(c1,c2,r,x,y,plot=FALSE) {
	out <- (c1-x)^2+ (c2-y)^2 < r^2
	if (plot) {
		rr <- 1.5
		plot(c1,c2,xlim=c(c1-r*rr,c1+r*rr),ylim=c(c2-r*rr,c2+r*rr),xlab='',ylab='')
		plotrix::draw.circle(c1,c2,r)
		points(x,y,pch=19)
	}
	return(out)
}

###########################################################################
#'  Print a message with current time
#'
#' @return
#' @export
print_with_time <- function(m) print(paste(m,Sys.time(),sep=' @ '))


###########################################################################
#' Convert a list to an array
#'
#'
#' @return
#' @export
list2array <- function(x) {
	l <- length(x)
	d <- dim(x[[1]])
	if (is.null(d)) {
		a <- array(0,c(l,length(x[[1]])))
		for (i in 1:l) {
			a[i,] <- x[[i]]
		}
	} else {
		a <- array(0,c(l,d[1],d[2]))
		for (i in 1:l) {
			a[i,,] <- x[[i]]
		}
	}
	return(a)
}


###########################################################################
#' Calculate the Euclidian distance between a point (centroid) and the locations of STW sites
#'
#'  Examples
#'
#' For a single point,
#'  point2points.distances(c(3,4),matrix(1:6,ncol=2))
#'
#' For multiple points,
#'  t(apply(cbind(rep(3,2),rep(4,2)),1,function(x) {point2points.distances(x,matrix(1:6,ncol=2))}))
#'
#' Run time for all England LSOA centroids takes around 20 minutes
#' n <- 100
#' system.time(
#'   m <- t(apply(cbind(rep(3,n),rep(4,n)),1,function(x) {point2points.distances(x,matrix(1:(2*271),ncol=2))}))
#' )
#'   user  system elapsed
#'  3.038   0.021   3.050
#'
#' @return
#' @export
point2points.distances <- function(coord.pt,coord.pts) {
	coord.pts <- t(coord.pts)
	if (!is.data.frame(coord.pts)) coord.pts <- as.data.frame(coord.pts)
	if (!is.matrix(coord.pt)) coord.pt <- as.matrix(coord.pt)
	ds <- sqrt(colSums((coord.pt - coord.pts)^2))
	return(ds)
}
