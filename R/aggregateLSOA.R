###########################################################################
#' format data for aggregation
#'
#' @description
#' This function carries out two operatoins: (a) constructs the lookup table between LSOA and the geography for which the aggregation is to be carried out; and (b) calculate population weights
#' @return
#' @export
get_lsoa_agg_popw <- function(geography,dat,aggregate.method) {
  #  select lookup between LSOA and the selected geography
	agg <- NULL
	colnms <- c('LSOA11CD','LSOA11NM',sapply(c('CD','NM'),function(x){paste0(geography,'21',x)}),'population')
	agg <- dat[colnms]
	save.file <- paste0(geography,'pred_',aggregate.method,'.RData')

  #  calculate population weights
	pw <- rep(0,length(agg$LSOA11CD))
  acd <- unique(agg[[3]]) #  the geography to aggregate to
	nacds <- length(acd)
	for (i in 1:nacds) {
	  ac <- acd[i]
		ids <- which(agg[[3]]==ac)  # LSOAs within this area
		#  population of the LSOAs within
		pop <- agg$population[ids]
		#  population weights = pop / sum(pop)
		pw[ids] <- pop/sum(pop)
		names(pw)[ids] <- agg$LSOA11CD[ids]
	}
	npw <- names(pw)
	ids <- sapply(agg$LSOA11CD,function(x){which(npw==x)})
	agg$popw <- pw[ids]
	return(list(agg=agg,save.file=save.file))
}

###########################################################################
#' aggregate LSOA concentration to the given geographical level for a given area at a given time point
#'
#' @description
#'
#' @return
#' @export
aggregate_lsoa <- function(area, m,  agg.index, pw=NULL, method='fast') {
	tmp <- rep(0,length(agg.index))
	ids <- which(agg.index==area)
	tmp[ids] <- 1
	if (!is.null(pw)) {
		tmp <- tmp * pw # weighted by population
	}
	if (method=='fast') {
		#  this roughly halves the time of the apply one below
		# it cames from https://stackoverflow.com/questions/51110216/how-to-multiply-each-column-by-each-scalar-in-r
		p <- t(t(m)* tmp)
		aa <- rowSums(p)
	} else {
		aa <- apply(m,1,function(x){sum(x*tmp)})
	}
	return(aa)
}
