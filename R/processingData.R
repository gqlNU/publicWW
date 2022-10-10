###########################################################################
#' Compose helper functions to make weekly data from public WW data
#'
#'
#'
#' @return
#' @export
make_weekly_data <- function(dat, analysis_unit, covars, 
                             LOD, syears = NULL, smonths = NULL,
                             construct_catchment_covariates) {

  #  transform the wide format data to long
  ww.long <- ww_wide_to_long(dat,analysis_unit$space,show.warning=TRUE)

  #  read in the locations of the sewage wastewater plants
  site_geo <- get_sites_geo(ww.long$site_uniq)

  #  combine the daily-site measurements with site locations
  ww.with.geo <-  merge_ww_with_site_geo(ww.long$ww.long,site_geo,by=analysis_unit$space)

  #  subset data by year and months if required
  ww <- select_region_time(ww.with.geo,sel.year=syears,sel.month=smonths)

  #  create continuous numeric time and site indices
  ww <- add_time_site_indices(ww,analysis_unit)

  #  replace LOD by L/2 and create an LOD indicator
  ww <- replace_LOD(ww,raw_gc_column='gc',LOD.means=LOD/2,LOD.sds=NULL,dist='Uniform',isim_LOD=NULL,LOD_holder='tLOD')

  #  aggregated to weekly-site data
  ww.weekly <- aggregate_time(ww,analysis_unit,outcome_col='gc',log_e_y=TRUE)

  if (!is.null(covars)) {
  	#  getting catchment-level time-invariant covariates (excluding the genomic variables)
	X_names <- covars$covariates[covars$aggregation!='none']
	op <- covars$aggregation[covars$aggregation!='none']
  	if (construct_catchment_covariates) {
  		#  get catchment lookup for all sites (including the missing ones in South West Water)
		catchments_lookup <- get_catchments_lookup(ww.weekly)
  		#  aggregate LSOA-level values to catchment using the catchment-LSOA lookup
    	covariates.data <- LSOA.to.catchment.covariates(X_names,op,catchments_lookup)
  	} else {
  		#  load pre-constructed catchment-level covariates (no genomic)
  		data_file <- system.file("extdata","catchment_covariates_no_genomic.rds",package = "publicWW")
  		covariates.data <- readRDS(data_file)
  	}
  	ww.weekly <- merge_covariates_ww_data(ww.weekly,covariates.data,X_names)

    #  genomic data
    covars_time_dependent <- covars[covars$aggregation=='none',]
    if (nrow(covars_time_dependent)!=0){
        genomics.weekly <- get_weekly_genomics_data(ww, syears, smonths, analysis_unit, 
                                                    covars_time_dependent$covariates)
        ww.weekly <- dplyr::left_join(ww.weekly, genomics.weekly, by = c("time_index"),all.x=TRUE)
    }
  }

  #  categorical regional effects
  ww.weekly$region <- as.factor(ww.weekly$RegionName)

  #  a numeric index for region (used when including region as random effects)
  ww.weekly$region_index <- add_region_index(ww.weekly$RegionName)

  if (!is.null(covars)) {
    # scale covariates
    ww.weekly <- scale_covariates(ww.weekly, covars$covariates)
  }

  return(ww.weekly)
}


###########################################################################
#' Transform the raw public wastewater data from wide to long format
#'
#'  week is the week of the year and is obtained using the ISOweek function
#'  in the package of the same name.
#'
#'  For 01/01/2022 (a Saturday), it is allocated to week 52 in 2021.
#'  However, day, indicating the day of the year, is 1 and month, the month of the year, is 1.
#'  If the day (or the month) column is
#'  used as time_index for the analysis, care needs to be taken!
#'  For different date definitions, see https://stackoverflow.com/questions/45549449/transform-year-week-to-date-object/45587644#45587644
#'
#' @param dat The dataset from the function read_public_data()
#' @return A data matrix
#' @export
ww_wide_to_long <- function(dat,site_column,show.warning=TRUE) {
	###  extract time invariant variables
	#usites <- unique(dat[c('uwwCode','uwwName','Population')])
	usites <- dat[c('RegionName','uwwCode','uwwName','Population')]

	###  extract measurements
	ids <- grep('/',colnames(dat))
	meas <- t(dat[ids])  #  transpose so rows=dates and columns=sites
	colnames(meas) <- usites$uwwCode
	meas <- as.data.frame(meas)
	meas$Date <- rownames(meas)

	###  wide to long via gather
	###    stacking the sites so that
	###    gc_11,...,gc_1T,gc_21,...,gc_2T,...,gc_ST (S=site, T=time)
  ### meas.long <- tidyr::gather(meas,site_column,'gc',-Date,na.rm=FALSE)
  ## use pivot_longer instead of gather because gather is practically deprecated and
  ## with this you can remove one "no visible binding" warning
  meas.long <- tidyr::pivot_longer(meas, 1:(ncol(meas)-1), names_to="site_column", values_to="gc")
  ## pivot_longer returns in a different order - sort by site name here to keep the sites
  ## together but the site ordering is different because the source order isn't alphabetical
  meas.long = meas.long[order(meas.long$site_column),]
  ## pivot_longer insists on returning a tibble so lets make it a data frame
  ## really `cdata` is a much better framework for reshaping data...
  meas.long = data.frame(meas.long)
  names(meas.long)[which(names(meas.long)=='site_column')] <- site_column

	###  getting temporal variables
	d <- as.character(strptime(as.character(meas.long$Date),"%d/%m/%Y"))
	yw <- ISOweek::ISOweek(as.Date(d))
	tyw <- (strsplit(yw,'-W'))
	meas.long$year <- as.numeric(unlist(lapply(tyw,function(x){x[1]})))  # year (all the same as those from lubridate::year apart from the first couple of days in 2022 which are assigned to 2021 using ISOweek::ISOweek)
	meas.long$week <- as.numeric(unlist(lapply(tyw,function(x){x[2]})))  # week of the year (this gives the same output from lubridate::isoweek)
	d <- strptime(as.character(meas.long$Date),"%d/%m/%Y")
	meas.long$month <- lubridate::month(d)  # month of the year
	meas.long$day <- lubridate::yday(d)     # day of the year
	if (FALSE) {
		#  previously used
		d <- strptime(as.character(meas.long$Date),"%d/%m/%Y")
		meas.long$year <- lubridate::year(d)
		meas.long$week <- lubridate::isoweek(meas.long$Date) # week of the year
	}
	if (show.warning) {
		warning(' It is fine to use week to construct time_index but care needs to be taken when using day or month to construct time_index (see ?ww_wide_to_long for detail) ')
	}
	return(list(ww.long=meas.long,site_uniq=usites))
}


###########################################################################
#' Merging site geo data with wastewater data
#'
#'
#'
#' @return
#' @export
merge_ww_with_site_geo <- function(ww.long.data,site_geo,by) {
	ww <- merge(ww.long.data,site_geo,by=by,all.x=TRUE,sort=FALSE)
	return(ww)
}


###########################################################################
#'
#' @return
#' @export
select_region_time <- function(dat,sel.year=NULL,sel.region=NULL,sel.month=NULL,sel.week=NULL) {
	###  select region(s) (if required)
	if (!is.null(sel.region)) {
        nr <- length(sel.region)
        temp <- NULL
        for (ir in 1:nr){
            temp <- rbind(temp,dat[which(dat$RegionName==sel.region[ir]),])
        }
        dat <- temp
    }
    ###  select year(s)
	if (!is.null(sel.year)) {
		ny <- length(sel.year)
		temp <- NULL
		for (iy in 1:ny){
        temp <- rbind(temp,dat[which(dat$year==sel.year[iy],),])
    }
		dat <- temp
	}
	###  select months (if required)
	if (!is.null(sel.month)) {
		nm <- length(sel.month)
		temp <- NULL
		for (im in 1:nm) {
			if (sel.month[im]==6) weeks <- c(22,26)  ## WARNING: for 2021 only!!!
			if (sel.month[im]==7) weeks <- c(27,30)
			if (sel.month[im]==8) weeks <- c(31,34)
			if (sel.month[im]==9) weeks <- c(35,39)
			temp <- rbind(temp,dat[which(dat$week>=weeks[1] & dat$week<=weeks[2]),])
		}
		dat <- temp
	}
	###  select weeks (if required)
	if (!is.null(sel.week)) {
		nw <- length(sel.week)
		temp <- NULL
		for (iw in 1:nw){
        temp <- rbind(temp,dat[which(dat$week==sel.week[iw]),])
    }
		dat <- temp
	}
	return(dat)
}


###########################################################################
#'
#' @return
#' @export
add_time_site_indices <- function(ww.agg,analysis_unit) {
  time_index_exists <- any(colnames(ww.agg)=='time_index')
  site_index_exists <- any(colnames(ww.agg)=='site_index')
  if (time_index_exists | site_index_exists) stop(' Error: time_index and/or site_index already exists in the dataset')
  #  create time index
  if (analysis_unit$time!='week') {
    stop('Error from add_time_site_indices: currently only week can be used to define time_index')
  }
  time_index <- create_consecutive_time_index(ww.agg)
  if (!is_consecutive(time_index,unique=TRUE)) stop('Error: time index used is not consecutive.')
  ww.agg$time_index <- time_index
  ww.agg$time_index_yw <- names(time_index)
  #  create numeric site index
  ww.agg$site_index <- as.numeric(as.factor(ww.agg[[analysis_unit$space]]))
  return(ww.agg)
}


###########################################################################
#' Construct time index
#'
#' @return
#' @export
create_consecutive_time_index <- function(dat) {
	yw <- paste(dat$year,dat$week,sep='w')
	years <- unique(dat$year)
	nyears <- length(years)
	wks <- lapply(years, function(x){unique(dat$week[which(dat$year==x)])})
	nwks <- unlist(lapply(wks,length))
	names(wks) <- names(nwks) <- years
	#  unique time indices
	inx <- 1:nwks[1]
	nms <- paste(years[1],wks[[1]],sep='w')
	if (nyears>1) {
		for (iy in 2:nyears) {
			inx <- c(inx,inx[length(inx)] + 1:nwks[iy])
			nms <- c(nms,paste(years[iy],wks[[iy]],sep='w'))
		}
	}
	names(inx) <- nms

	inx.time <- rep(0,length(dat$year))
	for (i in 1:length(inx)) {
		ids <- which(yw==names(inx)[i])
		inx.time[ids] <- inx[i]
		names(inx.time)[ids] <- names(inx)[i]
	}
	return(inx.time)
}




##########################################################################
#'  todo: L differs by type of measurement technique. Need to accommodate that
#'
#'  If you want to replace all LODs with a fixed value (L/2 say), set
#'    LOD.means <- L/2 and leave the rest as default
#'
#'  If you want to replace each LOD with a value randomly generated from U(0,L), set
#'    LOD.means <- L, LOD.sds <- 0 and leave dist as default
#' @return
#' @export
replace_LOD <- function(ww,raw_gc_column='gc',LOD.means,LOD.sds=NULL,dist='Uniform',isim_LOD=NULL,LOD_holder='tLOD') {
	if (!any(colnames(ww)==paste0(raw_gc_column,'_raw'))) {
		###  create the raw gc column with NA, tLOD and other observed values
		ww$gc_raw <- ww[[raw_gc_column]]
		###  add LOD indicator
		ww$lod <- rep(0,nrow(ww))
		ww$lod[which(ww$gc_raw==LOD_holder)] <- 1
	}
	temp <- ww[[raw_gc_column]]
	temp[which(ww$lod==1)] <- -99
	temp <- as.numeric(temp)

	lod.ids <- which(ww$lod==1)
	#  a fixed value replacement for LOD (no uncertainty)
	if (length(LOD.means)==1 & is.null(LOD.sds)) temp[lod.ids] <- LOD.means

	#
	if (length(LOD.means)==1 & !is.null(LOD.sds)) {
		if (dist=='Uniform') temp[lod.ids] <- stats::runif(length(lod.ids),0,LOD.means)
	}

	if (any(temp[which(!is.na(temp))]<0)) stop(' Error: there are still LOD entried not been replaced')

	if (is.null(isim_LOD)) {
		ww[[raw_gc_column]] <- temp
		msg <- paste0(length(lod.ids),' ',LOD_holder,' in column ',raw_gc_column,' have been replaced')
		cat(msg)
	} else {
		temp <- data.frame(a=temp)
		names(temp) <- paste0('gc_sim',isim_LOD)
		ww <- cbind(ww,temp)
		msg <- paste0(' A new column ',names(temp),' has been added to the dataset')
		cat(msg)
	}
	return(ww)
}


###########################################################################
#'
#' @return
#' @export
aggregate_time <- function(ww.dat,analysis_unit,outcome_col='gc',log_e_y=TRUE) {
	obs <- ww.dat[[outcome_col]]
	if (log_e_y) obs <- log(obs)
	tm <- ww.dat$time_index
	site <- ww.dat$site_index
	###  mean gc by week/month by site
	if (analysis_unit$time=='week') {
		agg <- stats::aggregate(obs,list(site_index=site,time_index=tm),function(x){mean(x,na.rm=TRUE)})
		###  add week index
		agg <- add_index(ww.dat,agg,by='time_index',index_to_add='week',type='numeric')
		###  add month index
		agg <- add_index(ww.dat,agg,by='time_index',index_to_add='month',type='numeric')
	}
	#if (analysis_unit$time=='month') {
	#	agg <- stats::aggregate(obs,list(site=site,month=tm),function(x){mean(x,na.rm=TRUE)})
	#	###  add month index
	#	agg <- add_index(ww.dat,agg,by='week',index_to_add='month',type='numeric')
	#}
	if (log_e_y) outcome_col <- paste0('log_e_',outcome_col)
	names(agg)[which(names(agg)=='x')] <- outcome_col
	agg[[outcome_col]][which(is.nan(agg[[outcome_col]]))] <- NA
	###  add year index
	agg <- add_index(ww.dat,agg,by='month',index_to_add='year',type='numeric')
	###  add first day of the week
	agg <- add_index(ww.dat,agg,by='week',index_to_add='day1week',type='character')
	###  add region index
	agg <- add_index(ww.dat,agg,by= 'site_index',index_to_add='RegionName',type='character')
	###  add eastings and northings
	agg <- add_index(ww.dat,agg,by= 'site_index',index_to_add='eastings',type='numeric')
	agg <- add_index(ww.dat,agg,by= 'site_index',index_to_add='northings',type='numeric')
	agg <- add_index(ww.dat,agg,by= 'site_index',index_to_add='uwwLatitude',type='numeric')
	agg <- add_index(ww.dat,agg,by= 'site_index',index_to_add='uwwLongitude',type='numeric')
	###  add site names and site code
	agg <- add_index(ww.dat,agg,by= 'site_index',index_to_add='uwwCode',type='character')
	agg <- add_index(ww.dat,agg,by= 'site_index',index_to_add='uwwName',type='character')
	return(agg)
}


###########################################################################
#'
#' @return
#' @export
add_index <- function(long_dat,agg_dat,by,index_to_add,type='character') {
	###  used in the aggregated.time function below
	n <- nrow(agg_dat)
	uby <- unique(agg_dat[[by]])
	nuby <- length(uby)
	if (type=='numeric') temp <- rep(0,n)
	if (type=='character') temp <- rep('',n)
	for (i in 1:nuby) {
		ids_agg <- which(agg_dat[[by]]==uby[i])
		if (index_to_add=='day1week') {
			temp[ids_agg] <- long_dat[['Date']][which(long_dat[[by]]==uby[i])][1]
		} else {
			temp[ids_agg] <- long_dat[[index_to_add]][which(long_dat[[by]]==uby[i])][1]
		}
	}
	agg_dat <- cbind(agg_dat,temp)
	names(agg_dat)[which(names(agg_dat)=='temp')] <- index_to_add
	return(agg_dat)
}


###########################################################################
#' get approximate catchment for all sites (including missing)
#'
#'
#'
#' @return
#' @export
#'
get_catchments_lookup <- function(ww.with.geo){
  lookup <- LSOA.to.catchment.lookup()  #  all but no South West
  #  identify STW catchments (from South West Water) that are not in the catchment-LSOA lookup table
  usites <- unique(ww.with.geo[c('uwwName','uwwCode','eastings','northings')])
  ids <- sapply(usites$uwwCode,function(x){which(lookup$uwwCode==x)})
  l <- unlist(lapply(ids,length))
  i <- which(l==0)
  #  the sites for which we need to define catchments
  st <- usites[i,]
  #  LSOA centroids
  pc <- read_LSOA_centroids()
  #  approximate the catchments using circles of varying radius
  out <- approximate.catchment.using.circle(st,pc,initial.radius=1,nlsoas.req=10,max.nsearches=40,r.inc=0.2,plot=TRUE)
  approx.catchments <- out$approx.catchments

  lookup <- rbind(lookup,approx.catchments)

  return(lookup)
}


###########################################################################
#'  Approximating catchments using circles with varying radius to achieve a required number of LSOA centroids to fall within
#'
#' @return
#' @export
approximate.catchment.using.circle <- function(site.coords,pc,initial.radius=1,nlsoas.req=10,max.nsearches=40,r.inc=0.2,plot=FALSE) {
	nusites <- nrow(site.coords)
	#  initialise
	lsoa.in <- as.list(1:nusites)
	lsoa.covered <-
	  R <-
	  uww.nm <-
	  uww.cd <- NULL
	met <- not.met <- 0
	if (plot) {
		plot(site.coords$eastings,site.coords$northings,
		     asp=1,pch=19,col=2,cex=0.2,xlab='Eastings',ylab='Northings')#,xlim=c(75,400),ylim=c(0,200))
		points(pc$eastings,pc$northings,pch=19,col='light grey',cex=0.2)
	}
	for (j in 1:nusites) {
		ic <- 0
		r <- initial.radius
		a <- unlist(site.coords[j,c('eastings','northings')])
		#  which LSOA centroids are in the circle
		tmp <- apply(as.matrix(pc[,c('eastings','northings')]),1,function(x){inCircle(a[1], a[2], r, x[1], x[2])})
		#  number of LSOA centroids in the circle
		ntmp <- length(which(tmp))
		while (ic<=max.nsearches & ntmp<= nlsoas.req) {
			#  increase the circle radius
			r <- r + r.inc
			tmp <- apply(as.matrix(pc[,c('eastings','northings')]),1,function(x){inCircle(a[1], a[2], r, x[1], x[2])})
			ntmp <- length(which(tmp))
			ic <- ic + 1  #  increase search counter
		}
		if (ntmp>=nlsoas.req) {
			met <- met + 1
		} else {
			not.met <- not.met + 1
		}
		R <- c(R,r)  #  store final radius of circle
		uww.nm <- c(uww.nm,rep(site.coords[j,'uwwName'],ntmp))  # STW site names for these LSOAs
		uww.cd <- c(uww.cd,rep(site.coords[j,'uwwCode'],ntmp))  # STW site codes for these LSOAs
		#  store LSOA within the final circle
		lsoa.in[[j]] <- pc[which(tmp),]
		lsoa.covered <- rbind(lsoa.covered,pc[which(tmp),])
		if (plot) {
			plotrix::draw.circle(a[1], a[2],r,col='transparent')
			points(lsoa.in[[j]]$eastings,lsoa.in[[j]]$northings,pch=19,cex=0.1,col=1)
			if (ntmp>=nlsoas.req) points(a[1], a[2],col=3,cex=0.3)
			if (j==nusites) {
				legend('topleft',legend=c(paste0(met,' STWs with >=',nlsoas.req,' LSOAs'),
				                          paste0(not.met,' STWs with <',nlsoas.req,' LSOAs'),
				                          'LSOAs covered',
				                          'LSOAs not covered','approx. catchment boundary'),
				                  pch=c(rep(19,4),-1),col=c(3,2,1,'light grey',1),
				                  bty='n',cex=0.7,lty=c(rep(0,4),1)
				       )
			}
		}
	}
	#  outputs
	out <- list()
	out$R <- R  #  radius at which the required number of LSOAs is achieved or the max number of searches is reached
	out$lsoa.in <- lsoa.in  #  detail on the LSOAs included
	out$lsoa.covered <- lsoa.covered  #  all LSOAs covered by all the sites
	#  in the format of the output from LSOA.to.catchment.lookup()
	out$approx.catchments <- data.frame(identifier=rep('circular',nrow(out$lsoa.covered)),
	                                    LSOA11CD=lsoa.covered$lsoa11cd,
	                                    intersection_area=rep(0,nrow(out$lsoa.covered)),
	                                    name=uww.nm,
	                                    uwwCode=uww.cd,
	                                    uwwName=uww.nm,
	                                    distance=rep(0,nrow(out$lsoa.covered))
	                                   )
	#  output messages
	m1 <- range(unlist(lapply(lsoa.in,nrow)))
	m1 <- paste0('**  Range of LSOAs covered is between ',m1[1],' (min) and ',m1[2],' (max)')

	m2a <- table(as.numeric(table(lsoa.covered$lsoa11cd)))
	m2b <- nrow(lsoa.covered)
	if (m2a==m2b) {
		m2 <- paste0('**  Each of the ',m2a,' LSOAs only appears in at most one circle')
	} else {
		m2 <- paste0('**  Some of the ',m2a,' LSOAs appear in two or more circles')
	}
	if (plot) {
		legend('bottomright',legend=c(m1,m2),bty='n',cex=0.6,pch=rep(-1,2))
	} else {
		print(m1)
		print(m2)
	}
	return(out)
}


###########################################################################
#'  get catchment level covariates from LSOA
#'
#' @return
#' @export
LSOA.to.catchment.covariates <- function(X_names=c('IMD_score','all_ages'),op=c('mean','sum')
                                        ,lsoa2catchment.lookup){
  lsoa.covas <- get.LSOA.covariates()
  lookup <- lsoa2catchment.lookup  #  this includes the approximated catchments
  #  merging the two
  all <- merge(lsoa.covas,lookup,by='LSOA11CD',all=FALSE)

  all_grouped <- all %>% group_by(uwwCode)
  #  aggregate by catchment
  for (ix in 1:length(X_names)) {
    x <- X_names[ix]
    if (op[ix]=='mean') tX <- unique((all_grouped %>% mutate(tmp = mean(!!as.name(x))))[c('uwwCode','tmp')])
    if (op[ix]=='sum') tX <- unique((all_grouped %>% mutate(tmp = sum(!!as.name(x))))[c('uwwCode','tmp')])
    if (op[ix]=='wmean') tX <- unique((all_grouped %>% mutate(tmp = weighted.mean(!!as.name(x), population)))[c('uwwCode','tmp')])

    names(tX)[2] <- x
    if (ix==1) X <- tX
    if (ix>1) X <- merge(X,tX,by='uwwCode',sort=FALSE)
  }
  return(X)
}


###########################################################################
#'  merge the catchment covariates with the wastewater data

#'
#' @return
#' @export
#'
merge_covariates_ww_data <- function(ww.weekly,covariates.data,covariate_list){

  nxs <- length(covariate_list)

  #  merge the catchment covariates with the wastewater data
  usite_codes <- unique(ww.weekly$uwwCode)
  for (ix in 1:nxs) {
    xn <- covariate_list[ix]
    mx <- mean(covariates.data[[xn]],na.rm=TRUE)
    tmp <- data.frame(tmp=rep(0,nrow(ww.weekly)))
    for (isite in 1:length(usite_codes)) {
      site <- usite_codes[isite]
      x.ids <- which(covariates.data$uwwCode==site)
      if (length(x.ids)>0) {
        tmp$tmp[which(ww.weekly$uwwCode==site)] <- covariates.data[[xn]][x.ids]
      } else {
        # catchments of some sites are not yet available, use the mean to replace them
        tmp$tmp[which(ww.weekly$uwwCode==site)] <- mx
      }
    }
    names(tmp) <- xn
    ww.weekly <- cbind(ww.weekly,tmp)
  }
  return (ww.weekly)
}



###########################################################################
#' Prepare the genomic data for model fitting
#'
#' @description
#' This function reads and formats the genomic data provided by UKHSA.
#' The genomic data cover the period from xx to 27/03/2022 across the sewage treatment works in England.
#' The variables of interest are num_snp and pct_genome_coverage. These two variables are aggregated
#' to the national level and entered the wastewater model as temporal-only covariates.
#' For the wastewater data available till 30/03/2022, the values of the genomic variables are
#' imputed using zoo::na.fill(x,'extend') that fills the NA by using rightmost non-NA value. A rolling average over
#' 3 weeks on each side is applied to the filled time series values, the resulting values from which enters the
#' wastewater model.
#'
#' @return
#' @export
get_weekly_genomics_data <- function(ww ,syears = NULL, smonths = NULL, analysis_unit, gen_covars) {
  data.file1 <- system.file("extdata", "meta_Imperial_SWT.csv", package = "publicWW")
  dat1 <- utils::read.csv(data.file1,header=TRUE)
  data.file2 <- system.file("extdata", "metaImperial2.csv", package = "publicWW")
  dat2 <- utils::read.csv(data.file2,header=TRUE)
  data.file3 <- system.file("extdata", "metaImperial3.csv", package = "publicWW")
  dat3 <- utils::read.csv(data.file3,header=TRUE)

  common_col_names <- intersect(names(dat1), names(dat2)) %>% intersect(., names(dat3))

  dat <- rbind(dat1[c(common_col_names)], dat2[c(common_col_names)], dat3[c(common_col_names)])

  dat <- dplyr::rename(dat, Date=date, uwwCode=sample_site_code)

  gendat <- dplyr::left_join(x = ww , y = dat, by = c('Date','uwwCode')) %>% rename(week = week.x)

  agg_df <- lapply(gen_covars, aggregate_time_simple, ww.dat=gendat, analysis_unit=analysis_unit, log_e_y=FALSE) %>%
    purrr::reduce(dplyr::left_join, by = c("site_index", "time_index","week","day1week","uwwCode","month","year","RegionName"))

  # summary of data availability at the weekly-site level
  ll <- NULL
  for (i in unique(agg_df$time_index)) {
	tmp <- subset(agg_df,agg_df$time_index==i)
	ll <- c(ll,length(which(is.na(tmp$num_snp)==FALSE)))
  }
  print("number of sites (out of 303 in total) with num_snp measurements across the 44 weeks")
  print(ll)

  ll <- NULL
  for (i in unique(agg_df$time_index)) {
	tmp <- subset(agg_df,agg_df$time_index==i)
	ll <- c(ll,length(which(is.na(tmp$pct_genome_coverage)==FALSE)))
  }
  print("number of sites (out of 303 in total) with pct_genome_coverage measurements across the 44 weeks")
  print(ll)
	
  # At the national level, take the average, linear interpolate missing values.
  #and assign them to the site-level
  agg_df <- agg_df %>% group_by(time_index)  %>%
    mutate(across(all_of(gen_covars), ~ mean(.x, na.rm = T), .names="{.col}")) %>%
    ungroup()

  ra <- agg_df %>%
    select(time_index, all_of(gen_covars)) %>%
    unique() %>%
    arrange(time_index) %>%
    mutate(across(all_of(gen_covars), ~ zoo::rollmean(zoo::na.fill(.x, fill='extend'), k=3, fill='extend'), .names="{.col}"))

  #  print warning if ww data exceeds April 2022
  d <- sapply(ww$Date,function(x){substr(x,4,nchar(x))})
  if (any(sapply(4:12,function(x){any(paste0('0',x,'/2022')==d)})==TRUE)) {
    warning('The genomic data are being extrapolated to future months. Make sure you are happy about that.')
  }
  return(ra)
}


###########################################################################
#'
#' @return
#' @export
aggregate_time_simple <- function(ww.dat,analysis_unit,outcome_col='gc',log_e_y=TRUE) {
  obs <- ww.dat[[outcome_col]]
  if (log_e_y) obs <- log(obs)
  tm <- ww.dat$time_index
  site <- ww.dat$site_index
  ###  mean gc by week/month by site
  if (analysis_unit$time=='week') {
    agg <- stats::aggregate(obs,list(site_index=site,time_index=tm),function(x){mean(x,na.rm=TRUE)})
    ###  add week index
    agg <- add_index(ww.dat,agg,by='time_index',index_to_add='week',type='numeric')
    ###  add month index
    agg <- add_index(ww.dat,agg,by='time_index',index_to_add='month',type='numeric')
  }
  if (log_e_y) outcome_col <- paste0('log_e_',outcome_col)
  names(agg)[which(names(agg)=='x')] <- outcome_col
  agg[[outcome_col]][which(is.nan(agg[[outcome_col]]))] <- NA
  ###  add year index
  agg <- add_index(ww.dat,agg,by='month',index_to_add='year',type='numeric')
  ###  add first day of the week
  agg <- add_index(ww.dat,agg,by='week',index_to_add='day1week',type='character')
  ###  add region index
  agg <- add_index(ww.dat,agg,by= 'site_index',index_to_add='uwwCode',type='character')
  agg <- add_index(ww.dat,agg,by= 'site_index',index_to_add='RegionName',type='character')

  return(agg)
}


###########################################################################
#' Assign a numerical index to regions in England
#'
#' @description
#' The order of the regions comes from the names of table(ww.weekly$RegionName).
#'
#' @return
#' @export
add_region_index <- function(rgn) {
  levels <- c("East Midlands","East of England","London",
              "North East","North West","South East",
              "South West","West Midlands","Yorkshire and The Humber")
  ff <- factor(rgn,levels=levels)
  out <- as.numeric(ff)
  names(out) <- rgn
  return(out)
}
