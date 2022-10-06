###########################################################################
#'  Load predictions of viral concentration 
#'
#' @return
#' @export
load_pred <- function(model,geography,aggregate.method) {
  file <- paste0(model,'_',geography,'pred_',aggregate.method,'.RData')
	load(file)
	return(agg.sims)
}

###########################################################################
#' Picking the first week of each month
#'
#' @return
#' @export
pick.1st.week <- function(week.start) {
	w <- week.start$week_start
	n <- length(w)
	m <- as.numeric(sapply(w,function(xx){substr(xx,6,7)}))
	ids <- c(1,which((m[1:(n-1)]-m[2:n])!=0)+1)
	names(ids) <- sapply(w[ids],function(x){sub('-','/',sub('-','/',x))})
	return(ids)
}


###########################################################################
#' posterior summary
#'
#' @return
#' @export
post_summary <- function(x) {
	out <- c(mean(x),sd(x),quantile(x,c(0.025,0.5,0.975)))
	names(out) <- c('mean','sd','lowCI95','q50','highCI95')
	return(out)
}


###########################################################################
#' Get the coordinate for the map inset
#'
#' @return
#' @export
get_inset_coords <- function(selected.region='London',ltla.shp) {
  lp <- get.geo.lookups()
  selected.ltlas <- unique(lp$lookup$LAD21CD[which(lp$lookup$RGN21NM== selected.region)])
  selected.coords <- sapply(selected.ltlas,function(x){
	(st_coordinates(ltla.shp[which(ltla.shp$LAD21CD==x),'geometry'])[,c('X','Y')])
  })
  lc <- NULL
  for (i in 1:length(selected.ltlas)) {lc <- rbind(lc, selected.coords[[i]])}
  zx <- range(lc[,1])*c(0.99,1.01)
  zy <- range(lc[,2])*c(0.99,1.01)
  return(list(zx=zx,zy=zy))
}



###########################################################################
#' produce various LTLA maps
#'
#' @return
#' @export
ltla_map <- function(ltla.shp,values,inset.region='London',colour.scheme=NULL) {
  shp <- merge(ltla.shp,values,by='LAD21CD',all.x=FALSE,all.y=TRUE)
  #  main map
  if (is.null(colour.scheme)) {
    main.map <- shp %>%
    ggplot() +
      geom_sf(aes(fill = mean),lwd = 0.001,colour = "light grey") +
	  scale_fill_distiller(palette = "RdYlBu",name='Viral conc. \nlog(gc/L)') + 
	  theme_void() +
	  theme(
        # legend.justification defines the edge of the legend that the legend.position coordinates refer to
        legend.justification = c(0, 1),
        # Set the legend flush with the left side of the plot, and just slightly below the top of the plot
        legend.position = c(0.1, 0.85),
        legend.text = element_text(size=14),
        legend.title = element_text(size=18),
        legend.key.size = unit(.7, 'cm')
      ) + 
      coord_sf(expand = FALSE)  # Prevent ggplot from slightly expanding the map limits beyond the bounding box of the spatial objects
  } else {
    main.map <- shp %>%
      ggplot() +
	  geom_sf(aes(fill = cut(mean,breaks= colour.scheme$breaks,labels= colour.scheme$labels)),
	                 lwd = 0.001,colour = "light grey") +
	  scale_fill_manual(values=colour.scheme$cols,limits=colour.scheme$labels) +
	  theme_void() +
	  theme(legend.position='none') +
	  coord_sf(expand = FALSE)  # Prevent ggplot from slightly expanding the map limits beyond the bounding box of the spatial objects
  }
  #  putting the main map and the inset together
  if (!is.null(inset.region)) {
    coords <- get_inset_coords(inset.region,ltla.shp)
    zx <- coords$zx
    zy <- coords$zy
    #  add a rectangle on the main map to indicate inset
    main.map <- main.map + geom_rect(xmin = zx[1],ymin = zy[1],xmax = zx[2],ymax = zy[2],
                                                                fill = NA,colour = "black",size = 0.6)
    final <- ggdraw(main.map) +
                draw_plot({main.map + coord_sf(xlim = zx,ylim = zy,expand = FALSE) +
                                   theme(legend.position = "none")
                                  },
                                  # The distance along a (0,1) x-axis to draw the left edge of the plot
                                  x = 0.02,
                                  # The distance along a (0,1) y-axis to draw the bottom edge of the plot
                                  y = 0.2,
                                  # The width and height of the plot expressed as proportion of the entire ggdraw object
                                  width = 0.35,height = 0.35)
  } else {#  no inset
  	final <- map.map
  }
  return(final)
}

###########################################################################
#'  Add legend to map
#'
#' @return
#' @export
add.map.legend <- function(breaks,cols,main='log(gc/L)',dis.bks=c('6'=6,'8'=8,'10'=10,'12'=12)) {
	#https://stackoverflow.com/questions/41861079/ggplot2-plotting-several-boxes-using-a-loop
	nbks <- length(breaks)
	n <- nbks - 1
	df <- data.frame(xmin = rep(0, n),
	                 xmax = rep(1,n),
	                 ymin = breaks[1:n],
	                 ymax = breaks[2:nbks],
	                 fill = cols)
	gp <- ggplot(df) +
	           geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill)) +
	           scale_fill_identity() +
	           theme_bw() +
	           theme(panel.border = element_blank(), panel.grid.major = element_blank(),
	                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"),
	                      axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()
	                     )  +
	           scale_y_continuous(breaks=dis.bks,labels=names(dis.bks)) +
	           ggtitle(main) +
	           theme(plot.title = element_text(hjust = 0.7))
	return(gp)
}
#add.map.legend(breaks,cols)

###########################################################################
#' Get simulated values for the debiased prevalence on the logit scale
#'
#' @return
#' @export
get_prev_sims <- function(nsims,lad21cd) {
  data.file <- system.file("extdata", "logit_moments.rds", package = "publicWW")
  prev <- readRDS(file=data.file)
  mp <- 'mean'
  sp <- 'sd'
  nts.prev <- max(prev$time_index_ww)
  nlads <- length(lad21cd)
  prev.sims <- array(0,c(nsims,nlads,ntimes))
  for (i in 1:nlads) {
    for (t in 1: nts.prev) {
	  id <- which(prev$LAD21CD==lad21cd[i] & prev$time_index_ww==t)
	  if (length(id)>0) {
	    prev.sims[,i,t] <- rnorm(nsims,prev[[mp]][id],prev[[sp]][id])
	  }
    }
  }
  return(list(prev.sims=prev.sims,mean=mean,sd=sd))
}

###########################################################################
#' Posterior probability
#'
#' @return
#' @export
ppp <- function(x,thres=0) {
	if (is.null(dim(x))) {
		nx <- length(x)
		#  rise/fall of one week
		out <- length(which(x>thres))/nx
	} else {
		nx <- dim(x)[[1]]
		if (dim(x)[[2]]==2) {
			#  rise/fall of 2 consecutive weeks
			out <- length(which(x[,1]>thres & x[,2]>thres))/nx
		}
		if (dim(x)[[2]]==3) {
			#  rise/fall of 3 consecutive weeks
			out <- length(which(x[,1]>thres & x[,2]>thres & x[,3]>thres))/nx
		}
		if (dim(x)[[2]]==4) {
			#  rise/fall of 4 consecutive weeks
			out <- length(which(x[,1]>thres & x[,2]>thres & x[,3]>thres & x[,4]>thres))/nx
		}
	}
	return(out)
}

###########################################################################
#'  Construct biscale
#'
#' @return
#' @export
bi.legend  <- function(cols,xlab,ylab,dis.bks=c(0,0.8,1),counts=NULL,size.lab,size.bks) {
	#https://stackoverflow.com/questions/41861079/ggplot2-plotting-several-boxes-using-a-loop
	nbs <- length(cols)
	coords <- fill <- NULL
	for (ib in 1:nbs) {
		nb <- as.numeric(unlist(strsplit(names(cols)[ib],'-')))
		xmin <- nb[1]-1
		xmax <- nb[1]
		ymin <- nb[2]-1
		ymax <- nb[2]
		tmp <- c(xmin,xmax,ymin,ymax)
		coords <- rbind(coords,tmp)
		fill <- c(fill,cols[ib])
	}
	coords <- as.data.frame(coords)
	names(coords) <- c('xmin','xmax','ymin','ymax')
	coords$fill <- fill
	gp <- ggplot(coords) +
	           geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill)) +
	           scale_fill_identity() +
	           theme_bw() +
	           theme(panel.border = element_blank(), panel.grid.major = element_blank(),
	                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"),
	                      axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
	                      axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()
	                     )

	                     l1 <- c(0.15,0.6,1.05)
	                     l2 <- rep(0.13,3)
    out <-
    ggdraw(xlim=c(0,1.1),ylim=c(0,1.1)) +
  			  draw_plot(gp, 0.1, 0.1, 1,1) +
  			  draw_text(dis.bks,x=l1,y=l2,size=size.bks) +
  			  draw_text(dis.bks[2:length(dis.bks)],x=l2[2:length(dis.bks)],y=l1[2:length(dis.bks)],size=size.bks) +
  			  draw_text(xlab,0.6,0.03,size=size.lab) +
  			  draw_text(ylab,0.01,0.6,angle=90,size=size.lab)
#  			  draw_text(xlab,0.6,0.05,size=size.lab) + draw_text('      \u2192',0.8,0.065,size=size.arrow) +
#  			  draw_text(ylab,0.025,0.6,angle=90,size=size.lab) + draw_text('      \u2191',-0.069,0.95,size=size.arrow)

  	if (!is.null(counts)) {
  		out <- out  + draw_text(counts,x=c(0.4,.83,0.4,.83),y=c(0.4,0.4,.83,.83),size=25)
  	}
	return(out)
}


###########################################################################
#' for producing a forest plot to visualise level correlation between WW and debiased prevalence
#'
#' @return
#' @export
forestplot <- function(d, xlab='week start', main="Correlation of level between log(gc/L) and debiased prevalence",ylab='Correlation',week.ids) {
  p <- ggplot(d, aes(x=index, y=mean, ymin=low, ymax=high)) +
         geom_pointrange() +
         geom_hline(aes(yintercept=0), lty=2) +
         ylab(ylab) +
         xlab(xlab) +
         ggtitle(main) +
         theme_bw() +
         theme(panel.border = element_blank(), 
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
                    scale_x_continuous(breaks = week.ids, labels = names(week.ids)) +
         theme(axis.text=element_text(size=13),
                    plot.title = element_text(size=18),
                    axis.title.y = element_text(size = 16),
                    axis.title.x = element_text(size = 16))
    return(p)
}

