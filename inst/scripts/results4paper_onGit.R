if (FALSE) {
	source(system.file("scripts", "results4paper_onGit.R", package = "publicWW"))
}

rm(list=ls())
library(publicWW)

library(ggplot2)
library(sf)
library(dplyr)
library(cowplot)
library(publicWW)
library(ggridges)
library(biscale)
library(inlabru)
library(gridExtra)

produce.plot <- !FALSE 

#  this is the folder where you store the three files of prediction (see https://github.com/gqlNU/publicWW#visualising-predictions)
workdir <- '/Volumes/WorkSpace/publicWW_results/'
setwd(workdir)

figdir <- '/Volumes/WorkSpace/publicWW_results/ww_figures/'

model <- 'ww_model_fit-20220801-195150'  #  regional random effect
aggregate.method <- 'pw'
pred <- load_pred(model,'LAD',aggregate.method)
lp <- get.geo.lookups()

######################################################################
###  figure 2: mesh
######################################################################
#  UK boundary
shp <- sf::st_read(system.file("extdata", "CTRY_DEC_2021_GB_BUC.shp", package = "publicWW"))
id <- grep('England',as.character(shp$CTRY21NM),ignore.case=TRUE)
shp <- shp[id,1]
shp <- as(shp,'Spatial')
n <- length(shp@polygons[[1]]@Polygons)
for (i in 1:n) {
	shp@polygons[[1]]@Polygons[[i]]@coords <- shp@polygons[[1]]@Polygons[[i]]@coords/1000
}

load(system.file("extdata", "ww_mesh-20220801-195150.RData", package = "publicWW"))
locs <- get_sites_geo((read_public_data(upto = '30Mar2022'))[c('uwwCode','uwwName')])

pt <-
      ggplot() +
      gg(mesh) +
      geom_point(data=data.frame(locs),aes(eastings,northings)) +
			gg(shp,alpha=0.15,lwd=0.3) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"),
            axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

if (produce.plot) {
  pdf(file=paste0(figdir,'mesh.pdf'),height=7,width=7)
  print(pt)
  dev.off()
}



################################################################
#  Fig3a: national trend
################################################################

ntimes <- dim(pred)[[3]]
week.start <- read.csv(system.file("extdata", "time_index.csv",package = "publicWW"))[1:ntimes,]
week.ids <- pick.1st.week(week.start)
##   1a: national trend  (use the nationally aggregated estimates when ready)
england.pred <- load_pred(model,'CTRY',aggregate.method)
national <- cbind(t(apply(england.pred[,1,],2,post_summary)),week.start)

ylab <- 'Viral concentration (log(gc/L))'

p1 <- ggplot(national, aes(x = time_index, y = mean)) +
  labs(title='',x='Week starting',y=ylab) +
  geom_line(col=1) +
  geom_point(col=1) +
  geom_ribbon(aes(ymin = lowCI95, ymax = highCI95), alpha = 0.3, linetype=0, color="grey") +
  theme_bw() +
  theme(panel.border = element_blank(),
             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks = week.ids, labels = names(week.ids)) +
  theme(axis.text=element_text(size=13),
             axis.title.y = element_text(size = 18),
             axis.title.x = element_text(size = 18))



################################################################
#  Fig3b: spatial distribution of mean log(gc/L)
################################################################
shpfile <- '/Volumes/WorkSpace/publicWW/inst/extdata/Local_Authority_Districts_(May_2021)_UK_BFE_V3/LAD_MAY_2021_UK_BFE_V2.shp'
ltla.shp <- sf::st_read(shpfile)

values <- data.frame(mean=apply(pred,2,mean),LAD21CD=dimnames(pred)[[2]])
p2 <- ltla_map(ltla.shp,values,inset.region='London')

################################################################
#  Fig3c: boxplot of mean log(gc/L) by region
################################################################
rgn.pred <- load_pred(model,'RGN',aggregate.method)

nsims <- dim(pred)[1]
regions <- unique(lp$lookup[c('RGN21NM','RGN21CD')])
region.cd <- dimnames(rgn.pred)[[2]]
nregions <- nrow(regions)
rgs <- rgsind <- NULL
for (ir in 1:nregions) {
	r <- regions[ir,'RGN21NM']
	dc <- regions[ir,'RGN21CD']
	mn <- apply(rgn.pred[,which(region.cd==dc),],1,mean)
	rgs <- c(rgs,mn)
	rgsind <- c(rgsind,rep(r,length(mn)))
}
rgs <- data.frame(region=rgs)
rgsind[which(rgsind=='Yorkshire and The Humber')] <- 'Yorkshire and \nThe Humber'
rgs$ind <- factor(rgsind)
mn <- tapply(rgs$region,rgs$ind,mean)
od <- names(mn)[order(mn,decreasing=FALSE)]
rgs$ordered <- factor(rgsind,levels=od)

p3 <- ggplot(rgs, aes(x = region, y = ordered)) +
	labs(x='Viral concentration (log(gc/L))',y='')+
	geom_density_ridges(scale = 0.96,
	                    rel_min_height = 0.0005,
	                    quantile_lines = TRUE, quantiles = c(0.025, 0.975)) +
	 theme_bw() +
theme(panel.border = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=14),
  axis.title.y = element_text(size = 18),
  axis.title.x = element_text(size = 18))
  
###  putting Figure 1 together
if (produce.plot) {
  pdf(file=paste0(figdir,'fig3_30March2022.pdf'),height=12,width=12)
  pt <- ggdraw() +
             draw_plot(p1, x = 0, y = .65, width=1,height=0.35) +
             draw_plot(p2, x = 0, y = 0.03, width=0.55,height=0.6) +
             draw_plot(p3, x = 0.55, y = 0, width=0.45,height=0.6) +
             draw_plot_label(label = c("A", "B", "C"), size = 35,
                                         x = c(rep(0.05,2) ,0.55), y = c(0.97, rep(0.63,2)))
  print(pt)
  dev.off()
}


################################################################
#  fig4 map sequence
################################################################
#week.ids <- pick.1st.week(week.start)
ids <- week.start[floor(seq(1,ntimes,length.out=12)),]
week.ids <- ids$time_index
names(week.ids) <- ids$week_start

#  define colour scheme
mm <- apply(pred,c(2,3),mean)
breaks <- seq(min(mm)*0.99,max(mm)*1.01,length.out=25)
nbks <- length(breaks)
labels <- c(paste0(breaks[1:(nbks-2)],'-'),paste(breaks[nbks-1],breaks[nbks],sep='-'))
labels <- 1:(nbks-1)
#  from https://colorbrewer2.org/#type=diverging&scheme=RdYlBu&n=11
hexc <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')[11:1]
palette <- colorRampPalette(hexc)
cols <- palette(nbks-1)  #RColorBrewer::brewer.pal(nbks-1,'RdYlBu')[(nbks-1):1]
names(cols) <- labels
colour.scheme <- list(breaks = breaks, labels = labels,cols=cols)

#  going through the map sequence
nweeks <- length(week.ids)
p <- as.list(1:nweeks)
for (iw in 1:nweeks) {
#  values <- data.frame(mean=apply(pred[,,week.ids[iw]],2,mean),LAD21CD=dimnames(pred)[[2]])
  values <- data.frame(mean=mm[,week.ids[iw]],LAD21CD=dimnames(pred)[[2]])
  p[[iw]] <- ltla_map(ltla.shp,values,inset.region='London',colour.scheme)
}

if (produce.plot) {
  ##  for a 3 x 4 grid
  w <- 0.32*1.15
  h <- 0.37*1.1
  xd <- 0.335
  sz <- 4
  pdf(file=paste0(figdir,'fig4_30March2022.pdf'),height=sz*4,width=sz*3)
  x <- c(0,0.3,0.6)*1.15
  y <- c(1,0.67,0.33,0)*1.5
  pt <- ggdraw(xlim=c(0,1),ylim=c(0,2)) +
    #  row 1
    draw_plot(p[[1]], x =x[1] , y = y[1], height=h,width=w) +
    draw_plot(p[[2]], x =x[2] , y = y[1], height=h,width=w) +
    draw_plot(p[[3]], x =x[3] , y = y[1], height=h,width=w) +
    #  row 2
    draw_plot(p[[4]], x =x[1] , y = y[2], height=h,width=w) +
    draw_plot(p[[5]], x =x[2] , y = y[2], height=h,width=w) +
    draw_plot(p[[6]], x =x[3] , y = y[2], height=h,width=w) +
    #  row 3
    draw_plot(p[[7]], x =x[1] , y = y[3], height=h,width=w) +
    draw_plot(p[[8]], x =x[2] , y = y[3], height=h,width=w) +
    draw_plot(p[[9]], x =x[3] , y = y[3], height=h,width=w) +
    #  row 4
    draw_plot(p[[10]], x =x[1] , y = y[4], height=h,width=w) +
    draw_plot(p[[11]], x =x[2] , y = y[4], height=h,width=w) +
    draw_plot(p[[12]], x =x[3] , y = y[4], height=h,width=w) +
    draw_plot_label(label = names(week.ids), size = 20,
                               x = rep(x,4)+0.07, y=rep(y+y[1]-y[2]-0.02,each=3)) +
    draw_plot(add.map.legend(breaks,cols), x = x[3]+0.05, y = 0.229, height=0.25,width=0.06)
  print(pt)
  dev.off()
}


#################################################
#  fig5
#################################################
#  getting draw-level prevalence 
lad21cd <- dimnames(pred)[[2]]
prev.sims <- get_prev_sims(nsims,lad21cd)$prev.sims

nts.prev <- dim(prev.sims)[3]
lag <- 0
cpv.sim <- prev.sims[,,2:nts.prev] - prev.sims[,,1:(nts.prev-1)]
cpw.sim <- pred[,,2:nts.prev] - pred[,,1:(nts.prev-1)]


hex <- c('#850344','#de77ae','#f1b6da','#fde0ef','#fbf1f6','#f7f7f7','#b8e186','#7fbc41','#285c12')
palette <- colorRampPalette(hex)
cols <- palette(9)

tcs <- as.list(1:2)
nts.start <- 1
nts.end <- 3
tcs[[1]] <- 1:2  #  Jun 2021
tcs[[2]] <- 19:20 # Oct 2021

mt <- c('Jun01 to Jun20','Jun7 to Jun27','Jun14 to Jul04')
mt <- c(mt, c('Oct04 to Oct24','Oct11 to Oct31','Oct18 to Nov07'))
# custom palette
custom_pal <- c("1-1" = '#dbf7ba', # low x, low y
		                   "2-1" = as.character(cols[5]), # high x, low y
                           "1-2" = as.character(cols[2]), # low x, high y
                           "2-2" = as.character(cols[9]) # high x, high y
		                  )
bp <- bi_pal(pal = custom_pal, dim = 2, preview = FALSE)
legend <- bi.legend(cols = bp,
                    				 xlab = "Pr(inc. in ww) ",
                     				 ylab = "Pr(inc. in prev) ",
                     				 size.lab=10,size.bks=9)

ic <- 0
final.plots <- as.list(1:length(mt))
for (tpp in 1:length(tcs)) {
	for (t in nts.start:nts.end) {
		ic <- ic + 1
		tp <- tcs[[tpp]]+(t-1)
		tw <- tcs[[tpp]]+(t-1)
		pdw <- data.frame(LAD21CD=lad21cd,pdw=apply(cpw.sim[,,tw],2,ppp))
		pdp <- data.frame(LAD21CD=lad21cd,pdp=apply(cpv.sim[,,tp],2,ppp))
		pwp <- merge(pdw,pdp)
		pwp$cdw <- cut(pwp$pdw,breaks=c(0,0.8,1),include.lowest=TRUE)
		pwp$cdp <- cut(pwp$pdp,breaks=c(0,0.8,1),include.lowest=TRUE)
		shp <- merge(ltla.shp,pwp,by='LAD21CD',all.x=FALSE,all.y=TRUE)
		bc <- biscale::bi_class(shp, x = cdw, y = cdp, style = "equal", dim = 2)

	    map <- ggplot() +
				geom_sf(data = bc, mapping = aes(fill = bi_class), color = "white", size = 0.01, show.legend = FALSE) +
				bi_scale_fill(pal = bp, dim = 2) +
				labs(title = '',subtitle = "") +
				theme_bw() +
	            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
	               panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"),
	               axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
	               axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()
	              )
		### combine map with legend
	    final.plots[[ic]] <- ggdraw() +
  		                              draw_plot(map, 0, 0, 1, 1) +
                                      draw_text(mt[ic],0.55,0.95,size=24) #+ draw_plot(legend, 0.1, 0.55, 0.2, 0.2)
    }
}
if (produce.plot) {
  p <- final.plots
  ##  for a 3 x 4 grid
  w <- 0.55
  h <- 0.55
  xd <- 0.335
  sz <- 4
  pdf(file=paste0(figdir,'fig5_30March2022.pdf'),height=sz*3,width=sz*3)
  x <- c(0,0.3,0.6)*1.15-0.2
  y <- c(1,0.67,0.33,0)*1.75
  pt <- ggdraw(xlim=c(-0.15,1),ylim=c(1,2.4)) +
    #  row 1
    draw_plot(p[[1]], x =x[1] , y = y[1], height=h,width=w) +
    draw_plot(p[[2]], x =x[2] , y = y[1], height=h,width=w) +
    draw_plot(p[[3]], x =x[3] , y = y[1], height=h,width=w) +
    #  row 2
    draw_plot(p[[4]], x =x[1] , y = y[2], height=h,width=w) +
    draw_plot(p[[5]], x =x[2] , y = y[2], height=h,width=w) +
    draw_plot(p[[6]], x =x[3] , y = y[2], height=h,width=w) +
    draw_text('-->',x=c(0.295,0.65),y=rep(1.4,2),size=25) +
    draw_text('-->',x=c(0.295,0.65),y=rep(2.0,2),size=25) +
    draw_text(c('A','B'),x=rep(-.08,2),y=c(2.28,1.7),size=35) +
    draw_plot(legend, 0.57, 1.45, 0.12, 0.13)
  print(pt)
  dev.off()
}


#################################################
#  sup fig4:  correlation of level between WW and prevalence
#################################################
###########  A: correlation of level over time
week.start <- read.csv(system.file("extdata", "time_index.csv",package = "publicWW"))[1:ntimes,]
week.ids <- pick.1st.week(week.start)
nlads <- length(lad21cd)
lad.ids <- (1:nlads)[-sapply(c('E06000053','E09000001'),function(x){which(lad21cd==x)})]
###  correlation of level over time
nts.prev <- length(which(prev.sims[1,lad.ids[1],]!=0)) #  the last week does not have prevalence
lag <- 0
cr.sim <- sapply(1:nsims,function(isim) {
	sapply(1:nts.prev,function(tt){cor(pred[isim,lad.ids,tt],prev.sims[isim,lad.ids,tt],method='spearman')})
})
crs <- t(apply(cr.sim,1,function(x){quantile(x,c(0.025,0.5,0.975))}))
colnames(crs) <- c('low','mean','high')
crs <- as.data.frame(crs)
crs$index <- week.start$time_index[1:nts.prev]
A <- forestplot(crs,main="Correlation of level between log(gc/L) and debiased prevalence",week.ids=week.ids)

###########  B: correlation of level over LADs
###  correlation between levels over LAD
lag <- 0
cr.sim <- sapply(1:nsims,function(isim) {
  sapply(lad.ids,function(x){
    cor(prev.sims[isim,x,(1+lag):nts.prev],pred[isim,x,1:(nts.prev-lag)],method='spearman')
  })
})
crs <- t(apply(cr.sim,1,function(x){c(mean(x),quantile(x,c(0.025,0.975)))}))
colnames(crs) <- c('mean','low','high')
crs <- as.data.frame(crs)
crs$index <- lad.ids
pk <- apply(cr.sim,1,function(x){length(which(x>0))/length(x)})
print(paste0('There are ',length(which(pk>0.95)),' LTLAs with non negative correlation in levels > 0.95 '))

#  colour scheme
breaks <- c(-0.4,seq(0,max(0.7,max(crs$mean)),length.out=40))
nbks <- length(breaks)
labels <- c(paste0(breaks[1:(nbks-2)],'-'),paste(breaks[nbks-1],breaks[nbks],sep='-'))
labels <- 1:(nbks-1)
hexc <- c('#f7fcf5','#74c476','#00441b')
palette <- colorRampPalette(hexc)
cols <- palette(nbks-1)
names(cols) <- labels
colour.scheme <- list(breaks = breaks, labels = labels,cols=cols)

values <- data.frame(mean=c(crs$mean),LAD21CD=dimnames(pred)[[2]][lad.ids])
B <- ltla_map(ltla.shp,values,inset.region='London',colour.scheme)
#   add legend
db <- c(round(min(values$mean),digits=1),0,.2,.4,.6)
names(db) <- db
B <- ggdraw(xlim=c(0,1),ylim=c(0,1)) +
        draw_plot(B, x =0 , y = 0, height=1,width=1) +
        draw_plot(add.map.legend(breaks,cols,main='Cor. between levels',dis.bks=db), x = 0.12, y = 0.55, height=0.4,width=0.15)

###########  C: Prob(correlation of level over LADs > 0)
breaks <- c(-0.001,0.05,0.95,1)
nbks <- length(breaks)
labels <- c(paste0(breaks[1:(nbks-2)],'-'),paste(breaks[nbks-1],breaks[nbks],sep='-'))
labels <- 1:(nbks-1)
hexc <- c('#dcf1fc','#92c5de','#4393c3')
palette <- colorRampPalette(hexc)
cols <- palette(nbks-1)
names(cols) <- labels
colour.scheme <- list(breaks = breaks, labels = labels,cols=cols)

crs <- t(apply(cr.sim,1,function(x){c(length(which(x>0))/length(x),length(x))}))
colnames(crs) <- c('mean','l')
crs <- as.data.frame(crs)
crs$index <- lad.ids
values <- data.frame(mean=c(crs$mean),LAD21CD=dimnames(pred)[[2]][lad.ids])
C <- ltla_map(ltla.shp,values,inset.region='London',colour.scheme)
db <- c(0,0.05,0.95,1)
names(db) <- db
C <- ggdraw(xlim=c(0,1),ylim=c(0,1)) +
        draw_plot(C, x =0 , y = 0, height=1,width=1) +
        draw_plot(add.map.legend(breaks,cols,main='Posterior prob. of \nCor. between levels > 0',dis.bks=db), x = 0.12, y = 0.55, height=0.4,width=0.15)

if (produce.plot) {
  pdf(file=paste0(figdir,'sup_fig4_30March2022.pdf'),height=12,width=12)
  pt <- ggdraw(xlim=c(-0.05,1.15),ylim=c(0,1.1)) +
             draw_plot(A, x = 0, y = .75, width=1.1,height=0.32) +
             draw_plot(B, x = 0, y = 0.03, width=0.51,height=0.65) +
             draw_plot(C, x = 0.56, y = 0.03, width=0.51,height=0.65)+
             draw_plot_label(label = c("", "A", "B","C"), size = 35, x = c(rep(0.00,3) ,0.54), y = c(1.45,1.1, rep(0.72,2)))
  print(pt)
  dev.off()
}

#################################################
#  sup fig5: nonlinear relationship over space and time
#################################################
#  note: City of London and Isles of Scilly are not in the prevalence data
lad21cd <- dimnames(pred)[[2]]
nlads <- length(lad21cd)
lad.ids <- (1:nlads)[-sapply(c('E06000053','E09000001'),function(x){which(lad21cd==x)})]
names(lad.ids) <- lad21cd[lad.ids]

#  plot by region
rgn <- rep(0,length(lad.ids))
for (i in 1:length(lad.ids)) rgn[i] <- lp$lookup$RGN21NM[which(lp$lookup$LAD21CD==names(lad.ids)[i])[1]]
rgn.ids <- add_region_index(rgn)
names(rgn.ids) <- rgn

#  space-time mean for WW and prev
nts.prev <- length(which(prev.sims[1,lad.ids[1],]!=0)) #  the last week does not have prevalence
mpw <- apply(pred[,lad.ids,1:nts.prev],c(2,3),mean)
mpv <- apply(prev.sims[,lad.ids,1:nts.prev],c(2,3),mean)

#  plotting
hexc <- c('#3288bd','#4d9221','#fee08b','#d53e4f')
palette <- colorRampPalette(hexc)
cols <- palette(10)
xl <- range(c(mpw))
yl <- range(c(mpv))
mg <- sapply(1:9,function(x){names(rgn.ids)[which(rgn.ids==x)[1]]})
g <- list(1:10)
for (ir in 1:9) {
  x <- mpw[which(rgn.ids==ir),]
  y <- mpv[which(rgn.ids==ir),]
  wk <- get.starting.day()
  w <- p <- week <- month <- NULL
  im0 <- as.numeric(substr(wk$week_start[1],6,7))
  im <- 1
  for (tt in 1:ncol(x)) {
    im1 <- as.numeric(substr(wk$week_start[tt],6,7))
    if (im1!=im0) {
      im0 <- im1
      im <- im + 1
    }
    w <- c(w,x[,tt])
    p <- c(p,y[,tt])
    week <- c(week,rep(tt,length(x[,tt])))
    month <- c(month,rep(im,length(x[,tt])))
  }
  dd <- data.frame(ww=w,prev=p,week=week,month=month)
  g[[ir]] <- ggplot(dd, aes(x=ww, y=prev)) +
    geom_point(aes(color = month), size = 1.5,alpha=0.6, stroke=0.01) +
      theme_bw() +
      scale_color_gradientn(colors = cols) + 
      theme(legend.position="none") +
      xlim(xl[1],xl[2]) + ylim(yl[1],yl[2]) +
      labs(x ="", y = "") +
      ggtitle(mg[ir]) + theme(plot.title=element_text(size=10))
}
db <- 0:9 + 0.5
names(db) <- c('Jun21','Jul21','Aug21','Sept21','Oct21','Nov21','Dec21','Jan22','Feb22','Mar22')
g[[10]] <- add.map.legend(0:10,cols,main='',dis.bks=db)

hlay <- rbind(c(1,1,2,2,3,3,NA),
              c(4,4,5,5,6,6,NA),
              c(7,7,8,8,9,9,10))
if (produce.plot) {
  pdf(file=paste0(figdir,'sup_fig5_30March2022.pdf'),height=8,width=8)
  grid.arrange(grobs = g, layout_matrix = hlay,
                       top="", bottom="Viral concentration in wastewater (log(gc/L))",
                       left="Debiased prevalence on logit scale", right="")
  dev.off()
}
