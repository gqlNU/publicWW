#  source('/Volumes/WorkSpace/publicWW/scrapyard/aggregate/summarise_aggregate.R')


rm(list=ls())

PAT <- "" # enter token here
devtools::install_github("gqlNU/publicWW",ref='INLApredict',auth_token = PAT,force=TRUE)

library(publicWW)

workdir <- '/Volumes/WorkSpace/publicWW_results/'
setwd(workdir)

library(publicWW)

posterior.summary <- function(x) {
	out <- c(mean(x),sd(x),quantile(x,c(0.025,0.5,0.975)))
	return(out)
}

##############################
#  user inputs
##############################
model <- 'Jun2021_to_Jan2022_allsites_boundary_and_site_mesh'

pw <- TRUE  #  population-weighted aggregation from LSOAs

geography <- 'RGN'  #  output geography
#geography <- 'CTRY'  #  output geography
#geography <- 'LAD'  #  output geography

################################################
#  extract predictions from two-stage 
################################################
file <- paste0(model,'_',geography,'pred.RData')
load(file)
ntimes <- dim(agg.sims)[[3]]
nareas <- dim(agg.sims)[[2]]

################################################
#  get area names
################################################
lookup <- get.geo.lookups()
code <- paste0(geography,'21')
name <- paste0(geography,'21NM')
nms <- lookup$code.name[[code]]
code <- paste0(geography,'21CD')

################################################
#  starting day of each week
################################################
week_start <- get.starting.day()[1:ntimes,]

################################################
#  output 1: posterior summaries per area-time cell
################################################
p <- apply(agg.sims,c(2,3),posterior.summary)
if (pw) {
	p <- apply(agg.sims.pw,c(2,3),posterior.summary)
}
area.codes <- dimnames(agg.sims)[[2]]
s <- NULL
for (iarea in 1:nareas) {
	out <- data.frame(area.code=area.codes[iarea],
	                  area.name=nms[[name]][which(nms[[code]]==area.codes[iarea])],
	                  time_index=1:ntimes,
	                  week_start=week_start$week_start)
	tmp <- t(p[,iarea,])
	colnames(tmp) <- c('mean','sd','q025','q50','q975')
	tmp <- as.data.frame(tmp)
	tmp <- cbind(out,tmp)
	s <- rbind(s,tmp)
}
names(s)[1:2] <- c(code,name)
#  save output
save.file <- gsub('.RData','_summary.csv',file)
if (pw) {
	save.file <- gsub('pred','pred_pw',save.file)
}
write.csv(s,file=save.file,quote=FALSE,row.names=FALSE)