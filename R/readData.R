###########################################################################
#' Read the publicly available data into R
#'
#' This function reads in the public wastewater data (already included in the package)
#'
#' The data can be donwloaded via https://www.gov.uk/government/publications/monitoring-of-sars-cov-2-rna-in-england-wastewater-monthly-statistics-15-july-2020-to-30-march-2022
#'
#' @return A wide-format wastewater data matrix (site by day)
#' @export
read_public_data <- function(upto='30Mar2022') {
	if (upto=='Jan2022') {
		#  read data between Jun2021 and Jan2022
		data.file <- system.file("extdata", "Final_EMHP_Wastewater_Data_January_2022.ods", package = "publicWW")
		dat <- readODS::read_ods(data.file,sheet=5,skip=4)
	}
	if (upto=='07Mar2022') {
		#  read data between Jun2021 and 07Mar2022
		data.file <- system.file("extdata", "EMHP_wastewater_data_March_-1.xlsx", package = "publicWW")
		dat <- readxl::read_excel(data.file,sheet=5,skip=4)
	}
	if (upto=='30Mar2022') {
		#  read data between Jun2021 and 30Mar2022
		data.file <- system.file("extdata", "EMHP_SARS-CoV-2_RNA_Concentration_May_2022.xlsx", package = "publicWW")
		dat <- readxl::read_excel(data.file,sheet=4,skip=4)
	}
	colnames(dat)[which(colnames(dat)=='Site code')] <- 'uwwCode'
	colnames(dat)[which(colnames(dat)=='Site name')] <- 'uwwName'
	colnames(dat)[which(colnames(dat)=='Region name')] <- 'RegionName'
	return(dat)
}


#
###########################################################################
#' Getting locations of the STW sites
#'
#' https://uwwtd.eu/United-Kingdom/download (choose csv under Title Urban Waste Water Treatment plants)
#'
#' @return
#' @export
get_sites_geo <- function(site_uniq) {
	data.file <- system.file("extdata", "UWWTD_United-Kingdom_UrbanWasteWaterTreatmentPlant.csv", package = "publicWW")
	geo <- utils::read.csv(data.file,header=TRUE)
	latlong <- geo[c('uwwCode','uwwLatitude','uwwLongitude')]

	###  merge locations with time invariant variables
	geo <- merge(site_uniq,latlong,by='uwwCode',all.x=TRUE,all.y=FALSE)

	### Setting existing coordinate as lat-long system
	cord.dec <- sp::SpatialPoints(cbind(geo$uwwLongitude, geo$uwwLatitude), proj4string = sp::CRS("+proj=longlat"))

	### Transforming coordinate to UTM using EPSG=32748 for WGS=84, UTM Zone=48M, Southern Hemisphere)
	cord.UTM <- sp::spTransform(cord.dec, sp::CRS("+init=epsg:27700"))@coords
	#par(mfrow=c(1,2))
	#plot(cord.UTM[,1],cord.UTM[,2])
	#plot(cord.dec@coords[,1],cord.dec@coords[,2])

	geo$eastings <- cord.UTM[,1]/1000  # converted to km
	geo$northings <- cord.UTM[,2]/1000
	return(geo)
}

###########################################################################
#' Get England boundary for mesh construction
#'
#'
#'
#' @return
#' @export
get.england.boundary <- function() {
	data.file <- system.file("extdata", "CTRY_DEC_2021_GB_BUC.shp", package = "publicWW")
	shp <- maptools::readShapePoly(data.file)

	id <- grep('England',as.character(shp@data$CTRY21NM),ignore.case=TRUE)
	p <- shp@polygons[[id]]@Polygons

	coords <- NULL
	for (i in 1:length(p)) {
		coords <- rbind(coords,p[[i]]@coords)
	}
	coords <- coords/1000  # convert to KM
	return(coords)
}


###########################################################################
#'  aggregating LSOA covariates to STW catchment using the open source
#'  lookup table https://github.com/alan-turing-institute/ukhsa-turing-rss-wastewater/tree/hoffmann-data/data/wastewater_catchment_areas_public
#'
#' @return
#' @export
LSOA.to.catchment.lookup <- function(){
  #  LSOA-catchment lookup (with STW identifier)
  #x <- RCurl::getURL('https://raw.githubusercontent.com/alan-turing-institute/ukhsa-turing-rss-wastewater/hoffmann-data/data/wastewater_catchment_areas_public/lsoa_catchment_lookup.csv')
  data.file <- system.file("extdata","lsoa_catchment_lookup.txt",package="publicWW")
  lookup <- read.csv(file=data.file)
  #  STW identifier - STW Name/code lookup
  #x <- RCurl::getURL('https://raw.githubusercontent.com/alan-turing-institute/ukhsa-turing-rss-wastewater/hoffmann-data/data/wastewater_catchment_areas_public/waterbase_catchment_lookup.csv')
  data.file <- system.file("extdata","waterbase_catchment_lookup.txt",package="publicWW")
  stw.ids <- read.csv(file=data.file)
  lsoa.catch.lookup <- merge(lookup,stw.ids,by='identifier')
  return(lsoa.catch.lookup)
}


###########################################################################
#'  get LSOA-level covariates from the alan-turing-institute/ukhsa-turing-rss-wastewater repro
#'
#' @return
#' @export
get.LSOA.covariates <- function() {
	#  IMD, population counts by ethnicity (2011 census) and population density (km2, 2019)
	#x <- RCurl::getURL('https://raw.githubusercontent.com/alan-turing-institute/ukhsa-turing-rss-wastewater/main/data/lsoa_covariates.csv')
	#d <- read.csv(text=x)
	data.file <- system.file("extdata","lsoa_covariates.txt",package="publicWW")
	d <- read.csv(file=data.file)
	#  population
	#x <- RCurl::getURL('https://raw.githubusercontent.com/alan-turing-institute/ukhsa-turing-rss-wastewater/main/data/lsoa_pop_estimates.csv')
	data.file <- system.file("extdata","lsoa_pop_estimates.txt",package="publicWW")
	p <- read.csv(file=data.file)
	colnames(p)[grep('16',colnames(p))] <- 'age<=16'
	colnames(p)[grep('75',colnames(p))] <- 'age>=75'

	p['young_prop'] =p['age<=16'] / p['all_ages']
	p['old_prop'] = p['age>=75'] / p['all_ages']
	p['population'] = p['all_ages']

	#  land cover: fraction of urban/industry/vegetation/other
	data.file <- system.file("extdata","lsoa_landcover.csv",package="publicWW")
	landcover <- read.csv(data.file)
	#  keep only LSOAs in England
	landcover <- landcover[grep('E',substr(landcover$LSOA11CD,1,1)),]

	#  merge the two
	cova <- merge(d,p,by=c('lsoa11_nm','lsoa11_cd'),all=FALSE)
	names(cova)[grep('lsoa11_nm',names(cova))] <- 'LSOA11NM'
	names(cova)[grep('lsoa11_cd',names(cova))] <- 'LSOA11CD'

	# covariates from lsoa_covariates
	cova['bame_proportion'] = 1-cova['white']/cova['all_ethnicities']

	# add geographical indicators for LAD/CCG/RGN/CTRY
	lookup <- get.geo.lookups()
	cova <- merge(cova,lookup$lookup,by=c('LSOA11CD','LSOA11NM'),all.x=TRUE,all.y=TRUE)

	# add land cover
	cova <- merge(cova,landcover,by=c('LSOA11CD','LSOA11NM'),all.x=TRUE,all.y=TRUE)

  # a numeric region indicator
	cova$region_index <- add_region_index(cova$RGN21NM)

	return(cova)
}

###########################################################################
#' Construct lookup tables between LSOA11 and various other admin geographies (LAD21, CCG21, RGN21 and CTRY21)
#'
#'
#'
#' @return
#' @export
get.geo.lookups <- function() {
	#  LSOA to CCG to LAD
	data.file <- system.file("extdata", "Lower_Layer_Super_Output_Area_(2011)_to_Clinical_Commissioning_Group_to_Local_Authority_District_(April_2021)_Lookup_in_England.csv", package = "publicWW")
	dat1 <- read.csv(data.file)
	#  Ward to LAD to Region to Country
	data.file <- system.file("extdata", "Ward_to_Local_Authority_District_to_County_to_Region_to_Country_(December_2021)_Lookup_in_United_Kingdom.csv", package = "publicWW")
	dat2 <- read.csv(data.file)
	#  retain only entries in England
	dat2 <- subset(dat2,dat2$CTRY21NM=='England')

	#  add regions to the LSOA lookup
	uregions <- unique(dat2[c('RGN21CD','RGN21NM')])
	nregions <- nrow(uregions)
	RGN21CD <- RGN21NM <- rep('',length(dat1))
	for (ir in 1:nregions) {
		ids <- which(dat2$RGN21CD==uregions$RGN21CD[ir])
		lads <- unique(dat2[ids,c('LAD21CD','LAD21NM')])
		ids <- unlist(sapply(lads$LAD21CD,function(x){which(dat1$LAD21CD==x)}))
		RGN21CD[ids] <- uregions$RGN21CD[ir]
		RGN21NM[ids] <- uregions$RGN21NM[ir]
	}
	dat1$RGN21CD <- RGN21CD
	dat1$RGN21NM <- RGN21NM

	#  add country indicator to the LSOA lookup
	dat1$CTRY21CD <- unique(dat2$CTRY21CD)
	dat1$CTRY21NM <- unique(dat2$CTRY21NM)

	#  code-name lookups
	code.name <- list()
	code.name$LSOA11 <- unique(dat1[c('LSOA11CD','LSOA11NM')])
	code.name$LAD21 <- unique(dat1[c('LAD21CD','LAD21NM')])
	code.name$CCG21 <- unique(dat1[c('CCG21CD','CCG21NM','CCG21CDH')])
	code.name$RGN21 <- unique(dat1[c('RGN21CD','RGN21NM')])
	code.name$CTRY21 <- unique(dat1[c('CTRY21CD','CTRY21NM')])

	#  remove a column named FID in the lookup table
	id <- grep('FID',names(dat1))
	out <- list(lookup=dat1[-c(id)],code.name=code.name)
	return(out)
}

###########################################################################
#' Read the publicly available population-weighted LSOA centroids for England
#'
#'
#' @return
#' @export
read_LSOA_centroids <- function(england.only=TRUE) {
	data.file <- system.file("extdata", "Lower_Layer_Super_Output_Areas_(December_2011)_Population_Weighted_Centroids.csv", package = "publicWW")
	all_c <- E_c <- read.csv(data.file)

	## select only England
	if (england.only) E_c <- all_c[grep('E',all_c$lsoa11cd),]

	## convert easting/northing to km
	E_c$X <- E_c$X/1000
	E_c$Y <- E_c$Y/1000
	## rename some columns
	names(E_c)[grep('X',names(E_c))] <- 'eastings'
	names(E_c)[grep('Y',names(E_c))] <- 'northings'
	return(E_c)
}



###########################################################################
#' A lookup table between time index and first day of the week
#'
#'
#'
#' @return
#' @export
get.starting.day <- function() {
	data.file <- system.file("extdata", "time_index.csv", package = "publicWW")
	d <- read.csv(data.file)
	return(d)
}


###########################################################################
#' Add regional indicator to each LSOA centriod
#'
#'
#' @return
#' @export
add.region.to.lsoas <- function(lsoa.x) {
	# https://data.gov.uk/dataset/ec39697d-e7f4-4419-a146-0b9c9c15ee06/output-area-to-lsoa-to-msoa-to-local-authority-district-december-2017-lookup-with-area-classifications-in-great-britain
	data.file <- system.file("extdata", "Output_Area_to_LSOA_to_MSOA_to_Local_Authority_District_(December_2017)_Lookup_with_Area_Classifications_in_Great_Britain.csv", package = "publicWW")
	tb <- read.csv(data.file)
	tbs <- tb[c('LSOA11NM','LSOA11CD','LAD17CD','LAD17NM','RGN11CD','RGN11NM')]
	out <- merge(lsoa.x,tbs,by=c('LSOA11NM','LSOA11CD'),all.x=TRUE,all.y=FALSE,sort=FALSE)
	return(out)
}
