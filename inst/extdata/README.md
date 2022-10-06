## Data sources

### Covariates and catchment-LSOA mapping

The covariate data files contained here were constructed from the following, publicly available data sources:
-  Census 2011 ethnicity by LSOA from [nomis](https://www.nomisweb.co.uk/query/construct/summary.asp?mode=construct&version=0&dataset=1087).
- Population estimates (2019) from [ONS](https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/lowersuperoutputareamidyearpopulationestimates)
- Population density (2019) from [ONS](https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/lowersuperoutputareapopulationdensity)

The above LSOA-level covariates are mapped onto the SWT catchments using the catchment-LSOA lookup obtained from [this public repo](https://github.com/tillahoffmann/wastewater-catchment-areas). The resulting catchment-level covariates are in `catchment_covariates_no_genomic.rds`.

Note that the 21 STWs in the South West are not in [the public repo](https://github.com/tillahoffmann/wastewater-catchment-areas). We approximated their catchments using circles that are centered at the STW locations and include the centroids of at least 10 LSOAs. See the function `publicWW::approximate.catchment.using.circle` for more detail.


### Data on geography
- [The Lower Layer Super Output Areas (December 2011) Population Weighted Centroids](https://geoportal.statistics.gov.uk/datasets/ons::lower-layer-super-output-areas-december-2011-population-weighted-centroids/explore?location=52.905447%2C-2.000000%2C7.21&showTable=true)
- [Lower Layer Super Output Area (2011) to Clinical Commissioning Group to Local Authority District (April 2021) Lookup in England](https://geoportal.statistics.gov.uk/datasets/ons::lower-layer-super-output-area-2011-to-clinical-commissioning-group-to-local-authority-district-april-2021-lookup-in-england/explore)
- [Local Authority Districts (May 2021) UK BFE boundary for mapping](https://geoportal.statistics.gov.uk/datasets/1119a90ec5f343678f044374392e6bda_0/explore?location=55.215451%2C-3.313875%2C6.48)
- [(Ward to) Local Authority District to County to Region to Country (December 2021) Lookup in United Kingdom](https://geoportal.statistics.gov.uk/datasets/7244b2542e4a43c0ae528ec52fea2cf5/explore)

The model predicts at the LSOA (2011 definition, the most up-to-date version available online) population weighted centroids. The LSOA-predictions can then be aggregated to one of the following:
- Lower Tier Local Authority
- Clinical Commissioning Group 
- Region
- England

All geographies above are based on their 2021 definitions.
