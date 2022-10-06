## Data sources

The covariate data files contained here were constructed from the following, publicly available data sources:
-  Census 2011 ethnicity by LSOA from [nomis](https://www.nomisweb.co.uk/query/construct/summary.asp?mode=construct&version=0&dataset=1087).
- Population estimates (2019) from [ONS](https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/lowersuperoutputareamidyearpopulationestimates)
- Population density (2019) from [ONS](https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/lowersuperoutputareapopulationdensity)

The above LSOA-level covariates are mapped onto the SWT catchments using the catchment-LSOA lookup obtained from [this public repo](https://github.com/tillahoffmann/wastewater-catchment-areas). The resulting catchment-level covariates are in `catchment_covariates_no_genomic.rds`.

Note that the 21 STWs in the South West are not in [the public repo](https://github.com/tillahoffmann/wastewater-catchment-areas). We approximated their catchments using circles that are centered at the STW locations and include the centroids of at least 10 LSOAs. See the function `publicWW::approximate.catchment.using.circle` for more detail.

