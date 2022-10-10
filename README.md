# publicWW

## Overview

The publicWW package implements a spatio-temporal modelling framework that produces probabilistic estimates of SARS-CoV-2 viral concentration in wastewater at high spatio-temporal resolution across England. The modelling framework utilises [the publicly available wastewater data](https://www.gov.uk/government/publications/monitoring-of-sars-cov-2-rna-in-england-wastewater-monthly-statistics-15-july-2020-to-30-march-2022) collected as [part of the Environmental Monitoring for Health Protection (EMHP) wastewater monitoring programme](https://www.gov.uk/government/publications/wastewater-testing-coverage-data-for-19-may-2021-emhp-programme/wastewater-testing-coverage-data-for-the-environmental-monitoring-for-health-protection-emhp-programme). The functionality of the package includes data management, model fitting, model evaluation and output visualisation. All data files are in the [inst/extdata/ directory](inst/extdata/README.md). Outputs of this modelling framework are visualised via [a dynamic and interactive dashboard](https://b-rowlingson.gitlab.io/wwatlas/).


##  Installing the package

You first need to obtain a personal access token (PAT) by following the instruction via this link:  https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token.

Enter the token as a string to the object `PAT` then run the following code to install the package:

```R
PAT <- " " # enter your PAT token here
devtools::install_github("gqlNU/publicWW", auth_token = PAT)
```

##  Steps in running the framework

### 1. Modelling the wastewater data using INLA

We use [INLA](https://www.r-inla.org/) to fit a spatio-temporal model to the wastewater data observed across the 303 sewage treatment works (STWs) in England over the period from 1 June, 2021 to 30 March 2022. The script to fit the model in INLA can be found via

[inst/scripts/fit_model.R](inst/scripts/fit_model.R)

The top of the script has a `USER INPUT` section where some of the parameters of the model can be changed as described in the script comments. By default the script includes the following covariates in the model

- IMD_score - Index of Multiple Deprivation for 2019
- old_prop - age structure of the population, defines as percentage of the adults older than 75 in the population (ONS 2019)
- young_prop - age structure of the population, defines as percentage of the young population aging less than 16 years old (ONS 2019)
- bame_proportion - Black, Asian and Minority Ethnic (BAME) proportion in each area (2011 Census)
- population_density - Estimated by ONS 2019.
- [f_industrial](https://land.copernicus.eu/pan-european/corine-land-cover)

:warning: There are also two genomic covariates included in the full model but the data are not included with this package as they cannot be released so the script it set to run without them in the model (on line 20 we set `include_genomic_covars <- FALSE`). :warning:

The resulting model fit will be saved in an `RData` file. This file also saves the auxiliary variables needed to run the predictions (such as the covariate list, the setting for regional effect and the SPDE mesh).

### 2. Making predictions at the LSOA level

The following script takes the model fit from the previous step to predict weekly viral concentration levels at the population-weighted centroids of all the 32844 LSOAs in England:

[inst/scripts/predict_lsoas.R](inst/scripts/predict_lsoas.R)

#### A note on memory:
In `predict_lsoas.R` the fitted `RData` file is loaded. Due to memory constraints predictions will, by default, be made in batches of iteration for all LSOAs.

The results of the predictions are saved in an `outputs` directory. These are two prediction outputs:

1. The full iterations of the batch.
2. The summary statistics for that batch.

:warning: The reason for the batches is to solve memory issues, if you are still facing this problems is best to reduce the size of the batches of LSOAs on the predictions or the simulations (`nsims`). The time to run all LSOAs sequentially is about 8 hours and it requires at least 8.2 GB of hard disk space to store the predictions at iteration level. :warning:

### 3. Aggregating LSOA-level weekly predictions to other spatial resolution levels
The LSOA-level weekly concentrations can be aggregated to other adminstrative levels in England, including
- the Lower Tier Local Authority (LTLA) level
- the Clinical commissioning groups (CCGs) level
- the regional level
- the national level

This spatial aggregation is carried out by using the following script
[inst/scripts/aggregate_LSOApred.R](inst/scripts/aggregate_LSOApred.R)

## Visualising predictions

The links below provide the iteration-level predictions at LTLA, CCG, region and England level.

- [LTLA predictions](https://livenorthumbriaac-my.sharepoint.com/:u:/g/personal/guangquan_li_northumbria_ac_uk/EQQWWWyinQFCvqCD9Z8E2iABk9L0a6AG8tqUDfmazvMi2Q?e=caaQpe)
- [CCG predictions](https://livenorthumbriaac-my.sharepoint.com/:u:/g/personal/guangquan_li_northumbria_ac_uk/ER_-A7oaHVZOpJxZFtv8g08BghnzAgPYzQBWIVG3idKj8g?e=0hkWDI)
- [Region predictions](https://livenorthumbriaac-my.sharepoint.com/:u:/g/personal/guangquan_li_northumbria_ac_uk/EawmudokolNCjMLDzmvyzmEBUYpkhYrO4Ih5ycp1eR9oig?e=KKC2P4)
- [Country predictions](https://livenorthumbriaac-my.sharepoint.com/:u:/g/personal/guangquan_li_northumbria_ac_uk/EcM65OsrhJdEqGrZSe98VG0Bxq0fkaaXgcYgt_K9oZmAvQ)

- [The England LTLA shapefile](https://livenorthumbriaac-my.sharepoint.com/:u:/g/personal/guangquan_li_northumbria_ac_uk/EWkFrGmxWGFNreJX8EZ34koBKiheoJKJtxuGWM3xqXKotA?e=hkCWUT)

Download the files above and use the following script to produce the figures in the paper
[inst/scripts/results4paper_onGit.R](inst/scripts/results4paper_onGit.R)



## Model evaluation via cross validation

A 10-fold cross validation can be run using the following script:

[inst/scripts/run_cross_validation.R](inst/scripts/run_cross_validation.R)

Similarly to the the fitting section, the script is ready to be run with the full dataset (1 June 2021 to 30 March 2022) and using all the covariates, but can be modified. Each fold of out-of-sample sites has been randomly selected without replacement. The list of sites to hold-out per fold is found in [inst/extdata/cv_sites_10folds_30Mar2022.csv](inst/extdata/cv_sites_10folds_30Mar2022.csv). It is also possible to create new folds if necessary.

### Other scripts

There are a few other scripts that are used to explore genomics variables, aggregate batches of predictions, aggregate predictions to the LTLA level in the directory [inst/scripts/](inst/scripts/).
