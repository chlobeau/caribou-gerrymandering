# How boundary delineation shapes population trend estimates: a boreal caribou case study

Chloe Beaupré, Allicia Kelly, James Hodson, Ashley A.D. McLaren, Ève Lamontagne, Eliezer Gurarie. 2026. Manuscript in review. [DOI when published]

## Overview

This repository holds the R code for evaluating how boundary schemes affect estimates of boreal caribou (*Rangifer tarandus caribou*) survival, recruitment, and population growth rate (lambda) in the southern Northwest Territories, Canada. 
We compare six boundary schemes (StudyArea, StudyArea2, Landscape, RangePlan, Harvest, Admin) using three analyses:

1. **Polygon fidelity**, the proportion of each animal's annual GPS locations that fall inside its assigned polygon
2. **Randomization tests** to assess spatial structure of vital rates
3. **Temporal structure in vital rates**, the share of variation attributable to year versus delienation

## Data availability

Code for all analyses and demographic data, with spatial coordinates removed and polygon delineation assignments retained, are available in this repository. 
Boreal caribou movement data are owned by the Government of the Northwest Territories. Researchers seeking access should contact the Department of Environment and Climate Change, which can grant access to the dataset on Movebank.

## Requirements

R (v4.4.0 or later) with the packages `sf`, `glmmTMB`, `survival`, `emmeans`, `performance`, `plyr`, `dplyr`, `tidyr`, `TuktuTools`, and `taigatools`. 

## Installing taigatools

taigatools is not on CRAN. Install it from R-universe:

    install.packages("taigatools",
                     repos = c("https://chlobeau.r-universe.dev",
                               "https://cloud.r-project.org"))

## Contact

Chloe Beaupre, cbeaupre@esf.edu
