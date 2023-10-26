# TailCalibration

This repository contains R code to evaluate the calibration of forecasts for extreme outcomes.

## Forecast calibration

Loosely speaking, probabilistic forecasts are _calibrated_ if they align statistically with the corresponding outcomes. Calibration is a necessary property for forecasts to be considered trustworthy.

While several notions of calibration exist, conventional forecast evaluation techniques are ill-equipped to assess forecasts made for extreme outcomes. In this work, we introduce several notions of tail calibration, which allow us to assess whether or not probabilistic forecasts issue reliable forecasts for extreme events.

This repository provides the functionality to implement these different notions of tail calibration in practice. Usage of the repository is thoroughly documented in the vignette.

## Installation and development

This package has not been submitted to CRAN, and can therefore be installed in R using devtools
```r
# install.packages("devtools")
library(devtools)
install_github("sallen12/TailCalibration")
```
The package is still in active development, and the vignette lists several possible extensions that could be implemented. Comments, suggestions, and input are more than welcome.
