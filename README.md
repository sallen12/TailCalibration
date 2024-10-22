# TailCalibration

This repository contains R code to evaluate the tail calibration of probabilistic forecasts, as introduced in the following pre-print.

> Allen, S., Koh, J., Segers, J. and Ziegel, J. (2024). 
> Tail calibration of probabilistic forecasts.
> ArXiv preprint.
> [arxiv.org/abs/2407.03167](https://arxiv.org/abs/2407.03167)

The results in this paper can be reproduced by sourcing the files in the scripts folder.

## Forecast calibration

Loosely speaking, probabilistic forecasts are _calibrated_ if they align statistically with the corresponding outcomes. Calibration is a necessary property for forecasts to be considered trustworthy. 

While several notions of calibration have been proposed, it is most common to evaluate the _probabilistic calibration_ of forecasts, which corresponds to uniformity of forecast probability integral transform (PIT) values, or the _marginal calibration_, which corresponds to equality between the average forecast distribution and the unconditional distribution of the outcomes. However, conventional forecast evaluation techniques are ill-equipped to assess forecasts made for extreme outcomes. Notions of tail calibration have therefore been proposed to assess whether or not probabilistic forecasts issue reliable forecasts for the occurrence and severity of extreme events.

This repository provides the functionality to evaluate different notions of tail calibration in practice.

## Installation and development

This package has not been submitted to CRAN, but can be installed in R using devtools
```r
# install.packages("devtools")
library(devtools)
install_github("sallen12/TailCalibration")
```
The package is still in active development. Comments, suggestions, and input are more than welcome.


