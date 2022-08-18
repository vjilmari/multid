## Update

This is an update (to version 0.7.0). Changes to this version can be found below.

## Changes

* Switched from sample to dplyr::sample_n for bootstrap example in the multivariate sex difference vignette
* Minor changes to style and text in the multivariate sex difference vignette
* Additional descriptive statistics to reliability_dms output
* Allow vpc_at for models with no intercept-slope covariation (conditional level-2 variances are same for all requested level-1 values)
* Added function qcc for quantile correlation coefficient
* Updated README

## Test environments
* local R installation, R 4.2.1
* win-builder (devel)
* Linux, using devtools::check_rhub()

## R CMD check results

0 errors | 0 warnings | 0 note
