## Update

This is an update (to version 0.8.0). Changes to this version can be found below.

## Changes

* Fixed a typo in D_regularized manual
* Added difference between two dependent correlations -function (diff_two_dep_cors) which enables simultaneous estimation and testing of Cohen's q under variable dependency
* Added possibility to use manually constructed regularization and estimation datasets, supplied as a list of two dataframes to "data"-argument in D_regularized
* Improved the output of diff_two_dep_cors and included an argument for missing data ("ML")
* Improved clarity in multivariate sex difference vignette
* Added value_correlation -function for testing and quantifying how ipsatizing values influences associations with other variables
* Added ddsc_sem -function for deconstructing difference score correlation with structural equation modeling

## Test environments
* local R installation, R 4.2.2
* win-builder (devel)
* Linux, using devtools::check_rhub()

## R CMD check results

0 errors | 0 warnings | 0 note
