## Update

This is an update (to version 0.7.1). Changes to this version can be found below.

## Changes

* Added na.rm to qcc bootstrap summary over tau-values
* Added main and interaction effects, and comparison of their absolute magnitudes to ml_dadas and sem_dadas outputs
* Added moderator/intercept difference estimates for dadas-functions
* Added abs_coef_diff_test in sem_dadas and ml_dadas to enable tests for slope difference that is not against null but a different numeric value
* Added an estimate of scaled difference of slopes for ml_dadas and a derived estimate of component correlation
* Fixed URLs and output in README

## Test environments
* local R installation, R 4.2.2
* win-builder (devel)
* Linux, using devtools::check_rhub()

## R CMD check results

0 errors | 0 warnings | 0 note
