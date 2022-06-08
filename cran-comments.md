## Update

This is an update (version 0.5.0.9000). Changes to this version can be found below.

## Changes

* Bug fixed in D_regularized_out
* Addend.data argument added to all D_regularized -functions. In _fold -functions, test-partition of the data is appended, else the entire data frame is added.
* Include ICC2 (group-mean reliability) for vpc_at. Enables calculation of sub-group mean-level reliabilities, in case the "at" is a group
* Include reliability_dms that calculates difference score reliability coefficient for data that is difference score between two mean values across some upper-level units (e.g., sex differences across countries)
* Vignette on estimation of multivariate sex differences with multid addede

## Test environments
* local R installation, R 4.2.0
* win-builder (devel)
* Linux, using devtools::check_rhub()

## R CMD check results

0 errors | 0 warnings | 0 note
