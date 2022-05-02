## Update

This is an update (version 0.5.0). Changes to this version can be found below.

## Changes

* Added option to obtain scaled estimates in ml_dadas. Scaling is done for both difference score components and the difference scores, based on random intercept SDs and random slope SD, respectively, in a reduced model without the predictor and the interaction between predictor and moderator
* Added option to test for random effect covariation with likelihood ratio test in ml_dadas from a reduced model without the predictor and the interaction between predictor and moderator
* Added option to include covariates in sem_dadas
* Added variance test with sem in sem_dadas
* Added variance test via parametric bootstrap in ml_dadas
* Added cvv_manual -function for calculation of coefficients of variance variation from manually inputted sample sizes and variances of multiple variables

## Test environments
* local R installation, R 4.2.0
* win-builder (devel)
* Linux, using devtools::check_rhub()

## R CMD check results

0 errors | 0 warnings | 0 note
