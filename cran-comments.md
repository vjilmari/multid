## Update

This is an update (to version 1.0.0). Changes to this version can be found below.

## Changes

* Added a possibility to run ddsc_ml with just two observations per upper-level unit
* Added a possibility to obtain bootstrap estimates and percentile confidence intervals for non-scaled parameter estimates in ddsc_ml results-table
* Added confidence intervals for ddsc_ml results table
* Added a possibility to bootstrap in ddsc_sem
* ml_dadas and sem_dadas deprecated (superceded by ddsc_ml and ddsc_sem)
* Added plot_ddsc function for directly plotting ddsc_sem results
* Removed ml_dadas and sem_dadas from README examples. Replaced with ddsc_ml and ddsc_sem

## Test environments
* local R installation, R 4.3.2
* win-builder (devel)
* Linux, using devtools::check_rhub()

## R CMD check results

0 errors | 0 warnings | 0 note
