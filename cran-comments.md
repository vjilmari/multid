## Update

This is an update (to version 1.0.1). Changes to this version can be found below.

## Changes

* Added coef_text_size to plot_ddsc
* Added x_scale and y_scale for plot_ddsc
* Fixed problem in checking for multiple observation data in ml_ddsc
* Defined sign difference in quantile correlation coefficient (qcc) resulting as rho_tau = zero, not NaN as previously
* Updated description and roxygenNote version to DESCRIPTION
* Fixed several typos and minor problems (e.g, location of deprecated_functions.R)
* Fixed a note related to "no visible binding for global variable" in D_regularized_out dplyr-pipes

## Test environments
* local R installation, R 4.5.1
* win-builder (devel)
* linux, m1-san, macos, macos-arm64, windows, using rhub::rhub_check()

## R CMD check results

0 errors | 0 warnings | 0 note
