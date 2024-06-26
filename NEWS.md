# multid 1.0.0.9000

* Development version
* Added coef_text_size to plot_ddsc
* Added x_scale and y_scale for plot_ddsc
* Fixed problem in checking for multiple observation data in ml_ddsc
* Defined sign difference in quantile correlation coefficient (qcc) resulting as rho_tau = zero, not NaN as previously

# multid 1.0.0

* Added a possibility to run ddsc_ml with just two observations per upper-level unit
* Added a possibility to obtain bootstrap estimates and percentile confidence intervals for non-scaled parameter estimates in ddsc_ml results-table
* Added confidence intervals for ddsc_ml results table
* Added a possibility to bootstrap in ddsc_sem
* ml_dadas and sem_dadas deprecated (superceded by ddsc_ml and ddsc_sem)
* Added plot_ddsc function for directly plotting ddsc_sem results
* Removed ml_dadas and sem_dadas from README examples. Replaced with ddsc_ml and ddsc_sem

# multid 0.9.0

* Renamed variance_test output in ddsc_sem
* Added ddsc_ml -function for deconstructing difference score correlation with multi-level modeling 

# multid 0.8.0

* Fixed a typo in D_regularized manual
* Added difference between two dependent correlations -function (diff_two_dep_cors) which enables simultaneous estimation and testing of Cohen's q under variable dependency
* Added possibility to use manually constructed regularization and estimation datasets, supplied as a list of two dataframes to "data"-argument in D_regularized
* Improved the output of diff_two_dep_cors and included an argument for missing data ("ML")
* Improved clarity in multivariate sex difference vignette
* Added value_correlation -function for testing and quantifying how ipsatizing values influences associations with other variables
* Added ddsc_sem -function for deconstructing difference score correlation with structural equation modeling 

# multid 0.7.1

* Added na.rm to qcc bootstrap summary over tau-values
* Added main and interaction effects, and comparison of their absolute magnitudes to ml_dadas and sem_dadas outputs
* Added moderator/intercept difference estimates for dadas-functions
* Added abs_coef_diff_test in sem_dadas and ml_dadas to enable tests for slope difference that is not against null but a different numeric value
* Added an estimate of scaled difference of slopes for ml_dadas and a derived estimate of component correlation
* Fixed URLs and output in README

# multid 0.7.0

* Switched from sample to dplyr::sample_n for bootstrap example in the multivariate sex difference vignette
* Minor changes to style and text in the multivariate sex difference vignette
* Additional descriptive statistics to reliability_dms output
* Allow vpc_at for models with no intercept-slope covariation (conditional level-2 variances are same for all requested level-1 values)
* Added function qcc for quantile correlation coefficient
* Updated README

# multid 0.6.0

* Bug fixed in D_regularized_out
* Addend.data argument added to all D_regularized -functions. In _fold -functions, test-partition of the data is appended, else the entire data frame is added.
* Include ICC2 (group-mean reliability) for vpc_at. Enables calculation of sub-group mean-level reliabilities, in case the "at" is a group
* Include reliability_dms that calculates difference score reliability coefficient for data that is difference score between two mean values across some upper-level units (e.g., sex differences across countries)
* Vignette on estimation of multivariate sex differences with multid added

# multid 0.5.0

* Added option to obtain scaled estimates in ml_dadas. Scaling is done for both difference score components and the difference scores, based on random intercept SDs and random slope SD, respectively, in a reduced model without the predictor and the interaction between predictor and moderator
* Added option to test for random effect covariation with likelihood ratio test in ml_dadas from a reduced model without the predictor and the interaction between predictor and moderator
* Added option to include covariates in sem_dadas
* Added variance test with sem in sem_dadas
* Added variance test via parametric bootstrap in ml_dadas
* Added cvv_manual -function for calculation of coefficients of variance variation from manually inputted sample sizes and variances of multiple variables

# multid 0.4.0

* Replaced two-sided tests in sem_dadas for absolute parameters with one-sided tests
* Added three variants of coefficient of variance variation in cvv -function (CVV=coefficient of variance variation, SVH=standardized variance heterogeneity, and VR=variance ration between the largest and the smallest variance)
* Added vpc_at -function for calculation of intercept variances and variance partition coefficients (VPCs) at selected values of level-1 predictors in two-level models

# multid 0.3.0

* Added sem_dadas and ml_dadas functions for predicting algebraic difference scores in structural equation (sem_dadas) and multilevel model (ml_dadas). DADAS acronym follows from the joint hypothesis test of Difference between Absolute Differences and Absolute Sums between (regression coefficients on difference score components).

# multid 0.2.0

* Fixed bug: renaming output in D_regularized_fold functions
* Fixed joining data frames by fold.var in D_regularized_fold functions
* Added statistical inference to d_pooled_sd
* Added probabilities of correct classification option to D_regularized_out and D_regularized_fold_out
* Added area under the receiver operating characteristics to D_regularized_out and D_regularized_fold_out
* Added probability classification tables to D_regularized_out and D_regularized_fold_out
* Added more examples to README file

# multid 0.1.0

* First submission of multid package
