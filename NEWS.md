# multid 0.4.0.9000

* Development version
* Added option to obtain scaled estimates in ml_dadas. Scaling is done for both difference score components and the difference scores, based on random intercept SDs and random slope SD, respectively, in a reduced model without the predictor and the interaction between predictor and moderator
* Added option to test for random effect covariation with likelihood ratio test in ml_dadas from a reduced model without the predictor and the interaction between predictor and moderator
* Added option to include covariates in sem_dadas
* Added variance test with sem in sem_dadas

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
