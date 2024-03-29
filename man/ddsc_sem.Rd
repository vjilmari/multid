% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ddsc_sem.R
\name{ddsc_sem}
\alias{ddsc_sem}
\title{Deconstructing difference score correlation with structural equation modeling}
\usage{
ddsc_sem(
  data,
  x,
  y1,
  y2,
  center_yvars = FALSE,
  covariates = NULL,
  estimator = "ML",
  level = 0.95,
  sampling.weights = NULL,
  q_sesoi = 0,
  min_cross_over_point_location = 0,
  boot_ci = FALSE,
  boot_n = 5000,
  boot_ci_type = "perc"
)
}
\arguments{
\item{data}{A data frame.}

\item{x}{Character string. Variable name of independent variable.}

\item{y1}{Character string. Variable name of first component score of difference score.}

\item{y2}{Character string. Variable name of second component score of difference score.}

\item{center_yvars}{Logical. Should y1 and y2 be centered around their grand mean? (Default FALSE)}

\item{covariates}{Character string or vector. Variable names of covariates (Default NULL).}

\item{estimator}{Character string. Estimator used in SEM (Default "ML").}

\item{level}{Numeric. The confidence level required for the result output (Default .95)}

\item{sampling.weights}{Character string. Name of sampling weights variable.}

\item{q_sesoi}{Numeric. The smallest effect size of interest for Cohen's q estimates (Default 0; See Lakens et al. 2018).}

\item{min_cross_over_point_location}{Numeric. Z-score for the minimal slope cross-over point of interest (Default 0).}

\item{boot_ci}{Logical. Calculate confidence intervals based on bootstrap (Default FALSE).}

\item{boot_n}{Numeric. How many bootstrap redraws (Default 5000).}

\item{boot_ci_type}{If bootstrapping was used, the type of interval required. The value should be one of "norm", "basic", "perc" (default), or "bca.simple".}
}
\value{
\item{descriptives}{Means, standard deviations, and intercorrelations.}
\item{parameter_estimates}{Parameter estimates from the structural equation model.}
\item{variance_test}{Variances and covariances of component scores.}
\item{data}{Data frame with original and scaled variables used in SEM.}
\item{results}{Summary of key results.}
}
\description{
Deconstructs a bivariate association between x and a difference score y1-y2 with SEM.
A difference score correlation is indicative that slopes for y1 as function of x and y2 as function of x are non-parallel.
Deconstructing the bivariate association to these slopes allows for understanding the pattern and magnitude of this non-parallelism.
}
\examples{
\dontrun{
set.seed(342356)
d <- data.frame(
  y1 = rnorm(50),
  y2 = rnorm(50),
  x = rnorm(50)
)
ddsc_sem(
  data = d, y1 = "y1", y2 = "y2",
  x = "x",
  q_sesoi = 0.20,
  min_cross_over_point_location = 1
)$results
}
}
\references{
Edwards, J. R. (1995). Alternatives to Difference Scores as Dependent Variables in the Study of Congruence in Organizational Research. Organizational Behavior and Human Decision Processes, 64(3), 307–324.

Lakens, D., Scheel, A. M., & Isager, P. M. (2018). Equivalence Testing for Psychological Research: A Tutorial. Advances in Methods and Practices in Psychological Science, 1(2), 259–269. https://doi.org/10.1177/2515245918770963
}
