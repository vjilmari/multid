% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vpc_at.R
\name{vpc_at}
\alias{vpc_at}
\title{Variance partition coefficient calculated at different level-1 values}
\usage{
vpc_at(model, lvl1.var, lvl1.values)
}
\arguments{
\item{model}{Two-level model fitted with lme4. Must include random intercept, slope, and their covariation.}

\item{lvl1.var}{Character string. Level 1 variable name to which random slope is also estimated.}

\item{lvl1.values}{Level 1 variable values.}
}
\value{
Data frame of level 2 variance and std.dev. estimates at level 1 variable values, respective VPCs (ICC1s) and group-mean reliabilities (ICC2s) (Bliese, 2000).
}
\description{
Calculates variance estimates (level-2 Intercept variance) and variance partition coefficients (i.e., intra-class correlation) at selected values of predictor values in two-level linear models with random effects (intercept, slope, and their covariation).
}
\examples{
fit <- lme4::lmer(Sepal.Length ~ Petal.Length +
  (Petal.Length | Species),
data = iris
)

lvl1.values <-
  c(
    mean(iris$Petal.Length) - stats::sd(iris$Petal.Length),
    mean(iris$Petal.Length),
    mean(iris$Petal.Length) + stats::sd(iris$Petal.Length)
  )

vpc_at(
  model = fit,
  lvl1.var = "Petal.Length",
  lvl1.values = lvl1.values
)
}
\references{
Goldstein, H., Browne, W., & Rasbash, J. (2002). Partitioning Variation in Multilevel Models. Understanding Statistics, 1(4), 223–231. https://doi.org/10.1207/S15328031US0104_02

Bliese, P. D. (2000). Within-group agreement, non-independence, and reliability: Implications for data aggregation and analysis. In K. J. Klein & S. W. J. Kozlowski (Eds.), Multilevel theory, research, and methods in organizations: Foundations, extensions, and new directions (pp. 349–381). Jossey-Bass.
}
