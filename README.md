
<!-- README.md is generated from README.Rmd. Please edit that file -->

# multid

<!-- badges: start -->

[![Travis build
status](https://app.travis-ci.com/vjilmari/multid.svg?branch=main)](https://app.travis-ci.com/vjilmari/multid/)
<!-- badges: end -->

The goal of multid is to provide tools for regularized measurement of
multivariate differences between two groups (e.g., sex differences).
Regularization via logistic regression variants enables inclusion of
large number of correlated variables in the multivariate set while
providing k-fold cross-validation and regularization to avoid
overfitting.

Studies in which the method has been used:

1.  [Lönnqvist & Ilmarinen (2021). Using a Continuous Measure of
    Genderedness to Assess Sex Differences in the Attitudes of the
    Political Elite. Political
    Behavior.](https://doi.org/10.1007/s11109-021-09681-2)
2.  [Ilmarinen et al. (2021). Is There a g-factor of Genderedness? Using
    a Continuous Measure of Genderedness to Assess Sex Differences in
    Personality, Values, Cognitive Ability, School Grades, and
    Educational Track. Manuscript in
    review.](https://doi.org/10.31234/osf.io/j59bs)

## Installation

You can install the released version of multid from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("multid")
```

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("vjilmari/multid")
```

## Example

This is a basic example which shows you how to measure standardized
multivariate (both Sepal and Petal dimensions, four variables in total)
distance between setosa and versicolor Species in iris dataset. Similar
example is provided with artificially generated multi group data which
are used as separate data folds in the procedure.

``` r
library(multid)
set.seed(91237)

D_regularized(
  data = iris[iris$Species == "setosa" | iris$Species == "versicolor", ],
  mv.vars = c(
    "Sepal.Length", "Sepal.Width",
    "Petal.Length", "Petal.Width"
  ),
  group.var = "Species",
  group.values = c("setosa", "versicolor")
)$D
#>      n.setosa n.versicolor m.setosa m.versicolor sd.setosa sd.versicolor
#> [1,]       50           50 8.348364       -9.379  1.386586      2.247504
#>      pooled.sd     diff        D
#> [1,]  1.867337 17.72736 9.493393
# Use different partitions of data for regularization and estimation
D_regularized(
  data = iris[iris$Species == "setosa" |
    iris$Species == "versicolor", ],
  mv.vars = c(
    "Sepal.Length", "Sepal.Width",
    "Petal.Length", "Petal.Width"
  ),
  group.var = "Species",
  group.values = c("setosa", "versicolor"),
  size = 35,
  out = TRUE
)$D
#>      n.setosa n.versicolor m.setosa m.versicolor sd.setosa sd.versicolor
#> [1,]       15           15 7.614128    -8.973335   1.79109      2.353035
#>      pooled.sd     diff        D
#> [1,]  2.091026 16.58746 7.932692
# multigroup model where the group variable is defined as fold.var, and output is produced separately for each group/fold. Groups are also used as folds in the k-fold cross-validation procedure.

# separate sample folds
# generate data for 10 groups
set.seed(34246)
n1 <- 100
n2 <- 10
d <-
  data.frame(
    sex = sample(c("male", "female"), n1 * n2, replace = TRUE),
    fold = sample(x = LETTERS[1:n2], size = n1 * n2, replace = TRUE),
    x1 = rnorm(n1 * n2),
    x2 = rnorm(n1 * n2),
    x3 = rnorm(n1 * n2)
  )
#'
# Fit and predict with same data
round(D_regularized(
  data = d,
  mv.vars = c("x1", "x2", "x3"),
  group.var = "sex",
  group.values = c("female", "male"),
  fold.var = "fold",
  fold = TRUE,
  rename.output = TRUE
)$D,2)
#>   n.female n.male m.female m.male sd.female sd.male pooled.sd  diff     D
#> A       53     48     0.01  -0.02      0.07    0.06      0.07  0.04  0.53
#> B       56     53     0.00  -0.02      0.07    0.06      0.07  0.02  0.27
#> C       39     56    -0.02  -0.01      0.08    0.06      0.07 -0.01 -0.18
#> D       53     58     0.00  -0.01      0.06    0.06      0.06  0.01  0.17
#> E       44     52     0.00   0.01      0.06    0.07      0.06 -0.01 -0.14
#> F       60     31    -0.01  -0.01      0.06    0.06      0.06 -0.01 -0.10
#> G       56     51    -0.01  -0.03      0.06    0.05      0.06  0.01  0.20
#> H       47     51     0.00  -0.01      0.07    0.07      0.07  0.01  0.14
#> I       42     45     0.02  -0.01      0.07    0.07      0.07  0.03  0.48
#> J       48     57     0.00  -0.01      0.05    0.07      0.06  0.01  0.10
#>   pooled.sd.female pooled.sd.male pooled.sd.total d.sd.total
#> A             0.07           0.06            0.07       0.56
#> B             0.07           0.06            0.07       0.26
#> C             0.07           0.06            0.07      -0.19
#> D             0.07           0.06            0.07       0.16
#> E             0.07           0.06            0.07      -0.14
#> F             0.07           0.06            0.07      -0.10
#> G             0.07           0.06            0.07       0.18
#> H             0.07           0.06            0.07       0.15
#> I             0.07           0.06            0.07       0.49
#> J             0.07           0.06            0.07       0.10
#'
# Different partitions for regularization and estimation for each data fold
round(D_regularized(
  data = d,
  mv.vars = c("x1", "x2", "x3"),
  group.var = "sex",
  group.values = c("female", "male"),
  fold.var = "fold",
  size = 17,
  out = TRUE,
  fold = TRUE,
  rename.output = TRUE
)$D,2)
#>   n.female n.male m.female m.male sd.female sd.male pooled.sd  diff     D
#> A       36     31     0.00  -0.01      0.18    0.15      0.17  0.00  0.03
#> B       39     36    -0.04  -0.04      0.16    0.15      0.15 -0.01 -0.04
#> C       22     39    -0.07   0.01      0.17    0.13      0.15 -0.08 -0.56
#> D       36     41     0.01   0.01      0.14    0.17      0.15  0.00 -0.01
#> E       27     35    -0.07   0.00      0.13    0.11      0.12 -0.07 -0.60
#> F       43     14    -0.04   0.04      0.17    0.17      0.17 -0.08 -0.47
#> G       39     34    -0.03  -0.04      0.12    0.15      0.14  0.01  0.07
#> H       30     34     0.04   0.01      0.16    0.15      0.16  0.03  0.18
#> I       25     28     0.02  -0.03      0.14    0.14      0.14  0.05  0.35
#> J       31     40     0.00   0.00      0.15    0.15      0.15  0.00  0.02
#>   pooled.sd.female pooled.sd.male pooled.sd.total d.sd.total
#> A             0.15           0.15            0.15       0.03
#> B             0.15           0.15            0.15      -0.04
#> C             0.15           0.15            0.15      -0.54
#> D             0.15           0.15            0.15      -0.01
#> E             0.15           0.15            0.15      -0.48
#> F             0.15           0.15            0.15      -0.54
#> G             0.15           0.15            0.15       0.06
#> H             0.15           0.15            0.15       0.19
#> I             0.15           0.15            0.15       0.32
#> J             0.15           0.15            0.15       0.02
```
