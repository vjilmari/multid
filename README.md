
<!-- README.md is generated from README.Rmd. Please edit that file -->

# multid

<!-- badges: start -->
<!-- badges: end -->

The goal of multid is to provide tools for regularized measurement of
multivariate differences between two groups (e.g., sex differences).
Regularization via logistic regression variants enables inclusion of
large number of correlated variables in the multivariate set while
avoiding overfitting.

<!--  ## Installation

You can install the released version of multid from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("multid")
```-->

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("vjilmari/multid")
```

## Example

This is a basic example which shows you how to measure standardized
multivariate (both Sepal and Petal dimensions, four variables in total)
distance between setosa and versicolor Species in iris dataset:

``` r
library(multid)
set.seed(91237)

D_regularized(
   data = iris[iris$Species == "setosa" | iris$Species == "versicolor", ],
   mv.vars = c("Sepal.Length", "Sepal.Width",
               "Petal.Length", "Petal.Width"),
   group.var = "Species",
   group.values = c("setosa", "versicolor")
 )$D
#>      n.setosa n.versicolor m.setosa m.versicolor sd.setosa sd.versicolor
#> [1,]       50           50 8.348364       -9.379  1.386586      2.247504
#>      pooled.sd     diff        D
#> [1,]  1.867337 17.72736 9.493393
D_regularized(
  data = iris[iris$Species=="setosa" |
                iris$Species=="versicolor",],
  mv.vars = c("Sepal.Length","Sepal.Width",
              "Petal.Length","Petal.Width"),
  group.var="Species",
  group.values=c("setosa","versicolor"),
  size=35,
  out=TRUE
)$D
#>      n.setosa n.versicolor m.setosa m.versicolor sd.setosa sd.versicolor
#> [1,]       15           15 7.614128    -8.973335   1.79109      2.353035
#>      pooled.sd     diff        D
#> [1,]  2.091026 16.58746 7.932692
```
