
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

Predictive approach as implemented with regularized methods also allows
for examination of group-membership probabilities and their
distributions across individuals. In the context of statistical
predictions of sex, these distributions are an updated variant to
gender-typicality distributions used in gender diagnosticity methodology
[(Lippa & Connelly, 1990)](https://doi.org/10.1037/0022-3514.59.5.1051).

Studies in which these methods have been used:

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

## Examples

### Single sample with two groups

This example shows how to measure standardized multivariate (both Sepal
and Petal dimensions, four variables in total) distance between setosa
and versicolor Species in iris dataset.

``` r
library(multid)
set.seed(91237)

D.iris<-
  D_regularized(
  data = iris[iris$Species == "setosa" | iris$Species == "versicolor", ],
  mv.vars = c(
    "Sepal.Length", "Sepal.Width",
    "Petal.Length", "Petal.Width"
  ),
  group.var = "Species",
  group.values = c("setosa", "versicolor")
)

round(D.iris$D,2)
#>      n.setosa n.versicolor m.setosa m.versicolor sd.setosa sd.versicolor
#> [1,]       50           50     8.35        -9.38      1.39          2.25
#>      pooled.sd  diff    D
#> [1,]      1.87 17.73 9.49
# Use different partitions of data for regularization and estimation
D.iris_out<-
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
  out = TRUE,
  pred.prob = TRUE,
  prob.cutoffs = seq(0,1,0.25)
  
)

# print group differences (D)
round(D.iris_out$D,2)
#>      n.setosa n.versicolor m.setosa m.versicolor sd.setosa sd.versicolor
#> [1,]       15           15     7.61        -8.97      1.79          2.35
#>      pooled.sd  diff    D
#> [1,]      2.09 16.59 7.93
# print table of predicted probabilities
D.iris_out$P.table
#>             
#>              [0,0.25) [0.25,0.5) [0.5,0.75) [0.75,1]
#>   setosa            0          0          0        1
#>   versicolor        1          0          0        0
```

### Multiple samples with two groups in each

This example first generates artificial multi-group data which are then
used as separate data folds in the regularization procedure following
separate predictions made for each fold.

``` r
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

# Different partitions for regularization and estimation for each data fold.

# Request probabilities of correct classification (pcc) and area under the receiver operating characteristics (auc) for the output.

round(D_regularized(
  data = d,
  mv.vars = c("x1", "x2", "x3"),
  group.var = "sex",
  group.values = c("female", "male"),
  fold.var = "fold",
  size = 17,
  out = TRUE,
  fold = TRUE,
  rename.output = TRUE,
  pcc = TRUE,
  auc = TRUE
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
#>   pcc.female pcc.male pcc.total  auc pooled.sd.female pooled.sd.male
#> A       0.56     0.39      0.48 0.49             0.15           0.15
#> B       0.36     0.61      0.48 0.50             0.15           0.15
#> C       0.32     0.46      0.41 0.34             0.15           0.15
#> D       0.53     0.49      0.51 0.49             0.15           0.15
#> E       0.37     0.49      0.44 0.37             0.15           0.15
#> F       0.44     0.36      0.42 0.39             0.15           0.15
#> G       0.38     0.59      0.48 0.53             0.15           0.15
#> H       0.73     0.47      0.59 0.57             0.15           0.15
#> I       0.48     0.64      0.57 0.60             0.15           0.15
#> J       0.55     0.45      0.49 0.52             0.15           0.15
#>   pooled.sd.total d.sd.total
#> A            0.15       0.03
#> B            0.15      -0.04
#> C            0.15      -0.54
#> D            0.15      -0.01
#> E            0.15      -0.48
#> F            0.15      -0.54
#> G            0.15       0.06
#> H            0.15       0.19
#> I            0.15       0.32
#> J            0.15       0.02
```

### Comparison of Mahalanobis’ D and Regularized D when Difference in Population Exists

This example compares a measure of standardized distance between group
centroids (Mahalanobis’ D) and a regularized variant provided in the
multid-package in small-sample scenario when the distance between group
centroids in the population is D = 1.

``` r
set.seed(8327482)
# generate data from sixteen correlated (r = .20) variables each with d = .50 difference (equals to Mahalanobis' D = 1)
k=16
r=0.2
d=0.5
n=200

# population correlation matrix
cor_mat<-matrix(ncol=k,nrow=k,rep(r,k*k))
diag(cor_mat)<-1

# population difference vector
d_vector<-rep(d,k)

# population Mahalanobis' D is exactly 1

sqrt(t(d_vector) %*% solve(cor_mat) %*% d_vector)
#>      [,1]
#> [1,]    1
# generate data
library(MASS)

male.dat<-
  data.frame(sex="male",
             mvrnorm(n = n/2,
                     mu = 0.5*d_vector,
                     Sigma = cor_mat,empirical = F))

female.dat<-
  data.frame(sex="female",
             mvrnorm(n = n/2,
                     mu = -0.5*d_vector,
                     Sigma = cor_mat,empirical = F))

dat<-rbind(male.dat,female.dat)

# sample Mahalanobis' D

# obtain mean differences

d_vector_sample<-rep(NA,k)

for (i in 1:k){
  d_vector_sample[i]<-mean(male.dat[,i+1]-female.dat[,i+1])
  
}

# sample pooled covariance matrix (use mean, because equal sample sizes)

cov_mat_sample<-
  (cov(male.dat[,2:17])+cov(female.dat[,2:17]))/2

# calculate sample Mahalanobis' D
sqrt(t(d_vector_sample) %*% solve(cov_mat_sample) %*% d_vector_sample)
#>          [,1]
#> [1,] 1.265318
# calculate elastic net D

D.ela<-
  D_regularized(data=dat,
              mv.vars=paste0("X",1:k),
              group.var = "sex",
              group.values = c("male","female"))

round(D.ela$D,2)
#>      n.male n.female m.male m.female sd.male sd.female pooled.sd diff    D
#> [1,]    100      100   0.54    -0.52    0.86      0.88      0.87 1.06 1.22
# use separate data for regularization and estimation

D.ela_out<-D_regularized(data=dat,
              mv.vars=paste0("X",1:k),
              group.var = "sex",
              group.values = c("male","female"),
              out=T,size = 50,pcc = T, auc=T,pred.prob = T)

round(D.ela_out$D,2)
#>      n.male n.female m.male m.female sd.male sd.female pooled.sd diff    D
#> [1,]     50       50   0.55    -0.33    1.04      0.92      0.99 0.88 0.89
#>      pcc.male pcc.female pcc.total  auc
#> [1,]     0.68       0.68      0.68 0.73
# Table of predicted probabilites
D.ela_out$P.table
#>         
#>          [0,0.2) [0.2,0.4) [0.4,0.6) [0.6,0.8) [0.8,1]
#>   female    0.10      0.34      0.32      0.20    0.04
#>   male      0.00      0.18      0.26      0.36    0.20
```

### Comparison of Mahalanobis’ D and Regularized D when Sex Difference in Population Does Not Exist

This example compares a measure of standardized distance between group
centroids (Mahalanobis’ D) and a regularized variant provided in the
multid-package in small-sample scenario when the group centroids in the
population is are at the same location, D = 0. In this sample,
Mahalanobis’ D is measured at D = 0.5, elastic net D with same data used
for regularization and estimation at D = 0.35, whereas elastic net D
with independent estimation data shows D = 0.

``` r
set.seed(8327482)
# generate data from sixteen correlated (r = .20) variables each with d = .00 difference (equals to Mahalanobis' D = 0)
k=16
r=0.2
d=0.0
n=200

# population correlation matrix
cor_mat<-matrix(ncol=k,nrow=k,rep(r,k*k))
diag(cor_mat)<-1

# population difference vector
d_vector<-rep(d,k)

# population Mahalanobis' D is exactly 1

sqrt(t(d_vector) %*% solve(cor_mat) %*% d_vector)
#>      [,1]
#> [1,]    0
# generate data

male.dat<-
  data.frame(sex="male",
             mvrnorm(n = n/2,
                     mu = 0.5*d_vector,
                     Sigma = cor_mat,empirical = F))

female.dat<-
  data.frame(sex="female",
             mvrnorm(n = n/2,
                     mu = -0.5*d_vector,
                     Sigma = cor_mat,empirical = F))

dat<-rbind(male.dat,female.dat)

# sample Mahalanobis' D

# obtain mean differences

d_vector_sample<-rep(NA,k)

for (i in 1:k){
  d_vector_sample[i]<-mean(male.dat[,i+1]-female.dat[,i+1])
  
}

# sample pooled covariance matrix (use mean, because equal sample sizes)

cov_mat_sample<-
  (cov(male.dat[,2:17])+cov(female.dat[,2:17]))/2

# calculate sample Mahalanobis' D
sqrt(t(d_vector_sample) %*% solve(cov_mat_sample) %*% d_vector_sample)
#>           [,1]
#> [1,] 0.5316555
# calculate elastic net D

D.ela.zero<-
  D_regularized(data=dat,
              mv.vars=paste0("X",1:k),
              group.var = "sex",
              group.values = c("male","female"))

round(D.ela.zero$D,2)
#>      n.male n.female m.male m.female sd.male sd.female pooled.sd diff    D
#> [1,]    100      100   0.02    -0.02     0.1       0.1       0.1 0.04 0.35
# use separate data for regularization and estimation

D.ela.zero_out<-
  D_regularized(data=dat,
              mv.vars=paste0("X",1:k),
              group.var = "sex",
              group.values = c("male","female"),
              out=T,size = 50,pcc = T, auc=T,pred.prob = T)

round(D.ela.zero_out$D,2)
#>      n.male n.female m.male m.female sd.male sd.female pooled.sd diff   D
#> [1,]     50       50      0        0       0         0         0    0 NaN
#>      pcc.male pcc.female pcc.total auc
#> [1,]        1          0       0.5 0.5
# Table of predicted probabilites
D.ela.zero_out$P.table
#>         
#>          [0,0.2) [0.2,0.4) [0.4,0.6) [0.6,0.8) [0.8,1]
#>   female       0         0         1         0       0
#>   male         0         0         1         0       0
```

### Distribution overlap

This example shows how the degree of overlap between the predicted
values across the two groups can be visualized and estimated.

For parametric variants, see [Del Giudice (in
press)](https://marcodgdotnet.files.wordpress.com/2019/04/delgiudice_measuring_sex-differences-similarities_pre.pdf).

For non-parametric variants, see [Pastore
(2018)](https://doi.org/10.21105/joss.01023) and [Pastore & Calcagnì
(2019)](https://doi.org/10.3389/fpsyg.2019.01089).

``` r
# Use predicted values from elastic net D (out) when difference in population exists-

library(ggplot2)

ggplot(D.ela_out$pred.dat,
       aes(x=pred,fill=group))+
  geom_density(alpha=0.5)+
  xlab("Predicted log odds of being male (FM-score)")
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

``` r
# parametric overlap 

## Proportion of overlap relative to a single distribution (OVL)

## obtain D first
(D<-unname(D.ela_out$D[,"D"]))
#> [1] 0.8884007
(OVL<-2*pnorm((-D/2)))
#> [1] 0.6568977
## Proportion of overlap relative to the joint distribution

(OVL2<-OVL/(2-OVL))
#> [1] 0.4890899
# non-parametric overlap

library(overlapping)

np.overlap<-
  overlap(x = list(D.ela_out$pred.dat[
  D.ela_out$pred.dat$group=="male","pred"],
  D.ela_out$pred.dat[
  D.ela_out$pred.dat$group=="female","pred"]),
  plot=T)
```

<img src="man/figures/README-unnamed-chunk-4-2.png" width="100%" />

``` r
# this corresponds to Proportion of overlap relative to the joint distribution (OVL2)
(np.OVL2<-unname(np.overlap$OV))
#> [1] 0.5412493
# from which Proportion of overlap relative to a single distribution (OVL) is approximated at
(np.OVL<-(2*np.OVL2)/(1+np.OVL2))
#> [1] 0.7023514
# compare overlaps

round(cbind(OVL,np.OVL,OVL2,np.OVL2),2)
#>       OVL np.OVL OVL2 np.OVL2
#> [1,] 0.66    0.7 0.49    0.54
```
