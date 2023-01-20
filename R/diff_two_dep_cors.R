#' Difference between two dependent Pearson's correlations (with common index)
#'
#' Calculates Cohen's q effect size statistic for difference between two correlations, r_yx1 and r_yx2.
#' Tests if Cohen's q is different from zero while accounting for dependency between the two correlations.
#'
#' @param data Data frame.
#' @param y Character. Variable name of the common index variable.
#' @param x1 Character. Variable name.
#' @param x2 Character. Variable name.
#' @param level Numeric. The confidence level required for the result output (Default .95)
#'
#' @return Parameter estimates from the fitted structural path model.
#' @export
#'
#' @examples
#' set.seed(3864)
#' d<-data.frame(y=rnorm(100),x=rnorm(100))
#' d$x1<-d$x+rnorm(100)
#' d$x2<-d$x+rnorm(100)
#' diff_two_dep_cors(data=d,y="y",x1="x1",x2="x2")
diff_two_dep_cors <-
  function(data, y, x1, x2, level = .95) {
    mod <- paste0(
      paste0(y, "~~var_y*", y), "\n",
      paste0(x1, "~~var_x1*", x1), "\n",
      paste0(x2, "~~var_x2*", x2), "\n",
      paste0(y, "~~cov_yx1*", x1), "\n",
      paste0(y, "~~cov_yx2*", x2), "\n",
      paste0(x1, "~~cov_x1x2*", x2), "\n",
      paste0("r_yx1:=cov_yx1/(sqrt(var_y)*sqrt(var_x1))"), "\n",
      paste0("r_yx2:=cov_yx2/(sqrt(var_y)*sqrt(var_x2))"), "\n",
      paste0("r_x1x2:=cov_x1x2/(sqrt(var_x1)*sqrt(var_x2))"), "\n",
      paste0("diff_r_yx1_vs_r_yx2:=r_yx1-r_yx2"), "\n",
      paste0("z_yx1:=(1/2)*log((1+r_yx1)/(1-r_yx1))"), "\n",
      paste0("z_yx2:=(1/2)*log((1+r_yx2)/(1-r_yx2))"), "\n",
      paste0("Cohens_q:=z_yx1-z_yx2")
    )

    output <-
      lavaan::parameterestimates(lavaan::sem(mod, data = data, estimator = "ML"),
        ci = TRUE, level = level, output = "data.frame"
      )
    rownames(output) <- output$label
    output <- output[, c(5:ncol(output))]

    return(output)
  }
