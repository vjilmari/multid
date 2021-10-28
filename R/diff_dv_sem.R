#' Algebraic difference scores as dependent variables in regression analysis conducted via structural equation modeling
#'
#' @param data A data frame.
#' @param var1 Character string. Variable name of first component of difference score (Y_1).
#' @param var2 Character string. Variable name of second component of difference score (Y_2).
#' @param predictor Character string. Variable name of independent variable predicting difference.
#' @param estimator Character string. Estimator used in SEM (Default "MLR").
#' @param sampling.weights Character string. Name of sampling weights variable.
#'
#' @return
#' \item{aic}{Fit statistics and model probabilities for all fitted models.}
#' \item{models}{List of fitted models.}
#' @references Edwards, J. R. (1995). Alternatives to Difference Scores as Dependent Variables in the Study of Congruence in Organizational Research. Organizational Behavior and Human Decision Processes, 64(3), 307–324. <doi:10.1006/obhd.1995.1108>.
#'
#' @export
#'
#' @examples set.seed(342356)
#' d <- data.frame(
#'   var1 = rnorm(50),
#'   var2 = rnorm(50),
#'   x = rnorm(50)
#' )
#'
#' dds <- diff_dv_sem(data = d, var1 = "var1", var2 = "var2", predictor = "x")
#' round(dds$aic, 2)
#' lavaan::summary(dds$models$opposite.ev,
#'   standardized = TRUE
#' )
diff_dv_sem <- function(data,
                        var1,
                        var2,
                        predictor,
                        estimator = "MLR",
                        sampling.weights = NULL) {
  model <- paste0(
    paste0(var1, "~b_11*", predictor), "\n",
    paste0(var2, "~b_21*", predictor), "\n",
    paste0(var1, "~b_10*1"), "\n",
    paste0(var2, "~b_20*1"), "\n",
    paste0(var1, "~~", var2), "\n",
    paste0(var1, "~~e_1*", var1), "\n",
    paste0(var2, "~~e_2*", var2), "\n",
    paste0("coef_diff:=b_21-b_11"), "\n",
    paste0("coef_abs_diff:=abs(b_21)-abs(b_11)"), "\n",
    paste0("coef_sum:=b_21+b_11")
  )

  fit <-
    lavaan::sem(
      model = model,
      data = data,
      estimator = estimator,
      sampling.weights = sampling.weights
    )

  addList <- list(
    all.free = "",
    both.null = "b_21==0\nb_11==0",
    b.null = "b_21==0",
    a.null = "b_11==0",
    equal = "b_21==b_11",
    opposite = "b_21==-1*b_11",
    all.free.ev = "e_1==e_2",
    both.null.ev = "b_21==0\nb_11==0\ne_1==e_2",
    b.null.ev = "b_21==0\ne_1==e_2",
    a.null.ev = "b_11==0\ne_1==e_2",
    equal.ev = "b_21==b_11\ne_1==e_2",
    opposite.ev = "b_21==-1*b_11\ne_1==e_2"
  )


  newmodels <- c(model)

  fit <- stats::update(fit, model = newmodels)

  newFit <- lapply(addList, function(x) {
    print(x)
    stats::update(fit, model = c(newmodels, x))
  })

  get_sem_fit_aic <-
    function(model) {
      if (estimator == "MLR") {
        fits <-
          lavaan::inspect(model, "fit")[c(
            "npar",
            "chisq.scaled", "df.scaled", "pvalue.scaled", "cfi.robust",
            "tli.robust", "rmsea.robust", "srmr", "logl", "aic"
          )]
      } else {
        fits <-
          lavaan::inspect(model, "fit")[c(
            "npar",
            "chisq", "df", "pvalue", "cfi",
            "tli", "rmsea", "srmr", "logl", "aic"
          )]
      }


      n <- lavaan::inspect(model, "nobs")

      c(n = n, fits)
    }

  aictab <-
    data.frame(t(data.frame(lapply(newFit, get_sem_fit_aic))))

  aictab$AICc <- aictab$aic + 2 * aictab$npar +
    (2 * aictab$npar * (aictab$npar + 1)) / (aictab$n - aictab$npar - 1)

  aictab$dAICc <-
    aictab$AICc - min(aictab$AICc)

  aictab$rel_likelihood <-
    exp((-0.5) * aictab$dAICc)

  aictab$norm_likelihood <-
    aictab$rel_likelihood / sum(aictab$rel_likelihood)

  aicresults <-
    aictab[order(-aictab$norm_likelihood), ]

  output <- list(aic = aicresults, models = newFit)

  return(output)
}
