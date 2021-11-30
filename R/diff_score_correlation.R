#' Test assumptions and infer about associations with algebraic difference scores
#'
#' @param data A data frame.
#' @param var1 Character string. Variable name of first component of difference score (Y_1).
#' @param var2 Character string. Variable name of second component of difference score (Y_2).
#' @param predictor Character string. Variable name of independent variable predicting difference.
#' @param estimator Character string. Estimator used in SEM (Default "MLR").
#' @param scale Logical. Are var1 and var2 scaled with their pooled sd? (Default FALSE)
#' @param center Logical. Are var1 and var2 centered around their grand mean? (Default FALSE)
#' @param sampling.weights Character string. Name of sampling weights variable.
#'
#' @return
#' \item{descriptives}{Means, standard deviations, and intercorrelations.}
#' \item{results}{Parameter estimates from the structural equation model.}
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
#' diff_score_correlation(
#'   data = d, var1 = "var1", var2 = "var2",
#'   predictor = "x", center = TRUE, scale = TRUE
#' )
diff_score_correlation <- function(data,
                                   var1,
                                   var2,
                                   center = FALSE,
                                   scale = FALSE,
                                   predictor,
                                   estimator = "MLR",
                                   sampling.weights = NULL) {
  vt <- stats::var.test(
    data[, var1], data[, var2]
  )

  F.test <-
    paste0(
      "F test: F(",
      unname(vt$parameter[1]), ", ",
      unname(vt$parameter[2]), ") = ",
      unname(as.character(round(vt$statistic, 5))),
      ", p = ",
      unname(as.character(round(vt$p.value, 8)))
    )

  fl <- stats::fligner.test(list(
    data[, var1],
    data[, var2]
  ))

  fl.test <-
    paste0(
      "Fligner-Killeen test: chi-squared(df = ",
      unname(fl$parameter), ") = ",
      unname(as.character(round(fl$statistic, 5))),
      ", p = ",
      unname(as.character(round(fl$p.value, 8)))
    )

  variances <- unname(rbind(F.test, fl.test))


  if (center) {
    pooled_mean <-
      mean(c(mean(data[, var1]), mean(data[, var2])))

    data[, var1] <- data[, var1] - pooled_mean
    data[, var2] <- data[, var2] - pooled_mean
  }

  if (scale) {
    pooled_sd <-
      sqrt(((nrow(data) - 1) * stats::sd(data[, var1])^2 +
        (nrow(data) - 1) * stats::sd(data[, var2])^2) /
        (nrow(data) * 2 - 2))

    data[, var1] <- data[, var1] / pooled_sd
    data[, var2] <- data[, var2] / pooled_sd
    data[, predictor] <-
      (data[, predictor] - mean(data[, predictor])) /
        stats::sd(data[, predictor])
  }


  data[, "diff"] <- data[, var1] - data[, var2]

  descriptives <-
    cbind(
      rbind(
        c(mean(data[, var1]), stats::sd(data[, var1])),
        c(mean(data[, var2]), stats::sd(data[, var2])),
        c(mean(data[, "diff"]), stats::sd(data[, "diff"])),
        c(mean(data[, predictor]), stats::sd(data[, predictor]))
      ),
      stats::cor(data[, c(var1, var2, "diff", predictor)])
    )
  colnames(descriptives) <-
    c("M", "SD", var1, var2, "diff", predictor)

  model <- paste0(
    paste0(var1, "~b_11*", predictor), "\n",
    paste0(var2, "~b_21*", predictor), "\n",
    paste0(var1, "~b_10*1"), "\n",
    paste0(var2, "~b_20*1"), "\n",
    paste0(var1, "~~rescov_12*", var2), "\n",
    paste0(var1, "~~e_1*", var1), "\n",
    paste0(var2, "~~e_2*", var2), "\n",
    paste0(predictor, "~pred_int*1"), "\n",
    paste0(predictor, "~~pred_var*", predictor), "\n",
    paste0("coef_diff:=b_11-b_21"), "\n",
    paste0(
      "coef_diff_std:=((b_11-b_21)*sqrt(pred_var))/",
      as.character(round(descriptives["diff", 2], 5))
    ), "\n",
    paste0("coef_sum:=b_21+b_11")
  )

  fit <-
    lavaan::sem(
      model = model,
      data = data,
      estimator = estimator,
      sampling.weights = sampling.weights
    )

  pretty.out <-
    lavaan::parameterestimates(fit, output = "pretty")


  output <- list(
    variances = variances,
    descriptives = descriptives,
    results = pretty.out
  )

  return(output)
}
