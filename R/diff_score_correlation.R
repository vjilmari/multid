#' Test assumptions and infer about associations with algebraic difference scores
#'
#' @param data A data frame.
#' @param var1 Character string. Variable name of first component of difference score (Y_1).
#' @param var2 Character string. Variable name of second component of difference score (Y_2).
#' @param predictor Character string. Variable name of independent variable predicting difference.
#' @param estimator Character string. Estimator used in SEM (Default "MLR").
#' @param scale Logical. Are var1 and var2 scaled with their pooled sd? (Default FALSE)
#' @param center Logical. Are var1 and var2 centered around their grand mean? (Default FALSE)
#' @param level Numeric. The confidence level required for the result output (Default .95)
#' @param bound_l Equivalence test lower bound based on SESOI (smallest effect size of interest; Default = -0.05).
#' @param bound_u Equivalence test upper bound based on SESOI (smallest effect size of interest; Default = 0.05).
#' @param sampling.weights Character string. Name of sampling weights variable.
#'
#' @return
#' \item{descriptives}{Means, standard deviations, and intercorrelations.}
#' \item{parameter_estimates}{Parameter estimates from the structural equation model.}
#' \item{variances}{F test and Fligner-Killeen test for homogeneity of variances between var1 and var2.}
#' \item{transformed_data}{Data frame with variables used in SEM.}
#' \item{coef_sum_equivalence}{Equivalence test for sum of regression coefficiens for var1 and var2.}
#' \item{abs_test_one_sided}{One sided abs_test for positivity of abs(b_11-b_21)-abs(b_11+b_21).}
#' \item{results}{Summary of key results.}
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
                                   level = .95,
                                   bound_l = -0.05,
                                   bound_u = 0.05,
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
    paste0("coef_sum:=b_11+b_21"), "\n",
    # paste0("sum_eq_lower:=(b_11+b_21)-", as.character(bound_l)), "\n",
    # paste0("sum_eq_upper:=(b_11+b_21)-", as.character(bound_u)), "\n",
    paste0("diff_abs_magnitude:=sqrt(b_11^2)-sqrt(b_21^2)"), "\n",
    paste0("abs_test_two_sided:=sqrt((b_11-b_21)^2)-sqrt((b_11+b_21)^2)")
  )

  fit <-
    lavaan::sem(
      model = model,
      data = data,
      estimator = estimator,
      sampling.weights = sampling.weights
    )

  pretty.out <-
    lavaan::parameterestimates(fit,
      output = "pretty",
      rsquare = TRUE, level = level
    )

  output.data <-
    data[, c(var1, var2, "diff", predictor)]

  pars <- data.frame(lavaan::parameterestimates(fit,
                                                rsquare = TRUE,
                                                level = level))


  coef_sum_equivalence <-
    c(
      Estimate = pars[pars$label == "coef_sum", "est"],
      "Std. Error" = pars[pars$label == "coef_sum", "se"]
    )

  coef_sum_equivalence["z_lower"] <-
    (coef_sum_equivalence["Estimate"] - (bound_l)) / coef_sum_equivalence["Std. Error"]
  coef_sum_equivalence["z_upper"] <-
    (coef_sum_equivalence["Estimate"] - (bound_u)) / coef_sum_equivalence["Std. Error"]
  coef_sum_equivalence["p_lower"] <-
    stats::pnorm(coef_sum_equivalence["z_lower"], lower.tail = F)
  coef_sum_equivalence["p_upper"] <-
    stats::pnorm(coef_sum_equivalence["z_upper"], lower.tail = T)


  abs_test_one_sided <- c(
    Estimate = pars[pars$label == "abs_test_two_sided", "est"],
    "Std. Error" = pars[pars$label == "abs_test_two_sided", "se"],
    z = pars[pars$label == "abs_test_two_sided", "z"]
  )

  abs_test_one_sided["p.pos"] <-
    stats::pnorm(abs_test_one_sided["z"], lower.tail = F)

  rsquared=pars[pars[,"op"]=="r2",c(1,5)]
  rownames(rsquared)<-NULL
  colnames(rsquared)<-c("component","r2")

  res.pars<-c("b_11","b_21","b_10","b_20","rescov_12",
              "coef_diff","coef_diff_std",
              "coef_sum","diff_abs_magnitude",
              "abs_test_two_sided")

  results<-pars[pars$label %in% res.pars,4:ncol(pars)]
  rownames(results)<-results$label
  results<-results[,2:ncol(results)]


  results["abs_test_two_sided","pvalue"]<-
    unname(abs_test_one_sided["p.pos"])
  rownames(results)<-c(rownames(results)[1:(nrow(results)-1)],
                       "abs_test_one_sided")

  output <- list(
    variances = variances,
    descriptives = descriptives,
    parameter_estimates = pars,
    transformed_data = output.data,
    coef_sum_equivalence = coef_sum_equivalence,
    abs_test_one_sided = abs_test_one_sided,
    rsquared=rsquared,
    results=results
  )

  return(output)
}
