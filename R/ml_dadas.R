#' Predicting algebraic difference scores in multilevel model
#'
#' @param model Multilevel model fitted with lmerTest.
#' @param predictor Character string. Variable name of independent variable predicting difference score.
#' @param diff_var Character string. A variable indicative of difference score components (two groups).
#' @param diff_var_values Vector. Values of the component score groups in diff_var.
#'
#' @return A data frame including regression coefficients for component scores and dadas (One sided dadas-test for positivity of abs(b_11-b_21)-abs(b_11+b_21)).
#' @export
#'
#' @examples
ml_dadas <- function(model,
                     predictor,
                     diff_var,
                     diff_var_values) {
  at.list <- list(diff_var_values)
  names(at.list) <- diff_var

  trends <- emmeans::emtrends(
    object = model,
    specs = diff_var,
    var = predictor,
    lmerTest.limit = nrow(model@frame),
    disable.pbkrtest = TRUE,
    infer = c(FALSE, TRUE),
    at = at.list
  )
  # obtain trend signs for correct contrasts

  df.trends <- data.frame(trends)
  trend.signs <- c(
    df.trends[1, 2] / abs(df.trends[1, 2]),
    df.trends[2, 2] / abs(df.trends[2, 2])
  )

  # define correct contrasts based on signs


  if (trend.signs[1] == -1 & trend.signs[2] == -1) {
    mlist <- list(
      abs_diff = c(-1, 1),
      abs_sum = c(-1, -1)
    )
  } else if (trend.signs[1] == 1 & trend.signs[2] == 1) {
    mlist <- list(
      abs_diff = c(-1, 1),
      abs_sum = c(1, 1)
    )
  } else if (trend.signs[1] == -1 & trend.signs[2] == 1) {
    mlist <- list(
      abs_diff = c(-1, 1),
      abs_sum = c(-1, -1)
    )
  } else if (trend.signs[1] == 1 & trend.signs[2] == -1) {
    mlist <- list(
      abs_diff = c(1, -1),
      abs_sum = c(-1, -1)
    )
  }

  temp.cont <-
    emmeans::contrast(trends,
      method = mlist
    )

  ml_abstest <-
    data.frame(emmeans::contrast(temp.cont,
      method = list(abs_test = c(1, -1)),
      side = ">"
    ))

  trends.df <- data.frame(trends)
  colnames(trends.df) <- colnames(data.frame(temp.cont))


  output <- rbind(
    trends.df,
    data.frame(temp.cont),
    ml_abstest
  )
  rownames(output) <- output$contrast
  output <- output[, 2:ncol(output)]

  return(output)
}
