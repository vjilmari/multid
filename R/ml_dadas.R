#' Predicting algebraic difference scores in multilevel model
#'
#' Decomposes difference score predictions to predictions of difference score components by probing simple effects at the levels of the binary moderator.
#'
#' @param model Multilevel model fitted with lmerTest.
#' @param predictor Character string. Variable name of independent variable predicting difference score.
#' @param diff_var Character string. A variable indicative of difference score components (two groups).
#' @param diff_var_values Vector. Values of the component score groups in diff_var.
#' @param scaled_estimates Logical. Are scaled estimates obtained? Does fit a reduced model for correct standard deviations. (Default FALSE)
#' @param re_cov_test Logical. Significance test for random effect covariation? Does fit a reduced model without the correlation. (Default FALSE)
#'
#' @return
#' \item{dadas}{A data frame including regression coefficients for component scores and dadas.}
#' \item{scaled_estimates}{Scaled regression coefficients for difference score components and difference score.}
#' \item{vpc_at_reduced}{Variance partition coefficients in the model without the predictor and interactions.}
#' \item{re_cov_test}{Likelihood ratio significance test for random effect covariation.}
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(95332)
#' n1 <- 10 # groups
#' n2 <- 10 # observations per group
#'
#' dat <- data.frame(
#'   group = rep(c(LETTERS[1:n1]), each = n2),
#'   x = sample(c(-0.5, 0.5), n1 * n2, replace = TRUE),
#'   w = rep(sample(1:5, n1, replace = TRUE), each = n2),
#'   y = sample(1:5, n1 * n2, replace = TRUE)
#' )
#' library(lmerTest)
#' fit <- lmerTest::lmer(y ~ x * w + (x | group),
#'   data = dat
#' )
#'
#' round(ml_dadas(fit,
#'   predictor = "w",
#'   diff_var = "x",
#'   diff_var_values = c(0.5, -0.5)
#' )$dadas, 3)
#' }
ml_dadas <- function(model,
                     predictor,
                     diff_var,
                     diff_var_values,
                     scaled_estimates = FALSE,
                     re_cov_test = FALSE) {

  # reorder diff_var_values
  diff_var_values <-
    diff_var_values[order(diff_var_values)]

  # format list for contrast values
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

  trend.diff <- df.trends[1, 2] - df.trends[2, 2]
  trend.sum <- df.trends[1, 2] + df.trends[2, 2]

  trend.strength <- df.trends[1, 2] - df.trends[2, 2]

  # define correct contrasts based on signs

  if (trend.diff > 0 & trend.sum < 0) {
    mlist <- list(
      abs_diff = c(1, -1),
      abs_sum = c(-1, -1)
    )
  } else if (trend.diff > 0 & trend.sum > 0) {
    mlist <- list(
      abs_diff = c(1, -1),
      abs_sum = c(1, 1)
    )
  } else if (trend.diff < 0 & trend.sum < 0) {
    mlist <- list(
      abs_diff = c(-1, 1),
      abs_sum = c(-1, -1)
    )
  } else if (trend.diff < 0 & trend.sum > 0) {
    mlist <- list(
      abs_diff = c(-1, 1),
      abs_sum = c(1, 1)
    )
  }

  temp.cont <-
    emmeans::contrast(trends,
      method = mlist,
      side = ">"
    )

  ml_abstest <-
    data.frame(emmeans::contrast(temp.cont,
      method = list(abs_test = c(1, -1)),
      side = ">"
    ))

  trends.df <- data.frame(trends)
  colnames(trends.df) <- colnames(data.frame(temp.cont))

  dadas <- rbind(
    trends.df,
    data.frame(temp.cont),
    ml_abstest
  )
  rownames(dadas) <- dadas$contrast
  dadas <- dadas[, 2:ncol(dadas)]

  rownames(dadas) <- c(
    rownames(dadas)[1:(nrow(dadas) - 1)],
    "dadas"
  )

  output <- list(dadas = dadas)

  if (scaled_estimates | re_cov_test) {
    # refit a reduced model to obtain scaled estimates
    # obtain fixed effects as character vector

    FEs <- attributes(stats::terms(model))$term.labels

    # obtain random effects and DV as a character vector

    DV <- as.character(stats::formula(model,
      random.only = T
    ))[2]
    RE <- as.character(stats::formula(model,
      random.only = T
    ))[3]

    # define the interaction term correctly

    if (paste0(predictor, ":", diff_var) %in% FEs) {
      int.term <- paste0(predictor, ":", diff_var)
    } else {
      int.term <- paste0(diff_var, ":", predictor)
    }

    # reformulate

    new.formula <-
      stats::reformulate(c(FEs[-c(
        which(FEs == predictor),
        which(FEs == int.term)
      )], RE), response = DV)

    # update to reduced model

    reduced_model <-
      stats::update(model, new.formula)

    # get component standard deviations

    vpc_at_reduced <-
      vpc_at(
        model = reduced_model,
        lvl1.var = diff_var,
        lvl1.values = diff_var_values
      )

    # get slope sd

    vc <-
      as.data.frame(lme4::VarCorr(reduced_model))

    slope_sd_reduced <-
      vc[vc$var1 == diff_var &
        !is.na(vc$var1) &
        is.na(vc$var2), "sdcor"]

    # compile to same frame

    scaled_estimates_df <-
      rbind(
        c(
          est = dadas[1, "estimate"],
          scaling_SD =
            vpc_at_reduced[
              vpc_at_reduced$lvl1.value ==
                rownames(dadas[1, ]),
              "Intercept.sd"
            ],
          scaled_est = dadas[1, "estimate"] /
            vpc_at_reduced[
              vpc_at_reduced$lvl1.value ==
                rownames(dadas[1, ]),
              "Intercept.sd"
            ]
        ),
        c(
          est = dadas[2, "estimate"],
          scaling_SD =
            vpc_at_reduced[
              vpc_at_reduced$lvl1.value ==
                rownames(dadas[2, ]),
              "Intercept.sd"
            ],
          scaled_est = dadas[2, "estimate"] /
            vpc_at_reduced[
              vpc_at_reduced$lvl1.value ==
                rownames(dadas[2, ]),
              "Intercept.sd"
            ]
        ),
        c(
          est = dadas[1, "estimate"] - dadas[2, "estimate"],
          scaling_SD = slope_sd_reduced,
          scaled_est = (dadas[1, "estimate"] - dadas[2, "estimate"]) /
            slope_sd_reduced
        )
      )

    rownames(scaled_estimates_df) <-
      c(diff_var_values, "difference")

    output <- list(
      dadas = dadas,
      scaled_estimates = scaled_estimates_df,
      vpc_at_reduced = vpc_at_reduced
    )
  }

  if (re_cov_test) {

    # drop the random effect correlation
    RE.no.cov <-
      gsub(
        x = RE, pattern = " | ",
        replacement = " || ",
        fixed = TRUE
      )

    # reformulate
    no.cov.formula <-
      stats::reformulate(c(FEs[-c(
        which(FEs == predictor),
        which(FEs == int.term)
      )], RE.no.cov), response = DV)

    # update to model without the covariance

    no.cov.mod <-
      stats::update(model, no.cov.formula)

    # test against the less reduced model

    re.cov.test <-
      stats::anova(
        no.cov.mod,
        reduced_model
      )

    # obtain important numbers

    re_cov_reduced <-
      vc[vc$var1 == "(Intercept)" &
        vc$var2 == diff_var &
        !is.na(vc$var1) &
        !is.na(vc$var2), c("vcov", "sdcor")]

    # compile to a frame

    re_cov_test_df <-
      c(
        RE_cov = re_cov_reduced[1, 1],
        RE_cor = re_cov_reduced[1, 2],
        Chisq = re.cov.test$Chisq[2],
        Df = re.cov.test$Df[2],
        p = re.cov.test$`Pr(>Chisq)`[2]
      )

    output <- list(
      dadas = dadas,
      scaled_estimates = scaled_estimates_df,
      vpc_at_reduced = vpc_at_reduced,
      re_cov_test = re_cov_test_df
    )
  }

  return(output)
}
