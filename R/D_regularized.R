#' Multivariate group difference estimation with regularized binomial regression
#'
#' @param data A data frame.
#' @param mv.vars Character vector. Variable names in the multivariate variable set.
#' @param group.var The name of the group variable.
#' @param group.values Vector of length 2, group values (e.g. c("male", "female) or c(0,1)).
#' @param alpha Alpha-value for penalizing function ranging from 0 to 1: 0 = ridge regression, 1 = lasso, 0.5 = elastic net (default).
#' @param nfolds Number of folds used for obtaining lambda (range from 3 to n-1, default 10).
#' @param s Which lambda value is used for predicted values? Either "lambda.min" (default) or "lambda.1se".
#' @param type.measure Which measure is used during cross-validation. Default "deviance".
#' @param rename.output Logical. Should the output values be renamed according to the group.values? Default TRUE.
#' @param out Logical. Should results and predictions be calculated on out-of-bad data set? (Default FALSE)
#' @param size Integer. Size of regularization data per each group. Default 1/4 of cases.
#'
#' @return
#' \item{D}{Multivariate descriptives and differences}
#' \item{pred.dat}{A data.frame with predicted values}
#' @export
#'
#' @examples D_regularized(
#'   data = iris[iris$Species == "setosa" | iris$Species == "versicolor", ],
#'   mv.vars = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"),
#'   group.var = "Species", group.values = c("setosa", "versicolor")
#' )$D
D_regularized <-
  function(data,
           mv.vars,
           group.var,
           group.values,
           alpha = 0.5,
           nfolds = 10,
           s = "lambda.min",
           type.measure = "deviance",
           rename.output = TRUE,
           out=FALSE,
           size=NULL) {
    data$group.var.num <-
      ifelse(data[, group.var] == group.values[1], 1,
        ifelse(data[, group.var] == group.values[2], 0,
          NA
        )
      )
    if (out){
      multid::D_regularized_out(data=data,
                                mv.vars=mv.vars,
                                group.var=group.var,
                                group.values=group.values,
                                alpha=alpha,
                                size=size,
                                rename.output=rename.output,
                                type.measure=type.measure,
                                s=s,
                                nfolds=nfolds)
    }

    else {

    cv.mod <-
      glmnet::cv.glmnet(
        x = as.matrix(data[, c(mv.vars)]),
        y = data$group.var.num,
        family = c("binomial"),
        nfolds = nfolds,
        type.measure = type.measure,
        alpha = alpha
      )

    preds <- data.frame(
      group = data[, group.var],
      pred = as.numeric(
        stats::predict(cv.mod,
          newx = as.matrix(data[, c(mv.vars)]),
          s = s
        )
      )
    )

    D <- multid::d_pooled_sd(
      data = preds,
      var = "pred",
      group.var = "group",
      group.values = group.values,
      rename.output = rename.output
    )

    comb.output <- list(D = D, pred.dat = preds)
    return(comb.output) }
  }
