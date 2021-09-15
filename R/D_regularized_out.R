#' Use out-of-bag predictions with out and size arguments.
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
#' @param size
#'
#' @return
#' @export
#'
#' @examples
D_regularized_out <-
  function(data,
           mv.vars,
           group.var,
           group.values,
           alpha = 0.5,
           nfolds = 10,
           s = "lambda.min",
           type.measure = "deviance",
           rename.output = TRUE,
           size=NULL) {
    data$group.var.num <-
      ifelse(data[, group.var] == group.values[1], 1,
             ifelse(data[, group.var] == group.values[2], 0,
                    NA
             )
      )

    data$row.nmbr<-rownames(data)

    data.grouped<-dplyr::group_by(data,group.var.num)

    train.data<-dplyr::sample_n(data.grouped,
                         size = size,
                         replace = F)

    test.data<-data[!(data$row.nmbr %in% train.data$row.nmbr),]
    train.data<-dplyr::ungroup(train.data)

    cv.mod <-
      glmnet::cv.glmnet(
        x = as.matrix(train.data[, c(mv.vars)]),
        y = train.data$group.var.num,
        family = c("binomial"),
        nfolds = nfolds,
        type.measure = type.measure,
        alpha = alpha
      )

    preds <- data.frame(
      group = test.data[, group.var],
      pred = as.numeric(
        stats::predict(cv.mod,
                       newx = as.matrix(test.data[, c(mv.vars)]),
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
    return(comb.output)
  }
