#' Quantile correlation coefficient
#'
#' For computation of tail dependence as correlations estimated at different variable quantiles (Choi & Shin, 2022; Lee et al., 2022) summarized across two quantile regression models where x and y switch roles as independent/dependent variable.
#'
#' @param data Data frame.
#' @param x Name of x variable. Character string.
#' @param y Name of y variable. Character string.
#' @param tau The quantile(s) to be estimated. A vector of values between 0 and 1, default c(.1,.5,.9). @seealso \code{\link[quantreg]{rq}}
#' @param method The algorithmic method used to compute the fit (default "br"). @seealso \code{\link[quantreg]{rq}}
#'
#' @return Returns bivariate correlations at different tau values (quantiles) of the variables and a linear Pearson's correlation estimate for comparison.
#' @references Choi, J.-E., & Shin, D. W. (2022). Quantile correlation coefficient: A new tail dependence measure. Statistical Papers, 63(4), 1075–1104. https://doi.org/10.1007/s00362-021-01268-7
#' @references Lee, J. A., Bardi, A., Gerrans, P., Sneddon, J., van Herk, H., Evers, U., & Schwartz, S. (2022). Are value–behavior relations stronger than previously thought? It depends on value importance. European Journal of Personality, 36(2), 133–148. https://doi.org/10.1177/08902070211002965
#' @export
#'
#' @examples
#' d <- data.frame(
#'   x = rnorm(100),
#'   y = rnorm(100)
#' )
#' qcc(x="x",y="y",data=d)
qcc<-function(x,y,tau=c(.1,.5,.9),data,method="br"){
  rq1<-quantreg::rq(stats::as.formula(paste(y,"~",x)),
                    data=data,tau = tau,method = method)

  b1<-rq1$coefficients[rownames(rq1$coefficients)==x,]

  rq2<-quantreg::rq(as.formula(paste(x,"~",y)),
                    data=data,tau = tau,method = method)

  b2<-rq2$coefficients[rownames(rq2$coefficients)==y,]

  rho_tau<-sign(b1)*sqrt(b1*b2)

  r<-cor(data[,x],data[,y],method="pearson")

  output=
    list(r=r,
         rho_tau=rho_tau)

  return(output)
}
