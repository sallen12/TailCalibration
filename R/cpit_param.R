#' Conditional PIT values for parametric distributions
#'
#' Calculate conditional probability integral transform (PIT) values for probabilistic
#' forecasts in the form of familiar parametric distributions.
#'
#' @param y vector of observations.
#' @param mean,sd,location,scale,rate,shape,shape1,shape2,df1,df2,df,ncp,meanlog,sdlog,nmeans,nranges,min,max vector of parameters corresponding
#'  to the parametric predictive distributions.
#' @param a,b numeric; lower and upper threshold defining the interval of outcomes that are of interest.
#'
#' @details
#'
#' Probability integral transform (PIT) histograms are a well-established diagnostic
#' tool with which to assess the calibration of probabilistic forecasts. Conditional
#' PIT histograms (cPIT) allow forecast calibration to be assessed whilst focusing on a particular range
#' of outcomes [\code{a}, \code{b}]. This is achieved by restricting attention to the observed values \code{y} that are between
#' \code{a} and \code{b}, and calculating the PIT values corresponding to the conditional
#' distribution given that the observation is in this range (see references for details).
#'
#' The arguments \code{a} and \code{b} are numerical values corresponding to the
#' lower and upper bound of the interval of interest. If we wish to assess forecast calibration
#' when predicting outcomes that exceed a threshold \code{t}, then we can set \code{a = t}
#' and \code{b = Inf}.
#'
#' A NA is returned for the entries of \code{y} that are not in the specified range,
#' otherwise these functions output the cPIT value associated with \code{y} and the
#' corresponding forecast distribution.
#'
#' The cPIT values can be obtained for all continuous parametric distributions
#' that are available in base \proglang{R}: beta, Cauchy, chi-square, exponential,
#' F, gamma, logistic, log-normal, normal, Student's t, Tukey, uniform, and Weibull.
#' The parameters used in these functions are the same as in the standard R
#' distribution handles, e.g. \code{pnorm}.
#'
#'
#' @return vector of conditional PIT values.
#'
#' @author Sam Allen
#'
#' @references
#' \emph{PIT histograms}
#'
#' Dawid, A. P. (1984):
#' `Present position and potential developments: Some personal views: Statistical theory: The prequential approach'.
#' \emph{Journal of the Royal Statistical Society: Series A} 147, 278-290.
#' \doi{10.2307/2981683}
#'
#' \emph{Conditional PIT histograms}
#'
#' Allen, S., Bhend, J., Martius, O. and Ziegel, J (2023):
#' `Weighted verification tools to evaluate univariate and multivariate probabilistic forecasts for high-impact weather events'.
#' \emph{Weather and Forecasting} 38, 499–516.
#' \doi{10.1175/WAF-D-22-0161.1}
#'
#'
#' @examples
#' mu <- rnorm(20, mean = 0, sd = 5)
#' y <- rnorm(20, mean = mu, sd = 1)
#'
#' cpit_norm(y = y, mean = mu, sd = 0.5)
#' pnorm(q = y, mean = mu, sd = 0.5) # pnorm() returns the same values if a = -Inf and b = Inf
#' cpit_norm(y = y, mean = mu, sd = 0.5, a = 0) # restrict attention to values that exceed 0
#' cpit_norm(y = y, mean = mu, sd = 0.5, a = -1, b = 1) # restrict attention to between -1 and 1
#'
#' mu <- rnorm(10000, mean = 0, sd = 5)
#' y <- rnorm(10000, mean = mu, sd = 1)
#'
#' pit_hist(cpit_norm(y = y, mean = mu, sd = 1), ranks = FALSE) # PIT hist
#'
#' pit_hist(cpit_norm(y = y, mean = mu, sd = 1, a = 0), ranks = FALSE) # cPIT hist
#' pit_hist(cpit_norm(y = y, mean = mu + 1, sd = 1, a = 0), ranks = FALSE) # positive bias
#' pit_hist(cpit_norm(y = y, mean = mu - 1, sd = 1, a = 0), ranks = FALSE) # negative bias
#' pit_hist(cpit_t(y = y, df = 1, a = 0), ranks = FALSE)
#'
#' @import stats
#' @name cpit_param
NULL

#' @rdname cpit_param
#' @export
cpit_beta <- function(y, shape1 = 1, shape2 = 1, a = -Inf, b = Inf, return_na = TRUE){
  cpit_dist(y, pbeta, a, b, return_na, shape1 = shape1, shape2 = shape2)
}

#' @rdname cpit_param
#' @export
cpit_cauchy <- function(y, location = 0, scale = 1, a = -Inf, b = Inf, return_na = TRUE){
  cpit_dist(y, pcauchy, a, b, return_na, location = location, scale = scale)
}

#' @rdname cpit_param
#' @export
cpit_chisq <- function(y, df, ncp = 0, a = -Inf, b = Inf, return_na = TRUE){
  cpit_dist(y, pchisq, a, b, return_na, df = df, ncp = ncp)
}

#' @rdname cpit_param
#' @export
cpit_exp <- function(y, rate = 1, a = -Inf, b = Inf, return_na = TRUE){
  cpit_dist(y, pexp, a, b, return_na, rate = rate)
}

#' @rdname cpit_param
#' @export
cpit_f <- function(y, df1, df2, ncp = 0, a = -Inf, b = Inf, return_na = TRUE){
  cpit_dist(y, pf, a, b, return_na, df1 = df1, df2 = df2, ncp = ncp)
}

#' @rdname cpit_param
#' @export
cpit_gamma <- function(y, shape, rate = 1, scale = 1/rate, a = -Inf, b = Inf, return_na = TRUE){
  cpit_dist(y, pgamma, a, b, return_na, shape = shape, rate = rate, scale = scale)
}

#' @rdname cpit_param
#' @export
cpit_logis <- function(y, location = 0, scale = 1, a = -Inf, b = Inf, return_na = TRUE){
  cpit_dist(y, plogis, a, b, return_na, location = location, scale = scale)
}

#' @rdname cpit_param
#' @export
cpit_lnorm <- function(y, meanlog = 0, sdlog = 1, ncp = 0, a = -Inf, b = Inf, return_na = TRUE){
  cpit_dist(y, plnorm, a, b, return_na, meanlog = meanlog, sdlog = sdlog)
}

#' @rdname cpit_param
#' @export
cpit_norm <- function(y, mean = 0, sd = 1, a = -Inf, b = Inf, return_na = TRUE){
  cpit_dist(y, pnorm, a, b, return_na, mean = mean, sd = sd)
}

#' @rdname cpit_param
#' @export
cpit_t <- function(y, df, ncp = 0, a = -Inf, b = Inf, return_na = TRUE){
  cpit_dist(y, pt, a, b, return_na, df = df, ncp = ncp)
}

#' @rdname cpit_param
#' @export
cpit_tukey <- function(y, nmeans, df, nranges = 1, a = -Inf, b = Inf, return_na = TRUE){
  cpit_dist(y, ptukey, a, b, return_na, nmeans = nmeans, df = df, nranges = nranges)
}

#' @rdname cpit_param
#' @export
cpit_unif <- function(y, min, max, a = -Inf, b = Inf, return_na = TRUE){
  cpit_dist(y, punif, a, b, return_na, min = min, max = max)
}

#' @rdname cpit_param
#' @export
cpit_weibull <- function(y, shape, scale = 1, a = -Inf, b = Inf, return_na = TRUE){
  cpit_dist(y, pweibull, a, b, return_na, shape = shape, scale = scale)
}

#' @rdname cpit_param
#' @export
cpit_dist <- function(y, F_x, a = -Inf, b = Inf, return_na = TRUE, ...) {
  if (sum(y >= a & y <= b) == 0) warning(paste("no values in y fall between a =", a, "and b =", b))
  p_y <- F_x(y, ...)
  if (a == -Inf) {
    p_a <- 0
  } else {
    p_a <- F_x(a, ...)
  }
  if (b == Inf) {
    p_b <- 1
  } else {
    p_b <- F_x(b, ...)
  }
  pit <- (p_y - p_a)/(p_b - p_a)
  if (a == -Inf) {
    pit[p_b == 0] <- 0
  } else if (b == Inf) {
    pit[p_a == 1] <- 1
  }
  pit[y <= a | y >= b] <- NA
  if (return_na) {
    return(pit)
  } else {
    return(na.omit(pit))
  }
}


