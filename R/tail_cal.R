#' Tail calibration of probabilistic forecasts
#'
#' Assess marginal and probabilistic tail calibration of probabilistic forecasts.
#' This allows forecast calibration to be assessed when interest is on values above
#' thresholds of interest.
#'
#' @param y vector of observations.
#' @param F_x forecast distribution function to be evaluated, or a vector or matrix of samples from the forecast distribution.
#' @param t vector of threshold(s) at which to evaluate tail calibration.
#' @param type string denoting whether to assess marginal ('marg') or probabilistic ('prob') tail calibration.
#' @param ratio ratio to return; one of 'com', 'sev', 'occ'.
#' @param u vector of values at which to assess tail calibration.
#' @param lower numeric value at which the forecast distributions are censored.
#' @param group vector specifying observations that should be grouped together.
#' @param sup logical specifying whether to quantify miscalibration using the supremum distance; default is \code{FALSE}.
#' @param qu logical specifying whether the thresholds should be renamed in terms of
#'  quantiles of \code{y}; default is \code{FALSE}.
#' @param subset logical vector of the same length as \code{y}, allowing only a subset
#'  of the forecasts and observations (where \code{TRUE}) to be assessed.
#' @param ... additional arguments to F_x.
#'
#'
#' @section Details:
#'
#' \code{tail_cal()} is a wrapper for \code{\link{tc_prob}} if \code{type == 'prob'} and
#' \code{\link{tc_marg}} if \code{type == 'marg'}. If \code{group} is not \code{NULL},
#' then \code{tail_cal()} is a wrapper for \code{\link{tc_cond}}.
#'
#' Forecast evaluation is typically performed using a sequence of forecasts \eqn{F_{i}}
#' and corresponding observations \eqn{y_{i}}, for \eqn{i = 1, \dots, n}. The forecasts
#' are said to be calibrated if the forecasts align statistically with the observations.
#' Several notions of forecast calibration exist. Forecast tail calibration allows
#' the calibration of forecasts to be assessed when interest is on values that exceed
#' a high threshold.
#'
#'
#' @inheritSection tc_prob Details
#' @inheritSection tc_marg Details
#' @inheritSection tc_cond Details
#'
#'
#' @inheritSection tc_prob Value
#'
#'
#' @seealso \code{\link{tc_prob}} \code{\link{tc_marg}}
#'
#'
#' @inheritSection tc_prob References
#' @inheritSection tc_marg References
#'
#'
#' @author Sam Allen
#'
#' @examples
#' n <- 1e5
#' y <- rnorm(n)
#' F_x <- pnorm
#' mu <- rnorm(n, sd = 0.5)
#'
#' # standard calibration is assessed by setting t to -Inf
#' ptc <- tail_cal(y, F_x, t = -Inf, mean = mu)
#' mtc <- tail_cal(y, F_x, t = -Inf, type = "marg", mean = mu)
#'
#' t <- c(-Inf, 0, 1, 2)
#' com <- tail_cal(y, F_x, t = t, sup = TRUE, mean = mu)
#' com2 <- tc_prob(y, F_x, t = t, sup = TRUE, mean = mu)
#' identical(com, com2)
#'
#' sev <- tail_cal(y, F_x, t = t, type = 'marg', ratio = 'sev', mean = mu)
#' sev2 <- tc_marg(y, F_x, t = t, ratio = 'sev', mean = mu)
#' identical(sev, sev2)
#'
#'
#' # ensemble forecast
#' M <- 100
#' F_x <- matrix(rnorm(n*M, mu), nrow = n, ncol = M)
#' ptc_ens <- tail_cal(y, F_x, t = -Inf, mean = mu)
#' plot_ptc(ptc)
#' plot_ptc(ptc_ens)
#'
#'
#' @name tail_cal
NULL


#' @rdname tail_cal
#' @export
tail_cal <- function(y, F_x, t, type = c('prob', 'marg'), ratio = c('com', 'sev', 'occ'),
                     u = seq(0.01, 0.99, 0.01), lower = -Inf, group = NULL, sup = FALSE, qu = FALSE,
                     subset = rep(TRUE, length(y)),  ...) {
  type <- match.arg(type)
  ratio <- match.arg(ratio)
  if (is.null(group)) {
    if (type == 'prob') {
      tc_prob(y, F_x, t, ratio = ratio, u = u, lower = lower, sup = sup, qu = qu, subset = subset, ...)
    } else if (type == 'marg') {
      tc_marg(y, F_x, t, ratio = ratio, u = u, sup = sup, qu = qu, subset = subset, ...)
    }
  } else {
    tc_cond(y, F_x, t, type = type, ratio = ratio, u = u, group = group, sup = sup, qu = qu, subset = subset, ...)
  }

}


check_tc_inputs <- function(y, F_x, t, u, group, sup, qu, subset) {
  if (!is.numeric(y) || !is.vector(y)) stop("'y' must be a numeric value or vector")
  if (any(is.na(y))) stop("'y' contains missing values")
  if (!is.numeric(t) || !is.vector(t)) stop("'t' must be a numeric value or vector")
  if (!is.numeric(u) || !is.vector(u)) stop("'u' must be a numeric value or vector")
  if (!is.null(group)) {
    if (!is.vector(group) || length(group) != length(y))
      stop("'group' must be a vector of the same length as 'y'")
    if (length(unique(group)) > length(y)/2)
      warning("the number of unique elements in 'group' is large relative to the number of elements in 'y'")
  }
  if (!is.logical(sup) || length(sup) > 1) stop("'sup' must be either TRUE or FALSE")
  if (!is.logical(qu) || length(qu) > 1) stop("'qu' must be either TRUE or FALSE")
  if (!is.logical(subset) || !is.vector(subset) || length(subset) != length(y))
    stop("'subset' must be a logical value or vector of the same length as 'y'")
  if (!is.function(F_x) && !is.matrix(F_x) && !is.vector(F_x)) stop("'F_x' must be a function or a numeric vector or matrix")
  if (is.matrix(F_x)) {
    if (nrow(F_x) != length(y)) stop("the dimensions of 'F_x' do not match the dimensions of 'y'")
  } else if (is.vector(F_x)) {
    if (length(y) > 1) stop("the dimensions of 'F_x' do not match the dimensions of 'y'")
  }
}

