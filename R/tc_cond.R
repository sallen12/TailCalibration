#' Conditional tail calibration
#'
#' Assess the tail calibration of probabilistic forecasts conditional on non-trivial
#' information sets.
#'
#' @inheritParams tail_cal
#'
#' @return
#' A list of data frames (for each group) or a list of a list of data frames (for each group
#' and threshold \code{t}) containing the requested ratio quantifying tail calibration.
#'
#'
#' @section Details:
#'
#' \code{tc_cond()} is a wrapper for \code{tc_cprob} if \code{type == 'prob'} and
#' \code{tc_cmarg} if \code{type == 'marg'}, which are themselves wrappers for
#' \code{\link{tc_prob}} and \code{\link{tc_marg}}, respectively. More details can be
#' found in their respective help pages.
#'
#' \code{tc_cprob} and \code{tc_cmarg} stratify the forecasts and observations into
#' different groups, according to the input vector \code{group}. \code{group} should
#' be a vector of the same length as \code{y}, containing a relatively small number of
#' unique entries. Observations \code{y} that have the same value of \code{group} are
#' grouped together, and the standard notions of probabilistic (\code{type == 'prob'})
#' or marginal tail calibration (\code{type == 'marg'}) are then assessed using the
#' functionality in \code{tc_prob} or \code{tc_marg}, respectively.
#'
#'
#' @inheritSection tc_prob References
#' @inheritSection tc_marg References
#'
#'
#' @author Sam Allen
#'
#'
#' @examples
#' n <- 1e5
#' y <- rnorm(n)
#' F_x <- pnorm
#' mu <- rnorm(n, sd = 0.5)
#'
#' # standard calibration is assessed by setting t to -Inf
#' ptc <- tc_prob(y, F_x, t = -Inf, mean = mu)
#' ptc_cond <- tc_cprob(y, F_x, t = -Inf, group = rep(1, length(y)), sup = FALSE, mean = mu)
#' ptc_cond2 <- tc_cond(y, F_x, t = -Inf, sup = FALSE, mean = mu)
#'
#' mtc <- tc_marg(y, F_x, t = -Inf, mean = mu)
#' mtc_cond <- tc_cmarg(y, F_x, t = -Inf, group = rep(1, length(y)), sup = FALSE, mean = mu)
#' mtc_cond2 <- tc_cond(y, F_x, t = -Inf, type = 'marg', sup = FALSE, mean = mu)
#'
#' # we can condition on arbitrary groupings of the observations
#' group <- sample(c('A', 'B'), length(y), replace = TRUE)
#' t <- c(-Inf, 0, 1, 2)
#' sev <- tc_cond(y, F_x, t = t, ratio = 'sev', group = group, mean = mu)
#' occ <- tc_cond(y, F_x, t = t, ratio = 'occ', group = group, mean = mu)
#'
#'
#' @name tc_cond
NULL


#' @rdname tc_cond
#' @export
tc_cond <- function(y, F_x, t, type = c('prob', 'marg'), ratio = c('com', 'sev', 'occ'),
                    u = seq(0.01, 0.99, 0.01), group = rep(1, length(y)), lower = -Inf,
                    sup = TRUE, qu = FALSE, subset = rep(TRUE, length(y)),  ...) {
  type <- match.arg(type)
  ratio <- match.arg(ratio)
  if (type == 'prob') {
    tc_cprob(y, F_x, t, group = group, ratio = ratio, u = u, lower = lower, sup = sup, qu = qu, subset = subset, ...)
  } else if (type == 'marg') {
    tc_cmarg(y, F_x, t, group = group, ratio = ratio, u = u, sup = sup, qu = qu, subset = subset, ...)
  }
}


#' @rdname tc_cond
#' @export
tc_cprob <- function(y, F_x, t, group, ratio = c("com", "sev", "occ"), u = seq(0, 1, 0.01),
                     lower = -Inf, sup = TRUE, qu = FALSE, subset = rep(TRUE, length(y)), ...) {
  check_tc_inputs(y, F_x, t, u = u, group = group, sup = sup, qu = qu, subset = subset)
  ratio <- match.arg(ratio)

  grps <- unique(group)

  R_g <- lapply(grps, function(g) {
    ind <- (group == g) & subset
    tc_prob(y, F_x, t, ratio = ratio, u = u, lower = lower, sup = sup, qu = qu, subset = ind, ...)
  })

  names(R_g) <- grps
  if (length(grps) == 1) R_g <- R_g[[1]]

  return(R_g)
}


#' @rdname tc_cond
#' @export
tc_cmarg <- function(y, F_x, t, group, ratio = c("com", "sev", "occ"), u = seq(0, 1, 0.01),
                     sup = TRUE, qu = FALSE, subset = rep(TRUE, length(y)), ...) {
  check_tc_inputs(y, F_x, t, u = u, group = group, sup = sup, qu = qu, subset = subset)
  ratio <- match.arg(ratio)

  grps <- unique(group)

  R_g <- lapply(grps, function(g) {
    ind <- (group == g) & subset
    tc_marg(y, F_x, t, ratio = ratio, u = u, sup = sup, qu = qu, subset = ind, ...)
  })

  names(R_g) <- grps
  if (length(grps) == 1) R_g <- R_g[[1]]

  return(R_g)
}

