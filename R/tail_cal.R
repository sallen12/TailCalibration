#' Tail calibration of probabilistic forecasts
#'
#' Calculate marginal and probabilistic tail calibration of probabilistic forecasts
#'
#' @param y vector of realized values.
#' @param F_x forecast distribution function to be evaluated.
#' @param t threshold(s) at which to evaluate tail calibration.
#' @param z exceedances at which to calculate marginal tail calibration.
#' @param type string denoting whether to assess marginal or probabilistic tail calibration.
#' @param ... additional arguments to F_x.
#'
#' @return
#' vector of conditional probability integral transform (PIT) values (for \code{type = "prob"}),
#' or vector of differences (for \code{type = "marg"}).
#'
#' @references
#' \emph{Standard notions of calibration:}
#'
#' Gneiting, T. and Resin, J. (2022):
#' `Regression diagnostics meets forecast evaluation: Conditional calibration, reliability diagrams, and coefficient of determination',
#' \emph{arXiv pre-print}
#' \doi{arxiv.org/abs/2108.03210}
#'
#'
#' @author Sam Allen
#'
#' @details
#' Details to be added here.
#'
#' @examples
#' Examples to be added here.
#'
#' @name tail_cal
NULL

#' @rdname tail_cal
#' @export
tail_cal <- function(y, F_x, t, z = seq(0, 10, 0.1), type = "prob", ...) {
  check_inputs(y = y, F_x = F_x, t = t, z = z, type = type, ...)
  if (type == "prob") {
    tail_prob_cal(y, F_x, t, ...)
  } else if (type == "marg") {
    tail_marg_cal(y, F_x, t, z, ...)
  }
}

#' @rdname tail_cal
#' @export
tail_marg_cal <- function(y, F_x, t, z = seq(0, 10, 0.1), ...) {
  check_inputs(y = y, F_x = F_x, t = t, z = z, type = "marg", ...)
  d <- sapply(t, function(tt) {
    ind <- y > tt
    Fhat_t <- sapply(z, function(zz) (F_x(zz + tt, ...) - F_x(tt, ...))/(1 - F_x(tt, ...)))
    Fhat_t <- colMeans(Fhat_t[ind, ])
    if (sum(ind) == 0) {
      Qhat_t <- NA
    } else if (sum(ind) == 1) {
      Qhat_t <- (y[ind] - tt) <= z
    } else {
      Qhat_t <- colMeans(outer(y[ind] - tt, z, FUN = function(x1, x2) x1 <= x2))
    }
    return(Fhat_t - Qhat_t)
  })
  d <- as.data.frame(d)
  colnames(d) <- paste0("t", 1:length(t))
  d <- cbind(z, d)
  return(d)
}

#' @rdname tail_cal
#' @export
tail_prob_cal <- function(y, F_x, t, ...) {
  check_inputs(y = y, F_x = F_x, t = t, z = NULL, type = "prob", ...)
  cpit <- lapply(t, function(tt) {
    ind <- y > tt
    if (sum(ind) == 0) {
      cpit_t <- NA
    } else {
      cpit_t <- (F_x(y, ...) - F_x(tt, ...))/(1 - F_x(tt, ...))
    }
    return(data.frame(y = y[ind], cpit = cpit_t[ind]))
  })
  if (length(t) == 1) {
    cpit <- cpit[[1]]
  } else {
    names(cpit) <- round(t, 2)
  }
  return(cpit)
}

check_inputs <- function(y, F_x, t, z, type, ...) {
  if (!is.numeric(y)) stop("'y' is not numeric")
  if (!is.vector(y)) stop("'y' is not a vector")
  if (!is.function(F_x)) stop("'F_x' is not a function")
  if (!is.numeric(t)) stop("'t' is not numeric")

  if (!(type %in% c("prob", "marg"))) stop("'type' must be one of 'prob' or 'marg'")
  if (type == "marg") {
    if (!is.numeric(z)) stop("'z' is not numeric")
  }
}







