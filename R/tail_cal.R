#' Tail calibration of probabilistic forecasts
#'
#' Calculate marginal and probabilistic tail calibration of probabilistic forecasts
#'
#' @param y vector of realized values.
#' @param F_x forecast distribution function to be evaluated.
#' @param t threshold(s) at which to evaluate tail calibration.
#' @param ratio ratio to return; one of 'comb', 'sev', 'occ'.
#' @param z exceedances at which to calculate marginal tail calibration.
#' @param type string denoting whether to assess marginal or probabilistic tail calibration.
#' @param sup logical specifying whether to take the supremum of the differences in marginal tail calibration.
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
  if (type == "prob") {
    tail_prob_cal(y, F_x, t, ...)
  } else if (type == "marg") {
    tail_marg_cal(y, F_x, t, z, ...)
  }
}


################################################################################
##### probabilistic tail calibration


#' @rdname tail_cal
#' @export
tail_prob_cal <- function(y, F_x, t, ratio = c("comb", "sev", "occ"), u = seq(0, 1, 0.01),
                          sup = FALSE, qu = FALSE, subset = NULL, ...) {
  check_inputs(y = y, F_x = F_x, t = t, z = NULL, type = "prob", ...)
  ratio <- match.arg(ratio)

  if (is.null(subset)) subset <- rep(TRUE, length(y))
  n <- sum(subset)

  if (ratio == "comb") {
    R <- lapply(t, function(tt) {
      exc_p <- 1 - F_x(tt, ...)
      if (length(exc_p) > 1) exc_p <- mean(exc_p[subset])
      cpit <- cpit_dist(y, F_x, a = tt, ...)
      cpit <- na.omit(cpit[subset])
      rat <- sapply(u, function(uu) sum(cpit <= uu)/(n*exc_p))
      data.frame(u = u, rat = rat)
    })
  } else if (ratio == "sev") {
    R <- lapply(t, function(tt) {
      cpit <- cpit_dist(y, F_x, a = tt, ...)
      cpit <- na.omit(cpit[subset])
      rat <- sapply(u, function(uu) mean(cpit <= uu))
      data.frame(u = u, rat = rat)
    })
  } else if (ratio == "occ") {
    R <- tail_cal_occ(y, F_x, t, type = "ratio", qu = qu, subset = subset, ...)
  }

  if (sup) {
    if (ratio %in% c("comb", "sev")) {
      R <- sapply(R, function(r) max(abs(r$rat - r$u)))
    } else {
      R <- abs(R$rat - 1)
    }
  }

  if (length(t) == 1) {
    R <- R[[1]]
  } else if (length(y) == 1) {
    R <- as.data.frame(t(sapply(seq_along(t), function(i) R[[i]][1, ])))
  } else if (ratio %in% c("comb", "sev")) {
    names(R) <- round(t, 2)
  }
  return(R)
}


#' @rdname tail_cal
#' @export
tail_cal_occ <- function(y, F_x, t, type = c("all", "ratio", "dif"), qu = FALSE, subset = NULL, ...) {
  type <- match.arg(type)
  if (is.null(subset)) subset <- rep(TRUE, length(y))
  G_t <- sapply(t, function(tt) mean(y[subset] > tt))
  F_t <- sapply(t, function(tt) 1 - F_x(tt, ...))
  if (is.matrix(F_t)) F_t <- colMeans(F_t[subset, ])
  if (qu) t <- 1 - G_t
  df <- data.frame(t = t, G_t = G_t, F_t = F_t, rat = G_t/F_t, dif = G_t - F_t)
  if (type == "all") {
    return(df)
  } else if (type == "ratio") {
    return(df[c("t", "rat")])
  } else if (type == "dif") {
    return(df[c("t", "dif")])
  }
}


################################################################################
##### conditional probabilistic tail calibration


#' @rdname tail_cal
#' @export
tail_prob_div <- function(y, F_x, t, group, ratio = c("comb", "sev", "occ"), u = seq(0, 1, 0.01), sup = TRUE, qu = FALSE, ...) {
  check_inputs(y = y, F_x = F_x, t = t, z = NULL, type = "prob", ...)
  ratio <- match.arg(ratio)

  grps <- unique(group)

  R_g <- lapply(grps, function(g) {
    ind <- group == g
    tail_prob_cal(y, F_x, t, ratio = ratio, u = u, sup = sup, qu = qu, subset = ind, ...)
  })

  if (sup) {
    if (qu) t <- sapply(t, function(tt) mean(y <= tt))
    R_g <- lapply(R_g, function(x) data.frame(t = t, d = x))
  }

  names(R_g) <- grps

  return(R_g)
}


################################################################################
##### marginal tail calibration


#' @rdname tail_cal
#' @export
tail_marg_cal <- function(y, F_x, t, ratio = c("comb", "sev", "occ"), z = seq(0, 10, 0.1), sup = FALSE, ...) {
  check_inputs(y = y, F_x = F_x, t = t, z = z, type = "marg", ...)

  if (ratio == "comb") {
    D <- lapply(t, function(tt) {
      F_t <- F_x(tt, ...)
      exc_p <- 1 - mean(F_tt)
      ind <- y > tt
      dif <- sapply(z, function(zz) {
        cprob <- mean(y > tt & y <= (tt + zz))
        Fhat_t <- (F_x(zz + tt, ...) - F_t)/(1 - F_t)
        Fhat_t <- mean(Fhat_t[ind])
        (cprob/exc_p) - Fhat_t
      })
      data.frame(x = z, dif = dif)
    })
  } else if (ratio == "sev") {
    D <- lapply(t, function(tt) {
      F_t <- F_x(tt, ...)
      ind <- y > tt
      dif <- sapply(z, function(zz) {
        cprob <- mean(y[ind] <= (tt + zz))
        Fhat_t <- (F_x(zz + tt, ...) - F_t)/(1 - F_t)
        Fhat_t <- mean(Fhat_t[ind])
        cprob - Fhat_t
      })
      data.frame(x = z, dif = dif)
    })
  } else if (ratio == "occ") {
    D <- tail_cal_occ(y, F_x, t, type = "dif", ...)
  }

  if (sup) {
    if (ratio %in% c("comb", "sev")) {
      D <- sapply(D, function(d) max(abs(d$dif)))
    } else {
      D <- max(abs(D$dif))
    }
  }

  if (length(t) == 1) {
    D <- D[[1]]
  } else if (length(y) == 1) {
    D <- as.data.frame(t(sapply(seq_along(t), function(i) D[[i]][1, ])))
  } else {
    names(D) <- round(t, 2)
  }
  return(D)
}


################################################################################
##### conditional marginal tail calibration


################################################################################
##### checks


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















