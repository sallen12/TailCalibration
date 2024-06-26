#' Marginal tail calibration
#'
#' Assess the marginal tail calibration of probabilistic forecasts. This allows marginal
#' forecast calibration to be assessed when interest is on values above thresholds of
#' interest.
#'
#' @inheritParams tail_cal
#'
#' @inheritSection tc_prob Value
#'
#'
#' @section Details:
#'
#' \emph{Marginal tail calibration}
#'
#' Forecasts \eqn{F_{i}} for observations \eqn{y_{i}} are \emph{marginally calibrated}
#' if average forecast probabilities are the same as unconditional observed
#' probabilities. That is, for all \eqn{u \in (-\infty, \infty)},
#' \deqn{\frac{1}{n} \sum_{i=1}^{n} 1 \{y_{i} \le u\} \approxeq \frac{1}{n} \sum_{i=1}^{n} F_{i}(u),}
#' where \eqn{1\{ \cdot \}} denotes the indicator function.
#' To assess this, we can plot the difference between these two terms as a function of \eqn{u}.
#'
#' Forecasts are said to be \emph{marginally tail calibrated} if, for all \eqn{u \in [0, \infty)},
#' \deqn{\frac{\sum_{i=1}^{n} 1 \{t < y_{i} \le u + t\} }{\sum_{i=1}^{n}(1 - F_{i}(t))} \approxeq \frac{1}{n_{t}} \sum_{i=1}^{n} F_{t,i}(u) 1 \{y_{i} > t\},}
#' as the threshold \eqn{t} tends to the upper end point of the distribution of the observations,
#' where \eqn{n_{t}} is the number of observations that exceed \eqn{t}, and
#' \deqn{F_{t,i}(u) = \frac{F_{i}(u + t) - F_{i}(t)}{1 - F_{i}(t)}}
#' for \eqn{u \ge 0} is the forecast excess distribution.
#' To assess marginal tail calibration, we can plot the difference between the left and right hand sides as a
#' function of \eqn{u}.
#'
#' We refer to the difference between these terms as the \emph{combined ratio}, since it can be
#' decomposed into two more intuitive conditions:
#' \deqn{\frac{\sum_{i=1}^{n} 1 \{t < y_{i} \le u + t\} }{\sum_{i=1}^{n} 1 \{y_{i} > t\}} \approxeq \frac{1}{n_{t}} \sum_{i=1}^{n} F_{t,i}(u) 1 \{y_{i} > t\},}
#' for all \eqn{u \in [0, \infty)}; and
#' \deqn{\frac{\sum_{i=1}^{n} 1 \{y_{i} > t\}}{\sum_{i=1}^{n}(1 - F_{i}(t))} \approxeq 1.}
#'
#' The first of these two conditions assesses whether the average forecast excess distribution
#' is equal to the unconditional excess distribution of the observations, conditional on
#' the outcomes \eqn{y_{i}} being larger than threshold \eqn{t}. We therefore refer to the
#' differences between these terms as the \emph{severity ratio}, since it only assesses
#' forecasts for the severity of the extreme event and does not consider the probability
#' assigned to being above the threshold.
#'
#' The second of these two conditions evaluates the forecast's ability to predict the
#' occurrence of a threshold exceedance, by comparing the average forecast probability
#' that the outcome will exceed the threshold with the unconditional observed probability
#' of the outcome exceeding the threshold. This ratio should be close to one if the
#' forecast is calibrated. We refer to the ratio in this condition as the
#' \emph{occurrence} ratio, since it only assesses the forecast probability for
#' extreme event occurrence and does not consider by how much the outcome exceeds
#' the threshold.
#'
#' The function \code{tc_marg()} takes the observations \eqn{y_{i}} and forecast
#' distributions \eqn{F_{i}} as inputs (via the arguments \code{y} and \code{F_x}),
#' and returns one of these three ratios: the combined ratio if \code{ratio == 'com'}
#' (the default), the severity ratio if \code{ratio == 'sev'}, and the occurrence
#' ratio if \code{ratio == 'occ'}. The vector \code{u} specifies the values
#' \eqn{u \in [0, \infty)} at which to assess the combined and severity ratios, and
#' the numeric \code{t} specifies the threshold(s) of interest. \code{t} can either be a
#' single value or a vector.
#'
#' For \code{ratio == 'com'} or \code{ratio == 'sev'}, \code{tc_marg()} returns a data frame
#' containing the vector \code{u} and the chosen ratio evaluated at each element in
#' this vector. If \code{t} is a vector, a list of data frames is returned, with each
#' element of the list corresponding to a different threshold. For \code{ratio == 'occ'},
#' \code{tc_marg()} returns a single value containing the occurrence ratio at this threshold.
#' If \code{t} is a vector, a data frame is instead returned that contains the thresholds \eqn{t}
#' and the corresponding occurrence ratios.
#'
#' If \code{sup = TRUE}, then, for each threshold, the maximum absolute distance between the
#' ratio and the expected output for a calibrated forecast is returned.
#' For example, for \code{ratio == 'com'}, \code{tc_marg()} calculates
#' \deqn{\max_{u} | \frac{\sum_{i=1}^{n} 1 \{t < y_{i} \le u + t\} }{\sum_{i=1}^{n}(1 - F_{i}(t))} - \frac{1}{n_{t}} \sum_{i=1}^{n} F_{t,i}(u) 1 \{y_{i} > t\} |,}
#' with the maximum taken over all elements in the argument \code{u}. This single value is
#' returned, or, if \code{t} is a vector, a data frame is returned containing the thresholds
#' \code{t} and the corresponding maximum absolute distance.
#'
#' If \code{qu = TRUE}, then the thresholds are re-expressed in terms of quantiles of the observation
#' vector \code{y}. This is only relevant for the presentation of the output.
#'
#' Additional arguments to \code{F_x} can be specified using \code{...}. \code{tc_marg()}
#' assumes that the distribution function is the same in all forecast cases, and that
#' \code{...} contains parameters corresponding to this distribution function.
#' If this is not the case, then a work around is to define \code{F_x} as a regime-switching
#' distribution, with an additional parameter specifying which distribution/regime to choose;
#' see examples.
#'
#'
#' @seealso \code{\link{tail_cal}} \code{\link{tc_prob}}
#'
#'
#' @section References:
#' \emph{Marginal calibration:}
#'
#' Gneiting, T., Balabdaoui, F., and Raftery, A. E. (2007):
#' `Probabilistic forecasts, calibration and sharpness'
#' \emph{Journal of the Royal Statistical Society Series B: Statistical Methodology}, 69, 243-268.
#' DOI: 10.1111/j.1467-9868.2007.00587.x
#'
#'
#' Gneiting, T. and Resin, J. (2022):
#' `Regression diagnostics meets forecast evaluation: Conditional calibration, reliability diagrams, and coefficient of determination',
#' \emph{Electronic Journal of Statistics}, 17, 3226-3286.
#' DOI: 10.1214/23-EJS2180
#'
#'
#' \emph{Marginal tail calibration:}
#'
#' Allen, S., Koh, J., Segers, J. and Ziegel, J. (2024+):
#' `Tail calibration of probabilistic forecasts'.
#' \emph{arXiv pre-print}.
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
#' # standard marginal calibration is assessed by setting t to -Inf
#' mtc <- tc_marg(y, F_x, t = -Inf, mean = mu)
#'
#' com <- tc_marg(y, F_x, t = c(-Inf, 0, 1, 2), mean = mu)
#' sev <- tc_marg(y, F_x, t = c(-Inf, 0, 1, 2), ratio = 'sev', mean = mu)
#' occ <- tc_marg(y, F_x, t = c(-Inf, 0, 1, 2), ratio = 'occ', mean = mu)
#'
#' sup <- tc_marg(y, F_x, t = c(-Inf, 0, 1, 2), sup = TRUE, qu = TRUE, mean = mu)
#' sevsup <- tc_marg(y, F_x, t = c(-Inf, 0, 1, 2), ratio = 'sev', sup = TRUE, qu = TRUE, mean = mu)
#'
#'
#' # mixture between a standard normal and logistic distribution
#' F_x <- function(x, mu1, mu2, lo) plogis(x, mu1)*lo + pnorm(x, mu2)*(1 - lo)
#'
#' mu1 <- rnorm(mu + 1, sd = 0.5)
#' mu2 <- mu
#' lo <- mu > 1
#'
#' mtc <- tc_marg(y, F_x, t = c(-2, 0, 2), mu1 = mu1, mu2 = mu2, lo = lo)
#'
#'
#' @name tc_marg
NULL


#' @rdname tc_marg
#' @export
tc_marg <- function(y, F_x, t, ratio = c('com', 'sev', 'occ'), u = seq(0, 10, 0.1),
                   sup = FALSE, qu = FALSE, subset = rep(TRUE, length(y)), ...) {
  check_tc_inputs(y, F_x, t, u = u, group = NULL, sup = sup, qu = qu, subset = subset)
  ratio <- match.arg(ratio)
  if (!is.function(F_x)) {
    dat <- F_x
    if (is.vector(dat)) {
      F_x <- function(x, ...) mean(dat <= x)
    } else {
      F_x <- function(x, ...) rowMeans(dat <= x)
    }
  }

  n <- sum(subset)

  if (ratio == 'com') {
    R <- lapply(t, function(tt) {
      F_t <- F_x(tt, ...)
      exc_p <- 1 - F_t
      if (length(exc_p) > 1) exc_p <- mean(exc_p[subset])
      ind <- (y > tt) & subset
      if (tt == -Inf) {
        dif <- sapply(u, function(uu) {
          cprob <- mean(y[subset] <= uu)
          Fhat_t <- F_x(uu, ...)
          Fhat_t <- mean(Fhat_t[ind])
          (cprob/exc_p) - Fhat_t
        })
      } else {
        dif <- sapply(u, function(uu) {
          cprob <- mean(y[subset] > tt & y[subset] <= (tt + uu))
          Fhat_t <- (F_x(uu + tt, ...) - F_t)/(1 - F_t)
          Fhat_t <- mean(Fhat_t[ind])
          (cprob/exc_p) - Fhat_t
        })
      }
      data.frame(u = u, rat = dif)
    })
  } else if (ratio == 'sev') {
    R <- lapply(t, function(tt) {
      F_t <- F_x(tt, ...)
      ind <- (y > tt) & subset
      if (tt == -Inf) {
        dif <- sapply(u, function(uu) {
          cprob <- mean(y[ind] <= uu)
          Fhat_t <- F_x(uu, ...)
          Fhat_t <- mean(Fhat_t[ind])
          cprob - Fhat_t
        })
      } else {
        dif <- sapply(u, function(uu) {
          cprob <- mean(y[ind] <= (tt + uu))
          Fhat_t <- (F_x(uu + tt, ...) - F_t)/(1 - F_t)
          Fhat_t <- mean(Fhat_t[ind])
          cprob - Fhat_t
        })
      }
      data.frame(u = u, rat = dif)
    })
  } else if (ratio == 'occ') {
    G_t <- sapply(t, function(tt) mean(y[subset] > tt))
    F_t <- sapply(t, function(tt) 1 - F_x(tt, ...))
    if (is.matrix(F_t)) F_t <- colMeans(F_t[subset, ])
    if (qu) t <- 1 - G_t
    R <- G_t/F_t
    if (length(t) > 1) R <- data.frame(t = t, rat = R)
  }

  if (qu) t <- sapply(t, function(tt) mean(y[subset] <= tt))

  if (sup) {
    if (ratio %in% c('com', 'sev')) {
      R <- sapply(R, function(r) max(abs(r$rat)))
    } else {
      R <- abs(R$rat - 1)
    }
    if (length(t) > 1) R <- data.frame(t = t, rat = R)
  } else if (length(t) == 1) {
    R <- R[[1]]
  } else if (ratio %in% c('com', 'sev')) {
    names(R) <- round(t, 2)
  }

  return(R)
}

