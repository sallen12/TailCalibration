#' Conditional PIT values for simulated forecast distributions
#'
#' Calculate conditional probability integral transform (PIT) values for probabilistic
#' forecasts in the form of samples or ensembles.
#'
#' @param y vector of observations.
#' @param dat matrix of simulation draws from forecast distribution.
#' @param a,b numeric; lower and upper threshold defining the interval of outcomes that are of interest.
#' @param bw vector of bandwidth parameters to use in kernel density estimation.
#'
#' @details
#' \code{cpit_sample()} calculates conditional PIT values for forecast distributions
#' in the form of an ensemble or predictive sample. See the documentation for
#' \code{\link{cpit_param}} for details on conditional PIT histograms.
#'
#' Calculating conditional PIT values requires the conditional predictive distribution
#' given that the outcome is in the range [\code{a}, \code{b}]. This is achieved here
#' by first smoothing the predictive sample using kernel density estimation.
#'
#' The argument \code{bw} specifies the bandwidth parameter to use within kernel
#' density estimation. If this is not provided, then it is selected using the core
#' function \code{\link{bw.nrd}}.
#'
#' The argument \code{dat} contains the sample from the forecast distribution. This
#' should be a matrix with the same number of rows as the length of \code{y}, with
#' the columns representing different samples or ensemble members.
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
#' \emph{Journal of the Royal Statistical Society: Series A (General)} 147, 278-290.
#' \doi{10.2307/2981683}
#'
#' \emph{Conditional PIT histograms}
#'
#' Allen, S., Bhend, J., Martius, O. and Ziegel, J (2023):
#' `Weighted verification tools to evaluate univariate and multivariate probabilistic forecasts for high-impact weather events'.
#' \emph{Weather and Forecasting} 38, 499–516.
#' \doi{10.1175/WAF-D-22-0161.1}
#'
#' @examples
#' mu <- rnorm(20, mean = 0, sd = 5)
#' y <- rnorm(20, mean = mu, sd = 1)
#' dat <- matrix(rnorm(2000, mean = mu, sd = 1), nrow = 20)
#'
#' cpit_sample(y = y, dat = dat)
#' cpit_norm(y = y, mean = mu, sd = 1)
#'
#' cpit_sample(y = y, dat = dat, a = 0)
#' cpit_norm(y = y, mean = mu, sd = 1, a = 0)
#'
#' @export
cpit_sample <- function(y, dat, a = -Inf, b = Inf, bw = NULL, return_na = TRUE){
  if (is.null(bw)) {
    message("Computing conditional PIT values using kernel density estimation")
    bw <- apply(dat, 1, bw.nrd)
    pit <- sapply(seq_along(y), function(i) pkde(y[i], m_vec = dat[i, ], s_vec = bw[i], a, b))
  } else {
    pit <- pens(y, dat, a, b)
  }
  pit[y <= a | y >= b] <- NA
  if (return_na) {
    return(pit)
  } else {
    return(na.omit(pit))
  }
}

# distribution function for normal kernel density estimation
pkde <- function(q, m_vec, s_vec, a = -Inf, b = Inf){
  p_q <- mean(pnorm(q, m_vec, s_vec))
  p_a <- mean(pnorm(a, m_vec, s_vec))
  p_b <- mean(pnorm(b, m_vec, s_vec))
  p <- (p_q - p_a)/(p_b - p_a)
  if (a == -Inf) {
    p[p_b == 0] <- 0
  } else if (b == Inf) {
    p[p_a == 1] <- 1
  }
  return(p)
}

# distribution function for samples
pens <- function(q, dat, a = -Inf, b = Inf){
  if (is.vector(dat)) {
    p_q <- mean(dat <= q)
    p_a <- mean(dat <= a)
    p_b <- mean(dat <= b)
  } else {
    p_q <- rowMeans(dat <= q)
    p_a <- rowMeans(dat <= a)
    p_b <- rowMeans(dat <= b)
  }
  p <- (p_q - p_a)/(p_b - p_a)
  if (a == -Inf) {
    p[p_b == 0] <- 0
  } else if (b == Inf) {
    p[p_a == 1] <- 1
  }
  return(p)
}
