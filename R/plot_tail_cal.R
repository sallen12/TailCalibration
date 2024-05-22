#' Tail calibration of probabilistic forecasts
#'
#' Plot marginal and probabilistic tail calibration of probabilistic forecasts
#'
#' @param cal data frame or list of data frames containing the data to be plotted.
#' @param names names to be used in legend.
#' @param t threshold(s) corresponding to the data in \code{cal}.
#' @param xlims,ylims lower and upper limits of the axes.
#' @param xlab,ylab axes labels
#' @param title plot title.
#' @param type type of divergence used to measure miscalibration.
#'
#' @return
#' ggplot object
#'
#' @references
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
#' @import ggplot2 WeightedForecastVerification
#' @name plot_tail_cal
NULL


#' @rdname plot_tail_cal
#' @export
plot_mtc <- function(cal, names = NULL, ylims = c(-0.2, 0.2), xlims = NULL, ylab = NULL, xlab = NULL, title = NULL) {

  if (is.data.frame(cal)) {
    df <- cal
    mtc <- ggplot(df) + geom_line(aes(x = z, y = d))
  } else {
    df <- do.call(rbind, cal)
    if (is.null(names)) {
      df$names <- rep(names(cal), each = nrow(cal[[1]]))
    } else {
      df$names <- rep(names, each = nrow(cal[[1]]))
    }
    mtc <- ggplot(df) + geom_line(aes(x = z, y = d, col = as.factor(names)))
  }

  if (is.null(ylab)) ylab <- expression("E[" ~ F[t] ~ "(y)] - Q(Y - t ≤ y | Y > t)")
  if (is.null(xlab)) xlab <- "y"

  mtc <- mtc + geom_hline(aes(yintercept = 0), lty = "dotted") +
    scale_x_continuous(name = xlab, limits = xlims, expand = c(0, 0)) +
    scale_y_continuous(name = ylab, limits = ylims, expand = c(0, 0)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.position = "bottom") +
    ggtitle(title)

  return(mtc)
}


#' @rdname plot_tail_cal
#' @export
plot_ptc <- function(cal, names = NULL, title = NULL) {
  if (is.data.frame(cal)) {
    cal <- cal$cpit
    ptc <- pit_reldiag(cal, resampling = T, title = title)
  } else {
    cal <- lapply(cal, function(z) z$cpit)
    if (!is.null(names)) names(cal) <- names
    ptc <- pit_reldiag(cal, resampling = F, title = title)
  }

  return(ptc)
}


#' @rdname plot_tail_cal
#' @export
plot_ptc_div <- function(cal, t, names = NULL, ylims = NULL, xlims = NULL,
                         ylab = "Miscalibration", xlab = "Threshold", title = NULL,
                         type = c("sup", "cramer")) {

  if (is.data.frame(cal[[1]])) {
    if (type == "cramer") {
      div <- sapply(seq_along(t), function(i) crps_div(cal[[i]]$cpit))
    } else {
      div <- sapply(seq_along(t), function(i) sup_div(cal[[i]]$cpit))
    }
    df <- data.frame(t = t, d = div)
    ptc_div <- ggplot(df) + geom_line(aes(x = t, y = d))
  } else {
    if (type == "cramer") {
      div <- sapply(cal, function(x) sapply(seq_along(t), function(i) crps_div(x[[i]]$cpit)))
    } else {
      div <- sapply(cal, function(x) sapply(seq_along(t), function(i) sup_div(x[[i]]$cpit)))
    }
    if (is.null(names)) names <- colnames(div)
    df <- data.frame(t = t, d = c(div),
                     mth = rep(names, each = length(t)))
    ptc_div <- ggplot(df) + geom_line(aes(x = t, y = d, col = mth))
  }

  ptc_div <- ptc_div + geom_hline(aes(yintercept = 0), lty = "dotted") +
    scale_x_continuous(name = xlab, limits = xlims, expand = c(0, 0)) +
    scale_y_continuous(name = ylab, limits = ylims) +
    theme_bw() + theme(panel.grid = element_blank(),
                       legend.title = element_blank(),
                       legend.position = "bottom") +
    ggtitle(title)

  return(ptc_div)
}


sup_div <- function(z, u = seq(0, 1, 0.01)) {
  F_x <- sapply(u, function(uu) mean(z <= uu))
  max(abs(F_x - u))
}
