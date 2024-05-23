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


################################################################################
##### probabilistic tail calibration


#' @rdname plot_tail_cal
#' @export
plot_ptc <- function(cal, ratio = c("com", "sev", "occ"), names = NULL,
                     xlab = NULL, ylab = NULL, xlims = NULL, ylims = NULL, title = NULL) {
  ratio <- match.arg(ratio)

  if (!is.data.frame(cal)) {
    if (!is.null(names)) names(cal) <- names
  }

  if (ratio == "com") {
    tc <- plot_ptc_com(cal, xlab, ylab, xlims, ylims, title)
  } else if (ratio == "sev") {
    tc <- plot_ptc_sev(cal, xlab, ylab, xlims, ylims, title)
  } else if (ratio == "occ") {
    tc <- plot_tc_occ(cal, xlab, ylab, xlims, ylims, title)
  }

  return(tc)
}


#' @rdname plot_tail_cal
#' @export
plot_ptc_com <- function(cal, names = NULL, xlab = NULL, ylab = NULL, xlims = NULL, ylims = NULL, title = NULL) {

  if (is.null(ylab)) ylab <- "Combined ratio"

  tc <- plot_ptc_sev(cal, names, xlab, ylab, xlims, ylims, title)

  return(tc)
}


#' @rdname plot_tail_cal
#' @export
plot_ptc_sev <- function(cal, names = NULL, xlab = NULL, ylab = NULL, xlims = NULL, ylims = NULL, title = NULL) {

  if (is.null(xlab)) xlab <- "u"
  if (is.null(ylab)) ylab <- "Severity ratio"
  if (is.null(xlims)) xlims <- c(0, 1)
  if (is.null(ylims)) ylims <- c(0, 1)

  if (is.data.frame(cal)) {
    tc <- ggplot(cal) + geom_line(aes(x = u, y = rat))
  } else {
    if (!is.null(names)) names(cal) <- names
    df <- do.call(rbind, cal)
    df$mth <- rep(names(cal), sapply(cal, nrow))
    tc <- ggplot(df) + geom_line(aes(x = u, y = rat, col = mth))
  }

  tc <- tc + geom_abline(aes(intercept = 0, slope = 1), linetype = "dotted") +
    scale_x_continuous(name = xlab, limits = xlims, expand = c(0, 0)) +
    scale_y_continuous(name = ylab, limits = ylims, expand = c(0, 0)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.justification = c(0, 1),
          legend.position = c(0.01, 0.99),
          plot.margin = margin(c(5.5, 10.5, 5.5, 5.5)))

  return(tc)
}


#' @rdname plot_tail_cal
#' @export
plot_tc_occ <- function(cal, names = NULL, xlab = "t", ylab = "Occurrence ratio", xlims = NULL, ylims = NULL, title = NULL) {

  if (is.null(xlab)) xlab <- "t"
  if (is.null(ylab)) ylab <- "Occurrence ratio"

  if (is.data.frame(cal)) {
    tc <- ggplot(cal) + geom_line(aes(x = t, y = rat))
  } else {
    if (!is.null(names)) names(cal) <- names
    df <- do.call(rbind, cal)
    df$mth <- rep(names(cal), sapply(cal, nrow))
    tc <- ggplot(cal) + geom_line(aes(x = t, y = rat, col = mth))
  }

  tc <- tc + geom_hline(aes(yintercept = 1), linetype = "dotted") +
    scale_x_continuous(name = xlab, limits = xlims, expand = c(0, 0)) +
    scale_y_continuous(name = ylab, limits = ylims, expand = c(0, 0)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.margin = margin(c(5.5, 10.5, 5.5, 5.5))) +
    ggtitle(title)

  return(tc)
}


################################################################################
##### conditional probabilistic tail calibration


#' @rdname plot_tail_cal
#' @export
plot_ptc_div <- function(cal,  names = NULL, xlab = NULL, ylab = NULL, xlims = NULL, ylims = NULL, title = NULL) {

  if (is.data.frame(cal)) {
    tc <- ggplot(cal) + geom_line(aes(x = t, y = d))
  } else {
    if (!is.null(names)) names(cal) <- names
    df <- do.call(rbind, cal)
    df$g <- rep(names(cal), sapply(cal, nrow))
    tc <- ggplot(df) + geom_line(aes(x = t, y = d, col = g))
  }

  tc <- tc + geom_hline(aes(yintercept = 0), lty = "dotted") +
    scale_x_continuous(name = xlab, limits = xlims, expand = c(0, 0)) +
    scale_y_continuous(name = ylab, limits = ylims, expand = c(0, 0)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.justification = c(0, 1),
          legend.position = c(0.01, 0.99),
          plot.margin = margin(c(5.5, 10.5, 5.5, 5.5))) +
    ggtitle(title)

  return(tc)
}


################################################################################
##### marginal tail calibration


#' @rdname plot_tail_cal
#' @export
plot_mtc <- function(cal, names = NULL, ylims = c(-0.2, 0.2), xlims = NULL, ylab = NULL, xlab = NULL, title = NULL) {

  if (!is.data.frame(cal)) {
    if (!is.null(names)) names(cal) <- names
    resampling <- FALSE
  }

  if (ratio == "comb") {
    tc <- plot_mtc_com(cal, resampling = resampling, title = title)
  } else if (ratio == "sev") {
    tc <- plot_mtc_sev(cal, resampling = resampling, title = title)
  } else if (ratio == "occ") {
    tc <- plot_tc_occ(cal, resampling = resampling, title = title)
  }

  return(tc)
}


#' @rdname plot_tail_cal
#' @export
plot_mtc_com <- function(cal, xlab = NULL, ylab = NULL, xlims = NULL, ylims = NULL, title = NULL) {


  return(tc)
}


#' @rdname plot_tail_cal
#' @export
plot_mtc_sev <- function(cal, xlab = NULL, ylab = NULL, xlims = NULL, ylims = NULL, title = NULL) {

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


################################################################################
##### conditional marginal tail calibration



















