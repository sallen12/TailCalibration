#' Probabilistic tail calibration diagnostic plots
#'
#' Diagnostic plots to assess the probabilistic tail calibration of probabilistic forecasts.
#'
#' @inheritParams plot_tc
#'
#' @section Details:
#'
#' \code{ratio} denotes the ratio of interest to be displayed. This must be one of
#' \code{'com'}, \code{'sev'}, and \code{'occ'}, which correspond to the combined,
#' severity, and occurrence ratios, respectively; see \code{\link{tail_cal}} for details.
#'
#' \code{cal} is a data frame containing the data to be plotted. For \code{ratio == 'com'}
#' or \code{ratio == 'sev'}, \code{cal} should contain two columns: one labelled \code{u},
#' containing the vector of values at which calibration has been assessed, and the
#' other labelled \code{rat}, containing the corresponding ratio measuring miscalibration.
#' If \code{ratio == 'occ'}, then rather than a column labelled \code{u}, one column should
#' be labelled \code{t}, containing the thresholds at which calibration has been assessed.
#' These data frames align with the output of \code{\link{tail_cal}}.
#' The output of \code{plot_tc()} is a \code{\link{ggplot2}} object displaying the vector \code{u}
#' or \code{t} against the corresponding measure of miscalibration in \code{rat}.
#'
#' \code{cal} can also be a list of such data frames, corresponding to different thresholds
#' or different methods. In this case, the different data frames are displayed
#' on one plot, with a different colour representing each element in the list. The labels
#' for each method will be taken to be the names of the list, unless the argument
#' \code{names} is specified. \code{names} should be a character string of the same length
#' as the list \code{cal}, with the first element of \code{names} the label of the first
#' data frame in \code{cal}, and so on.
#'
#' \code{xlims} and \code{ylims} are both numeric vectors of length two, specifying the
#' limits of the plots, while \code{xlab}, \code{ylab}, and \code{title} are character
#' strings specifying the label of the x-axis, y-axis, and title of the plot.
#'
#'
#' @seealso \code{\link{tc_prob}}
#'
#'
#' @return
#' A \code{\link{ggplot2}} object.
#'
#'
#' @inheritSection tc_prob References
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
#' t <- c(-Inf, 0, 1, 2)
#'
#' com <- tail_cal(y, F_x, t = t, mean = mu)
#' plot_ptc(com[[1]])
#' plot_ptc_com(com[[1]])
#' plot_ptc(com, names = c('t1', 't2', 't3', 't4'))
#'
#' sev <- tail_cal(y, F_x, t = t, ratio = 'sev', mean = mu)
#' plot_ptc(sev, ratio = 'sev')
#' plot_ptc_sev(sev)
#'
#' occ <- tail_cal(y, F_x, t = t, ratio = 'occ', mean = mu)
#' plot_ptc(occ, ratio = 'occ', ylims = c(0, 2))
#' plot_tc_occ(occ, ylims = c(0, 2))
#'
#'
#' @import ggplot2
#' @name plot_ptc
NULL


#' @rdname plot_ptc
#' @export
plot_ptc <- function(cal, ratio = c("com", "sev", "occ"), names = NULL,
                     xlab = NULL, ylab = NULL, xlims = NULL, ylims = NULL, title = NULL) {
  ratio <- match.arg(ratio)

  if (!is.data.frame(cal) && !is.null(names)) names(cal) <- names

  if (ratio == "com") {
    tc <- plot_ptc_com(cal, names, xlab, ylab, xlims, ylims, title)
  } else if (ratio == "sev") {
    tc <- plot_ptc_sev(cal, names, xlab, ylab, xlims, ylims, title)
  } else if (ratio == "occ") {
    tc <- plot_tc_occ(cal, names, xlab, ylab, xlims, ylims, title)
  }

  return(tc)
}


#' @rdname plot_ptc
#' @export
plot_ptc_com <- function(cal, names = NULL, xlab = NULL, ylab = NULL,
                         xlims = NULL, ylims = NULL, title = NULL) {

  if (is.null(ylab)) ylab <- "Combined ratio"

  tc <- plot_ptc_sev(cal, names, xlab, ylab, xlims, ylims, title)

  return(tc)
}


#' @rdname plot_ptc
#' @export
plot_ptc_sev <- function(cal, names = NULL, xlab = NULL, ylab = NULL,
                         xlims = NULL, ylims = NULL, title = NULL) {

  if (is.null(xlab)) xlab <- "u"
  if (is.null(ylab)) ylab <- "Severity ratio"
  if (is.null(xlims)) xlims <- c(0, 1)

  if (is.data.frame(cal)) {
    tc <- ggplot(cal) + geom_line(aes(x = u, y = rat))
  } else {
    df <- do.call(rbind, cal)
    if (!is.null(names)) {
       df$mth <- rep(names, sapply(cal, nrow))
    } else {
      df$mth <- rep(names(cal), sapply(cal, nrow))
    }
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
          plot.margin = margin(c(5.5, 10.5, 5.5, 5.5))) +
    ggtitle(title)

  return(tc)
}

