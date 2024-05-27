#' Tail calibration diagnostic plots
#'
#' Diagnostic plots to assess the tail calibration of probabilistic forecasts.
#'
#' @param cal data frame or list of data frames containing the data to be plotted.
#' @param type string denoting whether to assess marginal ('marg') or probabilistic ('prob') tail calibration.
#' @param ratio ratio to plot; one of 'com', 'sev', 'occ'.
#' @param names names to be used in legend, if \code{cal} is a list.
#' @param xlims,ylims lower and upper limits of the axes.
#' @param xlab,ylab axes labels.
#' @param title plot title.
#'
#'
#' @section Details:
#'
#' \code{plot_tc()} is a wrapper for \code{\link{plot_ptc}} if \code{type == 'prob'} and
#' \code{\link{plot_mtc}} if \code{type == 'marg'}. Further details about tail calibration
#' can be found in the help pages for \code{\link{tail_cal}}, \code{\link{tc_prob}}, and
#' \code{\link{tc_marg}}.
#'
#' @inheritSection plot_ptc Details
#'
#'
#' @seealso \code{\link{tail_cal}} \code{\link{tc_prob}} \code{\link{tc_marg}}
#'
#'
#' @return
#' A \code{\link{ggplot2}} object.
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
#' t <- c(-Inf, 0, 1, 2)
#'
#' # probabilistic tail calibration
#' com <- tail_cal(y, F_x, t = t, mean = mu)
#' plot_tc(com[[1]])
#' plot_tc(com, names = c('t1', 't2', 't3', 't4'))
#'
#' sev <- tail_cal(y, F_x, t = t, ratio = 'sev', mean = mu)
#' plot_tc(sev, ratio = 'sev')
#'
#' occ <- tail_cal(y, F_x, t = t, ratio = 'occ', mean = mu)
#' plot_tc(occ, ratio = 'occ', ylims = c(0, 2))
#'
#' # marginal tail calibration
#' com <- tail_cal(y, F_x, t = t, type = 'marg', mean = mu)
#' plot_tc(com[[1]], type = 'marg')
#' plot_tc(com, type = 'marg', names = c('t1', 't2', 't3', 't4'))
#'
#' sev <- tail_cal(y, F_x, t = t, ratio = 'sev', mean = mu)
#' plot_tc(sev, type = 'marg', ratio = 'sev')
#'
#' occ <- tail_cal(y, F_x, t = t, ratio = 'occ', mean = mu)
#' plot_tc(occ, type = 'marg', ratio = 'occ', ylims = c(0, 2))
#'
#'
#' @import ggplot2
#' @name plot_tc
NULL


#' @rdname plot_tc
#' @export
plot_tc <- function(cal, type = c('prob', 'marg'), ratio = c('com', 'sev', 'occ'),
                    names = NULL, xlab = NULL, ylab = NULL,
                    xlims = NULL, ylims = NULL, title = NULL) {
  type <- match.arg(type)
  ratio <- match.arg(ratio)
  if (type == 'prob') {
    tc <- plot_ptc(cal, ratio, names, xlab, ylab, xlims, ylims, title)
  } else if (type == 'marg') {
    tc <- plot_mtc(cal, ratio, names, xlab, ylab, xlims, ylims, title)
  }
  return(tc)
}


#' @rdname plot_tc
#' @export
plot_tc_occ <- function(cal, names = NULL, xlab = NULL, ylab = NULL,
                        xlims = NULL, ylims = NULL, title = NULL) {

  if (is.null(xlab)) xlab <- "t"
  if (is.null(ylab)) ylab <- "Occurrence ratio"

  if (is.data.frame(cal)) {
    tc <- ggplot(cal) + geom_line(aes(x = t, y = rat))
  } else {
    df <- do.call(rbind, cal)
    if (!is.null(names)) {
      df$mth <- rep(names, sapply(cal, nrow))
    } else {
      df$mth <- rep(names(cal), sapply(cal, nrow))
    }
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


check_tc_inputs <- function(y, F_x, t, u, sup, qu, subset, ...) {
  if (!is.numeric(y) || !is.vector(y)) stop("'y' must be a numeric value or vector")
  if (!is.numeric(t) || !is.vector(t)) stop("'t' must be a numeric value or vector")
  if (!is.numeric(u) || !is.vector(u)) stop("'u' must be a numeric value or vector")
  if (!is.logical(sup) || length(sup) > 1) stop("'sup' must be either TRUE or FALSE")
  if (!is.logical(qu) || length(qu) > 1) stop("'qu' must be either TRUE or FALSE")
  if (!is.logical(subset) || !is.vector(subset) || length(subset) != length(y))
    stop("'subset' must be a logical value or vector of the same length as 'y'")
  if (!is.function(F_x)) stop("'F_x' is not a function")
}

