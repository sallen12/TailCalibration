#' Marginal tail calibration diagnostic plots
#'
#' Diagnostic plots to assess the marginal tail calibration of probabilistic forecasts.
#'
#' @inheritParams plot_tc
#'
#' @inheritSection plot_ptc Details
#'
#'
#' @seealso \code{\link{tc_marg}}
#'
#'
#' @return
#' A \code{\link{ggplot2}} object.
#'
#'
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
#' com <- tail_cal(y, F_x, t = t, type = 'marg', mean = mu)
#' plot_mtc(com[[1]])
#' plot_mtc_com(com[[1]])
#' plot_mtc(com, names = c('t1', 't2', 't3', 't4'))
#'
#' sev <- tail_cal(y, F_x, t = t, ratio = 'sev', mean = mu)
#' plot_mtc(sev, ratio = 'sev')
#' plot_mtc_sev(sev)
#'
#' occ <- tail_cal(y, F_x, t = t, ratio = 'occ', mean = mu)
#' plot_mtc(occ, ratio = 'occ', ylims = c(0, 2))
#' plot_tc_occ(occ)
#'
#'
#' @import ggplot2
#' @name plot_mtc
NULL


#' @rdname plot_mtc
#' @export
plot_mtc <- function(cal, names = NULL, xlab = NULL, ylab = NULL,
                     xlims = NULL, ylims = NULL, title = NULL) {
  ratio <- match.arg(ratio)

  if (!is.data.frame(cal) && !is.null(names)) names(cal) <- names

  if (ratio == "com") {
    tc <- plot_mtc_com(cal, names, xlab, ylab, xlims, ylims, title)
  } else if (ratio == "sev") {
    tc <- plot_mtc_sev(cal, names, xlab, ylab, xlims, ylims, title)
  } else if (ratio == "occ") {
    tc <- plot_tc_occ(cal, names, xlab, ylab, xlims, ylims, title)
  }

  return(tc)
}


#' @rdname plot_mtc
#' @export
plot_mtc_com <- function(cal, names = NULL, xlab = NULL, ylab = NULL,
                         xlims = NULL, ylims = NULL, title = NULL) {

  if (is.null(ylab)) ylab <- "Combined ratio"

  tc <- plot_mtc_sev(cal, names, xlab, ylab, xlims, ylims, title)

  return(tc)
}


#' @rdname plot_mtc
#' @export
plot_mtc_sev <- function(cal, names = NULL, xlab = NULL, ylab = NULL,
                         xlims = NULL, ylims = NULL, title = NULL) {

  if (is.null(xlab)) xlab <- "x"
  if (is.null(ylab)) ylab <- "Severity ratio"
  if (is.null(ylims)) ylims <- c(-1, 1)

  if (is.data.frame(cal)) {
    tc <- ggplot(cal) + geom_line(aes(x = x, y = d))
  } else {
    df <- do.call(rbind, cal)
    if (!is.null(names)) {
      df$mth <- rep(names, sapply(cal, nrow))
    } else {
      df$mth <- rep(names(cal), sapply(cal, nrow))
    }
    tc <- ggplot(df) + geom_line(aes(x = x, y = d, col = mth))
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

