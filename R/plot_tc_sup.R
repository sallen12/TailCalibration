#' Conditional tail calibration diagnostic plots
#'
#' Diagnostic plots to assess the tail calibration of probabilistic forecasts
#' conditional on non-trivial information sets.
#'
#' @inheritParams plot_tc
#'
#' @details
#'
#' \code{cal} is a data frame containing the data to be plotted. \code{cal} should contain
#' two columns: one labelled \code{t}, containing the thresholds at which calibration has
#' been assessed, and the other labelled \code{rat}, containing the corresponding supremum
#' distance measuring miscalibration at this threshold.
#' These data frames align with the output of \code{\link{tail_cal}} and \code{\link{tc_cond}}.
#' The output of \code{plot_tc_sup()} is a \code{\link{ggplot2}} object displaying the vector
#' \code{t} against the corresponding measure of miscalibration in \code{rat}.
#'
#' \code{cal} can also be a list of such data frames, corresponding to different groupings
#' of the data; if the groupings are constructed to represent an information set, then this
#' allows conditional tail calibration to be assessed with respect to non-trivial information sets.
#' In this case, the different data frames are displayed on one plot, with a different colour
#' representing each element in the list. The labels for each method will be taken to be the
#' names of the list, unless the argument \code{names} is specified. \code{names} should be a
#' character string of the same length as the list \code{cal}, with the first element of
#' \code{names} the label of the first data frame in \code{cal}, and so on.
#'
#' \code{xlims} and \code{ylims} are both numeric vectors of length two, specifying the
#' limits of the plots, while \code{xlab}, \code{ylab}, and \code{title} are character
#' strings specifying the label of the x-axis, y-axis, and title of the plot.
#'
#' @seealso \code{\link{tc_cond}} \code{\link{plot_tc}}
#'
#'
#' @return
#' A \code{\link{ggplot2}} object.
#'
#'
#' @references
#' Allen, S., Koh, J., Segers, J. and Ziegel, J. (2024+):
#' `Tail calibration of probabilistic forecasts'.
#' \emph{arXiv pre-print}.
#'
#'
#' @author Sam Allen
#'
#' @examples
#' n <- 1e5
#' y <- rnorm(n)
#' F_x <- pnorm
#' mu <- rnorm(n, sd = 0.5)
#'
#' group <- sample(c('A', 'B'), length(y), replace = TRUE)
#' t <- c(-Inf, 0, 1, 2)
#' sev <- tc_cond(y, F_x, t = t, ratio = 'sev', group = group, mean = mu)
#'
#' plot_tc_sup(sev, ylims = c(-0.01, 0.05))
#'
#'
#' @import ggplot2
#' @name plot_tc_sup
NULL


#' @rdname plot_tc_sup
#' @export
plot_tc_sup <- function(cal, names = NULL, xlab = NULL, ylab = NULL,
                        xlims = NULL, ylims = NULL, title = NULL) {

  if (is.null(xlab)) xlab <- "t"
  if (is.null(ylab)) ylab <- "Miscalibration"

  if (is.data.frame(cal)) {
    tc <- ggplot(cal) + geom_line(aes(x = t, y = rat))
  } else {
    df <- do.call(rbind, cal)
    if (!is.null(names)) {
      if (is.numeric(names)) names <- as.factor(names)
      df$g <- rep(names, sapply(cal, nrow))
    } else {
      df$g <- rep(names(cal), sapply(cal, nrow))
    }
    tc <- ggplot(df) + geom_line(aes(x = t, y = rat, col = g))
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

