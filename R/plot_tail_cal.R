#' Tail calibration of probabilistic forecasts
#'
#' Plot marginal and probabilistic tail calibration of probabilistic forecasts
#'
#' @param cal data frame or list of data frames containing the data to be plotted.
#' @param t threshold(s) corresponding to the data in \code{cal}.
#' @param xlims,ylims lower and upper limits of the axes.
#' @param xlab,ylab axes labels
#' @param title plot title.
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
#' @name plot_tail_cal
NULL


#' @rdname plot_tail_cal
#' @export
plot_mtc <- function(cal, t = NULL, ylims = c(-0.2, 0.2), xlims = NULL, ylab = NULL, xlab = NULL, title = NULL) {
  if (is.data.frame(cal)) {
    df <- cal
    df$t <- NULL
    mtc <- ggplot(df) + geom_line(aes(x = z, y = d))
  } else {
    df <- do.call(rbind, cal)
    if (is.null(t)) {
      df$t <- rep(names(cal), each = nrow(cal[[1]]))
    } else {
      df$t <- rep(t, each = nrow(cal[[1]]))
    }
    mtc <- ggplot(df) + geom_line(aes(x = z, y = d, col = as.factor(t)))
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


