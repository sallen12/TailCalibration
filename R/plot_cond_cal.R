#' Conditional calibration of probabilistic forecasts
#'
#' Plot conditional probabilistic tail calibration of probabilistic forecasts
#'
#' @param cal data frame or list of data frames containing the data to be plotted.
#' @param t threshold(s) corresponding to the data in \code{cal}.
#' @param y vector of observations.
#' @param grp vector of values that are used to perform the conditioning.
#' @param n_grp number of groups to condition on.
#' @param names names to be used in legend.
#' @param xlims,ylims lower and upper limits of the axes.
#' @param xlab,ylab axes labels
#' @param title plot title.
#' @param quantile_sc logical specifying whether the x-axis of the plot should be
#'  on the quantile scale or threshold scale.
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
#' @name plot_cond_cal
NULL

#' @rdname plot_cond_cal
#' @export
plot_ptc_cond <- function(cal, t, y, grp, n_grp = 1, names = NULL, ylims = NULL, xlims = NULL,
                          ylab = "Miscalibration", xlab = "Threshold", title = NULL, quantile_sc = T) {

  quants <- seq(0, 1, length.out = n_grp + 1)

  cond_plots <- list()
  for (i_g in 1:n_grp) {

    # get category limits
    grp_ab <- quantile(grp, quants[c(i_g, i_g + 1)])

    # get calibration
    cal_grp <- lapply(cal, function(x) lapply(seq_along(t), function(i) {
      ind <- (grp[y > t[i]] > grp_ab[1]) & (grp[y > t[i]] <= grp_ab[2])
      return(x[[i]][ind, ])
    }))

    # convert thresholds to quantile scale
    if (quantile_sc) {
      q_vec <- sapply(t, function(tt) mean(y[(grp > grp_ab[1]) & (grp <= grp_ab[2])] <= tt))
    } else {
      q_vec <- t
    }

    # plot calibration
    cond_plots[[i_g]] <- plot_ptc_div(cal_grp, q_vec, names, ylims, xlims, ylab, xlab, title)
  }

  ptc_div <- do.call(gridExtra::grid.arrange, c(cond_plots, nrow = 1))

  return(ptc_div)
}
