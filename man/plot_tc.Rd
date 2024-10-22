% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_tc.R
\name{plot_tc}
\alias{plot_tc}
\alias{plot_tc_occ}
\title{Tail calibration diagnostic plots}
\usage{
plot_tc(
  cal,
  type = c("prob", "marg"),
  ratio = c("com", "sev", "occ"),
  names = NULL,
  xlab = NULL,
  ylab = NULL,
  xlims = NULL,
  ylims = NULL,
  title = NULL
)

plot_tc_occ(
  cal,
  names = NULL,
  xlab = NULL,
  ylab = NULL,
  xlims = NULL,
  ylims = NULL,
  title = NULL
)
}
\arguments{
\item{cal}{data frame or list of data frames containing the data to be plotted.}

\item{type}{string denoting whether to assess marginal ('marg') or probabilistic ('prob') tail calibration.}

\item{ratio}{ratio to plot; one of 'com', 'sev', 'occ'.}

\item{names}{names to be used in legend, if \code{cal} is a list.}

\item{xlab, ylab}{axes labels.}

\item{xlims, ylims}{lower and upper limits of the axes.}

\item{title}{plot title.}
}
\value{
A \code{\link{ggplot2}} object.
}
\description{
Diagnostic plots to assess the tail calibration of probabilistic forecasts.
}
\section{Details}{


\code{plot_tc()} is a wrapper for \code{\link{plot_ptc}} if \code{type == 'prob'} and
\code{\link{plot_mtc}} if \code{type == 'marg'}. Further details about tail calibration
can be found in the help pages for \code{\link{tail_cal}}, \code{\link{tc_prob}}, and
\code{\link{tc_marg}}.



\code{ratio} denotes the ratio of interest to be displayed. This must be one of
\code{'com'}, \code{'sev'}, and \code{'occ'}, which correspond to the combined,
severity, and occurrence ratios, respectively; see \code{\link{tail_cal}} for details.

\code{cal} is a data frame containing the data to be plotted. For \code{ratio == 'com'}
or \code{ratio == 'sev'}, \code{cal} should contain two columns: one labelled \code{u},
containing the vector of values at which calibration has been assessed, and the
other labelled \code{rat}, containing the corresponding ratio measuring miscalibration.
If \code{ratio == 'occ'}, then rather than a column labelled \code{u}, one column should
be labelled \code{t}, containing the thresholds at which calibration has been assessed.
These data frames align with the output of \code{\link{tail_cal}}.
The output of \code{plot_tc()} is a \code{\link{ggplot2}} object displaying the vector \code{u}
or \code{t} against the corresponding measure of miscalibration in \code{rat}.

\code{cal} can also be a list of such data frames, corresponding to different thresholds
or different methods. In this case, the different data frames are displayed
on one plot, with a different colour representing each element in the list. The labels
for each method will be taken to be the names of the list, unless the argument
\code{names} is specified. \code{names} should be a character string of the same length
as the list \code{cal}, with the first element of \code{names} the label of the first
data frame in \code{cal}, and so on.

\code{xlims} and \code{ylims} are both numeric vectors of length two, specifying the
limits of the plots, while \code{xlab}, \code{ylab}, and \code{title} are character
strings specifying the label of the x-axis, y-axis, and title of the plot.
}

\section{References}{

\emph{Probabilistic calibration:}

Diebold, F., Gunther, T.A., and Tay, A.S. (1998):
`Evaluating Density Forecasts with Applications to Financial Risk Management'
\emph{International Economic Review}, 39, 863-883.
DOI: 10.2307/2527342


Gneiting, T., Balabdaoui, F., and Raftery, A. E. (2007):
`Probabilistic forecasts, calibration and sharpness'
\emph{Journal of the Royal Statistical Society Series B: Statistical Methodology}, 69, 243-268.
DOI: 10.1111/j.1467-9868.2007.00587.x


Gneiting, T. and Resin, J. (2022):
`Regression diagnostics meets forecast evaluation: Conditional calibration, reliability diagrams, and coefficient of determination',
\emph{Electronic Journal of Statistics}, 17, 3226-3286.
DOI: 10.1214/23-EJS2180


\emph{Probabilistic tail calibration:}

Allen, S., Koh, J., Segers, J. and Ziegel, J. (2024+):
`Tail calibration of probabilistic forecasts'.
\emph{arXiv pre-print}.


\emph{Marginal calibration:}

Gneiting, T., Balabdaoui, F., and Raftery, A. E. (2007):
`Probabilistic forecasts, calibration and sharpness'
\emph{Journal of the Royal Statistical Society Series B: Statistical Methodology}, 69, 243-268.
DOI: 10.1111/j.1467-9868.2007.00587.x


Gneiting, T. and Resin, J. (2022):
`Regression diagnostics meets forecast evaluation: Conditional calibration, reliability diagrams, and coefficient of determination',
\emph{Electronic Journal of Statistics}, 17, 3226-3286.
DOI: 10.1214/23-EJS2180


\emph{Marginal tail calibration:}

Allen, S., Koh, J., Segers, J. and Ziegel, J. (2024+):
`Tail calibration of probabilistic forecasts'.
\emph{arXiv pre-print}.
}

\examples{
n <- 1e5
y <- rnorm(n)
F_x <- pnorm
mu <- rnorm(n, sd = 0.5)
t <- c(-Inf, 0, 1, 2)

# probabilistic tail calibration
com <- tail_cal(y, F_x, t = t, mean = mu)
plot_tc(com[[1]])
plot_tc(com, names = c('t1', 't2', 't3', 't4'))

sev <- tail_cal(y, F_x, t = t, ratio = 'sev', mean = mu)
plot_tc(sev, ratio = 'sev')

occ <- tail_cal(y, F_x, t = t, ratio = 'occ', mean = mu)
plot_tc(occ, ratio = 'occ', ylims = c(0, 2))

# marginal tail calibration
com <- tail_cal(y, F_x, t = t, type = 'marg', mean = mu)
plot_tc(com[[1]], type = 'marg')
plot_tc(com, type = 'marg', names = c('t1', 't2', 't3', 't4'))

sev <- tail_cal(y, F_x, t = t, ratio = 'sev', mean = mu)
plot_tc(sev, type = 'marg', ratio = 'sev')

occ <- tail_cal(y, F_x, t = t, ratio = 'occ', mean = mu)
plot_tc(occ, type = 'marg', ratio = 'occ', ylims = c(0, 2))


}
\seealso{
\code{\link{tail_cal}} \code{\link{tc_prob}} \code{\link{tc_marg}}
}
\author{
Sam Allen
}
