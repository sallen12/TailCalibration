% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tail_cal.R
\name{tail_cal}
\alias{tail_cal}
\title{Tail calibration of probabilistic forecasts}
\usage{
tail_cal(
  y,
  F_x,
  t,
  type = c("prob", "marg"),
  ratio = c("com", "sev", "occ"),
  u = seq(0.01, 0.99, 0.01),
  lower = -Inf,
  group = NULL,
  sup = FALSE,
  qu = FALSE,
  subset = rep(TRUE, length(y)),
  ...
)
}
\arguments{
\item{y}{vector of observations.}

\item{F_x}{forecast distribution function to be evaluated, or a vector or matrix of samples from the forecast distribution.}

\item{t}{vector of threshold(s) at which to evaluate tail calibration.}

\item{type}{string denoting whether to assess marginal ('marg') or probabilistic ('prob') tail calibration.}

\item{ratio}{ratio to return; one of 'com', 'sev', 'occ'.}

\item{u}{vector of values at which to assess tail calibration.}

\item{lower}{numeric value at which the forecast distributions are censored.}

\item{group}{vector specifying observations that should be grouped together.}

\item{sup}{logical specifying whether to quantify miscalibration using the supremum distance; default is \code{FALSE}.}

\item{qu}{logical specifying whether the thresholds should be renamed in terms of
quantiles of \code{y}; default is \code{FALSE}.}

\item{subset}{logical vector of the same length as \code{y}, allowing only a subset
of the forecasts and observations (where \code{TRUE}) to be assessed.}

\item{...}{additional arguments to F_x.}
}
\description{
Assess marginal and probabilistic tail calibration of probabilistic forecasts.
This allows forecast calibration to be assessed when interest is on values above
thresholds of interest.
}
\section{Details}{


\code{tail_cal()} is a wrapper for \code{\link{tc_prob}} if \code{type == 'prob'} and
\code{\link{tc_marg}} if \code{type == 'marg'}. If \code{group} is not \code{NULL},
then \code{tail_cal()} is a wrapper for \code{\link{tc_cond}}.

Forecast evaluation is typically performed using a sequence of forecasts \eqn{F_{i}}
and corresponding observations \eqn{y_{i}}, for \eqn{i = 1, \dots, n}. The forecasts
are said to be calibrated if the forecasts align statistically with the observations.
Several notions of forecast calibration exist. Forecast tail calibration allows
the calibration of forecasts to be assessed when interest is on values that exceed
a high threshold.



\emph{Probabilistic tail calibration}

Forecasts \eqn{F_{i}} for observations \eqn{y_{i}} are \emph{probabilistically calibrated}
if the probability integral transform (PIT) values \eqn{z_{i} = F_{i}(y_{i})} resemble a sample
from a standard uniform distribution. That is, for all \eqn{u \in [0, 1]},
\deqn{\frac{1}{n} \sum_{i=1}^{n} 1 \{z_{i} \le u\} \approxeq u}
where \eqn{1\{ \cdot \}} denotes the indicator function.
To assess this, we can plot the term on the left hand side as a function of \eqn{u}.

Forecasts are said to be \emph{probabilistically tail calibrated} if, for all \eqn{u \in [0, 1]},
\deqn{\frac{\sum_{i=1}^{n} 1 \{z_{t,i} \le u, y_{i} > t\}}{\sum_{i=1}^{n}(1 - F_{i}(t))} \approxeq u,}
as the threshold \eqn{t} tends to the upper end point of the distribution of the observations.
To assess this, we can choose increasingly large thresholds and plot the term on the left hand
side as a function of \eqn{u}.

We refer to the term of the left hand side as the \emph{combined ratio}, since it can be
decomposed into two more intuitive conditions:
\deqn{\frac{\sum_{i=1}^{n} 1 \{z_{t,i} \le u, y_{i} > t\}}{\sum_{i=1}^{n} 1 \{y_{i} > t\}} \approxeq u,}
for all \eqn{u \in [0, 1]}, where \eqn{z_{t,i} = F_{t,i}(y_{i} - t)} are conditional
PIT values, with
\deqn{F_{t,i}(x) = \frac{F_{i}(x + t) - F_{i}(t)}{1 - F_{i}(t)}}
for \eqn{x \ge 0} the forecast excess distribution; and
\deqn{\frac{\sum_{i=1}^{n} 1 \{y_{i} > t\}}{\sum_{i=1}^{n}(1 - F_{i}(t))} \approxeq 1.}

The first of these two conditions evaluates the conditional forecast distribution,
assessing whether the conditional PIT values \eqn{z_{t,i}} resemble a sample from
a standard uniform distribution, conditional on the outcomes \eqn{y_{i}} being
larger than threshold \eqn{t}. We therefore refer to the ratio in this condition
as the \emph{severity ratio}, since it only assesses forecasts for the severity of
the extreme event and does not consider the probability assigned to being above
the threshold.

The second of these two conditions evaluates the forecast's ability to predict the
occurrence of a threshold exceedance, by comparing the average forecast probability
that the outcome will exceed the threshold with the unconditional observed probability
of the outcome exceeding the threshold. This ratio should be close to one if the
forecast is calibrated. We refer to the ratio in this condition as the
\emph{occurrence} ratio, since it only assesses the forecast probability for
extreme event occurrence and does not consider by how much the outcome exceeds
the threshold.

The function \code{tc_prob()} takes the observations \eqn{y_{i}} and forecast
distributions \eqn{F_{i}} as inputs (via the arguments \code{y} and \code{F_x}),
and returns one of these three ratios: the combined ratio if \code{ratio == 'com'}
(the default), the severity ratio if \code{ratio == 'sev'}, and the occurrence
ratio if \code{ratio == 'occ'}. The vector \code{u} specifies the values
\eqn{u \in [0, 1]} at which to assess the combined and severity ratios, and
the numeric \code{t} specifies the threshold(s) of interest. \code{t} can either be a
single value or a vector.

For \code{ratio == 'com'} or \code{ratio == 'sev'}, \code{tc_prob()} returns a data frame
containing the vector \code{u} and the chosen ratio evaluated at each element in
this vector. If \code{t} is a vector, a list of data frames is returned, with each
element of the list corresponding to a different threshold. For \code{ratio == 'occ'},
\code{tc_prob()} returns a single value containing the occurrence ratio at this threshold.
If \code{t} is a vector, a data frame is instead returned that contains the thresholds \eqn{t}
and the corresponding occurrence ratios.

If \code{sup = TRUE}, then, for each threshold, the maximum absolute distance between the
ratio and the expected output for a calibrated forecast is returned.
For example, for \code{ratio == 'com'}, \code{tc_prob()} calculates
\deqn{\max_{u} | \frac{\sum_{i=1}^{n} 1 \{z_{t,i} \le u, y_{i} > t\}}{\sum_{i=1}^{n}(1 - F_{i}(t))} - u |,}
with the maximum taken over all elements in the argument \code{u}. This single value is
returned, or, if \code{t} is a vector, a data frame is returned containing the thresholds
\code{t} and the corresponding maximum absolute distance.

If \code{qu = TRUE}, then the thresholds are re-expressed in terms of quantiles of the observation
vector \code{y}. This is only relevant for the presentation of the output.

Additional arguments to \code{F_x} can be specified using \code{...}. \code{tc_prob()}
assumes that the distribution function is the same in all forecast cases, and that
\code{...} contains parameters corresponding to this distribution function.
If this is not the case, then a work around is to define \code{F_x} as a regime-switching
distribution, with an additional parameter specifying which distribution/regime to choose;
see examples.



\emph{Marginal tail calibration}

Forecasts \eqn{F_{i}} for observations \eqn{y_{i}} are \emph{marginally calibrated}
if average forecast probabilities are the same as unconditional observed
probabilities. That is, for all \eqn{u \in (-\infty, \infty)},
\deqn{\frac{1}{n} \sum_{i=1}^{n} 1 \{y_{i} \le u\} \approxeq \frac{1}{n} \sum_{i=1}^{n} F_{i}(u),}
where \eqn{1\{ \cdot \}} denotes the indicator function.
To assess this, we can plot the difference between these two terms as a function of \eqn{u}.

Forecasts are said to be \emph{marginally tail calibrated} if, for all \eqn{u \in [0, \infty)},
\deqn{\frac{\sum_{i=1}^{n} 1 \{t < y_{i} \le u + t\} }{\sum_{i=1}^{n}(1 - F_{i}(t))} \approxeq \frac{1}{n_{t}} \sum_{i=1}^{n} F_{t,i}(u) 1 \{y_{i} > t\},}
as the threshold \eqn{t} tends to the upper end point of the distribution of the observations,
where \eqn{n_{t}} is the number of observations that exceed \eqn{t}, and
\deqn{F_{t,i}(u) = \frac{F_{i}(u + t) - F_{i}(t)}{1 - F_{i}(t)}}
for \eqn{u \ge 0} is the forecast excess distribution.
To assess marginal tail calibration, we can plot the difference between the left and right hand sides as a
function of \eqn{u}.

We refer to the difference between these terms as the \emph{combined ratio}, since it can be
decomposed into two more intuitive conditions:
\deqn{\frac{\sum_{i=1}^{n} 1 \{t < y_{i} \le u + t\} }{\sum_{i=1}^{n} 1 \{y_{i} > t\}} \approxeq \frac{1}{n_{t}} \sum_{i=1}^{n} F_{t,i}(u) 1 \{y_{i} > t\},}
for all \eqn{u \in [0, \infty)}; and
\deqn{\frac{\sum_{i=1}^{n} 1 \{y_{i} > t\}}{\sum_{i=1}^{n}(1 - F_{i}(t))} \approxeq 1.}

The first of these two conditions assesses whether the average forecast excess distribution
is equal to the unconditional excess distribution of the observations, conditional on
the outcomes \eqn{y_{i}} being larger than threshold \eqn{t}. We therefore refer to the
differences between these terms as the \emph{severity ratio}, since it only assesses
forecasts for the severity of the extreme event and does not consider the probability
assigned to being above the threshold.

The second of these two conditions evaluates the forecast's ability to predict the
occurrence of a threshold exceedance, by comparing the average forecast probability
that the outcome will exceed the threshold with the unconditional observed probability
of the outcome exceeding the threshold. This ratio should be close to one if the
forecast is calibrated. We refer to the ratio in this condition as the
\emph{occurrence} ratio, since it only assesses the forecast probability for
extreme event occurrence and does not consider by how much the outcome exceeds
the threshold.

The function \code{tc_marg()} takes the observations \eqn{y_{i}} and forecast
distributions \eqn{F_{i}} as inputs (via the arguments \code{y} and \code{F_x}),
and returns one of these three ratios: the combined ratio if \code{ratio == 'com'}
(the default), the severity ratio if \code{ratio == 'sev'}, and the occurrence
ratio if \code{ratio == 'occ'}. The vector \code{u} specifies the values
\eqn{u \in [0, \infty)} at which to assess the combined and severity ratios, and
the numeric \code{t} specifies the threshold(s) of interest. \code{t} can either be a
single value or a vector.

For \code{ratio == 'com'} or \code{ratio == 'sev'}, \code{tc_marg()} returns a data frame
containing the vector \code{u} and the chosen ratio evaluated at each element in
this vector. If \code{t} is a vector, a list of data frames is returned, with each
element of the list corresponding to a different threshold. For \code{ratio == 'occ'},
\code{tc_marg()} returns a single value containing the occurrence ratio at this threshold.
If \code{t} is a vector, a data frame is instead returned that contains the thresholds \eqn{t}
and the corresponding occurrence ratios.

If \code{sup = TRUE}, then, for each threshold, the maximum absolute distance between the
ratio and the expected output for a calibrated forecast is returned.
For example, for \code{ratio == 'com'}, \code{tc_marg()} calculates
\deqn{\max_{u} | \frac{\sum_{i=1}^{n} 1 \{t < y_{i} \le u + t\} }{\sum_{i=1}^{n}(1 - F_{i}(t))} - \frac{1}{n_{t}} \sum_{i=1}^{n} F_{t,i}(u) 1 \{y_{i} > t\} |,}
with the maximum taken over all elements in the argument \code{u}. This single value is
returned, or, if \code{t} is a vector, a data frame is returned containing the thresholds
\code{t} and the corresponding maximum absolute distance.

If \code{qu = TRUE}, then the thresholds are re-expressed in terms of quantiles of the observation
vector \code{y}. This is only relevant for the presentation of the output.

Additional arguments to \code{F_x} can be specified using \code{...}. \code{tc_marg()}
assumes that the distribution function is the same in all forecast cases, and that
\code{...} contains parameters corresponding to this distribution function.
If this is not the case, then a work around is to define \code{F_x} as a regime-switching
distribution, with an additional parameter specifying which distribution/regime to choose;
see examples.



\code{tc_cond()} is a wrapper for \code{tc_cprob} if \code{type == 'prob'} and
\code{tc_cmarg} if \code{type == 'marg'}, which are themselves wrappers for
\code{\link{tc_prob}} and \code{\link{tc_marg}}, respectively. More details can be
found in their respective help pages.

\code{tc_cprob} and \code{tc_cmarg} stratify the forecasts and observations into
different groups, according to the input vector \code{group}. \code{group} should
be a vector of the same length as \code{y}, containing a relatively small number of
unique entries. Observations \code{y} that have the same value of \code{group} are
grouped together, and the standard notions of probabilistic (\code{type == 'prob'})
or marginal tail calibration (\code{type == 'marg'}) are then assessed using the
functionality in \code{tc_prob} or \code{tc_marg}, respectively.
}

\section{Value}{


Data frame or list of data frames (for each threshold \code{t}) containing the
requested ratio quantifying tail calibration.

If \code{sup = TRUE} or \code{ratio == 'occ'}, a data frame is returned that contains
the thresholds \code{t} along
with the corresponding measure of tail calibration at this threshold.
If \code{t} is a single numeric, then \code{sup = TRUE} returns a single numeric value.
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

# standard calibration is assessed by setting t to -Inf
ptc <- tail_cal(y, F_x, t = -Inf, mean = mu)
mtc <- tail_cal(y, F_x, t = -Inf, type = "marg", mean = mu)

t <- c(-Inf, 0, 1, 2)
com <- tail_cal(y, F_x, t = t, sup = TRUE, mean = mu)
com2 <- tc_prob(y, F_x, t = t, sup = TRUE, mean = mu)
identical(com, com2)

sev <- tail_cal(y, F_x, t = t, type = 'marg', ratio = 'sev', mean = mu)
sev2 <- tc_marg(y, F_x, t = t, ratio = 'sev', mean = mu)
identical(sev, sev2)


# ensemble forecast
M <- 100
F_x <- matrix(rnorm(n*M, mu), nrow = n, ncol = M)
ptc_ens <- tail_cal(y, F_x, t = -Inf, mean = mu)
plot_ptc(ptc)
plot_ptc(ptc_ens)


}
\seealso{
\code{\link{tc_prob}} \code{\link{tc_marg}}
}
\author{
Sam Allen
}
