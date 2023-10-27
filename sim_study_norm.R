################################################################################
### set up

#devtools::install_github("sallen12/WeightedForecastVerification")
library(TailCalibration)

set.seed(631)


################################################################################
### simulate data

N <- 1e6

mu <- rnorm(N)
y <- rnorm(N, mu)

tau <- sample(c(-1, 1), length(mu), replace = T)
F_uf <- function(x, m, ta) 0.5*pnorm(x, m) + 0.5*pnorm(x, m + ta)

rd_q <- c(0, 0.9, 0.95, 0.99, 0.999)
rd_vec <- quantile(y, rd_q)

t_vec <- seq(0, 10, 0.1)


################################################################################
### marginal tail calibration

mtc_cl <- tail_marg_cal(y, pnorm, t = rd_vec, z = t_vec, sd = sqrt(2))
mtc_id <- tail_marg_cal(y, pnorm, t = rd_vec, z = t_vec, mean = mu)
mtc_uf <- tail_marg_cal(y, F_uf, t = rd_vec, z = t_vec, m = mu, ta = tau)
mtc_sr <- tail_marg_cal(y, pnorm, t = rd_vec, z = t_vec, mean = -mu)

## marginal difference plots

plot_mtc(mtc_sr, names = rd_q, ylims = c(-0.2, 0.4))
#ggsave("plots/mtc_dif_1e6_norm_id.png", width = 3.2, height = 3)


################################################################################
### probabilistic tail calibration

ptc_cl <- tail_prob_cal(y, pnorm, t = rd_vec, sd = sqrt(2))
ptc_id <- tail_prob_cal(y, pnorm, t = rd_vec, mean = mu)
ptc_uf <- tail_prob_cal(y, F_uf, t = rd_vec, m = mu, ta = tau)
ptc_sr <- tail_prob_cal(y, pnorm, t = rd_vec, mean = -mu)


## cpit reliability diagrams

plot_ptc(ptc_id, names = rd_q)
#ggsave("plots/ptc_rh_1e6_norm_id.png", width = 3.2, height = 3)

