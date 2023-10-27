################################################################################
### set up

#devtools::install_github("sallen12/WeightedForecastVerification")
library(TailCalibration)
library(evmix)

set.seed(89)


################################################################################
### simulate data

gamma <- 1/4
N <- 1e5
v <- 1.4

delta <- rgamma(N, shape = 1/gamma, rate = 1/gamma)
y <- rexp(N, rate = delta)

rd_q <- c(0, 0.9, 0.95, 0.99, 0.999)
rd_vec <- quantile(y, rd_q)

t_vec <- seq(0, 100, 0.1)


################################################################################
### marginal tail calibration

mtc_cl <- tail_marg_cal(y, pgpd, t = rd_vec, z = t_vec, xi = gamma)
mtc_id <- tail_marg_cal(y, pexp, t = rd_vec, z = t_vec, rate = delta)
mtc_ex <- tail_marg_cal(y, pexp, t = rd_vec, z = t_vec, rate = delta/v)


## marginal difference plots

plot_mtc(mtc_id, names = rd_q)
#ggsave("plots/mtc_dif_1e6_id.png", width = 3.2, height = 3)


################################################################################
### probabilistic tail calibration

ptc_cl <- tail_prob_cal(y, pgpd, t = rd_vec, xi = gamma)
ptc_id <- tail_prob_cal(y, pexp, t = rd_vec, rate = delta)
ptc_ex <- tail_prob_cal(y, pexp, t = rd_vec, rate = delta/v)


## cpit reliability diagrams

plot_ptc(ptc_id, names = rd_q)
#ggsave("plots/ptc_rh_1e6_id.png", width = 3.2, height = 3)


################################################################################
### conditional probabilistic tail calibration

ptc_cl <- tail_prob_cal(y, pgpd, t = t_vec, xi = gamma)
ptc_id <- tail_prob_cal(y, pexp, t = t_vec, rate = delta)
ptc_ex <- tail_prob_cal(y, pexp, t = t_vec, rate = delta/v)


## cpit divergence plot

cond_plots <- plot_ptc_cond(list(cl = ptc_cl, id = ptc_id, ex = ptc_ex),
                            names = c("Clim.", "Ideal", "Extr."),
                            t = t_vec, y = y, grp = delta, n_grp = 4, ylim = c(-0.01, 0.1))
#ggsave(plot = cond_plots, "plots/ptc_div_1e6_cond.png", width = 10, height = 2.8)




