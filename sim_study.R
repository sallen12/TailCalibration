################################################################################
### set up

library(evmix)
library(ggplot2)
#devtools::install_github("sallen12/WeightedForecastVerification")
library(WeightedForecastVerification)
source("reldiag.R")

set.seed(89)


################################################################################
### simulate data

gamma <- 1/4
N <- 1e5
v <- 1.4

delta <- rgamma(N, shape = 1/gamma, rate = 1/gamma)
y <- rexp(N, rate = delta)

t_vec <- seq(0, 100, 0.1)

plot_calibration <- function(df, type = "", yref = 1, title = "", ylims = NULL, both = TRUE) {
  t_plot <- ggplot(df) + geom_line(aes(x = t, y = mtc, col = mth)) +
    geom_hline(aes(yintercept = yref), lty = "dotted") +
    geom_vline(aes(xintercept = max(y)), lty = "dotted") +
    scale_x_continuous(name = "Threshold", expand = c(0, 0)) +
    scale_y_continuous(name = type, limits = ylims) +
    theme_bw() + theme(panel.grid = element_blank(),
                       legend.title = element_blank(),
                       legend.position = "bottom") +
    ggtitle(title)
  q_plot <- ggplot(df) + geom_line(aes(x = q, y = mtc, col = mth)) +
    geom_hline(aes(yintercept = yref), lty = "dotted") +
    scale_x_continuous(name = "Threshold", expand = c(0, 0)) +
    scale_y_continuous(name = type, limits = ylims) +
    theme_bw() + theme(panel.grid = element_blank(),
                       legend.title = element_blank(),
                       legend.position = "bottom")
  if (both) {
    q_plot <- q_plot + ggtitle(title)
    gridExtra::grid.arrange(t_plot, q_plot)
  } else {
    q_plot <- q_plot + ggtitle("")
    q_plot
  }

}


################################################################################
### marginal tail calibration

rd_q <- c(0, 0.9, 0.95, 0.99, 0.999)
rd_vec <- quantile(y, rd_q)

mtc_cl <- tail_marg_cal(y, pgpd, t = rd_vec, z = t_vec, xi = gamma)
mtc_id <- tail_marg_cal(y, pexp, t = rd_vec, z = t_vec, rate = delta)
mtc_ex <- tail_marg_cal(y, pexp, t = rd_vec, z = t_vec, rate = delta/v)


## marginal difference plots

plot_mtc(mtc_id)
#ggsave("plots/mtc_dif_1e6_id.png", width = 3.2, height = 3)


################################################################################
### probabilistic tail calibration

ptc_cl <- tail_prob_cal(y, pgpd, t = rd_vec, xi = gamma)
ptc_id <- tail_prob_cal(y, pexp, t = rd_vec, rate = delta)
ptc_ex <- tail_prob_cal(y, pexp, t = rd_vec, rate = delta/v)


## cpit reliability diagrams

plot_ptc(ptc_ex, names = rd_q)
#ggsave("plots/ptc_rh_1e6_cl.png", width = 3.2, height = 3)


################################################################################
### conditional probabilistic tail calibration
Rcpp::sourceCpp("crps_div.cpp")

ptc_clim_del <- lapply(seq_along(t_vec), function(i) cbind(ptc_clim[[i]], delta[y > t_vec[i]]))
ptc_id_del <- lapply(seq_along(t_vec), function(i) cbind(ptc_id[[i]], delta[y > t_vec[i]]))
ptc_ex_del <- lapply(seq_along(t_vec), function(i) cbind(ptc_ex[[i]], delta[y > t_vec[i]]))

del_n <- 5
del_quants <- seq(0, 1, length.out = del_n + 1)

del_plots <- list()
for (del_i in 1:del_n) {

  # get category limits
  del_t <- quantile(delta, del_quants[c(del_i, del_i + 1)])
  if (del_i == 1) del_t[1] <- 0

  # plot calibration
  df <- data.frame(t = t_vec,
                   q = sapply(t_vec, function(t) mean(y[delta > del_t[1] & delta <= del_t[2]] <= t)),
                   mtc = c(sapply(ptc_clim_del, function(x)
                     tryCatch(crps_div(x[x[, 2] > del_t[1] & x[, 2] <= del_t[2], 1]), error = function(x) NA)),
                     sapply(ptc_id_del, function(x)
                       tryCatch(crps_div(x[x[, 2] > del_t[1] & x[, 2] <= del_t[2], 1]), error = function(x) NA)),
                     sapply(ptc_ex_del, function(x)
                       tryCatch(crps_div(x[x[, 2] > del_t[1] & x[, 2] <= del_t[2], 1]), error = function(x) NA))),
                   mth = rep(c("Clim", "Ideal", "Extremist"), each = length(t_vec)))

  del_plots[[del_i]] <-
    plot_calibration(df, type = "Probabilistic tail calibration", yref = 0, both = F, ylims = c(-0.02, 0.1),
                     title = paste0("Min: ", round(del_t[1], 2), ", Max: ", round(del_t[2], 2)))
}

cond_plots <- do.call(gridExtra::grid.arrange, c(del_plots, nrow = 1))
ggsave(plot = cond_plots, "plots/ptc_div_1e6_cond.png", width = 17, height = 3.9)










