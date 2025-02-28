################################################################################
# normal simulation study

################################################################################
### set up

library(TailCalibration)
library(ggplot2)


set.seed(631)

# wrapper to save plots
save_plots <- function(prefix, mth = c("cl", "id", "uf", "sr"), width = 3.5, height = 3.4) {
  for (m in mth) {
    plotname <- paste0(prefix, "_", m)
    filename <- paste0("plots/sim_norm_", prefix, "_1e6_", m, ".png")
    ggsave(filename, plot = get(plotname), width = width, height = height)
  }
}


################################################################################
### simulate data

N <- 1e6

mu <- rnorm(N)
y <- rnorm(N, mu)

tau <- sample(c(-1, 1), length(mu), replace = T)
F_uf <- function(x, m, ta) 0.5*pnorm(x, m) + 0.5*pnorm(x, m + ta)

rd_q <- c(0.9, 0.95, 0.99, 0.999)
rd_vec <- c(-Inf, quantile(y, rd_q))
names <- c(0, rd_q)

t_vec <- quantile(y, c(seq(0, 0.99, 0.01), 0.999))


################################################################################
### combined ratio

com_cl <- tc_prob(y, pnorm, t = rd_vec, sd = sqrt(2))
com_id <- tc_prob(y, pnorm, t = rd_vec, mean = mu)
com_uf <- tc_prob(y, F_uf, t = rd_vec, m = mu, ta = tau)
com_sr <- tc_prob(y, pnorm, t = rd_vec, mean = -mu)

com_cl <- com_cl |> plot_ptc(ratio = "com", names = names, ylims = c(0, 1.02), title = "Climatological")
com_id <- com_id |> plot_ptc(ratio = "com", names = names, ylims = c(0, 1.02), title = "Ideal")
com_uf <- com_uf |> plot_ptc(ratio = "com", names = names, ylims = c(0, 1.02), title = "Unfocused")
com_sr <- com_sr |> plot_ptc(ratio = "com", names = names, ylims = c(0, 1.02), title = "Sign-reversed")

save_plots("com")


################################################################################
### severity ratio

sev_cl <- tc_prob(y, pnorm, t = rd_vec, ratio = "sev", sd = sqrt(2))
sev_id <- tc_prob(y, pnorm, t = rd_vec, ratio = "sev", mean = mu)
sev_uf <- tc_prob(y, F_uf, t = rd_vec, ratio = "sev", m = mu, ta = tau)
sev_sr <- tc_prob(y, pnorm, t = rd_vec, ratio = "sev", mean = -mu)

sev_cl <- sev_cl |> plot_ptc(ratio = "sev", names = names, ylims = c(0, 1.02), title = "")
sev_id <- sev_id |> plot_ptc(ratio = "sev", names = names, ylims = c(0, 1.02), title = "")
sev_uf <- sev_uf |> plot_ptc(ratio = "sev", names = names, ylims = c(0, 1.02), title = "")
sev_sr <- sev_sr |> plot_ptc(ratio = "sev", names = names, ylims = c(0, 1.02), title = "")

save_plots("sev")


################################################################################
### occurrence ratio

occ_cl <- tc_prob(y, pnorm, t = t_vec, ratio = "occ", qu = T, sd = sqrt(2))
occ_id <- tc_prob(y, pnorm, t = t_vec, ratio = "occ", qu = T, mean = mu)
occ_uf <- tc_prob(y, F_uf, t = t_vec, ratio = "occ", qu = T, m = mu, ta = tau)
occ_sr <- tc_prob(y, pnorm, t = t_vec, ratio = "occ", qu = T, mean = -mu)

occ_cl <- occ_cl |> plot_ptc(ratio = "occ", ylims = c(0, 2), title = "")
occ_id <- occ_id |> plot_ptc(ratio = "occ", ylims = c(0, 2), title = "")
occ_uf <- occ_uf |> plot_ptc(ratio = "occ", ylims = c(0, 2), title = "")
occ_sr <- occ_sr |> plot_ptc(ratio = "occ", ylims = c(0, 2), title = "")

save_plots("occ")


