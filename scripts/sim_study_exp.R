################################################################################
# exponential simulation study


################################################################################
### set up

library(TailCalibration)
library(evd)
library(ggplot2)


set.seed(89)

# wrapper to save plots
save_plots <- function(prefix, mth = c("cl", "id", "ex"), width = 3.5, height = 3.4) {
  for (m in mth) {
    plotname <- paste0(prefix, "_", m)
    filename <- paste0("plots/sim_ex_", prefix, "_1e6_", m, ".png")
    ggsave(filename, plot = get(plotname), width = width, height = height)
  }
}


################################################################################
### simulate data

gamma <- 1/4
N <- 1e6
v <- 1.4

delta <- rgamma(N, shape = 1/gamma, rate = 1/gamma)
y <- rexp(N, rate = delta)

rd_q <- c(0.9, 0.95, 0.99, 0.999)
rd_vec <- c(0, quantile(y, rd_q))
names <- c(0, rd_q)

t_vec <- quantile(y, c(seq(0, 0.99, 0.01), 0.999))


################################################################################
### combined ratio

com_cl <- tc_prob(y, pgpd, t = rd_vec, shape = gamma)
com_id <- tc_prob(y, pexp, t = rd_vec, rate = delta)
com_ex <- tc_prob(y, pexp, t = rd_vec, rate = delta/v)

com_cl <- com_cl |> plot_ptc(ratio = "com", names = names, ylims = c(0, 1.02), title = "Climatological")
com_id <- com_id |> plot_ptc(ratio = "com", names = names, ylims = c(0, 1.02), title = "Ideal")
com_ex <- com_ex |> plot_ptc(ratio = "com", names = names, ylims = c(0, 1.02), title = "Extremist")

save_plots("com")


################################################################################
### severity ratio

sev_cl <- tc_prob(y, pgpd, t = rd_vec, ratio = "sev", shape = gamma)
sev_id <- tc_prob(y, pexp, t = rd_vec, ratio = "sev", rate = delta)
sev_ex <- tc_prob(y, pexp, t = rd_vec, ratio = "sev", rate = delta/v)

sev_cl <- sev_cl |> plot_ptc(ratio = "sev", names = names, ylims = c(0, 1.02), title = "")
sev_id <- sev_id |> plot_ptc(ratio = "sev", names = names, ylims = c(0, 1.02), title = "")
sev_ex <- sev_ex |> plot_ptc(ratio = "sev", names = names, ylims = c(0, 1.02), title = "")

save_plots("sev")


################################################################################
### occurrence ratio

occ_cl <- tc_prob(y, pgpd, t = t_vec, ratio = "occ", qu = T, shape = gamma)
occ_id <- tc_prob(y, pexp, t = t_vec, ratio = "occ", qu = T, rate = delta)
occ_ex <- tc_prob(y, pexp, t = t_vec, ratio = "occ", qu = T, rate = delta/v)

occ_cl <- occ_cl |> plot_ptc(ratio = "occ", ylims = c(0, 2), title = "")
occ_id <- occ_id |> plot_ptc(ratio = "occ", ylims = c(0, 2), title = "")
occ_ex <- occ_ex |> plot_ptc(ratio = "occ", ylims = c(0, 2), title = "")

save_plots("occ")


################################################################################
### conditional combined ratio

n_grp <- 3
group <- numeric(length(delta))
for (i in 1:n_grp) group[delta >= quantile(delta, (i - 1)/n_grp)] <- paste0("B", i)

com_cl <- tc_cprob(y, pgpd, t = t_vec, group = group, qu = T, shape = gamma)
com_id <- tc_cprob(y, pexp, t = t_vec, group = group, qu = T, rate = delta)
com_ex <- tc_cprob(y, pexp, t = t_vec, group = group, qu = T, rate = delta/v)

com_b1 <- list("Clim." = com_cl[["B1"]], "Ideal" = com_id[["B1"]], "Extr." = com_ex[["B1"]])
com_b2 <- list("Clim." = com_cl[["B2"]], "Ideal" = com_id[["B2"]], "Extr." = com_ex[["B2"]])
com_b3 <- list("Clim." = com_cl[["B3"]], "Ideal" = com_id[["B3"]], "Extr." = com_ex[["B3"]])

com_div_b1 <- com_b1 |> plot_tc_sup(ylab = NULL, ylims = c(-0.1, 2.4), title = expression(B[1]))
com_div_b2 <- com_b2 |> plot_tc_sup(ylab = NULL, ylims = c(-0.1, 2.4), title = expression(B[2]))
com_div_b3 <- com_b3 |> plot_tc_sup(ylab = NULL, ylims = c(-0.1, 2.4), title = expression(B[3]))

save_plots("com_div", mth = c("b1", "b2", "b3"), width = 3, height = 3)

