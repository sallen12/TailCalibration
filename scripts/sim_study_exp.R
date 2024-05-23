################################################################################
### set up

library(TailCalibration)
library(evmix)
library(ggplot2)


set.seed(89)

# wrapper to save plots
save_plots <- function(prefix, mth = c("cl", "id", "ex"), width = 3.5, height = 3.4) {
  for (m in mth) {
    plotname <- paste0(prefix, "_", m)
    filename <- paste0("plots/ptc_", prefix, "_1e6_", m, ".png")
    ggsave(filename, plot = get(plotname), width = width, height = height)
  }
}


################################################################################
### simulate data

gamma <- 1/4
N <- 1e5
v <- 1.4

delta <- rgamma(N, shape = 1/gamma, rate = 1/gamma)
y <- rexp(N, rate = delta)

rd_q <- c(0.9, 0.95, 0.99, 0.999)
rd_vec <- c(0, quantile(y, rd_q))
names <- c(0, rd_q)

t_vec <- quantile(y, c(seq(0, 0.99, 0.01), 0.999))


################################################################################
### occurrence ratio

occ_cl <- tail_prob_cal(y, pgpd, t = t_vec, ratio = "occ", qu = T, xi = gamma)
occ_id <- tail_prob_cal(y, pexp, t = t_vec, ratio = "occ", qu = T, rate = delta)
occ_ex <- tail_prob_cal(y, pexp, t = t_vec, ratio = "occ", qu = T, rate = delta/v)

occ_cl <- occ_cl |> plot_ptc(ratio = "occ", ylims = c(0, 2))
occ_id <- occ_id |> plot_ptc(ratio = "occ", ylims = c(0, 2))
occ_ex <- occ_ex |> plot_ptc(ratio = "occ", ylims = c(0, 2))

## unconditional exceedance ratio plots
save_plots("occ")


################################################################################
### severity ratio

sev_cl <- tail_prob_cal(y, pgpd, t = rd_vec, ratio = "sev", xi = gamma)
sev_id <- tail_prob_cal(y, pexp, t = rd_vec, ratio = "sev", rate = delta)
sev_ex <- tail_prob_cal(y, pexp, t = rd_vec, ratio = "sev", rate = delta/v)

sev_cl <- sev_cl |> plot_ptc(ratio = "sev", names = names)
sev_id <- sev_id |> plot_ptc(ratio = "sev", names = names)
sev_ex <- sev_ex |> plot_ptc(ratio = "sev", names = names)


## cpit reliability diagrams
save_plots("sev")


################################################################################
### combined ratio

com_cl <- tail_prob_cal(y, pgpd, t = rd_vec, xi = gamma)
com_id <- tail_prob_cal(y, pexp, t = rd_vec, rate = delta)
com_ex <- tail_prob_cal(y, pexp, t = rd_vec, rate = delta/v)

com_cl <- com_cl |> plot_ptc(ratio = "com", names = names, ylims = c(0, 1.02), title = "Climatological")
com_id <- com_id |> plot_ptc(ratio = "com", names = names, ylims = c(0, 1.02), title = "Ideal")
com_ex <- com_ex |> plot_ptc(ratio = "com", names = names, ylims = c(0, 1.02), title = "Extremist")


## combined diagnostic plot
save_plots("com")


################################################################################
### conditional combined ratio

n_grp <- 3
group <- numeric(length(delta))
for (i in 1:n_grp) group[delta >= quantile(delta, (i - 1)/n_grp)] <- paste0("B", i)

com_cl <- tail_prob_div(y, pgpd, t = t_vec, group = group, qu = T, xi = gamma)
com_id <- tail_prob_div(y, pexp, t = t_vec, group = group, qu = T, rate = delta)
com_ex <- tail_prob_div(y, pexp, t = t_vec, group = group, qu = T, rate = delta/v)

com_b1 <- list("Clim." = com_cl[["B1"]], "Ideal" = com_id[["B1"]], "Extr." = com_ex[["B1"]])
com_b2 <- list("Clim." = com_cl[["B2"]], "Ideal" = com_id[["B2"]], "Extr." = com_ex[["B2"]])
com_b3 <- list("Clim." = com_cl[["B3"]], "Ideal" = com_id[["B3"]], "Extr." = com_ex[["B3"]])

com_div_b1 <- com_b1 |> plot_ptc_div(ylab = NULL, ylims = c(-0.1, 2.4), title = expression(B[1]))
com_div_b2 <- com_b2 |> plot_ptc_div(ylab = NULL, ylims = c(-0.1, 2.4), title = expression(B[2]))
com_div_b3 <- com_b3 |> plot_ptc_div(ylab = NULL, ylims = c(-0.1, 2.4), title = expression(B[3]))


## combined divergence plots
save_plots("com_div", mth = c("b1", "b2", "b3"), width = 3, height = 3)

