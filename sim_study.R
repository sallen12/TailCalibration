################################################################################
### set up

library(evmix)
library(ggplot2)
#devtools::install_github("sallen12/WeightedForecastVerification")
library(WeightedForecastVerification)

set.seed(89)


################################################################################
### simulate data

gamma <- 1/4
N <- 1e4
v <- 1.4

delta <- rgamma(N, shape = 1/gamma, rate = 1/gamma)
y <- rexp(N, rate = delta)

t_vec <- seq(0, 100, 0.1)
q_vec <- sapply(t_vec, function(t) mean(y <= t))

plot_calibration <- function(df, type = "", yref = 1, title = "", ylims = NULL) {
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
    scale_x_continuous(name = "Quantile", expand = c(0, 0)) + 
    scale_y_continuous(name = type, limits = ylims) +
    theme_bw() + theme(panel.grid = element_blank(),
                       legend.title = element_blank(),
                       legend.position = "bottom") +
    ggtitle("")
  gridExtra::grid.arrange(t_plot, q_plot)
}


################################################################################
### marginal tail calibration

mtc_clim <- sapply(t_vec, function(t) (y > t)/(1 - pgpd(t, 0, 1, gamma)))
mtc_id <- sapply(t_vec, function(t) (y > t)/(1 - pexp(t, delta)))
mtc_ex <- sapply(t_vec, function(t) (y > t)/(1 - pexp(t, delta/v)))

df <- data.frame(t = t_vec, 
                 q = q_vec,
                 mtc = c(colMeans(mtc_clim), colMeans(mtc_id), colMeans(mtc_ex)), 
                 mth = rep(c("Clim", "Ideal", "Extremist"), each = length(t_vec)))

plot_calibration(df, type = "Marginal tail calibration")


## condition on delta

del_n <- 5
del_quants <- seq(0, 1, length.out = del_n + 1)

cal_clim <- matrix(NA, nrow = del_n, ncol = length(t_vec))
del_plots <- list()
for (del_i in 1:del_n) {
  print(del_i)
  
  # get category limits
  del_t <- quantile(delta, del_quants[c(del_i, del_i + 1)])
  if (del_i == 1) del_t[1] <- 0
  del_ind <- delta > del_t[1] & delta <= del_t[2]

  # plot calibration 
  df <- data.frame(t = t_vec, 
                 q = sapply(t_vec, function(t) mean(y[del_ind] <= t)),
                 mtc = c(colMeans((mtc_clim - 1)*del_ind), 
                         colMeans((mtc_id - 1)*del_ind), 
                         colMeans((mtc_ex - 1)*del_ind)), 
                 mth = rep(c("Clim", "Ideal", "Extremist"), each = length(t_vec)))
  
  del_plots[[del_i]] <- plot_calibration(df, 
                                         type = "Marginal tail calibration", 
                                         yref = 0, 
                                         title = paste0("Min: ", round(del_t[1], 2), ", Max: ", round(del_t[2], 2)), 
                                         ylims = c(-0.5, 2.1))
  
  # check that average over all delta is the overall marginal tail calibration
  cal_clim[del_i, ] <- colSums((mtc_clim - 1)*del_ind)
}

do.call(gridExtra::grid.arrange, c(del_plots, nrow = 1))

# check that average over all delta is the overall marginal tail calibration
plot(t_vec, colSums(cal_clim)/N)
lines(t_vec, colMeans(mtc_clim - 1), col = "blue")


################################################################################
### climatological marginal tail calibration

cmtc_clim <- sapply(t_vec, function(t) mean(y > t)/(1 - pgpd(t, 0, 1, gamma)))
cmtc_id <- sapply(t_vec, function(t) mean(y > t)/mean(1 - pexp(t, delta)))
cmtc_ex <- sapply(t_vec, function(t) mean(y > t)/mean(1 - pexp(t, delta/v)))

df <- data.frame(t = t_vec, 
                 q = q_vec,
                 mtc = c(cmtc_clim, cmtc_id, cmtc_ex), 
                 mth = rep(c("Clim", "Ideal", "Extremist"), each = length(t_vec)))

plot_calibration(df, type = "Climatological marginal tail calibration")


################################################################################
### probabilistic tail calibration

ptc_clim <- sapply(t_vec, function(t) {
  ind <- y > t
  if (sum(ind) == 0) {
    NA
  } else {
    (pgpd(y[ind], 0, 1, gamma) - pgpd(t, 0, 1, gamma))/(1 - pgpd(t, 0, 1, gamma))
  }
})
ptc_id <- sapply(t_vec, function(t) {
  ind <- y > t
  if (sum(ind) == 0) {
    NA
  } else {
    (pexp(y[ind], delta[ind]) - pexp(t, delta[ind]))/(1 - pexp(t, delta[ind]))
  }
})
ptc_ex <- sapply(t_vec, function(t) {
  ind <- y > t
  if (sum(ind) == 0) {
    NA
  } else {
    (pexp(y[ind], delta[ind]/v) - pexp(t, delta[ind]/v))/(1 - pexp(t, delta[ind]/v))
  }
})


## mean
ptc_clim_mn <- sapply(ptc_clim, mean)
ptc_id_mn <- sapply(ptc_id, mean)
ptc_ex_mn <- sapply(ptc_ex, mean)

df <- data.frame(t = t_vec, 
                 q = q_vec,
                 mtc = c(ptc_clim_mn, ptc_id_mn, ptc_ex_mn), 
                 mth = rep(c("Clim", "Ideal", "Extremist"), each = length(t_vec)))

plot_calibration(df, type = "Probabilistic tail calibration", 0.5)


## variance
ptc_clim_var <- sapply(ptc_clim, var)
ptc_id_var <- sapply(ptc_id, var)
ptc_ex_var <- sapply(ptc_ex, var)

df <- data.frame(t = t_vec, 
                 q = q_vec,
                 mtc = c(ptc_clim_var, ptc_id_var, ptc_ex_var), 
                 mth = rep(c("Clim", "Ideal", "Extremist"), each = length(t_vec)))

plot_calibration(df, type = "Probabilistic tail calibration", yref = 1/12)


## histogram
t <- round(quantile(y, 0.99), 1)
pitplot_clim <- pit_hist(ptc_clim[[which(t_vec == t)]], ranks = F, ymax = 0.25, title = "Climatological")
pitplot_id <- pit_hist(ptc_id[[which(t_vec == t)]], ranks = F, ymax = 0.25, title = "Ideal")
pitplot_ex <- pit_hist(ptc_ex[[which(t_vec == t)]], ranks = F, ymax = 0.25, title = "Extremist")
gridExtra::grid.arrange(pitplot_clim, pitplot_id, pitplot_ex)


## divergence
Rcpp::sourceCpp("crps_div.cpp")

ptc_clim_div <- sapply(ptc_clim, crps_div)
ptc_id_div <- sapply(ptc_id, crps_div)
ptc_ex_div <- sapply(ptc_ex, crps_div)

df <- data.frame(t = t_vec, 
                 q = q_vec,
                 mtc = c(ptc_clim_div, ptc_id_div, ptc_ex_div), 
                 mth = rep(c("Clim", "Ideal", "Extremist"), each = length(t_vec)))

plot_calibration(df, type = "Probabilistic tail calibration", yref = 0)


## condition on delta
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
    plot_calibration(df, type = "Probabilistic tail calibration", yref = 0, 
                     title = paste0("Min: ", round(del_t[1], 2), ", Max: ", round(del_t[2], 2)))
}

do.call(gridExtra::grid.arrange, c(del_plots, nrow = 1))


## reliability diagrams
rd_vec <- seq(0, 20, 5)

# clim
z <- ptc_clim[sapply(rd_vec, function(t) which(t_vec == t))]
names(z) <- sprintf("%02d", rd_vec)
pit_reldiag(z)

# id
z <- ptc_id[sapply(rd_vec, function(t) which(t_vec == t))]
names(z) <- sprintf("%02d", rd_vec)
pit_reldiag(z)

# ext
z <- ptc_ex[sapply(rd_vec, function(t) which(t_vec == t))]
names(z) <- sprintf("%02d", rd_vec)
pit_reldiag(z)


### climatological probabilistic tail calibration (EXPENSIVE)

# cptc_clim <- ptc_clim
# cptc_id <- sapply(t_vec, function(t) {
#   print(t)
#   ind <- y > t
#   if (sum(ind) == 0) {
#     NA
#   } else {
#     mean((sapply(y[ind], function(z) mean(pexp(z, delta))) - mean(pexp(t, delta)))/(1 - mean(pexp(t, delta))))
#   }
# })
# cptc_ex <- sapply(t_vec, function(t) {
#   ind <- y > t
#   if (sum(ind) == 0) {
#     NA
#   } else {
#     mean((sapply(y[ind], function(z) mean(pexp(z, delta/v))) - mean(pexp(t, delta/v)))/(1 - mean(pexp(t, delta/v))))
#   }
# })
# 
# df <- data.frame(t = t_vec, 
#                  q = q_vec,
#                  mtc = c(cptc_clim, cptc_id, cptc_ex), 
#                  mth = rep(c("Clim", "Ideal", "Extremist"), each = length(t_vec)))
# 
# plot_calibration(df, type = "Climatological probabilistic tail calibration")



### quantile calibration

source("reldiag.R")

N <- 1e4
delta <- rgamma(N, shape = 1/gamma, rate = 1/gamma)
y <- rexp(N, rate = delta)

alpha <- 0.99

F_q_clim <- rep(qgpd(alpha, 0, 1, gamma), N)
F_q_id <- qexp(alpha, delta)
F_q_ex <- qexp(alpha, delta/v)

qc_clim <- reldiag(x = F_q_clim, y = y, type = list("quantile", alpha), 
                   resampling = F, ext_decomp = F, inset_hist = F, lim = range(y))
qc_id <- reldiag(x = F_q_id, y = y, type = list("quantile", alpha), 
                 resampling = F, ext_decomp = F, inset_hist = F, lim = range(y))
qc_ex <- reldiag(x = F_q_ex, y = y, type = list("quantile", alpha), 
                 resampling = F, ext_decomp = F, inset_hist = F, lim = range(y))


### threshold calibration

t <- round(quantile(y, 0.99), 1)

F_t_clim <- list(F = function(x) rep(pgpd(x, 0, 1, gamma), N))
F_t_id <- list(F = function(x) pexp(x, delta))
F_t_ex <- list(F = function(x) pexp(x, delta/v))

tc_clim <- threshreldiag(fcast = F_t_clim, y = y, t = t, inset_hist = F)
tc_id <- threshreldiag(fcast = F_t_id, y = y, t = t, inset_hist = F)
tc_ex <- threshreldiag(fcast = F_t_ex, y = y, t = t, inset_hist = F)













