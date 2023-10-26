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
N <- 1e6
v <- 1.4

delta <- rgamma(N, shape = 1/gamma, rate = 1/gamma)
delta_mi <- rgamma(N, shape = 1/gamma, rate = 1/gamma)
y <- rexp(N, rate = delta)

t_vec <- seq(0, 100, 0.1)
q_vec <- sapply(t_vec, function(t) mean(y <= t))

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

mtc_clim <- lapply(rd_vec, function(t) {
  print(t)
  ind <- y > t
  if (sum(ind) == 0) {
    NA
  } else {
    Y <- t_vec
    Fhat_t <- (pgpd(Y + t, 0, 1, gamma) - pgpd(t, 0, 1, gamma))/(1 - pgpd(t, 0, 1, gamma))
    if (sum(ind) == 1) {
      Qhat_t <- (y[ind] - t) <= Y
    } else {
      Qhat_t <- colMeans(outer(y[ind] - t, Y, FUN = function(x, y) x <= y))
    }
    data.frame(Y = Y, d = Fhat_t - Qhat_t) # need Fhat - Qhat = 0 for all Y > t
  }
})
mtc_id <- lapply(rd_vec, function(t) {
  print(t)
  ind <- y > t
  Y <- t_vec
  Fhat_t <- colMeans(sapply(Y, function(z) (pexp(z + t, delta[ind]) - pexp(t, delta[ind]))/(1 - pexp(t, delta[ind]))))
  if (sum(ind) == 0) {
    Qhat_t <- NA
  } else if (sum(ind) == 1) {
    Qhat_t <- (y[ind] - t) <= Y
  } else {
    Qhat_t <- colMeans(outer(y[ind] - t, Y, FUN = function(x, y) x <= y))
  }
  data.frame(Y = Y, d = Fhat_t - Qhat_t) # need Fhat - Qhat = 0 for all Y > t
})
mtc_ex <- lapply(rd_vec, function(t) {
  print(t)
  ind <- y > t
  Y <- t_vec
  Fhat_t <- colMeans(sapply(Y, function(z) (pexp(z + t, delta[ind]/v) - pexp(t, delta[ind]/v))/(1 - pexp(t, delta[ind]/v))))
  if (sum(ind) == 0) {
    Qhat_t <- NA
  } else if (sum(ind) == 1) {
    Qhat_t <- (y[ind] - t) <= Y
  } else {
    Qhat_t <- colMeans(outer(y[ind] - t, Y, FUN = function(x, y) x <= y))
  }
  data.frame(Y = Y, d = Fhat_t - Qhat_t) # need Fhat - Qhat = 0 for all Y > t
})
mtc_mi <- lapply(rd_vec, function(t) {
  print(t)
  ind <- y > t
  Y <- t_vec
  Fhat_t <- colMeans(sapply(Y, function(z) (pexp(z + t, delta_mi[ind]) - pexp(t, delta_mi[ind]))/(1 - pexp(t, delta_mi[ind]))))
  if (sum(ind) == 0) {
    Qhat_t <- NA
  } else if (sum(ind) == 1) {
    Qhat_t <- (y[ind] - t) <= Y
  } else {
    Qhat_t <- colMeans(outer(y[ind] - t, Y, FUN = function(x, y) x <= y))
  }
  data.frame(Y = Y, d = Fhat_t - Qhat_t) # need Fhat - Qhat = 0 for all Y > t
})

## marginal difference plots
df <- do.call(rbind, mtc_id)
df$t <- rep(rd_q, each = length(t_vec))

ggplot(df) + geom_line(aes(x = Y, y = d, col = as.factor(t))) +
  geom_hline(aes(yintercept = 0), lty = "dotted") +
  scale_x_continuous(name = "y", expand = c(0, 0)) +
  scale_y_continuous(name = expression("E[" ~ F[t] ~ "(y)] - Q(Y - t ≤ y | Y > t)"),
                     limits = c(-0.2, 0.2), expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1, 1),
        legend.position = c(0.99, 0.99)) +
  guides(col = guide_legend(nrow = 2, byrow = TRUE))
ggsave("plots/mtc_dif_1e6_mi.png", width = 3.2, height = 3)


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
ptc_mi <- sapply(t_vec, function(t) {
  ind <- y > t
  if (sum(ind) == 0) {
    NA
  } else {
    (pexp(y[ind], delta_mi[ind]) - pexp(t, delta_mi[ind]))/(1 - pexp(t, delta_mi[ind]))
  }
})

## reliability diagrams
rd_q <- c(0, 0.9, 0.95, 0.99, 0.999)
rd_vec <- quantile(y, rd_q)
rd_vec <- sapply(rd_vec, function(z) t_vec[which.min(abs(z - t_vec))])

z <- ptc_clim[sapply(rd_vec, function(t) which(t_vec == t))]
names(z) <- rd_q
pit_reldiag(z, resampling = F)
ggsave("plots/ptc_rh_1e6_cl.png", width = 3.2, height = 3)


## condition divergence on delta
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


