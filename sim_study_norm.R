################################################################################
### set up

library(ggplot2)
#devtools::install_github("sallen12/WeightedForecastVerification")
library(WeightedForecastVerification)
source("reldiag.R")

set.seed(631)


################################################################################
### simulate data

N <- 1e6

mu <- rnorm(N)
y <- rnorm(N, mu)

t_vec <- seq(-10, 10, 0.1)


################################################################################
### marginal tail calibration

rd_q <- c(0, 0.9, 0.95, 0.99, 0.999)
rd_vec <- quantile(y, rd_q)

mtc_clim <- lapply(rd_vec, function(t) {
  print(t)
  ind <- y > t
  Y <- seq(0, 10, 0.1)
  Fhat_t <- (pnorm(Y + t, 0, sqrt(2)) - pnorm(t, 0, sqrt(2)))/(1 - pnorm(t, 0, sqrt(2)))
  if (sum(ind) == 0) {
    Qhat_t <- NA
  } else if (sum(ind) == 1) {
    Qhat_t <- (y[ind] - t) <= Y
  } else {
    Qhat_t <- colMeans(outer(y[ind] - t, Y, FUN = function(x, y) x <= y))
  }
  data.frame(Y = Y, d = Fhat_t - Qhat_t)
})
mtc_id <- lapply(rd_vec, function(t) {
  print(t)
  ind <- y > t
  Y <- seq(0, 10, 0.1)
  Fhat_t <- colMeans(sapply(Y, function(z) (pnorm(z + t, mu[ind]) - pnorm(t, mu[ind]))/(1 - pnorm(t, mu[ind]))))
  if (sum(ind) == 0) {
    Qhat_t <- NA
  } else if (sum(ind) == 1) {
    Qhat_t <- (y[ind] - t) <= Y
  } else {
    Qhat_t <- colMeans(outer(y[ind] - t, Y, FUN = function(x, y) x <= y))
  }
  data.frame(Y = Y, d = Fhat_t - Qhat_t)
})
mtc_uf <- lapply(rd_vec, function(t) {
  print(t)
  ind <- y > t
  Y <- seq(0, 10, 0.1)
  tau <- sample(c(-1, 1), length(mu), replace = T)
  F_x <- function(x, m, ta) 0.5*pnorm(x, m) + 0.5*pnorm(x, m + ta)
  Fhat_t <- colMeans(sapply(Y, function(z)
    (F_x(z + t, mu[ind], tau[ind]) - F_x(t, mu[ind], tau[ind]))/(1 - F_x(t, mu[ind], tau[ind]))))
  if (sum(ind) == 0) {
    Qhat_t <- NA
  } else if (sum(ind) == 1) {
    Qhat_t <- (y[ind] - t) <= Y
  } else {
    Qhat_t <- colMeans(outer(y[ind] - t, Y, FUN = function(x, y) x <= y))
  }
  data.frame(Y = Y, d = Fhat_t - Qhat_t)
})
mtc_sr <- lapply(rd_vec, function(t) {
  print(t)
  ind <- y > t
  Y <- seq(0, 10, 0.1)
  Fhat_t <- colMeans(sapply(Y, function(z) (pnorm(z + t, -mu[ind]) - pnorm(t, -mu[ind]))/(1 - pnorm(t, -mu[ind]))), na.rm = T)
  if (sum(ind) == 0) {
    Qhat_t <- NA
  } else if (sum(ind) == 1) {
    Qhat_t <- (y[ind] - t) <= Y
  } else {
    Qhat_t <- colMeans(outer(y[ind] - t, Y, FUN = function(x, y) x <= y))
  }
  data.frame(Y = Y, d = Fhat_t - Qhat_t)
})


## marginal difference plots
df <- do.call(rbind, mtc_id)
df$t <- as.vector(sapply(1:length(rd_vec), function(i) rep(rd_q[i], nrow(mtc_clim[[i]]))))

ggplot(df) + geom_line(aes(x = Y, y = d, col = as.factor(t))) +
  geom_hline(aes(yintercept = 0), lty = "dotted") +
  scale_x_continuous(name = "y", expand = c(0, 0)) +
  scale_y_continuous(name = expression("E[" ~ F[t] ~ "(y)] - Q(Y - t ≤ y | Y > t)"),
                     limits = c(-0.2, 0.4), expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1, 1),
        legend.position = c(0.99, 0.99)) +
  guides(col = guide_legend(nrow = 2, byrow = TRUE))
ggsave("plots/mtc_dif_1e6_norm_id.png", width = 3.2, height = 3)


################################################################################
### probabilistic tail calibration

ptc_clim <- sapply(t_vec, function(t) {
  ind <- y > t
  if (sum(ind) == 0) {
    NA
  } else {
    (pnorm(y[ind], 0, sqrt(2)) - pnorm(t, 0, sqrt(2)))/(1 - pnorm(t, 0, sqrt(2)))
  }
})
ptc_id <- sapply(t_vec, function(t) {
  ind <- y > t
  if (sum(ind) == 0) {
    NA
  } else {
    (pnorm(y[ind], mu[ind]) - pnorm(t, mu[ind]))/(1 - pnorm(t, mu[ind]))
  }
})
ptc_uf <- sapply(t_vec, function(t) {
  ind <- y > t
  if (sum(ind) == 0) {
    NA
  } else {
    F_x <- function(x, m) 0.5*pnorm(x, m) + 0.5*pnorm(x, m + sample(c(-1, 1), length(m), replace = T))
    (F_x(y[ind], mu[ind]) - F_x(t, mu[ind]))/(1 - F_x(t, mu[ind]))
  }
})
ptc_sr <- sapply(t_vec, function(t) {
  ind <- y > t
  if (sum(ind) == 0) {
    NA
  } else {
    (pnorm(y[ind], -mu[ind]) - pnorm(t, -mu[ind]))/(1 - pnorm(t, -mu[ind]))
  }
})


## reliability diagrams
rd_vec <- sapply(rd_vec, function(z) t_vec[which.min(abs(z - t_vec))])

z <- ptc_clim[sapply(rd_vec, function(t) which(t_vec == t))]
names(z) <- rd_q
plt <- pit_reldiag(z, resampling = F)
plt <- plt + guides(col = guide_legend(ncol = 2))
ggsave("plots/ptc_rh_1e6_norm_cl.png", plt, width = 3.2, height = 3)

