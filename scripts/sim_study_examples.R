################################################################################
###### examples

###### non-random

library(evmix)
library(ggplot2)

mu2 <- 0.5
sig2 <- 0.9
th <- mu2 / (1 - sig2)

Fcst <- function(x, lower.tail = TRUE) {
  Fcst <- pgpd(x, 0, 1, 1/4)
  Fcst[x < th] <- pgpd(x[x < th], mu2, sig2, 1/4)
  if (lower.tail) {
    return(Fcst)
  } else {
    return(1 - Fcst)
  }
}

x <- seq(0, 10, 0.001)
plot(x, pgpd(x, 0, 1, 1/4), type = "l")
lines(x, pgpd(x, mu2, sig2, 1/4), col = "red")
lines(x, Fcst(x), col = "blue")
abline(v = th)

t <- seq(0, 10, 0.01)

u <- seq(0.01, 0.99, 0.01)
denom <- Fcst(t, lower.tail = FALSE)
F_inv <- function(u) {
  u_th <- pgpd(th, 0, 1, 1/4)
  q <- qgpd(u, 0, 1, 1/4)
  if (any(u < u_th)) {
    q[u < u_th] <- qgpd(u[u < u_th], mu2, sig2, 1/4)
  }
  return(q)
}
F_t_inv <- function(u, t) F_inv(u * (1 - Fcst(t)) + Fcst(t)) - t
numer <- sapply(u, function(uu) pgpd(t + F_t_inv(uu, t), 0, 1, 1/4) - pgpd(t, 0, 1, 1/4))


## overall
ind <- c(1, 201, 401, 601)

rat <- numer/denom

df <- data.frame(r = as.vector(t(rat[ind, ])),
                 u = u,
                 t = as.factor(rep(t[ind], each = length(u))))
ggplot(df) + geom_line(aes(x = u, y = r, col = t)) +
  geom_abline(aes(slope = 1, intercept = 0), linetype = "dotted") +
  scale_x_continuous(name = "u", breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
  scale_y_continuous(name = "Combined ratio", breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(0, 1),
        legend.position = c(0.01, 0.99))
ggsave("plots/sim_nonrandom_all.png", width = 3.2, height = 3)


## occurrence
rat <- pgpd(t, 0, 1, 1/4, lower.tail = FALSE) / denom

df <- data.frame(t = t, r = rat)
ggplot(df) + geom_line(aes(x = t, y = r)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted") +
  scale_x_continuous(name = "t", breaks = seq(0, 10, 1), expand = c(0, 0)) +
  scale_y_continuous(name = "Occurrence ratio", limits = c(0, 2), expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave("plots/sim_nonrandom_occ.png", width = 3.2, height = 3)


## severity
rat <- numer / pgpd(t, 0, 1, 1/4, lower.tail = FALSE)

df <- data.frame(r = as.vector(t(rat[ind, ])),
                 u = u,
                 t = as.factor(rep(t[ind], each = length(u))))
ggplot(df) + geom_line(aes(x = u, y = r, col = t)) +
  geom_abline(aes(slope = 1, intercept = 0), linetype = "dotted") +
  scale_x_continuous(name = "u", breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
  scale_y_continuous(name = "Severity ratio", breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(0, 1),
        legend.position = c(0.01, 0.99))
ggsave("plots/sim_nonrandom_sev.png", width = 3.2, height = 3)




###### unfocused (uniform)

u <- seq(0, 1, 0.01)
t <- seq(-1, 0.99, 0.01)

E_F <- function(y) {
  F_y <- numeric(length(y))
  F_y[y <= -1] <- 0
  F_y[y >= 2] <- 1
  ind <- y > -1 & y <= 1
  F_y[ind] <- (y[ind] + 1)/4
  ind <- y > 1 & y < 2
  F_y[ind] <- (y[ind] + 2)/4
  return(F_y)
}
Q_t <- function(u, t) {
  if (t <= -1) {
    u
  } else if (t > -1 & t <= 0) {
    (pmin(2*u, 1) + pmax(u*(1 - t) + t, 0))/2
  } else if (t > 0) {
    (pmin(u*(2 - t)/(1 - t), 1) + u)/2
  } else {
    NA
  }
}

denom <- 1 - E_F(t)
G_t <- pmax(pmin(1 - t, 1), 0)
numer <- t(sapply(t, Q_t, u = u))


## overall
ind <- c(1, 51, 101, 196)

rat <- G_t * numer/denom

df <- data.frame(r = as.vector(t(rat[ind, ])),
                 u = u,
                 t = as.factor(rep(t[ind], each = length(u))))
ggplot(df) + geom_line(aes(x = u, y = r, col = t)) +
  geom_abline(aes(slope = 1, intercept = 0), linetype = "dotted") +
  scale_x_continuous(name = "u", breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
  scale_y_continuous(name = "Combined ratio", expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(0, 1),
        legend.position = c(0.01, 0.99))
ggsave("plots/sim_unfocused_all.png", width = 3.2, height = 3)


## occurrence
rat <- G_t / denom

df <- data.frame(t = t, r = rat)
ggplot(df) + geom_line(aes(x = t, y = r)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted") +
  scale_x_continuous(name = "t", limits = c(-1, 1), expand = c(0, 0)) +
  scale_y_continuous(name = "Occurrence ratio", limits = c(0, 2)) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave("plots/sim_unfocused_occ.png", width = 3.2, height = 3)


## severity
df <- data.frame(r = as.vector(t(numer[ind, ])),
                 u = u,
                 t = rep(as.factor(t[ind]), each = length(u)))
                 #t = rep(c("A", "B"), each = length(u)))
ggplot(df) + geom_line(aes(x = u, y = r, col = t)) +
  geom_abline(aes(slope = 1, intercept = 0), linetype = "dotted") +
  scale_x_continuous(name = "u", breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
  scale_y_continuous(name = "Severity ratio", breaks = seq(0, 1, 0.1), limits = c(0, 1), expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.99, 0.01))
ggsave("plots/sim_unfocused_sev.png", width = 3.2, height = 3)



###### stochastic

library(evmix)
library(TailCalibration)

n <- 1e6
gam <- 1/4

delta <- rgamma(n, shape = 1/gam, rate = 1/gam)
u <- runif(n)
x1 <- qexp(u, delta)
x2 <- qexp(u, delta/2)
L <- rgpd(n, xi = gam/2)

y <- apply(cbind(x1, x2, L), 1, function(z) sort(z)[2])


rd_q <- c(0, 0.9, 0.95, 0.99, 0.999)
rd_vec <- quantile(y, rd_q)

t_vec <- quantile(y, c(seq(0, 0.99, 0.01), 0.999))



## unconditional exceedance ratio plots
unc_exc_rat_plot <- function(y, F_x, t, ylims = c(0, 2), xlims = c(0, 1), qu = TRUE, ...) {

  y_exc_p <- sapply(t, function(tt) mean(y > tt))
  F_exc_p <- sapply(t, function(tt) mean(1 - F_x(tt, ...)))

  if (qu) {
    a <- sapply(t, function(tt) mean(y <= tt))
    xlab <- expression(alpha)
  } else {
    a <- t
    xlab <- "Threshold"
  }

  df <- data.frame(rat = y_exc_p/F_exc_p, t = a)
  plt <- ggplot(df) + geom_line(aes(x = t, y = rat)) +
    geom_hline(aes(yintercept = 1), linetype = "dotted") +
    scale_x_continuous(name = xlab, limits = xlims, expand = c(0, 0)) +
    #scale_y_continuous(name = "Q(Y > t)/E[1 - F(t)]", limits = ylims) +
    scale_y_continuous(name = "Occurrence ratio", limits = ylims) +
    theme_bw() +
    theme(panel.grid = element_blank())

  return(plt)
}

unc_exc_rat_plot(y, pgpd, t = t_vec, xi = gam)#, xlims = c(0.7, 1))
ggsave("plots/sim_stoch_occ.png", width = 3.2, height = 3)


## cpit reliability diagrams
ptc <- tail_prob_cal(y, pgpd, t = rd_vec, xi = gam)
plot_ptc(ptc, names = rd_q)
ggsave("plots/sim_stoch_sev.png", width = 3.2, height = 3)


## combined diagnostic plot
get_comb_ratio <- function(y, F_x, t, ...) {

  u <- seq(0, 1, 0.01)

  rat_list <- lapply(t, function(tt) {
    print(tt)
    pit <- (F_x(y, ...) - F_x(tt, ...)) / (1 - F_x(tt, ...))
    exc <- y > tt
    numer <- sapply(u, function(uu) mean(pit <= uu & exc, na.rm=T))
    denom <- mean(1 - F_x(tt, ...))
    df <- data.frame(cpit = numer/denom)
    return(df)
  })

  return(rat_list)
}
plot_comb_ratio <- function(comb, t, u = seq(0, 1, 0.01), ylims = c(0, 1), title = "") {
  df <- data.frame(u = u,
                   r = as.vector(sapply(comb, function(i) i$cpit)),
                   mth = as.factor(rep(t, each = length(u))))
  comb_plot <- ggplot(df) + geom_line(aes(x = u, y = r, col = mth)) +
    geom_abline(aes(intercept = 0, slope = 1), lty = "dotted") +
    scale_x_continuous(name = "u", expand = c(0, 0), limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(name = "Combined ratio", limits = ylims, expand = c(0, 0)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.justification = c(0, 1),
          legend.position = c(0.01, 0.99),
          plot.margin = margin(c(5.5, 10.5, 5.5, 5.5), "points")) +
    ggtitle(title)
  return(comb_plot)
}


comb <- get_comb_ratio(y, pgpd, t = rd_vec, xi = gam)
plot_comb_ratio(comb, t = rd_q, ylims = c(0, 1.3))
ggsave("plots/sim_stoch_all.png", width = 3, height = 3)



################################################################################
### combined ratio div plots

plot_comb_div <- function(div, t, u = seq(0, 1, 0.01), ylims = NULL,
                          xlims = range(u), ylab = NULL, xlab = "t", title = NULL) {

  df <- data.frame(t = t, r = div)
  div_plot <- ggplot(df) + geom_line(aes(x = t, y = r)) +
    geom_hline(aes(yintercept = 0), linetype = "dotted") +
    scale_x_continuous(name = xlab, limits = xlims, expand = c(0, 0)) +
    scale_y_continuous(name = ylab, limits = ylims) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.margin = margin(c(5.5, 10.5, 5.5, 5.5), "points")) +
    ggtitle(title)
  return(div_plot)
}

plot_comb_ratio_div <- function(y, gam, t, grp, n_grp = 1, u = seq(0, 1, 0.01), names = NULL, ylims = NULL, xlims = NULL,
                                ylab = NULL, xlab = NULL, title = NULL, quantile_sc = T) {

  quants <- seq(0, 1, length.out = n_grp + 1)

  cond_plots <- vector("list", n_grp)
  for (i_g in 1:n_grp) {

    # get category limits
    grp_ab <- quantile(grp, quants[c(i_g, i_g + 1)])

    ind <- (grp > grp_ab[1]) & (grp <= grp_ab[2])

    # get calibration
    comb <- get_comb_ratio(y[ind], pgpd, t = t, xi = gam)

    div <- sapply(comb, function(x) max(abs(x$cpit - u)))

    # convert thresholds to quantile scale
    if (quantile_sc) {
      q_vec <- sapply(t, function(tt) mean(y[(grp > grp_ab[1]) & (grp <= grp_ab[2])] <= tt))
    } else {
      q_vec <- t
    }

    # plot calibration
    cond_plots[[i_g]] <- plot_comb_div(div, q_vec, u, ylims, xlims, ylab, xlab, title)
  }

  return(cond_plots)
}



## combined divergence diagnostic plot
comb_div_plots <- plot_comb_ratio_div(y, gam, t = t_vec, grp = delta, n_grp = 3,
                                      ylims = c(0, 2.2), xlab = "t")
comb_div_plots[[3]] + ggtitle(expression(B[3]))
ggsave("plots/sim_stoch_cond3.png", width = 3, height = 3)

