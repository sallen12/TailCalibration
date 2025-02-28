################################################################################
###### examples

library(TailCalibration)
library(evd)
library(ggplot2)

set.seed(298301)


################################################################################
###### Example 15: Non-random forecaster

mu2 <- 0.5
sig2 <- 0.9
th <- mu2 / (1 - sig2)

F_x <- function(x, lower.tail = TRUE) {
  Fcst <- pgpd(x, 0, 1, 1/4)
  Fcst[x < th] <- pgpd(x[x < th], mu2, sig2, 1/4)
  if (lower.tail) {
    return(Fcst)
  } else {
    return(1 - Fcst)
  }
}

F_inv <- function(u) {
  u_th <- pgpd(th, 0, 1, 1/4)
  q <- qgpd(u, 0, 1, 1/4)
  if (any(u < u_th)) {
    q[u < u_th] <- qgpd(u[u < u_th], mu2, sig2, 1/4)
  }
  return(q)
}

t <- seq(0, 10, 0.01)
u <- seq(0.01, 0.99, 0.01)

denom <- F_x(t, lower.tail = FALSE)

F_t_inv <- function(u, t) F_inv(u * (1 - F_x(t)) + F_x(t)) - t
numer <- sapply(u, function(uu) pgpd(t + F_t_inv(uu, t), 0, 1, 1/4) - pgpd(t, 0, 1, 1/4))

ind <- c(1, 201, 401, 601)


##### combined ratio
rat <- numer[ind, ]/denom[ind]
cal <- lapply(1:length(ind), function(i) data.frame(u = u, rat = rat[i, ]))
names(cal) <- t[ind]

plot_ptc(cal)
ggsave("plots/ex_non_com.png", width = 3.2, height = 3)


##### severity ratio
rat <- numer[ind, ] / pgpd(t[ind], 0, 1, 1/4, lower.tail = FALSE)
cal <- lapply(1:length(ind), function(i) data.frame(u = u, rat = rat[i, ]))
names(cal) <- t[ind]

plot_ptc(cal, ratio = "sev")
ggsave("plots/ex_non_sev.png", width = 3.2, height = 3)


##### occurrence ratio
rat <- pgpd(t, 0, 1, 1/4, lower.tail = FALSE) / denom

df <- data.frame(t = t, r = rat)
plot_ptc(df, ratio = "occ", ylims = c(0, 2))
ggsave("plots/ex_non_occ.png", width = 3.2, height = 3)



################################################################################
###### Example 16: Uniform unfocused forecaster

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

ind <- c(1, 51, 101, 196)


##### combined ratio
rat <- G_t[ind] * numer[ind, ]/denom[ind]
cal <- lapply(1:length(ind), function(i) data.frame(u = u, rat = rat[i, ]))

plot_ptc(cal, ylims = c(0, 1.4), names = as.factor(t[ind]))
ggsave("plots/ex_uuf_com.png", width = 3.2, height = 3)


##### severity ratio
rat <- numer[ind, ]
cal <- lapply(1:length(ind), function(i) data.frame(u = u, rat = rat[i, ]))

plot_ptc(cal, ratio = "sev", names = as.factor(t[ind]))
ggsave("plots/ex_uuf_sev.png", width = 3.2, height = 3)


##### occurrence ratio
rat <- G_t / denom

df <- data.frame(t = t, r = rat)
plot_ptc(df, ratio = "occ", ylims = c(-0.2, 2.2))
ggsave("plots/ex_uuf_occ.png", width = 3.2, height = 3)



################################################################################
###### Example 18: Optimistic forecaster

N <- 1e6
gamma <- 1/4

delta <- rgamma(N, shape = 1/gamma, rate = 1/gamma)
u <- runif(N)
x1 <- qexp(u, delta)
x2 <- qexp(u, delta/2)
L <- rgpd(N, shape = gamma/2)

y <- pmax(pmin(x1, x2), pmin(pmax(x1, x2), L))

rd_q <- c(0.9, 0.95, 0.99, 0.999)
rd_vec <- c(0, quantile(y, rd_q))
names <- c(0, rd_q)

t_vec <- quantile(y, c(seq(0, 0.99, 0.01), 0.999))


##### unconditional combined ratio diagnostic plot

com <- tc_prob(y, pgpd, t = rd_vec, qu = T, shape = gamma)
occ <- tc_prob(y, pgpd, t = t_vec, ratio = "occ", qu = T, shape = gamma)
sev <- tc_prob(y, pgpd, t = rd_vec, ratio = "sev", qu = T, shape = gamma)

com <- com |> plot_ptc(names = names, ylims = c(0, 1.3), title = "")
occ <- occ |> plot_ptc(ratio = "occ", ylims = c(0, 2))
sev <- sev |> plot_ptc(ratio = "sev", names = names)

ggsave(plot = com, "plots/ex_opt_com_1e6.png", width = 3, height = 3)


##### conditional combined ratio diagnostic plots

n_grp <- 3
group <- numeric(length(delta))
for (i in 1:n_grp) group[delta >= quantile(delta, (i - 1)/n_grp)] <- paste0("B", i)

com_div <- tc_cprob(y, pgpd, t = t_vec, group = group, qu = T, shape = gamma)

com_div[["B1"]] |> plot_tc_sup(ylab = NULL, ylims = c(-0.1, 2.4), title = expression(B[1]))
ggsave("plots/ex_opt_com_div_1e6_B1.png", width = 3, height = 3)
com_div[["B2"]] |> plot_tc_sup(ylab = NULL, ylims = c(-0.1, 2.4), title = expression(B[2]))
ggsave("plots/ex_opt_com_div_1e6_B2.png", width = 3, height = 3)
com_div[["B3"]] |> plot_tc_sup(ylab = NULL, ylims = c(-0.1, 2.4), title = expression(B[3]))
ggsave("plots/ex_opt_com_div_1e6_B3.png", width = 3, height = 3)

