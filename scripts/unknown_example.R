################################################################################
### set up

library(TailCalibration)
library(evmix)
library(ggplot2)

set.seed(89)


################################################################################
### simulate data

gamma <- 1/4
gamma2 <- gamma*0.8
N <- 1e4
delta <- rgamma(N, shape = 1/gamma, rate = 1/gamma)

G_1 <- function(x, delta) pexp(x, rate = delta)
G_2 <- function(x, delta) pexp(x, rate = delta/2)
G_2bar <- function(x) 1 - (2/gamma)^(1/gamma) * ((2/gamma) + x)^(-1/gamma)
lambda <- function(x) pgpd(x, u = 0, sigmau = 1, xi = gamma2)

G <- function(x) lambda(x) * pgpd(x, u = 0, sigmau = 1, xi = gamma) + (1 - lambda(x)) * G_2bar(x)

xseq <- seq(0, 20, 0.001)
plot(xseq, G(xseq), type = "l")

# quantile function
Gq <- function(u) xseq[which.min(abs(G(xseq) - u))]

# sample from the quantile function
y <- sapply(1:N, function(i) {
  print(i)
  Gq(runif(1))
})




rd_q <- c(0, 0.9, 0.95, 0.99, 0.999)
rd_vec <- quantile(y, rd_q)

t_vec <- quantile(y, c(seq(0, 1, 0.01), 0.999))
z_vec <- seq(0, 20, 0.1)




# plot for F = G1

# ptc
ptc <- tail_prob_cal(y, G_1, t = rd_vec, delta = delta)
plot_ptc(ptc, names = rd_q)

# mtc
mtc <- tail_marg_cal(y, G_1, t = rd_vec, z = z_vec, delta = delta)
plot_mtc(mtc, xlab = expression(alpha), ylab = "Sup. difference")



