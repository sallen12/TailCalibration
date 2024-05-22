################################################################################
### set up

library(TailCalibration)
library(ggplot2)
library(evmix)

set.seed(412149)


################################################################################
### equiprobable forecaster

### plot
x <- seq(-4, 4, 0.01)
G_x <- pnorm(x, mean = 0, sd = 0.5)
H_x <- pnorm(x, mean = 0, sd = 2)
Fbar_x <- (G_x + H_x)/2

F1_x <- c(G_x[x <= 0], Fbar_x[x > 0])
F2_x <- c(H_x[x <= 0], Fbar_x[x > 0])

## unconditional distributions
df <- data.frame(x = x,
                 F_x = c(F1_x, F2_x, G_x, H_x),
                 mth = rep(c("F[1]", "F[2]", "Q(Y < y | F = F[1])", "Q(Y < y | F = F[2])"), each = length(x)))
ggplot(df) + geom_line(aes(x = x, y = F_x, col = mth, lty = mth)) +
  scale_color_manual(values = c("red", "blue", "red", "blue")) +
  scale_linetype_manual(values = c(1, 1, 2, 2)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.99, 0.01))


## excess distributions
t <- 1
Fbar_t <- (Fbar_x[x > t] - Fbar_x[x == t])/(1 - Fbar_x[x == t])
Ft_bar <- (G_x[x > t] - G_x[x == t])/(1 - G_x[x == t])/2 + (H_x[x > t] - H_x[x == t])/(1 - H_x[x == t])/2

df <- data.frame(x = x[x > t],
                 F_x = c(Fbar_t, Ft_bar),
                 mth = rep(c("Fbar_t", "Ft_bar"), each = sum(x > t)))
ggplot(df) + geom_line(aes(x = x, y = F_x, col = mth)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.99, 0.01))



### sample

N <- 1e6

F1 <- sample(c(T, F), size = N, replace = T)
y <- rnorm(N, mean = 0, sd = 0.5)
y[!F1] <- rnorm(sum(!F1), mean = 0, sd = 2)

F1_x <- function(x, sig1, sig2) {
  if (x <= 0) {
    pnorm(x, mean = 0, sd = sig1)
  } else {
    (pnorm(x, mean = 0, sd = sig1) + pnorm(x, mean = 0, sd = sig2))/2
  }
}

F2_x <- function(x, sig1, sig2) {
  if (x <= 0) {
    pnorm(x, mean = 0, sd = sig2)
  } else {
    (pnorm(x, mean = 0, sd = sig1) + pnorm(x, mean = 0, sd = sig2))/2
  }
}

F_x <- function(x, sig1, sig2, F1) {
  out <- numeric(length(F1))
  out[F1] <- F1_x(x, sig1, sig2)
  out[!F1] <- F2_x(x, sig1, sig2)
  return(out)
}

z_vec <- seq(0, 5, 0.01)
t_vec <- 1:5


### marginal tail calibration

mtc <- tail_marg_cal(y, F_x, t = t_vec, z = z_vec, sig1 = 0.5, sig2 = 2, F1 = F1)
plot_mtc(mtc)


### probabilistic tail calibration

ptc <- tail_prob_cal(y, F_x, t = t_vec, sig1 = 0.5, sig2 = 2, F1 = F1)
plot_mtc(mtc)


################################################################################
### Correct regime, wrong shape

N <- 1e8

y <- rgpd(N, sigmau = 1, xi = 0.5)

z_vec <- seq(0, 10, 1)
a_vec <- c(0.99, 0.999, 0.9999, 0.99999)
t_vec <- c(0, quantile(y, a_vec))

### marginal tail calibration

mtc <- tail_marg_cal(y, pgpd, t = t_vec, z = z_vec, sigmau = 2, xi = 0.5)
plot_mtc(mtc, names = c(0, a_vec), title = paste("scale =", 2, "shape =", 0.5))


### probabilistic tail calibration

ptc <- tail_prob_cal(y, pgpd, t = t_vec, sigmau = 2, xi = 0.7)
plot_ptc(ptc, names = c(0, a_vec), title = paste("scale =", 2, "shape =", 0.7))

