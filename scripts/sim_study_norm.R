################################################################################
### set up

#devtools::install_github("sallen12/WeightedForecastVerification")
library(TailCalibration)
library(ggplot2)

set.seed(631)


################################################################################
### simulate data

N <- 1e6

mu <- rnorm(N)
y <- rnorm(N, mu)

tau <- sample(c(-1, 1), length(mu), replace = T)
F_uf <- function(x, m, ta) 0.5*pnorm(x, m) + 0.5*pnorm(x, m + ta)

rd_q <- c(0, 0.9, 0.95, 0.99, 0.999)
rd_vec <- quantile(y, rd_q)

t_vec <- quantile(y, c(seq(0, 1, 0.01), 0.999))
z_vec <- seq(0, 10, 0.1)


################################################################################
### reliability diagrams

rel_diag <- function(y, F_x, t, names = NULL, ...) {
  d <- lapply(t, function(tt) {
    p <- 1 - F_x(tt, ...)
    if (length(p) == 1) p <- rep(p, length(y))
    o <- as.numeric(y > tt)
    out <- reliabilitydiag::reliabilitydiag(p, y = o)[[1]]$bins
    if (!is.null(names)) {
      name <- names[which(t == tt)]
    } else {
      name <- format(round(tt, 2), nsmall = 2)
    }
    return(data.frame(x = c(out$x_min, out$x_max[length(out$x_max)]),
                      y = c(out$CEP_pav, out$CEP_pav[length(out$CEP_pav)]),
                      t = name))
  })
  d <- do.call(rbind, d)
  rd_plot <- ggplot(d) +
    geom_step(aes(x = x, y = y, col = as.factor(t))) +
    geom_abline(aes(slope = 1, intercept = 0), linetype = "dotted") +
    scale_x_continuous(name = "p", limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(name = "CEP", limits = c(0, 1), expand = c(0, 0)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.justification = c(0, 1),
          legend.position = c(0.01, 0.99))

  return(rd_plot)
}

rd_cl <- rel_diag(y, pnorm, t = rd_vec, sd = sqrt(2), names = rd_q)
rd_id <- rel_diag(y, pnorm, t = rd_vec, mean = mu, names = rd_q)
rd_uf <- rel_diag(y, F_uf, t = rd_vec, m = mu, ta = tau, names = rd_q)
rd_sr <- rel_diag(y, pnorm, t = rd_vec, mean = -mu, names = rd_q)


## reliability diagrams

ggsave(plot = rd_uf, "plots/rd_1e6_norm_uf.png", width = 3.2, height = 3)


################################################################################
### unconditional exceedance ratio plots

unc_exc_rat_plot <- function(y, F_x, t, title = "", ylims = c(0, 2), xlims = c(0, 1), qu = TRUE, ...) {

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
    theme(panel.grid = element_blank(),
          plot.margin = margin(c(5.5, 10.5, 5.5, 5.5), "points")) +
    ggtitle(title)

  return(plt)
}

uer_cl <- unc_exc_rat_plot(y, pnorm, t = t_vec, sd = sqrt(2))#, xlims = c(0.7, 1))
uer_id <- unc_exc_rat_plot(y, pnorm, t = t_vec, mean = mu)#, xlims = c(0.7, 1))
uer_uf <- unc_exc_rat_plot(y, F_uf, t = t_vec, m = mu, ta = tau)#, xlims = c(0.7, 1))
uer_sr <- unc_exc_rat_plot(y, pnorm, t = t_vec, mean = -mu)#, xlims = c(0.7, 1))


## unconditional exceedance ratio plots

ggsave(plot = uer_uf, "plots/uer_1e6_norm_uf.png", width = 3.5, height = 3.4)


################################################################################
### marginal tail calibration

mtc_cl <- tail_marg_cal(y, pnorm, t = t_vec, z = z_vec, sd = sqrt(2))
mtc_id <- tail_marg_cal(y, pnorm, t = t_vec, z = z_vec, mean = mu)
mtc_uf <- tail_marg_cal(y, F_uf, t = t_vec, z = z_vec, m = mu, ta = tau)
mtc_sr <- tail_marg_cal(y, pnorm, t = t_vec, z = z_vec, mean = -mu)

## marginal difference plots

plot_mtc(mtc_sr, xlab = expression(alpha), ylab = "Sup. difference")
#ggsave("plots/mtc_dif_1e6_norm_id.png", width = 3.2, height = 3)


################################################################################
### probabilistic tail calibration

ptc_cl <- tail_prob_cal(y, pnorm, t = rd_vec, sd = sqrt(2))
ptc_id <- tail_prob_cal(y, pnorm, t = rd_vec, mean = mu)
ptc_uf <- tail_prob_cal(y, F_uf, t = rd_vec, m = mu, ta = tau)
ptc_sr <- tail_prob_cal(y, pnorm, t = rd_vec, mean = -mu)


## cpit reliability diagrams

plot_ptc(ptc_cl, names = rd_q)
#ggsave("plots/ptc_rh_1e6_norm_id.png", width = 3.5, height = 3.4)


################################################################################
### combined ratio

get_comb_ratio <- function(y, F_x, t, ...) {

  u <- seq(0, 1, 0.01)

  rat_list <- lapply(t, function(tt) {
    print(tt)
    numer <- sapply(u, function(uu) {
      pit <- (F_x(y, ...) - F_x(tt, ...)) / (1 - F_x(tt, ...))
      exc <- y > tt
      return(mean(pit <= uu & exc, na.rm=T))
    })
    denom <- mean(1 - F_x(tt, ...))
    df <- data.frame(cpit = numer/denom)
    return(df)
  })

  return(rat_list)
}

plot_comb_ratio <- function(comb, t, title = "", u = seq(0, 1, 0.01), ylims = c(0, 1)) {
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

comb_cl <- get_comb_ratio(y, pnorm, t = rd_vec, sd = sqrt(2))
comb_id <- get_comb_ratio(y, pnorm, t = rd_vec, mean = mu)
comb_uf <- get_comb_ratio(y, F_uf, t = rd_vec, m = mu, ta = tau)
comb_sr <- get_comb_ratio(y, pnorm, t = rd_vec, mean = -mu)

## combined diagnostic plot
plot_comb_ratio(comb_uf, t = rd_q, title = "Unfocused")
ggsave("plots/cmb_1e6_norm_uf.png", width = 3.5, height = 3.4)

