################################################################################
### set up

library(TailCalibration)
library(evmix)
library(ggplot2)

set.seed(89)


################################################################################
### simulate data

gamma <- 1/4
N <- 1e5
v <- 1.4

delta <- rgamma(N, shape = 1/gamma, rate = 1/gamma)
y <- rexp(N, rate = delta)

rd_q <- c(0, 0.9, 0.95, 0.99, 0.999)
rd_vec <- quantile(y, rd_q)

t_vec <- quantile(y, c(seq(0, 0.99, 0.01), 0.999))
z_vec <- seq(0, 20, 0.1)


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

rd_cl <- rel_diag(y, pgpd, t = rd_vec, xi = gamma, names = rd_q)
rd_id <- rel_diag(y, pexp, t = rd_vec, rate = delta, names = rd_q)
rd_ex <- rel_diag(y, pexp, t = rd_vec, rate = delta/v, names = rd_q)


## reliability diagrams

ggsave(plot = rd_ex, "plots/rd_1e6_ex.png", width = 3.2, height = 3)


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

uer_cl <- unc_exc_rat_plot(y, pgpd, t = t_vec, xi = gamma)#, xlims = c(0.7, 1))
uer_id <- unc_exc_rat_plot(y, pexp, t = t_vec, rate = delta)#, xlims = c(0.7, 1))
uer_ex <- unc_exc_rat_plot(y, pexp, t = t_vec, rate = delta/v)#, xlims = c(0.7, 1))


## unconditional exceedance ratio plots

ggsave(plot = uer_id, "plots/uer_1e6_id.png", width = 3.5, height = 3.4)


################################################################################
### marginal tail calibration

mtc_cl <- tail_marg_cal(y, pgpd, t = t_vec, z = z_vec, xi = gamma)
mtc_id <- tail_marg_cal(y, pexp, t = t_vec, z = z_vec, rate = delta)
mtc_ex <- tail_marg_cal(y, pexp, t = t_vec, z = z_vec, rate = delta/v)


## marginal difference plots

plot_mtc(mtc_id, xlab = expression(alpha), ylab = "Sup. difference")
ggsave("plots/mtc_sup_1e6_id.png", width = 3.2, height = 3)


################################################################################
### probabilistic tail calibration

ptc_cl <- tail_prob_cal(y, pgpd, t = rd_vec, xi = gamma)
ptc_id <- tail_prob_cal(y, pexp, t = rd_vec, rate = delta)
ptc_ex <- tail_prob_cal(y, pexp, t = rd_vec, rate = delta/v)


## cpit reliability diagrams

plot_ptc(ptc_id, names = rd_q)
ggsave("plots/ptc_rh_1e6_id.png", width = 3.5, height = 3.4)


################################################################################
### conditional probabilistic tail calibration

ptc_cl <- tail_prob_cal(y, pgpd, t = t_vec, xi = gamma)
ptc_id <- tail_prob_cal(y, pexp, t = t_vec, rate = delta)
ptc_ex <- tail_prob_cal(y, pexp, t = t_vec, rate = delta/v)


## cpit divergence plot

cond_plots <- plot_ptc_cond(list(cl = ptc_cl, id = ptc_id, ex = ptc_ex),
                            names = c("Clim.", "Ideal", "Extr."),
                            t = t_vec, y = y, grp = delta, n_grp = 4, ylim = c(-0.01, 0.1))
#ggsave(plot = cond_plots, "plots/ptc_div_1e6_cond.png", width = 10, height = 2.8)


################################################################################
### combined ratio

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

plot_comb_ratio <- function(comb, t, u = seq(0, 1, 0.01), ylims = c(0, 1), title = NULL) {
  df <- data.frame(u = u,
                   r = as.vector(sapply(comb, function(i) i$cpit)),
                   mth = as.factor(rep(t, each = length(u))))
  comb_plot <- ggplot(df) + geom_line(aes(x = u, y = r, col = mth)) +
    geom_abline(aes(intercept = 0, slope = 1), lty = "dotted") +
    scale_x_continuous(name = "u", expand = c(0, 0), limits = c(0, 1)) +
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


comb_cl <- get_comb_ratio(y, pgpd, t = rd_vec, xi = gamma)
comb_id <- get_comb_ratio(y, pexp, t = rd_vec, rate = delta)
comb_ex <- get_comb_ratio(y, pexp, t = rd_vec, rate = delta/v)

## combined diagnostic plot
plot_comb_ratio(comb_ex, t = rd_q, ylims = c(0, 1.02), title = "Extremist")
#ggsave("plots/ptc_cmb_1e6_id.png", width = 3.5, height = 3.4)


################################################################################
### combined ratio div plots

plot_comb_div <- function(div_list, t, u = seq(0, 1, 0.01), names = NULL, ylims = NULL,
                          xlims = range(u), ylab = NULL, xlab = "t", title = NULL) {

  if (is.null(names)) names <- names(div_list)

  df <- data.frame(t = t, r = unlist(div_list), mth = rep(names, each = length(t)))
  div_plot <- ggplot(df) + geom_line(aes(x = t, y = r, col = mth)) +
    geom_hline(aes(yintercept = 0), linetype = "dotted") +
    scale_x_continuous(name = xlab, limits = xlims, expand = c(0, 0)) +
    scale_y_continuous(name = ylab, limits = ylims) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.justification = c(0, 1),
          legend.position = c(0.01, 0.99),
          plot.margin = margin(c(5.5, 10.5, 5.5, 5.5), "points")) +
    ggtitle(title)
  return(div_plot)
}

plot_comb_ratio_div <- function(y, gamma, rate1, rate2, t, grp, n_grp = 1, u = seq(0, 1, 0.01), names = NULL, ylims = NULL, xlims = NULL,
                                ylab = NULL, xlab = NULL, title = NULL, quantile_sc = T) {

  quants <- seq(0, 1, length.out = n_grp + 1)

  cond_plots <- vector("list", n_grp)
  for (i_g in 1:n_grp) {

    # get category limits
    grp_ab <- quantile(grp, quants[c(i_g, i_g + 1)])

    ind <- (grp > grp_ab[1]) & (grp <= grp_ab[2])

    # get calibration
    comb_cl <- get_comb_ratio(y[ind], pgpd, t = t, xi = gamma)
    comb_id <- get_comb_ratio(y[ind], pexp, t = t, rate = rate1[ind])
    comb_ex <- get_comb_ratio(y[ind], pexp, t = t, rate = rate2[ind])

    div_cl <- sapply(comb_cl, function(x) max(abs(x$cpit - u)))
    div_id <- sapply(comb_id, function(x) max(abs(x$cpit - u)))
    div_ex <- sapply(comb_ex, function(x) max(abs(x$cpit - u)))

    # convert thresholds to quantile scale
    if (quantile_sc) {
      q_vec <- sapply(t, function(tt) mean(y[(grp > grp_ab[1]) & (grp <= grp_ab[2])] <= tt))
    } else {
      q_vec <- t
    }

    # plot calibration
    cond_plots[[i_g]] <- plot_comb_div(list(cl = div_cl, id = div_id, ex = div_ex),
                                       q_vec, u, names, ylims, xlims, ylab, xlab, title)
  }

  return(cond_plots)
}



## combined divergence diagnostic plot
comb_div_plots <- plot_comb_ratio_div(y, gamma, delta, delta/v, t = t_vec, grp = delta, n_grp = 3,
                                      ylims = c(0, 2.2), xlab = "t", names = c("Clim.", "Ideal", "Extr."))
comb_div_plots[[2]] + ggtitle(expression(B[2]))
ggsave("plots/ptc_cmb_1e6_div2.png", width = 3, height = 3)



