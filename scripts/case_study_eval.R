################################################################################
# script to evaluate the forecasts in the case study

################################################################################
## set up

set.seed(69249837)

library(TailCalibration)
library(crch)
library(evd)
library(ggplot2)

# load forecast data
fc_dat <- readRDS(here::here("scripts", "cs_fc_data.RDS"))


# remove missing data
y <- fc_dat$aux_data$obs
na_ind <- is.na(y)
y <- y[!na_ind]


# specify thresholds at which to evaluate tail calibration
y_vec <- seq(0, 20, 0.1)
t_vec <- c(-Inf, 5, 10, 15)


# specify plot dimensions
width <- 3.5
height <- 3.4


################################################################################
## evaluation


##### IFS

dat <- fc_dat$ifs$dat
dat <- matrix(dat, ncol = 51)[!na_ind, ]

## combined ratio
com <- tc_prob(y, dat, t = t_vec, lower = 0)
plot_ptc(com, names = c("  S", " 5", "10", "15"), ylims = c(0, 1.9), title = "IFS")
ggsave("plots/cs_com_ifs.png", width = width, height = height)


## severity ratio
sev <- tc_prob(y, dat, t = t_vec, ratio = "sev", lower = 0)
plot_ptc(sev, ratio = "sev", names = c("  S", " 5", "10", "15"), ylims = c(0, 1.02), title = "")
ggsave("plots/cs_sev_ifs.png", width = width, height = height)


## occurrence
occ <- tc_prob(y, dat, t = y_vec, ratio = "occ", lower = 0)
plot_ptc(occ, ratio = "occ", ylims = c(-0.2, 6), xlab = "t (mm)", title = "")
ggsave("plots/cs_occ_ifs.png", width = width, height = height)


## testing
sev <- tc_prob(y, dat, t = t_vec, ratio = "sev", lower = 0, test = TRUE)
occ <- tc_prob(y, dat, t = t_vec, ratio = "occ", lower = 0, test = TRUE)



##### smoothed IFS

F_x <- fc_dat$smooth$F_x
mu <- as.vector(fc_dat$smooth$location)[!na_ind]
sig <- as.vector(fc_dat$smooth$scale)[!na_ind]


## combined ratio
com <- tc_prob(y, F_x, t = t_vec, lower = 0, location = mu, scale = sig)
plot_ptc(com, names = c("  S", " 5", "10", "15"), ylims = c(0, 1.9), title = "Smoothed IFS")
ggsave("plots/cs_com_smo.png", width = width, height = height)


## severity ratio
sev <- tc_prob(y, F_x, t = t_vec, ratio = "sev", lower = 0, location = mu, scale = sig)
plot_ptc(sev, ratio = "sev", names = c("  S", " 5", "10", "15"), ylims = c(0, 1.02), title = "")
ggsave("plots/cs_sev_smo.png", width = width, height = height)


## occurrence
occ <- tc_prob(y, F_x, t = y_vec, ratio = "occ", lower = 0, location = mu, scale = sig)
plot_ptc(occ, ratio = "occ", names = c("  S", " 5", "10", "15"), ylims = c(-0.2, 6), xlab = "t (mm)", title = "")
ggsave("plots/cs_occ_smo.png", width = width, height = height)


## testing
sev <- tc_prob(y, F_x, t = t_vec, lower = 0, ratio = "sev", test = TRUE, location = mu, scale = sig)
occ <- tc_prob(y, F_x, t = t_vec, ratio = "occ", lower = 0, test = TRUE, location = mu, scale = sig)



##### post-processed (logistic)

F_x <- fc_dat$emos_cl$F_x
mu <- as.vector(fc_dat$emos_cl$location)[!na_ind]
sig <- as.vector(fc_dat$emos_cl$scale)[!na_ind]


## combined ratio
com <- tc_prob(y, F_x, t = t_vec, lower = 0, location = mu, scale = sig)
plot_ptc(com, names = c("  S", " 5", "10", "15"), ylims = c(0, 1.9), title = "Post-processed (Logistic)")
ggsave("plots/cs_com_pp.png", width = width, height = height)


## severity ratio
sev <- tc_prob(y, F_x, t = t_vec, lower = 0, ratio = "sev", location = mu, scale = sig)
plot_ptc(sev, ratio = "sev", names = c("  S", " 5", "10", "15"), ylims = c(0, 1.02), title = "")
ggsave("plots/cs_sev_pp.png", width = width, height = height)


## occurrence
occ <- tc_prob(y, F_x, t = y_vec, ratio = "occ", lower = 0, location = mu, scale = sig)
plot_ptc(occ, ratio = "occ", names = c("  S", " 5", "10", "15"), ylims = c(-0.2, 6), xlab = "t (mm)", title = "")
ggsave("plots/cs_occ_pp.png", width = width, height = height)


## testing
sev <- tc_prob(y, F_x, t = t_vec, lower = 0, ratio = "sev", test = TRUE, location = mu, scale = sig)
occ <- tc_prob(y, F_x, t = t_vec, ratio = "occ", lower = 0, test = TRUE, location = mu, scale = sig)


##### post-processed (GEV)

F_x <- fc_dat$emos_cgev$F_x
mu <- as.vector(fc_dat$emos_cgev$location)[!na_ind]
sig <- as.vector(fc_dat$emos_cgev$scale)[!na_ind]
shape <- as.vector(fc_dat$emos_cgev$shape)[!na_ind]


## combined ratio
com <- tc_prob(y, F_x, t = t_vec, lower = 0, location = mu, scale = sig, shape = shape)
plot_ptc(com, names = c("  S", " 5", "10", "15"), ylims = c(0, 1.9), title = "Post-processed (GEV)")
ggsave("plots/cs_com_ppgev.png", width = width, height = height)


## severity ratio
sev <- tc_prob(y, F_x, t = t_vec, lower = 0, ratio = "sev", location = mu, scale = sig, shape = shape)
plot_ptc(sev, ratio = "sev", names = c("  S", " 5", "10", "15"), ylims = c(0, 1.02), title = "")
ggsave("plots/cs_sev_ppgev.png", width = width, height = height)


## occurrence
occ <- tc_prob(y, F_x, t = y_vec, ratio = "occ", lower = 0, location = mu, scale = sig, shape = shape)
plot_ptc(occ, ratio = "occ", names = c("  S", " 5", "10", "15"), ylims = c(-0.2, 6), xlab = "t (mm)", title = "")
ggsave("plots/cs_occ_ppgev.png", width = width, height = height)


## testing
sev <- tc_prob(y, F_x, t = t_vec, lower = 0, ratio = "sev", test = TRUE, location = mu, scale = sig, shape = shape)
occ <- tc_prob(y, F_x, t = t_vec, ratio = "occ", lower = 0, test = TRUE, location = mu, scale = sig, shape = shape)



################################################################################
## evaluation using quantile-based thresholds

q_vec <- c(0.97, 0.99, 0.995)
q_vec_oc <- seq(0.7, 0.999, 0.001)
q_name <- c(" S", as.character(q_vec))

y <- fc_dat$aux_data$obs
a_vec <- sapply(q_vec, function(q) apply(y, 1, quantile, q, na.rm=T) |> rep(ncol(y)))
a_vec_oc <- sapply(q_vec_oc, function(q) apply(y, 1, quantile, q, na.rm=T) |> rep(ncol(y)))

y <- y[!na_ind]
a_vec <- a_vec[!na_ind, ]
a_vec_oc <- a_vec_oc[!na_ind, ]

a_vec <- cbind(-Inf, a_vec)
q_vec <- c(0, q_vec)



##### IFS

dat <- fc_dat$ifs$dat
dat <- matrix(dat, ncol = 51)[!na_ind, ]


## combined ratio
com <- lapply(seq_along(q_vec), function(q) tc_prob(y, dat, t = a_vec[, q], lower = 0, var_t = TRUE))
plot_ptc(com, names = q_name, ylims = c(0, 1.9), title = "IFS")
ggsave("plots/cs_com_ifs_qu.png", width = width, height = height)


## severity ratio
sev <- lapply(seq_along(q_vec), function(q) tc_prob(y, dat, t = a_vec[, q], ratio = "sev", lower = 0, var_t = TRUE))
plot_ptc(sev, ratio = "sev", names = q_name, ylims = c(0, 1.1), title = "")
ggsave("plots/cs_sev_ifs_qu.png", width = width, height = height)


## occurrence
occ <- sapply(seq_along(q_vec_oc), function(q) tc_prob(y, dat, t = a_vec_oc[, q], ratio = "occ", lower = 0, var_t = TRUE))
occ <- data.frame(t = q_vec_oc, rat = occ)
plot_ptc(occ, ratio = "occ", ylims = c(-0.2, 6), xlab = expression(alpha), title = "")
ggsave("plots/cs_occ_ifs_qu.png", width = width, height = height)



##### smoothed IFS

F_x <- fc_dat$smooth$F_x
mu <- as.vector(fc_dat$smooth$location)[!na_ind]
sig <- as.vector(fc_dat$smooth$scale)[!na_ind]


## combined ratio
com <- lapply(seq_along(q_vec), function(q) tc_prob(y, F_x, t = a_vec[, q], lower = 0, location = mu, scale = sig, var_t = TRUE))
plot_ptc(com, names = q_name, ylims = c(0, 1.9), title = "Smoothed IFS")
ggsave("plots/cs_com_smo_qu.png", width = width, height = height)


## severity ratio
sev <- lapply(seq_along(q_vec), function(q) tc_prob(y, F_x, t = a_vec[, q], ratio = "sev", lower = 0, location = mu, scale = sig, var_t = TRUE))
plot_ptc(sev, ratio = "sev", names = q_name, ylims = c(0, 1.1), title = "")
ggsave("plots/cs_sev_smo_qu.png", width = width, height = height)


## occurrence
occ <- sapply(seq_along(q_vec_oc), function(q) tc_prob(y, F_x, t = a_vec_oc[, q], ratio = "occ", lower = 0, location = mu, scale = sig, var_t = TRUE))
occ <- data.frame(t = q_vec_oc, rat = occ)
plot_ptc(occ, ratio = "occ", ylims = c(-0.2, 6), xlab = expression(alpha), title = "")
ggsave("plots/cs_occ_smo_qu.png", width = width, height = height)



##### post-processed (logistic)

F_x <- fc_dat$emos_cl$F_x
mu <- as.vector(fc_dat$emos_cl$location)[!na_ind]
sig <- as.vector(fc_dat$emos_cl$scale)[!na_ind]


## combined ratio
com <- lapply(seq_along(q_vec), function(q) tc_prob(y, F_x, t = a_vec[, q], lower = 0, location = mu, scale = sig, var_t = TRUE))
plot_ptc(com, names = q_name, ylims = c(0, 1.9), title = "Post-processed (Logistic)")
ggsave("plots/cs_com_pp_qu.png", width = width, height = height)


## severity ratio
sev <- lapply(seq_along(q_vec), function(q) tc_prob(y, F_x, t = a_vec[, q], ratio = "sev", lower = 0, location = mu, scale = sig, var_t = TRUE))
plot_ptc(sev, ratio = "sev", names = q_name, ylims = c(0, 1.1), title = "")
ggsave("plots/cs_sev_pp_qu.png", width = width, height = height)


## occurrence
occ <- sapply(seq_along(q_vec_oc), function(q) tc_prob(y, F_x, t = a_vec_oc[, q], ratio = "occ", lower = 0, location = mu, scale = sig, var_t = TRUE))
occ <- data.frame(t = q_vec_oc, rat = occ)
plot_ptc(occ, ratio = "occ", ylims = c(-0.2, 6), xlab = expression(alpha), title = "")
ggsave("plots/cs_occ_pp_qu.png", width = width, height = height)



##### post-processed (GEV)

F_x <- fc_dat$emos_cgev$F_x
mu <- as.vector(fc_dat$emos_cgev$location)[!na_ind]
sig <- as.vector(fc_dat$emos_cgev$scale)[!na_ind]
shape <- as.vector(fc_dat$emos_cgev$shape)[!na_ind]


## combined ratio
com <- lapply(seq_along(q_vec), function(q) tc_prob(y, F_x, t = a_vec[, q], lower = 0, location = mu, scale = sig, shape = shape, var_t = TRUE))
plot_ptc(com, names = q_name, ylims = c(0, 1.9), title = "Post-processed (GEV)")
ggsave("plots/cs_com_ppgev_qu.png", width = width, height = height)


## severity ratio
sev <- lapply(seq_along(q_vec), function(q) tc_prob(y, F_x, t = a_vec[, q], ratio = "sev", lower = 0, location = mu, scale = sig, shape = shape, var_t = TRUE))
plot_ptc(sev, ratio = "sev", names = q_name, ylims = c(0, 1.1), title = "")
ggsave("plots/cs_sev_ppgev_qu.png", width = width, height = height)


## occurrence
occ <- sapply(seq_along(q_vec_oc), function(q) tc_prob(y, F_x, t = a_vec_oc[, q], ratio = "occ", lower = 0, location = mu, scale = sig, shape = shape, var_t = TRUE))
occ <- data.frame(t = q_vec_oc, rat = occ)
plot_ptc(occ, ratio = "occ", ylims = c(-0.2, 6), xlab = expression(alpha), title = "")
ggsave("plots/cs_occ_ppgev_qu.png", width = width, height = height)


