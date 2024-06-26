################################################################################
# script to evaluate the forecasts in the case study

################################################################################
## set up

set.seed(69249837)

library(TailCalibration)
library(crch)

# load forecast data
fc_dat <- readRDS("scripts/cs_fc_data.RDS")


# remove missing data
y <- as.vector(fc_dat$aux_data$obs)
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
plot_ptc(occ, ratio = "occ", names = c("  S", " 5", "10", "15"), ylims = c(-0.2, 6), xlab = "t (mm)", title = "")
ggsave("plots/cs_occ_ifs.png", width = width, height = height)



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


##### post-processed

F_x <- fc_dat$emos$F_x
mu <- as.vector(fc_dat$emos$location)[!na_ind]
sig <- as.vector(fc_dat$emos$scale)[!na_ind]


## combined ratio
com <- tc_prob(y, F_x, t = t_vec, lower = 0, location = mu, scale = sig)
plot_ptc(com, names = c("  S", " 5", "10", "15"), ylims = c(0, 1.9), title = "Post-processed")
ggsave("plots/cs_com_pp.png", width = width, height = height)


## severity ratio
sev <- tc_prob(y, F_x, t = t_vec, lower = 0, ratio = "sev", location = mu, scale = sig)
plot_ptc(sev, ratio = "sev", names = c("  S", " 5", "10", "15"), ylims = c(0, 1.02), title = "")
ggsave("plots/cs_sev_pp.png", width = width, height = height)


## occurrence
occ <- tc_prob(y, F_x, t = y_vec, ratio = "occ", lower = 0, location = mu, scale = sig)
plot_ptc(occ, ratio = "occ", names = c("  S", " 5", "10", "15"), ylims = c(-0.2, 6), xlab = "t (mm)", title = "")
ggsave("plots/cs_occ_pp.png", width = width, height = height)


