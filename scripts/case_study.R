################################################################################
## set up

set.seed(91)

library(ncdf4)
library(ggplot2)
library(zoo)
library(crch)
library(scoringRules)
library(WeightedForecastVerification)
library(reliabilitydiag)


################################################################################
## load data

path <- "C:/Users/sa20i493/Documents/Data/EUMetNet/tp6_station_"
load_data <- function(path, na_prop = 10) {

  ## training data

  ## fcst
  fcst_file <- nc_open(paste0(path, "refo_fc.ncdf4"))
  train_stid_fc <- ncvar_get(fcst_file, varid = "station_id")
  train_year_fc <- ncvar_get(fcst_file, varid = "year")
  #train_time_fc <- as.POSIXct(ncvar_get(fcst_file, varid = "time")*(24*60*60), origin = '2017-01-02')
  train_time_fc <- ncvar_get(fcst_file, varid = "time")
  train_lt_fc <- ncvar_get(fcst_file, varid = "step")
  train_ens_fc <- ncvar_get(fcst_file, varid = "number")
  train_fc <<- ncvar_get(fcst_file, "tp6")*1000 # convert from m to mm
  train_fc <<- aperm(train_fc, c(5, 1, 2, 4, 3))
  nc_close(fcst_file)

  ## obs
  obs_file <- nc_open(paste0(path, "refo_obs.ncdf4"))
  train_stid_obs <- ncvar_get(obs_file, varid = "station_id")
  train_year_obs <- ncvar_get(obs_file, varid = "year")
  #train_time_obs <- as.POSIXct(ncvar_get(obs_file, varid = "time")*(24*60*60), origin = '2017-01-02')
  train_time_obs <- ncvar_get(obs_file, varid = "time")
  train_lt_obs <- ncvar_get(obs_file, varid = "step")
  train_obs <<- ncvar_get(obs_file, "tp6")*1000 # convert from m to mm
  lonlatalt <<- data.frame(lon = ncvar_get(obs_file, "longitude"),
                           lat = ncvar_get(obs_file, "latitude"),
                           alt = ncvar_get(obs_file, "altitude"),
                           id = ncvar_get(obs_file, "station_id"))
  nc_close(obs_file)



  ### test data

  ## fcst
  fcst_file <- nc_open(paste0(path, "1718_fc.ncdf4"))
  test_stid_fc <- ncvar_get(fcst_file, varid = "station_id")
  test_time_fc <- as.POSIXct(ncvar_get(fcst_file, varid = "time"), origin = '1970-01-01')
  test_lt_fc <- ncvar_get(fcst_file, varid = "step")
  test_ens_fc <- ncvar_get(fcst_file, varid = "number")
  test_fc <<- ncvar_get(fcst_file, varid = "tp6")*1000 # convert from m to mm
  test_fc <<- aperm(test_fc, c(4, 1, 2, 3))
  nc_close(fcst_file)

  ## obs
  obs_file <- nc_open(paste0(path, "1718_obs.ncdf4"))
  test_stid_obs <- ncvar_get(obs_file, varid = "station_id")
  test_time_obs <- as.POSIXct(ncvar_get(obs_file, varid = "time"), origin = '1970-01-01')
  test_lt_obs <- ncvar_get(obs_file, varid = "step")
  test_obs <<- ncvar_get(obs_file, varid = "tp6")*1000 # convert from m to mm
  nc_close(obs_file)


  ### checks
  if (identical(train_lt_fc, train_lt_obs) &
      identical(test_lt_fc, test_lt_obs) &
      identical(train_lt_fc, test_lt_fc)) {
    lead_times <<- test_lt_obs
  } else {
    stop("Lead times in training and test data do not match")
  }
  if (identical(train_stid_fc, train_stid_obs) &
      identical(test_stid_fc, test_stid_obs) &
      identical(train_stid_fc, test_stid_fc)) {
    stat_ids <<- test_stid_obs
  } else {
    stop("Station IDs in training and test data do not match")
  }
  if (identical(train_time_fc, train_time_obs) &
      identical(test_time_fc, test_time_obs)) {
    train_times <<- train_time_obs
    test_times <<- test_time_obs
  } else {
    stop("Forecast reference times in training and test data do not match")
  }
  if (identical(train_year_fc, train_year_obs)) {
    train_years <<- train_year_obs
  } else {
    stop("Years in training data do not match")
  }
  train_n_ens <<- length(train_ens_fc)
  test_n_ens <<- length(test_ens_fc)


  train_obs[is.nan(train_obs)] <<- NA
  test_obs[is.nan(test_obs)] <<- NA
  train_obs[train_obs < 0] <<- 0
  test_obs[test_obs < 0] <<- 0
  train_fc[train_fc < 0] <<- 0
  test_fc[test_fc < 0] <<- 0
  train_obs[train_obs > 100] <<- NA # needs choosing properly

  ### remove stations with missing data
  i <- 1
  while (i <= length(stat_ids)) {
    id <- stat_ids[i]
    train_na <- sapply(1:length(lead_times), function(lt) mean(is.na(train_obs[i, lt, , ])))
    test_na <- sapply(1:length(lead_times), function(lt) mean(is.na(test_obs[i, lt, ])))
    if (any(100*train_na > na_prop) || any(100*test_na > na_prop)) {
      train_obs <<- train_obs[-i, , , ]
      train_fc <<- train_fc[-i, , , , ]
      test_obs <<- test_obs[-i, , ]
      test_fc <<- test_fc[-i, , , ]
      stat_ids <<- stat_ids[-i]
      lonlatalt <<- lonlatalt[-i, ]
      print(paste("Station", id, "has been removed due to a high proportion of missing values"))
    } else {
      i <- i + 1
    }
  }



}
load_data(path)

plot_example <- function(test = sample(c(T, F), 1)) {
  s <- sample(seq_along(stat_ids), 1)
  if (test) {
    n_ens <- test_n_ens
    t <- sample(seq_along(test_times), 1)
    df <- data.frame(lt = lead_times,
                     y = c(as.vector(test_fc[s, , t, ]), test_obs[s, , t]),
                     m = as.factor(rep(0:n_ens, each = length(lead_times))))
  } else {
    n_ens <- train_n_ens
    t <- sample(seq_along(train_times), 1)
    y <- sample(seq_along(train_years), 1)
    df <- data.frame(lt = lead_times,
                     y = c(as.vector(train_fc[s, , y, t, ]), train_obs[s, , y, t]),
                     m = as.factor(rep(0:n_ens, each = length(lead_times))))
  }

  ggplot(df) + geom_line(aes(x = lt, y = y, col = m)) +
    scale_x_continuous(name = "Lead time (hours)", expand = c(0, 0)) +
    scale_y_continuous(name = "Precip. (mm)") +
    scale_color_manual(values = c(rep("grey", n_ens), "black")) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "none")
}
plot_example(F)

plot_map <- function(lons, lats, z, filename = NULL, ymin = 0, ymax = 15, title = NULL){
  if (is.matrix(z)) z <- as.vector(z)

  world <- map_data("world")
  ind <- (world$long >= min(lons) & world$long < max(lons)) & (world$lat >= min(lats) & world$lat <= max(lats))
  world <- world[ind, ]

  df <- data.frame(lat = lats, lon = lons, z = z)
  plot_obj <- ggplot() + geom_point(data = df, aes(lon, lat, fill = z), shape = 21, size = 2) +
    borders("world") +
    coord_fixed(ylim = range(lats), xlim = range(lons)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_gradient(low = "blue", high = "red", limits = c(ymin, ymax),
                        guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
    theme_void() + theme(legend.title = element_blank(), legend.position = "bottom",
                         panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
                         legend.key.width = unit(0.3, "in")) +
    ggtitle(title)


  if (!is.null(filename)) {
    ggsave(filename, plot_obj, width = 2, height = 2.5)
  }

  return(plot_obj)

}
plot_map(lonlatalt$lon, lonlatalt$lat, lonlatalt$alt, ymin = -5, ymax = 1600, filename = "altplot.png")


################################################################################
## post-processing

y_vec <- seq(0.1, 20, 0.1)
t_vec <- c(5, 10, 15)

pit <- crps <- list(raw = array(NA, dim(test_obs)),
                    emos = array(NA, dim(test_obs)),
                    raw_sm = array(NA, dim(test_obs)))
margd <- cmargd <- list(raw = array(NA, c(dim(test_obs)[1:2], length(y_vec))),
                        emos = array(NA, c(dim(test_obs)[1:2], length(y_vec))),
                        raw_sm = array(NA, c(dim(test_obs)[1:2], length(y_vec))))
cpit <- exc_p <- list(raw = array(NA, c(dim(test_obs), length(t_vec))),
                      emos = array(NA, c(dim(test_obs), length(t_vec))),
                      raw_sm = array(NA, c(dim(test_obs), length(t_vec))))

get_test_data <- function(test_obs, test_fc, test_times, stat_ids, lt, st) {

  obs <- test_obs[which(stat_ids == st), which(lead_times == lt), ]
  time <- test_times + lt * 3600
  obs <- cbind(test_times, lt, st, obs)
  obs <- zoo(obs, time)
  colnames(obs) <- c('init', 'lt', 'stat', 'obs')

  fc <- test_fc[which(stat_ids == st), which(lead_times == lt), , ]
  ens.mu <- apply(fc, 1, mean)
  ens.sd <- apply(fc, 1, sd)
  test <- cbind(fc, ens.mu, ens.sd)
  test <- zoo(test, time)
  test <- merge(obs, test)

  #test <- na.omit(test)
  yday <- as.POSIXlt(index(test))$yday
  test$sin.y1 <- sin(2 * pi * yday / 365)
  test$cos.y1 <- cos(2 * pi * yday / 365)
  test$sin.y2 <- sin(4 * pi * yday / 365)
  test$cos.y2 <- cos(4 * pi * yday / 365)

  return(test)
}

get_train_data <- function(train_obs, train_fc, train_times, train_years, stat_ids, lt, st) {

  trai.all <- NULL
  for (year in train_years) {
    obs <- train_obs[which(stat_ids == st), which(lead_times == lt), which(train_years == year), ]
    year.dummy <- as.POSIXct(paste0(2017 - 20 + year, '-01-02'))
    time <- year.dummy + train_times * 3600 * 24 + lt * 3600
    init_time <- year.dummy + train_times * 3600 * 24
    obs <- cbind(init_time, lt, st, obs)
    obs <- zoo(obs, time)
    colnames(obs) <- c('init', 'lt', 'stat', 'obs')

    fc <- train_fc[which(stat_ids == st), which(lead_times == lt), which(train_years == year), , ]
    ens.mu <- apply(fc, 1, mean)
    ens.sd <- apply(fc, 1, sd)
    trai <- cbind(ens.mu, ens.sd)
    trai <- zoo(trai, time)
    trai <- merge(obs, trai)

    trai.all <- rbind(trai.all, trai)
  }

  trai <- na.omit(trai.all)
  trai <- trai[which(index(trai) < as.POSIXct('2017-01-01')),]

  yday <- as.POSIXlt(index(trai))$yday
  trai$sin.y1 <- sin(2 * pi * yday / 365)
  trai$cos.y1 <- cos(2 * pi * yday / 365)
  trai$sin.y2 <- sin(4 * pi * yday / 365)
  trai$cos.y2 <- cos(4 * pi * yday / 365)

  return(trai)
}


fit_emos <- function(trai, test, sd_aug = 0.01, dist = "logistic") {
  pred <- tryCatch(
    {
      #fit <- crch(obs ~ ens.mu + sin.y1 + cos.y1 + sin.y2 + cos.y2 |
      #       ens.sd + sin.y1 + cos.y1 + sin.y2 + cos.y2,
      #     data = trai,
      #     dist = dist,
      #     left = 0,
      #     type = "crps")

      fit <- crch(obs ~ ens.mu | ens.sd,
           data = trai,
           link.scale = "quadratic",
           dist = dist,
           left = 0,
           type = "crps")

      mu <- predict(fit, newdata = test)
      sig <- predict(fit, newdata = test, type = 'scale')
      pred <- cbind(as.numeric(mu), as.numeric(sig))
      colnames(pred) <- c("mu", "sig")
      pred
    },
    error = function(cond) {
      message("Fit failed, returning ensemble mean and standard deviation")
      pred <- as.matrix(subset(test, select = c('ens.mu', 'ens.sd')))
      colnames(pred) <- c("mu", "sig")
      rownames(pred) <- NULL
      pred[pred[, 2] == 0, 2] <- sd_aug
      return(pred)
      })

  return(pred)
}

fit_emos <- function(trai, test) {

  loss <- function(par) {
    mu <- par[1] + par[2]*trai$ens.mu
    sig <- (par[3]^2) + (par[4]^2)*trai$ens.sd
    s <- crps_clogis(trai$obs, mu, sig, lower = 0)
    return(mean(s))
  }
  par <- optim(par = c(0, 1, 1, 1), loss)$par

  mu <- par[1] + par[2]*test$ens.mu
  sig <- (par[3]^2) + (par[4]^2)*test$ens.sd
  pred <- cbind(mu, sig)
  return(pred)
}

eval_emos <- function(preds, obs, t_vec, y_vec, dist = "logistic") {
  mu <- preds[, 1]
  sig <- preds[, 2]

  if (dist == "normal") {
    F_x <- pnorm
    crps <- crps_cnorm(as.numeric(obs), mu, sig, lower = 0)
  } else if (dist == "logistic") {
    F_x <- plogis
    crps <- crps_clogis(as.numeric(obs), mu, sig, lower = 0)
  }

  pit <- F_x(obs, mu, sig)
  zero <- obs == 0 & !is.na(obs)
  pit[zero] <- runif(sum(zero), 0, F_x(obs[zero], mu[zero], sig[zero]))
  cpit <- sapply(t_vec, function(t) {
    z <- (F_x(obs, mu, sig) - F_x(t, mu, sig))/(1 - F_x(t, mu, sig))
    z[obs < t] <- NA
    return(z)
    })
  margd <- sapply(y_vec, function(y) mean(F_x(y, mu, sig)) - mean(obs <= y, na.rm = T))
  cmargd <- sapply(y_vec, function(t) { # note y_vec instead of t_vec since we have a plot against threshold
    dd <- sapply(y_vec, function(y) {
      if (sum(obs > y, na.rm = T) == 0) {
        d <- NA
      } else {
        ind <- obs > t
        F_hat <- mean((F_x(y + t, mu[ind], sig[ind]) - F_x(t, mu[ind], sig[ind]))/(1 - F_x(t, mu[ind], sig[ind])), na.rm = T)
        Q_hat <- mean(obs[ind] - t <= y, na.rm = T)
        d <- F_hat - Q_hat
      }
      return(d)
    })
    return(max(abs(dd)))
  })
  exc_p <- sapply(t_vec, function(t) 1 - F_x(t, mu, sig))
  return(list(pit = pit, crps = crps, cpit = cpit, margd = margd, cmargd = cmargd, exc_p = exc_p))
}

eval_ens <- function(pred, obs, t_vec, y_vec, dist = "empirical") {

  mu <- rowMeans(pred)
  sig <- apply(pred, 1, sd)

  if (dist == "normal") {
    sig[sig == 0] <- 0.1
    F_x <- pnorm
    crps <- crps_cnorm(as.numeric(obs), mu, sig, lower = 0)
  } else if (dist == "logistic") {
    sig[sig == 0] <- 0.1
    F_x <- plogis
    crps <- crps_clogis(as.numeric(obs), mu, sig, lower = 0)
  } else {
    F_x <- function(x, pred) {
      if (is.vector(pred)) {
        mean(pred <= x)
      } else {
        rowMeans(pred <= x)
      }
    }
  }


  if (dist %in% c("normal", "logistic")) {

    pit <- F_x(obs, mu, sig)
    zero <- obs == 0 & !is.na(obs)
    pit[zero] <- runif(sum(zero), 0, F_x(obs[zero], mu[zero], sig[zero]))

    cpit <- sapply(t_vec, function(t) {
      z <- (F_x(obs, mu, sig) - F_x(t, mu, sig))/(1 - F_x(t, mu, sig))
      z[obs < t] <- NA
      return(z)
    })

    margd <- sapply(y_vec, function(y) mean(F_x(y, mu, sig)) - mean(obs <= y, na.rm = T))

    cmargd <- sapply(y_vec, function(t) { # note y_vec instead of t_vec since we have a plot against threshold
      dd <- sapply(y_vec, function(y) {
        if (sum(obs > y, na.rm = T) == 0) {
          d <- NA
        } else {
          ind <- obs > t
          F_hat <- mean((F_x(y + t, mu[ind], sig[ind]) - F_x(t, mu[ind], sig[ind]))/(1 - F_x(t, mu[ind], sig[ind])), na.rm = T)
          Q_hat <- mean(obs[ind] - t <= y, na.rm = T)
          d <- F_hat - Q_hat
        }
        return(d)
      })
      return(max(abs(dd)))
    })

    exc_p <- sapply(t_vec, function(t) 1 - F_x(t, mu, sig))

  } else {

    S <- cbind(obs, pred)
    pit <- apply(S, 1, rank, ties.method = "random")[1, ]

    cpit <- sapply(t_vec, function(t) {
      z <- (F_x(obs, pred) - F_x(t, pred))/(1 - F_x(t, pred))
      z[F_x(t, pred) == 1] <- 1
      z[obs < t] <- NA
      return(z)
    })

    margd <- sapply(y_vec, function(y) mean(F_x(y, pred)) - mean(obs <= y, na.rm = T))

    cmargd <- sapply(y_vec, function(t) { # note y_vec instead of t_vec since we have a plot against threshold
      dd <- sapply(y_vec, function(y) {
        if (sum(obs > y, na.rm = T) == 0) {
          d <- NA
        } else {
          ind <- obs > t
          F_hat <- mean((F_x(y + t, pred[ind, ]) - F_x(t, pred[ind, ]))/(1 - F_x(t, pred[ind, ])), na.rm = T)
          Q_hat <- mean(obs[ind] - t <= y, na.rm = T)
          d <- F_hat - Q_hat
        }
        return(d)
      })
      return(max(abs(dd)))
    })

    exc_p <- sapply(t_vec, function(t) 1 - F_x(t, pred))

    crps <- rep(NA, length(obs))
    na_ind <- is.na(obs)
    crps[!na_ind] <- crps_sample(y = obs[!na_ind], dat = pred[!na_ind, ])
  }

  return(list(pit = pit, crps = crps, cpit = cpit, margd = margd, cmargd = cmargd, exc_p = exc_p))
}

emos_all_preds <- array(NA, dim = c(730, 2, 112))
for (j in seq_along(stat_ids)) {
  st <- stat_ids[j]
  for (i in seq_along(lead_times)[4]) { # 4 restricts attention to 24 hours
    lt <- lead_times[i]
    print(paste0('Forecast at Station: ', st, ' (', j, ' from ', length(stat_ids), ') and Lead Time: ', lt))

    ### Get train data
    trai <- get_train_data(train_obs, train_fc, train_times, train_years, stat_ids, lt, st)

    ### Get test data
    test <- get_test_data(test_obs, test_fc, test_times, stat_ids, lt, st)

    ### Fit EMOS
    emos_preds <- fit_emos(trai, test)
    emos_all_preds[, , j] <- emos_preds

    ### Eval EMOS
    scores <- eval_emos(emos_preds, as.numeric(test$obs), t_vec, y_vec)
    pit[['emos']][j, i, ] <- scores$pit
    crps[['emos']][j, i, ] <- scores$crps
    cpit[['emos']][j, i, , ] <- scores$cpit
    margd[['emos']][j, i, ] <- scores$margd
    cmargd[['emos']][j, i, ] <- scores$cmargd
    exc_p[['emos']][j, i, , ] <- scores$exc_p

    ### Eval Raw
    ens <- as.matrix(test[, names(test) %in% as.character(1:51)])
    scores <- eval_ens(ens, as.numeric(test$obs), t_vec, y_vec, dist = "empirical")
    pit[['raw']][j, i, ] <- scores$pit
    crps[['raw']][j, i, ] <- scores$crps
    cpit[['raw']][j, i, , ] <- scores$cpit
    margd[['raw']][j, i, ] <- scores$margd
    cmargd[['raw']][j, i, ] <- scores$cmargd
    exc_p[['raw']][j, i, , ] <- scores$exc_p

    ### Eval smoothed Raw
    ens <- as.matrix(test[, names(test) %in% as.character(1:51)])
    scores <- eval_ens(ens, as.numeric(test$obs), t_vec, y_vec, dist = "logistic")
    pit[['raw_sm']][j, i, ] <- scores$pit
    crps[['raw_sm']][j, i, ] <- scores$crps
    cpit[['raw_sm']][j, i, , ] <- scores$cpit
    margd[['raw_sm']][j, i, ] <- scores$margd
    cmargd[['raw_sm']][j, i, ] <- scores$cmargd
    exc_p[['raw_sm']][j, i, , ] <- scores$exc_p

  }
}


#save(pit, crps, margd, cpit, cmargd, exc_p, file = "scripts/cs_all_scores.RData")
load("scripts/cs_all_scores.RData")

################################################################################
## evaluation


lt <- which(lead_times == 24)

## crps
1 - mean(crps[['emos']][, lt, ], na.rm = T)/mean(crps[['raw']][, lt, ], na.rm = T)
1 - mean(crps[['raw_sm']][, lt, ], na.rm = T)/mean(crps[['raw']][, lt, ], na.rm = T)


## reliability diagrams

rel_diag <- function(y, exc_p, t, names = NULL, ...) {
  d <- lapply(seq_along(t), function(i) {
    o <- as.numeric(as.vector(y > t[i]))
    na_ind <- is.na(o)
    o <- o[!na_ind]
    p <- as.vector(exc_p[, , i])[!na_ind]
    out <- reliabilitydiag::reliabilitydiag(p, y = o)[[1]]$bins
    if (!is.null(names)) {
      name <- names[i]
    } else {
      name <- round(t[i], 2)
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

mth <- 'emos'
rel_diag(test_obs[, lt, ], exc_p[[mth]][, lt, , ], t = t_vec)
ggsave(paste0("plots/cs_rd_", mth, ".png"), width = 3.3, height = 3)


## pit histograms
mth <- 'raw_sm'
if (mth == 'raw') {
  pl_pit1 <- pit_hist(pit[[mth]][, lt, ], ymax = 0.2, xticks = F, title = "PIT")
} else {
  pl_pit1 <- pit_hist(pit[[mth]][, lt, ], bins = test_n_ens + 1, ranks = F, ymax = 0.2, xticks = F, title = "PIT")
}
pl_pit2 <- pit_hist(cpit[[mth]][, lt,  , which(t_vec == 5)],
                    bins = 20, ranks = F, ymax = 0.6, xlab = NULL, xticks = F, title = "cPIT: t = 5mm")
pl_pit3 <- pit_hist(cpit[[mth]][, lt,  , which(t_vec == 10)],
                    bins = 20, ranks = F, ymax = 0.6, xlab = NULL, xticks = F, title = "cPIT: t = 10mm")
pl_pit4 <- pit_hist(cpit[[mth]][, lt,  , which(t_vec == 15)],
                    bins = 20, ranks = F, ymax = 0.6, xlab = NULL, xticks = F, title = "cPIT: t = 15mm")
pl_pit <- gridExtra::grid.arrange(pl_pit1, pl_pit2, pl_pit3, pl_pit4, ncol = 1)
ggsave(plot = pl_pit, paste0("plots/cs_pit_", mth, ".png"), width = 2.5, height = 9)


## pit reliability diagrams
mth <- 'raw_sm'

z <- list("  S" = pit[[mth]][, lt, ],
          " 5" = cpit[[mth]][, lt, , which(t_vec == 5)],
          "10" = cpit[[mth]][, lt, , which(t_vec == 10)],
          "15" = cpit[[mth]][, lt, , which(t_vec == 15)])

if (mth == 'raw') {
  pit_reldiag(z, ranks = c(T, F, F, F), resampling = F)
} else {
  pit_reldiag(z, ranks = F, resampling = F)
}
ggsave(paste0("plots/cs_pit_rd_", mth, ".png"), width = 3.5, height = 3.4)


## combined ratio
mth <- 'raw_sm'

u <- seq(0.01, 0.99, 0.01)
df <- sapply(u, function(u) sapply(t_vec, function(t) {
  Gbar <- mean(test_obs[, lt, ] > t, na.rm = T)
  Fbar <- mean(exc_p_Fx[, which(y_vec == t)])
  Z_F <- mean(cpit[[mth]][, lt, , which(t_vec == t)] <= u, na.rm=T)
  return(Z_F*Gbar/Fbar)
}))
df <- as.data.frame(t(df))
colnames(df) <- t_vec

if (mth == "raw") {
  z <- pit[[mth]][, lt, ]
  z <- (z + runif(length(z)) - 1)/max(z) # convert ranks to PIT values
  df$S <- sapply(u, function(uu) mean(z <= uu, na.rm=T))
} else {
  df$S <- sapply(u, function(uu) mean(pit[[mth]][, lt, ] <= uu, na.rm=T))
}


df <- data.frame(u = u,
                 r = c(df$S, df$`5`, df$`10`, df$`15`),
                 mth = rep(c("  S", " 5", "10", "15"), each = length(u)))
ggplot(df) + geom_line(aes(x = u, y = r, col = mth)) +
  geom_abline(aes(intercept = 0, slope = 1), lty = "dotted") +
  scale_x_continuous(name = "u", expand = c(0, 0), limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(name = "Combined ratio", limits = c(0, 1.6), expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(0, 1),
        legend.position = c(0.01, 0.99),
        plot.margin = margin(c(5.5, 10.5, 5.5, 5.5), "points")) +
  ggtitle("Smoothed IFS")

ggsave(paste0("plots/cs_cmb_rd_", mth, ".png"), width = 3.5, height = 3.4)


## marginal calibration plots
mth <- 'raw_sm'

# recover predicted exceedance probabilities
exc_p_Fx <- 1 - sapply(seq_along(y_vec), function(j) sapply(seq_along(stat_ids), function(i)
  margd[[mth]][i, lt, j] + mean(test_obs[i, lt, ] <= y_vec[j], na.rm = T)))

df <- data.frame(rat = sapply(y_vec, function(y) mean(test_obs[, lt, ] > y, na.rm = T))/colMeans(exc_p_Fx), t = y_vec)

ggplot(df) + geom_line(aes(x = t, y = rat)) +
  geom_hline(aes(yintercept = 1), linetype = "dotted") +
  scale_x_continuous(name = "t (mm)", limits = c(0, 20), expand = c(0, 0)) +
  scale_y_continuous(name = "Occurrence ratio", limits = c(0, 6)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ggtitle(" ")
ggsave(paste0("plots/cs_marg_rat_", mth, ".png"), width = 3.5, height = 3.4)


## marginal tail calibration plots
mth <- 'raw'

df <- data.frame(cal = c(apply(cmargd[[mth]][, lt, ], 2, mean, na.rm = T)),
                 t = y_vec)

#a_vec <- sapply(y_vec, function(x) mean(test_obs <= x, na.rm = T))

ggplot(df) + geom_line(aes(x = t, y = cal)) +
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  scale_x_continuous(name = "Threshold (mm)", expand = c(0, 0)) +
  scale_y_continuous(name = "Sup. difference", limits = c(-0.1, 1)) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(paste0("plots/cs_marg_sup_", mth, ".png"), width = 3.2, height = 3)




## marginal tail calibration plots (old)
# mth <- 'emos'
#
# z <- list("  S" = colMeans(margd[[mth]][, lt, ], na.rm = T),
#           " 5" = colMeans(cmargd[[mth]][, lt, , which(t_vec == 5)], na.rm = T),
#           "10" = colMeans(cmargd[[mth]][, lt, , which(t_vec == 10)], na.rm = T),
#           "15" = colMeans(cmargd[[mth]][, lt, , which(t_vec == 15)], na.rm = T))
#
# df <- data.frame(d = unlist(z),
#                  y = y_vec,
#                  m = rep(names(z), each = length(y_vec)))
# ggplot(df) + geom_line(aes(x = y, y = d, col = as.factor(m))) +
#   geom_hline(aes(yintercept = 0), lty = "dotted") +
#   scale_x_continuous(name = "y", limits = c(0, max(y_vec)), expand = c(0, 0)) +
#   scale_y_continuous(name = expression("E[" ~ F[t] ~ "(y)] - Q(Y - t ≤ y | Y > t)"),
#                      limits = c(-0.2, 0.4), expand = c(0, 0)) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         legend.title = element_blank(),
#         legend.justification = c(0.5, 1),
#         legend.position = c(0.5, 0.99)) +
#   guides(col = guide_legend(nrow = 1))
# ggsave(paste0("plots/cs_marg_", mth, ".png"), width = 3.3, height = 3)


