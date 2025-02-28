################################################################################
# script to extract raw data and and save forecast data for the case study

################################################################################
## set up

set.seed(91)

library(ncdf4)
library(zoo)
library(crch)
library(ggplot2)
library(evd)
library(pracma)


################################################################################
## load data

path <- here::here("data", "tp6_station_")
load_data <- function(path, na_prop = 0) {

  ## training data

  ## fcst
  fcst_file <- nc_open(paste0(path, "refo_fc.ncdf4"))
  train_stid_fc <- ncvar_get(fcst_file, varid = "station_id")
  train_year_fc <- ncvar_get(fcst_file, varid = "year")
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

  ### remove stations whose proportion of missing data is > na_prop
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
load_data(path, na_prop = 10) # remove stations with more than 10% missing values

# plot example of the ensemble forecast trajectories at a random time and station
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
    theme(panel.grid = element_blank(),
          legend.position = "none")
}
plot_example(test = F)


# plot map of stations with their corresponding altitudes
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
plot_map(lonlatalt$lon, lonlatalt$lat, lonlatalt$alt, ymin = -5, ymax = 1600)


# restrict attention to forecasts issued 24h in advance
lead <- 24
lt <- which(lead_times == lead)


################################################################################
## IFS ensemble

fc_ifs <- list(dat = test_fc[, lt, , ])


################################################################################
## smoothed IFS ensemble

# calculate ensemble mean and standard deviation
ens_mn <- apply(fc_ifs$dat, c(1, 2), mean)
ens_sd <- apply(fc_ifs$dat, c(1, 2), sd)
ens_sd[ens_sd == 0] <- 0.1 # augment ensemble standard deviation when all ensemble members are zero

F_x <- function(q, location, scale) pclogis(q, location, scale, left = 0)

fc_smooth <- list(F_x = F_x, location = ens_mn, scale = ens_sd)


################################################################################
## post-processing

# wrapper to extract training data at a particular lead time (lt) and station (st)
get_train_data <- function(train_obs, train_fc, train_times, train_years, stat_ids, lt, st) {

  trai.all <- NULL
  for (year in train_years) {
    obs <- train_obs[which(stat_ids == st), which(lead_times == lt), which(train_years == year), ]
    year.dummy <- as.POSIXct(paste0(2017 - 20 + year, '-01-02')) # first Monday in 2017 is 2nd Jan
    init_time <- year.dummy + train_times * 3600 * 24
    time <- init_time + lt * 3600
    obs <- cbind(lt, st, obs)
    obs <- zoo(obs, time)
    colnames(obs) <- c('lt', 'stat', 'obs')

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

  return(trai)
}

# wrapper to extract test data at a particular lead time (lt) and station (st)
get_test_data <- function(test_obs, test_fc, test_times, stat_ids, lt, st) {

  obs <- test_obs[which(stat_ids == st), which(lead_times == lt), ]
  obs <- cbind(lt, st, obs)
  obs <- zoo(obs, test_times)
  colnames(obs) <- c('lt', 'stat', 'obs')

  fc <- test_fc[which(stat_ids == st), which(lead_times == lt), , ]
  ens.mu <- apply(fc, 1, mean)
  ens.sd <- apply(fc, 1, sd)
  test <- cbind(ens.mu, ens.sd)
  test <- zoo(test, test_times)
  test <- merge(obs, test)

  return(test)
}

# wrapper to fit emos models using crch
fit_emos_crch <- function(trai, test, dist = "logistic") {
  fit <- tryCatch(
    {
      fit <- crch(obs ~ ens.mu | ens.sd,
                  data = trai,
                  link.scale = "quadratic",
                  dist = dist,
                  left = 0,
                  type = "crps")

    },
    error = function(cond) {
      message("Fit failed, refitting the model using maximum likelihood")
      fit <- crch(obs ~ ens.mu | ens.sd,
                  data = trai,
                  link.scale = "quadratic",
                  dist = dist,
                  left = 0)
    }
  )
  mu <- predict(fit, newdata = test)
  sig <- predict(fit, newdata = test, type = 'scale')
  pred <- cbind(as.numeric(mu), as.numeric(sig))
  colnames(pred) <- c("mu", "sig")
  pred
}

# wrapper to fit emos model with GEV
crps_cgev <- function(y, location, scale, shape, eps = 0.01){
  if (shape >= 1) {
    return(1000)
  } else if(abs(shape) > eps){
    p0 <- py <- NA*numeric(length(y))
    na_ind <- !is.na(y)
    p0[na_ind] <- pgev(0, loc = location[na_ind], scale = scale[na_ind], shape = shape)
    py[na_ind] <- pgev(y[na_ind], loc = location[na_ind], scale = scale[na_ind], shape = shape)
    p0[p0 < 1e-10] <- 1e-10
    py[py < 1e-10] <- 1e-10

    gamma_1 <- gamma_2 <- NA*numeric(length(y))
    gamma_1[na_ind] <- sapply(py[na_ind], function(x){gammainc(-log(x), 1 - shape)[1]})
    gamma_2[na_ind] <- sapply(p0[na_ind], function(x){gammainc(-2*log(x), 1 - shape)[1]})

    crps1 <- (location - y)*(1 - 2*py)
    crps2 <- location*(p0^2)
    crps3 <- 2*(scale/shape)*(1 - py - gamma_1)
    crps4 <- (scale/shape)*(1 - (p0^2) - (2^shape)*gamma_2)

    crps <- crps1 + crps2 - crps3 + crps4
  }else{
    shape_0 <- shape

    shape <- -eps
    p0 <- py <- NA*numeric(length(y))
    na_ind <- !is.na(y)
    p0[na_ind] <- pgev(0, loc = location[na_ind], scale = scale[na_ind], shape = shape)
    py[na_ind] <- pgev(y[na_ind], loc = location[na_ind], scale = scale[na_ind], shape = shape)
    p0[p0 < 1e-10] <- 1e-10
    py[py < 1e-10] <- 1e-10

    gamma_1 <- gamma_2 <- NA*numeric(length(y))
    gamma_1[na_ind] <- sapply(py[na_ind], function(x){gammainc(-log(x), 1 - shape)[1]})
    gamma_2[na_ind] <- sapply(p0[na_ind], function(x){gammainc(-2*log(x), 1 - shape)[1]})

    crps1 <- (location - y)*(1 - 2*py)
    crps2 <- location*(p0^2)
    crps3 <- 2*(scale/shape)*(1 - py - gamma_1)
    crps4 <- (scale/shape)*(1 - (p0^2) - (2^shape)*gamma_2)

    crps_neg <- crps1 + crps2 - crps3 + crps4


    shape <- eps
    p0 <- py <- NA*numeric(length(y))
    na_ind <- !is.na(y)
    p0[na_ind] <- pgev(0, loc = location[na_ind], scale = scale[na_ind], shape = shape)
    py[na_ind] <- pgev(y[na_ind], loc = location[na_ind], scale = scale[na_ind], shape = shape)
    p0[p0 < 1e-10] <- 1e-10
    py[py < 1e-10] <- 1e-10

    gamma_1 <- gamma_2 <- NA*numeric(length(y))
    gamma_1[na_ind] <- sapply(py[na_ind], function(x){gammainc(-log(x), 1 - shape)[1]})
    gamma_2[na_ind] <- sapply(p0[na_ind], function(x){gammainc(-2*log(x), 1 - shape)[1]})

    crps1 <- (location - y)*(1 - 2*py)
    crps2 <- location*(p0^2)
    crps3 <- 2*(scale/shape)*(1 - py - gamma_1)
    crps4 <- (scale/shape)*(1 - (p0^2) - (2^shape)*gamma_2)

    crps_pos <- crps1 + crps2 - crps3 + crps4


    crps <- ((eps - shape_0)*crps_neg + (eps + shape_0)*crps_pos)/(2*eps)
  }
  return(crps)
}
crps_gev_obj <- function(par, y, ens.mu, ens.sd) {
  mu <- par[1] + par[2]*ens.mu
  sig <- (par[3]^2) + (par[4]^2)*ens.sd
  shape <- par[5]^2

  crps_cgev(y, location = mu, scale = sig, shape = shape) |> mean()
}
fit_emos_cgev <- function(trai, test) {

  outpar <- optim(c(1, 1, 1, 1, 0.5), crps_gev_obj, y=trai$obs, ens.mu=trai$ens.mu, ens.sd=trai$ens.sd)$par
  mu <- outpar[1] + outpar[2]*test$ens.mu
  sig <- (outpar[3]^2) + (outpar[4]^2)*test$ens.sd
  shape <- outpar[5]^2

  pred <- cbind(as.numeric(mu), as.numeric(sig), as.numeric(shape))
  colnames(pred) <- c("location", "scale", "shape")
  pred
}


# estimate censored logistic emos parameters at each station
emos_preds_cl <- lapply(seq_along(stat_ids), function(j) {
  st <- stat_ids[j]
  print(paste0('Forecast at Station: ', st, ' (', j, ' from ', length(stat_ids), ')'))
  trai <- get_train_data(train_obs, train_fc, train_times, train_years, stat_ids, lead_times[lt], st)
  test <- get_test_data(test_obs, test_fc, test_times, stat_ids, lead_times[lt], st)
  fit_emos_crch(trai, test)
})
emos_preds_cl <- simplify2array(emos_preds_cl)

fc_emos_cl <- list(F_x = F_x, location = t(emos_preds_cl[, 1, ]), scale = t(emos_preds_cl[, 2, ]))


# estimate censored GEV emos parameters at each station
emos_preds_cgev <- lapply(seq_along(stat_ids), function(j) {
  st <- stat_ids[j]
  print(paste0('Forecast at Station: ', st, ' (', j, ' from ', length(stat_ids), ')'))
  trai <- get_train_data(train_obs, train_fc, train_times, train_years, stat_ids, lead_times[lt], st)
  test <- get_test_data(test_obs, test_fc, test_times, stat_ids, lead_times[lt], st)
  fit_emos_cgev(trai, test)
})
emos_preds_cgev <- simplify2array(emos_preds_cgev)

F_x <- function(q, location, scale, shape) {
  if (length(q) == 1) q <- rep(q, length(location))
  out <- sapply(seq_along(location), function(i) pgev(q[i], location[i], scale[i], shape[i]))
  out[q < 0] <- 0
  return(out)
}

fc_emos_cgev <- list(F_x = F_x,
                     location = t(emos_preds_cgev[, 1, ]),
                     scale = t(emos_preds_cgev[, 2, ]),
                     shape = t(emos_preds_cgev[, 3, ]))


################################################################################
## save forecasts

# save auxiliary data
aux_data <- list(obs = test_obs[, lt, ],
                 stat_id = stat_ids,
                 stat_info = lonlatalt,
                 time = test_times,
                 lead = lead_times[lt])

fc_dat <- list(ifs = fc_ifs, smooth = fc_smooth, emos_cl = fc_emos_cl, emos_cgev = fc_emos_cgev, aux_data = aux_data)

saveRDS(fc_dat, file = "scripts/cs_fc_data.RDS")

