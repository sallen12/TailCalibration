# from https://github.com/resinj/replication_GR21

# reliability diagram for mean and quantile calibration
reldiag = function(x,
                   y,
                   type = "mean",      # type = list("quantile", alpha = 0.1)
                   resampling = TRUE, n_resamples = 1000, replace = TRUE, region_level = 0.9,
                   digits = 3, inset_hist = TRUE, hist.breaks = 8, scatter_plot = TRUE, ext_decomp = TRUE,
                   lim = NULL, main = "", xlab = NULL, ylab = NULL, adj_xlab = NA) {

  if (type[[1]] == "mean") {
    pava = function(x, y) isoreg(x, y)$yf
    score = function(x, y) mean((x - y)^2) # canonical score
    marg = base::mean
    identif = function(x,y) x - y
    if(is.null(xlab)) xlab = expression(x == m[1](F))
    score_label = "MSE "
  } else {
    if (type[[1]] == "quantile") {
      require(isotone)
      alpha = type[[2]]
      pava = function(x, y){
        ranking = match(1:length(x), order(x, y, decreasing = c(FALSE, TRUE)))
        return(gpava(ranking, y, solver = weighted.fractile, p = alpha, ties = "secondary")$x)
      }
      score = function(x, y) mean(2*(as.numeric(x >= y) - alpha)*(x - y))
      marg = function(x) quantile(x, alpha, type = 1)
      identif = function(x, y) as.numeric(x > y) - alpha
      if(is.null(xlab)) xlab = bquote(x == q[.(alpha)](F))
      score_label = "QS "
    }
    else stop("type must be \"mean\" or list(\"quantile\",level)")
  }

  if (is.null(lim)) lim = range(x) + c(-1, 1)*max(abs(range(x)))*0.2
  if (is.null(ylab)) ylab = expression(hat(x)[rc])

  plot(NULL, xlim = lim, ylim = lim, main = main, xlab = "", ylab = ylab)
  mtext(xlab, side = 1, line = par()$mgp[1], cex = par()$cex, adj = adj_xlab)

  ord_x = order(x)
  x = x[ord_x]
  y = y[ord_x]

  x_rc = pava(x,y)

  res = y - x

  # score decomposition
  s = score(x, y)
  c_rc_ucond = optim(par = 0, fn = function(c) score(x + c, y),
                     method = "Brent", lower = min(res), upper = max(res))$par
  s_rc_ucond = score(x + c_rc_ucond, y)
  s_rc = score(x_rc, y)
  s_mg = score(marg(y), y)

  mcb = s - s_rc
  umcb = s - s_rc_ucond
  cmcb = s_rc_ucond - s_rc
  dsc = s_mg - s_rc
  unc = s_mg

  # compute p-value for hypothesis of unconditional calibration (t-test)
  v = identif(x, y)
  t = sqrt(length(v)) * mean(v)/sd(v)
  pval_ucond = 1 - abs(pt(t, length(v) - 1) - 0.5)*2

  if (resampling) {
    low = floor(n_resamples * (1 - region_level)/2)
    up = n_resamples - low

    resamples = sapply(1:n_resamples, function(i) x + sample(res, length(y), replace = replace))
    x_rc_resamples = apply(resamples, 2, function(y) pava(x, y))
    x_rc_resamples_sorted = apply(x_rc_resamples, 1, sort) - marg(res)

    ran_x = range(x)
    polygon(c(ran_x[1], x,ran_x[2], rev(x), ran_x[1]),
            c(ran_x[1], x_rc_resamples_sorted[up, ], ran_x[2], rev(x_rc_resamples_sorted[low, ]), ran_x[1]),
            border = NA, col = "lightblue1")
    points(x, x_rc_resamples_sorted[low, ], type = "l", lty = 1, col = "lightblue2")
    points(x, x_rc_resamples_sorted[up, ], type = "l", lty = 1, col = "lightblue2")
    box()

    # compute Monte-Carlo p-value (for hypothesis of conditional calibration)
    mcb_resamples = sapply(1:n_resamples,
                           function(i) score(x, resamples[ , i]) - score(x_rc_resamples[ , i], resamples[ , i]))
    mcb_bounds = sort(c(mcb, mcb_resamples))[c(low, up)]
    rank_obs = tail(rank(c(mcb_resamples, mcb)),1)
    pval = 1 - (rank_obs - 1)/(n_resamples + 1)
  }

  if (ext_decomp) {
    mcb_label = c(expression(MCB[u]),expression(MCB[c]))
    mcb_comp = c(umcb,cmcb)
  } else {
    mcb_label = "MCB"
    mcb_comp = mcb

  }
  s = sum(round(c(mcb_comp, -dsc, unc), digits))

  text_pos = legend("topleft", legend = parse(text = c(score_label, mcb_label, "DSC", "UNC")), plot = FALSE)
  text(x = lim[1] ,y = text_pos$text$y, parse(text = c(score_label, mcb_label, "DSC", "UNC")), adj = c(0,0.5))
  text(x = 0.8*lim[1] + 0.2*lim[2], y = text_pos$text$y,
       bquote(.(format(round(c(s, mcb_comp, dsc,unc), digits = digits)), nsmall = digits)),
       adj = c(0, 0.5))

  abline(a = 0, b = 1, col = "grey", lty = 2)
  points(x, x_rc, type = "l")

  if (scatter_plot) {
    points(x, y, pch = 20, col = adjustcolor("black", alpha = 0.25), cex = 0.5)
  }

  if (inset_hist) {
    a = par("usr")
    a = c(grconvertX(a[1:2], "user", "ndc"), grconvertY(a[3:4], "user", "ndc"))
    par.old = par(fig = c(0.3*a[1] + 0.7*a[2], 0.05*a[1] + 0.95*a[2],0.9*a[3] + 0.1*a[4], 0.65*a[3] + 0.35*a[4]),
                  pty = "m", mar = c(1, 0, 0, 0), mgp = c(1, 0.4, 0), tcl = -0.25, new = TRUE)
    plot(hist(x, breaks = hist.breaks, main = "", yaxt = "n", xlab = "", ylab = "", xlim = lim), add = TRUE)
    par(par.old)
  }

  return(list(x = x, y = y, res = res, x_rc = x_rc,
              decomp = c(umcb, cmcb, mcb, dsc, unc),
              pval_u = pval_ucond, pval_c = if (resampling) pval else NA))
}


# reliability diagram for threshold calibration
threshreldiag = function(fcast, y, t, resampling = TRUE, n_resamples = 1000, region_level = 0.9,
                         digits = 3, inset_hist = TRUE, main = "", xlab = NULL, ylab = NULL, ...){
  if (is.null(xlab)) xlab = bquote(x == F(.(t)))
  if (is.null(ylab)) ylab = expression(x[rc])

  # binarize problem
  x = fcast$F(t)
  y = as.numeric(y <= t)
  sample = function() as.numeric(runif(y) <= x)
  score = function(x, y) mean((x - y)^2)

  plot(NULL, xlim = c(0, 1), ylim = c(0, 1), main = main, xlab = xlab, ylab = ylab)

  x_rc = isoreg(x, y)$yf # values correspond to ORDERED forecast values!

  # compute score decomposition
  s = score(x, y)
  s_rc = score(x_rc, y[order(x)])
  s_mg = score(mean(y), y)
  mcb = s - s_rc
  dsc = s_mg - s_rc
  unc = s_mg

  if (resampling) {
    low = floor(n_resamples * (1 - region_level)/2)
    up = n_resamples - low
    pval_digits = ceiling(log(n_resamples, 10))

    resamples = sapply(1:n_resamples, function(i) sample())
    x_rc_resamples = apply(resamples, 2, function(y) isoreg(x, y)$yf)
    x_rc_resamples_sorted = apply(x_rc_resamples, 1, sort)

    ran_x = range(x)
    polygon(c(ran_x[1], sort(x), ran_x[2], rev(sort(x)), ran_x[1]),
            c(ran_x[1], x_rc_resamples_sorted[up, ], ran_x[2], rev(x_rc_resamples_sorted[low, ]), ran_x[1]),
            border = NA, col = "lightblue1")
    points(sort(x), x_rc_resamples_sorted[low, ], type = "l", lty = 1, col = "lightblue2")
    points(sort(x), x_rc_resamples_sorted[up, ], type = "l", lty = 1, col = "lightblue2")
    box()

    # compute Monte-Carlo p-value
    mcb_resamples = sapply(1:n_resamples, function(i) score(x, resamples[ , i]) -
                             score(x_rc_resamples[ , i], resamples[order(x), i]))
    rank_obs = tail(rank(c(mcb_resamples, mcb)), 1)
    pval = 1 - (rank_obs - 1)/(n_resamples + 1)
  }

  s = sum(round(c(mcb, -dsc, unc), digits))

  text(x = 0, y = 1, paste0(c("BS ", "MCB ", "DSC ", "UNC "), collapse = "\n"), adj = c(0, 1))
  text(x = 0.2, y = 1, paste0(bquote(.(format(round(c(s, mcb, dsc, unc), digits = digits)),
                                       nsmall = digits)), collapse = "\n"), adj = c(0, 1))

  abline(a = 0, b = 1, col = "grey", lty = 2)
  points(sort(x), x_rc, type = "l")

  if (inset_hist) {
    a = par("usr")
    a = c(grconvertX(a[1:2], "user", "ndc"), grconvertY(a[3:4], "user", "ndc"))
    par.old = par(fig = c(0.3*a[1] + 0.7*a[2], 0.05*a[1] + 0.95*a[2], 0.9*a[3] + 0.1*a[4], 0.65*a[3] + 0.35*a[4]),
                  pty = "m", mar = c(1, 0, 0, 0), mgp = c(1, 0.4, 0), tcl = -0.25, new = TRUE)
    plot(hist(x, main = "", yaxt = "n", xlab = "", ylab = ""), add = TRUE)
    par(par.old)
  }

  return(list(decomp = c(mcb, dsc, unc), pval = if (resampling) pval else NULL), x = data.frame(x = sort(x), x_rc = x_rc))
}
