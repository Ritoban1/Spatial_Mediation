summarize_parameters <- function(mle, var_par, ok, truth) {
  do.call(rbind, lapply(seq_along(truth), function(j) {
    valid <- ok & is.finite(mle[, j])
    x <- mle[valid, j]
    t <- truth[j]
    v <- var_par[, j]
    good_se <- valid & is.finite(v) & v > 0
    se <- sqrt(v[good_se])
    data.frame(
      parameter = colnames(mle)[j], truth = t,
      mean_est = if (length(x)) mean(x) else NA_real_,
      bias = if (length(x)) mean(x - t) else NA_real_,
      RMSE = if (length(x)) sqrt(mean((x - t)^2)) else NA_real_,
      ESE = if (length(x) > 1L) sd(x) else NA_real_,
      ASE = if (length(se)) mean(se) else NA_real_,
      coverage = if (length(se)) {
        mean(abs(mle[good_se, j] - t) <= qnorm(0.975) * se)
      } else NA_real_,
      n_fit = length(x), n_SE = length(se)
    )
  }))
}

summarize_effects <- function(effects, ok) {
  truth <- effects$truth
  est <- effects$estimates
  se <- effects$SE
  do.call(rbind, lapply(seq_along(truth), function(j) {
    has_fit <- ok & is.finite(est[, j])
    good_se <- has_fit & is.finite(se[, j]) & se[, j] > 0
    e <- est[has_fit, j]
    tr <- unname(truth[j])
    ss <- se[good_se, j]
    if (!length(e)) {
      return(data.frame(effect = names(truth)[j], truth = tr,
                        mean_est = NA_real_, bias = NA_real_,
                        RMSE = NA_real_, RRMSE = NA_real_,
                        ESE = NA_real_, ASE = NA_real_, coverage = NA_real_,
                        n_fit = 0L, n_SE = 0L))
    }
    RMSE <- sqrt(mean((e - tr)^2))
    data.frame(
      effect = names(truth)[j], truth = tr, mean_est = mean(e),
      bias = mean(e - tr), RMSE = RMSE,
      RRMSE = if (abs(tr) > 1e-12) RMSE / abs(tr) else NA_real_,
      ESE = if (length(e) > 1L) sd(e) else NA_real_,
      ASE = if (length(ss)) mean(ss) else NA_real_,
      coverage = if (length(ss)) {
        mean(abs(est[good_se, j] - tr) <= qnorm(0.975) * ss)
      } else NA_real_,
      n_fit = length(e), n_SE = length(ss)
    )
  }))
}
