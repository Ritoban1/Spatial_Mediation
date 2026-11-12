# ================================================================
# JRSSC spatial mediation simulation: Gaussian conditioning PRIMARY
# ================================================================
# Replaces the original data-generation/MLE/variance script.
# Requires only base R. If an old graph is unavailable or reuse is
# disabled, it also needs the original dist.rds distance matrix.
#
# Changes:
#   * Preserve the original ER + distance-weighted E and Gaussian
#     mediator/outcome equations and regression coefficients.
#   * Replace random, across-node standardization of A with fixed
#     population-scale centering/scaling. This preserves roughly the
#     original A scale while retaining independent Gaussian errors
#     conditional on C and the fixed graph.
#   * For each Gaussian block, profile out coefficients and the
#     residual variance, optimizing only ONE covariance ratio.
#   * Diagonalize K = E %*% t(E) only ONCE for each graph.
#   * Use coherent covariance parameters: sigma^2 I + sigma_b^2 K.
#   * Use conditional expected Gaussian information for mean
#     coefficients (no custom 2n-by-2n Hessian or pseudoinverse).
#   * Correct weighted tau3 and weighted Gaussian-conditioned tau6.
#     Under conditional independent equal-variance Gaussian exposures,
#     the conditional slope of S2_i on S1_i is calculated exactly.
#     A uniform-local shift is retained only as an optional comparison.
#   * Save MLEs, variances, all six Gaussian-primary effects, SEs,
#     and summary tables after each sample size.
#
# Interpretation: six FOCAL-LOCAL Gaussian-conditioned pathway
# contrasts, NOT an exhaustive effect of a global county intervention.
# Political propensity is bounded in reality; the unbounded Gaussian
# exposure allocation here is a mathematical working intervention.
# ================================================================

# From the repository root: Rscript Simulations/Spatial_Mediation_MLE_Gaussian_Primary.R
# Optional arguments: DIST_FILE OUTPUT_DIR
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) {
  stop("Usage: Rscript Simulations/Spatial_Mediation_MLE_Gaussian_Primary.R [DIST_FILE [OUTPUT_DIR]]")
}
DIST_FILE <- if (length(args) >= 1L) args[[1L]] else "Data_Analysis/Data/dist.rds"
OUT_DIR <- if (length(args) >= 2L) args[[2L]] else "Simulations/Results"
source("Simulations/Simulation_Summaries.R")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------- User settings ----------------
SAMPLE_SIZES <- c(200L, 500L)     # For one run, e.g. 500L
N_REPS <- 500L
SEED <- 100000000L
MEAN_DEGREE <- 10
REUSE_SAVED_GRAPH <- FALSE      # Regenerate the published graphs from the fixed seeds
USE_FIXED_A_STANDARDIZATION <- TRUE  # Keep original approximate scale, iid Gaussian errors
SAVE_UNIFORM_LOCAL_COMPARISON <- TRUE   # Optional second set of metrics

# Original coefficients and error SDs
beta_true  <- c(-2, 1.5, 1, 1.2, 0.4, 0.2, 0.3)
gamma_true <- c(-1, 2, 1, 0.5, 0.2)
delta_true <- c(-0.5, 0.5, 0)
sd_a <- 0.25
sd_m <- 0.5
sd_y <- 0.5
sd_bm <- 0.25
sd_by <- 0.25

# Same numerical local contrasts as the original analysis
reference <- c(own = 0.1, neighbor = 0.2)
target    <- c(own = 0.8, neighbor = 0.9)
d_own <- unname(target["own"] - reference["own"])
d_nbr <- unname(target["neighbor"] - reference["neighbor"])

# ---------------- County distances (loaded only if generating graphs) ----------------
load_county_distances <- function(path) {
  if (!file.exists(path)) stop("Missing distance file: ", path)
  distances <- as.matrix(readRDS(path))
  stopifnot(nrow(distances) == ncol(distances))
  bad <- which(!is.finite(distances[1L, ]))
  if (length(bad)) distances <- distances[-bad, -bad, drop = FALSE]
  if (any(!is.finite(distances))) {
    stop("dist.rds contains additional non-finite distances; inspect its entries.")
  }
  if (!isTRUE(all.equal(distances, t(distances), tolerance = 1e-8))) {
    stop("Distance matrix is not symmetric.")
  }
  if (max(SAMPLE_SIZES) > nrow(distances)) {
    stop("Not enough counties for the requested sample size.")
  }
  distances
}
dist_all <- NULL

# ---------------- Graph and effect coefficients ----------------
prepare_network <- function(E) {
  E <- as.matrix(E)
  stopifnot(nrow(E) == ncol(E), all(is.finite(E)))
  if (!isTRUE(all.equal(E, t(E), tolerance = 1e-8))) {
    stop("Spatial matrix E must be symmetric.")
  }
  W <- E
  diag(W) <- 0
  strength <- rowSums(W)
  W <- sweep(W, 1L, ifelse(strength > 0, strength, 1), "/")
  list(E = E, W = W, strength = strength)
}

make_network <- function(n, distances) {
  ids <- sample.int(nrow(distances), n, replace = FALSE)
  dsub <- distances[ids, ids, drop = FALSE]

  # Symmetric ER connectivity (original p = 10/n).
  B <- matrix(0, n, n)
  B[upper.tri(B)] <- rbinom(n * (n - 1L) / 2L, 1L, MEAN_DEGREE / n)
  B <- B + t(B)

  # Distance-weighted E; self-loading is 1, as in original code.
  E <- B * exp(-dsub / 1000)
  diag(E) <- 1

  prepare_network(E)
}

network_coefficients <- function(E, W, strength) {
  n <- nrow(W)
  W2 <- W %*% W
  adj <- E > 0
  diag(adj) <- FALSE

  # For the local own-A -> neighboring-M -> own-Y route.
  rho <- mean(diag(W2))

  # For neighbor-A -> neighboring-M -> own-Y, count only
  # first-degree-to-first-degree two-step paths (triangles).
  kappa <- sum(W2[adj]) / n

  # PRIMARY tau6: Gaussian conditioning. Conditional on the FULL
  # fixed covariate array and graph, A_k have independent Gaussian
  # errors of the SAME variance. For S1_i = sum_k W_ik A_k and
  # S2_i = sum_{k in N_i} (W^2)_ik A_k, the conditional slope
  # d E(S2_i | S1_i = s, H) / ds equals
  #      sum_k (W^2)_ik W_ik / sum_k W_ik^2.
  # Any heterogeneity in the conditional means cancels between
  # the two values of s; equal conditional variances are essential.
  den_gauss <- rowSums(W^2)
  num_gauss <- rowSums(W2 * W)
  ratio_gauss <- numeric(n)
  active_rows <- den_gauss > 0
  ratio_gauss[active_rows] <- num_gauss[active_rows] / den_gauss[active_rows]
  kappa_gaussian <- mean(ratio_gauss)

  # An isolated county has no neighboring summary or shifted
  # neighbor effect: all slope contributions for it are zero.
  active <- mean(strength > 0)

  list(rho = rho, kappa = kappa,
       kappa_gaussian = kappa_gaussian,
       active = active, isolates = sum(strength == 0L),
       edges = sum(adj) / 2)
}

# ---------------- Gaussian spatial profiled MLE ----------------
# Model: y | X ~ N(X b, v I + w K), K = E %*% t(E),
# with v = residual VARIANCE, w = random-effect VARIANCE.
# Parametrize lambda = w/v. Given lambda, estimate b by GLS and
# v = weighted residual SS / n, leaving a one-dimensional search.
# eig is the single eigendecomposition of K for the fixed graph.
fit_spatial_gaussian <- function(y, X, eig) {
  n <- length(y)
  X <- as.matrix(X)
  U <- eig$vectors
  kval <- pmax(eig$values, 0)
  Xstar <- crossprod(U, X)
  ystar <- drop(crossprod(U, y))

  gls_at_lambda <- function(lambda) {
    den <- 1 + lambda * kval
    sw <- 1 / sqrt(den)
    Xw <- Xstar * sw  # R recycles sw down each column
    yw <- ystar * sw
    XtX <- crossprod(Xw)
    Uchol <- tryCatch(chol(XtX), error = function(e) NULL)
    if (is.null(Uchol)) return(list(objective = Inf))

    rhs <- crossprod(Xw, yw)
    coef <- drop(backsolve(Uchol, forwardsolve(t(Uchol), rhs)))
    residual <- drop(ystar - Xstar %*% coef)
    rss <- sum(residual^2 / den)
    if (!is.finite(rss) || rss <= 0) return(list(objective = Inf))
    v <- rss / n
    # n/2 * (1 + log(2*pi)) is constant over lambda.
    objective <- 0.5 * (sum(log(den)) + n * log(v))
    list(objective = objective, coef = coef, v = v,
         lambda = lambda, XtX = XtX, Uchol = Uchol)
  }

  # Short coarse scan, then local refinement. Include lambda=0 as
  # a legitimate finite-sample boundary solution.
  tgrid <- seq(log(1e-7), log(1e4), length.out = 23L)
  grid_values <- vapply(tgrid, function(t) gls_at_lambda(exp(t))$objective,
                        numeric(1))
  if (!any(is.finite(grid_values))) stop("Singular mean design in spatial GLS.")
  best_index <- which.min(grid_values)
  left <- tgrid[max(1L, best_index - 1L)]
  right <- tgrid[min(length(tgrid), best_index + 1L)]
  local_opt <- optimize(function(t) gls_at_lambda(exp(t))$objective,
                        interval = c(left, right), tol = 1e-6)
  candidates <- c(0, exp(tgrid[best_index]), exp(local_opt$minimum))
  objs <- vapply(candidates, function(lambda) gls_at_lambda(lambda)$objective,
                 numeric(1))
  best <- gls_at_lambda(candidates[which.min(objs)])
  if (!is.finite(best$objective)) stop("Spatial GLS failed at the optimum.")

  v <- best$v
  w <- best$lambda * v
  sigma <- sqrt(v)
  sigma_b <- sqrt(w)

  # Conditional expected information for regression coefficients:
  # X^T Sigma^{-1} X = (1/v) Xstar^T diag(1/(1+lambda*k)) Xstar.
  vcov_coef <- v * chol2inv(best$Uchol)

  # Optional expected information for log v, log w, subsequently
  # converted to the variances of sigma=sqrt(v), sigma_b=sqrt(w).
  # At w=0 this variance-component inference is nonregular; do
  # not fabricate a standard error for sigma_b on the boundary.
  var_sigma <- c(NA_real_, NA_real_)
  if (w > 0) {
    den <- 1 + best$lambda * kval
    Ueta <- cbind(1 / den, best$lambda * kval / den)
    info_eta <- 0.5 * crossprod(Ueta)
    inv_eta <- tryCatch(chol2inv(chol(info_eta)), error = function(e) NULL)
    if (!is.null(inv_eta)) {
      var_sigma <- c((sigma / 2)^2 * inv_eta[1, 1],
                     (sigma_b / 2)^2 * inv_eta[2, 2])
    }
  }

  list(coef = best$coef, sigma = sigma, sigma_b = sigma_b,
       vcov_coef = vcov_coef, var_sigma = var_sigma,
       lambda = best$lambda, boundary = (best$lambda == 0))
}

# ---------------- Six Gaussian-conditioned effects + uncertainty ----------------
# mle ordering, same as original file:
#   1:7  beta0,...,beta6    8:9   SD(y), SD(b_y)
#  10:14 gamma0,...,gamma4  15:16 SD(m), SD(b_m)
# var_par, as in original file:
#   1:16 matching diagonal variances; column 17 = Cov(gamma1,gamma2)
compute_effects <- function(mle, var_par, net,
                            mode = c("gaussian_conditional", "uniform_local")) {
  mode <- match.arg(mode)
  b <- mle[, 1:7, drop = FALSE]
  g <- mle[, 10:14, drop = FALSE]
  V <- var_par
  h <- net$active
  r <- net$rho
  k <- if (mode == "uniform_local") net$kappa else net$kappa_gaussian

  est <- cbind(
    tau1 = d_own * b[, 2],
    tau2 = d_own * b[, 4] * g[, 2],
    tau3 = d_own * r * b[, 5] * g[, 3],
    tau4 = d_nbr * h * b[, 3],
    tau5 = d_nbr * h * b[, 4] * g[, 3],
    tau6 = d_nbr * b[, 5] * (h * g[, 2] + k * g[, 3])
  )

  # Delta-method effect variances. In particular, the Gaussian
  # policy coefficient k is a KNOWN graph summary (equal-variance
  # assumption), not an exposure nuisance-parameter estimate.
  # Outcome and mediator score
  # blocks are orthogonal in expected information under the
  # correctly specified conditional Gaussian structural model.
  vv <- cbind(
    tau1 = d_own^2 * V[, 2],
    tau2 = d_own^2 * (g[, 2]^2 * V[, 4] + b[, 4]^2 * V[, 11]),
    tau3 = (d_own * r)^2 * (g[, 3]^2 * V[, 5] + b[, 5]^2 * V[, 12]),
    tau4 = (d_nbr * h)^2 * V[, 3],
    tau5 = (d_nbr * h)^2 * (g[, 3]^2 * V[, 4] + b[, 4]^2 * V[, 12]),
    tau6 = d_nbr^2 * (
      (h * g[, 2] + k * g[, 3])^2 * V[, 5] +
      b[, 5]^2 * (h^2 * V[, 11] + k^2 * V[, 12] +
                     2 * h * k * V[, 17])
    )
  )
  se <- sqrt(pmax(vv, 0))
  se[!is.finite(vv) | vv <= 0] <- NA_real_

  truth <- c(
    tau1 = d_own * beta_true[2],
    tau2 = d_own * beta_true[4] * gamma_true[2],
    tau3 = d_own * r * beta_true[5] * gamma_true[3],
    tau4 = d_nbr * h * beta_true[3],
    tau5 = d_nbr * h * beta_true[4] * gamma_true[3],
    tau6 = d_nbr * beta_true[5] * (h * gamma_true[2] + k * gamma_true[3])
  )
  list(truth = truth, estimates = est, variances = vv,
       SE = se, definition = mode)
}

# ---------------- Run both sample sizes ----------------
for (n in SAMPLE_SIZES) {
  set.seed(SEED + n)
  cat("\n===== Spatial simulation: N =", n, ", repetitions =", N_REPS, "=====\n")
  original_rds <- file.path(OUT_DIR, paste0("New_simu_", n, ".rds"))
  graph_source <- "new ER distance-weighted graph"
  if (REUSE_SAVED_GRAPH && file.exists(original_rds)) {
    saved <- readRDS(original_rds)
    if (length(saved) < 1L || nrow(saved[[1]]) != n) {
      stop("Saved network dimension mismatch in ", original_rds)
    }
    graph <- prepare_network(saved[[1]])
    rm(saved)
    graph_source <- paste0("reused E from ", original_rds)
  } else {
    if (is.null(dist_all)) dist_all <- load_county_distances(DIST_FILE)
    graph <- make_network(n, dist_all)
  }
  E <- graph$E
  W <- graph$W
  net <- network_coefficients(E, W, graph$strength)
  if (!is.finite(net$kappa_gaussian) || !is.finite(net$rho)) {
    stop("Nonfinite Gaussian-conditioning or weighted tau3 graph coefficient.")
  }
  cat(sprintf("Network: rho = %.6f, kappa_G = %.6f, kappa_uniform = %.6f, isolates = %d\n",
              net$rho, net$kappa_gaussian, net$kappa, net$isolates))
  K <- tcrossprod(E)
  # Gaussian-conditioning requires exposure assignments with equal
  # conditional variances given the full C array and fixed graph.
  # Use a fixed (non-exposure-realization-dependent) A scale.
  # delta3=0 here, giving sqrt(delta2^2 + sd_a^2). The expression
  # below is a convenient approximate total scale if delta3 != 0;
  # the conditional error variance remains equal across counties.
  A_fixed_scale <- sqrt(delta_true[2]^2 + sd_a^2 +
                        delta_true[3]^2 * mean(rowSums(W^2)))
  eig <- eigen(K, symmetric = TRUE)   # ONCE per N, not inside optim

  mle <- matrix(NA_real_, N_REPS, 16L)
  var_par <- matrix(NA_real_, N_REPS, 17L)
  ok <- rep(FALSE, N_REPS)
  boundary_y <- boundary_m <- rep(NA, N_REPS)
  errors <- rep("", N_REPS)
  colnames(mle) <- c(paste0("beta", 0:6), "sd_y", "sd_by",
                     paste0("gamma", 0:4), "sd_m", "sd_bm")
  colnames(var_par) <- c(paste0("Var_", colnames(mle)),
                         "Cov_gamma1_gamma2")

  for (rep_id in seq_len(N_REPS)) {
    # Reproduce original covariate standardization. Conditional
    # likelihood treats this full covariate array as given.
    C <- as.numeric(scale(rnorm(n)))
    WC <- drop(W %*% C)

    # IMPORTANT: do NOT standardize A using its observed mean/SD.
    # Fixed affine centering/scaling preserves conditional Gaussian
    # independence and keeps the original approximate A mean 0, SD 1.
    A_raw <- delta_true[1] + delta_true[2] * C +
             delta_true[3] * WC + rnorm(n, sd = sd_a)
    A <- if (USE_FIXED_A_STANDARDIZATION) {
      (A_raw - delta_true[1]) / A_fixed_scale
    } else A_raw
    WA <- drop(W %*% A)

    M <- gamma_true[1] + gamma_true[2] * A + gamma_true[3] * WA +
         gamma_true[4] * C + gamma_true[5] * WC +
         drop(E %*% rnorm(n, sd = sd_bm)) + rnorm(n, sd = sd_m)
    WM <- drop(W %*% M)

    Y <- beta_true[1] + beta_true[2] * A + beta_true[3] * WA +
         beta_true[4] * M + beta_true[5] * WM +
         beta_true[6] * C + beta_true[7] * WC +
         drop(E %*% rnorm(n, sd = sd_by)) + rnorm(n, sd = sd_y)

    X_m <- cbind(1, A, WA, C, WC)
    X_y <- cbind(1, A, WA, M, WM, C, WC)

    attempt <- tryCatch({
      fit_m <- fit_spatial_gaussian(M, X_m, eig)
      fit_y <- fit_spatial_gaussian(Y, X_y, eig)
      list(m = fit_m, y = fit_y)
    }, error = function(e) e)

    if (inherits(attempt, "error")) {
      errors[rep_id] <- conditionMessage(attempt)
      next
    }

    fy <- attempt$y
    fm <- attempt$m
    mle[rep_id, ] <- c(fy$coef, fy$sigma, fy$sigma_b,
                       fm$coef, fm$sigma, fm$sigma_b)
    var_par[rep_id, ] <- c(diag(fy$vcov_coef), fy$var_sigma,
                           diag(fm$vcov_coef), fm$var_sigma,
                           fm$vcov_coef[2, 3])
    boundary_y[rep_id] <- fy$boundary
    boundary_m[rep_id] <- fm$boundary
    ok[rep_id] <- TRUE

    if (rep_id %% 50L == 0L) {
      cat("N =", n, ":", rep_id, "replications attempted;",
          sum(ok), "successful fits\n")
    }
  }

  # The default publication table uses Gaussian conditioning.
  effects <- compute_effects(mle, var_par, net, mode = "gaussian_conditional")
  metrics <- summarize_effects(effects, ok)
  effects_uniform_local <- if (SAVE_UNIFORM_LOCAL_COMPARISON) {
    compute_effects(mle, var_par, net, mode = "uniform_local")
  } else NULL
  metrics_uniform_local <- if (SAVE_UNIFORM_LOCAL_COMPARISON) {
    summarize_effects(effects_uniform_local, ok)
  } else NULL
  parameter_metrics <- summarize_parameters(
    mle, var_par, ok, c(beta_true, sd_y, sd_by, gamma_true, sd_m, sd_bm)
  )
  diagnostics <- data.frame(
    n = n, n_reps = N_REPS, n_fit = sum(ok), n_failed = sum(!ok),
    graph_source = graph_source,
    edges = net$edges, isolates = net$isolates,
    rho_tau3 = net$rho, kappa_tau6_local = net$kappa,
    kappa_tau6_gaussian = net$kappa_gaussian,
    frac_nonisolates = net$active,
    outcome_boundary_fits = sum(boundary_y %in% TRUE),
    mediator_boundary_fits = sum(boundary_m %in% TRUE)
  )

  # Backward compatibility: result[[1]], [[2]], [[3]] have the
  # same roles and parameter/variance column order as previously.
  result <- list(
    E = E, mle = mle, var_par = var_par,
    metrics = metrics, effects = effects,
    metrics_uniform_local = metrics_uniform_local,
    effects_uniform_local = effects_uniform_local,
    parameter_metrics = parameter_metrics,
    diagnostics = diagnostics,
    success = ok, errors = errors,
    config = list(n = n, n_reps = N_REPS, seed = SEED + n,
                  graph_source = graph_source,
                  beta = beta_true, gamma = gamma_true,
                  delta_raw = delta_true, sd_a_raw = sd_a,
                  fixed_A_standardization = USE_FIXED_A_STANDARDIZATION,
                  A_fixed_scale = A_fixed_scale,
                  sigma_m = sd_m, sigma_y = sd_y,
                  sigma_bm = sd_bm, sigma_by = sd_by,
                  reference = reference, target = target,
                  effect_interpretation = "Gaussian-conditioned focal-local contrast; second-degree sources marginalized under independent exposure law",
                  primary_effect = "gaussian_conditional",
                  uniform_comparison_saved = SAVE_UNIFORM_LOCAL_COMPARISON)
  )

  base <- file.path(OUT_DIR, paste0("New_simu_", n, "_gaussian_primary"))
  saveRDS(result, paste0(base, ".rds"))
  write.csv(metrics, paste0(base, "_metrics.csv"), row.names = FALSE)
  if (SAVE_UNIFORM_LOCAL_COMPARISON) {
    write.csv(metrics_uniform_local, paste0(base, "_metrics_uniform_local.csv"),
              row.names = FALSE)
  }
  write.csv(parameter_metrics, paste0(base, "_parameters.csv"),
            row.names = FALSE)
  write.csv(diagnostics, paste0(base, "_diagnostics.csv"), row.names = FALSE)
  cat("\nNetwork and fit diagnostics:\n")
  print(diagnostics, row.names = FALSE, digits = 5)
  cat("\nGAUSSIAN-CONDITIONING PRIMARY: all six effects (truth, bias, RMSE, RRMSE, ESE, ASE, CP):\n")
  print(metrics, row.names = FALSE, digits = 5)
  if (SAVE_UNIFORM_LOCAL_COMPARISON) {
    cat("\nSecondary comparison: uniform-local tau6 (not primary):\n")
    print(metrics_uniform_local[6L, , drop = FALSE], row.names = FALSE, digits = 5)
  }
  cat("\nSaved:", paste0(base, ".rds"), "and summary CSVs\n")
}
