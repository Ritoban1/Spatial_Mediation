# Uniform spatial mediation application with and without mobility adjustment.
# From the repository root: Rscript Data_Analysis/Main.R
# Optional arguments: DATA_DIR DIST_FILE OUTPUT_DIR [REGION]
# REGION: all (default), USA, Midwest, Northeast, Southeast, or West.
# One change at a time, relative to the main analysis:
# Main                 : state weighted, with mobility, original contrast.
# Alternative contrast : state weighted, with mobility, alternative contrast.
# Binary B             : binary B, with mobility, original contrast.
# Without mobility     : state weighted, without mobility, original contrast.
# Binary B keeps distance decay. Without mobility removes C2 and WC2 from
# both equations. All four analyses use the main complete-case county sample.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) {
  args <- c("Data_Analysis/Data", "Data_Analysis/Data/dist.rds",
            "Data_Analysis/Results")
}
if (!length(args) %in% c(3L, 4L)) {
  stop("Usage: Rscript Data_Analysis/Main.R DATA_DIR DIST_FILE OUTPUT_DIR [REGION]")
}
DATA_DIR <- path.expand(args[1L])
DIST_FILE <- path.expand(args[2L])
OUT_DIR <- path.expand(args[3L])
if (!file.exists(DIST_FILE)) stop("Distance matrix not found: ", DIST_FILE)
for (pkg in c("data.table", "sf", "spdep", "Matrix")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Package not installed: ", pkg)
}

# Preserve the study's regional definitions and the agreed Nebraska exclusion.
# A common eligibility rule is used everywhere; recorded zeros remain eligible.
# The USA is fitted jointly with its own graph, not pooled from regional fits.
CONUS_FIPS <- c("01", "04", "05", "06", "08", "09", "10", "11", "12", "13",
                "16", "17", "18", "19", "20", "21", "22", "23", "24", "25",
                "26", "27", "28", "29", "30", "31", "32", "33", "34", "35",
                "36", "37", "38", "39", "40", "41", "42", "44", "45", "46",
                "47", "48", "49", "50", "51", "53", "54", "55", "56")
NEW_ENGLAND_FIPS <- c("09", "23", "25", "33", "44", "50")
REGIONS <- list(
  USA = setdiff(CONUS_FIPS, NEW_ENGLAND_FIPS),
  Midwest = c("17", "18", "19", "20", "26", "27", "29", "31", "38", "39", "46", "55"),
  Northeast = c("09", "10", "23", "24", "25", "33", "34", "36", "42", "44", "50"),
  Southeast = c("01", "05", "12", "13", "21", "22", "28", "37", "45", "47", "51", "54"),
  West = c("04", "06", "08", "16", "30", "32", "35", "40", "41", "48", "49", "53", "56")
)
EXCLUDED_STATE_FIPS <- "31" # Nebraska: unavailable/frozen mortality in this source
REGIONS <- lapply(REGIONS, setdiff, y = EXCLUDED_STATE_FIPS)
REGION <- if (length(args) == 4L) args[4L] else "all"
if (tolower(REGION) != "all") {
  idx <- match(tolower(REGION), tolower(names(REGIONS)))
  if (is.na(idx)) stop("REGION must be all, USA, Midwest, Northeast, Southeast, or West.")
  REGIONS <- REGIONS[idx]
}
MOBILITY_START <- as.Date("2020-01-01")
MOBILITY_END <- as.Date("2021-05-31")
MORTALITY_START <- as.Date("2021-06-01")
MORTALITY_END <- as.Date("2022-02-01")
DECAY_DENOMINATOR <- 1000
MAX_DISTANCE_DROP_FRACTION <- 0.05
REFERENCE <- c(own = 0.1, neighbor = 0.2)
MAIN_TARGET <- c(own = 0.8, neighbor = 0.9)
ALTERNATIVE_TARGET <- c(own = 0.2, neighbor = 0.8)

fips5 <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.0+$", "", x)
  z <- suppressWarnings(as.integer(x))
  out <- rep(NA_character_, length(z))
  good <- !is.na(z) & z >= 0L & z <= 99999L
  out[good] <- sprintf("%05d", z[good])
  out
}

require_unique <- function(dt, id = "FIPS", label = "table") {
  if (anyNA(dt[[id]]) || anyDuplicated(dt[[id]])) {
    stop(label, " has missing or duplicate FIPS values; resolve before fitting.")
  }
}

# ---- 1. Read source tables and build complete county data ----
cat("Reading the county source files...\n")
# Read the original CSV or its losslessly compressed copy with the same parser.
read_survey_csv <- function(path) {
  if (file.exists(path)) return(data.table::fread(path))
  compressed <- paste0(path, ".gz")
  if (!file.exists(compressed)) stop("Survey file not found: ", path, " or ", compressed)

  unpacked <- tempfile(fileext = ".csv")
  on.exit(unlink(unpacked), add = TRUE)
  input <- gzfile(compressed, "rb")
  on.exit(close(input), add = TRUE)
  output <- file(unpacked, "wb")
  tryCatch({
    repeat {
      chunk <- readBin(input, what = "raw", n = 1024L * 1024L)
      if (!length(chunk)) break
      writeBin(chunk, output)
    }
  }, finally = close(output))
  data.table::fread(unpacked)
}

survey_raw <- read_survey_csv(file.path(DATA_DIR, "VaccineHesitancyData.csv"))
election_raw <- data.table::fread(file.path(DATA_DIR, "VaccineHesitancyData_v2.csv"))
# fread() normally preserves spaces in header names, whereas the original
# read.csv() script converted spaces to periods (via check.names = TRUE).
# Retain the exact original downstream field names by enabling that behavior.
trip_raw <- data.table::fread(
  file.path(DATA_DIR, "trip_mortality_county_weekly_data.csv"),
  check.names = TRUE
)

survey_fields <- c("FIPS Code", "Social Vulnerability Index (SVI)",
                   "Estimated hesitant", "Percent non-Hispanic White")
if (!all(survey_fields %in% names(survey_raw))) {
  stop("Missing survey columns: ",
       paste(setdiff(survey_fields, names(survey_raw)), collapse = ", "))
}
survey <- data.table::data.table(
  FIPS = fips5(survey_raw[["FIPS Code"]]),
  C3 = as.numeric(survey_raw[["Social Vulnerability Index (SVI)"]]),
  M_raw = as.numeric(survey_raw[["Estimated hesitant"]]),
  C1 = as.numeric(survey_raw[["Percent non-Hispanic White"]])
)
require_unique(survey, label = "Vaccine hesitancy survey")

# Original code renamed the second column of the election file as FIPS.
if (!("FIPS" %in% names(election_raw))) {
  data.table::setnames(election_raw, old = names(election_raw)[2], new = "FIPS")
}
if (!("Prop_Biden" %in% names(election_raw))) {
  stop("Expected election column Prop_Biden was not found.")
}
election <- data.table::data.table(
  FIPS = fips5(election_raw$FIPS), A = as.numeric(election_raw$Prop_Biden)
)
require_unique(election, label = "Election data")

trip_fields <- c("County.FIPS", "date", "Population.not.staying.at.Home",
                 "Total.weekly.population", "weekly.deaths.adj")
if (!all(trip_fields %in% names(trip_raw))) {
  stop("Missing mobility/mortality columns after normalizing CSV headers: ",
       paste(setdiff(trip_fields, names(trip_raw)), collapse = ", "),
       "\nDetected columns: ", paste(names(trip_raw), collapse = ", "),
       "\nPlease check the input CSV header; do not substitute columns by guesswork.")
}
trip_raw[, FIPS := fips5(County.FIPS)]
trip_raw[, date := as.Date(date)]
if (anyNA(trip_raw$date)) stop("Some trip dates cannot be interpreted as Date.")
trip_raw <- trip_raw[!is.na(FIPS)]

# Keep EXACTLY the same rate definitions as the original code, but close
# the previously open-ended mortality window at February 1, 2022.
trip_raw[, in_mob := date >= MOBILITY_START & date <= MOBILITY_END]
trip_raw[, in_death := date >= MORTALITY_START & date <= MORTALITY_END]
trip_raw[, mob_num := ifelse(in_mob, Population.not.staying.at.Home, NA_real_)]
trip_raw[, mob_den := ifelse(in_mob, Total.weekly.population, NA_real_)]
trip_raw[, deaths := ifelse(in_death, weekly.deaths.adj, NA_real_)]
trip_raw[, death_population := ifelse(in_death, Total.weekly.population, NA_real_)]

mobility <- trip_raw[, .(
  C2 = if (all(is.na(mob_num)) || all(is.na(mob_den))) NA_real_
  else sum(mob_num, na.rm = TRUE) / sum(mob_den, na.rm = TRUE),
  Y_raw = if (all(is.na(deaths)) || all(is.na(death_population))) NA_real_
  else sum(deaths, na.rm = TRUE) / mean(death_population, na.rm = TRUE)
), by = FIPS]
require_unique(mobility, label = "Aggregated mobility/mortality")

cat("Loading 2020 Census county polygons and matching inputs by FIPS...\n")
shapefile <- file.path(DATA_DIR, "county_boundaries", "cb_2020_us_county_500k.shp")
counties <- sf::st_read(shapefile, quiet = TRUE)
counties <- counties[counties$STATEFP %in% unique(unlist(REGIONS)), ]
counties$FIPS <- paste0(counties$STATEFP, counties$COUNTYFP)
require_unique(counties, label = "County polygons")
i_survey <- match(counties$FIPS, survey$FIPS)
i_elect <- match(counties$FIPS, election$FIPS)
i_mob <- match(counties$FIPS, mobility$FIPS)
counties$C3 <- survey$C3[i_survey]
counties$M_raw <- survey$M_raw[i_survey]
counties$C1 <- survey$C1[i_survey]
counties$A <- election$A[i_elect]
counties$C2 <- mobility$C2[i_mob]
counties$Y_raw <- mobility$Y_raw[i_mob]

D <- as.matrix(readRDS(DIST_FILE))
if (!is.numeric(D) || nrow(D) != ncol(D) || is.null(rownames(D)) || is.null(colnames(D))) {
  stop("dist.rds must be a square numeric matrix with FIPS row and column names.")
}
row_fips <- fips5(rownames(D))
col_fips <- fips5(colnames(D))
if (anyNA(row_fips) || anyNA(col_fips) || anyDuplicated(row_fips) ||
    anyDuplicated(col_fips) || !setequal(row_fips, col_fips)) {
  stop("Invalid or mismatched distance-matrix FIPS.")
}

# Freeze the cohort and Queen links once per geography. Both weightings
# require the same adjacent-pair distances, so their samples are identical.
prepare_graph <- function(states, region) {
  g <- counties[counties$STATEFP %in% states, ]
  eligible <- with(sf::st_drop_geometry(g), is.finite(A) & is.finite(C1) &
                     is.finite(C2) & is.finite(C3) & is.finite(M_raw) & is.finite(Y_raw) &
                     M_raw >= 0 & Y_raw >= 0 & FIPS %in% row_fips)
  audit <- data.frame(Region = region, FIPS = g$FIPS, State_FIPS = g$STATEFP,
                      County = g$NAME, Status = ifelse(eligible, "Retained", "Incomplete data or no distance ID"))
  g <- g[eligible, ]
  initial_n <- nrow(g)
  removed <- character()
  repeat {
    n <- nrow(g)
    if (n < 30L) stop(region, ": fewer than 30 counties remain.")
    nb <- spdep::poly2nb(g, queen = TRUE)
    nb <- lapply(nb, function(x) x[x > 0L])
    i <- rep.int(seq_len(n), lengths(nb))
    j <- unlist(nb, use.names = FALSE)
    if (!length(i)) stop(region, ": no Queen links found.")
    ri <- match(g$FIPS, row_fips)
    ci <- match(g$FIPS, col_fips)
    dij <- D[cbind(ri[i], ci[j])]
    dji <- D[cbind(ri[j], ci[i])]
    both <- is.finite(dij) & is.finite(dji)
    if (any(both & abs(dij - dji) > 1e-8 * pmax(1, abs(dij), abs(dji)))) {
      stop(region, ": conflicting distances in adjacent county pairs.")
    }
    reverse_ok <- !is.finite(dij) & is.finite(dji)
    dij[reverse_ok] <- dji[reverse_ok]
    bad <- !is.finite(dij) | dij < 0
    if (!any(bad)) break
    bad_pairs <- unique(cbind(pmin(i[bad], j[bad]), pmax(i[bad], j[bad])))
    degree <- tabulate(as.integer(bad_pairs), nbins = n)
    if ((length(removed) + 1L) / initial_n > MAX_DISTANCE_DROP_FRACTION) {
      stop(region, ": missing adjacent distances would remove more than 5% of counties.")
    }
    k <- which.max(degree)
    removed <- c(removed, g$FIPS[k])
    g <- g[-k, ]
  }
  audit$Status[audit$FIPS %in% removed] <- "Unavailable adjacent distance"
  dat <- sf::st_drop_geometry(g)
  dat$Y <- sqrt(100 * dat$Y_raw)
  dat$M <- sqrt(dat$M_raw)
  cat(sprintf("%s: %d counties, %d recorded zero outcomes, %d Queen links.\n",
              region, n, sum(dat$Y == 0), length(i)/2))
  list(data = dat, i = i, j = j, distance = dij, audit = audit, removed = removed)
}

fit_profile_gls <- function(y, X, U, lambda, label) {
  n <- length(y)
  k <- ncol(X)
  Z <- crossprod(U, X)
  z <- drop(crossprod(U, y))
  
  fit_at_ratio <- function(r) {
    den <- 1 + r * lambda
    sw <- 1 / sqrt(den)
    Zw <- Z * sw
    zw <- z * sw
    H <- crossprod(Zw)
    R <- tryCatch(chol(H), error = function(e) NULL)
    if (is.null(R)) return(list(objective = Inf))
    b <- drop(backsolve(R, forwardsolve(t(R), crossprod(Zw, zw))))
    residual <- drop(z - Z %*% b)
    v <- sum(residual^2 / den) / n
    if (!is.finite(v) || v <= 0) return(list(objective = Inf))
    list(objective = 0.5 * (n * log(v) + sum(log1p(r * lambda))),
         coef = b, v = v, w = r * v, r = r, R = R)
  }
  
  # Coarse scan plus local refinement; explicitly include r=0 boundary.
  log_grid <- seq(log(1e-8), log(1e5), length.out = 37)
  grid_obj <- vapply(log_grid, function(t) fit_at_ratio(exp(t))$objective, numeric(1))
  idx <- which.min(grid_obj)
  lo <- log_grid[max(1L, idx - 1L)]
  hi <- log_grid[min(length(log_grid), idx + 1L)]
  opt <- optimize(function(t) fit_at_ratio(exp(t))$objective,
                  interval = c(lo, hi), tol = 1e-7)
  candidates <- c(0, exp(log_grid[idx]), exp(opt$minimum))
  obj <- vapply(candidates, function(r) fit_at_ratio(r)$objective, numeric(1))
  best <- fit_at_ratio(candidates[which.min(obj)])
  if (!is.finite(best$objective)) stop(label, ": spatial Gaussian MLE failed.")
  
  # Covariance of regression coefficients from expected information:
  # Var(beta_hat | X) = (X' Sigma^{-1} X)^{-1}.
  vcov_b <- best$v * chol2inv(best$R)
  dimnames(vcov_b) <- list(colnames(X), colnames(X))
  names(best$coef) <- colnames(X)
  
  # Expected information for residual and random-effect variances (v,w).
  # At w=0, SD(random effect) inference is nonregular and set to NA.
  vd <- best$v + best$w * lambda
  Q <- cbind(1 / vd, lambda / vd)
  I_vw <- 0.5 * crossprod(Q)
  inv_I <- tryCatch(chol2inv(chol(I_vw)), error = function(e) NULL)
  sd_vcov <- matrix(NA_real_, 2, 2)
  if (!is.null(inv_I) && best$v > 0 && best$w > 0) {
    D_sd <- diag(c(1 / (2 * sqrt(best$v)), 1 / (2 * sqrt(best$w))))
    sd_vcov <- D_sd %*% inv_I %*% D_sd
  }
  
  cat(sprintf("%s: r=w/v=%.5g; sigma=%.5g; sigma_b=%.5g\n",
              label, best$r, sqrt(best$v), sqrt(best$w)))
  list(coef = best$coef, vcov = vcov_b,
       sigma = sqrt(best$v), sigma_b = sqrt(best$w),
       var_sigma = diag(sd_vcov), vcov_sigma = sd_vcov,
       r = best$r, nll_profile = best$objective)
}

fit_graph <- function(g, cross_state_multiplier, label, adjustments) {
  dat <- g$data
  n <- nrow(dat)
  same <- dat$STATEFP[g$i] == dat$STATEFP[g$j]
  edge_w <- exp(-g$distance / DECAY_DENOMINATOR) * ifelse(same, 1, cross_state_multiplier)
  E <- Matrix::sparseMatrix(i = c(seq_len(n), g$i), j = c(seq_len(n), g$j),
                            x = c(rep(1, n), edge_w), dims = c(n, n))
  E <- Matrix::forceSymmetric(E, uplo = "U")
  E0 <- E
  diag(E0) <- 0
  strength <- as.numeric(Matrix::rowSums(E0))
  active <- strength > 0
  W <- Matrix::Diagonal(x = ifelse(active, 1 / pmax(strength, .Machine$double.eps), 0)) %*% E0
  stopifnot(all(abs(as.numeric(Matrix::rowSums(W)) - as.numeric(active)) < 1e-12))
  W2 <- W %*% W
  rho <- sum(Matrix::diag(W2)) / n
  den <- as.numeric(Matrix::rowSums(W * W))
  num <- as.numeric(Matrix::rowSums(W2 * W))
  kappa <- numeric(n)
  kappa[den > 0] <- num[den > 0] / den[den > 0]
  lag_value <- function(x) as.numeric(W %*% x)
  X_m <- with(dat, cbind(Intercept = 1, A = A, WA = lag_value(A),
                         C1 = C1, WC1 = lag_value(C1), C2 = C2, WC2 = lag_value(C2), C3 = C3, WC3 = lag_value(C3)))
  X_y <- cbind(X_m[, 1:3, drop = FALSE], M = dat$M, WM = lag_value(dat$M), X_m[, 4:9, drop = FALSE])
  cat("Fitting ", label, " (cross-state multiplier ", cross_state_multiplier, ")...\n", sep = "")
  # The covariance basis depends on E and can be shared by both adjustments.
  eig <- eigen(as.matrix(E), symmetric = TRUE)
  lambda <- eig$values^2
  fits <- lapply(adjustments, function(adjustment) {
    use_mobility <- adjustment == "With mobility"
    omit <- c("C2", "WC2")
    dy <- if (use_mobility) X_y else X_y[, !colnames(X_y) %in% omit, drop = FALSE]
    dm <- if (use_mobility) X_m else X_m[, !colnames(X_m) %in% omit, drop = FALSE]
    fy <- fit_profile_gls(dat$Y, dy, eig$vectors, lambda, paste0("  Outcome (", adjustment, ")"))
    fm <- fit_profile_gls(dat$M, dm, eig$vectors, lambda, paste0("  Mediator (", adjustment, ")"))
    list(FIPS = dat$FIPS, E = E, W = W, X_y = dy, X_m = dm, fit_y = fy, fit_m = fm,
         h = mean(active), rho = rho, kappa_gaussian = mean(kappa), n_isolates = sum(!active),
         n_edges = length(g$i)/2, cross_state_multiplier = cross_state_multiplier)
  })
  setNames(fits, adjustments)
}

estimate_effects <- function(fit, target, region, scenario, adjustment) {
  da <- unname(target["own"] - REFERENCE["own"])
  ds <- unname(target["neighbor"] - REFERENCE["neighbor"])
  b <- fit$fit_y$coef
  g <- fit$fit_m$coef
  Vb <- fit$fit_y$vcov
  Vg <- fit$fit_m$vcov
  rho <- fit$rho
  h <- fit$h
  # Isolates stay in the fit and contribute zero to the neighbor pathways.
  ds_active <- h * ds
  kappa_gaussian <- fit$kappa_gaussian
  # b['A'], b['WA'], b['M'], b['WM']; g['A'], g['WA'].
  point <- c(
    tau1 = da * b["A"],
    tau2 = da * b["M"] * g["A"],
    tau3 = da * rho * b["WM"] * g["WA"],
    tau4 = ds_active * b["WA"],
    tau5 = ds_active * b["M"] * g["WA"],
    tau6 = ds * b["WM"] * (h * g["A"] + kappa_gaussian * g["WA"])
  )
  point <- setNames(as.numeric(point), paste0("tau", 1:6))
  
  # Gradient with respect to y coefficients and mediator coefficients.
  Gy <- matrix(0, nrow = 6, ncol = length(b),
               dimnames = list(names(point), names(b)))
  Gm <- matrix(0, nrow = 6, ncol = length(g),
               dimnames = list(names(point), names(g)))
  Gy["tau1", "A"] <- da
  Gy["tau2", "M"] <- da * g["A"]
  Gm["tau2", "A"] <- da * b["M"]
  Gy["tau3", "WM"] <- da * rho * g["WA"]
  Gm["tau3", "WA"] <- da * rho * b["WM"]
  Gy["tau4", "WA"] <- ds_active
  Gy["tau5", "M"] <- ds_active * g["WA"]
  Gm["tau5", "WA"] <- ds_active * b["M"]
  Gy["tau6", "WM"] <- ds * (h * g["A"] + kappa_gaussian * g["WA"])
  Gm["tau6", "A"] <- ds * b["WM"] * h
  Gm["tau6", "WA"] <- ds * b["WM"] * kappa_gaussian
  
  # Gaussian conditional outcome and mediator scores are orthogonal.
  # Combined first-order delta covariance; this retains the gamma_A,
  # gamma_WA covariance for tau6 (omitted by many shortcuts).
  Vtau <- Gy %*% Vb %*% t(Gy) + Gm %*% Vg %*% t(Gm)
  se <- sqrt(pmax(diag(Vtau), 0))
  z <- qnorm(0.975)
  results <- data.frame(Region = region, Mobility_Adjustment = adjustment,
                        Scenario = scenario, N = length(fit$FIPS),
                        Effect = names(point), Estimate = unname(point), SE = unname(se),
                        Lower_95 = unname(point - z * se), Upper_95 = unname(point + z * se),
                        Own_Reference = unname(REFERENCE["own"]), Neighbor_Reference = unname(REFERENCE["neighbor"]),
                        Own_Target = unname(target["own"]), Neighbor_Target = unname(target["neighbor"]),
                        Cross_State_Multiplier = fit$cross_state_multiplier, Distance_Decay_Denominator = DECAY_DENOMINATOR)
  list(results = results, effect_vcov = Vtau)
}

parameters <- function(fit, region, weighting, adjustment) {
  f <- list(Outcome = fit$fit_y, Mediator = fit$fit_m)
  do.call(rbind, lapply(names(f), function(model) {
    p <- f[[model]]
    data.frame(Region = region, Mobility_Adjustment = adjustment, Weighting = weighting, Model = model,
               Term = names(p$coef), Estimate = unname(p$coef), SE = sqrt(unname(diag(p$vcov))))
  }))
}

all_runs <- lapply(names(REGIONS), function(region) {
  g <- prepare_graph(REGIONS[[region]], region)
  fits <- list(
    State_weighted = fit_graph(g, 0.5, paste(region, "state weighted"),
                               c("With mobility", "Without mobility")),
    Binary_B = fit_graph(g, 1, paste(region, "binary B"), "With mobility")
  )
  main <- fits$State_weighted[["With mobility"]]
  binary <- fits$Binary_B[["With mobility"]]
  no_mobility <- fits$State_weighted[["Without mobility"]]
  stopifnot(identical(main$FIPS, g$data$FIPS), identical(binary$FIPS, g$data$FIPS),
            identical(no_mobility$FIPS, g$data$FIPS))
  effect_sets <- list(
    main = estimate_effects(main, MAIN_TARGET, region, "Main", "With mobility"),
    contrast = estimate_effects(main, ALTERNATIVE_TARGET, region, "Alternative contrast", "With mobility"),
    binary = estimate_effects(binary, MAIN_TARGET, region, "Binary B (distance decay retained)", "With mobility"),
    no_mobility = estimate_effects(no_mobility, MAIN_TARGET, region, "Without mobility", "Without mobility")
  )
  scale <- c(rep(1/7, 3), rep(6/7, 3))
  ci_cols <- c("Estimate", "SE", "Lower_95", "Upper_95")
  a <- as.matrix(effect_sets$main$results[ci_cols])
  b <- as.matrix(effect_sets$contrast$results[ci_cols])
  stopifnot(max(abs(b - scale * a)) < 1e-12)
  
  diagnostics <- do.call(rbind, lapply(names(fits), function(weighting) {
    do.call(rbind, lapply(names(fits[[weighting]]), function(adjustment) {
      f <- fits[[weighting]][[adjustment]]
      data.frame(Region = region, Mobility_Adjustment = adjustment, Weighting = weighting,
                 N = nrow(g$data), Recorded_Zero_Outcomes = sum(g$data$Y == 0),
                 Unavailable_Distance_Exclusions = length(g$removed),
                 Cross_State_Multiplier = f$cross_state_multiplier, Isolates = f$n_isolates,
                 Edges = f$n_edges, rho = f$rho, kappa_gaussian = f$kappa_gaussian,
                 Nonisolate_Fraction = f$h, Outcome_Sigma = f$fit_y$sigma, Outcome_Sigma_B = f$fit_y$sigma_b,
                 Mediator_Sigma = f$fit_m$sigma, Mediator_Sigma_B = f$fit_m$sigma_b)
    }))
  }))
  coefficients <- do.call(rbind, lapply(names(fits), function(weighting) {
    do.call(rbind, lapply(names(fits[[weighting]]), function(adjustment) {
      parameters(fits[[weighting]][[adjustment]], region, weighting, adjustment)
    }))
  }))
  list(fits = fits, data = g$data, audit = g$audit, effects = effect_sets,
       diagnostics = diagnostics, parameters = coefficients)
})
names(all_runs) <- names(REGIONS)

collect <- function(field) do.call(rbind, lapply(all_runs, `[[`, field))
all_effects <- do.call(rbind, lapply(all_runs, function(x) {
  do.call(rbind, lapply(x$effects, `[[`, "results"))
}))
stopifnot(nrow(all_effects) == length(REGIONS) * 24L,
          all(is.finite(as.matrix(all_effects[c("Estimate", "SE", "Lower_95", "Upper_95")]))))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(all_effects, file.path(OUT_DIR, "effects_all.csv"), row.names = FALSE)
write.csv(collect("diagnostics"), file.path(OUT_DIR, "fit_diagnostics.csv"), row.names = FALSE)
write.csv(collect("parameters"), file.path(OUT_DIR, "model_coefficients.csv"), row.names = FALSE)
write.csv(collect("audit"), file.path(OUT_DIR, "county_inclusion.csv"), row.names = FALSE)
saveRDS(all_runs, file.path(OUT_DIR, "application_fits.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "R_session_info.txt"))
cat(sprintf("Completed %d geographic samples, %d model-fit pairs, and %d effect estimates.\n",
            length(REGIONS), 3L * length(REGIONS), nrow(all_effects)))
