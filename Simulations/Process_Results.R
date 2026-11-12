# From the repository root: Rscript Simulations/Process_Results.R
# Optional arguments: RESULTS_DIR TABLE_DIR
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) {
  stop("Usage: Rscript Simulations/Process_Results.R [RESULTS_DIR [TABLE_DIR]]")
}
results_dir <- if (length(args) >= 1L) args[[1L]] else "Simulations/Results"
table_dir <- if (length(args) >= 2L) args[[2L]] else file.path(results_dir, "Reported")
source("Simulations/Simulation_Summaries.R")

effect_tables <- parameter_tables <- list()
for (n in c(200L, 500L)) {
  name <- paste0("New_simu_", n, "_gaussian_primary.rds")
  run <- readRDS(file.path(results_dir, name))
  config <- run$config
  truth <- with(config, c(beta, sigma_y, sigma_by, gamma, sigma_m, sigma_bm))
  effects <- summarize_effects(run$effects, run$success)
  parameters <- summarize_parameters(run$mle, run$var_par, run$success, truth)

  effect_tables[[as.character(n)]] <- data.frame(
    N = n, Effect = effects$effect,
    True = sprintf("%.4f", effects$truth),
    Bias = sprintf("%.4f", effects$bias),
    RRMSE = sprintf("%.4f", effects$RRMSE),
    ESE = sprintf("%.4f", effects$ESE),
    ASE = sprintf("%.4f", effects$ASE),
    CP = sprintf("%.3f", effects$coverage)
  )
  parameter_tables[[as.character(n)]] <- parameters
}

params200 <- parameter_tables[["200"]]
params500 <- parameter_tables[["500"]]
stopifnot(identical(params200$parameter, params500$parameter),
          identical(params200$truth, params500$truth))
parameter_table <- data.frame(
  Parameter = params200$parameter,
  True = sprintf("%.3f", params200$truth),
  Mean_N200 = sprintf("%.3f", params200$mean_est),
  Mean_N500 = sprintf("%.3f", params500$mean_est)
)

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, effect_tables), file.path(table_dir, "simulation_effects_table.csv"),
          row.names = FALSE)
write.csv(parameter_table, file.path(table_dir, "simulation_parameters_table.csv"),
          row.names = FALSE)
