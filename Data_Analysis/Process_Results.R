# From the repository root: Rscript Data_Analysis/Process_Results.R
# Optional arguments: FITS_RDS TABLE_DIR
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) {
  stop("Usage: Rscript Data_Analysis/Process_Results.R [FITS_RDS [TABLE_DIR]]")
}
fit_file <- if (length(args) >= 1L) args[[1L]] else "Data_Analysis/Results/application_fits.rds"
table_dir <- if (length(args) >= 2L) args[[2L]] else "Data_Analysis/Results/Reported"
fits <- readRDS(fit_file)
stopifnot(setequal(names(fits), c("USA", "Midwest", "Northeast", "Southeast", "West")))

format_value <- function(x, digits) {
  printed <- sprintf(paste0("%.", digits, "f"), x)
  near_zero <- x != 0 & abs(x) < 0.5 * 10^(-digits)
  printed[near_zero] <- sprintf("%.6f", x[near_zero])
  printed
}

effects <- do.call(rbind, lapply(fits, function(run) {
  do.call(rbind, lapply(run$effects, function(effect) effect$results))
}))
stopifnot(nrow(effects) == 120L)
effect_table <- effects[c("Region", "Scenario", "N", "Effect", "Estimate", "Lower_95", "Upper_95")]
for (field in c("Estimate", "Lower_95", "Upper_95")) {
  effect_table[[field]] <- format_value(effect_table[[field]], 4L)
}

parameters <- do.call(rbind, lapply(fits, function(run) run$parameters))
main <- parameters[parameters$Weighting == "State_weighted" &
                   parameters$Mobility_Adjustment == "With mobility", ]
labels <- c(A = "PP", M = "VH", C3 = "SVI", C1 = "Race", C2 = "Mobility")
region_order <- c("Midwest", "Northeast", "West", "Southeast", "USA")
own <- main[main$Term %in% names(labels), ]
own <- own[order(match(own$Region, region_order),
                 match(own$Model, c("Outcome", "Mediator")),
                 match(own$Term, names(labels))), ]
neighbor <- main[match(paste(own$Region, own$Model, paste0("W", own$Term)),
                       paste(main$Region, main$Model, main$Term)), ]
stopifnot(nrow(own) == 45L, all(is.finite(neighbor$Estimate)))
z <- qnorm(0.975)
coefficient_table <- data.frame(
  Region = own$Region, Model = own$Model, Covariate = unname(labels[own$Term]),
  Endogenous_Estimate = format_value(own$Estimate, 3L),
  Endogenous_Lower_95 = format_value(own$Estimate - z * own$SE, 3L),
  Endogenous_Upper_95 = format_value(own$Estimate + z * own$SE, 3L),
  Exogenous_Estimate = format_value(neighbor$Estimate, 3L),
  Exogenous_Lower_95 = format_value(neighbor$Estimate - z * neighbor$SE, 3L),
  Exogenous_Upper_95 = format_value(neighbor$Estimate + z * neighbor$SE, 3L)
)

diagnostics <- do.call(rbind, lapply(fits, function(run) run$diagnostics))
sample_counts <- diagnostics[diagnostics$Weighting == "State_weighted" &
                             diagnostics$Mobility_Adjustment == "With mobility",
                             c("Region", "N", "Recorded_Zero_Outcomes",
                               "Unavailable_Distance_Exclusions", "Isolates")]

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(effect_table, file.path(table_dir, "pathway_estimates_table.csv"), row.names = FALSE)
write.csv(coefficient_table, file.path(table_dir, "main_coefficient_table.csv"), row.names = FALSE)
write.csv(sample_counts, file.path(table_dir, "sample_counts.csv"), row.names = FALSE)
