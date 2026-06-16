# ============================================================================
# 07_ipm_run.R  —  IPM v3.2 RUN FILE
# ----------------------------------------------------------------------------
# Driver for fitting and post-processing the Bayesian IPM. Run blocks
# individually (not source()'d top to bottom). Object shape is consistent
# throughout: every result object exposes $fit and $sim_data.
#
# Pipeline assumed already run for real data:
#   source("Code/01_data_prep_6site.R")   # -> dat (18 x 29), years
#   source("Code/02_covariates_6site.R")  # -> cov_t_scaled (16 x 29)
#
# Stan tuning lives in the run_full_analysis_v3.2() call below (Block 2).
# ============================================================================


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║ BLOCK 0 — STARTUP                                                          ║
# ╚══════════════════════════════════════════════════════════════════════════╝

source("Code/00_load.R")                       # packages + dirs + banner
source("Code/05_ipm_model.R")             # Stan model + functions + orchestrator
source("Code/06_ipm_plots.R")             # plotting + load_seal_results()
filter <- dplyr::filter                   # guard against stats::filter masking


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║ BLOCK 1 — BUILD / RELOAD IPM INPUTS                                        ║
# ╚══════════════════════════════════════════════════════════════════════════╝
# Either (a) build inputs fresh from the prep scripts, or (b) reload a saved
# bundle. Inputs are: dat (18x29), cov_t_scaled (16x29), years.

# (a) Build fresh:
# source("Code/01_data_prep_6site.R")
# source("Code/02_covariates_6site.R")
# saveRDS(list(dat = dat, cov_t_scaled = cov_t_scaled, years = years),
#         "Output/ipm_input_data.rds")

# (b) Reload a saved bundle:
# input_data   <- readRDS("Output/ipm_input_data.rds")
# dat          <- input_data$dat
# cov_t_scaled <- input_data$cov_t_scaled
# years        <- input_data$years

stopifnot(exists("dat"), exists("cov_t_scaled"), exists("years"))
cat(sprintf("Inputs ready: dat %dx%d | cov %dx%d | years %d-%d\n",
            nrow(dat), ncol(dat), nrow(cov_t_scaled), ncol(cov_t_scaled),
            min(years), max(years)))


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║ BLOCK 2 — FIT REAL DATA                                                    ║
# ╚══════════════════════════════════════════════════════════════════════════╝
# STAN TUNING is set here. Defaults below (2000/2000, adapt_delta 0.97,
# treedepth 12) reproduce the 8000-draw fit reported in the manuscript.
# run_full_analysis_v3.2() saves the fit + input bundle and runs all plots.

out <- run_full_analysis_v3.2(
  use_real_data = TRUE,
  dat           = dat,
  cov_t_scaled  = cov_t_scaled,
  years         = years,
  seed          = 123,
  iter_warmup   = 2000,
  iter_sampling = 2000,
  adapt_delta   = 0.97,
  max_treedepth = 12,
  run_portfolio = TRUE,
  run_synchrony = TRUE
)
# out$fit       -> CmdStanFit
# out$sim_data  -> stan_data + metadata (NOT $data)
# out$portfolio, out$sync, out$forest, ... -> plot/table results

cat("Real-data fit complete:", format(Sys.time()), "\n")

# Quick sanity check
out$fit$summary("phi_juv_base")[, c("variable", "mean")]   # expect ~0.70


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║ BLOCK 3 — RELOAD + REPLOT (no re-sampling)                                 ║
# ╚══════════════════════════════════════════════════════════════════════════╝
# Use in later sessions to regenerate plots/tables from the saved fit.

out <- load_seal_results("IPM_v3.2_real")   # -> list(fit, sim_data)
filter <- dplyr::filter

cat("fit class:", class(out$fit)[1], "| draws() works:",
    is.function(out$fit$draws), "\n")

run_all_plots_v3.2(
  fit           = out$fit,
  sim_data      = out$sim_data,            # consistent shape ($sim_data)
  prefix        = "IPM_v3.2_real",
  run_portfolio = TRUE,
  run_synchrony = TRUE
)

# Single plot in the RStudio pane, e.g.:
# sa <- create_site_age_timeseries_v3.2(out$fit, out$sim_data,
#                                       save = TRUE, prefix = "IPM_v3.2_real")
# sa$by_age


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║ BLOCK 4 — SIMULATION FIT + PARAMETER RECOVERY                              ║
# ╚══════════════════════════════════════════════════════════════════════════╝

source("Code/00_load.R")
source("Code/05_ipm_model.R")      # writes the .stan, defines functions (with the edit)
source("Code/05b_recovery_sbc.R")
filter <- dplyr::filter


sbc <- run_recovery_check_v3.2(n_datasets = 10, base_seed = 1000)"
)

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║ BLOCK 5 — HOUSEKEEPING                                                     ║
# ╚══════════════════════════════════════════════════════════════════════════╝

# Delete CmdStan temp CSVs after the fit object is saved (frees disk).
# out$fit$output_files()          # show paths
# unlink(out$fit$output_files())  # delete

# Free memory between heavy steps.
# rm(list = setdiff(ls(), c("out", "filter"))); gc(); gc()

# Commit.
# system('git add -A')
# system('git commit -m "IPM v3.2 rerun"')
# system('git push origin main')