##2026-06-08 startup process


# 1. ONLY FOR NEW CONTAINERS   If packages missing after container restart:
source("restore_packages.R")   # restores + restarts automatically


# everytime
source("harbor_seal_setup.R")  # sets up git, CmdStan, .Rprofile

# 2. After restart, normal workflow:
source("00_load.R")
source("Code/harbor_seal_ipm_v3.2.R")
source("Code/harbor_seal_ipm_v3.2_plots.R")
out <- load_seal_results("IPM_v3.2_real")
filter <- dplyr::filter   # prevent stats::filter masking


source("restore_packages.R")   # installs packages + restarts R automatically
# after restart:
# 
# #only if cmdstan missing
# source("harbor_seal_setup.R")  # git config, CmdStan path, .Rprofile
# # then normal startup above


# LOAD (after restart)
invisible(lapply(readRDS("session_packages.rds"), library, character.only = TRUE))

# if just running post_model analyses/plots
Results.real <- readRDS("Output/harbor_seal_IPM_v3.2_real_fit.rds")


# Save IPM input data
# saveRDS(
#   list(dat = dat, cov_t_scaled = cov_t_scaled, years = years),
#   "Output/ipm_input_data.rds"
# )

# Reload
input_data <- readRDS("Output/ipm_input_data.rds")
dat <- input_data$dat
cov_t_scaled <- input_data$cov_t_scaled
years <- input_data$years



fit <- readRDS("Output/harbor_seal_IPM_v3.2_real_fit.rds")

#
## run sim

library(bayesplot)
results.sim <- run_full_analysis_v3.2(use_real_data = FALSE, #shorter runs and lower adapt delta
                                      iter_warmup = 1000,
                                      iter_sampling = 500,
                                      adapt_delta = 0.95,
                                      seed = 42)
saveRDS(results.sim, "Output/harbor_seal_IPM_v3.2_sim_results.rds")

out.sim <- load_seal_results("IPM_v3.2_sim")
out.sim$sim_data$true_params   # should now have values


## run real

library(bayesplot) 
Results.real <- run_full_analysis_v3.2(
use_real_data=TRUE, 
dat=dat, 
cov_t_scaled=cov_t_scaled, 
years=years,
seed = 123,
iter_warmup=3000, 
iter_sampling=1000, 
adapt_delta=0.98,
max_treedepth = 12,
store_warmup    = FALSE,    # explicit — don't store warmup draws
save_warmup     = 0,        # Stan-level: don't write warmup to CSV
refresh         = 200)


saveRDS(Results.real, "Output/harbor_seal_IPM_v3.2_real_results.rds")

# After saving fit object, delete temp CSV files
fit$save_object("Output/harbor_seal_IPM_v3.2_real_fit.rds")
fit$output_files()    # shows where the CSVs are
unlink(fit$output_files())    # deletes them


# Use fit directly for all functions:
create_effect_plots_v3.2(fit, prefix = "IPM_v3.2_real")



# ── Reload the fit using load_seal_results() ──────────────────────────────
out <- load_seal_results("IPM_v3.2_real")  # reads from Output/RDS files

# Quick check before running plots
cat("fit class:", class(out$fit), "\n")          # should show CmdStanFit
cat("draws method works:", is.function(out$fit$draws), "\n")  # should be TRUE






## run it all 

# ── Now run plots ─────────────────────────────────────────────────────────
source("Code/harbor_seal_ipm_v3.2.R")
source("Code/harbor_seal_ipm_v3.2_plots.R")

out <- load_seal_results("IPM_v3.2_real")
filter <- dplyr::filter

run_all_plots_v3.2(
  fit           = out$fit,
  sim_data      = out$sim_data,
  prefix        = "IPM_v3.2_real",
  run_portfolio = TRUE,
  run_synchrony = TRUE
)





# ── Step 1: fit model on SIMULATED data ──────────────────────────────────────
out.sim <- load_seal_results("IPM_v3.2_sim")
names(out.sim$sim_data)
out.sim$sim_data$true_params




results.sim <- run_full_analysis_v3.2(
  use_real_data  = FALSE,        # generates fake data from true_params
  iter_warmup    = 2000,
  iter_sampling  = 2000,
  adapt_delta    = 0.97
)

# ── Step 2: check recovery against known true values ─────────────────────────
rec <- check_parameter_recovery_v3.2(
  fit      = results.sim$fit,
  sim_data = results.sim$sim_data,   # contains true_params
  save     = TRUE,
  prefix   = "IPM_v3.2_sim"
)

# ── Step 3: separately reload and plot real data results ─────────────────────
out <- load_seal_results("IPM_v3.2_real")
filter <- dplyr::filter
run_all_plots_v3.2(
  fit           = out$fit,
  sim_data      = out$sim_data,
  prefix        = "IPM_v3.2_real",
  run_portfolio = TRUE,
  run_synchrony = TRUE
)







#before close: 
 # ── Clean up after heavy model processing ────────────────────────────────────
  # Remove large intermediate objects that accumulate during plotting
  rm(list = setdiff(ls(), c("out", "filter")))   # keep only what you need

# Force garbage collection twice — second pass frees memory from first
gc(); gc()


# Stage, commit, push — run these three lines
system('git add -A')
system('git commit -m "commit"')
system('git push origin main')

