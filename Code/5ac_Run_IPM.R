##2026-06-08 startup process


# 1. ONLY FOR NEW CONTAINERS   If packages missing after container restart:
source("restore_packages.R")   # restores + restarts automatically


#everytime
source("harbor_seal_setup.R")  # sets up git, CmdStan, .Rprofile

# 2. After restart, normal workflow:
source("00_load.R")
source("Code/harbor_seal_ipm_v3.2.R")
source("Code/harbor_seal_ipm_v3.2_plots.R")
out <- load_seal_results("IPM_v3.2_real")





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



sync_analysis <- create_synchrony_projections_v3.2(
  fit      = fit,
  sim_data = prepared,
  prefix   = "IPM_v3.2_real"
)





#
## run sim

library(bayesplot)
results.sim <- run_full_analysis_v3.2(use_real_data = FALSE, #shorter runs and lower adapt delta
                                      iter_warmup = 1000,
                                      iter_sampling = 500,
                                      adapt_delta = 0.95,
                                      seed = 42)
saveRDS(results.sim, "Output/harbor_seal_IPM_v3.2_sim_results.rds")




## run real

library(bayesplot) 
Results.real <- run_full_analysis_v3.2(
use_real_data=TRUE, 
dat=dat, 
cov_t_scaled=cov_t_scaled, 
years=years,
seed = 123,
iter_warmup=2000, 
iter_sampling=2000, 
adapt_delta=0.97,
max_treedepth = 12)


saveRDS(Results.real, "Output/harbor_seal_IPM_v3.2_real_results.rds")




# Use fit directly for all functions:
create_effect_plots_v3.2(fit, prefix = "IPM_v3.2_real")



# ── Reload the fit using load_seal_results() ──────────────────────────────
out <- load_seal_results("IPM_v3.2_real")  # reads from Output/RDS files

# Quick check before running plots
cat("fit class:", class(out$fit), "\n")          # should show CmdStanFit
cat("draws method works:", is.function(out$fit$draws), "\n")  # should be TRUE




sync_analysis <- create_synchrony_projections_v3.2(
  fit      = out$fit,
  sim_data = out$sim_data,
  prefix = "IPM_v3.2_real")


source("Code/harbor_seal_ipm_v3.2_plots.R", local=FALSE)



# ── Now run plots ─────────────────────────────────────────────────────────
run_all_plots_v3.2(
  fit      = out$fit,
  sim_data = out$sim_data,
  prefix   = "IPM_v3.2_real"
)


portfolio.real <- create_portfolio_analysis_v3.2(
  fit      = out$fit,
  sim_data = out$sim_data,
  prefix   = "IPM_v3.2_real"
)

source("Code/harbor_seal_ipm_v3.2_plots.R", local=FALSE)
eff.real <- create_effect_plots_v3.2(
  fit      = out$fit,
  #sim_data = out$sim_data,
  prefix="IPM_v3.2_real")

source("Code/harbor_seal_ipm_v3.2_plots.R", local=FALSE)
decomp <- create_covariate_decomposition_plots_v3.2(
  fit      = out$fit,
  sim_data = out$sim_data,
  prefix   = "IPM_v3.2_real"
)


proj <- create_projection_plots_v3.2(
  fit      = out$fit,
  sim_data = out$sim_data,
  prefix   = "IPM_v3.2_real"
)




## run it all 

source("Code/harbor_seal_ipm_v3.2.R")
source("Code/harbor_seal_ipm_v3.2_plots.R")

out <- load_seal_results("IPM_v3.2_real")

run_all_plots_v3.2(
  fit      = out$fit,
  sim_data = out$sim_data,
  prefix   = "IPM_v3.2_real"
)




