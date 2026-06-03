# SAVE (before restart)
saveRDS((.packages()), "session_packages.rds")

# LOAD (after restart)
invisible(lapply(readRDS("session_packages.rds"), library, character.only = TRUE))

# if just running post_model analyses/plots
results.real <- readRDS("Output/harbor_seal_IPM_v3.1_real_fit.rds")


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



fit <- readRDS("Output/harbor_seal_IPM_v3.1_real_fit.rds")




#
## run sim

library(bayesplot)
results.sim <- run_full_analysis_v3.1(use_real_data = FALSE, #shorter runs and lower adapt delta
                                      iter_warmup = 1000,
                                      iter_sampling = 500,
                                      adapt_delta = 0.95,
                                      seed = 42)
saveRDS(results.sim, "Output/harbor_seal_IPM_v3.1_sim_results.rds")

## run real

library(bayesplot) 
results.real <- run_full_analysis_v3.1(
  use_real_data = TRUE,
  dat = dat,
  cov_t_scaled = cov_t_scaled,
  years = years,
  iter_warmup = 2000,
  iter_sampling = 1000,
  adapt_delta = 0.995,
  max_treedepth = 15)







# Use fit directly for all functions:
create_effect_plots_v3.1(fit, prefix = "IPM_v3.1_real")


sync_analysis <- create_synchrony_projections_v3.1(fit, sim_data = results.real$data, prefix = "IPM_v3.1_real")




# Now run portfolio analysis
portfolio <- create_portfolio_analysis_v3.1(fit, sim_data = sim_data, prefix = "IPM_v3.1_real")

# Rebuild the full sim_data structure
data_list <- prepare_real_data_for_ipm_v3.1(dat, cov_t_scaled, years, T_proj = 10)
sim_data <- list(
  stan_data = data_list$stan_data,
  site_names = data_list$site_names,
  years = data_list$years,
  scenario_names = data_list$scenario_names
)


portfolio.real <- create_portfolio_analysis_v3.1(
  fit = fit,
  sim_data = sim_data,
  prefix = "IPM_v3.1_real"
)


# And synchrony analysis
sync_analysis <- create_synchrony_projections_v3.1(fit, sim_data = sim_data, prefix = "IPM_v3.1_real")

