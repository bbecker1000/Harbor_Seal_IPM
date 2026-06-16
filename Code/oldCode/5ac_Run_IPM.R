##2026-06-08 startup process

# 1
source("restore_packages.R")   # installs packages + restarts R automatically
# 2. After restart, normal workflow:
source("00_load.R")
source("Code/harbor_seal_ipm_v3.2.R")
source("Code/harbor_seal_ipm_v3.2_plots.R")
out <- load_seal_results("IPM_v3.2_real")
filter <- dplyr::filter   # prevent stats::filter masking



# after restart:
# 
# #only if cmdstan missing
# source("harbor_seal_setup.R")  # git config, CmdStan path, .Rprofile
# # then normal startup above



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


## run real

library(bayesplot) 
results.real <- run_full_analysis_v3.2(
  use_real_data  = TRUE,
  dat            = dat,
  cov_t_scaled   = cov_t_scaled,
  years          = years,
  seed           = 123,
  iter_warmup    = 3000,
  iter_sampling  = 1000,
  adapt_delta    = 0.97,
  max_treedepth  = 12
  # prefix         = "IPM_v3.2_real"
)

# Explicit save immediately after run completes
results.real$fit$save_object("Output/harbor_seal_IPM_v3.2_real_fit.rds")
cat("Saved at:", format(Sys.time()), "\n")


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

results.real <- out <- load_seal_results("IPM_v3.2_real")

out <- load_seal_results('IPM_v3.2_real')
run_all_plots_v3.2(out$fit, out$sim_data, prefix='IPM_v3.2_real')

filter <- dplyr::filter

run_all_plots_v3.2(
  fit           = out$fit,
  sim_data      = out$data,   # note: $data not $sim_data
  prefix        = "IPM_v3.2_real",
  run_portfolio = TRUE,
  run_synchrony = TRUE
)

# Save the new fit over the old RDS
results.real$fit$save_object("Output/harbor_seal_IPM_v3.2_real_fit.rds")

# Then future sessions can use load_seal_results normally
out <- load_seal_results("IPM_v3.2_real")
# Verify
out$fit$summary("phi_juv_base")[, c("variable","mean")]  # should be ~0.698




out <- load_seal_results("IPM_v3.2_real")

sa <- create_site_age_timeseries_v3.2(out$fit, out$sim_data, save=TRUE, prefix="IPM_v3.2_real")

sa$by_age   # view the plot in RStudio's plot pane










# ── Step 1: fit model on SIMULATED data ──────────────────────────────────────
out.sim <- load_seal_results("IPM_v3.2_sim")
names(out.sim$sim_data)
out.sim$sim_data$true_params



## run sim

library(bayesplot)
results.sim <- run_full_analysis_v3.2(use_real_data = FALSE, #shorter runs and lower adapt delta
                                      iter_warmup = 1000,
                                      iter_sampling = 500,
                                      adapt_delta = 0.9, #for speed
                                      seed = 42)
saveRDS(results.sim, "Output/harbor_seal_IPM_v3.2_sim_results.rds")

out.sim <- load_seal_results("IPM_v3.2_sim")
out.sim$sim_data$true_params   # should now have values

# ── Step 2: check recovery against known true values ─────────────────────────
filter <- dplyr::filter
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


# see recovery plot: 

# Access the results directly
results.sim$recovery$table
results.sim$recovery$plot
results.sim$recovery$manuscript_table

# View the plot
p.sim.recovery <- results.sim$recovery$plot

ggsave(paste0("Output/Plots/p.sim.recovery.jpeg"),
       p.sim.recovery, width=30, height=15, units="cm")



# Check coverage — already printed in subtitle:
# "Overall coverage: 30/30 (100%); mean |bias|: 23.8%"


results.sim$recovery$table |> 
  dplyr::arrange(desc(abs(rel_bias_pct))) |>
  dplyr::select(variable, true_value, mean, rel_bias_pct) |>
  print(n = 30)

#what drives bias?
results.sim$recovery$table |>
  dplyr::arrange(desc(abs(rel_bias_pct))) |>
  dplyr::select(variable, true_value, mean, rel_bias_pct, identifiability) |>
  print(n=30)







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

