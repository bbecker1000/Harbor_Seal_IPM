# ============================================================================
# HARBOR SEAL IPM v3.2 — PLOTS, TABLES & POST-PROCESSING
# ============================================================================
#
# COMPANION SCRIPT to harbor_seal_ipm_v3.2.R
#
# Can be used in two ways:
#
#   (A) Automatically sourced by run_full_analysis_v3.2() — no action needed.
#
#   (B) Independently, to replot or refine outputs from a saved fit:
#
#         source("harbor_seal_ipm_v3.2.R")       # load model + data functions
#         source("harbor_seal_ipm_v3.2_plots.R")  # load plot functions
#
#         out <- load_seal_results("IPM_v3.2_real")
#
#         # Re-run everything:
#         run_all_plots_v3.2(out$fit, out$sim_data, prefix="IPM_v3.2_real")
#
#         # Or refine a single plot:
#         ts <- create_timeseries_plots_v3.2(out$fit, out$sim_data,
#                                            prefix="IPM_v3.2_real")
#         ts$total   # view in RStudio
#         ts$phi_pup
#
#         # Portfolio / synchrony (not run by default):
#         port <- create_portfolio_analysis_v3.2(out$fit, out$sim_data)
#         sync <- create_synchrony_projections_v3.2(out$fit, out$sim_data)
#
# Functions in this file:
#   load_seal_results()               — load saved fit + reconstruct sim_data
#   run_all_plots_v3.2()              — orchestrate all plots and tables
#   check_diagnostics_v3.2()          — Part 4
#   create_trace_plots_v3.2()         — Part 5
#   check_parameter_recovery_v3.2()   — Part 6
#   create_ppc_plots_v3.2()           — Part 7
#   create_timeseries_plots_v3.2()    — Part 8
#   create_site_age_timeseries_v3.2() — Part 9
#   create_projection_plots_v3.2()    — Part 10
#   create_effect_plots_v3.2()        — Part 11
#   create_summary_table_v3.2()       — Part 12
#   save_model_output_v3.2()          — Part 13
#   create_portfolio_analysis_v3.2()  — Part 15
#   create_synchrony_projections_v3.2() — Part 16
#
# ============================================================================

library(tidyverse)
library(posterior)
library(bayesplot)
library(patchwork)

dir.create("Output",       showWarnings = FALSE)
dir.create("Output/Plots", showWarnings = FALSE)

# ============================================================================
# LOAD HELPER — reconstruct sim_data from saved input_data and fit
# ============================================================================

load_seal_results <- function(prefix         = "IPM_v3.2_real",
                              fit_path       = NULL,
                              input_data_path = NULL,
                              years          = 1997:2025,
                              T_proj         = 10) {
  # Locate fit object
  if (is.null(fit_path))
    fit_path <- paste0("Output/harbor_seal_", prefix, "_fit.rds")
  if (!file.exists(fit_path))
    stop("Fit file not found: ", fit_path,
         "\nRun run_full_analysis_v3.2() first, or pass fit_path explicitly.")
  cat("Loading fit from", fit_path, "...\n")
  fit <- readRDS(fit_path)
  
  # Reconstruct minimal sim_data for plotting
  # If the original stan input data was saved, reload it for full fidelity.
  # Otherwise construct a minimal shell with site/year metadata.
  # Auto-detect input_data_path if not supplied
  if (is.null(input_data_path)) {
    auto_path <- paste0("Output/harbor_seal_", prefix, "_input_data.rds")
    if (file.exists(auto_path)) input_data_path <- auto_path
  }
  if (!is.null(input_data_path) && file.exists(input_data_path)) {
    inp      <- readRDS(input_data_path)
    sim_data <- list(
      stan_data      = inp$stan_data,
      site_names     = c("BL","DE","DP","PRH","TB","TP"),
      years          = inp$years %||% years,
      scenario_names = c("Status Quo","Warm (MOCI +1)",
                         "Cool (MOCI -1)","Warm + High Coyote"),
      true_params    = inp$true_params %||% NULL
    )
  } else {
    cat("No input_data_path supplied — using metadata shell.\n")
    cat("  Pass input_data_path for PPC plots and parameter recovery.\n")
    T_val <- fit$metadata()$data$T %||% length(years)
    sim_data <- list(
      stan_data      = list(T=T_val, S=6,
                            T_proj=T_proj, N_scenarios=4,
                            y_adult_obs=matrix(1L,6,T_val),
                            y_pup_obs  =matrix(1L,6,T_val),
                            y_molt_obs =matrix(1L,6,T_val),
                            y_adult=matrix(0,6,T_val),
                            y_pup  =matrix(0,6,T_val),
                            y_molt =matrix(0,6,T_val)),
      site_names     = c("BL","DE","DP","PRH","TB","TP"),
      years          = years[seq_len(T_val)],
      scenario_names = c("Status Quo","Warm (MOCI +1)",
                         "Cool (MOCI -1)","Warm + High Coyote"),
      true_params    = NULL
    )
  }
  # Validate that the fit object has working methods
  fit_ok <- tryCatch(is.function(fit$draws), error = function(e) FALSE)
  if (!fit_ok) {
    warning(paste0(
      "The loaded fit object may not have working methods.\n",
      "This can happen when a CmdStanFit is saved inside a list.\n",
      "Try: fit <- readRDS(\"Output/harbor_seal_", prefix, "_fit.rds\")"))
  }
  cat("Ready. Call run_all_plots_v3.2(fit, sim_data) to regenerate all outputs.\n")
  list(fit=fit, sim_data=sim_data)
}

`%||%` <- function(x, y) if (!is.null(x)) x else y

# ============================================================================
# ORCHESTRATOR — run all plots and tables in sequence
# ============================================================================

run_all_plots_v3.2 <- function(fit,
                               sim_data,
                               prefix        = "IPM_v3.2",
                               save          = TRUE,
                               run_recovery  = FALSE,
                               run_portfolio = FALSE,
                               run_synchrony = FALSE) {
  
  # Each step is wrapped in tryCatch so a single plotting error does not
  # abort the whole pipeline. Failures print a clear recovery message.
  # The fit and sim_data are saved to RDS before this function is called
  # (in run_full_analysis_v3.2), so you can always reload and rerun:
  #   out <- load_seal_results(prefix)
  #   run_all_plots_v3.2(out$fit, out$sim_data, prefix=prefix)
  
  safe_run <- function(label, expr) {
    cat(sprintf("\n%s\n", label))
    # withCallingHandlers muffles warnings IN PLACE (no restart).
    # tryCatch catches errors after the fact and returns NULL.
    # Using tryCatch(warning=) here would restart the expression on every
    # warning, causing the entire step to execute twice.
    tryCatch(
      withCallingHandlers(
        expr,
        warning = function(w) invokeRestart("muffleWarning")
      ),
      error = function(e) {
        cat(sprintf("  !! FAILED: %s\n", conditionMessage(e)))
        NULL
      }
    )
  }

  # ── v3.2 parameter validation ──────────────────────────────────────────────
  # Checks that all parameters introduced in v3.2 are present in this fit.
  # A mismatch indicates the fit was produced by an older model version;
  # some plots will silently drop parameters or use incorrect baselines.
  v3.2_required <- c(
    "phi_pup_logit", "phi_adult_F_logit", "delta_adult", "p_male_breed",
    "beta_moci_ond_fecund", "beta_moci_ond_pup", "beta_moci_jfm_pup",
    "beta_moci_amj_pup", "beta_moci_jfm_juv", "beta_moci_jfm_adult",
    "detect_breed_logit", "detect_molt_logit"
  )
  avail_vars   <- tryCatch(fit$summary()$variable, error = function(e) character(0))
  missing_v3.2 <- v3.2_required[!v3.2_required %in% avail_vars]
  if (length(missing_v3.2) > 0) {
    cat(sprintf(paste0(
      "\n!! VERSION WARNING: fit is missing %d v3.2 parameter(s):\n",
      "     %s\n",
      "   This fit may have been produced by an older model version.\n",
      "   Some plots will silently omit these parameters.\n\n"),
      length(missing_v3.2), paste(missing_v3.2, collapse = ", ")))
  } else {
    cat("v3.2 parameter check passed.\n")
  }

  diag <- safe_run("── Diagnostics ───────────────────────────────────────────",
                   check_diagnostics_v3.2(fit))
  
  params_candidate <- c(
    "phi_pup_logit","phi_juv_base","phi_adult_F_logit","phi_adult_F_base","delta_adult",
    "fecund_primip","fecund_mature","prop_female","p_male_breed",
    "beta_coy[1]","beta_coy[2]","beta_coy[3]",
    paste0("beta_dist_surv[",1:6,"]"),
    paste0("beta_dist_detect[",1:6,"]"),
    "beta_moci_ond_fecund","beta_moci_amj_pup","beta_moci_jfm_juv",
    "beta_moci_jfm_adult","beta_eseal_pup",
    "detect_breed_logit","detect_molt_logit",
    "beta_moci_ond_pup","beta_moci_jfm_pup",
    "sigma_process","sigma_obs_adult","sigma_obs_pup","sigma_obs_molt","sigma_site"
  )
  # Filter to only params present in this model fit — prevents trace/summary
  # failures when plotting against a fit that pre-dates new parameter additions.
  available_vars <- tryCatch(
    fit$summary()$variable,
    error = function(e) params_candidate
  )
  params <- if (!is.null(diag) && !is.null(diag$params)) {
    diag$params[diag$params %in% available_vars]
  } else {
    params_candidate[params_candidate %in% available_vars]
  }
  
  traces <- safe_run("── Trace plots ───────────────────────────────────────────",
                     create_trace_plots_v3.2(fit, params, save=save, prefix=prefix))
  
  rec <- NULL
  if (run_recovery && !is.null(sim_data$true_params))
    rec <- safe_run("── Parameter recovery ────────────────────────────────────",
                    check_parameter_recovery_v3.2(fit, sim_data, save=save, prefix=prefix))
  else
    cat("\n── Parameter recovery: skipped (simulation only; set run_recovery=TRUE) ─\n")
  
  ppc  <- safe_run("── Posterior predictive checks ───────────────────────────",
                   create_ppc_plots_v3.2(fit, sim_data, save=save, prefix=prefix))
  
  ts   <- safe_run("── Population time series ────────────────────────────────",
                   create_timeseries_plots_v3.2(fit, sim_data, save=save, prefix=prefix))
  
  sa   <- safe_run("── Site × age class time series ──────────────────────────",
                   create_site_age_timeseries_v3.2(fit, sim_data, save=save, prefix=prefix))
  
  sasite <- safe_run("── Population by age class at each site ──────────────────",
                     create_site_panels_v3.2(fit, sim_data, save=save, prefix=prefix))
  
  proj <- safe_run("── Projection plots ──────────────────────────────────────",
                   create_projection_plots_v3.2(fit, sim_data, save=save, prefix=prefix))
  
  eff  <- safe_run("── Covariate effect plots ────────────────────────────────",
                   create_effect_plots_v3.2(fit, save=save, prefix=prefix))
  
  jveff <- safe_run("── Juvenile + adult survival effects ─────────────────────",
                    create_juv_adult_effect_plots_v3.2(fit, save=save, prefix=prefix))
  
  forest <- safe_run("── Coefficient forest plot ───────────────────────────────",
                     create_forest_plot_v3.2(fit, save=save, prefix=prefix))
  
  decomp <- safe_run("── Covariate decomposition stacked bars ──────────────────",
                     create_covariate_decomposition_plots_v3.2(fit, sim_data, save=save, prefix=prefix))
  
  tbl  <- safe_run("── Parameter summary table ───────────────────────────────",
                   create_summary_table_v3.2(fit, save=save, prefix=prefix))
  
  safe_run("── Key results printout ──────────────────────────────────",
           save_model_output_v3.2(fit, prefix=prefix))
  
  port <- sync <- NULL
  
  if (run_portfolio)
    port <- safe_run("── Portfolio analysis ────────────────────────────────────",
                     create_portfolio_analysis_v3.2(fit, sim_data, save=save, prefix=prefix))
  else
    cat("\n── Portfolio analysis: skipped (set run_portfolio=TRUE) ─────\n")
  
  if (run_synchrony)
    sync <- safe_run("── Synchrony projections ─────────────────────────────────",
                     create_synchrony_projections_v3.2(fit, sim_data, save=save, prefix=prefix))
  else
    cat("\n── Synchrony projections: skipped (set run_synchrony=TRUE) ──\n")
  
  # Summarise outcome
  results <- list(diagnostics=diag, traces=traces, recovery=rec,
                  ppc=ppc, ts=ts, site_age=sa, site_panels=sasite,
                  projections=proj, effects=eff, juv_adult_effects=jveff,
                  forest=forest, decomposition=decomp,
                  table=tbl, portfolio=port, sync=sync)
  # ── Tally results ─────────────────────────────────────────────────────────
  # Distinguish genuinely skipped steps from genuine failures.
  # Skipped steps are expected NULLs; failures are unexpected NULLs.
  skipped_keys <- c(
    if (!run_recovery || is.null(sim_data$true_params)) "recovery" else character(0),
    if (!run_portfolio)                                  "portfolio" else character(0),
    if (!run_synchrony)                                  "sync"      else character(0)
  )
  attempted <- results[!names(results) %in% skipped_keys]
  n_ok      <- sum(!sapply(attempted, is.null))
  n_fail    <- sum( sapply(attempted, is.null))
  failed_names <- names(attempted)[sapply(attempted, is.null)]
  
  cat(sprintf("\n── Pipeline complete: %d/%d steps succeeded ──────────────\n",
              n_ok, length(attempted)))
  cat(sprintf("   Plots  -> Output/Plots/%s_*\n", prefix))
  cat(sprintf("   Fit    -> Output/harbor_seal_%s_fit.rds\n", prefix))
  if (length(skipped_keys) > 0)
    cat(sprintf("   Skipped (%d): %s\n",
                length(skipped_keys), paste(skipped_keys, collapse=", ")))
  if (n_fail > 0) {
    cat(sprintf("   Failed  (%d): %s\n", n_fail, paste(failed_names, collapse=", ")))
    cat(sprintf(paste0(
      "   To retry failed steps, re-source and run:\n",
      "     out <- load_seal_results('%s')\n",
      "     run_all_plots_v3.2(out$fit, out$sim_data, prefix='%s')\n"
    ), prefix, prefix))
  }
  
  invisible(results)
}

# ============================================================================
# PLOT THEME + CREDIBLE INTERVAL CONSTANTS
# ============================================================================
# Unified theme applied to all ggplot outputs for consistent publication style.
# 89% credible intervals used throughout (Kruschke 2014; McElreath 2020):
#   - More honest about uncertainty than 95% CrI for Bayesian posteriors
#   - Corresponds to quantiles 0.055 (lower) and 0.945 (upper)

CI_LO    <- 0.055   # lower quantile: 89% CrI
CI_HI    <- 0.945   # upper quantile: 89% CrI
CI_LABEL <- "89% CrI"   # for paste0 / direct string use
CI_FMT   <- "89%% CrI"  # for sprintf (% must be escaped as %%)

# Colour palette (consistent across all plots)
SEAL_COLS <- list(
  pop      = "#2166AC",   # total population / adult lines
  pup      = "#1B7837",   # pup survival / pup counts
  juv      = "#762A83",   # juvenile
  adult_f  = "#B2182B",   # adult female
  adult_m  = "#D6604D",   # adult male
  molt     = "#8C510A",   # molt counts
  ribbon   = "#AECDE8",   # default ribbon fill
  neutral  = "#888888",   # neutral / historical
  warn     = "#FF7F00"    # warning / observed sex ratio
)

# Shared theme function — call instead of theme_minimal() throughout
theme_seal <- function(base_size = 16) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      # Panel
      panel.grid.major   = element_line(colour = "grey88", linewidth = 0.4),
      panel.grid.minor   = element_line(colour = "grey93", linewidth = 0.2),
      panel.border       = element_rect(colour = "grey70", fill = NA, linewidth = 0.5),
      # Axes
      axis.title         = element_text(size = rel(0.95), colour = "grey20"),
      axis.text          = element_text(size = rel(0.88), colour = "grey30"),
      axis.ticks         = element_line(colour = "grey70", linewidth = 0.3),
      # Legend
      legend.position    = "bottom",
      legend.title       = element_text(size = rel(0.88), face = "bold"),
      legend.text        = element_text(size = rel(0.82)),
      legend.key.width   = unit(1.6, "cm"),
      legend.background  = element_blank(),
      # Facet labels
      strip.text         = element_text(size = rel(0.90), face = "bold", colour = "grey20"),
      strip.background   = element_rect(fill = "grey94", colour = "grey80", linewidth = 0.3),
      # Plot titles
      plot.title         = element_text(size = rel(1.05), face = "bold", margin = margin(b=6)),
      plot.subtitle      = element_text(size = rel(0.88), colour = "grey40", margin = margin(b=8)),
      plot.caption       = element_text(size = rel(0.78), colour = "grey50", hjust = 1),
      plot.margin        = margin(10, 14, 10, 10)
    )
}

# Convenience: compute 89% CrI summary for a matrix of draws (cols = time/params)
ci89_matrix <- function(m) {
  list(
    mean = colMeans(m),
    lo   = apply(m, 2, quantile, CI_LO),
    hi   = apply(m, 2, quantile, CI_HI)
  )
}

# ============================================================================
# ============================================================================
# PART 4: DIAGNOSTICS
# ============================================================================
# 89% CrI summary helper for cmdstanr fit objects
seal_fit_summary <- function(fit, variables) {
  # Defensive implementation using only stable base-R operations.
  # Both fit$summary() and fit$draws() are wrapped in tryCatch so a
  # stale/RDS-reloaded fit gives a clear reload instruction rather than
  # a cryptic "attempt to apply non-function" error.
  
  ci_lo <- if (exists("CI_LO")) CI_LO else 0.055
  ci_hi <- if (exists("CI_HI")) CI_HI else 0.945
  
  stale_msg <- paste0(
    "fit object methods are unavailable (stale or RDS-reloaded fit).\n",
    "Reload with:  out <- load_seal_results(\"IPM_v3.2_real\")\n",
    "then pass:    out$fit  and  out$sim_data")
  
  s <- tryCatch(
    fit$summary(variables = variables),
    error = function(e) stop(stale_msg, "\nOriginal: ", conditionMessage(e))
  )
  d <- tryCatch(
    fit$draws(variables = variables, format = "matrix"),
    error = function(e) stop(stale_msg, "\nOriginal: ", conditionMessage(e))
  )
  
  s$q_lo <- apply(d[, s$variable, drop = FALSE], 2, quantile, ci_lo)
  s$q_hi <- apply(d[, s$variable, drop = FALSE], 2, quantile, ci_hi)
  
  # Newer posterior versions (>=1.6) may return rhat/ess columns as S3
  # class objects or list columns rather than plain doubles, which breaks
  # the > operator. Coerce element-by-element to handle both cases.
  .to_num <- function(x) {
    if (is.list(x))
      vapply(x, function(v) suppressWarnings(as.numeric(v)[1L]), numeric(1))
    else
      suppressWarnings(as.numeric(x))
  }
  for (col in c("rhat", "ess_bulk", "ess_tail"))
    if (col %in% names(s)) s[[col]] <- .to_num(s[[col]])
  
  s
}


check_diagnostics_v3.2 <- function(fit) {
  
  cat("\n============================================\n")
  cat("MODEL DIAGNOSTICS — IPM v3.2\n")
  cat("============================================\n\n")
  # fit$cmdstan_diagnose() reads original CSV files (temp dir — gone after restart).
  # Use fit$diagnostic_summary() which works on cached results, with cmdstan_diagnose
  # as a fallback if CSVs are still present.
  diag_ok <- tryCatch({
    ds <- fit$diagnostic_summary(quiet = TRUE)
    cat("Checking sampler transitions treedepth.\n")
    cat(if (all(ds$num_max_treedepth == 0))
      "Treedepth satisfactory for all transitions.\n"
      else sprintf("WARNING: %d transitions hit max treedepth.\n",
                   sum(ds$num_max_treedepth)))
    cat("\nChecking sampler transitions for divergences.\n")
    cat(if (all(ds$num_divergent == 0))
      "No divergent transitions found.\n"
      else sprintf("WARNING: %d divergent transitions.\n",
                   sum(ds$num_divergent)))
    cat("\nChecking E-BFMI.\n")
    cat(if (all(ds$ebfmi > 0.2))
      "E-BFMI satisfactory.\n"
      else sprintf("WARNING: low E-BFMI in %d chain(s).\n",
                   sum(ds$ebfmi <= 0.2)))
    cat("\nProcessing complete, no problems detected.\n")
    invisible(ds)
  }, error = function(e) {
    # CSV files unavailable — print a note and skip
    cat(sprintf("NOTE: Full sampler diagnostics unavailable (CSV files not found).\n"))
    cat(sprintf("      Reload fit with load_seal_results() for full diagnostics.\n"))
    cat(sprintf("      Original error: %s\n", conditionMessage(e)))
    NULL
  })
  
  params <- c(
    "phi_pup_logit","phi_juv_base","phi_adult_F_logit","phi_adult_F_base","delta_adult",
    "fecund_primip","fecund_mature",
    "prop_female","p_male_breed",
    "beta_coy[1]","beta_coy[2]","beta_coy[3]",
    "beta_dist_surv[1]","beta_dist_surv[2]","beta_dist_surv[3]",
    "beta_dist_surv[4]","beta_dist_surv[5]","beta_dist_surv[6]",
    "beta_moci_ond_fecund","beta_moci_ond_pup","beta_moci_amj_pup",
    "beta_moci_jfm_pup","beta_moci_jfm_juv",
    "beta_moci_jfm_adult","beta_eseal_pup",
    "detect_breed_logit","detect_molt_logit",
    "sigma_process","sigma_obs_adult","sigma_obs_pup","sigma_obs_molt","sigma_site"
  )
  
  s <- seal_fit_summary(fit, params)
  cat("\nParameter Summary:\n")
  print(s |> select(variable,mean,sd,q_lo,q_hi,rhat,ess_bulk))
  
  # Convenience: report phi_pup on probability scale
  pup_logit_draws <- fit$draws(variables="phi_pup_logit", format="df")$phi_pup_logit
  cat(sprintf("\nphi_pup_base (prob scale): median=%.3f, 90%% CrI=[%.3f, %.3f]\n",
              median(plogis(pup_logit_draws)),
              quantile(plogis(pup_logit_draws), CI_LO),
              quantile(plogis(pup_logit_draws), CI_HI)))
  
  # Force rhat and ess to plain double before comparison — immune to
  # S3-class columns, list columns, or dplyr dispatch issues.
  # unclass() strips S3, as.vector() strips all attributes, as.double() converts.
  rhat_v <- tryCatch(as.double(as.vector(unclass(s$rhat))),
                     error = function(e) rep(NA_real_, nrow(s)))
  ess_v  <- tryCatch(as.double(as.vector(unclass(s$ess_bulk))),
                     error = function(e) rep(NA_real_, nrow(s)))
  bad <- s[!is.na(rhat_v) & rhat_v > 1.05, , drop=FALSE]
  low <- s[!is.na(ess_v)  & ess_v  < 400,  , drop=FALSE]
  cat(paste0("\n", CI_LABEL, " used for all intervals\n"))
  if (nrow(bad)>0) { cat("\nWARNING: Rhat > 1.05:\n"); print(bad[, c('variable','rhat')]) }
  if (nrow(low)>0) { cat("\nWARNING: low ESS:\n");     print(low[, c('variable','ess_bulk')]) }
  
  list(params=params, summary=s)
}


# ============================================================================
# PART 5: TRACE PLOTS
# ============================================================================

create_trace_plots_v3.2 <- function(fit, params, save=TRUE, prefix="IPM_v3.2") {
  
  draws <- fit$draws(format="df")
  
  p1 <- mcmc_trace(draws,
                   pars=c("phi_pup_logit","phi_juv_base","phi_adult_F_logit",
                          "delta_adult","p_male_breed",
                          "detect_breed_logit","detect_molt_logit")) +
    labs(title="Trace: Survival + Observation Parameters")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_trace_survival.jpeg"),
                   p1, width=30, height=18, units="cm")
  
  p2 <- mcmc_trace(draws,
                   pars=c("fecund_primip","fecund_mature","prop_female")) +
    labs(title="Trace: Fecundity Parameters")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_trace_fecundity.jpeg"),
                   p2, width=30, height=18, units="cm")
  
  p3 <- mcmc_trace(draws, pars=c("beta_coy[1]","beta_coy[2]","beta_coy[3]")) +
    labs(title="Trace: Coyote Effects (BL, DE, DP)")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_trace_coyote.jpeg"),
                   p3, width=30, height=12, units="cm")
  
  p4 <- mcmc_trace(draws, pars=c(paste0("beta_dist_surv[",1:6,"]"),
                                  paste0("beta_dist_detect[",1:6,"]"))) +
    labs(title="Trace: Site-Specific Disturbance Effects")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_trace_disturbance.jpeg"),
                   p4, width=30, height=18, units="cm")
  
  p5 <- mcmc_trace(draws,
                   pars=c("beta_moci_ond_fecund","beta_moci_ond_pup",
                          "beta_moci_amj_pup","beta_moci_jfm_pup",
                          "beta_moci_jfm_juv","beta_moci_jfm_adult",
                          "beta_moci_amj_molt")) +
    labs(title="Trace: MOCI Effects")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_trace_moci.jpeg"),
                   p5, width=30, height=18, units="cm")
  
  p6 <- mcmc_trace(draws,
                   pars=c("sigma_process","sigma_obs_adult","sigma_obs_pup","sigma_obs_molt","sigma_site")) +
    labs(title="Trace: Error Terms")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_trace_errors.jpeg"),
                   p6, width=30, height=18, units="cm")
  
  list(survival=p1, fecundity=p2, coyote=p3, disturbance=p4, moci=p5, errors=p6)
}


# ============================================================================
# PART 6: PARAMETER RECOVERY
# ============================================================================

check_parameter_recovery_v3.2 <- function(fit, sim_data, save=TRUE, prefix="IPM_v3.2") {
  
  tp <- sim_data$true_params
  
  true_vals <- tibble(
    parameter  = c(
      "phi_pup_logit","phi_juv_base","phi_adult_F_logit","delta_adult",
      "fecund_primip","fecund_mature",
      "prop_female","p_male_breed",
      "beta_coy[1]","beta_coy[2]","beta_coy[3]",
      "beta_dist_surv[1]","beta_dist_surv[2]","beta_dist_surv[3]",
      "beta_dist_surv[4]","beta_dist_surv[5]","beta_dist_surv[6]",
      "beta_moci_ond_fecund","beta_moci_amj_pup","beta_moci_jfm_juv",
      "beta_moci_jfm_adult","beta_eseal_pup",
      "sigma_process","sigma_obs_adult","sigma_obs_pup","sigma_obs_molt"
    ),
    true_value = c(
      tp$phi_pup_logit, tp$phi_juv_base, tp$phi_adult_F_logit, tp$delta_adult,
      tp$fecund_primip, tp$fecund_mature,
      tp$prop_female, tp$p_male_breed,
      tp$beta_coy[1], tp$beta_coy[2], tp$beta_coy[3],
      tp$beta_dist_surv[1], tp$beta_dist_surv[2], tp$beta_dist_surv[3],
      tp$beta_dist_surv[4], tp$beta_dist_surv[5], tp$beta_dist_surv[6],
      tp$beta_moci_ond_fecund, tp$beta_moci_amj_pup, tp$beta_moci_jfm_juv,
      tp$beta_moci_jfm_adult, tp$beta_eseal_pup,
      tp$sigma_process, tp$sigma_obs_adult, tp$sigma_obs_pup, tp$sigma_obs_molt
    )
  )
  
  rec <- seal_fit_summary(fit, true_vals$parameter) |>
    left_join(true_vals, by=c("variable"="parameter")) |>
    mutate(recovered    = true_value>=q_lo & true_value<=q_hi,
           rel_bias_pct = (mean-true_value)/abs(true_value)*100)
  
  cat("Recovery rate:", sum(rec$recovered),"/",nrow(rec),
      sprintf("(%.1f%%)\n", 100*mean(rec$recovered)))
  print(rec |> select(variable,true_value,mean,q_lo,q_hi,recovered,rel_bias_pct) |>
          mutate(across(where(is.numeric),~round(.x,3))))
  
  p <- ggplot(rec, aes(x=true_value,y=mean)) +
    geom_abline(slope=1,intercept=0,linetype=2,color="gray50") +
    geom_pointrange(aes(ymin=q_lo,ymax=q_hi,color=recovered),size=0.8) +
    geom_text(aes(label=variable),hjust=-0.1,vjust=-0.3,size=2.5,check_overlap=TRUE) +
    scale_color_manual(values=c("TRUE"=SEAL_COLS$pop,"FALSE"=SEAL_COLS$adult_f)) +
    labs(x="True Value",y=paste0("Estimated (", CI_LABEL, ")"),title="Parameter Recovery: IPM v3.2") +
    theme_seal() + coord_equal()
  if (save) ggsave(paste0("Output/Plots/",prefix,"_parameter_recovery.jpeg"),
                   p, width=25, height=22, units="cm")
  
  list(table=rec, plot=p)
}


# ============================================================================
# PART 7: POSTERIOR PREDICTIVE CHECKS
# ============================================================================

create_ppc_plots_v3.2 <- function(fit, sim_data, save=TRUE, prefix="IPM_v3.2") {
  
  rep_a  <- fit$draws(variables="y_adult_rep", format="matrix")
  rep_p  <- fit$draws(variables="y_pup_rep",   format="matrix")
  rep_m  <- fit$draws(variables="y_molt_rep",  format="matrix")
  
  obs_a  <- as.vector(t(sim_data$stan_data$y_adult))
  obs_p  <- as.vector(t(sim_data$stan_data$y_pup))
  obs_m  <- as.vector(t(sim_data$stan_data$y_molt))
  
  ind_a  <- as.vector(t(sim_data$stan_data$y_adult_obs))==1
  ind_p  <- as.vector(t(sim_data$stan_data$y_pup_obs))  ==1
  ind_m  <- as.vector(t(sim_data$stan_data$y_molt_obs)) ==1
  
  p_comb <- (ppc_dens_overlay(obs_a[ind_a],rep_a[1:100,ind_a]) + labs(title="PPC: Adult")) /
    (ppc_dens_overlay(obs_p[ind_p],rep_p[1:100,ind_p]) + labs(title="PPC: Pup"))   /
    (ppc_dens_overlay(obs_m[ind_m],rep_m[1:100,ind_m]) + labs(title="PPC: Molt"))
  if (save) ggsave(paste0("Output/Plots/",prefix,"_ppc_density.jpeg"),
                   p_comb, width=25, height=30, units="cm")
  
  sg <- rep(sim_data$site_names, each=sim_data$stan_data$T)
  p_site <- ppc_stat_grouped(obs_a[ind_a], rep_a[1:100,ind_a],
                             group=sg[ind_a], stat="mean") +
    labs(title="PPC: Mean Adult by Site")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_ppc_by_site.jpeg"),
                   p_site, width=25, height=15, units="cm")
  
  list(density=p_comb, by_site=p_site)
}


# ============================================================================
# PART 8: TIME SERIES PLOTS
# ============================================================================

create_timeseries_plots_v3.2 <- function(fit, sim_data, save=TRUE, prefix="IPM_v3.2") {
  
  years <- sim_data$years; site_names <- sim_data$site_names
  S <- length(site_names); T <- length(years)
  
  # Total population
  Ntot <- fit$draws(variables="N_total_all", format="df") |> select(starts_with("N_total_all"))
  p_total <- ggplot(tibble(Year=years, mean=colMeans(Ntot),
                           lo=as.numeric(apply(Ntot,2,quantile,CI_LO)),
                           hi=as.numeric(apply(Ntot,2,quantile,CI_HI))),
                    aes(x=Year)) +
    geom_ribbon(aes(ymin=lo,ymax=hi),alpha=0.25,fill=SEAL_COLS$ribbon) +
    geom_line(aes(y=mean),linewidth=1.2,color=SEAL_COLS$pop) +
    labs(x="Year",y="Total Population",
         title="Estimated Total Harbor Seal Population") +
    theme_seal()
  if (save) ggsave(paste0("Output/Plots/",prefix,"_total_population.jpeg"),
                   p_total, width=25, height=15, units="cm")
  
  # Pup survival on probability scale
  pup_logit <- fit$draws(variables="phi_pup_logit",format="df")$phi_pup_logit
  pup_prob  <- plogis(pup_logit)
  p_phipup  <- ggplot(tibble(x=pup_prob), aes(x=x)) +
    geom_density(fill=SEAL_COLS$pup,alpha=0.5) +
    geom_vline(xintercept=median(pup_prob),linetype="dashed",color="darkgreen") +
    geom_vline(xintercept=0.5, linetype="dotted", color="gray50") +
    annotate("text",x=0.5,y=Inf,label=" Field data ~0.50",
             hjust=0,vjust=1.5,size=3.5,color="gray40") +
    labs(x="Pup survival probability", y="Density",
         title="Posterior: Pup Survival (sex-neutral)",
         subtitle=sprintf(paste0("Median = %.3f  (", CI_FMT, ": %.3f–%.3f)"),
                          median(pup_prob),
                          quantile(pup_prob,CI_LO),
                          quantile(pup_prob,CI_HI))) +
    theme_seal()
  if (save) ggsave(paste0("Output/Plots/",prefix,"_phi_pup_posterior.jpeg"),
                   p_phipup, width=20, height=12, units="cm")
  
  # Male haul-out fraction
  pmb <- fit$draws(variables="p_male_breed",format="df")$p_male_breed
  p_pmb <- ggplot(tibble(x=pmb),aes(x=x)) +
    geom_density(fill=SEAL_COLS$pop,alpha=0.5) +
    geom_vline(xintercept=median(pmb),linetype="dashed",color="#08306b") +
    labs(x="p_male_breed",y="Density",
         title="Posterior: Male Haul-out Fraction (Breeding Season)",
         subtitle=sprintf(paste0("Median=%.3f  (", CI_FMT, ": %.3f–%.3f)"),
                          median(pmb),quantile(pmb,CI_LO),quantile(pmb,CI_HI))) +
    theme_seal()
  if (save) ggsave(paste0("Output/Plots/",prefix,"_p_male_breed.jpeg"),
                   p_pmb, width=20, height=12, units="cm")
  
  # Sex ratio: true vs observed
  sr_t <- fit$draws(variables="sex_ratio_adult",   format="matrix")
  sr_o <- fit$draws(variables="sex_ratio_observed",format="matrix")
  
  make_sr <- function(mat) {
    m <- sapply(1:T, function(t) rowMeans(mat[,paste0("sex_ratio_adult[",1:S,",",t,"]"),drop=FALSE]))
    tibble(Year=years, mean=colMeans(m),
           lo=as.numeric(apply(m,2,quantile,CI_LO)), hi=as.numeric(apply(m,2,quantile,CI_HI)))
  }
  make_sr_obs <- function(mat) {
    m <- sapply(1:T, function(t) rowMeans(mat[,paste0("sex_ratio_observed[",1:S,",",t,"]"),drop=FALSE]))
    tibble(Year=years, mean=colMeans(m),
           lo=as.numeric(apply(m,2,quantile,CI_LO)), hi=as.numeric(apply(m,2,quantile,CI_HI)))
  }
  
  sr_true_df <- tryCatch(make_sr(sr_t),    error=function(e) NULL)
  sr_obs_df  <- tryCatch(make_sr_obs(sr_o),error=function(e) NULL)
  
  p_sr <- ggplot() +
    geom_hline(yintercept=0.5,linetype="dotted",color="gray50") +
    { if (!is.null(sr_obs_df))  geom_ribbon(data=sr_obs_df,  aes(x=Year,ymin=lo,ymax=hi),alpha=0.15,fill="orange") } +
    { if (!is.null(sr_obs_df))  geom_line(data=sr_obs_df,    aes(x=Year,y=mean,color="Observed (spring)"),linewidth=1.1,linetype="dashed") } +
    { if (!is.null(sr_true_df)) geom_ribbon(data=sr_true_df, aes(x=Year,ymin=lo,ymax=hi),alpha=0.25,fill="purple") } +
    { if (!is.null(sr_true_df)) geom_line(data=sr_true_df,   aes(x=Year,y=mean,color="True (population)"),linewidth=1.2) } +
    scale_color_manual(values=c("True (population)"=SEAL_COLS$juv,"Observed (spring)"=SEAL_COLS$warn)) +
    labs(x="Year",y="Proportion female",color=NULL,
         title="Adult Sex Ratio: True vs Spring Survey Observation",
         subtitle="Observed > true because most males remain in water during breeding") +
    ylim(0.45,0.80) + theme_seal() + theme(legend.position="bottom")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_sex_ratio.jpeg"),
                   p_sr, width=25, height=15, units="cm")
  
  # Mean pup survival over time (sex-neutral)
  mp <- fit$draws(variables="mean_phi_pup", format="matrix")
  mp_df <- tibble(Year=years, mean=colMeans(mp),
                  lo=as.numeric(apply(mp,2,quantile,CI_LO)), hi=as.numeric(apply(mp,2,quantile,CI_HI)))
  p_phipup_t <- ggplot(mp_df,aes(x=Year)) +
    geom_ribbon(aes(ymin=lo,ymax=hi),alpha=0.25,fill=SEAL_COLS$pup) +
    geom_line(aes(y=mean),linewidth=1.2,color=SEAL_COLS$pup) +
    geom_hline(yintercept=0.5,linetype="dotted",color="gray50") +
    labs(x="Year",y="Pup survival (sex-neutral)",
         title="Mean Pup Survival Over Time (across all sites)") +
    theme_seal()
  if (save) ggsave(paste0("Output/Plots/",prefix,"_phi_pup_timeseries.jpeg"),
                   p_phipup_t, width=22, height=12, units="cm")
  
  list(total=p_total, phi_pup=p_phipup, p_male_breed=p_pmb,
       sex_ratio=p_sr, phi_pup_time=p_phipup_t)
}


# ============================================================================
# PART 9: SITE-BY-AGE TIME SERIES
# ============================================================================

create_site_age_timeseries_v3.2 <- function(fit, sim_data, save=TRUE, prefix="IPM_v3.2") {
  
  years <- sim_data$years; site_names <- sim_data$site_names
  S <- length(site_names); T <- length(years)
  
  pull_st <- function(var) {
    d <- fit$draws(variables=var, format="matrix")
    map_dfr(1:S, function(s) map_dfr(1:T, function(t) {
      cn <- paste0(var,"[",s,",",t,"]")
      if (!cn %in% colnames(d)) return(NULL)
      tibble(Site=site_names[s], Year=years[t],
             mean=mean(d[,cn]), lo=as.numeric(quantile(d[,cn],CI_LO)), hi=as.numeric(quantile(d[,cn],CI_HI)))
    }))
  }
  
  all_sum <- bind_rows(
    pull_st("N_pup")         |> mutate(Age_Class="Pup"),
    pull_st("N_juv_total")   |> mutate(Age_Class="Juvenile"),
    pull_st("N_adult_total") |> mutate(Age_Class="Adult")
  ) |> mutate(Age_Class=factor(Age_Class,levels=c("Pup","Juvenile","Adult")))
  
  p <- ggplot(all_sum,aes(x=Year,y=mean,color=Site,fill=Site)) +
    geom_ribbon(aes(ymin=lo,ymax=hi),alpha=0.15,color=NA) +
    geom_line(linewidth=1) + facet_wrap(~Age_Class,scales="free_y",ncol=1) +
    labs(x="Year",y="Population Size",title="Population by Age Class Across Sites") +
    theme_seal() + theme(legend.position="bottom")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_age_class_timeseries.jpeg"),
                   p, width=30, height=35, units="cm")
  
  list(by_age=p, data=all_sum)
}

# ============================================================================
# PART 9b: POPULATION BY AGE CLASS — ONE PANEL PER SITE
# ============================================================================
# Companion to create_site_age_timeseries_v3.2 but transposed:
#   rows = sites, columns = age classes → easy comparison within each site.

create_site_panels_v3.2 <- function(fit, sim_data, save=TRUE, prefix="IPM_v3.2") {
  
  years      <- sim_data$years
  site_names <- sim_data$site_names
  S <- length(site_names); T <- length(years)
  
  age_cols <- c(Pup=SEAL_COLS$pup, Juvenile=SEAL_COLS$juv, Adult=SEAL_COLS$pop)
  
  pull_st <- function(var, label) {
    d <- fit$draws(variables=var, format="matrix")
    map_dfr(1:S, function(s) map_dfr(1:T, function(t) {
      cn <- paste0(var,"[",s,",",t,"]")
      if (!cn %in% colnames(d)) return(NULL)
      tibble(Site      = site_names[s],
             Year      = years[t],
             Age_Class = label,
             mean      = mean(d[, cn]),
             lo        = as.numeric(quantile(d[, cn], CI_LO)),
             hi        = as.numeric(quantile(d[, cn], CI_HI)))
    }))
  }
  
  all_df <- bind_rows(
    pull_st("N_pup",         "Pup"),
    pull_st("N_juv_total",   "Juvenile"),
    pull_st("N_adult_total", "Adult")
  ) |>
    mutate(Age_Class = factor(Age_Class, levels=c("Pup","Juvenile","Adult")),
           Site      = factor(Site, levels=site_names))
  
  # One column per age class, one row per site
  p <- ggplot(all_df, aes(x=Year, y=mean, colour=Age_Class, fill=Age_Class)) +
    geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.15, colour=NA) +
    geom_line(linewidth=0.9) +
    facet_grid(Site ~ Age_Class, scales="free_y") +
    scale_colour_manual(values=age_cols, guide="none") +
    scale_fill_manual(  values=age_cols, guide="none") +
    labs(x="Year", y="Population Size",
         title="Harbor Seal Population by Site and Age Class",
         subtitle=paste0("Posterior mean ± ", CI_LABEL)) +
    theme_seal(base_size=14) +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=8),
          strip.text.y = element_text(size=10, face="bold"),
          panel.spacing = unit(0.4, "lines"))
  
  if (save) ggsave(paste0("Output/Plots/",prefix,"_site_age_panels.jpeg"),
                   p, width=36, height=6*S, units="cm", dpi=200)
  
  list(plot=p, data=all_df)
}


# ============================================================================
# PART 11b: COEFFICIENT FOREST PLOT
# ============================================================================
# All covariate coefficients in one plot, ordered and colour-coded by group:
#   MOCI effects | Coyote effects | Disturbance effects | Elephant seal
# Point = posterior mean; thick bar = 50% CrI; thin bar = 89% CrI.
# Reference line at 0; significant effects (89% CrI excludes 0) filled.

create_forest_plot_v3.2 <- function(fit, save=TRUE, prefix="IPM_v3.2") {
  
  site_names <- c("BL","DE","DP","PRH","TB","TP")
  
  # Parameter list with display labels and group assignments
  params <- tribble(
    ~variable,               ~label,                              ~group,          ~stage,
    "beta_moci_ond_fecund",     "MOCI Fall (OND) → fecundity",  "MOCI",          "Fecundity",
    "beta_moci_ond_pup",     "MOCI Fall OND (t) → pup surv", "MOCI",          "Pup",
    "beta_moci_jfm_pup",     "MOCI Winter JFM (t) → pup surv","MOCI",         "Pup",
    "beta_moci_amj_pup",     "MOCI Spring AMJ (t-1) → pup",  "MOCI",          "Pup",
    "beta_moci_jfm_juv",     "MOCI Winter → juvenile",       "MOCI",          "Juvenile",
    "beta_moci_jfm_adult",   "MOCI Winter → adult",          "MOCI",          "Adult",
    "beta_moci_amj_molt",    "MOCI Spring → molt detect",     "MOCI",          "Observation",
    "detect_breed_logit",    "Detect breed baseline (logit)",  "Detection",     "Observation",
    "detect_molt_logit",     "Detect molt baseline (logit)",   "Detection",     "Observation",
    "beta_coy[1]",           "Coyote → pup (BL)",            "Coyote",        "Pup",
    "beta_coy[2]",           "Coyote → pup (DE)",            "Coyote",        "Pup",
    "beta_coy[3]",           "Coyote → pup (DP)",            "Coyote",        "Pup",
    "beta_dist_surv[1]",     "Disturbance → pup (BL)",       "Disturbance",   "Pup",
    "beta_dist_surv[2]",     "Disturbance → pup (DE)",       "Disturbance",   "Pup",
    "beta_dist_surv[3]",     "Disturbance → pup (DP)",       "Disturbance",   "Pup",
    "beta_dist_surv[4]",     "Disturbance → pup (PRH)",      "Disturbance",   "Pup",
    "beta_dist_surv[5]",     "Disturbance → pup (TB)",       "Disturbance",   "Pup",
    "beta_dist_surv[6]",     "Disturbance → pup (TP)",        "Disturbance",   "Pup",
    "beta_dist_detect[1]",   "Disturbance → detection (BL)", "Disturbance",   "Observation",
    "beta_dist_detect[2]",   "Disturbance → detection (DE)", "Disturbance",   "Observation",
    "beta_dist_detect[3]",   "Disturbance → detection (DP)", "Disturbance",   "Observation",
    "beta_dist_detect[4]",   "Disturbance → detection (PRH)","Disturbance",   "Observation",
    "beta_dist_detect[5]",   "Disturbance → detection (TB)", "Disturbance",   "Observation",
    "beta_dist_detect[6]",   "Disturbance → detection (TP)", "Disturbance",   "Observation",
    "beta_eseal_pup",        "Elephant seal → pup",          "Elephant seal", "Pup"
  )
  
  # Group colours
  grp_cols <- c(MOCI          = "#2166AC",
                Detection      = "#4D9221",
                Coyote        = "#B2182B",
                Disturbance   = "#8C510A",
                "Elephant seal" = "#762A83")
  
  all_draws <- tryCatch(
    fit$draws(format = "draws_df"),
    error = function(e) tryCatch(
      fit$draws(format = "matrix"),
      error = function(e2) stop(
        "Cannot extract draws from fit object.\n",
        "If fit was loaded from RDS use: out <- load_seal_results(\"IPM_v3.2_real\")\n",
        "then pass out$fit rather than a saved list element.\n",
        "Original error: ", conditionMessage(e))
    )
  )
  
  # Extract posterior summaries — 50% and 89% CrI
  post <- map_dfr(seq_len(nrow(params)), function(i) {
    v <- tryCatch(as.numeric(posterior::extract_variable(all_draws, params$variable[i])),
                  error = function(e) numeric(0))
    if (length(v) == 0L) return(NULL)
    tibble(
      variable = params$variable[i],
      mean     = mean(v),
      lo89     = as.numeric(quantile(v, CI_LO)),
      hi89     = as.numeric(quantile(v, CI_HI)),
      lo50     = as.numeric(quantile(v, 0.25)),
      hi50     = as.numeric(quantile(v, 0.75))
    )
  })
  
  df <- post |>
    left_join(params, by="variable") |>
    mutate(
      label     = factor(label, levels=rev(params$label)),
      group     = factor(group, levels=names(grp_cols)),
      stage     = factor(stage, levels=c("Pup","Juvenile","Adult","Observation")),
      sig       = (lo89 > 0) | (hi89 < 0),    # 89% CrI excludes zero
      dot_shape = ifelse(sig, 18, 16)           # filled diamond if significant
    )
  
  p <- ggplot(df, aes(y=label, colour=group)) +
    # 89% CrI — thin line
    geom_linerange(aes(xmin=lo89, xmax=hi89), linewidth=0.7, alpha=0.7) +
    # 50% CrI — thick line
    geom_linerange(aes(xmin=lo50, xmax=hi50), linewidth=2.2, alpha=0.9) +
    # Posterior mean — point (diamond if significant)
    geom_point(aes(x=mean, shape=sig), size=3) +
    scale_shape_manual(values=c("FALSE"=16, "TRUE"=18),
                       labels=c("FALSE"="Not significant","TRUE"="Significant (89% CrI)"),
                       name=NULL) +
    # Zero reference
    geom_vline(xintercept=0, linetype="dashed", colour="grey40", linewidth=0.5) +
    # Shaded group bands (alternating)
    scale_colour_manual(values=grp_cols, name="Covariate group") +
    facet_grid(group ~ ., scales="free_y", space="free_y") +
    labs(x="Coefficient (logit scale)", y=NULL,
         title="Covariate Effects — Coefficient Forest Plot",
         subtitle=paste0("Thick bar = 50% CrI; thin bar = ", CI_LABEL,
                         "; diamond = 89% CrI excludes zero")) +
    theme_seal(base_size=14) +
    theme(legend.position    = "bottom",
          strip.text.y       = element_text(angle=0, face="bold", size=11),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          axis.text.y        = element_text(size=10))
  
  if (save) ggsave(paste0("Output/Plots/",prefix,"_forest_plot.jpeg"),
                   p, width=28, height=26, units="cm", dpi=200)
  
  list(plot=p, data=df)
}



# ============================================================================
# PART 10: PROJECTION PLOTS
# ============================================================================

create_projection_plots_v3.2 <- function(fit, sim_data, save=TRUE, prefix="IPM_v3.2") {
  
  years <- sim_data$years; T <- length(years)
  T_proj <- sim_data$stan_data$T_proj
  N_sc   <- sim_data$stan_data$N_scenarios
  sc_nm  <- sim_data$scenario_names
  pyrs   <- (max(years)+1):(max(years)+T_proj)
  
  all_d <- fit$draws(format="matrix")
  
  proj_df <- map_dfr(1:N_sc, function(sc) {
    cols <- grep(paste0("^N_total_all_proj\\[",sc,","), colnames(all_d))
    if (!length(cols)) return(NULL)
    pm <- all_d[,cols]
    tibble(Scenario=sc_nm[sc], Year=pyrs,
           mean=colMeans(pm), lo=as.numeric(apply(pm,2,quantile,CI_LO)), hi=as.numeric(apply(pm,2,quantile,CI_HI)))
  })
  
  Ntot <- fit$draws(variables="N_total_all",format="df") |> select(starts_with("N_total_all"))
  hist <- tibble(Scenario="Historical",Year=years,mean=colMeans(Ntot),
                 lo=as.numeric(apply(Ntot,2,quantile,CI_LO)),hi=as.numeric(apply(Ntot,2,quantile,CI_HI)))
  
  full <- bind_rows(hist,proj_df) |>
    mutate(Period=ifelse(Scenario=="Historical","Historical","Projection"))
  
  p <- ggplot() +
    geom_ribbon(data=filter(full,Period=="Historical"),
                aes(x=Year,ymin=lo,ymax=hi),alpha=0.3,fill="gray50") +
    geom_line(data=filter(full,Period=="Historical"),
              aes(x=Year,y=mean),linewidth=1.2,color="black") +
    geom_ribbon(data=filter(full,Period=="Projection"),
                aes(x=Year,ymin=lo,ymax=hi,fill=Scenario),alpha=0.2) +
    geom_line(data=filter(full,Period=="Projection"),
              aes(x=Year,y=mean,color=Scenario),linewidth=1.2) +
    geom_vline(xintercept=max(years),linetype=2,color="red") +
    scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2") +
    labs(x="Year",y="Total Population",
         title="10-Year Projections",
         subtitle=paste0("Bands = ", CI_LABEL, "; dashed line = projection start")) +
    theme_seal() + theme(legend.position="bottom")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_projections.jpeg"),
                   p, width=30, height=20, units="cm")
  
  list(projection=p, data=proj_df)
}


# ============================================================================
# PART 11: COVARIATE EFFECT PLOTS
# ============================================================================

create_effect_plots_v3.2 <- function(fit, save=TRUE, prefix="IPM_v3.2") {
  
  site_names <- c("BL","DE","DP","PRH","TB","TP")
  
  # Pull all draws once for all make_plot calls
  all_draws <- tryCatch(
    fit$draws(format = "draws_df"),
    error = function(e) tryCatch(
      fit$draws(format = "matrix"),
      error = function(e2) stop(
        "Cannot extract draws from fit object.\n",
        "If fit was loaded from RDS use: out <- load_seal_results(\"IPM_v3.2_real\")\n",
        "then pass out$fit rather than a saved list element.\n",
        "Original error: ", conditionMessage(e))
    )
  )
  
  # ── Core plot builder ───────────────────────────────────────────────────────
  # stage controls the survival baseline:
  #   "pup"     — logit baseline from phi_pup_logit
  #   "juv"     — logit baseline from qlogis(phi_juv_base)
  #   "adult_F" — logit baseline from phi_adult_F_logit
  #   "adult_M" — logit baseline from qlogis(phi_adult_F_base - delta_adult)
  # Uses .data[[]] pronoun so ggplot resolves columns in the data frame at
  # render time, not in the (already-gone) make_plot closure environment.
  make_plot <- function(param, xlab, title,
                        stage = "pup",
                        xr    = seq(-2, 2, length.out = 100),
                        ylims = NULL) {
    
    beta_v <- tryCatch(
      as.numeric(posterior::extract_variable(all_draws, param)),
      error = function(e) { warning(sprintf("No draws for '%s'", param)); numeric(0) }
    )
    if (length(beta_v) == 0L)
      return(list(plot=ggplot()+labs(title=paste("MISSING:",param))+theme_seal(), yr=c(0,1)))
    
    # Stage-specific baseline on logit scale
    logit_v <- switch(stage,
                      pup     = as.numeric(posterior::extract_variable(all_draws, "phi_pup_logit")),
                      juv     = qlogis(as.numeric(posterior::extract_variable(all_draws, "phi_juv_base"))),
                      adult_F = as.numeric(posterior::extract_variable(all_draws, "phi_adult_F_logit")),
                      adult_M = {
                        aF  <- as.numeric(posterior::extract_variable(all_draws, "phi_adult_F_base"))
                        del <- as.numeric(posterior::extract_variable(all_draws, "delta_adult"))
                        qlogis(pmax(aF - del, 0.001))
                      },
                      fecund  = {
                        # Logit-scale baseline for avg_fecundity
                        # Shows how OND MOCI shifts fecundity (prob of pupping) from its posterior mean
                        af <- as.numeric(posterior::extract_variable(all_draws, "avg_fecundity"))
                        qlogis(pmax(pmin(af, 0.999), 0.001))
                      }
    )
    
    y_label <- switch(stage,
                      pup     = "Pup Survival",
                      juv     = "Juvenile Survival",
                      adult_F = "Adult Female Survival",
                      adult_M = "Adult Male Survival",
                      fecund  = "Fecundity (prob. pupping)"
    )
    
    idx <- sample(seq_along(beta_v), min(500L, length(beta_v)))
    df  <- do.call(rbind, lapply(idx, function(i)
      data.frame(cov_val  = xr,
                 survival = plogis(logit_v[i] + beta_v[i] * xr))))
    agg <- aggregate(survival ~ cov_val, data = df,
                     FUN = function(v) c(mean = mean(v),
                                         lo   = as.numeric(quantile(v, CI_LO)),
                                         hi   = as.numeric(quantile(v, CI_HI))))
    sm  <- data.frame(cov_val = agg$cov_val,
                      mn      = agg$survival[, "mean"],
                      lo      = agg$survival[, "lo"],
                      hi      = agg$survival[, "hi"])
    
    clr       <- ifelse(mean(beta_v) < 0, "red3", "blue3")
    base_surv <- mean(plogis(logit_v))
    
    p <- ggplot(sm, aes(x = .data[["cov_val"]])) +
      geom_ribbon(aes(ymin = .data[["lo"]], ymax = .data[["hi"]]),
                  alpha = .2, fill = clr) +
      geom_line(aes(y = .data[["mn"]]), linewidth = 1.2, color = clr) +
      geom_hline(yintercept = base_surv, linetype = 2, color = "gray50") +
      geom_vline(xintercept = 0,         linetype = 2, color = "gray50") +
      labs(x = xlab, y = y_label, title = title) +
      theme_seal()
    if (!is.null(ylims)) p <- p + coord_cartesian(ylim = ylims)
    list(plot = p, yr = c(min(sm$lo), max(sm$hi)))
  }
  
  # ── PUP: coyote effects ─────────────────────────────────────────────────────
  coy <- lapply(1:3, function(i) make_plot(paste0("beta_coy[",i,"]"),
                                           "Coyote (SD)", paste0("Coyote → Pup Survival (",c("BL","DE","DP")[i],")")))
  ylc <- c(max(0,min(sapply(coy,`[[`,"yr"))-0.03), min(1,max(sapply(coy,`[[`,"yr"))+0.03))
  p_coy <- wrap_plots(lapply(1:3, function(i)
    make_plot(paste0("beta_coy[",i,"]"),"Coyote (SD)",
              paste0("(",c("BL","DE","DP")[i],")"), ylims=ylc)$plot), ncol=3) +
    plot_annotation(title="Site-Specific Coyote Effects on Pup Survival")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_effects_coyote.jpeg"),
                   p_coy, width=36, height=12, units="cm")
  
  # ── PUP: disturbance effects ────────────────────────────────────────────────
  dst <- lapply(1:6, function(s) make_plot(paste0("beta_dist_surv[",s,"]"),
                                           "Disturbance (SD)", paste0("(",site_names[s],")")))
  yld <- c(max(0,min(sapply(dst,`[[`,"yr"))-0.03), min(1,max(sapply(dst,`[[`,"yr"))+0.03))
  p_dst <- wrap_plots(lapply(1:6, function(s)
    make_plot(paste0("beta_dist_surv[",s,"]"),"Disturbance (SD)",
              paste0("(",site_names[s],")"), ylims=yld)$plot), ncol=3) +
    plot_annotation(title="Site-Specific Disturbance Effects on Pup Survival")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_effects_disturbance.jpeg"),
                   p_dst, width=36, height=24, units="cm")
  
  # ── PUP: shared MOCI + elephant seal ───────────────────────────────────────
  mci <- lapply(list(c("beta_moci_ond_fecund","MOCI Fall (SD)","MOCI Fall (OND) → Fecundity","fecund"),
                     c("beta_moci_ond_pup",   "MOCI Fall (SD)","MOCI Fall OND → Pup Survival","pup"),
                     c("beta_moci_jfm_pup",   "MOCI Winter (SD)","MOCI Winter JFM → Pup Survival","pup"),
                     c("beta_moci_amj_pup","MOCI Spring AMJ (SD)","MOCI Spring AMJ (t-1) → Pup","pup"),
                     c("beta_eseal_pup","Elephant Seal (SD)","Elephant Seal → Pup Survival")),
                function(x) make_plot(x[1],x[2],x[3],stage=if(length(x)>=4)x[4] else "pup"))
  ylm <- c(max(0,min(sapply(mci,`[[`,"yr"))-0.03), min(1,max(sapply(mci,`[[`,"yr"))+0.03))
  p_mci <- wrap_plots(lapply(list(c("beta_moci_ond_fecund","MOCI Fall (SD)","OND → Fecundity","fecund"),
                                  c("beta_moci_ond_pup",   "MOCI Fall (SD)","OND → Pup","pup"),
                                  c("beta_moci_jfm_pup",   "MOCI Winter (SD)","JFM → Pup","pup"),
                                  c("beta_moci_amj_pup","MOCI Spring AMJ (SD)","AMJ → Pup","pup"),
                                  c("beta_eseal_pup","Elephant Seal (SD)","Eseal")),
                             function(x) make_plot(x[1],x[2],x[3],stage=if(length(x)>=4)x[4] else "pup",ylims=ylm)$plot), ncol=3) +
    plot_annotation(title="Shared Covariate Effects on Pup Survival")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_effects_moci.jpeg"),
                   p_mci, width=36, height=12, units="cm")
  
  # ── JUVENILE: MOCI winter effect ────────────────────────────────────────────
  juv_moci <- make_plot("beta_moci_jfm_juv", "MOCI Winter (SD)",
                        "MOCI Winter → Juvenile Survival", stage="juv")
  p_juv <- juv_moci$plot
  if (save) ggsave(paste0("Output/Plots/",prefix,"_effects_juv_moci.jpeg"),
                   p_juv, width=20, height=14, units="cm")
  
  # ── ADULT: MOCI winter effect — female and male ─────────────────────────────
  adF_moci <- make_plot("beta_moci_jfm_adult", "MOCI Winter (SD)",
                        "MOCI Winter → Adult Female Survival", stage="adult_F")
  adM_moci <- make_plot("beta_moci_jfm_adult", "MOCI Winter (SD)",
                        "MOCI Winter → Adult Male Survival",   stage="adult_M")
  
  # Shared y-axis range across adult panels
  yla <- c(max(0, min(adF_moci$yr[1], adM_moci$yr[1]) - 0.02),
           min(1, max(adF_moci$yr[2], adM_moci$yr[2]) + 0.02))
  p_adult <- wrap_plots(
    make_plot("beta_moci_jfm_adult","MOCI Winter (SD)",
              "Female", stage="adult_F", ylims=yla)$plot,
    make_plot("beta_moci_jfm_adult","MOCI Winter (SD)",
              "Male",   stage="adult_M", ylims=yla)$plot,
    ncol=2) +
    plot_annotation(title="MOCI Winter Effects on Adult Survival")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_effects_adult_moci.jpeg"),
                   p_adult, width=28, height=14, units="cm")
  
  # ── COMBINED: all MOCI effects across stages ─────────────────────────────────
  p_all_moci <- wrap_plots(
    make_plot("beta_moci_ond_fecund",   "MOCI Fall (SD)",    "Fecundity (OND)",     stage="fecund")$plot,
    make_plot("beta_moci_ond_pup",      "MOCI Fall (SD)",    "Pup surv (OND)",      stage="pup")$plot,
    make_plot("beta_moci_jfm_pup",      "MOCI Winter (SD)",  "Pup surv (JFM)",      stage="pup")$plot,
    make_plot("beta_moci_amj_pup",      "MOCI Spring (SD)",  "Pup surv (AMJ t-1)",  stage="pup")$plot,
    make_plot("beta_moci_jfm_juv",   "MOCI Winter (SD)", "Juvenile",          stage="juv")$plot,
    make_plot("beta_moci_jfm_adult", "MOCI Winter (SD)", "Adult Female",      stage="adult_F")$plot,
    make_plot("beta_moci_jfm_adult", "MOCI Winter (SD)", "Adult Male",        stage="adult_M")$plot,
    ncol=3) +
    plot_annotation(title="MOCI Effects Across All Life Stages")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_effects_moci_all_stages.jpeg"),
                   p_all_moci, width=36, height=24, units="cm")
  
  list(coyote        = p_coy,
       disturbance   = p_dst,
       moci_pup      = p_mci,
       juv_moci      = p_juv,
       adult_moci    = p_adult,
       moci_all      = p_all_moci)
}


# ============================================================================
# PART 11c: JUVENILE + ADULT SURVIVAL EFFECTS
# ============================================================================
# Dedicated function for stage-specific MOCI effects on juvenile and adult
# survival. Complements create_effect_plots_v3.2 (which focuses on pup stage).
# Produces three panels: juvenile, adult female, adult male — all on a shared
# y-axis range for direct cross-stage comparison.

create_juv_adult_effect_plots_v3.2 <- function(fit, save=TRUE, prefix="IPM_v3.2") {
  
  all_draws <- tryCatch(
    fit$draws(format = "draws_df"),
    error = function(e) tryCatch(
      fit$draws(format = "matrix"),
      error = function(e2) stop(
        "Cannot extract draws from fit object.\n",
        "If fit was loaded from RDS use: out <- load_seal_results(\"IPM_v3.2_real\")\n",
        "then pass out$fit rather than a saved list element.\n",
        "Original error: ", conditionMessage(e))
    )
  )
  
  # Reuse the same make_plot engine as create_effect_plots_v3.2
  # by defining a local copy with the same implementation
  make_plot_stage <- function(param, xlab, title, stage,
                              xr=seq(-2,2,length.out=100), ylims=NULL) {
    beta_v <- tryCatch(
      as.numeric(posterior::extract_variable(all_draws, param)),
      error = function(e) { warning(sprintf("No draws for '%s'", param)); numeric(0) }
    )
    if (length(beta_v) == 0L)
      return(list(plot=ggplot()+labs(title=paste("MISSING:",param))+theme_seal(), yr=c(0,1)))
    
    logit_v <- switch(stage,
                      juv     = qlogis(as.numeric(posterior::extract_variable(all_draws, "phi_juv_base"))),
                      adult_F = as.numeric(posterior::extract_variable(all_draws, "phi_adult_F_logit")),
                      adult_M = {
                        aF  <- as.numeric(posterior::extract_variable(all_draws, "phi_adult_F_base"))
                        del <- as.numeric(posterior::extract_variable(all_draws, "delta_adult"))
                        qlogis(pmax(aF - del, 0.001))
                      }
    )
    y_label <- switch(stage,
                      juv     = "Juvenile Survival",
                      adult_F = "Adult Female Survival",
                      adult_M = "Adult Male Survival"
    )
    
    idx <- sample(seq_along(beta_v), min(500L, length(beta_v)))
    df  <- do.call(rbind, lapply(idx, function(i)
      data.frame(cov_val  = xr,
                 survival = plogis(logit_v[i] + beta_v[i] * xr))))
    agg <- aggregate(survival ~ cov_val, data=df,
                     FUN=function(v) c(mean=mean(v),
                                       lo=as.numeric(quantile(v, CI_LO)),
                                       hi=as.numeric(quantile(v, CI_HI))))
    sm  <- data.frame(cov_val=agg$cov_val,
                      mn=agg$survival[,"mean"],
                      lo=agg$survival[,"lo"],
                      hi=agg$survival[,"hi"])
    
    clr       <- ifelse(mean(beta_v) < 0, "red3", "blue3")
    base_surv <- mean(plogis(logit_v))
    
    p <- ggplot(sm, aes(x=.data[["cov_val"]])) +
      geom_ribbon(aes(ymin=.data[["lo"]], ymax=.data[["hi"]]), alpha=.2, fill=clr) +
      geom_line(aes(y=.data[["mn"]]), linewidth=1.2, color=clr) +
      geom_hline(yintercept=base_surv, linetype=2, color="gray50") +
      geom_vline(xintercept=0,         linetype=2, color="gray50") +
      labs(x=xlab, y=y_label, title=title) + theme_seal()
    if (!is.null(ylims)) p <- p + coord_cartesian(ylim=ylims)
    list(plot=p, yr=c(min(sm$lo), max(sm$hi)))
  }
  
  # ── MOCI winter → juvenile survival ─────────────────────────────────────────
  j1 <- make_plot_stage("beta_moci_jfm_juv", "MOCI Winter JFM (SD)",
                        "MOCI Winter → Juvenile Survival", stage="juv")
  
  # ── MOCI winter → adult survival (female + male, shared y-axis) ──────────────
  a1 <- make_plot_stage("beta_moci_jfm_adult", "MOCI Winter JFM (SD)",
                        "Adult Female", stage="adult_F")
  a2 <- make_plot_stage("beta_moci_jfm_adult", "MOCI Winter JFM (SD)",
                        "Adult Male",   stage="adult_M")
  
  # Shared y-axis across adult panels
  yla <- c(max(0, min(a1$yr[1], a2$yr[1]) - 0.02),
           min(1, max(a1$yr[2], a2$yr[2]) + 0.02))
  a1p <- make_plot_stage("beta_moci_jfm_adult","MOCI Winter JFM (SD)",
                         "Adult Female",stage="adult_F",ylims=yla)$plot
  a2p <- make_plot_stage("beta_moci_jfm_adult","MOCI Winter JFM (SD)",
                         "Adult Male",  stage="adult_M",ylims=yla)$plot
  
  p_adult <- wrap_plots(a1p, a2p, ncol=2) +
    plot_annotation(title="MOCI Winter Effects on Adult Survival")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_effects_adult_moci_detail.jpeg"),
                   p_adult, width=28, height=14, units="cm")
  
  # ── Cross-stage comparison: MOCI winter on juvenile vs adult ─────────────────
  # Shared y-axis across all three stages
  all_yr <- c(j1$yr, a1$yr, a2$yr)
  yall   <- c(max(0, min(all_yr)-0.02), min(1, max(all_yr)+0.02))
  
  p_compare <- wrap_plots(
    make_plot_stage("beta_moci_jfm_juv",   "MOCI Winter JFM (SD)", "Juvenile",
                    stage="juv",     ylims=yall)$plot,
    make_plot_stage("beta_moci_jfm_adult", "MOCI Winter JFM (SD)", "Adult Female",
                    stage="adult_F", ylims=yall)$plot,
    make_plot_stage("beta_moci_jfm_adult", "MOCI Winter JFM (SD)", "Adult Male",
                    stage="adult_M", ylims=yall)$plot,
    ncol=3) +
    plot_annotation(
      title   = "MOCI Winter (JFM) Effects: Juvenile and Adult Survival",
      subtitle= paste0("Shared y-axis for direct cross-stage comparison; ",
                       "dashed line = posterior mean baseline; ", CI_LABEL))
  if (save) ggsave(paste0("Output/Plots/",prefix,"_effects_juv_adult_moci_comparison.jpeg"),
                   p_compare, width=36, height=14, units="cm")
  
  list(juv_moci    = j1$plot,
       adult_moci  = p_adult,
       comparison  = p_compare)
}



# ============================================================================
# PART 12: SUMMARY TABLE
# ============================================================================

create_summary_table_v3.2 <- function(fit, save=TRUE, prefix="IPM_v3.2") {
  
  params <- c(
    "phi_pup_logit","phi_juv_base","phi_adult_F_logit","phi_adult_F_base","delta_adult",
    "fecund_primip","fecund_mature",
    "prop_female","avg_fecundity","p_male_breed",
    "beta_coy[1]","beta_coy[2]","beta_coy[3]",
    paste0("beta_dist_surv[",1:6,"]"),
    "beta_moci_ond_fecund","beta_moci_ond_pup","beta_moci_amj_pup",
    "beta_moci_jfm_pup","beta_moci_jfm_juv",
    "beta_moci_jfm_adult","beta_eseal_pup","beta_moci_amj_molt",
    "detect_breed_logit","detect_molt_logit",
    "sigma_process","sigma_obs_adult","sigma_obs_pup","sigma_obs_molt","sigma_site"
  )
  
  # Derive phi_pup_base on probability scale for reporting
  pup_draws <- fit$draws(variables="phi_pup_logit",format="df")$phi_pup_logit
  pup_base_summary <- tibble(
    variable="phi_pup_base (prob)",
    mean=mean(plogis(pup_draws)), sd=sd(plogis(pup_draws)),
    q_lo=quantile(plogis(pup_draws),CI_LO), q_hi=quantile(plogis(pup_draws),CI_HI),
    rhat=NA_real_, ess_bulk=NA_real_
  )
  
  tbl <- seal_fit_summary(fit, params) |>
    bind_rows(pup_base_summary) |>
    mutate(
      Parameter = case_when(
        variable=="phi_pup_logit"        ~ "Pup survival, logit scale (sex-neutral)",
        variable=="phi_pup_base (prob)"  ~ "Pup survival, probability scale (derived)",
        variable=="phi_juv_base"         ~ "Juvenile survival (sex-neutral)",
        variable=="phi_adult_F_logit"    ~ "Adult female survival (logit scale)",
        variable=="phi_adult_F_base"     ~ "Adult female survival (probability scale, derived)",
        variable=="delta_adult"          ~ "Adult M-F survival difference",
        variable=="fecund_primip"        ~ "Fecundity: primiparous (4-5 yr)",
        variable=="fecund_mature"        ~ "Fecundity: mature (6+ yr, experienced breeders)",
        variable=="prop_female"          ~ "Proportion female at birth",
        variable=="avg_fecundity"        ~ "Average fecundity (weighted)",
        variable=="p_male_breed"         ~ "Male haul-out fraction (breeding season)",
        variable=="beta_coy[1]"          ~ "Coyote → pup survival (BL)",
        variable=="beta_coy[2]"          ~ "Coyote → pup survival (DE)",
        variable=="beta_coy[3]"          ~ "Coyote → pup survival (DP)",
        variable=="beta_dist_surv[1]"    ~ "Disturbance → pup survival (BL)",
        variable=="beta_dist_surv[2]"    ~ "Disturbance → pup survival (DE)",
        variable=="beta_dist_surv[3]"    ~ "Disturbance → pup survival (DP)",
        variable=="beta_dist_surv[4]"    ~ "Disturbance → pup survival (PRH)",
        variable=="beta_dist_surv[5]"    ~ "Disturbance → pup survival (TB)",
        variable=="beta_dist_surv[6]"    ~ "Disturbance → pup survival (TP)",
        variable=="beta_moci_ond_fecund"    ~ "MOCI Fall (OND) → fecundity",
        variable=="beta_moci_ond_pup"    ~ "MOCI fall OND (t) → pup survival",
        variable=="beta_moci_jfm_pup"    ~ "MOCI winter JFM (t) → pup survival",
        variable=="beta_moci_amj_pup"    ~ "MOCI spring AMJ (t-1) → pup survival",
        variable=="beta_moci_jfm_juv"    ~ "MOCI winter → juvenile survival",
        variable=="beta_moci_jfm_adult"  ~ "MOCI winter → adult survival",
        variable=="beta_eseal_pup"       ~ "Elephant seal → pup survival",
        variable=="beta_moci_amj_molt"   ~ "MOCI spring → molt detection",
        variable=="sigma_process"        ~ "Process error (σ)",
        variable=="sigma_obs_adult"      ~ "Obs error adult (σ)",
        variable=="sigma_obs_pup"        ~ "Obs error pup (σ)",
        variable=="sigma_obs_molt"       ~ "Obs error molt (σ)",
        variable=="detect_breed_logit"  ~ "Detect breed baseline (logit)",
        variable=="detect_molt_logit"   ~ "Detect molt baseline (logit)",
        variable=="sigma_site"           ~ "Site random effect (σ)",
        TRUE ~ variable
      ),
      Estimate = sprintf("%.3f (%.3f, %.3f)", mean, q_lo, q_hi),
      Category = case_when(
        str_detect(variable,"phi|delta")              ~ "Survival",
        str_detect(variable,"fecund|prop|avg")        ~ "Reproduction",
        variable %in% c("p_male_breed","phi_pup_base (prob)") ~ "Observation / Derived",
        str_detect(variable,"beta_coy")               ~ "Coyote (site-specific)",
        str_detect(variable,"beta_dist")              ~ "Disturbance (site-specific)",
        str_detect(variable,"beta_moci|beta_eseal")   ~ "Shared covariates",
        variable %in% c("detect_breed_logit",
                        "detect_molt_logit")          ~ "Observation / Derived",
        str_detect(variable,"sigma")                  ~ "Error terms"
      )
    ) |>
    select(Category, Parameter, Estimate, rhat, ess_bulk)
  
  cat("\n============================================\n")
  cat("PARAMETER SUMMARY — IPM v3.2\n")
  cat("============================================\n\n")
  print(tbl, n=nrow(tbl))
  
  if (save) write_csv(tbl, paste0("Output/",prefix,"_parameter_summary.csv"))
  tbl
}


# ============================================================================
# PART 13: SAVE OUTPUT
# ============================================================================

save_model_output_v3.2 <- function(fit, prefix="IPM_v3.2") {
  
  fit$save_object(paste0("Output/harbor_seal_",prefix,"_fit.rds"))
  cat(sprintf("Fit saved: Output/harbor_seal_%s_fit.rds\n", prefix))
  
  surv  <- seal_fit_summary(fit, c("phi_juv_base","phi_adult_F_base","delta_adult"))
  logit_aF <- fit$draws(variables="phi_adult_F_logit",format="df")$phi_adult_F_logit
  pup_l <- fit$draws(variables="phi_pup_logit",format="df")$phi_pup_logit
  pmb   <- seal_fit_summary(fit, "p_male_breed")
  coy   <- seal_fit_summary(fit, paste0("beta_coy[",1:3,"]"))
  dst   <- seal_fit_summary(fit, paste0("beta_dist_surv[",1:6,"]"))
  sns   <- c("BL","DE","DP","PRH","TB","TP")
  
  cat("\n============================================\n")
  cat("HARBOR SEAL IPM v3.2 — KEY RESULTS\n")
  cat("============================================\n")
  cat(sprintf("Pup survival (sex-neutral, prob):  %.3f (%.3f–%.3f)\n",
              median(plogis(pup_l)), quantile(plogis(pup_l),CI_LO), quantile(plogis(pup_l),CI_HI)))
  cat(sprintf("Juvenile survival (sex-neutral):   %.3f (%.3f–%.3f)\n",
              surv$mean[1],surv$q_lo[1],surv$q_hi[1]))
  cat(sprintf("Adult female survival (baseline):  %.3f (%.3f–%.3f)\n",
              surv$mean[2],surv$q_lo[2],surv$q_hi[2]))
  cat(sprintf("Adult male survival (F - delta):   %.3f (%.3f–%.3f)\n",
              surv$mean[2]-surv$mean[3], surv$q_lo[2]-surv$q_hi[3], surv$q_hi[2]-surv$q_lo[3]))
  cat(sprintf("Male breeding haul-out (p_male_breed): %.3f (%.3f–%.3f)\n",
              pmb$mean,pmb$q_lo,pmb$q_hi))
  cat("\nCoyote effects on pup survival:\n")
  for (i in 1:3) cat(sprintf("  %s: %.3f (%.3f–%.3f)\n",c("BL","DE","DP")[i],
                             coy$mean[i],coy$q_lo[i],coy$q_hi[i]))
  moci <- seal_fit_summary(fit, c("beta_moci_jfm_juv","beta_moci_jfm_adult",
                                  "beta_moci_ond_fecund","beta_moci_amj_pup"))
  
  cat("\nDisturbance effects on pup survival:\n")
  for (i in 1:6) cat(sprintf("  %s: %.3f (%.3f–%.3f)\n",sns[i],
                             dst$mean[i],dst$q_lo[i],dst$q_hi[i]))
  cat("\nMOCI effects by stage (logit scale; negative = warm MOCI reduces survival):\n")
  cat(sprintf("  Pup    OND lag:  %.3f (%.3f–%.3f)\n",
              moci$mean[3], moci$q_lo[3], moci$q_hi[3]))
  # Detection baselines — only present in updated model; skip gracefully if absent
  tryCatch({
    dtb <- seal_fit_summary(fit, c("detect_breed_logit","detect_molt_logit"))
    cat("\nDetection baselines (logit scale; plogis = detection probability):\n")
    cat(sprintf("  Breeding survey: %.3f (%.3f\u2013%.3f) -> p = %.2f\n",
                dtb$mean[1],dtb$q_lo[1],dtb$q_hi[1],plogis(dtb$mean[1])))
    cat(sprintf("  Molt survey:     %.3f (%.3f\u2013%.3f) -> p = %.2f\n",
                dtb$mean[2],dtb$q_lo[2],dtb$q_hi[2],plogis(dtb$mean[2])))
  }, error = function(e)
    cat("  NOTE: detect_breed_logit / detect_molt_logit not in this model fit\n"))
  cat("\nMOCI effects by stage:\n")
  cat(sprintf("  Pup    AMJ (t-1): %.3f (%.3f–%.3f)\n",
              moci$mean[4], moci$q_lo[4], moci$q_hi[4]))
  # New pup MOCI effects — only present in updated model; skip gracefully if absent
  tryCatch({
    pup_moci2 <- seal_fit_summary(fit, c("beta_moci_ond_pup", "beta_moci_jfm_pup"))
    cat(sprintf("  Pup    OND post-wean:  %.3f (%.3f–%.3f)\n",
                pup_moci2$mean[1], pup_moci2$q_lo[1], pup_moci2$q_hi[1]))
    cat(sprintf("  Pup    JFM 1st-winter: %.3f (%.3f–%.3f)\n",
                pup_moci2$mean[2], pup_moci2$q_lo[2], pup_moci2$q_hi[2]))
  }, error = function(e)
    cat("  NOTE: beta_moci_ond_pup / beta_moci_jfm_pup not in this model fit\n"))
  cat(sprintf("  Juv    JFM:      %.3f (%.3f–%.3f)\n",
              moci$mean[1], moci$q_lo[1], moci$q_hi[1]))
  cat(sprintf("  Adult  JFM:      %.3f (%.3f–%.3f)\n",
              moci$mean[2], moci$q_lo[2], moci$q_hi[2]))
  cat("============================================\n")
}


# ============================================================================
# NOTE: run_full_analysis_v3.2() is defined in harbor_seal_ipm_v3.2.R.
# Source that file before this one. A second definition here would silently
# overwrite the canonical version (which includes run_portfolio,
# run_synchrony, and safe_run tryCatch wrapping). Do not redefine it here.
# ============================================================================


# ============================================================================
# PART 15: PORTFOLIO ANALYSIS
# ============================================================================

create_portfolio_analysis_v3.2 <- function(fit, sim_data, save=TRUE, prefix="IPM_v3.2") {
  
  years <- sim_data$years; site_names <- sim_data$site_names
  S <- length(site_names); T <- length(years)
  
  ldraws <- fit$draws(variables="lambda", format="matrix")
  lmat   <- matrix(NA,S,T-1, dimnames=list(site_names,years[1:(T-1)]))
  for (s in 1:S) for (t in 1:(T-1)) {
    cn <- paste0("lambda[",s,",",t,"]")
    if (cn %in% colnames(ldraws)) lmat[s,t] <- mean(ldraws[,cn])
  }
  
  cv_meta  <- sd(colMeans(lmat,na.rm=TRUE)) / mean(colMeans(lmat,na.rm=TRUE))
  cv_sites <- mean(apply(lmat,1,function(x) sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)))
  per      <- cv_meta / cv_sites
  lcor     <- cor(t(lmat), use="pairwise.complete.obs")
  
  cat(sprintf("Portfolio Effect Ratio: %.3f  (< 1 = buffering)\n",per))
  cat(sprintf("Mean site correlation:  %.3f\n",mean(lcor[lower.tri(lcor)])))
  
  phi_d <- fit$draws(variables="phi_pup", format="matrix")
  phi_m <- matrix(NA,S,T, dimnames=list(site_names,years))
  for (s in 1:S) for (t in 1:T) {
    cn <- paste0("phi_pup[",s,",",t,"]")
    if (cn %in% colnames(phi_d)) phi_m[s,t] <- mean(phi_d[,cn])
  }
  
  ldf <- expand.grid(Site=site_names,Year=years[1:(T-1)]) |>
    mutate(Site=factor(Site,levels=site_names), lambda=as.vector(t(lmat)))
  
  p_heat <- ggplot(ldf,aes(x=Year,y=Site,fill=lambda)) +
    geom_tile(color="white",linewidth=0.5) +
    scale_fill_gradient2(low="red3",mid="white",high="darkgreen",
                         midpoint=1,limits=c(0.7,1.3),oob=scales::squish,name="λ") +
    geom_text(aes(label=sprintf("%.2f",lambda)),size=2.5) +
    labs(x="Year",y="Site",title="Site-Specific λ by Year") +
    theme_seal() +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=8),panel.grid=element_blank())
  if (save) ggsave(paste0("Output/Plots/",prefix,"_portfolio_lambda_heatmap.jpeg"),
                   p_heat, width=35, height=15, units="cm")
  
  phidf <- expand.grid(Site=site_names,Year=years) |>
    mutate(Site=factor(Site,levels=site_names), phi=as.vector(t(phi_m)))
  p_phi <- ggplot(phidf,aes(x=Year,y=Site,fill=phi)) +
    geom_tile(color="white",linewidth=0.5) +
    scale_fill_viridis_c(name="φ_pup",option="plasma",
                         limits=c(0.25,0.75),oob=scales::squish) +
    geom_text(aes(label=sprintf("%.2f",phi)),size=2.2) +
    labs(x="Year",y="Site",title="Site-Specific Pup Survival by Year (sex-neutral)") +
    theme_seal() +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=8),panel.grid=element_blank())
  if (save) ggsave(paste0("Output/Plots/",prefix,"_portfolio_pup_survival_heatmap.jpeg"),
                   p_phi, width=35, height=15, units="cm")
  
  p_async <- ggplot(ldf,aes(x=Year,y=lambda,color=Site,group=Site)) +
    geom_hline(yintercept=1,linetype=2,color="gray50") +
    geom_line(linewidth=1,alpha=0.7) + geom_point(size=2) +
    scale_color_brewer(palette="Dark2") +
    scale_fill_brewer(palette="Dark2") +
    labs(x="Year",y=expression(lambda),title="Site-Level Population Growth — Asynchrony") +
    theme_seal() + theme(legend.position="bottom")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_portfolio_asynchrony.jpeg"),
                   p_async, width=30, height=18, units="cm")
  
  bw <- tibble(Year=years[1:(T-1)],
               Best_Site =site_names[apply(lmat,2,which.max)],
               Worst_Site=site_names[apply(lmat,2,which.min)])
  p_bw <- bw |>
    pivot_longer(c(Best_Site,Worst_Site),names_to="Performance",values_to="Site") |>
    mutate(Performance=recode(Performance,Best_Site="Best",Worst_Site="Worst"),
           Site=factor(Site,levels=site_names)) |>
    ggplot(aes(x=Year,y=Site,fill=Performance)) +
    geom_tile(color="white",linewidth=1,alpha=0.8) +
    scale_fill_manual(values=c(Best="darkgreen",Worst="red3")) +
    labs(x="Year",y="Site",title="Best and Worst Performing Sites Each Year") +
    theme_seal() +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=8),panel.grid=element_blank())
  if (save) ggsave(paste0("Output/Plots/",prefix,"_portfolio_best_worst.jpeg"),
                   p_bw, width=35, height=12, units="cm")
  
  cordf <- expand.grid(Site1=site_names,Site2=site_names) |>
    mutate(r=as.vector(lcor))
  p_cor <- ggplot(cordf,aes(x=Site1,y=Site2,fill=r)) +
    geom_tile(color="white") + geom_text(aes(label=sprintf("%.2f",r)),size=4) +
    scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0,limits=c(-1,1)) +
    labs(x="",y="",title="Between-Site Correlation in λ") +
    theme_seal() + theme(panel.grid=element_blank()) + coord_fixed()
  if (save) ggsave(paste0("Output/Plots/",prefix,"_portfolio_correlation.jpeg"),
                   p_cor, width=20, height=18, units="cm")
  
  # ── Portfolio tables ─────────────────────────────────────────────────────────
  
  # Table 1: summary statistics
  port_summary_tbl <- tibble(
    Metric    = c("Portfolio Effect Ratio (CV_meta / CV_sites)",
                  "Mean between-site correlation (lambda)",
                  "Min between-site correlation (lambda)",
                  "Max between-site correlation (lambda)",
                  "Overall mean lambda (across sites and years)",
                  "SD of annual mean lambda"),
    Value     = c(round(per, 3),
                  round(mean(lcor[lower.tri(lcor)]), 3),
                  round(min(lcor[lower.tri(lcor)]), 3),
                  round(max(lcor[lower.tri(lcor)]), 3),
                  round(mean(lmat, na.rm=TRUE), 3),
                  round(sd(colMeans(lmat, na.rm=TRUE)), 3)),
    Interpretation = c(
      if (per < 1) "Portfolio buffering present (< 1)" else "No buffering (>= 1)",
      if (mean(lcor[lower.tri(lcor)]) < 0.3) "Low synchrony — strong portfolio effect"
      else if (mean(lcor[lower.tri(lcor)]) < 0.6) "Moderate synchrony"
      else "High synchrony — weak portfolio effect",
      "—", "—",
      if (mean(lmat,na.rm=TRUE) > 1) "Population growing on average"
      else "Population declining on average",
      "—"
    )
  )
  cat("
── Portfolio Summary Table ──────────────────────────────────────
")
  print(port_summary_tbl, n=nrow(port_summary_tbl))
  if (save) write_csv(port_summary_tbl,
                      paste0("Output/",prefix,"_portfolio_summary.csv"))
  
  # Table 2: mean lambda per site
  lambda_site_tbl <- as.data.frame(lmat) |>
    tibble::rownames_to_column("Site") |>
    mutate(
      Mean_lambda = round(rowMeans(across(where(is.numeric)), na.rm=TRUE), 3),
      SD_lambda   = round(apply(across(where(is.numeric)), 1,
                                function(x) sd(x, na.rm=TRUE)), 3),
      CV_lambda   = round(SD_lambda / Mean_lambda, 3),
      # table() must be subsetted BEFORE as.integer() — as.integer() strips names,
      # making named lookup impossible and returning NA → 0 for every site.
      Times_best  = as.integer(table(bw$Best_Site) [Site]) |> replace_na(0L),
      Times_worst = as.integer(table(bw$Worst_Site)[Site]) |> replace_na(0L)
    ) |>
    select(Site, Mean_lambda, SD_lambda, CV_lambda, Times_best, Times_worst)
  cat("
── Per-Site Lambda Table ────────────────────────────────────────
")
  print(lambda_site_tbl)
  if (save) write_csv(lambda_site_tbl,
                      paste0("Output/",prefix,"_portfolio_lambda_by_site.csv"))
  
  # Table 3: between-site correlation matrix as a tidy table
  cor_long <- as.data.frame(lcor) |>
    tibble::rownames_to_column("Site1") |>
    pivot_longer(cols=-all_of("Site1"), names_to="Site2", values_to="r") |>
    filter(Site1 < Site2) |>
    mutate(r = round(r, 3),
           Interpretation = case_when(
             abs(r) < 0.3 ~ "Low",
             abs(r) < 0.6 ~ "Moderate",
             TRUE          ~ "High"
           )) |>
    arrange(desc(abs(r)))
  if (save) write_csv(cor_long,
                      paste0("Output/",prefix,"_portfolio_correlations.csv"))
  
  # ── Additional portfolio plots ────────────────────────────────────────────────
  
  # CV by site (bar chart)
  cv_site <- lambda_site_tbl |>
    ggplot(aes(x=reorder(Site, CV_lambda), y=CV_lambda, fill=Mean_lambda)) +
    geom_col(alpha=0.85) +
    geom_text(aes(label=sprintf("%.3f", CV_lambda)), vjust=-0.3, size=4) +
    scale_fill_gradient2(low="red3", mid="white", high="darkgreen",
                         midpoint=1, name=expression(bar(lambda))) +
    labs(x="Site", y="CV of lambda",
         title="Among-Year Variability in Population Growth by Site",
         subtitle="Lower CV = more stable site") +
    theme_seal() + theme(legend.position="right")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_portfolio_cv_by_site.jpeg"),
                   cv_site, width=22, height=14, units="cm")
  
  # Pair-wise correlation scatter: lambda at each site pair
  ldf_wide <- as.data.frame(t(lmat)) |>
    setNames(site_names)
  site_pairs <- combn(site_names, 2, simplify=FALSE)
  pair_plots <- lapply(site_pairs, function(pair) {
    d <- data.frame(x_val=ldf_wide[[pair[1]]], y_val=ldf_wide[[pair[2]]])
    r <- cor(d$x_val, d$y_val, use="complete.obs")
    ggplot(d, aes(x=.data[["x_val"]], y=.data[["y_val"]])) +
      geom_point(alpha=0.6, colour=SEAL_COLS$pop, size=2) +
      geom_smooth(method="lm", se=FALSE, colour="gray40", linewidth=0.7) +
      geom_abline(slope=1, intercept=0, linetype=2, colour="gray70") +
      labs(x=paste0(pair[1]," lambda"), y=paste0(pair[2]," lambda"),
           title=sprintf("%s vs %s  (r=%.2f)", pair[1], pair[2], r)) +
      theme_seal(base_size=12)
  })
  p_pairs <- wrap_plots(pair_plots, ncol=3) +
    plot_annotation(title="Pairwise Site Synchrony in Population Growth (λ)")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_portfolio_pairwise_lambda.jpeg"),
                   p_pairs, width=36, height=30, units="cm")
  
  list(lambda_heatmap=p_heat, phi_heatmap=p_phi, asynchrony=p_async,
       best_worst=p_bw, correlation=p_cor,
       cv_by_site=cv_site, pairwise=p_pairs,
       lambda_matrix=lmat, phi_matrix=phi_m,
       summary_table=port_summary_tbl,
       lambda_site_table=lambda_site_tbl,
       correlation_table=cor_long,
       summary=list(portfolio_effect_ratio=per,
                    mean_correlation=mean(lcor[lower.tri(lcor)]),
                    times_best=table(bw$Best_Site),
                    times_worst=table(bw$Worst_Site)))
}


# ============================================================================
# PART 16: SYNCHRONY PROJECTIONS
# ============================================================================

create_synchrony_projections_v3.2 <- function(fit, sim_data, n_sims=500,
                                              T_proj=10, save=TRUE, prefix="IPM_v3.2") {
  
  years <- sim_data$years; site_names <- sim_data$site_names
  S <- length(site_names); T <- length(years)
  pyrs <- (max(years)+1):(max(years)+T_proj)
  
  draws <- fit$draws(format="df")
  idx   <- sample(seq_len(nrow(draws)), min(n_sims,nrow(draws)))
  
  pup_l    <- draws$phi_pup_logit[idx]
  phi_juv  <- draws$phi_juv_base[idx]
  phi_aF_logit <- draws$phi_adult_F_logit[idx]
  delta_a  <- draws$delta_adult[idx]
  pf       <- draws$prop_female[idx]
  avgf     <- draws$avg_fecundity[idx]
  b_ond    <- draws$beta_moci_ond_fecund[idx]     # OND → fecundity
  b_amj    <- draws$beta_moci_amj_pup[idx]         # AMJ → pup survival
  b_jfmJ   <- draws$beta_moci_jfm_juv[idx]
  b_jfmA   <- draws$beta_moci_jfm_adult[idx]
  # New pup MOCI parameters — may be absent in old model fits; default to 0
  b_ond_pup <- tryCatch(draws$beta_moci_ond_pup[idx],  error=function(e) rep(0, length(idx)))
  b_jfm_pup <- tryCatch(draws$beta_moci_jfm_pup[idx],  error=function(e) rep(0, length(idx)))
  
  se   <- sapply(1:S, function(s) draws[[paste0("site_effect[",s,"]")]][idx])
  bcoy <- sapply(1:3, function(k) draws[[paste0("beta_coy[",k,"]")]][idx])
  cidx <- c(1,2,3,0,0,0)
  
  NAF <- sapply(1:S, function(s) draws[[paste0("N_adult_F[",s,",",T,"]")]][idx])
  NAM <- sapply(1:S, function(s) draws[[paste0("N_adult_M[",s,",",T,"]")]][idx])
  NJF <- sapply(1:S, function(s) draws[[paste0("N_juv_F[",s,",",T,"]")]][idx])
  NJM <- sapply(1:S, function(s) draws[[paste0("N_juv_M[",s,",",T,"]")]][idx])
  NP  <- sapply(1:S, function(s) draws[[paste0("N_pup[",s,",",T,"]")]][idx])
  
  scenarios <- list(
    list(name="Status Quo",        moci=0,  coyote=0),
    list(name="Cool (MOCI -1)",    moci=-1, coyote=0),
    list(name="Warm (MOCI +1)",    moci=1,  coyote=0),
    list(name="Warm + High Coyote",moci=1,  coyote=1)
  )
  
  run_proj <- function(sync, sc, psd=0.15) {
    Nt <- matrix(NA, length(idx), T_proj)
    for (i in seq_along(idx)) {
      naf<-NAF[i,]; nam<-NAM[i,]; njf<-NJF[i,]; njm<-NJM[i,]; np<-NP[i,]
      Nt[i,1] <- sum(naf+nam+njf+njm+np)
      for (tp in 2:T_proj) {
        sh <- if (sync) rep(rnorm(1,0,psd),S) else rnorm(S,0,psd)
        for (s in 1:S) {
          ce <- if (cidx[s]>0) bcoy[i,cidx[s]]*sc$coyote else 0
          # OND MOCI → fecundity (maternal condition), NOT phi_pup
          fecund_t <- plogis(qlogis(pmax(pmin(avgf[i],0.999),0.001)) + b_ond[i]*sc$moci)
          # phi_pup: AMJ + new OND pup-survival + JFM first-winter terms
          pp  <- plogis(pup_l[i]+se[i,s]+ce +
                          b_amj[i]*sc$moci +
                          b_ond_pup[i]*sc$moci +
                          b_jfm_pup[i]*sc$moci)
          pj  <- plogis(qlogis(phi_juv[i])+se[i,s]*0.5+b_jfmJ[i]*sc$moci)
          paF <- plogis(phi_aF_logit[i] + se[i,s]*0.25 + b_jfmA[i]*sc$moci)
          paM <- plogis(qlogis(pmax(plogis(phi_aF_logit[i]) - delta_a[i], 0.001)) +
                          se[i,s]*0.25 + b_jfmA[i]*sc$moci)
          new_p <- naf[s]*fecund_t*exp(sh[s])   # fecundity now MOCI-modulated
          njF2  <- np[s]*pf[i]*pp;       njM2 <- np[s]*(1-pf[i])*pp
          jsF   <- njf[s]*pj*(2/3);      jsM  <- njm[s]*pj*(2/3)
          jaF   <- njf[s]*pj*(1/3);      jaM  <- njm[s]*pj*(1/3)
          np[s]  <- max(new_p,1); njf[s] <- max(njF2+jsF,.1); njm[s] <- max(njM2+jsM,.1)
          naf[s] <- max(naf[s]*paF+jaF,1); nam[s] <- max(nam[s]*paM+jaM,1)
        }
        Nt[i,tp] <- sum(naf+nam+njf+njm+np)
      }
    }
    Nt
  }
  
  res <- lapply(scenarios, function(sc) {
    cat(sprintf("Projecting: %s\n",sc$name))
    list(async=run_proj(FALSE,sc), sync=run_proj(TRUE,sc))
  })
  names(res) <- sapply(scenarios,`[[`,"name")
  
  comp_df <- map_dfr(names(res), function(sn)
    bind_rows(
      tibble(Scenario=sn,Year=pyrs,Synchrony="Asynchronous (current)",
             mean=colMeans(res[[sn]]$async),
             lo=as.numeric(apply(res[[sn]]$async,2,quantile,CI_LO)),
             hi=as.numeric(apply(res[[sn]]$async,2,quantile,CI_HI))),
      tibble(Scenario=sn,Year=pyrs,Synchrony="Synchronous (hypothetical)",
             mean=colMeans(res[[sn]]$sync),
             lo=as.numeric(apply(res[[sn]]$sync,2,quantile,CI_LO)),
             hi=as.numeric(apply(res[[sn]]$sync,2,quantile,CI_HI)))
    ))
  
  cv_df <- map_dfr(names(res), function(sn) {
    ca <- sd(res[[sn]]$async[,T_proj])/mean(res[[sn]]$async[,T_proj])
    cs <- sd(res[[sn]]$sync[,T_proj]) /mean(res[[sn]]$sync[,T_proj])
    tibble(Scenario=sn, CV_Async=ca, CV_Sync=cs, CV_Ratio=ca/cs,
           Buffering_Pct=(1-ca/cs)*100)
  })
  cat("\n--- PORTFOLIO BUFFERING ---\n"); print(cv_df)
  
  p_comp <- ggplot(comp_df,aes(x=Year,y=mean,color=Synchrony,fill=Synchrony,linetype=Synchrony)) +
    geom_ribbon(aes(ymin=lo,ymax=hi),alpha=0.15,color=NA) + geom_line(linewidth=1.2) +
    facet_wrap(~Scenario,ncol=2,scales="free_y") +
    scale_color_manual(values=c("Asynchronous (current)"=SEAL_COLS$pop,"Synchronous (hypothetical)"=SEAL_COLS$adult_f)) +
    scale_fill_manual( values=c("Asynchronous (current)"=SEAL_COLS$ribbon,"Synchronous (hypothetical)"="#FFCDD2")) +
    scale_linetype_manual(values=c("Asynchronous (current)"="solid","Synchronous (hypothetical)"="dashed")) +
    labs(x="Year",y="Total Population",title="Portfolio Buffering: Async vs Sync") +
    theme_seal() + theme(legend.position="bottom")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_synchrony_comparison.jpeg"),
                   p_comp, width=30, height=25, units="cm")
  
  p_cv <- cv_df |>
    pivot_longer(c(CV_Async,CV_Sync),names_to="Type",values_to="CV") |>
    mutate(Type=recode(Type,CV_Async="Asynchronous",CV_Sync="Synchronous")) |>
    ggplot(aes(x=Scenario,y=CV,fill=Type)) +
    geom_col(position="dodge",alpha=0.8) +
    scale_fill_manual(values=c(Asynchronous=SEAL_COLS$pop,Synchronous=SEAL_COLS$adult_f)) +
    labs(x="Scenario",y="CV of aggregate abundance",title="Portfolio Buffering by Scenario") +
    theme_seal() +
    theme(axis.text.x=element_text(angle=20,hjust=1),legend.position="bottom")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_synchrony_cv.jpeg"),
                   p_cv, width=25, height=18, units="cm")
  
  if (save) {
    write_csv(comp_df, paste0("Output/",prefix,"_synchrony_projections.csv"))
    write_csv(cv_df,   paste0("Output/",prefix,"_synchrony_cv_comparison.csv"))
  }
  
  list(comparison=p_comp, cv_plot=p_cv,
       projection_data=comp_df, cv_comparison=cv_df, raw_projections=res)
}


# ============================================================================
# USAGE
# ============================================================================

cat("\n")
cat("============================================\n")
cat("IPM v3.2 LOADED\n")

# ============================================================================
# PART 14: COVARIATE DECOMPOSITION — STACKED BAR BY SITE × YEAR
# ============================================================================
# Shows how each covariate contributes to vital rates at each site by year.
# y-axis = covariate effect in logit units (beta × covariate value).
# Positive bars push survival/fecundity up; negative bars push it down.
# Uses posterior MEANS for coefficients; actual observed covariate values.
#
# Three panels produced:
#   (a) Pup survival: coyote, disturbance, eseal, MOCI AMJ — all at t-1 (birth year)
#   (b) Fecundity:    MOCI OND (at t = OND of year t-1, pre-lagged)
#   (c) All combined: all 7 covariate effects across all vital rates

create_covariate_decomposition_plots_v3.2 <- function(fit, sim_data,
                                                      save=TRUE,
                                                      prefix="IPM_v3.2") {
  
  years      <- sim_data$years
  site_names <- sim_data$site_names
  S          <- length(site_names)
  nyr        <- length(years)             # avoid masking T=TRUE
  
  # Covariate arrays — find where moci_amj lives (nested or flat layout)
  sdat <- if (!is.null(sim_data$stan_data) &&
              !is.null(sim_data$stan_data$moci_amj)) {
    sim_data$stan_data
  } else if (!is.null(sim_data$moci_amj)) {
    sim_data
  } else {
    # Try one level deeper
    sub <- Filter(function(x) is.list(x) && !is.null(x$moci_amj), sim_data)
    if (length(sub) > 0) sub[[1]]
    else stop(paste(
      "Decomp: cannot find moci_amj in sim_data.",
      "Available names:", paste(names(sim_data), collapse=", "),
      "\nTry: out <- load_seal_results('IPM_v3.2_real') and use out$sim_data"))
  }
  
  # ── Posterior means via fit$summary() ────────────────────────────────────
  # fit$summary() uses cmdstanr's cached summaries and does NOT require
  # re-reading CSV files, so it works robustly after session reloads and
  # when called from functions appended after the script's top-level block.
  # Core parameters always present in the model
  beta_vars_core <- c(
    paste0("beta_coy[",       1:3, "]"),
    paste0("beta_dist_surv[", 1:S, "]"),
    "beta_eseal_pup",
    "beta_moci_amj_pup",
    "beta_moci_ond_fecund",
    "beta_moci_jfm_juv",
    "beta_moci_jfm_adult"
  )
  bsum <- fit$summary(variables = beta_vars_core)
  # gm0: returns 0 (not numeric(0)) when variable absent — prevents zero-row tibbles
  gm   <- function(v) {
    r <- bsum$mean[bsum$variable == v]
    if (length(r) == 0L) 0 else r
  }
  
  b_coy  <- sapply(1:3, function(k) gm(paste0("beta_coy[",       k, "]")))
  b_dst  <- sapply(1:S, function(k) gm(paste0("beta_dist_surv[", k, "]")))
  b_es   <- gm("beta_eseal_pup")
  b_amj  <- gm("beta_moci_amj_pup")
  b_ond  <- gm("beta_moci_ond_fecund")
  b_jfmJ <- gm("beta_moci_jfm_juv")
  b_jfmA <- gm("beta_moci_jfm_adult")
  
  # New pup MOCI parameters — may not exist if this is an older model fit.
  # Use a separate fit$summary() call with tryCatch so that a missing parameter
  # returns 0 (scalar) rather than numeric(0), which would collapse the tibble.
  get_new_param <- function(param) {
    tryCatch({
      v <- fit$summary(variables = param)$mean
      if (length(v) == 0L) { warning(paste("Decomp:", param, "not in model — set to 0")); 0 }
      else v
    }, error = function(e) {
      warning(paste("Decomp:", param, "not in model — set to 0"))
      0
    })
  }
  b_ond_pup <- get_new_param("beta_moci_ond_pup")
  b_jfm_pup <- get_new_param("beta_moci_jfm_pup")
  
  coy_idx <- c(1, 2, 3, 0, 0, 0)   # BL=1, DE=2, DP=3
  has_es  <- c(0, 1, 0, 1, 0, 0)   # DE=2, PRH=4
  
  # ── Compute contributions for every site × year ───────────────────────────
  contrib <- map_dfr(1:S, function(si) {
    map_dfr(seq_along(years), function(t) {
      
      tb <- if (t > 1L) t - 1L else 1L   # birth-year index (t-1)
      
      # Pup survival covariates — ALL at t_birth (birth year t-1)
      coy  <- if (coy_idx[si] > 0) b_coy[coy_idx[si]] * sdat$coyote[si, tb]        else 0
      dst  <- b_dst[si]                               * sdat$disturbance[si, tb]
      es   <- if (has_es[si]  > 0) b_es               * sdat$elephant_seal[si, tb]  else 0
      amj  <- b_amj  * sdat$moci_amj[tb]   # birth year t-1
      # New: OND and JFM MOCI on pup survival
      ond_pup <- b_ond_pup * sdat$moci_ond[t]    # post-weaning fall (OND of year t-1)
      jfm_pup <- b_jfm_pup * sdat$moci_jfm[t]    # first winter
      
      # Fecundity covariate (OND at t = OND of year t-1, pre-lagged in data)
      ond  <- b_ond  * sdat$moci_ond[t]
      
      # Juvenile and adult survival (JFM at t)
      jfmJ <- b_jfmJ * sdat$moci_jfm[t]
      jfmA <- b_jfmA * sdat$moci_jfm[t]
      
      tibble(
        Site         = site_names[si],
        Year         = years[t],
        `Coyote (t-1)`           = coy,
        `Disturbance (t-1)`      = dst,
        `Elephant seal (t-1)`    = es,
        `MOCI AMJ → pup (t-1)`   = amj,
        `MOCI OND → pup (t)`     = ond_pup,
        `MOCI JFM → pup (t)`     = jfm_pup,
        `MOCI OND → fecund (t)`  = ond,
        `MOCI JFM → juv (t)`     = jfmJ,
        `MOCI JFM → adult (t)`   = jfmA
      )
    })
  }) |> mutate(Site = factor(Site, levels = site_names))
  
  # ── Colour palette ─────────────────────────────────────────────────────────
  cov_cols <- c(
    `Coyote (t-1)`          = "#B2182B",
    `Disturbance (t-1)`     = "#8C510A",
    `Elephant seal (t-1)`   = "#762A83",
    `MOCI AMJ → pup (t-1)` = "#4DAC26",
    `MOCI OND → pup (t)`   = "#B8E186",
    `MOCI JFM → pup (t)`   = "#74C476",
    `MOCI OND → fecund (t)` = "#1B7837",
    `MOCI JFM → juv (t)`    = "#2166AC",
    `MOCI JFM → adult (t)`  = "#74ADD1"
  )
  
  # Shared theme additions
  bar_theme <- list(
    geom_hline(yintercept = 0, linewidth = 0.6, colour = "gray30"),
    facet_wrap(~Site, ncol = 2, scales = "fixed"),
    scale_x_continuous(breaks = seq(min(years), max(years), by = 4)),
    labs(x = "Year", y = "Covariate effect (logit scale)", fill = NULL),
    theme_seal(),
    theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 9),
          legend.position = "bottom",
          strip.text      = element_text(size = 11, face = "bold"),
          panel.spacing   = unit(0.5, "lines"))
  )
  
  # ── (a) Pup survival decomposition ────────────────────────────────────────
  pup_vars <- c("Coyote (t-1)", "Disturbance (t-1)",
                "Elephant seal (t-1)", "MOCI AMJ → pup (t-1)",
                "MOCI OND → pup (t)", "MOCI JFM → pup (t)")
  
  pup_long <- contrib |>
    select(Site, Year, all_of(pup_vars)) |>
    pivot_longer(cols=-all_of(c("Site","Year")), names_to="Covariate", values_to="Effect") |>
    filter(!is.na(Effect) & Effect != 0) |>
    mutate(Covariate = factor(Covariate, levels = pup_vars))
  
  p_pup <- ggplot(pup_long, aes(x = Year, y = Effect, fill = Covariate)) +
    geom_col(position = "stack", width = 0.75, alpha = 0.85) +
    scale_fill_manual(values = cov_cols[pup_vars]) +
    labs(title   = "Covariate Contributions to Pup Survival by Site and Year",
         subtitle = paste0("Logit-scale effects: stacked positive = higher φ_pup, ",
                           "negative = lower φ_pup; birth-year covariates at t−1")) +
    bar_theme
  if (save) ggsave(paste0("Output/Plots/", prefix,
                          "_decomp_pup_survival.jpeg"),
                   p_pup, width = 32, height = 30, units = "cm", dpi = 200)
  
  # ── (b) Fecundity decomposition ────────────────────────────────────────────
  fec_long <- contrib |>
    select(Site, Year, `MOCI OND → fecund (t)`) |>
    pivot_longer(cols=-all_of(c("Site","Year")), names_to="Covariate", values_to="Effect") |>
    mutate(Covariate = factor(Covariate))
  
  p_fec <- ggplot(fec_long, aes(x = Year, y = Effect, fill = Covariate)) +
    geom_col(width = 0.75, alpha = 0.85) +
    scale_fill_manual(values = cov_cols["MOCI OND → fecund (t)"]) +
    labs(title   = "MOCI Fall Contribution to Fecundity by Site and Year",
         subtitle = paste0("Logit-scale effect: moci_ond(t) = OND of year t−1 ",
                           "(pre-lagged); negative = warm fall → lower prob. pupping")) +
    bar_theme
  if (save) ggsave(paste0("Output/Plots/", prefix,
                          "_decomp_fecundity.jpeg"),
                   p_fec, width = 32, height = 30, units = "cm", dpi = 200)
  
  # ── (c) All covariates combined ────────────────────────────────────────────
  all_vars <- names(cov_cols)
  
  all_long <- contrib |>
    pivot_longer(cols=-all_of(c("Site","Year")), names_to="Covariate", values_to="Effect") |>
    filter(!is.na(Effect) & Effect != 0) |>
    mutate(Covariate = factor(Covariate, levels = all_vars))
  
  p_all <- ggplot(all_long, aes(x = Year, y = Effect, fill = Covariate)) +
    geom_col(position = "stack", width = 0.75, alpha = 0.82) +
    scale_fill_manual(values = cov_cols) +
    labs(title   = "All Covariate Contributions by Site and Year",
         subtitle = paste0("Additive logit-scale effects across φ_pup, fecundity, ",
                           "φ_juv, φ_adult; positive bars = beneficial conditions")) +
    bar_theme +
    guides(fill = guide_legend(nrow = 3))
  if (save) ggsave(paste0("Output/Plots/", prefix,
                          "_decomp_all_covariates.jpeg"),
                   p_all, width = 32, height = 32, units = "cm", dpi = 200)
  
  # ── (d) Net effect line overlay (sum of all effects by site × year) ────────
  net_df <- all_long |>
    group_by(Site, Year) |>
    summarise(Net = sum(Effect), .groups = "drop")
  
  p_net <- ggplot() +
    geom_col(data = all_long,
             aes(x = Year, y = Effect, fill = Covariate),
             position = "stack", width = 0.75, alpha = 0.78) +
    geom_line(data = net_df,
              aes(x = Year, y = Net),
              colour = "black", linewidth = 0.9, linetype = "solid") +
    geom_point(data = net_df,
               aes(x = Year, y = Net),
               colour = "black", size = 1.8) +
    geom_hline(yintercept = 0, linewidth = 0.6, colour = "gray30") +
    facet_wrap(~Site, ncol = 2, scales = "fixed") +
    scale_fill_manual(values = cov_cols, name = NULL) +
    scale_x_continuous(breaks = seq(min(years), max(years), by = 4)) +
    labs(x = "Year", y = "Covariate effect (logit scale)",
         title   = "Covariate Contributions with Net Effect by Site and Year",
         subtitle = "Stacked bars = individual effects; black line = net sum of all covariates") +
    theme_seal() +
    theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 9),
          legend.position = "bottom",
          strip.text      = element_text(size = 11, face = "bold"),
          panel.spacing   = unit(0.5, "lines")) +
    guides(fill = guide_legend(nrow = 3))
  if (save) ggsave(paste0("Output/Plots/", prefix,
                          "_decomp_net_effect.jpeg"),
                   p_net, width = 32, height = 32, units = "cm", dpi = 200)
  
  list(pup_survival = p_pup,
       fecundity    = p_fec,
       all_combined = p_all,
       net_effect   = p_net,
       data         = contrib)
}

cat("============================================\n")
cat("  Credible intervals: 89% (CI_LO=0.055, CI_HI=0.945)\n")
cat("  Plot theme: theme_seal(base_size=16)\n")
cat("\nKEY CHANGES FROM v3.1:\n")
cat("  (1) Pup + juvenile survival: sex-neutral (phi_pup_base, phi_juv_base)\n")
cat("      Adult survival only: sex-specific (phi_adult_F_logit, delta_adult)\n")
cat("  (2) Pup survival prior: logit-normal(0, 0.8) → mean p ≈ 0.50\n")
cat("      logit-normal Normal(0,0.4) → mean p ≈ 0.50; SD tightened from 0.8\n")
cat("  (3) Breeding adult counts: N_adult_F + N_adult_M * p_male_breed\n")
cat("      Molt counts: N_juvF+N_juvM+N_adultF+N_adultM (50:50, no extra param)\n")
cat("\nTo run with REAL data:\n")
cat("  results.real <- run_full_analysis_v3.2(\n")
cat("    use_real_data=TRUE, dat=dat, cov_t_scaled=cov_t_scaled, years=years,\n")
cat("    iter_warmup=3000, iter_sampling=1000, adapt_delta=0.995\n")
cat("  )\n")
cat("\nSite indexing: 1=BL 2=DE 3=DP 4=PRH 5=TB 6=TP\n")
cat("============================================\n")
