# ============================================================================
# 19_regional_run.R
# ----------------------------------------------------------------------------
# REGIONAL HARBOR SEAL IPM — RUN SCRIPT
#
# This is the top-level execution script. Source this file to run the
# complete regional IPM pipeline. Modify the settings in Part 0 to switch
# between simulated and real data, adjust MCMC settings, etc.
#
# Parallel to 07_ipm_run.R in the Marin IPM workflow.
#
# Pipeline:
#   15_regional_data_prep.R   — data cleaning + Stan data list
#   16_regional_simulate.R    — simulate data for model testing
#   17_regional_ipm_model.R   — Stan model + run_regional_ipm() orchestrator
#   18_regional_plots.R       — all plotting functions
#   19_regional_run.R         — THIS FILE: calls everything in sequence
#
# Recommended workflow:
#   1. Run with USE_REAL_DATA = FALSE first (simulated data, ~8 min)
#      Verify delta_moci_bay posteriors contain true values
#   2. Run with USE_REAL_DATA = TRUE for publication results (~9 min)
# ============================================================================

# ── Part 0: SETTINGS ──────────────────────────────────────────────────────────
USE_REAL_DATA  <- TRUE     # FALSE = simulated data for testing; TRUE = real data
RUN_SIMTEST    <- FALSE    # Set TRUE to run simulation test before real data
# Adds ~8 min but confirms model identifiability

ITER_WARMUP    <- 2000   # mass matrix adaptation for ~960 eps_adult_raw parameters
ITER_SAMPLING  <- 1000   # sufficient: worst-case ESS ~1000+ (detect_molt_logit)
ADAPT_DELTA    <- 0.95
MAX_TREEDEPTH  <- 12
SEED           <- 456
T_PROJ         <- 5    # 5-year projections: defensible near-term scenario analysis

# ── Output file naming ────────────────────────────────────────────────────────
# These prefixes are the single source of truth for ALL output file names
# across scripts 17, 18, and 19.  Every fit RDS, input RDS, CSV, and plot
# JPEG produced by the regional IPM pipeline will begin with one of these
# strings, clearly distinguishing regional outputs from the 6-site Marin IPM
# (which uses prefixes "IPM_v3.3_real", "IPM_v3.3_sim", etc.).
#
# File inventory produced by this script:
#   Output/harbor_seal_Regional_real_fit.rds
#   Output/harbor_seal_Regional_real_input_data.rds
#   Output/regional_ipm_input_data.rds          (from 15_regional_data_prep.R)
#   Output/Regional_real_parameter_summary.csv
#   Output/Regional_real_portfolio_summary.csv
#   Output/Plots/Regional_real_*.jpeg           (all plot files)
#   Output/Regional_sim_*                       (if RUN_SIMTEST = TRUE)

PREFIX_SIM  <- "Regional_sim"
PREFIX_REAL <- "Regional_real"

# ── Part 1: SOURCE ALL SCRIPTS ────────────────────────────────────────────────
cat("=================================================================\n")
cat("  REGIONAL HARBOR SEAL IPM\n")
cat("  Central and Northern California, 2005-2025\n")
cat("=================================================================\n\n")

source("Code/15_regional_data_prep.R")   # creates regional_stan, regional_sites, etc.
source("Code/17_regional_ipm_model.R")   # creates run_regional_ipm(), writes Stan file
source("Code/18_regional_plots.R")        # creates run_all_regional_plots(), etc.

# ── Part 2: OPTIONAL SIMULATION TEST ─────────────────────────────────────────
if (RUN_SIMTEST) {
  cat("\n── Running simulation test ──────────────────────────────────────────\n")
  cat("Goal: confirm delta_moci_bay parameters are recoverable\n\n")
  
  source("Code/16_regional_simulate.R")  # creates sim_regional
  
  out_sim <- run_regional_ipm(
    use_real_data = FALSE,
    sim_data      = sim_regional,
    prefix        = PREFIX_SIM,
    seed          = SEED,
    iter_warmup   = 1000,    # shorter for test run
    iter_sampling = 1000,
    adapt_delta   = ADAPT_DELTA,
    max_treedepth = MAX_TREEDEPTH
  )
  
  # ── Simulation recovery — 89% CrI method matching Marin IPM ────────────────
  # Uses draw-level quantiles (CI_LO=0.055, CI_HI=0.945) not CmdStan defaults,
  # and covers all key parameters including the new bay/coast modifiers.
  cat("\n── Simulation recovery check (89% CrI, Marin IPM method) ─────────\n")
  tp <- sim_regional$true_params
  
  true_vals_tbl <- tibble::tribble(
    ~parameter,                ~true_value,
    "phi_pup_logit",            tp$phi_pup_logit,
    "phi_juv_base",             tp$phi_juv_base,
    "phi_adult_F_logit",        tp$phi_adult_F_logit,
    "delta_adult",              tp$delta_adult,
    "avg_fecundity",            0.20*tp$fecund_primip + 0.80*tp$fecund_mature,
    "rho_pup",                  tp$rho_pup,
    "rho_juv_molt",             tp$rho_juv_molt,
    "beta_moci_ond_fecund",     tp$beta_moci_ond_fecund,
    "beta_moci_jfm_adult",      tp$beta_moci_jfm_adult,
    "delta_moci_mouth_fecund",  tp$delta_moci_mouth_fecund,
    "delta_moci_mouth_surv",    tp$delta_moci_mouth_surv,
    "delta_moci_south_fecund",  tp$delta_moci_south_fecund,
    "delta_moci_south_surv",    tp$delta_moci_south_surv,
    "delta_moci_marin_fecund",  tp$delta_moci_marin_fecund,
    "sigma_county",             tp$sigma_county,
    "sigma_site",               tp$sigma_site,
    "detect_breed_logit",       tp$detect_breed_logit,
    "detect_molt_logit",        tp$detect_molt_logit,
    "sigma_process",            tp$sigma_process,
    "sigma_obs_adult",          tp$sigma_obs_adult,
    "sigma_obs_molt",           tp$sigma_obs_molt
    # sigma_obs_pup removed: marginalised analytically (fixed at 0.148)
    # sigma_pup_eff = sqrt(sigma_process^2 + 0.148^2) is the effective pup SD
  )
  
  rec <- tryCatch({
    d <- out_sim$fit$draws(variables = true_vals_tbl$parameter, format = "matrix")
    s <- out_sim$fit$summary(variables = true_vals_tbl$parameter)
    s$q_lo <- apply(d[, s$variable, drop = FALSE], 2, quantile, 0.055)
    s$q_hi <- apply(d[, s$variable, drop = FALSE], 2, quantile, 0.945)
    s |>
      left_join(true_vals_tbl, by = c("variable" = "parameter")) |>
      mutate(
        recovered    = true_value >= q_lo & true_value <= q_hi,
        rel_bias_pct = (mean - true_value) / abs(true_value) * 100,
        identifiability = case_when(
          variable %in% c("phi_adult_F_logit","avg_fecundity",
                          "beta_moci_ond_fecund",
                          "sigma_obs_adult")                ~ "Well identified",
          variable %in% c("phi_pup_logit","phi_juv_base",
                          "beta_moci_jfm_adult",
                          "delta_moci_mouth_fecund",
                          "delta_moci_mouth_surv",
                          "delta_moci_south_fecund",
                          "detect_breed_logit",
                          "sigma_county","sigma_process")   ~ "Moderately identified",
          TRUE                                              ~ "Prior dominated")
      )
  }, error = function(e) { cat("Recovery check failed:", conditionMessage(e), "\n"); NULL })
  
  if (!is.null(rec)) {
    cat("\n")
    print(rec |> select(variable, true_value, mean, q_lo, q_hi,
                        recovered, rel_bias_pct, identifiability),
          n = nrow(rec))
    cat(sprintf("\nOverall recovery: %d/%d (%.0f%%)  mean |bias|: %.1f%%\n",
                sum(rec$recovered), nrow(rec),
                100 * mean(rec$recovered),
                mean(abs(rec$rel_bias_pct), na.rm = TRUE)))
    
    # Flag all four modifier parameters explicitly
    for (p in c("delta_moci_mouth_fecund","delta_moci_mouth_surv",
                "delta_moci_south_fecund","delta_moci_south_surv")) {
      r <- rec[rec$variable == p, ]
      if (nrow(r) > 0)
        cat(sprintf("  %-28s true=%6.3f  post=%6.3f  89%%CrI=[%6.3f,%6.3f]  %s\n",
                    p, r$true_value, r$mean, r$q_lo, r$q_hi,
                    ifelse(r$recovered, "RECOVERED", "MISSED")))
    }
    cat("  Expected: south > mouth > 0 gradient (more negative = coast-like)\n")
    
    if (mean(rec$recovered) < 0.70)
      warning("Recovery < 70% — check identifiability before fitting real data.")
  }
  
  # Quick simulation plots
  sim_plots <- run_all_regional_plots(
    fit        = out_sim$fit,
    model_data = out_sim$model_data,
    prefix     = PREFIX_SIM,
    save       = TRUE
  )
  cat("\nSimulation test complete. Check Output/Plots/ for sim diagnostics.\n")
  cat("Proceeding to real data run...\n\n")
}

# ── Part 3: FIT REAL DATA ─────────────────────────────────────────────────────
if (USE_REAL_DATA) {
  cat("\n── Fitting real data ────────────────────────────────────────────────\n")
  
  out_real <- run_regional_ipm(
    use_real_data  = TRUE,
    input_rds      = "Output/regional_ipm_input_data.rds",
    prefix         = PREFIX_REAL,
    T_proj         = T_PROJ,
    seed           = SEED,
    iter_warmup    = ITER_WARMUP,
    iter_sampling  = ITER_SAMPLING,
    adapt_delta    = ADAPT_DELTA,
    max_treedepth  = MAX_TREEDEPTH
  )
  
  # ── Part 4: GENERATE ALL PLOTS ──────────────────────────────────────────────
  cat("\n── Generating plots ─────────────────────────────────────────────────\n")
  
  real_plots <- run_all_regional_plots(
    fit        = out_real$fit,
    model_data = out_real$model_data,
    prefix     = PREFIX_REAL,
    save       = TRUE
  )
  
  # ── Part 5: PRINT KEY RESULTS SUMMARY ────────────────────────────────────────
  cat("\n=================================================================\n")
  cat("  KEY RESULTS SUMMARY — REGIONAL IPM\n")
  cat("=================================================================\n")
  
  fit <- out_real$fit
  
  # Vital rates
  vr_params <- c("phi_pup_logit", "phi_juv_base", "phi_adult_F_base",
                 "delta_adult", "avg_fecundity", "rho_pup")
  s_vr <- tryCatch(fit$summary(variables = vr_params), error = function(e) NULL)
  if (!is.null(s_vr)) {
    cat("\nVital rates:\n")
    pup_l <- fit$draws("phi_pup_logit", format = "df")$phi_pup_logit
    cat(sprintf("  phi_pup (prob): %.3f (89%% CrI: %.3f-%.3f)\n",
                median(plogis(pup_l)), quantile(plogis(pup_l), 0.055),
                quantile(plogis(pup_l), 0.945)))
    for (v in vr_params[-1]) {
      r <- s_vr[s_vr$variable == v, ]
      if (nrow(r) > 0)
        cat(sprintf("  %s: %.3f (q5=%.3f, q95=%.3f)\n",
                    v, r$mean, r$q5, r$q95))
    }
  }
  
  # Bay vs coast modifier — the key new result
  bay_params <- c("delta_moci_bay_fecund", "delta_moci_bay_surv")
  s_bay <- tryCatch(fit$summary(variables = bay_params), error = function(e) NULL)
  if (!is.null(s_bay)) {
    cat("\nBay vs. Coast MOCI modifier (key new parameters):\n")
    for (v in bay_params) {
      r   <- s_bay[s_bay$variable == v, ]
      d   <- as.numeric(fit$draws(v, format = "matrix"))
      p_pos <- mean(d > 0)
      if (nrow(r) > 0)
        cat(sprintf("  %s: %.3f (q5=%.3f, q95=%.3f) P(>0)=%.3f\n",
                    v, r$mean, r$q5, r$q95, p_pos))
    }
    cat("  Interpretation: positive values = SF Estuary less MOCI-sensitive than coast\n")
  }
  
  # County lambda summary
  T <- length(out_real$model_data$years)
  C <- out_real$model_data$stan_data$C %||% 6
  county_names <- out_real$model_data$county_names %||%
    c("Marin","Bay Mouth","South Bay","San Mateo","Sonoma","Mendocino")
  cat("\nMean annual lambda by county (2005-2024):\n")
  for (c in 1:C) {
    lam_vars <- paste0("lambda_county[", c, ",", 1:(T-1), "]")
    lam_d <- tryCatch(fit$draws(variables = lam_vars, format = "matrix"),
                      error = function(e) NULL)
    if (!is.null(lam_d)) {
      lam_mean <- mean(colMeans(lam_d))
      lam_lo   <- mean(apply(lam_d, 2, quantile, 0.055))
      lam_hi   <- mean(apply(lam_d, 2, quantile, 0.945))
      cat(sprintf("  %-14s: lambda_bar=%.3f (89%%CrI: %.3f-%.3f)  %s\n",
                  county_names[c], lam_mean, lam_lo, lam_hi,
                  if (lam_hi < 1) "DECLINING" else if (lam_lo > 1) "GROWING" else "~STABLE"))
    }
  }
  
  # Portfolio
  port <- real_plots$portfolio
  if (!is.null(port$portfolio_summary)) {
    cat("\nPortfolio analysis:\n")
    print(port$portfolio_summary)
  }
  
  cat("\n=================================================================\n")
  cat(sprintf("  COMPLETE — Regional IPM (%s)\n", PREFIX_REAL))
  cat("-----------------------------------------------------------------\n")
  cat(sprintf("  Fit RDS     : Output/harbor_seal_%s_fit.rds\n",       PREFIX_REAL))
  cat(sprintf("  Input RDS   : Output/harbor_seal_%s_input_data.rds\n", PREFIX_REAL))
  cat(sprintf("  Param CSV   : Output/%s_parameter_summary.csv\n",      PREFIX_REAL))
  cat(sprintf("  Portfolio   : Output/%s_portfolio_summary.csv\n",       PREFIX_REAL))
  cat(sprintf("  Plots       : Output/Plots/%s_*.jpeg\n",                PREFIX_REAL))
  cat("=================================================================\n\n")
  
} else {
  cat("\nUSE_REAL_DATA = FALSE. To fit real data, set USE_REAL_DATA <- TRUE\n")
  cat("and re-source this file.\n")
  if (!RUN_SIMTEST)
    cat("To run simulation test, set RUN_SIMTEST <- TRUE.\n")
}

# ============================================================================
# P_MALE SENSITIVITY ANALYSIS
# ----------------------------------------------------------------------------
# Mirrors Marin IPM Table S4.1. Three sensitivity runs fix p_male_fixed at
# 0.05, 0.10, and 0.20 and compare key posteriors against the base run
# (p_male_fixed = 0.057). Lambda estimates are expected to be robust;
# fecundity and detect_molt are expected to show systematic gradients.
#
# Can be run independently of USE_REAL_DATA — the base fit is reloaded
# from disk if out_real is not already in the R session:
#
#   RUN_SENSITIVITY <- TRUE
#   source("Code/19_regional_run.R")
#
# Set RUN_SENSITIVITY <- TRUE to execute (adds ~3 x chain time).
# ============================================================================

RUN_SENSITIVITY <- TRUE   # set TRUE to run; FALSE to skip

if (RUN_SENSITIVITY) {
  
  # ── Reload base fit from disk if not in session ──────────────────────────
  if (!exists("out_real") || is.null(out_real$fit)) {
    cat("Loading saved base fit from disk...\n")
    fit_path <- sprintf("Output/harbor_seal_%s_fit.rds", PREFIX_REAL)
    inp_path <- sprintf("Output/harbor_seal_%s_input_data.rds", PREFIX_REAL)
    if (!file.exists(fit_path))
      stop("Base fit not found at: ", fit_path,
           "\nRun USE_REAL_DATA <- TRUE first.")
    base_fit_obj  <- readRDS(fit_path)
    base_inp      <- readRDS(inp_path)
    # Reconstruct minimal out_real structure
    out_real <- list(
      fit        = base_fit_obj,
      model_data = list(
        stan_data    = base_inp$stan_data,
        years        = base_inp$years,
        county_names = c("Marin","Bay Mouth","South Bay",
                         "San Mateo","Sonoma","Mendocino")
      )
    )
    cat("Base fit loaded.\n")
  }
  
  cat("\n================================================================\n")
  cat("   P_MALE SENSITIVITY ANALYSIS (regional IPM)\n")
  cat("   Three fixed values: 0.05, 0.10, 0.20\n")
  cat("   Base model: p_male_fixed = 0.057\n")
  cat("================================================================\n\n")
  
  # ── Parameters to track ────────────────────────────────────────────────────
  SENS_PARAMS <- c(
    "phi_pup_logit", "phi_juv_base", "phi_adult_F_logit", "delta_adult",
    "avg_fecundity",
    "beta_moci_ond_fecund", "beta_moci_jfm_adult",
    "delta_moci_mouth_fecund", "delta_moci_mouth_surv",
    "delta_moci_marin_fecund",
    "detect_breed_logit", "detect_molt_logit",
    "sigma_process", "sigma_obs_adult", "sigma_obs_molt",
    "sigma_site", "sigma_county"
  )
  
  P_MALE_VALUES <- c(base = 0.057, s1 = 0.05, s2 = 0.10, s3 = 0.20)
  SENS_LABELS   <- c("Base (0.057)", "Fixed 0.05", "Fixed 0.10", "Fixed 0.20")
  
  # ── Helper: extract posterior summary for one fit ─────────────────────────
  extract_sens_summary <- function(fit, p_male_val, county_names) {
    s <- tryCatch(fit$summary(variables = SENS_PARAMS), error = function(e) NULL)
    if (is.null(s)) return(NULL)
    
    # Also extract lambda by county
    lam_rows <- map_dfr(seq_along(county_names), function(ci) {
      lam_vars <- paste0("lambda_county[", ci, ",", 1:20, "]")
      lam_d    <- tryCatch(fit$draws(variables = lam_vars, format = "matrix"),
                           error = function(e) NULL)
      if (is.null(lam_d)) return(NULL)
      tibble(variable = paste0("lambda_bar_", county_names[ci]),
             mean     = mean(colMeans(lam_d)),
             q5       = mean(apply(lam_d, 2, quantile, 0.055)),
             q95      = mean(apply(lam_d, 2, quantile, 0.945)))
    })
    
    bind_rows(
      s |> select(variable, mean, q5, q95),
      lam_rows
    ) |> mutate(p_male = p_male_val)
  }
  
  # ── Extract base run posteriors ───────────────────────────────────────────
  base_fit     <- out_real$fit
  county_names <- out_real$model_data$county_names %||%
    c("Marin","Bay Mouth","South Bay","San Mateo","Sonoma","Mendocino")
  inp_rds      <- sprintf("Output/harbor_seal_%s_input_data.rds", PREFIX_REAL)
  
  sens_results <- list()
  sens_results[["base"]] <- extract_sens_summary(base_fit, 0.057, county_names)
  for (nm in c("s1","s2","s3")) {
    p_val <- P_MALE_VALUES[[nm]]
    cat(sprintf("\n--- Sensitivity run: p_male_fixed = %.2f ---\n", p_val))
    
    sens_dat <- prepare_real_data_regional(inp_rds)
    sens_dat$p_male_fixed <- p_val
    
    sens_fit_obj <- tryCatch(
      run_regional_ipm(
        use_real_data  = TRUE,
        stan_data_list = sens_dat,
        prefix         = sprintf("Regional_pmale_%.2f", p_val),
        iter_warmup    = ITER_WARMUP,
        iter_sampling  = ITER_SAMPLING,
        adapt_delta    = ADAPT_DELTA,
        save_fit       = FALSE),
      error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL })
    
    if (!is.null(sens_fit_obj))
      sens_results[[nm]] <- extract_sens_summary(sens_fit_obj$fit, p_val, county_names)
  }
  
  # ── Compile sensitivity table ─────────────────────────────────────────────
  sens_df <- bind_rows(sens_results) |>
    mutate(
      label    = SENS_LABELS[match(p_male, P_MALE_VALUES)],
      Estimate = sprintf("%.3f (%.3f, %.3f)", mean, q5, q95)
    )
  
  # Pivot wide for comparison
  sens_wide <- sens_df |>
    select(variable, label, mean, q5, q95) |>
    pivot_wider(names_from = label,
                values_from = c(mean, q5, q95),
                names_glue = "{label}_{.value}")
  
  # ── Sensitivity flags ─────────────────────────────────────────────────────
  flag_sensitivity <- function(means, lo_base, hi_base, lo_high, hi_high) {
    # Check monotonic trend across p_male values
    is_monotone <- all(diff(means) > 0) || all(diff(means) < 0)
    range_span  <- max(means) - min(means)
    # Check if CrIs of extreme runs overlap
    overlap <- (lo_high < hi_base) && (lo_base < hi_high)
    if (!overlap)                             return("\u25b2\u25b2")   # ▲▲
    if (is_monotone && range_span > 0.02)     return("\u25b2")        # ▲
    if (is_monotone && range_span > 0.005)    return("~")
    return("")
  }
  
  run_labels <- c("Base (0.057)", "Fixed 0.05", "Fixed 0.10", "Fixed 0.20")
  sens_flags <- sens_wide |>
    rowwise() |>
    mutate(
      Flag = flag_sensitivity(
        means   = c(`Base (0.057)_mean`, `Fixed 0.05_mean`,
                    `Fixed 0.10_mean`,   `Fixed 0.20_mean`),
        lo_base = `Base (0.057)_q5`,   hi_base  = `Base (0.057)_q95`,
        lo_high = `Fixed 0.20_q5`,     hi_high  = `Fixed 0.20_q95`
      )
    ) |> ungroup()
  
  # ── Print summary table ───────────────────────────────────────────────────
  cat("\n=== P_MALE SENSITIVITY TABLE — REGIONAL IPM ===\n")
  cat(sprintf("%-32s  %-22s  %-22s  %-22s  %-22s  %s\n",
              "Parameter", "Base (p=0.057)", "Fixed p=0.05",
              "Fixed p=0.10", "Fixed p=0.20", "Flag"))
  cat(strrep("-", 130), "\n")
  
  for (i in seq_len(nrow(sens_flags))) {
    row <- sens_flags[i, ]
    fmt <- function(lbl) {
      m  <- row[[paste0(lbl, "_mean")]]
      lo <- row[[paste0(lbl, "_q5")]]
      hi <- row[[paste0(lbl, "_q95")]]
      if (any(is.null(c(m, lo, hi)) | is.na(c(m, lo, hi)))) return("—")
      sprintf("%.3f (%.3f, %.3f)", m, lo, hi)
    }
    cat(sprintf("%-32s  %-22s  %-22s  %-22s  %-22s  %s\n",
                row$variable,
                fmt("Base (0.057)"), fmt("Fixed 0.05"),
                fmt("Fixed 0.10"),   fmt("Fixed 0.20"),
                row$Flag))
  }
  
  # ── Save CSV ──────────────────────────────────────────────────────────────
  write_csv(sens_flags, "Output/Regional_pmale_sensitivity.csv")
  cat("\nSensitivity table saved: Output/Regional_pmale_sensitivity.csv\n")
  cat("\nKey parameters to flag:\n")
  flagged <- sens_flags |> filter(Flag != "")
  for (i in seq_len(nrow(flagged))) {
    cat(sprintf("  %-30s  %s\n", flagged$variable[i], flagged$Flag[i]))
  }
}

