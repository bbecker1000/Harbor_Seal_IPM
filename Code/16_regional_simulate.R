# ============================================================================
# 16_regional_simulate.R
# ----------------------------------------------------------------------------
# REGIONAL IPM — SIMULATE DATA FOR MODEL TESTING
#
# Generates a synthetic regional dataset under known parameters to:
#   (1) Verify that county-level population dynamics + site availability
#       parameters are identifiable from the combined observation streams
#   (2) Confirm that the Bay vs. Open Coast MOCI modifier (delta_moci_bay)
#       is recoverable from data — the key new parameter in this model
#   (3) Test that the 2020 latent-year handling propagates state uncertainty
#       correctly without biasing vital rate estimates
#   (4) Confirm rho_pup (pup molt attendance) remains orthogonal to vital rates
#
# The simulation mirrors the 6-county structure in 15_regional_data_prep.R.
# Key structural differences from the Marin IPM simulation:
#   - County-level (not site-level) Leslie matrix dynamics
#   - Site availability parameters (alpha_breed, alpha_pup, alpha_molt)
#   - Bay vs. Coast MOCI modifier (delta_moci_bay_fecund, delta_moci_bay_surv)
#   - No coyote, disturbance, or elephant seal covariates
#   - p_male_breed fixed (not estimated) at 0.057
#
# Prereq: source("Code/15_regional_data_prep.R")
# Next:   source("Code/17_regional_ipm_model.R")
# ============================================================================

library(tidyverse)

# ── Pull site metadata from data prep (or reload) ─────────────────────────────
if (!exists("regional_sites")) {
  inp <- readRDS("Output/regional_ipm_input_data.rds")
  regional_sites <- inp$site_meta
  years_regional <- inp$years
}

# ── SIMULATION FUNCTION ───────────────────────────────────────────────────────
simulate_regional_ipm <- function(T    = 21,
                                  S    = 24,
                                  C    = 6,
                                  T_proj = 10,
                                  seed = 456,
                                  site_meta = regional_sites,
                                  true_params = NULL) {
  set.seed(seed)
  
  if (is.null(true_params)) true_params <- list(
    # ── Vital rates ───────────────────────────────────────────────────────────
    phi_pup_logit     = qlogis(0.22),
    phi_juv_base      = 0.70,
    phi_adult_F_logit = qlogis(0.93),
    delta_adult       = 0.017,
    fecund_primip     = 0.60,
    fecund_mature     = 0.85,
    prop_female       = 0.50,
    rho_pup           = 0.18,
    p_male_fixed      = 0.057,
    
    # ── MOCI (open coast baseline) ────────────────────────────────────────────
    beta_moci_ond_fecund = -0.15,
    beta_moci_ond_pup    = -0.15,
    beta_moci_amj_pup    = -0.15,
    beta_moci_jfm_pup    = -0.15,
    beta_moci_jfm_juv    = -0.12,
    beta_moci_jfm_adult  = -0.10,
    beta_moci_amj_molt   =  0.05,
    
    # ── Bay modifiers ─────────────────────────────────────────────────────────
    delta_moci_mouth_fecund = 0.08,
    delta_moci_mouth_surv   = 0.06,
    delta_moci_south_fecund = 0.14,
    delta_moci_south_surv   = 0.10,
    
    # ── County random effects ─────────────────────────────────────────────────
    sigma_county  = 0.15,
    county_effect = NULL,
    
    # ── Site detection (replaces county alpha in site-level model) ─────────────
    detect_breed_logit = 1.20,   # inv_logit(1.20) ≈ 0.77
    detect_molt_logit  = 0.75,   # inv_logit(0.75) ≈ 0.68
    sigma_site         = 0.20,   # site-to-site detection SD
    site_detect        = NULL,   # drawn below
    
    # ── Error structure ───────────────────────────────────────────────────────
    sigma_process   = 0.15,
    sigma_obs_adult = 0.18,
    sigma_obs_pup   = 0.15,
    sigma_obs_molt  = 0.30,
    
    # ── Site initial populations (per-site, site-level model) ─────────────────
    mu_log_adult = 5.5,   # exp(5.5) ≈ 245 per site
    mu_log_juv   = 4.5,
    mu_log_pup   = 4.5,
    sigma_init   = 0.30
  )
  
  # ── Draw random effects if not supplied ──────────────────────────────────────
  if (is.null(true_params$county_effect))
    true_params$county_effect <- rnorm(C, 0, true_params$sigma_county)
  if (is.null(true_params$site_detect))
    true_params$site_detect <- rnorm(nrow(site_meta), 0, true_params$sigma_site)
  
  # ── Derived scalars ───────────────────────────────────────────────────────────
  phi_adult_M_base <- plogis(true_params$phi_adult_F_logit) - true_params$delta_adult
  avg_fecundity    <- 0.20 * true_params$fecund_primip + 0.80 * true_params$fecund_mature
  county_type      <- c(0L, 1L, 2L, 0L, 0L, 0L)  # coast/mouth/south/coast/coast/coast
  years_sim        <- 2005:2025
  
  # ── Site start indices ────────────────────────────────────────────────────────
  site_t1 <- vapply(seq_len(nrow(site_meta)), function(si) {
    which(years_sim == site_meta$starts[si])
  }, integer(1))
  
  # ── MOCI time series ──────────────────────────────────────────────────────────
  T_mid <- 10
  moci_base <- c(seq(-0.5, 0.5, length.out = T_mid),
                 seq( 0.5,-0.5, length.out = T - T_mid)) +
    as.vector(arima.sim(list(ar = 0.5), n = T)) * 0.5
  moci_jfm <- as.vector(scale(moci_base))
  moci_amj <- as.vector(scale(moci_base * 0.80 +
                                as.vector(arima.sim(list(ar = 0.3), n = T)) * 0.4))
  moci_ond <- as.vector(scale(moci_base * 0.70 +
                                as.vector(arima.sim(list(ar = 0.3), n = T)) * 0.5))
  
  # ── Site initial populations (site-level) ─────────────────────────────────────
  N_adult_F_init <- exp(rnorm(nrow(site_meta), true_params$mu_log_adult, true_params$sigma_init))
  N_adult_M_init <- N_adult_F_init * 0.9
  N_juv_F_init   <- exp(rnorm(nrow(site_meta), true_params$mu_log_juv,   true_params$sigma_init)) * 0.5
  N_juv_M_init   <- N_juv_F_init
  N_pup_init     <- exp(rnorm(nrow(site_meta), true_params$mu_log_pup,   true_params$sigma_init))
  
  # ── Site-level state arrays ───────────────────────────────────────────────────
  S_model <- nrow(site_meta)
  N_adult_F  <- N_adult_M  <- matrix(NA, S_model, T)
  N_juv_F    <- N_juv_M    <- matrix(NA, S_model, T)
  N_pup      <- matrix(NA, S_model, T)
  
  # ── Simulate Leslie matrix per site ──────────────────────────────────────────
  for (s in seq_len(S_model)) {
    ci      <- site_meta$county_id[s]
    bay_fec <- (county_type[ci] == 1) * true_params$delta_moci_mouth_fecund +
      (county_type[ci] == 2) * true_params$delta_moci_south_fecund
    bay_sur <- (county_type[ci] == 1) * true_params$delta_moci_mouth_surv +
      (county_type[ci] == 2) * true_params$delta_moci_south_surv
    
    for (t in seq_len(T)) {
      if (t <= site_t1[s]) {
        # Hold at initial values before site's first survey year
        N_adult_F[s,t] <- N_adult_F_init[s]
        N_adult_M[s,t] <- N_adult_M_init[s]
        N_juv_F[s,t]   <- N_juv_F_init[s]
        N_juv_M[s,t]   <- N_juv_M_init[s]
        N_pup[s,t]     <- N_pup_init[s]
      } else {
        t_birth <- t - 1
        phi_pup_val <- plogis(
          true_params$phi_pup_logit + true_params$county_effect[ci] +
            (true_params$beta_moci_amj_pup + bay_sur) * moci_amj[t_birth] +
            (true_params$beta_moci_ond_pup + bay_sur) * moci_ond[t] +
            (true_params$beta_moci_jfm_pup + bay_sur) * moci_jfm[t])
        phi_juv_val <- plogis(
          qlogis(true_params$phi_juv_base) + true_params$county_effect[ci] * 0.5 +
            (true_params$beta_moci_jfm_juv + bay_sur) * moci_jfm[t])
        phi_aF_val <- plogis(
          true_params$phi_adult_F_logit + true_params$county_effect[ci] * 0.25 +
            (true_params$beta_moci_jfm_adult + bay_sur) * moci_jfm[t])
        phi_aM_val <- plogis(
          qlogis(phi_adult_M_base) + true_params$county_effect[ci] * 0.25 +
            (true_params$beta_moci_jfm_adult + bay_sur) * moci_jfm[t])
        fecund_t <- plogis(qlogis(avg_fecundity) +
                             (true_params$beta_moci_ond_fecund + bay_fec) * moci_ond[t])
        
        ep  <- N_adult_F[s,t-1] * fecund_t
        njF <- N_pup[s,t-1] * true_params$prop_female     * phi_pup_val
        njM <- N_pup[s,t-1] * (1-true_params$prop_female) * phi_pup_val
        jsF <- N_juv_F[s,t-1] * phi_juv_val * (2/3)
        jsM <- N_juv_M[s,t-1] * phi_juv_val * (2/3)
        jaF <- N_juv_F[s,t-1] * phi_juv_val * (1/3)
        jaM <- N_juv_M[s,t-1] * phi_juv_val * (1/3)
        
        N_pup[s,t]     <- exp(rnorm(1, log(max(ep,           1)), true_params$sigma_process))
        N_juv_F[s,t]   <- exp(rnorm(1, log(max(njF+jsF, 0.1)), true_params$sigma_process*0.5))
        N_juv_M[s,t]   <- exp(rnorm(1, log(max(njM+jsM, 0.1)), true_params$sigma_process*0.5))
        N_adult_F[s,t] <- exp(rnorm(1, log(max(N_adult_F[s,t-1]*phi_aF_val+jaF, 1)),
                                    true_params$sigma_process*0.5))
        N_adult_M[s,t] <- exp(rnorm(1, log(max(N_adult_M[s,t-1]*phi_aM_val+jaM, 1)),
                                    true_params$sigma_process*0.5))
      }
    }
  }
  
  N_adult_total <- N_adult_F + N_adult_M
  N_juv_total   <- N_juv_F   + N_juv_M
  N_molt_true   <- N_juv_total + N_adult_total + true_params$rho_pup * N_pup
  
  # ── Generate site-level observations ─────────────────────────────────────────
  y_adult <- y_pup <- y_molt <- matrix(NA, S_model, T)
  
  for (s in seq_len(S_model)) {
    ci <- site_meta$county_id[s]
    
    detect_breed <- plogis(true_params$detect_breed_logit + true_params$site_detect[s])
    bay_sur <- (county_type[ci] == 1) * true_params$delta_moci_mouth_surv +
      (county_type[ci] == 2) * true_params$delta_moci_south_surv
    
    for (t in seq_len(T)) {
      if (years_sim[t] == 2020) next          # latent year
      if (years_sim[t] < site_meta$starts[s]) next  # before site open
      
      detect_molt <- plogis(true_params$detect_molt_logit +
                              true_params$beta_moci_amj_molt * moci_amj[t] +
                              true_params$site_detect[s])
      
      N_ao <- N_adult_F[s,t] + N_adult_M[s,t] * true_params$p_male_fixed
      if (N_ao > 0)
        y_adult[s,t] <- rnorm(1, log(N_ao * detect_breed), true_params$sigma_obs_adult)
      
      if (N_pup[s,t] > 0)
        y_pup[s,t] <- rnorm(1, log(N_pup[s,t] * detect_breed), true_params$sigma_obs_pup)
      
      if (N_molt_true[s,t] > 0)
        y_molt[s,t] <- rnorm(1, log(N_molt_true[s,t] * detect_molt), true_params$sigma_obs_molt)
    }
  }
  
  # Random 5% additional missingness
  set_na <- function(mat) {
    idx <- which(!is.na(mat))
    mat[sample(idx, round(0.05 * length(idx)))] <- NA
    mat
  }
  y_adult <- set_na(y_adult)
  y_pup   <- set_na(y_pup)
  y_molt  <- set_na(y_molt)
  
  y_adult_obs <- ifelse(!is.na(y_adult), 1L, 0L)
  y_pup_obs   <- ifelse(!is.na(y_pup),   1L, 0L)
  y_molt_obs  <- ifelse(!is.na(y_molt),  1L, 0L)
  y_adult[is.na(y_adult)] <- 0
  y_pup[is.na(y_pup)]     <- 0
  y_molt[is.na(y_molt)]   <- 0
  
  # ── Stan data ─────────────────────────────────────────────────────────────────
  N_scen    <- 3
  moci_proj <- matrix(c(0, 1, -1), nrow = N_scen, ncol = T_proj)
  
  stan_data <- list(
    T = T, S = S_model, C = C, T_proj = T_proj, N_scenarios = N_scen,
    y_adult = y_adult, y_pup = y_pup, y_molt = y_molt,
    y_adult_obs = y_adult_obs, y_pup_obs = y_pup_obs, y_molt_obs = y_molt_obs,
    county_id   = site_meta$county_id,
    site_t1     = site_t1,
    county_type = county_type,
    moci_jfm = moci_jfm, moci_amj = moci_amj, moci_ond = moci_ond,
    moci_proj    = moci_proj,
    p_male_fixed = 0.057
  )
  
  true_params$county_type  <- county_type
  true_params$site_t1      <- site_t1
  
  list(
    stan_data   = stan_data,
    true_params = true_params,
    true_states = list(
      N_adult_F = N_adult_F, N_adult_M = N_adult_M,
      N_juv_F = N_juv_F, N_juv_M = N_juv_M, N_pup = N_pup,
      N_adult_total = N_adult_total, N_juv_total = N_juv_total,
      N_molt_true = N_molt_true
    ),
    county_names   = c("Marin","Bay Mouth","South Bay","San Mateo","Sonoma","Mendocino"),
    site_meta      = site_meta,
    years          = years_regional,
    scenario_names = c("Status Quo","Warm (MOCI +1)","Cool (MOCI -1)"),
    moci           = list(jfm = moci_jfm, amj = moci_amj, ond = moci_ond)
  )
}

# ── RUN SIMULATION ────────────────────────────────────────────────────────────
cat("\nSimulating regional IPM data for model testing...\n")
sim_regional <- simulate_regional_ipm(seed = 456)

cat(sprintf(
  "Simulation complete:\n  %d counties | %d sites | %d years\n",
  sim_regional$stan_data$C,
  sim_regional$stan_data$S,
  sim_regional$stan_data$T))

cat(sprintf(
  "  y_adult obs: %d | y_pup obs: %d | y_molt obs: %d\n",
  sum(sim_regional$stan_data$y_adult_obs),
  sum(sim_regional$stan_data$y_pup_obs),
  sum(sim_regional$stan_data$y_molt_obs)))

cat(sprintf(
  "  True delta_moci_mouth_fecund: %.3f  delta_moci_south_fecund: %.3f\n",
  sim_regional$true_params$delta_moci_mouth_fecund,
  sim_regional$true_params$delta_moci_south_fecund))

cat(sprintf("  2020 observations (should be 0): adult=%d, pup=%d, molt=%d\n",
            sum(sim_regional$stan_data$y_adult_obs[, 16]),  # index 16 = 2020
            sum(sim_regional$stan_data$y_pup_obs[,   16]),
            sum(sim_regional$stan_data$y_molt_obs[,  16])))

# Quick sanity check: site totals look biologically plausible
cat("\nSite N_adult at t=1 (selected sites):\n")
check_sites <- c(1, 9, 12, 14, 20, 24)  # BL, Castro, Mowry, Fitzgerald, Jenner, PointArena
for (si in check_sites) {
  n1 <- round(sim_regional$true_states$N_adult_total[si, 1])
  n_last <- round(sim_regional$true_states$N_adult_total[si, T])
  cat(sprintf("  %s: %d → %d\n",
              sim_regional$site_meta$site_name[si], n1, n_last))
}

# County totals (sum across sites)
cat("\nCounty total N_adult at t=1 and t=21 (summed across sites):\n")
for (ci in seq_len(6)) {
  sites_in_c <- which(sim_regional$site_meta$county_id == ci)
  n1   <- round(sum(sim_regional$true_states$N_adult_total[sites_in_c, 1]))
  n21  <- round(sum(sim_regional$true_states$N_adult_total[sites_in_c, T]))
  cat(sprintf("  %s: %d → %d\n", sim_regional$county_names[ci], n1, n21))
}

cat("\nObjects created: sim_regional\n")
cat("Next: source(\"Code/17_regional_ipm_model.R\") then fit on simulated data\n")
