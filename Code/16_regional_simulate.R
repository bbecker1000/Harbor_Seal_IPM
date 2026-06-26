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
simulate_regional_ipm <- function(T    = 21,           # 2005–2025 incl. 2020 latent
                                  S    = 24,
                                  C    = 6,
                                  T_proj = 10,
                                  seed = 456,
                                  site_meta = regional_sites,
                                  true_params = NULL) {
  set.seed(seed)
  
  if (is.null(true_params)) true_params <- list(
    # ── Vital rates (shared across all counties) ────────────────────────────
    phi_pup_logit      = qlogis(0.22),       # ~22% pup survival
    phi_juv_base       = 0.70,
    phi_adult_F_logit  = qlogis(0.93),       # ~93% adult female survival
    delta_adult        = 0.017,              # female survival advantage
    fecund_primip      = 0.60,
    fecund_mature      = 0.85,
    prop_female        = 0.50,
    rho_pup            = 0.18,               # pup molt attendance fraction
    p_male_fixed       = 0.057,              # fixed; not estimated
    
    # ── MOCI effects (open coast baseline) ─────────────────────────────────
    beta_moci_ond_fecund = -0.15,
    beta_moci_ond_pup    = -0.15,
    beta_moci_amj_pup    = -0.15,
    beta_moci_jfm_pup    = -0.15,
    beta_moci_jfm_juv    = -0.12,
    beta_moci_jfm_adult  = -0.10,
    beta_moci_amj_molt   =  0.05,
    
    # ── Bay vs. Coast MOCI modifier (KEY TEST PARAMETER) ───────────────────
    # Positive = Bay counties less negatively affected by poor MOCI
    # If recovery fails, the model cannot distinguish bay from coast responses
    delta_moci_mouth_fecund = 0.08,   # Bay Mouth less MOCI-sensitive than coast
    delta_moci_mouth_surv   = 0.06,
    delta_moci_south_fecund = 0.14,   # South Bay most decoupled (expected > mouth)
    delta_moci_south_surv   = 0.10,
    
    # ── County random effects ───────────────────────────────────────────────
    sigma_county     = 0.15,
    county_effect    = NULL,         # drawn below
    
    # ── Site availability (log scale; county-specific means) ─────────────────
    # These are set per-county in the Stan model using -log(n_sites_county[c]).
    # For simulation we use county-mean values consistent with that prior.
    mu_log_alpha_breed = -1.5,   # average across counties (weighted)
    mu_log_alpha_pup   = -1.2,
    mu_log_alpha_molt  = -1.5,
    sigma_log_alpha    = 0.30,   # Fix 3: tightened from 0.40
    
    # ── Error structure ─────────────────────────────────────────────────────
    sigma_process   = 0.15,
    sigma_obs_adult = 0.18,
    sigma_obs_pup   = 0.15,
    sigma_obs_molt  = 0.30
  )
  
  # ── County random effects ─────────────────────────────────────────────────
  if (is.null(true_params$county_effect))
    true_params$county_effect <- rnorm(C, 0, true_params$sigma_county)
  
  # ── Derived quantities ────────────────────────────────────────────────────
  phi_adult_M_base <- plogis(true_params$phi_adult_F_logit) - true_params$delta_adult
  avg_fecundity    <- 0.20 * true_params$fecund_primip + 0.80 * true_params$fecund_mature
  
  # ── County-level oceanographic type ──────────────────────────────────────
  # 0=coast, 1=bay mouth, 2=south bay — matches Stan data county_type vector
  county_type <- c(0L, 1L, 2L, 0L, 0L, 0L)
  
  # Year vector — defined locally to avoid global dependency on years_regional
  years_sim <- 2005:2025
  
  # ── MOCI time series (simulate realistic structure) ───────────────────────
  T_mid  <- 10
  moci_base <- c(seq(-0.5, 0.5, length.out = T_mid),
                 seq( 0.5,-0.5, length.out = T - T_mid)) +
    as.vector(arima.sim(list(ar = 0.5), n = T)) * 0.5
  moci_jfm  <- as.vector(scale(moci_base))
  moci_amj  <- as.vector(scale(moci_base * 0.80 +
                                 as.vector(arima.sim(list(ar = 0.3), n = T)) * 0.4))
  moci_ond  <- as.vector(scale(moci_base * 0.70 +
                                 as.vector(arima.sim(list(ar = 0.3), n = T)) * 0.5))
  
  # ── Site availability parameters (per site, log scale) ───────────────────
  log_alpha_breed <- rnorm(nrow(site_meta), true_params$mu_log_alpha_breed,
                           true_params$sigma_log_alpha)
  log_alpha_pup   <- rnorm(nrow(site_meta), true_params$mu_log_alpha_pup,
                           true_params$sigma_log_alpha)
  log_alpha_molt  <- rnorm(nrow(site_meta), true_params$mu_log_alpha_molt,
                           true_params$sigma_log_alpha)
  
  # Type H sites: set alpha_pup to effectively 0 (will be excluded from likelihood)
  log_alpha_pup[site_meta$is_breeder == 0] <- -10
  
  # ── Initialize county populations ────────────────────────────────────────
  # County-scale populations calibrated to raw count totals:
  #   Marin molt counts sum to ~1,000-1,500/yr across 8 sites
  #   Bay Mouth breed counts sum to ~300-500/yr across 3 sites
  #   South Bay breed counts sum to ~150-250/yr across 2 sites
  #   San Mateo breed counts sum to ~400-700/yr across 6 sites
  #   Sonoma breed counts sum to ~500-800/yr across 4 sites
  #   Mendocino breed counts ~100-150/yr at 1 site
  N_adult_F_init <- c(1200, 400, 200, 700, 650, 120)  # Marin, BayMouth, SouthBay, SMat, Son, Men
  N_adult_M_init <- N_adult_F_init * 0.9
  N_juv_F_init   <- N_adult_F_init * 0.30
  N_juv_M_init   <- N_juv_F_init
  N_pup_init     <- N_adult_F_init * 0.40
  
  # ── County-level state arrays ─────────────────────────────────────────────
  N_adult_F  <- N_adult_M  <- matrix(NA, C, T)
  N_juv_F    <- N_juv_M    <- matrix(NA, C, T)
  N_pup      <- matrix(NA, C, T)
  
  N_adult_F[, 1] <- N_adult_F_init; N_adult_M[, 1] <- N_adult_M_init
  N_juv_F[, 1]   <- N_juv_F_init;   N_juv_M[, 1]   <- N_juv_M_init
  N_pup[, 1]     <- N_pup_init
  
  # Time-varying vital rates (county × time)
  phi_pup_ct     <- phi_juv_ct <- phi_adult_F_ct <- phi_adult_M_ct <- matrix(NA, C, T)
  detect_molt_ct <- matrix(NA, C, T)
  
  for (c in 1:C) {
    for (t in 1:T) {
      t_birth <- max(t - 1, 1)
      
      # Three-way MOCI modifier by county type
      bay_fec <- (county_type[c] == 1) * true_params$delta_moci_mouth_fecund +
        (county_type[c] == 2) * true_params$delta_moci_south_fecund
      bay_sur <- (county_type[c] == 1) * true_params$delta_moci_mouth_surv +
        (county_type[c] == 2) * true_params$delta_moci_south_surv
      
      phi_pup_ct[c, t] <- plogis(
        true_params$phi_pup_logit + true_params$county_effect[c] +
          (true_params$beta_moci_amj_pup + bay_sur) * moci_amj[t_birth] +
          (true_params$beta_moci_ond_pup + bay_sur) * moci_ond[t] +
          (true_params$beta_moci_jfm_pup + bay_sur) * moci_jfm[t])
      
      phi_juv_ct[c, t] <- plogis(
        qlogis(true_params$phi_juv_base) + true_params$county_effect[c] * 0.5 +
          (true_params$beta_moci_jfm_juv + bay_sur) * moci_jfm[t])
      
      phi_adult_F_ct[c, t] <- plogis(
        true_params$phi_adult_F_logit + true_params$county_effect[c] * 0.25 +
          (true_params$beta_moci_jfm_adult + bay_sur) * moci_jfm[t])
      
      phi_adult_M_ct[c, t] <- plogis(
        qlogis(phi_adult_M_base) + true_params$county_effect[c] * 0.25 +
          (true_params$beta_moci_jfm_adult + bay_sur) * moci_jfm[t])
      
      detect_molt_ct[c, t] <- plogis(
        0.25 + true_params$beta_moci_amj_molt * moci_amj[t])
      
      if (t > 1) {
        fecund_t <- plogis(qlogis(avg_fecundity) +
                             (true_params$beta_moci_ond_fecund + bay_fec) * moci_ond[t])
        
        ep  <- N_adult_F[c, t-1] * fecund_t
        njF <- N_pup[c, t-1] * true_params$prop_female       * phi_pup_ct[c, t]
        njM <- N_pup[c, t-1] * (1 - true_params$prop_female) * phi_pup_ct[c, t]
        jsF <- N_juv_F[c, t-1] * phi_juv_ct[c, t] * (2/3)
        jsM <- N_juv_M[c, t-1] * phi_juv_ct[c, t] * (2/3)
        jaF <- N_juv_F[c, t-1] * phi_juv_ct[c, t] * (1/3)
        jaM <- N_juv_M[c, t-1] * phi_juv_ct[c, t] * (1/3)
        
        N_pup[c, t]     <- exp(rnorm(1, log(max(ep, 1)),
                                     true_params$sigma_process))
        N_juv_F[c, t]   <- exp(rnorm(1, log(max(njF + jsF, 0.1)),
                                     true_params$sigma_process * 0.5))
        N_juv_M[c, t]   <- exp(rnorm(1, log(max(njM + jsM, 0.1)),
                                     true_params$sigma_process * 0.5))
        N_adult_F[c, t] <- exp(rnorm(1, log(max(N_adult_F[c, t-1] * phi_adult_F_ct[c, t] + jaF, 1)),
                                     true_params$sigma_process * 0.5))
        N_adult_M[c, t] <- exp(rnorm(1, log(max(N_adult_M[c, t-1] * phi_adult_M_ct[c, t] + jaM, 1)),
                                     true_params$sigma_process * 0.5))
      }
    }
  }
  
  N_adult_total <- N_adult_F + N_adult_M
  N_juv_total   <- N_juv_F   + N_juv_M
  N_molt_true   <- N_juv_total + N_adult_total + true_params$rho_pup * N_pup
  
  # ── Generate site-level observations ────────────────────────────────────
  S_model <- nrow(site_meta)
  y_adult <- y_pup <- y_molt <- matrix(NA, S_model, T)
  
  for (s in 1:S_model) {
    c   <- site_meta$county_id[s]
    bay <- site_meta$is_bay[s]
    br  <- site_meta$is_breeder[s]
    
    for (t in 1:T) {
      # 2020 latent year: skip observations
      if (years_sim[t] == 2020) next
      # Before site started: skip
      if (years_sim[t] < site_meta$starts[s]) next
      
      alpha_b <- exp(log_alpha_breed[s])
      alpha_p <- exp(log_alpha_pup[s])
      alpha_m <- exp(log_alpha_molt[s])
      
      # Adult/breeding count
      N_adult_obs <- (N_adult_F[c, t] + N_adult_M[c, t] * true_params$p_male_fixed)
      if (N_adult_obs > 0 & alpha_b > 0)
        y_adult[s, t] <- log(N_adult_obs * alpha_b) +
        rnorm(1, 0, true_params$sigma_obs_adult)
      
      # Pup count (breeding sites only)
      if (br == 1 && N_pup[c, t] > 0 && alpha_p > 0)
        y_pup[s, t] <- log(N_pup[c, t] * alpha_p) +
        rnorm(1, 0, true_params$sigma_obs_pup)
      
      # Molt count (all sites)
      if (N_molt_true[c, t] > 0 && alpha_m > 0)
        y_molt[s, t] <- log(N_molt_true[c, t] * alpha_m * detect_molt_ct[c, t]) +
        rnorm(1, 0, true_params$sigma_obs_molt)
    }
  }
  
  # Introduce ~5% random missingness beyond structural gaps
  n_obs <- S_model * T
  set_na <- function(mat) {
    mat[sample(which(!is.na(mat)), round(0.05 * sum(!is.na(mat))))] <- NA; mat
  }
  y_adult <- set_na(y_adult)
  y_pup   <- set_na(y_pup)
  y_molt  <- set_na(y_molt)
  
  # Build indicator matrices and replace NA with 0
  y_adult_obs <- ifelse(!is.na(y_adult), 1L, 0L)
  y_pup_obs   <- ifelse(!is.na(y_pup),   1L, 0L)
  y_molt_obs  <- ifelse(!is.na(y_molt),  1L, 0L)
  y_adult[is.na(y_adult)] <- 0
  y_pup[is.na(y_pup)]     <- 0
  y_molt[is.na(y_molt)]   <- 0
  
  # Enforce Type H: no pup obs
  y_pup_obs[site_meta$is_breeder == 0, ] <- 0L
  y_pup[site_meta$is_breeder == 0, ]     <- 0
  
  # ── Projection scenarios ─────────────────────────────────────────────────
  N_scen <- 3
  moci_proj <- matrix(c(0, 1, -1), nrow = N_scen, ncol = T_proj)
  
  n_sites_county <- as.integer(table(factor(site_meta$county_id, levels = 1:C)))
  county_t1 <- vapply(seq_len(C), function(ci) {
    min_start <- min(site_meta$starts[site_meta$county_id == ci])
    which(2005:2025 == min_start)
  }, integer(1))
  
  stan_data <- list(
    T = T, S = S_model, C = C, T_proj = T_proj, N_scenarios = N_scen,
    y_adult = y_adult, y_pup = y_pup, y_molt = y_molt,
    y_adult_obs = y_adult_obs, y_pup_obs = y_pup_obs, y_molt_obs = y_molt_obs,
    county_id      = site_meta$county_id,
    is_breeder     = site_meta$is_breeder,
    county_t1      = county_t1,
    n_sites_county = n_sites_county,
    county_type    = county_type,
    moci_jfm = moci_jfm, moci_amj = moci_amj, moci_ond = moci_ond,
    moci_proj = moci_proj,
    p_male_fixed = 0.057
  )
  
  true_params$county_effect <- true_params$county_effect  # ensure saved
  true_params$log_alpha_breed <- log_alpha_breed
  true_params$log_alpha_pup   <- log_alpha_pup
  true_params$log_alpha_molt  <- log_alpha_molt
  true_params$county_type     <- county_type
  
  list(
    stan_data   = stan_data,
    true_params = true_params,
    true_states = list(
      N_adult_F = N_adult_F, N_adult_M = N_adult_M,
      N_juv_F = N_juv_F, N_juv_M = N_juv_M, N_pup = N_pup,
      N_adult_total = N_adult_total, N_juv_total = N_juv_total,
      N_molt_true = N_molt_true,
      phi_pup = phi_pup_ct, phi_juv = phi_juv_ct,
      phi_adult_F = phi_adult_F_ct, phi_adult_M = phi_adult_M_ct
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

# Quick sanity check: county trajectories look biologically plausible
cat("\nCounty total N (adults) at t=1 and t=21:\n")
for (c in 1:6) {
  n1 <- round(sum(sim_regional$true_states$N_adult_total[c, 1]))
  n21 <- round(sum(sim_regional$true_states$N_adult_total[c, 21]))
  cat(sprintf("  %s: %d → %d\n",
              sim_regional$county_names[c], n1, n21))
}

cat("\nObjects created: sim_regional\n")
cat("Next: source(\"Code/17_regional_ipm_model.R\") then fit on simulated data\n")
