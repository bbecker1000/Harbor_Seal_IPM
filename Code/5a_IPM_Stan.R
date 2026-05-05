
# ============================================================================
# HARBOR SEAL INTEGRATED POPULATION MODEL v3.1
# Sex-Specific Survival, Age-Dependent Fecundity, 
# SITE-SPECIFIC Coyote & Disturbance Effects, and 10-Year Projections
# ============================================================================
# 
# CORRECTED: Site-specific covariate effects matching MARSS structure:
#   - beta_coy[3]: Separate coyote effects for BL, DE, DP
#   - beta_dist[6]: Separate disturbance effects for each site
#
# Key Features:
# - Sex-specific survival (females ~5-10% higher than males)
# - Age-specific fecundity (primiparous < multiparous < prime < senescent)
# - Three observation types: Adults (breeding), Pups, Molts
# - Age structure: Pup → Juvenile (1-3 yr) → Adult (4+ yr)
# - Site-specific covariates: Coyote (BL, DE, DP), Disturbance (all sites)
# - Shared covariates: MOCI, Elephant seals
# - 10-year population projections under climate scenarios
# ============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(patchwork)

# Create output directories
dir.create("Output", showWarnings = FALSE)
dir.create("Output/Plots", showWarnings = FALSE)

# ============================================================================
# SITE INDEXING REFERENCE
# ============================================================================
# Site indices: 1=BL, 2=DE, 3=DP, 4=PRH, 5=TB, 6=TP
#
# Coyote sites: BL(1), DE(2), DP(3) → indices 1, 2, 3
# Disturbance sites: All 6 sites have disturbance data
# Elephant seal sites: DE(2), PRH(4)
# ============================================================================

# ============================================================================
# PART 1: STAN MODEL CODE (v3.1 - Site-Specific Effects)
# ============================================================================

stan_code_v3.1 <- '
data {
  int<lower=1> T;                         // Number of years
  int<lower=1> S;                         // Number of sites (6)
  int<lower=1> N_coy;                     // Number of coyote sites (3: BL, DE, DP)
  
  // Observed log counts
  matrix[S, T] y_adult;
  matrix[S, T] y_pup;
  matrix[S, T] y_molt;
  
  // Observation indicators (1 = observed, 0 = missing)
  array[S, T] int<lower=0, upper=1> y_adult_obs;
  array[S, T] int<lower=0, upper=1> y_pup_obs;
  array[S, T] int<lower=0, upper=1> y_molt_obs;
  
  // Covariates - SITE-SPECIFIC
  matrix[S, T] coyote;                    // Coyote data (0 for non-coyote sites)
  matrix[S, T] disturbance;               // Disturbance data (all sites)
  matrix[S, T] elephant_seal;             // Elephant seal data (0 for non-eseal sites)
  
  // Covariates - SHARED (time-varying only)
  vector[T] moci_jfm;
  vector[T] moci_amj;
  vector[T] moci_ond;
  
  // Site indicators
  array[S] int<lower=0, upper=N_coy> coyote_idx;  // Maps site to coyote param (0 if no coyote)
  array[S] int<lower=0, upper=1> has_eseal;
  
  // Projection scenarios
  int<lower=0> T_proj;
  int<lower=1> N_scenarios;
  matrix[N_scenarios, T_proj] moci_proj;
  matrix[N_scenarios, T_proj] coyote_proj;  // Average coyote level for projections
}

parameters {
  // ===========================================
  // SEX-SPECIFIC BASELINE SURVIVAL RATES
  // ===========================================
  
  real<lower=0, upper=1> phi_pup_F_base;
  real<lower=0, upper=1> phi_juv_F_base;
  real<lower=0, upper=1> phi_adult_F_base;
  
  real<lower=0, upper=0.25> delta_pup;
  real<lower=0, upper=0.15> delta_juv;
  real<lower=0, upper=0.10> delta_adult;
  
  // ===========================================
  // AGE-SPECIFIC REPRODUCTION
  // ===========================================
  
  real<lower=0, upper=1> fecund_primip;
  real<lower=0, upper=1> fecund_young;
  real<lower=0, upper=1> fecund_prime;
  real<lower=0, upper=1> fecund_senior;
  real<lower=0.4, upper=0.6> prop_female;
  
  // ===========================================
  // SITE-SPECIFIC COVARIATE EFFECTS
  // ===========================================
  
  // Coyote effects - SITE SPECIFIC (BL, DE, DP)
  vector[N_coy] beta_coy;                 // Site-specific coyote effects on survival
  
  // Disturbance effects - SITE SPECIFIC (all 6 sites)
  vector[S] beta_dist_surv;               // Disturbance effect on pup survival
  vector[S] beta_dist_detect;             // Disturbance effect on detection
  
  // ===========================================
  // SHARED COVARIATE EFFECTS
  // ===========================================
  
  // MOCI effects (shared across sites)
  real beta_moci_ond_pup;
  real beta_moci_amj_pup;
  real beta_moci_jfm_juv;
  real beta_moci_jfm_adult;
  real beta_moci_amj_molt;
  
  // Elephant seal effect
  real beta_eseal_pup;
  
  // ===========================================
  // SITE RANDOM EFFECTS (non-centered)
  // ===========================================
  
  vector[S] site_effect_raw;
  real<lower=0.01, upper=1.5> sigma_site;  // Bounded to prevent funnel
  
  // ===========================================
  // ERROR TERMS - bounded to improve sampling
  // ===========================================
  
  real<lower=0.05, upper=0.5> sigma_process;    // Bounded process error
  real<lower=0.05, upper=0.4> sigma_obs_adult;  // Bounded obs error
  real<lower=0.02, upper=0.35> sigma_obs_pup;   // Tight bounds - was problematic
  real<lower=0.05, upper=0.6> sigma_obs_molt;   // Bounded obs error
  
  // ===========================================
  // INITIAL POPULATIONS (log scale, non-centered)
  // ===========================================
  
  vector[S] log_N_adult_F_init_raw;
  vector[S] log_N_adult_M_init_raw;
  vector[S] log_N_juv_init_raw;
  vector[S] log_N_pup_init_raw;
  
  real mu_log_adult;
  real mu_log_juv;
  real mu_log_pup;
  real<lower=0> sigma_init;
  
  // ===========================================
  // PROCESS ERRORS (non-centered)
  // ===========================================
  
  matrix[S, T-1] eps_adult_raw;
  matrix[S, T-1] eps_juv_raw;
  matrix[S, T-1] eps_pup_raw;
}

transformed parameters {
  // Site random effects
  vector[S] site_effect = sigma_site * site_effect_raw;
  
  // Derived sex-specific baseline survival
  real phi_pup_M_base = fmax(phi_pup_F_base - delta_pup, 0.01);
  real phi_juv_M_base = fmax(phi_juv_F_base - delta_juv, 0.01);
  real phi_adult_M_base = fmax(phi_adult_F_base - delta_adult, 0.01);
  
  // Initial populations (non-centered)
  vector<lower=0>[S] N_adult_F_init;
  vector<lower=0>[S] N_adult_M_init;
  vector<lower=0>[S] N_juv_F_init;
  vector<lower=0>[S] N_juv_M_init;
  vector<lower=0>[S] N_pup_init;
  
  for (s in 1:S) {
    N_adult_F_init[s] = exp(mu_log_adult + sigma_init * log_N_adult_F_init_raw[s]);
    N_adult_M_init[s] = exp(mu_log_adult + sigma_init * log_N_adult_M_init_raw[s]) * 0.9;
    N_juv_F_init[s] = exp(mu_log_juv + sigma_init * log_N_juv_init_raw[s]) * 0.5;
    N_juv_M_init[s] = exp(mu_log_juv + sigma_init * log_N_juv_init_raw[s]) * 0.5;
    N_pup_init[s] = exp(mu_log_pup + sigma_init * log_N_pup_init_raw[s]);
  }
  
  // Latent populations
  matrix<lower=0>[S, T] N_adult_F;
  matrix<lower=0>[S, T] N_adult_M;
  matrix<lower=0>[S, T] N_juv_F;
  matrix<lower=0>[S, T] N_juv_M;
  matrix<lower=0>[S, T] N_pup;
  
  matrix<lower=0>[S, T] N_adult_total;
  matrix<lower=0>[S, T] N_juv_total;
  matrix<lower=0>[S, T] N_molt_true;
  matrix<lower=0>[S, T] N_total;
  
  // Time-varying vital rates
  matrix<lower=0, upper=1>[S, T] phi_pup_F;
  matrix<lower=0, upper=1>[S, T] phi_pup_M;
  matrix<lower=0, upper=1>[S, T] phi_juv_F;
  matrix<lower=0, upper=1>[S, T] phi_juv_M;
  matrix<lower=0, upper=1>[S, T] phi_adult_F;
  matrix<lower=0, upper=1>[S, T] phi_adult_M;
  matrix<lower=0, upper=1>[S, T] detect_breed;
  matrix<lower=0, upper=1>[S, T] detect_molt;
  
  // Average fecundity
  real avg_fecundity = 0.30 * fecund_primip + 0.20 * fecund_young + 
                       0.40 * fecund_prime + 0.10 * fecund_senior;
  
  // ===========================================
  // CALCULATE TIME-VARYING VITAL RATES
  // ===========================================
  
  for (s in 1:S) {
    for (t in 1:T) {
      // -----------------------------------------
      // COYOTE EFFECT (site-specific)
      // -----------------------------------------
      real coyote_effect = 0;
      if (coyote_idx[s] > 0) {
        coyote_effect = beta_coy[coyote_idx[s]] * coyote[s, t];
      }
      
      // -----------------------------------------
      // DISTURBANCE EFFECT (site-specific)
      // -----------------------------------------
      real dist_surv_effect = beta_dist_surv[s] * disturbance[s, t];
      real dist_detect_effect = beta_dist_detect[s] * disturbance[s, t];
      
      // -----------------------------------------
      // FEMALE SURVIVAL RATES
      // -----------------------------------------
      
      // Pup survival (female)
      real logit_phi_pup_F = logit(phi_pup_F_base) + site_effect[s] +
        coyote_effect +
        beta_moci_ond_pup * moci_ond[t] +
        beta_moci_amj_pup * moci_amj[t] +
        has_eseal[s] * beta_eseal_pup * elephant_seal[s, t] +
        dist_surv_effect;
      phi_pup_F[s, t] = inv_logit(logit_phi_pup_F);
      
      // Juvenile survival (female) - no coyote/disturbance effect
      real logit_phi_juv_F = logit(phi_juv_F_base) + site_effect[s] * 0.5 +
        beta_moci_jfm_juv * moci_jfm[t];
      phi_juv_F[s, t] = inv_logit(logit_phi_juv_F);
      
      // Adult survival (female) - no coyote/disturbance effect
      real logit_phi_adult_F = logit(phi_adult_F_base) + site_effect[s] * 0.25 +
        beta_moci_jfm_adult * moci_jfm[t];
      phi_adult_F[s, t] = inv_logit(logit_phi_adult_F);
      
      // -----------------------------------------
      // MALE SURVIVAL RATES (same effects, lower baseline)
      // -----------------------------------------
      
      real logit_phi_pup_M = logit(phi_pup_M_base) + site_effect[s] +
        coyote_effect +
        beta_moci_ond_pup * moci_ond[t] +
        beta_moci_amj_pup * moci_amj[t] +
        has_eseal[s] * beta_eseal_pup * elephant_seal[s, t] +
        dist_surv_effect;
      phi_pup_M[s, t] = inv_logit(logit_phi_pup_M);
      
      real logit_phi_juv_M = logit(phi_juv_M_base) + site_effect[s] * 0.5 +
        beta_moci_jfm_juv * moci_jfm[t];
      phi_juv_M[s, t] = inv_logit(logit_phi_juv_M);
      
      real logit_phi_adult_M = logit(phi_adult_M_base) + site_effect[s] * 0.25 +
        beta_moci_jfm_adult * moci_jfm[t];
      phi_adult_M[s, t] = inv_logit(logit_phi_adult_M);
      
      // -----------------------------------------
      // DETECTION PROBABILITIES (site-specific disturbance)
      // -----------------------------------------
      
      detect_breed[s, t] = inv_logit(dist_detect_effect);
      detect_molt[s, t] = inv_logit(dist_detect_effect + beta_moci_amj_molt * moci_amj[t]);
    }
  }
  
  // ===========================================
  // POPULATION DYNAMICS
  // ===========================================
  
  for (s in 1:S) {
    N_adult_F[s, 1] = N_adult_F_init[s];
    N_adult_M[s, 1] = N_adult_M_init[s];
    N_juv_F[s, 1] = N_juv_F_init[s];
    N_juv_M[s, 1] = N_juv_M_init[s];
    N_pup[s, 1] = N_pup_init[s];
    
    N_adult_total[s, 1] = N_adult_F[s, 1] + N_adult_M[s, 1];
    N_juv_total[s, 1] = N_juv_F[s, 1] + N_juv_M[s, 1];
    N_molt_true[s, 1] = N_juv_total[s, 1] + N_adult_total[s, 1];
    N_total[s, 1] = N_pup[s, 1] + N_juv_total[s, 1] + N_adult_total[s, 1];
    
    for (t in 2:T) {
      real expected_pups = N_adult_F[s, t-1] * avg_fecundity;
      
      real new_juv_F = N_pup[s, t-1] * prop_female * phi_pup_F[s, t];
      real new_juv_M = N_pup[s, t-1] * (1 - prop_female) * phi_pup_M[s, t];
      
      real juv_stay_F = N_juv_F[s, t-1] * phi_juv_F[s, t] * (2.0/3.0);
      real juv_stay_M = N_juv_M[s, t-1] * phi_juv_M[s, t] * (2.0/3.0);
      
      real juv_to_adult_F = N_juv_F[s, t-1] * phi_juv_F[s, t] * (1.0/3.0);
      real juv_to_adult_M = N_juv_M[s, t-1] * phi_juv_M[s, t] * (1.0/3.0);
      
      real expected_juv_F = new_juv_F + juv_stay_F;
      real expected_juv_M = new_juv_M + juv_stay_M;
      
      real expected_adult_F = N_adult_F[s, t-1] * phi_adult_F[s, t] + juv_to_adult_F;
      real expected_adult_M = N_adult_M[s, t-1] * phi_adult_M[s, t] + juv_to_adult_M;
      
      N_pup[s, t] = exp(log(fmax(expected_pups, 1)) + 
                        sigma_process * eps_pup_raw[s, t-1]);
      N_juv_F[s, t] = exp(log(fmax(expected_juv_F, 0.1)) + 
                          sigma_process * eps_juv_raw[s, t-1] * 0.5);
      N_juv_M[s, t] = exp(log(fmax(expected_juv_M, 0.1)) + 
                          sigma_process * eps_juv_raw[s, t-1] * 0.5);
      N_adult_F[s, t] = exp(log(fmax(expected_adult_F, 1)) + 
                            sigma_process * eps_adult_raw[s, t-1] * 0.5);
      N_adult_M[s, t] = exp(log(fmax(expected_adult_M, 1)) + 
                            sigma_process * eps_adult_raw[s, t-1] * 0.5);
      
      N_adult_total[s, t] = N_adult_F[s, t] + N_adult_M[s, t];
      N_juv_total[s, t] = N_juv_F[s, t] + N_juv_M[s, t];
      N_molt_true[s, t] = N_juv_total[s, t] + N_adult_total[s, t];
      N_total[s, t] = N_pup[s, t] + N_juv_total[s, t] + N_adult_total[s, t];
    }
  }
}

model {
  // ===========================================
  // PRIORS
  // ===========================================
  
  // Female survival
  phi_pup_F_base ~ beta(16, 4);
  phi_juv_F_base ~ beta(17, 3);
  phi_adult_F_base ~ beta(18, 2);
  
  // Male survival difference
  delta_pup ~ normal(0.10, 0.025);
  delta_juv ~ normal(0.08, 0.020);
  delta_adult ~ normal(0.05, 0.015);
  
  // Fecundity
  fecund_primip ~ beta(14, 6);
  fecund_young ~ beta(17, 3);
  fecund_prime ~ beta(19, 1);
  fecund_senior ~ beta(16, 4);
  prop_female ~ beta(50, 50);
  
  // SITE-SPECIFIC COYOTE EFFECTS (informed by MARSS results)
  // MARSS found: BL = -3.3%, DE = -8.5%, DP = -0.9%
  beta_coy[1] ~ normal(-0.15, 0.15);    // BL - weak effect
  beta_coy[2] ~ normal(-0.40, 0.20);    // DE - strong effect
  beta_coy[3] ~ normal(-0.05, 0.15);    // DP - near zero
  
  // SITE-SPECIFIC DISTURBANCE EFFECTS ON SURVIVAL
  // MARSS found: BL = -4.3%, DE = -6.8%, DP = ?, PRH = ?, TB = -5.5%, TP = ?
  beta_dist_surv[1] ~ normal(-0.20, 0.15);   // BL
  beta_dist_surv[2] ~ normal(-0.30, 0.15);   // DE
  beta_dist_surv[3] ~ normal(-0.15, 0.15);   // DP
  beta_dist_surv[4] ~ normal(-0.10, 0.15);   // PRH (limited disturbance)
  beta_dist_surv[5] ~ normal(-0.25, 0.15);   // TB
  beta_dist_surv[6] ~ normal(-0.15, 0.15);   // TP
  
  // SITE-SPECIFIC DISTURBANCE EFFECTS ON DETECTION
  beta_dist_detect ~ normal(-0.20, 0.15);
  
  // SHARED COVARIATE EFFECTS
  beta_moci_ond_pup ~ normal(-0.25, 0.15);
  beta_moci_amj_pup ~ normal(-0.15, 0.15);
  beta_moci_jfm_juv ~ normal(-0.10, 0.15);
  beta_moci_jfm_adult ~ normal(-0.08, 0.10);
  beta_moci_amj_molt ~ normal(0.05, 0.15);
  beta_eseal_pup ~ normal(0.10, 0.20);
  
  // Random effects and error - tighter priors for better mixing
  sigma_site ~ normal(0.5, 0.3);         // Informative prior, centered on reasonable value
  site_effect_raw ~ std_normal();
  
  sigma_process ~ normal(0.15, 0.08);    // Informative prior ~15% process variation
  sigma_obs_adult ~ normal(0.18, 0.06);  // Based on previous estimate
  sigma_obs_pup ~ normal(0.12, 0.05);    // Informative - this was problematic
  sigma_obs_molt ~ normal(0.35, 0.10);   // Based on previous estimate
  
  // Initial populations
  mu_log_adult ~ normal(5, 0.5);
  mu_log_juv ~ normal(4, 0.5);
  mu_log_pup ~ normal(4, 0.5);
  sigma_init ~ exponential(3);
  
  log_N_adult_F_init_raw ~ std_normal();
  log_N_adult_M_init_raw ~ std_normal();
  log_N_juv_init_raw ~ std_normal();
  log_N_pup_init_raw ~ std_normal();
  
  // Process errors
  to_vector(eps_adult_raw) ~ std_normal();
  to_vector(eps_juv_raw) ~ std_normal();
  to_vector(eps_pup_raw) ~ std_normal();
  
  // ===========================================
  // LIKELIHOOD
  // ===========================================
  
  for (s in 1:S) {
    for (t in 1:T) {
      if (y_adult_obs[s, t] == 1) {
        y_adult[s, t] ~ normal(
          log(N_adult_total[s, t] * detect_breed[s, t]),
          sigma_obs_adult
        );
      }
      
      if (y_pup_obs[s, t] == 1) {
        y_pup[s, t] ~ normal(
          log(N_pup[s, t] * detect_breed[s, t]),
          sigma_obs_pup
        );
      }
      
      if (y_molt_obs[s, t] == 1) {
        y_molt[s, t] ~ normal(
          log(N_molt_true[s, t] * detect_molt[s, t]),
          sigma_obs_molt
        );
      }
    }
  }
}

generated quantities {
  // ===========================================
  // POSTERIOR PREDICTIVE
  // ===========================================
  
  matrix[S, T] y_adult_rep;
  matrix[S, T] y_pup_rep;
  matrix[S, T] y_molt_rep;
  
  // ===========================================
  // DERIVED QUANTITIES
  // ===========================================
  
  matrix[S, T] sex_ratio_adult;
  matrix[S, T-1] lambda;
  vector[T] N_total_all;
  
  vector[T] mean_phi_pup_F;
  vector[T] mean_phi_pup_M;
  vector[T] mean_phi_juv_F;
  vector[T] mean_phi_juv_M;
  vector[T] mean_phi_adult_F;
  vector[T] mean_phi_adult_M;
  
  // ===========================================
  // 10-YEAR POPULATION PROJECTIONS
  // ===========================================
  
  array[N_scenarios] matrix[S, T_proj] N_total_proj;
  array[N_scenarios] matrix[S, T_proj] N_pup_proj;
  array[N_scenarios] matrix[S, T_proj] N_adult_proj;
  array[N_scenarios] vector[T_proj] N_total_all_proj;
  array[N_scenarios] vector[T_proj-1] lambda_proj;
  
  // Generate replicated data
  for (s in 1:S) {
    for (t in 1:T) {
      y_adult_rep[s, t] = normal_rng(
        log(N_adult_total[s, t] * detect_breed[s, t]), 
        sigma_obs_adult
      );
      y_pup_rep[s, t] = normal_rng(
        log(N_pup[s, t] * detect_breed[s, t]), 
        sigma_obs_pup
      );
      y_molt_rep[s, t] = normal_rng(
        log(N_molt_true[s, t] * detect_molt[s, t]), 
        sigma_obs_molt
      );
      
      sex_ratio_adult[s, t] = N_adult_F[s, t] / N_adult_total[s, t];
    }
  }
  
  // Lambda
  for (s in 1:S) {
    for (t in 1:(T-1)) {
      lambda[s, t] = N_total[s, t+1] / N_total[s, t];
    }
  }
  
  // Totals and means
  for (t in 1:T) {
    N_total_all[t] = sum(col(N_total, t));
    mean_phi_pup_F[t] = mean(col(phi_pup_F, t));
    mean_phi_pup_M[t] = mean(col(phi_pup_M, t));
    mean_phi_juv_F[t] = mean(col(phi_juv_F, t));
    mean_phi_juv_M[t] = mean(col(phi_juv_M, t));
    mean_phi_adult_F[t] = mean(col(phi_adult_F, t));
    mean_phi_adult_M[t] = mean(col(phi_adult_M, t));
  }
  
  // ===========================================
  // POPULATION PROJECTIONS
  // ===========================================
  
  for (scen in 1:N_scenarios) {
    matrix[S, T_proj] proj_adult_F;
    matrix[S, T_proj] proj_adult_M;
    matrix[S, T_proj] proj_juv_F;
    matrix[S, T_proj] proj_juv_M;
    matrix[S, T_proj] proj_pup;
    
    for (s in 1:S) {
      proj_adult_F[s, 1] = N_adult_F[s, T];
      proj_adult_M[s, 1] = N_adult_M[s, T];
      proj_juv_F[s, 1] = N_juv_F[s, T];
      proj_juv_M[s, 1] = N_juv_M[s, T];
      proj_pup[s, 1] = N_pup[s, T];
      
      N_total_proj[scen][s, 1] = proj_pup[s, 1] + proj_juv_F[s, 1] + 
                                  proj_juv_M[s, 1] + proj_adult_F[s, 1] + proj_adult_M[s, 1];
      N_pup_proj[scen][s, 1] = proj_pup[s, 1];
      N_adult_proj[scen][s, 1] = proj_adult_F[s, 1] + proj_adult_M[s, 1];
      
      for (tp in 2:T_proj) {
        // Site-specific coyote effect for projections
        real proj_coy_effect = 0;
        if (coyote_idx[s] > 0) {
          proj_coy_effect = beta_coy[coyote_idx[s]] * coyote_proj[scen, tp];
        }
        
        // Calculate survival under scenario
        real proj_phi_pup_F = inv_logit(logit(phi_pup_F_base) + site_effect[s] +
          proj_coy_effect +
          beta_moci_ond_pup * moci_proj[scen, tp] +
          beta_moci_amj_pup * moci_proj[scen, tp]);
        
        real proj_phi_pup_M = inv_logit(logit(phi_pup_M_base) + site_effect[s] +
          proj_coy_effect +
          beta_moci_ond_pup * moci_proj[scen, tp] +
          beta_moci_amj_pup * moci_proj[scen, tp]);
        
        real proj_phi_juv_F = inv_logit(logit(phi_juv_F_base) + site_effect[s] * 0.5 +
          beta_moci_jfm_juv * moci_proj[scen, tp]);
        
        real proj_phi_juv_M = inv_logit(logit(phi_juv_M_base) + site_effect[s] * 0.5 +
          beta_moci_jfm_juv * moci_proj[scen, tp]);
        
        real proj_phi_adult_F = inv_logit(logit(phi_adult_F_base) + site_effect[s] * 0.25 +
          beta_moci_jfm_adult * moci_proj[scen, tp]);
        
        real proj_phi_adult_M = inv_logit(logit(phi_adult_M_base) + site_effect[s] * 0.25 +
          beta_moci_jfm_adult * moci_proj[scen, tp]);
        
        // Deterministic projection
        real new_pups = proj_adult_F[s, tp-1] * avg_fecundity;
        
        real new_juv_F = proj_pup[s, tp-1] * prop_female * proj_phi_pup_F;
        real new_juv_M = proj_pup[s, tp-1] * (1 - prop_female) * proj_phi_pup_M;
        
        real juv_stay_F = proj_juv_F[s, tp-1] * proj_phi_juv_F * (2.0/3.0);
        real juv_stay_M = proj_juv_M[s, tp-1] * proj_phi_juv_M * (2.0/3.0);
        
        real juv_to_adult_F = proj_juv_F[s, tp-1] * proj_phi_juv_F * (1.0/3.0);
        real juv_to_adult_M = proj_juv_M[s, tp-1] * proj_phi_juv_M * (1.0/3.0);
        
        proj_pup[s, tp] = new_pups;
        proj_juv_F[s, tp] = new_juv_F + juv_stay_F;
        proj_juv_M[s, tp] = new_juv_M + juv_stay_M;
        proj_adult_F[s, tp] = proj_adult_F[s, tp-1] * proj_phi_adult_F + juv_to_adult_F;
        proj_adult_M[s, tp] = proj_adult_M[s, tp-1] * proj_phi_adult_M + juv_to_adult_M;
        
        N_total_proj[scen][s, tp] = proj_pup[s, tp] + proj_juv_F[s, tp] + 
                                     proj_juv_M[s, tp] + proj_adult_F[s, tp] + proj_adult_M[s, tp];
        N_pup_proj[scen][s, tp] = proj_pup[s, tp];
        N_adult_proj[scen][s, tp] = proj_adult_F[s, tp] + proj_adult_M[s, tp];
      }
    }
    
    for (tp in 1:T_proj) {
      N_total_all_proj[scen][tp] = sum(col(N_total_proj[scen], tp));
    }
    
    for (tp in 1:(T_proj-1)) {
      lambda_proj[scen][tp] = N_total_all_proj[scen][tp+1] / N_total_all_proj[scen][tp];
    }
  }
}
'

# Write Stan model to file
write_lines(stan_code_v3.1, "harbor_seal_ipm_v3.1.stan")
cat("Stan model written to harbor_seal_ipm_v3.1.stan\n")


# ============================================================================
# PART 2: SIMULATE DATA WITH KNOWN PARAMETERS
# ============================================================================

simulate_seal_ipm_data_v3.1 <- function(
    T = 29,
    S = 6,
    T_proj = 10,
    seed = 123
) {
  
  set.seed(seed)
  
  site_names <- c("BL", "DE", "DP", "PRH", "TB", "TP")
  N_coy <- 3  # Number of coyote sites
  
  # -----------------------------------------
  # TRUE PARAMETER VALUES
  # -----------------------------------------
  
  true_params <- list(
    # Sex-specific survival
    phi_pup_F_base = 0.80,
    phi_juv_F_base = 0.87,
    phi_adult_F_base = 0.92,
    
    delta_pup = 0.10,
    delta_juv = 0.08,
    delta_adult = 0.05,
    
    phi_pup_M_base = 0.70,
    phi_juv_M_base = 0.79,
    phi_adult_M_base = 0.87,
    
    # Fecundity
    fecund_primip = 0.70,
    fecund_young = 0.85,
    fecund_prime = 0.95,
    fecund_senior = 0.80,
    avg_fecundity = 0.30 * 0.70 + 0.20 * 0.85 + 0.40 * 0.95 + 0.10 * 0.80,
    
    prop_female = 0.50,
    
    # SITE-SPECIFIC COYOTE EFFECTS (based on MARSS results)
    # BL = -3.3%, DE = -8.5%, DP = -0.9%
    beta_coy = c(-0.15, -0.40, -0.05),  # BL, DE, DP
    
    # SITE-SPECIFIC DISTURBANCE EFFECTS ON SURVIVAL
    # BL, DE, DP, PRH, TB, TP
    beta_dist_surv = c(-0.20, -0.30, -0.15, -0.10, -0.25, -0.15),
    
    # SITE-SPECIFIC DISTURBANCE EFFECTS ON DETECTION
    beta_dist_detect = c(-0.20, -0.20, -0.15, -0.10, -0.25, -0.15),
    
    # Shared effects
    beta_moci_ond_pup = -0.25,
    beta_moci_amj_pup = -0.15,
    beta_moci_jfm_juv = -0.10,
    beta_moci_jfm_adult = -0.08,
    beta_moci_amj_molt = 0.05,
    beta_eseal_pup = 0.10,
    
    # Error terms
    sigma_process = 0.06,
    sigma_obs_adult = 0.12,
    sigma_obs_pup = 0.15,
    sigma_obs_molt = 0.12,
    sigma_site = 0.10
  )
  
  # -----------------------------------------
  # SITE INDICATORS
  # -----------------------------------------
  
  # Maps site to coyote parameter index (0 if no coyote)
  coyote_idx <- c(1, 2, 3, 0, 0, 0)  # BL=1, DE=2, DP=3, others=0
  has_eseal <- c(0, 1, 0, 1, 0, 0)   # DE and PRH
  
  # -----------------------------------------
  # SIMULATE COVARIATES
  # -----------------------------------------
  
  # Coyote - site-specific time series
  coyote <- matrix(0, S, T)
  coyote[1, ] <- scale(cumsum(rnorm(T, 0.03, 0.3)))[,1]  # BL
  coyote[2, ] <- scale(cumsum(rnorm(T, 0.05, 0.3)))[,1]  # DE - increasing trend
  coyote[3, ] <- scale(cumsum(rnorm(T, 0.02, 0.3)))[,1]  # DP
  
  # Disturbance - all sites
  disturbance <- matrix(0, S, T)
  for (s in 1:S) {
    disturbance[s, ] <- scale(rnorm(T, 0, 1))[,1]
  }
  
  # MOCI - shared
  moci_jfm <- scale(arima.sim(list(ar = 0.5), n = T))[,1]
  moci_amj <- scale(arima.sim(list(ar = 0.5), n = T))[,1]
  moci_ond <- scale(arima.sim(list(ar = 0.5), n = T))[,1]
  
  # Elephant seals - DE and PRH only
  elephant_seal <- matrix(0, S, T)
  elephant_seal[2, ] <- scale(seq(0, 3, length.out = T) + rnorm(T, 0, 0.5))[,1]  # DE
  elephant_seal[4, ] <- scale(seq(0, 4, length.out = T) + rnorm(T, 0, 0.5))[,1]  # PRH
  
  # Site effects
  site_effect <- rnorm(S, 0, true_params$sigma_site)
  
  # -----------------------------------------
  # SIMULATE POPULATIONS
  # -----------------------------------------
  
  N_adult_F <- N_adult_M <- matrix(NA, S, T)
  N_juv_F <- N_juv_M <- matrix(NA, S, T)
  N_pup <- matrix(NA, S, T)
  
  phi_pup_F <- phi_pup_M <- matrix(NA, S, T)
  phi_juv_F <- phi_juv_M <- matrix(NA, S, T)
  phi_adult_F <- phi_adult_M <- matrix(NA, S, T)
  detect_breed <- detect_molt <- matrix(NA, S, T)
  
  # Initial populations
  init_adult_F <- c(180, 140, 70, 95, 115, 45)
  init_adult_M <- c(160, 125, 63, 85, 100, 40)
  init_juv_F <- c(55, 42, 21, 28, 35, 14)
  init_juv_M <- c(50, 38, 19, 25, 32, 12)
  init_pup <- c(130, 100, 50, 67, 82, 32)
  
  N_adult_F[, 1] <- init_adult_F
  N_adult_M[, 1] <- init_adult_M
  N_juv_F[, 1] <- init_juv_F
  N_juv_M[, 1] <- init_juv_M
  N_pup[, 1] <- init_pup
  
  for (s in 1:S) {
    for (t in 1:T) {
      # Site-specific coyote effect
      coy_effect <- 0
      if (coyote_idx[s] > 0) {
        coy_effect <- true_params$beta_coy[coyote_idx[s]] * coyote[s, t]
      }
      
      # Site-specific disturbance effect
      dist_surv_effect <- true_params$beta_dist_surv[s] * disturbance[s, t]
      dist_detect_effect <- true_params$beta_dist_detect[s] * disturbance[s, t]
      
      # Survival rates
      logit_phi_pup_F <- qlogis(true_params$phi_pup_F_base) + site_effect[s] +
        coy_effect +
        true_params$beta_moci_ond_pup * moci_ond[t] +
        true_params$beta_moci_amj_pup * moci_amj[t] +
        has_eseal[s] * true_params$beta_eseal_pup * elephant_seal[s, t] +
        dist_surv_effect
      phi_pup_F[s, t] <- plogis(logit_phi_pup_F)
      
      logit_phi_pup_M <- qlogis(true_params$phi_pup_M_base) + site_effect[s] +
        coy_effect +
        true_params$beta_moci_ond_pup * moci_ond[t] +
        true_params$beta_moci_amj_pup * moci_amj[t] +
        has_eseal[s] * true_params$beta_eseal_pup * elephant_seal[s, t] +
        dist_surv_effect
      phi_pup_M[s, t] <- plogis(logit_phi_pup_M)
      
      logit_phi_juv_F <- qlogis(true_params$phi_juv_F_base) + site_effect[s] * 0.5 +
        true_params$beta_moci_jfm_juv * moci_jfm[t]
      phi_juv_F[s, t] <- plogis(logit_phi_juv_F)
      
      logit_phi_juv_M <- qlogis(true_params$phi_juv_M_base) + site_effect[s] * 0.5 +
        true_params$beta_moci_jfm_juv * moci_jfm[t]
      phi_juv_M[s, t] <- plogis(logit_phi_juv_M)
      
      logit_phi_adult_F <- qlogis(true_params$phi_adult_F_base) + site_effect[s] * 0.25 +
        true_params$beta_moci_jfm_adult * moci_jfm[t]
      phi_adult_F[s, t] <- plogis(logit_phi_adult_F)
      
      logit_phi_adult_M <- qlogis(true_params$phi_adult_M_base) + site_effect[s] * 0.25 +
        true_params$beta_moci_jfm_adult * moci_jfm[t]
      phi_adult_M[s, t] <- plogis(logit_phi_adult_M)
      
      # Detection
      detect_breed[s, t] <- plogis(dist_detect_effect)
      detect_molt[s, t] <- plogis(dist_detect_effect + 
                                    true_params$beta_moci_amj_molt * moci_amj[t])
      
      if (t > 1) {
        expected_pup <- N_adult_F[s, t-1] * true_params$avg_fecundity
        
        new_juv_F <- N_pup[s, t-1] * true_params$prop_female * phi_pup_F[s, t]
        new_juv_M <- N_pup[s, t-1] * (1 - true_params$prop_female) * phi_pup_M[s, t]
        
        juv_stay_F <- N_juv_F[s, t-1] * phi_juv_F[s, t] * (2/3)
        juv_stay_M <- N_juv_M[s, t-1] * phi_juv_M[s, t] * (2/3)
        
        juv_to_adult_F <- N_juv_F[s, t-1] * phi_juv_F[s, t] * (1/3)
        juv_to_adult_M <- N_juv_M[s, t-1] * phi_juv_M[s, t] * (1/3)
        
        expected_juv_F <- new_juv_F + juv_stay_F
        expected_juv_M <- new_juv_M + juv_stay_M
        
        expected_adult_F <- N_adult_F[s, t-1] * phi_adult_F[s, t] + juv_to_adult_F
        expected_adult_M <- N_adult_M[s, t-1] * phi_adult_M[s, t] + juv_to_adult_M
        
        N_pup[s, t] <- exp(rnorm(1, log(max(expected_pup, 1)), true_params$sigma_process))
        N_juv_F[s, t] <- exp(rnorm(1, log(max(expected_juv_F, 0.1)), true_params$sigma_process * 0.5))
        N_juv_M[s, t] <- exp(rnorm(1, log(max(expected_juv_M, 0.1)), true_params$sigma_process * 0.5))
        N_adult_F[s, t] <- exp(rnorm(1, log(max(expected_adult_F, 1)), true_params$sigma_process * 0.5))
        N_adult_M[s, t] <- exp(rnorm(1, log(max(expected_adult_M, 1)), true_params$sigma_process * 0.5))
      }
    }
  }
  
  # -----------------------------------------
  # GENERATE OBSERVATIONS
  # -----------------------------------------
  
  N_adult_total <- N_adult_F + N_adult_M
  N_juv_total <- N_juv_F + N_juv_M
  N_molt_true <- N_juv_total + N_adult_total
  
  y_adult <- y_pup <- y_molt <- matrix(NA, S, T)
  
  for (s in 1:S) {
    for (t in 1:T) {
      y_adult[s, t] <- log(N_adult_total[s, t] * detect_breed[s, t]) + 
        rnorm(1, 0, true_params$sigma_obs_adult)
      y_pup[s, t] <- log(N_pup[s, t] * detect_breed[s, t]) + 
        rnorm(1, 0, true_params$sigma_obs_pup)
      y_molt[s, t] <- log(N_molt_true[s, t] * detect_molt[s, t]) + 
        rnorm(1, 0, true_params$sigma_obs_molt)
    }
  }
  
  # Add missing data
  n_obs <- S * T
  missing_adult <- sample(1:n_obs, size = round(0.05 * n_obs))
  missing_pup <- sample(1:n_obs, size = round(0.05 * n_obs))
  missing_molt <- sample(1:n_obs, size = round(0.05 * n_obs))
  
  y_adult[missing_adult] <- NA
  y_pup[missing_pup] <- NA
  y_molt[missing_molt] <- NA
  
  y_adult_obs <- ifelse(is.na(y_adult), 0, 1)
  y_pup_obs <- ifelse(is.na(y_pup), 0, 1)
  y_molt_obs <- ifelse(is.na(y_molt), 0, 1)
  
  y_adult[is.na(y_adult)] <- 0
  y_pup[is.na(y_pup)] <- 0
  y_molt[is.na(y_molt)] <- 0
  
  # -----------------------------------------
  # PROJECTION SCENARIOS
  # -----------------------------------------
  
  N_scenarios <- 4
  
  moci_proj <- matrix(0, N_scenarios, T_proj)
  coyote_proj <- matrix(0, N_scenarios, T_proj)
  
  moci_proj[1, ] <- 0                    # Status quo
  moci_proj[2, ] <- 1                    # Warm
  moci_proj[3, ] <- -1                   # Cool
  moci_proj[4, ] <- 1                    # Warm + high coyote
  
  coyote_proj[1, ] <- 0                  # Current level
  coyote_proj[2, ] <- 0                  # Current level
  coyote_proj[3, ] <- 0                  # Current level
  coyote_proj[4, ] <- 1                  # High coyote
  
  scenario_names <- c("Status Quo", "Warm (MOCI +1)", "Cool (MOCI -1)", "Warm + High Coyote")
  
  # -----------------------------------------
  # RETURN DATA
  # -----------------------------------------
  
  stan_data <- list(
    T = T,
    S = S,
    N_coy = N_coy,
    
    y_adult = y_adult,
    y_pup = y_pup,
    y_molt = y_molt,
    
    y_adult_obs = y_adult_obs,
    y_pup_obs = y_pup_obs,
    y_molt_obs = y_molt_obs,
    
    coyote = coyote,
    disturbance = disturbance,
    elephant_seal = elephant_seal,
    moci_jfm = as.vector(moci_jfm),
    moci_amj = as.vector(moci_amj),
    moci_ond = as.vector(moci_ond),
    
    coyote_idx = coyote_idx,
    has_eseal = has_eseal,
    
    T_proj = T_proj,
    N_scenarios = N_scenarios,
    moci_proj = moci_proj,
    coyote_proj = coyote_proj
  )
  
  return(list(
    stan_data = stan_data,
    true_params = true_params,
    true_states = list(
      N_adult_F = N_adult_F,
      N_adult_M = N_adult_M,
      N_juv_F = N_juv_F,
      N_juv_M = N_juv_M,
      N_pup = N_pup,
      N_adult_total = N_adult_total,
      N_molt_true = N_molt_true,
      phi_pup_F = phi_pup_F,
      phi_pup_M = phi_pup_M,
      phi_juv_F = phi_juv_F,
      phi_juv_M = phi_juv_M,
      phi_adult_F = phi_adult_F,
      phi_adult_M = phi_adult_M,
      detect_breed = detect_breed,
      detect_molt = detect_molt
    ),
    site_names = site_names,
    years = 1997:(1997 + T - 1),
    scenario_names = scenario_names
  ))
}


# ============================================================================
# PART 3: PREPARE REAL DATA FOR IPM
# ============================================================================

prepare_real_data_for_ipm_v3.1 <- function(dat, cov_t_scaled, years, T_proj = 10) {
  
  site_names <- c("BL", "DE", "DP", "PRH", "TB", "TP")
  S <- 6
  T <- length(years)
  N_coy <- 3
  
  # Extract counts
  adult_rows <- seq(1, 18, by = 3)
  molt_rows <- seq(2, 18, by = 3)
  pup_rows <- seq(3, 18, by = 3)
  
  y_adult <- as.matrix(dat[adult_rows, ])
  y_molt <- as.matrix(dat[molt_rows, ])
  y_pup <- as.matrix(dat[pup_rows, ])
  
  rownames(y_adult) <- site_names
  rownames(y_molt) <- site_names
  rownames(y_pup) <- site_names
  
  y_adult_obs <- ifelse(is.na(y_adult), 0, 1)
  y_molt_obs <- ifelse(is.na(y_molt), 0, 1)
  y_pup_obs <- ifelse(is.na(y_pup), 0, 1)
  
  y_adult[is.na(y_adult)] <- 0
  y_molt[is.na(y_molt)] <- 0
  y_pup[is.na(y_pup)] <- 0
  
  # Extract covariates from cov_t_scaled
  # Row mapping (based on your MARSS model):
  # 1-3: MOCI (JFM, AMJ, OND)
  # 4-9: Disturbance (BL, DE, DP, PRH, TB, TP)
  # 10-12: Coyote (BL, DE, DP)
  # 16: Elephant seal
  
  moci_jfm <- as.vector(cov_t_scaled[1, ])
  moci_amj <- as.vector(cov_t_scaled[2, ])
  moci_ond <- as.vector(cov_t_scaled[3, ])
  
  # Disturbance - SITE SPECIFIC (rows 4-9)
  disturbance <- matrix(0, S, T)
  disturbance[1, ] <- as.vector(cov_t_scaled[4, ])   # BL
  disturbance[2, ] <- as.vector(cov_t_scaled[5, ])   # DE
  disturbance[3, ] <- as.vector(cov_t_scaled[6, ])   # DP
  disturbance[4, ] <- as.vector(cov_t_scaled[7, ])   # PRH
  disturbance[5, ] <- as.vector(cov_t_scaled[8, ])   # TB
  disturbance[6, ] <- as.vector(cov_t_scaled[9, ])   # TP
  
  # Coyote - SITE SPECIFIC (rows 10-12, only BL/DE/DP)
  coyote <- matrix(0, S, T)
  coyote[1, ] <- as.vector(cov_t_scaled[10, ])   # BL
  coyote[2, ] <- as.vector(cov_t_scaled[11, ])   # DE
  coyote[3, ] <- as.vector(cov_t_scaled[12, ])   # DP
  # Sites 4-6 (PRH, TB, TP) stay at 0
  
  # Elephant seal - DE and PRH
  elephant_seal <- matrix(0, S, T)
  elephant_seal[2, ] <- as.vector(cov_t_scaled[16, ])   # DE
  elephant_seal[4, ] <- as.vector(cov_t_scaled[16, ])   # PRH (same data)
  
  # Site indicators
  coyote_idx <- c(1, 2, 3, 0, 0, 0)  # Maps site to coyote parameter (0 if none)
  has_eseal <- c(0, 1, 0, 1, 0, 0)   # DE and PRH
  
  # Projection scenarios
  N_scenarios <- 4
  
  moci_proj <- matrix(0, N_scenarios, T_proj)
  coyote_proj <- matrix(0, N_scenarios, T_proj)
  
  moci_proj[1, ] <- 0
  moci_proj[2, ] <- 1
  moci_proj[3, ] <- -1
  moci_proj[4, ] <- 1
  
  # Use recent coyote trend for projections (average across coyote sites)
  recent_coyote <- mean(c(
    mean(coyote[1, (T-4):T]),
    mean(coyote[2, (T-4):T]),
    mean(coyote[3, (T-4):T])
  ))
  coyote_proj[1, ] <- recent_coyote
  coyote_proj[2, ] <- recent_coyote
  coyote_proj[3, ] <- recent_coyote
  coyote_proj[4, ] <- recent_coyote + 1
  
  scenario_names <- c("Status Quo", "Warm (MOCI +1)", "Cool (MOCI -1)", "Warm + High Coyote")
  
  stan_data <- list(
    T = T,
    S = S,
    N_coy = N_coy,
    
    y_adult = y_adult,
    y_pup = y_pup,
    y_molt = y_molt,
    
    y_adult_obs = y_adult_obs,
    y_pup_obs = y_pup_obs,
    y_molt_obs = y_molt_obs,
    
    coyote = coyote,
    disturbance = disturbance,
    elephant_seal = elephant_seal,
    moci_jfm = moci_jfm,
    moci_amj = moci_amj,
    moci_ond = moci_ond,
    
    coyote_idx = coyote_idx,
    has_eseal = has_eseal,
    
    T_proj = T_proj,
    N_scenarios = N_scenarios,
    moci_proj = moci_proj,
    coyote_proj = coyote_proj
  )
  
  return(list(
    stan_data = stan_data,
    site_names = site_names,
    years = years,
    scenario_names = scenario_names,
    raw_counts = list(adult = y_adult, molt = y_molt, pup = y_pup)
  ))
}


# ============================================================================
# PART 4: MODEL DIAGNOSTICS
# ============================================================================

check_diagnostics_v3.1 <- function(fit) {
  
  cat("\n============================================\n")
  cat("MODEL DIAGNOSTICS\n")
  cat("============================================\n\n")
  
  fit$cmdstan_diagnose()
  
  params_to_check <- c(
    "phi_pup_F_base", "phi_juv_F_base", "phi_adult_F_base",
    "delta_pup", "delta_juv", "delta_adult",
    "fecund_primip", "fecund_young", "fecund_prime", "fecund_senior",
    "prop_female",
    # Site-specific coyote
    "beta_coy[1]", "beta_coy[2]", "beta_coy[3]",
    # Site-specific disturbance survival
    "beta_dist_surv[1]", "beta_dist_surv[2]", "beta_dist_surv[3]",
    "beta_dist_surv[4]", "beta_dist_surv[5]", "beta_dist_surv[6]",
    # Shared effects
    "beta_moci_ond_pup", "beta_moci_amj_pup", "beta_moci_jfm_juv", "beta_moci_jfm_adult",
    "beta_eseal_pup",
    "sigma_process", "sigma_obs_adult", "sigma_obs_pup", "sigma_obs_molt", "sigma_site"
  )
  
  diag_summary <- fit$summary(variables = params_to_check)
  
  cat("\nParameter Summary:\n")
  print(diag_summary %>% select(variable, mean, sd, q5, q95, rhat, ess_bulk))
  
  bad_rhat <- diag_summary %>% filter(rhat > 1.05)
  if (nrow(bad_rhat) > 0) {
    cat("\nWARNING: Parameters with R-hat > 1.05:\n")
    print(bad_rhat %>% select(variable, rhat))
  }
  
  low_ess <- diag_summary %>% filter(ess_bulk < 400)
  if (nrow(low_ess) > 0) {
    cat("\nWARNING: Parameters with low ESS (<400):\n")
    print(low_ess %>% select(variable, ess_bulk))
  }
  
  return(list(params = params_to_check, summary = diag_summary))
}


# ============================================================================
# PART 5: TRACE PLOTS
# ============================================================================

create_trace_plots_v3.1 <- function(fit, params, save = TRUE, prefix = "IPM_v3.1") {
  
  draws <- fit$draws(format = "df")
  
  # Survival
  p_trace_surv <- mcmc_trace(draws, 
                             pars = c("phi_pup_F_base", "phi_juv_F_base", "phi_adult_F_base",
                                      "delta_pup", "delta_juv", "delta_adult")) +
    labs(title = "Trace Plots: Survival Parameters")
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_trace_survival.jpeg"), p_trace_surv, 
                   width = 30, height = 20, units = "cm")
  
  # Fecundity
  p_trace_fec <- mcmc_trace(draws, 
                            pars = c("fecund_primip", "fecund_young", "fecund_prime", 
                                     "fecund_senior", "prop_female")) +
    labs(title = "Trace Plots: Fecundity Parameters")
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_trace_fecundity.jpeg"), p_trace_fec, 
                   width = 30, height = 20, units = "cm")
  
  # Site-specific coyote effects
  p_trace_coy <- mcmc_trace(draws, 
                            pars = c("beta_coy[1]", "beta_coy[2]", "beta_coy[3]")) +
    labs(title = "Trace Plots: Site-Specific Coyote Effects (BL, DE, DP)")
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_trace_coyote.jpeg"), p_trace_coy, 
                   width = 30, height = 15, units = "cm")
  
  # Site-specific disturbance effects
  p_trace_dist <- mcmc_trace(draws, 
                             pars = c("beta_dist_surv[1]", "beta_dist_surv[2]", "beta_dist_surv[3]",
                                      "beta_dist_surv[4]", "beta_dist_surv[5]", "beta_dist_surv[6]")) +
    labs(title = "Trace Plots: Site-Specific Disturbance Effects (BL, DE, DP, PRH, TB, TP)")
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_trace_disturbance.jpeg"), p_trace_dist, 
                   width = 30, height = 20, units = "cm")
  
  # MOCI effects
  p_trace_moci <- mcmc_trace(draws, 
                             pars = c("beta_moci_ond_pup", "beta_moci_amj_pup", 
                                      "beta_moci_jfm_juv", "beta_moci_jfm_adult",
                                      "beta_moci_amj_molt")) +
    labs(title = "Trace Plots: MOCI Effects")
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_trace_moci.jpeg"), p_trace_moci, 
                   width = 30, height = 20, units = "cm")
  
  # Error terms
  p_trace_err <- mcmc_trace(draws, 
                            pars = c("sigma_process", "sigma_obs_adult", 
                                     "sigma_obs_pup", "sigma_obs_molt", "sigma_site")) +
    labs(title = "Trace Plots: Error Terms")
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_trace_errors.jpeg"), p_trace_err, 
                   width = 30, height = 20, units = "cm")
  
  return(list(survival = p_trace_surv, fecundity = p_trace_fec, 
              coyote = p_trace_coy, disturbance = p_trace_dist,
              moci = p_trace_moci, errors = p_trace_err))
}


# ============================================================================
# PART 6: PARAMETER RECOVERY (SIMULATED DATA)
# ============================================================================

check_parameter_recovery_v3.1 <- function(fit, sim_data, save = TRUE, prefix = "IPM_v3.1") {
  
  cat("\n============================================\n")
  cat("PARAMETER RECOVERY CHECK\n")
  cat("============================================\n\n")
  
  tp <- sim_data$true_params
  
  true_vals <- tibble(
    parameter = c(
      "phi_pup_F_base", "phi_juv_F_base", "phi_adult_F_base",
      "delta_pup", "delta_juv", "delta_adult",
      "fecund_primip", "fecund_young", "fecund_prime", "fecund_senior",
      "prop_female",
      # Site-specific coyote
      "beta_coy[1]", "beta_coy[2]", "beta_coy[3]",
      # Site-specific disturbance
      "beta_dist_surv[1]", "beta_dist_surv[2]", "beta_dist_surv[3]",
      "beta_dist_surv[4]", "beta_dist_surv[5]", "beta_dist_surv[6]",
      # Shared effects
      "beta_moci_ond_pup", "beta_moci_amj_pup", "beta_moci_jfm_juv", "beta_moci_jfm_adult",
      "beta_eseal_pup",
      "sigma_process", "sigma_obs_adult", "sigma_obs_pup", "sigma_obs_molt"
    ),
    true_value = c(
      tp$phi_pup_F_base, tp$phi_juv_F_base, tp$phi_adult_F_base,
      tp$delta_pup, tp$delta_juv, tp$delta_adult,
      tp$fecund_primip, tp$fecund_young, tp$fecund_prime, tp$fecund_senior,
      tp$prop_female,
      tp$beta_coy[1], tp$beta_coy[2], tp$beta_coy[3],
      tp$beta_dist_surv[1], tp$beta_dist_surv[2], tp$beta_dist_surv[3],
      tp$beta_dist_surv[4], tp$beta_dist_surv[5], tp$beta_dist_surv[6],
      tp$beta_moci_ond_pup, tp$beta_moci_amj_pup, tp$beta_moci_jfm_juv, tp$beta_moci_jfm_adult,
      tp$beta_eseal_pup,
      tp$sigma_process, tp$sigma_obs_adult, tp$sigma_obs_pup, tp$sigma_obs_molt
    )
  )
  
  param_recovery <- fit$summary(variables = true_vals$parameter) %>%
    left_join(true_vals, by = c("variable" = "parameter")) %>%
    mutate(
      recovered = true_value >= q5 & true_value <= q95,
      bias = mean - true_value,
      rel_bias_pct = (mean - true_value) / abs(true_value) * 100
    )
  
  cat("Parameter Recovery Results:\n")
  print(param_recovery %>% 
          select(variable, true_value, mean, q5, q95, recovered, rel_bias_pct) %>%
          mutate(across(where(is.numeric), ~round(.x, 3))))
  
  cat(sprintf("\nRecovery rate: %d/%d (%.1f%%) parameters within 90%% CI\n",
              sum(param_recovery$recovered),
              nrow(param_recovery),
              100 * mean(param_recovery$recovered)))
  
  p_recovery <- ggplot(param_recovery, aes(x = true_value, y = mean)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray50") +
    geom_pointrange(aes(ymin = q5, ymax = q95, color = recovered), size = 0.8) +
    geom_text(aes(label = variable), hjust = -0.1, vjust = -0.3, size = 2.5, 
              check_overlap = TRUE) +
    scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "red")) +
    labs(
      x = "True Value",
      y = "Estimated Value (90% CI)",
      title = "Parameter Recovery: IPM v3.1 (Site-Specific Effects)",
      color = "Within 90% CI"
    ) +
    theme_minimal(base_size = 14) +
    coord_equal()
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_parameter_recovery.jpeg"), p_recovery, 
                   width = 25, height = 22, units = "cm")
  
  return(list(table = param_recovery, plot = p_recovery))
}


# ============================================================================
# PART 7: POSTERIOR PREDICTIVE CHECKS
# ============================================================================

create_ppc_plots_v3.1 <- function(fit, sim_data, save = TRUE, prefix = "IPM_v3.1") {
  
  cat("\n============================================\n")
  cat("POSTERIOR PREDICTIVE CHECKS\n")
  cat("============================================\n\n")
  
  y_adult_rep <- fit$draws(variables = "y_adult_rep", format = "matrix")
  y_pup_rep <- fit$draws(variables = "y_pup_rep", format = "matrix")
  y_molt_rep <- fit$draws(variables = "y_molt_rep", format = "matrix")
  
  y_adult_obs <- as.vector(t(sim_data$stan_data$y_adult))
  y_pup_obs <- as.vector(t(sim_data$stan_data$y_pup))
  y_molt_obs <- as.vector(t(sim_data$stan_data$y_molt))
  
  adult_indicator <- as.vector(t(sim_data$stan_data$y_adult_obs)) == 1
  pup_indicator <- as.vector(t(sim_data$stan_data$y_pup_obs)) == 1
  molt_indicator <- as.vector(t(sim_data$stan_data$y_molt_obs)) == 1
  
  p_ppc_adult <- ppc_dens_overlay(y_adult_obs[adult_indicator], 
                                  y_adult_rep[1:100, adult_indicator]) +
    labs(title = "PPC: Adult Counts (log scale)")
  
  p_ppc_pup <- ppc_dens_overlay(y_pup_obs[pup_indicator], 
                                y_pup_rep[1:100, pup_indicator]) +
    labs(title = "PPC: Pup Counts (log scale)")
  
  p_ppc_molt <- ppc_dens_overlay(y_molt_obs[molt_indicator], 
                                 y_molt_rep[1:100, molt_indicator]) +
    labs(title = "PPC: Molt Counts (log scale)")
  
  p_ppc_combined <- p_ppc_adult / p_ppc_pup / p_ppc_molt +
    plot_annotation(title = "Posterior Predictive Checks: IPM v3.1")
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_ppc_density.jpeg"), p_ppc_combined, 
                   width = 25, height = 30, units = "cm")
  
  site_groups <- rep(sim_data$site_names, each = sim_data$stan_data$T)
  
  p_ppc_site <- ppc_stat_grouped(
    y = y_adult_obs[adult_indicator],
    yrep = y_adult_rep[1:100, adult_indicator],
    group = site_groups[adult_indicator],
    stat = "mean"
  ) + labs(title = "PPC: Mean Adult Count by Site")
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_ppc_by_site.jpeg"), p_ppc_site, 
                   width = 25, height = 15, units = "cm")
  
  return(list(density = p_ppc_combined, by_site = p_ppc_site))
}


# ============================================================================
# PART 8: SITE-BY-AGE-CLASS TIME SERIES PLOTS
# ============================================================================

create_site_age_timeseries_v3.1 <- function(fit, sim_data, save = TRUE, prefix = "IPM_v3.1") {
  
  cat("\n============================================\n")
  cat("SITE-BY-AGE-CLASS TIME SERIES\n")
  cat("============================================\n\n")
  
  years <- sim_data$years
  site_names <- sim_data$site_names
  S <- length(site_names)
  T <- length(years)
  
  N_pup_draws <- fit$draws(variables = "N_pup", format = "matrix")
  N_juv_total_draws <- fit$draws(variables = "N_juv_total", format = "matrix")
  N_adult_total_draws <- fit$draws(variables = "N_adult_total", format = "matrix")
  
  create_site_summary <- function(draws, var_name) {
    result <- list()
    for (s in 1:S) {
      for (t in 1:T) {
        col_name <- paste0(var_name, "[", s, ",", t, "]")
        if (col_name %in% colnames(draws)) {
          result[[length(result) + 1]] <- tibble(
            Site = site_names[s],
            Year = years[t],
            mean = mean(draws[, col_name]),
            lo = quantile(draws[, col_name], 0.05),
            hi = quantile(draws[, col_name], 0.95)
          )
        }
      }
    }
    bind_rows(result)
  }
  
  pup_summary <- create_site_summary(N_pup_draws, "N_pup") %>% mutate(Age_Class = "Pup")
  juv_summary <- create_site_summary(N_juv_total_draws, "N_juv_total") %>% mutate(Age_Class = "Juvenile")
  adult_summary <- create_site_summary(N_adult_total_draws, "N_adult_total") %>% mutate(Age_Class = "Adult")
  
  all_summary <- bind_rows(pup_summary, juv_summary, adult_summary) %>%
    mutate(Age_Class = factor(Age_Class, levels = c("Pup", "Juvenile", "Adult")))
  
  # Plot: All sites, faceted by age class
  p_by_age <- ggplot(all_summary, aes(x = Year, y = mean, color = Site, fill = Site)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    facet_wrap(~Age_Class, scales = "free_y", ncol = 1) +
    labs(
      x = "Year",
      y = "Population Size",
      title = "Population by Age Class Across Sites"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_age_class_timeseries.jpeg"), p_by_age, 
                   width = 30, height = 35, units = "cm")
  
  # Plot: Each site separately
  site_plots <- list()
  
  for (site in site_names) {
    site_data <- all_summary %>% filter(Site == site)
    
    p <- ggplot(site_data, aes(x = Year, y = mean, color = Age_Class, fill = Age_Class)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, color = NA) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c("Pup" = "green3", "Juvenile" = "orange", "Adult" = "blue")) +
      scale_fill_manual(values = c("Pup" = "green", "Juvenile" = "orange", "Adult" = "blue")) +
      labs(
        x = "Year",
        y = "Population Size",
        title = site,
        color = "Age Class",
        fill = "Age Class"
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")
    
    site_plots[[site]] <- p
  }
  
  p_sites_combined <- wrap_plots(site_plots, ncol = 2) +
    plot_annotation(
      title = "Population Dynamics by Site and Age Class",
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_site_age_timeseries.jpeg"), p_sites_combined, 
                   width = 35, height = 40, units = "cm")
  
  # Stacked area chart per site
  stacked_plots <- list()
  
  for (site in site_names) {
    site_data <- all_summary %>% 
      filter(Site == site) %>%
      select(Year, Age_Class, mean)
    
    p <- ggplot(site_data, aes(x = Year, y = mean, fill = Age_Class)) +
      geom_area(alpha = 0.8) +
      scale_fill_manual(values = c("Pup" = "green3", "Juvenile" = "orange", "Adult" = "blue")) +
      labs(
        x = "Year",
        y = "Population Size",
        title = site
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")
    
    stacked_plots[[site]] <- p
  }
  
  p_stacked <- wrap_plots(stacked_plots, ncol = 2) +
    plot_annotation(
      title = "Population Composition by Site (Stacked)",
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_site_stacked.jpeg"), p_stacked, 
                   width = 35, height = 40, units = "cm")
  
  return(list(by_age = p_by_age, by_site = p_sites_combined, stacked = p_stacked,
              data = all_summary))
}


# ============================================================================
# PART 9: TIME SERIES PLOTS (TOTALS)
# ============================================================================

create_timeseries_plots_v3.1 <- function(fit, sim_data, save = TRUE, prefix = "IPM_v3.1") {
  
  cat("\n============================================\n")
  cat("TIME SERIES PLOTS\n")
  cat("============================================\n\n")
  
  years <- sim_data$years
  site_names <- sim_data$site_names
  S <- length(site_names)
  T <- length(years)
  
  # Total population
  N_total_draws <- fit$draws(variables = "N_total_all", format = "df") %>%
    select(starts_with("N_total_all"))
  
  N_total_summary <- tibble(
    Year = years,
    mean = colMeans(N_total_draws),
    lo = apply(N_total_draws, 2, quantile, 0.05),
    hi = apply(N_total_draws, 2, quantile, 0.95)
  )
  
  p_total <- ggplot(N_total_summary, aes(x = Year)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3, fill = "blue") +
    geom_line(aes(y = mean), linewidth = 1.2, color = "blue") +
    labs(
      x = "Year",
      y = "Total Population",
      title = "Estimated Total Harbor Seal Population (All Sites)"
    ) +
    theme_minimal(base_size = 14)
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_total_population.jpeg"), p_total, 
                   width = 25, height = 15, units = "cm")
  
  # Lambda
  lambda_draws <- fit$draws(variables = "lambda", format = "matrix")
  
  lambda_means <- matrix(NA, nrow = nrow(lambda_draws), ncol = T-1)
  for (t in 1:(T-1)) {
    cols <- paste0("lambda[", 1:S, ",", t, "]")
    lambda_means[, t] <- rowMeans(lambda_draws[, cols])
  }
  
  lambda_summary <- tibble(
    Year = years[1:(T-1)],
    mean = colMeans(lambda_means),
    lo = apply(lambda_means, 2, quantile, 0.05),
    hi = apply(lambda_means, 2, quantile, 0.95)
  )
  
  p_lambda <- ggplot(lambda_summary, aes(x = Year)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3, fill = "darkgreen") +
    geom_line(aes(y = mean), linewidth = 1.2, color = "darkgreen") +
    geom_hline(yintercept = 1, linetype = 2, color = "red") +
    labs(
      x = "Year",
      y = "Population Growth Rate (λ)",
      title = "Mean Population Growth Rate Over Time"
    ) +
    theme_minimal(base_size = 14)
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_lambda.jpeg"), p_lambda, 
                   width = 25, height = 15, units = "cm")
  
  # Sex ratio
  sex_ratio_draws <- fit$draws(variables = "sex_ratio_adult", format = "matrix")
  
  sex_ratio_means <- matrix(NA, nrow = nrow(sex_ratio_draws), ncol = T)
  for (t in 1:T) {
    cols <- paste0("sex_ratio_adult[", 1:S, ",", t, "]")
    sex_ratio_means[, t] <- rowMeans(sex_ratio_draws[, cols])
  }
  
  sex_ratio_summary <- tibble(
    Year = years,
    mean = colMeans(sex_ratio_means),
    lo = apply(sex_ratio_means, 2, quantile, 0.05),
    hi = apply(sex_ratio_means, 2, quantile, 0.95)
  )
  
  p_sex_ratio <- ggplot(sex_ratio_summary, aes(x = Year)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3, fill = "purple") +
    geom_line(aes(y = mean), linewidth = 1.2, color = "purple") +
    geom_hline(yintercept = 0.5, linetype = 2, color = "gray50") +
    labs(
      x = "Year",
      y = "Proportion Female",
      title = "Adult Sex Ratio Over Time (Derived from Literature Priors)"
    ) +
    ylim(0.45, 0.65) +
    theme_minimal(base_size = 14)
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_sex_ratio.jpeg"), p_sex_ratio, 
                   width = 25, height = 15, units = "cm")
  
  return(list(total = p_total, lambda = p_lambda, sex_ratio = p_sex_ratio))
}


# ============================================================================
# PART 10: POPULATION PROJECTIONS
# ============================================================================

create_projection_plots_v3.1 <- function(fit, sim_data, save = TRUE, prefix = "IPM_v3.1") {
  
  cat("\n============================================\n")
  cat("10-YEAR POPULATION PROJECTIONS\n")
  cat("============================================\n\n")
  
  years <- sim_data$years
  T <- length(years)
  T_proj <- sim_data$stan_data$T_proj
  N_scenarios <- sim_data$stan_data$N_scenarios
  scenario_names <- sim_data$scenario_names
  
  proj_years <- (max(years) + 1):(max(years) + T_proj)
  
  proj_summaries <- list()
  
  for (scen in 1:N_scenarios) {
    draws_var <- paste0("N_total_all_proj[", scen, ",")
    
    proj_draws <- fit$draws(format = "matrix")
    cols <- grep(draws_var, colnames(proj_draws), fixed = TRUE)
    
    if (length(cols) > 0) {
      proj_mat <- proj_draws[, cols]
      
      proj_summaries[[scen]] <- tibble(
        Scenario = scenario_names[scen],
        Year = proj_years,
        mean = colMeans(proj_mat),
        lo = apply(proj_mat, 2, quantile, 0.10),
        hi = apply(proj_mat, 2, quantile, 0.90)
      )
    }
  }
  
  proj_df <- bind_rows(proj_summaries)
  
  N_total_draws <- fit$draws(variables = "N_total_all", format = "df") %>%
    select(starts_with("N_total_all"))
  
  hist_summary <- tibble(
    Scenario = "Historical",
    Year = years,
    mean = colMeans(N_total_draws),
    lo = apply(N_total_draws, 2, quantile, 0.10),
    hi = apply(N_total_draws, 2, quantile, 0.90)
  )
  
  all_data <- bind_rows(hist_summary, proj_df) %>%
    mutate(Period = ifelse(Scenario == "Historical", "Historical", "Projection"))
  
  p_proj <- ggplot() +
    geom_ribbon(data = filter(all_data, Period == "Historical"),
                aes(x = Year, ymin = lo, ymax = hi), 
                alpha = 0.3, fill = "gray50") +
    geom_line(data = filter(all_data, Period == "Historical"),
              aes(x = Year, y = mean), linewidth = 1.2, color = "black") +
    geom_ribbon(data = filter(all_data, Period == "Projection"),
                aes(x = Year, ymin = lo, ymax = hi, fill = Scenario), 
                alpha = 0.2) +
    geom_line(data = filter(all_data, Period == "Projection"),
              aes(x = Year, y = mean, color = Scenario), linewidth = 1.2) +
    geom_vline(xintercept = max(years), linetype = 2, color = "red") +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(
      x = "Year",
      y = "Total Population",
      title = "10-Year Population Projections Under Different Scenarios",
      subtitle = "Vertical line marks start of projections; bands show 80% CI"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_projections.jpeg"), p_proj, 
                   width = 30, height = 20, units = "cm")
  
  baseline_pop <- proj_df %>% 
    filter(Scenario == scenario_names[1]) %>% 
    pull(mean)
  
  pct_change_df <- proj_df %>%
    group_by(Scenario) %>%
    mutate(
      pct_change = (mean - baseline_pop) / baseline_pop * 100
    ) %>%
    filter(Scenario != scenario_names[1])
  
  p_pct_change <- ggplot(pct_change_df, aes(x = Year, y = pct_change, color = Scenario)) +
    geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
    geom_line(linewidth = 1.2) +
    scale_color_brewer(palette = "Set1") +
    labs(
      x = "Year",
      y = "% Change from Status Quo",
      title = "Population Change Relative to Status Quo Scenario"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_projections_pct_change.jpeg"), p_pct_change, 
                   width = 25, height = 15, units = "cm")
  
  lambda_proj_summaries <- list()
  
  for (scen in 1:N_scenarios) {
    draws_var <- paste0("lambda_proj[", scen, ",")
    
    proj_draws <- fit$draws(format = "matrix")
    cols <- grep(draws_var, colnames(proj_draws), fixed = TRUE)
    
    if (length(cols) > 0) {
      proj_mat <- proj_draws[, cols]
      
      lambda_proj_summaries[[scen]] <- tibble(
        Scenario = scenario_names[scen],
        Year = proj_years[1:(T_proj-1)],
        mean = colMeans(proj_mat),
        lo = apply(proj_mat, 2, quantile, 0.10),
        hi = apply(proj_mat, 2, quantile, 0.90)
      )
    }
  }
  
  lambda_proj_df <- bind_rows(lambda_proj_summaries)
  
  p_lambda_proj <- ggplot(lambda_proj_df, aes(x = Year, y = mean, color = Scenario)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = Scenario), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_hline(yintercept = 1, linetype = 2, color = "red") +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(
      x = "Year",
      y = "Population Growth Rate (λ)",
      title = "Projected Population Growth Rate by Scenario"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_projections_lambda.jpeg"), p_lambda_proj, 
                   width = 25, height = 15, units = "cm")
  
  final_year_summary <- proj_df %>%
    filter(Year == max(Year)) %>%
    mutate(
      pct_of_baseline = mean / baseline_pop[T_proj] * 100
    )
  
  cat("\nProjection Summary (Final Year):\n")
  print(final_year_summary %>% select(Scenario, mean, lo, hi, pct_of_baseline))
  
  return(list(
    projection = p_proj, 
    pct_change = p_pct_change,
    lambda = p_lambda_proj,
    data = proj_df,
    final_summary = final_year_summary
  ))
}


# ============================================================================
# PART 11: SITE-SPECIFIC COVARIATE EFFECTS PLOTS
# ============================================================================

create_effect_plots_v3.1 <- function(fit, save = TRUE, prefix = "IPM_v3.1") {
  
  cat("\n============================================\n")
  cat("SITE-SPECIFIC COVARIATE EFFECTS PLOTS\n")
  cat("============================================\n\n")
  
  site_names <- c("BL", "DE", "DP", "PRH", "TB", "TP")
  
  # Function to create survival effect plot
  create_survival_effect_plot <- function(fit, param_name, baseline_param, 
                                          x_label, title,
                                          x_range = seq(-2, 2, length.out = 100),
                                          y_limits = NULL) {
    
    draws <- fit$draws(variables = param_name, format = "df")[[param_name]]
    baseline_draws <- fit$draws(variables = baseline_param, format = "df")[[baseline_param]]
    
    n_draws <- min(500, length(draws))
    idx <- sample(1:length(draws), n_draws)
    
    result_list <- lapply(idx, function(i) {
      logit_base <- qlogis(baseline_draws[i])
      logit_survival <- logit_base + draws[i] * x_range
      survival <- plogis(logit_survival)
      tibble(x = x_range, survival = survival, draw = i)
    })
    
    result_df <- bind_rows(result_list)
    
    summary_df <- result_df %>%
      group_by(x) %>%
      summarise(
        mean = mean(survival),
        lo = quantile(survival, 0.05),
        hi = quantile(survival, 0.95),
        .groups = "drop"
      )
    
    mean_effect <- mean(draws)
    baseline_mean <- mean(baseline_draws)
    line_color <- ifelse(mean_effect < 0, "red3", "blue3")
    fill_color <- ifelse(mean_effect < 0, "red", "blue")
    
    p <- ggplot(summary_df, aes(x = x)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = fill_color) +
      geom_line(aes(y = mean), linewidth = 1.2, color = line_color) +
      geom_hline(yintercept = baseline_mean, linetype = 2, color = "gray50") +
      geom_vline(xintercept = 0, linetype = 2, color = "gray50") +
      labs(x = x_label, y = "Pup Survival Probability", title = title) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(size = 11, face = "bold"))
    
    # Apply y-limits if provided, otherwise use data-driven limits
    if (!is.null(y_limits)) {
      p <- p + coord_cartesian(ylim = y_limits)
    }
    
    # Return both plot and data range for calculating common limits
    return(list(plot = p, y_range = c(min(summary_df$lo), max(summary_df$hi))))
  }
  
  # -----------------------------------------
  # First pass: get y-ranges for all coyote plots
  # -----------------------------------------
  
  coy_results <- list(
    create_survival_effect_plot(fit, "beta_coy[1]", "phi_pup_F_base",
                                "Coyote Index (SD)", "Coyote → Pup Survival (BL)"),
    create_survival_effect_plot(fit, "beta_coy[2]", "phi_pup_F_base",
                                "Coyote Index (SD)", "Coyote → Pup Survival (DE)"),
    create_survival_effect_plot(fit, "beta_coy[3]", "phi_pup_F_base",
                                "Coyote Index (SD)", "Coyote → Pup Survival (DP)")
  )
  
  # Calculate common y-limits for coyote plots
  coy_y_min <- min(sapply(coy_results, function(x) x$y_range[1]))
  coy_y_max <- max(sapply(coy_results, function(x) x$y_range[2]))
  coy_y_limits <- c(max(0, coy_y_min - 0.03), min(1, coy_y_max + 0.03))
  
  # Recreate with common limits
  p_coy_bl <- create_survival_effect_plot(fit, "beta_coy[1]", "phi_pup_F_base",
                                          "Coyote Index (SD)", "Coyote → Pup Survival (BL)",
                                          y_limits = coy_y_limits)$plot
  p_coy_de <- create_survival_effect_plot(fit, "beta_coy[2]", "phi_pup_F_base",
                                          "Coyote Index (SD)", "Coyote → Pup Survival (DE)",
                                          y_limits = coy_y_limits)$plot
  p_coy_dp <- create_survival_effect_plot(fit, "beta_coy[3]", "phi_pup_F_base",
                                          "Coyote Index (SD)", "Coyote → Pup Survival (DP)",
                                          y_limits = coy_y_limits)$plot
  
  coyote_plots <- p_coy_bl + p_coy_de + p_coy_dp +
    plot_layout(ncol = 3) +
    plot_annotation(title = "Site-Specific Coyote Effects on Pup Survival",
                    theme = theme(plot.title = element_text(size = 14, face = "bold")))
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_effects_coyote.jpeg"), coyote_plots, 
                   width = 36, height = 12, units = "cm")
  
  # -----------------------------------------
  # Site-specific disturbance effects - get common y-limits
  # -----------------------------------------
  
  dist_results <- list()
  for (s in 1:6) {
    dist_results[[s]] <- create_survival_effect_plot(
      fit, 
      paste0("beta_dist_surv[", s, "]"), 
      "phi_pup_F_base",
      "Disturbance Index (SD)", 
      paste0("Disturbance → Pup Survival (", site_names[s], ")")
    )
  }
  
  # Calculate common y-limits for disturbance plots
  dist_y_min <- min(sapply(dist_results, function(x) x$y_range[1]))
  dist_y_max <- max(sapply(dist_results, function(x) x$y_range[2]))
  dist_y_limits <- c(max(0, dist_y_min - 0.03), min(1, dist_y_max + 0.03))
  
  # Recreate with common limits
  dist_plots <- list()
  for (s in 1:6) {
    dist_plots[[s]] <- create_survival_effect_plot(
      fit, 
      paste0("beta_dist_surv[", s, "]"), 
      "phi_pup_F_base",
      "Disturbance Index (SD)", 
      paste0("Disturbance → Pup Survival (", site_names[s], ")"),
      y_limits = dist_y_limits
    )$plot
  }
  
  disturbance_plots <- wrap_plots(dist_plots, ncol = 3) +
    plot_annotation(title = "Site-Specific Disturbance Effects on Pup Survival",
                    theme = theme(plot.title = element_text(size = 14, face = "bold")))
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_effects_disturbance.jpeg"), disturbance_plots, 
                   width = 36, height = 24, units = "cm")
  
  # -----------------------------------------
  # MOCI effects (shared) - get common y-limits
  # -----------------------------------------
  
  moci_results <- list(
    create_survival_effect_plot(fit, "beta_moci_ond_pup", "phi_pup_F_base",
                                "MOCI Fall (SD)", "MOCI Fall → Pup Survival"),
    create_survival_effect_plot(fit, "beta_moci_amj_pup", "phi_pup_F_base",
                                "MOCI Spring (SD)", "MOCI Spring → Pup Survival"),
    create_survival_effect_plot(fit, "beta_eseal_pup", "phi_pup_F_base",
                                "Elephant Seal Index (SD)", "Elephant Seal → Pup Survival")
  )
  
  # Calculate common y-limits for MOCI plots
  moci_y_min <- min(sapply(moci_results, function(x) x$y_range[1]))
  moci_y_max <- max(sapply(moci_results, function(x) x$y_range[2]))
  moci_y_limits <- c(max(0, moci_y_min - 0.03), min(1, moci_y_max + 0.03))
  
  # Recreate with common limits
  p_moci_ond <- create_survival_effect_plot(fit, "beta_moci_ond_pup", "phi_pup_F_base",
                                            "MOCI Fall (SD)", "MOCI Fall → Pup Survival",
                                            y_limits = moci_y_limits)$plot
  p_moci_amj <- create_survival_effect_plot(fit, "beta_moci_amj_pup", "phi_pup_F_base",
                                            "MOCI Spring (SD)", "MOCI Spring → Pup Survival",
                                            y_limits = moci_y_limits)$plot
  p_eseal <- create_survival_effect_plot(fit, "beta_eseal_pup", "phi_pup_F_base",
                                         "Elephant Seal Index (SD)", "Elephant Seal → Pup Survival",
                                         y_limits = moci_y_limits)$plot
  
  moci_plots <- p_moci_ond + p_moci_amj + p_eseal +
    plot_layout(ncol = 3) +
    plot_annotation(title = "Shared Covariate Effects on Pup Survival",
                    theme = theme(plot.title = element_text(size = 14, face = "bold")))
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_effects_moci.jpeg"), moci_plots, 
                   width = 36, height = 12, units = "cm")
  
  # Combined all effects
  all_effects <- coyote_plots / disturbance_plots / moci_plots +
    plot_annotation(
      title = "IPM v3.1: All Covariate Effects on Pup Survival",
      subtitle = "Site-specific coyote and disturbance effects; shared MOCI and elephant seal effects",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12)
      )
    )
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_effects_all.jpeg"), all_effects, 
                   width = 36, height = 48, units = "cm")
  
  return(list(coyote = coyote_plots, disturbance = disturbance_plots, 
              moci = moci_plots, all = all_effects))
}


# ============================================================================
# PART 12: SUMMARY TABLE
# ============================================================================

create_summary_table_v3.1 <- function(fit, save = TRUE, prefix = "IPM_v3.1") {
  
  site_names <- c("BL", "DE", "DP", "PRH", "TB", "TP")
  
  params_to_check <- c(
    "phi_pup_F_base", "phi_juv_F_base", "phi_adult_F_base",
    "delta_pup", "delta_juv", "delta_adult",
    "fecund_primip", "fecund_young", "fecund_prime", "fecund_senior",
    "prop_female", "avg_fecundity",
    "beta_coy[1]", "beta_coy[2]", "beta_coy[3]",
    "beta_dist_surv[1]", "beta_dist_surv[2]", "beta_dist_surv[3]",
    "beta_dist_surv[4]", "beta_dist_surv[5]", "beta_dist_surv[6]",
    "beta_moci_ond_pup", "beta_moci_amj_pup", "beta_moci_jfm_juv", "beta_moci_jfm_adult",
    "beta_eseal_pup", "beta_moci_amj_molt",
    "sigma_process", "sigma_obs_adult", "sigma_obs_pup", "sigma_obs_molt", "sigma_site"
  )
  
  param_summary <- fit$summary(variables = params_to_check)
  
  param_table <- param_summary %>%
    mutate(
      Parameter = case_when(
        variable == "phi_pup_F_base" ~ "Female Pup Survival (baseline)",
        variable == "phi_juv_F_base" ~ "Female Juvenile Survival (baseline)",
        variable == "phi_adult_F_base" ~ "Female Adult Survival (baseline)",
        variable == "delta_pup" ~ "Male-Female Δ Pup Survival",
        variable == "delta_juv" ~ "Male-Female Δ Juvenile Survival",
        variable == "delta_adult" ~ "Male-Female Δ Adult Survival",
        variable == "fecund_primip" ~ "Fecundity: Primiparous (4-5 yr)",
        variable == "fecund_young" ~ "Fecundity: Young Multiparous (6-7 yr)",
        variable == "fecund_prime" ~ "Fecundity: Prime Age (8-25 yr)",
        variable == "fecund_senior" ~ "Fecundity: Senescent (25+ yr)",
        variable == "prop_female" ~ "Proportion Female at Birth",
        variable == "avg_fecundity" ~ "Average Fecundity (weighted)",
        variable == "beta_coy[1]" ~ "Coyote → Pup Survival (BL)",
        variable == "beta_coy[2]" ~ "Coyote → Pup Survival (DE)",
        variable == "beta_coy[3]" ~ "Coyote → Pup Survival (DP)",
        variable == "beta_dist_surv[1]" ~ "Disturbance → Pup Survival (BL)",
        variable == "beta_dist_surv[2]" ~ "Disturbance → Pup Survival (DE)",
        variable == "beta_dist_surv[3]" ~ "Disturbance → Pup Survival (DP)",
        variable == "beta_dist_surv[4]" ~ "Disturbance → Pup Survival (PRH)",
        variable == "beta_dist_surv[5]" ~ "Disturbance → Pup Survival (TB)",
        variable == "beta_dist_surv[6]" ~ "Disturbance → Pup Survival (TP)",
        variable == "beta_moci_ond_pup" ~ "MOCI Fall → Pup Survival",
        variable == "beta_moci_amj_pup" ~ "MOCI Spring → Pup Survival",
        variable == "beta_moci_jfm_juv" ~ "MOCI Winter → Juvenile Survival",
        variable == "beta_moci_jfm_adult" ~ "MOCI Winter → Adult Survival",
        variable == "beta_eseal_pup" ~ "Elephant Seal → Pup Survival",
        variable == "beta_moci_amj_molt" ~ "MOCI Spring → Molt Detection",
        variable == "sigma_process" ~ "Process Error (σ)",
        variable == "sigma_obs_adult" ~ "Obs Error Adult (σ)",
        variable == "sigma_obs_pup" ~ "Obs Error Pup (σ)",
        variable == "sigma_obs_molt" ~ "Obs Error Molt (σ)",
        variable == "sigma_site" ~ "Site Random Effect (σ)",
        TRUE ~ variable
      ),
      Estimate = sprintf("%.3f (%.3f, %.3f)", mean, q5, q95),
      Category = case_when(
        str_detect(variable, "phi|delta") ~ "Survival",
        str_detect(variable, "fecund|prop|avg") ~ "Reproduction",
        str_detect(variable, "beta_coy") ~ "Coyote (Site-Specific)",
        str_detect(variable, "beta_dist") ~ "Disturbance (Site-Specific)",
        str_detect(variable, "beta_moci|beta_eseal") ~ "Shared Covariates",
        str_detect(variable, "sigma") ~ "Error Terms"
      )
    ) %>%
    select(Category, Parameter, Estimate, rhat, ess_bulk)
  
  cat("\n============================================\n")
  cat("PARAMETER SUMMARY TABLE\n")
  cat("============================================\n\n")
  
  print(param_table, n = nrow(param_table))
  
  if (save) {
    write_csv(param_table, paste0("Output/", prefix, "_parameter_summary.csv"))
    cat(sprintf("\nSaved to Output/%s_parameter_summary.csv\n", prefix))
  }
  
  return(param_table)
}


# ============================================================================
# PART 13: SAVE MODEL OUTPUT
# ============================================================================

save_model_output_v3.1 <- function(fit, sim_data = NULL, save_draws = FALSE, prefix = "IPM_v3.1") {
  
  cat("\n============================================\n")
  cat("SAVING MODEL OUTPUT\n")
  cat("============================================\n\n")
  
  fit$save_object(paste0("Output/harbor_seal_", prefix, "_fit.rds"))
  cat(sprintf("Fit object saved to Output/harbor_seal_%s_fit.rds\n", prefix))
  
  if (save_draws) {
    draws <- fit$draws(format = "df")
    write_csv(as_tibble(draws), paste0("Output/", prefix, "_posterior_draws.csv"))
    cat(sprintf("Posterior draws saved to Output/%s_posterior_draws.csv\n", prefix))
  }
  
  # Print summary
  cat("\n============================================\n")
  cat("HARBOR SEAL IPM v3.1 - FINAL SUMMARY\n")
  cat("============================================\n\n")
  
  surv_params <- c("phi_pup_F_base", "phi_juv_F_base", "phi_adult_F_base")
  surv_summary <- fit$summary(variables = surv_params)
  
  cat("SURVIVAL RATES (Female baselines):\n")
  cat(sprintf("  Pup (0-1 yr): %.3f (%.3f - %.3f)\n",
              surv_summary$mean[1], surv_summary$q5[1], surv_summary$q95[1]))
  cat(sprintf("  Juvenile (1-3 yr): %.3f (%.3f - %.3f)\n",
              surv_summary$mean[2], surv_summary$q5[2], surv_summary$q95[2]))
  cat(sprintf("  Adult (4+ yr): %.3f (%.3f - %.3f)\n",
              surv_summary$mean[3], surv_summary$q5[3], surv_summary$q95[3]))
  
  # Site-specific coyote effects
  coy_params <- c("beta_coy[1]", "beta_coy[2]", "beta_coy[3]")
  coy_summary <- fit$summary(variables = coy_params)
  
  cat("\nSITE-SPECIFIC COYOTE EFFECTS:\n")
  cat(sprintf("  BL: %.3f (%.3f - %.3f)\n",
              coy_summary$mean[1], coy_summary$q5[1], coy_summary$q95[1]))
  cat(sprintf("  DE: %.3f (%.3f - %.3f)\n",
              coy_summary$mean[2], coy_summary$q5[2], coy_summary$q95[2]))
  cat(sprintf("  DP: %.3f (%.3f - %.3f)\n",
              coy_summary$mean[3], coy_summary$q5[3], coy_summary$q95[3]))
  
  # Site-specific disturbance effects
  dist_params <- paste0("beta_dist_surv[", 1:6, "]")
  dist_summary <- fit$summary(variables = dist_params)
  site_names <- c("BL", "DE", "DP", "PRH", "TB", "TP")
  
  cat("\nSITE-SPECIFIC DISTURBANCE EFFECTS:\n")
  for (i in 1:6) {
    cat(sprintf("  %s: %.3f (%.3f - %.3f)\n",
                site_names[i],
                dist_summary$mean[i], dist_summary$q5[i], dist_summary$q95[i]))
  }
  
  cat("\n============================================\n")
}


# ============================================================================
# PART 14: MAIN EXECUTION FUNCTION
# ============================================================================

run_full_analysis_v3.1 <- function(use_real_data = FALSE, 
                                   dat = NULL, 
                                   cov_t_scaled = NULL,
                                   years = NULL,
                                   T_proj = 10,
                                   seed = 42,
                                   iter_warmup = 3000,
                                   iter_sampling = 1000,
                                   adapt_delta = 0.995,
                                   max_treedepth = 15) {
  
  cat("\n")
  cat("================================================================\n")
  cat("   HARBOR SEAL INTEGRATED POPULATION MODEL v3.1\n")
  cat("   Site-Specific Coyote & Disturbance Effects\n")
  cat("================================================================\n\n")
  
  # Set file prefix based on data type
  prefix <- ifelse(use_real_data, "IPM_v3.1_real", "IPM_v3.1_sim")
  cat(sprintf("Output file prefix: %s\n", prefix))
  
  # -----------------------------------------
  # Step 1: Prepare data
  # -----------------------------------------
  
  if (use_real_data) {
    cat("Using REAL DATA\n")
    data_list <- prepare_real_data_for_ipm_v3.1(dat, cov_t_scaled, years, T_proj)
    sim_data <- list(
      stan_data = data_list$stan_data,
      site_names = data_list$site_names,
      years = data_list$years,
      scenario_names = data_list$scenario_names,
      true_params = NULL
    )
  } else {
    cat("Using SIMULATED DATA\n")
    sim_data <- simulate_seal_ipm_data_v3.1(T = 29, S = 6, T_proj = T_proj, seed = seed)
  }
  
  # -----------------------------------------
  # Step 2: Compile and run model
  # -----------------------------------------
  
  cat("\nCompiling Stan model...\n")
  model <- cmdstan_model("harbor_seal_ipm_v3.1.stan")
  
  cat(sprintf("\nRunning MCMC (warmup=%d, sampling=%d, adapt_delta=%.2f)...\n",
              iter_warmup, iter_sampling, adapt_delta))
  
  fit <- model$sample(
    data = sim_data$stan_data,
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = 200,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
  
  # SAVE IMMEDIATELY
  fit$save_object(paste0("Output/harbor_seal_", prefix, "_fit.rds"))
  cat(sprintf("Fit saved to Output/harbor_seal_%s_fit.rds\n", prefix))
  
  # -----------------------------------------
  # Step 3: Diagnostics
  # -----------------------------------------
  
  diag_results <- check_diagnostics_v3.1(fit)
  
  # -----------------------------------------
  # Step 4: Trace plots
  # -----------------------------------------
  
  trace_plots <- create_trace_plots_v3.1(fit, diag_results$params, prefix = prefix)
  
  # -----------------------------------------
  # Step 5: Parameter recovery (simulated only)
  # -----------------------------------------
  
  recovery <- NULL
  if (!use_real_data && !is.null(sim_data$true_params)) {
    recovery <- check_parameter_recovery_v3.1(fit, sim_data, prefix = prefix)
  }
  
  # -----------------------------------------
  # Step 6: Posterior predictive checks
  # -----------------------------------------
  
  ppc_plots <- create_ppc_plots_v3.1(fit, sim_data, prefix = prefix)
  
  # -----------------------------------------
  # Step 7: Site-by-age time series
  # -----------------------------------------
  
  site_age_plots <- create_site_age_timeseries_v3.1(fit, sim_data, prefix = prefix)
  
  # -----------------------------------------
  # Step 8: Total time series
  # -----------------------------------------
  
  ts_plots <- create_timeseries_plots_v3.1(fit, sim_data, prefix = prefix)
  
  # -----------------------------------------
  # Step 9: Projections
  # -----------------------------------------
  
  proj_plots <- create_projection_plots_v3.1(fit, sim_data, prefix = prefix)
  
  # -----------------------------------------
  # Step 10: Effects plots
  # -----------------------------------------
  
  effects_plot <- create_effect_plots_v3.1(fit, prefix = prefix)
  
  # -----------------------------------------
  # Step 11: Summary table
  # -----------------------------------------
  
  summary_table <- create_summary_table_v3.1(fit, prefix = prefix)
  
  # -----------------------------------------
  # Step 12: Save output
  # -----------------------------------------
  
  save_model_output_v3.1(fit, sim_data, prefix = prefix)
  
  cat("\n================================================================\n")
  cat("   ANALYSIS COMPLETE\n")
  cat(sprintf("   All plots saved to Output/Plots/%s_*.jpeg\n", prefix))
  cat(sprintf("   Model fit saved to Output/harbor_seal_%s_fit.rds\n", prefix))
  cat(sprintf("   Parameter summary saved to Output/%s_parameter_summary.csv\n", prefix))
  cat("================================================================\n\n")
  
  return(list(
    fit = fit,
    model = model,
    data = sim_data,
    diagnostics = diag_results,
    summary = summary_table,
    recovery = recovery,
    projections = proj_plots,
    prefix = prefix
  ))
}


# ============================================================================
# PART 15: PORTFOLIO EFFECT ANALYSIS
# ============================================================================

#' Analyze Portfolio Effect - Site-level asynchrony in population dynamics
#' 
#' This demonstrates how different sites have good/poor years at different times,
#' buffering the overall metapopulation from collapse.

create_portfolio_analysis_v3.1 <- function(fit, sim_data, save = TRUE, prefix = "IPM_v3.1") {
  
  cat("\n============================================\n")
  cat("PORTFOLIO EFFECT ANALYSIS\n")
  cat("============================================\n\n")
  
  years <- sim_data$years
  site_names <- sim_data$site_names
  S <- length(site_names)
  T <- length(years)
  
  # -----------------------------------------
  # 1. Extract site-specific lambda (growth rates)
  # -----------------------------------------
  
  lambda_draws <- fit$draws(variables = "lambda", format = "matrix")
  
  # Create site x year matrix of mean lambda
  lambda_site_year <- matrix(NA, S, T-1)
  lambda_site_year_lo <- matrix(NA, S, T-1)
  lambda_site_year_hi <- matrix(NA, S, T-1)
  
  for (s in 1:S) {
    for (t in 1:(T-1)) {
      col_name <- paste0("lambda[", s, ",", t, "]")
      if (col_name %in% colnames(lambda_draws)) {
        lambda_site_year[s, t] <- mean(lambda_draws[, col_name])
        lambda_site_year_lo[s, t] <- quantile(lambda_draws[, col_name], 0.10)
        lambda_site_year_hi[s, t] <- quantile(lambda_draws[, col_name], 0.90)
      }
    }
  }
  
  rownames(lambda_site_year) <- site_names
  colnames(lambda_site_year) <- years[1:(T-1)]
  
  # -----------------------------------------
  # 2. Calculate portfolio metrics
  # -----------------------------------------
  
  # Mean lambda across sites per year
  mean_lambda_year <- colMeans(lambda_site_year, na.rm = TRUE)
  
  # Variance across sites per year (asynchrony measure)
  var_lambda_year <- apply(lambda_site_year, 2, var, na.rm = TRUE)
  
  # CV across sites per year
  cv_lambda_year <- sqrt(var_lambda_year) / mean_lambda_year
  
  # Correlation matrix between sites
  lambda_cor <- cor(t(lambda_site_year), use = "pairwise.complete.obs")
  
  # Portfolio effect statistic: CV of sum vs sum of CVs
  # Lower ratio = stronger portfolio effect
  cv_metapop <- sd(mean_lambda_year) / mean(mean_lambda_year)
  mean_cv_sites <- mean(apply(lambda_site_year, 1, function(x) sd(x, na.rm=TRUE) / mean(x, na.rm=TRUE)))
  portfolio_effect_ratio <- cv_metapop / mean_cv_sites
  
  cat(sprintf("Portfolio Effect Ratio: %.3f\n", portfolio_effect_ratio))
  cat("  (Values < 1 indicate portfolio buffering; lower = stronger effect)\n\n")
  
  cat("Site-to-site lambda correlations:\n")
  print(round(lambda_cor, 2))
  
  # -----------------------------------------
  # 3. Site-specific pup survival over time
  # -----------------------------------------
  
  phi_pup_F_draws <- fit$draws(variables = "phi_pup_F", format = "matrix")
  
  phi_pup_site_year <- matrix(NA, S, T)
  for (s in 1:S) {
    for (t in 1:T) {
      col_name <- paste0("phi_pup_F[", s, ",", t, "]")
      if (col_name %in% colnames(phi_pup_F_draws)) {
        phi_pup_site_year[s, t] <- mean(phi_pup_F_draws[, col_name])
      }
    }
  }
  rownames(phi_pup_site_year) <- site_names
  colnames(phi_pup_site_year) <- years
  
  # -----------------------------------------
  # 4. Create visualization: Lambda heatmap
  # -----------------------------------------
  
  lambda_df <- expand.grid(Site = site_names, Year = years[1:(T-1)]) %>%
    mutate(
      lambda = as.vector(t(lambda_site_year)),
      status = case_when(
        lambda > 1.05 ~ "Growing",
        lambda < 0.95 ~ "Declining",
        TRUE ~ "Stable"
      )
    )
  
  p_lambda_heatmap <- ggplot(lambda_df, aes(x = Year, y = Site, fill = lambda)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(
      low = "red3", mid = "white", high = "darkgreen",
      midpoint = 1, 
      limits = c(0.7, 1.3),
      oob = scales::squish,
      name = "λ"
    ) +
    geom_text(aes(label = sprintf("%.2f", lambda)), size = 2.5) +
    labs(
      x = "Year",
      y = "Site",
      title = "Site-Specific Population Growth Rate (λ) by Year",
      subtitle = "Green = growing (λ > 1), Red = declining (λ < 1)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid = element_blank()
    )
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_portfolio_lambda_heatmap.jpeg"), 
                   p_lambda_heatmap, width = 35, height = 15, units = "cm")
  
  # -----------------------------------------
  # 5. Pup survival heatmap
  # -----------------------------------------
  
  phi_df <- expand.grid(Site = site_names, Year = years) %>%
    mutate(phi_pup = as.vector(t(phi_pup_site_year)))
  
  p_phi_heatmap <- ggplot(phi_df, aes(x = Year, y = Site, fill = phi_pup)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(
      low = "red3", mid = "yellow", high = "darkgreen",
      midpoint = 0.7,
      limits = c(0.4, 0.95),
      oob = scales::squish,
      name = "Pup\nSurvival"
    ) +
    geom_text(aes(label = sprintf("%.2f", phi_pup)), size = 2.2) +
    labs(
      x = "Year",
      y = "Site",
      title = "Site-Specific Pup Survival by Year",
      subtitle = "Shows how coyote, disturbance, and MOCI create spatial variation"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid = element_blank()
    )
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_portfolio_pup_survival_heatmap.jpeg"), 
                   p_phi_heatmap, width = 35, height = 15, units = "cm")
  
  # -----------------------------------------
  # 6. Asynchrony plot: site lambdas over time
  # -----------------------------------------
  
  p_async <- ggplot(lambda_df, aes(x = Year, y = lambda, color = Site, group = Site)) +
    geom_hline(yintercept = 1, linetype = 2, color = "gray50") +
    geom_line(linewidth = 1, alpha = 0.7) +
    geom_point(size = 2) +
    scale_color_brewer(palette = "Set1") +
    labs(
      x = "Year",
      y = "Population Growth Rate (λ)",
      title = "Asynchronous Population Dynamics Across Sites",
      subtitle = "Different sites thrive in different years → portfolio buffering"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_portfolio_asynchrony.jpeg"), 
                   p_async, width = 30, height = 18, units = "cm")
  
  # -----------------------------------------
  # 7. Best/worst site per year
  # -----------------------------------------
  
  best_worst_df <- tibble(
    Year = years[1:(T-1)],
    Best_Site = site_names[apply(lambda_site_year, 2, which.max)],
    Best_Lambda = apply(lambda_site_year, 2, max, na.rm = TRUE),
    Worst_Site = site_names[apply(lambda_site_year, 2, which.min)],
    Worst_Lambda = apply(lambda_site_year, 2, min, na.rm = TRUE),
    Mean_Lambda = mean_lambda_year,
    Range = Best_Lambda - Worst_Lambda
  )
  
  cat("\nBest and Worst Performing Sites by Year:\n")
  print(best_worst_df, n = 30)
  
  # Plot showing which site is best/worst
  best_worst_long <- best_worst_df %>%
    select(Year, Best_Site, Worst_Site) %>%
    pivot_longer(cols = c(Best_Site, Worst_Site), 
                 names_to = "Performance", values_to = "Site") %>%
    mutate(Performance = ifelse(Performance == "Best_Site", "Best", "Worst"))
  
  p_best_worst <- ggplot(best_worst_long, aes(x = Year, y = Site, fill = Performance)) +
    geom_tile(color = "white", linewidth = 1, alpha = 0.8) +
    scale_fill_manual(values = c("Best" = "darkgreen", "Worst" = "red3")) +
    labs(
      x = "Year",
      y = "Site",
      title = "Best and Worst Performing Sites Each Year",
      subtitle = "No single site is consistently best or worst → portfolio effect"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid = element_blank()
    )
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_portfolio_best_worst.jpeg"), 
                   p_best_worst, width = 35, height = 12, units = "cm")
  
  # -----------------------------------------
  # 8. Correlation heatmap between sites
  # -----------------------------------------
  
  cor_df <- expand.grid(Site1 = site_names, Site2 = site_names) %>%
    mutate(correlation = as.vector(lambda_cor))
  
  p_cor <- ggplot(cor_df, aes(x = Site1, y = Site2, fill = correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", correlation)), size = 4) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Correlation"
    ) +
    labs(
      x = "", y = "",
      title = "Correlation in Population Growth Between Sites",
      subtitle = "Low/negative correlations indicate asynchrony (stronger portfolio effect)"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank()) +
    coord_fixed()
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_portfolio_correlation.jpeg"), 
                   p_cor, width = 20, height = 18, units = "cm")
  
  # -----------------------------------------
  # 9. Decompose variance: what drives site differences?
  # -----------------------------------------
  
  # Get covariate data
  coyote <- sim_data$stan_data$coyote
  disturbance <- sim_data$stan_data$disturbance
  moci_ond <- sim_data$stan_data$moci_ond
  
  # Get effect estimates
  beta_coy <- fit$summary(variables = c("beta_coy[1]", "beta_coy[2]", "beta_coy[3]"))$mean
  beta_dist <- fit$summary(variables = paste0("beta_dist_surv[", 1:6, "]"))$mean
  beta_moci <- fit$summary(variables = "beta_moci_ond_pup")$mean
  
  # Calculate contribution to survival variation for each site-year
  contribution_df <- expand.grid(Site = site_names, Year = years) %>%
    mutate(
      site_idx = match(Site, site_names),
      year_idx = match(Year, years)
    ) %>%
    rowwise() %>%
    mutate(
      coyote_effect = ifelse(site_idx <= 3, 
                             beta_coy[site_idx] * coyote[site_idx, year_idx], 
                             0),
      disturbance_effect = beta_dist[site_idx] * disturbance[site_idx, year_idx],
      moci_effect = beta_moci * moci_ond[year_idx]
    ) %>%
    ungroup() %>%
    mutate(
      total_effect = coyote_effect + disturbance_effect + moci_effect
    )
  
  # Pivot for stacked plot
  contrib_long <- contribution_df %>%
    select(Site, Year, coyote_effect, disturbance_effect, moci_effect) %>%
    pivot_longer(cols = c(coyote_effect, disturbance_effect, moci_effect),
                 names_to = "Driver", values_to = "Effect") %>%
    mutate(Driver = case_when(
      Driver == "coyote_effect" ~ "Coyote",
      Driver == "disturbance_effect" ~ "Disturbance",
      Driver == "moci_effect" ~ "MOCI"
    ))
  
  p_drivers <- ggplot(contrib_long, aes(x = Year, y = Effect, fill = Driver)) +
    geom_bar(stat = "identity", position = "stack", alpha = 0.8) +
    geom_hline(yintercept = 0, color = "black") +
    facet_wrap(~Site, ncol = 2) +
    scale_fill_manual(values = c("Coyote" = "brown", "Disturbance" = "orange", "MOCI" = "steelblue")) +
    labs(
      x = "Year",
      y = "Effect on Pup Survival (logit scale)",
      title = "Drivers of Spatial Variation in Pup Survival",
      subtitle = "Different drivers dominate at different sites and times"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      legend.position = "bottom"
    )
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_portfolio_drivers.jpeg"), 
                   p_drivers, width = 30, height = 35, units = "cm")
  
  # -----------------------------------------
  # 10. Summary statistics
  # -----------------------------------------
  
  portfolio_summary <- list(
    portfolio_effect_ratio = portfolio_effect_ratio,
    cv_metapopulation = cv_metapop,
    mean_cv_sites = mean_cv_sites,
    site_correlations = lambda_cor,
    mean_correlation = mean(lambda_cor[lower.tri(lambda_cor)]),
    best_worst_table = best_worst_df,
    times_best = table(best_worst_df$Best_Site),
    times_worst = table(best_worst_df$Worst_Site)
  )
  
  cat("\n--- PORTFOLIO EFFECT SUMMARY ---\n")
  cat(sprintf("Portfolio Effect Ratio: %.3f (< 1 = buffering)\n", portfolio_effect_ratio))
  cat(sprintf("Mean site-to-site correlation: %.3f (lower = more asynchrony)\n", 
              portfolio_summary$mean_correlation))
  cat("\nTimes each site was BEST performer:\n")
  print(portfolio_summary$times_best)
  cat("\nTimes each site was WORST performer:\n")
  print(portfolio_summary$times_worst)
  
  if (save) {
    write_csv(best_worst_df, paste0("Output/", prefix, "_portfolio_best_worst.csv"))
    write_csv(contribution_df, paste0("Output/", prefix, "_portfolio_drivers.csv"))
  }
  
  return(list(
    lambda_heatmap = p_lambda_heatmap,
    phi_heatmap = p_phi_heatmap,
    asynchrony = p_async,
    best_worst = p_best_worst,
    correlation = p_cor,
    drivers = p_drivers,
    lambda_matrix = lambda_site_year,
    phi_matrix = phi_pup_site_year,
    summary = portfolio_summary,
    contributions = contribution_df
  ))
}


# ============================================================================
# PART 16: SYNCHRONY SCENARIO PROJECTIONS
# ============================================================================

#' Compare population projections under asynchronous vs synchronous dynamics
#' 
#' Shows how portfolio effect buffers the metapopulation by comparing:
#' - Current model: sites experience different conditions (asynchronous)
#' - Hypothetical: all sites experience identical conditions (synchronous)

create_synchrony_projections_v3.1 <- function(fit, sim_data, 
                                              n_sims = 500,
                                              T_proj = 10,
                                              save = TRUE, 
                                              prefix = "IPM_v3.1") {
  
  cat("\n============================================\n")
  cat("SYNCHRONY VS ASYNCHRONY PROJECTIONS\n")
  cat("============================================\n\n")
  
  years <- sim_data$years
  site_names <- sim_data$site_names
  S <- length(site_names)
  T <- length(years)
  
  proj_years <- (max(years) + 1):(max(years) + T_proj)
  
  # -----------------------------------------
  # Extract posterior samples of key parameters
  # -----------------------------------------
  
  draws <- fit$draws(format = "df")
  n_draws <- nrow(draws)
  sample_idx <- sample(1:n_draws, min(n_sims, n_draws))
  
  # Get parameter draws
  phi_pup_F_base <- draws$phi_pup_F_base[sample_idx]
  phi_juv_F_base <- draws$phi_juv_F_base[sample_idx]
  phi_adult_F_base <- draws$phi_adult_F_base[sample_idx]
  
  delta_pup <- draws$delta_pup[sample_idx]
  delta_juv <- draws$delta_juv[sample_idx]
  delta_adult <- draws$delta_adult[sample_idx]
  
  prop_female <- draws$prop_female[sample_idx]
  avg_fecundity <- draws$avg_fecundity[sample_idx]
  
  # Site effects
  site_effect <- matrix(NA, length(sample_idx), S)
  for (s in 1:S) {
    site_effect[, s] <- draws[[paste0("site_effect[", s, "]")]][sample_idx]
  }
  
  # Coyote effects (site-specific)
  beta_coy <- matrix(NA, length(sample_idx), 3)
  for (k in 1:3) {
    beta_coy[, k] <- draws[[paste0("beta_coy[", k, "]")]][sample_idx]
  }
  
  # MOCI effects
  beta_moci_ond <- draws$beta_moci_ond_pup[sample_idx]
  beta_moci_amj <- draws$beta_moci_amj_pup[sample_idx]
  beta_moci_jfm_juv <- draws$beta_moci_jfm_juv[sample_idx]
  beta_moci_jfm_adult <- draws$beta_moci_jfm_adult[sample_idx]
  
  # Get final year populations for initialization
  N_adult_F_final <- matrix(NA, length(sample_idx), S)
  N_adult_M_final <- matrix(NA, length(sample_idx), S)
  N_juv_F_final <- matrix(NA, length(sample_idx), S)
  N_juv_M_final <- matrix(NA, length(sample_idx), S)
  N_pup_final <- matrix(NA, length(sample_idx), S)
  
  for (s in 1:S) {
    N_adult_F_final[, s] <- draws[[paste0("N_adult_F[", s, ",", T, "]")]][sample_idx]
    N_adult_M_final[, s] <- draws[[paste0("N_adult_M[", s, ",", T, "]")]][sample_idx]
    N_juv_F_final[, s] <- draws[[paste0("N_juv_F[", s, ",", T, "]")]][sample_idx]
    N_juv_M_final[, s] <- draws[[paste0("N_juv_M[", s, ",", T, "]")]][sample_idx]
    N_pup_final[, s] <- draws[[paste0("N_pup[", s, ",", T, "]")]][sample_idx]
  }
  
  # Coyote site index
  coyote_idx <- c(1, 2, 3, 0, 0, 0)
  
  # -----------------------------------------
  # Define scenarios
  # -----------------------------------------
  
  scenarios <- list(
    list(name = "Status Quo", moci = 0, coyote = 0),
    list(name = "Cool (MOCI -1)", moci = -1, coyote = 0),
    list(name = "Warm (MOCI +1)", moci = 1, coyote = 0),
    list(name = "Warm + High Coyote", moci = 1, coyote = 1)
  )
  
  # -----------------------------------------
  # Projection function
  # -----------------------------------------
  
  run_projection <- function(synchronous = FALSE, scenario, process_sd = 0.15) {
    
    moci_val <- scenario$moci
    coyote_val <- scenario$coyote
    
    # Storage for all simulations
    N_total_all <- matrix(NA, length(sample_idx), T_proj)
    
    for (i in 1:length(sample_idx)) {
      
      # Initialize
      N_AF <- N_adult_F_final[i, ]
      N_AM <- N_adult_M_final[i, ]
      N_JF <- N_juv_F_final[i, ]
      N_JM <- N_juv_M_final[i, ]
      N_P <- N_pup_final[i, ]
      
      N_total_all[i, 1] <- sum(N_AF + N_AM + N_JF + N_JM + N_P)
      
      for (tp in 2:T_proj) {
        
        if (synchronous) {
          # SYNCHRONOUS: All sites get the same random shock
          common_shock <- rnorm(1, 0, process_sd)
          site_shocks <- rep(common_shock, S)
        } else {
          # ASYNCHRONOUS: Each site gets independent shock
          site_shocks <- rnorm(S, 0, process_sd)
        }
        
        # Project each site
        for (s in 1:S) {
          
          # Coyote effect (site-specific)
          coy_eff <- 0
          if (coyote_idx[s] > 0) {
            coy_eff <- beta_coy[i, coyote_idx[s]] * coyote_val
          }
          
          # Survival rates
          logit_phi_pup_F <- qlogis(phi_pup_F_base[i]) + site_effect[i, s] +
            coy_eff + beta_moci_ond[i] * moci_val + beta_moci_amj[i] * moci_val
          phi_pup_F <- plogis(logit_phi_pup_F)
          
          logit_phi_pup_M <- qlogis(phi_pup_F_base[i] - delta_pup[i]) + site_effect[i, s] +
            coy_eff + beta_moci_ond[i] * moci_val + beta_moci_amj[i] * moci_val
          phi_pup_M <- plogis(logit_phi_pup_M)
          
          phi_juv_F <- plogis(qlogis(phi_juv_F_base[i]) + 0.5 * site_effect[i, s] +
                                beta_moci_jfm_juv[i] * moci_val)
          phi_juv_M <- plogis(qlogis(phi_juv_F_base[i] - delta_juv[i]) + 0.5 * site_effect[i, s] +
                                beta_moci_jfm_juv[i] * moci_val)
          
          phi_adult_F <- plogis(qlogis(phi_adult_F_base[i]) + 0.25 * site_effect[i, s] +
                                  beta_moci_jfm_adult[i] * moci_val)
          phi_adult_M <- plogis(qlogis(phi_adult_F_base[i] - delta_adult[i]) + 0.25 * site_effect[i, s] +
                                  beta_moci_jfm_adult[i] * moci_val)
          
          # Population dynamics with site-specific (or common) process error
          new_pups <- N_AF[s] * avg_fecundity[i] * exp(site_shocks[s])
          
          new_juv_F <- N_P[s] * prop_female[i] * phi_pup_F
          new_juv_M <- N_P[s] * (1 - prop_female[i]) * phi_pup_M
          
          juv_stay_F <- N_JF[s] * phi_juv_F * (2/3)
          juv_stay_M <- N_JM[s] * phi_juv_M * (2/3)
          
          juv_to_adult_F <- N_JF[s] * phi_juv_F * (1/3)
          juv_to_adult_M <- N_JM[s] * phi_juv_M * (1/3)
          
          # Update populations
          N_P[s] <- max(new_pups, 1)
          N_JF[s] <- max(new_juv_F + juv_stay_F, 0.1)
          N_JM[s] <- max(new_juv_M + juv_stay_M, 0.1)
          N_AF[s] <- max(N_AF[s] * phi_adult_F + juv_to_adult_F, 1)
          N_AM[s] <- max(N_AM[s] * phi_adult_M + juv_to_adult_M, 1)
        }
        
        N_total_all[i, tp] <- sum(N_AF + N_AM + N_JF + N_JM + N_P)
      }
    }
    
    return(N_total_all)
  }
  
  # -----------------------------------------
  # Run projections for all scenarios
  # -----------------------------------------
  
  results <- list()
  
  for (sc in scenarios) {
    cat(sprintf("Running scenario: %s\n", sc$name))
    
    # Asynchronous (current model behavior)
    async_proj <- run_projection(synchronous = FALSE, scenario = sc)
    
    # Synchronous (hypothetical - all sites correlated)
    sync_proj <- run_projection(synchronous = TRUE, scenario = sc)
    
    results[[sc$name]] <- list(
      async = async_proj,
      sync = sync_proj
    )
  }
  
  # -----------------------------------------
  # Create comparison data frame
  # -----------------------------------------
  
  comparison_df <- tibble()
  
  for (sc_name in names(results)) {
    
    async_proj <- results[[sc_name]]$async
    sync_proj <- results[[sc_name]]$sync
    
    for (tp in 1:T_proj) {
      comparison_df <- bind_rows(comparison_df, tibble(
        Scenario = sc_name,
        Year = proj_years[tp],
        Synchrony = "Asynchronous (Current)",
        mean = mean(async_proj[, tp]),
        lo = quantile(async_proj[, tp], 0.10),
        hi = quantile(async_proj[, tp], 0.90)
      ))
      
      comparison_df <- bind_rows(comparison_df, tibble(
        Scenario = sc_name,
        Year = proj_years[tp],
        Synchrony = "Synchronous (Hypothetical)",
        mean = mean(sync_proj[, tp]),
        lo = quantile(sync_proj[, tp], 0.10),
        hi = quantile(sync_proj[, tp], 0.90)
      ))
    }
  }
  
  # -----------------------------------------
  # Calculate CV ratio (portfolio effect metric)
  # -----------------------------------------
  
  cv_comparison <- tibble()
  
  for (sc_name in names(results)) {
    async_proj <- results[[sc_name]]$async
    sync_proj <- results[[sc_name]]$sync
    
    # CV of final year populations
    cv_async <- sd(async_proj[, T_proj]) / mean(async_proj[, T_proj])
    cv_sync <- sd(sync_proj[, T_proj]) / mean(sync_proj[, T_proj])
    
    cv_comparison <- bind_rows(cv_comparison, tibble(
      Scenario = sc_name,
      CV_Async = cv_async,
      CV_Sync = cv_sync,
      CV_Ratio = cv_async / cv_sync,
      Buffering_Pct = (1 - cv_async / cv_sync) * 100
    ))
  }
  
  cat("\n--- PORTFOLIO BUFFERING BY SCENARIO ---\n")
  print(cv_comparison)
  cat("\nCV_Ratio < 1 means asynchrony reduces variability\n")
  cat("Buffering_Pct shows % reduction in variability due to portfolio effect\n")
  
  # -----------------------------------------
  # Plot 1: Side-by-side comparison
  # -----------------------------------------
  
  p_comparison <- ggplot(comparison_df, aes(x = Year, y = mean, 
                                            color = Synchrony, 
                                            fill = Synchrony,
                                            linetype = Synchrony)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.2) +
    facet_wrap(~Scenario, ncol = 2, scales = "free_y") +
    scale_color_manual(values = c("Asynchronous (Current)" = "blue3", 
                                  "Synchronous (Hypothetical)" = "red3")) +
    scale_fill_manual(values = c("Asynchronous (Current)" = "blue", 
                                 "Synchronous (Hypothetical)" = "red")) +
    scale_linetype_manual(values = c("Asynchronous (Current)" = "solid",
                                     "Synchronous (Hypothetical)" = "dashed")) +
    labs(
      x = "Year",
      y = "Total Population",
      title = "Portfolio Effect: Asynchronous vs Synchronous Dynamics",
      subtitle = "Synchronous = all sites experience identical environmental conditions"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_synchrony_comparison.jpeg"), 
                   p_comparison, width = 30, height = 25, units = "cm")
  
  # -----------------------------------------
  # Plot 2: Uncertainty bands comparison
  # -----------------------------------------
  
  # Focus on one scenario to show the difference clearly
  status_quo_df <- comparison_df %>% filter(Scenario == "Status Quo")
  
  p_uncertainty <- ggplot(status_quo_df, aes(x = Year)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = Synchrony), alpha = 0.3) +
    geom_line(aes(y = mean, color = Synchrony), linewidth = 1.5) +
    scale_color_manual(values = c("Asynchronous (Current)" = "blue3", 
                                  "Synchronous (Hypothetical)" = "red3")) +
    scale_fill_manual(values = c("Asynchronous (Current)" = "blue", 
                                 "Synchronous (Hypothetical)" = "red")) +
    labs(
      x = "Year",
      y = "Total Population",
      title = "Portfolio Effect on Population Uncertainty (Status Quo Scenario)",
      subtitle = "Wider bands under synchrony = loss of buffering from spatial heterogeneity"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_synchrony_uncertainty.jpeg"), 
                   p_uncertainty, width = 28, height = 18, units = "cm")
  
  # -----------------------------------------
  # Plot 3: CV comparison bar chart
  # -----------------------------------------
  
  cv_long <- cv_comparison %>%
    select(Scenario, CV_Async, CV_Sync) %>%
    pivot_longer(cols = c(CV_Async, CV_Sync), 
                 names_to = "Type", values_to = "CV") %>%
    mutate(Type = ifelse(Type == "CV_Async", "Asynchronous", "Synchronous"))
  
  p_cv <- ggplot(cv_long, aes(x = Scenario, y = CV, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = c("Asynchronous" = "blue3", "Synchronous" = "red3")) +
    labs(
      x = "Scenario",
      y = "Coefficient of Variation",
      title = "Population Variability: Asynchronous vs Synchronous",
      subtitle = "Lower CV = more stable population; difference shows portfolio buffering"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1),
      legend.position = "bottom"
    )
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_synchrony_cv.jpeg"), 
                   p_cv, width = 25, height = 18, units = "cm")
  
  # -----------------------------------------
  # Plot 4: All scenarios overlaid
  # -----------------------------------------
  
  # Get historical data for context
  N_total_draws <- fit$draws(variables = "N_total_all", format = "df") %>%
    select(starts_with("N_total_all"))
  
  hist_summary <- tibble(
    Year = years,
    mean = colMeans(N_total_draws),
    lo = apply(N_total_draws, 2, quantile, 0.10),
    hi = apply(N_total_draws, 2, quantile, 0.90)
  )
  
  p_all_scenarios <- ggplot() +
    # Historical
    geom_ribbon(data = hist_summary, aes(x = Year, ymin = lo, ymax = hi), 
                alpha = 0.2, fill = "gray50") +
    geom_line(data = hist_summary, aes(x = Year, y = mean), 
              linewidth = 1.2, color = "black") +
    # Projections
    geom_ribbon(data = comparison_df, 
                aes(x = Year, ymin = lo, ymax = hi, fill = interaction(Scenario, Synchrony)), 
                alpha = 0.1) +
    geom_line(data = comparison_df, 
              aes(x = Year, y = mean, color = Scenario, linetype = Synchrony), 
              linewidth = 1) +
    geom_vline(xintercept = max(years), linetype = 2, color = "red") +
    scale_linetype_manual(values = c("Asynchronous (Current)" = "solid",
                                     "Synchronous (Hypothetical)" = "dashed")) +
    scale_color_brewer(palette = "Set1") +
    guides(fill = "none") +
    labs(
      x = "Year",
      y = "Total Population",
      title = "10-Year Projections: Asynchronous (solid) vs Synchronous (dashed)",
      subtitle = "Loss of portfolio effect increases uncertainty under all scenarios"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  if (save) ggsave(paste0("Output/Plots/", prefix, "_synchrony_all_scenarios.jpeg"), 
                   p_all_scenarios, width = 32, height = 20, units = "cm")
  
  # -----------------------------------------
  # Save data
  # -----------------------------------------
  
  if (save) {
    write_csv(comparison_df, paste0("Output/", prefix, "_synchrony_projections.csv"))
    write_csv(cv_comparison, paste0("Output/", prefix, "_synchrony_cv_comparison.csv"))
  }
  
  cat("\n--- SUMMARY ---\n")
  cat("Under synchronous dynamics (all sites correlated):\n")
  cat("  - Population uncertainty increases\n")
  cat("  - Extreme outcomes (crashes/booms) become more likely\n")
  cat("  - Portfolio buffering is lost\n")
  
  return(list(
    comparison = p_comparison,
    uncertainty = p_uncertainty,
    cv_plot = p_cv,
    all_scenarios = p_all_scenarios,
    projection_data = comparison_df,
    cv_comparison = cv_comparison,
    raw_projections = results
  ))
}


# ============================================================================
# USAGE INSTRUCTIONS
# ============================================================================

cat("\n")
cat("============================================\n")
cat("IPM v3.1 (SITE-SPECIFIC) LOADED SUCCESSFULLY\n")
cat("============================================\n")
cat("\nThis version matches MARSS structure with:\n")
cat("  - Site-specific coyote effects: beta_coy[1-3] for BL, DE, DP\n")
cat("  - Site-specific disturbance effects: beta_dist_surv[1-6] for all sites\n")
cat("  - Shared MOCI and elephant seal effects\n")
cat("\nFILE NAMING:\n")
cat("  Simulated data: IPM_v3.1_sim_*.jpeg, IPM_v3.1_sim_*.csv\n")
cat("  Real data:      IPM_v3.1_real_*.jpeg, IPM_v3.1_real_*.csv\n")
cat("\nTo run with SIMULATED data:\n")
cat("  library(bayesplot)  # IMPORTANT!\n")
cat("  results.sim <- run_full_analysis_v3.1(use_real_data = FALSE, seed = 42)\n")
cat("\nTo run with REAL data:\n")
cat("  library(bayesplot)  # IMPORTANT!\n")
cat("  results.real <- run_full_analysis_v3.1(\n")
cat("    use_real_data = TRUE,\n")
cat("    dat = dat,\n")
cat("    cov_t_scaled = cov_t_scaled,\n")
cat("    years = years,\n")
cat("    iter_warmup = 3000,   # More warmup for variance params\n")
cat("    iter_sampling = 1000,\n")
cat("    adapt_delta = 0.995,  # Higher for fewer divergences\n")
cat("    max_treedepth = 15\n")
cat("  )\n")
cat("\nPORTFOLIO EFFECT ANALYSIS:\n")
cat("  # After running the model:\n")
cat("  portfolio <- create_portfolio_analysis_v3.1(\n")
cat("    fit = results.real$fit,\n")
cat("    sim_data = results.real$data,\n")
cat("    prefix = 'IPM_v3.1_real'\n")
cat("  )\n")
cat("  \n")
cat("  # Returns:\n")
cat("  #   - Lambda heatmap (site x year growth rates)\n")
cat("  #   - Pup survival heatmap\n")
cat("  #   - Asynchrony plot (site trajectories)\n")
cat("  #   - Best/worst site per year\n")
cat("  #   - Site correlation matrix\n")
cat("  #   - Driver decomposition (coyote vs disturbance vs MOCI)\n")
cat("  #   - Portfolio effect ratio (< 1 = buffering)\n")
cat("\nSYNCHRONY COMPARISON (What if all sites were correlated?):\n")
cat("  sync_analysis <- create_synchrony_projections_v3.1(\n")
cat("    fit = results.real$fit,\n")
cat("    sim_data = results.real$data,\n")
cat("    n_sims = 500,\n")
cat("    prefix = 'IPM_v3.1_real'\n")
cat("  )\n")
cat("  \n")
cat("  # Returns:\n")
cat("  #   - Comparison plot: async vs sync under 4 scenarios\n")
cat("  #   - Uncertainty plot: shows wider CI under synchrony\n")
cat("  #   - CV comparison: quantifies buffering effect\n")
cat("  #   - All scenarios overlaid with historical data\n")
cat("\nSite indexing:\n")
cat("  1 = BL (Bolinas Lagoon)\n")
cat("  2 = DE (Double Point)\n")
cat("  3 = DP (Drakes Estero)\n")
cat("  4 = PRH (Point Reyes Headlands)\n")
cat("  5 = TB (Tomales Bay)\n")
cat("  6 = TP (Tomales Point)\n")
cat("============================================\n")


