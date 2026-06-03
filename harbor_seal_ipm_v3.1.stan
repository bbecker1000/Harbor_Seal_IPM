
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

