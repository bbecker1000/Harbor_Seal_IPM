
data {
  int<lower=1> T;                         // Number of years
  int<lower=1> S;                         // Number of sites
  
  // Observed log counts
  matrix[S, T] y_adult;                   // Breeding adult counts (spring)
  matrix[S, T] y_pup;                     // Pup counts (spring)
  matrix[S, T] y_molt;                    // Molt counts (late summer/fall)
  
  // Observation indicators (1 = observed, 0 = missing)
  array[S, T] int<lower=0, upper=1> y_adult_obs;
  array[S, T] int<lower=0, upper=1> y_pup_obs;
  array[S, T] int<lower=0, upper=1> y_molt_obs;
  
  // Covariates
  matrix[S, T] coyote;
  matrix[S, T] disturbance;
  matrix[S, T] elephant_seal;
  vector[T] moci_jfm;
  vector[T] moci_amj;
  vector[T] moci_ond;
  
  // Site indicators
  array[S] int<lower=0, upper=1> has_coyote;
  array[S] int<lower=0, upper=1> has_eseal;
  array[S] int<lower=0, upper=1> has_disturbance;
}

parameters {
  // ===========================================
  // SEX-SPECIFIC BASELINE SURVIVAL RATES
  // ===========================================
  
  // Female survival (higher than males)
  real<lower=0, upper=1> phi_pup_F_base;      // First-year female survival ~0.82
  real<lower=0, upper=1> phi_juv_F_base;      // Juvenile female survival (1-3 yr) ~0.87
  real<lower=0, upper=1> phi_adult_F_base;    // Prime adult female survival (4+ yr) ~0.93
  
  // Male survival difference (males have lower survival)
  real<lower=0, upper=0.3> delta_pup;         // ~0.10 lower for male pups
  real<lower=0, upper=0.2> delta_juv;         // ~0.08 lower for male juveniles
  real<lower=0, upper=0.15> delta_adult;      // ~0.05 lower for male adults
  
  // ===========================================
  // AGE-SPECIFIC REPRODUCTION
  // ===========================================
  
  // Fecundity varies with breeding experience
  real<lower=0, upper=1> fecund_primip;       // First-time breeders (age 4-5): ~0.70
  real<lower=0, upper=1> fecund_young;        // Young multiparous (age 6-7): ~0.85
  real<lower=0, upper=1> fecund_prime;        // Prime age (age 8-25): ~0.95
  real<lower=0, upper=1> fecund_senior;       // Senescent (age 25+): ~0.80
  
  real<lower=0, upper=1> prop_female;         // ~0.50
  
  // ===========================================
  // COVARIATE EFFECTS
  // ===========================================
  
  // Effects on pup survival (strongest effects)
  real beta_coy_pup;
  real beta_moci_ond_pup;
  real beta_moci_amj_pup;
  real beta_eseal_pup;
  
  // Effects on juvenile survival
  real beta_coy_juv;
  real beta_moci_jfm_juv;
  
  // Effects on adult survival
  real beta_coy_adult;
  real beta_moci_jfm_adult;
  
  // Effects on detection
  real beta_dist_detect_breed;
  real beta_dist_detect_molt;
  real beta_moci_amj_molt;
  
  // ===========================================
  // SITE RANDOM EFFECTS
  // ===========================================
  
  vector[S] site_effect_pup_raw;
  vector[S] site_effect_juv_raw;
  vector[S] site_effect_adult_raw;
  real<lower=0> sigma_site;
  
  // ===========================================
  // ERROR TERMS
  // ===========================================
  
  real<lower=0> sigma_process;
  real<lower=0> sigma_obs_adult;
  real<lower=0> sigma_obs_pup;
  real<lower=0> sigma_obs_molt;
  
  // ===========================================
  // INITIAL POPULATIONS
  // ===========================================
  
  vector<lower=0>[S] N_adult_F_init;
  vector<lower=0>[S] N_adult_M_init;
  vector<lower=0>[S] N_juv_F_init;
  vector<lower=0>[S] N_juv_M_init;
  vector<lower=0>[S] N_pup_init;
  
  // ===========================================
  // PROCESS ERRORS (non-centered)
  // ===========================================
  
  matrix[S, T-1] eps_adult_raw;
  matrix[S, T-1] eps_juv_raw;
  matrix[S, T-1] eps_pup_raw;
}

transformed parameters {
  // Site random effects (non-centered parameterization)
  vector[S] site_effect_pup = sigma_site * site_effect_pup_raw;
  vector[S] site_effect_juv = sigma_site * site_effect_juv_raw;
  vector[S] site_effect_adult = sigma_site * site_effect_adult_raw;
  
  // ===========================================
  // DERIVED SEX-SPECIFIC BASELINE SURVIVAL
  // ===========================================
  
  real phi_pup_M_base = fmax(phi_pup_F_base - delta_pup, 0.01);
  real phi_juv_M_base = fmax(phi_juv_F_base - delta_juv, 0.01);
  real phi_adult_M_base = fmax(phi_adult_F_base - delta_adult, 0.01);
  
  // ===========================================
  // LATENT POPULATIONS (by sex and age class)
  // ===========================================
  
  matrix<lower=0>[S, T] N_adult_F;          // Adult females (4+ years)
  matrix<lower=0>[S, T] N_adult_M;          // Adult males (4+ years)
  matrix<lower=0>[S, T] N_juv_F;            // Juvenile females (1-3 years)
  matrix<lower=0>[S, T] N_juv_M;            // Juvenile males (1-3 years)
  matrix<lower=0>[S, T] N_pup;              // Pups (both sexes)
  
  // Derived totals
  matrix<lower=0>[S, T] N_adult_total;
  matrix<lower=0>[S, T] N_juv_total;
  matrix<lower=0>[S, T] N_molt_true;
  matrix<lower=0>[S, T] N_total;
  
  // Time-varying vital rates (sex-specific)
  matrix<lower=0, upper=1>[S, T] phi_pup_F;
  matrix<lower=0, upper=1>[S, T] phi_pup_M;
  matrix<lower=0, upper=1>[S, T] phi_juv_F;
  matrix<lower=0, upper=1>[S, T] phi_juv_M;
  matrix<lower=0, upper=1>[S, T] phi_adult_F;
  matrix<lower=0, upper=1>[S, T] phi_adult_M;
  matrix<lower=0, upper=1>[S, T] detect_breed;
  matrix<lower=0, upper=1>[S, T] detect_molt;
  
  // Average fecundity (weighted by approximate stable age distribution)
  // Weights: primip ~30%, young ~20%, prime ~40%, senior ~10%
  real avg_fecundity = 0.30 * fecund_primip + 0.20 * fecund_young + 
                       0.40 * fecund_prime + 0.10 * fecund_senior;
  
  // ===========================================
  // CALCULATE TIME-VARYING VITAL RATES
  // ===========================================
  
  for (s in 1:S) {
    for (t in 1:T) {
      // -----------------------------------------
      // FEMALE SURVIVAL RATES
      // -----------------------------------------
      
      // Pup survival (female) - affected by coyotes, MOCI, elephant seals
      real logit_phi_pup_F = logit(phi_pup_F_base) + site_effect_pup[s] +
        has_coyote[s] * beta_coy_pup * coyote[s, t] +
        beta_moci_ond_pup * moci_ond[t] +
        beta_moci_amj_pup * moci_amj[t] +
        has_eseal[s] * beta_eseal_pup * elephant_seal[s, t];
      phi_pup_F[s, t] = inv_logit(logit_phi_pup_F);
      
      // Juvenile survival (female) - affected by coyotes, winter MOCI
      real logit_phi_juv_F = logit(phi_juv_F_base) + site_effect_juv[s] +
        has_coyote[s] * beta_coy_juv * coyote[s, t] +
        beta_moci_jfm_juv * moci_jfm[t];
      phi_juv_F[s, t] = inv_logit(logit_phi_juv_F);
      
      // Adult survival (female) - affected by coyotes (weak), winter MOCI
      real logit_phi_adult_F = logit(phi_adult_F_base) + site_effect_adult[s] +
        has_coyote[s] * beta_coy_adult * coyote[s, t] +
        beta_moci_jfm_adult * moci_jfm[t];
      phi_adult_F[s, t] = inv_logit(logit_phi_adult_F);
      
      // -----------------------------------------
      // MALE SURVIVAL RATES (lower baseline, same covariate effects)
      // -----------------------------------------
      
      real logit_phi_pup_M = logit(phi_pup_M_base) + site_effect_pup[s] +
        has_coyote[s] * beta_coy_pup * coyote[s, t] +
        beta_moci_ond_pup * moci_ond[t] +
        beta_moci_amj_pup * moci_amj[t] +
        has_eseal[s] * beta_eseal_pup * elephant_seal[s, t];
      phi_pup_M[s, t] = inv_logit(logit_phi_pup_M);
      
      real logit_phi_juv_M = logit(phi_juv_M_base) + site_effect_juv[s] +
        has_coyote[s] * beta_coy_juv * coyote[s, t] +
        beta_moci_jfm_juv * moci_jfm[t];
      phi_juv_M[s, t] = inv_logit(logit_phi_juv_M);
      
      real logit_phi_adult_M = logit(phi_adult_M_base) + site_effect_adult[s] +
        has_coyote[s] * beta_coy_adult * coyote[s, t] +
        beta_moci_jfm_adult * moci_jfm[t];
      phi_adult_M[s, t] = inv_logit(logit_phi_adult_M);
      
      // -----------------------------------------
      // DETECTION PROBABILITIES
      // -----------------------------------------
      
      // Breeding season detection (affected by disturbance)
      detect_breed[s, t] = inv_logit(
        has_disturbance[s] * beta_dist_detect_breed * disturbance[s, t]
      );
      
      // Molt season detection (affected by disturbance and spring MOCI)
      detect_molt[s, t] = inv_logit(
        has_disturbance[s] * beta_dist_detect_molt * disturbance[s, t] +
        beta_moci_amj_molt * moci_amj[t]
      );
    }
  }
  
  // ===========================================
  // POPULATION DYNAMICS
  // ===========================================
  
  for (s in 1:S) {
    // -----------------------------------------
    // INITIAL CONDITIONS (year 1)
    // -----------------------------------------
    N_adult_F[s, 1] = N_adult_F_init[s];
    N_adult_M[s, 1] = N_adult_M_init[s];
    N_juv_F[s, 1] = N_juv_F_init[s];
    N_juv_M[s, 1] = N_juv_M_init[s];
    N_pup[s, 1] = N_pup_init[s];
    
    N_adult_total[s, 1] = N_adult_F[s, 1] + N_adult_M[s, 1];
    N_juv_total[s, 1] = N_juv_F[s, 1] + N_juv_M[s, 1];
    N_molt_true[s, 1] = N_juv_total[s, 1] + N_adult_total[s, 1];
    N_total[s, 1] = N_pup[s, 1] + N_juv_total[s, 1] + N_adult_total[s, 1];
    
    // -----------------------------------------
    // DYNAMICS (years 2 to T)
    // -----------------------------------------
    for (t in 2:T) {
      // PUPS: Born to adult females
      real expected_pups = N_adult_F[s, t-1] * avg_fecundity;
      
      // JUVENILES: Surviving pups enter juvenile class
      // Simplified 3-year juvenile stage: 2/3 stay juvenile, 1/3 recruit to adult
      real new_juv_F = N_pup[s, t-1] * prop_female * phi_pup_F[s, t];
      real new_juv_M = N_pup[s, t-1] * (1 - prop_female) * phi_pup_M[s, t];
      
      real juv_stay_F = N_juv_F[s, t-1] * phi_juv_F[s, t] * (2.0/3.0);
      real juv_stay_M = N_juv_M[s, t-1] * phi_juv_M[s, t] * (2.0/3.0);
      
      real juv_to_adult_F = N_juv_F[s, t-1] * phi_juv_F[s, t] * (1.0/3.0);
      real juv_to_adult_M = N_juv_M[s, t-1] * phi_juv_M[s, t] * (1.0/3.0);
      
      real expected_juv_F = new_juv_F + juv_stay_F;
      real expected_juv_M = new_juv_M + juv_stay_M;
      
      // ADULTS: Surviving adults + recruiting juveniles
      real expected_adult_F = N_adult_F[s, t-1] * phi_adult_F[s, t] + juv_to_adult_F;
      real expected_adult_M = N_adult_M[s, t-1] * phi_adult_M[s, t] + juv_to_adult_M;
      
      // Add process error (log scale, non-centered)
      N_pup[s, t] = exp(log(fmax(expected_pups, 0.1)) + 
                        sigma_process * eps_pup_raw[s, t-1]);
      N_juv_F[s, t] = exp(log(fmax(expected_juv_F, 0.1)) + 
                          sigma_process * eps_juv_raw[s, t-1] * 0.5);
      N_juv_M[s, t] = exp(log(fmax(expected_juv_M, 0.1)) + 
                          sigma_process * eps_juv_raw[s, t-1] * 0.5);
      N_adult_F[s, t] = exp(log(fmax(expected_adult_F, 0.1)) + 
                            sigma_process * eps_adult_raw[s, t-1] * 0.5);
      N_adult_M[s, t] = exp(log(fmax(expected_adult_M, 0.1)) + 
                            sigma_process * eps_adult_raw[s, t-1] * 0.5);
      
      // Derived totals
      N_adult_total[s, t] = N_adult_F[s, t] + N_adult_M[s, t];
      N_juv_total[s, t] = N_juv_F[s, t] + N_juv_M[s, t];
      N_molt_true[s, t] = N_juv_total[s, t] + N_adult_total[s, t];
      N_total[s, t] = N_pup[s, t] + N_juv_total[s, t] + N_adult_total[s, t];
    }
  }
}

model {
  // ===========================================
  // PRIORS - INFORMED BY LITERATURE
  // ===========================================
  
  // -----------------------------------------
  // Female survival (from Tugidak study, Hastings et al. 2012)
  // -----------------------------------------
  phi_pup_F_base ~ beta(16, 4);             // Mean ~0.80
  phi_juv_F_base ~ beta(17, 3);             // Mean ~0.85
  phi_adult_F_base ~ beta(18, 2);           // Mean ~0.90
  
  // Male survival difference (constrained positive)
  delta_pup ~ normal(0.10, 0.03);           // Males ~0.10 lower
  delta_juv ~ normal(0.08, 0.03);           // Males ~0.08 lower
  delta_adult ~ normal(0.05, 0.02);         // Males ~0.05 lower
  
  // -----------------------------------------
  // Age-specific fecundity (from Bigg 1969, Thompson & Wheeler 2008)
  // -----------------------------------------
  fecund_primip ~ beta(14, 6);              // Mean ~0.70 for first-time breeders
  fecund_young ~ beta(17, 3);               // Mean ~0.85 for young multiparous
  fecund_prime ~ beta(19, 1);               // Mean ~0.95 for prime age
  fecund_senior ~ beta(16, 4);              // Mean ~0.80 senescent
  
  prop_female ~ beta(10, 10);               // Mean 0.50
  
  // -----------------------------------------
  // Covariate effects (informed by MARSS analysis)
  // -----------------------------------------
  
  // Pup survival effects (strongest)
  beta_coy_pup ~ normal(-0.5, 0.3);
  beta_moci_ond_pup ~ normal(-0.3, 0.2);
  beta_moci_amj_pup ~ normal(-0.2, 0.2);
  beta_eseal_pup ~ normal(0.15, 0.3);
  
  // Juvenile survival effects
  beta_coy_juv ~ normal(-0.2, 0.2);
  beta_moci_jfm_juv ~ normal(-0.15, 0.2);
  
  // Adult survival effects (weakest)
  beta_coy_adult ~ normal(-0.05, 0.15);
  beta_moci_jfm_adult ~ normal(-0.1, 0.15);
  
  // Detection effects
  beta_dist_detect_breed ~ normal(-0.25, 0.2);
  beta_dist_detect_molt ~ normal(-0.25, 0.2);
  beta_moci_amj_molt ~ normal(0.1, 0.2);
  
  // -----------------------------------------
  // Random effects and error terms
  // -----------------------------------------
  sigma_site ~ exponential(5);
  site_effect_pup_raw ~ std_normal();
  site_effect_juv_raw ~ std_normal();
  site_effect_adult_raw ~ std_normal();
  
  sigma_process ~ exponential(10);
  sigma_obs_adult ~ exponential(5);
  sigma_obs_pup ~ exponential(5);
  sigma_obs_molt ~ exponential(5);
  
  // Initial populations (weakly informative)
  N_adult_F_init ~ lognormal(5, 1);
  N_adult_M_init ~ lognormal(5, 1);
  N_juv_F_init ~ lognormal(4, 1);
  N_juv_M_init ~ lognormal(4, 1);
  N_pup_init ~ lognormal(4, 1);
  
  // Process errors
  to_vector(eps_adult_raw) ~ std_normal();
  to_vector(eps_juv_raw) ~ std_normal();
  to_vector(eps_pup_raw) ~ std_normal();
  
  // ===========================================
  // LIKELIHOOD
  // ===========================================
  
  for (s in 1:S) {
    for (t in 1:T) {
      // Adult counts (breeding season, both sexes)
      if (y_adult_obs[s, t] == 1) {
        y_adult[s, t] ~ normal(
          log(N_adult_total[s, t] * detect_breed[s, t]),
          sigma_obs_adult
        );
      }
      
      // Pup counts (breeding season)
      if (y_pup_obs[s, t] == 1) {
        y_pup[s, t] ~ normal(
          log(N_pup[s, t] * detect_breed[s, t]),
          sigma_obs_pup
        );
      }
      
      // Molt counts (all non-pups: juveniles + adults)
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
  
  // Sex ratio (proportion female among adults)
  matrix[S, T] sex_ratio_adult;
  
  // Population growth rate
  matrix[S, T-1] lambda;
  
  // Total population across all sites
  vector[T] N_total_all;
  
  // Mean vital rates across sites (by sex)
  vector[T] mean_phi_pup_F;
  vector[T] mean_phi_pup_M;
  vector[T] mean_phi_juv_F;
  vector[T] mean_phi_juv_M;
  vector[T] mean_phi_adult_F;
  vector[T] mean_phi_adult_M;
  
  // Pup:Adult ratio
  matrix[S, T] pup_adult_ratio;
  
  // Log-likelihood for LOO-CV
  matrix[S, T] log_lik_adult;
  matrix[S, T] log_lik_pup;
  matrix[S, T] log_lik_molt;
  
  // Generate replicated data and log-likelihood
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
      pup_adult_ratio[s, t] = N_pup[s, t] / N_adult_total[s, t];
      
      // Log-likelihood
      log_lik_adult[s, t] = normal_lpdf(y_adult[s, t] | 
                                         log(N_adult_total[s, t] * detect_breed[s, t]),
                                         sigma_obs_adult);
      log_lik_pup[s, t] = normal_lpdf(y_pup[s, t] | 
                                       log(N_pup[s, t] * detect_breed[s, t]),
                                       sigma_obs_pup);
      log_lik_molt[s, t] = normal_lpdf(y_molt[s, t] | 
                                        log(N_molt_true[s, t] * detect_molt[s, t]),
                                        sigma_obs_molt);
    }
  }
  
  // Population growth rate (lambda)
  for (s in 1:S) {
    for (t in 1:(T-1)) {
      lambda[s, t] = N_total[s, t+1] / N_total[s, t];
    }
  }
  
  // Total population across sites and mean vital rates
  for (t in 1:T) {
    N_total_all[t] = sum(col(N_total, t));
    mean_phi_pup_F[t] = mean(col(phi_pup_F, t));
    mean_phi_pup_M[t] = mean(col(phi_pup_M, t));
    mean_phi_juv_F[t] = mean(col(phi_juv_F, t));
    mean_phi_juv_M[t] = mean(col(phi_juv_M, t));
    mean_phi_adult_F[t] = mean(col(phi_adult_F, t));
    mean_phi_adult_M[t] = mean(col(phi_adult_M, t));
  }
}

