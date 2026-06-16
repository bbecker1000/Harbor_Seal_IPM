#5A_ipm

#stan Code
# ============================================================================
# HARBOR SEAL IPM v3 - SEX-SPECIFIC SURVIVAL & AGE-DEPENDENT FECUNDITY
# ============================================================================

stan_code_v3 <- '
data {
  int<lower=1> T;                         // Number of years
  int<lower=1> S;                         // Number of sites
  int<lower=1> A_max;                     // Maximum age tracked (e.g., 30)
  
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
  real<lower=0, upper=1> phi_pup_F_base;      // First-year female survival
  real<lower=0, upper=1> phi_juv_F_base;      // Juvenile female survival (1-3 yr)
  real<lower=0, upper=1> phi_adult_F_base;    // Prime adult female survival (4-25 yr)
  real<lower=0, upper=1> phi_senior_F_base;   // Senior female survival (25+ yr)
  
  // Male survival difference (males have lower survival)
  real<lower=0> delta_pup;                    // ~0.10 lower for male pups
  real<lower=0> delta_juv;                    // ~0.08 lower for male juveniles
  real<lower=0> delta_adult;                  // ~0.05 lower for male adults
  
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
  // INITIAL AGE DISTRIBUTION (simplified)
  // ===========================================
  
  // Initial total females and males by broad age class
  vector<lower=0>[S] N_adult_F_init;
  vector<lower=0>[S] N_adult_M_init;
  vector<lower=0>[S] N_juv_F_init;
  vector<lower=0>[S] N_juv_M_init;
  vector<lower=0>[S] N_pup_init;
  
  // ===========================================
  // PROCESS ERRORS
  // ===========================================
  
  matrix[S, T-1] eps_adult_raw;
  matrix[S, T-1] eps_juv_raw;
  matrix[S, T-1] eps_pup_raw;
}

transformed parameters {
  // Site random effects
  vector[S] site_effect_pup = sigma_site * site_effect_pup_raw;
  vector[S] site_effect_juv = sigma_site * site_effect_juv_raw;
  vector[S] site_effect_adult = sigma_site * site_effect_adult_raw;
  
  // ===========================================
  // DERIVED SEX-SPECIFIC SURVIVAL
  // ===========================================
  
  real phi_pup_M_base = phi_pup_F_base - delta_pup;
  real phi_juv_M_base = phi_juv_F_base - delta_juv;
  real phi_adult_M_base = phi_adult_F_base - delta_adult;
  real phi_senior_M_base = phi_senior_F_base - delta_adult;  // Same delta for seniors
  
  // ===========================================
  // LATENT POPULATIONS (by sex and broad age class)
  // ===========================================
  
  matrix<lower=0>[S, T] N_adult_F;          // Breeding adult females (4+ years)
  matrix<lower=0>[S, T] N_adult_M;          // Adult males (4+ years)
  matrix<lower=0>[S, T] N_juv_F;            // Juvenile females (1-3 years)
  matrix<lower=0>[S, T] N_juv_M;            // Juvenile males (1-3 years)
  matrix<lower=0>[S, T] N_pup;              // Pups (both sexes)
  
  // Derived totals
  matrix<lower=0>[S, T] N_adult_total;
  matrix<lower=0>[S, T] N_molt_true;
  
  // Time-varying vital rates (sex-specific)
  matrix<lower=0, upper=1>[S, T] phi_pup_F;
  matrix<lower=0, upper=1>[S, T] phi_pup_M;
  matrix<lower=0, upper=1>[S, T] phi_juv_F;
  matrix<lower=0, upper=1>[S, T] phi_juv_M;
  matrix<lower=0, upper=1>[S, T] phi_adult_F;
  matrix<lower=0, upper=1>[S, T] phi_adult_M;
  matrix<lower=0, upper=1>[S, T] detect_breed;
  matrix<lower=0, upper=1>[S, T] detect_molt;
  
  // Average fecundity (weighted by age distribution - simplified)
  real avg_fecundity = 0.3 * fecund_primip + 0.2 * fecund_young + 
                       0.4 * fecund_prime + 0.1 * fecund_senior;
  
  // Calculate vital rates
  for (s in 1:S) {
    for (t in 1:T) {
      // -----------------------------------------
      // FEMALE SURVIVAL RATES
      // -----------------------------------------
      
      // Pup survival (female)
      real logit_phi_pup_F = logit(phi_pup_F_base) + site_effect_pup[s] +
        has_coyote[s] * beta_coy_pup * coyote[s, t] +
        beta_moci_ond_pup * moci_ond[t] +
        beta_moci_amj_pup * moci_amj[t] +
        has_eseal[s] * beta_eseal_pup * elephant_seal[s, t];
      phi_pup_F[s, t] = inv_logit(logit_phi_pup_F);
      
      // Juvenile survival (female)
      real logit_phi_juv_F = logit(phi_juv_F_base) + site_effect_juv[s] +
        has_coyote[s] * beta_coy_juv * coyote[s, t] +
        beta_moci_jfm_juv * moci_jfm[t];
      phi_juv_F[s, t] = inv_logit(logit_phi_juv_F);
      
      // Adult survival (female)
      real logit_phi_adult_F = logit(phi_adult_F_base) + site_effect_adult[s] +
        has_coyote[s] * beta_coy_adult * coyote[s, t] +
        beta_moci_jfm_adult * moci_jfm[t];
      phi_adult_F[s, t] = inv_logit(logit_phi_adult_F);
      
      // -----------------------------------------
      // MALE SURVIVAL RATES (lower than female)
      // -----------------------------------------
      
      // Apply same covariate effects but lower baseline
      real logit_phi_pup_M = logit(fmax(phi_pup_M_base, 0.01)) + site_effect_pup[s] +
        has_coyote[s] * beta_coy_pup * coyote[s, t] +
        beta_moci_ond_pup * moci_ond[t] +
        beta_moci_amj_pup * moci_amj[t] +
        has_eseal[s] * beta_eseal_pup * elephant_seal[s, t];
      phi_pup_M[s, t] = inv_logit(logit_phi_pup_M);
      
      real logit_phi_juv_M = logit(fmax(phi_juv_M_base, 0.01)) + site_effect_juv[s] +
        has_coyote[s] * beta_coy_juv * coyote[s, t] +
        beta_moci_jfm_juv * moci_jfm[t];
      phi_juv_M[s, t] = inv_logit(logit_phi_juv_M);
      
      real logit_phi_adult_M = logit(fmax(phi_adult_M_base, 0.01)) + site_effect_adult[s] +
        has_coyote[s] * beta_coy_adult * coyote[s, t] +
        beta_moci_jfm_adult * moci_jfm[t];
      phi_adult_M[s, t] = inv_logit(logit_phi_adult_M);
      
      // -----------------------------------------
      // DETECTION PROBABILITIES
      // -----------------------------------------
      
      detect_breed[s, t] = inv_logit(
        has_disturbance[s] * beta_dist_detect_breed * disturbance[s, t]
      );
      
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
    // Initial conditions
    N_adult_F[s, 1] = N_adult_F_init[s];
    N_adult_M[s, 1] = N_adult_M_init[s];
    N_juv_F[s, 1] = N_juv_F_init[s];
    N_juv_M[s, 1] = N_juv_M_init[s];
    N_pup[s, 1] = N_pup_init[s];
    
    N_adult_total[s, 1] = N_adult_F[s, 1] + N_adult_M[s, 1];
    N_molt_true[s, 1] = N_juv_F[s, 1] + N_juv_M[s, 1] + N_adult_total[s, 1];
    
    for (t in 2:T) {
      // -----------------------------------------
      // PUPS BORN THIS YEAR
      // Only adult females reproduce
      // -----------------------------------------
      real expected_pups = N_adult_F[s, t-1] * avg_fecundity;
      
      // -----------------------------------------
      // JUVENILE DYNAMICS (simplified: 3-year juvenile stage)
      // Pups → Juveniles, Juveniles → Adults
      // -----------------------------------------
      
      // New juveniles from surviving pups (split by sex)
      real new_juv_F = N_pup[s, t-1] * prop_female * phi_pup_F[s, t];
      real new_juv_M = N_pup[s, t-1] * (1 - prop_female) * phi_pup_M[s, t];
      
      // Juveniles that survive and stay juvenile (2/3 remain, 1/3 recruit to adult)
      real juv_stay_F = N_juv_F[s, t-1] * phi_juv_F[s, t] * (2.0/3.0);
      real juv_stay_M = N_juv_M[s, t-1] * phi_juv_M[s, t] * (2.0/3.0);
      
      // Juveniles that recruit to adult
      real juv_to_adult_F = N_juv_F[s, t-1] * phi_juv_F[s, t] * (1.0/3.0);
      real juv_to_adult_M = N_juv_M[s, t-1] * phi_juv_M[s, t] * (1.0/3.0);
      
      // Total juveniles
      real expected_juv_F = new_juv_F + juv_stay_F;
      real expected_juv_M = new_juv_M + juv_stay_M;
      
      // -----------------------------------------
      // ADULT DYNAMICS
      // Surviving adults + recruiting juveniles
      // -----------------------------------------
      real expected_adult_F = N_adult_F[s, t-1] * phi_adult_F[s, t] + juv_to_adult_F;
      real expected_adult_M = N_adult_M[s, t-1] * phi_adult_M[s, t] + juv_to_adult_M;
      
      // -----------------------------------------
      // ADD PROCESS ERROR
      // -----------------------------------------
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
      
      // Derived
      N_adult_total[s, t] = N_adult_F[s, t] + N_adult_M[s, t];
      N_molt_true[s, t] = N_juv_F[s, t] + N_juv_M[s, t] + N_adult_total[s, t];
    }
  }
}

model {
  // ===========================================
  // PRIORS - INFORMED BY LITERATURE
  // ===========================================
  
  // -----------------------------------------
  // Female survival (from Tugidak study)
  // -----------------------------------------
  phi_pup_F_base ~ beta(16, 4);             // Mean ~0.80, from literature
  phi_juv_F_base ~ beta(17, 3);             // Mean ~0.85, from literature  
  phi_adult_F_base ~ beta(18, 2);           // Mean ~0.90, from literature
  phi_senior_F_base ~ beta(16, 4);          // Mean ~0.80, senescence
  
  // Male survival difference (constrained to be positive)
  delta_pup ~ normal(0.10, 0.03);           // ~0.10 lower for males
  delta_juv ~ normal(0.08, 0.03);           // ~0.08 lower for males
  delta_adult ~ normal(0.05, 0.02);         // ~0.05 lower for males
  
  // -----------------------------------------
  // Age-specific fecundity (from literature)
  // -----------------------------------------
  fecund_primip ~ beta(14, 6);              // Mean ~0.70 for first-time breeders
  fecund_young ~ beta(17, 3);               // Mean ~0.85 for young multiparous
  fecund_prime ~ beta(19, 1);               // Mean ~0.95 for prime age
  fecund_senior ~ beta(16, 4);              // Mean ~0.80 senescent
  
  prop_female ~ beta(10, 10);               // Mean 0.50
  
  // -----------------------------------------
  // Covariate effects (informed by MARSS)
  // -----------------------------------------
  beta_coy_pup ~ normal(-0.5, 0.3);
  beta_coy_juv ~ normal(-0.2, 0.2);
  beta_coy_adult ~ normal(-0.05, 0.15);
  
  beta_moci_ond_pup ~ normal(-0.3, 0.2);
  beta_moci_amj_pup ~ normal(-0.2, 0.2);
  beta_moci_jfm_juv ~ normal(-0.15, 0.2);
  beta_moci_jfm_adult ~ normal(-0.1, 0.15);
  
  beta_eseal_pup ~ normal(0.15, 0.3);
  
  beta_dist_detect_breed ~ normal(-0.25, 0.2);
  beta_dist_detect_molt ~ normal(-0.25, 0.2);
  beta_moci_amj_molt ~ normal(0.1, 0.2);
  
  // -----------------------------------------
  // Random effects and error
  // -----------------------------------------
  sigma_site ~ exponential(5);
  site_effect_pup_raw ~ std_normal();
  site_effect_juv_raw ~ std_normal();
  site_effect_adult_raw ~ std_normal();
  
  sigma_process ~ exponential(10);
  sigma_obs_adult ~ exponential(5);
  sigma_obs_pup ~ exponential(5);
  sigma_obs_molt ~ exponential(5);
  
  // Initial populations
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
      // Adult counts (total adults, both sexes)
      if (y_adult_obs[s, t] == 1) {
        y_adult[s, t] ~ normal(
          log(N_adult_total[s, t] * detect_breed[s, t]),
          sigma_obs_adult
        );
      }
      
      // Pup counts
      if (y_pup_obs[s, t] == 1) {
        y_pup[s, t] ~ normal(
          log(N_pup[s, t] * detect_breed[s, t]),
          sigma_obs_pup
        );
      }
      
      // Molt counts (all non-pups)
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
  // Posterior predictive
  matrix[S, T] y_adult_rep;
  matrix[S, T] y_pup_rep;
  matrix[S, T] y_molt_rep;
  
  // Derived sex ratio over time
  matrix[S, T] sex_ratio_adult;  // Proportion female
  
  // Population growth rate
  matrix[S, T-1] lambda;
  
  // Total population
  vector[T] N_total_all;
  
  // Mean survival by sex
  vector[T] mean_phi_pup_F;
  vector[T] mean_phi_pup_M;
  vector[T] mean_phi_adult_F;
  vector[T] mean_phi_adult_M;
  
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
  
  for (s in 1:S) {
    for (t in 1:(T-1)) {
      real N_total_t = N_pup[s, t] + N_juv_F[s, t] + N_juv_M[s, t] + N_adult_total[s, t];
      real N_total_t1 = N_pup[s, t+1] + N_juv_F[s, t+1] + N_juv_M[s, t+1] + N_adult_total[s, t+1];
      lambda[s, t] = N_total_t1 / N_total_t;
    }
  }
  
  for (t in 1:T) {
    N_total_all[t] = 0;
    for (s in 1:S) {
      N_total_all[t] += N_pup[s, t] + N_juv_F[s, t] + N_juv_M[s, t] + N_adult_total[s, t];
    }
    mean_phi_pup_F[t] = mean(col(phi_pup_F, t));
    mean_phi_pup_M[t] = mean(col(phi_pup_M, t));
    mean_phi_adult_F[t] = mean(col(phi_adult_F, t));
    mean_phi_adult_M[t] = mean(col(phi_adult_M, t));
  }
}
'

# Write model
write_lines(stan_code_v3, "harbor_seal_ipm_v3.stan")