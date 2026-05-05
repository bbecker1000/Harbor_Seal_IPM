# ============================================================================
# HARBOR SEAL INTEGRATED POPULATION MODEL v3.2
# ============================================================================
#
# CHANGES FROM v3.1:
#
#  (1) SEX-SPECIFIC SURVIVAL STRUCTURE CORRECTED
#        - Pup survival:      sex-neutral  (phi_pup_base, same M and F)
#        - Juvenile survival: sex-neutral  (phi_juv_base, same M and F)
#        - Adult survival:    sex-specific (phi_adult_F_base; male = F - delta_adult)
#        - Removed: delta_pup, delta_juv
#
#  (2) PUP SURVIVAL PRIOR CORRECTED
#        - v3.1 used beta(16,4) → mean 0.80 (too high; biased estimates ~0.68)
#        - v3.2 uses beta(5,5)  → mean 0.50 (consistent with field data)
#        - Parameterised on logit scale for flexibility
#
#  (3) OBSERVATION MODEL — EXPLICIT SEX STRUCTURE
#        (a) BREEDING ADULT COUNTS (spring):
#              Predominantly female; p_male_breed ~ beta(2,18) fraction of
#              males present. Likelihood uses N_adult_F + N_adult_M * p_male_breed
#        (b) MOLT COUNTS (summer):
#              Both sexes haul out equally. Likelihood uses
#              N_molt_true = N_juv_F + N_juv_M + N_adult_F + N_adult_M
#              (50:50 by construction; no additional parameter needed)
#
# ============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(patchwork)

dir.create("Output",       showWarnings = FALSE)
dir.create("Output/Plots", showWarnings = FALSE)

# ============================================================================
# SITE INDEXING
# ============================================================================
# 1=BL  2=DE  3=DP  4=PRH  5=TB  6=TP
# Coyote:        BL(1), DE(2), DP(3)
# Disturbance:   all 6
# Elephant seal: DE(2), PRH(4)
# ============================================================================


# ============================================================================
# PART 1: STAN MODEL (v3.2)
# ============================================================================

stan_code_v3.2 <- '
data {
  int<lower=1> T;
  int<lower=1> S;
  int<lower=1> N_coy;

  matrix[S, T] y_adult;
  matrix[S, T] y_pup;
  matrix[S, T] y_molt;

  array[S, T] int<lower=0, upper=1> y_adult_obs;
  array[S, T] int<lower=0, upper=1> y_pup_obs;
  array[S, T] int<lower=0, upper=1> y_molt_obs;

  matrix[S, T] coyote;
  matrix[S, T] disturbance;
  matrix[S, T] elephant_seal;

  vector[T] moci_jfm;
  vector[T] moci_amj;
  vector[T] moci_ond;

  array[S] int<lower=0, upper=N_coy> coyote_idx;
  array[S] int<lower=0, upper=1>    has_eseal;

  int<lower=0> T_proj;
  int<lower=1> N_scenarios;
  matrix[N_scenarios, T_proj] moci_proj;
  matrix[N_scenarios, T_proj] coyote_proj;
}

parameters {
  // ── Survival ──────────────────────────────────────────────────────────────
  //
  // Pups and juveniles: sex-neutral survival.
  //   Empirical evidence for harbour seals does not support sex-specific
  //   survival in these stages; field data centre pup-to-year-1 survival
  //   at ~0.50 (prior: logit-normal centred at 0 ≈ p = 0.50).
  //
  // Adults: females survive at a higher rate than males (delta_adult > 0).
  //
  real phi_pup_logit;           // logit-scale; ~ Normal(0, 0.4) → 95% CrI on prob ≈ (0.31, 0.69)
  real<lower=0, upper=1> phi_juv_base;
  real phi_adult_F_logit;       // logit-scale; ~ Normal(2.44, 0.25) → 95% CrI on prob ≈ (0.87, 0.96)
  real<lower=0, upper=0.10> delta_adult;   // F survival advantage (probability scale)

  // ── Reproduction ──────────────────────────────────────────────────────────
  real<lower=0, upper=1>    fecund_primip;   // age 4-5 yr, first breeders
  real<lower=0, upper=1>    fecund_mature;   // age 6+ yr, all experienced breeders
  real<lower=0.4, upper=0.6> prop_female;

  // ── Observation: male haul-out during breeding season ─────────────────────
  // Harbour seal mating is aquatic; most adult males remain in the water
  // during the spring pupping / breeding surveys.  Only a small fraction
  // (p_male_breed) haul out and are counted alongside females.
  // Prior: beta(2,18) → mean ≈ 0.10, 95th percentile ≈ 0.23.
  // Molt counts are 50:50 by construction (both sexes haul out equally
  // in summer) and require no additional parameter.
  real<lower=0, upper=0.30> p_male_breed;

  // ── Site-specific covariate effects ───────────────────────────────────────
  vector[N_coy] beta_coy;
  vector[S]     beta_dist_surv;
  vector[S]     beta_dist_detect;

  // ── Shared covariate effects ───────────────────────────────────────────────
  real beta_moci_ond_pup;
  real beta_moci_amj_pup;
  real beta_moci_jfm_juv;
  real beta_moci_jfm_adult;
  real beta_moci_amj_molt;
  real beta_eseal_pup;

  // ── Site random effects ────────────────────────────────────────────────────
  vector[S] site_effect_raw;
  real<lower=0.01, upper=0.5> sigma_site;

  // ── Error terms ────────────────────────────────────────────────────────────
  real<lower=0.05, upper=0.5>  sigma_process;
  real<lower=0.05, upper=0.4>  sigma_obs_adult;
  real<lower=0.02, upper=0.35> sigma_obs_pup;
  real<lower=0.05, upper=0.6>  sigma_obs_molt;

  // ── Initial populations (non-centred) ─────────────────────────────────────
  vector[S] log_N_adult_F_init_raw;
  vector[S] log_N_adult_M_init_raw;
  vector[S] log_N_juv_init_raw;
  vector[S] log_N_pup_init_raw;
  real mu_log_adult;
  real mu_log_juv;
  real mu_log_pup;
  real<lower=0> sigma_init;

  // ── Process errors (non-centred) ──────────────────────────────────────────
  matrix[S, T-1] eps_adult_raw;
  matrix[S, T-1] eps_juv_raw;
  matrix[S, T-1] eps_pup_raw;
}

transformed parameters {
  vector[S] site_effect = sigma_site * site_effect_raw;

  // Baseline survival probabilities
  real phi_pup_base     = inv_logit(phi_pup_logit);     // sex-neutral
  real phi_adult_F_base = inv_logit(phi_adult_F_logit); // derived; constrains posterior to realistic range
  real phi_adult_M_base = fmax(phi_adult_F_base - delta_adult, 0.01);

  // Initial populations
  vector<lower=0>[S] N_adult_F_init;
  vector<lower=0>[S] N_adult_M_init;
  vector<lower=0>[S] N_juv_F_init;
  vector<lower=0>[S] N_juv_M_init;
  vector<lower=0>[S] N_pup_init;

  for (s in 1:S) {
    N_adult_F_init[s] = exp(mu_log_adult + sigma_init * log_N_adult_F_init_raw[s]);
    N_adult_M_init[s] = exp(mu_log_adult + sigma_init * log_N_adult_M_init_raw[s]) * 0.9;
    N_juv_F_init[s]   = exp(mu_log_juv   + sigma_init * log_N_juv_init_raw[s]) * 0.5;
    N_juv_M_init[s]   = exp(mu_log_juv   + sigma_init * log_N_juv_init_raw[s]) * 0.5;
    N_pup_init[s]     = exp(mu_log_pup   + sigma_init * log_N_pup_init_raw[s]);
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
  // Pups and juveniles: identical for M and F (single matrix each)
  matrix<lower=0, upper=1>[S, T] phi_pup;      // sex-neutral pup survival
  matrix<lower=0, upper=1>[S, T] phi_juv;      // sex-neutral juv survival
  matrix<lower=0, upper=1>[S, T] phi_adult_F;  // adult female survival
  matrix<lower=0, upper=1>[S, T] phi_adult_M;  // adult male survival
  matrix<lower=0, upper=1>[S, T] detect_breed;
  matrix<lower=0, upper=1>[S, T] detect_molt;

  // ~20% of adult females are primiparous in a stable population
  real avg_fecundity = 0.20 * fecund_primip + 0.80 * fecund_mature;

  // ── Vital rates ────────────────────────────────────────────────────────────
  for (s in 1:S) {
    for (t in 1:T) {
      real coyote_effect = 0;
      if (coyote_idx[s] > 0)
        coyote_effect = beta_coy[coyote_idx[s]] * coyote[s, t];

      real dist_surv_eff   = beta_dist_surv[s]   * disturbance[s, t];
      real dist_detect_eff = beta_dist_detect[s] * disturbance[s, t];

      // Pup survival — sex-neutral
      phi_pup[s, t] = inv_logit(
        phi_pup_logit + site_effect[s] + coyote_effect +
        beta_moci_ond_pup * moci_ond[t] +
        beta_moci_amj_pup * moci_amj[t] +
        has_eseal[s] * beta_eseal_pup * elephant_seal[s, t] +
        dist_surv_eff);

      // Juvenile survival — sex-neutral
      phi_juv[s, t] = inv_logit(
        logit(phi_juv_base) + site_effect[s] * 0.5 +
        beta_moci_jfm_juv * moci_jfm[t]);

      // Adult survival — sex-specific (female > male)
      phi_adult_F[s, t] = inv_logit(
        phi_adult_F_logit + site_effect[s] * 0.25 +
        beta_moci_jfm_adult * moci_jfm[t]);

      phi_adult_M[s, t] = inv_logit(
        logit(phi_adult_M_base) + site_effect[s] * 0.25 +
        beta_moci_jfm_adult * moci_jfm[t]);

      // Detection
      detect_breed[s, t] = inv_logit(dist_detect_eff);
      detect_molt[s, t]  = inv_logit(dist_detect_eff +
                                      beta_moci_amj_molt * moci_amj[t]);
    }
  }

  // ── Population dynamics ────────────────────────────────────────────────────
  for (s in 1:S) {
    N_adult_F[s, 1] = N_adult_F_init[s];
    N_adult_M[s, 1] = N_adult_M_init[s];
    N_juv_F[s, 1]   = N_juv_F_init[s];
    N_juv_M[s, 1]   = N_juv_M_init[s];
    N_pup[s, 1]     = N_pup_init[s];

    N_adult_total[s, 1] = N_adult_F[s, 1] + N_adult_M[s, 1];
    N_juv_total[s, 1]   = N_juv_F[s, 1]   + N_juv_M[s, 1];
    N_molt_true[s, 1]   = N_juv_total[s, 1] + N_adult_total[s, 1];
    N_total[s, 1]       = N_pup[s, 1] + N_juv_total[s, 1] + N_adult_total[s, 1];

    for (t in 2:T) {
      // Pups produced by adult females only
      real expected_pups = N_adult_F[s, t-1] * avg_fecundity;

      // Pup → juvenile transition: same phi_pup for both sexes
      real new_juv_F = N_pup[s, t-1] * prop_female       * phi_pup[s, t];
      real new_juv_M = N_pup[s, t-1] * (1 - prop_female) * phi_pup[s, t];

      // Juvenile dynamics: same phi_juv for both sexes
      real juv_stay_F     = N_juv_F[s, t-1] * phi_juv[s, t] * (2.0/3.0);
      real juv_stay_M     = N_juv_M[s, t-1] * phi_juv[s, t] * (2.0/3.0);
      real juv_to_adult_F = N_juv_F[s, t-1] * phi_juv[s, t] * (1.0/3.0);
      real juv_to_adult_M = N_juv_M[s, t-1] * phi_juv[s, t] * (1.0/3.0);

      N_pup[s, t]     = exp(log(fmax(expected_pups, 1)) +
                            sigma_process * eps_pup_raw[s, t-1]);
      N_juv_F[s, t]   = exp(log(fmax(new_juv_F + juv_stay_F, 0.1)) +
                            sigma_process * eps_juv_raw[s, t-1] * 0.5);
      N_juv_M[s, t]   = exp(log(fmax(new_juv_M + juv_stay_M, 0.1)) +
                            sigma_process * eps_juv_raw[s, t-1] * 0.5);
      // Adult survival IS sex-specific
      N_adult_F[s, t] = exp(log(fmax(N_adult_F[s, t-1] * phi_adult_F[s, t] +
                                     juv_to_adult_F, 1)) +
                            sigma_process * eps_adult_raw[s, t-1] * 0.5);
      N_adult_M[s, t] = exp(log(fmax(N_adult_M[s, t-1] * phi_adult_M[s, t] +
                                     juv_to_adult_M, 1)) +
                            sigma_process * eps_adult_raw[s, t-1] * 0.5);

      N_adult_total[s, t] = N_adult_F[s, t] + N_adult_M[s, t];
      N_juv_total[s, t]   = N_juv_F[s, t]   + N_juv_M[s, t];
      // Molt: all adults + all juveniles — 50:50 M:F by construction
      N_molt_true[s, t]   = N_juv_total[s, t] + N_adult_total[s, t];
      N_total[s, t]       = N_pup[s, t] + N_juv_total[s, t] + N_adult_total[s, t];
    }
  }
}

model {
  // ── Priors ─────────────────────────────────────────────────────────────────

  // Pup survival: logit-normal centred at 0 → p ≈ 0.50
  // SD = 0.4 on logit scale → 95% CrI on prob scale ≈ (0.31, 0.69)
  // Tighter than v3.1 to prevent drift to implausibly low values
  phi_pup_logit    ~ normal(0.0, 0.4);

  // Juvenile survival: beta prior on probability scale
  phi_juv_base     ~ beta(17, 3);    // mean = 0.85
  // Adult female survival: logit-normal prior; prevents drift to implausible values
  // Normal(2.44, 0.25): logit(0.92)=2.44; 95% CrI on prob scale ≈ (0.87, 0.96)
  phi_adult_F_logit ~ normal(2.44, 0.25);
  delta_adult       ~ normal(0.05, 0.015);   // F ≈ 5% higher than M (probability scale)

  // Reproduction — 2-class fecundity (collapses prior 4-class structure)
  fecund_primip ~ beta(12, 8);    // mean = 0.60; primiparous, first breeders
  fecund_mature ~ beta(17, 3);    // mean = 0.85; all experienced females (6+ yr)
  prop_female   ~ beta(50, 50);

  // Male breeding haul-out fraction
  p_male_breed ~ beta(2, 18);

  // Site-specific coyote effects (informed by MARSS: BL weak, DE strong, DP near-zero)
  beta_coy[1] ~ normal(-0.15, 0.15);
  beta_coy[2] ~ normal(-0.40, 0.20);
  beta_coy[3] ~ normal(-0.05, 0.15);

  // Site-specific disturbance effects
  beta_dist_surv[1] ~ normal(-0.20, 0.15);
  beta_dist_surv[2] ~ normal(-0.30, 0.15);
  beta_dist_surv[3] ~ normal(-0.15, 0.15);
  beta_dist_surv[4] ~ normal(-0.10, 0.15);
  beta_dist_surv[5] ~ normal(-0.25, 0.15);
  beta_dist_surv[6] ~ normal(-0.15, 0.15);
  beta_dist_detect  ~ normal(-0.20, 0.15);

  // Shared covariate effects
  beta_moci_ond_pup   ~ normal(-0.25, 0.15);
  beta_moci_amj_pup   ~ normal(-0.15, 0.15);
  beta_moci_jfm_juv   ~ normal(-0.10, 0.15);
  beta_moci_jfm_adult ~ normal(-0.08, 0.10);
  beta_moci_amj_molt  ~ normal( 0.05, 0.15);
  beta_eseal_pup      ~ normal( 0.10, 0.20);

  // Random effects and errors
  sigma_site      ~ normal(0.2,  0.1);   // tighter; was absorbing too much variation at 0.776
  site_effect_raw ~ std_normal();
  sigma_process   ~ normal(0.15, 0.08);
  sigma_obs_adult ~ normal(0.18, 0.06);
  sigma_obs_pup   ~ normal(0.15, 0.02);  // tightened to resolve sigma_process ridge
  sigma_obs_molt  ~ normal(0.35, 0.10);

  // Initial populations
  mu_log_adult ~ normal(5, 0.5);
  mu_log_juv   ~ normal(4, 0.5);
  mu_log_pup   ~ normal(4, 0.5);
  sigma_init   ~ exponential(3);
  log_N_adult_F_init_raw ~ std_normal();
  log_N_adult_M_init_raw ~ std_normal();
  log_N_juv_init_raw     ~ std_normal();
  log_N_pup_init_raw     ~ std_normal();

  // Process errors
  to_vector(eps_adult_raw) ~ std_normal();
  to_vector(eps_juv_raw)   ~ std_normal();
  to_vector(eps_pup_raw)   ~ std_normal();

  // ── Likelihood ─────────────────────────────────────────────────────────────
  for (s in 1:S) {
    for (t in 1:T) {

      // (a) BREEDING ADULT COUNT — predominantly female.
      //     Spring haul-out surveys capture all nursing females but only a
      //     fraction p_male_breed of adult males (aquatic mating system).
      if (y_adult_obs[s, t] == 1) {
        real N_adult_obs = N_adult_F[s, t] + N_adult_M[s, t] * p_male_breed;
        y_adult[s, t] ~ normal(
          log(N_adult_obs * detect_breed[s, t]),
          sigma_obs_adult);
      }

      // (b) PUP COUNT — all pups present on land with nursing females.
      if (y_pup_obs[s, t] == 1) {
        y_pup[s, t] ~ normal(
          log(N_pup[s, t] * detect_breed[s, t]),
          sigma_obs_pup);
      }

      // (c) MOLT COUNT — both sexes haul out equally in summer (50:50).
      //     N_molt_true = N_juv_F + N_juv_M + N_adult_F + N_adult_M.
      //     No additional male-haul-out parameter needed here.
      if (y_molt_obs[s, t] == 1) {
        y_molt[s, t] ~ normal(
          log(N_molt_true[s, t] * detect_molt[s, t]),
          sigma_obs_molt);
      }

    }
  }
}

generated quantities {
  // ── Posterior predictive ───────────────────────────────────────────────────
  matrix[S, T] y_adult_rep;
  matrix[S, T] y_pup_rep;
  matrix[S, T] y_molt_rep;

  // ── Derived quantities ─────────────────────────────────────────────────────
  matrix[S, T]   sex_ratio_adult;      // true sex ratio in population
  matrix[S, T]   sex_ratio_observed;   // apparent sex ratio in spring counts
  matrix[S, T-1] lambda;
  vector[T]      N_total_all;

  // Mean vital rates across sites
  vector[T] mean_phi_pup;        // sex-neutral
  vector[T] mean_phi_juv;        // sex-neutral
  vector[T] mean_phi_adult_F;
  vector[T] mean_phi_adult_M;

  // ── Projections ────────────────────────────────────────────────────────────
  array[N_scenarios] matrix[S, T_proj] N_total_proj;
  array[N_scenarios] matrix[S, T_proj] N_pup_proj;
  array[N_scenarios] matrix[S, T_proj] N_adult_proj;
  array[N_scenarios] vector[T_proj]    N_total_all_proj;
  array[N_scenarios] vector[T_proj-1]  lambda_proj;

  for (s in 1:S) {
    for (t in 1:T) {
      real N_adult_obs_rep = N_adult_F[s, t] + N_adult_M[s, t] * p_male_breed;
      y_adult_rep[s, t] = normal_rng(
        log(N_adult_obs_rep * detect_breed[s, t]), sigma_obs_adult);
      y_pup_rep[s, t] = normal_rng(
        log(N_pup[s, t] * detect_breed[s, t]), sigma_obs_pup);
      y_molt_rep[s, t] = normal_rng(
        log(N_molt_true[s, t] * detect_molt[s, t]), sigma_obs_molt);

      sex_ratio_adult[s, t]    = N_adult_F[s, t] / N_adult_total[s, t];
      sex_ratio_observed[s, t] = N_adult_F[s, t] /
        (N_adult_F[s, t] + N_adult_M[s, t] * p_male_breed);
    }
  }

  for (s in 1:S) {
    for (t in 1:(T-1))
      lambda[s, t] = N_total[s, t+1] / N_total[s, t];
  }

  for (t in 1:T) {
    N_total_all[t]    = sum(col(N_total, t));
    mean_phi_pup[t]   = mean(col(phi_pup,     t));
    mean_phi_juv[t]   = mean(col(phi_juv,     t));
    mean_phi_adult_F[t] = mean(col(phi_adult_F, t));
    mean_phi_adult_M[t] = mean(col(phi_adult_M, t));
  }

  // ── Scenario projections ───────────────────────────────────────────────────
  for (scen in 1:N_scenarios) {
    matrix[S, T_proj] pAF; matrix[S, T_proj] pAM;
    matrix[S, T_proj] pJF; matrix[S, T_proj] pJM;
    matrix[S, T_proj] pP;

    for (s in 1:S) {
      pAF[s,1]=N_adult_F[s,T]; pAM[s,1]=N_adult_M[s,T];
      pJF[s,1]=N_juv_F[s,T];   pJM[s,1]=N_juv_M[s,T];
      pP[s,1] =N_pup[s,T];

      N_total_proj[scen][s,1] = pP[s,1]+pJF[s,1]+pJM[s,1]+pAF[s,1]+pAM[s,1];
      N_pup_proj[scen][s,1]   = pP[s,1];
      N_adult_proj[scen][s,1] = pAF[s,1]+pAM[s,1];

      for (tp in 2:T_proj) {
        real ce = 0;
        if (coyote_idx[s] > 0)
          ce = beta_coy[coyote_idx[s]] * coyote_proj[scen, tp];

        // Pup and juv: sex-neutral in projections
        real pp = inv_logit(phi_pup_logit + site_effect[s] + ce +
                            beta_moci_ond_pup * moci_proj[scen,tp] +
                            beta_moci_amj_pup * moci_proj[scen,tp]);
        real pj = inv_logit(logit(phi_juv_base) + site_effect[s]*0.5 +
                            beta_moci_jfm_juv * moci_proj[scen,tp]);
        // Adult: sex-specific
        real paF = inv_logit(phi_adult_F_logit + site_effect[s]*0.25 +
                             beta_moci_jfm_adult * moci_proj[scen,tp]);
        real paM = inv_logit(logit(phi_adult_M_base) + site_effect[s]*0.25 +
                             beta_moci_jfm_adult * moci_proj[scen,tp]);

        real np  = pAF[s,tp-1] * avg_fecundity;
        real njF = pP[s,tp-1]  * prop_female       * pp;
        real njM = pP[s,tp-1]  * (1-prop_female)   * pp;
        real jsF = pJF[s,tp-1] * pj * (2.0/3.0);
        real jsM = pJM[s,tp-1] * pj * (2.0/3.0);
        real jaF = pJF[s,tp-1] * pj * (1.0/3.0);
        real jaM = pJM[s,tp-1] * pj * (1.0/3.0);

        pP[s,tp]  = np;
        pJF[s,tp] = njF + jsF;   pJM[s,tp] = njM + jsM;
        pAF[s,tp] = pAF[s,tp-1]*paF + jaF;
        pAM[s,tp] = pAM[s,tp-1]*paM + jaM;

        N_total_proj[scen][s,tp] = pP[s,tp]+pJF[s,tp]+pJM[s,tp]+pAF[s,tp]+pAM[s,tp];
        N_pup_proj[scen][s,tp]   = pP[s,tp];
        N_adult_proj[scen][s,tp] = pAF[s,tp]+pAM[s,tp];
      }
    }
    for (tp in 1:T_proj)
      N_total_all_proj[scen][tp] = sum(col(N_total_proj[scen], tp));
    for (tp in 1:(T_proj-1))
      lambda_proj[scen][tp] = N_total_all_proj[scen][tp+1] /
                               N_total_all_proj[scen][tp];
  }
}
'

write_lines(stan_code_v3.2, "harbor_seal_ipm_v3.2.stan")
cat("Stan model written to harbor_seal_ipm_v3.2.stan\n")


# ============================================================================
# PART 2: SIMULATE DATA
# ============================================================================

simulate_seal_ipm_data_v3.2 <- function(T=29, S=6, T_proj=10, seed=123) {
  
  set.seed(seed)
  site_names <- c("BL","DE","DP","PRH","TB","TP")
  N_coy <- 3
  
  true_params <- list(
    # Pup/juv: sex-neutral; adult: sex-specific
    phi_pup_base      = 0.50,            # centred on field-data expectation
    phi_pup_logit     = qlogis(0.50),    # = 0; prior Normal(0, 0.4)
    phi_juv_base      = 0.85,
    phi_adult_F_logit = qlogis(0.90),    # = 2.197; prior Normal(2.20, 0.18)
    phi_adult_F_base  = 0.90,            # kept for simulation loop
    phi_adult_M_base  = 0.85,            # = F_base - delta_adult
    delta_adult       = 0.05,
    
    # 2-class fecundity (collapses prior 4-class structure)
    fecund_primip = 0.60,
    fecund_mature = 0.85,
    avg_fecundity = 0.20*0.60 + 0.80*0.85,   # = 0.80
    prop_female   = 0.50,
    
    # ~10% of adult males hauled out during spring breeding survey
    p_male_breed = 0.10,
    
    beta_coy        = c(-0.15, -0.40, -0.05),
    beta_dist_surv  = c(-0.20, -0.30, -0.15, -0.10, -0.25, -0.15),
    beta_dist_detect = c(-0.20, -0.20, -0.15, -0.10, -0.25, -0.15),
    
    beta_moci_ond_pup   = -0.25, beta_moci_amj_pup   = -0.15,
    beta_moci_jfm_juv   = -0.10, beta_moci_jfm_adult = -0.08,
    beta_moci_amj_molt  =  0.05, beta_eseal_pup      =  0.10,
    
    sigma_process = 0.06, sigma_obs_adult = 0.12,
    sigma_obs_pup = 0.15, sigma_obs_molt  = 0.12,
    sigma_site    = 0.20       # matches new prior Normal(0.2, 0.1) and hard bound 0.5
  )
  
  coyote_idx <- c(1,2,3,0,0,0)
  has_eseal  <- c(0,1,0,1,0,0)
  
  # ── Trajectory design ────────────────────────────────────────────────────
  # Phase 1 (yrs 1–14, ~1997–2010): population grows
  #   MOCI predominantly cool/negative → higher pup survival
  #   Coyote pressure low at all sites
  # Phase 2 (yrs 15–29, ~2011–2025): population declines
  #   MOCI shifts warm/positive → reduced pup survival
  #   Coyote pressure increases at BL, DE, DP
  #   Disturbance increases at human-exposed sites
  # Driven through covariate effects only — no hardcoded survival shift.
  T_inflect <- 10   # peak year; adjust to move inflection point
  
  # ── MOCI: cool early, warm late, smooth 4-yr transition ──────────────────
  moci_mean <- c(rep(-0.6, T_inflect - 2),
                 seq(-0.6, 0.6, length.out = 4),
                 rep( 0.6, T - T_inflect - 2))
  moci_base <- moci_mean + as.vector(arima.sim(list(ar=0.45), n=T)) * 0.7
  moci_jfm  <- as.vector(scale(moci_base))
  moci_amj  <- as.vector(scale(moci_base * 0.85 + as.vector(arima.sim(list(ar=0.3),n=T))*0.4))
  moci_ond  <- as.vector(scale(moci_base * 0.75 + as.vector(arima.sim(list(ar=0.3),n=T))*0.5))
  
  # ── Coyote: low and flat Phase 1, increasing trend Phase 2 ───────────────
  coyote_trend <- c(rep(0, T_inflect),
                    seq(0, 1.2, length.out = T - T_inflect))
  coyote <- matrix(0, S, T)
  coyote[1,] <- as.vector(scale(coyote_trend       + as.vector(arima.sim(list(ar=0.4),n=T))*0.35))
  coyote[2,] <- as.vector(scale(coyote_trend * 1.3 + as.vector(arima.sim(list(ar=0.4),n=T))*0.35))
  coyote[3,] <- as.vector(scale(coyote_trend * 0.5 + as.vector(arima.sim(list(ar=0.4),n=T))*0.35))
  
  # ── Disturbance: moderate random variation with slight upward trend ───────
  dist_trend <- seq(0, 0.4, length.out = T)
  disturbance <- matrix(0, S, T)
  for (s in 1:S)
    disturbance[s,] <- as.vector(scale(dist_trend + as.vector(arima.sim(list(ar=0.3),n=T))*0.6))
  
  # ── Elephant seals: gradual increase at DE and PRH ────────────────────────
  elephant_seal <- matrix(0,S,T)
  elephant_seal[2,] <- as.vector(scale(seq(0,3,length.out=T)+rnorm(T,0,0.5)))
  elephant_seal[4,] <- as.vector(scale(seq(0,4,length.out=T)+rnorm(T,0,0.5)))
  
  site_effect <- rnorm(S, 0, true_params$sigma_site)
  
  # State arrays
  N_adult_F <- N_adult_M <- matrix(NA,S,T)
  N_juv_F   <- N_juv_M   <- matrix(NA,S,T)
  N_pup     <-              matrix(NA,S,T)
  phi_pup   <- phi_juv   <- matrix(NA,S,T)
  phi_adult_F <- phi_adult_M <- matrix(NA,S,T)
  detect_breed <- detect_molt <- matrix(NA,S,T)
  
  # ── Initial populations: moderate 1997 values, room to grow in Phase 1 ───
  N_adult_F[,1] <- c(120, 90, 45, 65, 75, 28)
  N_adult_M[,1] <- c(105, 80, 40, 58, 67, 25)
  N_juv_F[,1]   <- c( 38, 28, 14, 19, 24,  9)
  N_juv_M[,1]   <- c( 34, 25, 12, 17, 21,  8)
  N_pup[,1]     <- c( 90, 68, 34, 46, 56, 22)
  
  for (s in 1:S) {
    for (t in 1:T) {
      ce  <- if (coyote_idx[s]>0) true_params$beta_coy[coyote_idx[s]]*coyote[s,t] else 0
      dse <- true_params$beta_dist_surv[s]   * disturbance[s,t]
      dde <- true_params$beta_dist_detect[s] * disturbance[s,t]
      
      # Pup survival: sex-neutral
      phi_pup[s,t] <- plogis(true_params$phi_pup_logit + site_effect[s] + ce +
                               true_params$beta_moci_ond_pup*moci_ond[t] +
                               true_params$beta_moci_amj_pup*moci_amj[t] +
                               has_eseal[s]*true_params$beta_eseal_pup*elephant_seal[s,t] + dse)
      
      # Juvenile survival: sex-neutral
      phi_juv[s,t] <- plogis(qlogis(true_params$phi_juv_base) + site_effect[s]*0.5 +
                               true_params$beta_moci_jfm_juv*moci_jfm[t])
      
      # Adult survival: sex-specific
      phi_adult_F[s,t] <- plogis(qlogis(true_params$phi_adult_F_base) + site_effect[s]*0.25 +
                                   true_params$beta_moci_jfm_adult*moci_jfm[t])
      phi_adult_M[s,t] <- plogis(qlogis(true_params$phi_adult_M_base) + site_effect[s]*0.25 +
                                   true_params$beta_moci_jfm_adult*moci_jfm[t])
      
      detect_breed[s,t] <- plogis(dde)
      detect_molt[s,t]  <- plogis(dde + true_params$beta_moci_amj_molt*moci_amj[t])
      
      if (t > 1) {
        ep  <- N_adult_F[s,t-1] * true_params$avg_fecundity
        njF <- N_pup[s,t-1] * true_params$prop_female       * phi_pup[s,t]
        njM <- N_pup[s,t-1] * (1-true_params$prop_female)   * phi_pup[s,t]
        jsF <- N_juv_F[s,t-1] * phi_juv[s,t] * (2/3)
        jsM <- N_juv_M[s,t-1] * phi_juv[s,t] * (2/3)
        jaF <- N_juv_F[s,t-1] * phi_juv[s,t] * (1/3)
        jaM <- N_juv_M[s,t-1] * phi_juv[s,t] * (1/3)
        
        N_pup[s,t]     <- exp(rnorm(1,log(max(ep,1)),          true_params$sigma_process))
        N_juv_F[s,t]   <- exp(rnorm(1,log(max(njF+jsF,0.1)),   true_params$sigma_process*0.5))
        N_juv_M[s,t]   <- exp(rnorm(1,log(max(njM+jsM,0.1)),   true_params$sigma_process*0.5))
        N_adult_F[s,t] <- exp(rnorm(1,log(max(N_adult_F[s,t-1]*phi_adult_F[s,t]+jaF,1)),
                                    true_params$sigma_process*0.5))
        N_adult_M[s,t] <- exp(rnorm(1,log(max(N_adult_M[s,t-1]*phi_adult_M[s,t]+jaM,1)),
                                    true_params$sigma_process*0.5))
      }
    }
  }
  
  N_adult_total <- N_adult_F + N_adult_M
  N_juv_total   <- N_juv_F   + N_juv_M
  N_molt_true   <- N_juv_total + N_adult_total
  
  y_adult <- y_pup <- y_molt <- matrix(NA,S,T)
  for (s in 1:S) for (t in 1:T) {
    N_adult_obs  <- N_adult_F[s,t] + N_adult_M[s,t] * true_params$p_male_breed
    y_adult[s,t] <- log(N_adult_obs  * detect_breed[s,t]) + rnorm(1,0,true_params$sigma_obs_adult)
    y_pup[s,t]   <- log(N_pup[s,t]  * detect_breed[s,t]) + rnorm(1,0,true_params$sigma_obs_pup)
    y_molt[s,t]  <- log(N_molt_true[s,t]*detect_molt[s,t])+ rnorm(1,0,true_params$sigma_obs_molt)
  }
  
  n_obs <- S*T
  for (mat in list(y_adult, y_pup, y_molt)) {
    mat[sample(1:n_obs, round(0.05*n_obs))] <- NA
  }
  # (Re-apply: R passes by value)
  y_adult[sample(1:n_obs, round(0.05*n_obs))] <- NA
  y_pup[sample(1:n_obs,   round(0.05*n_obs))] <- NA
  y_molt[sample(1:n_obs,  round(0.05*n_obs))] <- NA
  
  y_adult_obs <- ifelse(is.na(y_adult),0L,1L)
  y_pup_obs   <- ifelse(is.na(y_pup),  0L,1L)
  y_molt_obs  <- ifelse(is.na(y_molt), 0L,1L)
  y_adult[is.na(y_adult)] <- 0
  y_pup[is.na(y_pup)]     <- 0
  y_molt[is.na(y_molt)]   <- 0
  
  N_scenarios    <- 4
  moci_proj      <- matrix(c(0,1,-1,1), N_scenarios, T_proj)
  coyote_proj    <- matrix(c(0,0, 0,1), N_scenarios, T_proj)
  scenario_names <- c("Status Quo","Warm (MOCI +1)","Cool (MOCI -1)","Warm + High Coyote")
  
  stan_data <- list(
    T=T, S=S, N_coy=N_coy,
    y_adult=y_adult, y_pup=y_pup, y_molt=y_molt,
    y_adult_obs=y_adult_obs, y_pup_obs=y_pup_obs, y_molt_obs=y_molt_obs,
    coyote=coyote, disturbance=disturbance, elephant_seal=elephant_seal,
    moci_jfm=as.vector(moci_jfm), moci_amj=as.vector(moci_amj),
    moci_ond=as.vector(moci_ond),
    coyote_idx=coyote_idx, has_eseal=has_eseal,
    T_proj=T_proj, N_scenarios=N_scenarios,
    moci_proj=moci_proj, coyote_proj=coyote_proj
  )
  
  list(stan_data=stan_data, true_params=true_params,
       true_states=list(N_adult_F=N_adult_F, N_adult_M=N_adult_M,
                        N_juv_F=N_juv_F, N_juv_M=N_juv_M, N_pup=N_pup,
                        N_adult_total=N_adult_total, N_molt_true=N_molt_true,
                        phi_pup=phi_pup, phi_juv=phi_juv,
                        phi_adult_F=phi_adult_F, phi_adult_M=phi_adult_M,
                        detect_breed=detect_breed, detect_molt=detect_molt),
       site_names=site_names, years=1997:(1997+T-1), scenario_names=scenario_names)
}


# ============================================================================
# PART 3: PREPARE REAL DATA
# ============================================================================

prepare_real_data_for_ipm_v3.2 <- function(dat, cov_t_scaled, years, T_proj=10) {
  
  site_names <- c("BL","DE","DP","PRH","TB","TP")
  S <- 6; T <- length(years); N_coy <- 3
  
  adult_rows <- seq(1,18,by=3); molt_rows <- seq(2,18,by=3); pup_rows <- seq(3,18,by=3)
  y_adult <- as.matrix(dat[adult_rows,]); rownames(y_adult) <- site_names
  y_molt  <- as.matrix(dat[molt_rows, ]); rownames(y_molt)  <- site_names
  y_pup   <- as.matrix(dat[pup_rows,  ]); rownames(y_pup)   <- site_names
  
  y_adult_obs <- ifelse(is.na(y_adult),0L,1L)
  y_molt_obs  <- ifelse(is.na(y_molt), 0L,1L)
  y_pup_obs   <- ifelse(is.na(y_pup),  0L,1L)
  y_adult[is.na(y_adult)] <- 0; y_molt[is.na(y_molt)] <- 0; y_pup[is.na(y_pup)] <- 0
  
  moci_jfm <- as.vector(cov_t_scaled[1,])
  moci_amj <- as.vector(cov_t_scaled[2,])
  moci_ond <- as.vector(cov_t_scaled[3,])
  
  disturbance <- matrix(0,S,T)
  for (s in 1:S) disturbance[s,] <- as.vector(cov_t_scaled[3+s,])
  
  coyote <- matrix(0,S,T)
  for (k in 1:3) coyote[k,] <- as.vector(cov_t_scaled[9+k,])
  
  elephant_seal <- matrix(0,S,T)
  elephant_seal[2,] <- as.vector(cov_t_scaled[16,])
  elephant_seal[4,] <- as.vector(cov_t_scaled[16,])
  
  coyote_idx <- c(1,2,3,0,0,0)
  has_eseal  <- c(0,1,0,1,0,0)
  
  N_scenarios    <- 4
  moci_proj      <- matrix(c(0,1,-1,1), N_scenarios, T_proj)
  recent_coyote  <- mean(c(mean(coyote[1,(T-4):T]),
                           mean(coyote[2,(T-4):T]),
                           mean(coyote[3,(T-4):T])))
  coyote_proj    <- matrix(recent_coyote, N_scenarios, T_proj)
  coyote_proj[4,] <- recent_coyote + 1
  scenario_names <- c("Status Quo","Warm (MOCI +1)","Cool (MOCI -1)","Warm + High Coyote")
  
  stan_data <- list(
    T=T, S=S, N_coy=N_coy,
    y_adult=y_adult, y_pup=y_pup, y_molt=y_molt,
    y_adult_obs=y_adult_obs, y_pup_obs=y_pup_obs, y_molt_obs=y_molt_obs,
    coyote=coyote, disturbance=disturbance, elephant_seal=elephant_seal,
    moci_jfm=moci_jfm, moci_amj=moci_amj, moci_ond=moci_ond,
    coyote_idx=coyote_idx, has_eseal=has_eseal,
    T_proj=T_proj, N_scenarios=N_scenarios,
    moci_proj=moci_proj, coyote_proj=coyote_proj
  )
  
  list(stan_data=stan_data, site_names=site_names, years=years,
       scenario_names=scenario_names,
       raw_counts=list(adult=y_adult, molt=y_molt, pup=y_pup))
}


# ============================================================================
# PART 4: DIAGNOSTICS
# ============================================================================

check_diagnostics_v3.2 <- function(fit) {
  
  cat("\n============================================\n")
  cat("MODEL DIAGNOSTICS — IPM v3.2\n")
  cat("============================================\n\n")
  fit$cmdstan_diagnose()
  
  params <- c(
    "phi_pup_logit","phi_juv_base","phi_adult_F_logit","phi_adult_F_base","delta_adult",
    "fecund_primip","fecund_mature",
    "prop_female","p_male_breed",
    "beta_coy[1]","beta_coy[2]","beta_coy[3]",
    "beta_dist_surv[1]","beta_dist_surv[2]","beta_dist_surv[3]",
    "beta_dist_surv[4]","beta_dist_surv[5]","beta_dist_surv[6]",
    "beta_moci_ond_pup","beta_moci_amj_pup","beta_moci_jfm_juv",
    "beta_moci_jfm_adult","beta_eseal_pup",
    "sigma_process","sigma_obs_adult","sigma_obs_pup","sigma_obs_molt","sigma_site"
  )
  
  s <- fit$summary(variables=params)
  cat("\nParameter Summary:\n")
  print(s |> select(variable,mean,sd,q5,q95,rhat,ess_bulk))
  
  # Convenience: report phi_pup on probability scale
  pup_logit_draws <- fit$draws(variables="phi_pup_logit", format="df")$phi_pup_logit
  cat(sprintf("\nphi_pup_base (prob scale): median=%.3f, 90%% CrI=[%.3f, %.3f]\n",
              median(plogis(pup_logit_draws)),
              quantile(plogis(pup_logit_draws), 0.05),
              quantile(plogis(pup_logit_draws), 0.95)))
  
  bad   <- s |> filter(rhat>1.05)
  low   <- s |> filter(ess_bulk<400)
  if (nrow(bad)>0) { cat("\nWARNING: Rhat > 1.05:\n"); print(bad |> select(variable,rhat)) }
  if (nrow(low)>0) { cat("\nWARNING: low ESS:\n");     print(low |> select(variable,ess_bulk)) }
  
  list(params=params, summary=s)
}


# ============================================================================
# PART 5: TRACE PLOTS
# ============================================================================

create_trace_plots_v3.2 <- function(fit, params, save=TRUE, prefix="IPM_v3.2") {
  
  draws <- fit$draws(format="df")
  
  p1 <- mcmc_trace(draws,
                   pars=c("phi_pup_logit","phi_juv_base","phi_adult_F_logit",
                          "delta_adult","p_male_breed")) +
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
  
  p4 <- mcmc_trace(draws, pars=paste0("beta_dist_surv[",1:6,"]")) +
    labs(title="Trace: Site-Specific Disturbance Effects")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_trace_disturbance.jpeg"),
                   p4, width=30, height=18, units="cm")
  
  p5 <- mcmc_trace(draws,
                   pars=c("beta_moci_ond_pup","beta_moci_amj_pup","beta_moci_jfm_juv",
                          "beta_moci_jfm_adult","beta_moci_amj_molt")) +
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
      "beta_moci_ond_pup","beta_moci_amj_pup","beta_moci_jfm_juv",
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
      tp$beta_moci_ond_pup, tp$beta_moci_amj_pup, tp$beta_moci_jfm_juv,
      tp$beta_moci_jfm_adult, tp$beta_eseal_pup,
      tp$sigma_process, tp$sigma_obs_adult, tp$sigma_obs_pup, tp$sigma_obs_molt
    )
  )
  
  rec <- fit$summary(variables=true_vals$parameter) |>
    left_join(true_vals, by=c("variable"="parameter")) |>
    mutate(recovered    = true_value>=q5 & true_value<=q95,
           rel_bias_pct = (mean-true_value)/abs(true_value)*100)
  
  cat("Recovery rate:", sum(rec$recovered),"/",nrow(rec),
      sprintf("(%.1f%%)\n", 100*mean(rec$recovered)))
  print(rec |> select(variable,true_value,mean,q5,q95,recovered,rel_bias_pct) |>
          mutate(across(where(is.numeric),~round(.x,3))))
  
  p <- ggplot(rec, aes(x=true_value,y=mean)) +
    geom_abline(slope=1,intercept=0,linetype=2,color="gray50") +
    geom_pointrange(aes(ymin=q5,ymax=q95,color=recovered),size=0.8) +
    geom_text(aes(label=variable),hjust=-0.1,vjust=-0.3,size=2.5,check_overlap=TRUE) +
    scale_color_manual(values=c("TRUE"="blue","FALSE"="red")) +
    labs(x="True Value",y="Estimated (90% CI)",title="Parameter Recovery: IPM v3.2") +
    theme_minimal(base_size=14) + coord_equal()
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
                           lo=apply(Ntot,2,quantile,0.05),
                           hi=apply(Ntot,2,quantile,0.95)),
                    aes(x=Year)) +
    geom_ribbon(aes(ymin=lo,ymax=hi),alpha=0.3,fill="blue") +
    geom_line(aes(y=mean),linewidth=1.2,color="blue") +
    labs(x="Year",y="Total Population",
         title="Estimated Total Harbor Seal Population") +
    theme_minimal(base_size=14)
  if (save) ggsave(paste0("Output/Plots/",prefix,"_total_population.jpeg"),
                   p_total, width=25, height=15, units="cm")
  
  # Pup survival on probability scale
  pup_logit <- fit$draws(variables="phi_pup_logit",format="df")$phi_pup_logit
  pup_prob  <- plogis(pup_logit)
  p_phipup  <- ggplot(tibble(x=pup_prob), aes(x=x)) +
    geom_density(fill="#2ca25f",alpha=0.6) +
    geom_vline(xintercept=median(pup_prob),linetype="dashed",color="darkgreen") +
    geom_vline(xintercept=0.5, linetype="dotted", color="gray50") +
    annotate("text",x=0.5,y=Inf,label=" Field data ~0.50",
             hjust=0,vjust=1.5,size=3.5,color="gray40") +
    labs(x="Pup survival probability", y="Density",
         title="Posterior: Pup Survival (sex-neutral)",
         subtitle=sprintf("Median = %.3f  (95%% CrI: %.3f–%.3f)",
                          median(pup_prob),
                          quantile(pup_prob,0.025),
                          quantile(pup_prob,0.975))) +
    theme_minimal(base_size=14)
  if (save) ggsave(paste0("Output/Plots/",prefix,"_phi_pup_posterior.jpeg"),
                   p_phipup, width=20, height=12, units="cm")
  
  # Male haul-out fraction
  pmb <- fit$draws(variables="p_male_breed",format="df")$p_male_breed
  p_pmb <- ggplot(tibble(x=pmb),aes(x=x)) +
    geom_density(fill="#4393c3",alpha=0.6) +
    geom_vline(xintercept=median(pmb),linetype="dashed",color="#08306b") +
    labs(x="p_male_breed",y="Density",
         title="Posterior: Male Haul-out Fraction (Breeding Season)",
         subtitle=sprintf("Median=%.3f  (95%% CrI: %.3f–%.3f)",
                          median(pmb),quantile(pmb,0.025),quantile(pmb,0.975))) +
    theme_minimal(base_size=14)
  if (save) ggsave(paste0("Output/Plots/",prefix,"_p_male_breed.jpeg"),
                   p_pmb, width=20, height=12, units="cm")
  
  # Sex ratio: true vs observed
  sr_t <- fit$draws(variables="sex_ratio_adult",   format="matrix")
  sr_o <- fit$draws(variables="sex_ratio_observed",format="matrix")
  
  make_sr <- function(mat) {
    m <- sapply(1:T, function(t) rowMeans(mat[,paste0("sex_ratio_adult[",1:S,",",t,"]"),drop=FALSE]))
    tibble(Year=years, mean=colMeans(m),
           lo=apply(m,2,quantile,0.05), hi=apply(m,2,quantile,0.95))
  }
  make_sr_obs <- function(mat) {
    m <- sapply(1:T, function(t) rowMeans(mat[,paste0("sex_ratio_observed[",1:S,",",t,"]"),drop=FALSE]))
    tibble(Year=years, mean=colMeans(m),
           lo=apply(m,2,quantile,0.05), hi=apply(m,2,quantile,0.95))
  }
  
  sr_true_df <- tryCatch(make_sr(sr_t),    error=function(e) NULL)
  sr_obs_df  <- tryCatch(make_sr_obs(sr_o),error=function(e) NULL)
  
  p_sr <- ggplot() +
    geom_hline(yintercept=0.5,linetype="dotted",color="gray50") +
    { if (!is.null(sr_obs_df))  geom_ribbon(data=sr_obs_df,  aes(x=Year,ymin=lo,ymax=hi),alpha=0.15,fill="orange") } +
    { if (!is.null(sr_obs_df))  geom_line(data=sr_obs_df,    aes(x=Year,y=mean,color="Observed (spring)"),linewidth=1.1,linetype="dashed") } +
    { if (!is.null(sr_true_df)) geom_ribbon(data=sr_true_df, aes(x=Year,ymin=lo,ymax=hi),alpha=0.25,fill="purple") } +
    { if (!is.null(sr_true_df)) geom_line(data=sr_true_df,   aes(x=Year,y=mean,color="True (population)"),linewidth=1.2) } +
    scale_color_manual(values=c("True (population)"="purple","Observed (spring)"="darkorange")) +
    labs(x="Year",y="Proportion female",color=NULL,
         title="Adult Sex Ratio: True vs Spring Survey Observation",
         subtitle="Observed > true because most males remain in water during breeding") +
    ylim(0.45,0.80) + theme_minimal(base_size=14) + theme(legend.position="bottom")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_sex_ratio.jpeg"),
                   p_sr, width=25, height=15, units="cm")
  
  # Mean pup survival over time (sex-neutral)
  mp <- fit$draws(variables="mean_phi_pup", format="matrix")
  mp_df <- tibble(Year=years, mean=colMeans(mp),
                  lo=apply(mp,2,quantile,0.05), hi=apply(mp,2,quantile,0.95))
  p_phipup_t <- ggplot(mp_df,aes(x=Year)) +
    geom_ribbon(aes(ymin=lo,ymax=hi),alpha=0.25,fill="#2ca25f") +
    geom_line(aes(y=mean),linewidth=1.0,color="darkgreen") +
    geom_hline(yintercept=0.5,linetype="dotted",color="gray50") +
    labs(x="Year",y="Pup survival (sex-neutral)",
         title="Mean Pup Survival Over Time (across all sites)") +
    theme_minimal(base_size=13)
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
             mean=mean(d[,cn]), lo=quantile(d[,cn],.05), hi=quantile(d[,cn],.95))
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
    theme_minimal(base_size=14) + theme(legend.position="bottom")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_age_class_timeseries.jpeg"),
                   p, width=30, height=35, units="cm")
  
  list(by_age=p, data=all_sum)
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
           mean=colMeans(pm), lo=apply(pm,2,quantile,.10), hi=apply(pm,2,quantile,.90))
  })
  
  Ntot <- fit$draws(variables="N_total_all",format="df") |> select(starts_with("N_total_all"))
  hist <- tibble(Scenario="Historical",Year=years,mean=colMeans(Ntot),
                 lo=apply(Ntot,2,quantile,.10),hi=apply(Ntot,2,quantile,.90))
  
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
    scale_color_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") +
    labs(x="Year",y="Total Population",
         title="10-Year Projections",
         subtitle="Bands = 80% CI; dashed line = projection start") +
    theme_minimal(base_size=14) + theme(legend.position="bottom")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_projections.jpeg"),
                   p, width=30, height=20, units="cm")
  
  list(projection=p, data=proj_df)
}


# ============================================================================
# PART 11: COVARIATE EFFECT PLOTS
# ============================================================================

create_effect_plots_v3.2 <- function(fit, save=TRUE, prefix="IPM_v3.2") {
  
  site_names <- c("BL","DE","DP","PRH","TB","TP")
  
  # For pup survival, baseline is on logit scale
  make_plot <- function(param, xlab, title, xr=seq(-2,2,length.out=100), ylims=NULL) {
    beta_v  <- fit$draws(variables=param,       format="df")[[param]]
    logit_v <- fit$draws(variables="phi_pup_logit", format="df")$phi_pup_logit
    idx     <- sample(seq_along(beta_v), min(500,length(beta_v)))
    df      <- map_dfr(idx, ~tibble(x=xr,
                                    survival=plogis(logit_v[.x]+beta_v[.x]*xr),draw=.x))
    sm <- df |> group_by(x) |>
      summarise(mean=mean(survival),lo=quantile(survival,.05),
                hi=quantile(survival,.95),.groups="drop")
    clr <- ifelse(mean(beta_v)<0,"red3","blue3")
    p <- ggplot(sm,aes(x=x)) +
      geom_ribbon(aes(ymin=lo,ymax=hi),alpha=.2,fill=clr) +
      geom_line(aes(y=mean),linewidth=1.2,color=clr) +
      geom_hline(yintercept=mean(plogis(logit_v)),linetype=2,color="gray50") +
      geom_vline(xintercept=0,linetype=2,color="gray50") +
      labs(x=xlab,y="Pup Survival",title=title) + theme_minimal(base_size=11)
    if (!is.null(ylims)) p <- p + coord_cartesian(ylim=ylims)
    list(plot=p, yr=c(min(sm$lo),max(sm$hi)))
  }
  
  coy <- lapply(1:3, function(i) make_plot(paste0("beta_coy[",i,"]"),
                                           "Coyote (SD)", paste0("Coyote → Pup Survival (",c("BL","DE","DP")[i],")")))
  ylc <- c(max(0,min(sapply(coy,`[[`,"yr"))-0.03), min(1,max(sapply(coy,`[[`,"yr"))+0.03))
  p_coy <- wrap_plots(lapply(1:3, function(i)
    make_plot(paste0("beta_coy[",i,"]"),"Coyote (SD)",
              paste0("(",c("BL","DE","DP")[i],")"), ylims=ylc)$plot), ncol=3) +
    plot_annotation(title="Site-Specific Coyote Effects on Pup Survival")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_effects_coyote.jpeg"),
                   p_coy, width=36, height=12, units="cm")
  
  dst <- lapply(1:6, function(s) make_plot(paste0("beta_dist_surv[",s,"]"),
                                           "Disturbance (SD)", paste0("(",site_names[s],")")))
  yld <- c(max(0,min(sapply(dst,`[[`,"yr"))-0.03), min(1,max(sapply(dst,`[[`,"yr"))+0.03))
  p_dst <- wrap_plots(lapply(1:6, function(s)
    make_plot(paste0("beta_dist_surv[",s,"]"),"Disturbance (SD)",
              paste0("(",site_names[s],")"), ylims=yld)$plot), ncol=3) +
    plot_annotation(title="Site-Specific Disturbance Effects on Pup Survival")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_effects_disturbance.jpeg"),
                   p_dst, width=36, height=24, units="cm")
  
  mci <- lapply(list(c("beta_moci_ond_pup","MOCI Fall (SD)","MOCI Fall → Pup Survival"),
                     c("beta_moci_amj_pup","MOCI Spring (SD)","MOCI Spring → Pup Survival"),
                     c("beta_eseal_pup","Elephant Seal (SD)","Elephant Seal → Pup Survival")),
                function(x) make_plot(x[1],x[2],x[3]))
  ylm <- c(max(0,min(sapply(mci,`[[`,"yr"))-0.03), min(1,max(sapply(mci,`[[`,"yr"))+0.03))
  p_mci <- wrap_plots(lapply(list(c("beta_moci_ond_pup","MOCI Fall (SD)","MOCI Fall"),
                                  c("beta_moci_amj_pup","MOCI Spring (SD)","MOCI Spring"),
                                  c("beta_eseal_pup","Elephant Seal (SD)","Eseal")),
                             function(x) make_plot(x[1],x[2],x[3],ylims=ylm)$plot), ncol=3) +
    plot_annotation(title="Shared Covariate Effects on Pup Survival")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_effects_moci.jpeg"),
                   p_mci, width=36, height=12, units="cm")
  
  list(coyote=p_coy, disturbance=p_dst, moci=p_mci)
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
    "beta_moci_ond_pup","beta_moci_amj_pup","beta_moci_jfm_juv",
    "beta_moci_jfm_adult","beta_eseal_pup","beta_moci_amj_molt",
    "sigma_process","sigma_obs_adult","sigma_obs_pup","sigma_obs_molt","sigma_site"
  )
  
  # Derive phi_pup_base on probability scale for reporting
  pup_draws <- fit$draws(variables="phi_pup_logit",format="df")$phi_pup_logit
  pup_base_summary <- tibble(
    variable="phi_pup_base (prob)",
    mean=mean(plogis(pup_draws)), sd=sd(plogis(pup_draws)),
    q5=quantile(plogis(pup_draws),.05), q95=quantile(plogis(pup_draws),.95),
    rhat=NA_real_, ess_bulk=NA_real_
  )
  
  tbl <- fit$summary(variables=params) |>
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
        variable=="beta_moci_ond_pup"    ~ "MOCI fall → pup survival",
        variable=="beta_moci_amj_pup"    ~ "MOCI spring → pup survival",
        variable=="beta_moci_jfm_juv"    ~ "MOCI winter → juvenile survival",
        variable=="beta_moci_jfm_adult"  ~ "MOCI winter → adult survival",
        variable=="beta_eseal_pup"       ~ "Elephant seal → pup survival",
        variable=="beta_moci_amj_molt"   ~ "MOCI spring → molt detection",
        variable=="sigma_process"        ~ "Process error (σ)",
        variable=="sigma_obs_adult"      ~ "Obs error adult (σ)",
        variable=="sigma_obs_pup"        ~ "Obs error pup (σ)",
        variable=="sigma_obs_molt"       ~ "Obs error molt (σ)",
        variable=="sigma_site"           ~ "Site random effect (σ)",
        TRUE ~ variable
      ),
      Estimate = sprintf("%.3f (%.3f, %.3f)", mean, q5, q95),
      Category = case_when(
        str_detect(variable,"phi|delta")              ~ "Survival",
        str_detect(variable,"fecund|prop|avg")        ~ "Reproduction",
        variable %in% c("p_male_breed","phi_pup_base (prob)") ~ "Observation / Derived",
        str_detect(variable,"beta_coy")               ~ "Coyote (site-specific)",
        str_detect(variable,"beta_dist")              ~ "Disturbance (site-specific)",
        str_detect(variable,"beta_moci|beta_eseal")   ~ "Shared covariates",
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
  
  surv  <- fit$summary(variables=c("phi_juv_base","phi_adult_F_base","delta_adult"))
  logit_aF <- fit$draws(variables="phi_adult_F_logit",format="df")$phi_adult_F_logit
  pup_l <- fit$draws(variables="phi_pup_logit",format="df")$phi_pup_logit
  pmb   <- fit$summary(variables="p_male_breed")
  coy   <- fit$summary(variables=paste0("beta_coy[",1:3,"]"))
  dst   <- fit$summary(variables=paste0("beta_dist_surv[",1:6,"]"))
  sns   <- c("BL","DE","DP","PRH","TB","TP")
  
  cat("\n============================================\n")
  cat("HARBOR SEAL IPM v3.2 — KEY RESULTS\n")
  cat("============================================\n")
  cat(sprintf("Pup survival (sex-neutral, prob):  %.3f (%.3f–%.3f)\n",
              median(plogis(pup_l)), quantile(plogis(pup_l),.025), quantile(plogis(pup_l),.975)))
  cat(sprintf("Juvenile survival (sex-neutral):   %.3f (%.3f–%.3f)\n",
              surv$mean[1],surv$q5[1],surv$q95[1]))
  cat(sprintf("Adult female survival (baseline):  %.3f (%.3f–%.3f)\n",
              surv$mean[2],surv$q5[2],surv$q95[2]))
  cat(sprintf("Adult male survival (F - delta):   %.3f–%.3f\n",
              surv$mean[2]-surv$mean[3]-2*surv$sd[3], surv$mean[2]-surv$mean[3]+2*surv$sd[3]))
  cat(sprintf("Male breeding haul-out (p_male_breed): %.3f (%.3f–%.3f)\n",
              pmb$mean,pmb$q5,pmb$q95))
  cat("\nCoyote effects on pup survival:\n")
  for (i in 1:3) cat(sprintf("  %s: %.3f (%.3f–%.3f)\n",c("BL","DE","DP")[i],
                             coy$mean[i],coy$q5[i],coy$q95[i]))
  cat("\nDisturbance effects on pup survival:\n")
  for (i in 1:6) cat(sprintf("  %s: %.3f (%.3f–%.3f)\n",sns[i],
                             dst$mean[i],dst$q5[i],dst$q95[i]))
  cat("============================================\n")
}


# ============================================================================
# PART 14: MAIN EXECUTION
# ============================================================================

run_full_analysis_v3.2 <- function(use_real_data  = FALSE,
                                   dat            = NULL,
                                   cov_t_scaled   = NULL,
                                   years          = NULL,
                                   T_proj         = 10,
                                   seed           = 42,
                                   iter_warmup    = 3000,
                                   iter_sampling  = 1000,
                                   adapt_delta    = 0.995,
                                   max_treedepth  = 15) {
  
  cat("\n================================================================\n")
  cat("   HARBOR SEAL IPM v3.2\n")
  cat("   Sex-neutral pup/juv survival | Corrected pup prior\n")
  cat("   Explicit observation sex structure\n")
  cat("================================================================\n\n")
  
  prefix <- ifelse(use_real_data,"IPM_v3.2_real","IPM_v3.2_sim")
  
  if (use_real_data) {
    dl       <- prepare_real_data_for_ipm_v3.2(dat,cov_t_scaled,years,T_proj)
    sim_data <- list(stan_data=dl$stan_data, site_names=dl$site_names,
                     years=dl$years, scenario_names=dl$scenario_names, true_params=NULL)
  } else {
    sim_data <- simulate_seal_ipm_data_v3.2(T=29,S=6,T_proj=T_proj,seed=seed)
  }
  
  cat("Compiling Stan model...\n")
  model <- cmdstan_model("harbor_seal_ipm_v3.2.stan")
  
  cat(sprintf("Running MCMC (warmup=%d, sampling=%d, adapt_delta=%.3f)...\n",
              iter_warmup,iter_sampling,adapt_delta))
  fit <- model$sample(
    data=sim_data$stan_data, seed=123, chains=4, parallel_chains=4,
    iter_warmup=iter_warmup, iter_sampling=iter_sampling,
    refresh=200, adapt_delta=adapt_delta, max_treedepth=max_treedepth
  )
  fit$save_object(paste0("Output/harbor_seal_",prefix,"_fit.rds"))
  
  diag   <- check_diagnostics_v3.2(fit)
  traces <- create_trace_plots_v3.2(fit, diag$params, prefix=prefix)
  rec    <- if (!use_real_data && !is.null(sim_data$true_params))
    check_parameter_recovery_v3.2(fit,sim_data,prefix=prefix) else NULL
  ppc    <- create_ppc_plots_v3.2(fit,sim_data,prefix=prefix)
  ts     <- create_timeseries_plots_v3.2(fit,sim_data,prefix=prefix)
  sa     <- create_site_age_timeseries_v3.2(fit,sim_data,prefix=prefix)
  proj   <- create_projection_plots_v3.2(fit,sim_data,prefix=prefix)
  eff    <- create_effect_plots_v3.2(fit,prefix=prefix)
  tbl    <- create_summary_table_v3.2(fit,prefix=prefix)
  save_model_output_v3.2(fit,prefix=prefix)
  
  cat("\n================================================================\n")
  cat("   COMPLETE — IPM v3.2\n")
  cat(sprintf("   Plots  → Output/Plots/%s_*.jpeg\n",prefix))
  cat(sprintf("   Fit    → Output/harbor_seal_%s_fit.rds\n",prefix))
  cat(sprintf("   Table  → Output/%s_parameter_summary.csv\n",prefix))
  cat("================================================================\n\n")
  
  list(fit=fit,model=model,data=sim_data,diagnostics=diag,
       summary=tbl,recovery=rec,projections=proj,prefix=prefix)
}


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
    theme_minimal(base_size=12) +
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
    theme_minimal(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=8),panel.grid=element_blank())
  if (save) ggsave(paste0("Output/Plots/",prefix,"_portfolio_pup_survival_heatmap.jpeg"),
                   p_phi, width=35, height=15, units="cm")
  
  p_async <- ggplot(ldf,aes(x=Year,y=lambda,color=Site,group=Site)) +
    geom_hline(yintercept=1,linetype=2,color="gray50") +
    geom_line(linewidth=1,alpha=0.7) + geom_point(size=2) +
    scale_color_brewer(palette="Set1") +
    labs(x="Year",y=expression(lambda),title="Site-Level Population Growth — Asynchrony") +
    theme_minimal(base_size=14) + theme(legend.position="bottom")
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
    theme_minimal(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=8),panel.grid=element_blank())
  if (save) ggsave(paste0("Output/Plots/",prefix,"_portfolio_best_worst.jpeg"),
                   p_bw, width=35, height=12, units="cm")
  
  cordf <- expand.grid(Site1=site_names,Site2=site_names) |>
    mutate(r=as.vector(lcor))
  p_cor <- ggplot(cordf,aes(x=Site1,y=Site2,fill=r)) +
    geom_tile(color="white") + geom_text(aes(label=sprintf("%.2f",r)),size=4) +
    scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0,limits=c(-1,1)) +
    labs(x="",y="",title="Between-Site Correlation in λ") +
    theme_minimal(base_size=12) + theme(panel.grid=element_blank()) + coord_fixed()
  if (save) ggsave(paste0("Output/Plots/",prefix,"_portfolio_correlation.jpeg"),
                   p_cor, width=20, height=18, units="cm")
  
  list(lambda_heatmap=p_heat, phi_heatmap=p_phi, asynchrony=p_async,
       best_worst=p_bw, correlation=p_cor,
       lambda_matrix=lmat, phi_matrix=phi_m,
       contributions=NULL,
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
  phi_aF_logit <- draws$phi_adult_F_logit[idx]   # logit scale — use directly in inv_logit()
  delta_a  <- draws$delta_adult[idx]
  pf       <- draws$prop_female[idx]
  avgf     <- draws$avg_fecundity[idx]
  b_ond    <- draws$beta_moci_ond_pup[idx]
  b_amj    <- draws$beta_moci_amj_pup[idx]
  b_jfmJ   <- draws$beta_moci_jfm_juv[idx]
  b_jfmA   <- draws$beta_moci_jfm_adult[idx]
  
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
          pp  <- plogis(pup_l[i]+se[i,s]+ce+b_ond[i]*sc$moci+b_amj[i]*sc$moci)
          pj  <- plogis(qlogis(phi_juv[i])+se[i,s]*0.5+b_jfmJ[i]*sc$moci)
          paF <- plogis(phi_aF_logit[i] + se[i,s]*0.25 + b_jfmA[i]*sc$moci)
          paM <- plogis(qlogis(plogis(phi_aF_logit[i]) - delta_a[i]) + se[i,s]*0.25 + b_jfmA[i]*sc$moci)
          new_p <- naf[s]*avgf[i]*exp(sh[s])
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
             lo=apply(res[[sn]]$async,2,quantile,.10),
             hi=apply(res[[sn]]$async,2,quantile,.90)),
      tibble(Scenario=sn,Year=pyrs,Synchrony="Synchronous (hypothetical)",
             mean=colMeans(res[[sn]]$sync),
             lo=apply(res[[sn]]$sync,2,quantile,.10),
             hi=apply(res[[sn]]$sync,2,quantile,.90))
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
    scale_color_manual(values=c("Asynchronous (current)"="blue3","Synchronous (hypothetical)"="red3")) +
    scale_fill_manual( values=c("Asynchronous (current)"="blue", "Synchronous (hypothetical)"="red")) +
    scale_linetype_manual(values=c("Asynchronous (current)"="solid","Synchronous (hypothetical)"="dashed")) +
    labs(x="Year",y="Total Population",title="Portfolio Buffering: Async vs Sync") +
    theme_minimal(base_size=12) + theme(legend.position="bottom")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_synchrony_comparison.jpeg"),
                   p_comp, width=30, height=25, units="cm")
  
  p_cv <- cv_df |>
    pivot_longer(c(CV_Async,CV_Sync),names_to="Type",values_to="CV") |>
    mutate(Type=recode(Type,CV_Async="Asynchronous",CV_Sync="Synchronous")) |>
    ggplot(aes(x=Scenario,y=CV,fill=Type)) +
    geom_col(position="dodge",alpha=0.8) +
    scale_fill_manual(values=c(Asynchronous="blue3",Synchronous="red3")) +
    labs(x="Scenario",y="CV of aggregate abundance",title="Portfolio Buffering by Scenario") +
    theme_minimal(base_size=12) +
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
cat("============================================\n")
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
