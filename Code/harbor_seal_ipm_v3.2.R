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
  real phi_pup_logit;           // logit-scale; ~ Normal(-1.2, 0.5) → median p ≈ 0.23
                                // Based on lit: Hansen 2013 6-mo=0.39; Oates 6-mo=0.48;
                                //   Lander 2002 CA 5-mo=0.20-0.60; Bigg 1969 1st-yr=35-80%
                                // Annual ≈ (6-mo survival)^2 ≈ 0.15-0.23 (constant hazard)
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
  // Free baseline detection probabilities (logit scale)
  // detect_breed_logit: fraction of breeding adults/pups detectable at spring survey
  // detect_molt_logit:  fraction of juv+adult detectable at molt survey
  // Priors informed by Pacific harbor seal haul-out studies:
  //   breeding: 0.70-0.85 typical → logit(0.77) ≈ 1.2
  //   molt:     0.60-0.75 typical → logit(0.68) ≈ 0.75
  real detect_breed_logit;
  real detect_molt_logit;

  // ── Shared covariate effects ───────────────────────────────────────────────
  // beta_moci_ond_fecund: effect of fall MOCI on fecundity (maternal condition at conception)
  // moci_ond is pre-lagged in data prep: moci_ond[t] = OND of year t-1
  real beta_moci_ond_fecund;   // OND MOCI → fecundity (maternal condition at conception)
  real beta_moci_ond_pup;      // OND MOCI → pup survival (post-weaning fall foraging, year t-1)
  real beta_moci_amj_pup;
  real beta_moci_jfm_pup;      // JFM MOCI → pup survival (first winter, year t)
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
      // t_birth: index for covariates that act during pup birth year (t-1).
      // phi_pup[s,1] is computed but never enters dynamics (loop starts at t=2),
      // so using index 1 at t=1 is harmless.
      int t_birth = (t > 1) ? t - 1 : 1;

      // Coyote and disturbance act during birth year (pup first summer/fall)
      real coyote_effect = 0;
      if (coyote_idx[s] > 0)
        coyote_effect = beta_coy[coyote_idx[s]] * coyote[s, t_birth];

      // dist_surv: birth-year predator/disturbance pressure on pup survival
      // dist_detect: detection at counting time (year t) — stays at t
      real dist_surv_eff   = beta_dist_surv[s]   * disturbance[s, t_birth];
      real dist_detect_eff = beta_dist_detect[s] * disturbance[s, t];

      // Pup survival — sex-neutral
      // All pup survival covariates now use t_birth (year t-1 = pup birth year):
      // coyote, disturbance, eseal, AND MOCI AMJ all act during nursing/post-weaning.
      // Pup survival covariates:
      //   t_birth = t-1: coyote, dist, eseal, MOCI AMJ (birth-year nursing period)
      //   t: MOCI OND = OND of year t-1 (post-weaning fall, via pre-lagged data)
      //      MOCI JFM = JFM of year t (first winter)
      phi_pup[s, t] = inv_logit(
        phi_pup_logit + site_effect[s] + coyote_effect +
        beta_moci_amj_pup * moci_amj[t_birth] +
        beta_moci_ond_pup  * moci_ond[t]       +
        beta_moci_jfm_pup  * moci_jfm[t]       +
        has_eseal[s] * beta_eseal_pup * elephant_seal[s, t_birth] +
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

      // Detection: free logit-scale baseline + covariate modifiers
      // detect_breed: breeding survey (spring) — adults + pups together
      // detect_molt:  molt survey (summer) — all juv + adult stages
      detect_breed[s, t] = inv_logit(detect_breed_logit + dist_detect_eff);
      detect_molt[s, t]  = inv_logit(detect_molt_logit  + dist_detect_eff +
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
      // Pups produced by adult females: fecundity modulated by fall MOCI
      // moci_ond[t] = OND of year t-1 (pre-lagged) — maternal energy at conception
      real fecund_t     = inv_logit(logit(avg_fecundity) +
                                    beta_moci_ond_fecund * moci_ond[t]);
      real expected_pups = N_adult_F[s, t-1] * fecund_t;

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

  // Pup survival prior updated from literature:
  //   Hansen 2013 (Scotland): P(survive 6 mo) = 0.390 (95% CI: 0.297-0.648)
  //   Oates thesis (CA): P(survive 6 mo) = 0.48 (SE = 0.10)
  //   Lander 2002 (central CA): P(survive 5 mo) = 0.20 (1995), 0.60 (1996)
  //   Bigg/Van Bemmel/Reijnders life tables: first-year 35-80%
  //   Annual ≈ 6-mo^2 under constant hazard: 0.15-0.23; life-table range 0.35-0.80
  //   Prior centered at annual ≈ 0.23 (logit = -1.2); SD = 0.5 spans 0.10-0.45
  phi_pup_logit    ~ normal(-1.2, 0.5);

  // Juvenile survival (ages 1-3): beta prior calibrated to conspecific mark-recapture
  //   Hastings et al. 2012 (P. v. richardii, Alaska): ages 1-3 female=0.865, male=0.782
  //   Sex-neutral average ≈ 0.82; Hansen 2013 deterministic model also requires ~0.82
  //   for 11%/yr decline when pup survival and fecundity are in observed range
  //   Prior: beta(16,4) → mean = 0.80, SD ≈ 0.089; spans 0.60–0.94 at ±2SD
  //   (vs. beta(14,6)→0.70 used previously; that was too pessimistic vs literature)
  phi_juv_base     ~ beta(16, 4);    // mean = 0.80, SD ≈ 0.089
  // Adult female survival: prior calibrated to PRNS-specific and conspecific estimates
  //   Manugian 2017 (Tomales Bay, adult females): 0.90 (95% CI 0.18-0.99) — most
  //     geographically relevant; Tomales Bay is one of our PRNS monitoring sites
  //   Härkönen & Heide-Jørgensen 1990 (East Atlantic): annual female = 0.902
  //   Mackey 2008, Cordes & Thompson 2014: 92-98% (other Pacific populations)
  //   Bigg 1969 life tables: 0.85
  //   Consensus: 0.90 best single estimate for PRNS; logit(0.90) = 2.197
  //   Prior: normal(2.20, 0.25) → median 0.900; ±1SD spans (0.865, 0.932)
  //   (vs. old normal(2.44,0.25)→0.920 which was slightly above PRNS evidence)
  phi_adult_F_logit ~ normal(2.20, 0.25);
  // delta_adult: literature is genuinely inconsistent on sex difference direction
  //   Bigg 1969: F=0.85 vs M=0.71 → Δ≈0.14 (large female advantage)
  //   Härkönen 1990: M=0.91 vs F=0.902 → slight male advantage
  //   No direct PRNS-specific male-vs-female comparison available
  //   Mean kept at 0.05 (small female advantage); SD widened to 0.025
  //   to reflect genuine uncertainty while keeping F≥M (lower=0 bound)
  delta_adult       ~ normal(0.05, 0.025);   // F ≈ 5% higher than M; wider SD = more uncertain

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
  // Detection baselines: logit-scale priors based on Pacific harbor seal haul-out rates
  // breeding survey: mean 0.77 (70-85% haul-out during peak breeding)
  // molt survey:     mean 0.68 (60-75% haul-out during peak molt)
  detect_breed_logit ~ normal(1.20, 0.50);  // plogis(1.20) ≈ 0.77
  detect_molt_logit  ~ normal(0.75, 0.50);  // plogis(0.75) ≈ 0.68

  // Shared covariate effects
  beta_moci_ond_fecund ~ normal(-0.25, 0.15);  // negative: warm fall → lower fecundity
  // OND → pup survival: post-weaning foraging (fall of birth year = moci_ond[t])
  // Separate from fecundity effect; affects survival of t-1 cohort, not production of t cohort
  beta_moci_ond_pup    ~ normal(-0.20, 0.15);  // negative: warm fall → poorer post-weaning prey
  beta_moci_amj_pup   ~ normal(-0.15, 0.15);
  // JFM → pup survival: first winter (moci_jfm[t])
  beta_moci_jfm_pup    ~ normal(-0.15, 0.15);  // negative: warm winter → lower prey availability
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

        // Pup survival: AMJ + OND (post-weaning) + JFM (first winter)
        // All MOCI effects use the same scenario MOCI value
        real pp = inv_logit(phi_pup_logit + site_effect[s] + ce +
                            beta_moci_amj_pup * moci_proj[scen,tp] +
                            beta_moci_ond_pup  * moci_proj[scen,tp] +
                            beta_moci_jfm_pup  * moci_proj[scen,tp]);
        real pj = inv_logit(logit(phi_juv_base) + site_effect[s]*0.5 +
                            beta_moci_jfm_juv * moci_proj[scen,tp]);
        // Adult: sex-specific
        real paF = inv_logit(phi_adult_F_logit + site_effect[s]*0.25 +
                             beta_moci_jfm_adult * moci_proj[scen,tp]);
        real paM = inv_logit(logit(phi_adult_M_base) + site_effect[s]*0.25 +
                             beta_moci_jfm_adult * moci_proj[scen,tp]);

        // Fecundity: modulated by fall MOCI (maternal energy at conception)
        real fecund_proj = inv_logit(logit(avg_fecundity) +
                                     beta_moci_ond_fecund * moci_proj[scen,tp]);
        real np = pAF[s,tp-1] * fecund_proj;
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
    # True values set to literature-informed estimates matching updated priors
    phi_pup_base      = 0.23,            # literature central: Hansen 6-mo²=0.15; Oates 6-mo²=0.23
    phi_pup_logit     = qlogis(0.23),    # ≈ -1.27; prior Normal(-1.2, 0.5)
    phi_juv_base      = 0.80,            # Hastings 2012 (P. v. richardii) sex-neutral ≈ 0.82; prior beta(16,4)
    phi_adult_F_logit = qlogis(0.90),    # = 2.197; Manugian 2017 Tomales Bay; prior Normal(2.20, 0.25)
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
    
    detect_breed_logit = 1.20, detect_molt_logit = 0.75,  # plogis ≈ 0.77, 0.68
    beta_moci_ond_fecund = -0.25,
    beta_moci_ond_pup    = -0.20, beta_moci_jfm_pup   = -0.15,
    beta_moci_amj_pup    = -0.15,
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
  T_inflect <- 14   # peak year; adjust to move inflection point
  
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
      # Covariates acting during pup birth year (t-1); use t for t=1 (harmless)
      t_birth <- if (t > 1) t - 1 else 1
      ce  <- if (coyote_idx[s]>0) true_params$beta_coy[coyote_idx[s]]*coyote[s,t_birth] else 0
      dse <- true_params$beta_dist_surv[s]   * disturbance[s,t_birth]
      dde <- true_params$beta_dist_detect[s] * disturbance[s,t]   # detection at counting year
      
      # Pup survival: birth-year covariates at t_birth, MOCI OND + JFM at t
      phi_pup[s,t] <- plogis(true_params$phi_pup_logit + site_effect[s] + ce +
                               true_params$beta_moci_amj_pup*moci_amj[t_birth] +
                               true_params$beta_moci_ond_pup*moci_ond[t] +
                               true_params$beta_moci_jfm_pup*moci_jfm[t] +
                               has_eseal[s]*true_params$beta_eseal_pup*elephant_seal[s,t_birth] + dse)
      
      # Juvenile survival: sex-neutral
      phi_juv[s,t] <- plogis(qlogis(true_params$phi_juv_base) + site_effect[s]*0.5 +
                               true_params$beta_moci_jfm_juv*moci_jfm[t])
      
      # Adult survival: sex-specific
      phi_adult_F[s,t] <- plogis(qlogis(true_params$phi_adult_F_base) + site_effect[s]*0.25 +
                                   true_params$beta_moci_jfm_adult*moci_jfm[t])
      phi_adult_M[s,t] <- plogis(qlogis(true_params$phi_adult_M_base) + site_effect[s]*0.25 +
                                   true_params$beta_moci_jfm_adult*moci_jfm[t])
      
      detect_breed[s,t] <- plogis(true_params$detect_breed_logit + dde)
      detect_molt[s,t]  <- plogis(true_params$detect_molt_logit  + dde +
                                    true_params$beta_moci_amj_molt*moci_amj[t])
      
      if (t > 1) {
        # Fecundity modulated by fall MOCI (maternal energy at conception)
        fecund_t <- plogis(qlogis(true_params$avg_fecundity) +
                             true_params$beta_moci_ond_fecund * moci_ond[t])
        ep  <- N_adult_F[s,t-1] * fecund_t
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
# PART 4: MAIN EXECUTION
# ============================================================================
# Plotting, tables, and post-processing functions live in the companion script.
# That script is sourced automatically here so run_full_analysis_v3.2() works
# as a single call, but each function can also be called individually after
# loading a saved fit — see harbor_seal_ipm_v3.2_plots.R for details.

# Pipe-safe null-coalescing operator (avoids rlang dependency)
`%||%` <- function(x, y) if (!is.null(x)) x else y

run_full_analysis_v3.2 <- function(use_real_data  = FALSE,
                                   dat            = NULL,
                                   cov_t_scaled   = NULL,
                                   years          = NULL,
                                   T_proj         = 10,
                                   seed           = 42,
                                   iter_warmup    = 3000,
                                   iter_sampling  = 1000,
                                   adapt_delta    = 0.995,
                                   max_treedepth  = 12,
                                   run_portfolio  = FALSE,
                                   run_synchrony  = FALSE) {
  
  # Source companion plots script if functions not yet loaded
  if (!exists("theme_seal", mode="function")) {
    plots_script <- file.path(dirname(sys.frame(1)$ofile %||% "."),
                              "harbor_seal_ipm_v3.2_plots.R")
    if (!file.exists(plots_script))
      plots_script <- "harbor_seal_ipm_v3.2_plots.R"
    source(plots_script)
    cat("Sourced: harbor_seal_ipm_v3.2_plots.R\n")
  }
  
  cat("\n================================================================\n")
  cat("   HARBOR SEAL IPM v3.2\n")
  cat("   Sex-neutral pup/juv survival | Corrected pup prior\n")
  cat("   Explicit observation sex structure\n")
  cat("================================================================\n\n")
  
  prefix <- ifelse(use_real_data, "IPM_v3.2_real", "IPM_v3.2_sim")
  
  if (use_real_data) {
    dl       <- prepare_real_data_for_ipm_v3.2(dat, cov_t_scaled, years, T_proj)
    sim_data <- list(stan_data      = dl$stan_data,
                     site_names     = dl$site_names,
                     years          = dl$years,
                     scenario_names = dl$scenario_names,
                     true_params    = NULL)
  } else {
    sim_data <- simulate_seal_ipm_data_v3.2(T=29, S=6, T_proj=T_proj, seed=seed)
  }
  
  cat("Compiling Stan model...\n")
  model <- cmdstan_model("harbor_seal_ipm_v3.2.stan")
  
  cat(sprintf("Running MCMC (warmup=%d, sampling=%d, adapt_delta=%.3f)...\n",
              iter_warmup, iter_sampling, adapt_delta))
  fit <- model$sample(
    data            = sim_data$stan_data,
    seed            = 123,
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = iter_warmup,
    iter_sampling   = iter_sampling,
    refresh         = 200,
    adapt_delta     = adapt_delta,
    max_treedepth   = max_treedepth
  )
  fit$save_object(paste0("Output/harbor_seal_", prefix, "_fit.rds"))
  
  # Save sim_data alongside fit so plots can be rerun without re-sampling
  saveRDS(sim_data, paste0("Output/harbor_seal_", prefix, "_input_data.rds"))
  cat(sprintf("Input data saved: Output/harbor_seal_%s_input_data.rds\n", prefix))
  
  # Run all plots and tables via companion script orchestrator
  results <- run_all_plots_v3.2(
    fit           = fit,
    sim_data      = sim_data,
    prefix        = prefix,
    run_recovery  = !use_real_data && !is.null(sim_data$true_params),
    run_portfolio = run_portfolio,
    run_synchrony = run_synchrony
  )
  
  cat("\n================================================================\n")
  cat("   COMPLETE — IPM v3.2\n")
  cat(sprintf("   Plots  -> Output/Plots/%s_*.jpeg\n", prefix))
  cat(sprintf("   Fit    -> Output/harbor_seal_%s_fit.rds\n", prefix))
  cat(sprintf("   Table  -> Output/%s_parameter_summary.csv\n", prefix))
  cat("================================================================\n\n")
  
  c(list(fit=fit, model=model, data=sim_data, prefix=prefix), results)
}


# ============================================================================
# USAGE
# ============================================================================

cat("\n")
cat("============================================\n")
cat("IPM v3.2 MODEL SCRIPT LOADED\n")
cat("============================================\n")
cat("  Parts 1-3: Stan model, simulate data, prepare real data\n")
cat("  Plots/tables: source harbor_seal_ipm_v3.2_plots.R\n")
cat("\nTo run full pipeline:\n")
cat("  results <- run_full_analysis_v3.2(\n")
cat("    use_real_data=TRUE, dat=dat, cov_t_scaled=cov_t_scaled, years=years,\n")
cat("    iter_warmup=5000, iter_sampling=1000, adapt_delta=0.99\n")
cat("  )\n")
cat("\nTo replot from a saved fit:\n")
cat("  source(\'harbor_seal_ipm_v3.2_plots.R\')\n")
cat("  out <- load_seal_results(\'IPM_v3.2_real\')\n")
cat("  run_all_plots_v3.2(out$fit, out$sim_data, prefix=\'IPM_v3.2_real\')\n")
cat("\nSite indexing: 1=BL 2=DE 3=DP 4=PRH 5=TB 6=TP\n")
cat("============================================\n")
