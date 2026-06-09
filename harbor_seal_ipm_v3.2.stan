
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

