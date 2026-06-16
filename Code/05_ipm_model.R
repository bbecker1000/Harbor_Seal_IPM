# ============================================================================
# 05_ipm_model.R
# ----------------------------------------------------------------------------
# HARBOR SEAL INTEGRATED POPULATION MODEL v3.2  —  STAN MODEL + DATA + RUNNER
#
# Contains:
#   PART 1  Stan model (v3.2, performance-optimized; results-identical)
#   PART 2  simulate_seal_ipm_data_v3.2()      — simulated data generator
#   PART 3  prepare_real_data_for_ipm_v3.2()   — real-data assembler (+ guards)
#   PART 4  run_full_analysis_v3.2()           — SINGLE orchestrator
#
# The plotting/orchestration helpers live in 06_ipm_plots.R, which this file
# sources automatically when run_full_analysis_v3.2() needs them.
#
# OPTIMIZATIONS vs the original Stan (all results-identical; see notes in chat):
#   (a) Loop-invariant logit() calls hoisted out of the s,t loops.
#   (b) Likelihood vectorized via observed-index arrays built in transformed
#       data (one normal_lpdf per class instead of S*T branched scalar calls).
#   (c) Output-size reduction: phi_juv / phi_adult_F / phi_adult_M are computed
#       in a LOCAL block (not written to every draw); their site-means are
#       exported as small vectors instead. Shrinks the fit RDS substantially.
#
# SCIENCE UNCHANGED: priors, fecundity parameterization, and simulation truth
# are exactly as before. Two flagged-but-unchanged items for the manuscript:
#   NOTE[P1] fecundity primip/mature are only jointly identified via avg; the
#            split is effectively prior-driven (see comment at fecund priors).
#   NOTE[P2] the simulation sets each true value = its prior mean, so parameter
#            "recovery" over-covers; revisit before claiming 100% recovery.
# ============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)

dir.create("Output",       showWarnings = FALSE)
dir.create("Output/Plots", showWarnings = FALSE)

# Pipe-safe null-coalescing operator (avoids rlang dependency)
`%||%` <- function(x, y) if (!is.null(x)) x else y

# ============================================================================
# SITE INDEXING
#   1=BL 2=DE 3=DP 4=PRH 5=TB 6=TP
#   Coyote: BL,DE,DP   Elephant seal: DE,PRH   Disturbance: all 6
# ============================================================================

# ============================================================================
# PART 1: STAN MODEL (v3.2 — optimized)
# ============================================================================

stan_code_v3.2 <- '
// ============================================================================
// HARBOR SEAL IPM v3.2 — STAN MODEL (performance-optimized, results-identical)
// ============================================================================

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

transformed data {
  // ── Precompute observed-cell index arrays for a vectorized likelihood ─────
  // (Optimization b) One pass to count, one to fill. Replaces S*T branch tests
  // and per-cell normal() calls in the model block.
  int N_oa = 0;
  int N_op = 0;
  int N_om = 0;
  for (s in 1:S) for (t in 1:T) {
    N_oa += y_adult_obs[s, t];
    N_op += y_pup_obs[s, t];
    N_om += y_molt_obs[s, t];
  }
  array[N_oa] int oa_s; array[N_oa] int oa_t; vector[N_oa] y_oa;
  array[N_op] int op_s; array[N_op] int op_t; vector[N_op] y_op;
  array[N_om] int om_s; array[N_om] int om_t; vector[N_om] y_om;
  {
    int ia = 0; int ip = 0; int im = 0;
    for (s in 1:S) for (t in 1:T) {
      if (y_adult_obs[s, t] == 1) { ia += 1; oa_s[ia]=s; oa_t[ia]=t; y_oa[ia]=y_adult[s,t]; }
      if (y_pup_obs[s, t]   == 1) { ip += 1; op_s[ip]=s; op_t[ip]=t; y_op[ip]=y_pup[s,t]; }
      if (y_molt_obs[s, t]  == 1) { im += 1; om_s[im]=s; om_t[im]=t; y_om[im]=y_molt[s,t]; }
    }
  }
}

parameters {
  // ── Survival ──────────────────────────────────────────────────────────────
  real phi_pup_logit;
  real<lower=0, upper=1> phi_juv_base;
  real phi_adult_F_logit;
  real<lower=0, upper=0.10> delta_adult;

  // ── Reproduction ──────────────────────────────────────────────────────────
  real<lower=0, upper=1>    fecund_primip;
  real<lower=0, upper=1>    fecund_mature;
  real<lower=0.4, upper=0.6> prop_female;

  // ── Observation: male haul-out during breeding season ─────────────────────
  real<lower=0, upper=0.30> p_male_breed;

  // ── Site-specific covariate effects ───────────────────────────────────────
  vector[N_coy] beta_coy;
  vector[S]     beta_dist_surv;
  vector[S]     beta_dist_detect;
  real detect_breed_logit;
  real detect_molt_logit;

  // ── Shared covariate effects ───────────────────────────────────────────────
  real beta_moci_ond_fecund;
  real beta_moci_ond_pup;
  real beta_moci_amj_pup;
  real beta_moci_jfm_pup;
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
  // ── Saved scalars ───────────────────────────────────────────────────────
  vector[S] site_effect = sigma_site * site_effect_raw;
  real phi_pup_base     = inv_logit(phi_pup_logit);
  real phi_adult_F_base = inv_logit(phi_adult_F_logit);
  real phi_adult_M_base = fmax(phi_adult_F_base - delta_adult, 0.01);
  real avg_fecundity    = 0.20 * fecund_primip + 0.80 * fecund_mature;

  // ── Saved initial-population vectors ──────────────────────────────────────
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

  // ── SAVED matrices (needed by model block, generated quantities, or plots) ─
  matrix<lower=0>[S, T] N_adult_F;
  matrix<lower=0>[S, T] N_adult_M;
  matrix<lower=0>[S, T] N_juv_F;
  matrix<lower=0>[S, T] N_juv_M;
  matrix<lower=0>[S, T] N_pup;
  matrix<lower=0>[S, T] N_adult_total;
  matrix<lower=0>[S, T] N_juv_total;
  matrix<lower=0>[S, T] N_molt_true;
  matrix<lower=0>[S, T] N_total;
  matrix<lower=0, upper=1>[S, T] phi_pup;        // used by portfolio + GQ
  matrix<lower=0, upper=1>[S, T] detect_breed;   // used by likelihood + GQ
  matrix<lower=0, upper=1>[S, T] detect_molt;    // used by likelihood + GQ

  // ── SAVED site-mean vital-rate vectors (small; replace full phi matrices) ──
  vector[T] mean_phi_pup;
  vector[T] mean_phi_juv;
  vector[T] mean_phi_adult_F;
  vector[T] mean_phi_adult_M;

  {
    // ── LOCAL (NOT saved): per-site/time juv & adult survival ──────────────
    // (Optimization c) These feed the dynamics and the means below, but are
    // never read by any downstream plot, so they are not written per draw.
    matrix[S, T] phi_juv;
    matrix[S, T] phi_adult_F;
    matrix[S, T] phi_adult_M;

    // ── Hoisted loop-invariants (Optimization a) ───────────────────────────
    real logit_phi_juv    = logit(phi_juv_base);
    real logit_phi_adultM = logit(phi_adult_M_base);
    real logit_avg_fec    = logit(avg_fecundity);

    // ── Vital rates ─────────────────────────────────────────────────────────
    for (s in 1:S) {
      for (t in 1:T) {
        int t_birth = (t > 1) ? t - 1 : 1;

        real coyote_effect = 0;
        if (coyote_idx[s] > 0)
          coyote_effect = beta_coy[coyote_idx[s]] * coyote[s, t_birth];

        real dist_surv_eff   = beta_dist_surv[s]   * disturbance[s, t_birth];
        real dist_detect_eff = beta_dist_detect[s] * disturbance[s, t];

        phi_pup[s, t] = inv_logit(
          phi_pup_logit + site_effect[s] + coyote_effect +
          beta_moci_amj_pup * moci_amj[t_birth] +
          beta_moci_ond_pup  * moci_ond[t]       +
          beta_moci_jfm_pup  * moci_jfm[t]       +
          has_eseal[s] * beta_eseal_pup * elephant_seal[s, t_birth] +
          dist_surv_eff);

        phi_juv[s, t] = inv_logit(
          logit_phi_juv + site_effect[s] * 0.5 +
          beta_moci_jfm_juv * moci_jfm[t]);

        phi_adult_F[s, t] = inv_logit(
          phi_adult_F_logit + site_effect[s] * 0.25 +
          beta_moci_jfm_adult * moci_jfm[t]);

        phi_adult_M[s, t] = inv_logit(
          logit_phi_adultM + site_effect[s] * 0.25 +
          beta_moci_jfm_adult * moci_jfm[t]);

        detect_breed[s, t] = inv_logit(detect_breed_logit + dist_detect_eff);
        detect_molt[s, t]  = inv_logit(detect_molt_logit  + dist_detect_eff +
                                        beta_moci_amj_molt * moci_amj[t]);
      }
    }

    // ── Population dynamics ───────────────────────────────────────────────────
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
        real fecund_t      = inv_logit(logit_avg_fec +
                                       beta_moci_ond_fecund * moci_ond[t]);
        real expected_pups = N_adult_F[s, t-1] * fecund_t;

        real new_juv_F = N_pup[s, t-1] * prop_female       * phi_pup[s, t];
        real new_juv_M = N_pup[s, t-1] * (1 - prop_female) * phi_pup[s, t];

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
        N_adult_F[s, t] = exp(log(fmax(N_adult_F[s, t-1] * phi_adult_F[s, t] +
                                       juv_to_adult_F, 1)) +
                              sigma_process * eps_adult_raw[s, t-1] * 0.5);
        N_adult_M[s, t] = exp(log(fmax(N_adult_M[s, t-1] * phi_adult_M[s, t] +
                                       juv_to_adult_M, 1)) +
                              sigma_process * eps_adult_raw[s, t-1] * 0.5);

        N_adult_total[s, t] = N_adult_F[s, t] + N_adult_M[s, t];
        N_juv_total[s, t]   = N_juv_F[s, t]   + N_juv_M[s, t];
        N_molt_true[s, t]   = N_juv_total[s, t] + N_adult_total[s, t];
        N_total[s, t]       = N_pup[s, t] + N_juv_total[s, t] + N_adult_total[s, t];
      }
    }

    // ── Site-mean vital rates (exported; full matrices stay local) ──────────
    for (t in 1:T) {
      mean_phi_pup[t]     = mean(col(phi_pup,     t));
      mean_phi_juv[t]     = mean(col(phi_juv,     t));
      mean_phi_adult_F[t] = mean(col(phi_adult_F, t));
      mean_phi_adult_M[t] = mean(col(phi_adult_M, t));
    }
  }
}

model {
  // ── Priors (UNCHANGED from the fitted model) ─────────────────────────────
  phi_pup_logit     ~ normal(-1.2, 0.5);
  phi_juv_base      ~ beta(14, 6);          // mean 0.70 (declining-pop literature)
  phi_adult_F_logit ~ normal(2.20, 0.25);
  delta_adult       ~ normal(0.05, 0.025);

  // NOTE[P1]: fecund_primip and fecund_mature are only jointly identified via
  // avg_fecundity (0.2*primip + 0.8*mature); the data pin the average and the
  // 80% weight pulls fecund_mature toward it, leaving fecund_primip prior-led.
  // Report avg_fecundity as the estimated quantity, or fix the primip:mature
  // ratio, before interpreting the two levels separately in the manuscript.
  fecund_primip ~ beta(12, 8);    // mean 0.60
  fecund_mature ~ beta(17, 3);    // mean 0.85
  prop_female   ~ beta(50, 50);

  p_male_breed ~ beta(2, 18);     // mean 0.10; aquatic mating system

  // Common weakly-informative covariate priors (data differentiate sites).
  beta_coy         ~ normal(-0.20, 0.20);
  beta_dist_surv   ~ normal(-0.15, 0.20);
  beta_dist_detect ~ normal(-0.15, 0.15);

  detect_breed_logit ~ normal(1.20, 0.50);
  detect_molt_logit  ~ normal(0.75, 0.50);

  beta_moci_ond_fecund ~ normal(-0.15, 0.20);
  beta_moci_ond_pup    ~ normal(-0.15, 0.20);
  beta_moci_amj_pup    ~ normal(-0.15, 0.20);
  beta_moci_jfm_pup    ~ normal(-0.15, 0.20);
  beta_moci_jfm_juv    ~ normal(-0.15, 0.15);
  beta_moci_jfm_adult  ~ normal(-0.10, 0.12);
  beta_moci_amj_molt   ~ normal(0.05, 0.15);

  beta_eseal_pup ~ normal(0.10, 0.20);

  sigma_site      ~ normal(0.2,  0.1);
  site_effect_raw ~ std_normal();
  sigma_process   ~ normal(0.15, 0.08);
  sigma_obs_adult ~ normal(0.18, 0.06);
  sigma_obs_pup   ~ normal(0.15, 0.02);   // tight: resolves sigma_proc ridge
  sigma_obs_molt  ~ normal(0.35, 0.10);

  mu_log_adult ~ normal(5, 0.5);
  mu_log_juv   ~ normal(4, 0.5);
  mu_log_pup   ~ normal(4, 0.5);
  sigma_init   ~ exponential(3);
  log_N_adult_F_init_raw ~ std_normal();
  log_N_adult_M_init_raw ~ std_normal();
  log_N_juv_init_raw     ~ std_normal();
  log_N_pup_init_raw     ~ std_normal();

  to_vector(eps_adult_raw) ~ std_normal();
  to_vector(eps_juv_raw)   ~ std_normal();
  to_vector(eps_pup_raw)   ~ std_normal();

  // ── Vectorized likelihood (Optimization b) ───────────────────────────────
  {
    vector[N_oa] mu_a;
    vector[N_op] mu_p;
    vector[N_om] mu_m;
    for (i in 1:N_oa)
      mu_a[i] = log((N_adult_F[oa_s[i], oa_t[i]]
                     + N_adult_M[oa_s[i], oa_t[i]] * p_male_breed)
                    * detect_breed[oa_s[i], oa_t[i]]);
    for (i in 1:N_op)
      mu_p[i] = log(N_pup[op_s[i], op_t[i]] * detect_breed[op_s[i], op_t[i]]);
    for (i in 1:N_om)
      mu_m[i] = log(N_molt_true[om_s[i], om_t[i]] * detect_molt[om_s[i], om_t[i]]);

    y_oa ~ normal(mu_a, sigma_obs_adult);
    y_op ~ normal(mu_p, sigma_obs_pup);
    y_om ~ normal(mu_m, sigma_obs_molt);
  }
}

generated quantities {
  // ── Posterior predictive (all cells, observed or not) ────────────────────
  matrix[S, T] y_adult_rep;
  matrix[S, T] y_pup_rep;
  matrix[S, T] y_molt_rep;

  // ── Derived quantities ─────────────────────────────────────────────────────
  matrix[S, T]   sex_ratio_adult;
  matrix[S, T]   sex_ratio_observed;
  matrix[S, T-1] lambda;
  vector[T]      N_total_all;

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

  for (s in 1:S)
    for (t in 1:(T-1))
      lambda[s, t] = N_total[s, t+1] / N_total[s, t];

  for (t in 1:T)
    N_total_all[t] = sum(col(N_total, t));

  // ── Scenario projections (recompute vital rates from parameters) ──────────
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

        // NOTE: scenario MOCI anomaly applied uniformly across the 3 pup
        // seasonal channels (sum of slopes) -> larger effective sensitivity
        // than any single coefficient. Document in Methods (manuscript P5).
        real pp = inv_logit(phi_pup_logit + site_effect[s] + ce +
                            beta_moci_amj_pup * moci_proj[scen,tp] +
                            beta_moci_ond_pup  * moci_proj[scen,tp] +
                            beta_moci_jfm_pup  * moci_proj[scen,tp]);
        real pj = inv_logit(logit(phi_juv_base) + site_effect[s]*0.5 +
                            beta_moci_jfm_juv * moci_proj[scen,tp]);
        real paF = inv_logit(phi_adult_F_logit + site_effect[s]*0.25 +
                             beta_moci_jfm_adult * moci_proj[scen,tp]);
        real paM = inv_logit(logit(phi_adult_M_base) + site_effect[s]*0.25 +
                             beta_moci_jfm_adult * moci_proj[scen,tp]);

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

write_lines(stan_code_v3.2, "Code/harbor_seal_ipm_v3.2.stan")
cat("Stan model written to Code/harbor_seal_ipm_v3.2.stan\n")


# ============================================================================
# PART 2: SIMULATE DATA
# ----------------------------------------------------------------------------
# NOTE[P2]: every true_params value below is set to its prior MEAN. This makes
# parameter "recovery" over-cover (truth always near the prior centre). For a
# genuine recovery / SBC check, draw true values FROM the priors (or offset
# them) instead. Left as-is here to reproduce existing results; revisit before
# claiming "100% recovery" in the manuscript.
# ============================================================================

simulate_seal_ipm_data_v3.2 <- function(T=29, S=6, T_proj=10, seed=123) {
  
  set.seed(seed)
  site_names <- c("BL","DE","DP","PRH","TB","TP")
  N_coy <- 3
  
  true_params <- list(
    phi_pup_logit     = qlogis(0.23),
    phi_juv_base      = 0.70,
    phi_adult_F_logit = qlogis(0.90),
    delta_adult       = 0.05,
    fecund_primip     = 0.60,
    fecund_mature     = 0.85,
    prop_female       = 0.50,
    p_male_breed      = 0.10,
    detect_breed_logit = 1.20,
    detect_molt_logit  = 0.75,
    beta_coy = c(-0.20, -0.20, -0.20),
    beta_dist_surv   = rep(-0.15, 6),
    beta_dist_detect = rep(-0.15, 6),
    beta_moci_ond_fecund = -0.15,
    beta_moci_ond_pup    = -0.15,
    beta_moci_amj_pup    = -0.15,
    beta_moci_jfm_pup    = -0.15,
    beta_moci_jfm_juv    = -0.15,
    beta_moci_jfm_adult  = -0.10,
    beta_moci_amj_molt   =  0.05,
    beta_eseal_pup = 0.10,
    sigma_process   = 0.15,
    sigma_obs_adult = 0.18,
    sigma_obs_pup   = 0.15,
    sigma_obs_molt  = 0.35,
    sigma_site      = 0.20
  )
  
  phi_adult_M_base <- plogis(true_params$phi_adult_F_logit) - true_params$delta_adult
  avg_fecundity    <- 0.20 * true_params$fecund_primip + 0.80 * true_params$fecund_mature
  
  coyote_idx <- c(1,2,3,0,0,0)
  has_eseal  <- c(0,1,0,1,0,0)
  
  T_inflect <- 14
  moci_mean <- c(rep(-0.6, T_inflect - 2),
                 seq(-0.6, 0.6, length.out = 4),
                 rep( 0.6, T - T_inflect - 2))
  moci_base <- moci_mean + as.vector(arima.sim(list(ar=0.45), n=T)) * 0.7
  moci_jfm  <- as.vector(scale(moci_base))
  moci_amj  <- as.vector(scale(moci_base * 0.85 +
                                 as.vector(arima.sim(list(ar=0.3), n=T)) * 0.4))
  moci_ond  <- as.vector(scale(moci_base * 0.75 +
                                 as.vector(arima.sim(list(ar=0.3), n=T)) * 0.5))
  
  coyote_trend <- c(rep(0, T_inflect),
                    seq(0, 1.2, length.out = T - T_inflect))
  coyote <- matrix(0, S, T)
  coyote[1,] <- as.vector(scale(coyote_trend +
                                  as.vector(arima.sim(list(ar=0.4), n=T)) * 0.35))
  coyote[2,] <- as.vector(scale(coyote_trend * 1.3 +
                                  as.vector(arima.sim(list(ar=0.4), n=T)) * 0.35))
  coyote[3,] <- as.vector(scale(coyote_trend * 0.5 +
                                  as.vector(arima.sim(list(ar=0.4), n=T)) * 0.35))
  
  dist_trend  <- seq(0, 0.4, length.out = T)
  disturbance <- matrix(0, S, T)
  for (s in 1:S)
    disturbance[s,] <- as.vector(scale(dist_trend +
                                         as.vector(arima.sim(list(ar=0.3), n=T)) * 0.6))
  
  elephant_seal    <- matrix(0, S, T)
  elephant_seal[2,] <- as.vector(scale(seq(0, 3, length.out=T) + rnorm(T, 0, 0.5)))
  elephant_seal[4,] <- as.vector(scale(seq(0, 4, length.out=T) + rnorm(T, 0, 0.5)))
  
  site_effect <- rnorm(S, 0, true_params$sigma_site)
  
  N_adult_F <- N_adult_M <- matrix(NA, S, T)
  N_juv_F   <- N_juv_M   <- matrix(NA, S, T)
  N_pup                  <- matrix(NA, S, T)
  phi_pup   <- phi_juv   <- matrix(NA, S, T)
  phi_adult_F <- phi_adult_M <- matrix(NA, S, T)
  detect_breed <- detect_molt <- matrix(NA, S, T)
  
  N_adult_F[,1] <- c(120, 90, 45, 65, 75, 28)
  N_adult_M[,1] <- c(105, 80, 40, 58, 67, 25)
  N_juv_F[,1]   <- c( 38, 28, 14, 19, 24,  9)
  N_juv_M[,1]   <- c( 34, 25, 12, 17, 21,  8)
  N_pup[,1]     <- c( 90, 68, 34, 46, 56, 22)
  
  for (s in 1:S) {
    for (t in 1:T) {
      t_birth <- if (t > 1) t - 1 else 1
      
      ce  <- if (coyote_idx[s] > 0)
        true_params$beta_coy[coyote_idx[s]] * coyote[s, t_birth] else 0
      dse <- true_params$beta_dist_surv[s]   * disturbance[s, t_birth]
      dde <- true_params$beta_dist_detect[s] * disturbance[s, t]
      
      phi_pup[s,t] <- plogis(
        true_params$phi_pup_logit + site_effect[s] + ce +
          true_params$beta_moci_amj_pup * moci_amj[t_birth] +
          true_params$beta_moci_ond_pup * moci_ond[t]       +
          true_params$beta_moci_jfm_pup * moci_jfm[t]       +
          has_eseal[s] * true_params$beta_eseal_pup * elephant_seal[s, t_birth] + dse)
      
      phi_juv[s,t] <- plogis(
        qlogis(true_params$phi_juv_base) + site_effect[s] * 0.5 +
          true_params$beta_moci_jfm_juv * moci_jfm[t])
      
      phi_adult_F[s,t] <- plogis(
        true_params$phi_adult_F_logit + site_effect[s] * 0.25 +
          true_params$beta_moci_jfm_adult * moci_jfm[t])
      
      phi_adult_M[s,t] <- plogis(
        qlogis(phi_adult_M_base) + site_effect[s] * 0.25 +
          true_params$beta_moci_jfm_adult * moci_jfm[t])
      
      detect_breed[s,t] <- plogis(true_params$detect_breed_logit + dde)
      detect_molt[s,t]  <- plogis(
        true_params$detect_molt_logit + dde +
          true_params$beta_moci_amj_molt * moci_amj[t])
      
      if (t > 1) {
        fecund_t <- plogis(
          qlogis(avg_fecundity) +
            true_params$beta_moci_ond_fecund * moci_ond[t])
        
        ep  <- N_adult_F[s, t-1] * fecund_t
        njF <- N_pup[s, t-1] * true_params$prop_female       * phi_pup[s, t]
        njM <- N_pup[s, t-1] * (1 - true_params$prop_female) * phi_pup[s, t]
        jsF <- N_juv_F[s, t-1] * phi_juv[s, t] * (2/3)
        jsM <- N_juv_M[s, t-1] * phi_juv[s, t] * (2/3)
        jaF <- N_juv_F[s, t-1] * phi_juv[s, t] * (1/3)
        jaM <- N_juv_M[s, t-1] * phi_juv[s, t] * (1/3)
        
        N_pup[s,t]     <- exp(rnorm(1, log(max(ep,         1.0)),
                                    true_params$sigma_process))
        N_juv_F[s,t]   <- exp(rnorm(1, log(max(njF + jsF,  0.1)),
                                    true_params$sigma_process * 0.5))
        N_juv_M[s,t]   <- exp(rnorm(1, log(max(njM + jsM,  0.1)),
                                    true_params$sigma_process * 0.5))
        N_adult_F[s,t] <- exp(rnorm(1,
                                    log(max(N_adult_F[s,t-1] * phi_adult_F[s,t] + jaF, 1)),
                                    true_params$sigma_process * 0.5))
        N_adult_M[s,t] <- exp(rnorm(1,
                                    log(max(N_adult_M[s,t-1] * phi_adult_M[s,t] + jaM, 1)),
                                    true_params$sigma_process * 0.5))
      }
    }
  }
  
  N_adult_total <- N_adult_F + N_adult_M
  N_juv_total   <- N_juv_F   + N_juv_M
  N_molt_true   <- N_juv_total + N_adult_total
  
  y_adult <- y_pup <- y_molt <- matrix(NA, S, T)
  for (s in 1:S) {
    for (t in 1:T) {
      N_adult_obs  <- N_adult_F[s,t] + N_adult_M[s,t] * true_params$p_male_breed
      y_adult[s,t] <- log(N_adult_obs      * detect_breed[s,t]) +
        rnorm(1, 0, true_params$sigma_obs_adult)
      y_pup[s,t]   <- log(N_pup[s,t]       * detect_breed[s,t]) +
        rnorm(1, 0, true_params$sigma_obs_pup)
      y_molt[s,t]  <- log(N_molt_true[s,t] * detect_molt[s,t])  +
        rnorm(1, 0, true_params$sigma_obs_molt)
    }
  }
  
  n_obs <- S * T
  y_adult[sample(1:n_obs, round(0.05 * n_obs))] <- NA
  y_pup[  sample(1:n_obs, round(0.05 * n_obs))] <- NA
  y_molt[ sample(1:n_obs, round(0.05 * n_obs))] <- NA
  
  y_adult_obs <- ifelse(is.na(y_adult), 0L, 1L)
  y_pup_obs   <- ifelse(is.na(y_pup),   0L, 1L)
  y_molt_obs  <- ifelse(is.na(y_molt),  0L, 1L)
  y_adult[is.na(y_adult)] <- 0
  y_pup[  is.na(y_pup)]   <- 0
  y_molt[ is.na(y_molt)]  <- 0
  
  N_scenarios    <- 4
  moci_proj      <- matrix(c( 0, 1, -1, 1), N_scenarios, T_proj)
  coyote_proj    <- matrix(c( 0, 0,  0, 1), N_scenarios, T_proj)
  scenario_names <- c("Status Quo", "Warm (MOCI +1)",
                      "Cool (MOCI -1)", "Warm + High Coyote")
  
  stan_data <- list(
    T = T, S = S, N_coy = N_coy,
    y_adult     = y_adult,     y_pup     = y_pup,     y_molt     = y_molt,
    y_adult_obs = y_adult_obs, y_pup_obs = y_pup_obs, y_molt_obs = y_molt_obs,
    coyote      = coyote,      disturbance = disturbance,
    elephant_seal = elephant_seal,
    moci_jfm    = as.vector(moci_jfm),
    moci_amj    = as.vector(moci_amj),
    moci_ond    = as.vector(moci_ond),
    coyote_idx  = coyote_idx,
    has_eseal   = has_eseal,
    T_proj      = T_proj,
    N_scenarios = N_scenarios,
    moci_proj   = moci_proj,
    coyote_proj = coyote_proj
  )
  
  list(
    stan_data   = stan_data,
    true_params = true_params,
    true_states = list(
      N_adult_F = N_adult_F, N_adult_M = N_adult_M,
      N_juv_F = N_juv_F, N_juv_M = N_juv_M, N_pup = N_pup,
      N_adult_total = N_adult_total, N_molt_true = N_molt_true,
      phi_pup = phi_pup, phi_juv = phi_juv,
      phi_adult_F = phi_adult_F, phi_adult_M = phi_adult_M,
      detect_breed = detect_breed, detect_molt = detect_molt
    ),
    site_names     = site_names,
    years          = 1997:(1997 + T - 1),
    scenario_names = scenario_names
  )
}

# ============================================================================
# PART 3: PREPARE REAL DATA  (+ canonical-order guards, fixes P9)
# ============================================================================

prepare_real_data_for_ipm_v3.2 <- function(dat, cov_t_scaled, years, T_proj=10) {
  
  site_names <- c("BL","DE","DP","PRH","TB","TP")
  S <- 6; T <- length(years); N_coy <- 3
  
  # ── Guard: shapes and canonical row order (fixes P9; matches 01/02 output) ─
  canon_state <- as.vector(t(outer(site_names, c("ADULT","MOLTING","PUP"),
                                   paste, sep = "_")))
  canon_cov   <- c("MOCI_JFM","MOCI_AMJ","MOCI_OND",
                   paste0("Dist_",   site_names),
                   paste0("Coyote_", site_names),
                   "eSeal_Sum_Imm_MaxCount")
  if (nrow(dat) != 18L)
    stop("`dat` must have 18 rows (6 sites x 3 classes); got ", nrow(dat),
         ". Did you pass an 8-site matrix?")
  if (nrow(cov_t_scaled) != 16L)
    stop("`cov_t_scaled` must have 16 rows; got ", nrow(cov_t_scaled),
         ". Did you pass the 20-row 8-site covariate matrix?")
  if (!is.null(rownames(dat)) && !identical(rownames(dat), canon_state))
    stop("`dat` row order is not canonical. Expected:\n  ",
         paste(canon_state, collapse = ", "))
  if (!is.null(rownames(cov_t_scaled)) &&
      !identical(rownames(cov_t_scaled), canon_cov))
    stop("`cov_t_scaled` row order is not canonical. Expected:\n  ",
         paste(canon_cov, collapse = ", "))
  if (ncol(dat) != T || ncol(cov_t_scaled) != T)
    stop("Column counts must equal length(years) = ", T)
  
  # ── Split count rows by class (canonical interleave: A,M,P per site) ──────
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
  for (s in 1:S) disturbance[s,] <- as.vector(cov_t_scaled[3+s,])   # rows 4-9
  
  coyote <- matrix(0,S,T)
  for (k in 1:3) coyote[k,] <- as.vector(cov_t_scaled[9+k,])        # rows 10-12
  
  elephant_seal <- matrix(0,S,T)
  elephant_seal[2,] <- as.vector(cov_t_scaled[16,])   # DE  (shared eSeal index)
  elephant_seal[4,] <- as.vector(cov_t_scaled[16,])   # PRH
  
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
# PART 4: MAIN ORCHESTRATOR  (single definition — fixes P6)
# ----------------------------------------------------------------------------
# Returns list(fit, model, sim_data, prefix, ...plot results). The element is
# `sim_data` (NOT `data`) so the return shape matches load_seal_results()
# (fixes P13): every downstream call uses out$fit and out$sim_data.
# ============================================================================

run_full_analysis_v3.2 <- function(use_real_data  = FALSE,
                                   dat            = NULL,
                                   cov_t_scaled   = NULL,
                                   years          = NULL,
                                   T_proj         = 10,
                                   seed           = 123,
                                   iter_warmup    = 2000,
                                   iter_sampling  = 2000,
                                   adapt_delta    = 0.97,
                                   max_treedepth  = 12,
                                   run_portfolio  = FALSE,
                                   run_synchrony  = FALSE) {
  
  # Source companion plots script if its functions aren't loaded yet.
  if (!exists("run_all_plots_v3.2", mode = "function")) {
    plots_script <- "Code/06_ipm_plots.R"
    if (!file.exists(plots_script)) plots_script <- "06_ipm_plots.R"
    source(plots_script)
    cat("Sourced:", plots_script, "\n")
  }
  
  cat("\n================================================================\n")
  cat("   HARBOR SEAL IPM v3.2\n")
  cat("================================================================\n\n")
  
  prefix <- ifelse(use_real_data, "IPM_v3.2_real", "IPM_v3.2_sim")
  
  if (use_real_data) {
    if (is.null(dat) || is.null(cov_t_scaled) || is.null(years))
      stop("use_real_data=TRUE requires dat, cov_t_scaled, and years.")
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
  model <- cmdstan_model("Code/harbor_seal_ipm_v3.2.stan")
  
  cat(sprintf("Running MCMC (warmup=%d, sampling=%d, adapt_delta=%.3f)...\n",
              iter_warmup, iter_sampling, adapt_delta))
  fit <- model$sample(
    data            = sim_data$stan_data,
    seed            = seed,
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = iter_warmup,
    iter_sampling   = iter_sampling,
    refresh         = 200,
    adapt_delta     = adapt_delta,
    max_treedepth   = max_treedepth
  )
  
  fit$save_object(paste0("Output/harbor_seal_", prefix, "_fit.rds"))
  
  # Save the model INPUTS (fixes P7: use sim_data$... not undefined bare names)
  saveRDS(list(stan_data   = sim_data$stan_data,
               years       = sim_data$years,
               true_params = sim_data$true_params),
          paste0("Output/harbor_seal_", prefix, "_input_data.rds"))
  cat(sprintf("Input data saved: Output/harbor_seal_%s_input_data.rds\n", prefix))
  
  # Run plots/tables via the companion orchestrator.
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
  cat(sprintf("   Fit -> Output/harbor_seal_%s_fit.rds\n", prefix))
  cat("================================================================\n\n")
  
  # Return shape exposes $fit + $sim_data (matches load_seal_results) — P13.
  c(list(fit = fit, model = model, sim_data = sim_data, prefix = prefix),
    results)
}

cat("\n05_ipm_model.R loaded: Stan model written, functions defined.\n")
cat("Run via 07_ipm_run.R, or call run_full_analysis_v3.2() directly.\n")
