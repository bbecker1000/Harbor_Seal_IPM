# ============================================================================
# 17_regional_ipm_model.R
# ----------------------------------------------------------------------------
# REGIONAL HARBOR SEAL IPM — STAN MODEL + DATA PREP + ORCHESTRATOR
#
# Key differences from Marin IPM (05_ipm_model.R):
#
# POPULATION DYNAMICS:
#   - Leslie matrix applied at the COUNTY level (C=6 county populations)
#   - County random effects on baseline vital rates (partial pooling)
#   - NO coyote, disturbance, or elephant seal covariates
#   - p_male_breed FIXED at 0.057 (Marin IPM posterior mean; not estimated)
#
# MOCI STRUCTURE:
#   - Same 3 seasonal pathways as Marin IPM (JFM, AMJ, OND)
#   - NEW: delta_moci_mouth_fecund, delta_moci_mouth_surv (Bay Mouth vs coast)
#          delta_moci_south_fecund, delta_moci_south_surv (South Bay vs coast)
#     for counties flagged as estuarine (county_is_bay = 1).
#     Positive delta = bay counties less negatively affected by poor upwelling.
#     Prior: N(0, 0.15) — weakly informative, allows learning from data.
#
# OBSERVATION MODEL:
#   - Each site s observes a fraction of county c(s)'s population
#   - Site availability: log_alpha_breed[s], log_alpha_pup[s], log_alpha_molt[s]
#     ~ N(mu_log_alpha[county], sigma_log_alpha) — partial pooling within county
#   - is_breeder[s] = 0 (Type H, Alcatraz): no pup likelihood
#   - rho_pup ~ Beta(2,5): pup molt attendance fraction (unchanged from Marin IPM)
#
# ZERO-COUNT CONVENTION (handled in 15_regional_data_prep.R):
#   - NA ("ND" in raw data) = not surveyed → y_obs = 0, excluded from likelihood
#   - Count = 0 = surveyed, no animals seen → y_obs = 1, y = log(0.5)
#   - log(0.5) lower-bound offset is biologically correct and does not bias
#     county N estimates since zero sites (Alcatraz molt, Point San Pedro pup)
#     sit within multi-site counties anchored by abundant positive observations.
#     The Stan likelihood sees log(0.5) as a regular real-valued observation;
#     no code changes required here relative to the positive-count case.
#
# 2020 LATENT YEAR:
#   - Year index t=16 (2020) has all y_*_obs = 0
#   - Leslie matrix transitions run normally; no likelihood contribution
#   - Uncertainty propagates correctly through the gap
#
# PRIORS: Identical to Marin IPM v3.3 for all shared parameters
#         (vital rates, MOCI effects, error terms, rho_pup)
#         Independent analysis — does NOT use Marin posteriors as priors.
#
# Prereq: source("Code/15_regional_data_prep.R")
#         source("Code/16_regional_simulate.R")  (to test on sim data first)
# ============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)

dir.create("Output",       showWarnings = FALSE)
dir.create("Output/Plots", showWarnings = FALSE)

`%||%` <- function(x, y) if (!is.null(x)) x else y

# ── STAN MODEL ────────────────────────────────────────────────────────────────
stan_code_regional <- '
// ============================================================================
// REGIONAL HARBOR SEAL IPM — SITE-LEVEL MODEL
// Parallel to Marin IPM v3.3 structure, extended to 24 sites with
// county-level vital rate pooling and Bay/Coast MOCI modifier.
//
// KEY CHANGE FROM COUNTY-LEVEL MODEL:
//   N state variables are at the SITE level (S=24), not county level (C=6).
//   This eliminates the N/alpha confound: each site N is directly anchored
//   by that site own counts without partitioning a shared county pool.
//   County totals are computed as sums in generated quantities.
//
// OBSERVATION MODEL (identical to Marin IPM v3.3):
//   y_adult[s,t] ~ Normal(log(N_adult_obs[s,t] * detect_breed[s,t]), sigma_obs_adult)
//   y_pup[s,t]   ~ Normal(log(N_pup[s,t]       * detect_breed[s,t]), sigma_obs_pup)
//   y_molt[s,t]  ~ Normal(log(N_molt_true[s,t]  * detect_molt[s,t]),  sigma_obs_molt)
//
// VITAL RATE POOLING:
//   phi_pup[s,t]     = f(phi_pup_logit + county_effect[county_id[s]] + MOCI)
//   phi_adult_F[s,t] = f(phi_adult_F_logit + county_effect[county_id[s]] * 0.25 + MOCI)
//   County random effects partially pool vital rates across sites within county.
//
// MOCI: same 3-season structure as Marin IPM + Bay/Coast modifier by county_type.
// site_t1[s]: Leslie matrix held flat for t < site_t1[s] (handles Point Arena).
// p_male_breed: FIXED at 0.057 (Marin IPM posterior).
// ============================================================================
data {
  int<lower=1> T;           // years (21: 2005-2025; index 16 = 2020 latent)
  int<lower=1> S;           // sites (24)
  int<lower=1> C;           // counties (6)
  int<lower=1> T_proj;
  int<lower=1> N_scenarios;

  // Observations (log-scale; 0 where unobserved; log(0.5) for true zeros)
  matrix[S, T] y_adult;
  matrix[S, T] y_pup;
  matrix[S, T] y_molt;
  array[S, T] int<lower=0,upper=1> y_adult_obs;
  array[S, T] int<lower=0,upper=1> y_pup_obs;
  array[S, T] int<lower=0,upper=1> y_molt_obs;

  // Site classifiers
  array[S] int<lower=1,upper=6>  county_id;  // county membership
  array[S] int<lower=1,upper=20> site_t1;    // first time index with data per site

  // County type: 0=coast, 1=bay mouth, 2=south bay (for MOCI modifier)
  array[C] int<lower=0,upper=2> county_type;

  // MOCI (z-scored)
  vector[T] moci_jfm;
  vector[T] moci_amj;
  vector[T] moci_ond;

  // Projections
  matrix[N_scenarios, T_proj] moci_proj;

  real<lower=0,upper=1> p_male_fixed;
}

transformed data {
  int N_oa = 0; int N_op = 0; int N_om = 0;
  for (s in 1:S) for (t in 1:T) {
    N_oa += y_adult_obs[s,t];
    N_op += y_pup_obs[s,t];
    N_om += y_molt_obs[s,t];
  }
  array[N_oa] int oa_s; array[N_oa] int oa_t; vector[N_oa] y_oa;
  array[N_op] int op_s; array[N_op] int op_t; vector[N_op] y_op;
  array[N_om] int om_s; array[N_om] int om_t; vector[N_om] y_om;
  {
    int ia = 0; int ip = 0; int im = 0;
    for (s in 1:S) for (t in 1:T) {
      if (y_adult_obs[s,t]) { ia += 1; oa_s[ia]=s; oa_t[ia]=t; y_oa[ia]=y_adult[s,t]; }
      if (y_pup_obs[s,t])   { ip += 1; op_s[ip]=s; op_t[ip]=t; y_op[ip]=y_pup[s,t]; }
      if (y_molt_obs[s,t])  { im += 1; om_s[im]=s; om_t[im]=t; y_om[im]=y_molt[s,t]; }
    }
  }
}

parameters {
  // ── Vital rates ───────────────────────────────────────────────────────────
  real phi_pup_logit;
  real<lower=0,upper=1> phi_juv_base;
  real phi_adult_F_logit;
  real<lower=0,upper=0.10> delta_adult;
  real<lower=0,upper=1> fecund_primip;
  real<lower=0,upper=1> fecund_mature;
  real<lower=0.4,upper=0.6> prop_female;
  real<lower=0,upper=1> rho_pup;
  real<lower=0,upper=1> rho_juv_molt;   // fraction of juveniles attending molt haul-outs

  // ── MOCI: open-coast baseline ─────────────────────────────────────────────
  real beta_moci_ond_fecund;
  real beta_moci_ond_pup;
  real beta_moci_amj_pup;
  real beta_moci_jfm_pup;
  real beta_moci_jfm_juv;
  real beta_moci_jfm_adult;
  real beta_moci_amj_molt;

  // ── Bay/Coast MOCI modifiers ──────────────────────────────────────────────
  real delta_moci_mouth_fecund;
  real delta_moci_mouth_surv;
  real delta_moci_south_fecund;
  real delta_moci_south_surv;

  // ── Marin core-site MOCI fecundity modifier ───────────────────────────────
  // Marin (county 1) is the primary pupping region with the strongest and most
  // consistent MOCI-fecundity link. The Marin IPM v3.3 estimated total
  // beta_moci_ond_fecund = -0.212 (89% CrI: -0.281, -0.141). The regional
  // baseline beta_moci_ond_fecund captures open-coast average sensitivity;
  // delta_moci_marin_fecund captures Marin excess sensitivity above baseline.
  // Prior N(-0.15, 0.10) informed by Marin IPM: expected Marin modifier
  // = -0.212 - baseline(~-0.05) = approximately -0.16.
  real delta_moci_marin_fecund;

  // ── County random effects on vital rates ──────────────────────────────────
  vector[C] county_effect_raw;
  real<lower=0.01,upper=0.5> sigma_county;

  // ── Site-level detection random effects ───────────────────────────────────
  // Single shared site random effect for both breed and molt.
  // Separate breed/molt effects inflated N by driving detect_breed down to 63%
  // (from 90%) while detect_molt remained pathologically low (~13%).
  // The actual fix for molt detection is removing juveniles from N_molt_true
  // (see transformed parameters) rather than separate detection effects.
  vector[S] site_detect_raw;
  real<lower=0.01,upper=0.75> sigma_site;
  real detect_breed_logit;
  real detect_molt_logit;

  // ── Error terms ───────────────────────────────────────────────────────────
  // sigma_obs_pup is NOT a free parameter — see model block for explanation.
  // sigma_process upper bound raised from 0.5 to 0.75.
  // Real harbor seal site counts have genuine ~40-50% year-to-year variability.
  // The 0.5 ceiling caused fecundity/phi_pup compensation when N_adult_M_init
  // was constrained. Prior widened to N(0.30, 0.10) consistent with real data.
  real<lower=0.05,upper=0.75> sigma_process;
  real<lower=0.05,upper=0.6>  sigma_obs_adult;
  real<lower=0.05,upper=0.7>  sigma_obs_molt;

  // ── Initial site populations (non-centred) ────────────────────────────────
  // log_N_adult_M_init_raw removed: N_adult_M_init is set as a fixed fraction
  // of N_adult_F_init (ratio 0.9, consistent with stable age distribution).
  // Separate initialization allowed N_adult_M to inflate unconstrained because
  // p_male_fixed=0.057 makes breed counts nearly insensitive to N_adult_M,
  // while the molt likelihood (which includes all adults) pulled N_adult_M up
  // to explain large molt counts, inflating total N by 3-5x.
  // After t=1 the Leslie matrix updates N_adult_M freely via survival and
  // juvenile recruitment, so the sex ratio can evolve naturally over time.
  vector[S] log_N_adult_F_init_raw;
  vector[S] log_N_juv_init_raw;
  vector[S] log_N_pup_init_raw;
  real      mu_log_adult;
  real      mu_log_juv;
  real      mu_log_pup;
  real<lower=0> sigma_init;

  // ── Site-level process errors ─────────────────────────────────────────────
  // eps_pup_raw removed: pup N is deterministic, observation+process variance
  // combined into sigma_pup_eff in the likelihood.
  // eps_juv_raw restored: removing it forced all process variance through
  // eps_adult_raw alone, driving sigma_process to its ceiling (0.49) and
  // causing fecundity/phi_pup compensation. With N_adult_M_init now fixed
  // at 0.9*N_adult_F_init, detect_molt is positive (~54%), meaning N_adult
  // already explains molt counts — eps_juv_raw has no pressure to inflate.
  matrix[S, T-1] eps_adult_raw;
  matrix[S, T-1] eps_juv_raw;
}

transformed parameters {
  vector[C] county_effect = sigma_county  * county_effect_raw;
  vector[S] site_detect   = sigma_site    * site_detect_raw;

  real phi_pup_base     = inv_logit(phi_pup_logit);
  real phi_adult_F_base = inv_logit(phi_adult_F_logit);
  real phi_adult_M_base = fmax(phi_adult_F_base - delta_adult, 0.01);
  real avg_fecundity    = 0.20 * fecund_primip + 0.80 * fecund_mature;

  // Site initial populations
  vector<lower=0>[S] N_adult_F_init;
  vector<lower=0>[S] N_adult_M_init;
  vector<lower=0>[S] N_juv_F_init;
  vector<lower=0>[S] N_juv_M_init;
  vector<lower=0>[S] N_pup_init;
  for (s in 1:S) {
    N_adult_F_init[s] = exp(mu_log_adult + sigma_init * log_N_adult_F_init_raw[s]);
    // N_adult_M_init linked to N_adult_F_init via stable sex ratio (0.9).
    // This prevents molt likelihood from inflating N_adult_M unconstrained.
    N_adult_M_init[s] = N_adult_F_init[s] * 0.9;
    N_juv_F_init[s]   = exp(mu_log_juv   + sigma_init * log_N_juv_init_raw[s]) * 0.5;
    N_juv_M_init[s]   = exp(mu_log_juv   + sigma_init * log_N_juv_init_raw[s]) * 0.5;
    N_pup_init[s]     = exp(mu_log_pup   + sigma_init * log_N_pup_init_raw[s]);
  }

  // Site-level state arrays
  matrix<lower=0>[S, T] N_adult_F;
  matrix<lower=0>[S, T] N_adult_M;
  matrix<lower=0>[S, T] N_juv_F;
  matrix<lower=0>[S, T] N_juv_M;
  matrix<lower=0>[S, T] N_pup;
  matrix<lower=0>[S, T] N_adult_total;
  matrix<lower=0>[S, T] N_juv_total;
  matrix<lower=0>[S, T] N_molt_true;
  matrix<lower=0>[S, T] N_total;
  matrix<lower=0,upper=1>[S, T] detect_breed_st;
  matrix<lower=0,upper=1>[S, T] detect_molt_st;

  {
    real logit_phi_juv    = logit(phi_juv_base);
    real logit_phi_adultM = logit(phi_adult_M_base);
    real logit_avg_fec    = logit(avg_fecundity);

    for (s in 1:S) {
      int c = county_id[s];
      // Bay/coast modifiers (additive on logit scale)
      real bay_f = (county_type[c] == 1) * delta_moci_mouth_fecund +
                   (county_type[c] == 2) * delta_moci_south_fecund;
      real bay_s = (county_type[c] == 1) * delta_moci_mouth_surv +
                   (county_type[c] == 2) * delta_moci_south_surv;
      // Marin core-site fecundity modifier (county 1 = Marin)
      real marin_f = (c == 1) ? delta_moci_marin_fecund : 0.0;

      // Vital rates and detection for all t
      for (t in 1:T) {
        int t_birth = (t > 1) ? t - 1 : 1;

        real phi_pup_val = inv_logit(
          phi_pup_logit + county_effect[c] +
          (beta_moci_amj_pup + bay_s) * moci_amj[t_birth] +
          (beta_moci_ond_pup + bay_s) * moci_ond[t] +
          (beta_moci_jfm_pup + bay_s) * moci_jfm[t]);

        real phi_juv_val = inv_logit(
          logit_phi_juv + county_effect[c] * 0.5 +
          (beta_moci_jfm_juv + bay_s) * moci_jfm[t]);

        real phi_aF_val = inv_logit(
          phi_adult_F_logit + county_effect[c] * 0.25 +
          (beta_moci_jfm_adult + bay_s) * moci_jfm[t]);

        real phi_aM_val = inv_logit(
          logit_phi_adultM + county_effect[c] * 0.25 +
          (beta_moci_jfm_adult + bay_s) * moci_jfm[t]);

        detect_breed_st[s,t] = inv_logit(detect_breed_logit + site_detect[s]);
        detect_molt_st[s,t]  = inv_logit(detect_molt_logit +
                                          beta_moci_amj_molt * moci_amj[t] +
                                          site_detect[s]);

        // ── Initialise or run Leslie matrix ──────────────────────────────
        // For t <= site_t1[s]: hold N at N_init (no dynamics, no likelihood)
        // For t >  site_t1[s]: Leslie matrix transitions
        if (t <= site_t1[s]) {
          N_adult_F[s,t] = N_adult_F_init[s];
          N_adult_M[s,t] = N_adult_M_init[s];
          N_juv_F[s,t]   = N_juv_F_init[s];
          N_juv_M[s,t]   = N_juv_M_init[s];
          N_pup[s,t]     = N_pup_init[s];
        } else {
          real fecund_t = inv_logit(logit_avg_fec +
                                    (beta_moci_ond_fecund + bay_f + marin_f) * moci_ond[t]);
          real ep  = N_adult_F[s,t-1] * fecund_t;
          real njF = N_pup[s,t-1] * prop_female       * phi_pup_val;
          real njM = N_pup[s,t-1] * (1-prop_female)   * phi_pup_val;
          real jsF = N_juv_F[s,t-1] * phi_juv_val * (2.0/3.0);
          real jsM = N_juv_M[s,t-1] * phi_juv_val * (2.0/3.0);
          real jaF = N_juv_F[s,t-1] * phi_juv_val * (1.0/3.0);
          real jaM = N_juv_M[s,t-1] * phi_juv_val * (1.0/3.0);

          // N_pup is deterministic from the Leslie matrix.
          // Process + observation variance for pups is handled jointly
          // in the likelihood via sigma_pup_eff (see model block).
          N_pup[s,t]     = fmax(ep, 1.0);
          N_juv_F[s,t]   = exp(log(fmax(njF+jsF, 0.1)) + sigma_process*eps_juv_raw[s,t-1]*0.5);
          N_juv_M[s,t]   = exp(log(fmax(njM+jsM, 0.1)) + sigma_process*eps_juv_raw[s,t-1]*0.5);
          N_adult_F[s,t] = exp(log(fmax(N_adult_F[s,t-1]*phi_aF_val+jaF, 1))
                               + sigma_process*eps_adult_raw[s,t-1]*0.5);
          N_adult_M[s,t] = exp(log(fmax(N_adult_M[s,t-1]*phi_aM_val+jaM, 1))
                               + sigma_process*eps_adult_raw[s,t-1]*0.5);
        }

        N_adult_total[s,t] = N_adult_F[s,t] + N_adult_M[s,t];
        N_juv_total[s,t]   = N_juv_F[s,t]   + N_juv_M[s,t];
        N_molt_true[s,t]   = N_adult_total[s,t]
                             + rho_juv_molt * N_juv_total[s,t]
                             + rho_pup * N_pup[s,t];
        N_total[s,t]       = N_pup[s,t] + N_juv_total[s,t] + N_adult_total[s,t];
      }
    }
  }
}

model {
  // ── Priors: vital rates (identical to Marin IPM v3.3) ────────────────────
  phi_pup_logit     ~ normal(-1.2, 0.5);
  phi_juv_base      ~ beta(14, 6);
  phi_adult_F_logit ~ normal(2.50, 0.25);   // prior updated from 2.20: Marin IPM v3.3 posterior mean
  delta_adult       ~ normal(0.05, 0.025);
  fecund_primip     ~ beta(12, 8);
  fecund_mature     ~ beta(17, 3);
  prop_female       ~ beta(50, 50);
  rho_pup           ~ beta(2, 5);
  // rho_juv_molt: juvenile molt attendance fraction.
  // Beta(8,14): mean=0.36, 89% CrI approx [0.21, 0.53].
  // Reflects that roughly a third of molt haul-out animals are juveniles,
  // with adults dominating but juveniles a meaningful presence.
  // Tighter than rho_pup prior to prevent confound with detect_molt_logit.
  rho_juv_molt      ~ beta(8, 14);

  // ── Priors: MOCI (identical to Marin IPM v3.3) ───────────────────────────
  beta_moci_ond_fecund ~ normal(-0.15, 0.20);
  beta_moci_ond_pup    ~ normal(-0.15, 0.20);
  beta_moci_amj_pup    ~ normal(-0.15, 0.20);
  beta_moci_jfm_pup    ~ normal(-0.15, 0.20);
  beta_moci_jfm_juv    ~ normal(-0.15, 0.15);
  beta_moci_jfm_adult  ~ normal(-0.10, 0.12);
  beta_moci_amj_molt   ~ normal(0.05, 0.15);

  // ── Priors: Bay modifiers ─────────────────────────────────────────────────
  delta_moci_mouth_fecund ~ normal(0, 0.15);
  delta_moci_mouth_surv   ~ normal(0, 0.15);
  delta_moci_south_fecund ~ normal(0, 0.20);
  delta_moci_south_surv   ~ normal(0, 0.20);
  // Marin core-site fecundity modifier: informed by Marin IPM v3.3.
  // Marin IPM total effect = -0.212; regional baseline ~-0.05;
  // implied Marin-specific component ~-0.16. Prior N(-0.15, 0.10).
  delta_moci_marin_fecund ~ normal(-0.15, 0.10);

  // ── Priors: county and site random effects ────────────────────────────────
  sigma_county       ~ normal(0.15, 0.10);
  county_effect_raw  ~ std_normal();
  sigma_site         ~ normal(0.20, 0.10);
  site_detect_raw    ~ std_normal();
  detect_breed_logit ~ normal(1.20, 0.50);
  detect_molt_logit  ~ normal(0.75, 0.50);

  // ── Priors: error terms ───────────────────────────────────────────────────
  sigma_process   ~ normal(0.30, 0.10);   // widened from N(0.20,0.08): real data shows ~0.49
  sigma_obs_adult ~ normal(0.18, 0.06);
  sigma_obs_molt  ~ normal(0.35, 0.10);
  // sigma_obs_pup is not a free parameter. Instead, pup process variance
  // (sigma_process^2) and observation variance (sigma_obs_pup_fixed^2) are
  // combined into a single effective sigma on the log scale:
  //   sigma_pup_eff = sqrt(sigma_process^2 + sigma_obs_pup_fixed^2)
  // sigma_obs_pup_fixed = 0.148 is the Marin IPM v3.3 posterior mean,
  // estimated from 28 years of intensive pup counts at 6 PRNS sites.
  // This marginalisation removes eps_pup_raw (480 parameters) and eliminates
  // the sigma_process/sigma_obs_pup geometric confound that caused E-BFMI < 0.3.

  // ── Priors: initial populations (site-level) ──────────────────────────────
  mu_log_adult ~ normal(5.5, 0.5);
  mu_log_juv   ~ normal(4.5, 0.5);
  mu_log_pup   ~ normal(4.5, 0.5);
  sigma_init   ~ exponential(3);
  log_N_adult_F_init_raw ~ std_normal();
  log_N_juv_init_raw     ~ std_normal();
  log_N_pup_init_raw     ~ std_normal();

  to_vector(eps_adult_raw) ~ std_normal();
  to_vector(eps_juv_raw)   ~ std_normal();

  // ── Likelihood ────────────────────────────────────────────────────────────
  {
    // Effective pup sigma: combines process error (sigma_process) and
    // observation error (sigma_obs_pup_fixed from Marin IPM) in quadrature.
    // On the log scale these add as variances: var_total = var_proc + var_obs.
    real sigma_obs_pup_fixed = 0.148;
    real sigma_pup_eff = sqrt(sigma_process^2 + sigma_obs_pup_fixed^2);

    vector[N_oa] mu_a;
    vector[N_op] mu_p;
    vector[N_om] mu_m;

    for (i in 1:N_oa) {
      int s = oa_s[i]; int t = oa_t[i];
      mu_a[i] = log((N_adult_F[s,t] + N_adult_M[s,t]*p_male_fixed) * detect_breed_st[s,t]);
    }
    for (i in 1:N_op) {
      int s = op_s[i]; int t = op_t[i];
      // Expected log pup count = log(N_pup_deterministic * detect)
      // Variance = sigma_pup_eff^2 (process + observation combined)
      mu_p[i] = log(N_pup[s,t] * detect_breed_st[s,t]);
    }
    for (i in 1:N_om) {
      int s = om_s[i]; int t = om_t[i];
      mu_m[i] = log(N_molt_true[s,t] * detect_molt_st[s,t]);
    }

    y_oa ~ normal(mu_a, sigma_obs_adult);
    y_op ~ normal(mu_p, sigma_pup_eff);
    y_om ~ normal(mu_m, sigma_obs_molt);
  }
}

generated quantities {
  // County totals: sum site N within each county
  matrix[C, T] N_total_county;
  vector[T]    N_total_regional;
  matrix[C, T-1] lambda_county;

  // PPC at site level
  matrix[S, T] y_adult_rep;
  matrix[S, T] y_pup_rep;
  matrix[S, T] y_molt_rep;

  {
    for (c in 1:C) for (t in 1:T) N_total_county[c,t] = 0;
    for (s in 1:S) for (t in 1:T)
      N_total_county[county_id[s],t] += N_total[s,t];
  }
  for (t in 1:T) N_total_regional[t] = sum(col(N_total_county, t));
  for (c in 1:C) for (t in 1:(T-1))
    lambda_county[c,t] = N_total_county[c,t+1] / N_total_county[c,t];

  for (s in 1:S) for (t in 1:T) {
    real N_ao = N_adult_F[s,t] + N_adult_M[s,t]*p_male_fixed;
    real sigma_obs_pup_fixed = 0.148;
    real sigma_pup_eff = sqrt(sigma_process^2 + sigma_obs_pup_fixed^2);
    y_adult_rep[s,t] = normal_rng(log(N_ao * detect_breed_st[s,t]), sigma_obs_adult);
    y_pup_rep[s,t]   = normal_rng(log(N_pup[s,t] * detect_breed_st[s,t]), sigma_pup_eff);
    y_molt_rep[s,t]  = normal_rng(log(N_molt_true[s,t] * detect_molt_st[s,t]), sigma_obs_molt);
  }

  // ── Projections ──────────────────────────────────────────────────────────
  // Forward projections by scenario: status quo (MOCI=0), warm (MOCI=+1),
  // cool (MOCI=-1). Starting from fitted N at T (last historical year).
  // N_total_county_proj[c, scenario, tp]: county totals by year and scenario.
  array[C, N_scenarios, T_proj] real N_total_county_proj;
  {
    real logit_avg_fec    = logit(avg_fecundity);
    real logit_phi_juv    = logit(phi_juv_base);
    real logit_phi_adultM = logit(phi_adult_M_base);

    for (sc in 1:N_scenarios) {
      // Site-level state at start of projection (= last historical year T)
      array[S] real pN_adult_F;
      array[S] real pN_adult_M;
      array[S] real pN_juv_F;
      array[S] real pN_juv_M;
      array[S] real pN_pup;
      for (s in 1:S) {
        pN_adult_F[s] = N_adult_F[s,T];
        pN_adult_M[s] = N_adult_M[s,T];
        pN_juv_F[s]   = N_juv_F[s,T];
        pN_juv_M[s]   = N_juv_M[s,T];
        pN_pup[s]     = N_pup[s,T];
      }
      // Initialise county totals
      for (c in 1:C) for (tp in 1:T_proj) N_total_county_proj[c,sc,tp] = 0;

      for (tp in 1:T_proj) {
        array[S] real nN_adult_F;
        array[S] real nN_adult_M;
        array[S] real nN_juv_F;
        array[S] real nN_juv_M;
        array[S] real nN_pup;

        for (s in 1:S) {
          int c   = county_id[s];
          real mf = moci_proj[sc, tp];   // MOCI value for this scenario and year
          real bay_f = (county_type[c] == 1) * delta_moci_mouth_fecund +
                       (county_type[c] == 2) * delta_moci_south_fecund;
          real bay_s = (county_type[c] == 1) * delta_moci_mouth_surv +
                       (county_type[c] == 2) * delta_moci_south_surv;
          real mar_f = (c == 1) ? delta_moci_marin_fecund : 0.0;

          real phi_pup_proj   = inv_logit(phi_pup_logit     + county_effect[c] +
                                          (beta_moci_amj_pup + bay_s) * mf +
                                          (beta_moci_ond_pup + bay_s) * mf +
                                          (beta_moci_jfm_pup + bay_s) * mf);
          real phi_juv_proj   = inv_logit(logit_phi_juv     + county_effect[c]*0.5 +
                                          (beta_moci_jfm_juv + bay_s) * mf);
          real phi_aF_proj    = inv_logit(phi_adult_F_logit + county_effect[c]*0.25 +
                                          (beta_moci_jfm_adult + bay_s) * mf);
          real phi_aM_proj    = inv_logit(logit_phi_adultM  + county_effect[c]*0.25 +
                                          (beta_moci_jfm_adult + bay_s) * mf);
          real fecund_proj    = inv_logit(logit_avg_fec +
                                          (beta_moci_ond_fecund + bay_f + mar_f) * mf);

          real ep  = pN_adult_F[s] * fecund_proj;
          real njF = pN_pup[s]   * prop_female     * phi_pup_proj;
          real njM = pN_pup[s]   * (1-prop_female) * phi_pup_proj;
          real jsF = pN_juv_F[s] * phi_juv_proj * (2.0/3.0);
          real jsM = pN_juv_M[s] * phi_juv_proj * (2.0/3.0);
          real jaF = pN_juv_F[s] * phi_juv_proj * (1.0/3.0);
          real jaM = pN_juv_M[s] * phi_juv_proj * (1.0/3.0);

          // Deterministic Leslie matrix projection (no process error)
          nN_pup[s]     = fmax(ep, 0.01);
          nN_juv_F[s]   = fmax(njF + jsF, 0.01);
          nN_juv_M[s]   = fmax(njM + jsM, 0.01);
          nN_adult_F[s] = fmax(pN_adult_F[s]*phi_aF_proj + jaF, 0.01);
          nN_adult_M[s] = fmax(pN_adult_M[s]*phi_aM_proj + jaM, 0.01);

          real N_site_proj = nN_pup[s] + nN_juv_F[s] + nN_juv_M[s] +
                             nN_adult_F[s] + nN_adult_M[s];
          N_total_county_proj[c,sc,tp] += N_site_proj;
        }
        // Advance state
        for (s in 1:S) {
          pN_adult_F[s] = nN_adult_F[s];
          pN_adult_M[s] = nN_adult_M[s];
          pN_juv_F[s]   = nN_juv_F[s];
          pN_juv_M[s]   = nN_juv_M[s];
          pN_pup[s]     = nN_pup[s];
        }
      }
    }
  }
}
'

cat("Stan model written to Code/harbor_seal_regional_ipm.stan\n")

# ── PREPARE REAL DATA FOR STAN ────────────────────────────────────────────────
prepare_real_data_regional <- function(input_rds = "Output/regional_ipm_input_data.rds",
                                       T_proj = 5) {   # 5-year projections
  inp <- readRDS(input_rds)
  sd  <- inp$stan_data
  
  # T_proj, N_scenarios and moci_proj are already saved in the RDS by
  # 15_regional_data_prep.R — do NOT re-add them or Stan will see duplicates.
  # Only override if they differ from what was saved (e.g. user changes T_proj).
  if (!is.null(sd$T_proj) && sd$T_proj == T_proj) return(sd)
  
  # T_proj changed: rebuild projection matrix and update
  sd$T_proj      <- T_proj
  sd$N_scenarios <- 3L
  sd$moci_proj   <- matrix(c(0, 1, -1), nrow = 3L, ncol = T_proj)
  sd
}

# ── ORCHESTRATOR ──────────────────────────────────────────────────────────────
run_regional_ipm <- function(use_real_data  = FALSE,
                             sim_data       = NULL,
                             stan_data_list = NULL,   # override Stan data (sensitivity runs)
                             input_rds      = "Output/regional_ipm_input_data.rds",
                             prefix         = NULL,
                             T_proj         = 10,
                             seed           = 456,
                             iter_warmup    = 2000,
                             iter_sampling  = 2000,
                             adapt_delta    = 0.97,
                             max_treedepth  = 12,
                             save_fit       = TRUE) {
  
  cat("\n================================================================\n")
  cat("   REGIONAL HARBOR SEAL IPM (Bay/Coast MOCI + County Dynamics)\n")
  cat("================================================================\n\n")
  
  # Prefix controls ALL output file names (fit RDS, input RDS, CSVs, plots).
  # Set explicitly from 19_regional_run.R to ensure consistency across scripts.
  if (is.null(prefix))
    prefix <- ifelse(use_real_data, "Regional_real", "Regional_sim")
  cat(sprintf("Output prefix: %s\n", prefix))
  
  if (use_real_data) {
    stan_dat <- if (!is.null(stan_data_list)) stan_data_list
    else prepare_real_data_regional(input_rds, T_proj)
    site_meta <- readRDS(input_rds)$site_meta
    years     <- readRDS(input_rds)$years
    scenario_names <- c("Status Quo","Warm (MOCI +1)","Cool (MOCI -1)")
    model_data <- list(stan_data      = stan_dat,
                       site_meta      = site_meta,
                       years          = years,
                       county_names   = c("Marin","Bay Mouth","South Bay",
                                          "San Mateo","Sonoma","Mendocino"),
                       scenario_names = scenario_names,
                       true_params    = NULL)
  } else {
    if (is.null(sim_data)) {
      if (!exists("sim_regional")) source("Code/16_regional_simulate.R")
      sim_data <- sim_regional
    }
    # Add projection fields if not present
    if (is.null(sim_data$stan_data$T_proj)) {
      sim_data$stan_data$T_proj      <- T_proj
      sim_data$stan_data$N_scenarios <- 3
      sim_data$stan_data$moci_proj   <- matrix(c(0,1,-1), 3, T_proj)
    }
    model_data <- sim_data
  }
  
  cat("Compiling Stan model...\n")
  write_lines(stan_code_regional, "Code/harbor_seal_regional_ipm.stan")
  # Delete compiled binary before compiling. CmdStan 2.38.0 uses incremental
  # compilation — stanc transpiles the .stan to C++ (~20 sec) and relinks
  # against cached object files. Deleting the binary ensures the final
  # executable is always built from the current .stan file.
  binary_path <- tools::file_path_sans_ext("Code/harbor_seal_regional_ipm.stan")
  if (file.exists(binary_path)) file.remove(binary_path)
  model <- cmdstan_model("Code/harbor_seal_regional_ipm.stan")
  
  cat(sprintf("Running MCMC (warmup=%d, sampling=%d, adapt_delta=%.3f)...\n",
              iter_warmup, iter_sampling, adapt_delta))
  fit <- model$sample(
    data            = model_data$stan_data,
    seed            = seed,
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = iter_warmup,
    iter_sampling   = iter_sampling,
    refresh         = 200,
    adapt_delta     = adapt_delta,
    max_treedepth   = max_treedepth
  )
  
  if (save_fit) {
    fit$save_object(paste0("Output/harbor_seal_", prefix, "_fit.rds"))
    saveRDS(list(stan_data = model_data$stan_data, years = model_data$years,
                 true_params = model_data$true_params %||% NULL),
            paste0("Output/harbor_seal_", prefix, "_input_data.rds"))
    cat(sprintf("Fit saved: Output/harbor_seal_%s_fit.rds\n", prefix))
  } else {
    cat(sprintf("Fit not saved to disk (save_fit=FALSE): %s\n", prefix))
  }
  
  # ── Quick diagnostics ──────────────────────────────────────────────────────
  cat("\n=== MODEL DIAGNOSTICS ===\n")
  tryCatch({
    ds <- fit$diagnostic_summary(quiet = TRUE)
    cat(if (all(ds$num_max_treedepth == 0)) "Treedepth OK.\n"
        else sprintf("WARN: %d treedepth hits.\n", sum(ds$num_max_treedepth)))
    cat(if (all(ds$num_divergent == 0)) "No divergences.\n"
        else sprintf("WARN: %d divergences.\n", sum(ds$num_divergent)))
    cat(if (all(ds$ebfmi > 0.2)) "E-BFMI OK.\n"
        else sprintf("WARN: low E-BFMI in %d chains.\n", sum(ds$ebfmi <= 0.2)))
  }, error = function(e) cat("Diagnostics unavailable (CSVs gone).\n"))
  
  # ── Key results ────────────────────────────────────────────────────────────
  key_params <- c("phi_pup_logit","phi_juv_base","phi_adult_F_logit","delta_adult",
                  "avg_fecundity","rho_pup","rho_juv_molt",
                  "beta_moci_ond_fecund","beta_moci_jfm_adult",
                  "delta_moci_mouth_fecund","delta_moci_mouth_surv",
                  "delta_moci_south_fecund","delta_moci_south_surv",
                  "delta_moci_marin_fecund",
                  "sigma_county","sigma_site",
                  "detect_breed_logit","detect_molt_logit",
                  "sigma_process","sigma_obs_adult","sigma_obs_molt")
  s <- tryCatch(fit$summary(variables = key_params), error = function(e) NULL)
  if (!is.null(s)) {
    cat("\n=== KEY PARAMETER ESTIMATES ===\n")
    print(s |> select(variable, mean, sd, q5, q95, rhat, ess_bulk), n = 30)
  }
  
  # Lambda by county
  cat("\n=== COUNTY LAMBDA (mean annual 2005-2024) ===\n")
  county_names <- model_data$county_names %||%
    c("Marin","Bay Mouth","South Bay","San Mateo","Sonoma","Mendocino")
  for (c in seq_along(county_names)) {
    lam_vars <- paste0("lambda_county[",c,",",1:(model_data$stan_data$T-1),"]")
    lam_draws <- tryCatch(fit$draws(variables = lam_vars, format = "matrix"),
                          error = function(e) NULL)
    if (!is.null(lam_draws)) {
      lam_mean <- mean(colMeans(lam_draws))
      cat(sprintf("  %s: lambda_bar = %.3f\n", county_names[c], lam_mean))
    }
  }
  
  # MOCI modifiers by county type — the key new results
  mod_params <- c("delta_moci_mouth_fecund","delta_moci_mouth_surv",
                  "delta_moci_south_fecund","delta_moci_south_surv")
  mod_summ <- tryCatch(fit$summary(variables = mod_params), error = function(e) NULL)
  if (!is.null(mod_summ)) {
    cat("\n=== MOCI MODIFIERS BY COUNTY TYPE ===\n")
    print(mod_summ |> select(variable, mean, sd, q5, q95))
    for (p in mod_params) {
      d <- as.numeric(fit$draws(p, format = "matrix"))
      cat(sprintf("  %-28s P(>0) = %.3f\n", p, mean(d > 0)))
    }
    cat("  Expected gradient: delta_south > delta_mouth > 0\n")
  }
  
  cat("\n================================================================\n")
  cat(sprintf("   COMPLETE — Regional IPM\n"))
  cat(sprintf("   Fit -> Output/harbor_seal_%s_fit.rds\n", prefix))
  cat("================================================================\n\n")
  
  list(fit         = fit,
       model       = model,
       model_data  = model_data,
       prefix      = prefix)
}

cat("\n17_regional_ipm_model.R loaded — SITE-LEVEL MODEL.\n")
cat("Key features:\n")
cat("  - Site-level N state variables (S=24): eliminates N/alpha confound\n")
cat("  - County-level vital rate pooling via county_effect[county_id[s]]\n")
cat("  - Site detection random effects (sigma_site)\n")
cat("  - MOCI modifiers by county type (mouth/south fecundity + survival)\n")
cat("  - site_t1[s]: Leslie matrix held flat before site's first survey year\n")
cat("  - County totals computed as sums in generated quantities\n")
cat("  - p_male_breed FIXED at 0.057 (Marin IPM posterior)\n")
