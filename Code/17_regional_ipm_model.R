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
// REGIONAL HARBOR SEAL IPM
// County-level dynamics + site availability + Three-way MOCI modifier
// (Open Coast / Bay Mouth / South Bay)
// ============================================================================
data {
  int<lower=1> T;          // years (21 for 2005-2025; index 16 = 2020 latent)
  int<lower=1> S;          // sites (24)
  int<lower=1> C;          // counties (6)
  int<lower=1> T_proj;
  int<lower=1> N_scenarios;

  // ── Observations ──────────────────────────────────────────────────────────
  // Log-scale counts; 0.0 where unobserved (y_*_obs = 0).
  // True zero counts (surveyed, no animals seen) are encoded as log(0.5)
  // in 15_regional_data_prep.R and enter the likelihood normally as real
  // observations — no special handling needed here.
  // NA (not surveyed) → y_*_obs = 0 → excluded from likelihood.
  matrix[S, T] y_adult;
  matrix[S, T] y_pup;
  matrix[S, T] y_molt;
  array[S, T] int<lower=0,upper=1> y_adult_obs;   // 1 = surveyed
  array[S, T] int<lower=0,upper=1> y_pup_obs;
  array[S, T] int<lower=0,upper=1> y_molt_obs;

  // ── Site-level classifiers ───────────────────────────────────────────────
  array[S] int<lower=1,upper=6> county_id;    // which county each site belongs to
  array[S] int<lower=0,upper=1> is_breeder;   // 1 = established pup production here

  // ── County-level oceanographic type ──────────────────────────────────────
  // 0 = open coast  (Marin, San Mateo, Sonoma, Mendocino)
  // 1 = bay mouth   (Castro Rocks, Yerba Buena, Alcatraz)
  // 2 = south bay   (Mowry Slough, Newark Slough)
  array[C] int<lower=0,upper=2> county_type;

  // ── County-level start index ──────────────────────────────────────────────
  // county_t1[c] = first time index with survey data for county c.
  // For counties with sites starting after 2005 (currently Mendocino, t1=15),
  // N is held at N_init for t < county_t1[c] — no Leslie transitions — then
  // the matrix runs normally from county_t1[c]+1 onward.  This prevents the
  // population from drifting to near-zero over unobserved years before data
  // appear.  For all other counties county_t1[c]=1 (behaviour unchanged).
  array[C] int<lower=1,upper=20> county_t1;

  // ── Number of sites per county ────────────────────────────────────────────
  array[C] int<lower=1> n_sites_county;

  // ── MOCI covariates (z-scored; length T each) ────────────────────────────
  vector[T] moci_jfm;
  vector[T] moci_amj;
  vector[T] moci_ond;

  // ── Projections ───────────────────────────────────────────────────────────
  matrix[N_scenarios, T_proj] moci_proj;

  // ── Fixed male detection fraction (from Marin IPM posterior) ─────────────
  real<lower=0,upper=1> p_male_fixed;
}

transformed data {
  // Flatten observed site-year-class triplets for vectorized likelihood
  int N_oa = 0; int N_op = 0; int N_om = 0;
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
      if (y_adult_obs[s,t]) { ia += 1; oa_s[ia]=s; oa_t[ia]=t; y_oa[ia]=y_adult[s,t]; }
      if (y_pup_obs[s,t])   { ip += 1; op_s[ip]=s; op_t[ip]=t; y_op[ip]=y_pup[s,t]; }
      if (y_molt_obs[s,t])  { im += 1; om_s[im]=s; om_t[im]=t; y_om[im]=y_molt[s,t]; }
    }
  }
}

parameters {
  // ── Vital rates (shared across all counties; county effect adds variation) ─
  real phi_pup_logit;
  real<lower=0,upper=1> phi_juv_base;
  real phi_adult_F_logit;
  real<lower=0,upper=0.10> delta_adult;
  real<lower=0,upper=1> fecund_primip;
  real<lower=0,upper=1> fecund_mature;
  real<lower=0.4,upper=0.6> prop_female;

  // ── Pup molt attendance (v3.3, unchanged from Marin IPM) ─────────────────
  real<lower=0,upper=1> rho_pup;

  // ── MOCI effects: open-coast baseline ────────────────────────────────────
  real beta_moci_ond_fecund;
  real beta_moci_ond_pup;
  real beta_moci_amj_pup;
  real beta_moci_jfm_pup;
  real beta_moci_jfm_juv;
  real beta_moci_jfm_adult;
  real beta_moci_amj_molt;

  // ── MOCI modifiers by county type (KEY NEW PARAMETERS) ──────────────────
  // Additive on logit scale relative to open-coast baseline.
  // Positive = that county type less negatively affected by poor upwelling.
  // Biological expectation: delta_south > delta_mouth > 0 (gradient of
  // MOCI decoupling from open coast → bay mouth → deep estuary).
  // Separate fecundity (OND implantation window) vs survival pathways.
  real delta_moci_mouth_fecund;   // Bay Mouth vs coast, fecundity
  real delta_moci_mouth_surv;     // Bay Mouth vs coast, survival
  real delta_moci_south_fecund;   // South Bay vs coast, fecundity
  real delta_moci_south_surv;     // South Bay vs coast, survival

  // ── County random effects (partial pooling across 6 counties) ────────────
  vector[C] county_effect_raw;      // non-centred parameterisation
  real<lower=0.01,upper=0.5> sigma_county;

  // ── Site availability parameters (partial pooling within county) ──────────
  // log(alpha) = log(site_availability × detection); one per site per survey type
  vector[S] log_alpha_breed_raw;
  vector[S] log_alpha_pup_raw;
  vector[S] log_alpha_molt_raw;
  // Hyperparameters: county-level mean availability
  vector[C] mu_log_alpha_breed;
  vector[C] mu_log_alpha_pup;
  vector[C] mu_log_alpha_molt;
  real<lower=0.05,upper=1.5> sigma_log_alpha;   // upper raised then tightened; real data showed hitting 1.0

  // ── Molt detection MOCI response ─────────────────────────────────────────
  // Separate from availability: how MOCI (AMJ) shifts detection during molt survey
  // (same structure as Marin IPM)
  real detect_molt_baseline;        // logit-scale baseline

  // ── Error terms ───────────────────────────────────────────────────────────
  real<lower=0.05,upper=0.5>  sigma_process;
  real<lower=0.05,upper=0.6>  sigma_obs_adult;
  real<lower=0.02,upper=0.5>  sigma_obs_pup;
  real<lower=0.05,upper=0.7>  sigma_obs_molt;

  // ── Initial county populations (non-centred) ──────────────────────────────
  vector[C] log_N_adult_F_init_raw;
  vector[C] log_N_adult_M_init_raw;
  vector[C] log_N_juv_init_raw;
  vector[C] log_N_pup_init_raw;
  real      mu_log_adult;
  real      mu_log_juv;
  real      mu_log_pup;
  real<lower=0> sigma_init;

  // ── County-level process errors (non-centred) ─────────────────────────────
  matrix[C, T-1] eps_adult_raw;
  matrix[C, T-1] eps_juv_raw;
  matrix[C, T-1] eps_pup_raw;
}

transformed parameters {
  // ── County effects and site availability ──────────────────────────────────
  vector[C] county_effect = sigma_county * county_effect_raw;

  // Site availability: county mean + site deviation
  vector[S] log_alpha_breed;
  vector[S] log_alpha_pup;
  vector[S] log_alpha_molt;
  for (s in 1:S) {
    log_alpha_breed[s] = mu_log_alpha_breed[county_id[s]] +
                         sigma_log_alpha * log_alpha_breed_raw[s];
    log_alpha_pup[s]   = mu_log_alpha_pup[county_id[s]]   +
                         sigma_log_alpha * log_alpha_pup_raw[s];
    log_alpha_molt[s]  = mu_log_alpha_molt[county_id[s]]  +
                         sigma_log_alpha * log_alpha_molt_raw[s];
  }

  // ── Derived vital rate scalars ─────────────────────────────────────────────
  real phi_pup_base     = inv_logit(phi_pup_logit);
  real phi_adult_F_base = inv_logit(phi_adult_F_logit);
  real phi_adult_M_base = fmax(phi_adult_F_base - delta_adult, 0.01);
  real avg_fecundity    = 0.20 * fecund_primip + 0.80 * fecund_mature;

  // ── County initial populations ────────────────────────────────────────────
  vector<lower=0>[C] N_adult_F_init;
  vector<lower=0>[C] N_adult_M_init;
  vector<lower=0>[C] N_juv_F_init;
  vector<lower=0>[C] N_juv_M_init;
  vector<lower=0>[C] N_pup_init;
  for (c in 1:C) {
    N_adult_F_init[c] = exp(mu_log_adult + sigma_init * log_N_adult_F_init_raw[c]);
    N_adult_M_init[c] = exp(mu_log_adult + sigma_init * log_N_adult_M_init_raw[c]) * 0.9;
    N_juv_F_init[c]   = exp(mu_log_juv   + sigma_init * log_N_juv_init_raw[c]) * 0.5;
    N_juv_M_init[c]   = exp(mu_log_juv   + sigma_init * log_N_juv_init_raw[c]) * 0.5;
    N_pup_init[c]     = exp(mu_log_pup   + sigma_init * log_N_pup_init_raw[c]);
  }

  // ── County-level population dynamics (Leslie matrix) ──────────────────────
  matrix<lower=0>[C, T] N_adult_F;
  matrix<lower=0>[C, T] N_adult_M;
  matrix<lower=0>[C, T] N_juv_F;
  matrix<lower=0>[C, T] N_juv_M;
  matrix<lower=0>[C, T] N_pup;
  matrix<lower=0>[C, T] N_adult_total;
  matrix<lower=0>[C, T] N_juv_total;
  matrix<lower=0>[C, T] N_molt_true;   // includes rho_pup * N_pup
  matrix<lower=0>[C, T] N_total;

  // Time-varying vital rates for each county
  matrix<lower=0,upper=1>[C, T] phi_pup_ct;
  matrix<lower=0,upper=1>[C, T] phi_juv_ct;
  matrix<lower=0,upper=1>[C, T] phi_adult_F_ct;
  matrix<lower=0,upper=1>[C, T] phi_adult_M_ct;
  matrix<lower=0,upper=1>[C, T] detect_molt_ct;

  {
    real logit_phi_juv    = logit(phi_juv_base);
    real logit_phi_adultM = logit(phi_adult_M_base);
    real logit_avg_fec    = logit(avg_fecundity);

    for (c in 1:C) {
      // MOCI modifier by county type:
      //   type 0 (coast):     no modifier
      //   type 1 (bay mouth): delta_moci_mouth_*
      //   type 2 (south bay): delta_moci_south_*
      real bay_f = (county_type[c] == 1) * delta_moci_mouth_fecund +
                   (county_type[c] == 2) * delta_moci_south_fecund;
      real bay_s = (county_type[c] == 1) * delta_moci_mouth_surv +
                   (county_type[c] == 2) * delta_moci_south_surv;

      for (t in 1:T) {
        int t_birth = (t > 1) ? t - 1 : 1;

        phi_pup_ct[c, t] = inv_logit(
          phi_pup_logit + county_effect[c] +
          (beta_moci_amj_pup + bay_s) * moci_amj[t_birth] +
          (beta_moci_ond_pup + bay_s) * moci_ond[t] +
          (beta_moci_jfm_pup + bay_s) * moci_jfm[t]);

        phi_juv_ct[c, t] = inv_logit(
          logit_phi_juv + county_effect[c] * 0.5 +
          (beta_moci_jfm_juv + bay_s) * moci_jfm[t]);

        phi_adult_F_ct[c, t] = inv_logit(
          phi_adult_F_logit + county_effect[c] * 0.25 +
          (beta_moci_jfm_adult + bay_s) * moci_jfm[t]);

        phi_adult_M_ct[c, t] = inv_logit(
          logit_phi_adultM + county_effect[c] * 0.25 +
          (beta_moci_jfm_adult + bay_s) * moci_jfm[t]);

        detect_molt_ct[c, t] = inv_logit(
          detect_molt_baseline + beta_moci_amj_molt * moci_amj[t]);
      }

      // ── Initialise and run Leslie matrix ─────────────────────────────────
      // For t <= county_t1[c]: hold N at N_init (no dynamics).
      //   - No likelihood fires (y_obs=0 for all sites in this county pre-start).
      //   - Prevents unobserved drift to zero for late-starting counties.
      // For t > county_t1[c]: run Leslie matrix normally.
      // For all counties except Mendocino, county_t1[c]=1 — behaviour unchanged.
      for (t in 1:county_t1[c]) {
        N_adult_F[c, t] = N_adult_F_init[c];
        N_adult_M[c, t] = N_adult_M_init[c];
        N_juv_F[c, t]   = N_juv_F_init[c];
        N_juv_M[c, t]   = N_juv_M_init[c];
        N_pup[c, t]     = N_pup_init[c];
        N_adult_total[c, t] = N_adult_F_init[c] + N_adult_M_init[c];
        N_juv_total[c, t]   = N_juv_F_init[c]   + N_juv_M_init[c];
        N_molt_true[c, t]   = N_juv_total[c,t] + N_adult_total[c,t]
                              + rho_pup * N_pup_init[c];
        N_total[c, t]       = N_pup_init[c] + N_juv_total[c,t] + N_adult_total[c,t];
      }

      for (t in (county_t1[c]+1):T) {
        real fecund_t = inv_logit(logit_avg_fec +
                                  (beta_moci_ond_fecund + bay_f) * moci_ond[t]);
        real ep  = N_adult_F[c, t-1] * fecund_t;
        real njF = N_pup[c,t-1] * prop_female       * phi_pup_ct[c,t];
        real njM = N_pup[c,t-1] * (1-prop_female)   * phi_pup_ct[c,t];
        real jsF = N_juv_F[c,t-1] * phi_juv_ct[c,t] * (2.0/3.0);
        real jsM = N_juv_M[c,t-1] * phi_juv_ct[c,t] * (2.0/3.0);
        real jaF = N_juv_F[c,t-1] * phi_juv_ct[c,t] * (1.0/3.0);
        real jaM = N_juv_M[c,t-1] * phi_juv_ct[c,t] * (1.0/3.0);

        N_pup[c,t]     = exp(log(fmax(ep,           1)) + sigma_process*eps_pup_raw[c,t-1]);
        N_juv_F[c,t]   = exp(log(fmax(njF+jsF,   0.1)) + sigma_process*eps_juv_raw[c,t-1]*0.5);
        N_juv_M[c,t]   = exp(log(fmax(njM+jsM,   0.1)) + sigma_process*eps_juv_raw[c,t-1]*0.5);
        N_adult_F[c,t] = exp(log(fmax(N_adult_F[c,t-1]*phi_adult_F_ct[c,t]+jaF, 1))
                             + sigma_process*eps_adult_raw[c,t-1]*0.5);
        N_adult_M[c,t] = exp(log(fmax(N_adult_M[c,t-1]*phi_adult_M_ct[c,t]+jaM, 1))
                             + sigma_process*eps_adult_raw[c,t-1]*0.5);

        N_adult_total[c,t] = N_adult_F[c,t] + N_adult_M[c,t];
        N_juv_total[c,t]   = N_juv_F[c,t]   + N_juv_M[c,t];
        N_molt_true[c,t]   = N_juv_total[c,t] + N_adult_total[c,t] + rho_pup*N_pup[c,t];
        N_total[c,t]       = N_pup[c,t] + N_juv_total[c,t] + N_adult_total[c,t];
      }
    }
  }
}

model {
  // ── Priors: vital rates (identical to Marin IPM v3.3) ───────────────────
  phi_pup_logit     ~ normal(-1.2, 0.5);
  phi_juv_base      ~ beta(14, 6);
  phi_adult_F_logit ~ normal(2.20, 0.25);
  delta_adult       ~ normal(0.05, 0.025);
  fecund_primip     ~ beta(12, 8);
  fecund_mature     ~ beta(17, 3);
  prop_female       ~ beta(50, 50);
  rho_pup           ~ beta(2, 5);

  // ── Priors: MOCI effects (identical to Marin IPM v3.3) ──────────────────
  beta_moci_ond_fecund ~ normal(-0.15, 0.20);
  beta_moci_ond_pup    ~ normal(-0.15, 0.20);
  beta_moci_amj_pup    ~ normal(-0.15, 0.20);
  beta_moci_jfm_pup    ~ normal(-0.15, 0.20);
  beta_moci_jfm_juv    ~ normal(-0.15, 0.15);
  beta_moci_jfm_adult  ~ normal(-0.10, 0.12);
  beta_moci_amj_molt   ~ normal(0.05, 0.15);

  // ── Priors: MOCI modifiers by county type ────────────────────────────────
  // N(0, 0.15): weakly informative; shrinks toward zero (no difference from
  // coast) if data are insufficient. Positive values biologically expected.
  // South Bay prior slightly wider to reflect greater expected decoupling.
  delta_moci_mouth_fecund ~ normal(0, 0.15);
  delta_moci_mouth_surv   ~ normal(0, 0.15);
  delta_moci_south_fecund ~ normal(0, 0.20);
  delta_moci_south_surv   ~ normal(0, 0.20);

  // ── Priors: county and site random effects ────────────────────────────────
  sigma_county ~ normal(0.15, 0.10);
  county_effect_raw ~ std_normal();

  // Fix 1: mu_log_adult calibrated to county-level population scale.
  // exp(7.0) ≈ 1,100 adults per county — consistent with raw molt counts
  // of 500-1,500 per county across multiple sites in a single survey period.
  // Previous prior N(5.5, 0.5) (≈245/county) was calibrated to site level,
  // causing the model to underestimate county N by ~5-10x.
  mu_log_adult ~ normal(7.0, 0.5);   // exp(7.0) ≈ 1,096
  mu_log_juv   ~ normal(6.0, 0.5);   // exp(6.0) ≈   403
  mu_log_pup   ~ normal(6.0, 0.5);   // exp(6.0) ≈   403
  sigma_init   ~ exponential(3);
  log_N_adult_F_init_raw ~ std_normal();
  log_N_adult_M_init_raw ~ std_normal();
  log_N_juv_init_raw     ~ std_normal();
  log_N_pup_init_raw     ~ std_normal();

  // Fix 2: mu_log_alpha calibrated per county by number of sites.
  // Prior centers each site at observing 1/n_sites_county of the county pool.
  // Counties with more sites get lower per-site alpha (more sharing).
  //   Marin (8 sites): N(-log(8), 0.5) = N(-2.08, 0.5) → each site ~12%
  //   Bay Mouth (3):   N(-log(3), 0.5) = N(-1.10, 0.5) → each site ~33%
  //   South Bay (2):   N(-log(2), 0.5) = N(-0.69, 0.5) → each site ~50%
  //   San Mateo (6):   N(-log(6), 0.5) = N(-1.79, 0.5) → each site ~17%
  //   Sonoma (4):      N(-log(4), 0.5) = N(-1.39, 0.5) → each site ~25%
  //   Mendocino (1):   N(-log(1), 0.5) = N( 0.00, 0.5) → site sees ~100%
  for (c in 1:C) {
    mu_log_alpha_breed[c] ~ normal(-log(n_sites_county[c]), 0.5);
    mu_log_alpha_pup[c]   ~ normal(-log(n_sites_county[c]) + 0.3, 0.5);
    mu_log_alpha_molt[c]  ~ normal(-log(n_sites_county[c]), 0.5);
  }

  // Fix 3: sigma_log_alpha prior tightened to prevent alpha absorbing all
  // N/alpha confound. Previous N(0.4, 0.2) allowed too much site-to-site
  // variation, letting the model fit data with implausibly small county N.
  sigma_log_alpha ~ normal(0.5, 0.25); //0.3, 0.15);
  log_alpha_breed_raw ~ std_normal();
  log_alpha_pup_raw   ~ std_normal();
  log_alpha_molt_raw  ~ std_normal();

  detect_molt_baseline ~ normal(0.25, 0.50);

  // ── Priors: error terms ───────────────────────────────────────────────────
  // sigma_obs widened slightly from Marin IPM: real data spans 24 sites with
  // more ecological diversity, so more unexplained count variation is expected.
  sigma_process   ~ normal(0.15, 0.08);
  sigma_obs_adult ~ normal(0.25, 0.10);
  sigma_obs_pup   ~ normal(0.20, 0.08);
  sigma_obs_molt  ~ normal(0.35, 0.10);

  to_vector(eps_adult_raw) ~ std_normal();
  to_vector(eps_juv_raw)   ~ std_normal();
  to_vector(eps_pup_raw)   ~ std_normal();

  // ── Vectorized likelihood ────────────────────────────────────────────────
  {
    vector[N_oa] mu_a;
    vector[N_op] mu_p;
    vector[N_om] mu_m;

    for (i in 1:N_oa) {
      int s = oa_s[i]; int t = oa_t[i]; int c = county_id[s];
      real N_adult_obs = N_adult_F[c,t] + N_adult_M[c,t] * p_male_fixed;
      mu_a[i] = log(N_adult_obs) + log_alpha_breed[s];
    }
    for (i in 1:N_op) {
      int s = op_s[i]; int t = op_t[i]; int c = county_id[s];
      mu_p[i] = log(N_pup[c,t]) + log_alpha_pup[s];
    }
    for (i in 1:N_om) {
      int s = om_s[i]; int t = om_t[i]; int c = county_id[s];
      mu_m[i] = log(N_molt_true[c,t]) + log_alpha_molt[s] +
                log(detect_molt_ct[c,t]);
    }

    y_oa ~ normal(mu_a, sigma_obs_adult);
    y_op ~ normal(mu_p, sigma_obs_pup);
    y_om ~ normal(mu_m, sigma_obs_molt);
  }
}

generated quantities {
  // ── County-level lambda (annual growth rate) ─────────────────────────────
  matrix[C, T-1] lambda_county;
  // ── Regional total (sum across all counties) ─────────────────────────────
  vector[T] N_total_regional;
  // ── County totals ─────────────────────────────────────────────────────────
  matrix[C, T] N_total_county;
  // ── Posterior predictive (for PPC) ────────────────────────────────────────
  matrix[S, T] y_adult_rep;
  matrix[S, T] y_pup_rep;
  matrix[S, T] y_molt_rep;

  for (c in 1:C) {
    for (t in 1:T) {
      N_total_county[c, t] = N_total[c, t];
    }
    for (t in 1:(T-1))
      lambda_county[c, t] = N_total[c, t+1] / N_total[c, t];
  }
  for (t in 1:T)
    N_total_regional[t] = sum(col(N_total_county, t));

  // PPCs at site level
  for (s in 1:S) {
    int c = county_id[s];
    for (t in 1:T) {
      real N_adult_obs = N_adult_F[c,t] + N_adult_M[c,t] * p_male_fixed;
      y_adult_rep[s,t] = normal_rng(log(N_adult_obs) + log_alpha_breed[s],
                                     sigma_obs_adult);
      y_pup_rep[s,t]   = normal_rng(log(N_pup[c,t]) + log_alpha_pup[s],
                                     sigma_obs_pup);
      y_molt_rep[s,t]  = normal_rng(log(N_molt_true[c,t]) + log_alpha_molt[s] +
                                     log(detect_molt_ct[c,t]),
                                     sigma_obs_molt);
    }
  }
}
'

write_lines(stan_code_regional, "Code/harbor_seal_regional_ipm.stan")
cat("Stan model written to Code/harbor_seal_regional_ipm.stan\n")

# ── PREPARE REAL DATA FOR STAN ────────────────────────────────────────────────
prepare_real_data_regional <- function(input_rds = "Output/regional_ipm_input_data.rds",
                                       T_proj = 10) {
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
                             input_rds      = "Output/regional_ipm_input_data.rds",
                             prefix         = NULL,   # set by 19_regional_run.R;
                             # defaults to Regional_real/sim
                             T_proj         = 10,
                             seed           = 456,
                             iter_warmup    = 2000,
                             iter_sampling  = 2000,
                             adapt_delta    = 0.97,
                             max_treedepth  = 12) {
  
  cat("\n================================================================\n")
  cat("   REGIONAL HARBOR SEAL IPM (Bay/Coast MOCI + County Dynamics)\n")
  cat("================================================================\n\n")
  
  # Prefix controls ALL output file names (fit RDS, input RDS, CSVs, plots).
  # Set explicitly from 19_regional_run.R to ensure consistency across scripts.
  if (is.null(prefix))
    prefix <- ifelse(use_real_data, "Regional_real", "Regional_sim")
  cat(sprintf("Output prefix: %s\n", prefix))
  
  if (use_real_data) {
    stan_dat <- prepare_real_data_regional(input_rds, T_proj)
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
  
  fit$save_object(paste0("Output/harbor_seal_", prefix, "_fit.rds"))
  saveRDS(list(stan_data = model_data$stan_data, years = model_data$years,
               true_params = model_data$true_params %||% NULL),
          paste0("Output/harbor_seal_", prefix, "_input_data.rds"))
  cat(sprintf("Fit saved: Output/harbor_seal_%s_fit.rds\n", prefix))
  
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
                  "avg_fecundity","rho_pup",
                  "beta_moci_ond_fecund","beta_moci_jfm_adult",
                  "delta_moci_mouth_fecund","delta_moci_mouth_surv",
                  "delta_moci_south_fecund","delta_moci_south_surv",
                  "sigma_county","sigma_log_alpha",
                  "sigma_process","sigma_obs_adult","sigma_obs_pup","sigma_obs_molt")
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

cat("\n17_regional_ipm_model.R loaded.\n")
cat("Key additions vs Marin IPM:\n")
cat("  - County-level population dynamics (C=6)\n")
cat("  - Site availability parameters (alpha_breed, alpha_pup, alpha_molt)\n")
cat("  - MOCI modifiers by county type (mouth/south fecundity + survival)\n")
cat("  - p_male_breed FIXED at 0.057 (Marin IPM posterior)\n")
cat("  - 2020 latent year: state dynamics run, no observations\n")
cat("\nTest on simulated data first:\n")
cat("  source(\"Code/16_regional_simulate.R\")\n")
cat("  out_sim <- run_regional_ipm(use_real_data=FALSE, sim_data=sim_regional)\n")
cat("\nFit real data:\n")
cat("  source(\"Code/15_regional_data_prep.R\")\n")
cat("  out_real <- run_regional_ipm(use_real_data=TRUE)\n")
