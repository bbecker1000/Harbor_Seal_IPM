
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

