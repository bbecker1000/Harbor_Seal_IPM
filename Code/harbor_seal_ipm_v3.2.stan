
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

