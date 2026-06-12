# MOCI Correlations 2026-06-11


# ── 1. Empirical correlation among MOCI series (raw covariates) ─────────────
# ── 1. Empirical correlation among MOCI series (raw covariates) ─────────────
sd <- out$sim_data$stan_data   # try this first

# If sd$moci_jfm is NULL, try the alternative storage location:
if (is.null(sd$moci_jfm)) {
  res <- readRDS("Output/harbor_seal_IPM_v3.2_real_input_data.rds")
  sd  <- res$stan_data
}

# Confirm it's there
str(sd$moci_jfm)
str(sd$moci_amj)
str(sd$moci_ond)
# More reliably, pull directly from the input data used to fit:
sd <- out$sim_data$stan_data  # or res$data$stan_data depending on object

T <- length(sd$moci_jfm)

# Build the relevant lagged/unlagged series used across the 7 pup/juv/adult/fecund MOCI terms
moci_df <- data.frame(
  jfm_t      = sd$moci_jfm,                 # JFM(t) -> pup, juv, adult survival
  amj_t      = sd$moci_amj,                 # AMJ(t) -> molt detection
  amj_tm1    = c(NA, sd$moci_amj[1:(T-1)]), # AMJ(t-1) -> pup survival (birth year)
  ond_t      = sd$moci_ond,                 # OND(t) -> fecundity, pup survival
  ond_tm1    = c(NA, sd$moci_ond[1:(T-1)])  # OND(t-1), if relevant to any birth-year term
)

cat("── Correlation matrix among MOCI predictor series ──────────────\n")
print(round(cor(moci_df, use="pairwise.complete.obs"), 3))


# ── 2. Posterior correlation among the corresponding beta parameters ────────
beta_vars <- c("beta_moci_jfm_pup", "beta_moci_amj_pup", "beta_moci_ond_pup",
               "beta_moci_jfm_juv", "beta_moci_jfm_adult",
               "beta_moci_ond_fecund", "beta_moci_amj_molt")

beta_draws <- out$fit$draws(variables = beta_vars, format = "df") |>
  dplyr::select(dplyr::all_of(beta_vars))

cat("\n── Posterior correlation matrix among MOCI beta parameters ─────\n")
print(round(cor(beta_draws), 3))

# ── 3. Aggregate "any-pathway" pup-survival MOCI signal ──────────────────────
# Sum of the three pup-survival MOCI betas per draw -- tests whether the
# combined seasonal signal is better identified than any individual term
pup_moci_sum <- beta_draws$beta_moci_jfm_pup +
  beta_draws$beta_moci_amj_pup +
  beta_draws$beta_moci_ond_pup

cat("\n── Aggregate pup-survival MOCI effect (sum of 3 seasonal terms) ─\n")
cat(sprintf("Median: %.3f, 89%% CI: [%.3f, %.3f]\n",
            median(pup_moci_sum),
            quantile(pup_moci_sum, 0.055),
            quantile(pup_moci_sum, 0.945)))