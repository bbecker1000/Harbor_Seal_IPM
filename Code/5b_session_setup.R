# =============================================================================
# HARBOR SEAL IPM v3.1 — SESSION SETUP
# Run this script first in every new R session before using
# harbor_seal_ipm_v3.1_standalone_plots.R
#
# What this script builds:
#   fit          — CmdStan fitted model object (loaded from .rds)
#   sim_data     — list with $stan_data, $site_names, $years, $scenario_names
#   portfolio    — list from create_portfolio_analysis_v3.1()
#   sync         — list from create_synchrony_projections_v3.1()
#   proj_results — list of scenario projections (used in plot Section 7)
#   draws        — posterior draws data frame (used by all plot sections)
#   YEARS, SITE_NAMES, SITE_COLORS, S_, T_ — shared constants
# =============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(patchwork)
library(scales)

# ─── Paths ──────────────────────────────────────────────────────────────────
PATH_SCRIPT     <- "Code/5a_IPM_Stan.R"
PATH_FIT        <- "Output/harbor_seal_IPM_v3.1_real_fit.rds"
PATH_DATA       <- "Output/ipm_input_data.rds"
PATH_PORTFOLIO  <- "Output/harbor_seal_IPM_v3.1_real_portfolio.rds"
PATH_SYNC       <- "Output/harbor_seal_IPM_v3.1_real_synchrony.rds"
# ────────────────────────────────────────────────────────────────────────────


# =============================================================================
# STEP 1: Source the IPM script to load all function definitions
# =============================================================================
message("Sourcing IPM v3.1 functions from ", PATH_SCRIPT, "...")
source(PATH_SCRIPT)
message("  ✓ Functions loaded")


# =============================================================================
# STEP 2: Load the fitted model
# =============================================================================
message("Loading fitted model from ", PATH_FIT, "...")
stopifnot(file.exists(PATH_FIT))
fit <- readRDS(PATH_FIT)
stopifnot(!is.null(fit), inherits(fit, "CmdStanMCMC"))
message(sprintf("  ✓ fit loaded  (class: %s)", class(fit)[1]))


# =============================================================================
# STEP 3: Build sim_data
# ipm_input_data.rds contains raw inputs: $dat, $cov_t_scaled, $years
# Handles two cases:
#   (a) file already contains processed list ($stan_data, $site_names, etc.)
#   (b) file contains raw inputs ($dat, $cov_t_scaled, $years) — runs prepare
# =============================================================================
message("Building sim_data from ", PATH_DATA, "...")
stopifnot(file.exists(PATH_DATA))

input_data <- readRDS(PATH_DATA)

if (!is.list(input_data)) {
  stop("Expected a list in ", PATH_DATA, " but got: ", class(input_data))
}

if ("stan_data" %in% names(input_data)) {
  # Already fully processed
  sim_data <- input_data
  message("  ✓ sim_data loaded directly (already processed)")
  
} else if (all(c("dat", "cov_t_scaled", "years") %in% names(input_data))) {
  # Raw inputs — run prepare function
  message("  Raw inputs found ($dat, $cov_t_scaled, $years)")
  message("  Running prepare_real_data_for_ipm_v3.1()...")
  
  prepared <- prepare_real_data_for_ipm_v3.1(
    dat          = input_data$dat,
    cov_t_scaled = input_data$cov_t_scaled,
    years        = input_data$years,
    T_proj       = 10
  )
  
  sim_data <- list(
    stan_data      = prepared$stan_data,
    site_names     = prepared$site_names,
    years          = prepared$years,
    scenario_names = prepared$scenario_names,
    raw_counts     = prepared$raw_counts
  )
  message("  ✓ prepare_real_data_for_ipm_v3.1() complete")
  
} else {
  stop(
    "Unrecognized structure in ", PATH_DATA, "\n",
    "Found names: ", paste(names(input_data), collapse = ", "), "\n",
    "Expected either '$stan_data' (processed) or '$dat/$cov_t_scaled/$years' (raw)."
  )
}

# Verify
stopifnot(
  !is.null(sim_data$stan_data),
  !is.null(sim_data$years),
  length(sim_data$years) > 0
)
message(sprintf("  ✓ sim_data ready: %d sites × %d years (%d–%d)",
                length(sim_data$site_names),
                length(sim_data$years),
                min(sim_data$years),
                max(sim_data$years)))


# =============================================================================
# STEP 4: Shared constants (used by every plot section)
# =============================================================================
SITE_NAMES  <- sim_data$site_names      # c("BL","DE","DP","PRH","TB","TP")
YEARS       <- sim_data$years           # numeric year vector
T_          <- length(YEARS)
S_          <- length(SITE_NAMES)

SITE_COLORS <- c(
  "BL"  = "#1b7837",
  "DE"  = "#762a83",
  "DP"  = "#d01c8b",
  "PRH" = "#f1b6da",
  "TB"  = "#4393c3",
  "TP"  = "#e08214"
)

message(sprintf("  ✓ Constants: S=%d, T=%d", S_, T_))


# =============================================================================
# STEP 5: Extract posterior draws (reused by all plot sections)
# =============================================================================
message("Extracting posterior draws...")
draws <- fit$draws(format = "df")
message(sprintf("  ✓ draws: %d rows × %d columns", nrow(draws), ncol(draws)))


# =============================================================================
# STEP 6: Portfolio analysis
# Loads from cache if available; otherwise runs and caches.
# =============================================================================
if (file.exists(PATH_PORTFOLIO)) {
  message("Loading cached portfolio analysis...")
  portfolio <- readRDS(PATH_PORTFOLIO)
  message("  ✓ portfolio loaded from cache")
} else {
  message("Running portfolio analysis (~1–2 min)...")
  portfolio <- create_portfolio_analysis_v3.1(
    fit      = fit,
    sim_data = sim_data,
    prefix   = "IPM_v3.1_real",
    save     = TRUE
  )
  saveRDS(portfolio, PATH_PORTFOLIO)
  message("  ✓ portfolio complete — saved to ", PATH_PORTFOLIO)
}


# =============================================================================
# STEP 7: Synchrony projections
# Loads from cache if available; otherwise runs and caches.
# =============================================================================
if (file.exists(PATH_SYNC)) {
  message("Loading cached synchrony projections...")
  sync <- readRDS(PATH_SYNC)
  message("  ✓ sync loaded from cache")
} else {
  message("Running synchrony projections (~2–3 min)...")
  sync <- create_synchrony_projections_v3.1(
    fit      = fit,
    sim_data = sim_data,
    n_sims   = 500,
    prefix   = "IPM_v3.1_real",
    save     = TRUE
  )
  saveRDS(sync, PATH_SYNC)
  message("  ✓ sync complete — saved to ", PATH_SYNC)
}


# =============================================================================
# STEP 8: Scenario projections for plot Section 7
# Extracts N_proj[scenario, site, year_idx] from Stan generated quantities
# and reshapes into a named list of per-scenario summary data frames.
# =============================================================================
message("Extracting scenario projections from draws...")

T_proj     <- sim_data$stan_data$T_proj %||% 10L
proj_years <- max(YEARS) + seq_len(T_proj)
n_scen     <- sim_data$stan_data$N_scenarios %||% 4L
scen_names <- sim_data$scenario_names

proj_list <- vector("list", n_scen)

for (sc in seq_len(n_scen)) {
  expected  <- paste0(
    "N_proj[", sc, ",",
    rep(seq_len(S_), each = T_proj), ",",
    rep(seq_len(T_proj), S_), "]"
  )
  params_sc <- intersect(expected, colnames(draws))
  
  if (length(params_sc) == 0) {
    message(sprintf("  ⚠  N_proj[%d,...] not found in draws — skipping", sc))
    next
  }
  
  proj_list[[sc]] <- draws |>
    select(all_of(params_sc)) |>
    mutate(.draw = row_number()) |>
    pivot_longer(-.draw, names_to = "param", values_to = "N") |>
    mutate(
      year_idx = as.integer(str_extract(param, "(?<=,)\\d+(?=\\])")),
      Year     = proj_years[year_idx]
    ) |>
    group_by(.draw, Year) |>
    summarise(N_total = sum(N), .groups = "drop") |>
    group_by(Year) |>
    summarise(
      median = median(N_total),
      lo95   = quantile(N_total, 0.025),
      hi95   = quantile(N_total, 0.975),
      lo50   = quantile(N_total, 0.25),
      hi50   = quantile(N_total, 0.75),
      .groups = "drop"
    )
}

names(proj_list) <- scen_names

proj_results <- list(
  status_quo     = proj_list[[1]],
  coyote_removal = proj_list[[2]],
  warm_moci      = proj_list[[3]],
  combined       = proj_list[[4]]
)

n_ok <- sum(!sapply(proj_list, is.null))
message(sprintf("  ✓ proj_results: %d/%d scenarios extracted", n_ok, n_scen))


# =============================================================================
# STEP 9: Parameter name crosswalk
# Your Stan model v3.1 uses the names on the LEFT.
# The plot script expects the aliases on the RIGHT.
# This block adds alias columns to draws so both names work.
#
#  Actual v3.1 Stan name       →  Plot script alias
# ─────────────────────────────────────────────────────
#  phi_pup_F_base              →  phi_pup_base
#  phi_juv_F_base              →  phi_juv_base
#  phi_adult_F_base            →  phi_adult_base
#  delta_adult                 →  delta_sex
#  fecund_primip               →  fecundity_prima
#  fecund_young                →  fecundity_multi
#  fecund_prime                →  fecundity_prime
#  sigma_obs_adult             →  sigma_obs
#  beta_coy[1]                 →  beta_coy_BL
#  beta_coy[2]                 →  beta_coy_DE
#  beta_coy[3]                 →  beta_coy_DP
#  beta_dist_surv[1..6]        →  beta_dist_surv_1..6
# =============================================================================
message("Applying parameter name crosswalk...")

PARAM_CROSSWALK <- c(
  "phi_pup_F_base"   = "phi_pup_base",
  "phi_juv_F_base"   = "phi_juv_base",
  "phi_adult_F_base" = "phi_adult_base",
  "delta_adult"      = "delta_sex",
  "fecund_primip"    = "fecundity_prima",
  "fecund_young"     = "fecundity_multi",
  "fecund_prime"     = "fecundity_prime",
  "sigma_obs_adult"  = "sigma_obs"
)

for (stan_nm in names(PARAM_CROSSWALK)) {
  plot_nm <- PARAM_CROSSWALK[[stan_nm]]
  if (stan_nm %in% colnames(draws) && !plot_nm %in% colnames(draws))
    draws[[plot_nm]] <- draws[[stan_nm]]
}

# beta_coy[n] → beta_coy_SITE
for (i in seq_along(c("BL", "DE", "DP"))) {
  old <- paste0("beta_coy[", i, "]")
  new <- paste0("beta_coy_", c("BL", "DE", "DP")[i])
  if (old %in% colnames(draws) && !new %in% colnames(draws))
    draws[[new]] <- draws[[old]]
}

# beta_dist_surv[n] → beta_dist_surv_n
for (i in seq_len(S_)) {
  old <- paste0("beta_dist_surv[", i, "]")
  new <- paste0("beta_dist_surv_", i)
  if (old %in% colnames(draws) && !new %in% colnames(draws))
    draws[[new]] <- draws[[old]]
}

message("  ✓ Crosswalk applied")


# =============================================================================
# STEP 10: Final verification
# =============================================================================
required <- c("fit", "sim_data", "portfolio", "sync",
              "proj_results", "draws",
              "YEARS", "SITE_NAMES", "SITE_COLORS", "S_", "T_")

missing <- required[!sapply(required, exists, envir = .GlobalEnv)]

if (length(missing) == 0) {
  message("\n============================================")
  message("SESSION SETUP COMPLETE — all objects ready")
  message("============================================")
  message(sprintf("  fit          CmdStanMCMC  (%d chains)", fit$num_chains()))
  message(sprintf("  sim_data     %d sites × %d years (%d–%d)",
                  S_, T_, min(YEARS), max(YEARS)))
  message(sprintf("  draws        %d posterior draws × %d params",
                  nrow(draws), ncol(draws)))
  message(sprintf("  portfolio    %d list elements", length(portfolio)))
  message(sprintf("  sync         %d list elements", length(sync)))
  message(sprintf("  proj_results %d scenarios", length(proj_results)))
  message("\nNext step:")
  message("  source('harbor_seal_ipm_v3.1_standalone_plots.R')")
  message("  — or run any individual plot section.\n")
} else {
  warning("The following required objects are missing:\n  ",
          paste(missing, collapse = "\n  "))
}