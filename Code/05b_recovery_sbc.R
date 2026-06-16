# ============================================================================
# 05b_recovery_sbc.R  —  PROPER PARAMETER-RECOVERY CHECK (Option B, fixes P2)
# ----------------------------------------------------------------------------
# Replaces the old "truth = prior mean" recovery (which over-covered by
# construction) with an SBC-lite design:
#   for each of n_datasets:
#     1. draw true parameter values FROM the priors
#     2. simulate data under those truths
#     3. fit the IPM
#     4. record whether each truth falls in its 89% CrI (+ posterior rank)
#   aggregate -> per-parameter coverage rate (should be ~0.89 if calibrated).
#
# Requires the 2-line edit to simulate_seal_ipm_data_v3.2() (accepts
# true_params=NULL). Source order:
#   source("Code/05_ipm_model.R"); source("Code/05b_recovery_sbc.R")
#
# RUNTIME: one fit ~ a few minutes; n_datasets fits run sequentially. Start
# small (n_datasets=10) to gauge, then scale up overnight for tighter coverage.
# ============================================================================

library(tidyverse)
library(cmdstanr)
library(posterior)

`%||%` <- function(x, y) if (!is.null(x)) x else y
.CI_LO <- 0.055; .CI_HI <- 0.945

# ── Draw a true_params list FROM the Stan priors (with the model's bounds) ───
draw_true_params_from_priors <- function(seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  rtnorm <- function(m, s, lo=-Inf, hi=Inf) {
    repeat { x <- rnorm(1, m, s); if (x >= lo && x <= hi) return(x) } }
  rtbeta <- function(a, b, lo=0, hi=1) {
    repeat { x <- rbeta(1, a, b); if (x >= lo && x <= hi) return(x) } }
  list(
    phi_pup_logit      = rnorm(1, -1.2, 0.5),
    phi_juv_base       = rbeta(1, 14, 6),
    phi_adult_F_logit  = rnorm(1, 2.20, 0.25),
    delta_adult        = rtnorm(0.05, 0.025, 0, 0.10),
    fecund_primip      = rbeta(1, 12, 8),
    fecund_mature      = rbeta(1, 17, 3),
    prop_female        = rtbeta(50, 50, 0.4, 0.6),
    p_male_breed       = rtbeta(2, 18, 0, 0.30),
    detect_breed_logit = rnorm(1, 1.20, 0.50),
    detect_molt_logit  = rnorm(1, 0.75, 0.50),
    beta_coy           = rnorm(3, -0.20, 0.20),
    beta_dist_surv     = rnorm(6, -0.15, 0.20),
    beta_dist_detect   = rnorm(6, -0.15, 0.15),
    beta_moci_ond_fecund = rnorm(1, -0.15, 0.20),
    beta_moci_ond_pup    = rnorm(1, -0.15, 0.20),
    beta_moci_amj_pup    = rnorm(1, -0.15, 0.20),
    beta_moci_jfm_pup    = rnorm(1, -0.15, 0.20),
    beta_moci_jfm_juv    = rnorm(1, -0.15, 0.15),
    beta_moci_jfm_adult  = rnorm(1, -0.10, 0.12),
    beta_moci_amj_molt   = rnorm(1,  0.05, 0.15),
    beta_eseal_pup       = rnorm(1,  0.10, 0.20),
    sigma_process   = rtnorm(0.15, 0.08, 0.05, 0.50),
    sigma_obs_adult = rtnorm(0.18, 0.06, 0.05, 0.40),
    sigma_obs_pup   = rtnorm(0.15, 0.02, 0.02, 0.35),
    sigma_obs_molt  = rtnorm(0.35, 0.10, 0.05, 0.60),
    sigma_site      = rtnorm(0.2,  0.1,  0.01, 0.50)
  )
}

# ── Flatten a true_params list to a (variable, true_value) tibble ───────────
# Order/names match the 30 monitored parameters in check_parameter_recovery.
.flatten_true <- function(tp) {
  tibble(
    variable = c(
      "phi_pup_logit","phi_juv_base","phi_adult_F_logit","delta_adult",
      "fecund_primip","fecund_mature","prop_female","p_male_breed",
      paste0("beta_coy[",1:3,"]"),
      paste0("beta_dist_surv[",1:6,"]"),
      "beta_moci_ond_fecund","beta_moci_ond_pup","beta_moci_jfm_pup",
      "beta_moci_amj_pup","beta_moci_jfm_juv","beta_moci_jfm_adult",
      "beta_eseal_pup","detect_breed_logit","detect_molt_logit",
      "sigma_process","sigma_obs_adult","sigma_obs_pup","sigma_obs_molt"),
    true_value = c(
      tp$phi_pup_logit, tp$phi_juv_base, tp$phi_adult_F_logit, tp$delta_adult,
      tp$fecund_primip, tp$fecund_mature, tp$prop_female, tp$p_male_breed,
      tp$beta_coy[1], tp$beta_coy[2], tp$beta_coy[3],
      tp$beta_dist_surv[1], tp$beta_dist_surv[2], tp$beta_dist_surv[3],
      tp$beta_dist_surv[4], tp$beta_dist_surv[5], tp$beta_dist_surv[6],
      tp$beta_moci_ond_fecund, tp$beta_moci_ond_pup, tp$beta_moci_jfm_pup,
      tp$beta_moci_amj_pup, tp$beta_moci_jfm_juv, tp$beta_moci_jfm_adult,
      tp$beta_eseal_pup, tp$detect_breed_logit, tp$detect_molt_logit,
      tp$sigma_process, tp$sigma_obs_adult, tp$sigma_obs_pup, tp$sigma_obs_molt)
  )
}

# Wilson score interval for a coverage proportion
.wilson <- function(k, n, z = 1.96) {
  if (n == 0) return(c(lo = NA, hi = NA))
  p <- k / n; d <- 1 + z^2/n
  c <- (p + z^2/(2*n)) / d
  h <- z * sqrt(p*(1-p)/n + z^2/(4*n^2)) / d
  c(lo = max(0, c - h), hi = min(1, c + h))
}

# ── Main SBC-lite runner ────────────────────────────────────────────────────
run_recovery_check_v3.2 <- function(n_datasets    = 10,
                                    iter_warmup   = 750,
                                    iter_sampling = 500,
                                    adapt_delta   = 0.90,
                                    max_treedepth = 12,
                                    base_seed     = 1000,
                                    rhat_max      = 1.05,
                                    save          = TRUE,
                                    prefix        = "IPM_v3.2_sbc") {
  
  if (!exists("simulate_seal_ipm_data_v3.2", mode = "function"))
    stop("Source 05_ipm_model.R first (need simulate_seal_ipm_data_v3.2 with ",
         "the true_params=NULL edit).")
  
  stan_path <- if (file.exists("Code/harbor_seal_ipm_v3.2.stan"))
    "Code/harbor_seal_ipm_v3.2.stan" else "harbor_seal_ipm_v3.2.stan"
  model <- cmdstan_model(stan_path)
  
  cat(sprintf("\nSBC-lite recovery: %d datasets x (warmup %d + sampling %d)\n",
              n_datasets, iter_warmup, iter_sampling))
  cat("Truths drawn FROM priors (not prior means) — expect coverage ~89%.\n\n")
  
  per_dataset <- vector("list", n_datasets)
  rhat_by_ds  <- numeric(n_datasets)
  
  for (d in seq_len(n_datasets)) {
    sd_d <- base_seed + d
    tp   <- draw_true_params_from_priors(seed = sd_d)
    sim  <- simulate_seal_ipm_data_v3.2(seed = sd_d + 5000, true_params = tp)
    tv   <- .flatten_true(tp)
    
    fit <- model$sample(
      data = sim$stan_data, seed = sd_d, chains = 4, parallel_chains = 4,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      refresh = 0, adapt_delta = adapt_delta, max_treedepth = max_treedepth)
    
    drw <- fit$draws(variables = tv$variable, format = "matrix")
    drw <- drw[, tv$variable, drop = FALSE]   # enforce column order
    qlo <- apply(drw, 2, quantile, .CI_LO)
    qhi <- apply(drw, 2, quantile, .CI_HI)
    pm  <- colMeans(drw)
    # posterior rank of the truth (for optional SBC uniformity at large N)
    rank <- vapply(seq_along(tv$variable),
                   function(j) mean(drw[, j] < tv$true_value[j]), numeric(1))
    rh   <- tryCatch(suppressWarnings(max(fit$summary(tv$variable)$rhat,
                                          na.rm = TRUE)), error = function(e) NA)
    rhat_by_ds[d] <- rh
    
    per_dataset[[d]] <- tibble(
      dataset = d, variable = tv$variable, true_value = tv$true_value,
      post_mean = pm, q_lo = qlo, q_hi = qhi,
      covered = tv$true_value >= qlo & tv$true_value <= qhi,
      rank = rank, max_rhat = rh)
    
    cat(sprintf("  dataset %2d/%d  seed %d  coverage %2d/%d  maxRhat %.3f\n",
                d, n_datasets, sd_d, sum(per_dataset[[d]]$covered),
                nrow(tv), rh))
  }
  
  detail <- bind_rows(per_dataset)
  n_conv <- sum(rhat_by_ds <= rhat_max, na.rm = TRUE)
  cat(sprintf("\nConverged datasets (maxRhat <= %.2f): %d/%d\n",
              rhat_max, n_conv, n_datasets))
  
  # ── Per-parameter coverage (all datasets, and converged-only) ─────────────
  conv_ids <- which(rhat_by_ds <= rhat_max)
  cover_tbl <- detail |>
    group_by(variable) |>
    summarise(n = n(), k = sum(covered),
              coverage = k / n,
              n_conv = sum(dataset %in% conv_ids),
              k_conv = sum(covered & dataset %in% conv_ids),
              coverage_conv = ifelse(n_conv > 0, k_conv / n_conv, NA),
              .groups = "drop") |>
    rowwise() |>
    mutate(ci_lo = .wilson(k, n)["lo"], ci_hi = .wilson(k, n)["hi"]) |>
    ungroup() |>
    arrange(coverage)
  
  overall <- mean(detail$covered)
  overall_conv <- mean(detail$covered[detail$dataset %in% conv_ids])
  
  cat(sprintf("\nOverall coverage: %.1f%% (all)  |  %.1f%% (converged)\n",
              100*overall, 100*overall_conv))
  cat("Nominal target: 89%\n\n")
  print(cover_tbl |> select(variable, n, coverage, coverage_conv, ci_lo, ci_hi),
        n = nrow(cover_tbl))
  
  # ── Coverage plot ──────────────────────────────────────────────────────────
  p_cov <- ggplot(cover_tbl,
                  aes(x = reorder(variable, coverage), y = coverage)) +
    geom_hline(yintercept = 0.89, linetype = 2, colour = "grey40") +
    geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi,
                        colour = (ci_lo <= 0.89 & ci_hi >= 0.89))) +
    scale_colour_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
                        name = "89% nominal\nin Wilson CI") +
    coord_flip(ylim = c(0, 1)) +
    labs(x = NULL, y = "Coverage of 89% CrI across simulated datasets",
         title = "Parameter Recovery — SBC-lite (truths drawn from priors)",
         subtitle = sprintf("%d datasets; overall coverage %.1f%% (nominal 89%%). Dashed = target.",
                            n_datasets, 100 * overall)) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.major.y = element_blank())
  
  if (save) {
    dir.create("Output/Plots", showWarnings = FALSE)
    write_csv(detail,    paste0("Output/", prefix, "_recovery_detail.csv"))
    write_csv(cover_tbl, paste0("Output/", prefix, "_recovery_coverage.csv"))
    ggsave(paste0("Output/Plots/", prefix, "_recovery_coverage.jpeg"),
           p_cov, width = 24, height = 22, units = "cm", dpi = 200)
  }
  
  list(detail = detail, coverage = cover_tbl, plot = p_cov,
       overall_coverage = overall, overall_coverage_converged = overall_conv,
       n_converged = n_conv, rhat_by_dataset = rhat_by_ds)
}

cat("05b_recovery_sbc.R loaded — run_recovery_check_v3.2(n_datasets=10)\n")