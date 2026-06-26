# ============================================================================
# 18_regional_plots.R
# ----------------------------------------------------------------------------
# REGIONAL HARBOR SEAL IPM — PLOTS, TABLES & POST-PROCESSING
#
# Parallel to 06_ipm_plots.R but adapted for the regional model structure:
#   - County-level (not site-level) population trajectories
#   - Bay vs. Coast MOCI modifier posteriors (the key new inferential target)
#   - Cross-county lambda correlation and portfolio analysis
#   - MOCI response curves stratified by coast type
#   - Projection scenarios by county and coast type
#
# Functions:
#   load_regional_results()               — reload fit + reconstruct model_data
#   run_all_regional_plots()              — orchestrator
#   check_diagnostics_regional()          — convergence + key parameters
#   create_county_trajectory_plots()      — Q1: is decline regional?
#   create_bay_coast_moci_plots()         — Q2: Bay vs Coast MOCI response
#   create_regional_projection_plots()    — Q3: projections by county/coast
#   create_regional_portfolio_plots()     — Q4: cross-county lambda correlation
#   create_regional_forest_plot()         — all MOCI + bay modifier posteriors
#   create_site_availability_plots()      — alpha parameter distributions
#   create_regional_summary_table()       — parameter summary CSV/print
#
# Usage:
#   source("Code/18_regional_plots.R")
#   out <- load_regional_results("Regional_real")
#   run_all_regional_plots(out$fit, out$model_data, prefix="Regional_real")
# ============================================================================

library(tidyverse)
library(posterior)
library(bayesplot)
library(patchwork)

dir.create("Output",       showWarnings = FALSE)
dir.create("Output/Plots", showWarnings = FALSE)

`%||%` <- function(x, y) if (!is.null(x)) x else y

CI_LO    <- 0.055
CI_HI    <- 0.945
CI_LABEL <- "89% CrI"

COUNTY_NAMES  <- c("Marin", "Bay Mouth", "South Bay", "San Mateo", "Sonoma", "Mendocino")
COUNTY_TYPE   <- c(0L, 1L, 2L, 0L, 0L, 0L)   # 0=coast, 1=bay mouth, 2=south bay

COUNTY_COLS <- c(
  "Marin"      = "#08519C",
  "Bay Mouth"  = "#F16913",   # medium orange — bay mouth (intermediate)
  "South Bay"  = "#D94801",   # deep orange  — most estuarine
  "San Mateo"  = "#6BAED6",
  "Sonoma"     = "#41AB5D",
  "Mendocino"  = "#74C476"
)
# Three coast types for MOCI response plots
COAST_TYPE_COLS <- c(
  "Open Coast"  = "#2166AC",
  "Bay Mouth"   = "#F16913",
  "South Bay"   = "#D94801"
)

theme_seal <- function(base_size = 14) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.grid.major   = element_line(colour = "grey88", linewidth = 0.4),
      panel.grid.minor   = element_line(colour = "grey93", linewidth = 0.2),
      panel.border       = element_rect(colour = "grey70", fill = NA, linewidth = 0.5),
      axis.title         = element_text(size = rel(0.95), colour = "grey20"),
      axis.text          = element_text(size = rel(0.88), colour = "grey30"),
      legend.position    = "bottom",
      legend.title       = element_text(size = rel(0.88), face = "bold"),
      legend.text        = element_text(size = rel(0.82)),
      strip.text         = element_text(size = rel(0.88), face = "bold", colour = "grey20"),
      strip.background   = element_rect(fill = "grey94", colour = "grey80"),
      plot.title         = element_text(size = rel(1.05), face = "bold", margin = margin(b = 6)),
      plot.subtitle      = element_text(size = rel(0.88), colour = "grey40", margin = margin(b = 8)),
      plot.margin        = margin(10, 14, 10, 10)
    )
}

# ── LOAD HELPER ───────────────────────────────────────────────────────────────
load_regional_results <- function(prefix        = "Regional_real",
                                  fit_path      = NULL,
                                  input_rds     = NULL,
                                  years         = 2005:2025) {
  if (is.null(fit_path))
    fit_path <- paste0("Output/harbor_seal_", prefix, "_fit.rds")
  if (!file.exists(fit_path)) stop("Fit file not found: ", fit_path)
  cat("Loading fit from", fit_path, "...\n")
  fit <- readRDS(fit_path)
  
  if (is.null(input_rds))
    input_rds <- paste0("Output/harbor_seal_", prefix, "_input_data.rds")
  
  if (file.exists(input_rds)) {
    inp <- readRDS(input_rds)
    model_data <- list(
      stan_data      = inp$stan_data,
      county_names   = COUNTY_NAMES,
      years          = inp$years %||% years,
      scenario_names = c("Status Quo", "Warm (MOCI +1)", "Cool (MOCI -1)"),
      true_params    = inp$true_params %||% NULL
    )
    # Try to reload site metadata
    reg_rds <- "Output/regional_ipm_input_data.rds"
    if (file.exists(reg_rds))
      model_data$site_meta <- readRDS(reg_rds)$site_meta
  } else {
    cat("No input_rds found — using minimal metadata shell.\n")
    model_data <- list(
      stan_data      = list(T = length(years), C = 6, S = 22),
      county_names   = COUNTY_NAMES,
      years          = years,
      scenario_names = c("Status Quo", "Warm (MOCI +1)", "Cool (MOCI -1)"),
      true_params    = NULL
    )
  }
  cat("Ready. run_all_regional_plots(out$fit, out$model_data) to generate outputs.\n")
  list(fit = fit, model_data = model_data)
}

# ── FIT SUMMARY HELPER ────────────────────────────────────────────────────────
regional_fit_summary <- function(fit, variables) {
  s <- fit$summary(variables = variables)
  d <- fit$draws(variables = variables, format = "matrix")
  s$q_lo <- apply(d[, s$variable, drop = FALSE], 2, quantile, CI_LO)
  s$q_hi <- apply(d[, s$variable, drop = FALSE], 2, quantile, CI_HI)
  s
}

# ============================================================================
# ORCHESTRATOR
# ============================================================================
run_all_regional_plots <- function(fit, model_data, prefix = "Regional",
                                   save = TRUE) {
  safe <- function(label, expr) {
    cat(sprintf("\n%s\n", label))
    tryCatch(
      withCallingHandlers(expr, warning = function(w) invokeRestart("muffleWarning")),
      error = function(e) { cat(sprintf("  !! FAILED: %s\n", conditionMessage(e))); NULL })
  }
  
  diag  <- safe("── Diagnostics ──",           check_diagnostics_regional(fit))
  traj  <- safe("── County trajectories ──",   create_county_trajectory_plots(fit, model_data, save, prefix))
  bay   <- safe("── Bay/Coast MOCI ──",        create_bay_coast_moci_plots(fit, model_data, save, prefix))
  proj  <- safe("── Projections ──",           create_regional_projection_plots(fit, model_data, save, prefix))
  port  <- safe("── Portfolio ──",             create_regional_portfolio_plots(fit, model_data, save, prefix))
  fplot <- safe("── Forest plot ──",           create_regional_forest_plot(fit, save, prefix))
  alpha <- safe("── Site availability ──",     create_site_availability_plots(fit, model_data, save, prefix))
  tbl   <- safe("── Summary table ──",         create_regional_summary_table(fit, save, prefix))
  
  results <- list(diagnostics = diag, trajectories = traj, bay_coast = bay,
                  projections = proj, portfolio = port, forest = fplot,
                  availability = alpha, table = tbl)
  n_ok <- sum(!sapply(results, is.null))
  cat(sprintf("\n── Pipeline complete: %d/%d steps produced output ──\n",
              n_ok, length(results)))
  invisible(results)
}

# ============================================================================
# DIAGNOSTICS
# ============================================================================
check_diagnostics_regional <- function(fit) {
  cat("\n=== MODEL DIAGNOSTICS — REGIONAL IPM ===\n")
  tryCatch({
    ds <- fit$diagnostic_summary(quiet = TRUE)
    cat(if (all(ds$num_max_treedepth == 0)) "Treedepth OK.\n"
        else sprintf("WARN: %d treedepth hits.\n", sum(ds$num_max_treedepth)))
    cat(if (all(ds$num_divergent == 0)) "No divergences.\n"
        else sprintf("WARN: %d divergences.\n", sum(ds$num_divergent)))
    cat(if (all(ds$ebfmi > 0.2)) "E-BFMI OK.\n"
        else sprintf("WARN: low E-BFMI in %d chains.\n", sum(ds$ebfmi <= 0.2)))
  }, error = function(e) cat("Diagnostics unavailable (CSVs gone).\n"))
  
  params <- c("phi_pup_logit", "phi_juv_base", "phi_adult_F_logit", "delta_adult",
              "avg_fecundity", "rho_pup",
              "beta_moci_ond_fecund", "beta_moci_jfm_adult",
              "delta_moci_mouth_fecund", "delta_moci_mouth_surv",
              "delta_moci_south_fecund", "delta_moci_south_surv",
              "sigma_county", "sigma_site", "detect_breed_logit", "detect_molt_logit",
              "sigma_process", "sigma_obs_adult", "sigma_obs_pup", "sigma_obs_molt")
  
  s <- tryCatch(regional_fit_summary(fit, params), error = function(e) NULL)
  if (!is.null(s)) {
    cat("\nParameter Summary:\n")
    print(s |> select(variable, mean, sd, q_lo, q_hi, rhat, ess_bulk), n = nrow(s))
    
    # Flag all four modifier parameters explicitly
    for (p in c("delta_moci_mouth_fecund","delta_moci_mouth_surv",
                "delta_moci_south_fecund","delta_moci_south_surv")) {
      r <- s[s$variable == p, ]
      if (nrow(r) > 0)
        cat(sprintf("%s: %.3f (%s: %.3f, %.3f)\n",
                    p, r$mean, CI_LABEL, r$q_lo, r$q_hi))
    }
  }
  list(params = params, summary = s)
}

# ============================================================================
# Q1: COUNTY TRAJECTORY PLOTS
# Addresses: "Is the regional population also declining like Marin?"
# ============================================================================
create_county_trajectory_plots <- function(fit, model_data, save = TRUE, prefix = "Regional") {
  years <- model_data$years
  T <- length(years); C <- model_data$stan_data$C %||% 6
  county_names <- model_data$county_names %||% COUNTY_NAMES
  # Guard: if county_names length differs from C, use generic labels
  if (length(county_names) != C)
    county_names <- paste0("County_", seq_len(C))
  
  # ── (a) County N_total trajectories ────────────────────────────────────────
  draws_all <- fit$draws(format = "matrix")
  
  traj_df <- map_dfr(1:C, function(c) {
    cols <- paste0("N_total_county[", c, ",", 1:T, "]")
    cols_ok <- cols[cols %in% colnames(draws_all)]
    if (!length(cols_ok)) return(NULL)
    mat <- draws_all[, cols_ok, drop = FALSE]
    tibble(
      County = county_names[c],
      Year   = years[seq_along(cols_ok)],
      mean   = colMeans(mat),
      lo     = apply(mat, 2, quantile, CI_LO),
      hi     = apply(mat, 2, quantile, CI_HI),
      coast_type = case_when(COUNTY_TYPE[c] == 1 ~ "Bay Mouth",
                             COUNTY_TYPE[c] == 2 ~ "South Bay",
                             TRUE                ~ "Open Coast")
    )
  }) |> mutate(County = factor(County, levels = county_names))
  
  # ── Survey coverage rug: which years had data in each county ────────────────
  y_adult_obs_mat <- model_data$stan_data$y_adult_obs   # [S, T]
  county_id_vec   <- model_data$stan_data$county_id     # [S]
  rug_df <- map_dfr(seq_len(C), function(ci) {
    sites_in_c <- which(county_id_vec == ci)
    has_survey  <- apply(y_adult_obs_mat[sites_in_c, , drop = FALSE], 2, any)
    tibble(County = county_names[ci],
           Year   = years[has_survey])
  }) |> mutate(County = factor(County, levels = county_names))
  
  p_traj <- ggplot(traj_df, aes(x = Year, y = mean, colour = County, fill = County)) +
    geom_vline(xintercept = 2020, linetype = "dotted", colour = "grey50", linewidth = 0.8) +
    annotate("text", x = 2020.2, y = Inf, label = "2020\n(latent)", vjust = 1.3,
             hjust = 0, size = 3, colour = "grey50") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1.1) +
    geom_rug(data = rug_df, aes(x = Year), sides = "b",
             colour = "grey35", linewidth = 0.5, length = unit(0.04, "npc"),
             inherit.aes = FALSE) +
    scale_colour_manual(values = COUNTY_COLS) +
    scale_fill_manual(values = COUNTY_COLS) +
    scale_y_continuous(labels = scales::comma, expand = c(0.04, 0)) +
    facet_wrap(~ County, scales = "free_y", ncol = 3) +
    labs(x = "Year", y = "Total population (all ages)",
         title = "County Population Trajectories (Historical)",
         subtitle = paste0("Posterior mean + ", CI_LABEL,
                           "; tick marks = years with survey data")) +
    theme_seal() + guides(colour = "none", fill = "none")
  if (save) ggsave(paste0("Output/Plots/", prefix, "_county_trajectories.jpeg"),
                   p_traj, width = 34, height = 22, units = "cm", dpi = 200)
  
  # ── (b) Regional total (all counties combined) ─────────────────────────────
  reg_cols <- paste0("N_total_regional[", 1:T, "]")
  reg_cols_ok <- reg_cols[reg_cols %in% colnames(draws_all)]
  if (length(reg_cols_ok) > 0) {
    reg_mat <- draws_all[, reg_cols_ok, drop = FALSE]
    reg_df  <- tibble(Year = years[seq_along(reg_cols_ok)],
                      mean = colMeans(reg_mat),
                      lo   = apply(reg_mat, 2, quantile, CI_LO),
                      hi   = apply(reg_mat, 2, quantile, CI_HI))
    p_total <- ggplot(reg_df, aes(x = Year)) +
      geom_vline(xintercept = 2020, linetype = "dotted", colour = "grey50") +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25, fill = "#2166AC") +
      geom_line(aes(y = mean), linewidth = 1.3, colour = "#2166AC") +
      scale_y_continuous(labels = scales::comma) +
      labs(x = "Year", y = "Total population (all counties)",
           title = "Regional Harbor Seal Population (Central-Northern California)",
           subtitle = paste0("Sum across 6 county groups; ", CI_LABEL)) +
      theme_seal()
    if (save) ggsave(paste0("Output/Plots/", prefix, "_regional_total.jpeg"),
                     p_total, width = 26, height = 14, units = "cm")
  }
  
  # ── (c) County lambda heatmap ──────────────────────────────────────────────
  lam_df <- map_dfr(1:C, function(c) {
    cols <- paste0("lambda_county[", c, ",", 1:(T-1), "]")
    cols_ok <- cols[cols %in% colnames(draws_all)]
    if (!length(cols_ok)) return(NULL)
    mat <- draws_all[, cols_ok, drop = FALSE]
    tibble(County = county_names[c],
           Year   = years[seq_along(cols_ok)],
           lambda = colMeans(mat))
  }) |> mutate(County = factor(County, levels = rev(county_names)))
  
  p_lam_heat <- ggplot(lam_df, aes(x = Year, y = County, fill = lambda)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", lambda)), size = 2.8, colour = "grey10") +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#1A7837",
                         midpoint = 1, limits = c(0.75, 1.25), oob = scales::squish,
                         name = expression(lambda)) +
    geom_vline(xintercept = 2020.5, linetype = "dotted", colour = "grey40") +
    labs(x = "Year", y = NULL,
         title = "Annual Population Growth Rate (λ) by County",
         subtitle = "Posterior mean; red < 1 < green; dotted = 2020 latent year") +
    theme_seal() + theme(panel.grid = element_blank(),
                         axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  if (save) ggsave(paste0("Output/Plots/", prefix, "_county_lambda_heatmap.jpeg"),
                   p_lam_heat, width = 34, height = 14, units = "cm", dpi = 200)
  
  list(trajectories = p_traj, lambda_heatmap = p_lam_heat,
       data = traj_df, lambda_data = lam_df)
}

# ============================================================================
# Q2: BAY vs. COAST MOCI MODIFIER PLOTS
# Addresses: "Does MOCI explain variation the same way everywhere?"
# ============================================================================
create_bay_coast_moci_plots <- function(fit, model_data, save = TRUE, prefix = "Regional") {
  draws <- fit$draws(format = "df")
  
  # ── (a) Prior vs posterior for all four modifiers ─────────────────────────
  mod_params <- list(
    list(nm="delta_moci_mouth_fecund", lbl="Bay Mouth: fecundity",  sd=0.15),
    list(nm="delta_moci_mouth_surv",   lbl="Bay Mouth: survival",   sd=0.15),
    list(nm="delta_moci_south_fecund", lbl="South Bay: fecundity",  sd=0.20),
    list(nm="delta_moci_south_surv",   lbl="South Bay: survival",   sd=0.20)
  )
  mod_df <- map_dfr(mod_params, function(mp) {
    d <- as.numeric(draws[[mp$nm]])
    bind_rows(tibble(draw = d,                        Distribution = "Posterior",
                     Parameter = mp$lbl),
              tibble(draw = rnorm(length(d), 0, mp$sd), Distribution = "Prior",
                     Parameter = mp$lbl))
  }) |> mutate(Parameter = factor(Parameter, levels = sapply(mod_params, `[[`, "lbl")))
  
  p_mod <- ggplot(mod_df, aes(x = draw, fill = Distribution, colour = Distribution)) +
    geom_density(alpha = 0.35, linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values   = c("Posterior" = "#D94801", "Prior" = "grey60")) +
    scale_colour_manual(values = c("Posterior" = "#D94801", "Prior" = "grey40")) +
    facet_wrap(~ Parameter, ncol = 2) +
    labs(x = "Coefficient value (logit scale)", y = "Density",
         title = "MOCI Modifiers by County Type: Prior vs Posterior",
         subtitle = "Expected gradient: South Bay > Bay Mouth > 0 (open coast baseline)") +
    theme_seal() + theme(legend.position = "top")
  if (save) ggsave(paste0("Output/Plots/", prefix, "_bay_coast_modifier.jpeg"),
                   p_mod, width = 28, height = 14, units = "cm")
  
  # ── (b) MOCI response curves: three county types, OND→fecundity ───────────
  xr  <- seq(-2, 2, length.out = 100)
  idx <- sample(seq_len(nrow(draws)), min(500L, nrow(draws)))
  
  beta_ond_f    <- as.numeric(draws$beta_moci_ond_fecund[idx])
  d_mouth_f     <- as.numeric(draws$delta_moci_mouth_fecund[idx])
  d_south_f     <- as.numeric(draws$delta_moci_south_fecund[idx])
  avg_fec       <- as.numeric(draws$avg_fecundity[idx])
  
  fec_df <- do.call(rbind, lapply(xr, function(x) {
    base  <- qlogis(pmax(pmin(avg_fec, 0.999), 0.001))
    fc    <- plogis(base + beta_ond_f * x)
    fm    <- plogis(base + (beta_ond_f + d_mouth_f) * x)
    fs    <- plogis(base + (beta_ond_f + d_south_f) * x)
    data.frame(moci = x,
               mn_coast = mean(fc), lo_coast = quantile(fc, CI_LO), hi_coast = quantile(fc, CI_HI),
               mn_mouth = mean(fm), lo_mouth = quantile(fm, CI_LO), hi_mouth = quantile(fm, CI_HI),
               mn_south = mean(fs), lo_south = quantile(fs, CI_LO), hi_south = quantile(fs, CI_HI))
  }))
  
  p_fec_resp <- ggplot(fec_df, aes(x = moci)) +
    geom_ribbon(aes(ymin=lo_coast,ymax=hi_coast),alpha=0.15,fill="#2166AC")+
    geom_line(aes(y=mn_coast,colour="Open Coast"),linewidth=1.2)+
    geom_ribbon(aes(ymin=lo_mouth,ymax=hi_mouth),alpha=0.15,fill="#F16913")+
    geom_line(aes(y=mn_mouth,colour="Bay Mouth"),linewidth=1.2)+
    geom_ribbon(aes(ymin=lo_south,ymax=hi_south),alpha=0.15,fill="#D94801")+
    geom_line(aes(y=mn_south,colour="South Bay"),linewidth=1.2)+
    geom_vline(xintercept=0,linetype="dashed",colour="grey50")+
    scale_colour_manual(values=COAST_TYPE_COLS,name="County type")+
    labs(x="MOCI OND (SD units)",y="Fecundity (prob. of pupping)",
         title="OND MOCI → Fecundity by County Type",
         subtitle="Expected: South Bay flattest, Open Coast steepest")+
    coord_cartesian(ylim=c(0,1))+theme_seal()
  if (save) ggsave(paste0("Output/Plots/",prefix,"_moci_fecundity_county_type.jpeg"),
                   p_fec_resp,width=22,height=12,units="cm")
  
  # ── (c) JFM → adult survival: three county types ──────────────────────────
  beta_jfm_a  <- as.numeric(draws$beta_moci_jfm_adult[idx])
  d_mouth_s   <- as.numeric(draws$delta_moci_mouth_surv[idx])
  d_south_s   <- as.numeric(draws$delta_moci_south_surv[idx])
  phi_aF      <- as.numeric(draws$phi_adult_F_logit[idx])
  
  surv_df <- do.call(rbind, lapply(xr, function(x) {
    sc <- plogis(phi_aF + beta_jfm_a * x)
    sm <- plogis(phi_aF + (beta_jfm_a + d_mouth_s) * x)
    ss <- plogis(phi_aF + (beta_jfm_a + d_south_s) * x)
    data.frame(moci = x,
               mn_coast=mean(sc),lo_coast=quantile(sc,CI_LO),hi_coast=quantile(sc,CI_HI),
               mn_mouth=mean(sm),lo_mouth=quantile(sm,CI_LO),hi_mouth=quantile(sm,CI_HI),
               mn_south=mean(ss),lo_south=quantile(ss,CI_LO),hi_south=quantile(ss,CI_HI))
  }))
  
  p_surv_resp <- ggplot(surv_df, aes(x = moci)) +
    geom_ribbon(aes(ymin=lo_coast,ymax=hi_coast),alpha=0.15,fill="#2166AC")+
    geom_line(aes(y=mn_coast,colour="Open Coast"),linewidth=1.2)+
    geom_ribbon(aes(ymin=lo_mouth,ymax=hi_mouth),alpha=0.15,fill="#F16913")+
    geom_line(aes(y=mn_mouth,colour="Bay Mouth"),linewidth=1.2)+
    geom_ribbon(aes(ymin=lo_south,ymax=hi_south),alpha=0.15,fill="#D94801")+
    geom_line(aes(y=mn_south,colour="South Bay"),linewidth=1.2)+
    geom_vline(xintercept=0,linetype="dashed",colour="grey50")+
    scale_colour_manual(values=COAST_TYPE_COLS,name="County type")+
    labs(x="MOCI JFM (SD units)",y="Adult female survival",
         title="JFM MOCI → Adult Survival by County Type",
         subtitle=paste0("Expected: South Bay flattest; ",CI_LABEL))+
    coord_cartesian(ylim=c(0.80,1))+theme_seal()
  if (save) ggsave(paste0("Output/Plots/",prefix,"_moci_adult_surv_county_type.jpeg"),
                   p_surv_resp,width=22,height=12,units="cm")
  
  # Print key numbers for all four modifiers
  cat("\n")
  for (mp in mod_params) {
    d    <- as.numeric(draws[[mp$nm]])
    q_lo <- quantile(d, CI_LO); q_hi <- quantile(d, CI_HI)
    cat(sprintf("  %-28s mean=%6.3f  89%%CrI=[%6.3f, %6.3f]  P(>0)=%.3f\n",
                mp$nm, mean(d), q_lo, q_hi, mean(d > 0)))
  }
  
  list(modifier_posteriors = p_mod,
       fecundity_response  = p_fec_resp,
       survival_response   = p_surv_resp)
}

# ============================================================================
# Q3: PROJECTION PLOTS BY COUNTY AND COAST TYPE
# ============================================================================
create_regional_projection_plots <- function(fit, model_data, save = TRUE, prefix = "Regional") {
  years <- model_data$years; T <- length(years)
  C <- model_data$stan_data$C %||% 6
  county_names   <- model_data$county_names   %||% COUNTY_NAMES
  scenario_names <- model_data$scenario_names %||% c("Status Quo","Warm (MOCI +1)","Cool (MOCI -1)")
  T_proj   <- model_data$stan_data$T_proj %||% 10
  N_scen   <- length(scenario_names)
  pyrs     <- (max(years) + 1):(max(years) + T_proj)
  draws_all <- fit$draws(format = "matrix")
  
  # ── County-level projection draws ─────────────────────────────────────────
  # Note: projections stored as N_total_county[c,t] extended forward
  # If not in generated quantities, reconstruct from regional total
  # For now build from county trajectories + simple Leslie matrix forward pass
  # using posterior draws — delegate to in-model generated quantities if available
  
  # Regional total projections (if in generated quantities)
  # These would need to be added in a future Stan model iteration;
  # for now compute from county totals at T (last observed year)
  cat("  (Projection plots require generated quantity N_total_county_proj;")
  cat(" placeholder produced — re-generate after Stan model update.)\n")
  
  # ── Historical county totals for context ──────────────────────────────────
  hist_df <- map_dfr(1:C, function(c) {
    cols <- paste0("N_total_county[", c, ",", 1:T, "]")
    cols_ok <- cols[cols %in% colnames(draws_all)]
    if (!length(cols_ok)) return(NULL)
    mat <- draws_all[, cols_ok, drop = FALSE]
    tibble(County = county_names[c], Year = years[seq_along(cols_ok)],
           mean = colMeans(mat),
           lo = apply(mat, 2, quantile, CI_LO), hi = apply(mat, 2, quantile, CI_HI),
           Period = "Historical")
  }) |> mutate(County = factor(County, levels = county_names))
  
  p_hist <- ggplot(hist_df, aes(x = Year, y = mean, colour = County, fill = County)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1.1) +
    geom_vline(xintercept = 2020, linetype = "dotted", colour = "grey50") +
    scale_colour_manual(values = COUNTY_COLS) +
    scale_fill_manual(values = COUNTY_COLS) +
    scale_y_continuous(labels = scales::comma) +
    facet_wrap(~ County, scales = "free_y", ncol = 3) +
    labs(x = "Year", y = "Total population",
         title = "County Population Trajectories (Historical)",
         subtitle = paste0("Posterior mean + ", CI_LABEL,
                           "; use these as baselines for projection comparison")) +
    theme_seal() + guides(colour = "none", fill = "none")
  if (save) ggsave(paste0("Output/Plots/", prefix, "_county_historical.jpeg"),
                   p_hist, width = 34, height = 22, units = "cm", dpi = 200)
  
  list(historical = p_hist, data = hist_df)
}

# ============================================================================
# Q4: REGIONAL PORTFOLIO ANALYSIS
# Addresses: "Does Bay/Coast asynchrony provide portfolio insurance?"
# ============================================================================
create_regional_portfolio_plots <- function(fit, model_data, save = TRUE, prefix = "Regional") {
  years <- model_data$years; T <- length(years)
  C <- model_data$stan_data$C %||% 6
  county_names <- model_data$county_names %||% COUNTY_NAMES
  # Guard: if county_names length doesn't match C, use generic labels
  if (length(county_names) != C)
    county_names <- paste0("County_", seq_len(C))
  draws_all <- fit$draws(format = "matrix")
  
  # ── County lambda matrix (posterior mean per county per year) ──────────────
  lam_mat <- matrix(NA, C, T - 1, dimnames = list(county_names, head(years, -1)))
  for (c in 1:C) {
    for (t in 1:(T-1)) {
      cn <- paste0("lambda_county[", c, ",", t, "]")
      if (cn %in% colnames(draws_all)) lam_mat[c, t] <- mean(draws_all[, cn])
    }
  }
  
  # ── Cross-county lambda correlation matrix ─────────────────────────────────
  lcor <- cor(t(lam_mat), use = "pairwise.complete.obs")
  cor_df <- expand.grid(County1 = county_names, County2 = county_names) |>
    mutate(r = as.vector(lcor),
           County1 = factor(County1, levels = county_names),
           County2 = factor(County2, levels = rev(county_names)))
  
  p_cor <- ggplot(cor_df, aes(x = County1, y = County2, fill = r)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", r)), size = 3.5) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), name = "r") +
    labs(x = NULL, y = NULL,
         title = "Cross-County Lambda Correlation",
         subtitle = "Negative or low correlation between SF Bay and open-coast counties = portfolio insurance") +
    theme_seal() + theme(panel.grid = element_blank(),
                         axis.text.x = element_text(angle = 30, hjust = 1)) +
    coord_fixed()
  if (save) ggsave(paste0("Output/Plots/", prefix, "_county_lambda_correlation.jpeg"),
                   p_cor, width = 20, height = 18, units = "cm")
  
  # ── Portfolio effect ratio: regional vs individual county ──────────────────
  # Compare CV of regional total lambda vs mean CV of county lambdas
  regional_lam <- colMeans(lam_mat, na.rm = TRUE)   # mean across counties per year
  cv_regional  <- sd(regional_lam, na.rm = TRUE) / mean(regional_lam, na.rm = TRUE)
  cv_counties  <- mean(apply(lam_mat, 1, function(x)
    sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)))
  per <- cv_regional / cv_counties
  
  # Bay vs. coast comparison
  bay_idx   <- which(COUNTY_IS_BAY == 1)
  coast_idx <- which(COUNTY_IS_BAY == 0)
  cor_bay_coast <- if (length(bay_idx) > 0 && length(coast_idx) > 0) {
    bay_lam   <- colMeans(lam_mat[bay_idx,   , drop = FALSE], na.rm = TRUE)
    coast_lam <- colMeans(lam_mat[coast_idx, , drop = FALSE], na.rm = TRUE)
    cor(bay_lam, coast_lam, use = "complete.obs")
  } else NA_real_
  
  port_summary <- tibble(
    Metric = c("Portfolio Effect Ratio (CVregional/CVcounty)",
               "Mean between-county corr (lambda)",
               "Bay vs Coast lambda correlation",
               "Mean regional lambda"),
    Value  = c(round(per, 3),
               round(mean(lcor[lower.tri(lcor)], na.rm = TRUE), 3),
               round(cor_bay_coast, 3),
               round(mean(regional_lam, na.rm = TRUE), 3))
  )
  cat("\n── Regional Portfolio Summary ──\n")
  print(port_summary)
  if (save) write_csv(port_summary, paste0("Output/", prefix, "_portfolio_summary.csv"))
  
  # ── Lambda trajectory by coast type ────────────────────────────────────────
  lam_long <- as_tibble(t(lam_mat), rownames = "Year_chr") |>
    mutate(Year = as.integer(Year_chr)) |>
    pivot_longer(-c(Year, Year_chr), names_to = "County", values_to = "lambda") |>
    mutate(CoastType = case_when(County == "Bay Mouth"  ~ "Bay Mouth",
                                 County == "South Bay"  ~ "South Bay",
                                 TRUE                   ~ "Open Coast"),
           County    = factor(County, levels = county_names))
  
  # Mean lambda by coast type
  lam_ct <- lam_long |>
    group_by(Year, CoastType) |>
    summarise(mean_lam = mean(lambda, na.rm = TRUE), .groups = "drop")
  
  p_lam_ct <- ggplot(lam_ct, aes(x = Year, y = mean_lam, colour = CoastType)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_line(linewidth = 1.3) +
    geom_point(size = 2) +
    scale_colour_manual(values = COAST_TYPE_COLS, name = NULL) +
    labs(x = "Year", y = expression(bar(lambda)~"(mean across counties)"),
         title = "Annual Lambda by Coast Type: Open Coast vs SF Bay",
         subtitle = "Divergence in poor upwelling years supports portfolio insurance hypothesis") +
    theme_seal()
  if (save) ggsave(paste0("Output/Plots/", prefix, "_lambda_coast_type.jpeg"),
                   p_lam_ct, width = 24, height = 12, units = "cm")
  
  list(correlation_plot = p_cor, lambda_coast_type = p_lam_ct,
       portfolio_summary = port_summary,
       lambda_matrix = lam_mat, cor_matrix = lcor)
}

# ============================================================================
# FOREST PLOT (MOCI + BAY MODIFIER)
# ============================================================================
create_regional_forest_plot <- function(fit, save = TRUE, prefix = "Regional") {
  params <- tribble(
    ~variable,                  ~label,                              ~group,
    "beta_moci_ond_fecund",     "MOCI OND → fecundity (baseline)",  "MOCI: open coast",
    "beta_moci_ond_pup",        "MOCI OND → pup survival",          "MOCI: open coast",
    "beta_moci_jfm_pup",        "MOCI JFM → pup survival",          "MOCI: open coast",
    "beta_moci_amj_pup",        "MOCI AMJ → pup survival",          "MOCI: open coast",
    "beta_moci_jfm_juv",        "MOCI JFM → juvenile survival",     "MOCI: open coast",
    "beta_moci_jfm_adult",      "MOCI JFM → adult survival",        "MOCI: open coast",
    "beta_moci_amj_molt",       "MOCI AMJ → molt detection",        "MOCI: open coast",
    "delta_moci_mouth_fecund",  "Bay Mouth: fecundity modifier",   "Bay modifier",
    "delta_moci_mouth_surv",    "Bay Mouth: survival modifier",    "Bay modifier",
    "delta_moci_south_fecund",  "South Bay: fecundity modifier",   "Bay modifier",
    "delta_moci_south_surv",    "South Bay: survival modifier",    "Bay modifier"
  )
  grp_cols <- c("MOCI: open coast" = "#2166AC", "Bay modifier" = "#D94801")
  
  draws_all <- tryCatch(fit$draws(format = "draws_df"),
                        error = function(e) fit$draws(format = "matrix"))
  post <- map_dfr(seq_len(nrow(params)), function(i) {
    v <- tryCatch(as.numeric(posterior::extract_variable(draws_all, params$variable[i])),
                  error = function(e) numeric(0))
    if (!length(v)) return(NULL)
    tibble(variable = params$variable[i], mean = mean(v),
           lo89 = quantile(v, CI_LO),  hi89 = quantile(v, CI_HI),
           lo50 = quantile(v, 0.25),   hi50 = quantile(v, 0.75))
  })
  df <- post |> left_join(params, by = "variable") |>
    mutate(label = factor(label, levels = rev(params$label)),
           group = factor(group, levels = names(grp_cols)),
           sig   = (lo89 > 0) | (hi89 < 0))
  
  p <- ggplot(df, aes(y = label, colour = group)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_linerange(aes(xmin = lo89, xmax = hi89), linewidth = 0.8, alpha = 0.7) +
    geom_linerange(aes(xmin = lo50, xmax = hi50), linewidth = 2.5, alpha = 0.9) +
    geom_point(aes(x = mean, shape = sig), size = 3.5) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), name = NULL,
                       labels = c("FALSE" = "CrI spans zero",
                                  "TRUE"  = "Excludes zero (89% CrI)")) +
    scale_colour_manual(values = grp_cols, name = "Parameter group") +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Coefficient (logit scale)", y = NULL,
         title = "Regional IPM: MOCI Effect Forest Plot",
         subtitle = paste0("Open-coast baseline (blue) + Bay vs. Coast modifier (orange); thick=50% CrI, thin=89% CrI")) +
    theme_seal(base_size = 13) +
    theme(panel.grid.major.y = element_blank(),
          strip.text.y = element_text(angle = 0, face = "bold"))
  if (save) ggsave(paste0("Output/Plots/", prefix, "_forest_plot.jpeg"),
                   p, width = 28, height = 20, units = "cm", dpi = 200)
  list(plot = p, data = df)
}

# ============================================================================
# SITE AVAILABILITY PARAMETER PLOTS
# ============================================================================
create_site_availability_plots <- function(fit, model_data, save = TRUE, prefix = "Regional") {
  site_meta <- model_data$site_meta
  if (is.null(site_meta)) {
    cat("  site_meta not available — skipping availability plots.\n")
    return(NULL)
  }
  S <- nrow(site_meta)
  
  alpha_summ <- regional_fit_summary(
    fit, c(paste0("log_alpha_breed[", 1:S, "]"),
           paste0("log_alpha_pup[",   1:S, "]"),
           paste0("log_alpha_molt[",  1:S, "]")))
  
  alpha_df <- alpha_summ |>
    mutate(
      param_type = case_when(
        str_detect(variable, "breed") ~ "Breeding (adult count)",
        str_detect(variable, "pup")   ~ "Pup count",
        str_detect(variable, "molt")  ~ "Molt count"),
      s_idx = as.integer(str_extract(variable, "\\d+")),
    ) |>
    left_join(site_meta |> select(site_id, site_name, county), by = c("s_idx" = "site_id")) |>
    filter(!is.na(site_name)) |>
    mutate(site_name = factor(site_name, levels = rev(site_meta$site_name)),
           alpha = exp(mean))   # back-transform for interpretability
  
  p_alpha <- ggplot(alpha_df, aes(y = site_name, colour = county)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_linerange(aes(xmin = q_lo, xmax = q_hi), linewidth = 0.8, alpha = 0.7) +
    geom_point(aes(x = mean), size = 2.5) +
    scale_colour_manual(values = COUNTY_COLS, name = "County") +
    facet_wrap(~ param_type, ncol = 3) +
    labs(x = "log(alpha): site availability × detection (log scale)", y = NULL,
         title = "Site Availability Parameters by Survey Type",
         subtitle = "Higher = more of county population visible at this site; bars = 89% CrI") +
    theme_seal(base_size = 12) +
    theme(axis.text.y = element_text(size = 8))
  if (save) ggsave(paste0("Output/Plots/", prefix, "_site_availability.jpeg"),
                   p_alpha, width = 34, height = 24, units = "cm", dpi = 200)
  list(plot = p_alpha, data = alpha_df)
}

# ============================================================================
# SUMMARY TABLE
# ============================================================================
create_regional_summary_table <- function(fit, save = TRUE, prefix = "Regional") {
  params <- c(
    "phi_pup_logit", "phi_juv_base", "phi_adult_F_logit", "phi_adult_F_base",
    "delta_adult", "fecund_primip", "fecund_mature", "prop_female", "avg_fecundity",
    "rho_pup",
    "beta_moci_ond_fecund", "beta_moci_ond_pup", "beta_moci_amj_pup",
    "beta_moci_jfm_pup", "beta_moci_jfm_juv", "beta_moci_jfm_adult",
    "beta_moci_amj_molt",
    "delta_moci_mouth_fecund", "delta_moci_mouth_surv",
    "delta_moci_south_fecund", "delta_moci_south_surv",
    "sigma_county", "sigma_site", "detect_breed_logit", "detect_molt_logit",
    "sigma_process", "sigma_obs_adult", "sigma_obs_pup", "sigma_obs_molt"
  )
  s <- regional_fit_summary(fit, params) |>
    mutate(
      Estimate = sprintf("%.3f (%.3f, %.3f)", mean, q_lo, q_hi),
      Category = case_when(
        str_detect(variable, "phi|delta_adult")                            ~ "Survival",
        str_detect(variable, "delta_moci_mouth|delta_moci_south")         ~ "MOCI: Bay modifier",
        str_detect(variable, "fecund|prop|avg")                           ~ "Reproduction",
        variable == "rho_pup"                                             ~ "Molt attendance",
        str_detect(variable, "beta_moci") & !str_detect(variable, "bay") ~ "MOCI: open coast",
        str_detect(variable, "detect|sigma_site")                         ~ "Detection",
        str_detect(variable, "sigma_county")                              ~ "Random effects",
        str_detect(variable, "sigma")                                     ~ "Error terms",
        TRUE                                                              ~ "Other"
      )
    ) |>
    select(Category, variable, Estimate, rhat, ess_bulk)
  
  cat("\n=== PARAMETER SUMMARY — REGIONAL IPM ===\n")
  print(s, n = nrow(s))
  if (save) write_csv(s, paste0("Output/", prefix, "_parameter_summary.csv"))
  s
}

cat("\n18_regional_plots.R loaded — regional plotting functions ready.\n")
cat("Usage:\n")
cat("  out <- load_regional_results(\"Regional_real\")\n")
cat("  run_all_regional_plots(out$fit, out$model_data, prefix=\"Regional_real\")\n")
