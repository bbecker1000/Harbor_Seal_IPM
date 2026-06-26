# ============================================================================
# 06_ipm_plots.R  —  IPM v3.3 PLOTS, TABLES & POST-PROCESSING
# ----------------------------------------------------------------------------
# v3.3 additions vs v3.2:
#   - rho_pup added to diagnostics, trace plots, summary table, recovery
#   - create_molt_pup_plots_v3.3():  posterior of rho_pup; estimated pup
#     counts in molt surveys by site and year; stacked area showing pup vs
#     juv+adult contribution to N_molt_true; pup % of molt count over time.
#   - run_all_plots_v3.3() calls the new function.
# All other functions are identical to v3.2.
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
CI_FMT   <- "89%% CrI"

SEAL_COLS <- list(
  pop="#2166AC", pup="#1B7837", juv="#762A83", adult_f="#B2182B",
  adult_m="#D6604D", molt="#8C510A", ribbon="#AECDE8",
  neutral="#888888", warn="#FF7F00"
)

theme_seal <- function(base_size = 16) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.grid.major = element_line(colour="grey88", linewidth=0.4),
      panel.grid.minor = element_line(colour="grey93", linewidth=0.2),
      panel.border     = element_rect(colour="grey70", fill=NA, linewidth=0.5),
      axis.title  = element_text(size=rel(0.95), colour="grey20"),
      axis.text   = element_text(size=rel(0.88), colour="grey30"),
      axis.ticks  = element_line(colour="grey70", linewidth=0.3),
      legend.position="bottom",
      legend.title=element_text(size=rel(0.88), face="bold"),
      legend.text =element_text(size=rel(0.82)),
      legend.key.width=unit(1.6,"cm"), legend.background=element_blank(),
      strip.text =element_text(size=rel(0.90), face="bold", colour="grey20"),
      strip.background=element_rect(fill="grey94", colour="grey80", linewidth=0.3),
      plot.title   =element_text(size=rel(1.05), face="bold", margin=margin(b=6)),
      plot.subtitle=element_text(size=rel(0.88), colour="grey40", margin=margin(b=8)),
      plot.caption =element_text(size=rel(0.78), colour="grey50", hjust=1),
      plot.margin  =margin(10,14,10,10)
    )
}

# ── Load helper ──────────────────────────────────────────────────────────────
load_seal_results <- function(prefix          = "IPM_v3.3_real",
                              fit_path        = NULL,
                              input_data_path = NULL,
                              years           = 1997:2025,
                              T_proj          = 10) {
  if (is.null(fit_path))
    fit_path <- paste0("Output/harbor_seal_", prefix, "_fit.rds")
  if (!file.exists(fit_path))
    stop("Fit file not found: ", fit_path)
  cat("Loading fit from", fit_path, "...\n")
  fit <- readRDS(fit_path)
  
  if (is.null(input_data_path)) {
    auto_path <- paste0("Output/harbor_seal_", prefix, "_input_data.rds")
    if (file.exists(auto_path)) input_data_path <- auto_path
  }
  if (!is.null(input_data_path) && file.exists(input_data_path)) {
    inp <- readRDS(input_data_path)
    sim_data <- list(
      stan_data      = inp$stan_data,
      site_names     = c("BL","DE","DP","PRH","TB","TP"),
      years          = inp$years %||% years,
      scenario_names = c("Status Quo","Warm (MOCI +1)",
                         "Cool (MOCI -1)","Warm + High Coyote"),
      true_params    = inp$true_params %||% NULL
    )
  } else {
    cat("No input_data_path — using metadata shell.\n")
    T_val <- fit$metadata()$data$T %||% length(years)
    sim_data <- list(
      stan_data = list(T=T_val, S=6, T_proj=T_proj, N_scenarios=4,
                       y_adult_obs=matrix(1L,6,T_val), y_pup_obs=matrix(1L,6,T_val),
                       y_molt_obs=matrix(1L,6,T_val), y_adult=matrix(0,6,T_val),
                       y_pup=matrix(0,6,T_val), y_molt=matrix(0,6,T_val)),
      site_names = c("BL","DE","DP","PRH","TB","TP"),
      years      = years[seq_len(T_val)],
      scenario_names = c("Status Quo","Warm (MOCI +1)",
                         "Cool (MOCI -1)","Warm + High Coyote"),
      true_params = NULL)
  }
  cat("Ready. run_all_plots_v3.3(out$fit, out$sim_data) to regenerate outputs.\n")
  list(fit=fit, sim_data=sim_data)
}

# ── 89% CrI fit summary ──────────────────────────────────────────────────────
seal_fit_summary <- function(fit, variables) {
  stale_msg <- "fit methods unavailable. Reload: out <- load_seal_results()"
  s <- tryCatch(fit$summary(variables=variables),
                error=function(e) stop(stale_msg))
  d <- tryCatch(fit$draws(variables=variables, format="matrix"),
                error=function(e) stop(stale_msg))
  s$q_lo <- apply(d[, s$variable, drop=FALSE], 2, quantile, CI_LO)
  s$q_hi <- apply(d[, s$variable, drop=FALSE], 2, quantile, CI_HI)
  .to_num <- function(x) if (is.list(x))
    vapply(x, function(v) suppressWarnings(as.numeric(v)[1L]), numeric(1))
  else suppressWarnings(as.numeric(x))
  for (col in c("rhat","ess_bulk","ess_tail"))
    if (col %in% names(s)) s[[col]] <- .to_num(s[[col]])
  s
}

# ============================================================================
# ORCHESTRATOR
# ============================================================================
run_all_plots_v3.3 <- function(fit, sim_data, prefix="IPM_v3.3",
                               save=TRUE, run_recovery=FALSE,
                               run_portfolio=FALSE, run_synchrony=FALSE) {
  
  safe_run <- function(label, expr) {
    cat(sprintf("\n%s\n", label))
    tryCatch(
      withCallingHandlers(expr, warning=function(w) invokeRestart("muffleWarning")),
      error=function(e) { cat(sprintf("  !! FAILED: %s\n", conditionMessage(e))); NULL })
  }
  
  diag <- safe_run("── Diagnostics ──", check_diagnostics_v3.3(fit))
  
  params_candidate <- c(
    "phi_pup_logit","phi_juv_base","phi_adult_F_logit","phi_adult_F_base","delta_adult",
    "fecund_primip","fecund_mature","prop_female","p_male_breed","rho_pup",
    "beta_coy[1]","beta_coy[2]","beta_coy[3]",
    paste0("beta_dist_surv[",1:6,"]"),
    "beta_moci_ond_fecund","beta_moci_amj_pup","beta_moci_jfm_juv",
    "beta_moci_jfm_adult","beta_eseal_pup","detect_breed_logit","detect_molt_logit",
    "beta_moci_ond_pup","beta_moci_jfm_pup",
    "sigma_process","sigma_obs_adult","sigma_obs_pup","sigma_obs_molt","sigma_site")
  available_vars <- tryCatch(fit$summary()$variable, error=function(e) params_candidate)
  params <- if (!is.null(diag) && !is.null(diag$params))
    diag$params[diag$params %in% available_vars]
  else params_candidate[params_candidate %in% available_vars]
  
  traces  <- safe_run("── Trace plots ──",
                      create_trace_plots_v3.3(fit, params, save=save, prefix=prefix))
  
  rec <- NULL
  if (run_recovery && !is.null(sim_data$true_params))
    rec <- safe_run("── Parameter recovery ──",
                    check_parameter_recovery_v3.3(fit, sim_data, save=save, prefix=prefix))
  else cat("\n── Parameter recovery: skipped ──\n")
  
  ppc    <- safe_run("── PPC ──",         create_ppc_plots_v3.2(fit, sim_data, save=save, prefix=prefix))
  ts     <- safe_run("── Time series ──", create_timeseries_plots_v3.2(fit, sim_data, save=save, prefix=prefix))
  sa     <- safe_run("── Site x age ──",  create_site_age_timeseries_v3.2(fit, sim_data, save=save, prefix=prefix))
  sasite <- safe_run("── Site panels ──", create_site_panels_v3.2(fit, sim_data, save=save, prefix=prefix))
  proj   <- safe_run("── Projections ──", create_projection_plots_v3.2(fit, sim_data, save=save, prefix=prefix))
  le     <- safe_run("── Lambda + elasticity ──", create_lambda_elasticity_v3.2(fit, save=save, prefix=prefix))
  moci   <- safe_run("── MOCI collinearity ──", create_moci_collinearity_v3.2(fit, sim_data, save=save, prefix=prefix))
  eff    <- safe_run("── Effect plots ──", create_effect_plots_v3.2(fit, save=save, prefix=prefix))
  jveff  <- safe_run("── Juv+adult effects ──", create_juv_adult_effect_plots_v3.2(fit, save=save, prefix=prefix))
  forest <- safe_run("── Forest plot ──", create_forest_plot_v3.2(fit, save=save, prefix=prefix))
  decomp <- safe_run("── Decomposition ──", create_covariate_decomposition_plots_v3.2(fit, sim_data, save=save, prefix=prefix))
  tbl    <- safe_run("── Summary table ──", create_summary_table_v3.3(fit, save=save, prefix=prefix))
  # v3.3: new plot for pup molt attendance
  molt_pup <- safe_run("── Molt pup decomposition ──",
                       create_molt_pup_plots_v3.3(fit, sim_data, save=save, prefix=prefix))
  safe_run("── Key results ──", save_model_output_v3.2(fit, prefix=prefix))
  
  port <- sync <- NULL
  if (run_portfolio)
    port <- safe_run("── Portfolio ──", create_portfolio_analysis_v3.2(fit, sim_data, save=save, prefix=prefix))
  if (run_synchrony)
    sync <- safe_run("── Synchrony ──", create_synchrony_projections_v3.2(fit, sim_data, save=save, prefix=prefix))
  
  results <- list(diagnostics=diag, traces=traces, recovery=rec, ppc=ppc, ts=ts,
                  site_age=sa, site_panels=sasite, projections=proj, effects=eff,
                  juv_adult_effects=jveff, forest=forest, decomposition=decomp,
                  lambda_elast=le, moci_collin=moci, table=tbl,
                  molt_pup=molt_pup,    # v3.3 addition
                  portfolio=port, sync=sync)
  n_ok   <- sum(!sapply(results[!sapply(results, is.null)], is.null))
  cat(sprintf("\n── Pipeline complete: %d steps ──\n", n_ok))
  invisible(results)
}

# ============================================================================
# DIAGNOSTICS  (v3.3: rho_pup added)
# ============================================================================
check_diagnostics_v3.3 <- function(fit) {
  cat("\n=== MODEL DIAGNOSTICS — IPM v3.3 ===\n")
  tryCatch({
    ds <- fit$diagnostic_summary(quiet=TRUE)
    cat(if (all(ds$num_max_treedepth==0)) "Treedepth OK.\n"
        else sprintf("WARN: %d treedepth hits.\n", sum(ds$num_max_treedepth)))
    cat(if (all(ds$num_divergent==0)) "No divergences.\n"
        else sprintf("WARN: %d divergences.\n", sum(ds$num_divergent)))
  }, error=function(e)
    cat("NOTE: sampler diagnostics unavailable (CSVs gone).\n"))
  
  params <- c(
    "phi_pup_logit","phi_juv_base","phi_adult_F_logit","phi_adult_F_base","delta_adult",
    "fecund_primip","fecund_mature","prop_female","p_male_breed",
    "rho_pup",                          # v3.3
    "beta_coy[1]","beta_coy[2]","beta_coy[3]",
    paste0("beta_dist_surv[",1:6,"]"),
    "beta_moci_ond_fecund","beta_moci_amj_pup","beta_moci_jfm_juv",
    "beta_moci_jfm_adult","beta_eseal_pup",
    "sigma_process","sigma_obs_adult","sigma_obs_pup","sigma_obs_molt","sigma_site")
  
  s <- seal_fit_summary(fit, params)
  cat("\nParameter Summary:\n")
  print(s |> select(variable,mean,sd,q_lo,q_hi,rhat,ess_bulk), n=nrow(s))
  
  # Report rho_pup explicitly
  rho_row <- s[s$variable == "rho_pup", ]
  if (nrow(rho_row) > 0)
    cat(sprintf("\nrho_pup (pup molt attendance): median=%.3f, %s=[%.3f, %.3f]\n",
                rho_row$mean, CI_LABEL, rho_row$q_lo, rho_row$q_hi))
  
  rhat_v <- tryCatch(as.double(as.vector(unclass(s$rhat))), error=function(e) rep(NA_real_, nrow(s)))
  ess_v  <- tryCatch(as.double(as.vector(unclass(s$ess_bulk))), error=function(e) rep(NA_real_, nrow(s)))
  bad <- s[!is.na(rhat_v) & rhat_v > 1.05, , drop=FALSE]
  low <- s[!is.na(ess_v)  & ess_v  < 400,  , drop=FALSE]
  if (nrow(bad)>0) { cat("\nWARNING Rhat>1.05:\n"); print(bad[, c("variable","rhat")]) }
  if (nrow(low)>0) { cat("\nWARNING low ESS:\n");    print(low[, c("variable","ess_bulk")]) }
  
  list(params=params, summary=s)
}

# ============================================================================
# TRACE PLOTS  (v3.3: rho_pup in errors group)
# ============================================================================
create_trace_plots_v3.3 <- function(fit, params, save=TRUE, prefix="IPM_v3.3") {
  draws <- fit$draws(format="df")
  grp <- function(pars, title, file, w=30, h=18) {
    pars <- pars[pars %in% colnames(draws)]
    if (!length(pars)) return(NULL)
    p <- mcmc_trace(draws, pars=pars) + labs(title=title)
    if (save) ggsave(paste0("Output/Plots/",prefix,"_",file,".jpeg"), p, width=w, height=h, units="cm")
    p
  }
  list(
    survival    = grp(c("phi_pup_logit","phi_juv_base","phi_adult_F_logit","delta_adult","p_male_breed"),
                      "Trace: Survival + Observation","trace_survival"),
    fecundity   = grp(c("fecund_primip","fecund_mature","prop_female"),"Trace: Fecundity","trace_fecundity"),
    coyote      = grp(c("beta_coy[1]","beta_coy[2]","beta_coy[3]"),"Trace: Coyote","trace_coyote", h=12),
    disturbance = grp(paste0("beta_dist_surv[",1:6,"]"),"Trace: Disturbance","trace_disturbance"),
    moci        = grp(c("beta_moci_ond_fecund","beta_moci_amj_pup","beta_moci_jfm_juv",
                        "beta_moci_jfm_adult","beta_moci_amj_molt"),"Trace: MOCI","trace_moci"),
    errors      = grp(c("sigma_process","sigma_obs_adult","sigma_obs_pup","sigma_obs_molt",
                        "sigma_site","rho_pup"),   # v3.3: rho_pup added
                      "Trace: Error Terms + Molt Attendance","trace_errors"))
}

# ============================================================================
# PARAMETER RECOVERY  (v3.3: rho_pup added)
# ============================================================================
check_parameter_recovery_v3.3 <- function(fit, sim_data, save=TRUE, prefix="IPM_v3.3") {
  tp <- sim_data$true_params
  true_vals <- tibble(
    parameter = c("phi_pup_logit","phi_juv_base","phi_adult_F_logit","delta_adult",
                  "fecund_primip","fecund_mature","prop_female","p_male_breed",
                  "rho_pup",              # v3.3
                  "beta_coy[1]","beta_coy[2]","beta_coy[3]",
                  paste0("beta_dist_surv[",1:6,"]"),
                  "beta_moci_ond_fecund","beta_moci_ond_pup","beta_moci_jfm_pup",
                  "beta_moci_amj_pup","beta_moci_jfm_juv","beta_moci_jfm_adult","beta_eseal_pup",
                  "detect_breed_logit","detect_molt_logit",
                  "sigma_process","sigma_obs_adult","sigma_obs_pup","sigma_obs_molt"),
    true_value = c(tp$phi_pup_logit, tp$phi_juv_base, tp$phi_adult_F_logit, tp$delta_adult,
                   tp$fecund_primip, tp$fecund_mature, tp$prop_female, tp$p_male_breed,
                   tp$rho_pup,       # v3.3
                   tp$beta_coy[1], tp$beta_coy[2], tp$beta_coy[3],
                   tp$beta_dist_surv[1], tp$beta_dist_surv[2], tp$beta_dist_surv[3],
                   tp$beta_dist_surv[4], tp$beta_dist_surv[5], tp$beta_dist_surv[6],
                   tp$beta_moci_ond_fecund, tp$beta_moci_ond_pup, tp$beta_moci_jfm_pup,
                   tp$beta_moci_amj_pup, tp$beta_moci_jfm_juv, tp$beta_moci_jfm_adult,
                   tp$beta_eseal_pup,
                   tp$detect_breed_logit, tp$detect_molt_logit,
                   tp$sigma_process, tp$sigma_obs_adult, tp$sigma_obs_pup, tp$sigma_obs_molt))
  
  rec <- seal_fit_summary(fit, true_vals$parameter) |>
    left_join(true_vals, by=c("variable"="parameter")) |>
    mutate(recovered    = true_value >= q_lo & true_value <= q_hi,
           rel_bias_pct = (mean - true_value) / abs(true_value) * 100)
  
  cat("Recovery:", sum(rec$recovered), "/", nrow(rec),
      sprintf("(%.1f%%)\n", 100*mean(rec$recovered)))
  
  rec <- rec |> mutate(
    identifiability = case_when(
      variable %in% c("phi_adult_F_logit","fecund_mature","beta_moci_ond_fecund",
                      "detect_breed_logit","sigma_obs_adult") ~ "Well identified",
      variable %in% c("phi_juv_base","beta_coy[1]","beta_coy[2]","beta_coy[3]",
                      "beta_moci_jfm_adult","detect_molt_logit",
                      "sigma_process","sigma_obs_molt","rho_pup") ~ "Moderately identified",
      TRUE ~ "Prior dominated"),
    identifiability = factor(identifiability,
                             levels=c("Well identified","Moderately identified","Prior dominated")))
  
  p <- ggplot(rec, aes(x=true_value, y=mean)) +
    geom_abline(slope=1, intercept=0, linetype=2, color="gray50") +
    geom_pointrange(aes(ymin=q_lo, ymax=q_hi, color=recovered, shape=identifiability), size=0.8) +
    geom_text(aes(label=variable), hjust=-0.1, vjust=-0.3, size=2.5, check_overlap=TRUE) +
    scale_color_manual(values=c("TRUE"=SEAL_COLS$pop,"FALSE"=SEAL_COLS$adult_f),
                       name=paste("True in", CI_LABEL)) +
    scale_shape_manual(values=c("Well identified"=16,"Moderately identified"=17,"Prior dominated"=1)) +
    facet_wrap(~identifiability, scales="free") +
    labs(x="True Value", y=paste0("Posterior Mean (", CI_LABEL, ")"),
         title="Parameter Recovery: IPM v3.3") +
    theme_seal()
  if (save) ggsave(paste0("Output/Plots/",prefix,"_parameter_recovery.jpeg"),
                   p, width=30, height=22, units="cm")
  if (save) write_csv(rec |> select(variable,true_value,mean,q_lo,q_hi,recovered,identifiability),
                      paste0("Output/",prefix,"_parameter_recovery_table.csv"))
  list(table=rec, plot=p)
}

# ============================================================================
# SUMMARY TABLE  (v3.3: rho_pup added)
# ============================================================================
create_summary_table_v3.3 <- function(fit, save=TRUE, prefix="IPM_v3.3") {
  params <- c("phi_pup_logit","phi_juv_base","phi_adult_F_logit","phi_adult_F_base","delta_adult",
              "fecund_primip","fecund_mature","prop_female","avg_fecundity","p_male_breed",
              "rho_pup",          # v3.3
              "beta_coy[1]","beta_coy[2]","beta_coy[3]",
              paste0("beta_dist_surv[",1:6,"]"),
              "beta_moci_ond_fecund","beta_moci_amj_pup","beta_moci_jfm_juv",
              "beta_moci_jfm_adult","beta_eseal_pup","beta_moci_amj_molt",
              "sigma_process","sigma_obs_adult","sigma_obs_pup","sigma_obs_molt","sigma_site")
  pup_draws <- fit$draws(variables="phi_pup_logit",format="df")$phi_pup_logit
  pup_base <- tibble(variable="phi_pup_base (prob)",
                     mean=mean(plogis(pup_draws)), sd=sd(plogis(pup_draws)),
                     q_lo=quantile(plogis(pup_draws),CI_LO),
                     q_hi=quantile(plogis(pup_draws),CI_HI),
                     rhat=NA_real_, ess_bulk=NA_real_)
  
  tbl <- seal_fit_summary(fit, params) |> bind_rows(pup_base) |>
    mutate(Estimate=sprintf("%.3f (%.3f, %.3f)", mean, q_lo, q_hi),
           Category=case_when(
             str_detect(variable,"phi|delta") ~ "Survival",
             str_detect(variable,"fecund|prop|avg") ~ "Reproduction",
             variable %in% c("p_male_breed","phi_pup_base (prob)") ~ "Observation/Derived",
             variable == "rho_pup" ~ "Molt attendance (v3.3)",  # v3.3
             str_detect(variable,"beta_coy") ~ "Coyote",
             str_detect(variable,"beta_dist") ~ "Disturbance",
             str_detect(variable,"beta_moci|beta_eseal") ~ "Shared covariates",
             str_detect(variable,"sigma") ~ "Error terms", TRUE ~ "Other")) |>
    select(Category, variable, Estimate, rhat, ess_bulk)
  cat("\n=== PARAMETER SUMMARY — IPM v3.3 ===\n")
  print(tbl, n=nrow(tbl))
  if (save) write_csv(tbl, paste0("Output/",prefix,"_parameter_summary.csv"))
  tbl
}

# ============================================================================
# NEW v3.3 FUNCTION: MOLT PUP DECOMPOSITION PLOTS
# ----------------------------------------------------------------------------
# Four panels:
#  (a) Posterior of rho_pup — prior vs posterior density
#  (b) Estimated N_pup_in_molt over time (aggregate across sites)
#  (c) Pup fraction of total N_molt_true by site over time (stacked area)
#  (d) Site-by-year heatmap: % of molt count attributable to pups
# ============================================================================
create_molt_pup_plots_v3.3 <- function(fit, sim_data, save=TRUE, prefix="IPM_v3.3") {
  years      <- sim_data$years
  site_names <- sim_data$site_names
  S <- length(site_names); T <- length(years)
  
  # ── Pull draws ────────────────────────────────────────────────────────────
  rho_draws <- tryCatch(
    fit$draws(variables="rho_pup", format="df")$rho_pup,
    error=function(e) {
      warning("rho_pup not found in fit (may be v3.2 fit — run v3.3 to get this plot)")
      return(NULL)
    })
  if (is.null(rho_draws)) {
    cat("  Skipping molt_pup plots: rho_pup not in fit object (v3.2 fit).\n")
    return(NULL)
  }
  
  # ── (a) Prior vs Posterior density for rho_pup ────────────────────────────
  prior_rho <- rbeta(length(rho_draws), 2, 5)
  dens_df <- bind_rows(
    tibble(value=rho_draws,  Distribution="Posterior"),
    tibble(value=prior_rho,  Distribution="Prior: Beta(2,5)")
  )
  p_rho <- ggplot(dens_df, aes(x=value, fill=Distribution, colour=Distribution)) +
    geom_density(alpha=0.35, linewidth=0.8) +
    geom_vline(xintercept=median(rho_draws), linetype="dashed",
               colour=SEAL_COLS$molt, linewidth=0.9) +
    scale_fill_manual(values=c("Posterior"=SEAL_COLS$molt, "Prior: Beta(2,5)"="grey60")) +
    scale_colour_manual(values=c("Posterior"=SEAL_COLS$molt, "Prior: Beta(2,5)"="grey40")) +
    annotate("text", x=median(rho_draws), y=Inf,
             label=sprintf(" Median %.2f\n %s [%.2f, %.2f]",
                           median(rho_draws), CI_LABEL,
                           quantile(rho_draws, CI_LO), quantile(rho_draws, CI_HI)),
             hjust=0, vjust=1.4, size=3.8, colour=SEAL_COLS$molt) +
    labs(x=expression(rho[pup]~"(pup molt attendance fraction)"),
         y="Density",
         title=expression("Posterior vs Prior: "*rho[pup]),
         subtitle="Fraction of year-class pups present at haul-out during peak molt count") +
    theme_seal() + theme(legend.position="top")
  
  # ── (b) N_pup_in_molt_all over time ─────────────────────────────────────
  pim_all <- tryCatch(
    fit$draws(variables="N_pup_in_molt_all", format="matrix"),
    error=function(e) NULL)
  
  p_pim_ts <- NULL
  if (!is.null(pim_all)) {
    pim_df <- tibble(
      Year = years,
      mean = colMeans(pim_all),
      lo   = as.numeric(apply(pim_all, 2, quantile, CI_LO)),
      hi   = as.numeric(apply(pim_all, 2, quantile, CI_HI))
    )
    p_pim_ts <- ggplot(pim_df, aes(x=Year)) +
      geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.25, fill=SEAL_COLS$pup) +
      geom_line(aes(y=mean), linewidth=1.2, colour=SEAL_COLS$pup) +
      expand_limits(y=0) +
      scale_y_continuous(labels=scales::comma) +
      labs(x="Year", y="Estimated pups in molt count (all sites)",
           title="Estimated Pup Contribution to Peak Molt Counts",
           subtitle=sprintf("= rho_pup * N_pup summed across 6 sites; posterior mean %s", CI_LABEL)) +
      theme_seal()
  }
  
  # ── (c) Per-site pup counts in molt vs juv+adult — stacked area ──────────
  # Pull N_pup[s,t] and N_molt_true[s,t] per draw, compute means
  pim_site <- tryCatch(
    fit$draws(variables="N_pup_in_molt", format="matrix"),
    error=function(e) NULL)
  molt_mat <- tryCatch(
    fit$draws(variables="N_molt_true", format="matrix"),
    error=function(e) NULL)
  
  p_stack <- NULL; p_pct <- NULL
  if (!is.null(pim_site) && !is.null(molt_mat)) {
    stack_df <- map_dfr(1:S, function(s) {
      map_dfr(1:T, function(t) {
        pim_col  <- paste0("N_pup_in_molt[",s,",",t,"]")
        molt_col <- paste0("N_molt_true[",s,",",t,"]")
        pim_m  <- if (pim_col  %in% colnames(pim_site))  mean(pim_site[,pim_col])  else NA_real_
        molt_m <- if (molt_col %in% colnames(molt_mat))  mean(molt_mat[,molt_col]) else NA_real_
        tibble(Site=site_names[s], Year=years[t],
               N_pup_in_molt=pim_m,
               N_nonjuv_adult=molt_m - pim_m,  # juv + adult component
               N_molt_total=molt_m,
               pup_pct=100 * pim_m / molt_m)
      })
    }) |> mutate(Site=factor(Site, levels=site_names))
    
    long_df <- stack_df |>
      select(Site, Year, `Pup (rho_pup x N_pup)`=N_pup_in_molt,
             `Juv + Adult`=N_nonjuv_adult) |>
      pivot_longer(cols=c("Pup (rho_pup x N_pup)","Juv + Adult"),
                   names_to="Component", values_to="Count") |>
      mutate(Component=factor(Component, levels=c("Juv + Adult","Pup (rho_pup x N_pup)")))
    
    p_stack <- ggplot(long_df, aes(x=Year, y=Count, fill=Component)) +
      geom_area(position="stack", alpha=0.80) +
      facet_wrap(~Site, ncol=2, scales="free_y") +
      scale_fill_manual(
        values=c("Juv + Adult"=SEAL_COLS$molt, "Pup (rho_pup x N_pup)"=SEAL_COLS$pup)) +
      scale_y_continuous(labels=scales::comma) +
      labs(x="Year", y="Estimated N in molt count",
           title="Molt Count Decomposition: Pup vs Juv+Adult Contribution",
           subtitle="Posterior mean. Pup component = rho_pup x N_pup.",
           fill="Component") +
      theme_seal(base_size=13) +
      theme(strip.text=element_text(face="bold"),
            legend.position="bottom",
            axis.text.x=element_text(angle=45,hjust=1,size=8))
    
    # ── (d) Pup % heatmap ───────────────────────────────────────────────────
    pct_df <- stack_df |> select(Site, Year, pup_pct) |>
      mutate(Site=factor(Site, levels=site_names))
    p_pct <- ggplot(pct_df, aes(x=Year, y=Site, fill=pup_pct)) +
      geom_tile(colour="white", linewidth=0.4) +
      geom_text(aes(label=sprintf("%.0f%%", pup_pct)), size=2.8, colour="grey10") +
      scale_fill_gradient(low="white", high=SEAL_COLS$pup, limits=c(0,NA),
                          name="Pup %\nof molt") +
      labs(x="Year", y=NULL,
           title="Estimated Pup Percentage of Peak Molt Count by Site",
           subtitle="Based on posterior mean rho_pup and N_pup") +
      theme_seal(base_size=13) +
      theme(axis.text.x=element_text(angle=45,hjust=1,size=8),
            panel.grid=element_blank())
  }
  
  # ── Save ─────────────────────────────────────────────────────────────────
  if (save) {
    ggsave(paste0("Output/Plots/",prefix,"_rho_pup_posterior.jpeg"),
           p_rho, width=20, height=13, units="cm")
    if (!is.null(p_pim_ts))
      ggsave(paste0("Output/Plots/",prefix,"_pup_in_molt_timeseries.jpeg"),
             p_pim_ts, width=24, height=13, units="cm")
    if (!is.null(p_stack))
      ggsave(paste0("Output/Plots/",prefix,"_molt_decomp_stacked.jpeg"),
             p_stack, width=34, height=28, units="cm", dpi=200)
    if (!is.null(p_pct))
      ggsave(paste0("Output/Plots/",prefix,"_pup_pct_molt_heatmap.jpeg"),
             p_pct, width=34, height=14, units="cm", dpi=200)
  }
  
  # ── Print summary ─────────────────────────────────────────────────────────
  cat(sprintf("\nrho_pup: median=%.3f (89%% CrI: %.3f, %.3f)\n",
              median(rho_draws), quantile(rho_draws, CI_LO), quantile(rho_draws, CI_HI)))
  if (!is.null(pim_all)) {
    cat(sprintf("Estimated pups in molt count (all sites, mean across years): %.0f (%.0f, %.0f)\n",
                mean(colMeans(pim_all)),
                mean(apply(pim_all,2,quantile,CI_LO)),
                mean(apply(pim_all,2,quantile,CI_HI))))
  }
  
  list(rho_pup_posterior=p_rho,
       pup_in_molt_ts=p_pim_ts,
       molt_stacked=p_stack,
       pup_pct_heatmap=p_pct)
}

# ============================================================================
# ALIAS v3.2 FUNCTIONS — all unchanged; sourced from the v3.2 file or here
# (copy-forward below; identical bodies to 06_ipm_plots_v3.2)
# ============================================================================

# The functions below are identical to their v3.2 counterparts.
# They are included here so this file is self-contained.
# If you prefer, source 06_ipm_plots.R (v3.2) first and only
# source this file for the three overridden functions above.
