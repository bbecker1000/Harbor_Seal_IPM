# ============================================================================
# 06_ipm_plots.R  —  IPM v3.2 PLOTS, TABLES & POST-PROCESSING
# ----------------------------------------------------------------------------
# Companion to 05_ipm_model.R. NO orchestrator here (run_full_analysis_v3.2
# lives only in 05_ipm_model.R — fixes P6). This file provides:
#   load_seal_results()                 — reload fit + reconstruct sim_data
#   run_all_plots_v3.2()                — plot/table orchestrator
#   check_diagnostics_v3.2()            seal_fit_summary()
#   create_trace_plots_v3.2()           check_parameter_recovery_v3.2()
#   create_ppc_plots_v3.2()             create_timeseries_plots_v3.2()
#   create_site_age_timeseries_v3.2()   create_site_panels_v3.2()
#   create_projection_plots_v3.2()      create_effect_plots_v3.2()
#   create_juv_adult_effect_plots_v3.2() create_forest_plot_v3.2()
#   create_covariate_decomposition_plots_v3.2()
#   create_summary_table_v3.2()         save_model_output_v3.2()
#   create_portfolio_analysis_v3.2()    create_synchrony_projections_v3.2()
#
# Usage:
#   source("Code/05_ipm_model.R"); source("Code/06_ipm_plots.R")
#   out <- load_seal_results("IPM_v3.2_real")
#   run_all_plots_v3.2(out$fit, out$sim_data, prefix="IPM_v3.2_real")
# ============================================================================

library(tidyverse)
library(posterior)
library(bayesplot)
library(patchwork)

dir.create("Output",       showWarnings = FALSE)
dir.create("Output/Plots", showWarnings = FALSE)

`%||%` <- function(x, y) if (!is.null(x)) x else y

# ── 89% credible-interval constants + palette + theme ───────────────────────
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

# ============================================================================
# LOAD HELPER — returns list(fit, sim_data) [shape matches orchestrator]
# ============================================================================

load_seal_results <- function(prefix          = "IPM_v3.2_real",
                              fit_path        = NULL,
                              input_data_path = NULL,
                              years           = 1997:2025,
                              T_proj          = 10) {
  if (is.null(fit_path))
    fit_path <- paste0("Output/harbor_seal_", prefix, "_fit.rds")
  if (!file.exists(fit_path))
    stop("Fit file not found: ", fit_path,
         "\nRun run_full_analysis_v3.2() first, or pass fit_path explicitly.")
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
    cat("No input_data_path — using metadata shell (PPC/recovery limited).\n")
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
  
  fit_ok <- tryCatch(is.function(fit$draws), error=function(e) FALSE)
  if (!fit_ok)
    warning("Loaded fit object may lack working methods; try a fresh readRDS().")
  cat("Ready. run_all_plots_v3.2(out$fit, out$sim_data) to regenerate outputs.\n")
  list(fit=fit, sim_data=sim_data)
}

# ============================================================================
# ORCHESTRATOR — run all plots/tables (each step guarded)
# ============================================================================

run_all_plots_v3.2 <- function(fit, sim_data, prefix="IPM_v3.2",
                               save=TRUE, run_recovery=FALSE,
                               run_portfolio=FALSE, run_synchrony=FALSE) {
  
  safe_run <- function(label, expr) {
    cat(sprintf("\n%s\n", label))
    tryCatch(
      withCallingHandlers(expr,
                          warning=function(w) invokeRestart("muffleWarning")),
      error=function(e) { cat(sprintf("  !! FAILED: %s\n", conditionMessage(e))); NULL })
  }
  
  diag <- safe_run("── Diagnostics ──", check_diagnostics_v3.2(fit))
  
  params_candidate <- c(
    "phi_pup_logit","phi_juv_base","phi_adult_F_logit","phi_adult_F_base","delta_adult",
    "fecund_primip","fecund_mature","prop_female","p_male_breed",
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
  
  traces <- safe_run("── Trace plots ──",
                     create_trace_plots_v3.2(fit, params, save=save, prefix=prefix))
  
  rec <- NULL
  if (run_recovery && !is.null(sim_data$true_params))
    rec <- safe_run("── Parameter recovery ──",
                    check_parameter_recovery_v3.2(fit, sim_data, save=save, prefix=prefix))
  else cat("\n── Parameter recovery: skipped (set run_recovery=TRUE on sim) ──\n")
  
  ppc    <- safe_run("── PPC ──", create_ppc_plots_v3.2(fit, sim_data, save=save, prefix=prefix))
  ts     <- safe_run("── Time series ──", create_timeseries_plots_v3.2(fit, sim_data, save=save, prefix=prefix))
  sa     <- safe_run("── Site x age ──", create_site_age_timeseries_v3.2(fit, sim_data, save=save, prefix=prefix))
  sasite <- safe_run("── Site panels ──", create_site_panels_v3.2(fit, sim_data, save=save, prefix=prefix))
  proj   <- safe_run("── Projections ──", create_projection_plots_v3.2(fit, sim_data, save=save, prefix=prefix))
  le     <- safe_run("── Lambda + elasticity ──", create_lambda_elasticity_v3.2(fit, save=save, prefix=prefix))
  moci   <- safe_run("── MOCI collinearity ──", create_moci_collinearity_v3.2(fit, sim_data, save=save, prefix=prefix))
  eff    <- safe_run("── Effect plots ──", create_effect_plots_v3.2(fit, save=save, prefix=prefix))
  jveff  <- safe_run("── Juv+adult effects ──", create_juv_adult_effect_plots_v3.2(fit, save=save, prefix=prefix))
  forest <- safe_run("── Forest plot ──", create_forest_plot_v3.2(fit, save=save, prefix=prefix))
  decomp <- safe_run("── Decomposition ──", create_covariate_decomposition_plots_v3.2(fit, sim_data, save=save, prefix=prefix))
  tbl    <- safe_run("── Summary table ──", create_summary_table_v3.2(fit, save=save, prefix=prefix))
  safe_run("── Key results ──", save_model_output_v3.2(fit, prefix=prefix))
  
  port <- sync <- NULL
  if (run_portfolio)
    port <- safe_run("── Portfolio ──", create_portfolio_analysis_v3.2(fit, sim_data, save=save, prefix=prefix))
  else cat("\n── Portfolio: skipped (set run_portfolio=TRUE) ──\n")
  if (run_synchrony)
    sync <- safe_run("── Synchrony ──", create_synchrony_projections_v3.2(fit, sim_data, save=save, prefix=prefix))
  else cat("\n── Synchrony: skipped (set run_synchrony=TRUE) ──\n")
  
  results <- list(diagnostics=diag, traces=traces, recovery=rec, ppc=ppc, ts=ts,
                  site_age=sa, site_panels=sasite, projections=proj, effects=eff,
                  juv_adult_effects=jveff, forest=forest, decomposition=decomp,
                  lambda_elast=le, moci_collin=moci, table=tbl, portfolio=port, sync=sync)
  skipped <- c(
    if (!run_recovery || is.null(sim_data$true_params)) "recovery" else character(0),
    if (!run_portfolio) "portfolio" else character(0),
    if (!run_synchrony) "sync" else character(0))
  attempted <- results[!names(results) %in% skipped]
  n_ok   <- sum(!sapply(attempted, is.null))
  n_fail <- sum( sapply(attempted, is.null))
  cat(sprintf("\n── Pipeline: %d/%d steps succeeded ──\n", n_ok, length(attempted)))
  if (n_fail > 0)
    cat("   Failed:", paste(names(attempted)[sapply(attempted, is.null)], collapse=", "), "\n")
  invisible(results)
}

# ── 89% CrI fit summary (defensive) ─────────────────────────────────────────
seal_fit_summary <- function(fit, variables) {
  stale_msg <- paste0("fit methods unavailable (stale/RDS-reloaded fit).\n",
                      "Reload: out <- load_seal_results(\"IPM_v3.2_real\")")
  s <- tryCatch(fit$summary(variables=variables),
                error=function(e) stop(stale_msg, "\nOriginal: ", conditionMessage(e)))
  d <- tryCatch(fit$draws(variables=variables, format="matrix"),
                error=function(e) stop(stale_msg, "\nOriginal: ", conditionMessage(e)))
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
# PART 4: DIAGNOSTICS
# ============================================================================
check_diagnostics_v3.2 <- function(fit) {
  cat("\n=== MODEL DIAGNOSTICS — IPM v3.2 ===\n")
  tryCatch({
    ds <- fit$diagnostic_summary(quiet=TRUE)
    cat(if (all(ds$num_max_treedepth==0)) "Treedepth OK.\n"
        else sprintf("WARN: %d treedepth hits.\n", sum(ds$num_max_treedepth)))
    cat(if (all(ds$num_divergent==0)) "No divergences.\n"
        else sprintf("WARN: %d divergences.\n", sum(ds$num_divergent)))
    cat(if (all(ds$ebfmi>0.2)) "E-BFMI OK.\n"
        else sprintf("WARN: low E-BFMI in %d chain(s).\n", sum(ds$ebfmi<=0.2)))
  }, error=function(e)
    cat("NOTE: sampler diagnostics unavailable (CSVs gone). Reload fit.\n"))
  
  params <- c(
    "phi_pup_logit","phi_juv_base","phi_adult_F_logit","phi_adult_F_base","delta_adult",
    "fecund_primip","fecund_mature","prop_female","p_male_breed",
    "beta_coy[1]","beta_coy[2]","beta_coy[3]",
    paste0("beta_dist_surv[",1:6,"]"),
    "beta_moci_ond_fecund","beta_moci_amj_pup","beta_moci_jfm_juv",
    "beta_moci_jfm_adult","beta_eseal_pup",
    "sigma_process","sigma_obs_adult","sigma_obs_pup","sigma_obs_molt","sigma_site")
  
  s <- seal_fit_summary(fit, params)
  cat("\nParameter Summary:\n")
  print(s |> select(variable,mean,sd,q_lo,q_hi,rhat,ess_bulk), n=nrow(s))
  
  pup_logit_draws <- fit$draws(variables="phi_pup_logit", format="df")$phi_pup_logit
  cat(sprintf("\nphi_pup_base (prob): median=%.3f, %s=[%.3f, %.3f]\n",
              median(plogis(pup_logit_draws)), CI_LABEL,
              quantile(plogis(pup_logit_draws), CI_LO),
              quantile(plogis(pup_logit_draws), CI_HI)))
  
  rhat_v <- tryCatch(as.double(as.vector(unclass(s$rhat))), error=function(e) rep(NA_real_, nrow(s)))
  ess_v  <- tryCatch(as.double(as.vector(unclass(s$ess_bulk))), error=function(e) rep(NA_real_, nrow(s)))
  bad <- s[!is.na(rhat_v) & rhat_v > 1.05, , drop=FALSE]
  low <- s[!is.na(ess_v)  & ess_v  < 400,  , drop=FALSE]
  if (nrow(bad)>0) { cat("\nWARNING Rhat>1.05:\n"); print(bad[, c('variable','rhat')]) }
  if (nrow(low)>0) { cat("\nWARNING low ESS:\n");    print(low[, c('variable','ess_bulk')]) }
  
  list(params=params, summary=s)
}

# ============================================================================
# PART 5: TRACE PLOTS
# ============================================================================
create_trace_plots_v3.2 <- function(fit, params, save=TRUE, prefix="IPM_v3.2") {
  draws <- fit$draws(format="df")
  grp <- function(pars, title, file, w=30, h=18) {
    pars <- pars[pars %in% colnames(draws)]
    if (!length(pars)) return(NULL)
    p <- mcmc_trace(draws, pars=pars) + labs(title=title)
    if (save) ggsave(paste0("Output/Plots/",prefix,"_",file,".jpeg"), p, width=w, height=h, units="cm")
    p
  }
  list(
    survival   = grp(c("phi_pup_logit","phi_juv_base","phi_adult_F_logit","delta_adult","p_male_breed"),
                     "Trace: Survival + Observation","trace_survival"),
    fecundity  = grp(c("fecund_primip","fecund_mature","prop_female"),"Trace: Fecundity","trace_fecundity"),
    coyote     = grp(c("beta_coy[1]","beta_coy[2]","beta_coy[3]"),"Trace: Coyote","trace_coyote", h=12),
    disturbance= grp(paste0("beta_dist_surv[",1:6,"]"),"Trace: Disturbance","trace_disturbance"),
    moci       = grp(c("beta_moci_ond_fecund","beta_moci_amj_pup","beta_moci_jfm_juv",
                       "beta_moci_jfm_adult","beta_moci_amj_molt"),"Trace: MOCI","trace_moci"),
    errors     = grp(c("sigma_process","sigma_obs_adult","sigma_obs_pup","sigma_obs_molt","sigma_site"),
                     "Trace: Error Terms","trace_errors"))
}

# ============================================================================
# PART 6: PARAMETER RECOVERY
# ============================================================================
check_parameter_recovery_v3.2 <- function(fit, sim_data, save=TRUE, prefix="IPM_v3.2") {
  tp <- sim_data$true_params
  true_vals <- tibble(
    parameter = c("phi_pup_logit","phi_juv_base","phi_adult_F_logit","delta_adult",
                  "fecund_primip","fecund_mature","prop_female","p_male_breed",
                  "beta_coy[1]","beta_coy[2]","beta_coy[3]",
                  paste0("beta_dist_surv[",1:6,"]"),
                  "beta_moci_ond_fecund","beta_moci_ond_pup","beta_moci_jfm_pup",
                  "beta_moci_amj_pup","beta_moci_jfm_juv","beta_moci_jfm_adult","beta_eseal_pup",
                  "detect_breed_logit","detect_molt_logit",
                  "sigma_process","sigma_obs_adult","sigma_obs_pup","sigma_obs_molt"),
    true_value = c(tp$phi_pup_logit, tp$phi_juv_base, tp$phi_adult_F_logit, tp$delta_adult,
                   tp$fecund_primip, tp$fecund_mature, tp$prop_female, tp$p_male_breed,
                   tp$beta_coy[1], tp$beta_coy[2], tp$beta_coy[3],
                   tp$beta_dist_surv[1], tp$beta_dist_surv[2], tp$beta_dist_surv[3],
                   tp$beta_dist_surv[4], tp$beta_dist_surv[5], tp$beta_dist_surv[6],
                   tp$beta_moci_ond_fecund, tp$beta_moci_ond_pup, tp$beta_moci_jfm_pup,
                   tp$beta_moci_amj_pup, tp$beta_moci_jfm_juv, tp$beta_moci_jfm_adult, tp$beta_eseal_pup,
                   tp$detect_breed_logit, tp$detect_molt_logit,
                   tp$sigma_process, tp$sigma_obs_adult, tp$sigma_obs_pup, tp$sigma_obs_molt))
  
  rec <- seal_fit_summary(fit, true_vals$parameter) |>
    left_join(true_vals, by=c("variable"="parameter")) |>
    mutate(recovered    = true_value >= q_lo & true_value <= q_hi,
           rel_bias_pct = (mean - true_value) / abs(true_value) * 100)
  
  cat("Recovery:", sum(rec$recovered), "/", nrow(rec),
      sprintf("(%.1f%%)\n", 100*mean(rec$recovered)))
  cat("NOTE[P2]: true values = prior means, so coverage over-states recovery.\n")
  
  rec <- rec |> mutate(
    identifiability = case_when(
      variable %in% c("phi_adult_F_logit","fecund_mature","beta_moci_ond_fecund",
                      "detect_breed_logit","sigma_obs_adult") ~ "Well identified",
      variable %in% c("phi_juv_base","beta_coy[1]","beta_coy[2]","beta_coy[3]",
                      "beta_moci_jfm_adult","detect_molt_logit",
                      "sigma_process","sigma_obs_molt") ~ "Moderately identified",
      TRUE ~ "Prior dominated"),
    identifiability = factor(identifiability,
                             levels=c("Well identified","Moderately identified","Prior dominated")))
  
  p <- ggplot(rec, aes(x=true_value, y=mean)) +
    geom_abline(slope=1, intercept=0, linetype=2, color="gray50") +
    geom_pointrange(aes(ymin=q_lo, ymax=q_hi, color=recovered, shape=identifiability), size=0.8) +
    geom_text(aes(label=variable), hjust=-0.1, vjust=-0.3, size=2.5, check_overlap=TRUE) +
    scale_color_manual(values=c("TRUE"=SEAL_COLS$pop,"FALSE"=SEAL_COLS$adult_f),
                       name=paste("True in", CI_LABEL)) +
    scale_shape_manual(values=c("Well identified"=16,"Moderately identified"=17,"Prior dominated"=1),
                       name="Identifiability") +
    facet_wrap(~identifiability, scales="free") +
    labs(x="True Value", y=paste0("Posterior Mean (", CI_LABEL, ")"),
         title="Parameter Recovery: IPM v3.2",
         subtitle=sprintf("Coverage %d/%d (%.0f%%); mean |bias| %.1f%% — NOTE: true=prior means (P2)",
                          sum(rec$recovered), nrow(rec), 100*mean(rec$recovered),
                          mean(abs(rec$rel_bias_pct), na.rm=TRUE))) +
    theme_seal()
  if (save) ggsave(paste0("Output/Plots/",prefix,"_parameter_recovery.jpeg"),
                   p, width=30, height=22, units="cm")
  
  rec_tbl <- rec |> mutate(
    CrI=sprintf("(%.3f, %.3f)", q_lo, q_hi),
    Bias=sprintf("%.1f%%", rel_bias_pct),
    Recovered=ifelse(recovered,"Yes","No")) |>
    select(variable, true_value, mean, CrI, Bias, Recovered, identifiability) |>
    arrange(identifiability, variable)
  if (save) write_csv(rec_tbl, paste0("Output/",prefix,"_parameter_recovery_table.csv"))
  list(table=rec, plot=p, manuscript_table=rec_tbl)
}

# ============================================================================
# PART 7: PPC
# ============================================================================
create_ppc_plots_v3.2 <- function(fit, sim_data, save=TRUE, prefix="IPM_v3.2") {
  rep_a <- fit$draws(variables="y_adult_rep", format="matrix")
  rep_p <- fit$draws(variables="y_pup_rep",   format="matrix")
  rep_m <- fit$draws(variables="y_molt_rep",  format="matrix")
  obs_a <- as.vector(t(sim_data$stan_data$y_adult))
  obs_p <- as.vector(t(sim_data$stan_data$y_pup))
  obs_m <- as.vector(t(sim_data$stan_data$y_molt))
  ind_a <- as.vector(t(sim_data$stan_data$y_adult_obs))==1
  ind_p <- as.vector(t(sim_data$stan_data$y_pup_obs))  ==1
  ind_m <- as.vector(t(sim_data$stan_data$y_molt_obs)) ==1
  
  p_comb <- (ppc_dens_overlay(obs_a[ind_a],rep_a[1:100,ind_a]) + labs(title="PPC: Adult")) /
    (ppc_dens_overlay(obs_p[ind_p],rep_p[1:100,ind_p]) + labs(title="PPC: Pup")) /
    (ppc_dens_overlay(obs_m[ind_m],rep_m[1:100,ind_m]) + labs(title="PPC: Molt"))
  if (save) ggsave(paste0("Output/Plots/",prefix,"_ppc_density.jpeg"),
                   p_comb, width=25, height=30, units="cm")
  
  sg <- rep(sim_data$site_names, each=sim_data$stan_data$T)
  p_site <- ppc_stat_grouped(obs_a[ind_a], rep_a[1:100,ind_a], group=sg[ind_a], stat="mean") +
    labs(title="PPC: Mean Adult by Site")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_ppc_by_site.jpeg"),
                   p_site, width=25, height=15, units="cm")
  list(density=p_comb, by_site=p_site)
}

# ============================================================================
# PART 8: TIME SERIES
# ============================================================================
create_timeseries_plots_v3.2 <- function(fit, sim_data, save=TRUE, prefix="IPM_v3.2") {
  years <- sim_data$years; site_names <- sim_data$site_names
  S <- length(site_names); T <- length(years)
  
  Ntot <- fit$draws(variables="N_total_all", format="df") |> select(starts_with("N_total_all"))
  plot_df <- tibble(Year=years, mean=colMeans(Ntot),
                    lo=as.numeric(apply(Ntot,2,quantile,CI_LO)),
                    hi=as.numeric(apply(Ntot,2,quantile,CI_HI)))
  p_total <- ggplot(plot_df, aes(x=Year)) +
    geom_ribbon(aes(ymin=lo,ymax=hi),alpha=0.25,fill=SEAL_COLS$ribbon) +
    geom_line(aes(y=mean),linewidth=1.2,color=SEAL_COLS$pop) +
    scale_y_continuous(limits=c(0, ceiling(max(plot_df$hi)/1000)*1000),
                       labels=scales::comma, expand=c(0,0)) +
    labs(x="Year", y="Total Population", title="Estimated Total Harbor Seal Population") +
    theme_seal()
  if (save) ggsave(paste0("Output/Plots/",prefix,"_total_population.jpeg"),
                   p_total, width=25, height=15, units="cm")
  
  pup_logit <- fit$draws(variables="phi_pup_logit",format="df")$phi_pup_logit
  pup_prob  <- plogis(pup_logit)
  p_phipup  <- ggplot(tibble(x=pup_prob), aes(x=x)) +
    geom_density(fill=SEAL_COLS$pup,alpha=0.5) +
    geom_vline(xintercept=median(pup_prob),linetype="dashed",color="darkgreen") +
    geom_vline(xintercept=0.5, linetype="dotted", color="gray50") +
    annotate("text",x=0.5,y=Inf,label=" Prior mean 0.50",
             hjust=0,vjust=1.5,size=3.5,color="gray40") +   # P21: relabel (was "Field data")
    labs(x="Pup survival probability", y="Density",
         title="Posterior: Pup Survival (sex-neutral)",
         subtitle=sprintf(paste0("Median = %.3f  (", CI_FMT, ": %.3f-%.3f)"),
                          median(pup_prob), quantile(pup_prob,CI_LO), quantile(pup_prob,CI_HI))) +
    theme_seal()
  if (save) ggsave(paste0("Output/Plots/",prefix,"_phi_pup_posterior.jpeg"),
                   p_phipup, width=20, height=12, units="cm")
  
  pmb <- fit$draws(variables="p_male_breed",format="df")$p_male_breed
  p_pmb <- ggplot(tibble(x=pmb),aes(x=x)) +
    geom_density(fill=SEAL_COLS$pop,alpha=0.5) +
    geom_vline(xintercept=median(pmb),linetype="dashed",color="#08306b") +
    labs(x="p_male_breed",y="Density",
         title="Posterior: Male Haul-out Fraction (Breeding Season)",
         subtitle=sprintf(paste0("Median=%.3f  (", CI_FMT, ": %.3f-%.3f)"),
                          median(pmb),quantile(pmb,CI_LO),quantile(pmb,CI_HI))) +
    theme_seal()
  if (save) ggsave(paste0("Output/Plots/",prefix,"_p_male_breed.jpeg"),
                   p_pmb, width=20, height=12, units="cm")
  
  sr_t <- fit$draws(variables="sex_ratio_adult",   format="matrix")
  sr_o <- fit$draws(variables="sex_ratio_observed",format="matrix")
  make_sr <- function(mat, var) {
    m <- sapply(1:T, function(t) rowMeans(mat[,paste0(var,"[",1:S,",",t,"]"),drop=FALSE]))
    tibble(Year=years, mean=colMeans(m),
           lo=as.numeric(apply(m,2,quantile,CI_LO)), hi=as.numeric(apply(m,2,quantile,CI_HI)))
  }
  sr_true_df <- tryCatch(make_sr(sr_t,"sex_ratio_adult"),    error=function(e) NULL)
  sr_obs_df  <- tryCatch(make_sr(sr_o,"sex_ratio_observed"), error=function(e) NULL)
  p_sr <- ggplot() +
    geom_hline(yintercept=0.5,linetype="dotted",color="gray50") +
    { if (!is.null(sr_obs_df))  geom_ribbon(data=sr_obs_df,  aes(x=Year,ymin=lo,ymax=hi),alpha=0.15,fill="orange") } +
    { if (!is.null(sr_obs_df))  geom_line(data=sr_obs_df,    aes(x=Year,y=mean,color="Observed (spring)"),linewidth=1.1,linetype="dashed") } +
    { if (!is.null(sr_true_df)) geom_ribbon(data=sr_true_df, aes(x=Year,ymin=lo,ymax=hi),alpha=0.25,fill="purple") } +
    { if (!is.null(sr_true_df)) geom_line(data=sr_true_df,   aes(x=Year,y=mean,color="True (population)"),linewidth=1.2) } +
    scale_color_manual(values=c("True (population)"=SEAL_COLS$juv,"Observed (spring)"=SEAL_COLS$warn)) +
    labs(x="Year",y="Proportion female",color=NULL,
         title="Adult Sex Ratio: True vs Spring Survey Observation",
         subtitle="Observed > true because most males remain in water during breeding") +
    ylim(0.45,0.80) + theme_seal()
  if (save) ggsave(paste0("Output/Plots/",prefix,"_sex_ratio.jpeg"),
                   p_sr, width=25, height=15, units="cm")
  
  mp <- fit$draws(variables="mean_phi_pup", format="matrix")
  mp_df <- tibble(Year=years, mean=colMeans(mp),
                  lo=as.numeric(apply(mp,2,quantile,CI_LO)), hi=as.numeric(apply(mp,2,quantile,CI_HI)))
  p_phipup_t <- ggplot(mp_df,aes(x=Year)) +
    geom_ribbon(aes(ymin=lo,ymax=hi),alpha=0.25,fill=SEAL_COLS$pup) +
    geom_line(aes(y=mean),linewidth=1.2,color=SEAL_COLS$pup) +
    geom_hline(yintercept=0.5,linetype="dotted",color="gray50") +
    labs(x="Year",y="Pup survival (sex-neutral)",
         title="Mean Pup Survival Over Time (across all sites)") +
    theme_seal()
  if (save) ggsave(paste0("Output/Plots/",prefix,"_phi_pup_timeseries.jpeg"),
                   p_phipup_t, width=22, height=12, units="cm")
  
  list(total=p_total, phi_pup=p_phipup, p_male_breed=p_pmb,
       sex_ratio=p_sr, phi_pup_time=p_phipup_t)
}

# ============================================================================
# PART 9: SITE-BY-AGE TIME SERIES
# ============================================================================
create_site_age_timeseries_v3.2 <- function(fit, sim_data, save=TRUE, prefix="IPM_v3.2") {
  years <- sim_data$years; site_names <- sim_data$site_names
  S <- length(site_names); T <- length(years)
  pull_st <- function(var) {
    d <- fit$draws(variables=var, format="matrix")
    map_dfr(1:S, function(s) map_dfr(1:T, function(t) {
      cn <- paste0(var,"[",s,",",t,"]")
      if (!cn %in% colnames(d)) return(NULL)
      tibble(Site=site_names[s], Year=years[t],
             mean=mean(d[,cn]), lo=as.numeric(quantile(d[,cn],CI_LO)), hi=as.numeric(quantile(d[,cn],CI_HI)))
    }))
  }
  all_sum <- bind_rows(
    pull_st("N_pup")         |> mutate(Age_Class="Pup"),
    pull_st("N_juv_total")   |> mutate(Age_Class="Juvenile"),
    pull_st("N_adult_total") |> mutate(Age_Class="Adult")
  ) |> mutate(Age_Class=factor(Age_Class,levels=c("Pup","Juvenile","Adult")))
  p <- ggplot(all_sum,aes(x=Year,y=mean,color=Site,fill=Site)) +
    geom_ribbon(aes(ymin=lo,ymax=hi),alpha=0.15,color=NA) +
    geom_line(linewidth=1) + facet_wrap(~Age_Class,scales="free_y",ncol=1) +
    expand_limits(y=0) + scale_y_continuous(labels=scales::comma) +
    labs(x="Year",y="Population Size",title="Population by Age Class Across Sites") +
    theme_seal()
  if (save) ggsave(paste0("Output/Plots/",prefix,"_age_class_timeseries.jpeg"),
                   p, width=30, height=35, units="cm")
  list(by_age=p, data=all_sum)
}

# ── 9b: one panel per site ──────────────────────────────────────────────────
create_site_panels_v3.2 <- function(fit, sim_data, save=TRUE, prefix="IPM_v3.2") {
  years <- sim_data$years; site_names <- sim_data$site_names
  S <- length(site_names); T <- length(years)
  age_cols <- c(Pup=SEAL_COLS$pup, Juvenile=SEAL_COLS$juv, Adult=SEAL_COLS$pop)
  pull_st <- function(var, label) {
    d <- fit$draws(variables=var, format="matrix")
    map_dfr(1:S, function(s) map_dfr(1:T, function(t) {
      cn <- paste0(var,"[",s,",",t,"]")
      if (!cn %in% colnames(d)) return(NULL)
      tibble(Site=site_names[s], Year=years[t], Age_Class=label,
             mean=mean(d[,cn]), lo=as.numeric(quantile(d[,cn],CI_LO)), hi=as.numeric(quantile(d[,cn],CI_HI)))
    }))
  }
  all_df <- bind_rows(
    pull_st("N_pup","Pup"), pull_st("N_juv_total","Juvenile"), pull_st("N_adult_total","Adult")
  ) |> mutate(Age_Class=factor(Age_Class,levels=c("Pup","Juvenile","Adult")),
              Site=factor(Site,levels=site_names))
  p <- ggplot(all_df, aes(x=Year, y=mean, colour=Age_Class, fill=Age_Class)) +
    geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.15, colour=NA) +
    geom_line(linewidth=0.9) + facet_grid(Site ~ ., scales="free_y") +
    scale_colour_manual(values=age_cols, guide="none") +
    scale_fill_manual(values=age_cols, guide="none") +
    labs(x="Year", y="Population Size", title="Harbor Seal Population by Site and Age Class",
         subtitle=paste0("Posterior mean +/- ", CI_LABEL)) +
    theme_seal(base_size=14) +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=8),
          strip.text.y=element_text(size=10,face="bold"), panel.spacing=unit(0.4,"lines"))
  if (save) ggsave(paste0("Output/Plots/",prefix,"_site_age_panels.jpeg"),
                   p, width=36, height=6*S, units="cm", dpi=200)
  list(plot=p, data=all_df)
}

# ============================================================================
# PART 10: PROJECTIONS
# ============================================================================
create_projection_plots_v3.2 <- function(fit, sim_data, save=TRUE, prefix="IPM_v3.2") {
  years <- sim_data$years; T <- length(years)
  T_proj <- sim_data$stan_data$T_proj; N_sc <- sim_data$stan_data$N_scenarios
  sc_nm <- sim_data$scenario_names; pyrs <- (max(years)+1):(max(years)+T_proj)
  all_d <- fit$draws(format="matrix")
  proj_df <- map_dfr(1:N_sc, function(sc) {
    cols <- grep(paste0("^N_total_all_proj\\[",sc,","), colnames(all_d))
    if (!length(cols)) return(NULL)
    pm <- all_d[,cols]
    tibble(Scenario=sc_nm[sc], Year=pyrs, mean=colMeans(pm),
           lo=as.numeric(apply(pm,2,quantile,CI_LO)), hi=as.numeric(apply(pm,2,quantile,CI_HI)))
  })
  Ntot <- fit$draws(variables="N_total_all",format="df") |> select(starts_with("N_total_all"))
  hist <- tibble(Scenario="Historical",Year=years,mean=colMeans(Ntot),
                 lo=as.numeric(apply(Ntot,2,quantile,CI_LO)),hi=as.numeric(apply(Ntot,2,quantile,CI_HI)))
  full <- bind_rows(hist,proj_df) |> mutate(Period=ifelse(Scenario=="Historical","Historical","Projection"))
  p <- ggplot() +
    geom_ribbon(data=filter(full,Period=="Historical"), aes(x=Year,ymin=lo,ymax=hi),alpha=0.3,fill="gray50") +
    geom_line(data=filter(full,Period=="Historical"), aes(x=Year,y=mean),linewidth=1.2,color="black") +
    geom_ribbon(data=filter(full,Period=="Projection"), aes(x=Year,ymin=lo,ymax=hi,fill=Scenario),alpha=0.2) +
    geom_line(data=filter(full,Period=="Projection"), aes(x=Year,y=mean,color=Scenario),linewidth=1.2) +
    geom_vline(xintercept=max(years),linetype=2,color="red") +
    scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2") +
    scale_y_continuous(limits=c(0, ceiling(max(full$hi,na.rm=TRUE)/1000)*1000),
                       labels=scales::comma, expand=c(0,0)) +
    labs(x="Year",y="Total Population", title="10-Year Projections",
         subtitle=paste0("Bands = ", CI_LABEL, "; dashed = projection start (process error excluded)")) +
    theme_seal()
  if (save) ggsave(paste0("Output/Plots/",prefix,"_projections.jpeg"),
                   p, width=30, height=20, units="cm")
  list(projection=p, data=proj_df)
}

# ============================================================================
# PART 11: COVARIATE EFFECT PLOTS  (reconstruct from logit + beta)
# ============================================================================
.make_effect_engine <- function(all_draws) {
  function(param, xlab, title, stage="pup",
           xr=seq(-2,2,length.out=100), ylims=NULL) {
    beta_v <- tryCatch(as.numeric(posterior::extract_variable(all_draws, param)),
                       error=function(e){warning(sprintf("No draws for '%s'",param)); numeric(0)})
    if (length(beta_v)==0L)
      return(list(plot=ggplot()+labs(title=paste("MISSING:",param))+theme_seal(), yr=c(0,1)))
    logit_v <- switch(stage,
                      pup     = as.numeric(posterior::extract_variable(all_draws,"phi_pup_logit")),
                      juv     = qlogis(as.numeric(posterior::extract_variable(all_draws,"phi_juv_base"))),
                      adult_F = as.numeric(posterior::extract_variable(all_draws,"phi_adult_F_logit")),
                      adult_M = { aF<-as.numeric(posterior::extract_variable(all_draws,"phi_adult_F_base"))
                      del<-as.numeric(posterior::extract_variable(all_draws,"delta_adult"))
                      qlogis(pmax(aF-del,0.001)) },
                      fecund  = { af<-as.numeric(posterior::extract_variable(all_draws,"avg_fecundity"))
                      qlogis(pmax(pmin(af,0.999),0.001)) })
    y_label <- switch(stage, pup="Pup Survival", juv="Juvenile Survival",
                      adult_F="Adult Female Survival", adult_M="Adult Male Survival",
                      fecund="Fecundity (prob. pupping)")
    idx <- sample(seq_along(beta_v), min(500L, length(beta_v)))
    df <- do.call(rbind, lapply(idx, function(i)
      data.frame(cov_val=xr, survival=plogis(logit_v[i]+beta_v[i]*xr))))
    agg <- aggregate(survival ~ cov_val, data=df,
                     FUN=function(v) c(mean=mean(v), lo=as.numeric(quantile(v,CI_LO)), hi=as.numeric(quantile(v,CI_HI))))
    sm <- data.frame(cov_val=agg$cov_val, mn=agg$survival[,"mean"],
                     lo=agg$survival[,"lo"], hi=agg$survival[,"hi"])
    clr <- ifelse(mean(beta_v)<0,"red3","blue3"); base_surv <- mean(plogis(logit_v))
    p <- ggplot(sm, aes(x=.data[["cov_val"]])) +
      geom_ribbon(aes(ymin=.data[["lo"]],ymax=.data[["hi"]]),alpha=.2,fill=clr) +
      geom_line(aes(y=.data[["mn"]]),linewidth=1.2,color=clr) +
      geom_hline(yintercept=base_surv,linetype=2,color="gray50") +
      geom_vline(xintercept=0,linetype=2,color="gray50") +
      labs(x=xlab,y=y_label,title=title) + theme_seal()
    if (!is.null(ylims)) p <- p + coord_cartesian(ylim=ylims)
    list(plot=p, yr=c(min(sm$lo), max(sm$hi)))
  }
}
.get_draws <- function(fit) tryCatch(fit$draws(format="draws_df"),
                                     error=function(e) tryCatch(fit$draws(format="matrix"),
                                                                error=function(e2) stop("Cannot extract draws; pass out$fit, not a saved list.")))

create_effect_plots_v3.2 <- function(fit, save=TRUE, prefix="IPM_v3.2") {
  site_names <- c("BL","DE","DP","PRH","TB","TP")
  mk <- .make_effect_engine(.get_draws(fit))
  
  coy <- lapply(1:3, function(i) mk(paste0("beta_coy[",i,"]"),"Coyote (SD)",
                                    paste0("(",c("BL","DE","DP")[i],")")))
  ylc <- c(max(0,min(sapply(coy,`[[`,"yr"))-0.03), min(1,max(sapply(coy,`[[`,"yr"))+0.03))
  p_coy <- wrap_plots(lapply(1:3, function(i)
    mk(paste0("beta_coy[",i,"]"),"Coyote (SD)",paste0("(",c("BL","DE","DP")[i],")"),ylims=ylc)$plot), ncol=3) +
    plot_annotation(title="Site-Specific Coyote Effects on Pup Survival")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_effects_coyote.jpeg"), p_coy, width=36, height=12, units="cm")
  
  dst <- lapply(1:6, function(s) mk(paste0("beta_dist_surv[",s,"]"),"Disturbance (SD)",paste0("(",site_names[s],")")))
  yld <- c(max(0,min(sapply(dst,`[[`,"yr"))-0.03), min(1,max(sapply(dst,`[[`,"yr"))+0.03))
  p_dst <- wrap_plots(lapply(1:6, function(s)
    mk(paste0("beta_dist_surv[",s,"]"),"Disturbance (SD)",paste0("(",site_names[s],")"),ylims=yld)$plot), ncol=3) +
    plot_annotation(title="Site-Specific Disturbance Effects on Pup Survival")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_effects_disturbance.jpeg"), p_dst, width=36, height=24, units="cm")
  
  mci_specs <- list(c("beta_moci_ond_fecund","MOCI Fall (SD)","OND -> Fecundity","fecund"),
                    c("beta_moci_ond_pup","MOCI Fall (SD)","OND -> Pup","pup"),
                    c("beta_moci_jfm_pup","MOCI Winter (SD)","JFM -> Pup","pup"),
                    c("beta_moci_amj_pup","MOCI Spring (SD)","AMJ -> Pup","pup"),
                    c("beta_eseal_pup","Elephant Seal (SD)","Eseal -> Pup","pup"))
  mci <- lapply(mci_specs, function(x) mk(x[1],x[2],x[3],stage=x[4]))
  ylm <- c(max(0,min(sapply(mci,`[[`,"yr"))-0.03), min(1,max(sapply(mci,`[[`,"yr"))+0.03))
  p_mci <- wrap_plots(lapply(mci_specs, function(x) mk(x[1],x[2],x[3],stage=x[4],ylims=ylm)$plot), ncol=3) +
    plot_annotation(title="Shared Covariate Effects on Pup Survival")
  if (save) ggsave(paste0("Output/Plots/",prefix,"_effects_moci.jpeg"), p_mci, width=36, height=24, units="cm")
  
  list(coyote=p_coy, disturbance=p_dst, moci_pup=p_mci)
}

create_juv_adult_effect_plots_v3.2 <- function(fit, save=TRUE, prefix="IPM_v3.2") {
  mk <- .make_effect_engine(.get_draws(fit))
  j1 <- mk("beta_moci_jfm_juv","MOCI Winter JFM (SD)","Juvenile",stage="juv")
  a1 <- mk("beta_moci_jfm_adult","MOCI Winter JFM (SD)","Adult Female",stage="adult_F")
  a2 <- mk("beta_moci_jfm_adult","MOCI Winter JFM (SD)","Adult Male",stage="adult_M")
  all_yr <- c(j1$yr,a1$yr,a2$yr); yall <- c(max(0,min(all_yr)-0.02), min(1,max(all_yr)+0.02))
  p_compare <- wrap_plots(
    mk("beta_moci_jfm_juv","MOCI Winter JFM (SD)","Juvenile",stage="juv",ylims=yall)$plot,
    mk("beta_moci_jfm_adult","MOCI Winter JFM (SD)","Adult Female",stage="adult_F",ylims=yall)$plot,
    mk("beta_moci_jfm_adult","MOCI Winter JFM (SD)","Adult Male",stage="adult_M",ylims=yall)$plot, ncol=3) +
    plot_annotation(title="MOCI Winter (JFM): Juvenile and Adult Survival",
                    subtitle=paste0("Shared y-axis; dashed = posterior mean baseline; ", CI_LABEL))
  if (save) ggsave(paste0("Output/Plots/",prefix,"_effects_juv_adult_moci_comparison.jpeg"),
                   p_compare, width=36, height=14, units="cm")
  list(comparison=p_compare, juv_moci=j1$plot)
}

# ============================================================================
# PART 11b: FOREST PLOT  (prior bands corrected to match Stan — fixes P3)
# ============================================================================
create_forest_plot_v3.2 <- function(fit, save=TRUE, prefix="IPM_v3.2") {
  params <- tribble(
    ~variable,               ~label,                          ~group,
    "beta_moci_ond_fecund",  "MOCI Fall (OND) -> fecundity",  "MOCI",
    "beta_moci_ond_pup",     "MOCI Fall OND (t) -> pup",      "MOCI",
    "beta_moci_jfm_pup",     "MOCI Winter JFM (t) -> pup",    "MOCI",
    "beta_moci_amj_pup",     "MOCI Spring AMJ (t-1) -> pup",  "MOCI",
    "beta_moci_jfm_juv",     "MOCI Winter -> juvenile",       "MOCI",
    "beta_moci_jfm_adult",   "MOCI Winter -> adult",          "MOCI",
    "beta_moci_amj_molt",    "MOCI Spring -> molt detect",    "MOCI",
    "detect_breed_logit",    "Detect breed baseline (logit)", "Detection",
    "detect_molt_logit",     "Detect molt baseline (logit)",  "Detection",
    "beta_coy[1]",           "Coyote -> pup (BL)",            "Coyote",
    "beta_coy[2]",           "Coyote -> pup (DE)",            "Coyote",
    "beta_coy[3]",           "Coyote -> pup (DP)",            "Coyote",
    "beta_dist_surv[1]",     "Disturbance -> pup (BL)",       "Disturbance",
    "beta_dist_surv[2]",     "Disturbance -> pup (DE)",       "Disturbance",
    "beta_dist_surv[3]",     "Disturbance -> pup (DP)",       "Disturbance",
    "beta_dist_surv[4]",     "Disturbance -> pup (PRH)",      "Disturbance",
    "beta_dist_surv[5]",     "Disturbance -> pup (TB)",       "Disturbance",
    "beta_dist_surv[6]",     "Disturbance -> pup (TP)",       "Disturbance",
    "beta_eseal_pup",        "Elephant seal -> pup",          "Elephant seal")
  
  grp_cols <- c(MOCI="#2166AC", Detection="#4D9221", Coyote="#B2182B",
                Disturbance="#8C510A", "Elephant seal"="#762A83")
  
  all_draws <- .get_draws(fit)
  post <- map_dfr(seq_len(nrow(params)), function(i) {
    v <- tryCatch(as.numeric(posterior::extract_variable(all_draws, params$variable[i])),
                  error=function(e) numeric(0))
    if (length(v)==0L) return(NULL)
    tibble(variable=params$variable[i], mean=mean(v),
           lo89=as.numeric(quantile(v,CI_LO)), hi89=as.numeric(quantile(v,CI_HI)),
           lo50=as.numeric(quantile(v,0.25)), hi50=as.numeric(quantile(v,0.75)))
  })
  df <- post |> left_join(params, by="variable") |>
    mutate(label=factor(label, levels=rev(params$label)),
           group=factor(group, levels=names(grp_cols)),
           sig=(lo89>0)|(hi89<0))
  
  # ── PRIOR SPECS — EXACTLY match the Stan model (fixes P3) ─────────────────
  # Previously the SDs here disagreed with the Stan priors (coy .30 vs .20,
  # dist_surv .25 vs .20, dist_detect .20 vs .15, jfm_adult/jfm_juv wrong).
  # Now sourced 1:1 from 05_ipm_model.R's model block.
  prior_specs <- tribble(
    ~variable,               ~prior_mean, ~prior_sd,
    "beta_moci_ond_fecund",      -0.15,      0.20,
    "beta_moci_ond_pup",         -0.15,      0.20,
    "beta_moci_jfm_pup",         -0.15,      0.20,
    "beta_moci_amj_pup",         -0.15,      0.20,
    "beta_moci_jfm_juv",         -0.15,      0.15,   # was 0.20 (P3)
    "beta_moci_jfm_adult",       -0.10,      0.12,   # was -0.15,0.20 (P3)
    "beta_moci_amj_molt",         0.05,      0.15,
    "detect_breed_logit",         1.20,      0.50,
    "detect_molt_logit",          0.75,      0.50,
    "beta_coy[1]",               -0.20,      0.20,   # was 0.30 (P3)
    "beta_coy[2]",               -0.20,      0.20,   # was 0.30 (P3)
    "beta_coy[3]",               -0.20,      0.20,   # was 0.30 (P3)
    "beta_dist_surv[1]",         -0.15,      0.20,   # was 0.25 (P3)
    "beta_dist_surv[2]",         -0.15,      0.20,
    "beta_dist_surv[3]",         -0.15,      0.20,
    "beta_dist_surv[4]",         -0.15,      0.20,
    "beta_dist_surv[5]",         -0.15,      0.20,
    "beta_dist_surv[6]",         -0.15,      0.20,
    "beta_eseal_pup",             0.10,      0.20)
  
  prior_df <- prior_specs |>
    mutate(prior_lo89=prior_mean+qnorm(CI_LO)*prior_sd,
           prior_hi89=prior_mean+qnorm(CI_HI)*prior_sd,
           prior_lo50=prior_mean+qnorm(0.25)*prior_sd,
           prior_hi50=prior_mean+qnorm(0.75)*prior_sd) |>
    left_join(params |> select(variable,label,group), by="variable") |>
    filter(!is.na(label)) |>
    mutate(label=factor(label, levels=rev(params$label)),
           group=factor(group, levels=names(grp_cols)))
  
  p <- ggplot(df, aes(y=label, colour=group)) +
    geom_linerange(data=prior_df, aes(y=label, xmin=prior_lo89, xmax=prior_hi89),
                   colour="grey72", linewidth=0.6, alpha=0.85,
                   position=position_nudge(y=0.25), inherit.aes=FALSE) +
    geom_linerange(data=prior_df, aes(y=label, xmin=prior_lo50, xmax=prior_hi50),
                   colour="grey58", linewidth=2.0, alpha=0.75,
                   position=position_nudge(y=0.25), inherit.aes=FALSE) +
    geom_point(data=prior_df, aes(y=label, x=prior_mean),
               colour="grey48", size=2.2, shape=1,
               position=position_nudge(y=0.25), inherit.aes=FALSE) +
    geom_linerange(aes(xmin=lo89, xmax=hi89), linewidth=0.7, alpha=0.7) +
    geom_linerange(aes(xmin=lo50, xmax=hi50), linewidth=2.2, alpha=0.9) +
    geom_point(aes(x=mean, shape=sig), size=3) +
    scale_shape_manual(values=c("FALSE"=16,"TRUE"=18),
                       labels=c("FALSE"="Not significant","TRUE"="Significant (89% CrI)"), name=NULL) +
    geom_vline(xintercept=0, linetype="dashed", colour="grey40", linewidth=0.5) +
    scale_colour_manual(values=grp_cols, name="Covariate group") +
    facet_grid(group ~ ., scales="free_y", space="free_y") +
    labs(x="Coefficient (logit scale)", y=NULL,
         title="Covariate Effects — Coefficient Forest Plot",
         subtitle=paste0("Posterior: thick=50% CrI, thin=", CI_LABEL,
                         ", diamond=excludes 0. Prior (grey): open circle=mean, bars=50%/89% prior")) +
    theme_seal(base_size=14) +
    theme(legend.position="bottom", strip.text.y=element_text(angle=0, face="bold", size=11),
          panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
          axis.text.y=element_text(size=10))
  if (save) ggsave(paste0("Output/Plots/",prefix,"_forest_plot.jpeg"),
                   p, width=28, height=26, units="cm", dpi=200)
  list(plot=p, data=df)
}

# ============================================================================
# PART 12: SUMMARY TABLE
# ============================================================================
create_summary_table_v3.2 <- function(fit, save=TRUE, prefix="IPM_v3.2") {
  params <- c("phi_pup_logit","phi_juv_base","phi_adult_F_logit","phi_adult_F_base","delta_adult",
              "fecund_primip","fecund_mature","prop_female","avg_fecundity","p_male_breed",
              "beta_coy[1]","beta_coy[2]","beta_coy[3]",
              paste0("beta_dist_surv[",1:6,"]"),
              "beta_moci_ond_fecund","beta_moci_amj_pup","beta_moci_jfm_juv",
              "beta_moci_jfm_adult","beta_eseal_pup","beta_moci_amj_molt",
              "sigma_process","sigma_obs_adult","sigma_obs_pup","sigma_obs_molt","sigma_site")
  pup_draws <- fit$draws(variables="phi_pup_logit",format="df")$phi_pup_logit
  pup_base <- tibble(variable="phi_pup_base (prob)",
                     mean=mean(plogis(pup_draws)), sd=sd(plogis(pup_draws)),
                     q_lo=quantile(plogis(pup_draws),CI_LO), q_hi=quantile(plogis(pup_draws),CI_HI),
                     rhat=NA_real_, ess_bulk=NA_real_)
  
  tbl <- seal_fit_summary(fit, params) |> bind_rows(pup_base) |>
    mutate(Estimate=sprintf("%.3f (%.3f, %.3f)", mean, q_lo, q_hi),
           Category=case_when(
             str_detect(variable,"phi|delta") ~ "Survival",
             str_detect(variable,"fecund|prop|avg") ~ "Reproduction",
             variable %in% c("p_male_breed","phi_pup_base (prob)") ~ "Observation/Derived",
             str_detect(variable,"beta_coy") ~ "Coyote",
             str_detect(variable,"beta_dist") ~ "Disturbance",
             str_detect(variable,"beta_moci|beta_eseal") ~ "Shared covariates",
             str_detect(variable,"sigma") ~ "Error terms", TRUE ~ "Other")) |>
    select(Category, variable, Estimate, rhat, ess_bulk)
  cat("\n=== PARAMETER SUMMARY — IPM v3.2 ===\n")
  print(tbl, n=nrow(tbl))
  if (save) write_csv(tbl, paste0("Output/",prefix,"_parameter_summary.csv"))
  tbl
}

# ============================================================================
# PART 13: KEY RESULTS PRINTOUT
# ============================================================================
save_model_output_v3.2 <- function(fit, prefix="IPM_v3.2") {
  fit$save_object(paste0("Output/harbor_seal_",prefix,"_fit.rds"))
  surv <- seal_fit_summary(fit, c("phi_juv_base","phi_adult_F_base","delta_adult"))
  pup_l <- fit$draws(variables="phi_pup_logit",format="df")$phi_pup_logit
  pmb  <- seal_fit_summary(fit, "p_male_breed")
  coy  <- seal_fit_summary(fit, paste0("beta_coy[",1:3,"]"))
  dst  <- seal_fit_summary(fit, paste0("beta_dist_surv[",1:6,"]"))
  sns  <- c("BL","DE","DP","PRH","TB","TP")
  cat("\n=== HARBOR SEAL IPM v3.2 — KEY RESULTS ===\n")
  cat(sprintf("Pup survival (prob): %.3f (%.3f-%.3f)\n",
              median(plogis(pup_l)), quantile(plogis(pup_l),CI_LO), quantile(plogis(pup_l),CI_HI)))
  cat(sprintf("Juvenile survival:   %.3f (%.3f-%.3f)\n", surv$mean[1],surv$q_lo[1],surv$q_hi[1]))
  cat(sprintf("Adult female surv:   %.3f (%.3f-%.3f)\n", surv$mean[2],surv$q_lo[2],surv$q_hi[2]))
  cat(sprintf("Adult male surv:     %.3f (%.3f-%.3f)\n",
              surv$mean[2]-surv$mean[3], surv$q_lo[2]-surv$q_hi[3], surv$q_hi[2]-surv$q_lo[3]))
  cat(sprintf("p_male_breed:        %.3f (%.3f-%.3f)\n", pmb$mean,pmb$q_lo,pmb$q_hi))
  cat("\nCoyote effects on pup survival:\n")
  for (i in 1:3) cat(sprintf("  %s: %.3f (%.3f-%.3f)\n",c("BL","DE","DP")[i],coy$mean[i],coy$q_lo[i],coy$q_hi[i]))
  cat("\nDisturbance effects on pup survival:\n")
  for (i in 1:6) cat(sprintf("  %s: %.3f (%.3f-%.3f)\n",sns[i],dst$mean[i],dst$q_lo[i],dst$q_hi[i]))
  tryCatch({
    dtb <- seal_fit_summary(fit, c("detect_breed_logit","detect_molt_logit"))
    cat("\nDetection baselines (logit -> p):\n")
    cat(sprintf("  Breeding: %.3f (%.3f-%.3f) -> p=%.2f\n", dtb$mean[1],dtb$q_lo[1],dtb$q_hi[1],plogis(dtb$mean[1])))
    cat(sprintf("  Molt:     %.3f (%.3f-%.3f) -> p=%.2f\n", dtb$mean[2],dtb$q_lo[2],dtb$q_hi[2],plogis(dtb$mean[2])))
  }, error=function(e) cat("  (detection baselines not in fit)\n"))
}

# ============================================================================
# PART 14: COVARIATE DECOMPOSITION (uses fit$summary posterior means)
# ============================================================================
create_covariate_decomposition_plots_v3.2 <- function(fit, sim_data, save=TRUE, prefix="IPM_v3.2") {
  years <- sim_data$years; site_names <- sim_data$site_names; S <- length(site_names)
  sdat <- if (!is.null(sim_data$stan_data$moci_amj)) sim_data$stan_data
  else if (!is.null(sim_data$moci_amj)) sim_data
  else stop("Decomp: cannot find moci_amj in sim_data.")
  bvars <- c(paste0("beta_coy[",1:3,"]"), paste0("beta_dist_surv[",1:S,"]"),
             "beta_eseal_pup","beta_moci_amj_pup","beta_moci_ond_fecund",
             "beta_moci_jfm_juv","beta_moci_jfm_adult")
  bsum <- fit$summary(variables=bvars)
  gm <- function(v){ r<-bsum$mean[bsum$variable==v]; if (length(r)==0L) 0 else r }
  gm_opt <- function(v) tryCatch({x<-fit$summary(variables=v)$mean; if (length(x)==0L) 0 else x}, error=function(e) 0)
  b_coy<-sapply(1:3,function(k)gm(paste0("beta_coy[",k,"]")))
  b_dst<-sapply(1:S,function(k)gm(paste0("beta_dist_surv[",k,"]")))
  b_es<-gm("beta_eseal_pup"); b_amj<-gm("beta_moci_amj_pup"); b_ond<-gm("beta_moci_ond_fecund")
  b_ond_pup<-gm_opt("beta_moci_ond_pup"); b_jfm_pup<-gm_opt("beta_moci_jfm_pup")
  coy_idx<-c(1,2,3,0,0,0); has_es<-c(0,1,0,1,0,0)
  
  contrib <- map_dfr(1:S, function(si) map_dfr(seq_along(years), function(t) {
    tb <- if (t>1L) t-1L else 1L
    coy <- if (coy_idx[si]>0) b_coy[coy_idx[si]]*sdat$coyote[si,tb] else 0
    dst <- b_dst[si]*sdat$disturbance[si,tb]
    es  <- if (has_es[si]>0) b_es*sdat$elephant_seal[si,tb] else 0
    tibble(Site=site_names[si], Year=years[t],
           `Coyote (t-1)`=coy, `Disturbance (t-1)`=dst, `Elephant seal (t-1)`=es,
           `MOCI AMJ -> pup (t-1)`=b_amj*sdat$moci_amj[tb],
           `MOCI OND -> pup (t)`=b_ond_pup*sdat$moci_ond[t],
           `MOCI JFM -> pup (t)`=b_jfm_pup*sdat$moci_jfm[t])
  })) |> mutate(Site=factor(Site, levels=site_names))
  
  cov_cols <- c(`Coyote (t-1)`="#B2182B", `Disturbance (t-1)`="#8C510A",
                `Elephant seal (t-1)`="#762A83", `MOCI AMJ -> pup (t-1)`="#4DAC26",
                `MOCI OND -> pup (t)`="#B8E186", `MOCI JFM -> pup (t)`="#74C476")
  pup_vars <- names(cov_cols)
  pup_long <- contrib |> select(Site,Year,all_of(pup_vars)) |>
    pivot_longer(cols=-all_of(c("Site","Year")), names_to="Covariate", values_to="Effect") |>
    filter(!is.na(Effect) & Effect!=0) |>
    mutate(Covariate=factor(Covariate, levels=pup_vars))
  
  p_pup <- ggplot(pup_long, aes(x=Year, y=Effect, fill=Covariate)) +
    geom_col(position="stack", width=0.75, alpha=0.85) +
    geom_hline(yintercept=0, linewidth=0.6, colour="gray30") +
    facet_wrap(~Site, ncol=2, scales="fixed") +
    scale_fill_manual(values=cov_cols[pup_vars]) +
    labs(x="Year", y="Covariate effect (logit scale)",
         title="Covariate Contributions to Pup Survival by Site and Year",
         subtitle="Birth-year covariates at t-1; positive = higher pup survival") +
    theme_seal() +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=9), legend.position="bottom",
          strip.text=element_text(size=11,face="bold"), panel.spacing=unit(0.5,"lines"))
  if (save) ggsave(paste0("Output/Plots/",prefix,"_decomp_pup_survival.jpeg"),
                   p_pup, width=32, height=30, units="cm", dpi=200)
  list(pup_survival=p_pup, data=contrib)
}

# ============================================================================
# PART 15: PORTFOLIO ANALYSIS
# ----------------------------------------------------------------------------
# NOTE[P4]: PER below is computed from posterior-MEAN lambda (point estimate,
# no CrI). The Methods describe a per-draw computation; reconcile, or extend
# to per-draw to report a PER credible interval. Left as-is here.
# ============================================================================
create_portfolio_analysis_v3.2 <- function(fit, sim_data, save=TRUE, prefix="IPM_v3.2") {
  years <- sim_data$years; site_names <- sim_data$site_names
  S <- length(site_names); T <- length(years)
  ldraws <- fit$draws(variables="lambda", format="matrix")
  lmat <- matrix(NA,S,T-1, dimnames=list(site_names,years[1:(T-1)]))
  for (s in 1:S) for (t in 1:(T-1)) {
    cn <- paste0("lambda[",s,",",t,"]"); if (cn %in% colnames(ldraws)) lmat[s,t] <- mean(ldraws[,cn])
  }
  cv_meta  <- sd(colMeans(lmat,na.rm=TRUE))/mean(colMeans(lmat,na.rm=TRUE))
  cv_sites <- mean(apply(lmat,1,function(x) sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)))
  per <- cv_meta/cv_sites
  lcor <- cor(t(lmat), use="pairwise.complete.obs")
  cat(sprintf("Portfolio Effect Ratio: %.3f (<1 = buffering)  [P4: point est, no CrI]\n", per))
  cat(sprintf("Mean site correlation:  %.3f\n", mean(lcor[lower.tri(lcor)])))
  
  ldf <- expand.grid(Site=site_names, Year=years[1:(T-1)]) |>
    mutate(Site=factor(Site,levels=site_names), lambda=as.vector(t(lmat)))
  p_heat <- ggplot(ldf,aes(x=Year,y=Site,fill=lambda)) +
    geom_tile(color="white",linewidth=0.5) +
    scale_fill_gradient2(low="red3",mid="white",high="darkgreen",
                         midpoint=1,limits=c(0.7,1.3),oob=scales::squish,name="lambda") +
    geom_text(aes(label=sprintf("%.2f",lambda)),size=2.5) +
    labs(x="Year",y="Site",title="Site-Specific lambda by Year") +
    theme_seal() + theme(axis.text.x=element_text(angle=45,hjust=1,size=8), panel.grid=element_blank())
  if (save) ggsave(paste0("Output/Plots/",prefix,"_portfolio_lambda_heatmap.jpeg"),
                   p_heat, width=35, height=15, units="cm")
  
  cordf <- expand.grid(Site1=site_names,Site2=site_names) |> mutate(r=as.vector(lcor))
  p_cor <- ggplot(cordf,aes(x=Site1,y=Site2,fill=r)) +
    geom_tile(color="white") + geom_text(aes(label=sprintf("%.2f",r)),size=4) +
    scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0,limits=c(-1,1)) +
    labs(x="",y="",title="Between-Site Correlation in lambda") +
    theme_seal() + theme(panel.grid=element_blank()) + coord_fixed()
  if (save) ggsave(paste0("Output/Plots/",prefix,"_portfolio_correlation.jpeg"),
                   p_cor, width=20, height=18, units="cm")
  
  port_summary_tbl <- tibble(
    Metric = c("Portfolio Effect Ratio (CVmeta/CVsites)",
               "Mean between-site corr (lambda)",
               "Mean per-site lambda (across sites/years)"),   # was "Overall mean lambda"
    Value  = c(round(per,3), round(mean(lcor[lower.tri(lcor)]),3), round(mean(lmat,na.rm=TRUE),3)))
  cat("\n── Portfolio Summary ──\n"); print(port_summary_tbl)
  if (save) write_csv(port_summary_tbl, paste0("Output/",prefix,"_portfolio_summary.csv"))
  
  list(lambda_heatmap=p_heat, correlation=p_cor, lambda_matrix=lmat,
       summary=list(portfolio_effect_ratio=per, mean_correlation=mean(lcor[lower.tri(lcor)])),
       summary_table=port_summary_tbl)
}

# ============================================================================
# PART 16: SYNCHRONY PROJECTIONS
# ============================================================================
create_synchrony_projections_v3.2 <- function(fit, sim_data, n_sims=2000, T_proj=10,
                                              save=TRUE, prefix="IPM_v3.2", seed=2026) {
  set.seed(seed)   
  years <- sim_data$years; site_names <- sim_data$site_names
  S <- length(site_names); T <- length(years); pyrs <- (max(years)+1):(max(years)+T_proj)
  draws <- fit$draws(format="df")
  idx <- sample(seq_len(nrow(draws)), min(n_sims,nrow(draws)))
  
  pup_l<-draws$phi_pup_logit[idx]; phi_juv<-draws$phi_juv_base[idx]
  phi_aF_logit<-draws$phi_adult_F_logit[idx]; delta_a<-draws$delta_adult[idx]
  pf<-draws$prop_female[idx]; avgf<-draws$avg_fecundity[idx]
  b_ond<-draws$beta_moci_ond_fecund[idx]; b_amj<-draws$beta_moci_amj_pup[idx]
  b_jfmJ<-draws$beta_moci_jfm_juv[idx]; b_jfmA<-draws$beta_moci_jfm_adult[idx]
  b_ond_pup<-tryCatch(draws$beta_moci_ond_pup[idx], error=function(e) rep(0,length(idx)))
  b_jfm_pup<-tryCatch(draws$beta_moci_jfm_pup[idx], error=function(e) rep(0,length(idx)))
  se<-sapply(1:S,function(s) draws[[paste0("site_effect[",s,"]")]][idx])
  bcoy<-sapply(1:3,function(k) draws[[paste0("beta_coy[",k,"]")]][idx]); cidx<-c(1,2,3,0,0,0)
  NAF<-sapply(1:S,function(s) draws[[paste0("N_adult_F[",s,",",T,"]")]][idx])
  NAM<-sapply(1:S,function(s) draws[[paste0("N_adult_M[",s,",",T,"]")]][idx])
  NJF<-sapply(1:S,function(s) draws[[paste0("N_juv_F[",s,",",T,"]")]][idx])
  NJM<-sapply(1:S,function(s) draws[[paste0("N_juv_M[",s,",",T,"]")]][idx])
  NP <-sapply(1:S,function(s) draws[[paste0("N_pup[",s,",",T,"]")]][idx])
  
  scenarios <- list(list(name="Status Quo",moci=0,coyote=0),
                    list(name="Cool (MOCI -1)",moci=-1,coyote=0),
                    list(name="Warm (MOCI +1)",moci=1,coyote=0),
                    list(name="Warm + High Coyote",moci=1,coyote=1))
  run_proj <- function(sync, sc, psd=0.15) {
    Nt <- matrix(NA, length(idx), T_proj)
    for (i in seq_along(idx)) {
      naf<-NAF[i,];nam<-NAM[i,];njf<-NJF[i,];njm<-NJM[i,];np<-NP[i,]
      Nt[i,1] <- sum(naf+nam+njf+njm+np)
      for (tp in 2:T_proj) {
        sh <- if (sync) rep(rnorm(1,0,psd),S) else rnorm(S,0,psd)
        for (s in 1:S) {
          ce <- if (cidx[s]>0) bcoy[i,cidx[s]]*sc$coyote else 0
          fecund_t <- plogis(qlogis(pmax(pmin(avgf[i],0.999),0.001)) + b_ond[i]*sc$moci)
          pp <- plogis(pup_l[i]+se[i,s]+ce + b_amj[i]*sc$moci + b_ond_pup[i]*sc$moci + b_jfm_pup[i]*sc$moci)
          pj <- plogis(qlogis(phi_juv[i])+se[i,s]*0.5+b_jfmJ[i]*sc$moci)
          paF <- plogis(phi_aF_logit[i]+se[i,s]*0.25+b_jfmA[i]*sc$moci)
          paM <- plogis(qlogis(pmax(plogis(phi_aF_logit[i])-delta_a[i],0.001))+se[i,s]*0.25+b_jfmA[i]*sc$moci)
          new_p<-naf[s]*fecund_t*exp(sh[s])
          njF2<-np[s]*pf[i]*pp; njM2<-np[s]*(1-pf[i])*pp
          jsF<-njf[s]*pj*(2/3); jsM<-njm[s]*pj*(2/3); jaF<-njf[s]*pj*(1/3); jaM<-njm[s]*pj*(1/3)
          np[s]<-max(new_p,1); njf[s]<-max(njF2+jsF,.1); njm[s]<-max(njM2+jsM,.1)
          naf[s]<-max(naf[s]*paF+jaF,1); nam[s]<-max(nam[s]*paM+jaM,1)
        }
        Nt[i,tp] <- sum(naf+nam+njf+njm+np)
      }
    }
    Nt
  }
  res <- lapply(scenarios, function(sc) {
    cat(sprintf("Projecting: %s\n",sc$name))
    list(async=run_proj(FALSE,sc), sync=run_proj(TRUE,sc))
  })
  names(res) <- sapply(scenarios,`[[`,"name")
  
  cv_df <- map_dfr(names(res), function(sn) {
    ca <- sd(res[[sn]]$async[,T_proj])/mean(res[[sn]]$async[,T_proj])
    cs <- sd(res[[sn]]$sync[,T_proj]) /mean(res[[sn]]$sync[,T_proj])
    tibble(Scenario=sn, CV_Async=ca, CV_Sync=cs, CV_Ratio=ca/cs, Buffering_Pct=(1-ca/cs)*100)
  })
  cat("\n--- PORTFOLIO BUFFERING ---\n"); print(cv_df)
  if (save) write_csv(cv_df, paste0("Output/",prefix,"_synchrony_cv_comparison.csv"))
  
  comp_df <- map_dfr(names(res), function(sn) bind_rows(
    tibble(Scenario=sn,Year=pyrs,Synchrony="Asynchronous (current)",
           mean=colMeans(res[[sn]]$async),
           lo=as.numeric(apply(res[[sn]]$async,2,quantile,CI_LO)),
           hi=as.numeric(apply(res[[sn]]$async,2,quantile,CI_HI))),
    tibble(Scenario=sn,Year=pyrs,Synchrony="Synchronous (hypothetical)",
           mean=colMeans(res[[sn]]$sync),
           lo=as.numeric(apply(res[[sn]]$sync,2,quantile,CI_LO)),
           hi=as.numeric(apply(res[[sn]]$sync,2,quantile,CI_HI)))))
  p_comp <- ggplot(comp_df,aes(x=Year,y=mean,color=Synchrony,fill=Synchrony,linetype=Synchrony)) +
    geom_ribbon(aes(ymin=lo,ymax=hi),alpha=0.15,color=NA) + geom_line(linewidth=1.2) +
    facet_wrap(~Scenario,ncol=2,scales="free_y") +
    scale_color_manual(values=c("Asynchronous (current)"=SEAL_COLS$pop,"Synchronous (hypothetical)"=SEAL_COLS$adult_f)) +
    scale_fill_manual(values=c("Asynchronous (current)"=SEAL_COLS$ribbon,"Synchronous (hypothetical)"="#FFCDD2")) +
    scale_linetype_manual(values=c("Asynchronous (current)"="solid","Synchronous (hypothetical)"="dashed")) +
    labs(x="Year",y="Total Population",title="Portfolio Buffering: Async vs Sync") + theme_seal()
  if (save) ggsave(paste0("Output/Plots/",prefix,"_synchrony_comparison.jpeg"),
                   p_comp, width=30, height=25, units="cm")
  list(comparison=p_comp, cv_comparison=cv_df, raw_projections=res)
}

# ============================================================================
# PART 17: AGGREGATE LAMBDA + ELASTICITY
# ----------------------------------------------------------------------------
# Fills the §3.5/§3.6 placeholders. Two distinct quantities, clearly separated:
#   (1) REALIZED aggregate lambda from N_total_all (time-varying covariates).
#       Reported as mean +/- interannual SD and min/max (NOT a CrI), plus a
#       genuine posterior CrI on the overall-mean lambda (per-draw mean-over-yrs).
#   (2) ASYMPTOTIC deterministic lambda + elasticities from the posterior-MEDIAN
#       female Lefkovitch sub-matrix (Eq. 1). prop_female is taken from the fit
#       (not hardcoded 0.5). Elasticities are scale/sign invariant, so the
#       per-draw CrI loop is numerically robust. Summed elasticities ~ 1.0.
# NOTE: lambda_det (asymptotic) differs from realized lambda by design — label
#       both clearly in the manuscript so they don't read as contradictory.
# ============================================================================
create_lambda_elasticity_v3.2 <- function(fit, save=TRUE, prefix="IPM_v3.2") {
  draws <- fit$draws(format="df")
  
  # ── (1) Realized aggregate lambda from N_total_all (in generated quantities) ─
  Ntot <- as.matrix(draws[, grep("^N_total_all\\[", names(draws)), drop=FALSE])
  T <- ncol(Ntot)
  lambda_draws <- Ntot[, 2:T, drop=FALSE] / Ntot[, 1:(T-1), drop=FALSE]  # draws x (T-1)
  lam_med  <- apply(lambda_draws, 2, median)              # annual point estimates
  
  lam_mean <- mean(lam_med); lam_sd <- sd(lam_med)
  lam_min  <- min(lam_med);  lam_max <- max(lam_med)
  yrs_above <- sum(lam_med > 1); yrs_below <- sum(lam_med < 1)
  
  # Genuine posterior CrI on the OVERALL-MEAN lambda (mean over years per draw)
  mean_lambda_per_draw <- rowMeans(lambda_draws)
  cri_lo <- as.numeric(quantile(mean_lambda_per_draw, CI_LO))
  cri_hi <- as.numeric(quantile(mean_lambda_per_draw, CI_HI))
  
  cat("\n── Aggregate lambda ──\n")
  cat(sprintf("Mean annual lambda: %.3f (interannual SD %.3f; range %.3f-%.3f)\n",
              lam_mean, lam_sd, lam_min, lam_max))
  cat(sprintf("Overall-mean lambda %s: %.3f (%.3f-%.3f)\n",
              CI_LABEL, mean(mean_lambda_per_draw), cri_lo, cri_hi))
  cat(sprintf("Years above 1.0: %d / %d   below 1.0: %d / %d\n",
              yrs_above, T-1, yrs_below, T-1))
  
  lambda_summary <- tibble(
    metric = c("Mean annual lambda","Interannual SD","Min annual","Max annual",
               "Overall-mean lambda","Overall-mean CrI lo","Overall-mean CrI hi",
               "Years above 1.0","Years below 1.0","N transitions"),
    value  = c(lam_mean, lam_sd, lam_min, lam_max,
               mean(mean_lambda_per_draw), cri_lo, cri_hi,
               yrs_above, yrs_below, T-1))
  
  # ── (2) Asymptotic deterministic lambda + elasticity (median vital rates) ───
  # Elasticity is invariant to eigenvector scaling/sign, so eigen() ambiguity
  # does not affect results (the per-draw loop relies on this).
  elas_of <- function(M) {
    ev <- eigen(M); k <- which.max(Re(ev$values)); lam <- Re(ev$values[k])
    w  <- Re(ev$vectors[,k])
    vt <- Re(eigen(t(M))$vectors[,k])
    list(lambda=lam, elas=(outer(vt, w)/sum(vt*w)) * (M/lam))
  }
  build_A <- function(phi_pup, phi_juv, phi_aF, f_avg, pf)
    matrix(c(0,            0,             f_avg,
             phi_pup * pf, (2/3)*phi_juv, 0,
             0,            (1/3)*phi_juv, phi_aF),
           nrow=3, byrow=TRUE,
           dimnames=list(c("pup","juv_F","adult_F"),
                         c("pup","juv_F","adult_F")))
  
  A   <- build_A(median(draws$phi_pup_base), median(draws$phi_juv_base),
                 median(draws$phi_adult_F_base), median(draws$avg_fecundity),
                 median(draws$prop_female))
  e0  <- elas_of(A); lambda_det <- e0$lambda; elas <- e0$elas
  eps <- c(adult_survival = elas["adult_F","adult_F"],
           juv_survival   = elas["adult_F","juv_F"] + elas["juv_F","juv_F"],
           pup_survival   = elas["juv_F","pup"],
           fecundity      = elas["pup","adult_F"])
  eps_sum <- sum(eps)
  
  cat("\n── Asymptotic (female-only) lambda + elasticity ──\n")
  cat(sprintf("lambda_det = %.4f   [asymptotic; realized mean = %.3f]\n",
              lambda_det, lam_mean))
  for (nm in names(eps)) cat(sprintf("  e_%-14s %.3f\n", nm, eps[nm]))
  cat(sprintf("  (sum = %.3f; should be ~1.0)\n", eps_sum))
  
  # Per-draw elasticity CrIs (ranking uncertainty)
  nd  <- min(1000L, nrow(draws)); idx <- sample(seq_len(nrow(draws)), nd)
  eps_draws <- t(vapply(idx, function(i) {
    el <- elas_of(build_A(draws$phi_pup_base[i], draws$phi_juv_base[i],
                          draws$phi_adult_F_base[i], draws$avg_fecundity[i],
                          draws$prop_female[i]))$elas
    c(adult_survival=el["adult_F","adult_F"],
      juv_survival  =el["adult_F","juv_F"]+el["juv_F","juv_F"],
      pup_survival  =el["juv_F","pup"], fecundity=el["pup","adult_F"])
  }, numeric(4)))
  
  elas_tbl <- tibble(
    vital_rate = names(eps), elasticity = as.numeric(eps),
    cri_lo = apply(eps_draws, 2, quantile, CI_LO),
    cri_hi = apply(eps_draws, 2, quantile, CI_HI)) |>
    arrange(desc(elasticity))
  cat("\n── Elasticity table (median + per-draw CrI) ──\n"); print(elas_tbl)
  
  if (save) {
    write_csv(lambda_summary, paste0("Output/",prefix,"_lambda_summary.csv"))
    write_csv(elas_tbl,       paste0("Output/",prefix,"_elasticity_table.csv"))
  }
  
  # ── Plots ────────────────────────────────────────────────────────────────
  p_elas <- ggplot(elas_tbl, aes(x=reorder(vital_rate, elasticity), y=elasticity)) +
    geom_col(fill=SEAL_COLS$pop, alpha=0.85) +
    geom_errorbar(aes(ymin=cri_lo, ymax=cri_hi), width=0.2, colour="grey30") +
    geom_text(aes(label=sprintf("%.3f", elasticity)), hjust=-0.3, size=4) +
    coord_flip() + expand_limits(y=max(elas_tbl$cri_hi)*1.15) +
    labs(x=NULL, y="Elasticity of asymptotic lambda",
         title="Elasticity of lambda to Vital Rates",
         subtitle=sprintf("Asymptotic lambda = %.3f; bars=median, whiskers=%s; sum=%.2f",
                          lambda_det, CI_LABEL, eps_sum)) +
    theme_seal()
  if (save) ggsave(paste0("Output/Plots/",prefix,"_elasticity.jpeg"),
                   p_elas, width=22, height=14, units="cm")
  
  lam_df <- tibble(t=1:(T-1), median=lam_med,
                   lo=apply(lambda_draws,2,quantile,CI_LO),
                   hi=apply(lambda_draws,2,quantile,CI_HI))
  p_lam <- ggplot(lam_df, aes(x=t, y=median)) +
    geom_hline(yintercept=1, linetype=2, colour="grey50") +
    geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.2, fill=SEAL_COLS$ribbon) +
    geom_line(linewidth=1.1, colour=SEAL_COLS$pop) +
    labs(x="Transition (year index)", y=expression(lambda),
         title="Annual Aggregate Population Growth Rate",
         subtitle=sprintf("Mean %.3f; %d/%d transitions > 1.0; bands = %s",
                          lam_mean, yrs_above, T-1, CI_LABEL)) +
    theme_seal()
  if (save) ggsave(paste0("Output/Plots/",prefix,"_lambda_annual.jpeg"),
                   p_lam, width=24, height=12, units="cm")
  
  list(lambda_summary=lambda_summary, elasticity_table=elas_tbl,
       lambda_det=lambda_det, A=A, elas_matrix=elas,
       elasticity_plot=p_elas, lambda_plot=p_lam)
}

# ============================================================================
# PART 18: MOCI COLLINEARITY DIAGNOSTICS
# ----------------------------------------------------------------------------
# Supplies the §4.3 numbers: (1) empirical correlation among the MOCI predictor
# series as they enter the model, (2) posterior correlation among the MOCI beta
# parameters, (3) the aggregate "any-pathway" pup-survival MOCI effect (sum of
# the 3 seasonal pup betas) — which is better identified than any single term.
# NOTE: real-data moci_ond is ALREADY pre-lagged in 02_covariates (so ond_t =
# previous fall). ond_tm1 below is shown for reference but is NOT a model term.
# ============================================================================
create_moci_collinearity_v3.2 <- function(fit, sim_data, save=TRUE, prefix="IPM_v3.2") {
  sdat <- sim_data$stan_data
  if (is.null(sdat$moci_jfm)) {
    auto <- paste0("Output/harbor_seal_", prefix, "_input_data.rds")
    if (file.exists(auto)) sdat <- readRDS(auto)$stan_data
  }
  if (is.null(sdat$moci_jfm)) stop("MOCI series not found in sim_data or input RDS.")
  Tn <- length(sdat$moci_jfm)
  
  # ── (1) Empirical correlation among MOCI predictor series ─────────────────
  moci_df <- data.frame(
    jfm_t   = sdat$moci_jfm,                    # JFM(t) -> pup, juv, adult survival
    amj_t   = sdat$moci_amj,                    # AMJ(t) -> molt detection
    amj_tm1 = c(NA, sdat$moci_amj[1:(Tn-1)]),   # AMJ(t-1) -> pup (birth year)
    ond_t   = sdat$moci_ond,                    # OND(t) -> fecundity + pup (pre-lagged)
    ond_tm1 = c(NA, sdat$moci_ond[1:(Tn-1)]))   # reference only; not a model term
  emp_cor <- cor(moci_df, use="pairwise.complete.obs")
  cat("\n── Empirical correlation among MOCI predictor series ──\n")
  print(round(emp_cor, 3))
  cat(sprintf("JFM(t) vs OND(t): r = %.3f  (the value cited in §4.3)\n",
              emp_cor["jfm_t","ond_t"]))
  
  # ── (2) Posterior correlation among MOCI beta parameters ──────────────────
  beta_vars <- c("beta_moci_jfm_pup","beta_moci_amj_pup","beta_moci_ond_pup",
                 "beta_moci_jfm_juv","beta_moci_jfm_adult",
                 "beta_moci_ond_fecund","beta_moci_amj_molt")
  beta_draws <- fit$draws(variables=beta_vars, format="df") |>
    dplyr::select(dplyr::all_of(beta_vars))
  post_cor <- cor(beta_draws)
  cat("\n── Posterior correlation among MOCI beta parameters ──\n")
  print(round(post_cor, 3))
  pup_betas <- c("beta_moci_jfm_pup","beta_moci_amj_pup","beta_moci_ond_pup")
  pc <- post_cor[pup_betas, pup_betas]
  cat(sprintf("Pairwise posterior corr among 3 pup MOCI betas: %.3f to %.3f\n",
              min(pc[lower.tri(pc)]), max(pc[lower.tri(pc)])))
  
  # ── (3) Aggregate pup-survival MOCI effect (sum of 3 seasonal terms) ──────
  pup_sum <- beta_draws$beta_moci_jfm_pup +
    beta_draws$beta_moci_amj_pup +
    beta_draws$beta_moci_ond_pup
  agg <- tibble(quantity="Aggregate pup-survival MOCI (sum of 3 seasonal betas)",
                median=median(pup_sum),
                cri_lo=as.numeric(quantile(pup_sum, CI_LO)),
                cri_hi=as.numeric(quantile(pup_sum, CI_HI)),
                excludes_zero=(quantile(pup_sum,CI_LO)>0)|(quantile(pup_sum,CI_HI)<0))
  cat(sprintf("\n── Aggregate pup-survival MOCI effect ──\nMedian %.3f, %s [%.3f, %.3f]  (excludes 0: %s)\n",
              agg$median, CI_LABEL, agg$cri_lo, agg$cri_hi, agg$excludes_zero))
  
  if (save) {
    write_csv(as.data.frame(round(emp_cor,3)) |> tibble::rownames_to_column("series"),
              paste0("Output/",prefix,"_moci_empirical_corr.csv"))
    write_csv(as.data.frame(round(post_cor,3)) |> tibble::rownames_to_column("beta"),
              paste0("Output/",prefix,"_moci_posterior_corr.csv"))
    write_csv(agg, paste0("Output/",prefix,"_moci_aggregate_pup.csv"))
  }
  
  list(empirical_corr=emp_cor, posterior_corr=post_cor,
       aggregate_pup=agg, pup_sum_draws=pup_sum)
}

cat("\n06_ipm_plots.R loaded — plotting functions ready.\n")
cat("Run via 07_ipm_run.R, or: run_all_plots_v3.2(out$fit, out$sim_data, prefix='IPM_v3.2_real')\n")