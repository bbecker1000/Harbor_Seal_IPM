# ============================================================================
# HARBOR SEAL IPM v3.2 — PORTFOLIO PROJECTION ANALYSIS
# 10-year projections under varying MOCI, disturbance, and coyote
# (a) By site    (b) Total population
# Run after model has completed.
# ============================================================================
#
# SETUP — choose one:
#
#   Option A: results.real exists in session
#     fit      <- Results.real$fit
#     sim_data <- Results.real$data
#
#   Option B: fresh session, load from saved RDS
#     source("harbor_seal_ipm_v3.2.R")
#     source("Code/harbor_seal_ipm_v3.2_plots.R")
#     out      <- load_seal_results("IPM_v3.2_real")
#     fit      <- out$fit
#     sim_data <- out$sim_data
#
# ============================================================================

library(tidyverse)
library(patchwork)

dir.create("Output",       showWarnings=FALSE)
dir.create("Output/Plots", showWarnings=FALSE)

prefix     <- "IPM_v3.2_real"
site_names <- c("BL","DE","DP","PRH","TB","TP")
S          <- length(site_names)
T          <- length(sim_data$years)
T_proj     <- 10
proj_years <- (max(sim_data$years)+1):(max(sim_data$years)+T_proj)
coyote_idx <- c(1,2,3,0,0,0)   # BL=1, DE=2, DP=3, others=0 (no coyote)
n_sims     <- 500               # posterior draws to use
psd        <- 0.15              # process SD for stochastic projection

# ── Constants (sourced from plots script; define fallbacks if not loaded) ────
if (!exists("CI_LO"))    CI_LO    <- 0.055
if (!exists("CI_HI"))    CI_HI    <- 0.945
if (!exists("CI_LABEL")) CI_LABEL <- "89% CrI"
if (!exists("SEAL_COLS")) SEAL_COLS <- list(
  pop="navy", pup="#1B7837", juv="#762A83",
  adult_f="#B2182B", ribbon="#AECDE8", neutral="gray60"
)
if (!exists("theme_seal")) theme_seal <- function(base_size=16) {
  theme_minimal(base_size=base_size) %+replace%
    theme(panel.border=element_rect(colour="grey70",fill=NA,linewidth=0.5),
          strip.text=element_text(face="bold"),
          plot.title=element_text(face="bold"),
          legend.position="bottom")
}

# ============================================================================
# SCENARIO DEFINITIONS
# Each scenario specifies standardised covariate shifts (in SD units):
#   moci_ond: OND MOCI → fecundity (maternal condition)
#   moci_amj: AMJ MOCI → pup survival
#   moci_jfm: JFM MOCI → juvenile and adult survival
#   dist:     disturbance → pup survival (all sites scaled equally)
#   coyote:   coyote pressure → pup survival (BL, DE, DP only)
# ============================================================================

scenarios <- list(
  list(name="Status Quo",
       moci_ond=0,  moci_amj=0,  moci_jfm=0,  dist=0,  coyote=0),
  list(name="Warm MOCI (+1 SD)",
       moci_ond=1,  moci_amj=1,  moci_jfm=1,  dist=0,  coyote=0),
  list(name="Cool MOCI (-1 SD)",
       moci_ond=-1, moci_amj=-1, moci_jfm=-1, dist=0,  coyote=0),
  list(name="High Disturbance (+1 SD)",
       moci_ond=0,  moci_amj=0,  moci_jfm=0,  dist=1,  coyote=0),
  list(name="High Coyote (+1 SD)",
       moci_ond=0,  moci_amj=0,  moci_jfm=0,  dist=0,  coyote=1),
  list(name="Warm + High Stress",
       moci_ond=1,  moci_amj=1,  moci_jfm=1,  dist=1,  coyote=1),
  list(name="Cool + Low Stress",
       moci_ond=-1, moci_amj=-1, moci_jfm=-1, dist=-1, coyote=-1)
)

scen_names <- sapply(scenarios, `[[`, "name")

# ============================================================================
# EXTRACT POSTERIOR DRAWS
# ============================================================================
cat("Extracting posterior draws...\n")

draws_df <- fit$draws(format="df")
set.seed(42)
idx <- sample(seq_len(nrow(draws_df)), min(n_sims, nrow(draws_df)))

# Survival parameters
pup_l    <- draws_df$phi_pup_logit[idx]
phi_juv  <- draws_df$phi_juv_base[idx]
phi_aF_l <- draws_df$phi_adult_F_logit[idx]
delta_a  <- draws_df$delta_adult[idx]
pf       <- draws_df$prop_female[idx]
avgf     <- draws_df$avg_fecundity[idx]

# MOCI coefficients
b_ond  <- draws_df$beta_moci_ond_fecund[idx]   # OND → fecundity
b_amj  <- draws_df$beta_moci_amj_pup[idx]      # AMJ → pup survival
b_jfmJ <- draws_df$beta_moci_jfm_juv[idx]      # JFM → juvenile survival
b_jfmA <- draws_df$beta_moci_jfm_adult[idx]    # JFM → adult survival

# Site-specific coefficients (n_sims × S or n_sims × 3)
se   <- sapply(1:S, function(s)
  draws_df[[paste0("site_effect[",s,"]")]][idx])
bcoy <- sapply(1:3, function(k)
  draws_df[[paste0("beta_coy[",k,"]")]][idx])
bdst <- sapply(1:S, function(s)
  draws_df[[paste0("beta_dist_surv[",s,"]")]][idx])

# Starting populations (end of observed time series)
NAF <- sapply(1:S, function(s) draws_df[[paste0("N_adult_F[",s,",",T,"]")]][idx])
NAM <- sapply(1:S, function(s) draws_df[[paste0("N_adult_M[",s,",",T,"]")]][idx])
NJF <- sapply(1:S, function(s) draws_df[[paste0("N_juv_F[",s,",",T,"]")]][idx])
NJM <- sapply(1:S, function(s) draws_df[[paste0("N_juv_M[",s,",",T,"]")]][idx])
NP  <- sapply(1:S, function(s) draws_df[[paste0("N_pup[",s,",",T,"]")]][idx])

cat(sprintf("  %d draws extracted from %d total\n", length(idx), nrow(draws_df)))

# ============================================================================
# PROJECTION ENGINE
# sync=TRUE:  all sites share the same annual process shock (hypothetical)
# sync=FALSE: each site gets an independent shock (current spatial structure)
# Returns array [n_sims × T_proj × S] (by-site) and matrix [n_sims × T_proj] (total)
# ============================================================================

run_projection <- function(sc, sync=FALSE) {
  N_site  <- array(NA_real_, dim=c(length(idx), T_proj, S))
  N_total <- matrix(NA_real_, nrow=length(idx), ncol=T_proj)
  
  for (i in seq_along(idx)) {
    naf <- NAF[i,]; nam <- NAM[i,]
    njf <- NJF[i,]; njm <- NJM[i,]
    np  <- NP[i,]
    
    # Year 1: starting state (last observed year)
    N_site[i,1,]  <- naf + nam + njf + njm + np
    N_total[i,1]  <- sum(N_site[i,1,])
    
    for (tp in 2:T_proj) {
      # Shared vs independent process shocks
      sh <- if (sync) rep(rnorm(1, 0, psd), S) else rnorm(S, 0, psd)
      
      for (s in 1:S) {
        # Coyote effect on pup survival (uses current projection year)
        ce  <- if (coyote_idx[s] > 0) bcoy[i, coyote_idx[s]] * sc$coyote else 0
        
        # Disturbance effect on pup survival
        dse <- bdst[i, s] * sc$dist
        
        # Vital rates under scenario
        # Fecundity: OND MOCI acts on maternal condition
        fecund_t <- plogis(qlogis(pmax(pmin(avgf[i], 0.999), 0.001)) +
                             b_ond[i] * sc$moci_ond)
        
        # Pup survival: AMJ MOCI + coyote + disturbance + site effect
        pp <- plogis(pup_l[i] + se[i,s] + ce + dse +
                       b_amj[i] * sc$moci_amj)
        
        # Juvenile survival: JFM MOCI + half site effect
        pj <- plogis(qlogis(phi_juv[i]) + se[i,s]*0.5 +
                       b_jfmJ[i] * sc$moci_jfm)
        
        # Adult survival: JFM MOCI + quarter site effect (sex-specific)
        paF <- plogis(phi_aF_l[i] + se[i,s]*0.25 + b_jfmA[i] * sc$moci_jfm)
        paM <- plogis(qlogis(pmax(plogis(phi_aF_l[i]) - delta_a[i], 0.001)) +
                        se[i,s]*0.25 + b_jfmA[i] * sc$moci_jfm)
        
        # Population transitions
        new_p <- naf[s] * fecund_t * exp(sh[s])        # pups produced
        njF2  <- np[s]  * pf[i] * pp                   # pup → juv female
        njM2  <- np[s]  * (1-pf[i]) * pp               # pup → juv male
        jsF   <- njf[s] * pj * (2/3)                   # juv stay female
        jsM   <- njm[s] * pj * (2/3)                   # juv stay male
        jaF   <- njf[s] * pj * (1/3)                   # juv → adult female
        jaM   <- njm[s] * pj * (1/3)                   # juv → adult male
        
        np[s]  <- max(new_p,    1)
        njf[s] <- max(njF2+jsF, 0.1)
        njm[s] <- max(njM2+jsM, 0.1)
        naf[s] <- max(naf[s]*paF + jaF, 1)
        nam[s] <- max(nam[s]*paM + jaM, 1)
      }
      
      N_site[i,tp,]  <- naf + nam + njf + njm + np
      N_total[i,tp]  <- sum(N_site[i,tp,])
    }
  }
  
  list(by_site=N_site, total=N_total)
}

# ============================================================================
# RUN ALL SCENARIOS (async and sync for portfolio buffering)
# ============================================================================
cat("Running projections across", length(scenarios), "scenarios...\n")

proj_results <- lapply(scenarios, function(sc) {
  cat(sprintf("  %s\n", sc$name))
  list(
    async = run_projection(sc, sync=FALSE),
    sync  = run_projection(sc, sync=TRUE)
  )
})
names(proj_results) <- scen_names

# ============================================================================
# PORTFOLIO STATISTICS PER SCENARIO
# ============================================================================

port_stats <- map_dfr(scen_names, function(sn) {
  pr <- proj_results[[sn]]
  
  # Total population CV at end of projection
  cv_async <- sd(pr$async$total[,T_proj]) / mean(pr$async$total[,T_proj])
  cv_sync  <- sd(pr$sync$total[,T_proj])  / mean(pr$sync$total[,T_proj])
  buffering <- (1 - cv_async/cv_sync) * 100
  
  # Mean final population
  n_final_mean <- mean(pr$async$total[,T_proj])
  n_final_lo   <- quantile(pr$async$total[,T_proj], CI_LO)
  n_final_hi   <- quantile(pr$async$total[,T_proj], CI_HI)
  
  # Lambda_proj (final year / initial year total)
  lambda_proj <- pr$async$total[,T_proj] / pr$async$total[,1]
  
  tibble(
    Scenario          = sn,
    CV_Async          = round(cv_async,  4),
    CV_Sync           = round(cv_sync,   4),
    Portfolio_Ratio   = round(cv_async/cv_sync, 4),
    Buffering_Pct     = round(buffering, 2),
    N_final_mean      = round(n_final_mean, 0),
    N_final_lo        = round(n_final_lo,   0),
    N_final_hi        = round(n_final_hi,   0),
    Lambda_proj_mean  = round(mean(lambda_proj), 3),
    Lambda_proj_lo    = round(quantile(lambda_proj, CI_LO), 3),
    Lambda_proj_hi    = round(quantile(lambda_proj, CI_HI), 3)
  )
})

cat("\n── Portfolio Statistics by Scenario ────────────────────────────────────\n")
print(port_stats |> select(Scenario, CV_Async, CV_Sync, Portfolio_Ratio,
                           Buffering_Pct, Lambda_proj_mean), n=nrow(port_stats))
write_csv(port_stats, paste0("Output/",prefix,"_projection_portfolio_stats.csv"))

# Per-site CV table
site_cv_tbl <- map_dfr(scen_names, function(sn) {
  pr <- proj_results[[sn]]$async$by_site
  map_dfr(1:S, function(s) {
    vals <- pr[,T_proj,s]
    tibble(Scenario=sn, Site=site_names[s],
           Mean_N   = round(mean(vals), 1),
           CV       = round(sd(vals)/mean(vals), 4),
           N_lo     = round(quantile(vals, CI_LO), 1),
           N_hi     = round(quantile(vals, CI_HI), 1))
  })
})
cat("\n── Per-Site Final N and CV by Scenario ─────────────────────────────────\n")
print(site_cv_tbl |> select(Scenario, Site, Mean_N, CV), n=nrow(site_cv_tbl))
write_csv(site_cv_tbl, paste0("Output/",prefix,"_projection_site_cv.csv"))

# ============================================================================
# EXTRACT HISTORICAL TRAJECTORIES (1997–2025)
# ============================================================================
cat("Extracting historical population trajectories...
")

# Total population: N_total_all[t] from generated quantities
Ntot_hist_mat <- fit$draws(variables="N_total_all", format="matrix")
hist_total_df <- tibble(
  Period = "Historical",
  Year   = sim_data$years,
  mean   = colMeans(Ntot_hist_mat),
  lo     = as.numeric(apply(Ntot_hist_mat, 2, quantile, CI_LO)),
  hi     = as.numeric(apply(Ntot_hist_mat, 2, quantile, CI_HI))
)

# By-site total: N_total[s,t] from transformed parameters
# N_total = N_pup + N_juv_total + N_adult_total at each site
cat("  Extracting site-level historical N...
")
Ntot_site_mat <- fit$draws(variables="N_total", format="matrix")

hist_site_df <- map_dfr(1:S, function(s) {
  map_dfr(seq_along(sim_data$years), function(t) {
    cn <- paste0("N_total[",s,",",t,"]")
    if (!cn %in% colnames(Ntot_site_mat)) return(NULL)
    v <- Ntot_site_mat[, cn]
    tibble(Period="Historical",
           Site=site_names[s], Year=sim_data$years[t],
           mean=mean(v),
           lo=as.numeric(quantile(v, CI_LO)),
           hi=as.numeric(quantile(v, CI_HI)))
  })
}) |> mutate(Site=factor(Site, levels=site_names))

# Boundary year for vertical separator line
hist_boundary <- max(sim_data$years)   # 2025
cat(sprintf("  Historical: %d–%d | Projection: %d–%d
",
            min(sim_data$years), hist_boundary,
            min(proj_years), max(proj_years)))

# ============================================================================
# BUILD TIDY PLOTTING DATA FRAMES
# ============================================================================

# (a) Total population trajectories (async)
total_df <- map_dfr(scen_names, function(sn) {
  m <- proj_results[[sn]]$async$total
  tibble(
    Scenario = sn,
    Year     = proj_years,
    mean     = colMeans(m),
    lo       = as.numeric(apply(m, 2, quantile, CI_LO)),
    hi       = as.numeric(apply(m, 2, quantile, CI_HI))
  )
}) |> mutate(Scenario=factor(Scenario, levels=scen_names))

# (b) By-site trajectories (async)
site_df <- map_dfr(scen_names, function(sn) {
  arr <- proj_results[[sn]]$async$by_site
  map_dfr(1:S, function(s) {
    m <- arr[,,s]
    tibble(Scenario=sn, Site=site_names[s], Year=proj_years,
           mean = colMeans(m),
           lo   = as.numeric(apply(m, 2, quantile, CI_LO)),
           hi   = as.numeric(apply(m, 2, quantile, CI_HI)))
  })
}) |> mutate(Scenario=factor(Scenario, levels=scen_names),
             Site=factor(Site, levels=site_names))

# Palette: 7 scenarios
scen_cols <- setNames(
  c("#2166AC","#B2182B","#4DAC26","#8C510A","#762A83","#D6604D","#1B7837"),
  scen_names
)

# ============================================================================
# PLOTS — (b) TOTAL POPULATION
# ============================================================================

# Plot 1: Total trajectories — historical + all scenarios overlaid
p_total_all <- ggplot() +
  # Historical period: single gray ribbon + line (same for all scenarios)
  geom_ribbon(data=hist_total_df,
              aes(x=Year, ymin=lo, ymax=hi),
              fill="gray70", alpha=0.35, colour=NA) +
  geom_line(data=hist_total_df,
            aes(x=Year, y=mean),
            colour="gray25", linewidth=1.3) +
  # Projection period: coloured by scenario
  geom_ribbon(data=total_df,
              aes(x=Year, ymin=lo, ymax=hi, fill=Scenario),
              alpha=0.12, colour=NA) +
  geom_line(data=total_df,
            aes(x=Year, y=mean, colour=Scenario),
            linewidth=1.1) +
  # Boundary between observed and projected
  geom_vline(xintercept=hist_boundary + 0.5,
             linetype="dashed", colour="gray40", linewidth=0.7) +
  annotate("text", x=hist_boundary - 0.5, y=Inf,
           label="Observed", hjust=1, vjust=1.5,
           size=3.5, colour="gray40") +
  annotate("text", x=hist_boundary + 1.5, y=Inf,
           label="Projected", hjust=0, vjust=1.5,
           size=3.5, colour="gray40") +
  scale_colour_manual(values=scen_cols) +
  scale_fill_manual(  values=scen_cols) +
  labs(x="Year", y="Total Population",
       title="Harbor Seal Population: Observed (1997–2025) and 10-Year Projection (2026–2035)",
       subtitle=paste0("Posterior mean ± ", CI_LABEL, " | Asynchronous site dynamics"),
       colour=NULL, fill=NULL) +
  theme_seal() +
  theme(legend.position="bottom",
        legend.text=element_text(size=10)) +
  guides(colour=guide_legend(nrow=3), fill=guide_legend(nrow=3))
ggsave(paste0("Output/Plots/",prefix,"_proj_portfolio_total.jpeg"),
       p_total_all, width=34, height=18, units="cm")

# Plot 2: Total trajectories — faceted by scenario (historical repeated in each panel)
# hist_total_df has no "Scenario" column → ggplot repeats it across all facets
p_total_facet <- ggplot(total_df, aes(x=Year, y=mean, colour=Scenario)) +
  # Historical (in every panel — inherit.aes=FALSE so Scenario aes not required)
  geom_ribbon(data=hist_total_df, inherit.aes=FALSE,
              aes(x=Year, ymin=lo, ymax=hi),
              fill="gray70", alpha=0.35, colour=NA) +
  geom_line(data=hist_total_df, inherit.aes=FALSE,
            aes(x=Year, y=mean),
            colour="gray25", linewidth=1.2) +
  # Projection
  geom_ribbon(aes(ymin=lo, ymax=hi, fill=Scenario), alpha=0.18, colour=NA) +
  geom_line(linewidth=1.1) +
  # Boundary
  geom_vline(xintercept=hist_boundary + 0.5,
             linetype="dashed", colour="gray40", linewidth=0.5) +
  facet_wrap(~Scenario, ncol=3, scales="free_y") +
  scale_colour_manual(values=scen_cols, guide="none") +
  scale_fill_manual(  values=scen_cols, guide="none") +
  labs(x="Year", y="Total Population",
       title="Harbor Seal Population by Scenario: Observed + 10-Year Projection",
       subtitle="Gray = observed 1997–2025; coloured = projected 2026–2035; dashed = projection start") +
  theme_seal() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8))
ggsave(paste0("Output/Plots/",prefix,"_proj_portfolio_total_facet.jpeg"),
       p_total_facet, width=36, height=22, units="cm")

# Plot 3: Portfolio buffering ratio by scenario (bar chart)
p_buffer <- port_stats |>
  mutate(Scenario=factor(Scenario, levels=scen_names)) |>
  ggplot(aes(x=Scenario, y=Buffering_Pct, fill=Scenario)) +
  geom_col(alpha=0.85) +
  geom_text(aes(label=sprintf("%.1f%%", Buffering_Pct)),
            vjust=-0.3, size=4, fontface="bold") +
  geom_hline(yintercept=0, linetype=2, colour="gray40") +
  scale_fill_manual(values=scen_cols, guide="none") +
  labs(x=NULL, y="Portfolio Buffering (%)",
       title="Portfolio Buffering by Scenario",
       subtitle="% reduction in CV of total N from spatial asynchrony\n(positive = asynchronous sites buffer aggregate variability)") +
  theme_seal() +
  theme(axis.text.x=element_text(angle=25, hjust=1, size=10))
ggsave(paste0("Output/Plots/",prefix,"_proj_portfolio_buffering.jpeg"),
       p_buffer, width=28, height=16, units="cm")

# Plot 4: Final-year N with uncertainty — all scenarios
p_final_n <- port_stats |>
  mutate(Scenario=factor(Scenario, levels=rev(scen_names))) |>
  ggplot(aes(y=Scenario, colour=Scenario)) +
  geom_linerange(aes(xmin=N_final_lo, xmax=N_final_hi), linewidth=2, alpha=0.7) +
  geom_point(aes(x=N_final_mean), size=4) +
  geom_vline(xintercept=port_stats$N_final_mean[port_stats$Scenario=="Status Quo"],
             linetype=2, colour="gray40") +
  scale_colour_manual(values=scen_cols, guide="none") +
  labs(x=paste0("Projected Total N in year ",max(proj_years)),
       y=NULL,
       title=paste0("Projected Total Population in Year +",T_proj),
       subtitle=paste0(CI_LABEL, " | dashed line = Status Quo mean")) +
  theme_seal() +
  theme(panel.grid.major.y=element_blank())
ggsave(paste0("Output/Plots/",prefix,"_proj_portfolio_final_n.jpeg"),
       p_final_n, width=26, height=16, units="cm")

# Plot 5: Async vs Sync CV comparison (portfolio effect visualisation)
p_cv_comp <- port_stats |>
  select(Scenario, CV_Async, CV_Sync) |>
  pivot_longer(c(CV_Async,CV_Sync), names_to="Type", values_to="CV") |>
  mutate(Type=recode(Type,
                     CV_Async="Asynchronous (current)",
                     CV_Sync ="Synchronous (hypothetical)"),
         Scenario=factor(Scenario, levels=scen_names)) |>
  ggplot(aes(x=Scenario, y=CV, fill=Type, alpha=Type)) +
  geom_col(position="dodge") +
  scale_fill_manual(values=c("Asynchronous (current)"=SEAL_COLS$pop,
                             "Synchronous (hypothetical)"=SEAL_COLS$adult_f),
                    name=NULL) +
  scale_alpha_manual(values=c("Asynchronous (current)"=0.9,
                              "Synchronous (hypothetical)"=0.6),
                     guide="none") +
  labs(x=NULL, y="CV of total population",
       title="Spatial Asynchrony Buffers Population Variability",
       subtitle="CV(async) < CV(sync) confirms portfolio effect; gap = buffering magnitude") +
  theme_seal() +
  theme(axis.text.x=element_text(angle=25, hjust=1, size=10),
        legend.position="top")
ggsave(paste0("Output/Plots/",prefix,"_proj_portfolio_cv_comparison.jpeg"),
       p_cv_comp, width=28, height=16, units="cm")

# ============================================================================
# PLOTS — (a) BY SITE
# ============================================================================

# Plot 6: By-site trajectories — historical + all scenarios, one panel per site
# hist_site_df has Site but no Scenario → repeated across all scenarios in each site panel
p_site_all <- ggplot(site_df, aes(x=Year, y=mean, colour=Scenario, fill=Scenario)) +
  # Historical per site (no Scenario → shows in every scenario layer; facet=Site matches)
  geom_ribbon(data=hist_site_df, inherit.aes=FALSE,
              aes(x=Year, ymin=lo, ymax=hi, group=Site),
              fill="gray70", alpha=0.35, colour=NA) +
  geom_line(data=hist_site_df, inherit.aes=FALSE,
            aes(x=Year, y=mean, group=Site),
            colour="gray25", linewidth=1.1) +
  # Projection
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.10, colour=NA) +
  geom_line(linewidth=0.9) +
  # Boundary
  geom_vline(xintercept=hist_boundary + 0.5,
             linetype="dashed", colour="gray40", linewidth=0.5) +
  facet_wrap(~Site, scales="free_y", ncol=2) +
  scale_colour_manual(values=scen_cols) +
  scale_fill_manual(  values=scen_cols) +
  labs(x="Year", y="Site Population",
       title="Harbor Seal Population by Site: Observed + 10-Year Projection",
       subtitle=paste0("Gray = observed 1997–2025; coloured = projected 2026–2035 | ",
                       CI_LABEL),
       colour=NULL, fill=NULL) +
  theme_seal() +
  theme(legend.position="bottom",
        legend.text=element_text(size=9),
        axis.text.x=element_text(angle=45, hjust=1, size=8)) +
  guides(colour=guide_legend(nrow=3), fill=guide_legend(nrow=3))
ggsave(paste0("Output/Plots/",prefix,"_proj_portfolio_by_site.jpeg"),
       p_site_all, width=30, height=36, units="cm")

# Plot 7: By-site × key-scenario grid with historical overlay
key_scens   <- c("Status Quo","Warm + High Stress","Cool + Low Stress")
key_scen_df <- site_df |>
  filter(Scenario %in% key_scens) |>
  mutate(Scenario=factor(Scenario, levels=key_scens))

p_site_key <- ggplot(key_scen_df,
                     aes(x=Year, y=mean, colour=Scenario, fill=Scenario)) +
  # Historical per site — hist_site_df has Site but no Scenario, so it appears
  # in every column of the facet_grid (one column per scenario)
  geom_ribbon(data=hist_site_df, inherit.aes=FALSE,
              aes(x=Year, ymin=lo, ymax=hi, group=Site),
              fill="gray70", alpha=0.35, colour=NA) +
  geom_line(data=hist_site_df, inherit.aes=FALSE,
            aes(x=Year, y=mean, group=Site),
            colour="gray25", linewidth=1.1) +
  # Projection
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.15, colour=NA) +
  geom_line(linewidth=1.0) +
  # Boundary
  geom_vline(xintercept=hist_boundary + 0.5,
             linetype="dashed", colour="gray40", linewidth=0.5) +
  facet_grid(Site~Scenario, scales="free_y") +
  scale_colour_manual(values=scen_cols[key_scens], guide="none") +
  scale_fill_manual(  values=scen_cols[key_scens], guide="none") +
  labs(x="Year", y="Site Population",
       title="Harbor Seal Population by Site: Observed + 3 Key Scenarios",
       subtitle="Gray = observed 1997–2025; coloured = projected 2026–2035; y-axis free within each site") +
  theme_seal(base_size=13) +
  theme(strip.text.y=element_text(size=10, face="bold"),
        strip.text.x=element_text(size=10, face="bold"),
        axis.text.x=element_text(angle=45, hjust=1, size=8))
ggsave(paste0("Output/Plots/",prefix,"_proj_portfolio_site_key_scenarios.jpeg"),
       p_site_key, width=34, height=38, units="cm")

# Plot 8: Per-site CV heatmap across scenarios
p_site_cv <- site_cv_tbl |>
  mutate(Scenario=factor(Scenario, levels=scen_names),
         Site=factor(Site, levels=site_names)) |>
  ggplot(aes(x=Scenario, y=Site, fill=CV)) +
  geom_tile(colour="white", linewidth=0.5) +
  geom_text(aes(label=sprintf("%.3f", CV)), size=3.5, fontface="bold") +
  scale_fill_gradient2(low="darkgreen", mid="lightyellow", high="red3",
                       midpoint=median(site_cv_tbl$CV),
                       name="CV of N") +
  labs(x=NULL, y="Site",
       title=paste0("Site-Level Population Variability (CV) at Year +",T_proj),
       subtitle="Darker red = higher variability; values are CV of projected N across posterior draws") +
  theme_seal() +
  theme(axis.text.x=element_text(angle=30, hjust=1, size=10),
        panel.grid=element_blank())
ggsave(paste0("Output/Plots/",prefix,"_proj_portfolio_site_cv_heatmap.jpeg"),
       p_site_cv, width=30, height=16, units="cm")

# Plot 9: Stacked area — composition of total pop by site under each scenario
# (Shows which sites drive population under each future)
site_final_df <- site_cv_tbl |>
  mutate(Scenario=factor(Scenario, levels=scen_names),
         Site=factor(Site, levels=site_names))

p_stack <- ggplot(site_final_df,
                  aes(x=Scenario, y=Mean_N, fill=Site)) +
  geom_col(position="stack", alpha=0.85) +
  scale_fill_brewer(palette="Dark2") +
  labs(x=NULL, y=paste0("Projected Total N in year +",T_proj),
       title="Site Contributions to Projected Total Population",
       subtitle="Stacked = total population; segment height = site contribution") +
  theme_seal() +
  theme(axis.text.x=element_text(angle=25, hjust=1, size=10))
ggsave(paste0("Output/Plots/",prefix,"_proj_portfolio_site_stack.jpeg"),
       p_stack, width=28, height=16, units="cm")

# ============================================================================
# COMBINED SUMMARY FIGURE
# ============================================================================
p_summary <- (p_total_all + theme(legend.position="right")) /
  (p_buffer + p_cv_comp) +
  plot_annotation(
    title   = "Harbor Seal Portfolio Projection — 10-Year Summary",
    subtitle= paste0(T_proj," year projection from ", max(sim_data$years),
                     "; ", n_sims, " posterior draws"),
    theme   = theme(plot.title=element_text(size=16, face="bold"))
  ) +
  plot_layout(heights=c(1.4, 1))
ggsave(paste0("Output/Plots/",prefix,"_proj_portfolio_summary.jpeg"),
       p_summary, width=36, height=30, units="cm")

# ============================================================================
# COLLECT RESULTS
# ============================================================================
portfolio_proj <- list(
  # Plots
  total_all       = p_total_all,
  total_facet     = p_total_facet,
  buffering       = p_buffer,
  final_n         = p_final_n,
  cv_comparison   = p_cv_comp,
  by_site         = p_site_all,
  site_key        = p_site_key,
  site_cv_heatmap = p_site_cv,
  site_stack      = p_stack,
  summary_figure  = p_summary,
  # Tables
  portfolio_stats = port_stats,
  site_cv_table   = site_cv_tbl,
  # Raw projections (array: n_sims × T_proj × S for by_site; matrix for total)
  raw             = proj_results
)

cat("\n── Portfolio projection analysis complete ───────────────────────\n")
cat(sprintf("   Scenarios: %d\n", length(scenarios)))
cat(sprintf("   Draws:     %d\n", n_sims))
cat(sprintf("   Horizon:   %d years (%d–%d)\n",
            T_proj, min(proj_years), max(proj_years)))
cat(sprintf("   Plots (9)  -> Output/Plots/%s_proj_portfolio_*.jpeg\n", prefix))
cat(sprintf("   Tables (2) -> Output/%s_projection_portfolio_*.csv\n", prefix))
cat("\n── Portfolio Effect Ratio by Scenario ───────────────────────────\n")
for (i in seq_len(nrow(port_stats))) {
  cat(sprintf("   %-30s  ratio=%.3f  buffering=%.1f%%  lambda=%.3f\n",
              port_stats$Scenario[i],
              port_stats$Portfolio_Ratio[i],
              port_stats$Buffering_Pct[i],
              port_stats$Lambda_proj_mean[i]))
}
