# ============================================================================
# 14_marss_legacy_plots.R   (was 7c_plots_legacy_MARSS.R)
# ----------------------------------------------------------------------------
# PLOTS — LONG-TERM MARSS ANALYSIS (1975–2025). Style matches 11_marss8_plots.R.
#
# Source first:
#   source("Code/12_data_prep_legacy.R")
#   source("Code/13_marss_legacy.R")
#   OR load saved objects (auto-loaded below if missing):
#   load("Output/m.LA_legacy.RData"); load("Output/CIs_legacy.RData")
#
# FIX (P8): d_obs is now defined BEFORE the blocks that consume it
# (first_obs_year / d_obs_diff / d_obs_BL_DP). The original referenced d_obs
# above its own definition and crashed on first run. Duplicate d_obs_diff
# definition also removed.
# ============================================================================

library(MARSS)
library(tidyverse)
library(patchwork)
library(broom)   # tidy() for CIs_legacy

# ── Constants ───────────────────────────────────────────────────────────────
CI_LO <- 0.055; CI_HI <- 0.945; CI_LABEL <- "89% CI"
SITE_COLS_LEGACY <- c(BL="#E41A1C", DE="#FF7F00", DP="#4DAF4A",
                      PRH="#377EB8", TB="#A65628", TP="#F781BF")
bp_year <- 2004
dir.create("Output/Plots", showWarnings = FALSE)

theme_seal <- function(base_size = 14) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.grid.major = element_line(colour="grey88", linewidth=0.4),
      panel.border     = element_rect(colour="grey70", fill=NA, linewidth=0.5),
      axis.title       = element_text(size=rel(0.95)),
      legend.position  = "bottom",
      strip.text       = element_text(size=rel(0.90), face="bold"),
      strip.background = element_rect(fill="grey94", colour="grey80"),
      plot.title       = element_text(size=rel(1.05), face="bold"),
      plot.caption     = element_text(size=rel(0.78), colour="grey50", hjust=1),
      plot.margin      = margin(10,14,10,10))
}

# ── Load model + CIs if not in memory ───────────────────────────────────────
if (!exists("m.LA"))        { load("Output/m.LA_legacy.RData");  cat("Loaded m.LA_legacy.RData\n") }
if (!exists("CIs_legacy"))  { load("Output/CIs_legacy.RData");   cat("Loaded CIs_legacy.RData\n") }
if (!exists("dat_legacy"))  stop("dat_legacy not found — run 12_data_prep_legacy.R.")
if (!exists("years_legacy")) stop("years_legacy not found — run 12_data_prep_legacy.R.")

BESTMODEL_L <- if (exists("BEST_LEGACY")) BEST_LEGACY else m.LA

# ── Helper: long-form smoothed state data ───────────────────────────────────
build_legacy_df <- function(model, state_names, yrs) {
  d    <- as_tibble(t(model$states))
  d.se <- as_tibble(t(model$states.se))
  names(d) <- names(d.se) <- state_names
  d$Year <- yrs
  d_long <- d %>%
    pivot_longer(-Year, names_to="State", values_to="log_est") %>%
    separate(State, into=c("Site","Class"), sep="_", remove=FALSE, extra="merge")
  d_se_long <- d.se %>%
    pivot_longer(everything(), names_to="State", values_to="log_se")
  bind_cols(d_long, d_se_long %>% select(log_se)) %>%
    mutate(lo89=log_est+qnorm(CI_LO)*log_se, hi89=log_est+qnorm(CI_HI)*log_se,
           estimate=exp(log_est), lo89_n=exp(lo89), hi89_n=exp(hi89),
           Class=factor(Class, levels=c("PUP","ADULT","MOLTING"),
                        labels=c("Pup","Adult","Molt")))
}
d_legacy <- build_legacy_df(BESTMODEL_L, state_names_legacy, years_legacy)

# ── Observed data points (DEFINED FIRST — fixes P8) ─────────────────────────
# dat_legacy holds log counts; NA = unobserved. Exponentiate for plotting.
d_obs <- as.data.frame(t(dat_legacy)) %>%
  mutate(Year = years_legacy) %>%
  pivot_longer(-Year, names_to="State", values_to="log_obs") %>%
  filter(!is.na(log_obs)) %>%
  separate(State, into=c("Site","Class"), sep="_", remove=FALSE, extra="merge") %>%
  mutate(obs=exp(log_obs),
         Class=factor(Class, levels=c("PUP","ADULT","MOLTING"),
                      labels=c("Pup","Adult","Molt")))
cat("Observed data points:", nrow(d_obs), "across",
    length(unique(d_obs$State)), "state-classes\n")

# ── Now derived objects that consume d_obs ──────────────────────────────────
first_obs_year <- d_obs %>%
  group_by(State) %>%
  summarise(first_obs_year = min(Year), .groups="drop")

# Don't plot smoothed estimates before first real observation per state.
d_legacy <- d_legacy %>%
  left_join(first_obs_year, by="State") %>%
  filter(Year >= first_obs_year)

d_obs_BL_DP <- d_obs %>% filter(Site %in% c("BL","DP"))

# ── PLOT L1: smoothed trajectories + observed points (real scale) ───────────
p.L1 <- ggplot(d_legacy, aes(Year, estimate, colour=Site, fill=Site, group=Site)) +
  geom_vline(xintercept=bp_year, linetype=2, colour="grey40", linewidth=0.6) +
  geom_ribbon(aes(ymin=lo89_n, ymax=hi89_n), alpha=0.12, colour=NA) +
  geom_line(linewidth=1) +
  geom_point(data=d_obs, aes(Year, obs, colour=Site, group=Site),
             shape=16, size=1.8, alpha=0.7, inherit.aes=FALSE) +
  expand_limits(y=0) +
  scale_colour_manual(values=SITE_COLS_LEGACY) +
  scale_fill_manual(values=SITE_COLS_LEGACY) +
  scale_y_continuous(labels=scales::comma) +
  facet_wrap(~Class, ncol=1, scales="free_y") +
  labs(x="Year", y="Estimated abundance",
       title="Long-Term Harbor Seal Population Trends — PRNS (1975–2025)",
       subtitle=paste0("Lines = MARSS smoothed; bands = ", CI_LABEL,
                       "; points = observed annual maxima"),
       caption=paste0("Vertical dashed = ", bp_year,
                      " breakpoint. Pre-1997 from legacy survey records.")) +
  theme_seal() + theme(legend.position="right")
ggsave("Output/Plots/legacy_smoothed_trajectories.jpeg", p.L1,
       width=24, height=30, units="cm", dpi=200)

# ── PLOT L2: percent change relative to first observation + observed ─────────
d_diff_legacy <- d_legacy %>%
  group_by(State) %>%
  mutate(first_obs=first(log_est[!is.na(log_est)]),
         log_diff=log_est-first_obs,
         lo89_d=log_diff+qnorm(CI_LO)*log_se,
         hi89_d=log_diff+qnorm(CI_HI)*log_se,
         pct_change=(exp(log_diff)-1)*100,
         pct_lo=(exp(lo89_d)-1)*100, pct_hi=(exp(hi89_d)-1)*100) %>%
  ungroup()

d_obs_diff <- d_obs %>%
  group_by(State) %>%
  mutate(first_obs_log=first(log_obs),
         log_diff_obs=log_obs-first_obs_log,
         pct_change_obs=(exp(log_diff_obs)-1)*100) %>%
  ungroup()

p.L2 <- ggplot(d_diff_legacy, aes(Year, pct_change, colour=Site, group=Site)) +
  geom_hline(yintercept=0, linetype=2, colour="grey40") +
  geom_vline(xintercept=bp_year, linetype=2, colour="grey40", linewidth=0.6) +
  geom_ribbon(aes(ymin=pct_lo, ymax=pct_hi, fill=Site), alpha=0.12, colour=NA) +
  geom_line(linewidth=1) +
  geom_point(data=d_obs_diff, aes(Year, pct_change_obs, colour=Site, group=Site),
             shape=16, size=1.8, alpha=0.7, inherit.aes=FALSE) +
  scale_colour_manual(values=SITE_COLS_LEGACY) +
  scale_fill_manual(values=SITE_COLS_LEGACY) +
  scale_y_continuous(labels=function(x) paste0(x,"%")) +
  facet_wrap(~Class, ncol=1, scales="free_y") +
  labs(x="Year", y="% change from first observation",
       title="Long-Term Population Trends — Percent Change by Site",
       subtitle=paste0("Lines = MARSS smoothed; points = observed; bands = ", CI_LABEL),
       caption=paste0("Vertical dashed = ", bp_year, " breakpoint.")) +
  theme_seal() + theme(legend.position="right")
ggsave("Output/Plots/legacy_pct_change.jpeg", p.L2, width=24, height=30, units="cm", dpi=200)

# ── PLOT L3: total abundance summed across 6 sites + observed sum ────────────
d_total_legacy <- d_legacy %>%
  group_by(Year, Class) %>%
  summarise(total=sum(estimate, na.rm=TRUE), total_lo=sum(lo89_n, na.rm=TRUE),
            total_hi=sum(hi89_n, na.rm=TRUE), n_sites=sum(!is.na(estimate)),
            .groups="drop")
d_obs_total <- d_obs %>%
  group_by(Year, Class) %>%
  summarise(obs_total=sum(obs, na.rm=TRUE), n_sites=n(), .groups="drop")
p.L3 <- ggplot(d_total_legacy, aes(Year, total)) +
  geom_vline(xintercept=bp_year, linetype=2, colour="grey40", linewidth=0.6) +
  geom_ribbon(aes(ymin=total_lo, ymax=total_hi), alpha=0.20, fill="#2166AC", colour=NA) +
  geom_line(colour="#2166AC", linewidth=1.2) +
  geom_point(data=d_obs_total, aes(Year, obs_total, size=n_sites),
             colour="#2166AC", shape=16, alpha=0.75, inherit.aes=FALSE) +
  scale_size_continuous(name="Sites observed", range=c(1.5,4), breaks=1:6) +
  annotate("rect", xmin=2019.5, xmax=2020.5, ymin=-Inf, ymax=Inf, fill="grey60", alpha=0.15) +
  annotate("text", x=2020, y=Inf, label="COVID\n(incomplete)", vjust=1.4, hjust=0.5,
           size=3, colour="grey40") +
  expand_limits(y=0) + scale_y_continuous(labels=scales::comma) +
  facet_wrap(~Class, ncol=1, scales="free_y") +
  labs(x="Year", y="Total estimated abundance (6 sites)",
       title="Total Estimated Harbor Seal Abundance — All 6 Breeding Sites",
       subtitle=paste0("Line = MARSS sum; points = observed sum; bands = ", CI_LABEL),
       caption=paste0("Vertical dashed = ", bp_year,
                      " breakpoint. Point size = # sites observed.")) +
  theme_seal()
ggsave("Output/Plots/legacy_total_abundance.jpeg", p.L3, width=22, height=30, units="cm", dpi=200)

# ── PLOT L4: growth-rate estimates (bar chart) ──────────────────────────────
u_df <- tidy(CIs_legacy) %>%
  tibble::as_tibble() %>%
  filter(grepl("^U\\.", term)) %>%
  mutate(
    Phase = case_when(grepl("t1_", term) ~ paste0("Phase 1 (pre-", bp_year, ")"),
                      grepl("t2_", term) ~ paste0("Phase 2 (", bp_year, "+)"),
                      TRUE ~ "Constant"),
    Class = case_when(grepl("_A", term)~"Adult", grepl("_M", term)~"Molt",
                      grepl("_P", term)~"Pup", TRUE~"Unknown"),
    Class = factor(Class, levels=c("Pup","Adult","Molt")),
    direction = ifelse(estimate > 0, "Positive (growth)", "Negative (decline)"))
p.L4 <- ggplot(u_df, aes(Phase, estimate, fill=direction, colour=direction)) +
  geom_hline(yintercept=0, linetype=2, colour="grey40") +
  geom_col(alpha=0.7, width=0.55) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.up), width=0.2, linewidth=0.8, colour="grey20") +
  geom_text(aes(label=sprintf("%.3f", estimate),
                vjust=ifelse(estimate >= 0, -0.5, 1.4)), size=3.5, colour="grey20") +
  scale_fill_manual(values=c("Positive (growth)"="#2166AC","Negative (decline)"="#B2182B"), guide="none") +
  scale_colour_manual(values=c("Positive (growth)"="#2166AC","Negative (decline)"="#B2182B"), guide="none") +
  facet_wrap(~Class, ncol=3) +
  labs(x=NULL, y="Log-scale growth rate (per year)",
       title="Estimated Population Growth Rates by Phase and Class",
       subtitle=paste0("Error bars = ", CI_LABEL)) +
  theme_seal() + theme(axis.text.x=element_text(angle=20, hjust=1))
ggsave("Output/Plots/legacy_growth_rates.jpeg", p.L4, width=26, height=14, units="cm", dpi=200)

# ── PLOT L5: BL & DP detail + observed points ───────────────────────────────
d_long_sites <- d_legacy %>% filter(Site %in% c("BL","DP"))
p.L5 <- ggplot(d_long_sites, aes(Year, estimate, colour=Site, fill=Site)) +
  geom_vline(xintercept=bp_year, linetype=2, colour="grey40", linewidth=0.6) +
  geom_ribbon(aes(ymin=lo89_n, ymax=hi89_n), alpha=0.20, colour=NA) +
  geom_line(linewidth=1.3) +
  geom_point(data=d_obs_BL_DP, aes(Year, obs, colour=Site, group=Site),
             shape=16, size=2.5, alpha=0.8, inherit.aes=FALSE) +
  expand_limits(y=0) +
  scale_colour_manual(values=c(BL="#E41A1C", DP="#4DAF4A")) +
  scale_fill_manual(values=c(BL="#E41A1C", DP="#4DAF4A")) +
  scale_y_continuous(labels=scales::comma) +
  facet_wrap(~Class, ncol=1, scales="free_y") +
  labs(x="Year", y="Estimated abundance",
       title="Bolinas Lagoon and Drakes Estero — Full 50-Year Record",
       subtitle=paste0("Lines = MARSS smoothed; points = observed; bands = ", CI_LABEL),
       caption="BL = Bolinas Lagoon; DP = Drakes Estero Pup Beach") +
  theme_seal() + theme(legend.position="bottom")
ggsave("Output/Plots/legacy_BL_DP_detail.jpeg", p.L5, width=22, height=30, units="cm", dpi=200)

# ── Summary ─────────────────────────────────────────────────────────────────
cat("\n── Long-term trend summary ─────────────────────────────────────────\n")
cat("Year range:", min(years_legacy), "-", max(years_legacy),
    "(", length(years_legacy), "years)\n")
cat("\nPhase 1 adult growth (log scale):",
    round(BESTMODEL_L$par$U["t1_A", 1], 4), "\n")
if ("t2_A" %in% rownames(BESTMODEL_L$par$U))
  cat("Phase 2 adult change (log scale):",
      round(BESTMODEL_L$par$U["t2_A", 1], 4), "\n")
cat("\nAll plots -> Output/Plots/legacy_*\n")