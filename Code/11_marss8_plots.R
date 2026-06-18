# ============================================================================
# 11_marss8_plots.R
# ----------------------------------------------------------------------------
# 8-SITE MARSS — PLOTS (Appendix A)
# Visualises results from the three candidate models (A, B, C) fit in
# 10_marss8_models.R. Generates figures for Appendix A in the manuscript.
#
# Source first (either):
#   source("Code/08_data_prep_8site.R")
#   source("Code/09_covariates_8site.R")
#   source("Code/10_marss8_models.R")       # runs models + saves RData
# OR auto-loads saved objects if they exist:
#   load("Output/m.A_8site.RData"); load("Output/CIs_8site.RData"); etc.
#
# Produces:
#   Output/Plots/marss8_smoothed_trajectories.jpeg  (Figure A1)
#   Output/Plots/marss8_summed_totals.jpeg          (Figure A2)
#   Output/Plots/marss8_covariate_effects.jpeg      (Figure A3)
#   Output/model_comparison_8site_AIC.csv           (Table S2 data)
# ============================================================================

library(MARSS)
library(tidyverse)
library(patchwork)
library(broom)

dir.create("Output/Plots", showWarnings = FALSE)

CI_LO <- 0.055; CI_HI <- 0.945; CI_LABEL <- "89% CI"
SITES_8 <- c("BL","DE","DP","DR","PB","PRH","TB","TP")
bp_year <- 2004

# ── Colour palette ────────────────────────────────────────────────────────────
# Breeding sites match IPM palette; haul-out sites use distinct colours.
SITE_COLS_8 <- c(
  BL  = "#E41A1C",   # red
  DE  = "#FF7F00",   # orange
  DP  = "#4DAF4A",   # green
  PRH = "#377EB8",   # blue
  TB  = "#A65628",   # brown
  TP  = "#F781BF",   # pink
  DR  = "#984EA3",   # purple  (haul-out)
  PB  = "#00CED1"    # cyan    (haul-out)
)

theme_seal <- function(base_size = 14) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.grid.major = element_line(colour = "grey88", linewidth = 0.4),
      panel.border     = element_rect(colour = "grey70", fill = NA, linewidth = 0.5),
      axis.title       = element_text(size = rel(0.95)),
      legend.position  = "right",
      strip.text       = element_text(size = rel(0.90), face = "bold"),
      strip.background = element_rect(fill = "grey94", colour = "grey80"),
      plot.title       = element_text(size = rel(1.05), face = "bold"),
      plot.caption     = element_text(size = rel(0.78), colour = "grey50", hjust = 1),
      plot.margin      = margin(10, 14, 10, 10))
}

# ── Auto-load saved objects if not in memory ──────────────────────────────────
if (!exists("m.A_8site")) {
  load("Output/m.A_8site.RData"); cat("Loaded m.A_8site\n")
}
if (!exists("CIs_8site")) {
  load("Output/CIs_8site.RData"); cat("Loaded CIs_8site\n")
}
if (!exists("df_aic_8site")) {
  if (file.exists("Output/df_aic_8site.RData")) {
    load("Output/df_aic_8site.RData")
  } else if (file.exists("Output/model_comparison_8site_AIC.csv")) {
    df_aic_8site <- read.csv("Output/model_comparison_8site_AIC.csv")
  }
}
if (!exists("dat_8site"))
  stop("dat_8site not found — run 08_data_prep_8site.R first.")
if (!exists("years_8site"))
  stop("years_8site not found — run 08_data_prep_8site.R first.")
if (!exists("cov_t_scaled_8site"))
  stop("cov_t_scaled_8site not found — run 09_covariates_8site.R first.")

BESTMODEL_8 <- m.A_8site

# ── Build long-form smoothed state data ───────────────────────────────────────
state_names_8 <- rownames(dat_8site)
TT <- ncol(dat_8site)

# IMPORTANT: MARSS labels states from the Z factor levels ("1".."24"), NOT
# from rownames(dat_8site), even though Z = factor(1:24) was built in the
# same order as state_names_8. Without this explicit rename, every state
# name fails to split on "_" below (silent NA fill, not an error) — fixed
# here the same way 14_marss_legacy_plots.R does it for the legacy model.
states_mat    <- BESTMODEL_8$states
states_se_mat <- BESTMODEL_8$states.se
stopifnot(nrow(states_mat) == length(state_names_8))   # order must match
rownames(states_mat)    <- state_names_8
rownames(states_se_mat) <- state_names_8

d8 <- as.data.frame(t(states_mat)) %>%
  mutate(Year = years_8site) %>%
  pivot_longer(-Year, names_to = "State", values_to = "log_est") %>%
  bind_cols(
    as.data.frame(t(states_se_mat)) %>%
      pivot_longer(everything(), names_to = "State2", values_to = "log_se") %>%
      select(log_se)
  ) %>%
  separate(State, into = c("Site","Class"), sep = "_", remove = FALSE, extra = "merge") %>%
  mutate(
    lo89 = log_est + qnorm(CI_LO) * log_se,
    hi89 = log_est + qnorm(CI_HI) * log_se,
    estimate  = exp(log_est),
    lo89_n    = exp(lo89),
    hi89_n    = exp(hi89),
    Site      = factor(Site, levels = SITES_8),
    Class     = factor(Class, levels = c("ADULT","MOLTING","PUP"),
                       labels = c("Adult","Molt","Pup")),
    SiteType  = ifelse(Site %in% c("DR","PB"), "Haul-out", "Breeding")
  )

# Observed data points
d8_obs <- as.data.frame(t(dat_8site)) %>%
  mutate(Year = years_8site) %>%
  pivot_longer(-Year, names_to = "State", values_to = "log_obs") %>%
  filter(!is.na(log_obs)) %>%
  separate(State, into = c("Site","Class"), sep = "_", remove = FALSE, extra = "merge") %>%
  mutate(
    obs   = exp(log_obs),
    Site  = factor(Site, levels = SITES_8),
    Class = factor(Class, levels = c("ADULT","MOLTING","PUP"),
                   labels = c("Adult","Molt","Pup")),
    SiteType = ifelse(Site %in% c("DR","PB"), "Haul-out", "Breeding")
  )

# Filter smoothed lines to first observed year per state
first_obs <- d8_obs %>%
  group_by(State) %>%
  summarise(first_obs_year = min(Year), .groups = "drop")

d8 <- d8 %>%
  left_join(first_obs, by = "State") %>%
  filter(Year >= first_obs_year)

# ── FIGURE A1: Smoothed trajectories — all 8 sites ───────────────────────────
p.A1 <- ggplot(d8 %>% filter(Class != "Pup" | SiteType == "Breeding"),
               aes(x = Year, y = estimate, colour = Site, fill = Site,
                   linetype = SiteType, group = Site)) +
  geom_vline(xintercept = bp_year, linetype = 2, colour = "grey40", linewidth = 0.6) +
  geom_ribbon(aes(ymin = lo89_n, ymax = hi89_n), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(data = d8_obs %>% filter(Class != "Pup" | SiteType == "Breeding"),
             aes(x = Year, y = obs), shape = 16, size = 1.5, alpha = 0.65,
             inherit.aes = TRUE) +
  scale_colour_manual(values = SITE_COLS_8) +
  scale_fill_manual(values = SITE_COLS_8) +
  scale_linetype_manual(values = c("Breeding" = "solid", "Haul-out" = "dashed"),
                        name = "Site type") +
  expand_limits(y = 0) +
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(~ Class, ncol = 1, scales = "free_y") +
  labs(x = "Year", y = "Estimated abundance",
       title = "MARSS Smoothed Trajectories — All 8 Sites (1997–2025)",
       subtitle = paste0("Lines = posterior-smoothed; bands = ", CI_LABEL,
                         "; points = observed annual maxima"),
       caption = paste0("Dashed lines = haul-out-only sites (DR, PB); pup states not shown for DR/PB.",
                        " Vertical dashed = 2004 breakpoint.")) +
  theme_seal()

ggsave("Output/Plots/marss8_smoothed_trajectories.jpeg", p.A1,
       width = 28, height = 32, units = "cm", dpi = 200)
cat("Saved: Figure A1 — smoothed trajectories\n")

# ── FIGURE A2: Summed totals across all 8 sites ───────────────────────────────
d8_total <- d8 %>%
  filter(!is.na(estimate)) %>%
  group_by(Year, Class) %>%
  summarise(total      = sum(estimate, na.rm = TRUE),
            total_lo   = sum(lo89_n, na.rm = TRUE),
            total_hi   = sum(hi89_n, na.rm = TRUE),
            .groups = "drop")

d8_obs_total <- d8_obs %>%
  group_by(Year, Class) %>%
  summarise(obs_total = sum(obs, na.rm = TRUE),
            n_sites   = n_distinct(Site), .groups = "drop")

p.A2 <- ggplot(d8_total, aes(x = Year, y = total)) +
  geom_vline(xintercept = bp_year, linetype = 2, colour = "grey40", linewidth = 0.6) +
  geom_ribbon(aes(ymin = total_lo, ymax = total_hi),
              alpha = 0.20, fill = "#2166AC", colour = NA) +
  geom_line(colour = "#2166AC", linewidth = 1.2) +
  geom_point(data = d8_obs_total,
             aes(x = Year, y = obs_total, size = n_sites),
             colour = "#2166AC", shape = 16, alpha = 0.70,
             inherit.aes = FALSE) +
  scale_size_continuous(name = "Sites observed", range = c(1.5, 4), breaks = 1:8) +
  expand_limits(y = 0) +
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(~ Class, ncol = 1, scales = "free_y") +
  labs(x = "Year", y = "Summed estimated abundance (8 sites)",
       title = "Total Estimated Harbor Seal Abundance — All 8 Sites",
       subtitle = paste0("Line = MARSS sum; points = sum of observed annual maxima; bands = ", CI_LABEL),
       caption = "Point size = number of sites observed. Vertical dashed = 2004 breakpoint.") +
  theme_seal() +
  theme(legend.position = "bottom")

ggsave("Output/Plots/marss8_summed_totals.jpeg", p.A2,
       width = 24, height = 32, units = "cm", dpi = 200)
cat("Saved: Figure A2 — summed totals\n")

# ── FIGURE A3: Covariate effects (CIs_8site, Model A) ───────────────────────
C_df <- tryCatch(
  tidy(CIs_8site) %>%
    tibble::as_tibble() %>%
    dplyr::filter(str_detect(term, "^C\\.")) %>%
    mutate(coef_name = str_remove(term, "^C\\.")),
  error = function(e)
    stop("Could not tidy CIs_8site: ", conditionMessage(e))
)

# Annotate with site and covariate type
C_df <- C_df %>%
  mutate(
    significant = (conf.low > 0) | (conf.up < 0),
    site = str_extract(coef_name, "BL|DE|DP|PRH|TB|TP|DR|PB"),
    cov_type = case_when(
      str_detect(coef_name, "MOCI")  ~ "MOCI",
      str_detect(coef_name, "Coy")   ~ "Coyote",
      str_detect(coef_name, "Dist")  ~ "Disturbance",
      str_detect(coef_name, "ES-")   ~ "Elephant seal",
      TRUE ~ "Other"),
    site_type = ifelse(site %in% c("DR","PB"), "Haul-out", "Breeding")
  ) %>%
  filter(!is.na(site))

# Forest plot grouped by covariate type
p.A3 <- ggplot(C_df, aes(x = estimate, y = reorder(coef_name, estimate),
                         colour = site, shape = significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_linerange(aes(xmin = conf.low, xmax = conf.up),
                 linewidth = 0.7, alpha = 0.7) +
  geom_point(size = 2.5) +
  scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16),
                     labels = c("FALSE" = "Includes zero",
                                "TRUE"  = "Excludes zero"),
                     name = NULL) +
  scale_colour_manual(values = SITE_COLS_8, name = "Site") +
  facet_wrap(~ cov_type, ncol = 2, scales = "free") +
  labs(x = "Coefficient estimate (log scale)", y = NULL,
       title = "MARSS Covariate Effects — All 8 Sites",
       subtitle = paste0("Model A (24 independent states); bars = 89% CI; ",
                         "filled = excludes zero"),
       caption = paste0("Haul-out sites (DR, PB) have MOCI and disturbance terms only ",
                        "(coyote/eSeal structurally absent).")) +
  theme_seal(base_size = 12) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = rel(0.95), face = "bold"))

ggsave("Output/Plots/marss8_covariate_effects.jpeg", p.A3,
       width = 30, height = 32, units = "cm", dpi = 200)
cat("Saved: Figure A3 — covariate effects\n")

# ── MODEL COMPARISON TABLE (printed for Appendix A, Table S2) ─────────────────
if (exists("df_aic_8site")) {
  cat("\n── 8-Site Model Comparison (Table S2) ─────────────────────────────\n")
  print(df_aic_8site)
  cat("\nInterpretation:\n")
  cat("  Model A vs C: tests whether Estuary/Coastal grouping explains site dynamics\n")
  cat("  Model A vs B: tests haul-out vs breeding distinction (confounded with covariate pooling)\n")
  cat("  Large deltaAICc for B and C confirms sites are better treated as independent.\n")
}

# ── U parameter trends from best model ───────────────────────────────────────
cat("\n── Phase-specific growth rates (Model A, 89% CI) ──────────────────\n")
u_ci <- tryCatch(
  tidy(CIs_8site) %>%
    filter(grepl("^U\\.", term)) %>%
    distinct(term, .keep_all = TRUE) %>%
    arrange(term),
  error = function(e) NULL
)
if (!is.null(u_ci)) print(u_ci)

cat("\nAll plots -> Output/Plots/marss8_*.jpeg\n")
cat("Next: manuscript Appendix A figures inserted from these outputs.\n")
