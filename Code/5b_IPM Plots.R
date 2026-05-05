##IPM plots

# =============================================================================
# HARBOR SEAL IPM v3.1 — STANDALONE PLOT COMMANDS
# =============================================================================
# Each section is self-contained. Assumes the following objects exist:
#   fit        — CmdStan fitted model object
#   sim_data   — list from prepare_ipm_data_v3.1(), containing $stan_data
#   portfolio  — list from create_portfolio_analysis_v3.1()
#   sync       — list from create_synchrony_projections_v3.1()
#
# Site index:  1=BL, 2=DE, 3=DP, 4=PRH, 5=TB, 6=TP
# To adjust any plot: modify the THEME ADJUSTMENTS block at the bottom
#   of each section, then re-run from "p <- ggplot(...)" onward.
# =============================================================================

library(tidyverse)
library(posterior)
library(bayesplot)
library(patchwork)
library(scales)

# Shared constants — adjust once here, applies everywhere
SITE_NAMES  <- c("BL", "DE", "DP", "PRH", "TB", "TP")
SITE_COLORS <- c("BL"  = "#1b7837",
                 "DE"  = "#762a83",
                 "DP"  = "#d01c8b",
                 "PRH" = "#f1b6da",
                 "TB"  = "#4393c3",
                 "TP"  = "#e08214")
YEARS       <- sim_data$stan_data$years   # numeric year vector
T_          <- length(YEARS)
S_          <- 6

# Pull posterior draws once; reused by all sections
draws <- fit$draws(format = "df")


# =============================================================================
# SECTION 1: POPULATION TRAJECTORIES (observed vs. fitted, by site & class)
# =============================================================================
# Produces one panel per site; rows = Adults / Pups / Molts

# --- Extract ---
N_adult_post <- draws |> select(starts_with("N_adult[")) |>
  pivot_longer(everything(), names_to = "param", values_to = "value") |>
  mutate(site = as.integer(str_extract(param, "(?<=\\[)\\d+")),
         year = as.integer(str_extract(param, "\\d+(?=\\])")),
         Year = YEARS[year], Site = SITE_NAMES[site], Class = "Adult")

N_pup_post <- draws |> select(starts_with("N_pup[")) |>
  pivot_longer(everything(), names_to = "param", values_to = "value") |>
  mutate(site = as.integer(str_extract(param, "(?<=\\[)\\d+")),
         year = as.integer(str_extract(param, "\\d+(?=\\])")),
         Year = YEARS[year], Site = SITE_NAMES[site], Class = "Pup")

N_molt_post <- draws |> select(starts_with("N_molt[")) |>
  pivot_longer(everything(), names_to = "param", values_to = "value") |>
  mutate(site = as.integer(str_extract(param, "(?<=\\[)\\d+")),
         year = as.integer(str_extract(param, "\\d+(?=\\])")),
         Year = YEARS[year], Site = SITE_NAMES[site], Class = "Molt")

traj_post <- bind_rows(N_adult_post, N_pup_post, N_molt_post) |>
  group_by(Site, Year, Class) |>
  summarise(median = median(value),
            lo95   = quantile(value, 0.025),
            hi95   = quantile(value, 0.975),
            lo50   = quantile(value, 0.25),
            hi50   = quantile(value, 0.75), .groups = "drop") |>
  mutate(Class = factor(Class, levels = c("Adult", "Pup", "Molt")),
         Site  = factor(Site,  levels = SITE_NAMES))

# Observed counts (log scale in Stan; back-transform)
obs_adult <- sim_data$stan_data$y_adult
obs_pup   <- sim_data$stan_data$y_pup
obs_molt  <- sim_data$stan_data$y_molt

obs_df <- bind_rows(
  expand.grid(site = 1:S_, year_idx = 1:T_) |>
    mutate(Year = YEARS[year_idx], Site = SITE_NAMES[site],
           Class = "Adult", Observed = exp(obs_adult[cbind(site, year_idx)])),
  expand.grid(site = 1:S_, year_idx = 1:T_) |>
    mutate(Year = YEARS[year_idx], Site = SITE_NAMES[site],
           Class = "Pup",   Observed = exp(obs_pup[cbind(site, year_idx)])),
  expand.grid(site = 1:S_, year_idx = 1:T_) |>
    mutate(Year = YEARS[year_idx], Site = SITE_NAMES[site],
           Class = "Molt",  Observed = exp(obs_molt[cbind(site, year_idx)]))
) |>
  filter(!is.na(Observed)) |>
  mutate(Class = factor(Class, levels = c("Adult", "Pup", "Molt")),
         Site  = factor(Site,  levels = SITE_NAMES))

# --- Plot ---
p_trajectories <- ggplot(traj_post, aes(x = Year)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95, fill = Class), alpha = 0.15) +
  geom_ribbon(aes(ymin = lo50, ymax = hi50, fill = Class), alpha = 0.30) +
  geom_line(aes(y = median, color = Class), linewidth = 0.9) +
  geom_point(data = obs_df, aes(y = Observed, color = Class),
             shape = 21, size = 1.8, fill = "white", stroke = 0.8) +
  facet_grid(Class ~ Site, scales = "free_y") +
  scale_color_manual(values = c(Adult = "#2166ac", Pup = "#d73027", Molt = "#4dac26")) +
  scale_fill_manual( values = c(Adult = "#2166ac", Pup = "#d73027", Molt = "#4dac26")) +
  # ── THEME ADJUSTMENTS ────────────────────────────────────────────────────
  labs(x = "Year", y = "Count", title = "Population Trajectories by Site and Class",
       caption = "Ribbon = 50% and 95% CrI; points = observed counts") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey92"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
# ─────────────────────────────────────────────────────────────────────────

p_trajectories
ggsave("Output/Plots/P01_population_trajectories.jpeg",
       p_trajectories, width = 40, height = 22, units = "cm", dpi = 300)


# =============================================================================
# SECTION 2: VITAL RATE POSTERIORS (baseline, all sites pooled)
# =============================================================================
# Violin / half-eye plot of phi_pup, phi_juv, phi_adult_f, phi_adult_m,
# fecundity by age class, sigma_proc, sigma_obs

vr_params <- c(
  "phi_pup_base"    = "φ[pup]",
  "phi_juv_base"    = "φ[juv]",
  "phi_adult_base"  = "φ[adult, male]",
  "delta_sex"       = "δ[sex]  (f − m logit)",
  "fecundity_prima" = "f[primiparous]",
  "fecundity_multi" = "f[multiparous]",
  "fecundity_prime" = "f[prime]",
  "sigma_process"   = "σ[proc]",
  "sigma_obs"       = "σ[obs]"
)

# phi_adult_female derived: plogis(qlogis(phi_adult_base) + delta_sex)
vr_draws <- draws |>
  select(any_of(names(vr_params))) |>
  mutate(phi_adult_female = plogis(qlogis(phi_adult_base) + delta_sex)) |>
  pivot_longer(everything(), names_to = "param", values_to = "value") |>
  mutate(label = recode(param, !!!vr_params,
                        phi_adult_female = "φ[adult, female]"),
         label = factor(label,
                        levels = c("φ[pup]","φ[juv]",
                                   "φ[adult, male]","φ[adult, female]",
                                   "δ[sex]  (f − m logit)",
                                   "f[primiparous]","f[multiparous]","f[prime]",
                                   "σ[proc]","σ[obs]")))

p_vital_rates <- ggplot(vr_draws, aes(x = value, y = label)) +
  stat_halfeye(point_interval = median_qi, .width = c(0.5, 0.95),
               fill = "#4393c3", color = "#08306b", alpha = 0.8,
               normalize = "panels") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  # ── THEME ADJUSTMENTS ────────────────────────────────────────────────────
  labs(x = "Parameter value", y = NULL,
       title = "IPM v3.1 — Baseline Vital Rate Posteriors",
       subtitle = "Point = median; thick bar = 50% CrI; thin bar = 95% CrI") +
  theme_bw(base_size = 13) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.y = element_text(size = 11))
# ─────────────────────────────────────────────────────────────────────────

p_vital_rates
ggsave("Output/Plots/P02_vital_rate_posteriors.jpeg",
       p_vital_rates, width = 22, height = 18, units = "cm", dpi = 300)


# =============================================================================
# SECTION 3: COVARIATE EFFECTS — FOREST PLOT (all betas, by vital rate)
# =============================================================================

cov_params <- tribble(
  ~param,                  ~label,                         ~vital_rate,
  "beta_coy_BL",           "Coyote — BL",                  "Pup survival",
  "beta_coy_DE",           "Coyote — DE",                  "Pup survival",
  "beta_coy_DP",           "Coyote — DP",                  "Pup survival",
  "beta_dist_surv_1",      "Disturbance — BL",             "Pup survival",
  "beta_dist_surv_2",      "Disturbance — DE",             "Pup survival",
  "beta_dist_surv_3",      "Disturbance — DP",             "Pup survival",
  "beta_dist_surv_4",      "Disturbance — PRH",            "Pup survival",
  "beta_dist_surv_5",      "Disturbance — TB",             "Pup survival",
  "beta_dist_surv_6",      "Disturbance — TP",             "Pup survival",
  "beta_moci_amj_pup",     "MOCI AMJ lag",                 "Pup survival",
  "beta_moci_ond_pup",     "MOCI OND lag",                 "Pup survival",
  "beta_eseal_pup",        "Elephant seal",                "Pup survival",
  "beta_moci_jfm_adult",   "MOCI JFM lag",                 "Adult survival",
  "beta_eseal_adult",      "Elephant seal",                "Adult survival",
  "beta_moci_amj_fec",     "MOCI AMJ lag",                 "Fecundity"
)

cov_draws <- draws |>
  select(any_of(cov_params$param)) |>
  pivot_longer(everything(), names_to = "param", values_to = "value") |>
  left_join(cov_params, by = "param") |>
  group_by(param, label, vital_rate) |>
  summarise(median = median(value),
            lo95   = quantile(value, 0.025),
            hi95   = quantile(value, 0.975),
            lo50   = quantile(value, 0.25),
            hi50   = quantile(value, 0.75),
            p_neg  = mean(value < 0), .groups = "drop") |>
  mutate(sig = ifelse((lo95 < 0 & hi95 < 0) | (lo95 > 0 & hi95 > 0),
                      "95% CI ≠ 0", "CI spans 0"),
         vital_rate = factor(vital_rate,
                             levels = c("Pup survival","Adult survival","Fecundity")))

p_covariate_effects <- ggplot(cov_draws,
                              aes(x = median, y = reorder(label, median),
                                  color = sig, shape = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = lo95, xmax = hi95), height = 0, linewidth = 0.6) +
  geom_errorbarh(aes(xmin = lo50, xmax = hi50), height = 0, linewidth = 1.5) +
  geom_point(size = 2.8) +
  facet_wrap(~ vital_rate, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("95% CI ≠ 0" = "#b2182b", "CI spans 0" = "#4393c3")) +
  scale_shape_manual(values = c("95% CI ≠ 0" = 16, "CI spans 0" = 21)) +
  # ── THEME ADJUSTMENTS ────────────────────────────────────────────────────
  labs(x = "Coefficient (logit or log scale)", y = NULL,
       color = NULL, shape = NULL,
       title = "IPM v3.1 — Covariate Effects on Vital Rates",
       subtitle = "Thick bar = 50% CrI; thin bar = 95% CrI") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top",
        strip.background = element_rect(fill = "grey92"),
        panel.grid.major.x = element_line(color = "grey90"),
        panel.grid.minor   = element_blank())
# ─────────────────────────────────────────────────────────────────────────

p_covariate_effects
ggsave("Output/Plots/P03_covariate_effects.jpeg",
       p_covariate_effects, width = 20, height = 28, units = "cm", dpi = 300)


# =============================================================================
# SECTION 4: LAMBDA OVER TIME — AGGREGATE AND SITE-SPECIFIC
# =============================================================================

lambda_post <- draws |> select(starts_with("lambda[")) |>
  pivot_longer(everything(), names_to = "param", values_to = "value") |>
  mutate(site     = as.integer(str_extract(param, "(?<=\\[)\\d+")),
         year_idx = as.integer(str_extract(param, "\\d+(?=\\])")),
         Year     = YEARS[year_idx],
         Site     = SITE_NAMES[site]) |>
  group_by(Site, Year) |>
  summarise(median = median(value),
            lo95   = quantile(value, 0.025),
            hi95   = quantile(value, 0.975),
            lo50   = quantile(value, 0.25),
            hi50   = quantile(value, 0.75), .groups = "drop") |>
  mutate(Site = factor(Site, levels = SITE_NAMES))

lambda_agg <- draws |> select(starts_with("lambda[")) |>
  pivot_longer(everything(), names_to = "param", values_to = "value") |>
  mutate(year_idx = as.integer(str_extract(param, "\\d+(?=\\])")),
         Year     = YEARS[year_idx]) |>
  group_by(Year, .draw = .draw) |>
  summarise(agg_lambda = mean(value), .groups = "drop") |>
  group_by(Year) |>
  summarise(median = median(agg_lambda),
            lo95   = quantile(agg_lambda, 0.025),
            hi95   = quantile(agg_lambda, 0.975),
            lo50   = quantile(agg_lambda, 0.25),
            hi50   = quantile(agg_lambda, 0.75), .groups = "drop")

# 4a — Site-specific lambda
p_lambda_site <- ggplot(lambda_post, aes(x = Year)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_ribbon(aes(ymin = lo95, ymax = hi95, fill = Site), alpha = 0.15) +
  geom_ribbon(aes(ymin = lo50, ymax = hi50, fill = Site), alpha = 0.30) +
  geom_line(aes(y = median, color = Site), linewidth = 0.9) +
  facet_wrap(~ Site, ncol = 3) +
  scale_color_manual(values = SITE_COLORS) +
  scale_fill_manual( values = SITE_COLORS) +
  # ── THEME ADJUSTMENTS ────────────────────────────────────────────────────
  labs(x = "Year", y = expression(lambda ~ "(population growth rate)"),
       title = "Site-Specific Population Growth Rate Over Time") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey92"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
# ─────────────────────────────────────────────────────────────────────────

p_lambda_site
ggsave("Output/Plots/P04a_lambda_by_site.jpeg",
       p_lambda_site, width = 30, height = 18, units = "cm", dpi = 300)

# 4b — Aggregate lambda
p_lambda_agg <- ggplot(lambda_agg, aes(x = Year)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "#4393c3", alpha = 0.20) +
  geom_ribbon(aes(ymin = lo50, ymax = hi50), fill = "#4393c3", alpha = 0.40) +
  geom_line(aes(y = median), color = "#08306b", linewidth = 1.1) +
  # ── THEME ADJUSTMENTS ────────────────────────────────────────────────────
  labs(x = "Year", y = expression(lambda ~ "(aggregate)"),
       title = "Aggregate Population Growth Rate",
       subtitle = "Mean across all 6 sites; ribbon = 50% and 95% CrI") +
  scale_y_continuous(limits = c(NA, NA), oob = scales::squish) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())
# ─────────────────────────────────────────────────────────────────────────

p_lambda_agg
ggsave("Output/Plots/P04b_lambda_aggregate.jpeg",
       p_lambda_agg, width = 22, height = 13, units = "cm", dpi = 300)


# =============================================================================
# SECTION 5: ELASTICITY ANALYSIS
# =============================================================================
# Assumes fit$summary() includes "elasticity_*" generated quantities;
# if not, compute from posterior lambda and vital rate perturbation.

elast_params <- c(
  "elasticity_phi_pup"   = "φ[pup]",
  "elasticity_phi_juv"   = "φ[juv]",
  "elasticity_phi_adult" = "φ[adult]",
  "elasticity_fecundity" = "Fecundity"
)

elast_draws <- draws |>
  select(any_of(names(elast_params))) |>
  pivot_longer(everything(), names_to = "param", values_to = "value") |>
  mutate(label = recode(param, !!!elast_params),
         label = factor(label, levels = rev(elast_params))) |>
  group_by(label) |>
  summarise(median = median(value),
            lo95   = quantile(value, 0.025),
            hi95   = quantile(value, 0.975),
            lo50   = quantile(value, 0.25),
            hi50   = quantile(value, 0.75), .groups = "drop")

p_elasticity <- ggplot(elast_draws, aes(x = median, y = label)) +
  geom_errorbarh(aes(xmin = lo95, xmax = hi95), height = 0, linewidth = 0.7) +
  geom_errorbarh(aes(xmin = lo50, xmax = hi50), height = 0, linewidth = 2.0,
                 color = "#2166ac") +
  geom_point(size = 3.5, color = "#08306b") +
  # ── THEME ADJUSTMENTS ────────────────────────────────────────────────────
  labs(x = "Elasticity (proportional change in λ per proportional change in vital rate)",
       y = NULL, title = "Elasticity of λ to Vital Rates") +
  theme_bw(base_size = 13) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank())
# ─────────────────────────────────────────────────────────────────────────

p_elasticity
ggsave("Output/Plots/P05_elasticity.jpeg",
       p_elasticity, width = 22, height = 12, units = "cm", dpi = 300)


# =============================================================================
# SECTION 6: SEX RATIO TREND OVER TIME
# =============================================================================
# Derived from phi_adult_base + delta_sex; shows convergence to stable ratio.

sex_ratio_post <- draws |>
  select(phi_adult_base, delta_sex, starts_with("N_adult[")) |>
  # Quick computation: stable sex ratio = phi_f / (phi_f + phi_m)
  mutate(phi_female  = plogis(qlogis(phi_adult_base) + delta_sex),
         phi_male    = phi_adult_base,   # already on probability scale from Stan
         stable_ratio = phi_female / (phi_female + phi_male))

stable_ratio_summary <- sex_ratio_post |>
  summarise(median = median(stable_ratio),
            lo95   = quantile(stable_ratio, 0.025),
            hi95   = quantile(stable_ratio, 0.975))

# Observed sex ratio from count data (requires sex-disaggregated counts if available;
# otherwise use sim_data covariate or external field data):
obs_sex_ratio <- sim_data$stan_data$sex_ratio_obs   # matrix [S x T] or data frame
# If not in stan_data, load from your field data separately:
# obs_sex_ratio <- read_csv("Data/harbor_seal_sex_ratio.csv")

obs_sr_df <- if (exists("obs_sex_ratio") && !is.null(obs_sex_ratio)) {
  as.data.frame(obs_sex_ratio) |>
    rownames_to_column("site_idx") |>
    pivot_longer(-site_idx, names_to = "year_idx", values_to = "prop_female") |>
    mutate(site     = as.integer(site_idx),
           year_idx = as.integer(str_extract(year_idx, "\\d+")),
           Year     = YEARS[year_idx],
           Site     = SITE_NAMES[site]) |>
    filter(!is.na(prop_female))
} else {
  # Build from counts directly: n_female / (n_female + n_male)
  # Placeholder — replace with actual sex-disaggregated count data
  tibble(Year = YEARS, prop_female = NA_real_, Site = "All sites")
}

# Aggregate across sites, plot trend + stable ratio reference line
p_sex_ratio <- ggplot(obs_sr_df |> filter(!is.na(prop_female)),
                      aes(x = Year, y = prop_female)) +
  # Stable sex ratio posterior band
  geom_hline(yintercept = stable_ratio_summary$median,
             linetype = "dashed", color = "#b2182b", linewidth = 0.9) +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = stable_ratio_summary$lo95,
           ymax = stable_ratio_summary$hi95,
           fill = "#b2182b", alpha = 0.08) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey50") +
  geom_line(aes(color = Site), linewidth = 0.8) +
  geom_point(aes(color = Site), size = 2) +
  scale_color_manual(values = SITE_COLORS) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0.35, 0.70)) +
  # ── THEME ADJUSTMENTS ────────────────────────────────────────────────────
  labs(x = "Year", y = "Proportion female (adults)",
       color = "Site",
       title = "Adult Sex Ratio Trend — Convergence to Stable Ratio",
       subtitle = "Dashed line = IPM-estimated stable sex ratio (median + 95% CrI shaded);\ndotted = 50:50 reference") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "right")
# ─────────────────────────────────────────────────────────────────────────

p_sex_ratio
ggsave("Output/Plots/P06_sex_ratio_trend.jpeg",
       p_sex_ratio, width = 24, height = 14, units = "cm", dpi = 300)


# =============================================================================
# SECTION 7: POPULATION PROJECTIONS — SCENARIO COMPARISON
# =============================================================================
# Uses projection output stored in results from 10-year forward simulation.
# Assumes: proj_results = list with $status_quo, $coyote_removal, $warm_moci, $combined
# Each element: data frame with columns Year, median, lo50, hi50, lo95, hi95

proj_scenarios <- bind_rows(
  proj_results$status_quo      |> mutate(Scenario = "Status Quo"),
  proj_results$coyote_removal  |> mutate(Scenario = "Coyote Removal"),
  proj_results$warm_moci       |> mutate(Scenario = "Warm MOCI (+1 SD)"),
  proj_results$combined        |> mutate(Scenario = "Combined")
) |>
  mutate(Scenario = factor(Scenario, levels = c("Status Quo","Coyote Removal",
                                                "Warm MOCI (+1 SD)","Combined")))

scenario_colors <- c("Status Quo"        = "#4393c3",
                     "Coyote Removal"    = "#2ca25f",
                     "Warm MOCI (+1 SD)" = "#f03b20",
                     "Combined"          = "#9e4f8a")

proj_hist <- traj_post |>
  filter(Class == "Adult") |>
  group_by(Year) |>
  summarise(median = sum(median), .groups = "drop")

p_projections <- ggplot(proj_scenarios, aes(x = Year, color = Scenario, fill = Scenario)) +
  geom_line(data = proj_hist, aes(y = median), color = "grey30",
            linewidth = 0.9, linetype = "solid", inherit.aes = FALSE) +
  geom_vline(xintercept = max(YEARS), linetype = "dotted", color = "grey50") +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.10, color = NA) +
  geom_ribbon(aes(ymin = lo50, ymax = hi50), alpha = 0.25, color = NA) +
  geom_line(aes(y = median), linewidth = 1.0) +
  scale_color_manual(values = scenario_colors) +
  scale_fill_manual( values = scenario_colors) +
  # ── THEME ADJUSTMENTS ────────────────────────────────────────────────────
  labs(x = "Year", y = "Total adult abundance (all sites)",
       color = "Scenario", fill = "Scenario",
       title = "10-Year Population Projections Under Management Scenarios",
       subtitle = "Grey line = historical posterior median; ribbon = 50% and 95% CrI") +
  theme_bw(base_size = 13) +
  theme(legend.position  = "right",
        panel.grid.minor = element_blank())
# ─────────────────────────────────────────────────────────────────────────

p_projections
ggsave("Output/Plots/P07_projections_scenarios.jpeg",
       p_projections, width = 28, height = 16, units = "cm", dpi = 300)


# =============================================================================
# SECTION 8: POSTERIOR PREDICTIVE CHECK (observed vs. replicated counts)
# =============================================================================

y_adult_rep <- draws |> select(starts_with("y_adult_rep[")) |>
  pivot_longer(everything(), values_to = "value") |>
  mutate(site     = as.integer(str_extract(name, "(?<=\\[)\\d+")),
         year_idx = as.integer(str_extract(name, "\\d+(?=\\])")),
         Year     = YEARS[year_idx],
         Site     = SITE_NAMES[site])

y_adult_rep_summary <- y_adult_rep |>
  group_by(Site, Year) |>
  summarise(median = median(value),
            lo95   = quantile(value, 0.025),
            hi95   = quantile(value, 0.975), .groups = "drop")

obs_adult_df <- obs_df |> filter(Class == "Adult") |>
  mutate(log_obs = log(Observed))

p_ppc_adult <- ggplot(y_adult_rep_summary, aes(x = Year)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "#4393c3", alpha = 0.25) +
  geom_line(aes(y = median), color = "#2166ac", linewidth = 0.8) +
  geom_point(data = obs_adult_df, aes(y = log_obs),
             color = "#b2182b", size = 2, shape = 16) +
  facet_wrap(~ Site, ncol = 3, scales = "free_y") +
  # ── THEME ADJUSTMENTS ────────────────────────────────────────────────────
  labs(x = "Year", y = "log(Adult count)",
       title = "Posterior Predictive Check — Adult Counts",
       subtitle = "Blue = replicated (median + 95% CrI); red points = observed") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey92"),
        panel.grid.minor  = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
# ─────────────────────────────────────────────────────────────────────────

p_ppc_adult
ggsave("Output/Plots/P08_ppc_adults.jpeg",
       p_ppc_adult, width = 30, height = 18, units = "cm", dpi = 300)


# =============================================================================
# SECTION 9: PORTFOLIO — LAMBDA HEATMAP (site × year)
# =============================================================================
# Replot from portfolio$lambda_matrix (already computed inside portfolio function)

lambda_heat_df <- expand.grid(site = 1:S_, year_idx = 1:(T_-1)) |>
  mutate(Site     = factor(SITE_NAMES[site], levels = SITE_NAMES),
         Year     = YEARS[year_idx],
         Lambda   = as.vector(portfolio$lambda_matrix))

p_lambda_heatmap <- ggplot(lambda_heat_df, aes(x = Year, y = Site, fill = Lambda)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = "#d73027", mid = "white", high = "#1a9641",
                       midpoint = 1, name = expression(lambda),
                       limits = c(0.7, 1.4), oob = scales::squish) +
  # ── THEME ADJUSTMENTS ────────────────────────────────────────────────────
  labs(x = "Year", y = "Site",
       title = expression("Site × Year Population Growth Rate (" * lambda * ")"),
       subtitle = "Green > 1 (growing); red < 1 (declining)") +
  theme_bw(base_size = 12) +
  theme(panel.grid   = element_blank(),
        axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "right")
# ─────────────────────────────────────────────────────────────────────────

p_lambda_heatmap
ggsave("Output/Plots/P09_portfolio_lambda_heatmap.jpeg",
       p_lambda_heatmap, width = 35, height = 12, units = "cm", dpi = 300)


# =============================================================================
# SECTION 10: PORTFOLIO — PUP SURVIVAL HEATMAP (site × year)
# =============================================================================

phi_heat_df <- expand.grid(site = 1:S_, year_idx = 1:T_) |>
  mutate(Site       = factor(SITE_NAMES[site], levels = SITE_NAMES),
         Year       = YEARS[year_idx],
         Phi_pup    = as.vector(portfolio$phi_matrix))

p_phi_heatmap <- ggplot(phi_heat_df, aes(x = Year, y = Site, fill = Phi_pup)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_viridis_c(name = expression(phi[pup]), option = "plasma",
                       limits = c(0.1, 0.9), oob = scales::squish) +
  # ── THEME ADJUSTMENTS ────────────────────────────────────────────────────
  labs(x = "Year", y = "Site",
       title = "Site × Year Pup Survival Probability",
       subtitle = "Darker = lower pup survival") +
  theme_bw(base_size = 12) +
  theme(panel.grid   = element_blank(),
        axis.text.x  = element_text(angle = 45, hjust = 1, size = 9))
# ─────────────────────────────────────────────────────────────────────────

p_phi_heatmap
ggsave("Output/Plots/P10_portfolio_phi_heatmap.jpeg",
       p_phi_heatmap, width = 35, height = 12, units = "cm", dpi = 300)


# =============================================================================
# SECTION 11: PORTFOLIO — SITE TRAJECTORIES (asynchrony)
# =============================================================================

async_df <- lambda_heat_df |>
  group_by(Year) |>
  mutate(Mean_Lambda = mean(Lambda, na.rm = TRUE),
         Deviation   = Lambda - Mean_Lambda)

p_asynchrony <- ggplot(lambda_heat_df, aes(x = Year, y = Lambda,
                                           color = Site, group = Site)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  geom_point(size = 1.4) +
  scale_color_manual(values = SITE_COLORS) +
  # ── THEME ADJUSTMENTS ────────────────────────────────────────────────────
  labs(x = "Year", y = expression(lambda), color = "Site",
       title = "Site-Specific Population Growth Rates — Portfolio Asynchrony",
       subtitle = "Lines crossing and diverging from each other = portfolio effect") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))
# ─────────────────────────────────────────────────────────────────────────

p_asynchrony
ggsave("Output/Plots/P11_portfolio_asynchrony.jpeg",
       p_asynchrony, width = 28, height = 14, units = "cm", dpi = 300)


# =============================================================================
# SECTION 12: PORTFOLIO — BEST / WORST SITE PER YEAR
# =============================================================================
# Directly reuse portfolio$best_worst if available, or rebuild:

best_worst_df <- lambda_heat_df |>
  group_by(Year) |>
  summarise(
    Best_Site    = Site[which.max(Lambda)],
    Best_Lambda  = max(Lambda, na.rm = TRUE),
    Worst_Site   = Site[which.min(Lambda)],
    Worst_Lambda = min(Lambda, na.rm = TRUE), .groups = "drop"
  ) |>
  pivot_longer(cols = c(Best_Site, Worst_Site),
               names_to = "Performance", values_to = "Site") |>
  mutate(Performance = recode(Performance,
                              Best_Site  = "Best",
                              Worst_Site = "Worst"),
         Site = factor(Site, levels = SITE_NAMES))

p_best_worst <- ggplot(best_worst_df, aes(x = Year, y = Site, fill = Performance)) +
  geom_tile(color = "white", linewidth = 0.8, alpha = 0.85) +
  scale_fill_manual(values = c("Best" = "#2ca25f", "Worst" = "#de2d26")) +
  # ── THEME ADJUSTMENTS ────────────────────────────────────────────────────
  labs(x = "Year", y = "Site", fill = NULL,
       title = "Best and Worst Performing Sites Each Year",
       subtitle = "No single site consistently worst → portfolio insurance") +
  theme_bw(base_size = 12) +
  theme(panel.grid  = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "top")
# ─────────────────────────────────────────────────────────────────────────

p_best_worst
ggsave("Output/Plots/P12_portfolio_best_worst.jpeg",
       p_best_worst, width = 35, height = 12, units = "cm", dpi = 300)


# =============================================================================
# SECTION 13: PORTFOLIO — SITE CORRELATION MATRIX
# =============================================================================

lambda_wide <- lambda_heat_df |>
  select(Site, Year, Lambda) |>
  pivot_wider(names_from = Site, values_from = Lambda) |>
  select(-Year)

lambda_cor <- cor(lambda_wide, use = "pairwise.complete.obs")

cor_df <- as.data.frame(as.table(lambda_cor)) |>
  rename(Site1 = Var1, Site2 = Var2, correlation = Freq) |>
  mutate(Site1 = factor(Site1, levels = SITE_NAMES),
         Site2 = factor(Site2, levels = rev(SITE_NAMES)))

p_correlation <- ggplot(cor_df, aes(x = Site1, y = Site2, fill = correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 4.5) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-1, 1),
                       name = "r") +
  coord_fixed() +
  # ── THEME ADJUSTMENTS ────────────────────────────────────────────────────
  labs(x = NULL, y = NULL,
       title = "Between-Site Correlation in Population Growth (λ)",
       subtitle = "Low / negative correlations = asynchrony = stronger portfolio effect") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        axis.text  = element_text(size = 12))
# ─────────────────────────────────────────────────────────────────────────

p_correlation
ggsave("Output/Plots/P13_portfolio_correlation.jpeg",
       p_correlation, width = 20, height = 18, units = "cm", dpi = 300)


# =============================================================================
# SECTION 14: PORTFOLIO — DRIVER DECOMPOSITION
# =============================================================================
# portfolio$contributions has: Site, Driver, Contribution (% of total variance)

p_drivers <- ggplot(portfolio$contributions,
                    aes(x = Site, y = Contribution, fill = Driver)) +
  geom_col(position = "stack", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("MOCI"        = "#d73027",
                               "Coyote"      = "#9970ab",
                               "Disturbance" = "#4393c3",
                               "Residual"    = "grey70")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  # ── THEME ADJUSTMENTS ────────────────────────────────────────────────────
  labs(x = "Site", y = "% contribution to between-site λ variance",
       fill = "Driver",
       title = "Driver Decomposition of Site-to-Site Asynchrony",
       subtitle = "MOCI synchronizes; coyote and disturbance desynchronize") +
  theme_bw(base_size = 13) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank())
# ─────────────────────────────────────────────────────────────────────────

p_drivers
ggsave("Output/Plots/P14_portfolio_driver_decomposition.jpeg",
       p_drivers, width = 22, height = 14, units = "cm", dpi = 300)


# =============================================================================
# SECTION 15: SYNCHRONY PROJECTIONS — ASYNC VS SYNC COMPARISON
# =============================================================================
# sync$comparison_data has: Year, Scenario, Type (Async/Sync), median, lo95, hi95

p_sync_comparison <- ggplot(sync$projection_data,
                            aes(x = Year, color = Scenario, linetype = Type)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95, fill = Scenario),
              alpha = 0.10, color = NA) +
  geom_line(aes(y = median), linewidth = 0.9) +
  scale_linetype_manual(values = c("Async" = "solid", "Sync" = "dashed")) +
  scale_color_manual(values = scenario_colors) +
  scale_fill_manual( values = scenario_colors) +
  facet_wrap(~ Scenario, ncol = 2) +
  # ── THEME ADJUSTMENTS ────────────────────────────────────────────────────
  labs(x = "Year", y = "Total abundance",
       color = "Scenario", linetype = "Dynamics",
       title = "Portfolio Buffering: Asynchronous vs. Synchronous Site Dynamics",
       subtitle = "Solid = observed inter-site correlation; dashed = full synchrony (φ = 1)") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        strip.background = element_rect(fill = "grey92"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
# ─────────────────────────────────────────────────────────────────────────

p_sync_comparison
ggsave("Output/Plots/P15_synchrony_comparison.jpeg",
       p_sync_comparison, width = 30, height = 22, units = "cm", dpi = 300)


# =============================================================================
# SECTION 16: SYNCHRONY PROJECTIONS — CV BUFFERING BAR CHART
# =============================================================================
# sync$cv_comparison has: Scenario, CV_Async, CV_Sync, CV_Ratio, Buffering_Pct

cv_plot_df <- sync$cv_comparison |>
  pivot_longer(cols = c(CV_Async, CV_Sync),
               names_to = "Type", values_to = "CV") |>
  mutate(Type     = recode(Type, CV_Async = "Asynchronous", CV_Sync = "Synchronous"),
         Scenario = factor(Scenario, levels = c("Status Quo","Coyote Removal",
                                                "Warm MOCI","Combined")))

p_cv_buffering <- ggplot(cv_plot_df, aes(x = Scenario, y = CV, fill = Type)) +
  geom_col(position = "dodge", color = "white", linewidth = 0.3, width = 0.65) +
  geom_text(data = sync$cv_comparison,
            aes(x = Scenario, y = pmax(CV_Async, CV_Sync) + 0.005,
                label = sprintf("Buffer: %.1f%%", Buffering_Pct)),
            inherit.aes = FALSE, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("Asynchronous" = "#4393c3", "Synchronous" = "#d73027")) +
  # ── THEME ADJUSTMENTS ────────────────────────────────────────────────────
  labs(x = "Scenario", y = "CV of aggregate abundance",
       fill = "Site dynamics",
       title = "Portfolio Buffering Effect by Scenario",
       subtitle = "Lower CV = more stable; buffering % = reduction from synchrony to asynchrony") +
  theme_bw(base_size = 13) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1))
# ─────────────────────────────────────────────────────────────────────────

p_cv_buffering
ggsave("Output/Plots/P16_cv_buffering_by_scenario.jpeg",
       p_cv_buffering, width = 22, height = 14, units = "cm", dpi = 300)


# =============================================================================
# SECTION 17: MCMC DIAGNOSTICS — TRACE PLOTS (key parameters)
# =============================================================================
# Using bayesplot; easy to add/remove parameters from the vector below.

diag_pars_vital <- c("phi_pup_base","phi_juv_base","phi_adult_base","delta_sex",
                     "fecundity_base","sigma_process","sigma_obs")

diag_pars_cov   <- c("beta_coy_BL","beta_coy_DE","beta_coy_DP",
                     "beta_moci_amj_pup","beta_moci_ond_pup","beta_moci_jfm_adult",
                     "beta_eseal_pup","beta_eseal_adult")

draws_array <- fit$draws()  # 3D array for bayesplot

# 17a — Vital rates traces
p_trace_vital <- mcmc_trace(draws_array, pars = diag_pars_vital) +
  # ── THEME ADJUSTMENTS ────────────────────────────────────────────────────
  labs(title = "MCMC Trace Plots — Vital Rate Parameters") +
  theme_bw(base_size = 10)
# ─────────────────────────────────────────────────────────────────────────
p_trace_vital
ggsave("Output/Plots/P17a_trace_vital_rates.jpeg",
       p_trace_vital, width = 30, height = 22, units = "cm", dpi = 200)

# 17b — Covariate traces
p_trace_cov <- mcmc_trace(draws_array, pars = diag_pars_cov) +
  labs(title = "MCMC Trace Plots — Covariate Parameters") +
  theme_bw(base_size = 10)
p_trace_cov
ggsave("Output/Plots/P17b_trace_covariates.jpeg",
       p_trace_cov, width = 30, height = 22, units = "cm", dpi = 200)


# =============================================================================
# SECTION 18: PRIOR vs. POSTERIOR COMPARISON (key parameters)
# =============================================================================

prior_post_params <- c("phi_pup_base","phi_adult_base","delta_sex",
                       "beta_coy_BL","beta_coy_DE","beta_moci_amj_pup","sigma_process")

# Sample from priors (must match Stan prior specs)
set.seed(42)
n_prior <- 4000
prior_draws <- tibble(
  phi_pup_base     = plogis(rnorm(n_prior, 0,    1)),
  phi_adult_base   = plogis(rnorm(n_prior, 2.5,  0.5)),
  delta_sex        = rnorm(n_prior,  0.3,  0.2),
  beta_coy_BL      = rnorm(n_prior, -0.5,  0.3),
  beta_coy_DE      = rnorm(n_prior, -0.5,  0.3),
  beta_moci_amj_pup= rnorm(n_prior, -0.3,  0.2),
  sigma_process    = abs(rnorm(n_prior, 0,   0.3))
) |> pivot_longer(everything(), names_to = "param", values_to = "value") |>
  mutate(distribution = "Prior")

post_draws_pp <- draws |>
  select(all_of(prior_post_params)) |>
  pivot_longer(everything(), names_to = "param", values_to = "value") |>
  mutate(distribution = "Posterior")

pp_combined <- bind_rows(prior_draws, post_draws_pp) |>
  mutate(distribution = factor(distribution, levels = c("Prior","Posterior")))

p_prior_post <- ggplot(pp_combined, aes(x = value, fill = distribution,
                                        color = distribution)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  facet_wrap(~ param, scales = "free", ncol = 3) +
  scale_fill_manual( values = c("Prior" = "grey60", "Posterior" = "#2166ac")) +
  scale_color_manual(values = c("Prior" = "grey40", "Posterior" = "#08306b")) +
  # ── THEME ADJUSTMENTS ────────────────────────────────────────────────────
  labs(x = "Parameter value", y = "Density",
       fill = NULL, color = NULL,
       title = "Prior vs. Posterior — Key Parameters") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "top",
        strip.background = element_rect(fill = "grey92"),
        panel.grid.minor = element_blank())
# ─────────────────────────────────────────────────────────────────────────

p_prior_post
ggsave("Output/Plots/P18_prior_posterior.jpeg",
       p_prior_post, width = 30, height = 22, units = "cm", dpi = 300)


# =============================================================================
# SECTION 19: PATCHWORK COMPOSITES (publication-quality multi-panel figures)
# =============================================================================
# Mix and match individual plots defined above.

# --- Figure 1 (main text): Trajectories + Lambda ----
fig1 <- p_trajectories / p_lambda_agg +
  plot_annotation(tag_levels = "A",
                  title = "Harbor Seal Population Dynamics, PRNS 1997–2023") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave("Output/Plots/FIG1_trajectories_lambda.jpeg",
       fig1, width = 40, height = 30, units = "cm", dpi = 300)

# --- Figure 2 (main text): Covariate effects + Vital rates ----
fig2 <- (p_vital_rates | p_covariate_effects) +
  plot_layout(widths = c(1, 1.4)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave("Output/Plots/FIG2_vital_rates_covariates.jpeg",
       fig2, width = 40, height = 22, units = "cm", dpi = 300)

# --- Figure 3 (main text): Lambda + Elasticity + Sex ratio ----
fig3 <- (p_lambda_site / (p_elasticity | p_sex_ratio)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave("Output/Plots/FIG3_lambda_elasticity_sexratio.jpeg",
       fig3, width = 40, height = 28, units = "cm", dpi = 300)

# --- Figure 4 (main text): Projections ----
ggsave("Output/Plots/FIG4_projections.jpeg",
       p_projections, width = 28, height = 16, units = "cm", dpi = 300)

# --- Figure 5 (main text): Portfolio 4-panel ----
fig5 <- (p_lambda_heatmap / p_best_worst) |
  (p_correlation / p_drivers) +
  plot_annotation(tag_levels = "A",
                  title = "Portfolio Analysis — Spatial Insurance across PRNS Haul-out Sites") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave("Output/Plots/FIG5_portfolio.jpeg",
       fig5, width = 45, height = 28, units = "cm", dpi = 300)

# --- Figure 6 (main text): Synchrony ----
fig6 <- (p_sync_comparison / p_cv_buffering) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave("Output/Plots/FIG6_synchrony.jpeg",
       fig6, width = 30, height = 28, units = "cm", dpi = 300)

message("All plots saved to Output/Plots/")