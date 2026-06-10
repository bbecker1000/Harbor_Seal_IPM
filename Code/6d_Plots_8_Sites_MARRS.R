# ============================================================================
# PLOTS — 8-SITE MARSS ANALYSIS
# Aligned with IPM v3.2 style (89% CI, consistent palette, y-axis from 0)
#
# Source first:
#   source("4a_data_prep_8sites_MARSS.R")
#   source("4b_covariates_8sites_MARSS.R")
#   source("4cc_models_8sites_MARSS.R")   # or load saved objects:
#   load("Output/m.A_8site.RData")
#   load("Output/CIs_8site.RData")
# ============================================================================

library(MARSS)
library(tidyverse)
library(ggplot2)
library(patchwork)

# ── Constants matching IPM v3.2 ───────────────────────────────────────────────
CI_LO    <- 0.055    # 89% CI lower
CI_HI    <- 0.945    # 89% CI upper
CI_LABEL <- "89% CI"

SITE_TYPE <- c(BL = "Breeding", DE = "Breeding", DP = "Breeding",
               DR = "Haul-out only", PB = "Haul-out only",
               PRH = "Breeding", TB = "Breeding", TP = "Breeding")

SITE_COLS <- c(
  BL = "#E41A1C", DE = "#FF7F00", DP = "#4DAF4A",
  DR = "#984EA3", PB = "#00CED1",              # haul-out sites: distinct colours
  PRH = "#377EB8", TB = "#A65628", TP = "#F781BF"
)

SITE_LTY <- c(BL=1, DE=1, DP=1, DR=2, PB=2, PRH=1, TB=1, TP=1)  # dashed for haul-out

theme_seal <- function(base_size = 14) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.grid.major   = element_line(colour = "grey88", linewidth = 0.4),
      panel.border       = element_rect(colour = "grey70", fill = NA, linewidth = 0.5),
      axis.title         = element_text(size = rel(0.95)),
      legend.position    = "bottom",
      strip.text         = element_text(size = rel(0.90), face = "bold"),
      strip.background   = element_rect(fill = "grey94", colour = "grey80"),
      plot.title         = element_text(size = rel(1.05), face = "bold"),
      plot.caption       = element_text(size = rel(0.78), colour = "grey50", hjust = 1),
      plot.margin        = margin(10, 14, 10, 10)
    )
}

dir.create("Output/Plots", showWarnings = FALSE)

BESTMODEL  <- m.A_8site
years      <- years_8site
state_names_8site <- c(
  "BL_Adult",  "BL_Molt",  "BL_Pup",
  "DE_Adult",  "DE_Molt",  "DE_Pup",
  "DP_Adult",  "DP_Molt",  "DP_Pup",
  "DR_Adult",  "DR_Molt",  "DR_Pup",
  "PB_Adult",  "PB_Molt",  "PB_Pup",
  "PRH_Adult", "PRH_Molt", "PRH_Pup",
  "TB_Adult",  "TB_Molt",  "TB_Pup",
  "TP_Adult",  "TP_Molt",  "TP_Pup"
)

# ── Helper: build long-form state data ────────────────────────────────────────
build_state_df <- function(model, state_names, yrs) {
  d    <- as_tibble(t(model$states))
  d.se <- as_tibble(t(model$states.se))
  names(d) <- names(d.se) <- state_names
  d$Year <- yrs
  
  d_long <- d %>%
    pivot_longer(-Year, names_to = "Site_Class", values_to = "log_est") %>%
    separate(Site_Class, into = c("Site", "Class"), sep = "_", remove = FALSE)
  
  d_se_long <- d.se %>%
    pivot_longer(everything(), names_to = "Site_Class", values_to = "log_se")
  
  bind_cols(d_long, d_se_long %>% select(log_se)) %>%
    mutate(
      SiteType = SITE_TYPE[Site],
      lo89 = log_est + qnorm(CI_LO) * log_se,
      hi89 = log_est + qnorm(CI_HI) * log_se
    )
}

d_plot <- build_state_df(BESTMODEL, state_names_8site, years)

# ── PLOT 1: Log-abundance index — all 8 sites, faceted by class ───────────────
# (replicate of 6-site plot, now including DR and PB with dashed lines)
d_diff <- d_plot %>%
  group_by(Site_Class) %>%
  mutate(
    log_diff = log_est - first(log_est),
    lo89_d   = log_diff + qnorm(CI_LO) * log_se,
    hi89_d   = log_diff + qnorm(CI_HI) * log_se
  ) %>%
  ungroup()

p1 <- ggplot(d_diff, aes(x = Year, y = log_diff,
                         colour = Site, linetype = Site, group = Site)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_ribbon(aes(ymin = lo89_d, ymax = hi89_d, fill = Site),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = SITE_COLS) +
  scale_fill_manual(values = SITE_COLS) +
  scale_linetype_manual(values = SITE_LTY,
                        labels = c(paste0(names(SITE_LTY)[SITE_LTY==1]," (breeding)"),
                                   paste0(names(SITE_LTY)[SITE_LTY==2]," (haul-out)"))) +
  facet_wrap(~Class, ncol = 1) +
  labs(x = "Year", y = "Log-abundance index (relative to 1997)",
       title = "8-Site Population Trends: Breeding vs Haul-out Sites",
       subtitle = paste0("Bands = ", CI_LABEL,
                         "; dashed lines = DR and PB (haul-out only)"),
       caption = "DR and PB pup states unobserved (NA) — excluded from pup panel") +
  theme_seal() +
  theme(legend.position = "right")

ggsave("Output/Plots/8site_log_abundance.jpeg", p1, width = 22, height = 30, units = "cm")

# ── PLOT 2: Haul-out sites vs breeding site mean ──────────────────────────────
# KEY NEW PLOT: Do DR and PB track the breeding sites?
# Tests ecological hypothesis: same population, different use pattern

d_comparison <- d_diff %>%
  dplyr::filter(Class %in% c("Adult", "Molt")) %>%
  mutate(Group = case_when(
    Site %in% c("DR", "PB") ~ paste0(Site, " (haul-out)"),
    TRUE ~ "Breeding site mean"
  ))

breed_mean <- d_diff %>%
  dplyr::filter(Class %in% c("Adult", "Molt"),
                !Site %in% c("DR", "PB")) %>%
  group_by(Year, Class) %>%
  summarise(
    log_diff = mean(log_diff),
    lo89_d   = mean(lo89_d),
    hi89_d   = mean(hi89_d),
    .groups  = "drop"
  ) %>%
  mutate(Group = "Breeding site mean", Site = "Breeding")

haulout_sites <- d_diff %>%
  dplyr::filter(Site %in% c("DR","PB"), Class %in% c("Adult","Molt")) %>%
  mutate(Group = paste0(Site, " (haul-out)"))

compare_df <- bind_rows(
  breed_mean %>% select(Year, Class, log_diff, lo89_d, hi89_d, Group),
  haulout_sites %>% select(Year, Class, log_diff, lo89_d, hi89_d, Group)
)

p2 <- ggplot(compare_df, aes(x = Year, y = log_diff,
                             colour = Group, fill = Group, linetype = Group)) +
  geom_ribbon(aes(ymin = lo89_d, ymax = hi89_d), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  scale_colour_manual(values = c("Breeding site mean" = "grey30",
                                 "DR (haul-out)" = SITE_COLS["DR"],
                                 "PB (haul-out)" = SITE_COLS["PB"])) +
  scale_fill_manual(  values = c("Breeding site mean" = "grey50",
                                 "DR (haul-out)" = SITE_COLS["DR"],
                                 "PB (haul-out)" = SITE_COLS["PB"])) +
  scale_linetype_manual(values = c("Breeding site mean" = 1,
                                   "DR (haul-out)" = 2,
                                   "PB (haul-out)" = 2)) +
  facet_wrap(~Class, ncol = 1) +
  labs(x = "Year", y = "Log-abundance index relative to 1997",
       title = "Haul-out Sites vs Breeding Site Mean Trend",
       subtitle = paste0("Do DR and PB track the breeding population? Bands = ", CI_LABEL),
       colour = NULL, fill = NULL, linetype = NULL) +
  theme_seal()

ggsave("Output/Plots/8site_haulout_vs_breeding.jpeg", p2,
       width = 22, height = 24, units = "cm")

# ── PLOT 3: Between-site correlation (adult class) — all 8 sites ─────────────
# Extends IPM portfolio correlation matrix to include DR and PB
adult_states <- d_plot %>%
  dplyr::filter(Class == "Adult") %>%
  select(Year, Site, log_est) %>%
  pivot_wider(names_from = Site, values_from = log_est) %>%
  select(-Year)

cor_mat_8site <- cor(adult_states, use = "pairwise.complete.obs")

cor_long_8site <- as.data.frame(cor_mat_8site) %>%
  tibble::rownames_to_column("Site1") %>%
  pivot_longer(-Site1, names_to = "Site2", values_to = "r") %>%
  mutate(
    Type1 = SITE_TYPE[Site1],
    Type2 = SITE_TYPE[Site2],
    PairType = case_when(
      Type1 == "Breeding"      & Type2 == "Breeding"      ~ "Breeding-Breeding",
      Type1 == "Haul-out only" & Type2 == "Haul-out only" ~ "Haul-out-Haul-out",
      TRUE                                                  ~ "Breeding-Haul-out"
    )
  )

# Order sites: breeding first, haul-out last
site_order <- c("BL","DE","DP","PRH","TB","TP","DR","PB")
cor_long_8site <- cor_long_8site %>%
  mutate(Site1 = factor(Site1, levels = site_order),
         Site2 = factor(Site2, levels = site_order))

p3 <- ggplot(cor_long_8site, aes(x = Site1, y = fct_rev(Site2), fill = r)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", r), colour = abs(r) > 0.5), size = 3.5) +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "grey30"), guide = "none") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  geom_hline(yintercept = 2.5, linewidth = 1.2, colour = "grey40") +
  geom_vline(xintercept = 6.5, linewidth = 1.2, colour = "grey40") +
  annotate("text", x = 3.5, y = 8.7, label = "Breeding sites", size = 3.5) +
  annotate("text", x = 7.5, y = 8.7, label = "Haul-out", size = 3.5) +
  labs(x = NULL, y = NULL,
       title = "Between-site Synchrony: Adult λ Correlation (All 8 Sites)",
       subtitle = "Grey line separates breeding sites (left/bottom) from haul-out-only sites (right/top)",
       caption = "High correlation between DR/PB and breeding sites suggests metapopulation integration") +
  coord_fixed() +
  theme_seal() +
  theme(axis.text = element_text(size = 11),
        panel.grid = element_blank(),
        legend.position = "right")

ggsave("Output/Plots/8site_correlation_matrix.jpeg", p3,
       width = 22, height = 20, units = "cm")

# ── PLOT 4: Proportion of total estimated abundance at haul-out sites ─────────
# Addresses: what fraction of PRNS harbor seals are at non-breeding sites?

d_abund <- d_plot %>%
  mutate(est_count = exp(log_est)) %>%
  dplyr::filter(Class != "Pup") %>%   # Pup NA for DR/PB; focus on adult+molt
  group_by(Year) %>%
  mutate(total_all = sum(est_count, na.rm = TRUE)) %>%
  group_by(Year, Site) %>%
  summarise(site_total = sum(est_count, na.rm = TRUE),
            total_all  = first(total_all), .groups = "drop") %>%
  mutate(
    prop_haulout = site_total / total_all,
    SiteType = SITE_TYPE[Site]
  )

haulout_prop <- d_abund %>%
  dplyr::filter(SiteType == "Haul-out only") %>%
  group_by(Year) %>%
  summarise(haulout_total = sum(site_total),
            total_all     = first(total_all), .groups = "drop") %>%
  mutate(prop = haulout_total / total_all)

p4 <- ggplot(haulout_prop, aes(x = Year, y = prop * 100)) +
  geom_area(fill = "#984EA3", alpha = 0.4) +
  geom_line(colour = "#984EA3", linewidth = 1.2) +
  scale_y_continuous(limits = c(0, NA),
                     labels = function(x) paste0(x, "%")) +
  labs(x = "Year",
       y = "% of total adult + molt abundance",
       title = "Proportion of Estimated Abundance at Haul-out-Only Sites (DR + PB)",
       subtitle = "Adult and molting counts only (pup counts not available at DR/PB)",
       caption = "Derived from MARSS latent state estimates; does not account for IPM demographic structure") +
  theme_seal()

ggsave("Output/Plots/8site_haulout_proportion.jpeg", p4,
       width = 22, height = 14, units = "cm")

# ── PLOT 5: Covariate effects — 8-site forest plot ────────────────────────────
# Aligned with IPM forest plot style (89% CI, grouped by covariate type)

coef_df <- tidy(CIs_8site) %>%
  dplyr::filter(str_detect(term, "^C\\."))

coef_df <- coef_df %>%
  mutate(
    Class = case_when(
      str_detect(term, "_A$|_A,") ~ "Adult",
      str_detect(term, "_M$|_M,") ~ "Molting",
      str_detect(term, "_P$|_P,") ~ "Pup",
      TRUE ~ "All"
    ),
    Significant = case_when(
      conf.up < 0 ~ "Negative",
      conf.low > 0 ~ "Positive",
      TRUE ~ "Neutral"
    ),
    Group = case_when(
      str_detect(term, "MOCI")     ~ "MOCI",
      str_detect(term, "Dist")     ~ "Disturbance",
      str_detect(term, "Coy")      ~ "Coyote",
      str_detect(term, "ES|eSeal") ~ "Elephant seal",
      TRUE ~ "Other"
    ),
    Label = term %>%
      str_remove("^C\\.") %>%
      str_replace_all("_", " ") %>%
      str_replace_all("AMJ", "Spring AMJ") %>%
      str_replace_all("JFM", "Winter JFM") %>%
      str_replace_all("OND", "Fall OND") %>%
      str_replace_all("\\bDist\\b", "Disturbance") %>%
      str_replace_all("\\bCoy\\b", "Coyote") %>%
      str_to_title() %>%
      str_replace_all("\\bBl\\b","BL") %>% str_replace_all("\\bDe\\b","DE") %>%
      str_replace_all("\\bDp\\b","DP") %>% str_replace_all("\\bDr\\b","DR") %>%
      str_replace_all("\\bPb\\b","PB") %>% str_replace_all("\\bPrh\\b","PRH") %>%
      str_replace_all("\\bTb\\b","TB") %>% str_replace_all("\\bTp\\b","TP") %>%
      str_replace_all("\\bMoci\\b","MOCI") %>%
      str_replace_all("\\bEs-De\\b","Elephant Seal DE") %>%
      str_replace_all("\\bEs-Prh\\b","Elephant Seal PRH"),
    Site = str_extract(term, "BL|DE|DP|DR|PB|PRH|TB|TP"),
    SiteType = case_when(
      Site %in% c("DR","PB") ~ "Haul-out",
      !is.na(Site)           ~ "Breeding",
      TRUE                   ~ "Shared"
    )
  ) %>%
  arrange(Group, Class, Label) %>%
  mutate(Label = fct_inorder(Label))

grp_cols <- c(MOCI = "#2166AC", Disturbance = "#8C510A",
              Coyote = "#B2182B", "Elephant seal" = "#762A83")

p5 <- ggplot(coef_df, aes(y = Label, x = estimate)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_linerange(aes(xmin = conf.low, xmax = conf.up,
                     colour = Group, alpha = SiteType), linewidth = 0.8) +
  geom_point(aes(colour = Group, shape = Significant,
                 size = SiteType == "Haul-out")) +
  scale_colour_manual(values = grp_cols) +
  scale_shape_manual(values = c("Negative" = 25, "Neutral" = 16, "Positive" = 24),
                     name = NULL) +
  scale_alpha_manual(values = c("Breeding" = 0.8, "Haul-out" = 1.0, "Shared" = 0.7),
                     name = "Site type") +
  scale_size_manual(values = c("TRUE" = 3.5, "FALSE" = 2.5), guide = "none") +
  facet_grid(Group ~ Class, scales = "free_y", space = "free_y") +
  labs(x = "Coefficient estimate", y = NULL,
       title = "MARSS Covariate Effects — 8-Site Analysis",
       subtitle = paste0(CI_LABEL,
                         "; triangles = significant; larger points = haul-out sites (DR/PB)"),
       colour = "Covariate group") +
  theme_seal(base_size = 12) +
  theme(strip.text.y = element_text(angle = 0, face = "bold"),
        panel.grid.major.y = element_blank())

ggsave("Output/Plots/8site_covariate_forest.jpeg", p5,
       width = 32, height = 36, units = "cm", dpi = 200)

# ── PLOT 6: MOCI response comparison — haul-out vs breeding sites ─────────────
# Do DR and PB respond to MOCI similarly to breeding sites?
# This tests whether ocean forcing acts uniformly across the population.

# Extract MOCI covariate effects by class
moci_coefs <- coef_df %>%
  dplyr::filter(Group == "MOCI") %>%
  mutate(
    Season = case_when(
      str_detect(term, "JFM") ~ "Winter JFM",
      str_detect(term, "AMJ") ~ "Spring AMJ",
      str_detect(term, "OND") ~ "Fall OND"
    )
  )

p6 <- ggplot(moci_coefs, aes(x = Season, y = estimate,
                             colour = SiteType, shape = Significant)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.up),
                  position = position_dodge(width = 0.4), size = 0.7) +
  scale_colour_manual(values = c("Breeding" = "#2166AC",
                                 "Haul-out" = "#984EA3",
                                 "Shared"   = "grey50")) +
  scale_shape_manual(values = c("Negative" = 25, "Neutral" = 16, "Positive" = 24),
                     name = NULL) +
  facet_wrap(~Class, ncol = 1) +
  labs(x = "MOCI Season", y = "Coefficient estimate",
       title = "MOCI Effects: Breeding Sites vs Haul-out Sites",
       subtitle = paste0("Similar MOCI responses between DR/PB and breeding sites\n",
                         "suggest shared prey field; divergence suggests site-specific foraging. ",
                         CI_LABEL),
       colour = "Site type") +
  theme_seal()

ggsave("Output/Plots/8site_MOCI_comparison.jpeg", p6,
       width = 22, height = 28, units = "cm")

# ── PLOT 7: Total population — breeding vs haul-out separate ──────────────────
# Analogous to IPM total population plot but disaggregated by site type

d_abund_class <- d_plot %>%
  mutate(est_count = exp(log_est),
         lo_count  = exp(lo89),
         hi_count  = exp(hi89)) %>%
  group_by(Year, SiteType, Class) %>%
  summarise(
    total    = sum(est_count, na.rm = TRUE),
    total_lo = sum(lo_count,  na.rm = TRUE),
    total_hi = sum(hi_count,  na.rm = TRUE),
    .groups  = "drop"
  )

p7 <- ggplot(d_abund_class %>% dplyr::filter(Class != "Pup" | SiteType == "Breeding"),
             aes(x = Year, y = total, colour = SiteType, fill = SiteType)) +
  geom_ribbon(aes(ymin = total_lo, ymax = total_hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.1) +
  scale_y_continuous(limits = c(0, NA), labels = scales::comma,
                     expand = c(0, 0)) +
  scale_colour_manual(values = c("Breeding" = "#2166AC", "Haul-out only" = "#984EA3")) +
  scale_fill_manual(values   = c("Breeding" = "#AECDE8", "Haul-out only" = "#D9B3E8")) +
  facet_wrap(~Class, ncol = 1, scales = "free_y") +
  labs(x = "Year", y = "Estimated abundance",
       title = "Total Estimated Abundance by Site Type and Class",
       subtitle = paste0("Breeding sites (BL/DE/DP/PRH/TB/TP) vs Haul-out only (DR/PB). ",
                         CI_LABEL),
       colour = NULL, fill = NULL,
       caption = "DR and PB pup counts not available; pup panel shows breeding sites only") +
  theme_seal()

ggsave("Output/Plots/8site_total_by_type.jpeg", p7, width = 22, height = 30, units = "cm")

# ── PLOT 8: % change by site (pct change plot from 6-site, extended) ──────────
d_pct <- d_diff %>%
  mutate(
    pct_change = (exp(log_diff) - 1) * 100,
    pct_lo     = (exp(lo89_d) - 1) * 100,
    pct_hi     = (exp(hi89_d) - 1) * 100
  )

p8 <- ggplot(d_pct, aes(x = Year, y = pct_change, colour = Site, group = Site,
                        linetype = Site)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_ribbon(aes(ymin = pct_lo, ymax = pct_hi, fill = Site),
              alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = SITE_COLS) +
  scale_fill_manual(  values = SITE_COLS) +
  scale_linetype_manual(values = SITE_LTY) +
  facet_wrap(~Class, ncol = 1, scales = "free_y") +
  labs(x = "Year", y = "% change from 1997",
       title = "Percentage Change in Estimated Abundance from 1997 — All 8 Sites",
       subtitle = paste0("Dashed lines = DR and PB (haul-out only). Bands = ", CI_LABEL)) +
  theme_seal() +
  theme(legend.position = "right")

ggsave("Output/Plots/8site_pct_change.jpeg", p8, width = 22, height = 28, units = "cm")

# ── Summary: key statistics ───────────────────────────────────────────────────
cat("\n── Summary: Haul-out site statistics ────────────────────────────────────\n")
cat("Mean correlation of DR/PB with breeding sites (adult class):\n")
dr_pb_corr <- cor_mat_8site[c("DR","PB"), c("BL","DE","DP","PRH","TB","TP")]
cat("  DR with breeding sites:", round(mean(dr_pb_corr["DR",]), 3), "\n")
cat("  PB with breeding sites:", round(mean(dr_pb_corr["PB",]), 3), "\n")
cat("  Mean correlation DR/PB vs breeding:", round(mean(dr_pb_corr), 3), "\n\n")

cat("Mean % haul-out site abundance as fraction of total:\n")
cat(" ", round(mean(haulout_prop$prop) * 100, 1), "% of adult + molt counts\n")
cat("Range:", round(min(haulout_prop$prop)*100, 1), "-",
    round(max(haulout_prop$prop)*100, 1), "%\n")

cat("\n── All plots saved to Output/Plots/8site_* ──────────────────────────────\n")
