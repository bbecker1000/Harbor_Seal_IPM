#Effects plots

# ============================================================================
# COVARIATE EFFECTS PLOTS
# ============================================================================

library(tidyverse)
library(patchwork)

# Get coefficient estimates from the model
coefs <- coef(BESTMODEL)

# Extract C matrix coefficients
C_coefs <- coefs[grep("^C\\.", names(coefs))]
names(C_coefs)

# ----------------------------------------------------------------------------
# FUNCTION TO CREATE EFFECT PLOT
# ----------------------------------------------------------------------------

create_effect_plot <- function(coef_name, coef_value, 
                               covariate_name, 
                               x_label, 
                               title,
                               class_type = "All") {
  
  # Create range of covariate values (scaled: -2 to +2 SD)
  x_range <- seq(-2, 2, length.out = 100)
  
  # Calculate effect on log scale
  log_effect <- coef_value * x_range
  
  # Convert to multiplicative effect (percent change from mean)
  mult_effect <- (exp(log_effect) - 1) * 100
  
  # Create data frame
  plot_df <- tibble(
    x = x_range,
    effect = mult_effect
  )
  
  # Get CI from CIs object
  ci_row <- tidy(CIs) %>% filter(term == paste0("C.", coef_name))
  
  if (nrow(ci_row) > 0) {
    lo_effect <- (exp(ci_row$conf.low * x_range) - 1) * 100
    hi_effect <- (exp(ci_row$conf.up * x_range) - 1) * 100
    plot_df$lo <- lo_effect
    plot_df$hi <- hi_effect
  }
  
  # Determine color based on coefficient sign
  line_color <- ifelse(coef_value < 0, "red3", "blue3")
  fill_color <- ifelse(coef_value < 0, "red", "blue")
  
  p <- ggplot(plot_df, aes(x = x, y = effect)) +
    geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
    geom_vline(xintercept = 0, linetype = 2, color = "gray50")
  
  if ("lo" %in% names(plot_df)) {
    p <- p + geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = fill_color)
  }
  
  p <- p +
    geom_line(linewidth = 1.2, color = line_color) +
    labs(
      x = x_label,
      y = "% Change in Abundance",
      title = title
    ) +
    ylim(-30, 50) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 10)
    )
  
  return(p)
}

# ----------------------------------------------------------------------------
# EXTRACT SIGNIFICANT COEFFICIENTS
# ----------------------------------------------------------------------------

# Get coefficient values from model
C_df <- tidy(CIs) %>% 
  filter(str_detect(term, "^C\\.")) %>%
  mutate(coef_name = str_remove(term, "^C\\."))

# View all C coefficients
print(C_df %>% select(term, estimate, conf.low, conf.up))

# ----------------------------------------------------------------------------
# CREATE INDIVIDUAL EFFECT PLOTS
# ----------------------------------------------------------------------------

# COYOTE EFFECTS
p_coy_bl <- create_effect_plot(
  coef_name = "Coy_BL",
  coef_value = C_df %>% filter(coef_name == "Coy_BL") %>% pull(estimate),
  covariate_name = "Coyote_BL",
  x_label = "Coyote Sighting Rate (SD)",
  title = "Coyote - BL (All Classes)"
)

p_coy_de <- create_effect_plot(
  coef_name = "Coy_DE",
  coef_value = C_df %>% filter(coef_name == "Coy_DE") %>% pull(estimate),
  covariate_name = "Coyote_DE",
  x_label = "Coyote Sighting Rate (SD)",
  title = "Coyote - DE (All Classes)"
)

p_coy_dp <- create_effect_plot(
  coef_name = "Coy_DP",
  coef_value = C_df %>% filter(coef_name == "Coy_DP") %>% pull(estimate),
  covariate_name = "Coyote_DP",
  x_label = "Coyote Sighting Rate (SD)",
  title = "Coyote - DP (All Classes)"
)

# DISTURBANCE EFFECTS
p_dist_bl <- create_effect_plot(
  coef_name = "Dist_BL",
  coef_value = C_df %>% filter(coef_name == "Dist_BL") %>% pull(estimate),
  covariate_name = "Dist_BL",
  x_label = "Disturbance Rate (SD)",
  title = "Disturbance - BL (All Classes)"
)

p_dist_de <- create_effect_plot(
  coef_name = "Dist_DE",
  coef_value = C_df %>% filter(coef_name == "Dist_DE") %>% pull(estimate),
  covariate_name = "Dist_DE",
  x_label = "Disturbance Rate (SD)",
  title = "Disturbance - DE (All Classes)"
)

p_dist_tb <- create_effect_plot(
  coef_name = "Dist_TB",
  coef_value = C_df %>% filter(coef_name == "Dist_TB") %>% pull(estimate),
  covariate_name = "Dist_TB",
  x_label = "Disturbance Rate (SD)",
  title = "Disturbance - TB (All Classes)"
)

# MOCI EFFECTS - PUPS
p_moci_ond_p <- create_effect_plot(
  coef_name = "MOCI_OND_P",
  coef_value = C_df %>% filter(coef_name == "MOCI_OND_P") %>% pull(estimate),
  covariate_name = "MOCI_OND",
  x_label = "MOCI Prior Fall (SD)",
  title = "MOCI Prior Fall (Pup)"
)

p_moci_amj_p <- create_effect_plot(
  coef_name = "MOCI_AMJ_P",
  coef_value = C_df %>% filter(coef_name == "MOCI_AMJ_P") %>% pull(estimate),
  covariate_name = "MOCI_AMJ",
  x_label = "MOCI Spring (SD)",
  title = "MOCI Spring (Pup)"
)

# MOCI EFFECTS - MOLTING
p_moci_amj_m <- create_effect_plot(
  coef_name = "MOCI_AMJ_M",
  coef_value = C_df %>% filter(coef_name == "MOCI_AMJ_M") %>% pull(estimate),
  covariate_name = "MOCI_AMJ",
  x_label = "MOCI Spring (SD)",
  title = "MOCI Spring (Molting)"
)

p_moci_jfm_m <- create_effect_plot(
  coef_name = "MOCI_JFM_M",
  coef_value = C_df %>% filter(coef_name == "MOCI_JFM_M") %>% pull(estimate),
  covariate_name = "MOCI_JFM",
  x_label = "MOCI Winter (SD)",
  title = "MOCI Winter (Molting)"
)

# ELEPHANT SEALS
p_es_prh <- create_effect_plot(
  coef_name = "ES-PRH",
  coef_value = C_df %>% filter(coef_name == "ES-PRH") %>% pull(estimate),
  covariate_name = "eSeal",
  x_label = "Elephant Seal Abundance (SD)",
  title = "Elephant Seals - PRH (All Classes)"
)

# ----------------------------------------------------------------------------
# COMBINE PLOTS
# ----------------------------------------------------------------------------

# Coyote effects
coyote_plots <- p_coy_bl + p_coy_de + p_coy_dp +
  plot_layout(ncol = 3) +
  plot_annotation(title = "Coyote Effects on Harbor Seal Abundance",
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))

ggsave("Output/Plots/Effects_Coyote.jpeg", coyote_plots, width = 30, height = 10, units = "cm")

# Disturbance effects
disturbance_plots <- p_dist_bl + p_dist_de + p_dist_tb +
  plot_layout(ncol = 3) +
  plot_annotation(title = "Disturbance Effects on Harbor Seal Abundance",
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))

ggsave("Output/Plots/Effects_Disturbance.jpeg", disturbance_plots, width = 30, height = 10, units = "cm")

# MOCI effects
moci_plots <- (p_moci_ond_p + p_moci_amj_p) / (p_moci_amj_m + p_moci_jfm_m) +
  plot_annotation(title = "MOCI (Ocean Conditions) Effects on Harbor Seal Abundance",
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))

ggsave("Output/Plots/Effects_MOCI.jpeg", moci_plots, width = 24, height = 20, units = "cm")

# Elephant seal effects
ggsave("Output/Plots/Effects_ElephantSeals.jpeg", p_es_prh, width = 15, height = 12, units = "cm")

# ----------------------------------------------------------------------------
# ALL SIGNIFICANT EFFECTS IN ONE PLOT
# ----------------------------------------------------------------------------

all_effects <- wrap_plots(
  p_coy_bl, p_coy_de, p_coy_dp,
  p_dist_bl, p_dist_de, p_dist_tb,
  p_moci_ond_p, p_moci_amj_p, p_es_prh,
  p_moci_amj_m, p_moci_jfm_m, plot_spacer(),
  ncol = 3
) +
  plot_annotation(
    title = "Significant Covariate Effects on Harbor Seal Abundance",
    subtitle = "Effect shown as % change in abundance per SD change in covariate",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

ggsave("Output/Plots/Effects_All_Significant.jpeg", all_effects,
       width = 36, height = 40, units = "cm")


# ----------------------------------------------------------------------------
# SUMMARY TABLE OF EFFECTS
# ----------------------------------------------------------------------------

effects_summary <- C_df %>%
  filter(coef_name %in% c("Coy_BL", "Coy_DE", "Coy_DP",
                          "Dist_BL", "Dist_DE", "Dist_TB",
                          "MOCI_OND_P", "MOCI_AMJ_P", 
                          "MOCI_AMJ_M", "MOCI_JFM_M",
                          "ES-PRH")) %>%
  mutate(
    pct_change_per_SD = (exp(estimate) - 1) * 100,
    pct_change_lo = (exp(conf.low) - 1) * 100,
    pct_change_hi = (exp(conf.up) - 1) * 100
  ) %>%
  select(coef_name, estimate, conf.low, conf.up, 
         pct_change_per_SD, pct_change_lo, pct_change_hi) %>%
  arrange(pct_change_per_SD)

print(effects_summary)

write.csv(effects_summary, "Output/significant_effects_summary.csv", row.names = FALSE)