# ============================================================================
# 04_marss6_plots.R
# ----------------------------------------------------------------------------
# 6-SITE MARSS — COVARIATE EFFECT PLOTS
# Renders the C-matrix covariate effects from the best framework model (Model A)
# as % change in abundance per +1 SD covariate, with 89% CIs.
#
# Prereqs (either):
#   source("Code/03_marss6_models.R")        # creates BESTMODEL, CIs in memory
#   OR  load("Output/marss6_best.RData")     # auto-loaded below if missing
#
# Saves: Output/Plots/marss6_effects_*.jpeg
#        Output/marss6_significant_effects_summary.csv
# ============================================================================

library(MARSS)
library(tidyverse)
library(patchwork)
library(broom)   # tidy() for marssMLE / marssParamCIs

dir.create("Output/Plots", showWarnings = FALSE)

# ── Load best model + CIs if not already in memory ──────────────────────────
if (!exists("BESTMODEL") || !exists("CIs")) {
  if (file.exists("Output/marss6_best.RData")) {
    load("Output/marss6_best.RData")   # BESTMODEL, CIs, df_aic_6site
    cat("Loaded Output/marss6_best.RData\n")
  } else {
    stop("BESTMODEL/CIs not found. Run 03_marss6_models.R first.")
  }
}

# ── Constants matching IPM v3.2 style ───────────────────────────────────────
CI_LABEL <- "89% CI"
EFFECT_YLIM <- c(-30, 50)   # % change axis range for all effect panels

theme_seal <- function(base_size = 14) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.grid.major = element_line(colour = "grey88", linewidth = 0.4),
      panel.border     = element_rect(colour = "grey70", fill = NA, linewidth = 0.5),
      axis.title       = element_text(size = rel(0.95)),
      legend.position  = "bottom",
      strip.text       = element_text(size = rel(0.90), face = "bold"),
      strip.background = element_rect(fill = "grey94", colour = "grey80"),
      plot.title       = element_text(size = rel(1.05), face = "bold"),
      plot.subtitle    = element_text(size = rel(0.88), colour = "grey40"),
      plot.margin      = margin(10, 14, 10, 10)
    )
}

# ── Tidy C-matrix coefficients once (guarded) ───────────────────────────────
C_df <- tryCatch(
  tidy(CIs) %>%
    tibble::as_tibble() %>%                       # tidy.marssMLE returns a data.frame
    dplyr::filter(str_detect(term, "^C\\.")) %>%
    dplyr::mutate(coef_name = str_remove(term, "^C\\.")),
  error = function(e)
    stop("Could not tidy() the CIs object: ", conditionMessage(e),
         "\n  Ensure `broom` is installed and CIs came from MARSSparamCIs().")
)

cat("\n── C-matrix coefficients (89% CI) ──────────────────────────────────\n")
print(C_df %>% select(coef_name, estimate, conf.low, conf.up), n = nrow(C_df))

# ── Guarded effect-plot builder ─────────────────────────────────────────────
# Returns NULL (with a warning) if the coefficient isn't in the model, so a
# single renamed/absent term can't abort the whole script.
make_effect_plot <- function(coef_name, x_label, title,
                             x_range = seq(-2, 2, length.out = 100),
                             ylim = EFFECT_YLIM) {
  row <- C_df %>% dplyr::filter(coef_name == !!coef_name)
  if (nrow(row) == 0) {
    warning(sprintf("Coefficient '%s' not found — skipping panel.", coef_name))
    return(NULL)
  }
  est <- row$estimate[1]; lo <- row$conf.low[1]; hi <- row$conf.up[1]
  
  # MARSS C effects are additive on log-abundance -> multiplicative % change
  df <- tibble(
    x      = x_range,
    effect = (exp(est * x_range) - 1) * 100,
    lo     = (exp(lo  * x_range) - 1) * 100,
    hi     = (exp(hi  * x_range) - 1) * 100
  )
  clr  <- if (est < 0) "red3" else "blue3"
  fill <- if (est < 0) "red"  else "blue"
  
  ggplot(df, aes(x = x, y = effect)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "gray50") +
    geom_vline(xintercept = 0, linetype = 2, colour = "gray50") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = fill) +
    geom_line(linewidth = 1.2, colour = clr) +
    labs(x = x_label, y = "% change in abundance", title = title) +
    coord_cartesian(ylim = ylim) +
    theme_seal(base_size = 13) +
    theme(plot.title = element_text(size = rel(0.95)))
}

# Helper: wrap only the panels that exist (drops NULLs)
wrap_existing <- function(plots, ncol, title, subtitle = NULL) {
  plots <- Filter(Negate(is.null), plots)
  if (length(plots) == 0) return(NULL)
  wrap_plots(plots, ncol = ncol) +
    plot_annotation(title = title, subtitle = subtitle,
                    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                                  plot.subtitle = element_text(size = 11)))
}

# ── Coyote effects (BL, DE, DP) ─────────────────────────────────────────────
coy <- list(
  make_effect_plot("Coy_BL", "Coyote rate (SD)", "Coyote - BL"),
  make_effect_plot("Coy_DE", "Coyote rate (SD)", "Coyote - DE"),
  make_effect_plot("Coy_DP", "Coyote rate (SD)", "Coyote - DP")
)
p_coy <- wrap_existing(coy, 3, "Coyote Effects on Harbor Seal Abundance")
if (!is.null(p_coy))
  ggsave("Output/Plots/marss6_effects_coyote.jpeg", p_coy,
         width = 30, height = 10, units = "cm")

# ── Disturbance effects (BL, DE, TB shown; add others as needed) ────────────
dst <- list(
  make_effect_plot("Dist_BL", "Disturbance rate (SD)", "Disturbance - BL"),
  make_effect_plot("Dist_DE", "Disturbance rate (SD)", "Disturbance - DE"),
  make_effect_plot("Dist_TB", "Disturbance rate (SD)", "Disturbance - TB")
)
p_dst <- wrap_existing(dst, 3, "Disturbance Effects on Harbor Seal Abundance")
if (!is.null(p_dst))
  ggsave("Output/Plots/marss6_effects_disturbance.jpeg", p_dst,
         width = 30, height = 10, units = "cm")

# ── MOCI effects (pup + molt classes) ───────────────────────────────────────
mci <- list(
  make_effect_plot("MOCI_OND_P", "MOCI prior fall (SD)",  "MOCI Fall (Pup)"),
  make_effect_plot("MOCI_AMJ_P", "MOCI spring (SD)",      "MOCI Spring (Pup)"),
  make_effect_plot("MOCI_AMJ_M", "MOCI spring (SD)",      "MOCI Spring (Molt)"),
  make_effect_plot("MOCI_JFM_M", "MOCI winter (SD)",      "MOCI Winter (Molt)")
)
p_mci <- wrap_existing(mci, 2, "MOCI (Ocean Climate) Effects on Abundance")
if (!is.null(p_mci))
  ggsave("Output/Plots/marss6_effects_moci.jpeg", p_mci,
         width = 24, height = 20, units = "cm")

# ── Elephant seal effects (DE, PRH) ─────────────────────────────────────────
es <- list(
  make_effect_plot("ES-DE",  "Elephant seal abundance (SD)", "Elephant Seal - DE"),
  make_effect_plot("ES-PRH", "Elephant seal abundance (SD)", "Elephant Seal - PRH")
)
p_es <- wrap_existing(es, 2, "Elephant Seal Effects on Harbor Seal Abundance")
if (!is.null(p_es))
  ggsave("Output/Plots/marss6_effects_eseal.jpeg", p_es,
         width = 20, height = 10, units = "cm")

# ── Combined panel of the key effects ───────────────────────────────────────
all_panels <- list(
  make_effect_plot("Coy_BL", "Coyote (SD)",       "Coyote - BL"),
  make_effect_plot("Coy_DE", "Coyote (SD)",       "Coyote - DE"),
  make_effect_plot("Coy_DP", "Coyote (SD)",       "Coyote - DP"),
  make_effect_plot("Dist_BL","Disturbance (SD)",  "Disturbance - BL"),
  make_effect_plot("Dist_DE","Disturbance (SD)",  "Disturbance - DE"),
  make_effect_plot("Dist_TB","Disturbance (SD)",  "Disturbance - TB"),
  make_effect_plot("MOCI_OND_P","MOCI Fall (SD)", "MOCI Fall (Pup)"),
  make_effect_plot("MOCI_AMJ_P","MOCI Spring (SD)","MOCI Spring (Pup)"),
  make_effect_plot("ES-PRH",  "Elephant Seal (SD)","Elephant Seal - PRH"),
  make_effect_plot("MOCI_AMJ_M","MOCI Spring (SD)","MOCI Spring (Molt)"),
  make_effect_plot("MOCI_JFM_M","MOCI Winter (SD)","MOCI Winter (Molt)")
)
p_all <- wrap_existing(all_panels, 3,
                       "Covariate Effects on Harbor Seal Abundance (6-site MARSS)",
                       "Effect = % change in abundance per +1 SD covariate; bands = 89% CI")
if (!is.null(p_all))
  ggsave("Output/Plots/marss6_effects_all.jpeg", p_all,
         width = 36, height = 40, units = "cm")

# ── Summary table of effects (% change per SD) ──────────────────────────────
wanted <- c("Coy_BL","Coy_DE","Coy_DP","Dist_BL","Dist_DE","Dist_TB",
            "MOCI_OND_P","MOCI_AMJ_P","MOCI_AMJ_M","MOCI_JFM_M","ES-DE","ES-PRH")
effects_summary <- C_df %>%
  dplyr::filter(coef_name %in% wanted) %>%
  mutate(
    pct_change_per_SD = (exp(estimate) - 1) * 100,
    pct_change_lo     = (exp(conf.low) - 1) * 100,
    pct_change_hi     = (exp(conf.up)  - 1) * 100,
    significant       = (conf.low > 0) | (conf.up < 0)
  ) %>%
  select(coef_name, estimate, conf.low, conf.up,
         pct_change_per_SD, pct_change_lo, pct_change_hi, significant) %>%
  arrange(pct_change_per_SD)

cat("\n── Significant-effect summary (% change per SD) ────────────────────\n")
print(effects_summary, n = nrow(effects_summary))
write_csv(effects_summary, "Output/marss6_significant_effects_summary.csv")

cat("\nPlots -> Output/Plots/marss6_effects_*.jpeg\n")
cat("Table -> Output/marss6_significant_effects_summary.csv\n")
cat("Next: source(\"Code/05_ipm_model.R\")\n")