# ============================================================================
# 09_covariates_8site.R
# ----------------------------------------------------------------------------
# COVARIATE PREP — 8 SITES incl. DR & PB (1997–2025). Companion to the 8-site
# MARSS (Appendix A). Does NOT feed the IPM.
#
# Produces:
#   cov_t_scaled_8site : 20 x 29 numeric matrix (covariates x years), z-scored
#
# CANONICAL 20-row order (asserted; the 24-state C matrix indexes these rows):
#    1  MOCI_JFM     2  MOCI_AMJ     3  MOCI_OND
#    4  Dist_BL  5 Dist_DE  6 Dist_DP  7 Dist_DR  8 Dist_PB  9 Dist_PRH 10 Dist_TB 11 Dist_TP
#   12  Coyote_BL 13 Coyote_DE 14 Coyote_DP 15 Coyote_DR(0) 16 Coyote_PB(0)
#   17  Coyote_PRH(0) 18 Coyote_TB(0) 19 Coyote_TP(0)
#   20  eSeal_Sum_Imm_MaxCount
# DR/PB/PRH/TB/TP coyote rows are structurally zero (no sightings); kept for
# completeness and left unscaled (zero variance).
#
# Prereq: source("Code/08_data_prep_8site.R")  (years_8site alignment check)
# Next:   source("Code/10_marss8_models.R")
# ============================================================================

library(tidyverse)
library(readxl)
library(lubridate)

dir.create("Output", showWarnings = FALSE)

SITES_8 <- c("BL", "DE", "DP", "DR", "PB", "PRH", "TB", "TP")
CANON_COV_ROWS <- c(
  "MOCI_JFM", "MOCI_AMJ", "MOCI_OND",
  paste0("Dist_",   SITES_8),
  paste0("Coyote_", SITES_8),
  "eSeal_Sum_Imm_MaxCount"
)

# ── Coyote (3-yr weighted lag, per site) ────────────────────────────────────
coyote <- read_excel("Data/CoyoteSightings_2025.xlsx")
coyote_rate <- coyote %>%
  mutate(rate = `Number of days with coyote sightings` /
           `Total number of monitoring surveys`) %>%
  arrange(Site, Year) %>%
  group_by(Site) %>%
  mutate(
    rate_t1 = lag(rate, 1), rate_t2 = lag(rate, 2),
    has_t1 = !is.na(rate_t1), has_t2 = !is.na(rate_t2),
    weighted_rate = case_when(
      has_t1 & has_t2  ~ 0.5*rate + 0.3*rate_t1 + 0.2*rate_t2,
      has_t1 & !has_t2 ~ 0.7*rate + 0.3*rate_t1,
      TRUE             ~ 1.0*rate),
    weighted_rate = case_when(
      is.nan(weighted_rate) & has_t1 & has_t2 ~ (rate_t1 + rate_t2)/2,
      is.nan(weighted_rate) & has_t1           ~ rate_t1,
      is.nan(weighted_rate)                    ~ NA_real_,
      TRUE                                     ~ weighted_rate)) %>%
  ungroup() %>%
  select(Year, Site, weighted_rate)

coyote_wide <- coyote_rate %>%
  pivot_wider(names_from = Site, values_from = weighted_rate)

# Pre-survey years (1996–1999): all sites zero. new_rows fixes column order;
# DR & PB stay zero throughout (confirmed no coyote activity).
new_rows <- tibble(Year = 1996:1999,
                   BL = 0, DE = 0, DP = 0, DR = 0, PB = 0,
                   PRH = 0, TB = 0, TP = 0)
coyote_wide <- bind_rows(new_rows, coyote_wide) %>%
  arrange(Year) %>%
  rename_with(~ paste0("Coyote_", .), -Year)

# ── Human disturbance (per site; DR & PB RETAINED for the 8-site analysis) ──
HumanDisturbance <- read_excel("Data/HumanDisturbanceRate_1996To2025.xlsx")
HumanDisturbance <- HumanDisturbance[, -5]
HumanDisturbance <- HumanDisturbance %>%
  mutate(DistRate = SumOfDisturbanceCount / NSurveys) %>%
  select(-SumOfDisturbanceCount, -NSurveys)
HumanDisturbance.wide <- HumanDisturbance %>%
  pivot_wider(names_from = SiteCode, values_from = DistRate) %>%
  dplyr::filter(Year > 1996) %>%
  rename_with(~ paste0("Dist_", .), -Year)

# 2020 DP closure (COVID) -> 0, by name.
if ("Dist_DP" %in% names(HumanDisturbance.wide))
  HumanDisturbance.wide$Dist_DP[HumanDisturbance.wide$Year == 2020] <- 0

# ── MOCI (OND lead-shifted) ─────────────────────────────────────────────────
MOCI <- read_csv("Data/CaliforniaMOCI.csv", show_col_types = FALSE)
MOCI.dat <- MOCI %>%
  mutate(mean_value = (`North California (38-42N)` +
                         `Central California (34.5-38N)`) / 2) %>%
  select(Year, Season, mean_value) %>%
  mutate(Year = ifelse(Season == "OND", Year + 1, Year)) %>%
  dplyr::filter(Season %in% c("JFM", "AMJ", "OND")) %>%
  pivot_wider(names_from = Season, values_from = mean_value) %>%
  arrange(Year) %>%
  select(Year, JFM, AMJ, OND) %>%
  rename_with(~ paste0("MOCI_", .), -Year)

# ── Elephant seal (annual summed max, all PR sites) ─────────────────────────
eSeal <- read_excel("Data/Eseal_1981-2025_BySubsite.xlsx")
eSeal_max_imm <- eSeal %>%
  mutate(Year = year(StartDate)) %>%
  group_by(Year, SubSiteName, MatureCode) %>%
  summarise(MaxCount = max(Count, na.rm = TRUE), .groups = "drop") %>%
  group_by(Year) %>%
  summarise(eSeal_Sum_Imm_MaxCount = sum(MaxCount), .groups = "drop")

# ── Combine and restrict to study years ─────────────────────────────────────
covariates <- MOCI.dat %>%
  left_join(HumanDisturbance.wide, by = "Year") %>%
  left_join(coyote_wide,           by = "Year") %>%
  left_join(eSeal_max_imm,         by = "Year") %>%
  dplyr::filter(Year > 1996, Year < 2026) %>%
  arrange(Year)
cov_years <- covariates$Year

# ── Assemble rows in EXPLICIT canonical order ───────────────────────────────
pull_row <- function(name) {
  if (!name %in% names(covariates))
    stop("Covariate column not found: ", name,
         "\n  Available: ", paste(names(covariates), collapse = ", "))
  as.numeric(covariates[[name]])
}
cov_t <- do.call(rbind, setNames(lapply(CANON_COV_ROWS, pull_row), CANON_COV_ROWS))
colnames(cov_t) <- cov_years

# ── Single scaling pass (zero-variance rows left unscaled) ──────────────────
cov_t[is.na(cov_t)] <- 0
cov_t_scaled_8site <- cov_t
rows_to_scale <- apply(cov_t, 1, sd) > 0
cov_t_scaled_8site[rows_to_scale, ] <- t(scale(t(cov_t[rows_to_scale, ])))

# ── Assertions: shape, order, zero-coyote sanity, year alignment ────────────
stopifnot(
  identical(rownames(cov_t_scaled_8site), CANON_COV_ROWS),
  nrow(cov_t_scaled_8site) == 20L
)
for (z in c("Coyote_DR","Coyote_PB","Coyote_PRH","Coyote_TB","Coyote_TP"))
  if (any(cov_t_scaled_8site[z, ] != 0))
    warning(z, " expected all-zero but has non-zero values.")
if (exists("years_8site") &&
    !identical(as.numeric(cov_years), as.numeric(years_8site)))
  stop("Covariate years do not match years_8site from 08_data_prep_8site.R")

# ── Diagnostics ─────────────────────────────────────────────────────────────
cat("\n── 8-site covariate matrix ─────────────────────────────────────────\n")
cat(sprintf("cov_t_scaled_8site: %d covariates x %d years (%d-%d)\n",
            nrow(cov_t_scaled_8site), ncol(cov_t_scaled_8site),
            min(cov_years), max(cov_years)))
for (i in seq_len(nrow(cov_t_scaled_8site)))
  cat(sprintf("  %2d  %s\n", i, rownames(cov_t_scaled_8site)[i]))
cat("\nRow SDs (1 = scaled; 0 = structurally-zero coyote rows):\n")
print(round(apply(cov_t_scaled_8site, 1, sd), 2))

cat("\nObjects created: cov_t_scaled_8site\n")
cat("Next: source(\"Code/10_marss8_models.R\")\n")