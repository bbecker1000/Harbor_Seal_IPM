# ============================================================================
# COVARIATE PREP FOR 8-SITE MARSS ANALYSIS (1997-2025)
#
# IMPORTANT: This file is a COMPANION to the 6-site IPM analysis.
# Do NOT replace 4b_covariates.R — that file feeds the 6-site IPM.
#
# Key differences from 6-site version:
#   (1) Disturbance data kept for DR and PB (Dist_DR, Dist_PB added)
#   (2) Coyote_DR and Coyote_PB are NOT removed — they remain as zero series.
#       DR and PB have zero coyote sightings over the time series; these columns
#       are present for structural consistency but held at 0 in the C matrix.
#   (3) No elephant seal effect at DR or PB.
#
# Produces: cov_t_scaled_8site (20 × 29 matrix)
# Covariate order:
#   1:  MOCI_JFM
#   2:  MOCI_AMJ
#   3:  MOCI_OND
#   4:  Dist_BL
#   5:  Dist_DE
#   6:  Dist_DP
#   7:  Dist_DR   ← new
#   8:  Dist_PB   ← new
#   9:  Dist_PRH
#   10: Dist_TB
#   11: Dist_TP
#   12: Coyote_BL
#   13: Coyote_DE
#   14: Coyote_DP
#   15: Coyote_DR  (all zeros — confirmed no coyote activity)
#   16: Coyote_PB  (all zeros — confirmed no coyote activity)
#   17: Coyote_PRH (all zeros)
#   18: Coyote_TB  (all zeros)
#   19: Coyote_TP  (all zeros)
#   20: eSeal_Sum_Imm_MaxCount
# ============================================================================

library(readxl)
library(tidyverse)
library(dplyr)
library(lubridate)

# ── Coyote data ───────────────────────────────────────────────────────────────
coyote <- read_excel("Data/CoyoteSightings_2025.xlsx")

coyote_rate <- coyote %>%
  mutate(rate = `Number of days with coyote sightings` / `Total number of monitoring surveys`) %>%
  arrange(Site, Year) %>%
  group_by(Site) %>%
  mutate(
    rate_t1 = lag(rate, 1),
    rate_t2 = lag(rate, 2),
    has_t1  = !is.na(rate_t1),
    has_t2  = !is.na(rate_t2)
  ) %>%
  mutate(
    weighted_rate = case_when(
      has_t1 & has_t2  ~ 0.5 * rate + 0.3 * rate_t1 + 0.2 * rate_t2,
      has_t1 & !has_t2 ~ 0.7 * rate + 0.3 * rate_t1,
      TRUE             ~ 1.0 * rate
    ),
    weighted_rate = case_when(
      is.nan(weighted_rate) & has_t1 & has_t2 ~ (rate_t1 + rate_t2) / 2,
      is.nan(weighted_rate) & has_t1           ~ rate_t1,
      is.nan(weighted_rate)                    ~ NA_real_,
      TRUE                                     ~ weighted_rate
    )
  ) %>%
  select(Year, Site, weighted_rate) %>%
  ungroup() %>%
  arrange(Year, Site)

coyote_rate$weighted_rate[is.nan(coyote_rate$weighted_rate)] <- NA

coyote_wide <- coyote_rate %>%
  pivot_wider(names_from = Site, values_from = weighted_rate)

# Pre-1997 rows — DR and PB confirmed zero coyotes
new_rows <- tibble(
  Year = 1996:1999,
  BL = 0, DE = 0, DP = 0, DR = 0, PB = 0, PRH = 0, TB = 0, TP = 0
)

coyote_wide <- bind_rows(new_rows, coyote_wide) %>%
  arrange(Year) %>%
  rename_with(~ paste0("Coyote_", .), -Year)

# ── Human disturbance data ────────────────────────────────────────────────────
HumanDisturbance <- read_excel("Data/HumanDisturbanceRate_1996To2025.xlsx")
HumanDisturbance <- HumanDisturbance[, -5]
HumanDisturbance$DistRate <- HumanDisturbance$SumOfDisturbanceCount / HumanDisturbance$NSurveys
HumanDisturbance <- HumanDisturbance[, -c(3:4)]

HumanDisturbance.wide <- HumanDisturbance %>%
  pivot_wider(names_from = SiteCode, values_from = DistRate)

# KEY DIFFERENCE: Do NOT remove DR and PB — keep disturbance data for 8 sites
HumanDisturbance.wide <- HumanDisturbance.wide %>%
  dplyr::filter(Year > 1996)
# Note: NOT doing select(-c(DR, PB)) as in the 6-site version

# Fix any NA in disturbance (e.g. 2020 COVID closure)
HumanDisturbance.wide[HumanDisturbance.wide$Year == 2020, "DP"] <- 0
HumanDisturbance.wide[is.na(HumanDisturbance.wide)] <- 0

HumanDisturbance.wide <- HumanDisturbance.wide %>%
  rename_with(~ paste0("Dist_", .), -Year)

cat("Disturbance columns (should include DR and PB):\n")
print(names(HumanDisturbance.wide))

# ── MOCI data ─────────────────────────────────────────────────────────────────
MOCI <- read_csv("Data/CaliforniaMOCI.csv")

MOCI.dat <- MOCI %>%
  mutate(mean_value = (`North California (38-42N)` + `Central California (34.5-38N)`) / 2) %>%
  dplyr::select(Year, Season, mean_value) %>%
  mutate(Year = ifelse(Season == "OND", Year + 1, Year)) %>%
  dplyr::filter(Season %in% c("JFM", "AMJ", "OND")) %>%
  pivot_wider(names_from = Season, values_from = mean_value) %>%
  arrange(Year) %>%
  select(Year, JFM, AMJ, OND) %>%
  rename_with(~ paste0("MOCI_", .), -Year)

# ── Elephant seal data ────────────────────────────────────────────────────────
eSeal <- read_excel("Data/Eseal_1981-2025_BySubsite.xlsx")

eSeal_max_imm <- eSeal %>%
  mutate(Year = lubridate::year(StartDate)) %>%
  dplyr::group_by(Year, SubSiteName, MatureCode) %>%
  dplyr::summarise(MaxCount = max(Count, na.rm = TRUE), .groups = "drop") %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(eSeal_Sum_Imm_MaxCount = sum(MaxCount))

# ── Combine all covariates ─────────────────────────────────────────────────────
covariates_8site <- MOCI.dat %>%
  left_join(HumanDisturbance.wide, by = "Year") %>%
  left_join(coyote_wide, by = "Year") %>%
  left_join(eSeal_max_imm, by = "Year") %>%
  dplyr::filter(Year > 1996, Year < 2026)

cat("\nCovariate matrix dimensions:", nrow(covariates_8site), "×",
    ncol(covariates_8site) - 1, "covariates\n")

# ── Transpose and standardise ─────────────────────────────────────────────────
cov_t_8site <- t(covariates_8site)
cov_t_8site <- cov_t_8site[-1, ]   # remove Year row

# KEY DIFFERENCE: Do NOT remove Coyote_DR or Coyote_PB
# (They are all zeros, confirmed, but kept for structural completeness)
# Verify they are zero:
coy_dr <- cov_t_8site["Coyote_DR", ]
coy_pb <- cov_t_8site["Coyote_PB", ]
cat("\nCoyote_DR values (should all be zero):", unique(as.numeric(coy_dr)), "\n")
cat("Coyote_PB values (should all be zero):", unique(as.numeric(coy_pb)), "\n")

# Replace any remaining NAs with 0
cov_t_8site[is.na(cov_t_8site)] <- 0

# Standardise (z-score) rows with variance > 0
cov_t_scaled_8site <- cov_t_8site
rows_to_scale <- apply(cov_t_8site, 1, sd) > 0
cov_t_scaled_8site[rows_to_scale, ] <- t(scale(t(cov_t_8site[rows_to_scale, ])))

cat("\nScaled covariate rows (", sum(rows_to_scale), "of",
    nrow(cov_t_scaled_8site), "):\n")
cat("Row means (should be ~0 for scaled rows):\n")
print(round(rowMeans(cov_t_scaled_8site), 3))

cat("\nFinal covariate matrix: ", nrow(cov_t_scaled_8site), "×",
    ncol(cov_t_scaled_8site), "\n")
cat("Covariate order (columns in C matrix):\n")
for (i in seq_len(nrow(cov_t_scaled_8site))) {
  cat(sprintf("  %2d: %s\n", i, rownames(cov_t_scaled_8site)[i]))
}
