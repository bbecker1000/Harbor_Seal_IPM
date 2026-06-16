# ============================================================================
# 02_covariates_6site.R
# ----------------------------------------------------------------------------
# COVARIATE PREP — 6 BREEDING SITES (1997–2025)
# Feeds BOTH the 6-site MARSS framework (03_marss6_models.R) AND the IPM
# (07_ipm_run.R -> prepare_real_data_for_ipm_v3.2()).
#
# Produces in the global environment:
#   cov_t_scaled : 16 x 29 numeric matrix (covariates x years), z-scored
#
# Row order is CANONICAL and asserted (the IPM and MARSS C-matrix both index
# these rows positionally, so order is enforced here, not assumed downstream):
#    1  MOCI_JFM
#    2  MOCI_AMJ
#    3  MOCI_OND
#    4  Dist_BL     5  Dist_DE     6  Dist_DP
#    7  Dist_PRH    8  Dist_TB     9  Dist_TP
#   10  Coyote_BL  11  Coyote_DE  12  Coyote_DP
#   13  Coyote_PRH 14  Coyote_TB  15  Coyote_TP
#   16  eSeal_Sum_Imm_MaxCount
#
# Prereq: source("Code/01_data_prep_6site.R")  (for `years` alignment check)
# Next:   source("Code/03_marss6_models.R")
# ============================================================================

library(tidyverse)
library(readxl)
library(lubridate)

dir.create("Output", showWarnings = FALSE)

SITES_6 <- c("BL", "DE", "DP", "PRH", "TB", "TP")

CANON_COV_ROWS <- c(
  "MOCI_JFM", "MOCI_AMJ", "MOCI_OND",
  paste0("Dist_",   SITES_6),
  paste0("Coyote_", SITES_6),
  "eSeal_Sum_Imm_MaxCount"
)

# ── Coyote: survey-rate with a 3-year weighted lag, per site ────────────────
coyote <- read_excel("Data/CoyoteSightings_2025.xlsx")

coyote_rate <- coyote %>%
  mutate(rate = `Number of days with coyote sightings` /
           `Total number of monitoring surveys`) %>%
  arrange(Site, Year) %>%
  group_by(Site) %>%
  mutate(
    rate_t1 = lag(rate, 1),
    rate_t2 = lag(rate, 2),
    has_t1  = !is.na(rate_t1),
    has_t2  = !is.na(rate_t2),
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
  ungroup() %>%
  select(Year, Site, weighted_rate)

coyote_wide <- coyote_rate %>%
  pivot_wider(names_from = Site, values_from = weighted_rate)

# Pre-survey years (1996–1999) confirmed zero coyote at all sites.
# new_rows fixes the SITE COLUMN ORDER; bind_rows then aligns coyote_wide to it.
new_rows <- tibble(Year = 1996:1999,
                   BL = 0, DE = 0, DP = 0, DR = 0, PB = 0,
                   PRH = 0, TB = 0, TP = 0)

coyote_wide <- bind_rows(new_rows, coyote_wide) %>%
  arrange(Year) %>%
  rename_with(~ paste0("Coyote_", .), -Year)

# ── Human disturbance: events / surveys, per site ───────────────────────────
HumanDisturbance <- read_excel("Data/HumanDisturbanceRate_1996To2025.xlsx")
HumanDisturbance <- HumanDisturbance[, -5]
HumanDisturbance <- HumanDisturbance %>%
  mutate(DistRate = SumOfDisturbanceCount / NSurveys) %>%
  select(-SumOfDisturbanceCount, -NSurveys)

HumanDisturbance.wide <- HumanDisturbance %>%
  pivot_wider(names_from = SiteCode, values_from = DistRate) %>%
  dplyr::filter(Year > 1996) %>%
  rename_with(~ paste0("Dist_", .), -Year)

# 2020 DP closure (COVID): no surveys -> set rate to 0, BY NAME (fixes P12).
# (The blanket NA->0 below would also catch it; doing it explicitly documents
#  the decision and is robust to row/column reordering.)
if ("Dist_DP" %in% names(HumanDisturbance.wide))
  HumanDisturbance.wide$Dist_DP[HumanDisturbance.wide$Year == 2020] <- 0

# ── MOCI: seasonal index; OND lead-shifted to align with following pup year ──
MOCI <- read_csv("Data/CaliforniaMOCI.csv", show_col_types = FALSE)
MOCI.dat <- MOCI %>%
  mutate(mean_value = (`North California (38-42N)` +
                         `Central California (34.5-38N)`) / 2) %>%
  select(Year, Season, mean_value) %>%
  mutate(Year = ifelse(Season == "OND", Year + 1, Year)) %>%   # OND -> next yr
  dplyr::filter(Season %in% c("JFM", "AMJ", "OND")) %>%
  pivot_wider(names_from = Season, values_from = mean_value) %>%
  arrange(Year) %>%
  select(Year, JFM, AMJ, OND) %>%
  rename_with(~ paste0("MOCI_", .), -Year)

# ── Elephant seal: annual summed max immature+mature count (all PR sites) ────
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

# ── Assemble rows in EXPLICIT canonical order (fixes P10) ───────────────────
# Pull each row by name rather than trusting join/pivot column order. Any
# missing name errors out here instead of silently shifting a row.
pull_row <- function(name) {
  if (!name %in% names(covariates))
    stop("Covariate column not found: ", name,
         "\n  Available: ", paste(names(covariates), collapse = ", "))
  as.numeric(covariates[[name]])
}

cov_t <- rbind(
  MOCI_JFM   = pull_row("MOCI_JFM"),
  MOCI_AMJ   = pull_row("MOCI_AMJ"),
  MOCI_OND   = pull_row("MOCI_OND"),
  Dist_BL    = pull_row("Dist_BL"),
  Dist_DE    = pull_row("Dist_DE"),
  Dist_DP    = pull_row("Dist_DP"),
  Dist_PRH   = pull_row("Dist_PRH"),
  Dist_TB    = pull_row("Dist_TB"),
  Dist_TP    = pull_row("Dist_TP"),
  Coyote_BL  = pull_row("Coyote_BL"),
  Coyote_DE  = pull_row("Coyote_DE"),
  Coyote_DP  = pull_row("Coyote_DP"),
  Coyote_PRH = pull_row("Coyote_PRH"),
  Coyote_TB  = pull_row("Coyote_TB"),
  Coyote_TP  = pull_row("Coyote_TP"),
  eSeal_Sum_Imm_MaxCount = pull_row("eSeal_Sum_Imm_MaxCount")
)
colnames(cov_t) <- cov_years

# ── Single scaling pass (fixes the duplicate-scale bug) ─────────────────────
# Replace NAs with 0 on the RAW scale (so they land at the row mean after
# z-scoring), then z-score only rows with non-zero variance. Coyote rows that
# are all-zero (none here for the 6 sites) would be left unscaled by design.
cov_t[is.na(cov_t)] <- 0
cov_t_scaled <- cov_t
rows_to_scale <- apply(cov_t, 1, sd) > 0
cov_t_scaled[rows_to_scale, ] <- t(scale(t(cov_t[rows_to_scale, ])))

# ── Assertions: shape, order, year alignment (fixes P9/P10 at source) ───────
stopifnot(
  identical(rownames(cov_t_scaled), CANON_COV_ROWS),
  nrow(cov_t_scaled) == 16L
)
if (exists("years")) {
  if (!identical(as.numeric(cov_years), as.numeric(years)))
    stop("Covariate years do not match data years from 01_data_prep_6site.R\n",
         "  cov: ", paste(range(cov_years), collapse = "-"),
         " | dat: ", paste(range(years), collapse = "-"))
} else {
  message("Note: `years` not in environment — run 01_data_prep_6site.R to ",
          "enable the year-alignment check.")
}

# ── Diagnostics ─────────────────────────────────────────────────────────────
cat("\n── 6-site covariate matrix ─────────────────────────────────────────\n")
cat(sprintf("cov_t_scaled: %d covariates x %d years (%d-%d)\n",
            nrow(cov_t_scaled), ncol(cov_t_scaled),
            min(cov_years), max(cov_years)))
cat("Row order (canonical):\n")
for (i in seq_len(nrow(cov_t_scaled)))
  cat(sprintf("  %2d  %s\n", i, rownames(cov_t_scaled)[i]))
cat("\nRow means (scaled rows ~0):\n")
print(round(rowMeans(cov_t_scaled), 3))

# Optional: persist alongside dat (07_ipm_run.R also bundles these).
# saveRDS(cov_t_scaled, "Output/cov_t_scaled_6site.rds")

cat("\nObjects created: cov_t_scaled\n")
cat("Next: source(\"Code/03_marss6_models.R\")\n")