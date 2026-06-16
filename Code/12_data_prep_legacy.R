# ============================================================================
# 7a_data_prep_legacy_MARSS.R
# ============================================================================
# DATA PREP — LEGACY + MODERN COMBINED TIME SERIES (1975–2025)
# No-covariate MARSS for long-term trend characterization
#
# Sites: BL, DE, DP, PRH, TB, TP (6 IPM breeding sites)
#        DR and PB excluded — no legacy count data before 1997
#
# Produces: dat_legacy, years_legacy, state_names_legacy
# ============================================================================

library(MARSS)
library(readxl)
library(dplyr)
library(tidyverse)
library(lubridate)

dir.create("Output", showWarnings = FALSE)

LEGACY_SITES <- c("BL", "DE", "DP", "PRH", "TB", "TP")

# ── 1. Read and reshape legacy data (1975–1999 max counts, wide format) ───────
hseal_legacy <- read_excel("Data/1975-1999_max_counts.xlsx")

# Reclassify adults surveyed during molting season
hseal_legacy$Age <- ifelse(hseal_legacy$Season == "MOLTING", "MOLTING", hseal_legacy$Age)

hseal_legacy_long <- hseal_legacy %>%
  pivot_longer(cols      = -c(Year, Age, Season),
               names_to  = "Subsite",
               values_to = "Count") %>%
  filter(Year < 1996,                    # 1996+ comes from modern data
         Subsite %in% LEGACY_SITES)

cat("Legacy records:", nrow(hseal_legacy_long),
    "| Years:", min(hseal_legacy_long$Year), "–", max(hseal_legacy_long$Year), "\n")

# ── 2. Read and process modern data (1996–2025) ───────────────────────────────
Phoca <- read_excel("Data/1996_2025_Phocadata.xls")
Phoca <- Phoca[-c(3:4, 8:10)]

Phoca$Date2  <- ymd(Phoca$Date)
Phoca$Year   <- year(Phoca$Date2)
Phoca$Julian <- yday(Phoca$Date2)

Phoca$Age    <- ifelse(Phoca$Julian > 155, "MOLTING", Phoca$Age)
Phoca$Season <- ifelse(Phoca$Julian <= 140 & Phoca$Julian >= 105, "PUPPING", "MOLTING")

Phoca <- Phoca %>%
  filter(Age != "DEADPUP", Age != "DEADADULT") %>%
  mutate(Age     = ifelse(Age == "HPUP", "PUP", Age),
         Subsite = ifelse(Subsite == "PR", "PRH", Subsite)) %>%
  filter(Subsite %in% LEGACY_SITES)

# ── 3. Combine and take annual maximum per site × age class ──────────────────
all_long <- bind_rows(
  hseal_legacy_long %>% select(Year, Subsite, Age, Season, Count),
  Phoca             %>% select(Year, Subsite, Age, Season, Count)
)

top1_legacy <- all_long %>%
  group_by(Year, Subsite, Age) %>%
  slice_max(Count, n = 1, na_rm = TRUE) %>%
  filter(!(Age == "PUP" & Season == "MOLTING")) %>%  # remove pups in molt season
  ungroup() %>%
  distinct()

# ── 4. Log-transform; zeros and NAs → NA ─────────────────────────────────────
top1_legacy <- top1_legacy %>%
  mutate(Count = case_when(
    is.na(Count) | Count <= 0 ~ NA_real_,
    TRUE                      ~ log(Count)
  ))

top1_legacy$Subsite_Age <- paste0(top1_legacy$Subsite, "_", top1_legacy$Age)

# ── 5. BL 1997 Molt imputation (same as 4a) ──────────────────────────────────
if (!any(top1_legacy$Subsite_Age == "BL_MOLTING" &
         top1_legacy$Year == 1997 & !is.na(top1_legacy$Count))) {
  top1_legacy <- bind_rows(
    top1_legacy,
    data.frame(Year = 1997, Subsite = "BL", Age = "MOLTING",
               Season = "MOLTING", Count = 5.7745515,
               Subsite_Age = "BL_MOLTING", stringsAsFactors = FALSE)
  )
  cat("BL 1997 Molt imputation added (log scale: 5.775)\n")
}

# ── 6. Determine first usable year (≥ 1 site observed) ───────────────────────
first_year <- top1_legacy %>%
  filter(!is.na(Count)) %>%
  pull(Year) %>%
  min()
cat("First year with any observation:", first_year, "\n")

top1_legacy <- top1_legacy %>% filter(Year >= first_year)

# Coverage summary
coverage <- top1_legacy %>%
  filter(!is.na(Count)) %>%
  group_by(Year) %>%
  summarise(n_sites  = n_distinct(Subsite),
            n_states = n(), .groups = "drop")
cat("\nSite coverage over time:\n")
print(coverage, n = 20)

# ── 7. Reshape to MARSS wide format (states × years) ─────────────────────────
legacy_wide <- top1_legacy %>%
  select(Year, Subsite_Age, Count) %>%
  group_by(Year, Subsite_Age) %>%
  slice(1) %>%
  ungroup() %>%
  pivot_wider(names_from = Subsite_Age, values_from = Count) %>%
  arrange(Year)

dat_legacy   <- t(legacy_wide)
years_legacy <- as.numeric(dat_legacy[1, ])
dat_legacy   <- dat_legacy[2:nrow(dat_legacy), ]
dat_legacy   <- matrix(as.numeric(dat_legacy),
                       nrow      = nrow(dat_legacy),
                       ncol      = ncol(dat_legacy),
                       dimnames  = dimnames(dat_legacy))

state_names_legacy <- rownames(dat_legacy)
TT_legacy <- ncol(dat_legacy)

# ── 8. Diagnostics ───────────────────────────────────────────────────────────
cat("\n── Data matrix: ", nrow(dat_legacy), "states ×", TT_legacy, "years",
    "(", first_year, "–", max(years_legacy), ") ────\n")
cat("States:\n"); print(state_names_legacy)
cat("\nNon-NA observations per state:\n")
print(rowSums(!is.na(dat_legacy)))
cat("\nNA counts per state:\n")
print(rowSums(is.na(dat_legacy)))

# Quick coverage heatmap
coverage_long <- as.data.frame(dat_legacy) %>%
  tibble::rownames_to_column("State") %>%
  pivot_longer(-State, names_to = "Year_col", values_to = "val") %>%
  mutate(Year     = years_legacy[as.integer(gsub("V","",Year_col))],
         Observed = !is.na(val)) %>%
  separate(State, into = c("Site", "Class"), sep = "_", remove = FALSE)

cat("\n── Data coverage summary ────────────────────────────────────────────────\n")
print(coverage_long %>%
        group_by(Site, Class) %>%
        summarise(First_year = min(Year[Observed], na.rm=TRUE),
                  Last_year  = max(Year[Observed], na.rm=TRUE),
                  n_years    = sum(Observed), .groups = "drop"))

cat("\nObjects created: dat_legacy, years_legacy, state_names_legacy\n")
cat("Source 7b_models_legacy_MARSS.R to run analysis.\n")