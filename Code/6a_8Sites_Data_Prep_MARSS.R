# ============================================================================
# DATA PREP FOR 8-SITE MARSS ANALYSIS (1997-2025)
# Includes DR (Drakes) and PB (Point Bonita) — haul-out-only sites
#
# IMPORTANT: This file is a COMPANION to the 6-site IPM analysis.
# Do NOT replace 4a_data_prep.R — that file feeds the 6-site IPM.
# This file produces: dat_8site, years_8site for use in 4cc_models_8sites.R
# ============================================================================

library(MARSS)
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyverse)
library(lubridate)
library(plyr)

# ── Read raw data ─────────────────────────────────────────────────────────────
Phoca <- read_excel("Data/1996_2025_Phocadata.xls")
Phoca <- Phoca[-c(3:4, 8:10)]

Phoca$Date2  <- ymd(Phoca$Date)
Phoca$Year   <- year(Phoca$Date2)
Phoca$Julian <- yday(Phoca$Date2)

# Classify season
Phoca$Age    <- ifelse(Phoca$Julian > 155, "MOLTING", Phoca$Age)
Phoca$Season <- ifelse(Phoca$Julian <= 140 & Phoca$Julian >= 105, "PUPPING", "MOLTING")

# Clean data
Phoca <- Phoca %>% dplyr::filter(Age != "DEADPUP" & Age != "DEADADULT")
Phoca$Age[Phoca$Age == "HPUP"] <- "PUP"
Phoca$Yearf <- as.factor(Phoca$Year)

all_data <- Phoca
all_data$Subsite <- ifelse(all_data$Subsite == "PR", "PRH", all_data$Subsite)

# ── Maximum count per year/site/age (same as 6-site analysis) ─────────────────
top1_all_data <- all_data %>%
  dplyr::group_by(Year, Subsite, Age) %>%
  slice_max(Count, n = 1) %>%
  dplyr::filter(Age != "PUP" | Season != "MOLTING")

# ── MARSS data matrix ─────────────────────────────────────────────────────────
all_data.MARSS <- top1_all_data

all_data.MARSS <- all_data.MARSS[, c(7, 2, 4, 5)]

# Log-transform counts
# Key difference from 6-site: DR and PB pup counts are near-zero or zero.
# Setting zero pup counts to NA avoids log(0) = -Inf and correctly tells
# MARSS that the pup state is unobserved at those sites.
all_data.MARSS <- all_data.MARSS %>%
  mutate(Count = case_when(
    Count == 0 ~ NA_real_,   # log(0) undefined — treat as missing
    TRUE       ~ log(Count)
  ))

all_data.MARSS <- distinct(all_data.MARSS)

all_data.MARSS$Subsite_Age <- paste0(all_data.MARSS$Subsite, "_", all_data.MARSS$Age)

# Filter to 1997 onwards
all_data.MARSS <- all_data.MARSS %>% dplyr::filter(Year > 1996)

# ── Missing data imputation ───────────────────────────────────────────────────
# BL 1997 Molt missing — use 1999 value (stable non-ENSO year, same as 6-site analysis)
addBL97molt <- data.frame(Year = 1997, Subsite = "BL", Age = "MOLTING",
                          Count = 5.7745515, Subsite_Age = "BL_MOLTING",
                          stringsAsFactors = FALSE)

# Check if DR and PB have 1997 data — add NA rows if missing so matrix is complete
all_sites <- c("BL", "DE", "DP", "DR", "PB", "PRH", "TB", "TP")
all_ages  <- c("ADULT", "MOLTING", "PUP")
all_combos <- expand.grid(Year = 1997, Subsite = all_sites, Age = all_ages,
                          stringsAsFactors = FALSE) %>%
  mutate(Subsite_Age = paste0(Subsite, "_", Age), Count = NA_real_)

existing_1997 <- all_data.MARSS %>% dplyr::filter(Year == 1997)
missing_1997  <- all_combos %>%
  dplyr::filter(!Subsite_Age %in% existing_1997$Subsite_Age) %>%
  dplyr::filter(!Subsite_Age == "BL_MOLTING")  # handled above

all_data.MARSS <- bind_rows(addBL97molt, missing_1997, all_data.MARSS) %>%
  arrange(Year, Subsite, Age)

# ── NOTE: DR and PB are haul-out-only sites (no pups) ─────────────────────────
# Their _PUP rows will be NA throughout — MARSS treats NA as missing data,
# so those rows inform nothing about pup parameters. This is correct behavior.
# Ecological rationale: DR and PB are included to characterize adult and
# molting seal abundance trends at non-breeding sites.
cat("\nNA pup rows for DR and PB (expected):\n")
all_data.MARSS %>%
  dplyr::filter(Subsite %in% c("DR", "PB"), Age == "PUP") %>%
  summarise(n_NA = sum(is.na(Count)), n_total = n()) %>%
  print()

# ── Reshape to wide format ────────────────────────────────────────────────────
all_data.MARSS <- all_data.MARSS[, -c(2:3)]
all_data.MARSS$Count <- as.numeric(all_data.MARSS$Count)

all_data.MARSS.wide_8site <- all_data.MARSS %>%
  pivot_wider(names_from = Subsite_Age, values_from = Count)

# MARSS requires time across columns
dat_8site <- t(all_data.MARSS.wide_8site)
years_8site <- as.numeric(dat_8site[1, ])
dat_8site   <- dat_8site[2:nrow(dat_8site), ]
dat_8site   <- matrix(as.numeric(dat_8site),
                      nrow = nrow(dat_8site),
                      ncol = ncol(dat_8site),
                      dimnames = dimnames(dat_8site))

cat("\nData matrix dimensions (should be 24 × 29):\n")
print(dim(dat_8site))

cat("\nSite-class rows in matrix:\n")
print(rownames(dat_8site))

cat("\nNA counts per row (DR_PUP and PB_PUP should be all NA):\n")
print(rowSums(is.na(dat_8site)))

# ── Summary table ─────────────────────────────────────────────────────────────
cat("\n── Site summary: adult and molt counts only ──────────────────────────────\n")
all_data.MARSS.wide_8site %>%
  dplyr::select(Year, contains("DR"), contains("PB")) %>%
  print(n = 5)
