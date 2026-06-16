# ============================================================================
# 01_data_prep_6site.R
# ----------------------------------------------------------------------------
# DATA PREP — 6 BREEDING SITES (1997–2025)
# Feeds BOTH the 6-site MARSS framework analysis (03_marss6_models.R) AND the
# Bayesian IPM (via 07_ipm_run.R -> prepare_real_data_for_ipm_v3.2()).
#
# Sites:   BL, DE, DP, PRH, TB, TP  (the 6 pupping/breeding sites)
# Classes: ADULT, MOLTING, PUP
#
# Produces in the global environment:
#   dat   : 18 x 29 numeric matrix (states x years), log-scale counts, NA = missing
#   years : numeric vector 1997..2025
#
# Row order is CANONICAL and asserted (site-major; ADULT, MOLTING, PUP):
#   BL_ADULT, BL_MOLTING, BL_PUP, DE_ADULT, ... , TP_PUP
# Downstream code (MARSS C-matrix rows; IPM seq(1,18,by=3) extraction) depends
# on this exact order, so it is enforced here rather than assumed later.
#
# Next: source("Code/02_covariates_6site.R")
# ============================================================================

library(tidyverse)
library(readxl)
library(lubridate)

dir.create("Output", showWarnings = FALSE)

# ── Configuration ───────────────────────────────────────────────────────────
SITES_6      <- c("BL", "DE", "DP", "PRH", "TB", "TP")
CLASS_ORDER  <- c("ADULT", "MOLTING", "PUP")
BL97_MOLT_LOG <- 5.7745515   # log-scale imputation for missing BL 1997 molt
# (1999 value; stable non-ENSO year)

# Canonical 18-row order: site-major, class order ADULT, MOLTING, PUP.
# outer() builds a 6x3 site×class grid; t() then as.vector() reads it as
# BL_ADULT, BL_MOLTING, BL_PUP, DE_ADULT, ...  (this is what the IPM expects).
canonical_states <- as.vector(t(outer(SITES_6, CLASS_ORDER, paste, sep = "_")))

# ── Read raw PRNS survey data ───────────────────────────────────────────────
# NOTE: the column drop below is positional and matches the known layout of
# 1996_2025_Phocadata.xls. If the source spreadsheet layout changes, update
# this line (keeps Date, Subsite, Age, Count).
Phoca <- read_excel("Data/1996_2025_Phocadata.xls")
Phoca <- Phoca[, -c(3:4, 8:10)]

Phoca <- Phoca %>%
  mutate(
    Date2   = ymd(Date),
    Year    = year(Date2),
    Julian  = yday(Date2),
    # adults surveyed after ~June 5 are reclassified as molting
    Age     = ifelse(Julian > 155, "MOLTING", Age),
    Season  = ifelse(Julian <= 140 & Julian >= 105, "PUPPING", "MOLTING"),
    Subsite = ifelse(Subsite == "PR", "PRH", Subsite)       # name fix
  ) %>%
  dplyr::filter(Age != "DEADPUP", Age != "DEADADULT") %>%
  mutate(Age = ifelse(Age == "HPUP", "PUP", Age)) %>%       # HPUP typo -> PUP
  dplyr::filter(Subsite %in% SITES_6)

# ── Annual maximum count per Year × Subsite × Age ───────────────────────────
top1 <- Phoca %>%
  dplyr::group_by(Year, Subsite, Age) %>%
  dplyr::slice_max(Count, n = 1, with_ties = FALSE) %>%     # one row per group
  dplyr::ungroup() %>%
  dplyr::filter(!(Age == "PUP" & Season == "MOLTING")) %>%  # drop molt-season pups
  dplyr::filter(Year > 1996) %>%
  dplyr::distinct(Year, Subsite, Age, .keep_all = TRUE)

# ── Log-transform with a zero/NA guard (fixes P12) ──────────────────────────
# Counts of 0 (or NA) become NA so log() never produces -Inf, and MARSS/IPM
# correctly treat them as unobserved rather than as an impossible observation.
top1 <- top1 %>%
  mutate(
    Count       = ifelse(is.na(Count) | Count <= 0, NA_real_, log(Count)),
    Subsite_Age = paste0(Subsite, "_", Age)
  )

# ── BL 1997 MOLTING imputation — numeric, log scale, added once (fixes P12) ──
has_bl97 <- any(top1$Subsite_Age == "BL_MOLTING" &
                  top1$Year == 1997 & !is.na(top1$Count))
if (!has_bl97) {
  top1 <- bind_rows(
    top1,
    tibble(Year = 1997, Subsite = "BL", Age = "MOLTING", Season = "MOLTING",
           Count = BL97_MOLT_LOG, Subsite_Age = "BL_MOLTING")
  )
  message("Imputed BL 1997 MOLTING = ", BL97_MOLT_LOG, " (log scale).")
}

# ── Reshape to wide (years × states), then transpose to MARSS/IPM layout ────
wide <- top1 %>%
  dplyr::select(Year, Subsite_Age, Count) %>%
  dplyr::group_by(Year, Subsite_Age) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(names_from = Subsite_Age, values_from = Count) %>%
  dplyr::arrange(Year)

years <- wide$Year
dat   <- t(as.matrix(dplyr::select(wide, -Year)))   # states across rows

# ── Enforce canonical row order + assert completeness (fixes P11) ───────────
missing_states <- setdiff(canonical_states, rownames(dat))
if (length(missing_states) > 0)
  stop("Missing state rows after reshape: ",
       paste(missing_states, collapse = ", "),
       "\n  (a site-class was never observed; cannot build the 18-row matrix)")

dat <- dat[canonical_states, , drop = FALSE]
storage.mode(dat) <- "numeric"

stopifnot(
  identical(rownames(dat), canonical_states),
  nrow(dat) == 18L,
  ncol(dat) == length(years)
)

# ── Diagnostics ─────────────────────────────────────────────────────────────
cat("\n── 6-site data matrix ──────────────────────────────────────────────\n")
cat(sprintf("dat: %d states x %d years (%d-%d)\n",
            nrow(dat), ncol(dat), min(years), max(years)))
cat("Row order (canonical):\n  ", paste(rownames(dat), collapse = ", "), "\n")
cat("\nNon-NA observations per state:\n")
print(rowSums(!is.na(dat)))
cat("\nNA counts per state:\n")
print(rowSums(is.na(dat)))

# Optional: persist the IPM/MARSS input bundle (07_ipm_run.R also saves this).
# saveRDS(list(dat = dat, years = years), "Output/ipm_dat_6site.rds")

cat("\nObjects created: dat, years\n")
cat("Next: source(\"Code/02_covariates_6site.R\")\n")

