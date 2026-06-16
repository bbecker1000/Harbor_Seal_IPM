# ============================================================================
# 08_data_prep_8site.R
# ----------------------------------------------------------------------------
# DATA PREP — 8 SITES incl. haul-out-only DR & PB (1997–2025)
# Supplemental MARSS analysis (Appendix A). Companion to the 6-site IPM;
# does NOT feed the IPM. DR/PB have no meaningful pup production -> pup states
# are carried as NA throughout (MARSS treats NA as missing).
#
# Produces:
#   dat_8site   : 24 x 29 numeric matrix (states x years), log-scale, NA=missing
#   years_8site : numeric vector 1997..2025
#
# CANONICAL 24-row order (site-major; ADULT, MOLTING, PUP), asserted:
#   BL_*, DE_*, DP_*, DR_*, PB_*, PRH_*, TB_*, TP_*
# The 24-state C and U matrices in 10_marss8_models.R depend on this order.
#
# Next: source("Code/09_covariates_8site.R")
# ============================================================================

library(tidyverse)
library(readxl)
library(lubridate)

dir.create("Output", showWarnings = FALSE)

SITES_8       <- c("BL", "DE", "DP", "DR", "PB", "PRH", "TB", "TP")
CLASS_ORDER   <- c("ADULT", "MOLTING", "PUP")
BL97_MOLT_LOG <- 5.7745515
canonical_states <- as.vector(t(outer(SITES_8, CLASS_ORDER, paste, sep = "_")))

# ── Read raw PRNS survey data ───────────────────────────────────────────────
Phoca <- read_excel("Data/1996_2025_Phocadata.xls")
Phoca <- Phoca[, -c(3:4, 8:10)]

Phoca <- Phoca %>%
  mutate(
    Date2   = ymd(Date),
    Year    = year(Date2),
    Julian  = yday(Date2),
    Age     = ifelse(Julian > 155, "MOLTING", Age),
    Season  = ifelse(Julian <= 140 & Julian >= 105, "PUPPING", "MOLTING"),
    Subsite = ifelse(Subsite == "PR", "PRH", Subsite)
  ) %>%
  dplyr::filter(Age != "DEADPUP", Age != "DEADADULT") %>%
  mutate(Age = ifelse(Age == "HPUP", "PUP", Age))

# ── Annual maximum per Year × Subsite × Age ─────────────────────────────────
top1 <- Phoca %>%
  dplyr::group_by(Year, Subsite, Age) %>%
  dplyr::slice_max(Count, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!(Age == "PUP" & Season == "MOLTING")) %>%
  dplyr::select(Year, Subsite, Age, Count)

# ── Log-transform with zero/NA guard (0 pup counts at DR/PB -> NA) ──────────
top1 <- top1 %>%
  mutate(
    Count       = ifelse(is.na(Count) | Count <= 0, NA_real_, log(Count)),
    Subsite_Age = paste0(Subsite, "_", Age)
  ) %>%
  dplyr::filter(Year > 1996) %>%
  dplyr::distinct(Year, Subsite, Age, .keep_all = TRUE)

# ── BL 1997 MOLTING imputation (numeric, log scale, once) ───────────────────
if (!any(top1$Subsite_Age == "BL_MOLTING" &
         top1$Year == 1997 & !is.na(top1$Count))) {
  top1 <- bind_rows(
    top1,
    tibble(Year = 1997, Subsite = "BL", Age = "MOLTING",
           Count = BL97_MOLT_LOG, Subsite_Age = "BL_MOLTING")
  )
  message("Imputed BL 1997 MOLTING = ", BL97_MOLT_LOG, " (log scale).")
}

# ── Ensure all 24 site-class combos exist in 1997 (fill missing as NA) ──────
# Guarantees a complete 24-row matrix even where a site-class was unobserved
# in the first year (e.g. DR/PB pups, or sites whose surveys begin later).
all_combos_1997 <- expand_grid(Year = 1997, Subsite = SITES_8, Age = CLASS_ORDER) %>%
  mutate(Subsite_Age = paste0(Subsite, "_", Age), Count = NA_real_)
have_1997 <- top1 %>% dplyr::filter(Year == 1997) %>% pull(Subsite_Age)
top1 <- bind_rows(
  top1,
  all_combos_1997 %>% dplyr::filter(!Subsite_Age %in% have_1997)
) %>% arrange(Year, Subsite, Age)

# ── Reshape to wide (years × states); transpose to MARSS layout ─────────────
wide <- top1 %>%
  dplyr::select(Year, Subsite_Age, Count) %>%
  dplyr::group_by(Year, Subsite_Age) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(names_from = Subsite_Age, values_from = Count) %>%
  dplyr::arrange(Year)

years_8site <- wide$Year
dat_8site   <- t(as.matrix(dplyr::select(wide, -Year)))

# ── Enforce canonical order + assert completeness ───────────────────────────
missing_states <- setdiff(canonical_states, rownames(dat_8site))
if (length(missing_states) > 0)
  stop("Missing state rows after reshape: ",
       paste(missing_states, collapse = ", "))

dat_8site <- dat_8site[canonical_states, , drop = FALSE]
storage.mode(dat_8site) <- "numeric"

stopifnot(
  identical(rownames(dat_8site), canonical_states),
  nrow(dat_8site) == 24L,
  ncol(dat_8site) == length(years_8site)
)

# ── Diagnostics ─────────────────────────────────────────────────────────────
cat("\n── 8-site data matrix ──────────────────────────────────────────────\n")
cat(sprintf("dat_8site: %d states x %d years (%d-%d)\n",
            nrow(dat_8site), ncol(dat_8site),
            min(years_8site), max(years_8site)))
cat("\nNA counts per state (DR_PUP & PB_PUP should be all-NA = 29):\n")
print(rowSums(is.na(dat_8site)))

cat("\nObjects created: dat_8site, years_8site\n")
cat("Next: source(\"Code/09_covariates_8site.R\")\n")