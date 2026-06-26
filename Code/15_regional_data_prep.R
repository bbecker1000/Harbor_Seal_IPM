# ============================================================================
# 15_regional_data_prep.R
# ----------------------------------------------------------------------------
# REGIONAL HARBOR SEAL IPM — DATA PREPARATION (2005–2025)
#
# Combines:
#   (1) Marin county sites extracted from 1996_2025_Phocadata.xls
#       (same max logic as 01_data_prep_6site.R)
#   (2) Regional non-Marin max counts from RegionalPhoca_2005To2025_ForBecker.xlsx
#   (3) MOCI seasonal indices from CaliforniaMOCI.csv
#       (average of North + Central California, matching Marin IPM convention)
#
# COUNTY GROUPINGS (ecological, not administrative):
#   C1 = Marin      : BL, DE, DP, PRH, TB, TP           is_bay = 0
#   C2 = SF Coast   : Castro Rocks, Yerba Buena, Alcatraz  is_bay = 0
#   C3 = SF Estuary : Mowry Slough, Newark Slough          is_bay = 1
#   C4 = San Mateo  : Fitzgerald, Cowell Ranch, Pebble Beach,
#                     Pescadero, Purisima Creek, Point San Pedro  is_bay = 0
#   C5 = Sonoma     : Jenner, Sea Ranch, Fort Ross, Bodega MR  is_bay = 0
#   C6 = Mendocino  : Point Arena Lighthouse              is_bay = 0
#
# SITE TYPES:
#   is_breeder = 1 (Type B): full likelihood — adult + pup + molt
#   is_breeder = 0 (Type H): adult + molt only (no pup production observed)
#   Alcatraz is Type H. All other sites are Type B.
#
# 2020 GAP: No surveys conducted (COVID). Year 2020 is retained in the time
#   series as a LATENT year — Leslie matrix transitions occur normally but
#   all observation indicators are 0 (no likelihood contribution). This
#   correctly propagates state uncertainty through the gap without discarding
#   2019 or 2021 data.
#
# Produces in the global environment:
#   regional_sites  : tibble of site metadata (22 sites)
#   regional_counts : long-format tibble (site x year x class)
#   regional_moci   : tibble of MOCI covariates (2004–2025)
#   regional_stan   : Stan data list ready for 17_regional_ipm_model.R
#   years_regional  : integer vector 2005:2025 (T=21, index 16 = latent 2020)
#
# Next: source("Code/16_regional_simulate.R") to test model on simulated data,
#       then source("Code/17_regional_ipm_model.R") to fit real data.
# ============================================================================

library(tidyverse)
library(readxl)
library(lubridate)

dir.create("Output", showWarnings = FALSE)

# ── Constants ─────────────────────────────────────────────────────────────────
YEARS_REGIONAL <- 2005:2025          # T = 21; index 16 = 2020 (latent)
T_PROJ_REGIONAL <- 10
LATENT_YEAR <- 2020                  # no surveys; kept as unobserved state year

# ── SITE METADATA TABLE ───────────────────────────────────────────────────────
# Complete catalogue: 22 sites across 6 ecological county groups.
# county_id: integer 1–6 matching county groupings above.
# is_bay   : 1 = estuarine/bay-influenced prey; 0 = open coast
# is_breeder: 1 = Type B (pup likelihood active); 0 = Type H (adults/molt only)
# source   : "marin" (from raw XLS) | "regional" (from pre-processed XLSX)

regional_sites <- tribble(
  ~site_id, ~site_name,              ~county,     ~county_id, ~is_bay, ~is_breeder, ~source,     ~starts,
  # ── Marin (C=1, county_type=0: open coast) ───────────────────────────────
  1L,  "BL",                     "Marin",     1L,  0L,  1L,  "marin",    2005L,
  2L,  "DE",                     "Marin",     1L,  0L,  1L,  "marin",    2005L,
  3L,  "DP",                     "Marin",     1L,  0L,  1L,  "marin",    2005L,
  4L,  "PRH",                    "Marin",     1L,  0L,  1L,  "marin",    2005L,
  5L,  "TB",                     "Marin",     1L,  0L,  1L,  "marin",    2005L,
  6L,  "TP",                     "Marin",     1L,  0L,  1L,  "marin",    2005L,
  7L,  "DR",                     "Marin",     1L,  0L,  0L,  "marin",    2005L,
  8L,  "PB",                     "Marin",     1L,  0L,  0L,  "marin",    2005L,
  # ── Bay Mouth (C=2, county_type=1) ───────────────────────────────────────
  # Bay-mouth sites: strong tidal mixing but within SF Bay system; animals
  # exposed to Bay freshwater conditions during Nov-Jan implantation window.
  9L,  "Castro Rocks",           "Bay Mouth", 2L,  1L,  1L,  "regional", 2005L,
  10L, "Yerba Buena Island",     "Bay Mouth", 2L,  1L,  1L,  "regional", 2005L,
  11L, "Alcatraz",               "Bay Mouth", 2L,  1L,  0L,  "regional", 2005L,
  # ── South Bay (C=3, county_type=2) ───────────────────────────────────────
  # Deep estuarine sites: low winter salinity, freshwater-dominated prey;
  # most decoupled from coastal upwelling. Surveys ended after 2018.
  12L, "Mowry Slough",           "South Bay", 3L,  1L,  1L,  "regional", 2005L,
  13L, "Newark Slough",          "South Bay", 3L,  1L,  1L,  "regional", 2005L,
  # ── San Mateo (C=4, county_type=0: open coast) ───────────────────────────
  14L, "Fitzgerald",             "San Mateo", 4L,  0L,  1L,  "regional", 2005L,
  15L, "Cowell Ranch",           "San Mateo", 4L,  0L,  1L,  "regional", 2005L,
  16L, "Pebble Beach",           "San Mateo", 4L,  0L,  1L,  "regional", 2005L,
  17L, "Pescadero",              "San Mateo", 4L,  0L,  1L,  "regional", 2005L,
  18L, "Purisima Creek",         "San Mateo", 4L,  0L,  1L,  "regional", 2012L,
  19L, "Point San Pedro",        "San Mateo", 4L,  0L,  1L,  "regional", 2005L,
  # ── Sonoma (C=5, county_type=0: open coast) ──────────────────────────────
  20L, "Jenner",                 "Sonoma",    5L,  0L,  1L,  "regional", 2005L,
  21L, "Sea Ranch",              "Sonoma",    5L,  0L,  1L,  "regional", 2005L,
  22L, "Fort Ross",              "Sonoma",    5L,  0L,  1L,  "regional", 2006L,
  23L, "Bodega Marine Reserve",  "Sonoma",    5L,  0L,  1L,  "regional", 2010L,
  # ── Mendocino (C=6, county_type=0: open coast) ───────────────────────────
  24L, "Point Arena Lighthouse", "Mendocino", 6L,  0L,  1L,  "regional", 2019L
) |>
  mutate(county_label = factor(county,
                               levels = c("Marin","Bay Mouth","South Bay",
                                          "San Mateo","Sonoma","Mendocino")))

# County-level oceanographic type vector (indexed by county_id 1:C)
# 0 = open coast, 1 = bay mouth, 2 = south bay
COUNTY_TYPE <- c(0L, 1L, 2L, 0L, 0L, 0L)

S <- nrow(regional_sites)   # 24
C <- n_distinct(regional_sites$county_id)   # 6
T <- length(YEARS_REGIONAL)

cat(sprintf("Regional model: S=%d sites, C=%d counties, T=%d years (2005-2025)\n", S, C, T))
cat(sprintf("County types: %d coast, %d bay mouth, %d south bay\n",
            sum(COUNTY_TYPE == 0), sum(COUNTY_TYPE == 1), sum(COUNTY_TYPE == 2)))

# ── PART 1: MARIN COUNTS (from raw survey XLS) ────────────────────────────────
# Apply same max-count logic as 01_data_prep_6site.R.
# ADULT surveys: Julian ≤ 155 (breeding season, spring)
# MOLTING surveys: Julian > 155 (molt season, summer) — reclassified from ADULT
# PUP surveys: all seasons, but drop molt-season pups

marin_raw <- read_excel("Data/1996_2025_Phocadata.xls")
marin_raw <- marin_raw[, c("Date", "Subsite", "Age", "Count")]

MARIN_SITES <- c("BL", "DE", "DP", "PRH", "TB", "TP", "DR", "PB")

marin_processed <- marin_raw |>
  mutate(
    Date2   = ymd(Date),
    Year    = year(Date2),
    Julian  = yday(Date2),
    Age     = ifelse(Julian > 155, "MOLTING", Age),
    Season  = ifelse(Julian <= 140 & Julian >= 105, "PUPPING", "MOLTING"),
    Subsite = ifelse(Subsite == "PR", "PRH", Subsite)
  ) |>
  filter(
    Age     %in% c("ADULT", "MOLTING", "PUP"),
    Subsite %in% MARIN_SITES,
    Year    >= 2005
  ) |>
  filter(!(Age == "PUP" & Season == "MOLTING")) |>   # drop molt-season pups
  mutate(Count = ifelse(is.na(Count) | Count <= 0, NA_real_, Count))

# Annual maximum per site-year-class
marin_max <- marin_processed |>
  group_by(Year, Subsite, Age) |>
  slice_max(Count, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(Year, site_name = Subsite, age_class = Age, count = Count) |>
  mutate(age_class = recode(age_class,
                            "ADULT"   = "Breed",
                            "MOLTING" = "Molt",
                            "PUP"     = "Pup"))

cat(sprintf("Marin data: %d site-year-class records (2005-%d)\n",
            nrow(marin_max), max(marin_max$Year)))

# ── PART 2: REGIONAL NON-MARIN COUNTS (pre-processed max from XLSX) ───────────
regional_raw <- read_excel("Data/RegionalPhoca_2005To2025_ForBecker.xlsx")

# Recode AgeClass and handle "ND" as NA
regional_cleaned <- regional_raw |>
  mutate(
    count = suppressWarnings(as.numeric(Count)),   # "ND" → NA
    age_class = AgeClass                           # Breed, Molt, Pup — already correct
  ) |>
  select(Year, site_name = SiteName, age_class, count) |>
  filter(!site_name %in% MARIN_SITES)  # safety: exclude all Marin sites

# Verify all regional sites are in our metadata
check_sites <- setdiff(unique(regional_cleaned$site_name), regional_sites$site_name)
if (length(check_sites) > 0)
  warning("Sites in XLSX not in metadata: ", paste(check_sites, collapse=", "))

cat(sprintf("Regional data: %d site-year-class records (2005-%d)\n",
            nrow(regional_cleaned), max(regional_cleaned$Year)))

# ── PART 3: COMBINE ALL COUNTS ────────────────────────────────────────────────
regional_counts <- bind_rows(marin_max, regional_cleaned) |>
  left_join(regional_sites |> select(site_id, site_name, county, county_id,
                                     is_bay, is_breeder, starts),
            by = "site_name") |>
  arrange(site_id, Year, age_class)

# Validate: all site names matched
unmatched <- filter(regional_counts, is.na(site_id))
if (nrow(unmatched) > 0) stop("Unmatched sites: ", paste(unique(unmatched$site_name), collapse=", "))

# ── PART 4: MOCI (North + Central CA average; OND shifted +1 year) ────────────
moci_raw <- read_csv("Data/CaliforniaMOCI.csv", show_col_types = FALSE)

moci_proc <- moci_raw |>
  mutate(moci_avg = (`North California (38-42N)` + `Central California (34.5-38N)`) / 2) |>
  filter(Season %in% c("JFM", "AMJ", "OND")) |>
  mutate(Year_model = ifelse(Season == "OND", Year + 1, Year)) |>  # OND pre-shift
  select(Year_model, Season, moci_avg) |>
  pivot_wider(names_from = Season, values_from = moci_avg, names_prefix = "moci_") |>
  arrange(Year_model) |>
  rename(moci_jfm = moci_JFM, moci_amj = moci_AMJ, moci_ond = moci_OND)

# Restrict to study years + one prior year for OND shift
regional_moci <- moci_proc |>
  filter(Year_model %in% YEARS_REGIONAL) |>
  arrange(Year_model)

if (nrow(regional_moci) != T)
  stop(sprintf("MOCI has %d rows; expected %d (2005:2025)", nrow(regional_moci), T))

# Z-score MOCI (column-wise; within study period)
moci_scaled <- regional_moci |>
  mutate(across(starts_with("moci_"), ~ as.vector(scale(.x))))

cat(sprintf("MOCI data: %d years; JFM mean=%.3f, AMJ mean=%.3f, OND mean=%.3f\n",
            nrow(moci_scaled),
            mean(moci_scaled$moci_jfm), mean(moci_scaled$moci_amj), mean(moci_scaled$moci_ond)))

# ── PART 5: BUILD OBSERVATION ARRAYS ─────────────────────────────────────────
# Arrays: y_adult[S, T], y_pup[S, T], y_molt[S, T]
# Observation indicators: y_*_obs[S, T] = 1 if surveyed, 0 if not
#
# Three reasons for obs = 0:
#   (a) 2020 latent year — no surveys anywhere
#   (b) Site not yet started (before first year in metadata)
#   (c) Explicit ND (not surveyed) in a given year
#   (d) Type H site: y_pup_obs always 0 (is_breeder == 0)

# Initialize matrices
year_idx <- setNames(seq_along(YEARS_REGIONAL), YEARS_REGIONAL)

make_matrix <- function(fill = NA_real_) matrix(fill, nrow = S, ncol = T,
                                                dimnames = list(regional_sites$site_name,
                                                                YEARS_REGIONAL))
y_adult <- make_matrix(); y_pup <- make_matrix(); y_molt <- make_matrix()

# Fill observed counts (log scale for Normal likelihood, as in Marin IPM)
#
# ZERO CONVENTION: count = 0 means the site WAS surveyed but no animals were
# seen — this is real information and must be kept as an observation.
# NA (from "ND" in the XLSX or missing Marin years) means no survey occurred.
#
# log(0) is undefined under the log-normal model, so true zero counts use a
# lower-bound offset: log(0.5). This is equivalent to saying "we observed
# fewer than 1 animal" and constrains N * alpha * d to be near zero for that
# site-year while keeping the observation in the likelihood.

LOG_ZERO_OFFSET <- log(0.5)   # offset for true zero counts

n_zeros <- 0L  # track how many zero-count observations are present

for (i in seq_len(nrow(regional_counts))) {
  row <- regional_counts[i, ]
  s <- row$site_id
  t <- year_idx[as.character(row$Year)]
  if (is.na(t)) next
  if (row$Year == LATENT_YEAR) next   # 2020: latent year, no observations
  if (row$Year < row$starts)   next   # before this site began monitoring
  
  val <- if      (is.na(row$count))   NA_real_        # ND = not surveyed
  else if (row$count == 0)    { n_zeros <- n_zeros + 1L; LOG_ZERO_OFFSET }
  else                          log(row$count)  # positive count
  
  if (row$age_class == "Breed") y_adult[s, t] <- val
  if (row$age_class == "Pup")   y_pup[s, t]   <- val
  if (row$age_class == "Molt")  y_molt[s, t]  <- val
}

cat(sprintf("Zero-count observations retained (log(0.5) offset): %d\n", n_zeros))

# Observation indicators (set before any Type H suppression)
y_adult_obs <- ifelse(!is.na(y_adult), 1L, 0L)
y_pup_obs   <- ifelse(!is.na(y_pup),   1L, 0L)
y_molt_obs  <- ifelse(!is.na(y_molt),  1L, 0L)

# Type H sites (is_breeder=0): DR, PB, Alcatraz.
# These sites DO conduct pup surveys but are not established breeding colonies.
# Convention:
#   - Surveyed years with zero pups (log(0.5)): KEEP in likelihood — the zero
#     constrains alpha_pup[s] near zero without inflating county N_pup.
#   - Surveyed years with positive pups: KEEP — incidental pup presence is real
#     data; alpha_pup[s] absorbs the signal without attributing production to
#     the county breeding pool.
#   - Unsurveyed years (NA): already excluded by y_pup_obs=0 from the fill loop.
# No forced suppression needed — the fill loop already handled ND→NA correctly.
# is_breeder=0 is used only by the Stan model to label site type for reporting;
# it does NOT suppress the pup likelihood here.

# Diagnostic: confirm Type H pup obs come from real surveys not imputation
type_h_sites <- which(regional_sites$is_breeder == 0)
cat(sprintf("Type H sites with pup survey data (is_breeder=0 but pup obs present):\n"))
for (s in type_h_sites) {
  n_pup_obs <- sum(y_pup_obs[s, ])
  n_pup_pos <- sum(y_pup[s, ] > 0 & y_pup_obs[s, ] == 1)
  n_pup_zero <- n_pup_obs - n_pup_pos
  cat(sprintf("  %s: %d pup surveys (%d positive, %d zero)\n",
              regional_sites$site_name[s], n_pup_obs, n_pup_pos, n_pup_zero))
}

# Replace NA with 0 in count matrices (Stan requires numeric, not NA)
y_adult[is.na(y_adult)] <- 0
y_pup[is.na(y_pup)]     <- 0
y_molt[is.na(y_molt)]   <- 0

# Diagnostics
cat("\n── Observation coverage by site ──────────────────────────────────────\n")
obs_summary <- regional_sites |>
  mutate(
    n_adult = rowSums(y_adult_obs),
    n_pup   = rowSums(y_pup_obs),
    n_molt  = rowSums(y_molt_obs)
  ) |>
  select(site_id, site_name, county, n_adult, n_pup, n_molt)
print(obs_summary)

# ── PART 6: SCENARIO PROJECTIONS SETUP ───────────────────────────────────────
N_scenarios    <- 3
scenario_names <- c("Status Quo (MOCI 0)", "Warm (MOCI +1)", "Cool (MOCI -1)")
# Row per scenario, col per projection year; MOCI deviation from mean
moci_proj      <- matrix(c(0, 1, -1), nrow = N_scenarios, ncol = T_PROJ_REGIONAL)

# ── PART 7: BUILD STAN DATA LIST ──────────────────────────────────────────────
county_id_vec  <- regional_sites$county_id
is_breeder_vec <- regional_sites$is_breeder
n_sites_county <- as.integer(table(factor(county_id_vec, levels = 1:C)))

# First time index with survey data per county — used to initialize N at the
# correct year rather than 2005 for late-starting counties (Mendocino = 2019).
county_t1 <- vapply(seq_len(C), function(ci) {
  min_start <- min(regional_sites$starts[regional_sites$county_id == ci])
  which(YEARS_REGIONAL == min_start)
}, integer(1))

cat(sprintf("Sites per county: %s\n",
            paste(paste0(c("Marin","Bay Mouth","South Bay","San Mateo","Sonoma","Mendocino"),
                         "=", n_sites_county), collapse = ", ")))
cat(sprintf("County start indices (t1): %s\n",
            paste(paste0(c("Marin","Bay Mouth","South Bay","San Mateo","Sonoma","Mendocino"),
                         "=", county_t1), collapse = ", ")))

regional_stan <- list(
  # Dimensions
  T           = T,
  S           = S,
  C           = C,
  T_proj      = T_PROJ_REGIONAL,
  N_scenarios = N_scenarios,
  
  # Count observations (log scale; 0 where unobserved)
  y_adult     = y_adult,
  y_pup       = y_pup,
  y_molt      = y_molt,
  
  # Observation indicators (1 = observed, 0 = not)
  y_adult_obs = y_adult_obs,
  y_pup_obs   = y_pup_obs,
  y_molt_obs  = y_molt_obs,
  
  # Site-level classifiers
  county_id      = county_id_vec,
  is_breeder     = is_breeder_vec,
  county_t1      = county_t1,        # first obs year index per county
  n_sites_county = n_sites_county,   # Fix 2: calibrates mu_log_alpha per county
  
  # County-level oceanographic type (0=coast, 1=bay mouth, 2=south bay)
  county_type = COUNTY_TYPE,
  
  # MOCI covariates (scaled; length T each)
  moci_jfm    = moci_scaled$moci_jfm,
  moci_amj    = moci_scaled$moci_amj,
  moci_ond    = moci_scaled$moci_ond,
  
  # Projections
  moci_proj   = moci_proj,
  
  # p_male_breed fixed at Marin IPM posterior mean
  p_male_fixed = 0.057
)

# ── PART 8: SAVE OUTPUT ───────────────────────────────────────────────────────
years_regional <- YEARS_REGIONAL

saveRDS(list(stan_data    = regional_stan,
             site_meta    = regional_sites,
             moci         = moci_scaled,
             counts_long  = regional_counts,
             years        = years_regional,
             scenario_names = scenario_names),
        "Output/regional_ipm_input_data.rds")

cat("\n── Final data dimensions ─────────────────────────────────────────────\n")
cat(sprintf("y_adult obs: %d / %d site-years\n", sum(y_adult_obs), S * T))
cat(sprintf("y_pup   obs: %d / %d site-years\n", sum(y_pup_obs),   S * T))
cat(sprintf("y_molt  obs: %d / %d site-years\n", sum(y_molt_obs),  S * T))
cat(sprintf("Latent year 2020: all obs = 0 (confirmed: %s)\n",
            all(y_adult_obs[, year_idx["2020"]] == 0)))

cat("\nObjects created: regional_sites, regional_counts, regional_moci,")
cat("\n                 regional_stan, years_regional")
cat("\nSaved: Output/regional_ipm_input_data.rds")
cat("\nNext: source(\"Code/16_regional_simulate.R\")\n")
