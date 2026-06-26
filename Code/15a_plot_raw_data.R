# ============================================================================
# 15a_plot_raw_data.R
# ----------------------------------------------------------------------------
# RAW SURVEY DATA PLOTS — ALL SITES BY COUNTY
#
# Reads both raw data files and plots observed maximum annual counts for
# each site, faceted by county, for all three age classes (Breed/Adult,
# Pup, Molt). Shows the actual count data that drives the IPM — useful for
# checking data quality, identifying outliers, and understanding coverage
# gaps before modelling.
#
# Output: Output/Plots/Raw_data_*.jpeg (one file per county, one combined)
#
# Prereq: none (reads data files directly; does not require script 15)
# ============================================================================

library(tidyverse)
library(readxl)
library(lubridate)
library(patchwork)

dir.create("Output",       showWarnings = FALSE)
dir.create("Output/Plots", showWarnings = FALSE)

# ── Constants ─────────────────────────────────────────────────────────────────
YEARS_PLOT   <- 2005:2025
MARIN_SITES  <- c("BL","DE","DP","PRH","TB","TP","DR","PB")

# Site display names for Marin codes
MARIN_NAMES <- c(
  BL  = "Bolinas Lagoon",
  DE  = "Drakes Estero",
  DP  = "Double Point",
  PRH = "Pt Reyes Headlands",
  TB  = "Tomales Bay",
  TP  = "Tomales Point",
  DR  = "Duxbury Reef",
  PB  = "Point Bonita"
)

# County colour palette (matches IPM scripts)
COUNTY_COLS <- c(
  "Marin"      = "#08519C",
  "Bay Mouth"  = "#F16913",
  "South Bay"  = "#D94801",
  "San Mateo"  = "#6BAED6",
  "Sonoma"     = "#41AB5D",
  "Mendocino"  = "#74C476"
)

AGE_COLS <- c(Breed = "#2166AC", Pup = "#1B7837", Molt = "#8C510A")

theme_raw <- function(base_size = 12) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(colour = "grey70", fill = NA),
      strip.text         = element_text(face = "bold", size = rel(0.9)),
      strip.background   = element_rect(fill = "grey94", colour = "grey80"),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
      legend.position    = "bottom",
      plot.title         = element_text(face = "bold", size = rel(1.0)),
      plot.subtitle      = element_text(colour = "grey40", size = rel(0.85)),
      plot.margin        = margin(8, 12, 8, 8)
    )
}

# ── PART 1: READ MARIN RAW DATA ───────────────────────────────────────────────
cat("Reading Marin data...\n")
marin_raw <- read_excel("Data/1996_2025_Phocadata.xls")
marin_raw <- marin_raw[, c("Date","Subsite","Age","Count")]

marin_processed <- marin_raw |>
  mutate(
    Date2   = ymd(Date),
    Year    = year(Date2),
    Julian  = yday(Date2),
    Age     = ifelse(Julian > 155, "MOLTING", Age),
    Season  = ifelse(Julian <= 140 & Julian >= 105, "PUPPING", "MOLTING"),
    Subsite = ifelse(Subsite == "PR", "PRH", Subsite),
    Count   = suppressWarnings(as.numeric(Count))
  ) |>
  filter(
    Age     %in% c("ADULT","MOLTING","PUP"),
    Subsite %in% MARIN_SITES,
    Year    %in% YEARS_PLOT
  ) |>
  filter(!(Age == "PUP" & Season == "MOLTING")) |>
  filter(!is.na(Count), Count >= 0)

# Annual maximum per site-year-class
marin_max <- marin_processed |>
  group_by(Year, Subsite, Age) |>
  slice_max(Count, n = 1, with_ties = FALSE) |>
  ungroup() |>
  transmute(
    Year,
    site_name  = MARIN_NAMES[Subsite],
    county     = "Marin",
    age_class  = recode(Age, "ADULT" = "Breed", "MOLTING" = "Molt", "PUP" = "Pup"),
    count      = Count
  )

cat(sprintf("Marin: %d site-year-class records\n", nrow(marin_max)))

# ── PART 2: READ REGIONAL DATA ────────────────────────────────────────────────
cat("Reading regional data...\n")
regional_raw <- read_excel("Data/RegionalPhoca_2005To2025_ForBecker.xlsx")

regional_cleaned <- regional_raw |>
  mutate(count = suppressWarnings(as.numeric(Count))) |>
  filter(!is.na(count) | Count == "0") |>  # keep explicit zeros
  mutate(count = suppressWarnings(as.numeric(Count))) |>
  select(Year, site_name = SiteName, County, age_class = AgeClass, count) |>
  mutate(
    county = case_when(
      site_name %in% c("Castro Rocks","Yerba Buena Island","Alcatraz") ~ "Bay Mouth",
      site_name %in% c("Mowry Slough","Newark Slough")                 ~ "South Bay",
      County == "San Mateo"                                            ~ "San Mateo",
      County == "Sonoma"                                               ~ "Sonoma",
      County == "Mendocino"                                            ~ "Mendocino",
      TRUE                                                             ~ County
    )
  ) |>
  select(-County) |>
  filter(Year %in% YEARS_PLOT)

cat(sprintf("Regional: %d site-year-class records\n", nrow(regional_cleaned)))

# ── PART 3: COMBINE ───────────────────────────────────────────────────────────
all_counts <- bind_rows(marin_max, regional_cleaned) |>
  mutate(
    county    = factor(county, levels = c("Marin","Bay Mouth","South Bay",
                                          "San Mateo","Sonoma","Mendocino")),
    age_class = factor(age_class, levels = c("Breed","Pup","Molt"))
  ) |>
  arrange(county, site_name, Year, age_class)

cat(sprintf("Combined: %d records, %d sites, %d counties\n",
            nrow(all_counts),
            n_distinct(all_counts$site_name),
            n_distinct(all_counts$county)))

# ── PART 4: PLOT FUNCTIONS ────────────────────────────────────────────────────

plot_county_raw <- function(county_name, data = all_counts,
                            save = TRUE, prefix = "Raw_data") {
  df <- data |> filter(county == county_name)
  if (nrow(df) == 0) { cat("No data for", county_name, "\n"); return(NULL) }
  
  col <- COUNTY_COLS[county_name]
  sites_ordered <- sort(unique(df$site_name))
  df <- df |> mutate(site_name = factor(site_name, levels = sites_ordered))
  
  p <- ggplot(df, aes(x = Year, y = count, colour = age_class, shape = age_class)) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    geom_point(size = 2.2, na.rm = TRUE) +
    scale_colour_manual(values = AGE_COLS, name = "Survey type") +
    scale_shape_manual(values  = c(Breed = 16, Pup = 17, Molt = 15), name = "Survey type") +
    scale_x_continuous(breaks = seq(2005, 2025, 5)) +
    scale_y_continuous(labels = scales::comma, limits = c(0, NA), expand = c(0.02, 0)) +
    facet_wrap(~ site_name, ncol = 2, scales = "free_y") +
    labs(
      x        = "Year",
      y        = "Maximum annual count",
      title    = sprintf("%s County — Raw Survey Counts by Site", county_name),
      subtitle = "Annual maximum count per survey type; gaps = no survey (ND); zeros shown explicitly"
    ) +
    theme_raw()
  
  if (save) {
    n_sites <- n_distinct(df$site_name)
    h <- max(12, n_sites * 5)
    ggsave(paste0("Output/Plots/", prefix, "_", gsub(" ", "_", county_name), ".jpeg"),
           p, width = 26, height = h, units = "cm", dpi = 200)
    cat(sprintf("  Saved: %s\n", county_name))
  }
  p
}

# ── PART 5: COUNTY-LEVEL TOTALS ───────────────────────────────────────────────
plot_county_totals <- function(data = all_counts, save = TRUE, prefix = "Raw_data") {
  totals <- data |>
    group_by(county, Year, age_class) |>
    summarise(total = sum(count, na.rm = TRUE), .groups = "drop")
  
  p <- ggplot(totals, aes(x = Year, y = total, colour = county, linetype = age_class)) +
    geom_line(linewidth = 1.0) +
    geom_point(aes(shape = age_class), size = 2) +
    scale_colour_manual(values = COUNTY_COLS, name = "County") +
    scale_linetype_manual(values = c(Breed = "solid", Pup = "dashed", Molt = "dotted"),
                          name = "Survey type") +
    scale_shape_manual(values = c(Breed = 16, Pup = 17, Molt = 15), name = "Survey type") +
    scale_x_continuous(breaks = seq(2005, 2025, 5)) +
    scale_y_continuous(labels = scales::comma) +
    facet_wrap(~ age_class, ncol = 1, scales = "free_y") +
    labs(x = "Year", y = "Sum of maximum counts (all sites)",
         title = "County-Level Survey Count Totals",
         subtitle = "Sum across all sites per county per year; reflects raw count totals, not unique animals") +
    theme_raw()
  
  if (save)
    ggsave(paste0("Output/Plots/", prefix, "_county_totals.jpeg"),
           p, width = 28, height = 30, units = "cm", dpi = 200)
  p
}

# ── PART 6: ALL-SITE HEATMAP ─────────────────────────────────────────────────
plot_coverage_heatmap <- function(data = all_counts, save = TRUE, prefix = "Raw_data") {
  # Show which site-years have data for each survey type
  cov <- data |>
    group_by(county, site_name, Year) |>
    summarise(
      has_breed = any(age_class == "Breed" & !is.na(count)),
      has_pup   = any(age_class == "Pup"   & !is.na(count)),
      has_molt  = any(age_class == "Molt"  & !is.na(count)),
      n_types   = has_breed + has_pup + has_molt,
      .groups = "drop"
    ) |>
    mutate(
      site_county = paste0(site_name, " (", as.character(county), ")"),
      site_county = factor(site_county,
                           levels = rev(unique(site_county[order(county, site_name)])))
    )
  
  p <- ggplot(cov, aes(x = Year, y = site_county, fill = factor(n_types))) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_manual(
      values = c("0" = "grey90", "1" = "#FEE0D2", "2" = "#FC9272", "3" = "#DE2D26"),
      name = "Survey types\nconducted",
      labels = c("0" = "None", "1" = "1 type", "2" = "2 types", "3" = "All 3"),
      drop = FALSE
    ) +
    scale_x_continuous(breaks = seq(2005, 2025, 2)) +
    labs(x = "Year", y = NULL,
         title = "Survey Coverage Heatmap — All Sites",
         subtitle = "Red = all 3 survey types completed; grey = no survey (ND); 2020 = COVID gap") +
    theme_raw() +
    theme(axis.text.y = element_text(size = 7.5),
          panel.grid  = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  
  if (save)
    ggsave(paste0("Output/Plots/", prefix, "_coverage_heatmap.jpeg"),
           p, width = 30, height = 32, units = "cm", dpi = 200)
  p
}

# ── PART 7: RUN ALL PLOTS ─────────────────────────────────────────────────────
cat("\nGenerating raw data plots...\n")

county_plots <- map(
  levels(all_counts$county),
  ~ { cat(sprintf("  Plotting %s...\n", .x)); plot_county_raw(.x) }
)
names(county_plots) <- levels(all_counts$county)

cat("  Plotting county totals...\n")
p_totals <- plot_county_totals()

cat("  Plotting coverage heatmap...\n")
p_heatmap <- plot_coverage_heatmap()

cat("\nRaw data plots complete.\n")
cat("Output files:\n")
cat("  Output/Plots/Raw_data_<County>.jpeg  — one per county\n")
cat("  Output/Plots/Raw_data_county_totals.jpeg\n")
cat("  Output/Plots/Raw_data_coverage_heatmap.jpeg\n")

# Return for interactive inspection
invisible(list(
  county_plots = county_plots,
  totals       = p_totals,
  heatmap      = p_heatmap,
  data         = all_counts
))

