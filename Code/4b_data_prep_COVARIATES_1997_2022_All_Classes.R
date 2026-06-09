
## 2026-02-06
## simpler covariate setup

library(readxl)
library(tidyverse)
library(dplyr)

#rownames(small_c_Coyote_3yr_MOCI_MOCI_Dist_lag)
#we need: 
# "Coyote_BL_3yr"   "Coyote_DE_3yr"   "Coyote_DP_3yr"   "Coyote_PRH_3yr"  "Coyote_TB_3yr"  
# [6] "Coyote_TP_3yr"   "MOCI_JFM_NC"     "MOCI_LAG_AMJ_NC" "MOCI_LAG_OND_NC" "eSeal_IMM"      
# [11] "BL"              "DE"              "DP"              "PRH"             "TB"             
# [16] "TP" 

coyote <- read_excel("Data/CoyoteSightings_2025.xlsx")

coyote_rate <- coyote %>%
  mutate(rate = `Number of days with coyote sightings` / `Total number of monitoring surveys`) %>%
  arrange(Site, Year) %>%
  group_by(Site) %>%
  mutate(
    rate_t1 = lag(rate, 1),
    rate_t2 = lag(rate, 2),
    has_t1 = !is.na(rate_t1),
    has_t2 = !is.na(rate_t2)
  ) %>%
  mutate(
    weighted_rate = case_when(
      has_t1 & has_t2 ~ 0.5 * rate + 0.3 * rate_t1 + 0.2 * rate_t2,
      has_t1 & !has_t2 ~ 0.7 * rate + 0.3 * rate_t1,
      TRUE ~ 1.0 * rate
    )
  ) %>%
  # Replace NaN with mean of previous 2 years
  mutate(
    weighted_rate = case_when(
      is.nan(weighted_rate) & has_t1 & has_t2 ~ (rate_t1 + rate_t2) / 2,
      is.nan(weighted_rate) & has_t1 ~ rate_t1,
      is.nan(weighted_rate) ~ NA_real_,
      TRUE ~ weighted_rate
    )
  ) %>%
  select(Year, Site, rate, weighted_rate) %>%
  ungroup() %>%
  arrange(Year, Site)

# View results
print(coyote_rate)

# If you want just the final weighted rate column
coyote_rate <- coyote_rate %>%
  select(Year, Site, weighted_rate)

coyote_rate$weighted_rate[is.nan(coyote_rate$weighted_rate)] <- NA

print(coyote_rate)

coyote_wide <- coyote_rate %>%
  pivot_wider(names_from = Site, values_from = weighted_rate)

# Create rows for 1996-1999 with zeros
new_rows <- tibble(
  Year = 1996:1999,
  BL = 0, DE = 0, DP = 0, DR = 0, PB = 0, PRH = 0, TB = 0, TP = 0
)

# Bind and arrange by year
coyote_wide <- bind_rows(new_rows, coyote_wide) %>%
  arrange(Year)

#fix names
coyote_wide <- coyote_wide %>%
  rename_with(~ paste0("Coyote_", .), -Year)

## now human disturbance data---------
HumanDisturbance <- 
  read_excel("Data/HumanDisturbanceRate_1996To2025.xlsx")
HumanDisturbance <- 
  HumanDisturbance[,-5] #delete and give better name
HumanDisturbance$DistRate <- 
  HumanDisturbance$SumOfDisturbanceCount / HumanDisturbance$NSurveys
#remove raw data, keep rate
HumanDisturbance <- HumanDisturbance[,-c(3:4)]

#pivot wide
HumanDisturbance.wide <- 
  HumanDisturbance %>% pivot_wider(names_from = SiteCode,
                                   values_from = DistRate)
#remove DR and PB and 1996
HumanDisturbance.wide <-
  HumanDisturbance.wide %>%
  dplyr::filter(Year > 1996) %>%
  select(-c(DR, PB))

(HumanDisturbance.wide)

#assign the 2020 NA a zero
HumanDisturbance.wide[24,4] <- 0
HumanDisturbance.wide

#fix names
HumanDisturbance.wide <- HumanDisturbance.wide %>%
  rename_with(~ paste0("Dist_", .), -Year)


## now MOCI data--------
MOCI <- read_csv("Data/CaliforniaMOCI.csv")
# Data wrangling
MOCI.dat <- MOCI %>%
  # Calculate mean of North and Central California
  mutate(mean_value = (`North California (38-42N)` + `Central California (34.5-38N)`) / 2) %>%
  # Select only needed columns
  dplyr::select(Year, Season, mean_value) %>%
  # Create a lead year for OND values
  mutate(Year = ifelse(Season == "OND", Year + 1, Year)) %>%
  # Filter to keep only JFM, AMJ, and OND
  dplyr::filter(Season %in% c("JFM", "AMJ", "OND")) %>%
  # Pivot wider to create columns for each season
  pivot_wider(
    names_from = Season,
    values_from = mean_value
  ) %>%
  # Arrange by Year
  arrange(Year) %>%
  # Reorder columns
  select(Year, JFM, AMJ, OND)

# View the result
(MOCI.dat)

#fix names
MOCI.dat <- MOCI.dat %>%
  rename_with(~ paste0("MOCI_", .), -Year)


#now the eSeal data------
# Instructions from 2024:
# eSeal data-> SUM the IMMATURE MOLT and the WEANED PUP count as for all PR sites (total pop) covariate.

eSeal <- 
  read_excel("Data/Eseal_1981-2025_BySubsite.xlsx")

unique(eSeal$MatureCode)

# 2026-04-03
# !!!!! UPDATE THIS TO EXTRACT APRIL-MAY-JUNE MAX MOLT COUNTS SUM OF ALL AGE CLASSES
# !!!! 2026-04-08  no AMJ data in the codde data

# Extract year and calculate annual max by SubSiteName, then sum
library(lubridate)


eSeal_max_imm <- eSeal %>%
  dplyr::mutate(Year = lubridate::year(StartDate)) %>%
  # dplyr::filter(MatureCode %in% c("IMM", "WNR")) %>%  # use all age classes per SA 2026-04-27
  dplyr::group_by(Year, SubSiteName, MatureCode) %>%
  dplyr::summarise(MaxCount = max(Count, na.rm = TRUE), .groups = "drop") %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(eSeal_Sum_Imm_MaxCount = sum(MaxCount))

eSeal_max_imm

ggplot(eSeal_max_imm, aes(Year, eSeal_Sum_Imm_MaxCount)) +
  geom_line() +
  theme_gray(base_size = 16)

#####

## Now put all the covariates together

coyote_wide
HumanDisturbance.wide
MOCI.dat
eSeal_max_imm

covariates <- MOCI.dat %>%
  left_join(HumanDisturbance.wide, by = "Year") %>%
  left_join(coyote_wide, by = "Year") %>%
  left_join(eSeal_max_imm, by = "Year") %>%
  dplyr::filter(Year>1996) %>% ## seal data from 1997 - present
  dplyr::filter(Year<2026)  #some leftover 2026 data



#transform
cov_t <- t(covariates)
#remove 2 sites coyote data 
cov_t <- cov_t[!rownames(cov_t) %in% c("Coyote_PB", "Coyote_DR"), ]

### now bring in the seal data from 4a
dat

# Remove Year row from covariates
cov_t_marss <- cov_t[-1, ]  # Remove first row (Year)

# Standardize covariates (z-score: mean=0, sd=1)
cov_t_scaled <- t(scale(t(cov_t_marss)))

# Check
rowMeans(cov_t_scaled)  # Should be ~0
apply(cov_t_scaled, 1, sd)  # Should be ~1


cov_t_marss[is.na(cov_t_marss)] <- 0

# Scale only rows with variance
cov_t_scaled <- cov_t_marss
rows_to_scale <- apply(cov_t_marss, 1, sd) > 0
cov_t_scaled[rows_to_scale, ] <- t(scale(t(cov_t_marss[rows_to_scale, ])))


# dat appears ready - just confirm it's a matrix
dat_marss <- as.matrix(dat)

# Verify dimensions match (same number of columns = years)
ncol(cov_t_marss)  # Should match
ncol(dat_marss)    # Should match

dim(dat)          # Should be 18 × 29
dim(cov_t_marss)  # Should be 16 × 29

#remove NAs and replace with 0
cov_t_marss[is.na(cov_t_marss)] <- 0








