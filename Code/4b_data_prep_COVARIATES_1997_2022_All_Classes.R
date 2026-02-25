#4b covariate prep for 1997-2023
# all classes

# env data prep for all years
library(readxl)

############# env data --------------------------------

# MEI file has several ocean indices and coyote data
MEI <- read_excel("Data/MEI.xlsx", col_types = c("numeric", 
                                                 "numeric", "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric"))
show(MEI)

# 2026-02-06 looks ok

#add human disturbance data 
HumanDisturbance <- 
  read_excel("Data/HumanDisturbanceRate_1996To2023.xlsx")
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
  filter(Year > 1996) %>%
  select(-c(DR, PB))
#ready to join with MEI below

#assign the 2020 NA a zero
HumanDisturbance.wide[24,4] <- 0


### Use data only > 1996? -----------------
# MEI <- MEI %>%
#   filter(Year>1997)


#View(MEI)
##cut to sealData time series

MEI <- MEI %>% filter(Year < 2026 & Year > 1996)

MEI <- left_join(MEI, HumanDisturbance.wide, by = "Year")
#new HumanDist data is cols 29-34


#ESeal at DP in 2003 and 2004.  Count as a coyote disturbance
#MEI$Coyote_DP <- ifelse(MEI$Year == 2003, 1, MEI$Coyote_DP)
#MEI$Coyote_DP <- ifelse(MEI$Year == 2004, 1, MEI$Coyote_DP)

#Convert the COyote rate values to 3 year weighted means for memory effect.
# 20:30:50 for years t-2, t-1 and t




#
Coyote_01 <- MEI[,c(7:11)]  
small_c_Coyote_01 <- as.matrix(t(Coyote_01))

Coyote_Rate <- MEI[,c(12:16)] 
small_c_Coyote_Rate <- as.matrix(t(Coyote_Rate))

Coyote_3yr <- MEI[,c(20:24)] #added 2024-05-26 #scale for comparison)
small_c_Coyote_3yr <- as.matrix(t(Coyote_3yr))


#Matrix with coyote_3yr and UI and UI-lag

small_c_Coyote_3yr_UI_UI_lag <-tibble(MEI[,c(20:22, 18:19)]) #scale for the covariate plots!
small_c_Coyote_3yr_UI_UI_lag <- scale(small_c_Coyote_3yr_UI_UI_lag)
small_c_Coyote_3yr_UI_UI_lag <- as_tibble(small_c_Coyote_3yr_UI_UI_lag)

## zeros for DR, PB, PRH
small_c_Coyote_3yr_UI_UI_lag$Coyote_3yr_PRH <- 0
small_c_Coyote_3yr_UI_UI_lag$Coyote_3yr_TB <- 0
small_c_Coyote_3yr_UI_UI_lag$Coyote_3yr_TP <- 0


small_c_Coyote_3yr_UI_UI_lag <- small_c_Coyote_3yr_UI_UI_lag[,c(1:3, 6:8, 4:5)]
small_c_Coyote_3yr_UI_UI_lag <- as.matrix(t(small_c_Coyote_3yr_UI_UI_lag)) # 



####Matrix with coyote_3yr and MOCI and MOCI-lags

small_c_Coyote_3yr_MOCI_MOCI_lag <-tibble(MEI[,c(7:9, 5,26, 27, 28)]) #scale for the covariate plots!
small_c_Coyote_3yr_MOCI_MOCI_lag <- scale(small_c_Coyote_3yr_MOCI_MOCI_lag)
small_c_Coyote_3yr_MOCI_MOCI_lag <- as_tibble(small_c_Coyote_3yr_MOCI_MOCI_lag)

## zeros for DR, PB, PRH
small_c_Coyote_3yr_MOCI_MOCI_lag$Coyote_PRH_3yr <- 0
small_c_Coyote_3yr_MOCI_MOCI_lag$Coyote_TB_3yr <- 0
small_c_Coyote_3yr_MOCI_MOCI_lag$Coyote_TP_3yr <- 0


small_c_Coyote_3yr_MOCI_MOCI_lag <- small_c_Coyote_3yr_MOCI_MOCI_lag[,c(1:3, 8:10, 4:5, 7, 6)]
small_c_Coyote_3yr_MOCI_MOCI_lag <- as.matrix(t(small_c_Coyote_3yr_MOCI_MOCI_lag)) # 

####Matrix with coyote_3yr Human Disturbance and MOCI and MOCI-lags

small_c_Coyote_3yr_MOCI_MOCI_Dist_lag.table <-tibble(MEI[,c(7:9, 5,26, 27, 28, c(29:34))]) #scale for the covariate plots!
small_c_Coyote_3yr_MOCI_MOCI_Dist_lag <- scale(small_c_Coyote_3yr_MOCI_MOCI_Dist_lag.table)
small_c_Coyote_3yr_MOCI_MOCI_Dist_lag <- as_tibble(small_c_Coyote_3yr_MOCI_MOCI_Dist_lag)

## zeros for DR, PB, PRH
small_c_Coyote_3yr_MOCI_MOCI_Dist_lag$Coyote_PRH_3yr <- 0
small_c_Coyote_3yr_MOCI_MOCI_Dist_lag$Coyote_TB_3yr <- 0
small_c_Coyote_3yr_MOCI_MOCI_Dist_lag$Coyote_TP_3yr <- 0


small_c_Coyote_3yr_MOCI_MOCI_Dist_lag <- small_c_Coyote_3yr_MOCI_MOCI_Dist_lag[,c(1:3, 14:16, 4:5, 7, 6, 8:13)]
small_c_Coyote_3yr_MOCI_MOCI_Dist_lag <- as.matrix(t(small_c_Coyote_3yr_MOCI_MOCI_Dist_lag)) # 

#2024-08-21
#drop the eSeal_IMM from the model per Codde and Allen

#small_c_Coyote_3yr_MOCI_MOCI_Dist_lag <- small_c_Coyote_3yr_MOCI_MOCI_Dist_lag[-10,]


# Summary Table
summary(small_c_Coyote_3yr_MOCI_MOCI_Dist_lag.table)
library(vtable)
sumtable(small_c_Coyote_3yr_MOCI_MOCI_Dist_lag.table,
         out="return",
         group.long=TRUE)
st(small_c_Coyote_3yr_MOCI_MOCI_Dist_lag.table)

##MOCI plot
names(MEI)
MOCI.plot <- MEI[c(5, 26, 28, 27)]

MOCI.plot$year <- 1997:2023

p1 <- ggplot(MOCI.plot, aes(x=year, y = MOCI_JFM_NC)) +
  geom_line() +
  ylab("MOCI Pre \n pupping (JFM)") +
  xlab(NULL) + 
  ylim(-10, 10) +
  geom_hline(yintercept=0, color = "blue", linetype =2)+
  theme_minimal(base_size = 18)

p2 <- ggplot(MOCI.plot, aes(x=year, y = MOCI_LAG_OND_NC)) +
  geom_line() +
  ylab("MOCI Egg \n Implant (OND)") +
  xlab(NULL) +
  ylim(-10, 10) +
  geom_hline(yintercept=0, color = "blue", linetype =2)+
  theme_minimal(base_size = 18)

p3 <- ggplot(MOCI.plot, aes(x=year, y = MOCI_LAG_AMJ_NC)) +
  geom_line() +
  ylab("MOCI Prior \n Spring (AMJ)") +
  xlab(NULL) + 
  ylim(-10, 10) +
  geom_hline(yintercept=0, color = "blue", linetype =2)+
  theme_minimal(base_size = 18)

p4 <- ggplot(MOCI.plot, aes(x=year, y = eSeal_IMM)) +
  geom_line() +
  ylab("Immature \n elephant seals") +
  ylim(0,1500) +
  theme_minimal(base_size = 18)

cowplot::plot_grid(p1, p2, p3, p4, ncol = 1, labels="auto")

ggsave("Output/Plots/MOCI-Eseal.jpeg", width = 20, height = 30, units = "cm")

##################### ----------------------------------------------------

## 2026-02-06
## simpler covariate setup

library(readxl)

rownames(small_c_Coyote_3yr_MOCI_MOCI_Dist_lag)
#we need: 
# "Coyote_BL_3yr"   "Coyote_DE_3yr"   "Coyote_DP_3yr"   "Coyote_PRH_3yr"  "Coyote_TB_3yr"  
# [6] "Coyote_TP_3yr"   "MOCI_JFM_NC"     "MOCI_LAG_AMJ_NC" "MOCI_LAG_OND_NC" "eSeal_IMM"      
# [11] "BL"              "DE"              "DP"              "PRH"             "TB"             
# [16] "TP" 

#import  
coyote <- read_excel("Data/CoyoteSightings_2025.xlsx")
#for each year weight rate by 0.5, 0.3. 0.2 for t, t-1, and t-2

# Calculate weighted rates
coyote_rate <- coyote %>%
  # First calculate the raw rate for each year-site
  mutate(rate = `Number of days with coyote sightings` / `Total number of monitoring surveys`) %>%
  # Arrange by site and year to ensure proper ordering
  arrange(Site, Year) %>%
  # Group by site to calculate lagged values
  group_by(Site) %>%
  # Get lagged rates and check if they exist
  mutate(
    rate_t1 = lag(rate, 1),  # t-1 rate
    rate_t2 = lag(rate, 2),  # t-2 rate
    has_t1 = !is.na(rate_t1),
    has_t2 = !is.na(rate_t2)
  ) %>%
  # Calculate weighted rate based on available data
  mutate(
    weighted_rate = case_when(
      # If both t-1 and t-2 exist: 0.5 current + 0.3 t-1 + 0.2 t-2
      has_t1 & has_t2 ~ 0.5 * rate + 0.3 * rate_t1 + 0.2 * rate_t2,
      # If only t-1 exists: 0.7 current + 0.3 t-1
      has_t1 & !has_t2 ~ 0.7 * rate + 0.3 * rate_t1,
      # If neither exist: 1.0 current
      TRUE ~ 1.0 * rate
    )
  ) %>%
  # Select final columns
  select(Year, Site, rate, weighted_rate) %>%
  ungroup()  %>%
  # Sort by Year, then Site
  arrange(Year, Site)

# View results
print(coyote_rate)

# If you want just the final weighted rate column
coyote_rate <- coyote_rate %>%
  select(Year, Site, weighted_rate)

print(coyote_rate)

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
  filter(Year > 1996) %>%
  select(-c(DR, PB))

View(HumanDisturbance.wide)

#assign the 2020 NA a zero
HumanDisturbance.wide[24,4] <- 0
HumanDisturbance.wide

## now MOCI data--------
MOCI <- read_csv("Data/CaliforniaMOCI.csv")
# Data wrangling
MOCI.dat <- MOCI %>%
  # Calculate mean of North and Central California
  mutate(mean_value = (`North California (38-42N)` + `Central California (34.5-38N)`) / 2) %>%
  # Select only needed columns
  select(Year, Season, mean_value) %>%
  # Create a lead year for OND values
  mutate(Year = ifelse(Season == "OND", Year + 1, Year)) %>%
  # Filter to keep only JFM, AMJ, and OND
  filter(Season %in% c("JFM", "AMJ", "OND")) %>%
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

#now the eSeal data------
# Instructions from 2024:
# eSeal data-> SUM the IMMATURE MOLT and the WEANED PUP count as for all PR sites (total pop) covariate.

eSeal <- 
  read_excel("Data/Eseal_1981-2025_BySubsite.xlsx")

library(dplyr)
library(lubridate)
library(tidyverse)


unique(eSeal$MatureCode)

# Extract year and calculate annual max by SubSiteName, then sum
eSeal_max_imm <- eSeal %>%
  mutate(Year = year(StartDate)) %>%
  filter(MatureCode %in% c("IMM", "WNR")) %>%
  group_by(Year, SubSiteName, MatureCode) %>%
  summarise(MaxCount = max(Count, na.rm = TRUE), .groups = "drop") %>%
  group_by(Year) %>%
  summarise(TotalImmMaxCount = sum(MaxCount))

eSeal_max_imm

ggplot(eSeal_max_imm, aes(Year, TotalImmMaxCount)) +
  geom_line() +
  theme_gray(base_size = 16)

#####

## Now put all the covariates together










