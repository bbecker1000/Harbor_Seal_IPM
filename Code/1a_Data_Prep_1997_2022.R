library("MARSS")
library(ggplot2)
library("readxl")
library(dplyr)
library(tidyverse)
library(sjPlot) #used at end for 1 plot theme

##%######################################################%##
#                                                          #
####                   data wrangling                   ####
#                                                          #
##%######################################################%##

#2023-12-21
#ALL SITES ALL YEARS !!



## get PRNS data
Phoca <- read_excel("Data/1997_2022_Phocadata.xls")
Phoca <- Phoca[-c(3:5, 8:10)]

## need some year and julian date fields to plot by year and day of year
library(lubridate)
## tell lubridate the format
Phoca$Date2 <- ymd(Phoca$Date) #separates out month,day and year
## add year column for faceting plots
Phoca$Year <- year(Phoca$Date2)
## get Julian Date (yday) for within year dates
Phoca$Julian <- yday(Phoca$Date2)

#rename adults during molting season "MOLTING"
Phoca$Age <- ifelse(Phoca$Julian > 135, "MOLTING", Phoca$Age)
Phoca$Season <- ifelse(Phoca$Julian <= 135, "PUPPING", "MOLTING")

## Looks like need to remove Dead adult and deadpup
library(plyr)
library(dplyr)

Phoca <- filter(Phoca, Age != "DEADPUP" & Age != "DEADADULT")
unique(Phoca$Age)

## and need to convert HPUP (typo in database!) to PUP
Phoca$Age[Phoca$Age == "HPUP"] <- "PUP"
Phoca$Yearf <- as.factor(Phoca$Year)

Phoca_1997_2022 <- Phoca

# Change PR to PRH
unique(all_data$Subsite)

Phoca_1997_2022$Subsite <- ifelse(Phoca_1997_2022$Subsite == "PR", "PRH", Phoca_1997_2022$Subsite)
unique(Phoca_1997_2022$Subsite)

#top1 master file (all data)
top1_Phoca_1997_2022 <- Phoca_1997_2022 %>% 
  group_by(Year, Subsite, Age) %>%
  slice_max(Count, n = 1) %>%
  filter(Age != "PUP" | Season != "MOLTING") #remove pups counted in early molting season
top1_Phoca_1997_2022

top.plot <- ggplot(top1_Phoca_1997_2022,
                   aes(Year, Count)) + 
  geom_point() + 
  geom_smooth() + 
  #geom_vline(xintercept = c(1988, 1994), lty = 3) +
  ylab("Max Count") +
  #scale_y_continuous(trans='log10') +
 # geom_vline(xintercept = 1994, lty = 3) +
  theme_gray(base_size = 14)
  
top.plot + facet_grid(Age ~ Subsite)



## now within year plots to look at breeding and molting season peaks

#Phoca.Adult <- dplyr::filter(Phoca, Age == "ADULT")

seasonal.plot <- ggplot(Phoca_1997_2022, aes(Julian, Count, group = Year)) + 
  geom_point(alpha = 0.05) + 
  geom_smooth(se = FALSE)
seasonal.plot + facet_grid(.~ Age)

#remove Point Bonita and Duxbury PUP data since no or few pups
# new .MARSS file
Phoca_1997_2022.MARSS <- subset(top1_Phoca_1997_2022, Subsite != "DR" & Subsite != "PB")

#check
unique(Phoca_1997_2022.MARSS$Subsite)
unique(Phoca_1997_2022.MARSS$Season)

# remove extra columns
Phoca_1997_2022.MARSS <- Phoca_1997_2022.MARSS[,c(6, 2, 3, 4, 8, 9)]

## log all the counts, # log = ln in R
Phoca_1997_2022.MARSS$Count <- log(Phoca_1997_2022.MARSS$Count) 

# check no duplicate rows
Phoca_1997_2022.MARSS <- distinct(Phoca_1997_2022.MARSS)

#make new age-subsite column
Phoca_1997_2022.MARSS$Subsite_Age <- paste0(Phoca_1997_2022.MARSS$Subsite, "_", Phoca_1997_2022.MARSS$Age)

#remove extra columns
# Phoca_1997_2022.MARSS <- Phoca_1997_2022.MARSS[,-c(2:3)]

# #use only data post 1996 ? ------------------------------------
# Phoca_1997_2022.MARSS <- Phoca_1997_2022.MARSS %>%
#                       filter(Year>1997) # if filter 1996, then no BL molt so order off.
# #--------------------------------------------------------------

## now need to make the data wide for MARSS
Phoca_1997_2022.MARSS.wide <- Phoca_1997_2022.MARSS %>% 
  pivot_wider(names_from = Subsite_Age, values_from = Count)




dat2 <- Phoca_1997_2022.MARSS.wide ## for use in SealPopStructure Script

## now use names in the basic MARSS code 
dat <- t(Phoca_1997_2022.MARSS.wide)
## OR FOR PUPS
## dat <- t(top1.pup.breed.spread)

#Transpose since MARSS needs time ACROSS columns
years = dat[1,]     # remove years
n = nrow(dat)-1
dat = dat[2:nrow(dat),]
legendnames = (unlist(dimnames(dat)[1]))



############# env data --------------------------------

# MEI file has several ocean indices and coyote data
MEI <- read_excel("Data/MEI.xlsx")
show(MEI)



### Use data only > 1996? -----------------
# MEI <- MEI %>%
#   filter(Year>1997)

####--------------------------------------
MEI$MOCI_JFM_NC <- as.numeric(MEI$MOCI_JFM_NC)
MEI$BEUTI_FEB_APR_37N_39N <- as.numeric(MEI$BEUTI_FEB_APR_37N_39N)

#View(MEI)
##cut to sealData time series

MEI <- MEI %>% filter(Year < 2023 & Year > 1981)

#ESeal at DP in 2003 and 2004.  Count as a coyote disturbance
MEI$Coyote_DP <- ifelse(MEI$Year == 2003, 1, MEI$Coyote_DP)
MEI$Coyote_DP <- ifelse(MEI$Year == 2004, 1, MEI$Coyote_DP)

#if NPGO
#NPGO <- MEI[,c(1,3)]
## make NPGO a vector
#NPGO <- NPGO[,2]
#if Coyote


# plot(MEI$NPGO_MAR, MEI$NPGO_JAN )
# plot(MEI$NPGO_MAR, MEI$MOCI_JFM_NC )
# plot(MEI$PDO_MAR, MEI$MOCI_JFM_NC )
# plot(MEI$Year, MEI$PDO_MAR)
# plot(MEI$MEIv2_DEC_JAN, MEI$NPGO_MAR)
# cor.test(MEI$PDO_MAR, as.numeric(MEI$MOCI_JFM_NC))
# # cor = 0.79 (0.6-0.9) P < 0.001
# # So let's use PDO_MAR as our env covariate.'
# plot(MEI$PDO_MAR, MEI$NPGO_MAR  )
# plot(MEI$MEIv2_DEC_JAN, MEI$PDO_MAR  )
# plot(MEI$MOCI_JFM_NC, MEI$BEUTI_FEB_APR_37N_39N  )
# plot(MEI$Year, MEI$MOCI_JFM_NC)

# MARSS recommends start with more covariates and go smaller.
# create covariates for 
# - coyote Y/N + PDO_MAR
# - coyote count + PDO_MAR
# - PDO_MAR
# - whichever Coyote measure does better
# - adding the anthropogenic positive and negatives?

#
Coyote_01 <- MEI[,c(7:11)]  
small_c_Coyote_01 <- as.matrix(t(Coyote_01))

Coyote_Rate <- MEI[,c(12:16)] 
small_c_Coyote_Rate <- as.matrix(t(Coyote_Rate))


#Matrix with coyote and PDO_MAR

Coyote_01_PDO_MAR <- MEI[,c(7:11, 6)]
## zeros for DR, PB, PRH
#Coyote_01_PDO_MAR$Coyote_DR <- 0
#Coyote_01_PDO_MAR$Coyote_PB <- 0
Coyote_01_PDO_MAR$Coyote_PRH <- 0

Coyote_01_PDO_MAR <- Coyote_01_PDO_MAR[,c(1:3, 7, 4:6)]
small_c_Coyote_01_PDO_MAR <- as.matrix(t(Coyote_01_PDO_MAR))

PDO_MAR <- MEI[,6]
small_c_PDO_MAR <- as.matrix(t(PDO_MAR))


#Matrix with coyote and MOCI

Coyote_01_MOCI <- MEI[,c(7:11, 5)]
small_c_Coyote_01_MOCI <- as.matrix(t(Coyote_01_MOCI))  ## coyote 01 generally much better than coyote rate 

Coyote_rate_MOCI <- MEI[,c(12:16, 5)]
small_c_Coyote_rate_MOCI <- as.matrix(t(Coyote_rate_MOCI))

Coyote_01_BEUTI <- MEI[,c(7:11, 17)]
small_c_Coyote_01_BEUTI <- as.matrix(t(Coyote_01_BEUTI)) 



MOCI <- MEI[,5]
small_c_MOCI <- as.matrix(t(MOCI))





