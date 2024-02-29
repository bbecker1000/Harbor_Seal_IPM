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

Phoca <- Phoca[-c(3:4, 8:10)]

#older data
hseal_max_1975_1999 <- read_excel("Data/1975-1999_max_counts.xlsx")

hseal_max_1975_1999$Age <- ifelse(hseal_max_1975_1999$Season =="MOLTING", "MOLTING", hseal_max_1975_1999$Age)

hseal_max_1975_1999_long <- hseal_max_1975_1999 %>% 
  pivot_longer(!c(Year, Age, Season), names_to = "Subsite", values_to = "Count")
hseal_max_1975_1999_long

# remove pre-1996 data since in the modern data file
hseal_max_1975_1995_long <- hseal_max_1975_1999_long %>%
  filter(Year < 1996)



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


#now put files together
Phoca
hseal_max_1975_1995_long

all_data <- bind_rows(hseal_max_1975_1995_long, Phoca)



all_data$Yearf <- as.factor(all_data$Year)
all_data

# Change PR to PRH
unique(all_data$Subsite)

all_data$Subsite <- ifelse(all_data$Subsite == "PR", "PRH", all_data$Subsite)
unique(all_data$Subsite)

#top1 master file (all data)
top1_all_data <- all_data %>% 
  group_by(Year, Subsite, Age) %>%
  slice_max(Count, n = 1) %>%
  filter(Age != "PUP" | Season != "MOLTING") #remove pups counted in early molting season
top1_all_data

top.plot <- ggplot(top1_all_data,
                   aes(Year, Count)) + 
  geom_point() + 
  geom_smooth() + 
  geom_vline(xintercept = c(1988, 1994), lty = 3) +
  ylab("Max Count") +
  #scale_y_continuous(trans='log10') +
  geom_vline(xintercept = 1994, lty = 3) +
  theme_gray(base_size = 14)
  
top.plot + facet_grid(Age ~ Subsite)



## now within year plots to look at breeding and molting season peaks

#Phoca.Adult <- dplyr::filter(Phoca, Age == "ADULT")

seasonal.plot <- ggplot(all_data, aes(Julian, Count, color = Yearf)) + 
  geom_point() + 
  geom_smooth(aes(colour = Yearf))
seasonal.plot + facet_grid(Subsite ~ Age)

#remove Point Bonita and Duxbury PUP data since no or few pups
all_data.MARSS <- subset(top1_all_data, Age != "PUP") # | Subsite != "DR" & Subsite != "PB" & Subsite != "PRH")

#remove Point Bonita and Duxbury PUP data since no or few pups and shorter time series
all_data.MARSS <- subset(all_data.MARSS, Subsite != "DR" & Subsite != "PB")

#Year 1 only has bolinas data so delete 
#all sites start in 1982...so use that
all_data.MARSS <- subset(all_data.MARSS, Year > 1981)



# unique(all_data.MARSS$Subsite)

# remove extra columns
all_data.MARSS <- all_data.MARSS[,c(1, 2, 4,5)]

## log all the counts, # log = ln in R
all_data.MARSS$Count <- log(all_data.MARSS$Count) 

# check no duplicate rows
all_data.MARSS <- distinct(all_data.MARSS)

#make new age-subsite column
all_data.MARSS$Subsite_Age <- paste0(all_data.MARSS$Subsite, "_", all_data.MARSS$Age)

#remove extra columns
all_data.MARSS <- all_data.MARSS[,-c(2:3)]

# #use only data post 1997 ? ------------------------------------
# all_data.MARSS <- all_data.MARSS %>%
#                       filter(Year>1997) # if filter 1996, then no BL molt so order off.
# #--------------------------------------------------------------

## now need to make the data wide for MARSS
all_data.MARSS.wide <- all_data.MARSS %>% 
  pivot_wider(names_from = Subsite_Age, values_from = Count)




dat2 <- all_data.MARSS.wide ## for use in SealPopStructure Script

## now use names in the basic MARSS code 
dat <- t(all_data.MARSS.wide)
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





