# data prep for 1997-2023 years adults and molt
                                            

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
Phoca <- read_excel("Data/1997_2023_Phocadata.xls")
#Phoca <- read_excel("Data/1997_2022_Phocadata.xls")

Phoca <- Phoca[-c(3:4, 8:10)]

## need some year and julian date fields to plot by year and day of year
library(lubridate)
## tell lubridate the format
Phoca$Date2 <- ymd(Phoca$Date) #separates out month,day and year
## add year column for faceting plots
Phoca$Year <- year(Phoca$Date2)
## get Julian Date (yday) for within year dates
Phoca$Julian <- yday(Phoca$Date2)

## Codde 2024-06-18
# April 15 - May 20 for Breeding. For molt, it's just after June 10. 

#rename adults during molting season "MOLTING"
Phoca$Age <- ifelse(Phoca$Julian > 155, "MOLTING", Phoca$Age) # ~June 5 (if June 10 get zero molt for TB in 1997)
Phoca$Season <- ifelse(Phoca$Julian <= 140 & Phoca$Julian >= 105, "PUPPING", "MOLTING")



## Looks like need to remove Dead adult and deadpup
library(plyr)
library(dplyr)

Phoca <- filter(Phoca, Age != "DEADPUP" & Age != "DEADADULT")
unique(Phoca$Age)


## and need to convert HPUP (typo in database!) to PUP
Phoca$Age[Phoca$Age == "HPUP"] <- "PUP"
Phoca$Yearf <- as.factor(Phoca$Year)

all_data <- Phoca

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
  theme_gray(base_size = 14) + 
  scale_x_continuous(breaks=seq(1980,2025,20))

top.plot + facet_grid(Age ~ Subsite)



## now within year plots to look at breeding and molting season peaks

#Phoca.Adult <- dplyr::filter(Phoca, Age == "ADULT")

# seasonal.plot <- ggplot(all_data, aes(Julian, Count, color = Yearf)) + 
#   geom_point() + 
#   geom_smooth(aes(colour = Yearf))
# seasonal.plot + facet_grid(Subsite ~ Age)

#Year 1 only has bolinas data so delete 
#all sites start in 1982...so use that
all_data.MARSS <- top1_all_data



# unique(all_data.MARSS$Subsite)

# remove extra columns
all_data.MARSS <- all_data.MARSS[,c(7, 2, 4,5)]

## log all the counts, # log = ln in R
all_data.MARSS$Count <- log(all_data.MARSS$Count) 

# check no duplicate rows
all_data.MARSS <- distinct(all_data.MARSS)

#make new age-subsite column
all_data.MARSS$Subsite_Age <- paste0(all_data.MARSS$Subsite, "_", all_data.MARSS$Age)

#remove extra columns
all_data.MARSS <- all_data.MARSS %>% filter(Year > 1996)

#bolinas missing 1997 Molt data, so us ethe 1999 (non-ENSO) data...seems pretty steady from 1999
addBL97molt<- data.frame(list(1997, "BL", "MOLTING", 5.7745515, "BL_MOLTING"))


addBL97molt <- data.frame(Year=1997,
                          Subsite="BL",
                          Age="MOLTING",
                          Count="5.7745515",
                          Subsite_Age="BL_MOLTING",
                          stringsAsFactors=FALSE)
all_data.MARSS <- rbind(addBL97molt,all_data.MARSS)

all_data.MARSS <- all_data.MARSS %>%
  arrange(Year,
          Subsite,
          Age)

# Remove PB and DR since not pupping sites and sometimes have zero pups
all_data.MARSS <- all_data.MARSS %>% 
  filter(Subsite != "DR" & Subsite != "PB" )


all_data.MARSS <- all_data.MARSS[,-c(2:3)]

#make sure counts are number
all_data.MARSS$Count <- as.numeric(all_data.MARSS$Count)

## now need to make the data wide for MARSS
all_data.MARSS.wide <- all_data.MARSS %>% 
  pivot_wider(names_from = Subsite_Age, values_from = Count)

## now use names in the basic MARSS code 
dat <- t(all_data.MARSS.wide)
## OR FOR PUPS
## dat <- t(top1.pup.breed.spread)

#Transpose since MARSS needs time ACROSS columns
years = dat[1,]     # remove years
n = nrow(dat)-1
dat = dat[2:nrow(dat),]
legendnames = (unlist(dimnames(dat)[1]))


