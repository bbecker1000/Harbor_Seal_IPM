#test comment
library("readxl")
library(dplyr)
library(tidyverse)
hseals <- read_excel("Data/1997_2019_Phocadata.xls")
# drop extra cols, explained below
hseals <- hseals[-c(3:4, 8:10)]

# ended up with some extra cols:
# "...10" had N/A's and "* Not sure if TRB was skipped or just no seals" for 	1999-04-12; TP
# "...8" has the names of the people conducting the surveys (?)
# "...9" has notes for the data, copied below

# Notes on the data:
# "\"Pup\" age classification = Live pups only"                                                                                                                                       
# "1. For 1996-99 data, I did not include surveys that had comments about poor visibility or difficulty in identifying pups"                                                          
# "2. Tomales Bay and Tomales Point surveys in 1998 often had only one count"                                                                                                         
# "3. Bolinas, Drakes Estero, and Double Point start in 1996; Pt Reyes Headlands, Tomales Bay, and Tomales Point start in 1997; Duxbury data starts in 1999; Pt Bonita starts in 2000"
# "4. Pt Reyes data includes Sea Lion Overlook"                                                                                                                                       
# "5. TB and TP have only one molting count in 2001"                                                                                                                                  
# "6. No molting counts for any sites in 1996 except one at DE."                                                                                                                      
# "*Another survey on this day with higher counts, but at a higher tide height" for "1999-05-04; PR

# unique(hseals$Subsite)
# abbrev.:    "DP"         "DE"         "BL"      "TB"          "TP"          "PR"       "DR"       "PB"
# name:  Double Point  Drakes Estero  Bolinas  Tomales Bay  Tomales Point  Point Reyes  Duxbury  Point Bonita

get_year <- function(date) {
  if (date$mon > 10) { 
    return(date$year+1901)
  } else {
    return(date$year+1900)
  }
}

# function that does get_year on every item in a list
get_year_all <- function(lst) {
  vec <- numeric()
  for (i in 1:length(lst)) {
    vec <- c(vec, get_year(lst[[i]]))
  }
  return(vec)
}

hseals_date <- strptime(hseals$Date, format="%Y-%m-%d")
hseals$Year <- get_year_all(hseals_date)
hseals$Month_day <- format(hseals$Date, format="%m-%d")
hseals$Month <- format(hseals$Date, format="%m")

# two types of observations: during breeding season (March, April, and May) or during molting season (June July)
# add cols indicating which one for the observation

hseals$Breeding <- hseals$Month == "03" | hseals$Month == "04" | hseals$Month == "05"
hseals$Molting <- hseals$Month == "06" | hseals$Month == "07"

# only want to look at breeding between April 15 and May 15
peakBreed <- format(seq(as.Date("1981-04-15"), as.Date("1981-05-15"), by="days"), format="%m-%d")
hseals$Peak_breeding <- hseals $Month_day %in% peakBreed

# need to get subsite time in better format
hseals_time <- strptime(hseals$SubsiteTime, format="%Y-%m-%d %H:%M:%S")
hseals$Time <- format(hseals_time, "%H:%M")

# drop cols we don't need anymore
hseals <- hseals[-c(3, 7:9)]
# reorder cols
hseals <- hseals[, c(1, 8, 2, 3, 4, 6, 7, 5)]

# unique(hseals$Age)
# "ADULT"     "PUP"       "DEADPUP"   "DEADADULT" "HPUP" 

# hpup <- subset(hseals, Age == "HPUP")
# PRpup1819 <- subset(hseals, Age == "PUP" & Subsite == "PR" ) #& (Year == 2018 | Year == 2019))
# 12 obs in 2018 and 2019 labeled HPUP, all for site == PR; PR has no pup obs for these years so relabel to PUP
hseals$Age[hseals$Age == "HPUP"] <- "PUP"

# don't need to worry about deadpup and deadadult for right now
hseals <- subset(hseals, Age == "PUP" | Age == "ADULT")

# function to get the max count for each year.
# args are location (string) and table 
get_max <- function(tbl, location) {
  tbl <- subset(tbl, Subsite == location)
  m <- tbl %>% group_by(Year) %>% slice(which.max(Count))
  return(m)
}

#######################
### BREEDING SEASON ###
#######################

hseals_breed <- subset(hseals, Peak_breeding == TRUE)

# unique(hseals_breed$Age) gives "ADULT" and "PUP"  

adult_b <- subset(hseals_breed, Age == "ADULT")

DP_adult_b <- get_max(adult_b, "DP") #not missing any years
DE_adult_b <- get_max(adult_b, "DE") #not missing any years
BL_adult_b <- get_max(adult_b, "BL") #not missing any years 
TB_adult_b <- get_max(adult_b, "TB") #not missing any years 
TP_adult_b <- get_max(adult_b, "TP") #not missing any years 
PR_adult_b <- get_max(adult_b, "PR") #not missing any years 
DR_adult_b <- get_max(adult_b, "DR") #not missing any years 
PB_adult_b <- get_max(adult_b, "PB") #not missing any years 

# put it all back into one table
adult_breeding <- data.frame(Date=double(), Time=character(), Subsite=character(), Age=character(),
                             Count=numeric(), Molting=logical(), Peak_breeding=logical(), Year=numeric())

adult_breeding <- rbind(adult_breeding, DP_adult_b)
adult_breeding <- rbind(adult_breeding, DE_adult_b)
adult_breeding <- rbind(adult_breeding, BL_adult_b)
adult_breeding <- rbind(adult_breeding, TB_adult_b)
adult_breeding <- rbind(adult_breeding, TP_adult_b)
adult_breeding <- rbind(adult_breeding, PR_adult_b)
adult_breeding <- rbind(adult_breeding, DR_adult_b)
adult_breeding <- rbind(adult_breeding, PB_adult_b)

pup_b <- subset(hseals_breed, Age == "PUP")

DP_pup_b <- get_max(pup_b, "DP") #not missing any years
DE_pup_b <- get_max(pup_b, "DE") #not missing any years
BL_pup_b <- get_max(pup_b, "BL") #not missing any years
TB_pup_b <- get_max(pup_b, "TB") #not missing any years
TP_pup_b <- get_max(pup_b, "TP") #not missing any years
PR_pup_b <- get_max(pup_b, "PR") #not missing any years
DR_pup_b <- get_max(pup_b, "DR") #not missing any years
PB_pup_b <- get_max(pup_b, "PB") #not missing any years

# put it all back in one table for plotting
pup_breeding <- data.frame(Date=double(), Time=character(), Subsite=character(), Age=character(),
                             Count=numeric(), Molting=logical(), Peak_breeding=logical(), Year=numeric())

pup_breeding <- rbind(pup_breeding, DP_pup_b)
pup_breeding <- rbind(pup_breeding, DE_pup_b)
pup_breeding <- rbind(pup_breeding, BL_pup_b)
pup_breeding <- rbind(pup_breeding, TB_pup_b)
pup_breeding <- rbind(pup_breeding, TP_pup_b)
pup_breeding <- rbind(pup_breeding, PR_pup_b)
pup_breeding <- rbind(pup_breeding, DR_pup_b)
pup_breeding <- rbind(pup_breeding, PB_pup_b)

######################
### MOLTING SEASON ###
######################

hseals_molt <- subset(hseals, Molting == TRUE)

# drop breeding column

# unique(hseals_molt$Age) gives "ADULT" only

DP_adult_m <- get_max(hseals_molt, "DP") #not missing any years
DE_adult_m <- get_max(hseals_molt, "DE") #not missing any years
BL_adult_m <- get_max(hseals_molt, "BL") #not missing any years
TB_adult_m <- get_max(hseals_molt, "TB") #not missing any years
TP_adult_m <- get_max(hseals_molt, "TP") #not missing any years
PR_adult_m <- get_max(hseals_molt, "PR") #not missing any years
DR_adult_m <- get_max(hseals_molt, "DR") #not missing any years
PB_adult_m <- get_max(hseals_molt, "PB") #not missing any years

# put it all back in one table for plotting
molting <- data.frame(Date=double(), Time=character(), Subsite=character(), Age=character(),
                           Count=numeric(), Molting=logical(), Peak_breeding=logical(), Year=numeric())

molting <- rbind(molting, DP_adult_m)
molting <- rbind(molting, DE_adult_m)
molting <- rbind(molting, BL_adult_m)
molting <- rbind(molting, TB_adult_m)
molting <- rbind(molting, TP_adult_m)
molting <- rbind(molting, PR_adult_m)
molting <- rbind(molting, DR_adult_m)
molting <- rbind(molting, PB_adult_m)



################
### PLOTTING ###
################

plot_gg <- function(DATA){
  ggplot(DATA, aes(Year, Count)) +
    geom_line(color = "azure4") +
    geom_point(size = 1, colour="salmon") +
    theme_gray(base_size = 14) +
    facet_grid(Subsite ~ Age) #try facet_wrap ?
}

pup_breeding$Age[pup_breeding$Age == "PUP"] <- "Pup Breeding Season"
adult_breeding$Age[adult_breeding$Age == "ADULT"] <- "Adult Breeding Season"
molting$Age[molting$Age == "ADULT"] <- "Adult Molting Season"

plot_gg(pup_breeding)
plot_gg(adult_breeding)
plot_gg(molting)

