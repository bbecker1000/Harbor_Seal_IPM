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

top.plot <- ggplot(top1_all_data, aes(Year, Count)) + 
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

seasonal.plot <- plot.breeding <- ggplot(all_data, aes(Julian, Count, colour = Year)) + 
  geom_point() + 
  geom_smooth(aes(colour = Year))
seasonal.plot + facet_grid(Year ~ Age)

## remove Point Bonita, PRH, and Duxbury since no or few pups
all_data.MARSS <- subset(top1_all_data, Subsite != "DR" & Subsite != "PB" & Subsite != "PRH")
unique(all_data.MARSS$Subsite)

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

MEI <- MEI %>% filter(Year < 2023)

#if NPGO
#NPGO <- MEI[,c(1,3)]
## make NPGO a vector
#NPGO <- NPGO[,2]
#if Coyote


plot(MEI$NPGO_MAR, MEI$NPGO_JAN )
plot(MEI$NPGO_MAR, MEI$MOCI_JFM_NC )
plot(MEI$PDO_MAR, MEI$MOCI_JFM_NC )
plot(MEI$Year, MEI$PDO_MAR)
plot(MEI$MEIv2_DEC_JAN, MEI$NPGO_MAR)
cor.test(MEI$PDO_MAR, as.numeric(MEI$MOCI_JFM_NC))
# cor = 0.79 (0.6-0.9) P < 0.001
# So let's use PDO_MAR as our env covariate.'
plot(MEI$PDO_MAR, MEI$NPGO_MAR  )
plot(MEI$MEIv2_DEC_JAN, MEI$PDO_MAR  )
plot(MEI$MOCI_JFM_NC, MEI$BEUTI_FEB_APR_37N_39N  )
plot(MEI$Year, MEI$MOCI_JFM_NC)

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


##########################
#Two approaches/analyses
#  1. use 5 main sites from 1998 - 2022 for THREE classes (ADULT, MOLT, PUP)
#  2. use 5 main sites back to 1970s for TWO age classes (ADULT, PUP) since no MOLT data pre ~1997
#  3. consider running model for all 8 sites for either partial (3 age classes) of full (2 age classes) time series...


########################################################################

## Parameter key ########
## Z = design matrix = Spatial population structure
## R = observation errors          equal
## U = growth parameter            (un)equal, or matrix [Z x 1]
## Q = hidden state process        diagonal and (un)equal
#                                  equalvarcov or unconstrained
## B = effect of column on row     unequal (these are the interactions)

#####################################################################

## B in -------------------------------------------
#  terms of the interaction strengths between species; bij
#  equals dfi /dXj, the change in the log population growth
#  rate of species i with respect to changes in the log
#  population abundance of species j. T

#Data structure----------------------------------
# "BL_ADULT" "BL_MOLT"  "BL_PUP"   "DE_ADULT" "DE_MOLT"  
# "DE_PUP"   "DP_ADULT" "DP_MOLT"  "DP_PUP"   "TB_ADULT"
# "TB_MOLT"  "TB_PUP"   "TP_ADULT" "TP_MOLT"  "TP_PUP" 

# Z models/Hypotheses ---------------------------
#H01 INDEPENDENT pops and classes -> factor(c(1:15))
#H02 By POPS                      -> factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5))
#H03 By CLASSES                   -> factor(c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3)) #no coyote covariate
#H04 Molt Independent             -> factor(c(1,2,1,1,2,1,1,2,1,1,2,1,1,2,1)) #no coyote covariate
#H05 Molt Independent by Site
#H06 Pup Independent by Site      -> factor(c(1,1,2,3,3,4,5,5,6,7,7,8,9,9,10)) 
#H07 BL Independent               -> factor(c(1,1,1,2,2,2,2,2,2,2,2,2,2,2,2))
#H08 North/South/BL               -> factor(c(1,1,1,2,2,2,2,2,2,3,3,3,3,3,3))

#H09 3 Pops pup ind               -> factor(c(1,1,2,3,3,4,3,3,4,5,5,6,5,5,6))
#H10 3 Pops molt ind              -> factor(c(1,2,1,3,4,3,3,4,3,5,6,5,5,6,5))

#####################################
## 2023-12-15 continue here #########
#####################################
#H11 Ocean/Bay                    -> factor(c(1,1,1,1,1,1,2,2,2,1,1,1,2,2,2))
#H12 North/South                  -> factor(c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2))

## MARSS Models----------------------------------------
# AIC table setup
df_aic <- data.frame(model=character(), aic=integer())

#H1 - All sites and classes Independent-------------------
#RT 1 min
Z.model=factor(c(1:15))
R.model="diagonal and equal"
U.model="unequal"
Q.model="diagonal and equal"
B.model="unconstrained"  # > 90 min
C.model=matrix(              list("BL_A",0,0,0,0,"MOCI_A",
                                  "BL_M",0,0,0,0,"MOCI_M",
                                  "BL_P",0,0,0,0,"MOCI_P",

                                  0,"DE_A",0,0,0,"MOCI_A",
                                  0,"DE_M",0,0,0,"MOCI_M",
                                  0,"DE_P",0,0,0,"MOCI_P",

                                  0,0,"DP_A",0,0,"MOCI_A",
                                  0,0,"DP_M",0,0,"MOCI_M",
                                  0,0,"DP_P",0,0,"MOCI_P",

                                  0,0,0,0,0,"MOCI_A",
                                  0,0,0,0,0,"MOCI_M",
                                  0,0,0,0,0,"MOCI_P",

                                  0,0,0,0,0,"MOCI_A",
                                  0,0,0,0,0,"MOCI_M",
                                  0,0,0,0,0,"MOCI_P"),
                            
                                  nrow = 15, ncol = 6,
                                  byrow = TRUE)
m.5pop_Coyote_PDO_B_unc=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                            C = C.model,
                                            c = small_c_Coyote_01_PDO_MAR),
                                            control=list(maxit=5000, safe=TRUE, trace = 0)) 

df_aic <- df_aic %>% add_row(model = "H01B - m.5pop_Coyote_PDO_B_unc", aic = m.5pop_Coyote_PDO_B_unc$AIC)
df_aic

# END H1 ------------------------------------------------------------

#H1B
# one pop, all ages same, unconstrained movement between sites.
#RT 81min
Z.model=factor(c(1:15)) # (rep(1, times = 15))
R.model="diagonal and equal"
U.model="equal"  # One big population  =  growth rate overall
Q.model="equalvarcov"
B.model="unconstrained"  #want to make the two groups of sites N and S more likely to interact?  break up by ages too?



# B.model= matrix(list("BL_BL_AA","BL_BL_MA","BL_BL_PA","DE_BL_AA","DE_BL_MA","DE_BL_PA","DP_BL_AA","DP_BL_MA","DP_BL_PA","TB_BL_AA","TB_BL_MA","TB_BL_PA","TP_BL_AA","TP_BL_MA","TP_BL_PA",
#                      "BL_BL_AM","BL_BL_MM","BL_BL_AP","DE_AA","DE_MA","DE_PA","DP_AA","DP_MA","DP_PA","TB_AA","TB_MA","TB_PA","TP_AA","TP_MA","TP_PA",
                              # 
                              # nrow = 15, ncol = 15,
                              # byrow = TRUE)
C.model=matrix(              list("BL_A",0,0,0,0,"PDO_A",
                                  "BL_M",0,0,0,0,"PDO_M",
                                  "BL_P",0,0,0,0,"PDO_P",
                                  
                                  0,"DE_A",0,0,0,"PDO_A",
                                  0,"DE_M",0,0,0,"PDO_M",
                                  0,"DE_P",0,0,0,"PDO_P",
                                  
                                  0,0,"DP_A",0,0,"PDO_A",
                                  0,0,"DP_M",0,0,"PDO_M",
                                  0,0,"DP_P",0,0,"PDO_P",
                                  
                                  0,0,0,0,0,"PDO_A",
                                  0,0,0,0,0,"PDO_M",
                                  0,0,0,0,0,"PDO_P",
                                  
                                  0,0,0,0,0,"PDO_A",
                                  0,0,0,0,0,"PDO_M",
                                  0,0,0,0,0,"PDO_P"),
                             
                             nrow = 15, ncol = 6,
                             byrow = TRUE)
m.1pop_Coyote_PDO_B_unc=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                              C = C.model,
                                              c = small_c_Coyote_01_PDO_MAR),
                              control=list(maxit=5000, safe=TRUE, trace = 0)) 

df_aic <- df_aic %>% add_row(model = "H01B - m.m.1pop_Coyote_PDO_B_unc", aic = m.1pop_Coyote_PDO_B_unc$AIC)
df_aic

beepr::beep()

#################
#H1C
# Molt and pup/Adult pops
# Z = 15 sites
# R = equal
# U = "equal" or by age
# Q = whatever
# B = "unconstrained" or molt vs pup/adult

Z.model=factor(c(1:15)) # (rep(1, times = 15))
R.model="diagonal and equal"
U.model="equal"  # One big population  =  growth rate overall
Q.model="equalvarcov"
B.model="unconstrained"  #want to make the two groups of sites N and S more likely to interact?  break up by ages too?

# B.matrix
#BL, DE, DP, TB, TP
# A / M / P  pup =  zero, molt to molt, adult to adult

                        
# B.model= matrix(list("BL_BL_AA","BL_BL_MA","BL_BL_PA","DE_BL_AA","DE_BL_MA","DE_BL_PA","DP_BL_AA","DP_BL_MA","DP_BL_PA","TB_BL_AA","TB_BL_MA","TB_BL_PA","TP_BL_AA","TP_BL_MA","TP_BL_PA",
# #                    BL_BL_AM","BL_BL_MM","BL_BL_PA","DE_BL_AA","DE_BL_MA","DE_BL_PA","DP_BL_AA","DP_BL_MA","DP_BL_PA","TB_BL_AA","TB_BL_MA","TB_BL_PA","TP_BL_AA","TP_BL_MA","TP_BL_PA",
# #                    BL_BL_AM","BL_BL_AA","BL_BL_PP","DE_BL_AA","DE_BL_MA","DE_BL_PA","DP_BL_AA","DP_BL_MA","DP_BL_PA","TB_BL_AA","TB_BL_MA","TB_BL_PA","TP_BL_AA","TP_BL_MA","TP_BL_PA",

#                           BL,                  DE,                 DP,      TB,      TP
                            # A M P

matrix(              factor(list("1","16",0,   "18","19",0,   0,0,0,   0,0,0,   0,0,0,     #A
                          0,"2",0,   0,0,0,   0,0,0,   0,0,0,   0,0,0,   #BL        #M
                          0,0,"3",   0,0,0,   0,0,0,   0,0,0,   0,0,0,              #P
                          
                          0,0,0,   "4",0,0,   0,0,0,   0,0,0,   0,0,0,
                          0,0,0,   0,"5",0,   0,0,0,   0,0,0,   0,0,0,   #DE
                          0,0,0,   0,0,"6",   0,0,0,   0,0,0,   0,0,0,
                          
                          0,0,0,   0,0,0,   "7",0,0,   0,0,0,   0,0,0,
                          0,0,0,   0,0,0,   0,"8",0,   0,0,0,   0,0,0,   #DP
                          0,0,0,   0,0,0,   0,0,"9",   0,0,0,   0,0,0,
                          
                          0,0,0,   0,0,0,   0,0,0,   "10",0,0,   0,0,0,
                          0,0,0,   0,0,0,   0,0,0,   0,"11",0,   0,0,0,   #TB
                          0,0,0,   0,0,0,   0,0,0,   0,0,"12",   0,0,0,
                          
                          0,0,0,   0,0,0,   0,0,0,   0,0,0,   "13",0,0,
                          0,0,0,   0,0,0,   0,0,0,   0,0,0,   0,"14",0,   #TP
                          0,0,0,   0,0,0,   0,0,0,   0,0,0,   0,0,"15"),
                     
                     nrow = 15, ncol = 6,
                     byrow = TRUE))

f <- (LETTERS(1:225))
B.model = matrix(f, nrow = 15, ncol = 15, byrow = TRUE)
B.model[3,3] <- as.integer(0)
B.model[6,6] <- as.integer(0)
B.model[9,9] <- as.integer(0)
B.model[12,12] <- as.integer(0)
B.model[15,15] <- as.integer(0)

                     
                     






# nrow = 15, ncol = 15,
# byrow = TRUE)
C.model=matrix(              list("BL_A",0,0,0,0,"PDO_A",
                                  "BL_M",0,0,0,0,"PDO_M",
                                  "BL_P",0,0,0,0,"PDO_P",
                                  
                                  0,"DE_A",0,0,0,"PDO_A",
                                  0,"DE_M",0,0,0,"PDO_M",
                                  0,"DE_P",0,0,0,"PDO_P",
                                  
                                  0,0,"DP_A",0,0,"PDO_A",
                                  0,0,"DP_M",0,0,"PDO_M",
                                  0,0,"DP_P",0,0,"PDO_P",
                                  
                                  0,0,0,0,0,"PDO_A",
                                  0,0,0,0,0,"PDO_M",
                                  0,0,0,0,0,"PDO_P",
                                  
                                  0,0,0,0,0,"PDO_A",
                                  0,0,0,0,0,"PDO_M",
                                  0,0,0,0,0,"PDO_P"),
                             
                             nrow = 15, ncol = 6,
                             byrow = TRUE)
m.1pop_Coyote_PDO_B_unc=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                              C = C.model,
                                              c = small_c_Coyote_01_PDO_MAR),
                              control=list(maxit=5000, safe=TRUE, trace = 0)) 

df_aic <- df_aic %>% add_row(model = "H01B - m.m.1pop_Coyote_PDO_B_unc", aic = m.1pop_Coyote_PDO_B_unc$AIC)
df_aic

beepr::beep()




## H2 3 pops 
Z.model=factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5))
R.model="diagonal and equal"
U.model="unequal"
Q.model="equalvarcov"
B.model="diagonal and equal"

C.model=matrix(list("BL_1",0,0,0,0,"MOCI",
                 0,"DE_1",0,0,0,"MOCI",
                 0,0,"DP_1",0,0,"MOCI",
                 0,0,0,0,0,"MOCI",
                 0,0,0,0,0,"MOCI"),
               nrow = 5, ncol = 6,
               byrow = TRUE)

m.3pop_Coyote_PDO=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                       C = C.model,
                                       c = small_c_Coyote_01_PDO_MAR),
                       control=list(maxit=3000, safe=TRUE, trace = 1)) 

df_aic <- df_aic %>% add_row(model = "H02 - m.3pop_Coyote_PDO", aic = m.3pop_Coyote_PDO$AIC)
df_aic


#H03 By CLASSES --------------
# cannot run logically with site level coyote covariates...just MOCI

Z.model=factor(c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3))
R.model="diagonal and equal"
U.model="unequal"
Q.model="equalvarcov" #diagonal and equal"
B.model="diagonal and equal"
C.model=matrix(list(0,0,0,0,0,"MOCI_ADULT",
                    0,0,0,0,0,"MOCI_MOLT",
                    0,0,0,0,0,"MOCI_PUP"),
               nrow = 3, ncol = 6,
               byrow = TRUE)
m.3Class_PDO=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                       C = C.model,
                                       c = small_c_Coyote_01_PDO_MAR),
                       control=list(maxit=3000, safe=TRUE, trace = 0)) 

df_aic <- df_aic %>% add_row(model = "H03 - m.3Class_PDO", aic = m.3Class_PDO$AIC)
df_aic



#H04 By CLASSES --------------

#H04 Molt Independent _no coyote can't use         

Z.model=factor(c(1,2,1,1,2,1,1,2,1,1,2,1,1,2,1))
R.model="diagonal and equal"
U.model="unequal"
Q.model="equalvarcov"
B.model="diagonal and equal"
C.model=matrix(list(0,0,0,0,0,"MOCI_ADULT_PUP",
                    0,0,0,0,0,"MOCI_MOLT"),
               nrow = 2, ncol = 6,
               byrow = TRUE)

m.MOLT_IND_PDO=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                  C = C.model,
                                                  c = small_c_Coyote_01_PDO_MAR),
                                  
                                  control=list(maxit=5000, safe=TRUE, trace = 1)) 

df_aic <- df_aic %>% add_row(model = "H04 - m.MOLT_IND_PDO", aic = m.MOLT_IND_PDO$AIC)
df_aic


##------------------------------------------------------------------

######'
######
######
# H05A Q_equalvarcov
# molt unrelated to breeding season at each site
Z.model=factor(c(1,2,1,3,4,3,5,6,5,7,8,7,9,10,9))
R.model="diagonal and equal"
U.model="unequal"
Q.model="equalvarcov" #correlated process errors for all pops
B.model="diagonal and equal"#"unequal"
C.model=matrix(list("BL_B",0,0,0,0, "MOCI_B", #1
                    "BL_M",0,0,0,0, "MOCI_M",#2
                    
                    0,"DE_B",0,0,0, "MOCI_B",#3
                    0,"DE_M",0,0,0, "MOCI_M",#4
                    
                    0,0,"DP_B",0,0, "MOCI_B",#5
                    0,0,"DP_M",0,0, "MOCI_M",#6
                    
                    0,0,0,0,0,"MOCI_B",  #7
                    0,0,0,0,0,"MOCI_M",  #8
                    
                    0,0,0,0,0,"MOCI_B",  #9
                    0,0,0,0,0,"MOCI_M"), #10
               nrow = 10, ncol = 6,
               byrow = TRUE)
m.MOLT_IND_5SITE_COYOTE_PDO=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                   C = C.model,
                                   c = small_c_Coyote_01_PDO_MAR),
                   control=list(maxit=5000, safe=TRUE)) 

df_aic <- df_aic %>% add_row(model = "H05A - m.MOLT_IND_5SITE_COYOTE_PDO", aic = m.MOLT_IND_5SITE_COYOTE_PDO$AIC)
df_aic



###----------------------------------------
# H05B Q_unconstrained
# molt unrelated to breeding season at each site
Z.model=factor(c(1,2,1,3,4,3,5,6,5,7,8,7,9,10,9))
R.model="diagonal and equal"
U.model="unequal"
Q.model="unconstrained"  #correlated process errors that vary between pops
B.model="diagonal and equal"#"unequal"
C.model=matrix(list("BL_B",0,0,0,0, "MOCI_B", #1
                    "BL_M",0,0,0,0, "MOCI_M",#2
                    
                    0,"DE_B",0,0,0, "MOCI_B",#3
                    0,"DE_M",0,0,0, "MOCI_M",#4
                    
                    0,0,"DP_B",0,0, "MOCI_B",#5
                    0,0,"DP_M",0,0, "MOCI_M",#6
                    
                    0,0,0,0,0,"MOCI_B",  #7
                    0,0,0,0,0,"MOCI_M",  #8
                    
                    0,0,0,0,0,"MOCI_B",  #9
                    0,0,0,0,0,"MOCI_M"), #10
               nrow = 10, ncol = 6,
               byrow = TRUE)
m.MOLT_IND_SITE_PDO_Q_unc=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                           C = C.model,
                                           c = small_c_Coyote_01_PDO_MAR),
                           control=list(maxit=5000, safe=TRUE)) 

CIs <- MARSSparamCIs(m.MOLT_IND_SITE_PDO_Q_unc, alpha = 0.11)
CIs

beepr::beep(0)
m.MOLT_IND_SITE_PDO_Q_unc$AIC 

df_aic <- df_aic %>% add_row(model = "H05B - m.MOLT_IND_SITE_PDO_Q_unc", aic = m.MOLT_IND_SITE_PDO_Q_unc$AIC)
df_aic


###############
# H05C Q_equalvarcov / B. diagonal and unequal
# runtime ~ 7 min
# molt unrelated to breeding season at each site
Z.model=factor(c(1,2,1,3,4,3,5,6,5,7,8,7,9,10,9))
R.model="diagonal and equal"
U.model="unequal"
Q.model="equalvarcov"  #correlated process errors that vary between pops
B.model="diagonal and unequal"#"unequal"
C.model=matrix(list("BL_B",0,0,0,0, "MOCI_B",#1
                    "BL_M",0,0,0,0, "MOCI_M",#2
                    
                    0,"DE_B",0,0,0, "MOCI_B",#3
                    0,"DE_M",0,0,0, "MOCI_M",#4
                    
                    0,0,"DP_B",0,0, "MOCI_B",#5
                    0,0,"DP_M",0,0, "MOCI_M",#6
                    
                    0,0,0,0,0,"MOCI_B",  #7
                    0,0,0,0,0,"MOCI_M",  #8
                    
                    0,0,0,0,0,"MOCI_B",  #9
                    0,0,0,0,0,"MOCI_M"), #10
               nrow = 10, ncol = 6,
               byrow = TRUE)
m.MOLT_IND_SITE_PDO.equal_B_diag_uneq=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                 C = C.model,
                                                 c = small_c_Coyote_01_PDO_MAR),
                                 control=list(maxit=5000, safe=TRUE)) 

df_aic <- df_aic %>% add_row(model = "H05C - m.MOLT_IND_SITE_PDO.equal_B_diag_uneq", aic = m.MOLT_IND_SITE_PDO.equal_B_diag_uneq$AIC)
df_aic

###############

# H05D Q_equalvarcov / B. diagonal and unequal  No covariates
# runtime ~ 7 min
# molt unrelated to breeding season at each site
Z.model=factor(c(1,2,1,3,4,3,5,6,5,7,8,7,9,10,9))
R.model="diagonal and equal"
U.model="unequal"
Q.BL=matrix(list("offdiag"), 10, 10)
        diag(Q.BL) <- "q"
        Q.BL[, c(1,2)] <- 0
        Q.BL[c(1,2), ] <- 0
        Q.BL[1,1]<- "q.BL"
        Q.BL[2,2]<- "q.BL"

Q.BL
Q.model = Q.BL

#correlated process errors that vary between pops
B.model="diagonal and unequal"#"unequal"

H05D_m.MOLT_IND_SITE.equal_B_diag_uneq_Q.BL=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model),
                                                            # C = C.model,
                                                            # c = small_c_Coyote_01_MOCI),
                                             control=list(maxit=5000, safe=TRUE)) 

df_aic <- df_aic %>% add_row(model = "H05D - m.MOLT_IND_SITE.equal_B_diag_uneq_Q.BLq", 
                             aic = H05D_m.MOLT_IND_SITE.equal_B_diag_uneq_Q.BL$AIC)
df_aic

###############
# H05E Q_equalvarcov / B. diagonal and unequal  remove MOCI !
# runtime ~ 7 min
# molt unrelated to breeding season at each site
Z.model=factor(c(1,2,1,3,4,3,5,6,5,7,8,7,9,10,9))
R.model="diagonal and equal"
U.model="unequal"
Q.model="equalvarcov"  #correlated process errors that vary between pops
B.model="diagonal and unequal"#"unequal"
C.model=matrix(list("BL_B",0,0,0,0, #1
                    "BL_M",0,0,0,0,#2
                    
                    0,"DE_B",0,0,0,#3
                    0,"DE_M",0,0,0,#4
                    
                    0,0,"DP_B",0,0,#5
                    0,0,"DP_M",0,0,#6
                    
                    0,0,0,0,0,  #7
                    0,0,0,0,0,  #8
                    
                    0,0,0,0,0,  #9
                    0,0,0,0,0), #10
               nrow = 10, ncol = 5,
               byrow = TRUE)
m.MOLT_IND_SITE.equal_B_diag_uneq=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                             C = C.model,
                                                             c = small_c_Coyote_01),
                                             control=list(maxit=5000, safe=TRUE)) 

df_aic <- df_aic %>% add_row(model = "H05E - m.MOLT_IND_SITE.equal_B_diag_uneq", aic = m.MOLT_IND_SITE.equal_B_diag_uneq$AIC)
df_aic

###############
###############
# H05F Q_equalvarcov / B. diagonal and unequal  NO COVARIATES!
# runtime ~ start 1124
# molt unrelated to breeding season at each site
Z.model=factor(c(1,2,1,3,4,3,5,6,5,7,8,7,9,10,9))
R.model="diagonal and equal"
U.model="unequal"
Q.model="equalvarcov"  #correlated process errors that vary between pops
B.model="diagonal and unequal"#"unequal"

m.MOLT_IND_SITE.equal_B_diag_uneq_NO_COVARIATES=MARSS(dat, model=list(Z=Z.model, U=U.model, 
                                                                      Q=Q.model, R=R.model, B=B.model))

df_aic <- df_aic %>% add_row(model = "H05F - m.MOLT_IND_SITE.equal_B_diag_uneq_NO_COVARIATES", aic = m.MOLT_IND_SITE.equal_B_diag_uneq_NO_COVARIATES$AIC)
df_aic


###############
# H05G Q_equalvarcov / B. diagonal and unequal     small_c_Coyote_01_BEAUTI
# runtime ~ 7 min
# molt unrelated to breeding season at each site
Z.model=factor(c(1,2,1,3,4,3,5,6,5,7,8,7,9,10,9))
R.model="diagonal and equal"
U.model="unequal"
Q.model="equalvarcov"  #correlated process errors that vary between pops
B.model="diagonal and unequal"#"unequal"
C.model=matrix(list("BL_B",0,0,0,0, "BEUTI_B", #1
                    "BL_M",0,0,0,0, "BEUTI_M",#2
                    
                    0,"DE_B",0,0,0, "BEUTI_B",#3
                    0,"DE_M",0,0,0, "BEUTI_M",#4
                    
                    0,0,"DP_B",0,0, "BEUTI_B",#5
                    0,0,"DP_M",0,0, "BEUTI_M",#6
                    
                    0,0,0,0,0,"BEUTI_B",  #7
                    0,0,0,0,0,"BEUTI_M",  #8
                    
                    0,0,0,0,0,"BEUTI_B",  #9
                    0,0,0,0,0,"BEUTI_M"), #10
               nrow = 10, ncol = 6,
               byrow = TRUE)
m.MOLT_IND_SITE_BEUTI.equal_B_diag_uneq=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                             C = C.model,
                                                             c = small_c_Coyote_01_BEUTI),
                                             control=list(maxit=5000, safe=TRUE)) 

df_aic <- df_aic %>% add_row(model = "H05G - m.MOLT_IND_SITE_BEUTI.equal_B_diag_uneq", aic = m.MOLT_IND_SITE_BEUTI.equal_B_diag_uneq$AIC)
df_aic

###############
# H05H Q_equalvarcov / B. diagonal and unequal  coyote_rate and MOCI.
# runtime ~ 7 min
# molt unrelated to breeding season at each site
Z.model=factor(c(1,2,1,3,4,3,5,6,5,7,8,7,9,10,9))
R.model="diagonal and equal"
U.model="unequal"
Q.model="equalvarcov"  #correlated process errors that vary between pops
B.model="diagonal and unequal"#"unequal"
C.model=matrix(list("BL_B",0,0,0,0, "MOCI_B", #1
                    "BL_M",0,0,0,0, "MOCI_M",#2
                    
                    0,"DE_B",0,0,0, "MOCI_B",#3
                    0,"DE_M",0,0,0, "MOCI_M",#4
                    
                    0,0,"DP_B",0,0, "MOCI_B",#5
                    0,0,"DP_M",0,0, "MOCI_M",#6
                    
                    0,0,0,0,0,"MOCI_B",  #7
                    0,0,0,0,0,"MOCI_M",  #8
                    
                    0,0,0,0,0,"MOCI_B",  #9
                    0,0,0,0,0,"MOCI_M"), #10
               nrow = 10, ncol = 6,
               byrow = TRUE)
m.MOLT_IND_SITE_Coy_rate_PDO.equal_B_diag_uneq=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                             C = C.model,
                                                             c = small_c_Coyote_01_PDO_MAR),
                                             control=list(maxit=5000, safe=TRUE)) 

df_aic <- df_aic %>% add_row(model = "H05H - m.MOLT_IND_SITE_Coy_rate_PDO.equal_B_diag_uneq", aic = m.MOLT_IND_SITE_Coy_rate_PDO.equal_B_diag_uneq$AIC)
df_aic

I
###############
# H05I Q_equalvarcov / B. diagonal and unequal  remove MOCI !
# runtime ~ 7 min
# molt unrelated to breeding season at each site
Z.model=factor(c(1,2,1,3,4,3,5,6,5,7,8,7,9,10,9))
R.model="diagonal and equal"
U.model="unequal"
Q.model="unconstrained"  #correlated process errors that vary between pops
B.model="diagonal and unequal"#"unequal"
C.model=matrix(list("BL_B",0,0,0,0, #1
                    "BL_M",0,0,0,0,#2
                    
                    0,"DE_B",0,0,0,#3
                    0,"DE_M",0,0,0,#4
                    
                    0,0,"DP_B",0,0,#5
                    0,0,"DP_M",0,0,#6
                    
                    0,0,0,0,0,  #7
                    0,0,0,0,0,  #8
                    
                    0,0,0,0,0,  #9
                    0,0,0,0,0), #10
               nrow = 10, ncol = 5,
               byrow = TRUE)
m.MOLT_IND_SITE_B_diag_uneq_q_unc=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                        C = C.model,
                                                        c = small_c_Coyote_01),
                                        control=list(maxit=5000, safe=TRUE)) 

df_aic <- df_aic %>% add_row(model = "H05I - m.MOLT_IND_SITE_B_diag_uneq_q_unc", aic = m.MOLT_IND_SITE_B_diag_uneq_q_unc$AIC)
df_aic

beepr::beep()

###############

###############
# H05J Q_equalvarcov / B. diagonal and unequal  B unconstrained !
# runtime ~ 7 min
# molt unrelated to breeding season at each site
Z.model=factor(c(1,2,1,3,4,3,5,6,5,7,8,7,9,10,9))
R.model="diagonal and equal"
U.model="unequal"
Q.model="equalvarcov"  #correlated process errors that vary between pops
B.model="unconstrained"#"unequal"
C.model=matrix(list("BL_B",0,0,0,0, #1
                    "BL_M",0,0,0,0,#2
                    
                    0,"DE_B",0,0,0,#3
                    0,"DE_M",0,0,0,#4
                    
                    0,0,"DP_B",0,0,#5
                    0,0,"DP_M",0,0,#6
                    
                    0,0,0,0,0,  #7
                    0,0,0,0,0,  #8
                    
                    0,0,0,0,0,  #9
                    0,0,0,0,0), #10
               nrow = 10, ncol = 5,
               byrow = TRUE)
m.MOLT_IND_SITE.equal_B_unc=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                        C = C.model,
                                                        c = small_c_Coyote_01),
                                        control=list(maxit=5000, safe=TRUE)) 

df_aic <- df_aic %>% add_row(model = "H05J - m.MOLT_IND_SITE.equal_B_unc", aic = m.MOLT_IND_SITE.equal_B_unc$AIC)
df_aic



######-----------
# Pups unrelated to breeding adults and molt season at each site
# coyote
# H06
#RT = 1 min
Z.model=factor(c(1,1,2,3,3,4,5,5,6,7,7,8,9,9,10))
R.model="diagonal and equal"
U.model="unequal"
Q.model="equalvarcov"
B.model="diagonal and equal"#"unequal"
C.model=matrix(list("BL_ADULT",0,0,0,0, "MOCI_ADULT", #1
                    "BL_PUP",0,0,0,0, "MOCI_PUP",#2
                    
                    0,"DE_ADULT",0,0,0, "MOCI_ADULT",#3
                    0,"DE_PUP",0,0,0, "MOCI_PUP",#4
                    
                    0,0,"DP_ADULT",0,0, "MOCI_ADULT",#5
                    0,0,"DP_PUP",0,0, "MOCI_PUP",#6
                    
                    0,0,0,0,0,"MOCI_ADULT",  #7
                    0,0,0,0,0,"MOCI_PUP",  #8
                    
                    0,0,0,0,0,"MOCI_ADULT",  #9
                    0,0,0,0,0,"MOCI_PUP"), #10
               nrow = 10, ncol = 6,
               byrow = TRUE)
m.PUP_IND_SITE_PDO=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                           C = C.model,
                                           c = small_c_Coyote_01_PDO_MAR),
                           control=list(maxit=5000, safe=TRUE)) 

df_aic <- df_aic %>% add_row(model = "H06 - m.PUP_IND_SITE_PDO", aic = m.PUP_IND_SITE_PDO$AIC)
df_aic

###############
#H07 to do


##############------------
#H08 North/South/BL               -> 
Z.model=factor(c(1,1,1,2,2,2,2,2,2,3,3,3,3,3,3))
R.model="diagonal and equal"
U.model="unequal"
Q.model="equalvarcov"
B.model="diagonal and equal"#"unequal"
C.model=matrix(list(0,0,0,0,0, "MOCI", #1
                    0,0,0,0,0, "MOCI",#2
                    0,0,0,0,0, "MOCI"), #10
               nrow = 3, ncol = 6,
               byrow = TRUE)
m.BL_SOUTH_NORTH_PDO=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                          C = C.model,
                                          c = small_c_Coyote_01_PDO_MAR),
                          control=list(maxit=5000, safe=TRUE)) 

df_aic <- df_aic %>% add_row(model = "H08 - m.BL_SOUTH_NORTH_PDO", aic = m.BL_SOUTH_NORTH_PDO$AIC)
df_aic

###-------------------

#H09 3 Pops pup ind               -
Z.model=factor(c(1,1,2,3,3,4,3,3,4,5,5,6,5,5,6))
R.model="diagonal and equal"
U.model="unequal"
Q.model="equalvarcov"
B.model="diagonal and equal"#"unequal"
C.model=matrix(list("BL_ADULT",0,0,0,0, "MOCI_ADULT", #1
                    "BL_PUP",0,0,0,0, "MOCI_PUP",#2
                    
                    0,0,0,0,0, "MOCI_ADULT",#3
                    0,0,0,0,0, "MOCI_PUP",#4
                    
                    0,0,0,0,0, "MOCI_ADULT",#5
                    0,0,0,0,0, "MOCI_PUP"),#6
                    nrow = 6, ncol = 6,
                    byrow = TRUE)
m.BL_SOUTH_NORTH_PUP_IND_PDO=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                            C = C.model,
                                            c = small_c_Coyote_01_PDO_MAR),
                            control=list(maxit=5000, safe=TRUE)) 

df_aic <- df_aic %>% add_row(model = "H09 - m.BL_SOUTH_NORTH_PUP_IND_PDO", aic = m.BL_SOUTH_NORTH_PUP_IND_PDO$AIC)
df_aic

#####----------------------------------

#H10 3 Pops molt ind              -
Z.model=factor(c(1,2,1,3,4,3,3,4,3,5,6,5,5,6,5))
R.model="diagonal and equal"
U.model="unequal"
Q.model="equalvarcov"
B.model="diagonal and equal"#"unequal"
C.model=matrix(list("BL_ADULT",0,0,0,0, "MOCI_Breed", #1
                    "BL_PUP",0,0,0,0, "MOCI_Molt",#2
                    
                    0,0,0,0,0, "MOCI_Breed",#3
                    0,0,0,0,0, "MOCI_Molt",#4
                    
                    0,0,0,0,0, "MOCI_Breed",#5
                    0,0,0,0,0, "MOCI_Molt"),#6
               nrow = 6, ncol = 6,
               byrow = TRUE)
m.BL_SOUTH_NORTH_MOLT_IND_PDO=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                    C = C.model,
                                                    c = small_c_Coyote_01_PDO_MAR),
                                    control=list(maxit=5000, safe=TRUE)) 

df_aic <- df_aic %>% add_row(model = "H10 - m.BL_SOUTH_NORTH_MOLT_IND_PDO", aic = m.BL_SOUTH_NORTH_MOLT_IND_PDO$AIC)
df_aic


################

##############------------
#H11 2 molt ind   
# start 903
Z.model=factor(c(1,2,1,3,4,3,3,4,3,3,4,3,3,4,3))
R.model="diagonal and equal"
U.model="unequal"
Q.model="equalvarcov"  #correlated process errors that vary between pops
B.model="diagonal and unequal"#"unequal"
C.model=matrix(list("BL_ADULT",0,0,0,0, "MOCI_Breed", #1
                    "BL_PUP",0,0,0,0, "MOCI_Molt",#2
                    
                    0,0,0,0,0, "MOCI_Breed",#3
                    0,0,0,0,0, "MOCI_Molt"),#4
               nrow = 4, ncol = 6,
               byrow = TRUE)
m.2POP_MOLT_IND_PDO=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                     C = C.model,
                                                     c = small_c_Coyote_01_PDO_MAR),
                                     control=list(maxit=5000, safe=TRUE)) 

df_aic <- df_aic %>% add_row(model = "H11 - m.2POP_MOLT_IND_PDO", aic = m.2POP_MOLT_IND_PDO$AIC)
df_aic


beepr::beep()


####-----------
## assume all 5 sites part of one population and that there is both growth AND movement between sites based on coyotes, etc.
## molt and pup sepaprate?
####


Z.model=factor(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
R.model="diagonal and equal"
U.model="unequal"
Q.model="equalvarcov"  #correlated process errors that vary between pops
B.model="diagonal and unequal"#"unequal"
C.model=matrix(list("BL_B",0,0,0,0, "PDO_B", #1
                    "BL_M",0,0,0,0, "PDO_M",#2
                    
                    0,"DE_B",0,0,0, "PDO_B",#3
                    0,"DE_M",0,0,0, "_M",#4
                    
                    0,0,"DP_B",0,0, "MOCI_B",#5
                    0,0,"DP_M",0,0, "MOCI_M",#6
                    
                    0,0,0,0,0,"MOCI_B",  #7
                    0,0,0,0,0,"MOCI_M",  #8
                    
                    0,0,0,0,0,"MOCI_B",  #9
                    0,0,0,0,0,"MOCI_M"), #10
               nrow = 10, ncol = 6,
               byrow = TRUE)
m.ONE_POP_MOLT_IND_SITE_Coy_PDO.equal_B_diag_uneq=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                                     C = C.model,
                                                                     c = small_c_Coyote_01_PDO_MAR),
                                                     control=list(maxit=5000, safe=TRUE)) 






#############################
# Model plots -------------------------------------
#############################
BESTMODEL <-  m.1pop_Coyote_PDO_B_unc #   m.1pop_Coyote_PDO_B_unc   m.5pop_Coyote_PDO_B_unc



autoplot(BESTMODEL)

autoplot(BESTMODEL, plot.type = "fitted.ytT") + # xtT
  ylim(3.7,7.6) # xtT

autoplot(BESTMODEL, plot.type = "xtT") + # xtT
  ylim(3.1,7.6)

# plot temporal change from starting values

#BESTMODEL <- m.MOLT_IND_SITE_MOCI.equal_B_diag_uneq
CIs <- MARSSparamCIs(BESTMODEL, alpha = 0.11)
CIs <- MARSSparamCIs(BESTMODEL, alpha = 0.11, hessian.fun = "fdHess")

CIs

CIs$states.se



d <- as_tibble(t(BESTMODEL$states-BESTMODEL$states[,1]))
d.se <- as_tibble(t(BESTMODEL$states.se))
#add header names
names(d)[c(1:10)] <- c("BL_Breed", "BL_Molt", "DE_Breed", "DE_Molt", 
                       "DP_Breed", "DP_Molt", "TB_Breed", "TB_Molt", 
                       "TP_Breed", "TP_Molt")
names(d.se)[c(1:10)] <- c("BL_Breed", "BL_Molt", "DE_Breed", "DE_Molt", 
                       "DP_Breed", "DP_Molt", "TB_Breed", "TB_Molt", 
                       "TP_Breed", "TP_Molt")


d <- as_tibble(cbind(years,d))

d2 <- d %>% pivot_longer(cols = c(2:16), names_to = "Subsite_Season", values_to = "log_est")
d2.se <- d.se %>% pivot_longer(cols = c(1:15), names_to = "Subsite_Season", values_to = "log_est")
d2.se <- d2.se[,2]
colnames(d2.se)[1] = "log_se"

d4 <- as_tibble(cbind(d2,d2.se))

#add groupings for site and molt
d5 <- d4 %>% separate_wider_delim(Subsite_Season, "_", names = c("Subsite", "Season"))
d5$Subsite_Season <- paste0(d5$Subsite, "_", d5$Season)
d5

ggplot(d5, aes(x = years, y = log_est, group = Subsite_Season, color = Subsite)) +
  #geom_line(aes(linetype = Season), size = 1.25) +
  geom_point(size = 2) + 
  geom_line(size = 1.1) +
  geom_ribbon(aes(ymin = log_est-log_se, ymax = log_est+log_se),  
              alpha = 0.15, colour = NA) +
  geom_hline(yintercept = c(-1,0,1), lty = 2) +
  xlim(1997, 2022) +
  ylim(-1.5, 1.8) + 
  theme_gray(base_size = 20) +
  ylab("index of log abundance (N-No)") +
  xlab("Year") +
  facet_wrap(.~Season)

####-----plot covariate effects

coef.data <- tidy(CIs)

#get just the coefficient terms
coef.data <- coef.data %>%
  filter(str_detect(term, "^C."))

coef.data$Season <- ifelse(str_detect(coef.data$term, "_B"), "Breeding", "Molting")
#coef.data$Covariate <- ifelse(str_detect(coef.data$term, "MOCI"), "MOCI", "Coyote")
coef.data$Covariate <- ifelse(str_detect(coef.data$term, "PDO"), "PDO", "Coyote")
  
ggplot(coef.data, aes(term, estimate, color = Season, shape = Season)) +

  geom_pointrange(aes(ymin = conf.low, ymax = conf.up)) +
  coord_flip() +
  theme(legend.position = c(0.2, 0.2)) + 

  # scale_x_discrete(name ="Covariate",
  #                  labels=c("Coyote BL","Coyote BL","Coyote DE", "Coyote DE",
  #                           "Coyote DP","Coyote DP", "MOCI", "MOCI"), 
  #                  limits=rev(levels(coef.data$term))) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_sjplot(base_size = 18)


########################################################################

## Parameter key ########
## Z = design matrix = Spatial population structure
## R = observation errors          equal
## U = growth parameter            (un)equal, or matrix [Z x 1]
## Q = hidden state process        diagonal and (un)equal
#                                  equalvarcov or unconstrained
## B = effect of column on row     unequal (these are the interactions)

#####################################################################
## some coefs
coef(BESTMODEL, type="matrix")$R  # observation errors  
coef(BESTMODEL, type="matrix")$Q  # hidden state process
coef(BESTMODEL, type="matrix")$U  # growth
coef(m.MOLT_IND_SITE.equal_B_unc, type="matrix")$B  # effect of column on row


m.MOLT_IND_SITE_B_diag_uneq_q_unc

## Process correlation between sites --------
#recover q.matrix p 115-116
# use m.MOLT_IND_SITE_B_diag_uneq_q_unc to show all correlations for Q



Q.unc <- coef(m.MOLT_IND_SITE.equal_B_unc, type="matrix")$Q
h <- diag(1 / sqrt(diag(Q.unc)))
Q.corr <- h %*% Q.unc %*% h

## get Z.model from appropriate model
Z.model=factor(c(1,2,1,3,4,3,5,6,5,7,8,7,9,10,9)) 

rownames(Q.corr) <- unique(Z.model)
colnames(Q.corr) <- unique(Z.model)
# Process Correlations between sites
Q.corr






