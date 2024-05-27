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
  theme_gray(base_size = 14) + 
  scale_x_continuous(breaks=seq(1980,2025,20))
  
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
small_c_Coyote_3yr_UI_UI_lag <- as.matrix(t(small_c_Coyote_3yr_UI_UI_lag)) # ready for model 2024-05-26

PDO_MAR <- MEI[,6]
small_c_PDO_MAR <- as.matrix(t(PDO_MAR))


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
#ignore PUPS this run !!!!
# all sites

# Z models/Hypotheses ---------------------------


## MARSS Models----------------------------------------
# AIC table setup
df_aic <- data.frame(model=character(), aic=integer())


Sys.time()
#Model 1A - 12 indiv pops 6 sites and Adult and molt -------------------
# years = 1982-2022

Z.model=factor(c(1:12))  # 6 sites x 2 age classes
R.model="diagonal and equal"  #zero since focusing on B  else  "diagonal and equal"
#U.model="equal"
U.model=matrix(list("UA", "UM","UA", "UM","UA", "UM","UA", "UM","UA", "UM","UA", "UM"), 
               nrow = 12, ncol = 1,
               byrow = TRUE)# separate growth rates for molt and adult
Q.model="diagonal and equal"   #near zero add control
B.model="unconstrained"  # > 90 min
C.model=matrix(              list("BL_A",0,0,0,0,0,"PDO_A",
                                  "BL_M",0,0,0,0,0,"PDO_M",
                                #  "BL_P",0,0,0,0," PDO_P",

                                  0,"DE_A",0,0,0,0,"PDO_A",
                                  0,"DE_M",0,0,0,0,"PDO_M",
                               #   0,"DE_P",0,0,0," PDO_P",

                                  0,0,"DP_A",0,0,0,"PDO_A",
                                  0,0,"DP_M",0,0,0,"PDO_M",
                               #   0,0,"DP_P",0,0," PDO_P",

                                  # 0,0,0,0,0,0,0,0,"PDO_A",
                                  # 0,0,0,0,0,0,0,0,"PDO_M",       #DR
                               #   0,0,0,0,0," PDO_P",

                                  # 0,0,0,0,0,0,0,0,"PDO_A",       #PB
                                  # 0,0,0,0,0,0,0,0,"PDO_M",
                               #   0,0,0,0,0," PDO_P"
                               
                               0,0,0,0,0,0,"PDO_A",         #PRH
                               0,0,0,0,0,0,"PDO_M",
                               
                               0,0,0,0,0,0,"PDO_A",         #TB
                               0,0,0,0,0,0,"PDO_M",
                               
                               0,0,0,0,0,0,"PDO_A",         #TP
                               0,0,0,0,0,0,"PDO_M"
                               ),
                                  nrow = 12, ncol = 7,
                                  byrow = TRUE)
m.Ind_Molt_Adult_Coyote_PDO_B_unc_2U=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                            C = C.model,
                                            #x0 = x0.model, 
                                            tinitx=1, 
                                            c = small_c_Coyote_01_PDO_MAR),
                                            control=list(maxit=5000, safe=TRUE, 
                                                         trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "Model_1_A - m.Ind_Molt_Adult_Coyote_PDO_B_unc_2U", aic = m.2pop_Molt_Adult_Coyote_PDO_B_unc_2U$AIC)
df_aic
beepr::beep()
Sys.time()

save(m.2pop_Molt_Adult_Coyote_PDO_B_unc_2U, file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_unc_2U.RData")
#load(file = "Output/m.2pop_Molt_Adult_Coyote_PDO_B_unc_2U.RData")




# END H1 ------------------------------------------------------------

#Model 1B - 12 individual pops 6 sites Adult and molt have same U -------------------
# years = 1982-2022
# RT 25 min
t0 <- Sys.time()

Z.model=factor(c(1:12))  # 6 sites x 2 age classes
R.model="diagonal and equal"  #zero since focusing on B  else  "diagonal and equal"
U.model="equal"
Q.model="diagonal and equal"   #near zero add control
B.model="unconstrained"  # > 90 min
# x0.model=matrix(c(3.6, 4.4,# 1.7,  #BL  # estimated from first data point in time series at each site_age combo
#                   5.5, 6, # 3,       #DE
#                   5.3, 6, # 4.5,    #DP
#                   #2.5, 2.5,# 2,  # DR
#                   #3, 3, # 2,      # PB
#                   3, 3, # 2,      #PRH
#                   5, 5.5, # 4.5,   #TB
#                   6.1, 6.4), # 4.8), #TP
#                 nrow = 12, ncol = 1,
#                 byrow = TRUE)

C.model=matrix(              list("BL_A",0,0,0,0,0,"PDO_A",
                                  "BL_M",0,0,0,0,0,"PDO_M",
                                  #  "BL_P",0,0,0,0," PDO_P",
                                  
                                  0,"DE_A",0,0,0,0,"PDO_A",
                                  0,"DE_M",0,0,0,0,"PDO_M",
                                  #   0,"DE_P",0,0,0," PDO_P",
                                  
                                  0,0,"DP_A",0,0,0,"PDO_A",
                                  0,0,"DP_M",0,0,0,"PDO_M",
                                  #   0,0,"DP_P",0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",
                                  # 0,0,0,0,0,0,0,0,"PDO_M",       #DR
                                  #   0,0,0,0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",       #PB
                                  # 0,0,0,0,0,0,0,0,"PDO_M",
                                  #   0,0,0,0,0," PDO_P"
                                  
                                  0,0,0,0,0,0,"PDO_A",         #PRH
                                  0,0,0,0,0,0,"PDO_M",
                                  
                                  0,0,0,0,0,0,"PDO_A",         #TB
                                  0,0,0,0,0,0,"PDO_M",
                                  
                                  0,0,0,0,0,0,"PDO_A",         #TP
                                  0,0,0,0,0,0,"PDO_M"
),
nrow = 12, ncol = 7,
byrow = TRUE)
m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                             C = C.model,
                                                             #x0 = x0.model, 
                                                             tinitx=1, 
                                                             c = small_c_Coyote_01_PDO_MAR),
                                             control=list(maxit=5000, safe=TRUE, 
                                                          trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "Model_1_B - m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U", aic = m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U$AIC)
df_aic
beepr::beep()
Sys.time()
save(m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U, file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U.RData")
#load(file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U.RData")

#end time
t1 <- Sys.time()
#Runtime
t1-t0

# END H1 ------------------------------------------------------------


#Model 1C - 12 individual pops 6 sites Adult and molt have same U 
# custom B matrix   -------------------
# years = 1982-2022

#Start Time
t0 <- Sys.time()

Z.model=factor(c(1:12))  # 6 sites x 2 age classes
R.model="zero"  #zero since focusing on B  else  "diagonal and equal"
U.model="equal"
Q.model="diagonal and equal"   #near zero add control
B.model=matrix(list(0,0,"K",0,   "U",0,"e",0,   "p",0,"z",0,  #BL
                    0,0,0,"P",     0,"Z",0,"k",   0,"u",0,"ee",
                    
                    "A",0,0,0,   "V",0,"f",0,   "q",0,"aa",0, #DE
                    0,"F",0,0,   0,"a",0,"l",   0,"v",0,"ff",
                    
                    "B",0,"L",0,   0,0,"g",0,   "r",0,"bb",0,  #DP
                    0,"G",0,"Q",   0,0,0,"m",   0,"w",0,"gg",
                    
                    "C",0,"M",0, "W",0,0,0,   "s",0,"cc",0,  #PRH
                    0,"H",0,"R",   0,"b",0,0,   0,"x",0,"hh",
                    
                    "D",0,"N",0, "X",0,"i",0,   0,0,"dd",0,  #TB
                    0,"I",0,"S",   0,"c",0,"n",   0,0,0,"ii",
                    
                    "E",0,"O",0, "Y",0,"j",0,   "t",0,0,0,  #TP
                    0,"J",0,"T",   0,"d",0,"o",   0,"y",0,0),
               
                    nrow = 12, ncol = 12,
                    byrow = TRUE)
 
x0.model=matrix(c(3.6, 4.4,# 1.7,  #BL  # estimated from first data point in time series at each site_age combo
                  5.5, 6, # 3,       #DE
                  5.3, 6, # 4.5,    #DP
                  #2.5, 2.5,# 2,  # DR
                  #3, 3, # 2,      # PB
                  3, 3, # 2,      #PRH
                  5, 5.5, # 4.5,   #TB
                  6.1, 6.4), # 4.8), #TP
                nrow = 12, ncol = 1,
                byrow = TRUE)
  

C.model=matrix(              list("BL_A",0,0,0,0,0,"PDO_A",
                                  "BL_M",0,0,0,0,0,"PDO_M",
                                  #  "BL_P",0,0,0,0," PDO_P",
                                  
                                  0,"DE_A",0,0,0,0,"PDO_A",
                                  0,"DE_M",0,0,0,0,"PDO_M",
                                  #   0,"DE_P",0,0,0," PDO_P",
                                  
                                  0,0,"DP_A",0,0,0,"PDO_A",
                                  0,0,"DP_M",0,0,0,"PDO_M",
                                  #   0,0,"DP_P",0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",
                                  # 0,0,0,0,0,0,0,0,"PDO_M",       #DR
                                  #   0,0,0,0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",       #PB
                                  # 0,0,0,0,0,0,0,0,"PDO_M",
                                  #   0,0,0,0,0," PDO_P"
                                  
                                  0,0,0,0,0,0,"PDO_A",         #PRH
                                  0,0,0,0,0,0,"PDO_M",
                                  
                                  0,0,0,0,0,0,"PDO_A",         #TB
                                  0,0,0,0,0,0,"PDO_M",
                                  
                                  0,0,0,0,0,0,"PDO_A",         #TP
                                  0,0,0,0,0,0,"PDO_M"
),
nrow = 12, ncol = 7,
byrow = TRUE)
m.1pop_Molt_Adult_Coyote_PDO_B_custom_1U_R0=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                            C = C.model,
                                                            x0 = x0.model, 
                                                            tinitx=1, 
                                                            #method = "BFGS",
                                                            c = small_c_Coyote_01_PDO_MAR),
                                            control=list(maxit=5000, safe=TRUE, 
                                                         trace = 0, allow.degen=TRUE
                                                           ))

df_aic <- df_aic %>% add_row(model = "Model_1_C - m.1pop_Molt_Adult_Coyote_PDO_B_custom_1U_R0", aic = m.1pop_Molt_Adult_Coyote_PDO_B_custom_1U_R0$AIC)
df_aic
beepr::beep()
Sys.time()
save(m.1pop_Molt_Adult_Coyote_PDO_B_custom_1U, file = "Output/m.1pop_Molt_Adult_Coyote_PDO_B_custom_1U_R0.RData")
#load(file = "Output/m.1pop_Molt_Adult_Coyote_PDO_B_custom_1U_R0.RData")

#end time
t1 <- Sys.time()
#Runtime
t1-t0

# END H1 ------------------------------------------------------------


t0 <- Sys.time()

#Model 1C - 2 pops by age 6 sites and Adult and molt -------------------
# years = 1982-2022
# RT 4.7 min
Z.model=factor(c(1,2,1,2,1,2,1,2,1,2,1,2))  # 6 sites x 2 age classes 
## How get covariates for site???!!

Z.model <- matrix(0, 12, 2)
Z.model[c(1:8 ), 1] <- 1 # which elements in col 1 are 1
Z.model[c(9:12), 2] <- 1 # which elements in col 2 are 1
Z.model




R.model="diagonal and equal"  #zero since focusing on B  else  "diagonal and equal"
U.model="unequal"
# U.model=matrix(list("UA", "UM","UA", "UM","UA", "UM","UA", "UM","UA", "UM","UA", "UM"), 
#                nrow = 12, ncol = 1,
#                byrow = TRUE)# separate growth rates for molt and adult
Q.model="diagonal and equal"   #near zero add control
B.model="diagonal and equal"  # > 90 min
C.model=matrix(              list(0,0,0,0,0,0,"PDO_A",
                                  0,0,0,0,0,0,"PDO_M"                     
),
nrow = 2, ncol = 7,
byrow = TRUE)
m.age_pop_Molt_Adult_Coyote_PDO_B_unc_2U=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                             C = C.model,
                                                             #x0 = x0.model, 
                                                             tinitx=1, 
                                                             c = small_c_Coyote_01_PDO_MAR),
                                             control=list(maxit=5000, safe=TRUE, 
                                                          trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "Model_1_C - m.age_pop_Molt_Adult_Coyote_PDO_B_unc_2U", aic = m.age_pop_Molt_Adult_Coyote_PDO_B_unc_2U$AIC)
df_aic
beepr::beep()


save(m.age_pop_Molt_Adult_Coyote_PDO_B_unc_2U, file = "Output/m.age_pop_Molt_Adult_Coyote_PDO_B_unc_2U.RData")
#load(file = "Output/m.age_pop_Molt_Adult_Coyote_PDO_B_unc_2U.RData")


#end time
t1 <- Sys.time()
#Runtime
t1-t0


# END H1B ------------------------------------------------------------

#Model 1D - 2 pops by site (BL/DP/DE/PRH vs. TB/TP) 6 sites and Adult and molt -------------------
# years = 1982-202
#start 11 min 1U
#statrt 2U 11 min
t0 <- Sys.time()

Z.model=factor(c(1,1,1,1,1,1,1,1,2,2,2,2))  # 6 sites x 2 age classes
R.model="diagonal and equal"  #zero since focusing on B  else  "diagonal and equal"
U.model="unequal"
Q.model="diagonal and equal"   #near zero add control
B.model="unconstrained"  # > 90 min

C.model=matrix(              list(0,0,0,0,0,0,"PDO_S",
                                  0,0,0,0,0,0,"PDO_N"                     
),
nrow = 2, ncol = 7,
byrow = TRUE)
m.N_S_pop_Molt_Adult_Coyote_PDO_B_unc_2U=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                               C = C.model,
                                                               #x0 = x0.model, 
                                                               tinitx=1, 
                                                               c = small_c_Coyote_01_PDO_MAR),
                                               control=list(maxit=5000, safe=TRUE, 
                                                            trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "Model_1_C - m.N_S_pop_Molt_Adult_Coyote_PDO_B_unc_2U", aic = m.N_S_pop_Molt_Adult_Coyote_PDO_B_unc_2U$AIC)
df_aic
beepr::beep()
#end time
t1 <- Sys.time()
#Runtime
t1-t0

save(m.N_S_pop_Molt_Adult_Coyote_PDO_B_unc_2U, file = "Output/m.N_S_pop_Molt_Adult_Coyote_PDO_B_unc_2U.RData")
#load(file = "Output/m.N_S_pop_Molt_Adult_Coyote_PDO_B_unc_2U.RData")



# END H1B ------------------------------------------------------------


#Model 1E - 12 individual pops 6 sites Adult and molt all DIFFERENT U's -------------------
# years = 1982-2022
# RT 25 min
t0 <- Sys.time()

Z.model=factor(c(1:12))  # 6 sites x 2 age classes
Q.model="diagonal and equal"   #near zero add control
B.model="unconstrained"  # maybe instead of unconstrained
R.model="diagonal and equal"  #zero since focusing on B  else  "diagonal and equal"
U.model="unequal"
               Q.model="diagonal and equal"   #near zero add control
               B.model="unconstrained"
               C.model=matrix(              list("BL_A",0,0,0,0,0,"PDO_A",
                                                 "BL_M",0,0,0,0,0,"PDO_M",
                                                 #  "BL_P",0,0,0,0," PDO_P",
                                                 
                                                 0,"DE_A",0,0,0,0,"PDO_A",
                                                 0,"DE_M",0,0,0,0,"PDO_M",
                                                 #   0,"DE_P",0,0,0," PDO_P",
                                                 
                                                 0,0,"DP_A",0,0,0,"PDO_A",
                                                 0,0,"DP_M",0,0,0,"PDO_M",
                                                 #   0,0,"DP_P",0,0," PDO_P",
                                                 
                                                 # 0,0,0,0,0,0,0,0,"PDO_A",
                                                 # 0,0,0,0,0,0,0,0,"PDO_M",       #DR
                                                 #   0,0,0,0,0," PDO_P",
                                                 
                                                 # 0,0,0,0,0,0,0,0,"PDO_A",       #PB
                                                 # 0,0,0,0,0,0,0,0,"PDO_M",
                                                 #   0,0,0,0,0," PDO_P"
                                                 
                                                 0,0,0,0,0,0,"PDO_A",         #PRH
                                                 0,0,0,0,0,0,"PDO_M",
                                                 
                                                 0,0,0,0,0,0,"PDO_A",         #TB
                                                 0,0,0,0,0,0,"PDO_M",
                                                 
                                                 0,0,0,0,0,0,"PDO_A",         #TP
                                                 0,0,0,0,0,0,"PDO_M"
               ),
               nrow = 12, ncol = 7,
               byrow = TRUE)
               m.Ind_Molt_Adult_Coyote_PDO_B_unc_12U=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                                          C = C.model,
                                                                          #x0 = x0.model, 
                                                                          tinitx=1, 
                                                                          c = small_c_Coyote_01_PDO_MAR),
                                                          control=list(maxit=5000, safe=TRUE, 
                                                                       trace = 0, allow.degen=TRUE)) 
               
               df_aic <- df_aic %>% add_row(model = "Model_1_E - m.Ind_Molt_Adult_Coyote_PDO_B_unc_12U", aic = m.1pop_Molt_Adult_Coyote_PDO_B_unc_2U$AIC)
               df_aic
               beepr::beep()
               Sys.time()
               save(m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U, file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_unc_12U.RData")
               #load(file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U.RData")
               
               #end time
               t1 <- Sys.time()
               #Runtime
               t1-t0
               



#Model 1E - 12 individual pops 6 sites Adult and molt all DIFFERENT U's 
# B identity -------------------
# years = 1982-2022
# RT 25 min
t0 <- Sys.time()
               
               Z.model=factor(c(1:12))  # 6 sites x 2 age classes
               Q.model="diagonal and equal"   #near zero add control
               B.model="diagonal and equal"  # maybe instead of unconstrained
               R.model="diagonal and equal"  #zero since focusing on B  else  "diagonal and equal"
               U.model="unequal"
               C.model=matrix(              list("BL_A",0,0,0,0,0,"PDO_A",
                                                 "BL_M",0,0,0,0,0,"PDO_M",
                                                 #  "BL_P",0,0,0,0," PDO_P",
                                                 
                                                 0,"DE_A",0,0,0,0,"PDO_A",
                                                 0,"DE_M",0,0,0,0,"PDO_M",
                                                 #   0,"DE_P",0,0,0," PDO_P",
                                                 
                                                 0,0,"DP_A",0,0,0,"PDO_A",
                                                 0,0,"DP_M",0,0,0,"PDO_M",
                                                 #   0,0,"DP_P",0,0," PDO_P",
                                                 
                                                 # 0,0,0,0,0,0,0,0,"PDO_A",
                                                 # 0,0,0,0,0,0,0,0,"PDO_M",       #DR
                                                 #   0,0,0,0,0," PDO_P",
                                                 
                                                 # 0,0,0,0,0,0,0,0,"PDO_A",       #PB
                                                 # 0,0,0,0,0,0,0,0,"PDO_M",
                                                 #   0,0,0,0,0," PDO_P"
                                                 
                                                 0,0,0,0,0,0,"PDO_A",         #PRH
                                                 0,0,0,0,0,0,"PDO_M",
                                                 
                                                 0,0,0,0,0,0,"PDO_A",         #TB
                                                 0,0,0,0,0,0,"PDO_M",
                                                 
                                                 0,0,0,0,0,0,"PDO_A",         #TP
                                                 0,0,0,0,0,0,"PDO_M"
               ),
               nrow = 12, ncol = 7,
               byrow = TRUE)
               m.Ind_Molt_Adult_Coyote_PDO_B_ident_12U=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                                           C = C.model,
                                                                           #x0 = x0.model, 
                                                                           tinitx=1, 
                                                                           c = small_c_Coyote_01_PDO_MAR),
                                                           control=list(maxit=5000, safe=TRUE, 
                                                                        trace = 0, allow.degen=TRUE)) 
               
               df_aic <- df_aic %>% add_row(model = "Model_1_E - m.Ind_Molt_Adult_Coyote_PDO_B_ident_12U", aic = m.Ind_Molt_Adult_Coyote_PDO_B_ident_12U$AIC)
               df_aic
               beepr::beep()
               Sys.time()
               save(m.Ind_Molt_Adult_Coyote_PDO_B_ident_12U, file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_ident_12U.RData")
               #load(file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_ident_12U.RData")
               
               #end time
               t1 <- Sys.time()
               #Runtime
               t1-t0
               
###########-----------------
# 2024-02-06  
               
#Model 1B - 12 individual pops 6 sites Adult and molt have same U -------------------
# years = 1982-2022
# RT 25 min
t0 <- Sys.time()

Z.model=factor(c(1:12))  # 6 sites x 2 age classes
R.model="diagonal and equal"  #zero since focusing on B  else  "diagonal and equal"
U.model="equal"
Q.model="diagonal and equal"   #near zero add control
# nearby site B.model matrix
# assume same processes for adult and molt.
B.model=matrix(list(0,0,    "A",0,   "B",0,   0,0,   0,0,  0,0,  #BL
                    0,0,     0,"A",   0,"B",  0,0,   0,0,  0,0,
                    
                    "C",0,  0,0,   "D",0,  "E",0,    0,0,  0,0, #DE
                    0,"C",  0,0,   0,"D",  0,"E",    0,0,  0,0,
                    
                    "F",0,  "G",0,   0,0,  "H",0,   0,0,  0,0,  #DP
                    0,"F",  0,"G",   0,0,  0,"H",   0,0,  0,0,
                    
                    "I",0,  "J",0, "K",0,    0,0,   "L",0,  "M",0,  #PRH
                    0,"I",  0,"J",   0,"K",  0,0,   0,"L",  0,"M",
                    
                    0,0,   0,0,     0,0,   "N",0,   0,0,   "O",0,  #TB
                    0,0,   0,0,     0,0,   0,"N",   0,0,   0,"O",
                    
                    0,0,   0,0,     0,0,   0,0,    "P",0,   0,0,  #TP
                    0,0,   0,0,     0,0,   0,0,    0,"P",   0,0),
               
               nrow = 12, ncol = 12,
               byrow = TRUE)

C.model=matrix(              list("BL_A",0,0,0,0,0,"PDO_A",
                                  "BL_M",0,0,0,0,0,"PDO_M",
                                  #  "BL_P",0,0,0,0," PDO_P",
                                  
                                  0,"DE_A",0,0,0,0,"PDO_A",
                                  0,"DE_M",0,0,0,0,"PDO_M",
                                  #   0,"DE_P",0,0,0," PDO_P",
                                  
                                  0,0,"DP_A",0,0,0,"PDO_A",
                                  0,0,"DP_M",0,0,0,"PDO_M",
                                  #   0,0,"DP_P",0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",
                                  # 0,0,0,0,0,0,0,0,"PDO_M",       #DR
                                  #   0,0,0,0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",       #PB
                                  # 0,0,0,0,0,0,0,0,"PDO_M",
                                  #   0,0,0,0,0," PDO_P"
                                  
                                  0,0,0,0,0,0,"PDO_A",         #PRH
                                  0,0,0,0,0,0,"PDO_M",
                                  
                                  0,0,0,0,0,0,"PDO_A",         #TB
                                  0,0,0,0,0,0,"PDO_M",
                                  
                                  0,0,0,0,0,0,"PDO_A",         #TP
                                  0,0,0,0,0,0,"PDO_M"
),
nrow = 12, ncol = 7,
byrow = TRUE)
m.Ind_Molt_Adult_Coyote_PDO_B_Close_1U=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                             C = C.model,
                                                             #x0 = x0.model, 
                                                             tinitx=1, 
                                                             c = small_c_Coyote_01_PDO_MAR),
                                             control=list(maxit=5000, safe=TRUE, 
                                                          trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "Model_1_F - m.Ind_Molt_Adult_Coyote_PDO_B_Close_1U", aic = m.Ind_Molt_Adult_Coyote_PDO_B_Close_1U$AIC)
df_aic
beepr::beep()
Sys.time()
save(m.Ind_Molt_Adult_Coyote_PDO_B_Close_1U, file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_Close_1U.RData")
#load(file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_Close_1U.RData")

#end time
t1 <- Sys.time()
#Runtime
t1-t0

###########-----------------
# 2024-02-06  

#Model 1B - 12 individual pops 6 sites Adult and molt have same U -------------------
# years = 1982-2022
# RT 16 min
t0 <- Sys.time()

Z.model=factor(c(1:12))  # 6 sites x 2 age classes
R.model="diagonal and equal" # #zero since (need x0) focusing on B  else  "diagonal and equal"
<<<<<<< HEAD
=======
# x0.model=matrix(c(3.6, 4.4,# 1.7,  #BL  # estimated from first data point in time series at each site_age combo
#                   5.5, 6, # 3,       #DE
#                   5.3, 6, # 4.5,    #DP
#                   #2.5, 2.5,# 2,  # DR
#                   #3, 3, # 2,      # PB
#                   3, 3, # 2,      #PRH
#                   5, 5.5, # 4.5,   #TB
#                   6.1, 6.4), # 4.8), #TP
#                 nrow = 12, ncol = 1,
#                 byrow = TRUE)

>>>>>>> 9f3ae972c3d6bdc8c30283f3e165a6ef9ef929a0
U.model="equal"
Q.model="diagonal and equal"   #near zero add control

# nearby site B.model matrix
# assume same processes for adult and molt.
B.model=matrix(list(0,0,    "A",0,   "B",0,   0,0,   0,0,  0,0,  #BL
                    0,0,     0,"A",   0,"B",  0,0,   0,0,  0,0,
                    
                    "C",0,  0,0,   "D",0,  "E",0,    0,0,  0,0, #DE
                    0,"C",  0,0,   0,"D",  0,"E",    0,0,  0,0,
                    
                    "F",0,  "G",0,   0,0,  "H",0,   0,0,  0,0,  #DP
                    0,"F",  0,"G",   0,0,  0,"H",   0,0,  0,0,
                    
                    "I",0,  "J",0, "K",0,    0,0,   "L",0,  "M",0,  #PRH
                    0,"I",  0,"J",   0,"K",  0,0,   0,"L",  0,"M",
                    
                    0,0,   0,0,     0,0,   "N",0,   0,0,   "O",0,  #TB
                    0,0,   0,0,     0,0,   0,"N",   0,0,   0,"O",
                    
                    0,0,   0,0,     0,0,   0,0,    "P",0,   0,0,  #TP
                    0,0,   0,0,     0,0,   0,0,    0,"P",   0,0),
               
               nrow = 12, ncol = 12,
               byrow = TRUE)

C.model=matrix(              list("BL_A",0,0,0,0,0,"PDO_A",
                                  "BL_M",0,0,0,0,0,"PDO_M",
                                  #  "BL_P",0,0,0,0," PDO_P",
                                  
                                  0,"DE_A",0,0,0,0,"PDO_A",
                                  0,"DE_M",0,0,0,0,"PDO_M",
                                  #   0,"DE_P",0,0,0," PDO_P",
                                  
                                  0,0,"DP_A",0,0,0,"PDO_A",
                                  0,0,"DP_M",0,0,0,"PDO_M",
                                  #   0,0,"DP_P",0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",
                                  # 0,0,0,0,0,0,0,0,"PDO_M",       #DR
                                  #   0,0,0,0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",       #PB
                                  # 0,0,0,0,0,0,0,0,"PDO_M",
                                  #   0,0,0,0,0," PDO_P"
                                  
                                  0,0,0,0,0,0,"PDO_A",         #PRH
                                  0,0,0,0,0,0,"PDO_M",
                                  
                                  0,0,0,0,0,0,"PDO_A",         #TB
                                  0,0,0,0,0,0,"PDO_M",
                                  
                                  0,0,0,0,0,0,"PDO_A",         #TP
                                  0,0,0,0,0,0,"PDO_M"
),
nrow = 12, ncol = 7,
byrow = TRUE)
m.Ind_Molt_Adult_Coyote_PDO_B_Close_12U=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                             C = C.model,
                                                             #x0 = x0.model, 
                                                             tinitx=1, 
                                                             c = small_c_Coyote_01_PDO_MAR),
                                             control=list(maxit=5000, safe=TRUE, 
                                                          trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "Model_1_G - m.Ind_Molt_Adult_Coyote_PDO_B_Close_12U", aic = m.Ind_Molt_Adult_Coyote_PDO_B_Close_12U$AIC)
df_aic
beepr::beep()
Sys.time()
save(m.Ind_Molt_Adult_Coyote_PDO_B_Close_12U, file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_Close_12U.RData")
#load(file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_Close_12U.RData")

#end time
t1 <- Sys.time()
#Runtime
t1-t0

               
               
               
#Model 1H - 12 individual pops 6 sites Adult and molt have same U -------------------
# years = 1982-2022
# RT 25 min
t0 <- Sys.time()

Z.model=factor(c(1:12))  # 6 sites x 2 age classes
R.model="diagonal and equal"  #zero since focusing on B  else  "diagonal and equal"
U.model="equal"
Q.model="equalvarcov"   #near zero add control
B.model="equalvarcov"  # > 90 min
<<<<<<< HEAD
=======
# x0.model=matrix(c(3.6, 4.4,# 1.7,  #BL  # estimated from first data point in time series at each site_age combo
#                   5.5, 6, # 3,       #DE
#                   5.3, 6, # 4.5,    #DP
#                   #2.5, 2.5,# 2,  # DR
#                   #3, 3, # 2,      # PB
#                   3, 3, # 2,      #PRH
#                   5, 5.5, # 4.5,   #TB
#                   6.1, 6.4), # 4.8), #TP
#                 nrow = 12, ncol = 1,
#                 byrow = TRUE)

>>>>>>> 9f3ae972c3d6bdc8c30283f3e165a6ef9ef929a0
C.model=matrix(              list("BL_A",0,0,0,0,0,"PDO_A",
                                  "BL_M",0,0,0,0,0,"PDO_M",
                                  #  "BL_P",0,0,0,0," PDO_P",
                                  
                                  0,"DE_A",0,0,0,0,"PDO_A",
                                  0,"DE_M",0,0,0,0,"PDO_M",
                                  #   0,"DE_P",0,0,0," PDO_P",
                                  
                                  0,0,"DP_A",0,0,0,"PDO_A",
                                  0,0,"DP_M",0,0,0,"PDO_M",
                                  #   0,0,"DP_P",0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",
                                  # 0,0,0,0,0,0,0,0,"PDO_M",       #DR
                                  #   0,0,0,0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",       #PB
                                  # 0,0,0,0,0,0,0,0,"PDO_M",
                                  #   0,0,0,0,0," PDO_P"
                                  
                                  0,0,0,0,0,0,"PDO_A",         #PRH
                                  0,0,0,0,0,0,"PDO_M",
                                  
                                  0,0,0,0,0,0,"PDO_A",         #TB
                                  0,0,0,0,0,0,"PDO_M",
                                  
                                  0,0,0,0,0,0,"PDO_A",         #TP
                                  0,0,0,0,0,0,"PDO_M"
),
nrow = 12, ncol = 7,
byrow = TRUE)
m.Ind_Molt_Adult_CoyoteSame_PDO_B_eqcov_1U=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                           C = C.model,
                                                           #x0 = x0.model, 
                                                           tinitx=1, 
                                                           c = small_c_Coyote_01_PDO_MAR),
                                           control=list(maxit=5000, safe=TRUE, 
                                                        trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "Model_1_H - m.Ind_Molt_Adult_CoyoteSame_PDO_B_eqcov_1U", aic = m.Ind_Molt_Adult_CoyoteSame_PDO_B_eqcov_1U$AIC)
df_aic
beepr::beep()
Sys.time()
save(m.Ind_Molt_Adult_CoyoteSame_PDO_B_eqcov_1U, file = "Output/m.Ind_Molt_Adult_CoyoteSame_PDO_B_eqcov_1U.RData")
#load(file = "m.Ind_Molt_Adult_CoyoteSame_PDO_B_eqcov_1U.RData")

#end time
t1 <- Sys.time()
#Runtime
t1-t0

# END H1 ------------------------------------------------------------               
               

               
                  
               
#timevarying u
# all sites have same growth rates that can change around 2003

U1 <- matrix("t1", 12, 1)
U2 <- matrix("t2", 12, 1)
Ut <- array(U2, dim = c(dim(U1), dim(dat)[2]))
TT <- dim(dat)[2]
Ut[, , 1:ceiling(TT / 2)] <- U1


## TIME VARYING MODEL same all

#Model 1B - 12 individual pops 6 sites Adult and molt have same U -------------------
# years = 1982-2022
# RT 34 min
t0 <- Sys.time()

Z.model=factor(c(1:12))  # 6 sites x 2 age classes

R.model <- diag(0.01, 12) #known observation error variance

R.model <- diag(0.1, 12) #known observation error variance

#R.model="equal"  #zero since focusing on B  else  "diagonal and equal"
U.model=Ut  #Select from above
Q.model="equalvarcov" #"diagonal and equal"   #near zero add control
B.model="unconstrained"  # > 90 min

# x0.model=matrix(c(3.6, 4.4,# 1.7,  #BL  # estimated from first data point in time series at each site_age combo
#                   5.5, 6, # 3,       #DE
#                   5.3, 6, # 4.5,    #DP
#                   #2.5, 2.5,# 2,  # DR
#                   #3, 3, # 2,      # PB
#                   3, 3, # 2,      #PRH
#                   5, 5.5, # 4.5,   #TB
#                   6.1, 6.4), # 4.8), #TP
#                 nrow = 12, ncol = 1,
#                 byrow = TRUE)


C.model=matrix(              list("BL_A",0,0,0,0,0,"PDO_A",
                                  "BL_M",0,0,0,0,0,"PDO_M",
                                  #  "BL_P",0,0,0,0," PDO_P",
                                  
                                  0,"DE_A",0,0,0,0,"PDO_A",
                                  0,"DE_M",0,0,0,0,"PDO_M",
                                  #   0,"DE_P",0,0,0," PDO_P",
                                  
                                  0,0,"DP_A",0,0,0,"PDO_A",
                                  0,0,"DP_M",0,0,0,"PDO_M",
                                  #   0,0,"DP_P",0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",
                                  # 0,0,0,0,0,0,0,0,"PDO_M",       #DR
                                  #   0,0,0,0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",       #PB
                                  # 0,0,0,0,0,0,0,0,"PDO_M",
                                  #   0,0,0,0,0," PDO_P"
                                  
                                  0,0,0,0,0,0,"PDO_A",         #PRH
                                  0,0,0,0,0,0,"PDO_M",
                                  
                                  0,0,0,0,0,0,"PDO_A",         #TB
                                  0,0,0,0,0,0,"PDO_M",
                                  
                                  0,0,0,0,0,0,"PDO_A",         #TP
                                  0,0,0,0,0,0,"PDO_M"
),
nrow = 12, ncol = 7,
byrow = TRUE)
m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_R_fix=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, 
                                                                    R=R.model, B=B.model,
                                                                    C = C.model,
                                                           #x0 = x0.model, 
                                                                    tinitx=1, 
                                                                    c = small_c_Coyote_01_PDO_MAR),
                                           control=list(maxit=5000, safe=TRUE, 
                                                        trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "Model_1_I - m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_R_fix", aic = m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_R_fix$AIC)
df_aic
beepr::beep()

save(m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_R_fix, file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_R_fix.RData")

#load(file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_R_fix.RData")

#end time
t1 <- Sys.time()
#Runtime
t1-t0

# END H1 ------------------------------------------------------------


#timevarying u by site
U1 <- matrix(c("t1_BL_A","t1_BL_M","t1_DE_A","t1_DE_M","t1_DP_A","t1_DP_M",
               "t1_PRH_A","t1_PRH_M","t1_TB_A","t1_TB_M","t1_TP_A","t1_TP_M"),12, 1)
U2 <- matrix(c("t2_BL_A","t2_BL_M","t2_DE_A","t2_DE_M","t2_DP_A","t2_DP_M",
               "t2_PRH_A","t2_PRH_M","t2_TB_A","t2_TB_M","t2_TP_A","t2_TP_M"),12, 1)
Ut <- array(U2, dim = c(dim(U1), dim(dat)[2]))
TT <- dim(dat)[2]
Ut[, , 1:ceiling(TT / 2)] <- U1

## TIME VARYING MODEL by site

#Model 1B - 12 individual pops 6 sites Adult and molt have same U -------------------
# years = 1982-2022
# RT 43 min
t0 <- Sys.time()

Z.model=factor(c(1:12))  # 6 sites x 2 age classes
<<<<<<< HEAD
#Z.model=factor(rep(1, 12))
R.model <- diag(0.05, 12) #known observation error variance
R.model="diagonal and equal"
U.model=Ut  #Select from above
Q.model="equalvarcov"#"diagonal and equal"   #near zero add control
B.model="equalvarcov"#unconstrained"  # > 90 min
=======
R.model <- diag(0.1, 12) #known observation error variance
#R.model="diagonal and equal"  #zero since focusing on B  else  "diagonal and equal"
U.model=Ut  #Select from above
Q.model="diagonal and equal"   #near zero add control
B.model="unconstrained"  # > 90 min
# x0.model=matrix(c(3.6, 4.4,# 1.7,  #BL  # estimated from first data point in time series at each site_age combo
#                   5.5, 6, # 3,       #DE
#                   5.3, 6, # 4.5,    #DP
#                   #2.5, 2.5,# 2,  # DR
#                   #3, 3, # 2,      # PB
#                   3, 3, # 2,      #PRH
#                   5, 5.5, # 4.5,   #TB
#                   6.1, 6.4), # 4.8), #TP
#                 nrow = 12, ncol = 1,
#                 byrow = TRUE)

>>>>>>> 9f3ae972c3d6bdc8c30283f3e165a6ef9ef929a0
C.model=matrix(              list("BL_A",0,0,0,0,0,"PDO_A",
                                  "BL_M",0,0,0,0,0,"PDO_M",
                                  #  "BL_P",0,0,0,0," PDO_P",
                                  
                                  0,"DE_A",0,0,0,0,"PDO_A",
                                  0,"DE_M",0,0,0,0,"PDO_M",
                                  #   0,"DE_P",0,0,0," PDO_P",
                                  
                                  0,0,"DP_A",0,0,0,"PDO_A",
                                  0,0,"DP_M",0,0,0,"PDO_M",
                                  #   0,0,"DP_P",0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",
                                  # 0,0,0,0,0,0,0,0,"PDO_M",       #DR
                                  #   0,0,0,0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",       #PB
                                  # 0,0,0,0,0,0,0,0,"PDO_M",
                                  #   0,0,0,0,0," PDO_P"
                                  
                                  0,0,0,0,0,0,"PDO_A",         #PRH
                                  0,0,0,0,0,0,"PDO_M",
                                  
                                  0,0,0,0,0,0,"PDO_A",         #TB
                                  0,0,0,0,0,0,"PDO_M",
                                  
                                  0,0,0,0,0,0,"PDO_A",         #TP
                                  0,0,0,0,0,0,"PDO_M"
),
nrow = 12, ncol = 7,
byrow = TRUE)
<<<<<<< HEAD
m.Ind_Molt_Adult_Coyote_PDO_B_equalcov_1U_TV_Site=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
=======
m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_Site=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
>>>>>>> 9f3ae972c3d6bdc8c30283f3e165a6ef9ef929a0
                                                              C = C.model,
                                                              #x0 = x0.model, 
                                                              tinitx=1, 
                                                              c = small_c_Coyote_01_PDO_MAR),
                                              control=list(maxit=5000, safe=TRUE, 
                                                           trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "Model_1_J - m.Ind_Molt_Adult_Coyote_PDO_B_equalcov_1U_TV_Site", aic = m.Ind_Molt_Adult_Coyote_PDO_B_equalcov_1U_TV_Site$AIC)
df_aic
beepr::beep()
save(m.Ind_Molt_Adult_Coyote_PDO_B_equalcov_1U_TV_Site, file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_equalcov_1U_TV_Site.RData")
#load(file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_equalcov_1U_TV_Site.RData")
=======
df_aic <- df_aic %>% add_row(model = "Model_1_J - m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_Site", aic = m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_Site$AIC)
df_aic
beepr::beep()
save(m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_Site, file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_Site.RData")
#load(file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_Site.RData")
>>>>>>> 9f3ae972c3d6bdc8c30283f3e165a6ef9ef929a0

#end time
t1 <- Sys.time()
#Runtime
t1-t0


#####__________________




#timevarying u by site
U1 <- matrix(c("t1_BL_A","t1_BL_M","t1_DE_A","t1_DE_M","t1_DP_A","t1_DP_M",
               "t1_PRH_A","t1_PRH_M","t1_TB_A","t1_TB_M","t1_TP_A","t1_TP_M"),12, 1)
U2 <- matrix(c("t2_BL_A","t2_BL_M","t2_DE_A","t2_DE_M","t2_DP_A","t2_DP_M",
               "t2_PRH_A","t2_PRH_M","t2_TB_A","t2_TB_M","t2_TP_A","t2_TP_M"),12, 1)
Ut <- array(U2, dim = c(dim(U1), dim(dat)[2]))
TT <- dim(dat)[2]
Ut[, , 1:ceiling(TT / 2)] <- U1

t0 <- Sys.time()

Z.model=factor(c(1:12))  # 6 sites x 2 age classes
#R.model <- diag(0.05, 12) #known observation error variance
R.model="diagonal and equal"  #zero since focusing on B  else  "diagonal and equal"
U.model=Ut  #Select from above
Q.model="equalvarcov"#"diagonal and equal"   #near zero add control
B.model=matrix(list(0,0,    "A",0,   "B",0,   0,0,   0,0,  0,0,  #BL
                    0,0,     0,"A",   0,"B",  0,0,   0,0,  0,0,
                    
                    "C",0,  0,0,   "D",0,  "E",0,    0,0,  0,0, #DE
                    0,"C",  0,0,   0,"D",  0,"E",    0,0,  0,0,
                    
                    "F",0,  "G",0,   0,0,  "H",0,   0,0,  0,0,  #DP
                    0,"F",  0,"G",   0,0,  0,"H",   0,0,  0,0,
                    
                    "I",0,  "J",0, "K",0,    0,0,   "L",0,  "M",0,  #PRH
                    0,"I",  0,"J",   0,"K",  0,0,   0,"L",  0,"M",
                    
                    0,0,   0,0,     0,0,   "N",0,   0,0,   "O",0,  #TB
                    0,0,   0,0,     0,0,   0,"N",   0,0,   0,"O",
                    
                    0,0,   0,0,     0,0,   0,0,    "P",0,   0,0,  #TP
                    0,0,   0,0,     0,0,   0,0,    0,"P",   0,0),
               
               nrow = 12, ncol = 12,
               byrow = TRUE)

C.model=matrix(              list("BL_A",0,0,0,0,0,"PDO_A",
                                  "BL_M",0,0,0,0,0,"PDO_M",
                                  #  "BL_P",0,0,0,0," PDO_P",
                                  
                                  0,"DE_A",0,0,0,0,"PDO_A",
                                  0,"DE_M",0,0,0,0,"PDO_M",
                                  #   0,"DE_P",0,0,0," PDO_P",
                                  
                                  0,0,"DP_A",0,0,0,"PDO_A",
                                  0,0,"DP_M",0,0,0,"PDO_M",
                                  #   0,0,"DP_P",0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",
                                  # 0,0,0,0,0,0,0,0,"PDO_M",       #DR
                                  #   0,0,0,0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",       #PB
                                  # 0,0,0,0,0,0,0,0,"PDO_M",
                                  #   0,0,0,0,0," PDO_P"
                                  
                                  0,0,0,0,0,0,"PDO_A",         #PRH
                                  0,0,0,0,0,0,"PDO_M",
                                  
                                  0,0,0,0,0,0,"PDO_A",         #TB
                                  0,0,0,0,0,0,"PDO_M",
                                  
                                  0,0,0,0,0,0,"PDO_A",         #TP
                                  0,0,0,0,0,0,"PDO_M"
),
nrow = 12, ncol = 7,
byrow = TRUE)
m.Ind_Molt_Adult_Coyote_PDO_B_SET_1U_TV_Site=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                                   C = C.model,
                                                                   #x0 = x0.model, 
                                                                   tinitx=1, 
                                                                   c = small_c_Coyote_01_PDO_MAR),
                                                   control=list(maxit=5000, safe=TRUE, 
                                                                trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "Model_1_K - m.Ind_Molt_Adult_Coyote_PDO_B_SET_1U_TV_Site", aic = m.Ind_Molt_Adult_Coyote_PDO_B_SET_1U_TV_Site$AIC)
df_aic
beepr::beep()
save(m.Ind_Molt_Adult_Coyote_PDO_B_SET_1U_TV_Site, file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_SET_1U_TV_Site.RData")
#load(file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_SET_1U_TV_Site.RData")

#end time
t1 <- Sys.time()
#Runtime
t1-t0




###########-----
#timevarying u by site
U1 <- matrix(c("t1_BL_A","t1_BL_M","t1_DE_A","t1_DE_M","t1_DP_A","t1_DP_M",
               "t1_PRH_A","t1_PRH_M","t1_TB_A","t1_TB_M","t1_TP_A","t1_TP_M"),12, 1)
U2 <- matrix(c("t2_BL_A","t2_BL_M","t2_DE_A","t2_DE_M","t2_DP_A","t2_DP_M",
               "t2_PRH_A","t2_PRH_M","t2_TB_A","t2_TB_M","t2_TP_A","t2_TP_M"),12, 1)
Ut <- array(U2, dim = c(dim(U1), dim(dat)[2]))
TT <- dim(dat)[2]
Ut[, , 1:ceiling(TT / 2)] <- U1

## TIME VARYING MODEL by site

#Model 1B - 12 individual pops 6 sites Adult and molt have same U -------------------
# years = 1982-2022
# RT 43 min
t0 <- Sys.time()

Z.model=factor(c(1:12))  # 6 sites x 2 age classes
#R.model <- diag(0.05, 12) #known observation error variance
R.model="diagonal and equal"  #zero since focusing on B  else  "diagonal and equal"
U.model=Ut  #Select from above
Q.model="equalvarcov"#"diagonal and equal"   #near zero add control
B.model="equalvarcov"#identity"#unconstrained"  # > 90 min

C.model=matrix(              list("BL_A",0,0,0,0,0,"PDO_A",
                                  "BL_M",0,0,0,0,0,"PDO_M",
                                  #  "BL_P",0,0,0,0," PDO_P",
                                  
                                  0,"DE_A",0,0,0,0,"PDO_A",
                                  0,"DE_M",0,0,0,0,"PDO_M",
                                  #   0,"DE_P",0,0,0," PDO_P",
                                  
                                  0,0,"DP_A",0,0,0,"PDO_A",
                                  0,0,"DP_M",0,0,0,"PDO_M",
                                  #   0,0,"DP_P",0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",
                                  # 0,0,0,0,0,0,0,0,"PDO_M",       #DR
                                  #   0,0,0,0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",       #PB
                                  # 0,0,0,0,0,0,0,0,"PDO_M",
                                  #   0,0,0,0,0," PDO_P"
                                  
                                  0,0,0,0,0,0,"PDO_A",         #PRH
                                  0,0,0,0,0,0,"PDO_M",
                                  
                                  0,0,0,0,0,0,"PDO_A",         #TB
                                  0,0,0,0,0,0,"PDO_M",
                                  
                                  0,0,0,0,0,0,"PDO_A",         #TP
                                  0,0,0,0,0,0,"PDO_M"
),
nrow = 12, ncol = 7,
byrow = TRUE)
m.Ind_Molt_Adult_Coyote_PDO_B_EQUALVAR_1U_TV_Site=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                                   C = C.model,
                                                                   #x0 = x0.model, 
                                                                   tinitx=1, 
                                                                   c = small_c_Coyote_01_PDO_MAR),
                                                   control=list(maxit=5000, safe=TRUE, 
                                                                trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "Model_1_J - m.Ind_Molt_Adult_Coyote_PDO_B_EQUALVAR_1U_TV_Site", aic = m.Ind_Molt_Adult_Coyote_PDO_B_EQUALVAR_1U_TV_Site$AIC)
df_aic
beepr::beep()
save(m.Ind_Molt_Adult_Coyote_PDO_B_EQUALVAR_1U_TV_Site, file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_EQUALVAR_1U_TV_Site.RData")
#load(file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_EQUALVAR_1U_TV_Site.RData")

#end time
t1 <- Sys.time()
#Runtime
t1-t0











# END H1 ------------------------------------------------------------


#Many warnings  !!!

#timevarying u by 3 sites
U1 <- matrix(c("t1_BL_A","t1_BL_M","t1_South_A","t1_South_M","t1_South_A","t1_South_A",
               "t1_PRH_A","t1_PRH_M","t1_North_A","t1_North_A","t1_North_A","t1_North_A"),12, 1)
U2 <- matrix(c("t2_BL_A","t2_BL_M","t2_South_A","t2_South_M","t2_South_A","t2_South_A",
               "t2_PRH_A","t2_PRH_M","t2_North_A","t2_North_A","t2_North_A","t2_North_A"),12, 1)
Ut <- array(U2, dim = c(dim(U1), dim(dat)[2]))
TT <- dim(dat)[2]
Ut[, , 1:ceiling(TT / 2)] <- U1

#Model 1B - 12 individual pops 6 sites Adult and molt have same U -------------------
# years = 1982-2022
# RT 48 min
t0 <- Sys.time()

Z.model=factor(c(1:12))  # 6 sites x 2 age classes
#R.model <- diag(0.05, 12) #known observation error variance
R.model="diagonal and equal"  #zero since focusing on B  else  "diagonal and equal"
U.model=Ut  #Select from above
Q.model="equalvarcov"#"diagonal and equal"   #near zero add control
B.model="identity"#unconstrained"  # > 90 min

# x0.model=matrix(c(3.6, 4.4,# 1.7,  #BL  # estimated from first data point in time series at each site_age combo
#                   5.5, 6, # 3,       #DE
#                   5.3, 6, # 4.5,    #DP
#                   #2.5, 2.5,# 2,  # DR
#                   #3, 3, # 2,      # PB
#                   3, 3, # 2,      #PRH
#                   5, 5.5, # 4.5,   #TB
#                   6.1, 6.4), # 4.8), #TP
#                 nrow = 12, ncol = 1,
#                 byrow = TRUE)

C.model=matrix(              list("BL_A",0,0,0,0,0,"PDO_A",
                                  "BL_M",0,0,0,0,0,"PDO_M",
                                  #  "BL_P",0,0,0,0," PDO_P",
                                  
                                  0,"DE_A",0,0,0,0,"PDO_A",
                                  0,"DE_M",0,0,0,0,"PDO_M",
                                  #   0,"DE_P",0,0,0," PDO_P",
                                  
                                  0,0,"DP_A",0,0,0,"PDO_A",
                                  0,0,"DP_M",0,0,0,"PDO_M",
                                  #   0,0,"DP_P",0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",
                                  # 0,0,0,0,0,0,0,0,"PDO_M",       #DR
                                  #   0,0,0,0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",       #PB
                                  # 0,0,0,0,0,0,0,0,"PDO_M",
                                  #   0,0,0,0,0," PDO_P"
                                  
                                  0,0,0,0,0,0,"PDO_A",         #PRH
                                  0,0,0,0,0,0,"PDO_M",
                                  
                                  0,0,0,0,0,0,"PDO_A",         #TB
                                  0,0,0,0,0,0,"PDO_M",
                                  
                                  0,0,0,0,0,0,"PDO_A",         #TP
                                  0,0,0,0,0,0,"PDO_M"
),
nrow = 12, ncol = 7,
byrow = TRUE)
m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_3Site=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                                   C = C.model,
                                                                   #x0 = x0.model, 
                                                                   tinitx=1, 
                                                                   c = small_c_Coyote_01_PDO_MAR),

                                                   control=list(maxit=5000, safe=TRUE, 
                                                   control=list(maxit=7000, safe=TRUE, 
                                                                trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "Model_1_K - m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_3Site", aic = m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_3Site$AIC)
df_aic
beepr::beep()
save(m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_3Site, file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_3Site.Rdata")
#load(file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_3Site.RData")

#end time
t1 <- Sys.time()
#Runtime
t1-t0





#timevarying u by 3 sites Adult and Molt combined
# causes R to crash!!!! 
U1 <- matrix(c("t1_BL","t1_BL","t1_South","t1_South","t1_South","t1_South",
               "t1_PRH","t1_PRH","t1_North","t1_North","t1_North","t1_North"),12, 1)
U2 <- matrix(c("t2_BL","t2_BL","t2_South","t2_South","t2_South","t2_South",
               "t2_PRH","t2_PRH","t2_North","t2_North","t2_North","t2_North"),12, 1)
Ut <- array(U2, dim = c(dim(U1), dim(dat)[2]))
TT <- dim(dat)[2]
Ut[, , 1:ceiling(TT / 2)] <- U1


#Model 1B - 12 individual pops 6 sites Adult and molt have same U -------------------
# years = 1982-2022
# RT 25 min
t0 <- Sys.time()

Z.model=factor(c(1:12))  # 6 sites x 2 age classes

#R.model <- diag(0.05, 12) #known observation error variance
R.model="diagonal and equal"  #zero since focusing on B  else  "diagonal and equal"
U.model=Ut  #Select from above
Q.model="equalvarcov"#"diagonal and equal"   #near zero add control
B.model="identity"#unconstrained"  # > 90 min

R.model="diagonal and equal"  #zero since focusing on B  else  "diagonal and equal"
U.model=Ut  #Select from above
Q.model="diagonal and equal"   #near zero add control
B.model="unconstrained"  # > 90 min

# x0.model=matrix(c(3.6, 4.4,# 1.7,  #BL  # estimated from first data point in time series at each site_age combo
#                   5.5, 6, # 3,       #DE
#                   5.3, 6, # 4.5,    #DP
#                   #2.5, 2.5,# 2,  # DR
#                   #3, 3, # 2,      # PB
#                   3, 3, # 2,      #PRH
#                   5, 5.5, # 4.5,   #TB
#                   6.1, 6.4), # 4.8), #TP
#                 nrow = 12, ncol = 1,
#                 byrow = TRUE)

C.model=matrix(              list("BL_A",0,0,0,0,0,"PDO_A",
                                  "BL_M",0,0,0,0,0,"PDO_M",
                                  #  "BL_P",0,0,0,0," PDO_P",
                                  
                                  0,"DE_A",0,0,0,0,"PDO_A",
                                  0,"DE_M",0,0,0,0,"PDO_M",
                                  #   0,"DE_P",0,0,0," PDO_P",
                                  
                                  0,0,"DP_A",0,0,0,"PDO_A",
                                  0,0,"DP_M",0,0,0,"PDO_M",
                                  #   0,0,"DP_P",0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",
                                  # 0,0,0,0,0,0,0,0,"PDO_M",       #DR
                                  #   0,0,0,0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",       #PB
                                  # 0,0,0,0,0,0,0,0,"PDO_M",
                                  #   0,0,0,0,0," PDO_P"
                                  
                                  0,0,0,0,0,0,"PDO_A",         #PRH
                                  0,0,0,0,0,0,"PDO_M",
                                  
                                  0,0,0,0,0,0,"PDO_A",         #TB
                                  0,0,0,0,0,0,"PDO_M",
                                  
                                  0,0,0,0,0,0,"PDO_A",         #TP
                                  0,0,0,0,0,0,"PDO_M"
),
nrow = 12, ncol = 7,
byrow = TRUE)
m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_3Site_1U=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                                    C = C.model,
                                                                    #x0 = x0.model, 
                                                                    tinitx=1, 
                                                                    c = small_c_Coyote_01_PDO_MAR),
                                                    control=list(maxit=7000, safe=TRUE, 
                                                                 trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "Model_1_I - m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_3Site_1U", aic = m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_3Site_1U$AIC)
df_aic
beepr::beep()
save(m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_3Site_1U, file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_3Site_1U.Rdata")
#load(file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_3Site_1U.RData")

#end time
t1 <- Sys.time()
#Runtime
t1-t0


###--------------------
#2024-02-12

#one U changing through time
# but can't do covariates by site
#timevarying u by site
U1 <- matrix(c("U1"),1, 1)
U2 <- matrix(c("U2"),1, 1)
Ut <- array(U2, dim = c(dim(U1), dim(dat)[2]))
TT <- dim(dat)[2]
Ut[, , 1:ceiling(TT / 2)] <- U1

## TIME VARYING MODEL by site
t0 <- Sys.time()

#Z.model=factor(c(1:12))  # 6 sites x 2 age classes
Z.model=factor(rep(1, 12))
R.model="diagonal and equal"
U.model=Ut  #Select from above
Q.model="equalvarcov"#"diagonal and equal"   #near zero add control
B.model="unconstrained"
C.model = matrix(list("PDO"),
                 nrow = 1,
                 ncol = 1,
                 byrow = TRUE)
m.Ind_Molt_Adult_Coyote_PDO_B_equalcov_1U_TV_1pop=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                                        C = C.model,
                                                                        #x0 = x0.model, 
                                                                        tinitx=1, 
                                                                        c = small_c_PDO_MAR),
                                                        control=list(maxit=5000, safe=TRUE, 
                                                                     trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "Model_1_J - m.Ind_Molt_Adult_Coyote_PDO_B_equalcov_1U_TV_1pop", aic = m.Ind_Molt_Adult_Coyote_PDO_B_equalcov_1U_TV_1pop$AIC)
df_aic
beepr::beep()
save(m.Ind_Molt_Adult_Coyote_PDO_B_equalcov_1U_TV_Site, file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_equalcov_1U_TV_1pop.RData")
#load(file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_equalcov_1U_TV_1pop.RData")

#end time
t1 <- Sys.time()
#Runtime
t1-t0

# END H1 ------------------------------------------------------------   

########### 2024-05-26
########### using  UI and UI lag
## updated cumulative Coyote 3 yr
## TV at 2003

          
#timevarying u
# all sites have same growth rates that can change around 2003

U1 <- matrix("t1", 12, 1)
U2 <- matrix("t2", 12, 1)
Ut <- array(U2, dim = c(dim(U1), dim(dat)[2]))
TT <- dim(dat)[2]
Ut[, , 1:ceiling(TT / 2)] <- U1


## TIME VARYING MODEL same all

#Model 1B - 12 individual pops 6 sites Adult and molt have same U -------------------
# years = 1982-2022
# RT 34 min
t0 <- Sys.time()

Z.model=factor(c(1:12))  # 6 sites x 2 age classes

R.model <- diag(0.01, 12) #known observation error variance

# R.model <- diag(0.1, 12) #known observation error variance

#R.model="equal"  #zero since focusing on B  else  "diagonal and equal"
U.model=Ut  #Select from above
Q.model="equalvarcov" #"diagonal and equal"   #near zero add control
B.model="identity" # unconstrained"  > 90 min

# x0.model=matrix(c(3.6, 4.4,# 1.7,  #BL  # estimated from first data point in time series at each site_age combo
#                   5.5, 6, # 3,       #DE
#                   5.3, 6, # 4.5,    #DP
#                   #2.5, 2.5,# 2,  # DR
#                   #3, 3, # 2,      # PB
#                   3, 3, # 2,      #PRH
#                   5, 5.5, # 4.5,   #TB
#                   6.1, 6.4), # 4.8), #TP
#                 nrow = 12, ncol = 1,
#                 byrow = TRUE)


C.model=matrix(              list("BL_A",0,0,0,0,0,"UI_A","UI_lag_A",  #adding UI and UI_lag
                                  "BL_M",0,0,0,0,0,"UI_M","UI_lag_M",
                                  #  "BL_P",0,0,0,0," PDO_P",
                                  
                                  0,"DE_A",0,0,0,0,"UI_A","UI_lag_A",
                                  0,"DE_M",0,0,0,0,"UI_M","UI_lag_M",
                                  #   0,"DE_P",0,0,0," PDO_P",
                                  
                                  0,0,"DP_A",0,0,0,"UI_A","UI_lag_A",
                                  0,0,"DP_M",0,0,0,"UI_M","UI_lag_M",
                                  #   0,0,"DP_P",0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",
                                  # 0,0,0,0,0,0,0,0,"PDO_M",       #DR
                                  #   0,0,0,0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",       #PB
                                  # 0,0,0,0,0,0,0,0,"PDO_M",
                                  #   0,0,0,0,0," PDO_P"
                                  
                                  0,0,0,0,0,0,"UI_A", "UI_lag_A",        #PRH
                                  0,0,0,0,0,0,"UI_M", "UI_lag_M",
                                  
                                  0,0,0,0,0,0,"UI_A", "UI_lag_A",       #TB
                                  0,0,0,0,0,0,"UI_M", "UI_lag_M",
                                  
                                  0,0,0,0,0,0,"UI_A", "UI_lag_A",       #TP
                                  0,0,0,0,0,0,"UI_M", "UI_lag_M"
),
nrow = 12, ncol = 8,
byrow = TRUE)
m.Ind_Molt_Adult_Coyote_UI_B_ident_1U_TV_R_fix=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, 
                                                                    R=R.model, B=B.model,
                                                                    C = C.model,
                                                                    #x0 = x0.model, 
                                                                    tinitx=1, 
                                                                    c = small_c_Coyote_3yr_UI_UI_lag),
                                                    control=list(maxit=5000, safe=TRUE, 
                                                                 trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "Model_1_N - m.Ind_Molt_Adult_Coyote_UI_B_ident_1U_TV_R_fix", 
                             aic =m.Ind_Molt_Adult_Coyote_UI_B_ident_1U_TV_R_fix$AIC)
df_aic
beepr::beep()

save(m.Ind_Molt_Adult_Coyote_UI_B_ident_1U_TV_R_fix, file = "Output/m.Ind_Molt_Adult_Coyote_UI_B_ident_1U_TV_R_fix.Rdata")

#load(file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_ident_1U_TV_R_fix.RData")

#end time
t1 <- Sys.time()
#Runtime
t1-t0



########### 2024-05-26
########### using  UI and UI lag
## updated cumulative Coyote 3 yr
## TV at 2003


#timevarying u
# all sites have same growth rates that can change around 2003

U1 <- matrix("t1", 12, 1)
U2 <- matrix("t2", 12, 1)
Ut <- array(U2, dim = c(dim(U1), dim(dat)[2]))
TT <- dim(dat)[2]
Ut[, , 1:ceiling(TT / 2)] <- U1


## TIME VARYING MODEL same all

#Model 1B - 12 individual pops 6 sites Adult and molt have same U -------------------
# years = 1982-2022
# RT 34 min
t0 <- Sys.time()

Z.model=factor(c(1:12))  # 6 sites x 2 age classes

R.model <-  "diagonal and unequal" #known observation error variance

# R.model <- diag(0.01, 12) #known observation error variance

#R.model="equal"  #zero since focusing on B  else  "diagonal and equal"
U.model=Ut  #Select from above
Q.model="diagonal and equal" #"diagonal and equal"   #near zero add control
#B.model="diagonal and equal" # unconstrained"  > 90 min
#b.model give "AA" for same site correlation,
#               "BB" for molt-adult correlation
# other sites show potential relationships
# hypothesized "zeros" get an "XX"
B.model=matrix(list("AA","BB",    "A",0,   "A",0,   0,0,   0,0,  0,0,  #BL
                    "BB","AA",     0,"A",   0,"A",  0,0,   0,0,  0,0,
                    
                    "A",0,  "AA","BB",   "A",0,  "A",0,    0,0,  0,0, #DE
                    0,"A",  "BB","AA",   0,"A",  0,"A",    0,0,  0,0,
                    
                    "A",0,  "A",0,   "AA","BB",  "A",0,   0,0,  0,0,  #DP
                    0,"A",  0,"A",   "BB","AA",  0,"A",   0,0,  0,0,
                    
                    "B",0,  "B",0, "B",0,    "AA","BB",   "C",0,  "C",0,  #PRH
                    0,"B",  0,"B",   0,"B",  "BB","AA",   0,"C",  0,"C",
                    
                    0,0,   0,0,     0,0,   "D",0,   "AA","BB",   "D",0,  #TB
                    0,0,   0,0,     0,0,   0,"D",   "BB","AA",   0,"D",
                    
                    0,0,   0,0,     0,0,   0,0,    "D",0,   "AA","BB",  #TP
                    0,0,   0,0,     0,0,   0,0,    0,"D",   "BB","AA"),
               
               nrow = 12, ncol = 12,
               byrow = TRUE)

# x0.model=matrix(c(3.6, 4.4,# 1.7,  #BL  # estimated from first data point in time series at each site_age combo
#                   5.5, 6, # 3,       #DE
#                   5.3, 6, # 4.5,    #DP
#                   #2.5, 2.5,# 2,  # DR
#                   #3, 3, # 2,      # PB
#                   3, 3, # 2,      #PRH
#                   5, 5.5, # 4.5,   #TB
#                   6.1, 6.4), # 4.8), #TP
#                 nrow = 12, ncol = 1,
#                 byrow = TRUE)


C.model=matrix(              list("BL_A",0,0,0,0,0,"UI_A","UI_lag_A",  #adding UI and UI_lag
                                  "BL_M",0,0,0,0,0,"UI_M","UI_lag_M",
                                  #  "BL_P",0,0,0,0," PDO_P",
                                  
                                  0,"DE_A",0,0,0,0,"UI_A","UI_lag_A",
                                  0,"DE_M",0,0,0,0,"UI_M","UI_lag_M",
                                  #   0,"DE_P",0,0,0," PDO_P",
                                  
                                  0,0,"DP_A",0,0,0,"UI_A","UI_lag_A",
                                  0,0,"DP_M",0,0,0,"UI_M","UI_lag_M",
                                  #   0,0,"DP_P",0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",
                                  # 0,0,0,0,0,0,0,0,"PDO_M",       #DR
                                  #   0,0,0,0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",       #PB
                                  # 0,0,0,0,0,0,0,0,"PDO_M",
                                  #   0,0,0,0,0," PDO_P"
                                  
                                  0,0,0,0,0,0,"UI_A", "UI_lag_A",        #PRH
                                  0,0,0,0,0,0,"UI_M", "UI_lag_M",
                                  
                                  0,0,0,0,0,0,"UI_A", "UI_lag_A",       #TB
                                  0,0,0,0,0,0,"UI_M", "UI_lag_M",
                                  
                                  0,0,0,0,0,0,"UI_A", "UI_lag_A",       #TP
                                  0,0,0,0,0,0,"UI_M", "UI_lag_M"
),
nrow = 12, ncol = 8,
byrow = TRUE)
m.Ind_Molt_Adult_Coyote_UI_B_PRH_1U_TV_R_fix=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, 
                                                                     R=R.model, B=B.model,
                                                                     C = C.model,
                                                                     #x0 = x0.model, 
                                                                     tinitx=1, 
                                                                     c = small_c_Coyote_3yr_UI_UI_lag),
                                                     control=list(maxit=5000, safe=TRUE, 
                                                                  trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "Model_1_N - m.Ind_Molt_Adult_Coyote_UI_B_PRH_1U_TV_R_fix", 
                             aic =m.Ind_Molt_Adult_Coyote_UI_B_PRH_1U_TV_R_fix$AIC)
df_aic
beepr::beep()

save(m.Ind_Molt_Adult_Coyote_UI_B_PRH_1U_TV_R_fix, file = "Output/m.Ind_Molt_Adult_Coyote_UI_B_PRH_1U_TV_R_fix.Rdata")

#load(file = "Output/m.Ind_Molt_Adult_Coyote_UI_B_PRH_1U_TV_R_fix.Rdata")

#end time
t1 <- Sys.time()
#Runtime
t1-t0


#############################
# Model plots -------------------------------------
#############################
BESTMODEL <- m.Ind_Molt_Adult_Coyote_UI_B_PRH_1U_TV_R_fix #m.Ind_Molt_Adult_Coyote_PDO_B_equalcov_1U_TV_Site # m.1pop_Molt_Adult_Coyote_PDO_B_custom_1U_R0 # # m.1pop_Molt_Adult_Coyote_PDO_B_custom_1U_R0 # # m.1pop_Molt_Adult_Coyote_PDO_B_custom_1U# m.1pop_Molt_Adult_Coyote_PDO_B_unc_2U #m.1pop_Coyote_PDO_Xo_fixed_B_unc_tinitx_1 #m.1pop_Coyote_PDO_Xo_fixed_B_unc# m.1pop_Coyote_PDO_B_unc #   m.1pop_Coyote_PDO_B_unc   m.5pop_Coyote_PDO_B_unc

#############################
# Model plots -------------------------------------
#############################
BESTMODEL <- m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_R_fix # m.1pop_Molt_Adult_Coyote_PDO_B_custom_1U_R0 # # m.1pop_Molt_Adult_Coyote_PDO_B_custom_1U_R0 # # m.1pop_Molt_Adult_Coyote_PDO_B_custom_1U# m.1pop_Molt_Adult_Coyote_PDO_B_unc_2U #m.1pop_Coyote_PDO_Xo_fixed_B_unc_tinitx_1 #m.1pop_Coyote_PDO_Xo_fixed_B_unc# m.1pop_Coyote_PDO_B_unc #   m.1pop_Coyote_PDO_B_unc   m.5pop_Coyote_PDO_B_unc


autoplot(BESTMODEL)

autoplot(BESTMODEL, plot.type = "fitted.ytT") + # xtT
  ylim(4,8) # xtT

autoplot(BESTMODEL, plot.type = "xtT") + # xtT
  ylim(3.1,7.6)

# plot temporal change from starting values

#BESTMODEL <- m.MOLT_IND_SITE_MOCI.equal_B_diag_uneq

CIs <- MARSSparamCIs(BESTMODEL, alpha = 0.11)  #crashes if in results # start 420

#CIs <- MARSSparamCIs(BESTMODEL, alpha = 0.11, hessian.fun = "fdHess")  #runs if few NAs, about 60 min

CIs

CIs$states.se  # change this to the CIs rather than SEs...  use 89%CIs for everything?



d <- as_tibble(t(BESTMODEL$states-BESTMODEL$states[,1]))
d.se <- as_tibble(t(BESTMODEL$states.se))
#add header names
names(d)[c(1:12)] <- c("BL_Breed", "BL_Molt", #"BL_Pup",
                       "DE_Breed", "DE_Molt", #"DE_Pup",
                       "DP_Breed", "DP_Molt", #"DP_Pup",
                       "PRH_Breed", "PRH_Molt",# "PRH_Pup",
                       "TB_Breed", "TB_Molt", #"TB_Pup",
                       "TP_Breed", "TP_Molt" #"TP_Pup",
                      # "DR_Breed", "DR_Molt", #"DR_Pup",
                      # "PB_Breed", "PB_Molt" #"PB_Pup"
                       )
names(d.se)[c(1:12)] <- c("BL_Breed", "BL_Molt", #"BL_Pup",
                          "DE_Breed", "DE_Molt", #"DE_Pup",
                          "DP_Breed", "DP_Molt", #"DP_Pup",
                          "PRH_Breed", "PRH_Molt",# "PRH_Pup",
                          "TB_Breed", "TB_Molt", #"TB_Pup",
                          "TP_Breed", "TP_Molt" #"TP_Pup",
                         # "DR_Breed", "DR_Molt", #"DR_Pup",
                        #  "PB_Breed", "PB_Molt" #"PB_Pup"
)


d <- as_tibble(cbind(years,d))

d2 <- d %>% pivot_longer(cols = c(2:13), names_to = "Subsite_Season", values_to = "log_est")
d2.se <- d.se %>% pivot_longer(cols = c(1:12), names_to = "Subsite_Season", values_to = "log_est")
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
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = log_est-log_se, ymax = log_est+log_se),  
              alpha = 0.2, colour = NA) +
  geom_hline(yintercept = c(-1,0,1), lty = 2) +
  xlim(1982, 2022) +
  ylim(-1.5, 2) + 
  theme_minimal(base_size = 20) +
  ylab("index of log abundance (N-No)") +
  xlab("Year") +
  facet_wrap(.~Season, ncol = 1)


####----plot overall pop size for each age class summing sites

d.tot <- as_tibble(t(exp(BESTMODEL$states)))  #exp to original scale
d.tot.se <- as_tibble(t(exp(BESTMODEL$states.se))) #exp to original scale
#add header names

names(d.tot)[c(1:12)] <- c("BL_Breed", "BL_Molt", #"BL_Pup",
                       "DE_Breed", "DE_Molt", #"DE_Pup",
                       "DP_Breed", "DP_Molt", #"DP_Pup",
                       "PRH_Breed", "PRH_Molt",# "PRH_Pup",
                       "TB_Breed", "TB_Molt", #"TB_Pup",
                       "TP_Breed", "TP_Molt" #"TP_Pup",
                       # "DR_Breed", "DR_Molt", #"DR_Pup",
                       # "PB_Breed", "PB_Molt" #"PB_Pup"
)
names(d.tot.se)[c(1:12)] <- c("BL_Breed", "BL_Molt", #"BL_Pup",
                          "DE_Breed", "DE_Molt", #"DE_Pup",
                          "DP_Breed", "DP_Molt", #"DP_Pup",
                          "PRH_Breed", "PRH_Molt",# "PRH_Pup",
                          "TB_Breed", "TB_Molt", #"TB_Pup",
                          "TP_Breed", "TP_Molt" #"TP_Pup",
                          # "DR_Breed", "DR_Molt", #"DR_Pup",
                          #  "PB_Breed", "PB_Molt" #"PB_Pup"
)

#sum all molt and breeding sites by year
Breed.tot <- d.tot %>%
  dplyr::reframe(Breed = rowSums(across(c(BL_Breed, DE_Breed, DP_Breed, PRH_Breed, TB_Breed, TP_Breed))))
Molt.tot <- d.tot %>%
  dplyr::reframe(Molt = rowSums(across(c(BL_Molt, DE_Molt, DP_Molt, PRH_Molt, TB_Molt, TP_Molt))))

Breed.tot.se <- d.tot.se %>%
  dplyr::reframe(Breed.se = rowSums(across(c(BL_Breed, DE_Breed, DP_Breed, PRH_Breed, TB_Breed, TP_Breed))))
Molt.tot.se <- d.tot.se %>%
  dplyr::reframe(Molt.se = rowSums(across(c(BL_Molt, DE_Molt, DP_Molt, PRH_Molt, TB_Molt, TP_Molt))))


#put in a single table
d.Breed.Molt.tot <- as_tibble(cbind(years,Breed.tot, Molt.tot))


d.Breed.Molt.tot <- d.Breed.Molt.tot %>% pivot_longer(cols = c(2:3), names_to = "Season", values_to = "estimate")

d.Breed.Molt.se <- as_tibble(cbind(years,Breed.tot.se, Molt.tot.se))



d.Breed.Molt.se <- d.Breed.Molt.se %>% pivot_longer(cols = c(2:3), names_to = "Season", values_to = "SE")
d.Breed.Molt.se <- d.Breed.Molt.se[,3] #remove season and year

d4 <- as_tibble(cbind(d.Breed.Molt.tot,d.Breed.Molt.se))

#add groupings for site and molt

ggplot(d4, aes(x = years, y = estimate, color = Season)) +
  #geom_line(aes(linetype = Season), size = 1.25) +
  geom_point(size = 2) + 
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = estimate-SE, ymax = estimate+SE),  
              alpha = 0.2, colour = NA) +
  #geom_hline(yintercept = c(-1,0,1), lty = 2) +
  xlim(1982, 2022) +
  ylim(0, 5000) + 
  theme_minimal(base_size = 20) +
  ylab("Estimated abundance") +
  xlab("Year") 








####-----plot covariate effects

coef.data <- tidy(CIs)

#get just the coefficient terms
coef.data <- coef.data %>%
  filter(str_detect(term, "^C."))

coef.data$Season <- ifelse(str_detect(coef.data$term, "_A"), "Breeding Adults", 
                           ifelse(str_detect(coef.data$term, "_M"), "Molting", "Pups"))
#coef.data$Covariate <- ifelse(str_detect(coef.data$term, "MOCI"), "MOCI", "Coyote")
coef.data$Covariate <- ifelse(str_detect(coef.data$term, "UI"), "UI", "Coyote")
  
ggplot(coef.data, aes(term, estimate, color = Season, shape = Season)) +

  geom_pointrange(aes(ymin = conf.low, ymax = conf.up), size = 0.8) +
  coord_flip() +
  theme(legend.position = c(0.2, 0.2)) + 
  #ylim(-1.25, 0.25) +

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
(coef(BESTMODEL, type="matrix")$U)  # growth | su

exp(0.2389865)
exp(0.205)

B <- coef(BESTMODEL, type="matrix")$B  # effect of column on row

library("corrplot")

#Scale to get to -1 to +1
max(B)
min(B)
B <- B/abs(min(B))

#or scale with floor ceiling.

#B <- ifelse(B > 3, 3, B)
#B <- ifelse(B < -3, -3, B)



dimnames(B) <- list(c("BL Adult", "BL Molt", #"BL Pup",
                       "DE Adult", "DE Molt", #"DE Pup",
                       "DP Adult", "DP Molt", #"DP Pup",
                      #"DR Adult", "DR Molt", #"DP Pup",
                      #"DP Adult", "PB Molt", #"DP Pup",
                      "PRH Adult", "PRH Molt",
                       "TB Adult", "TB Molt", #"TB Pup",
                       "TP Adult", "TP Molt"), #"TP Pup"),
                      c("BL Adult", "BL Molt", #"BL Pup",
                        "DE Adult", "DE Molt", #"DE Pup",
                        "DP Adult", "DP Molt", #"DP Pup",
                       # "DR Adult", "DR Molt", #"DP Pup",
                        #"DP Adult", "PB Molt", #"DP Pup",
                        "PRH Adult", "PRH Molt",
                        "TB Adult", "TB Molt", #"TB Pup",
                        "TP Adult", "TP Molt")) #"TP Pup"),

par(mfrow = c(1,1))
corrplot(B, method="color", is.corr = FALSE)

## Process correlation between sites --------
#recover q.matrix p 115-116
# use m.MOLT_IND_SITE_B_diag_uneq_q_unc to show all correlations for Q



Q.unc <- coef(BESTMODEL, type="matrix")$Q
h <- diag(1 / sqrt(diag(Q.unc)))
Q.corr <- h %*% Q.unc %*% h

## get Z.model from appropriate model
Z.model=factor(1:12) 

rownames(Q.corr) <- unique(Z.model)
colnames(Q.corr) <- unique(Z.model)
# Process Correlations between sites
Q.corr

### PREDICT -----------

##foreward 5 years
# must add c data
# lets say coyote at DE and DP and warm PDO
#  
# example A
c_forecast=matrix(c(0,0,0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,0, #1,1,1,1,1,  Contrast with Y/N coyote
                    0,0,0,0,0,0,0,0,0,0, #1,1,1,1,1,
             0,0,0,0,0,0,0,0,0,0, #1,1,1,1,1,
             0,0,0,0,0,0,0,0,0,0, #1,1,1,1,1,
             0,0,0,0,0,0,0,0,0,0, #1,1,1,1,1,
             1,1,1,1,1,1,1,1,1,1), #PDO
nrow = 7, ncol = 10,
byrow = TRUE)

# example B
c_forecast=matrix(c(0,0,0,0,0,1,1,1,1,1,  #BL
                    1,1,1,1,1,1,1,1,1,1,  #DE
                    0,1,1,0,1,1,0,0,1,1,  #DP
                    0,0,0,0,0,0,0,0,0,0,  #PRH
                    1,1,1,1,1,1,1,1,1,1,  #TB
                    0,0,0,0,0,1,1,1,1,1,  #TP 
                    2,1,1,1,0,0,-1,1,2,1), #PDO
                  nrow = 7, ncol = 10,
                  byrow = TRUE)


c_new <- cbind(small_c_Coyote_01_PDO_MAR, c_forecast)

forecast1 <- predict(BESTMODEL, n.ahead = 10, interval = "prediction", 
                     nsim = 100,
                     newdata = list(c = c_forecast))


forecast1_plot_data <- forecast1$pred

d.1 <- forecast1_plot_data[,1:2]
d.2 <- exp(forecast1_plot_data[,3:9])
forecast1_plot_data <- tibble(d.1, d.2)

names(forecast1_plot_data)

ggplot(forecast1_plot_data, aes(t+1981, estimate)) +
  geom_line() +
  geom_ribbon(aes(ymin = `Lo 80`, ymax = `Hi 80`), alpha =0.2) + 
  geom_point(aes(t+1981, y), color = "blue4") + 
  geom_vline(xintercept = 2022.5, linetype = 2, color = "blue4") + 
  geom_vline(xintercept = 2003.5, linetype = 2) + #show TV timepoint
  xlab("Year") + 
  ylab("Seals") +
  facet_wrap(.~.rownames)


autoplot(forecast1) +
  theme_grey(base_size = 18) +
  #ylim(3,7.6) #+
  scale_y_continuous(trans="exp")
  

