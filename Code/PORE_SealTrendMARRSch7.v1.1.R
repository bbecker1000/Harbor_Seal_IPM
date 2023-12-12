##edited code from MARSS UserGuide chapter 7

library("MARSS")
library(ggplot2)
library("readxl")
library(dplyr)
library(tidyverse)

## B in
#  terms of the interaction strengths between species; bij
#  equals dfi /dXj, the change in the log population growth
#  rate of species i with respect to changes in the log
#  population abundance of species j. T

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
Phoca$Age <- ifelse(Phoca$Julian > 135, "MOLT", Phoca$Age)


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

all_data$Season <- ifelse(all_data$Julian <= 135 | !is.na(all_data$Season), "PUPPING", "MOLT")
all_data

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
  filter(Age != "PUP" | Season != "MOLT") #remove pups counted in early molting season
top1_all_data

top.plot <- ggplot(top1_all_data, aes(Year, Count)) + 
  geom_point() + 
  geom_smooth() + 
  geom_vline(xintercept = c(1988, 1994), lty = 3) +
  scale_y_continuous(trans='log10') 
top.plot + facet_grid(Age ~ Subsite)





## now within year plots to look at breeding and molting season peaks

#Phoca.Adult <- dplyr::filter(Phoca, Age == "ADULT")

seasonal.plot <- plot.breeding <- ggplot(all_data, aes(Julian, Count, colour = Year)) + 
  geom_point() + 
  geom_smooth(aes(colour = Year))
seasonal.plot + facet_grid(Age ~ Subsite)

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


## now need to make the data wide for MARSS
all_data.MARSS.wide <- all_data.MARSS %>% 
  pivot_wider(names_from = Subsite_Age, values_from = Count)




dat2 <- all_data.MARSS.wide ## for use in SealPopStructure Script

## now use names in the basic MARSS code 
dat <- t(all_data.MARSS.wide)
## OR FOR PUPS
## dat <- t(top1.pup.breed.spread)


############# env data --------------------------------

# MEI file has several ocean indices and coyote data
MEI <- read_excel("Data/MEI.xlsx")
show(MEI)
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
plot(MEI$MEIv2_DEC_JAN, MEI$NPGO_MAR)
cor.test(MEI$PDO_MAR, as.numeric(MEI$MOCI_JFM_NC))
# cor = 0.79 (0.6-0.9) P < 0.001
# So let's use PDO_MAR as our env covariate.'
plot(MEI$PDO_MAR, MEI$NPGO_MAR  )
plot(MEI$MEIv2_DEC_JAN, MEI$PDO_MAR  )

# MARSS recommends start with more covariates and go smaller.
# create covariates for 
# - coyote Y/N + PDO_MAR
# - coyote count + PDO_MAR
# - PDO_MAR
# - whichever Coyote measure does better
# - adding the anthropogenic positive and negatives?


# no coyote sightings at TB or TP so remove
Coyote_01 <- MEI[,c(7:11)]  
small_c_Coyote_01 <- as.matrix(t(Coyote_01))

Coyote_Rate <- MEI[,c(12:16)] 
small_c_Coyote_Rate <- as.matrix(t(Coyote_Rate))







##########################################################################
### 8.3 Single population model with independent and (non)identical errors
##########################################################################
#Code to fit the single population model with i.i.d. errors
#Read in data
#Transpose since MARSS needs time ACROSS columns
years = dat[1,]     # remove years
n = nrow(dat)-1
dat = dat[2:nrow(dat),]
legendnames = (unlist(dimnames(dat)[1]))

########################################################################

## Parameter key ########
## R = observation errors          equal
## U = growth parameter            unequal  
## Z = design matrix 
## Q = hidden state process        diagonal and equal
## B = effect of column on row     unequal (these are the interactions)

#####################################################################

df_aic <- data.frame(model=character(), aic=integer())

###NOW TO MODELS ##########
#estimate parameters
Z.model = factor(c(rep(1, 15)))
R.model = "diagonal and unequal" #error variance not the same
kem_onepop = MARSS(dat, model=list(Z=Z.model, R=R.model))

#make figure
par(mfrow=c(1,1))
matplot(years, t(dat),xlab="",ylab="index of log abundance",
        pch=c("1","2","3","4","5"), ylim=c(3,9), bty="L")
lines(years,kem_onepop$states-1.96*kem_onepop$states.se,type="l",
      lwd=1,lty=2,col="red")
lines(years,kem_onepop$states+1.96*kem_onepop$states.se,type="l",
      lwd=1,lty=2,col="red")
lines(years,kem_onepop$states,type="l",lwd=2)
title("Observations and total population estimate",cex.main=.9)

coef(kem_onepop, type="vector")  #show the estimated parameter elements as a vector
coef(kem_onepop, type="matrix")$R
#show estimated elements for each parameter matrix as a list
#coef(kem_onepop) 

kem_onepop$logLik   #show the log-likelihood
kem_onepop$AIC  #show the AIC

df_aic <- df_aic %>% add_row(model = "one pop - r unequal", aic = kem_onepop$AIC)

#plot residuals, from 8.3 (page 99)
plotdat = t(dat)
matrix.of.biases = matrix(coef(kem_onepop, type="matrix")$A,
                          nrow=nrow(plotdat),ncol=ncol(plotdat),byrow=T)
xs = matrix(kem_onepop$states,
            nrow=dim(plotdat)[1],ncol=dim(plotdat)[2],byrow=F)
resids = plotdat-matrix.of.biases-xs
par(mfrow=c(3, 5))
for(i in 1:n){
  plot(resids[!is.na(resids[,i]),i],ylab="residuals")
  title(paste("One Population", legendnames[i]))
}
#bad residuals: BL, DP
#good-ish residuals: DE, TB, TP(?)


##########################################################
### 8.2 A single well-mixed population with i.i.d. errors
##########################################################
#Code to fit the single population model with independent and equal errors
#(silas) model is linear if R. model = "diagonal and equal" changed to unequal 6/6


Z.model = factor(c(rep(1, 15)))
R.model = "diagonal and equal" 
kem2 = MARSS(dat, model=list(Z=Z.model, R=R.model))

coef(kem2) #the estimated parameter elements
kem2$logLik #log likelihood
kem2$AIC  #AICs

df_aic <- df_aic %>% add_row(model = "one pop - r equal", aic = kem2$AIC)


#plot residuals from page 99
plotdat = t(dat)
matrix.of.biases = matrix(coef(kem2, type="matrix")$A,
                          nrow=nrow(plotdat),ncol=ncol(plotdat),byrow=T)
xs = matrix(kem2$states,
            nrow=dim(plotdat)[1],ncol=dim(plotdat)[2],byrow=F)
resids = plotdat-matrix.of.biases-xs
par(mfrow=c(3,5))
for(i in 1:n){
  plot(resids[!is.na(resids[,i]),i],ylab="residuals")
  title(paste("One Population", legendnames[i]))
}
 ## BL and DP residuals look a bit problematic
 ## (Silas) I think most look somewhat problematic, maybe TB and TP are ok

par(mfrow=c(1,1))

###################################################
### code chunk number 23: Cs2_Code4
### Each population is independent

#(silas) BL and DE clearly nonlinear but DP, TB, and TP aren't
#(silas) is this because of the convergence problems??
#(Silas) tried running with maxit=10000 and it didn't change the plots (except the confidence estimates for TB??)


###
# 5 independent populations with unequal observation variances, unequal growth, equal hidden process variance, not too sure about B being identity
###
Z.model=factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5))
R.model="diagonal and equal"
U.model="unequal"
Q.model="diagonal and equal"
B.model="diagonal and equal"
kem4_ind1=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model),
                    control=list(maxit=1000, safe=TRUE)) 
kem4_ind1$AIC 

df_aic <- df_aic %>% add_row(model = "all ind 1", aic = kem4_ind1$AIC)

#same as above but now B.Model is "unequal instead of identity"
#DOESN'T REACH CONVERGENCE
Z.model=factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5))
R.model="diagonal and unequal"
U.model="unequal"
Q.model="diagonal and equal"
B.model="diagonal and equal"#"unequal"
kem4_ind2=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model),
               control=list(maxit=5000, safe=TRUE)) 
beepr::beep(0)
kem4_ind2$AIC 

df_aic <- df_aic %>% add_row(model = "Site ind 2", aic = kem4_ind2$AIC)


#by site, unequal R, U, and B
#RT about 5 min
Z.model=factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5))
R.model="diagonal and unequal"
U.model="unequal"
Q.model="diagonal and equal"
B.model="diagonal and unequal"#"unequal"   "LOWEST AIC SO FAR"
kem4_ind3=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model),
                control=list(maxit=5000, safe=TRUE)) 
beepr::beep(0)
kem4_ind3$AIC 

df_aic <- df_aic %>% add_row(model = "Site ind 3", aic = kem4_ind3$AIC)


plotdat = t(dat)
matrix.of.biases = matrix(coef(kem4_ind3, type="matrix")$A,
                          nrow=nrow(plotdat),ncol=ncol(plotdat),byrow=T)
xs = matrix(kem4_ind3$states,
            nrow=dim(plotdat)[1],ncol=dim(plotdat)[2],byrow=F)
resids = plotdat-matrix.of.biases-xs
par(mfrow=c(3, 5))
for(i in 1:n){
  plot(resids[!is.na(resids[,i]),i],ylab="residuals")
  title(paste("One Population", legendnames[i]))
}


#Molt not related to breeding season
#crashing R
Z.model=factor(c(1,2,1,1,2,1,1,2,1,1,2,1,1,2,1))
R.model="diagonal and unequal"
U.model="unequal"
Q.model="diagonal and equal"
B.model="diagonal and unequal"#"unequal"
kem4_MoltInd=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model),
                control=list(maxit=5000, safe=TRUE)) 
beepr::beep(0)
kem4_MoltInd$AIC 

df_aic <- df_aic %>% add_row(model = "Molt ind", aic = kem4_MoltInd$AIC)


# 3 sites
Z.model=factor(c(1,1,1,2,2,2,2,2,2,3,3,3,3,3,3))
R.model="diagonal and unequal"
U.model="unequal"
Q.model="diagonal and equal"
B.model="diagonal and equal"#"unequal"
kem4_3_sites_unequal=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model),
                control=list(maxit=5000, safe=TRUE)) 
beepr::beep(0)
kem4_3_sites_unequal$AIC 

df_aic <- df_aic %>% add_row(model = "3 Sites_unequal", aic = kem4_3_sites$AIC)
df_aic


## IND + NPGO
Z.model=factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5))
R.model="diagonal and unequal"
U.model="unequal"
Q.model="diagonal and equal"
B.model="diagonal and equal"#"unequal"
kem4_ind2_NPGO=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                c = small_c),
                control=list(maxit=5000, safe=TRUE)) 
beepr::beep(0)
kem4_ind2_NPGO$AIC 

df_aic <- df_aic %>% add_row(model = "Site ind 2 + NPGO", aic = kem4_ind2_NPGO$AIC)


#by site, unequal R, U, and B
#RT 10+ min
Z.model=factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5))
R.model="diagonal and unequal"
U.model="unequal"
Q.model="diagonal and equal"
B.model="diagonal and unequal"#"unequal"   "LOWEST AIC SO FAR"
kem4_ind3_NPGO=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                     c = small_c),
                control=list(maxit=5000, safe=TRUE)) 
beepr::beep(0)
kem4_ind3_NPGO$AIC 

df_aic <- df_aic %>% add_row(model = "Site ind 3 + NPGO", aic = kem4_ind3_NPGO$AIC)
df_aic


#by site, unequal R, U, and B and Coyote
#RT 10+ min
#works, but not matching covariates to sites correctly.
Z.model=factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5))
R.model="diagonal and unequal"
U.model="unequal"
Q.model="diagonal and equal"
B.model="diagonal and unequal"
# if Z is 1:15
# C.model=matrix(c("BL_1","OTH","OTH","OTH","OTH",
#                  "BL_1","OTH","OTH","OTH","OTH",
#                  "BL_1","OTH","OTH","OTH","OTH",
#                  
#                  "OTH","DE_1","OTH","OTH","OTH",
#                  "OTH","DE_1","OTH","OTH","OTH",
#                  "OTH","DE_1","OTH","OTH","OTH",
#                  
#                  "OTH","OTH","DP_1","OTH","OTH",
#                  "OTH","OTH","DP_1","OTH","OTH",
#                  "OTH","OTH","DP_1","OTH","OTH",
#                  
#                  "OTH","OTH","OTH","OTH","OTH",
#                  "OTH","OTH","OTH","OTH","OTH",
#                  "OTH","OTH","OTH","OTH","OTH",
#                  
#                  "OTH","OTH","OTH","OTH","OTH",
#                  "OTH","OTH","OTH","OTH","OTH",
#                  "OTH","OTH","OTH","OTH","OTH"),
#               nrow = 15, ncol = 5,
#               byrow = TRUE)

# if Z is rep(1:5), times = 3
# 6 min RT
C.model=matrix(c("BL_1","OTH","OTH","OTH","OTH",
                 "OTH","DE_1","OTH","OTH","OTH",
                 "OTH","OTH","DP_1","OTH","OTH",
                 "OTH","OTH","OTH","OTH","OTH",
                 "OTH","OTH","OTH","OTH","OTH"),
               nrow = 5, ncol = 5,
               byrow = TRUE)

kem4_ind3_Coyote=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                       C = C.model,
                                       c = small_c_Coyote_01),
                     control=list(maxit=5000, safe=TRUE, trace = 1)) 
beepr::beep(0)
kem4_ind3_Coyote$AIC 

plot(kem4_ind3_Coyote, plot.type="model.resids.ytT")





CIs <- MARSSparamCIs(kem4_ind3_Coyote)
CIs

coef(kem4_ind3_Coyote, type="matrix")$R
coef(kem4_ind3_Coyote, type="matrix")$Q
coef(kem4_ind3_Coyote, type="matrix")$U
coef(kem4_ind3_Coyote, type="matrix")$B

ggplot2::autoplot(kem4_ind3_Coyote, plot.type = "xtT") 


#Plot showing abundance change since time series start
par(mfrow = c(1,1))
matplot(years, 
        t(kem4_ind3_Coyote$states-kem4_ind3_Coyote$states[,1]), 
        xlab="Year", ylab="abundance index",
        type="l",lwd=2,col="black")
legend("topleft",c("BL","DE", "DP","TB", "TP"),lwd=2,lty=c(1:3),bty="n")



tidy(kem4_ind3_Coyote)

df_aic <- df_aic %>% add_row(model = "Site ind 5 pops + Coyote", aic = kem4_ind3_Coyote$AIC)
df_aic

## coyote continuous
Z.model=factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5))
R.model="diagonal and unequal"
U.model="unequal"
Q.model="diagonal and equal"
B.model="diagonal and unequal"
# if Z is 1:15
# C.model=matrix(c("BL_1","OTH","OTH","OTH","OTH",
#                  "BL_1","OTH","OTH","OTH","OTH",
#                  "BL_1","OTH","OTH","OTH","OTH",
#                  
#                  "OTH","DE_1","OTH","OTH","OTH",
#                  "OTH","DE_1","OTH","OTH","OTH",
#                  "OTH","DE_1","OTH","OTH","OTH",
#                  
#                  "OTH","OTH","DP_1","OTH","OTH",
#                  "OTH","OTH","DP_1","OTH","OTH",
#                  "OTH","OTH","DP_1","OTH","OTH",
#                  
#                  "OTH","OTH","OTH","OTH","OTH",
#                  "OTH","OTH","OTH","OTH","OTH",
#                  "OTH","OTH","OTH","OTH","OTH",
#                  
#                  "OTH","OTH","OTH","OTH","OTH",
#                  "OTH","OTH","OTH","OTH","OTH",
#                  "OTH","OTH","OTH","OTH","OTH"),
#               nrow = 15, ncol = 5,
#               byrow = TRUE)

# if Z is rep(1:5), times = 3
# 6 min RT
C.model=matrix(c("BL_1","OTH","OTH","OTH","OTH",
                 "OTH","DE_1","OTH","OTH","OTH",
                 "OTH","OTH","DP_1","OTH","OTH",
                 "OTH","OTH","OTH","OTH","OTH",
                 "OTH","OTH","OTH","OTH","OTH"),
               nrow = 5, ncol = 5,
               byrow = TRUE)

kem4_ind3_Coyote_rate=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                       C = C.model,
                                       c = small_c_Coyote_Rate),
                       control=list(maxit=5000, safe=TRUE, trace = 1)) 
beepr::beep(0)
kem4_ind3_Coyote_rate$AIC 

plot(kem4_ind3_Coyote_rate), plot.type="model.resids.ytT")

best.fit=fits[min.AICc][[1]]
par(mfrow = c(1,1))
matplot(years, t(best.fit$states-best.fit$states[,1]), 
        xlab="Year", ylab="abundance index",
        type="l",lwd=2,col="black")
legend("topleft",c("Bolinas Lagoon","South (DE/DP)","North (TB/TP)"),lwd=2,lty=c(1:3),bty="n")


CIs <- MARSSparamCIs(kem4_ind3_Coyote_rate)
CIs






df_aic <- df_aic %>% add_row(model = "Site ind 5 pops + Coyote_rate", aic = kem4_ind3_Coyote_rate$AIC)
df_aic








#each ts different
Z.model=factor(c(1:15))
R.model="diagonal and unequal"
U.model="unequal"
Q.model="diagonal and equal"
B.model="diagonal and unequal"
# if Z is 1:15
C.model=matrix(c("BL_1","OTH","OTH","OTH","OTH",
                 "BL_1","OTH","OTH","OTH","OTH",
                 "BL_1","OTH","OTH","OTH","OTH",

                 "OTH","DE_1","OTH","OTH","OTH",
                 "OTH","DE_1","OTH","OTH","OTH",
                 "OTH","DE_1","OTH","OTH","OTH",

                 "OTH","OTH","DP_1","OTH","OTH",
                 "OTH","OTH","DP_1","OTH","OTH",
                 "OTH","OTH","DP_1","OTH","OTH",

                 "OTH","OTH","OTH","OTH","OTH",
                 "OTH","OTH","OTH","OTH","OTH",
                 "OTH","OTH","OTH","OTH","OTH",

                 "OTH","OTH","OTH","OTH","OTH",
                 "OTH","OTH","OTH","OTH","OTH",
                 "OTH","OTH","OTH","OTH","OTH"),
              nrow = 15, ncol = 5,
              byrow = TRUE)

# if Z is rep(1:5), times = 3
# 28 min RT !
# Q did not converge

kem4_ind5_Coyote=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                       C = C.model,
                                       c = small_c),
                       control=list(maxit=5000, safe=TRUE, trace = 1)) 
beepr::beep(0)
kem4_ind5_Coyote$AIC 



# 15 pops is 2 AIC higher than five pops.

CIs <- MARSSparamCIs(kem4_ind5_Coyote)
CIs

#ggplot2::autoplot(kem4_ind5_Coyote) 



df_aic <- df_aic %>% add_row(model = "Site ind all ages + pops + Coyote", aic = kem4_ind5_Coyote$AIC)
df_aic

######
# molt unrelated to breeding season at each site
# RT = 10 min
Z.model=factor(c(1,2,1,3,4,3,5,6,5,7,8,7,9,10,9))
R.model="diagonal and unequal"
U.model="unequal"
Q.model="diagonal and equal"
B.model="diagonal and unequal"#"unequal"
C.model=matrix(c("BL_B","OTH","OTH","OTH","OTH",
                 "BL_M","OTH","OTH","OTH","OTH",
                 
                 "OTH","DE_B","OTH","OTH","OTH",
                 "OTH","DE_M","OTH","OTH","OTH",
                 
                 "OTH","OTH","DP_B","OTH","OTH",
                 "OTH","OTH","DP_M","OTH","OTH",
                 
                 "OTH","OTH","OTH","OTH","OTH",
                 "OTH","OTH","OTH","OTH","OTH",
                 
                 "OTH","OTH","OTH","OTH","OTH",
                 "OTH","OTH","OTH","OTH","OTH"),
               nrow = 10, ncol = 5,
               byrow = TRUE)
kem4_MoltInd.coyote=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                   C = C.model,
                                   c = small_c),
                   control=list(maxit=5000, safe=TRUE)) 

CIs <- MARSSparamCIs(kem4_MoltInd.coyote)
CIs

beepr::beep(0)
kem4_MoltInd.coyote$AIC 

df_aic <- df_aic %>% add_row(model = "Molt ind coyote", aic = kem4_MoltInd.coyote$AIC)


###############

######
# Pups unrelated to breeding adults and molt season at each site
# RT = 12 min
Z.model=factor(c(1,1,2,3,3,4,5,5,6,7,7,8,9,9,10))
R.model="diagonal and unequal"
U.model="unequal"
Q.model="diagonal and equal"
B.model="diagonal and unequal"#"unequal"
C.model=matrix(c("BL_A","OTH","OTH","OTH","OTH",
                 "BL_P","OTH","OTH","OTH","OTH",
                 
                 "OTH","DE_A","OTH","OTH","OTH",
                 "OTH","DE_P","OTH","OTH","OTH",
                 
                 "OTH","OTH","DP_A","OTH","OTH",
                 "OTH","OTH","DP_P","OTH","OTH",
                 
                 "OTH","OTH","OTH","OTH","OTH",
                 "OTH","OTH","OTH","OTH","OTH",
                 
                 "OTH","OTH","OTH","OTH","OTH",
                 "OTH","OTH","OTH","OTH","OTH"),
               nrow = 10, ncol = 5,
               byrow = TRUE)
kem4_Pup_Ind.coyote=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                          C = C.model,
                                          c = small_c),
                          control=list(maxit=5000, safe=TRUE)) 

CIs <- MARSSparamCIs(kem4_Pup_Ind.coyote)
CIs

beepr::beep(0)
kem4_Pup_Ind.coyote$AIC 

df_aic <- df_aic %>% add_row(model = "Pup ind coyote", aic = kem4_Pup_Ind.coyote$AIC)
df_aic

###############







#same as first but this time R is equal
Z.model=factor(c(1,2,3,4,5))
R.model="diagonal and equal"
U.model="unequal"
Q.model="diagonal and equal"
B.model="identity"
kem4_ind7=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model),
                control=list(maxit=1000, safe=TRUE)) 
kem4_ind7$AIC # 4.137058

df_aic <- df_aic %>% add_row(model = "all ind 7", aic = kem4_ind7$AIC)

## Ocean vs Bay
Z.model=factor(c(1,1,1,1,1,1,2,2,2,1,1,1,2,2,2))
R.model="diagonal and equal"
U.model="unequal"
Q.model="diagonal and equal"
B.model="diagonal and equal"#"unequal" 
kemOB2 = MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model),
                    control=list(maxit=500, safe=TRUE)) 
kemOB2$AIC #35.93709

df_aic <- df_aic %>% add_row(model = "ocean vs bay", aic = kemOB2$AIC)

## adult_molt_pup diff
Z.model=factor(c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3))
R.model="diagonal and equal"
U.model="unequal"
Q.model="diagonal and equal"
B.model="diagonal and equal"#"unequal" 
kem_adult_molt_pup = MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model),
               control=list(maxit=500, safe=TRUE)) 
kem_adult_molt_pup$AIC #35.93709

df_aic <- df_aic %>% add_row(model = "adult_molt_pup diff", aic = kem_adult_molt_pup$AIC)
df_aic


## adult_molt_pup diff
Z.model=factor(c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3))
R.model="diagonal and unequal"
U.model="equal"
Q.model="diagonal and equal"
B.model="diagonal and equal"#"unequal" 
kem_adult_molt_pup_u_equal = MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model),
                           control=list(maxit=500, safe=TRUE)) 
kem_adult_molt_pup_u_equal$AIC #35.93709

df_aic <- df_aic %>% add_row(model = "adult_molt_pup diff_U_unequal", aic = kem_adult_molt_pup_u_equal$AIC)
df_aic


##############

