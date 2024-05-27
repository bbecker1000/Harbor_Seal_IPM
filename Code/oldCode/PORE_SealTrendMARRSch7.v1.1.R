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

top.plot <- ggplot(top1_all_data, aes(Year, log(Count))) + 
  geom_point() + 
  geom_smooth() + 
  geom_vline(xintercept = c(1988, 1994), lty = 3) +
  #scale_y_continuous(trans='log10') +
  geom_vline(xintercept = 1994, lty = 3) 
  
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

#use only data post 1997 ? ------------------------------------
all_data.MARSS <- all_data.MARSS %>%
                      filter(Year>1996)
#--------------------------------------------------------------

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

### Use data only > 1996? -----------------
MEI <- MEI %>%
  filter(Year>1996)

####--------------------------------------



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



##########################################################################
### 8.3 Single population model with independent and (non)identical errors
##########################################################################
#Code to fit the single population model with i.i.d. errors
#Read in data


########################################################################

## Parameter key ########
## Z = design matrix = Spatial population structure
## R = observation errors          equal
## U = growth parameter            (un)equal, or matrix [Z x 1]
## Q = hidden state process        diagonal and (un)equal
#                                  equalvalcov or unconstrained
## B = effect of column on row     unequal (these are the interactions)

#####################################################################











df_aic <- data.frame(model=character(), aic=integer())



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


#by site, Coyote_01
#RT 10+ min
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

ggplot2::autoplot(kem4_ind3_Coyote)#, plot.type = "xtT") 
  
  
plot(kem4_ind3_Coyote, plot.type="fitted.ytT") +
  scale_x_continuous()

#Plot showing abundance change since time series start
#convert to ggplot
par(mfrow = c(1,1))
matplot(years, 
        t(kem4_ind3_Coyote$states-kem4_ind3_Coyote$states[,1]), 
        xlab="Year", ylab="abundance index",
        type="l",lwd=2,col="black")
legend("topleft",c("BL","DE", "DP","TB", "TP"),lwd=2,lty=c(1:5),bty="n")



tidy(kem4_ind3_Coyote)
CIs <- MARSSparamCIs(kem4_ind3_Coyote)
CIs


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

plot(kem4_ind3_Coyote, plot.type="fitted.ytT")

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


################ small_c_Coyote_01_PDO_MAR
#by site, Coyote_01
#RT 6 min
Z.model=factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5))
R.model="diagonal and unequal"
U.model="unequal"
Q.model="equalvarcov" #diagonal and equal"
B.model="diagonal and unequal"

# 6 min RT

C.model=matrix(list("BL_1",0,0,0,0,
                    0,"DE_1",0,0,0,
                    0,0,"DP_1",0,0,
                    0,0,0,0,0,
                    0,0,0,0,0,
                    "PDO","PDO","PDO","PDO","PDO"),
               nrow = 6, ncol = 5,
               byrow = TRUE)





kem4_ind3_Coyote_01_PDO_MAR=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                       C = C.model,
                                       c = small_c_Coyote_01_PDO_MAR),
                       
                       control=list(maxit=5000, safe=TRUE, trace = 1)) 
beepr::beep(0)
kem4_ind3_Coyote_01_PDO_MAR$AIC 

autoplot(kem4_ind3_Coyote_01_PDO_MAR)
plot(kem4_ind3_Coyote_01_PDO_MAR, plot.type="fitted.ytT")





CIs <- MARSSparamCIs(kem4_ind3_Coyote_01_PDO_MAR, alpha = 0.11)
CIs

coef(kem4_ind3_Coyote_01_PDO_MAR, type="matrix")$R
coef(kem4_ind3_Coyote_01_PDO_MAR, type="matrix")$Q
coef(kem4_ind3_Coyote_01_PDO_MAR, type="matrix")$U
coef(kem4_ind3_Coyote_01_PDO_MAR, type="matrix")$B

ggplot2::autoplot(kem4_ind3_Coyote_01_PDO_MAR)#, plot.type = "xtT") 


plot(kem4_ind3_Coyote_01_PDO_MAR, plot.type="fitted.ytT") +
  scale_x_continuous()

#Plot showing abundance change since time series start
#convert to ggplot
par(mfrow = c(1,1))
matplot(years, 
        t(kem4_ind3_Coyote_01_PDO_MAR$states-kem4_ind3_Coyote_01_PDO_MAR$states[,1]), 
        xlab="Year", ylab="abundance index",
        type="l",lwd=2,col="black")
legend("topleft",c("BL","DE", "DP","TB", "TP"),lwd=2,lty=c(1:5),bty="n")



tidy(kem4_ind3_Coyote_01_PDO_MAR)


df_aic <- df_aic %>% add_row(model = "Site ind 5 pops + Coyote + PDO", aic = kem4_ind3_Coyote_01_PDO_MAR$AIC)
df_aic



###############-------------
################ small_c_PDO_MAR
#by site, 
#RT 6 min
Z.model=factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5))
R.model="diagonal and unequal"
U.model="unequal"
Q.model="diagonal and equal"
B.model="diagonal and unequal"
C.model=matrix(c(#"BL_1","OTH","OTH","OTH","OTH",
                 #"OTH","DE_1","OTH","OTH","OTH",
                 #"OTH","OTH","DP_1","OTH","OTH",
                 #"OTH","OTH","OTH","OTH","OTH",
                 #"OTH","OTH","OTH","OTH","OTH",
                 "PDO","PDO","PDO","PDO","PDO"),
               nrow = 1, ncol = 5,
               byrow = TRUE)

kem4_ind3_PDO_MAR=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                  C = C.model,
                                                  c = small_c_PDO_MAR),
                                  
                                  control=list(maxit=5000, safe=TRUE, trace = 1)) 
beepr::beep(0)

kem4_ind3_PDO_MAR$AIC 

plot(kem4_ind3_PDO_MAR, plot.type="fitted.ytT")

CIs <- MARSSparamCIs(kem4_ind3_PDO_MAR)
CIs

coef(kem4_ind3_PDO_MAR, type="matrix")$R
coef(kem4_ind3_PDO_MAR, type="matrix")$Q
coef(kem4_ind3_PDO_MAR, type="matrix")$U
coef(kem4_ind3_PDO_MAR, type="matrix")$B

ggplot2::autoplot(kem4_ind3_PDO_MAR)#, plot.type = "xtT") 


plot(kem4_ind3_PDO_MAR, plot.type="fitted.ytT")

#Plot showing abundance change since time series start
#convert to ggplot
par(mfrow = c(1,1))
matplot(years, 
        t(kem4_ind3_PDO_MAR$states)-kem4_ind3_PDO_MAR$states[,1]), 
        xlab="Year", ylab="abundance index",
        type="l",lwd=2,col="black")
legend("topleft",c("BL","DE", "DP","TB", "TP"),lwd=2,lty=c(1:5),bty="n")

# try with ggplot

d <- as_tibble(t(kem4_ind3_PDO_MAR$states-kem4_ind3_PDO_MAR$states[,8]))
#add header names
names(d)[c(1:5)] <- c("BL", "DE", "DP", "TB", "TP")
d <- as_tibble(cbind(years,d))

d2 <- d %>% pivot_longer(cols = c(2:6), names_to = "Subsite", values_to = "log_est")
  
  

ggplot(d2, aes(x = years, y = log_est, color = Subsite)) +
  geom_line(linewidth = 1) +
  xlim(1982, 2022) +
  ylim(-1, 1.7)




tidy(kem4_ind3_PDO_MAR)


df_aic <- df_aic %>% add_row(model = "Site ind 5 pops + PDO", aic = kem4_ind3_PDO_MAR$AIC)
df_aic



###############-------------






#each ts different
Z.model=factor(c(1:15))
R.model="diagonal and unequal"
U.model="unequal"
Q.model="diagonal and equal"
B.model="diagonal and unequal"
# if Z is 1:15
C.model=matrix(list("BL_1",0,0,0,0,
                 "BL_1",0,0,0,0,
                 "BL_1",0,0,0,0,

                 0,"DE_1",0,0,0,
                 0,"DE_1",0,0,0,
                 0,"DE_1",0,0,0,

                 0,0,"DP_1",0,0,
                 0,0,"DP_1",0,0,
                 0,0,"DP_1",0,0,

                 0,0,0,0,0,
                 0,0,0,0,0,
                 0,0,0,0,0,

                 0,0,0,0,0,
                 0,0,0,0,0,
                 0,0,0,0,0),
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
C.model=matrix(c("BL_B","OTH","OTH","OTH","OTH", #1
                 "BL_M","OTH","OTH","OTH","OTH", #2
                 
                 "OTH","DE_B","OTH","OTH","OTH", #3
                 "OTH","DE_M","OTH","OTH","OTH", #4
                 
                 "OTH","OTH","DP_B","OTH","OTH", #5
                 "OTH","OTH","DP_M","OTH","OTH", #6
                 
                 "OTH","OTH","OTH","OTH","OTH",  #7
                 "OTH","OTH","OTH","OTH","OTH",  #8
                 
                 "OTH","OTH","OTH","OTH","OTH",  #9
                 "OTH","OTH","OTH","OTH","OTH"), #10
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
autoplot(kem4_MoltInd.coyote)

###############

######
# Pups unrelated to breeding adults and molt season at each site
# coyote
# RT = 12 min
Z.model=factor(c(1,1,2,3,3,4,5,5,6,7,7,8,9,9,10))
R.model="diagonal and equal"
U.model="unequal"   #growth rates same for each z
#U.model=matrix(c("u1", "u2", "u3", "u4", "u5", "u6", "u7", "u8", "u9", "u10"
#                 ), 10, 1)
Q.model="equalvarcov"
B.model="diagonal and unequal"#"unequal"
C.model=matrix(list("BL_A",0,0,0,0,
                 "BL_P",0,0,0,0,
                 
                 0,"DE_A",0,0,0,
                 0,"DE_P",0,0,0,
                 
                 0,0,"DP_A",0,0,
                 0,0,"DP_P",0,0,
                 
                 0,0,0,0,0,
                 0,0,0,0,0,
                 
                 0,0,0,0,0,
                 0,0,0,0,0),
               nrow = 10, ncol = 5,
               byrow = TRUE)
kem4_Pup_Ind.coyote.1997_2022=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                          C = C.model,
                                          c = small_c_Coyote_01),
                          control=list(maxit=5000, safe=TRUE)) 

CIs <- MARSSparamCIs(kem4_Pup_Ind.coyote.1997_2022)
CIs

beepr::beep(0)
kem4_Pup_Ind.coyote.1997_2022$AIC 

autoplot(kem4_Pup_Ind.coyote.1997_2022, plot.type = "fitted.ytT") + ylim(3.5,7.5)

plot(kem4_Pup_Ind.coyote.1997_2022, plot.type = "fitted.ytT")


df_aic <- df_aic %>% add_row(model = "Pup ind coyote", aic = kem4_Pup_Ind.coyote.1997_2022$AIC)
df_aic

###############


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

## PDO only effects pups 
## R: diagonal and unequal
## COyote effects both age classes
## RT = 13 min
Z.model=factor(c(1,1,2,3,3,4,5,5,6,7,7,8,9,9,10))
R.model="diagonal and unequal"  # Observation errors   "equal" better than "unequal"
U.model="unequal"               # Growth Parameter       
Q.model="diagonal and equal"    # hidden state process   #equalvarcov,  unconstrained
B.model="diagonal and unequal"  # effect of column on row 

C.model=matrix(list("BL_A",0,0,0,0,"PDO_A", 
                 "BL_P",0,0,0,0,"PDO_P",
                 
                 0,"DE_A",0,0,0,"PDO_A",
                 0,"DE_P",0,0,0,"PDO_P",
                 
                 0,0,"DP_A",0,0,"PDO_A",
                 0,0,"DP_P",0,0,"PDO_P",
                 
                 0,0,0,0,0,"PDO_A",
                 0,0,0,0,0,"PDO_P",
                 
                 0,0,0,0,0,"PDO_A",
                 0,0,0,0,0,"PDO_P"
                 
                 ),
               nrow = 10, ncol = 6,
               byrow = TRUE)

kem4_Pup_Ind.coyote_PDO=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                          C = C.model,
                                          c = small_c_Coyote_01_PDO_MAR),
                          control=list(maxit=5000, safe=TRUE)) 

CIs <- MARSSparamCIs(kem4_Pup_Ind.coyote_PDO)
CIs

beepr::beep(0)
kem4_Pup_Ind.coyote_PDO$AIC 

df_aic <- df_aic %>% add_row(model = "Pup ind coyote + PDO + R.unequal", aic = kem4_Pup_Ind.coyote_PDO$AIC)
df_aic

# try with ggplot

ggplot2::autoplot(kem4_Pup_Ind.coyote_PDO)#, plot.type = "xtT")


d <- as_tibble(t(kem4_ind3_PDO_MAR$states-kem4_ind3_PDO_MAR$states[,8]))
#add header names
names(d)[c(1:5)] <- c("BL", "DE", "DP", "TB", "TP")
d <- as_tibble(cbind(years,d))

d2 <- d %>% pivot_longer(cols = c(2:6), names_to = "Subsite", values_to = "log_est")



ggplot(d2, aes(x = years, y = log_est, color = Subsite)) +
  geom_line(linewidth = 1) +
  xlim(1982, 2022) +
  ylim(-1, 1.7)




## PDO only effects pups 
## R: diagonal and EQUAL
## COyote effects both age classes
## RT = >12 min
Z.model=factor(c(1,1,2,3,3,4,5,5,6,7,7,8,9,9,10))
R.model="diagonal and equal"  # Observation errors   diagonal and unequal
U.model="unequal"               # Growth Parameter       
Q.model="diagonal and equal"    # hidden state process   #equalvarcov,  unconstrained
B.model="diagonal and unequal"  # effect of column on row 
# C.model=matrix(c("BL_A","OTH","OTH","OTH","OTH","PDO_A", 
#                  "BL_P","OTH","OTH","OTH","OTH","PDO_P",
#                  
#                  "OTH","DE_A","OTH","OTH","OTH","PDO_A",
#                  "OTH","DE_P","OTH","OTH","OTH","PDO_P",
#                  
#                  "OTH","OTH","DP_A","OTH","OTH","PDO_A",
#                  "OTH","OTH","DP_P","OTH","OTH","PDO_P",
#                  
#                  "OTH","OTH","OTH","OTH","OTH","PDO_A",
#                  "OTH","OTH","OTH","OTH","OTH","PDO_P",
#                  
#                  "OTH","OTH","OTH","OTH","OTH","PDO_A",
#                  "OTH","OTH","OTH","OTH","OTH","PDO_P"),
#                nrow = 10, ncol = 6,
#                byrow = TRUE)

C.model=matrix(c("BL_A",0,0,0,0, 
                 "BL_P",0,0,0,0,
                 
                 0,"DE_A",0,0,0,
                 0,"DE_P",0,0,0,
                 
                 0,0,"DP_A",0,0,
                 0,0,"DP_P",0,0,
                 
                 0,0,0,0,0,
                 0,0,0,0,0,
                 
                 0,0,0,0,0,
                 0,0,0,0,0,
                 
                 "PDO_A","PDO_A","PDO_A","PDO_A","PDO_A",
                 "PDO_P","PDO_P","PDO_P","PDO_P","PDO_P"
                 ),
               nrow = 12, ncol = 5,
               byrow = TRUE)

kem4_Pup_Ind.coyote_PDO_r.equal=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
                                                      C = C.model,
                                                      c = small_c_Coyote_01_PDO_MAR),
                                      control=list(maxit=5000, safe=TRUE)) 

CIs <- MARSSparamCIs(kem4_Pup_Ind.coyote_PDO_r.equal)
CIs

beepr::beep(0)
kem4_Pup_Ind.coyote_PDO_r.equal$AIC 

df_aic <- df_aic %>% add_row(model = "Pup ind coyote + PDO + R.equal", aic = kem4_Pup_Ind.coyote_PDO_r.equal$AIC)
df_aic

ggplot2::autoplot(kem4_Pup_Ind.coyote_PDO_r.equal)




