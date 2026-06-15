# env data prep for all years


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
