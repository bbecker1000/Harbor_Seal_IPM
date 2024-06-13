#4b covariate prep for 1997-2023
# all classes

# env data prep for all years


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



### Use data only > 1996? -----------------
# MEI <- MEI %>%
#   filter(Year>1997)


#View(MEI)
##cut to sealData time series

MEI <- MEI %>% filter(Year < 2024 & Year > 1996)

#ESeal at DP in 2003 and 2004.  Count as a coyote disturbance
MEI$Coyote_DP <- ifelse(MEI$Year == 2003, 1, MEI$Coyote_DP)
MEI$Coyote_DP <- ifelse(MEI$Year == 2004, 1, MEI$Coyote_DP)

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



###Matrix with coyote_3yr and BEAUTI and BEAUTI-lag

small_c_Coyote_3yr_BEUTI_BEUTI_lag <-tibble(MEI[,c(20:22, 17, 25)]) #scale for the covariate plots!
small_c_Coyote_3yr_BEUTI_BEUTI_lag <- scale(small_c_Coyote_3yr_BEUTI_BEUTI_lag)
small_c_Coyote_3yr_BEUTI_BEUTI_lag <- as_tibble(small_c_Coyote_3yr_BEUTI_BEUTI_lag)

## zeros for DR, PB, PRH
small_c_Coyote_3yr_BEUTI_BEUTI_lag$Coyote_3yr_PRH <- 0
small_c_Coyote_3yr_BEUTI_BEUTI_lag$Coyote_3yr_TB <- 0
small_c_Coyote_3yr_BEUTI_BEUTI_lag$Coyote_3yr_TP <- 0


small_c_Coyote_3yr_BEUTI_BEUTI_lag <- small_c_Coyote_3yr_BEUTI_BEUTI_lag[,c(1:3, 6:8, 4:5)]
small_c_Coyote_3yr_BEUTI_BEUTI_lag <- as.matrix(t(small_c_Coyote_3yr_BEUTI_BEUTI_lag)) # 









