#adults and molts all years models



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
Z.model = factor(rep(c(1:3), times = 8))
#R.model = "diagonal and equal"
R.model=diag(0.025, 24)
Q.model=diag(0.025, 3)

AllYearAllClass.Z.class=MARSS(dat, model=list(
  Z=Z.model, 
  R=R.model,
  Q=Q.model))
  #x0 = x0.model, 
  #tinitx=1, 


Z.model = factor(rep(c(1:3), each = 8))
R.model=diag(0.025, 24)
Q.model=diag(0.025, 3)
AllYearAllClassZ.Site=MARSS(dat, model=list(
  Z=Z.model, 
  R=R.model,
  Q=Q.model))
#x0 = x0.model, 
#tinitx=1, 

Z.model = factor(1:24)
R.model=diag(0.025, 24)
Q.model=diag(0.025, 24)
AllYearAllClass=MARSS(dat, model=list(
  Z=Z.model, 
  R=R.model,
  Q=Q.model))

#timevarying u by site and class
#timevarying u by site and class
U1 <- matrix(c("t1_BL_A","t1_BL_M","t1_BL_P",
               "t1_DE_A","t1_DE_M","t1_DE_P",
               "t1_DP_A","t1_DP_M","t1_DP_P",
               "t1_DR_A","t1_DR_M","t1_DR_P",
               "t1_PB_A","t1_PB_M","t1_PB_P",
               "t1_PRH_A","t1_PRH_M","t1_PRH_P",
               "t1_TB_A","t1_TB_M","t1_TB_P",
               "t1_TP_A","t1_TP_M","t1_TP_P"),24, 1)
U2 <- matrix(c("t2_BL_A","t2_BL_M","t2_BL_P",
               "t2_DE_A","t2_DE_M","t2_DE_P",
               "t2_DP_A","t2_DP_M","t2_DP_P",
               "t2_DR_A","t2_DR_M","t2_DR_P",
               "t2_PB_A","t2_PB_M","t2_PB_P",
               "t2_PRH_A","t2_PRH_M","t2_PRH_P",
               "t2_TB_A","t2_TB_M","t2_TB_P",
               "t2_TP_A","t2_TP_M","t2_TP_P"),24, 1)
Ut.Site.Class <- array(U2, dim = c(dim(U1), dim(dat)[2]))
TT <- dim(dat)[2]
Ut.Site.Class[, , 1:ceiling(TT / 2)] <- U1


Z.model = factor(1:24)
R.model=diag(0.025, 24)
Q.model=diag(0.025, 24)
cntl.list <- list(maxit = 3000, allow.degen = FALSE)
AllYearAllClass.Ut=MARSS(dat, model=list(
  Z=Z.model, 
  R=R.model,
  Q=Q.model,
  U=Ut.Site.Class),
  control = cntl.list)

beepr::beep(0)




AllYearAllClass$AIC
AllYearAllClass.Z.class$AIC
AllYearAllClassZ.Site$AIC
AllYearAllClass.Ut$AIC

#use autoplot shortcut to gather data for plotting
p1<-autoplot(AllYearAllClass.Ut, plot.type = "fitted.ytT") +
  ylim(-2,8) +
  #scale_y_continuous(trans = exp(AllYearAllClass$ytT)) +
  ylab("estimate (log count)")# xtT



#now access the p1 object data
rownames <-as.character(p1[["data"]][[".rownames"]])
t <- p1[["data"]][["t"]]
y <- p1[["data"]][["y"]]
fitted <- p1[["data"]][[".fitted"]]
se <- p1[["data"]][[".se"]]
conf.low <- p1[["data"]][[".conf.low"]]
conf.up <- p1[["data"]][[".conf.up"]]

df <- data.frame(rownames, t, y, fitted, se, conf.low, conf.up)
head(df)

df$Site <- sub("_.*", "", df$rownames)
unique(df$Site)
df$Class <- sub(".*_", "", df$rownames)
head(df)
df$Year <- df$t+1974

#upper CIs for DR molt are unreasonable and go off chart.
#fix them to zero
head(df)

#df <- ifelse(df$rownames == "DR_MOLTING", df$conf.up == df$fitted, df$conf.up == df$conf.up)
# df <- within(df, conf.up[rownames == 'DR_MOLTING'] <- df$fitted)

#remove predictions for DR and PB prior to 1997 since no data
df$y <- ifelse(df$Site == "DR" & df$Year <= 1997, NA, df$y)
df$y <- ifelse(df$Site == "PB" & df$Year <= 1997, NA, df$y)
df$fitted <- ifelse(df$Site == "DR" & df$Year <= 1997, NA, df$fitted)
df$fitted <- ifelse(df$Site == "PB" & df$Year <= 1997, NA, df$fitted)
df$conf.low <- ifelse(df$Site == "DR" & df$Year <= 1997, NA, df$conf.low)
df$conf.low <- ifelse(df$Site == "PB" & df$Year <= 1997, NA, df$conf.low)
df$conf.up <- ifelse(df$Site == "DR" & df$Year <= 1997, NA, df$conf.up)
df$conf.up <- ifelse(df$Site == "PB" & df$Year <= 1997, NA, df$conf.up)

#plot MARSS predictions
ggplot(df, aes(Year, exp(y), fill = Class, color = Class)) +
  geom_ribbon(aes(ymin = exp(conf.low), ymax = exp(conf.up), fill = Class), 
              alpha = 0.3, colour = NA) + 
  geom_point(alpha = 0.75) +
  geom_line(aes(x = Year, y = exp(fitted)), size = 1) +
  ylim(0,2000) +
  ylab("Seals") + 
  scale_x_continuous(breaks=seq(1980,2020,20)) +
  theme_gray(base_size = 18) +
  facet_wrap(.~Site, ncol = 4)


#total pop plots.
#exponentiate the estimates
df.total <- df
head(df.total)
df.total$y <- exp(df.total$y)
df.total$fitted <- exp(df.total$fitted)
df.total$conf.low <- exp(df.total$conf.low)
df.total$conf.up <- exp(df.total$conf.up)

#sum by year and class  
df.total.sum <- 
  df.total %>%
    dplyr::group_by(Class, Year) %>%
    dplyr::summarize(y_sum=sum(y, na.rm = T),
              fitted_sum=sum(fitted, na.rm = T),
              conf.low_sum=sum(conf.low, na.rm = T),
              conf.up_sum=sum(conf.up, na.rm = T)) 


#plot - recall that errors for DR and PB < 1997 have been removed
ggplot(df.total.sum, aes(Year, y_sum, color = Class)) +
  geom_ribbon(aes(ymin = conf.low_sum, ymax = conf.up_sum, fill = Class), 
              alpha = 0.3, colour = NA) +
  #geom_point() +
  geom_line(aes(x = Year, y = fitted_sum), size = 1) +
  ylim(0, 6000) +
  ylab("Seals") + 
  scale_x_continuous(breaks=seq(1980,2020,20)) +
  theme_gray(base_size = 18)
  


    


#Plot the U's from best model
AllYearAllClass.Ut
CIs <- MARSSparamCIs(AllYearAllClass.Ut, alpha = 0.89)
CIs

CIs <- tidy(CIs)
U <- CIs %>% filter(str_detect(term, "^U"))
U$time <- rep(c(1:2), each = 24)
#get the class from term
U$Class <- sub(".*_", "", U$term)
U$Class2 <- ifelse(U$Class == "A", "Adult",
                  ifelse(U$Class == "M", "Molt", "Pup" ))

#get the site from term
U$Site = sapply(strsplit(U$term, "_"), function(x) x[2])

#make a forest plot.

Time_Periods <- c(
  `1` = "1975-2004",
  `2` = "2005-2023")


#parameter plot of U'ss
ggplot(U, aes(Site, estimate, color = Class2)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.up), 
                  position = position_dodge(width = 0.9), size = 0.75) +
  coord_flip() + 
  geom_hline(yintercept = 0, linetype = 2) +
  ylab("Annual rate of change") + 
  theme_gray(base_size = 18) +
  ylim(-0.1, 0.3) +
  facet_wrap(.~time, labeller = as_labeller(Time_Periods))

#total pop plots.




ytT <- as_tibble(t(AllYearAllClass[["ytT"]]))
ytT.se <- as_tibble(t(AllYearAllClass[["ytT.se"]]))

ytT.se.tidy <- ytT.se %>% pivot_longer(
  everything(),
  values_to = "se"
)


ytT.tidy <- ytT %>% pivot_longer(
  everything(),
  values_to = "count"
)

ytT.df <- cbind(ytT.se.tidy, ytT.tidy)
ytT.df$Year <- rep(c(1975:2023), each = 24)
ytT.df <- ytT.df[,-3]
ytT.df <- ytT.df %>%
  separate(name, "_")
ytT.df$Class <- rep(c("Adult", "Molting", "Pup"))
colnames(ytT.df)[1] <- "Site"
ytT.df$count.pred <- exp(ytT.df$count)
ytT.df$count.pred.se <- exp(ytT.df$c)
head(ytT.df)

ggplot(ytT.df, aes(Year, count.pred, color = Class)) +
         
         geom_ribbon(aes(ymin = count.pred-count.pred*se, ymax = count.pred+count.pred*se)) + 
        #             alpha = 0.5) +
         #geom_point(alpha = 0.6) +
         #geom_smooth(method = "loess") +
  geom_line() +
  theme_(base_size = 18) +
         facet_grid(.~Site)
  






## 2024-05-26 HYPOTHESES for Breed and Molt

## Z hypotheses

# 12 Independent pops since want to understand growth rates...know different...
Z.model=factor(c(1:12))  # 6 sites x 2 age classes
 

# R1: diagonal and equal ...no reason to assume otherwise, standardized methods, repeat visits
R1.model="diagonal and equal"

# Q1: diagonal and equal ... no reason to assume otherwise
Q1.model="diagonal and equal"

# B1 design matrix based on predicted north and south colony interactions (Sarah Allen)
# B.models
#breed then molt, so only breed effect on molt "BB"
#breeding fidelity = AA
#molting fidelity = CC
#A / B / C / D are within zone interactions by breed or molt

B.model=matrix(list("AA",0,    "A",0,   "A",0,   0,0,   0,0,  0,0,  #BL
                    "BB","CC",     0,"A",   0,"A",  0,0,   0,0,  0,0,
                    
                    "A",0,  "AA",0,   "A",0,  "A",0,    0,0,  0,0, #DE
                    0,"A",  "BB","CC",   0,"A",  0,"A",    0,0,  0,0,
                    
                    "A",0,  "A",0,   "AA",0,  "A",0,   0,0,  0,0,  #DP
                    0,"A",  0,"A",   "BB","CC",  0,"A",   0,0,  0,0,
                    
                    "B",0,  "B",0, "B",0,    "AA",0,   "C",0,  "C",0,  #PRH
                    0,"B",  0,"B",   0,"B",  "BB","CC",   0,"C",  0,"C",
                    
                    0,0,   0,0,     0,0,   "D",0,   "AA",0,   "D",0,  #TB
                    0,0,   0,0,     0,0,   0,"D",   "BB","CC",   0,"D",
                    
                    0,0,   0,0,     0,0,   0,0,    "D",0,   "AA",0,  #TP
                    0,0,   0,0,     0,0,   0,0,    0,"D",   "BB","CC"),
               
               nrow = 12, ncol = 12,
               byrow = TRUE)

## COVARIATES

# C1A: Spring UI
# C1B: Spring UI Lag (egg implantation)
# C2: MEI
# C3: NPGO_JAN (early for food web and egg implantation)
# C4: Coyote 3yr sum by site


# Initial models show R.diag and Q.diag between 0.01 and 0.04.  So setting both to 0.025


#U constant
U.Equal="equal"

#constant growth rates for Adult and Molt
U.Class=matrix(list("UA", "UM","UA", "UM", "UA", "UM", "UA", "UM", "UA", "UM", "UA", "UM"), 
               nrow = 12, ncol = 1,
               byrow = TRUE)# separate growth rates for molt and adult


#constant growth rates for Site and Class
U.Site=matrix(list("t1_BL","t1_BL","t1_DE","t1_DE","t1_DP","t1_DP",
                         "t1_PRH","t1_PRH","t1_TB","t1_TB","t1_TP","t1_TP"), 
               nrow = 12, ncol = 1,
               byrow = TRUE)# separate growth rates for molt and adult

#constant growth rates for BL, S, PRH, N. Adult and Molt
U.4Site.Class=matrix(list("t1_BL_A","t1_BL_M","t1_S_A","t1_S_M","t1_S_A","t1_S_M",
                         "t1_PRH_A","t1_PRH_M","t1_N_A","t1_N_M","t1_N_A","t1_N_M"), 
                    nrow = 12, ncol = 1,
                    byrow = TRUE)# separate growth rates for molt and adult

#constant growth rates for BL, S, PRH, N. Adult and Molt
U.4Site=matrix(list("t1_BL","t1_BL","t1_S","t1_S","t1_S","t1_S",
                          "t1_PRH","t1_PRH","t1_N","t1_N","t1_N","t1_N"), 
                     nrow = 12, ncol = 1,
                     byrow = TRUE)# separate growth rates for molt and adult



## some time varying U ---------------

# all sites have same growth rates that can change around 2003
U1 <- matrix("t1", 12, 1)
U2 <- matrix("t2", 12, 1)
Ut.All <- array(U2, dim = c(dim(U1), dim(dat)[2]))
TT <- dim(dat)[2]
Ut.All[, , 1:ceiling(TT / 2)] <- U1

#timevarying u by site
U1 <- matrix(c("t1_BL","t1_BL","t1_DE","t1_DE","t1_DP","t1_DP",
               "t1_PRH","t1_PRH","t1_TB","t1_TB","t1_TP","t1_TP"),12, 1)
U2 <- matrix(c("t2_BL","t2_BL","t2_DE","t2_DE","t2_DP","t2_DP",
               "t2_PRH","t2_PRH","t2_TB","t2_TB","t2_TP","t2_TP"),12, 1)
Ut.Site <- array(U2, dim = c(dim(U1), dim(dat)[2]))
TT <- dim(dat)[2]
Ut.Site[, , 1:ceiling(TT / 2)] <- U1

#timevarying u by site and class
U1 <- matrix(c("t1_BL_A","t1_BL_M","t1_DE_A","t1_DE_M","t1_DP_A","t1_DP_M",
               "t1_PRH_A","t1_PRH_M","t1_TB_A","t1_TB_M","t1_TP_A","t1_TP_M"),12, 1)
U2 <- matrix(c("t2_BL_A","t2_BL_M","t2_DE_A","t2_DE_M","t2_DP_A","t2_DP_M",
               "t2_PRH_A","t2_PRH_M","t2_TB_A","t2_TB_M","t2_TP_A","t2_TP_M"),12, 1)
Ut.Site.Class <- array(U2, dim = c(dim(U1), dim(dat)[2]))
TT <- dim(dat)[2]
Ut.Site.Class[, , 1:ceiling(TT / 2)] <- U1

#timevarying u by 4 sites
#Many warnings on early runs
U1 <- matrix(c("t1_BL_A","t1_BL_M","t1_South_A","t1_South_M","t1_South_A","t1_South_M",
               "t1_PRH_A","t1_PRH_M","t1_North_A","t1_North_M","t1_North_A","t1_North_M"),12, 1)
U2 <- matrix(c("t2_BL_A","t2_BL_M","t2_South_A","t2_South_M","t2_South_A","t2_South_M",
               "t2_PRH_A","t2_PRH_M","t2_North_A","t2_North_M","t2_North_A","t2_North_M"),12, 1)
Ut.4.Site.Class <- array(U2, dim = c(dim(U1), dim(dat)[2]))
TT <- dim(dat)[2]
Ut.4.Site.Class[, , 1:ceiling(TT / 2)] <- U1

# timevarying u by 4 sites Adult and Molt combined
# causes R to crash!!!! 
U1 <- matrix(c("t1_BL","t1_BL","t1_South","t1_South","t1_South","t1_South",
               "t1_PRH","t1_PRH","t1_North","t1_North","t1_North","t1_North"),12, 1)
U2 <- matrix(c("t2_BL","t2_BL","t2_South","t2_South","t2_South","t2_South",
               "t2_PRH","t2_PRH","t2_North","t2_North","t2_North","t2_North"),12, 1)
Ut.4.Site <- array(U2, dim = c(dim(U1), dim(dat)[2]))
TT <- dim(dat)[2]
Ut.4.Site[, , 1:ceiling(TT / 2)] <- U1

## R model options
R.model <- diag(0.05, 12) #known observation error variance
R.model="diagonal and equal"  #zero since focusing on B  else  "diagonal and equal"

## C model matrix for UI
C.model.UI=matrix(list("BL_A",0,0,0,0,0,  "UI_A","UI_A_lag",  #adding UI and UI_lag
                    "BL_M",0,0,0,0,0,  "UI_M","UI_M_lag",
                    #  "BL_P",0,0,0,0," PDO_P",
                                  
                    0,"DE_A",0,0,0,0,  "UI_A","UI_A_lag",
                    0,"DE_M",0,0,0,0,  "UI_M","UI_M_lag",
                    #   0,"DE_P",0,0,0," PDO_P",
                          
                    0,0,"DP_A",0,0,0,  "UI_A","UI_A_lag",
                    0,0,"DP_M",0,0,0,  "UI_M","UI_M_lag",
                    #   0,0,"DP_P",0,0," PDO_P",
                                  
                    # 0,0,0,0,0,0,0,0,"PDO_A",
                    # 0,0,0,0,0,0,0,0,"PDO_M",       #DR
                    #   0,0,0,0,0," PDO_P",
                                  
                    # 0,0,0,0,0,0,0,0,"PDO_A",       #PB
                    # 0,0,0,0,0,0,0,0,"PDO_M",
                    #   0,0,0,0,0," PDO_P"
                                  
                    0,0,0,0,0,0,  "UI_A", "UI_A_lag",        #PRH
                    0,0,0,0,0,0,  "UI_M", "UI_M_lag",
                                
                    0,0,0,0,0,0,  "UI_A", "UI_A_lag",       #TB
                    0,0,0,0,0,0,  "UI_M", "UI_M_lag",
                                  
                    0,0,0,0,0,0,  "UI_A", "UI_A_lag",       #TP
                    0,0,0,0,0,0,  "UI_M", "UI_M_lag"
),
nrow = 12, ncol = 8,
byrow = TRUE)

# C model for PDO (no lag)
C.model.PDO=matrix(list("BL_A",0,0,0,0,0,"PDO_A",
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

#2024-05-27  SETUP MODELS NOW BASED ON ABOVE parameters ---------------------------

## MARSS Models----------------------------------------
# AIC table setup
df_aic <- data.frame(model=character(), aic=integer())


## U models static
# equal 
# U.Class
# U.Site
# U.4Site.Class
# U.4Site


##m01 - U equal
t0 <- Sys.time()
m01b.1u.QR.fix=MARSS(dat, model=list(
                       Z=factor(c(1:12)), 
                       U="equal", 
                       R=diag(0.025, 12),
                       Q="diagonal and equal",
                       B=B.model,
                       C=C.model.UI,
                       #x0 = x0.model, 
                       tinitx=1, 
                       c = small_c_Coyote_3yr_UI_UI_lag),
                       control=list(maxit=5000, safe=TRUE, trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "m01b.1u.QR.fix", aic = m01b.1u.QR.fix$AIC)
df_aic
beepr::beep()
t1 <- Sys.time()
t1-t0 #run time #1.6 min
save(m01b.1u.QR.fix, file = "Output/m01b.1u.QR.fix.RData")
#load(file = "Output/m01b.1u.QR.fix.RData")


## m02 - U class
t0 <- Sys.time()
m02.u.Class=MARSS(dat, model=list(
  Z=factor(c(1:12)), 
  U=U.Class, 
  R=diag(0.025, 12),
  Q="diagonal and equal",
  B=B.model,
  C=C.model.UI,
  #x0 = x0.model, 
  tinitx=1, 
  c = small_c_Coyote_3yr_UI_UI_lag),
  control=list(maxit=5000, safe=TRUE, trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "m02.u.Class", aic = m02.u.Class$AIC)
df_aic
beepr::beep()
t1 <- Sys.time()
t1-t0 #run time 1.7 min
save(m02.u.Class, file = "Output/m02.u.Class.RData")
#load(file = "Output/m02.u.Class.RData")


## m03 - U Site 
## not converged after 5000 iter.
t0 <- Sys.time()
m03.u.Site=MARSS(dat, model=list(
  Z=factor(c(1:12)), 
  U=U.Site, 
  R=diag(0.025, 12),
  Q="diagonal and equal",
  B=B.model,
  C=C.model.UI,
  #x0 = x0.model, 
  tinitx=1, 
  c = small_c_Coyote_3yr_UI_UI_lag),
  control=list(maxit=5000, safe=TRUE, trace = 0, allow.degen=TRUE)) 
df_aic <- df_aic %>% add_row(model = "m03.u.Site", aic = m03.u.Site$AIC)
df_aic
beepr::beep()
t1 <- Sys.time()
t1-t0 #run time  9 min
save(m03.u.Site, file = "Output/m03.u.Site.RData")
#load(file = "Output/m03.u.Site.RData")

summary(m03.u.Site)


## m04 u 4site / class
## not converged after 5000 iter
t0 <- Sys.time()
m04.u.4Site.Class=MARSS(dat, model=list(
  Z=factor(c(1:12)), 
  U=U.4Site.Class, 
  R=diag(0.025, 12),
  Q="diagonal and equal",
  B=B.model,
  C=C.model.UI,
  #x0 = x0.model, 
  tinitx=1, 
  c = small_c_Coyote_3yr_UI_UI_lag),
  control=list(maxit=5000, safe=TRUE, trace = 0, allow.degen=TRUE)) 
df_aic <- df_aic %>% add_row(model = "m04.u.4Site.Class", aic = m04.u.4Site.Class$AIC)
df_aic
beepr::beep()
t1 <- Sys.time()
t1-t0 #run time 8.9 min
save(m04.u.4Site.Class, file = "Output/m04.u.4Site.Class.RData")
#load(file = "Output/m04.u.4Site.Class.RData")


## m05 u 4site
## not converged after 5000 iter
t0 <- Sys.time()
m05.u.4Site=MARSS(dat, model=list(
  Z=factor(c(1:12)), 
  U=U.4Site, 
  R=diag(0.025, 12),
  Q="diagonal and equal",
  B=B.model,
  C=C.model.UI,
  #x0 = x0.model, 
  tinitx=1, 
  c = small_c_Coyote_3yr_UI_UI_lag),
  control=list(maxit=5000, safe=TRUE, trace = 0, allow.degen=TRUE)) 
df_aic <- df_aic %>% add_row(model = "m05.u.4Site", aic = m05.u.4Site$AIC)
df_aic
beepr::beep()
t1 <- Sys.time()
t1-t0 #run time 8.9  min
save(m05.u.4Site, file = "Output/m05.u.4Site.RData")
#load(file = "Output/m05.u.4Site.RData")


## TV models
# Ut.All
# Ut.Site
# Ut.Site.Class
# Ut.4.Site.Class
# Ut.4.Site

##m06 - Ut.All
t0 <- Sys.time()
m06.ut.All=MARSS(dat, model=list(
  Z=factor(c(1:12)), 
  U=Ut.All, 
  R=diag(0.025, 12),
  Q="diagonal and equal",
  B=B.model,
  C=C.model.UI,
  #x0 = x0.model, 
  tinitx=1, 
  c = small_c_Coyote_3yr_UI_UI_lag),
  control=list(maxit=5000, safe=TRUE, trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "m06.ut.All", aic = m06.ut.All$AIC)
df_aic
beepr::beep()
t1 <- Sys.time()
t1-t0 #run time 1.6 min
save(m06.ut.All, file = "Output/m06.ut.All.RData")
#load(file = "Output/m06.ut.All.RData")


##m07 - Ut.Site
## warnings about Q or R going to zero.  See page 305.  Only effects loglik and AIC, but model still robust-ish
t0 <- Sys.time()
m07.Ut.Site=MARSS(dat, model=list(
  Z=factor(c(1:12)), 
  U=Ut.Site, 
  R=diag(0.025, 12),
  Q="diagonal and equal",
  B=B.model,
  C=C.model.UI,
  #x0 = x0.model, 
  tinitx=1, 
  c = small_c_Coyote_3yr_UI_UI_lag),
  control=list(maxit=5000, safe=TRUE, trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "m07.Ut.Site", aic = m07.Ut.Site$AIC)
df_aic
beepr::beep()
t1 <- Sys.time()
t1-t0 #run time 11  min
save(m07.Ut.Site, file = "Output/m07.Ut.Site.RData")
#load(file = "Output/m06.ut.All.RData")


##m08 - Ut.Site.Class
t0 <- Sys.time()
m08.Ut.Site.Class=MARSS(dat, model=list(
  Z=factor(c(1:12)), 
  U=Ut.Site.Class, 
  R=diag(0.025, 12),
  Q="diagonal and equal",
  B=B.model,
  C=C.model.UI,
  #x0 = x0.model, 
  tinitx=1, 
  c = small_c_Coyote_3yr_UI_UI_lag),
  control=list(maxit=5000, safe=TRUE, trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "m08.Ut.Site.Class", aic = m08.Ut.Site.Class$AIC)
df_aic
beepr::beep()
t1 <- Sys.time()
t1-t0 #run time #1.6 min
save(m08.Ut.Site.Class, file = "Output/m08.Ut.Site.Class.RData")
#load(file = "Output/m06.ut.All.RData")


##m09 - Ut.Site.Class
t0 <- Sys.time()
m09.Ut.4.Site.Class=MARSS(dat, model=list(
  Z=factor(c(1:12)), 
  U=Ut.4.Site.Class, 
  R=diag(0.025, 12),
  Q="diagonal and equal",
  B=B.model,
  C=C.model.UI,
  #x0 = x0.model, 
  tinitx=1, 
  c = small_c_Coyote_3yr_UI_UI_lag),
  control=list(maxit=5000, safe=TRUE, trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "m09.Ut.4.Site.Class", aic = m09.Ut.4.Site.Class$AIC)
df_aic
beepr::beep()
t1 <- Sys.time()
t1-t0 #run time #1.6 min
save(m09.Ut.4.Site.Class, file = "Output/m09.Ut.4.Site.Class.RData")
#load(file = "Output/m06.ut.All.RData")




























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
m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U=MARSS(dat, model=list(Z=Z.model, 
                                                           U=U.model,
                                                           Q=Q.model, 
                                                           R=R.model, 
                                                           B=B.model,
                                                           C=C.model,
                                                           #x0 = x0.model, 
                                                           tinitx=1, 
                                                           c=small_c_Coyote_01_PDO_MAR),
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
m.Ind_Molt_Adult_Coyote_PDO_B_equalcov_1U_TV_Site=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,

m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_Site=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
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

df_aic <- df_aic %>% add_row(model = "Model_1_J - m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_Site", aic = m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_Site$AIC)
df_aic
beepr::beep()
save(m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_Site, file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_Site.RData")
#load(file = "Output/m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_Site.RData")


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


C.model=matrix(              list("BL_A",0,0,0,0,0,"UI_A","UI_A_lag",  #adding UI and UI_lag
                                  "BL_M",0,0,0,0,0,"UI_M","UI_M_lag",
                                  #  "BL_P",0,0,0,0," PDO_P",
                                  
                                  0,"DE_A",0,0,0,0,"UI_A","UI_A_lag",
                                  0,"DE_M",0,0,0,0,"UI_M","UI_M_lag",
                                  #   0,"DE_P",0,0,0," PDO_P",
                                  
                                  0,0,"DP_A",0,0,0,"UI_A","UI_A_lag",
                                  0,0,"DP_M",0,0,0,"UI_M","UI_M_lag",
                                  #   0,0,"DP_P",0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",
                                  # 0,0,0,0,0,0,0,0,"PDO_M",       #DR
                                  #   0,0,0,0,0," PDO_P",
                                  
                                  # 0,0,0,0,0,0,0,0,"PDO_A",       #PB
                                  # 0,0,0,0,0,0,0,0,"PDO_M",
                                  #   0,0,0,0,0," PDO_P"
                                  
                                  0,0,0,0,0,0,"UI_A", "UI_A_lag",        #PRH
                                  0,0,0,0,0,0,"UI_M", "UI_M_lag",
                                  
                                  0,0,0,0,0,0,"UI_A", "UI_A_lag",       #TB
                                  0,0,0,0,0,0,"UI_M", "UI_M_lag",
                                  
                                  0,0,0,0,0,0,"UI_A", "UI_A_lag",       #TP
                                  0,0,0,0,0,0,"UI_M", "UI_M_lag"
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

#load(file = "Output/m.Ind_Molt_Adult_Coyote_UI_B_ident_1U_TV_R_fix.RData")

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

R.model <-  "diagonal and unequal" # observation error variance

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

