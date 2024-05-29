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




## 2024-05-26 HYPOTHESES for Breed and Molt

## Z hypotheses

# 12 Independent pops since want to understand growth rates...know different...
Z.model=factor(c(1:18))  # 6 sites x 2 age classes
 

# R1: diagonal and equal ...no reason to assume otherwise, standardized methods, repeat visits
R1.model="diagonal and equal"

# Q1: diagonal and equal ... no reason to assume otherwise
Q1.model="diagonal and equal"

# B1 design matrix based on predicted north and south colony interactions (Sarah Allen)
# B.models
#breed then molt, so only breed effect on molt "BB"
#breeding fidelity = AA
#molting fidelity = MM
#pup fidelity = PP
#adult effect pup = AP
#adult effect molt = AM
#A / B / C / D are within zone interactions by breed or molt

                     #A M P        #A M P   #A M P  #A  M  P #A  M  P  #A  M  P

B.model=matrix(list("AA",0,0,         "A",0,0,   "A",0,0,  0,0,0,   0,0,0,  0,0,0,  #A #BL
                    "AM","MM",0,      0,"A",0,   0,"A",0,  0,0,0,   0,0,0,  0,0,0,  #M
                    "AP",0,"PP",      0,0,"A",   0,0,"A",  0,0,0,   0,0,0,  0,0,0,  #P
                    
                    "A",0,0,  "AA",0,0,      "A",0,0,  "A",0,0,    0,0,0,  0,0,0, #DE
                    0,"A",0,  "AM","MM",0,   0,"A",0,  0,"A",0,    0,0,0,  0,0,0,
                    0,0,"A",  "AP",0,"PP",      0,0,"A",  0,0,"A",    0,0,0,  0,0,0,
                    
                    "A",0,0,  "A",0,0,   "AA",0,0,     "A",0,0,   0,0,0,  0,0,0,  #DP
                    0,"A",0,  0,"A",0,   "AM","MM",0,  0,"A",0,   0,0,0,  0,0,0,
                    0,0,"A",  0,0,"A",   "AP",0,"PP",  0,0,"A",    0,0,0,  0,0,0,
                    
                    "B",0,0,  "B",0,0,   "B",0,0,  "AA",0,0,      "C",0,0,  "C",0,0,  #PRH
                    0,"B",0,  0,"B",0,    0,"B",0,  "AM","MM",0,   0,"C",0,  0,"C",0,
                    0,0,"B",   0,0,"B",   0,0,"B",  "AP",0,"PP",   0,0,"C",  0,0,"C",
                    
                    0,0,0,   0,0,0,     0,0,0,   "D",0,0,   "AA",0,0,      "D",0,0,  #TB
                    0,0,0,   0,0,0,     0,0,0,   0,"D",0,   "AM","MM",0,   0,"D",0,
                    0,0,0,   0,0,0,     0,0,0,   0,0,"D",   "AP",0,"PP",   0,0,"D",  
                    
                    0,0,0,   0,0,0,     0,0,0,   0,0,0,    "D",0,0,   "AA",0,0,  #TP
                    0,0,0,   0,0,0,     0,0,0,   0,0,0,    0,"D",0,   "AM","MM",0,
                    0,0,0,   0,0,0,     0,0,0,   0,0,0,    0,0,"D",   "AP",0,"PP"),
               
               nrow = 18, ncol = 18,
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
U.Class=matrix(list("UA", "UM", "UP",
                    "UA", "UM", "UP",
                    "UA", "UM", "UP",
                    "UA", "UM", "UP",
                    "UA", "UM", "UP",
                    "UA", "UM", "UP"), 
               nrow = 18, ncol = 1,
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

# all sites have same growth rates that can change around 2003-4
U1 <- matrix("t1", 18, 1)
U2 <- matrix("t2", 18, 1)
Ut.All <- array(U2, dim = c(dim(U1), dim(dat)[2]))
TT <- dim(dat)[2]
Ut.All[, , 1:ceiling(TT / 4)] <- U1   #t1 through 2004

#timevarying u by site
U1 <- matrix(c("t1_BL","t1_BL","t1_BL",
               "t1_DE","t1_DE","t1_DE",
               "t1_DP","t1_DP","t1_DP",
               "t1_PRH","t1_PRH","t1_PRH",
               "t1_TB","t1_TB","t1_TB",
               "t1_TP","t1_TP","t1_TP"),18, 1)
U2 <- matrix(c("t2_BL","t2_BL","t2_BL",
               "t2_DE","t2_DE","t2_DE",
               "t2_DP","t2_DP","t2_DP",
               "t2_PRH","t2_PRH","t2_PRH",
               "t2_TB","t2_TB","t2_TB",
               "t2_TP","t2_TP","t2_TP"),18, 1)
Ut.Site <- array(U2, dim = c(dim(U1), dim(dat)[2]))
TT <- dim(dat)[2]
Ut.Site[, , 1:ceiling(TT / 4)] <- U1    #t1 through 2004

#timevarying u by site and class
U1 <- matrix(c("t1_BL_A","t1_BL_M","t1_BL_P",
               "t1_DE_A","t1_DE_M","t1_DE_P",
               "t1_DP_A","t1_DP_M","t1_DP_P",
               "t1_PRH_A","t1_PRH_M","t1_PRH_P",
               "t1_TB_A","t1_TB_M","t1_TB_P",
               "t1_TP_A","t1_TP_M","t1_TP_P"),18, 1)
U2 <- matrix(c("t2_BL_A","t2_BL_M","t2_BL_P",
               "t2_DE_A","t2_DE_M","t2_DE_P",
               "t2_DP_A","t2_DP_M","t2_DP_P",
               "t2_PRH_A","t2_PRH_M","t2_PRH_P",
               "t2_TB_A","t2_TB_M","t2_TB_P",
               "t2_TP_A","t2_TP_M","t2_TP_P"),18, 1)
Ut.Site.Class <- array(U2, dim = c(dim(U1), dim(dat)[2]))
TT <- dim(dat)[2]
Ut.Site.Class[, , 1:ceiling(TT / 4)] <- U1 #t1 through 2003

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
                       "BL_P",0,0,0,0,0,  "UI_P","UI_P_lag",
                                  
                       0,"DE_A",0,0,0,0,  "UI_A","UI_A_lag",
                       0,"DE_M",0,0,0,0,  "UI_M","UI_M_lag",
                       0,"DE_P",0,0,0,0,  "UI_P","UI_P_lag",
                          
                       0,0,"DP_A",0,0,0,  "UI_A","UI_A_lag",
                       0,0,"DP_M",0,0,0,  "UI_M","UI_M_lag",
                       0,0,"DP_P",0,0,0,  "UI_P","UI_P_lag",
                                  
                    # 0,0,0,0,0,0,0,0,"PDO_A",
                    # 0,0,0,0,0,0,0,0,"PDO_M",       #DR
                    #   0,0,0,0,0," PDO_P",
                                  
                    # 0,0,0,0,0,0,0,0,"PDO_A",       #PB
                    # 0,0,0,0,0,0,0,0,"PDO_M",
                    #   0,0,0,0,0," PDO_P"
                                  
                       0,0,0,0,0,0,  "UI_A", "UI_A_lag",       #PRH
                       0,0,0,0,0,0,  "UI_M", "UI_M_lag",
                       0,0,0,0,0,0,  "UI_P", "UI_P_lag",
                                
                       0,0,0,0,0,0,  "UI_A", "UI_A_lag",       #TB
                       0,0,0,0,0,0,  "UI_M", "UI_M_lag",
                       0,0,0,0,0,0,  "UI_P", "UI_P_lag",
                                  
                       0,0,0,0,0,0,  "UI_A", "UI_A_lag",       #TP
                       0,0,0,0,0,0,  "UI_M", "UI_M_lag",
                       0,0,0,0,0,0,  "UI_P", "UI_P_lag"
),
nrow = 18, ncol = 8,
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

.1997.2022.
##m.1997.2022.01 - U equal
t0 <- Sys.time()
m.1997.2022.01b.1u.QR.fix=MARSS(dat, model=list(
                       Z=factor(c(1:18)), 
                       U="equal", 
                       R=diag(0.025, 18),
                       Q="diagonal and equal",
                       B=B.model,
                       C=C.model.UI,
                       #x0 = x0.model, 
                       tinitx=1, 
                       c = small_c_Coyote_3yr_UI_UI_lag),
                       control=list(maxit=5000, safe=TRUE, trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "m.1997.2022.01b.1u.QR.fix", aic = m.1997.2022.01b.1u.QR.fix$AIC)
df_aic
beepr::beep()
t1 <- Sys.time()
t1-t0 #run time #1.6 min
save(m.1997.2022.01b.1u.QR.fix, file = "Output/m.1997.2022.01b.1u.QR.fix.RData")
#load(file = "Output/m.1997.2022.01b.1u.QR.fix.RData")


## m.1997.2022.02 - U class
t0 <- Sys.time()
m.1997.2022.02.u.Class=MARSS(dat, model=list(
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

df_aic <- df_aic %>% add_row(model = "m.1997.2022.02.u.Class", aic = m.1997.2022.02.u.Class$AIC)
df_aic
beepr::beep()
t1 <- Sys.time()
t1-t0 #run time 1.7 min
save(m.1997.2022.02.u.Class, file = "Output/m.1997.2022.02.u.Class.RData")
#load(file = "Output/m.1997.2022.02.u.Class.RData")


## m.1997.2022.03 - U Site 
## not converged after 5000 iter.
t0 <- Sys.time()
m.1997.2022.03.u.Site=MARSS(dat, model=list(
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
df_aic <- df_aic %>% add_row(model = "m.1997.2022.03.u.Site", aic = m.1997.2022.03.u.Site$AIC)
df_aic
beepr::beep()
t1 <- Sys.time()
t1-t0 #run time  9 min
save(m.1997.2022.03.u.Site, file = "Output/m.1997.2022.03.u.Site.RData")
#load(file = "Output/m.1997.2022.03.u.Site.RData")

summary(m.1997.2022.03.u.Site)


## m.1997.2022.04 u 4site / class
## not converged after 5000 iter
t0 <- Sys.time()
m.1997.2022.04.u.4Site.Class=MARSS(dat, model=list(
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
df_aic <- df_aic %>% add_row(model = "m.1997.2022.04.u.4Site.Class", aic = m.1997.2022.04.u.4Site.Class$AIC)
df_aic
beepr::beep()
t1 <- Sys.time()
t1-t0 #run time 8.9 min
save(m.1997.2022.04.u.4Site.Class, file = "Output/m.1997.2022.04.u.4Site.Class.RData")
#load(file = "Output/m.1997.2022.04.u.4Site.Class.RData")


## m.1997.2022.05 u 4site
## not converged after 5000 iter
t0 <- Sys.time()
m.1997.2022.05.u.4Site=MARSS(dat, model=list(
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
df_aic <- df_aic %>% add_row(model = "m.1997.2022.05.u.4Site", aic = m.1997.2022.05.u.4Site$AIC)
df_aic
beepr::beep()
t1 <- Sys.time()
t1-t0 #run time 8.9  min
save(m.1997.2022.05.u.4Site, file = "Output/m.1997.2022.05.u.4Site.RData")
#load(file = "Output/m.1997.2022.05.u.4Site.RData")


## TV models
# Ut.All
# Ut.Site
# Ut.Site.Class
# Ut.4.Site.Class
# Ut.4.Site

##m.1997.2022.06 - Ut.All
t0 <- Sys.time()
m.1997.2022.06.ut.All=MARSS(dat, model=list(
  Z=factor(c(1:18)), 
  U=Ut.All, 
  R=diag(0.025, 18),
  Q="diagonal and equal",
  B=B.model,
  C=C.model.UI,
  #x0 = x0.model, 
  tinitx=1, 
  c = small_c_Coyote_3yr_UI_UI_lag),
  control=list(maxit=5000, safe=TRUE, trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "m.1997.2022.06.ut.All", aic = m.1997.2022.06.ut.All$AIC)
df_aic
beepr::beep()
t1 <- Sys.time()
t1-t0 #run time 1.6 min
save(m.1997.2022.06.ut.All, file = "Output/m.1997.2022.06.ut.All.RData")
#load(file = "Output/m.1997.2022.06.ut.All.RData")


##m.1997.2022.07 - Ut.Site
## warnings about Q or R going to zero.  See page 305.  Only effects loglik and AIC, but model still robust-ish
t0 <- Sys.time()
m.1997.2022.07.Ut.Site=MARSS(dat, model=list(
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

df_aic <- df_aic %>% add_row(model = "m.1997.2022.07.Ut.Site", aic = m.1997.2022.07.Ut.Site$AIC)
df_aic
beepr::beep()
t1 <- Sys.time()
t1-t0 #run time 11  min
save(m.1997.2022.07.Ut.Site, file = "Output/m.1997.2022.07.Ut.Site.RData")
#load(file = "Output/m.1997.2022.06.ut.All.RData")


##m.1997.2022.08 - Ut.Site.Class
t0 <- Sys.time()
m.1997.2022.08.Ut.Site.Class=MARSS(dat, model=list(
  Z=factor(c(1:18)), 
  U=Ut.Site.Class, 
  R=diag(0.025, 18),
  Q="diagonal and equal",
  B=B.model,
  C=C.model.UI,
  #x0 = x0.model, 
  tinitx=1, 
  c = small_c_Coyote_3yr_UI_UI_lag),
  control=list(maxit=5000, safe=TRUE, trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "m.1997.2022.08.Ut.Site.Class", aic = m.1997.2022.08.Ut.Site.Class$AIC)
df_aic
beepr::beep()
t1 <- Sys.time()
t1-t0 #run time #1.6 min
save(m.1997.2022.08.Ut.Site.Class, file = "Output/m.1997.2022.08.Ut.Site.Class.RData")
#load(file = "Output/m.1997.2022.06.ut.All.RData")


##m.1997.2022.09 - Ut.Site.Class
t0 <- Sys.time()
m.1997.2022.09.Ut.4.Site.Class=MARSS(dat, model=list(
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

df_aic <- df_aic %>% add_row(model = "m.1997.2022.09.Ut.4.Site.Class", aic = m.1997.2022.09.Ut.4.Site.Class$AIC)
df_aic
beepr::beep()
t1 <- Sys.time()
t1-t0 #run time #1.6 min
save(m.1997.2022.09.Ut.4.Site.Class, file = "Output/m.1997.2022.09.Ut.4.Site.Class.RData")
#load(file = "Output/m.1997.2022.06.ut.All.RData")












