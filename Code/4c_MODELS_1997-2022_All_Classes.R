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

# 2024-06-04
# Add to the network matrix.
#    Hypotheses to add to matrix.
#      A. TB-TP  "D"--> DONE
#      B. DE-DP  "G" --> DONE
#      C. BL-DP   "F"  --> but try rerunning with "F" also in the BL--> DP effect

                     #A M P        #A M P   #A M P  #A  M  P #A  M  P  #A  M  P
# fixed DP-BL matric 2024-06-07

B.model=matrix(list("AA",0,0,     "A",0,0,    "F",0,0,  0,0,0,   0,0,0,  0,0,0,  #A #BL
                    "AM","MM",0,  0,"A",0,    0,"F",0,  0,0,0,   0,0,0,  0,0,0,  #M
                    "AP",0,"PP",  0,0,"A",    0,0,"F",  0,0,0,   0,0,0,  0,0,0,  #P
                    
                    "A",0,0,      "AA",0,0,        "G",0,0,  0,0,0,    0,0,0,  0,0,0, #DE
                    0,"A",0,      "AM","MM",0,     0,"G",0,  0,0,0,    0,0,0,  0,0,0,
                    0,0,"A",      "AP",0,"PP",      0,0,"G",  0,0,0,    0,0,0,  0,0,0,
                    
                    0,0,0,      "G",0,0,    "AA",0,0,     "A",0,0,   0,0,0,  0,0,0,  #DP
                    0,0,0,      0,"G",0,    "AM","MM",0,  0,"A",0,   0,0,0,  0,0,0,
                    0,0,0,      0,0,"G",    "AP",0,"PP",  0,0,"A",    0,0,0,  0,0,0,
                    
                    "B",0,0,      "B",0,0,   "B",0,0,  "AA",0,0,      "E",0,0,  "C",0,0,  #PRH
                    0,"B",0,      0,"B",0,    0,"B",0,  "AM","MM",0,   0,"E",0,  0,"C",0,
                    0,0,"B",      0,0,"B",   0,0,"B",  "AP",0,"PP",   0,0,"E",  0,0,"C",
                    
                    0,0,0,       0,0,0,     0,0,0,   "E",0,0,   "AA",0,0,      "D",0,0,  #TB
                    0,0,0,       0,0,0,     0,0,0,   0,"E",0,   "AM","MM",0,   0,"D",0,
                    0,0,0,       0,0,0,     0,0,0,   0,0,"E",   "AP",0,"PP",   0,0,"D",  
                    
                    0,0,0,       0,0,0,     0,0,0,   0,0,0,    "D",0,0,   "AA",0,0,  #TP
                    0,0,0,       0,0,0,     0,0,0,   0,0,0,    0,"D",0,   "AM","MM",0,
                    0,0,0,       0,0,0,     0,0,0,   0,0,0,    0,0,"D",   "AP",0,"PP"),
               
               nrow = 18, ncol = 18,
               byrow = TRUE)

## COVARIATES



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

#timevarying u by class
#adults and molt time vary at 2004, pups do not vary
U1 <- matrix(c("t1_A","t1_M","t1_P",
               "t1_A","t1_M","t1_P",
               "t1_A","t1_M","t1_P",
               "t1_A","t1_M","t1_P",
               "t1_A","t1_M","t1_P",
               "t1_A","t1_M","t1_P"),18, 1)
U2 <- matrix(c("t2_A","t2_M","t1_P",
               "t2_A","t2_M","t1_P",
               "t2_A","t2_M","t1_P",
               "t2_A","t2_M","t1_P",
               "t2_A","t2_M","t1_P",
               "t2_A","t2_M","t1_P"),18, 1)
Ut.Class <- array(U2, dim = c(dim(U1), dim(dat)[2]))
TT <- dim(dat)[2]
Ut.Class[, , 1:ceiling(TT / 4)] <- U1 #t1 through 2003


#timevarying u by class
#adults and molt time vary at 2004, pups do not vary
U1 <- matrix(c("t1_BL_A","t1_BL_M","t1_BL_P",
               "t1_DE_A","t1_DE_M","t1_DE_P",
               "t1_DP_A","t1_DP_M","t1_DP_P",
               "t1_PRH_A","t1_PRH_M","t1_PRH_P",
               "t1_TB_A","t1_TB_M","t1_TB_P",
               "t1_TP_A","t1_TP_M","t1_TP_P"),18, 1)
U2 <- matrix(c("t2_BL_A","t2_BL_M","t1_BL_P",
               "t2_DE_A","t2_DE_M","t1_DE_P",
               "t2_DP_A","t2_DP_M","t1_DP_P",
               "t2_PRH_A","t2_PRH_M","t1_PRH_P",
               "t2_TB_A","t2_TB_M","t1_TB_P",
               "t2_TP_A","t2_TP_M","t1_TP_P"),18, 1)
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
C.model.MOCI=matrix(
  list("BL_A",0,0,0,0,0,"MOCI_A","MOCI_AMJ_A_lag", "MOCI_OND_A_lag",   0,     "BL_Dist",0,0,0,0,0,
       "BL_M",0,0,0,0,0,  "MOCI_M","MOCI_AMJ_M_lag", "MOCI_OND_M_lag", 0,     "BL_Dist",0,0,0,0,0,
       "BL_P",0,0,0,0,0,  "MOCI_P","MOCI_AMJ_P_lag", "MOCI_OND_P_lag", 0,     "BL_Dist",0,0,0,0,0,
                                  
       0,"DE_A",0,0,0,0,  "MOCI_A","MOCI_AMJ_A_lag", "MOCI_OND_A_lag", "ES",  0,"DE_Dist",0,0,0,0, 
       0,"DE_M",0,0,0,0,  "MOCI_M","MOCI_AMJ_M_lag", "MOCI_OND_M_lag", "ES",  0,"DE_Dist",0,0,0,0,
       0,"DE_P",0,0,0,0,  "MOCI_P","MOCI_AMJ_P_lag", "MOCI_OND_P_lag", "ES",  0,"DE_Dist",0,0,0,0,
                        
       0,0,"DP_A",0,0,0,  "MOCI_A","MOCI_AMJ_A_lag", "MOCI_OND_A_lag", "ES",  0,0,"DP_Dist",0,0,0,
       0,0,"DP_M",0,0,0,  "MOCI_M","MOCI_AMJ_M_lag", "MOCI_OND_M_lag", "ES",  0,0,"DP_Dist",0,0,0,
       0,0,"DP_P",0,0,0,  "MOCI_P","MOCI_AMJ_P_lag", "MOCI_OND_P_lag", "ES",  0,0,"DP_Dist",0,0,0,
                                  
                    # 0,0,0,0,0,0,0,0,"PDO_A",
                    # 0,0,0,0,0,0,0,0,"PDO_M",       #DR
                    #   0,0,0,0,0," PDO_P",
                                  
                    # 0,0,0,0,0,0,0,0,"PDO_A",       #PB
                    # 0,0,0,0,0,0,0,0,"PDO_M",
                    #   0,0,0,0,0," PDO_P"
                                  
       0,0,0,0,0,0,  "MOCI_A", "MOCI_AMJ_A_lag", "MOCI_OND_A_lag", "ES",     0,0,0,"PRH_Dist",0,0, #PRH
       0,0,0,0,0,0,  "MOCI_M", "MOCI_AMJ_M_lag", "MOCI_OND_M_lag", "ES",     0,0,0,"PRH_Dist",0,0,
       0,0,0,0,0,0,  "MOCI_P", "MOCI_AMJ_P_lag", "MOCI_OND_P_lag", "ES",     0,0,0,"PRH_Dist",0,0,
                                
       0,0,0,0,0,0,  "MOCI_A", "MOCI_AMJ_A_lag", "MOCI_OND_A_lag", 0,        0,0,0,0,"TB_Dist",0,#TB
       0,0,0,0,0,0,  "MOCI_M", "MOCI_AMJ_M_lag", "MOCI_OND_M_lag", 0,        0,0,0,0,"TB_Dist",0,
       0,0,0,0,0,0,  "MOCI_P", "MOCI_AMJ_P_lag", "MOCI_OND_P_lag", 0,        0,0,0,0,"TB_Dist",0,
                                  
       0,0,0,0,0,0,  "MOCI_A", "MOCI_AMJ_A_lag", "MOCI_OND_A_lag", 0,        0,0,0,0,0,"TP_Dist",#TP
       0,0,0,0,0,0,  "MOCI_M", "MOCI_AMJ_M_lag", "MOCI_OND_M_lag", 0,        0,0,0,0,0,"TP_Dist",
       0,0,0,0,0,0,  "MOCI_P", "MOCI_AMJ_P_lag", "MOCI_OND_P_lag", 0,        0,0,0,0,0,"TP_Dist"
),
nrow = 18, ncol = 16,
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

####----
## adding MOCI MODEL

##m.1997.2022.06 - Ut.All
t0 <- Sys.time()
m.1997.2022.06.ut.All.MOCI=MARSS(dat, model=list(
  Z=factor(c(1:18)), 
  U=Ut.All, 
  R=diag(0.025, 18),
  Q="diagonal and equal",
  B=B.model,
  C=C.model.UI,
  #x0 = x0.model, 
  tinitx=1, 
  c = small_c_Coyote_3yr_MOCI_MOCI_lag),
  control=list(maxit=5000, safe=TRUE, trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "m.1997.2022.06.ut.All.MOCI", aic = m.1997.2022.06.ut.All.MOCI$AIC)
df_aic
beepr::beep()
t1 <- Sys.time()
t1-t0 #run time 11 min sec.
save(m.1997.2022.06.ut.All.MOCI, file = "Output/m.1997.2022.06.ut.All.MOCI.RData")
#load(file = "Output/m.1997.2022.06.ut.All.MOCI.RData")


##m.1997.2023.06.ut.All.MOCI.ES  4.8 sec !!
## time vary all classes same
t0 <- Sys.time()
m.1997.2023.06.ut.All.MOCI.ES=MARSS(dat, model=list(
  Z=factor(c(1:18)), 
  U=Ut.All, 
  R=diag(0.025, 18),
  Q="diagonal and equal",
  B=B.model,
  C=C.model.MOCI,
  #x0 = x0.model, 
  tinitx=1, 
  c = small_c_Coyote_3yr_MOCI_MOCI_lag),
  control=list(maxit=5000, safe=TRUE, trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "m.1997.2023.06.ut.All.MOCI.ES", aic = m.1997.2023.06.ut.All.MOCI.ES$AIC)
df_aic
beepr::beep()
t1 <- Sys.time()
t1-t0 #run time 60 sec.
save(m.1997.2023.06.ut.All.MOCI.ES, file = "Output/m.1997.2023.06.ut.All.MOCI.ES.RData")
#load(file = "Output/m.1997.2023.06.ut.All.MOCI.ES.RData")


##m.1997.2023.06.ut.All.MOCI.ES  
# time vary adult and molt,  not pup
# 31 min
t0 <- Sys.time()
m.1997.2023.06.ut.Site.Class.MOCI.ES=MARSS(dat, model=list(
  Z=factor(c(1:18)), 
  U=Ut.Site.Class, 
  R=diag(0.025, 18),
  Q="diagonal and equal",
  B=B.model,
  C=C.model.MOCI,
  #x0 = x0.model, 
  tinitx=1, 
  c = small_c_Coyote_3yr_MOCI_MOCI_lag),
  control=list(maxit=7500, safe=TRUE, trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "m.1997.2023.06.ut.Site.Class.MOCI.ES", aic = m.1997.2023.06.ut.Site.Class.MOCI.ES$AIC)
df_aic
beepr::beep()
t1 <- Sys.time()
t1-t0
save(m.1997.2023.06.ut.Site.Class.MOCI.ES, file = "Output/m.1997.2023.06.ut.Site.Class.MOCI.ES.RData")
#load(file = "Output/m.1997.2023.06.ut.Site.Class.MOCI.ES.RData")



##m.1997.2023.06.ut.All.MOCI.ES  
# time vary adult and molt,  not pup
# no site u diffs
# added human disturbance by site 2024-06-23
t0 <- Sys.time()
m.1997.2023.06.ut.Class.MOCI.ES.dist=MARSS(dat, model=list(
  Z=factor(c(1:18)), 
  U=Ut.Class, 
  R=diag(0.025, 18),
  Q="diagonal and equal",
  B=B.model,
  C=C.model.MOCI,
  #x0 = x0.model, 
  tinitx=1, 
  c = small_c_Coyote_3yr_MOCI_MOCI_Dist_lag),
  control=list(maxit=1000, safe=TRUE, trace = 0, allow.degen=TRUE)) 

df_aic <- df_aic %>% add_row(model = "m.1997.2023.06.ut.Class.MOCI.ES.dist", aic = m.1997.2023.06.ut.Class.MOCI.ES.dist$AIC)
df_aic
beepr::beep()
t1 <- Sys.time()
t1-t0
save(m.1997.2023.06.ut.Class.MOCI.ES.dist, file = "Output/m.1997.2023.06.ut.Class.MOCI.ES.dist.RData")
#load(file = "Output/m.1997.2023.06.ut.Class.MOCI.ES.dist.RData")













