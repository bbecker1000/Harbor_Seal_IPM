library("MARSS")
library(ggplot2)
library("readxl")
library(dplyr)
library(tidyverse)
library(sjPlot) #used at end for 1 plot theme

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
#R.model="equal"  #zero since focusing on B  else  "diagonal and equal"
U.model=Ut  #Select from above
Q.model="equalvarcov" #"diagonal and equal"   #near zero add control
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
m.Ind_Molt_Adult_Coyote_PDO_B_unc_1U_TV_R_fix=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model,
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
#Z.model=factor(rep(1, 12))
R.model <- diag(0.05, 12) #known observation error variance
R.model="diagonal and equal"
U.model=Ut  #Select from above
Q.model="equalvarcov"#"diagonal and equal"   #near zero add control
B.model="equalvarcov"#unconstrained"  # > 90 min
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



#############################
# Model plots -------------------------------------
#############################
BESTMODEL <- m.Ind_Molt_Adult_Coyote_PDO_B_equalcov_1U_TV_1pop # m.1pop_Molt_Adult_Coyote_PDO_B_custom_1U_R0 # # m.1pop_Molt_Adult_Coyote_PDO_B_custom_1U_R0 # # m.1pop_Molt_Adult_Coyote_PDO_B_custom_1U# m.1pop_Molt_Adult_Coyote_PDO_B_unc_2U #m.1pop_Coyote_PDO_Xo_fixed_B_unc_tinitx_1 #m.1pop_Coyote_PDO_Xo_fixed_B_unc# m.1pop_Coyote_PDO_B_unc #   m.1pop_Coyote_PDO_B_unc   m.5pop_Coyote_PDO_B_unc

autoplot(BESTMODEL)

autoplot(BESTMODEL, plot.type = "fitted.ytT") + # xtT
  ylim(4,8) # xtT

autoplot(BESTMODEL, plot.type = "xtT") + # xtT
  ylim(3.1,7.6)

# plot temporal change from starting values

#BESTMODEL <- m.MOLT_IND_SITE_MOCI.equal_B_diag_uneq
CIs <- MARSSparamCIs(BESTMODEL, alpha = 0.11)  #\
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
  theme_gray(base_size = 20) +
  ylab("index of log abundance (N-No)") +
  xlab("Year") +
  facet_wrap(.~Season, ncol = 1)

####-----plot covariate effects

coef.data <- tidy(CIs)

#get just the coefficient terms
coef.data <- coef.data %>%
  filter(str_detect(term, "^C."))

coef.data$Season <- ifelse(str_detect(coef.data$term, "_A"), "Breeding Adults", 
                           ifelse(str_detect(coef.data$term, "_M"), "Molting", "Pups"))
#coef.data$Covariate <- ifelse(str_detect(coef.data$term, "MOCI"), "MOCI", "Coyote")
coef.data$Covariate <- ifelse(str_detect(coef.data$term, "PDO"), "PDO", "Coyote")
  
ggplot(coef.data, aes(term, estimate, color = Season, shape = Season)) +

  geom_pointrange(aes(ymin = conf.low, ymax = conf.up), size = 0.8) +
  coord_flip() +
  theme(legend.position = c(0.2, 0.2)) + 
  ylim(-1.25, 0.25) +

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
forecast1

autoplot(forecast1) +
  theme_grey(base_size = 18) +
  ylim(3,7.6)

