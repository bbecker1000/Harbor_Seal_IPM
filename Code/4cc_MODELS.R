#4cc_MODELS

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



# 1:  MOCI_JFM
# 2:  MOCI_AMJ
# 3:  MOCI_OND
# 4:  Dist_BL
# 5:  Dist_DE
# 6:  Dist_DP
# 7:  Dist_PRH
# 8:  Dist_TB
# 9:  Dist_TP
# 10: Coyote_BL
# 11: Coyote_DE
# 12: Coyote_DP
# 13: Coyote_PRH
# 14: Coyote_TB
# 15: Coyote_TP
# 16: eSeal_Sum_Imm_MaxCount



# C matrix: 18 rows (states) × 16 columns (covariates)
# Rows: BL_A, BL_M, BL_P, DE_A, DE_M, DE_P, DP_A, DP_M, DP_P, PRH_A, PRH_M, PRH_P, TB_A, TB_M, TB_P, TP_A, TP_M, TP_P

# ============================================================================
# MARSS MODEL WITH UPDATED DATA (2025)
# ============================================================================

library(MARSS)


#scaled covariates
cov_t_marss <- cov_t_scaled


# ----------------------------------------------------------------------------
# C MATRIX: 18 rows (states) × 16 columns (covariates)
# Class-specific effects, pooled across sites
# ----------------------------------------------------------------------------

# Covariate order in cov_t_marss:
# 1:MOCI_JFM, 2:MOCI_AMJ, 3:MOCI_OND, 4:Dist_BL, 5:Dist_DE, 6:Dist_DP, 
# 7:Dist_PRH, 8:Dist_TB, 9:Dist_TP, 10:Coyote_BL, 11:Coyote_DE, 12:Coyote_DP,
# 13:Coyote_PRH, 14:Coyote_TB, 15:Coyote_TP, 16:eSeal

C.model.new <- matrix(
  list(
    #MOCI_JFM       MOCI_AMJ       MOCI_OND       Dist_BL    Dist_DE    Dist_DP    Dist_PRH    Dist_TB    Dist_TP    Coy_BL     Coy_DE     Coy_DP     Coy_PRH    Coy_TB     Coy_TP     eSeal
    "MOCI_JFM_A",  "MOCI_AMJ_A",  "MOCI_OND_A",  "Dist_BL", 0,         0,         0,          0,         0,         "Coy_BL",  0,         0,         0,         0,         0,         0,        #BL_A
    "MOCI_JFM_M",  "MOCI_AMJ_M",  "MOCI_OND_M",  "Dist_BL", 0,         0,         0,          0,         0,         "Coy_BL",  0,         0,         0,         0,         0,         0,        #BL_M
    "MOCI_JFM_P",  "MOCI_AMJ_P",  "MOCI_OND_P",  "Dist_BL", 0,         0,         0,          0,         0,         "Coy_BL",  0,         0,         0,         0,         0,         0,        #BL_P
    
    "MOCI_JFM_A",  "MOCI_AMJ_A",  "MOCI_OND_A",  0,         "Dist_DE", 0,         0,          0,         0,         0,         "Coy_DE",  0,         0,         0,         0,         "ES-DE",     #DE_A
    "MOCI_JFM_M",  "MOCI_AMJ_M",  "MOCI_OND_M",  0,         "Dist_DE", 0,         0,          0,         0,         0,         "Coy_DE",  0,         0,         0,         0,         "ES-DE",     #DE_M
    "MOCI_JFM_P",  "MOCI_AMJ_P",  "MOCI_OND_P",  0,         "Dist_DE", 0,         0,          0,         0,         0,         "Coy_DE",  0,         0,         0,         0,         "ES-DE",     #DE_P
    
    "MOCI_JFM_A",  "MOCI_AMJ_A",  "MOCI_OND_A",  0,         0,         "Dist_DP", 0,          0,         0,         0,         0,         "Coy_DP",  0,         0,         0,         0,     #DP_A    #removed ES DP per SA on 2026-04-27
    "MOCI_JFM_M",  "MOCI_AMJ_M",  "MOCI_OND_M",  0,         0,         "Dist_DP", 0,          0,         0,         0,         0,         "Coy_DP",  0,         0,         0,         0,     #DP_M
    "MOCI_JFM_P",  "MOCI_AMJ_P",  "MOCI_OND_P",  0,         0,         "Dist_DP", 0,          0,         0,         0,         0,         "Coy_DP",  0,         0,         0,         0,     #DP_P
    
    "MOCI_JFM_A",  "MOCI_AMJ_A",  "MOCI_OND_A",  0,         0,         0,         "Dist_PRH", 0,         0,         0,         0,         0,         0,         0,         0,         "ES-PRH",     #PRH_A    #adding site COV for ES per SA on 2026-04-27
    "MOCI_JFM_M",  "MOCI_AMJ_M",  "MOCI_OND_M",  0,         0,         0,         "Dist_PRH", 0,         0,         0,         0,         0,         0,         0,         0,         "ES-PRH",     #PRH_M
    "MOCI_JFM_P",  "MOCI_AMJ_P",  "MOCI_OND_P",  0,         0,         0,         "Dist_PRH", 0,         0,         0,         0,         0,         0,         0,         0,         "ES-PRH",     #PRH_P
    
    "MOCI_JFM_A",  "MOCI_AMJ_A",  "MOCI_OND_A",  0,         0,         0,         0,          "Dist_TB", 0,         0,         0,         0,         0,         0,         0,         0,        #TB_A
    "MOCI_JFM_M",  "MOCI_AMJ_M",  "MOCI_OND_M",  0,         0,         0,         0,          "Dist_TB", 0,         0,         0,         0,         0,         0,         0,         0,        #TB_M
    "MOCI_JFM_P",  "MOCI_AMJ_P",  "MOCI_OND_P",  0,         0,         0,         0,          "Dist_TB", 0,         0,         0,         0,         0,         0,         0,         0,        #TB_P
    
    "MOCI_JFM_A",  "MOCI_AMJ_A",  "MOCI_OND_A",  0,         0,         0,         0,          0,         "Dist_TP", 0,         0,         0,         0,         0,         0,         0,        #TP_A
    "MOCI_JFM_M",  "MOCI_AMJ_M",  "MOCI_OND_M",  0,         0,         0,         0,          0,         "Dist_TP", 0,         0,         0,         0,         0,         0,         0,        #TP_M
    "MOCI_JFM_P",  "MOCI_AMJ_P",  "MOCI_OND_P",  0,         0,         0,         0,          0,         "Dist_TP", 0,         0,         0,         0,         0,         0,         0         #TP_P
  ),
  nrow = 18, ncol = 16, byrow = TRUE
)
# ----------------------------------------------------------------------------
# TIME-VARYING U MATRIX
# Adults and Molt time-vary at 2004, Pups do not vary
# ----------------------------------------------------------------------------

U1 <- matrix(c("t1_A", "t1_M", "t1_P",
               "t1_A", "t1_M", "t1_P",
               "t1_A", "t1_M", "t1_P",
               "t1_A", "t1_M", "t1_P",
               "t1_A", "t1_M", "t1_P",
               "t1_A", "t1_M", "t1_P"), 18, 1)

U2 <- matrix(c("t2_A", "t2_M", "t1_P",
               "t2_A", "t2_M", "t1_P",
               "t2_A", "t2_M", "t1_P",
               "t2_A", "t2_M", "t1_P",
               "t2_A", "t2_M", "t1_P",
               "t2_A", "t2_M", "t1_P"), 18, 1)

Ut.Class <- array(U2, dim = c(dim(U1), ncol(dat)))
TT <- ncol(dat)
Ut.Class[, , 1:ceiling(TT / 4)] <- U1  # t1 through ~2004

# Check breakpoint year
years <- 1997:(1997 + TT - 1)
breakpoint <- ceiling(TT / 4)
cat("Breakpoint at column", breakpoint, "= year", years[breakpoint], "\n")

# ----------------------------------------------------------------------------
# VERIFY DIMENSIONS
# ----------------------------------------------------------------------------

cat("dat:", dim(dat), "\n")
cat("cov_t_marss:", dim(cov_t_marss), "\n")
cat("C.model.new:", dim(C.model.new), "\n")
cat("Ut.Class:", dim(Ut.Class), "\n")

# ----------------------------------------------------------------------------
# RUN MODEL
# ----------------------------------------------------------------------------

t0 <- Sys.time()
m.2025.ut.Class.MOCI.ES <- MARSS(dat, model = list(
  Z = factor(1:18),
  U = Ut.Class,
  R = diag(0.025, 18),
  Q = "diagonal and equal",
  B = "identity",
  C = C.model.new,
  c = cov_t_marss,
  tinitx = 1
), control = list(maxit = 5000, safe = TRUE, trace = 0, allow.degen = TRUE))

t1 <- Sys.time()
cat("Run time:", t1 - t0, "\n")

# ----------------------------------------------------------------------------
# RESULTS
# ----------------------------------------------------------------------------

m.2025.ut.Class.MOCI.ES$AICc
m.2025.ut.Class.MOCI.ES

# Save
save(m.2025.ut.Class.MOCI.ES, file = "Output/m.2025.ut.Class.MOCI.ES.RData")
summary(m.2025.ut.Class.MOCI.ES)


# ============================================================================
# MODEL B: ESTUARINE vs OUTER COAST GROUPING
# ============================================================================

# Z matrix groups states into subpopulations
# Estuarine: BL, DE, TB
# Outer Coast: DP, PRH, TP

# Model A: All 18 states independent
Z.indep <- factor(1:18)

# Model B: Group by location type AND class (6 groups)
# Row order: BL_A, BL_M, BL_P, DE_A, DE_M, DE_P, DP_A, DP_M, DP_P, 
#            PRH_A, PRH_M, PRH_P, TB_A, TB_M, TB_P, TP_A, TP_M, TP_P

Z.estuary <- factor(c(
  "Est_A", "Est_M", "Est_P",   # BL  - Estuarine
  "Est_A", "Est_M", "Est_P",   # DE  - Estuarine
  "Out_A", "Out_M", "Out_P",   # DP  - Outer Coast
  "Out_A", "Out_M", "Out_P",   # PRH - Outer Coast
  "Est_A", "Est_M", "Est_P",   # TB  - Estuarine
  "Out_A", "Out_M", "Out_P"    # TP  - Outer Coast
))

# ----------------------------------------------------------------------------
# MODEL A: Independent populations (18 states)
# ----------------------------------------------------------------------------

t0 <- Sys.time()
m.A.indep <- MARSS(dat, model = list(
  Z = Z.indep,
  U = Ut.Class,
  R = diag(0.025, 18),
  Q = "diagonal and equal",
  B = "identity",
  C = C.model.new,
  c = cov_t_marss,
  tinitx = 1
), control = list(maxit = 5000, safe = TRUE, trace = 0, allow.degen = TRUE))
t1 <- Sys.time()
cat("Model A run time:", t1 - t0, "\n")

# Z matrix: 6 hidden states observed through 18 observations
# ============================================================================
# MODEL B: ESTUARINE vs OUTER COAST GROUPING
# C MATRIX: 6 rows (hidden states) × 16 columns (covariates)
# ============================================================================

C.model.estuary <- matrix(
  list(
    #MOCI_JFM       MOCI_AMJ       MOCI_OND       Dist_BL       Dist_DE       Dist_DP       Dist_PRH       Dist_TB       Dist_TP       Coy_BL        Coy_DE        Coy_DP        Coy_PRH       Coy_TB        Coy_TP        eSeal
    "MOCI_JFM_A",  "MOCI_AMJ_A",  "MOCI_OND_A",  "Dist_Est_A", "Dist_Est_A", 0,            0,             "Dist_Est_A", 0,            "Coy_Est_A",  "Coy_Est_A",  0,            0,            "Coy_Est_A",  0,            "ES_A",       #Est_A (BL, DE, TB)
    "MOCI_JFM_M",  "MOCI_AMJ_M",  "MOCI_OND_M",  "Dist_Est_M", "Dist_Est_M", 0,            0,             "Dist_Est_M", 0,            "Coy_Est_M",  "Coy_Est_M",  0,            0,            "Coy_Est_M",  0,            "ES_M",       #Est_M
    "MOCI_JFM_P",  "MOCI_AMJ_P",  "MOCI_OND_P",  "Dist_Est_P", "Dist_Est_P", 0,            0,             "Dist_Est_P", 0,            "Coy_Est_P",  "Coy_Est_P",  0,            0,            "Coy_Est_P",  0,            "ES_P",       #Est_P
    
    "MOCI_JFM_A",  "MOCI_AMJ_A",  "MOCI_OND_A",  0,            0,            "Dist_Out_A", "Dist_Out_A",  0,            "Dist_Out_A", 0,            0,            "Coy_Out_A",  "Coy_Out_A",  0,            "Coy_Out_A",  "ES_A",       #Out_A (DP, PRH, TP)
    "MOCI_JFM_M",  "MOCI_AMJ_M",  "MOCI_OND_M",  0,            0,            "Dist_Out_M", "Dist_Out_M",  0,            "Dist_Out_M", 0,            0,            "Coy_Out_M",  "Coy_Out_M",  0,            "Coy_Out_M",  "ES_M",       #Out_M
    "MOCI_JFM_P",  "MOCI_AMJ_P",  "MOCI_OND_P",  0,            0,            "Dist_Out_P", "Dist_Out_P",  0,            "Dist_Out_P", 0,            0,            "Coy_Out_P",  "Coy_Out_P",  0,            "Coy_Out_P",  "ES_P"        #Out_P
  ),
  nrow = 6, ncol = 16, byrow = TRUE
)

# ----------------------------------------------------------------------------
# U MATRIX: 6 rows (hidden states) × 1 × 29 years
# ----------------------------------------------------------------------------

U1.6 <- matrix(c("t1_Est_A", "t1_Est_M", "t1_Est_P",
                 "t1_Out_A", "t1_Out_M", "t1_Out_P"), 6, 1)

U2.6 <- matrix(c("t2_Est_A", "t2_Est_M", "t1_Est_P",
                 "t2_Out_A", "t2_Out_M", "t1_Out_P"), 6, 1)

Ut.estuary <- array(U2.6, dim = c(dim(U1.6), ncol(dat)))
Ut.estuary[, , 1:ceiling(ncol(dat) / 4)] <- U1.6

# ----------------------------------------------------------------------------
# C MATRIX: 6 rows (hidden states) × 16 columns (covariates)
# ----------------------------------------------------------------------------

# Hidden state order: Est_A, Est_M, Est_P, Out_A, Out_M, Out_P
# Need to decide how covariates map to grouped states

# Option: Use site-averaged covariates OR pick representative covariates
# For simplicity, use class-specific effects only (MOCI, pooled Dist, pooled Coy, ES)

C.model.estuary <- matrix(
  list(
    #MOCI_JFM   MOCI_AMJ   MOCI_OND   Dist_BL   Dist_DE   Dist_DP   Dist_PRH   Dist_TB   Dist_TP   Coy_BL    Coy_DE    Coy_DP    Coy_PRH   Coy_TB    Coy_TP    eSeal
    "MOCI_A",  "MOCI_A",  "MOCI_A",  "Dist_A", "Dist_A", 0,        0,         "Dist_A", 0,        "Coy_A",  "Coy_A",  0,        0,        "Coy_A",  0,        "ES_A",   #Est_A (BL, DE, TB)
    "MOCI_M",  "MOCI_M",  "MOCI_M",  "Dist_M", "Dist_M", 0,        0,         "Dist_M", 0,        "Coy_M",  "Coy_M",  0,        0,        "Coy_M",  0,        "ES_M",   #Est_M
    "MOCI_P",  "MOCI_P",  "MOCI_P",  "Dist_P", "Dist_P", 0,        0,         "Dist_P", 0,        "Coy_P",  "Coy_P",  0,        0,        "Coy_P",  0,        "ES_P",   #Est_P
    
    "MOCI_A",  "MOCI_A",  "MOCI_A",  0,        0,        "Dist_A", "Dist_A",  0,        "Dist_A", 0,        0,        "Coy_A",  "Coy_A",  0,        "Coy_A",  "ES_A",   #Out_A (DP, PRH, TP)
    "MOCI_M",  "MOCI_M",  "MOCI_M",  0,        0,        "Dist_M", "Dist_M",  0,        "Dist_M", 0,        0,        "Coy_M",  "Coy_M",  0,        "Coy_M",  "ES_M",   #Out_M
    "MOCI_P",  "MOCI_P",  "MOCI_P",  0,        0,        "Dist_P", "Dist_P",  0,        "Dist_P", 0,        0,        "Coy_P",  "Coy_P",  0,        "Coy_P",  "ES_P"    #Out_P
  ),
  nrow = 6, ncol = 16, byrow = TRUE
)

# ----------------------------------------------------------------------------
# VERIFY DIMENSIONS
# ----------------------------------------------------------------------------

cat("Z.estuary levels:", nlevels(Z.estuary), "\n")
cat("Ut.estuary:", dim(Ut.estuary), "\n")
cat("C.model.estuary:", dim(C.model.estuary), "\n")

# ----------------------------------------------------------------------------
# RUN MODEL B
# ----------------------------------------------------------------------------

t0 <- Sys.time()
m.B.estuary <- MARSS(dat, model = list(
  Z = Z.estuary,
  U = Ut.estuary,
  R = diag(0.025, 18),
  Q = "diagonal and equal",
  B = "identity",
  C = C.model.estuary,
  c = cov_t_marss,
  tinitx = 1
), control = list(maxit = 5000, safe = TRUE, trace = 0, allow.degen = TRUE))
t1 <- Sys.time()
cat("Model B run time:", t1 - t0, "\n")

# ----------------------------------------------------------------------------
# COMPARE MODELS
# ----------------------------------------------------------------------------

df_aic <- data.frame(
  model = c("A: Independent (18 states)", "B: Estuarine vs Outer (6 states)"),
  AICc = c(m.A.indep$AICc, m.B.estuary$AICc),
  num_params = c(m.A.indep$num.params, m.B.estuary$num.params),
  logLik = c(m.A.indep$logLik, m.B.estuary$logLik)
)

df_aic$deltaAICc <- df_aic$AICc - min(df_aic$AICc)
df_aic <- df_aic[order(df_aic$AICc), ]

print(df_aic)

# ============================================================================
# MODEL C: ALL ONE POPULATION (3 hidden states by class only)
# ============================================================================

Z.one <- factor(c(
  "A", "M", "P",   # BL
  "A", "M", "P",   # DE
  "A", "M", "P",   # DP
  "A", "M", "P",   # PRH
  "A", "M", "P",   # TB
  "A", "M", "P"    # TP
))

# U matrix: 3 rows × 1 × 29 years
U1.3 <- matrix(c("t1_A", "t1_M", "t1_P"), 3, 1)
U2.3 <- matrix(c("t2_A", "t2_M", "t1_P"), 3, 1)

Ut.one <- array(U2.3, dim = c(dim(U1.3), ncol(dat)))
Ut.one[, , 1:ceiling(ncol(dat) / 4)] <- U1.3

# ============================================================================
# MODEL C: ALL ONE POPULATION (3 hidden states by class only)
# C MATRIX: 3 rows × 16 columns
# ============================================================================

C.model.one <- matrix(
  list(
    #MOCI_JFM       MOCI_AMJ       MOCI_OND       Dist_BL    Dist_DE    Dist_DP    Dist_PRH    Dist_TB    Dist_TP    Coy_BL     Coy_DE     Coy_DP     Coy_PRH    Coy_TB     Coy_TP     eSeal
    "MOCI_JFM_A",  "MOCI_AMJ_A",  "MOCI_OND_A",  "Dist_A",  "Dist_A",  "Dist_A",  "Dist_A",   "Dist_A",  "Dist_A",  "Coy_A",   "Coy_A",   "Coy_A",   "Coy_A",   "Coy_A",   "Coy_A",   "ES_A",   #A
    "MOCI_JFM_M",  "MOCI_AMJ_M",  "MOCI_OND_M",  "Dist_M",  "Dist_M",  "Dist_M",  "Dist_M",   "Dist_M",  "Dist_M",  "Coy_M",   "Coy_M",   "Coy_M",   "Coy_M",   "Coy_M",   "Coy_M",   "ES_M",   #M
    "MOCI_JFM_P",  "MOCI_AMJ_P",  "MOCI_OND_P",  "Dist_P",  "Dist_P",  "Dist_P",  "Dist_P",   "Dist_P",  "Dist_P",  "Coy_P",   "Coy_P",   "Coy_P",   "Coy_P",   "Coy_P",   "Coy_P",   "ES_P"    #P
  ),
  nrow = 3, ncol = 16, byrow = TRUE
)
# Verify
cat("Model C dimensions:\n")
cat("Z.one levels:", nlevels(Z.one), "\n")
cat("Ut.one:", dim(Ut.one), "\n")
cat("C.model.one:", dim(C.model.one), "\n")

# Run Model C
t0 <- Sys.time()
m.C.one <- MARSS(dat, model = list(
  Z = Z.one,
  U = Ut.one,
  R = diag(0.025, 18),
  Q = "diagonal and equal",
  B = "identity",
  C = C.model.one,
  c = cov_t_marss,
  tinitx = 1
), control = list(maxit = 5000, safe = TRUE, trace = 0, allow.degen = TRUE))
t1 <- Sys.time()
cat("Model C run time:", t1 - t0, "\n")


# ============================================================================
# MODEL D: MOLT SEPARATE FROM ADULTS/PUPS (2 populations × class structure)
# ============================================================================

# Hypothesis: Molting seals behave differently (different site fidelity, timing)
# Group 1: Adults + Pups (breeding population)
# Group 2: Molting animals

Z.molt <- factor(c(
  "Breed_A", "Molt", "Breed_P",   # BL
  "Breed_A", "Molt", "Breed_P",   # DE
  "Breed_A", "Molt", "Breed_P",   # DP
  "Breed_A", "Molt", "Breed_P",   # PRH
  "Breed_A", "Molt", "Breed_P",   # TB
  "Breed_A", "Molt", "Breed_P"    # TP
))

# U matrix: 3 rows (Breed_A, Molt, Breed_P) × 1 × 29 years
U1.molt <- matrix(c("t1_Breed_A", "t1_Molt", "t1_Breed_P"), 3, 1)
U2.molt <- matrix(c("t2_Breed_A", "t2_Molt", "t1_Breed_P"), 3, 1)  # Pups constant

Ut.molt <- array(U2.molt, dim = c(dim(U1.molt), ncol(dat)))
Ut.molt[, , 1:ceiling(ncol(dat) / 4)] <- U1.molt

# C matrix: 3 rows × 16 columns
C.model.molt <- matrix(
  list(
    #MOCI_JFM          MOCI_AMJ          MOCI_OND          Dist_BL        Dist_DE        Dist_DP        Dist_PRH        Dist_TB        Dist_TP        Coy_BL         Coy_DE         Coy_DP         Coy_PRH        Coy_TB         Coy_TP         eSeal
    "MOCI_JFM_Breed",  "MOCI_AMJ_Breed", "MOCI_OND_Breed", "Dist_Breed",  "Dist_Breed",  "Dist_Breed",  "Dist_Breed",   "Dist_Breed",  "Dist_Breed",  "Coy_Breed",   "Coy_Breed",   "Coy_Breed",   "Coy_Breed",   "Coy_Breed",   "Coy_Breed",   "ES_Breed",   #Breed_A
    "MOCI_JFM_Molt",   "MOCI_AMJ_Molt",  "MOCI_OND_Molt",  "Dist_Molt",   "Dist_Molt",   "Dist_Molt",   "Dist_Molt",    "Dist_Molt",   "Dist_Molt",   "Coy_Molt",    "Coy_Molt",    "Coy_Molt",    "Coy_Molt",    "Coy_Molt",    "Coy_Molt",    "ES_Molt",    #Molt
    "MOCI_JFM_Breed",  "MOCI_AMJ_Breed", "MOCI_OND_Breed", "Dist_Breed",  "Dist_Breed",  "Dist_Breed",  "Dist_Breed",   "Dist_Breed",  "Dist_Breed",  "Coy_Breed",   "Coy_Breed",   "Coy_Breed",   "Coy_Breed",   "Coy_Breed",   "Coy_Breed",   "ES_Breed"    #Breed_P
  ),
  nrow = 3, ncol = 16, byrow = TRUE
)

# Verify
cat("\nModel D dimensions:\n")
cat("Z.molt levels:", nlevels(Z.molt), "\n")
cat("Ut.molt:", dim(Ut.molt), "\n")
cat("C.model.molt:", dim(C.model.molt), "\n")

# Run Model D
t0 <- Sys.time()
m.D.molt <- MARSS(dat, model = list(
  Z = Z.molt,
  U = Ut.molt,
  R = diag(0.025, 18),
  Q = "diagonal and equal",
  B = "identity",
  C = C.model.molt,
  c = cov_t_marss,
  tinitx = 1
), control = list(maxit = 5000, safe = TRUE, trace = 0, allow.degen = TRUE))
t1 <- Sys.time()
cat("Model D run time:", t1 - t0, "\n")


# ============================================================================
# COMPARE ALL MODELS
# ============================================================================

df_aic <- data.frame(
  model = c(
    "A: Independent (18 states)",
    "B: Estuarine vs Outer (6 states)",
    "C: One population (3 states)",
    "D: Molt separate (3 states)"
  ),
  AICc = c(m.A.indep$AICc, m.B.estuary$AICc, m.C.one$AICc, m.D.molt$AICc),
  num_params = c(m.A.indep$num.params, m.B.estuary$num.params, m.C.one$num.params, m.D.molt$num.params),
  logLik = c(m.A.indep$logLik, m.B.estuary$logLik, m.C.one$logLik, m.D.molt$logLik)
)

df_aic$deltaAICc <- df_aic$AICc - min(df_aic$AICc)
df_aic <- df_aic[order(df_aic$AICc), ]

print(df_aic)

# Save all models
save(m.C.one, file = "Output/m.C.one.RData")
save(m.D.molt, file = "Output/m.D.molt.RData")



