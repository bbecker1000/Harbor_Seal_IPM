# ============================================================================
# MARSS MODELS — 8-SITE ANALYSIS (INCLUDING DR AND PB HAUL-OUT SITES)
# 
# Purpose: Supplemental analysis to characterize trend and covariate dynamics
# at non-breeding haul-out sites relative to the 6 IPM breeding sites.
# These MARSS results are NOT used as starting values for the IPM.
#
# Source files first:
#   source("4a_data_prep_8sites_MARSS.R")   # produces dat_8site, years_8site
#   source("4b_covariates_8sites_MARSS.R")  # produces cov_t_scaled_8site
#
# State ordering (24 rows, alphabetical by site then class):
#   BL_A BL_M BL_P | DE_A DE_M DE_P | DP_A DP_M DP_P |
#   DR_A DR_M DR_P | PB_A PB_M PB_P | PRH_A PRH_M PRH_P |
#   TB_A TB_M TB_P | TP_A TP_M TP_P
#
# DR_P and PB_P rows in dat_8site are NA throughout (haul-out-only sites).
# MARSS correctly treats NA as missing data.
#
# Covariate order (20 columns in C):
#   1:MOCI_JFM  2:MOCI_AMJ  3:MOCI_OND
#   4:Dist_BL  5:Dist_DE  6:Dist_DP  7:Dist_DR  8:Dist_PB
#   9:Dist_PRH  10:Dist_TB  11:Dist_TP
#   12:Coy_BL  13:Coy_DE  14:Coy_DP
#   15:Coy_DR(0)  16:Coy_PB(0)  17:Coy_PRH(0)  18:Coy_TB(0)  19:Coy_TP(0)
#   20:eSeal
# ============================================================================

library(MARSS)

# ── Verify source data are loaded ─────────────────────────────────────────────
if (!exists("dat_8site"))          stop("Run 4a_data_prep_8sites_MARSS.R first")
if (!exists("cov_t_scaled_8site")) stop("Run 4b_covariates_8sites_MARSS.R first")

cov <- cov_t_scaled_8site   # shorthand
cat("Data:      ", nrow(dat_8site), "×", ncol(dat_8site), "\n")
cat("Covariates:", nrow(cov), "×", ncol(cov), "\n")

# ── Time-varying U matrix (24 states × 1 × 29 years) ─────────────────────────
# Adults and Molt: growth rate shifts at ~2004 (column 8, year 2004)
# Pups: no temporal shift (constant U across time series)
# DR and PB: same U structure as other sites (shared metapopulation dynamics)

U1 <- matrix(c(
  "t1_A", "t1_M", "t1_P",   # BL
  "t1_A", "t1_M", "t1_P",   # DE
  "t1_A", "t1_M", "t1_P",   # DP
  "t1_A", "t1_M", "t1_P",   # DR ← haul-out site, same U class
  "t1_A", "t1_M", "t1_P",   # PB ← haul-out site, same U class
  "t1_A", "t1_M", "t1_P",   # PRH
  "t1_A", "t1_M", "t1_P",   # TB
  "t1_A", "t1_M", "t1_P"    # TP
), 24, 1)

U2 <- matrix(c(
  "t2_A", "t2_M", "t1_P",   # BL
  "t2_A", "t2_M", "t1_P",   # DE
  "t2_A", "t2_M", "t1_P",   # DP
  "t2_A", "t2_M", "t1_P",   # DR
  "t2_A", "t2_M", "t1_P",   # PB
  "t2_A", "t2_M", "t1_P",   # PRH
  "t2_A", "t2_M", "t1_P",   # TB
  "t2_A", "t2_M", "t1_P"    # TP
), 24, 1)

TT        <- ncol(dat_8site)
Ut_8site  <- array(U2, dim = c(dim(U1), TT))
breakpoint <- ceiling(TT / 4)
Ut_8site[, , 1:breakpoint] <- U1

cat("U breakpoint at column", breakpoint, "= year", years_8site[breakpoint], "\n")

# ── C matrix (24 rows × 20 columns) ──────────────────────────────────────────
# Column annotation (see header):
#  JFM AMJ OND  | D_BL D_DE D_DP D_DR D_PB D_PRH D_TB D_TP |
#  C_BL C_DE C_DP C_DR C_PB C_PRH C_TB C_TP | eSeal

C.8site <- matrix(list(
  # ── MOCI_JFM  MOCI_AMJ  MOCI_OND | D_BL   D_DE   D_DP   D_DR   D_PB   D_PRH  D_TB   D_TP | C_BL    C_DE    C_DP    C_DR  C_PB  C_PRH  C_TB   C_TP  | eSeal
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A", "Dist_BL",0,      0,     0,     0,    0,    0,    0,    "Coy_BL",0,      0,      0,    0,    0,    0,    0,    0,       #BL_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M", "Dist_BL",0,      0,     0,     0,    0,    0,    0,    "Coy_BL",0,      0,      0,    0,    0,    0,    0,    0,       #BL_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P", "Dist_BL",0,      0,     0,     0,    0,    0,    0,    "Coy_BL",0,      0,      0,    0,    0,    0,    0,    0,       #BL_P
  
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A", 0,"Dist_DE",      0,     0,     0,    0,    0,    0,    0,"Coy_DE",       0,      0,    0,    0,    0,    0,    "ES-DE", #DE_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M", 0,"Dist_DE",      0,     0,     0,    0,    0,    0,    0,"Coy_DE",       0,      0,    0,    0,    0,    0,    "ES-DE", #DE_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P", 0,"Dist_DE",      0,     0,     0,    0,    0,    0,    0,"Coy_DE",       0,      0,    0,    0,    0,    0,    "ES-DE", #DE_P
  
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A", 0,0,      "Dist_DP",     0,     0,    0,    0,    0,    0,0,      "Coy_DP",       0,    0,    0,    0,    0,    0,       #DP_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M", 0,0,      "Dist_DP",     0,     0,    0,    0,    0,    0,0,      "Coy_DP",       0,    0,    0,    0,    0,    0,       #DP_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P", 0,0,      "Dist_DP",     0,     0,    0,    0,    0,    0,0,      "Coy_DP",       0,    0,    0,    0,    0,    0,       #DP_P
  
  # DR: haul-out only — MOCI + Dist_DR; no coyote; no eSeal
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A", 0,0,      0,     "Dist_DR",0,    0,    0,    0,    0,0,      0,      0,    0,    0,    0,    0,    0,                   #DR_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M", 0,0,      0,     "Dist_DR",0,    0,    0,    0,    0,0,      0,      0,    0,    0,    0,    0,    0,                   #DR_M
  0,           0,           0,            0,0,      0,     0,        0,    0,    0,    0,    0,0,      0,      0,    0,    0,    0,    0,    0,                   #DR_P (NA row, no covariates)
  
  # PB: haul-out only — MOCI + Dist_PB; no coyote; no eSeal
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A", 0,0,      0,     0,        "Dist_PB",0,0,  0,    0,0,      0,      0,    0,    0,    0,    0,    0,                   #PB_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M", 0,0,      0,     0,        "Dist_PB",0,0,  0,    0,0,      0,      0,    0,    0,    0,    0,    0,                   #PB_M
  0,           0,           0,            0,0,      0,     0,        0,        0,0,  0,    0,0,      0,      0,    0,    0,    0,    0,    0,                   #PB_P (NA row)
  
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A", 0,0,      0,     0,     0,    "Dist_PRH",0, 0,    0,0,      0,      0,    0,    0,    0,    0,    "ES-PRH",            #PRH_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M", 0,0,      0,     0,     0,    "Dist_PRH",0, 0,    0,0,      0,      0,    0,    0,    0,    0,    "ES-PRH",            #PRH_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P", 0,0,      0,     0,     0,    "Dist_PRH",0, 0,    0,0,      0,      0,    0,    0,    0,    0,    "ES-PRH",            #PRH_P
  
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A", 0,0,      0,     0,     0,    0,    "Dist_TB",0,  0,0,      0,      0,    0,    0,    0,    0,    0,                   #TB_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M", 0,0,      0,     0,     0,    0,    "Dist_TB",0,  0,0,      0,      0,    0,    0,    0,    0,    0,                   #TB_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P", 0,0,      0,     0,     0,    0,    "Dist_TB",0,  0,0,      0,      0,    0,    0,    0,    0,    0,                   #TB_P
  
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A", 0,0,      0,     0,     0,    0,    0,    "Dist_TP",0,0,      0,      0,    0,    0,    0,    0,    0,                 #TP_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M", 0,0,      0,     0,     0,    0,    0,    "Dist_TP",0,0,      0,      0,    0,    0,    0,    0,    0,                 #TP_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P", 0,0,      0,     0,     0,    0,    0,    "Dist_TP",0,0,      0,      0,    0,    0,    0,    0,    0                  #TP_P
), nrow = 24, ncol = 20, byrow = TRUE)

cat("C matrix dimensions: ", nrow(C.8site), "×", ncol(C.8site),
    " (should be 24×20)\n")

# ── Model A: All 24 independent states ───────────────────────────────────────
cat("\n── Running Model A: 24 independent states ──────────────────────────────\n")
t0 <- Sys.time()
m.A_8site <- MARSS(dat_8site, model = list(
  Z = factor(1:24),
  U = Ut_8site,
  R = diag(0.025, 24),
  Q = "diagonal and equal",
  B = "identity",
  C = C.8site,
  c = cov,
  tinitx = 1
), control = list(maxit = 5000, safe = TRUE, trace = 0, allow.degen = TRUE))
cat("Model A run time:", difftime(Sys.time(), t0, units = "secs"), "sec\n")
cat("AICc:", m.A_8site$AICc, "\n")

save(m.A_8site, file = "Output/m.A_8site.RData")

# ── Model B: Breeding sites vs haul-out sites (ecological hypothesis) ─────────
# This directly tests whether DR and PB dynamics differ from breeding sites
# Breeding sites: BL, DE, DP, PRH, TB, TP (6 sites × 3 classes = 18 states)
# Haul-out sites: DR, PB (2 sites × 2 classes observed = 4 states, pup NA)
# 8 groups: Breed_A, Breed_M, Breed_P, HaulOut_A, HaulOut_M,
#           + site-specific offsets via a matrix (A parameter)

# For this comparison, restrict to Adult and Molt only at DR and PB
# (i.e., the 22-observation version without NA pup rows for DR/PB)
# Grouped Z: 5 hidden states (Breed_A, Breed_M, Breed_P, HaulOut_A, HaulOut_M)
Z.breed_haulout <- factor(c(
  "Breed_A", "Breed_M", "Breed_P",   # BL
  "Breed_A", "Breed_M", "Breed_P",   # DE
  "Breed_A", "Breed_M", "Breed_P",   # DP
  "HaulOut_A","HaulOut_M","HaulOut_P",# DR  (HaulOut_P will be NA → unidentified)
  "HaulOut_A","HaulOut_M","HaulOut_P",# PB  (HaulOut_P will be NA → unidentified)
  "Breed_A", "Breed_M", "Breed_P",   # PRH
  "Breed_A", "Breed_M", "Breed_P",   # TB
  "Breed_A", "Breed_M", "Breed_P"    # TP
))

U1.bh <- matrix(c("t1_BreedA", "t1_BreedM", "t1_BreedP",
                  "t1_HaulA",  "t1_HaulM",  "t1_HaulP"),  6, 1)
U2.bh <- matrix(c("t2_BreedA", "t2_BreedM", "t1_BreedP",
                  "t2_HaulA",  "t2_HaulM",  "t1_HaulP"),  6, 1)
Ut.bh <- array(U2.bh, dim = c(dim(U1.bh), TT))
Ut.bh[, , 1:breakpoint] <- U1.bh

# Simplified C for grouped model: class-specific MOCI + grouped effects
# Breeding and haul-out sites get same MOCI but different Dist/Coy effects
C.bh <- matrix(list(
  # JFM      AMJ      OND      Dist_Breed   Dist_HaulOut   Coy_Breed    eSeal
  "MOCI_A","MOCI_A","MOCI_A", "Dist_Breed", 0,            "Coy_Breed", "ES_A",  #Breed_A
  "MOCI_M","MOCI_M","MOCI_M", "Dist_Breed", 0,            "Coy_Breed", "ES_M",  #Breed_M
  "MOCI_P","MOCI_P","MOCI_P", "Dist_Breed", 0,            "Coy_Breed", "ES_P",  #Breed_P
  "MOCI_A","MOCI_A","MOCI_A", 0,            "Dist_Haulout",0,           0,       #HaulOut_A
  "MOCI_M","MOCI_M","MOCI_M", 0,            "Dist_Haulout",0,           0,       #HaulOut_M
  "MOCI_P","MOCI_P","MOCI_P", 0,            0,             0,           0        #HaulOut_P (unobserved)
), nrow = 6, ncol = 7, byrow = TRUE)

# For this grouped C matrix, need grouped covariates
# Use mean breeding-site disturbance and mean haul-out disturbance
# Simple approach: use the first 7 covariates subset from cov
# Build a 7-row grouped covariate matrix:
cov_bh <- rbind(
  cov["MOCI_JFM", ],   # 1: MOCI_JFM
  cov["MOCI_AMJ", ],   # 2: MOCI_AMJ
  cov["MOCI_OND", ],   # 3: MOCI_OND
  colMeans(cov[c("Dist_BL","Dist_DE","Dist_DP","Dist_PRH","Dist_TB","Dist_TP"), ]),  # 4: mean Breed Dist
  colMeans(cov[c("Dist_DR","Dist_PB"), ]),                                           # 5: mean HaulOut Dist
  colMeans(cov[c("Coyote_BL","Coyote_DE","Coyote_DP"), ]),                          # 6: mean Breed Coy
  cov["eSeal_Sum_Imm_MaxCount", ]                                                    # 7: eSeal
)
rownames(cov_bh) <- c("MOCI_A","MOCI_AMJ","MOCI_OND",
                      "Dist_Breed","Dist_Haulout","Coy_Breed","eSeal")

# Note: This is a simplified covariate structure for the grouped model
# For a more rigorous comparison, use site-specific covariates as in Model A

cat("\n── Running Model B: Breeding vs Haul-out (6 hidden states) ─────────────\n")
t0 <- Sys.time()
m.B_8site <- MARSS(dat_8site, model = list(
  Z = Z.breed_haulout,
  U = Ut.bh,
  R = diag(0.025, 24),
  Q = "diagonal and equal",
  B = "identity",
  C = C.bh,
  c = cov_bh,
  tinitx = 1
), control = list(maxit = 5000, safe = TRUE, trace = 0, allow.degen = TRUE))
cat("Model B run time:", difftime(Sys.time(), t0, units = "secs"), "sec\n")
cat("AICc:", m.B_8site$AICc, "\n")

save(m.B_8site, file = "Output/m.B_8site.RData")

# ── Model comparison ──────────────────────────────────────────────────────────
df_aic_8site <- data.frame(
  model = c(
    "A: All 24 independent states",
    "B: Breeding (18) vs Haul-out (6) [5 groups]"
  ),
  AICc       = c(m.A_8site$AICc, m.B_8site$AICc),
  num_params = c(m.A_8site$num.params, m.B_8site$num.params),
  logLik     = c(m.A_8site$logLik, m.B_8site$logLik)
)
df_aic_8site$deltaAICc <- df_aic_8site$AICc - min(df_aic_8site$AICc)
df_aic_8site <- df_aic_8site[order(df_aic_8site$AICc), ]

cat("\n── Model Comparison ──────────────────────────────────────────────────────\n")
print(df_aic_8site)

# ── CI estimation for best model ─────────────────────────────────────────────
cat("\n── Computing 89% CIs for best model ────────────────────────────────────\n")
CIs_8site <- MARSSparamCIs(m.A_8site, alpha = 0.11)  # alpha=0.11 → 89% CI
save(CIs_8site, file = "Output/CIs_8site.RData")
print(CIs_8site)

write.csv(df_aic_8site, "Output/model_comparison_8site_AIC.csv", row.names = FALSE)

