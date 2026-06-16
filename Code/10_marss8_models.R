# ============================================================================
# 10_marss8_models.R
# ----------------------------------------------------------------------------
# 8-SITE MARSS — SUPPLEMENTAL ANALYSIS (Appendix A)
# Characterizes trend + covariate dynamics at all 8 sites, incl. haul-out-only
# DR & PB, relative to the 6 IPM breeding sites. NOT used to seed the IPM.
#
# Prereqs (run in order):
#   source("Code/08_data_prep_8site.R")   # -> dat_8site (24 x 29), years_8site
#   source("Code/09_covariates_8site.R")  # -> cov_t_scaled_8site (20 x 29)
#
# State order (24, canonical site-major A,M,P): BL_*, DE_*, DP_*, DR_*, PB_*,
#   PRH_*, TB_*, TP_*   (DR_PUP, PB_PUP are all-NA → uninformative).
# Covariate order (20): MOCI x3, Dist x8, Coyote x8, eSeal.
#
# Models:
#   A  24 independent states (expected best)
#   B  Breeding vs Haul-out grouping — Z has 6 levels but only 5 are
#      IDENTIFIABLE: HaulOut_P maps only to DR_PUP & PB_PUP, both all-NA, so
#      that latent state carries no data. (Earlier comments said "5 states";
#      it is 6 levels / 5 identifiable.)
#
# Produces: m.A_8site, m.B_8site, CIs_8site -> Output/*.RData
# Next: source("Code/11_marss8_plots.R")
# ============================================================================

library(MARSS)

if (!exists("dat_8site"))          stop("Run 08_data_prep_8site.R first.")
if (!exists("cov_t_scaled_8site")) stop("Run 09_covariates_8site.R first.")

SITES_8 <- c("BL","DE","DP","DR","PB","PRH","TB","TP")
CANON_STATE_8 <- as.vector(t(outer(SITES_8, c("ADULT","MOLTING","PUP"),
                                   paste, sep = "_")))
CANON_COV_8 <- c("MOCI_JFM","MOCI_AMJ","MOCI_OND",
                 paste0("Dist_",   SITES_8),
                 paste0("Coyote_", SITES_8),
                 "eSeal_Sum_Imm_MaxCount")

stopifnot(
  nrow(dat_8site) == 24L, nrow(cov_t_scaled_8site) == 20L,
  ncol(dat_8site) == ncol(cov_t_scaled_8site),
  identical(rownames(dat_8site), CANON_STATE_8),
  identical(rownames(cov_t_scaled_8site), CANON_COV_8)
)

cov <- cov_t_scaled_8site
TT  <- ncol(dat_8site)
cat(sprintf("Data %d x %d | Covariates %d x %d\n",
            nrow(dat_8site), TT, nrow(cov), ncol(cov)))

# ── Time-varying U: adult & molt shift at 2004; pups constant ───────────────
U1 <- matrix(rep(c("t1_A","t1_M","t1_P"), 8), 24, 1)
U2 <- matrix(rep(c("t2_A","t2_M","t1_P"), 8), 24, 1)
Ut_8site <- array(U2, dim = c(dim(U1), TT))
breakpoint <- ceiling(TT / 4)
Ut_8site[, , 1:breakpoint] <- U1
cat("U breakpoint at column", breakpoint, "= year", years_8site[breakpoint], "\n")

# ── C matrix (24 states × 20 covariates), column order = CANON_COV_8 ────────
# Class-specific MOCI; site-specific Dist & Coy; eSeal at DE & PRH.
# DR/PB: MOCI + own Dist only (no coyote, no eSeal); their PUP rows = 0.
C.8site <- matrix(list(
  #            JFM          AMJ          OND         BL        DE        DP        DR        PB        PRH        TB        TP        cBL      cDE      cDP      cDR cPB cPRH cTB cTP eSeal
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A","Dist_BL", 0,        0,        0,        0,        0,         0,        0,        "Coy_BL",0,       0,       0,  0,  0,   0,  0,  0,        #BL_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M","Dist_BL", 0,        0,        0,        0,        0,         0,        0,        "Coy_BL",0,       0,       0,  0,  0,   0,  0,  0,        #BL_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P","Dist_BL", 0,        0,        0,        0,        0,         0,        0,        "Coy_BL",0,       0,       0,  0,  0,   0,  0,  0,        #BL_P
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A",0,        "Dist_DE", 0,        0,        0,        0,         0,        0,        0,       "Coy_DE",0,       0,  0,  0,   0,  0,  "ES-DE",  #DE_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M",0,        "Dist_DE", 0,        0,        0,        0,         0,        0,        0,       "Coy_DE",0,       0,  0,  0,   0,  0,  "ES-DE",  #DE_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P",0,        "Dist_DE", 0,        0,        0,        0,         0,        0,        0,       "Coy_DE",0,       0,  0,  0,   0,  0,  "ES-DE",  #DE_P
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A",0,        0,        "Dist_DP", 0,        0,        0,         0,        0,        0,       0,       "Coy_DP",0,  0,  0,   0,  0,  0,        #DP_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M",0,        0,        "Dist_DP", 0,        0,        0,         0,        0,        0,       0,       "Coy_DP",0,  0,  0,   0,  0,  0,        #DP_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P",0,        0,        "Dist_DP", 0,        0,        0,         0,        0,        0,       0,       "Coy_DP",0,  0,  0,   0,  0,  0,        #DP_P
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A",0,        0,        0,        "Dist_DR", 0,        0,         0,        0,        0,       0,       0,       0,  0,  0,   0,  0,  0,        #DR_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M",0,        0,        0,        "Dist_DR", 0,        0,         0,        0,        0,       0,       0,       0,  0,  0,   0,  0,  0,        #DR_M
  0,           0,           0,           0,        0,        0,        0,        0,        0,         0,        0,        0,       0,       0,       0,  0,  0,   0,  0,  0,        #DR_P (all-NA)
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A",0,        0,        0,        0,        "Dist_PB", 0,         0,        0,        0,       0,       0,       0,  0,  0,   0,  0,  0,        #PB_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M",0,        0,        0,        0,        "Dist_PB", 0,         0,        0,        0,       0,       0,       0,  0,  0,   0,  0,  0,        #PB_M
  0,           0,           0,           0,        0,        0,        0,        0,        0,         0,        0,        0,       0,       0,       0,  0,  0,   0,  0,  0,        #PB_P (all-NA)
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A",0,        0,        0,        0,        0,        "Dist_PRH", 0,        0,        0,       0,       0,       0,  0,  0,   0,  0,  "ES-PRH", #PRH_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M",0,        0,        0,        0,        0,        "Dist_PRH", 0,        0,        0,       0,       0,       0,  0,  0,   0,  0,  "ES-PRH", #PRH_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P",0,        0,        0,        0,        0,        "Dist_PRH", 0,        0,        0,       0,       0,       0,  0,  0,   0,  0,  "ES-PRH", #PRH_P
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A",0,        0,        0,        0,        0,        0,         "Dist_TB", 0,        0,       0,       0,       0,  0,  0,   0,  0,  0,        #TB_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M",0,        0,        0,        0,        0,        0,         "Dist_TB", 0,        0,       0,       0,       0,  0,  0,   0,  0,  0,        #TB_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P",0,        0,        0,        0,        0,        0,         "Dist_TB", 0,        0,       0,       0,       0,  0,  0,   0,  0,  0,        #TB_P
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A",0,        0,        0,        0,        0,        0,         0,        "Dist_TP", 0,       0,       0,       0,  0,  0,   0,  0,  0,        #TP_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M",0,        0,        0,        0,        0,        0,         0,        "Dist_TP", 0,       0,       0,       0,  0,  0,   0,  0,  0,        #TP_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P",0,        0,        0,        0,        0,        0,         0,        "Dist_TP", 0,       0,       0,       0,  0,  0,   0,  0,  0         #TP_P
), nrow = 24, ncol = 20, byrow = TRUE)

stopifnot(dim(C.8site) == c(24L, 20L))

R_fixed <- diag(0.025, 24)
ctrl    <- list(maxit = 5000, safe = TRUE, trace = 0, allow.degen = TRUE)

run_marss8 <- function(label, Z, U, C, cmat) {
  cat(sprintf("\n── %s ─────────────────────────────────\n", label))
  t0 <- Sys.time()
  m <- MARSS(dat_8site, model = list(Z = Z, U = U, R = R_fixed,
                                     Q = "diagonal and equal", B = "identity",
                                     C = C, c = cmat, tinitx = 1),
             control = ctrl)
  cat(sprintf("  %s | %.1fs | AICc=%.2f | params=%d\n",
              label, as.numeric(difftime(Sys.time(), t0, units="secs")),
              m$AICc, m$num.params))
  m
}

# ── Model A: 24 independent states ──────────────────────────────────────────
m.A_8site <- run_marss8("A: 24 independent states",
                        factor(1:24), Ut_8site, C.8site, cov)
save(m.A_8site, file = "Output/m.A_8site.RData")

# ── Model B: Breeding vs Haul-out (6 levels, 5 identifiable) ─────────────────
Z.breed_haulout <- factor(c(
  "Breed_A","Breed_M","Breed_P",      # BL
  "Breed_A","Breed_M","Breed_P",      # DE
  "Breed_A","Breed_M","Breed_P",      # DP
  "HaulOut_A","HaulOut_M","HaulOut_P",# DR  (HaulOut_P unidentified: DR/PB pups all-NA)
  "HaulOut_A","HaulOut_M","HaulOut_P",# PB
  "Breed_A","Breed_M","Breed_P",      # PRH
  "Breed_A","Breed_M","Breed_P",      # TB
  "Breed_A","Breed_M","Breed_P"))     # TP

U1.bh <- matrix(c("t1_BreedA","t1_BreedM","t1_BreedP",
                  "t1_HaulA", "t1_HaulM", "t1_HaulP"), 6, 1)
U2.bh <- matrix(c("t2_BreedA","t2_BreedM","t1_BreedP",
                  "t2_HaulA", "t2_HaulM", "t1_HaulP"), 6, 1)
Ut.bh <- array(U2.bh, dim = c(dim(U1.bh), TT))
Ut.bh[, , 1:breakpoint] <- U1.bh

# Grouped covariates: class MOCI + pooled Breed/HaulOut Dist + Breed Coy + eSeal
C.bh <- matrix(list(
  #         JFM       AMJ       OND       Dist_Breed   Dist_HaulOut  Coy_Breed   eSeal
  "MOCI_A","MOCI_A","MOCI_A", "Dist_Breed", 0,           "Coy_Breed", "ES_A",  #Breed_A
  "MOCI_M","MOCI_M","MOCI_M", "Dist_Breed", 0,           "Coy_Breed", "ES_M",  #Breed_M
  "MOCI_P","MOCI_P","MOCI_P", "Dist_Breed", 0,           "Coy_Breed", "ES_P",  #Breed_P
  "MOCI_A","MOCI_A","MOCI_A", 0,            "Dist_Haulout",0,          0,       #HaulOut_A
  "MOCI_M","MOCI_M","MOCI_M", 0,            "Dist_Haulout",0,          0,       #HaulOut_M
  "MOCI_P","MOCI_P","MOCI_P", 0,            0,            0,          0         #HaulOut_P (unidentified)
), nrow = 6, ncol = 7, byrow = TRUE)

cov_bh <- rbind(
  cov["MOCI_JFM", ],
  cov["MOCI_AMJ", ],
  cov["MOCI_OND", ],
  colMeans(cov[c("Dist_BL","Dist_DE","Dist_DP","Dist_PRH","Dist_TB","Dist_TP"), ]),
  colMeans(cov[c("Dist_DR","Dist_PB"), ]),
  colMeans(cov[c("Coyote_BL","Coyote_DE","Coyote_DP"), ]),
  cov["eSeal_Sum_Imm_MaxCount", ]
)
rownames(cov_bh) <- c("MOCI_A","MOCI_AMJ","MOCI_OND",
                      "Dist_Breed","Dist_Haulout","Coy_Breed","eSeal")

m.B_8site <- run_marss8("B: Breeding vs Haul-out (6 levels / 5 identifiable)",
                        Z.breed_haulout, Ut.bh, C.bh, cov_bh)
save(m.B_8site, file = "Output/m.B_8site.RData")

# ── Comparison ──────────────────────────────────────────────────────────────
# NOTE: B also collapses 20 site-specific covariates to 7 grouped ones, so the
# A-vs-B AICc gap reflects BOTH state structure AND covariate detail. The
# cleaner "sites independent" test is the 6-site model selection (03), where
# all candidates share one covariate structure. (Manuscript P25.)
df_aic_8site <- data.frame(
  model      = c("A: 24 independent states",
                 "B: Breeding(18) vs Haul-out(6) [5 identifiable]"),
  AICc       = c(m.A_8site$AICc, m.B_8site$AICc),
  num_params = c(m.A_8site$num.params, m.B_8site$num.params),
  logLik     = c(m.A_8site$logLik, m.B_8site$logLik)
)
df_aic_8site$deltaAICc <- df_aic_8site$AICc - min(df_aic_8site$AICc)
df_aic_8site <- df_aic_8site[order(df_aic_8site$AICc), ]
cat("\n── 8-site Model Comparison ─────────────────────────────────────────\n")
print(df_aic_8site)
write.csv(df_aic_8site, "Output/model_comparison_8site_AIC.csv", row.names = FALSE)

# ── 89% CIs for the best (independent) model ────────────────────────────────
cat("\n── Computing 89% CIs for m.A_8site (alpha = 0.11) ──────────────────\n")
CIs_8site <- MARSSparamCIs(m.A_8site, alpha = 0.11)
save(CIs_8site, file = "Output/CIs_8site.RData")

cat("\nObjects created: m.A_8site, m.B_8site, CIs_8site, df_aic_8site\n")
cat("Next: source(\"Code/11_marss8_plots.R\")\n")