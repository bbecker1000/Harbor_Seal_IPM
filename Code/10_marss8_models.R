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
# STATE ORDER (24, canonical site-major A,M,P):
#   BL_*, DE_*, DP_*, DR_*, PB_*, PRH_*, TB_*, TP_*
#   (DR_PUP, PB_PUP are all-NA throughout → uninformative for pup state)
# COVARIATE ORDER (20):
#   MOCI×3, Dist×8 (all sites incl. DR/PB), Coyote×8 (0 for DR/PB/PRH/TB/TP),
#   eSeal (shared regional index; modelled at DE and PRH only)
#
# MODELS compared via AICc:
#   A  24 independent states — expected best; justifies IPM site independence
#   B  Breeding (18) vs Haul-out (6) — ecotype grouping based on pup presence
#   C  Estuary/Bay (BL,DE,TB) vs Coastal/Exposed (DP,DR,PB,PRH,TP) — NEW
#      Tests whether sheltered-lagoon vs exposed-headland driver regimes create
#      distinct dynamics beyond the haul-out/breeding distinction.
#      Ecologically: Estuary sites are tidally sheltered; Coastal sites face
#      direct California Current upwelling variability. DR and PB are rocky
#      coastal haul-outs and slot naturally into the Coastal group.
#
# Produces: m.A_8site, m.B_8site, m.C_8site, CIs_8site -> Output/*.RData
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

# ── Helpers ─────────────────────────────────────────────────────────────────
run_marss8 <- function(label, Z, U, C, cmat) {
  cat(sprintf("\n── %s ─────────────────────────────────\n", label))
  t0 <- Sys.time()
  m <- MARSS(dat_8site,
             model = list(Z = Z, U = U, R = diag(0.025, nrow(dat_8site)),
                          Q = "diagonal and equal", B = "identity",
                          C = C, c = cmat, tinitx = 1),
             control = list(maxit = 5000, safe = TRUE, trace = 0,
                            allow.degen = TRUE))
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  %s | %.1fs | AICc=%.2f | params=%d | conv=%d\n",
              label, dt, m$AICc, m$num.params, m$convergence))
  m
}

# ── Time-varying U: adult & molt shift at 2004; pups constant ───────────────
# All three models share this structure (pup trend held constant per
# 8-site analysis finding of no significant pup shift in that window;
# legacy analysis estimated pup Phase 2 separately given the longer record).
U1 <- matrix(rep(c("t1_A","t1_M","t1_P"), 8), 24, 1)
U2 <- matrix(rep(c("t2_A","t2_M","t1_P"), 8), 24, 1)
breakpoint <- which(years_8site == 2004)
if (length(breakpoint) == 0) stop("2004 not found in years_8site")
Ut_8site <- array(U2, dim = c(24, 1, TT))
Ut_8site[, , 1:breakpoint] <- U1
cat("U breakpoint at column", breakpoint, "= year", years_8site[breakpoint], "\n")

# ── C matrix (24 states × 20 covariates): full site-specific structure ───────
# Column order exactly matches CANON_COV_8 (rows of cov_t_scaled_8site).
# DR/PB rows: MOCI and Dist only (no coyote [zero], no eSeal [absent]).
# DR_PUP / PB_PUP rows: all zero (states are all-NA; C has no effect).
C.8site <- matrix(list(
  #            JFM           AMJ           OND          BL        DE        DP        DR        PB        PRH        TB        TP        cBL      cDE      cDP  cDR cPB cPRH cTB cTP eSeal
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A","Dist_BL", 0,        0,        0,        0,        0,         0,        0,        "Coy_BL",0,       0,       0,  0,  0,   0,  0,  0,         #BL_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M","Dist_BL", 0,        0,        0,        0,        0,         0,        0,        "Coy_BL",0,       0,       0,  0,  0,   0,  0,  0,         #BL_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P","Dist_BL", 0,        0,        0,        0,        0,         0,        0,        "Coy_BL",0,       0,       0,  0,  0,   0,  0,  0,         #BL_P
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A",0,        "Dist_DE", 0,        0,        0,        0,         0,        0,        0,       "Coy_DE",0,       0,  0,  0,   0,  0,  "ES-DE",   #DE_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M",0,        "Dist_DE", 0,        0,        0,        0,         0,        0,        0,       "Coy_DE",0,       0,  0,  0,   0,  0,  "ES-DE",   #DE_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P",0,        "Dist_DE", 0,        0,        0,        0,         0,        0,        0,       "Coy_DE",0,       0,  0,  0,   0,  0,  "ES-DE",   #DE_P
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A",0,        0,        "Dist_DP", 0,        0,        0,         0,        0,        0,       0,       "Coy_DP",0,  0,  0,   0,  0,  0,         #DP_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M",0,        0,        "Dist_DP", 0,        0,        0,         0,        0,        0,       0,       "Coy_DP",0,  0,  0,   0,  0,  0,         #DP_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P",0,        0,        "Dist_DP", 0,        0,        0,         0,        0,        0,       0,       "Coy_DP",0,  0,  0,   0,  0,  0,         #DP_P
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A",0,        0,        0,        "Dist_DR", 0,        0,         0,        0,        0,       0,       0,       0,  0,  0,   0,  0,  0,         #DR_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M",0,        0,        0,        "Dist_DR", 0,        0,         0,        0,        0,       0,       0,       0,  0,  0,   0,  0,  0,         #DR_M
  0,           0,           0,           0,        0,        0,        0,        0,        0,         0,        0,        0,       0,       0,       0,  0,  0,   0,  0,  0,         #DR_P (all-NA)
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A",0,        0,        0,        0,        "Dist_PB", 0,         0,        0,        0,       0,       0,       0,  0,  0,   0,  0,  0,         #PB_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M",0,        0,        0,        0,        "Dist_PB", 0,         0,        0,        0,       0,       0,       0,  0,  0,   0,  0,  0,         #PB_M
  0,           0,           0,           0,        0,        0,        0,        0,        0,         0,        0,        0,       0,       0,       0,  0,  0,   0,  0,  0,         #PB_P (all-NA)
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A",0,        0,        0,        0,        0,        "Dist_PRH", 0,        0,        0,       0,       0,       0,  0,  0,   0,  0,  "ES-PRH",  #PRH_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M",0,        0,        0,        0,        0,        "Dist_PRH", 0,        0,        0,       0,       0,       0,  0,  0,   0,  0,  "ES-PRH",  #PRH_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P",0,        0,        0,        0,        0,        "Dist_PRH", 0,        0,        0,       0,       0,       0,  0,  0,   0,  0,  "ES-PRH",  #PRH_P
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A",0,        0,        0,        0,        0,        0,         "Dist_TB", 0,        0,       0,       0,       0,  0,  0,   0,  0,  0,         #TB_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M",0,        0,        0,        0,        0,        0,         "Dist_TB", 0,        0,       0,       0,       0,  0,  0,   0,  0,  0,         #TB_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P",0,        0,        0,        0,        0,        0,         "Dist_TB", 0,        0,       0,       0,       0,  0,  0,   0,  0,  0,         #TB_P
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A",0,        0,        0,        0,        0,        0,         0,        "Dist_TP", 0,       0,       0,       0,  0,  0,   0,  0,  0,         #TP_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M",0,        0,        0,        0,        0,        0,         0,        "Dist_TP", 0,       0,       0,       0,  0,  0,   0,  0,  0,         #TP_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P",0,        0,        0,        0,        0,        0,         0,        "Dist_TP", 0,       0,       0,       0,  0,  0,   0,  0,  0          #TP_P
), nrow = 24, ncol = 20, byrow = TRUE)

stopifnot(dim(C.8site) == c(24L, 20L))

# ── MODEL A: 24 independent states (expected best) ───────────────────────────
m.A_8site <- run_marss8("A: 24 independent states",
                        factor(1:24), Ut_8site, C.8site, cov)
save(m.A_8site, file = "Output/m.A_8site.RData")

# ── MODEL B: Breeding (18 states) vs Haul-out (6 states) ─────────────────────
# Note: HaulOut_P is unidentifiable (DR_PUP, PB_PUP both all-NA).
# AICc gap between A and B also reflects collapsing 20 site-specific
# covariates to 7 grouped covariates, so gap overstates pure state-structure
# effect. Model C (below) provides the fairer structural comparison.
Z.breed_haulout <- factor(c(
  "Breed_A","Breed_M","Breed_P",       # BL
  "Breed_A","Breed_M","Breed_P",       # DE
  "Breed_A","Breed_M","Breed_P",       # DP
  "HaulOut_A","HaulOut_M","HaulOut_P", # DR
  "HaulOut_A","HaulOut_M","HaulOut_P", # PB
  "Breed_A","Breed_M","Breed_P",       # PRH
  "Breed_A","Breed_M","Breed_P",       # TB
  "Breed_A","Breed_M","Breed_P"))      # TP

U1.bh <- matrix(c("t1_BreedA","t1_BreedM","t1_BreedP",
                  "t1_HaulA", "t1_HaulM", "t1_HaulP"), 6, 1)
U2.bh <- matrix(c("t2_BreedA","t2_BreedM","t1_BreedP",
                  "t2_HaulA", "t2_HaulM", "t1_HaulP"), 6, 1)
Ut.bh <- array(U2.bh, dim = c(dim(U1.bh), TT))
Ut.bh[, , 1:breakpoint] <- U1.bh

cov_bh <- rbind(
  cov["MOCI_JFM", ],
  cov["MOCI_AMJ", ],
  cov["MOCI_OND", ],
  colMeans(cov[c("Dist_BL","Dist_DE","Dist_DP","Dist_PRH","Dist_TB","Dist_TP"), ]),
  colMeans(cov[c("Dist_DR","Dist_PB"), ]),
  colMeans(cov[c("Coyote_BL","Coyote_DE","Coyote_DP"), ]),
  cov["eSeal_Sum_Imm_MaxCount", ]
)
rownames(cov_bh) <- c("MOCI_JFM","MOCI_AMJ","MOCI_OND",
                      "Dist_Breed","Dist_Haulout","Coy_Breed","eSeal")

C.bh <- matrix(list(
  "MOCI_JFM","MOCI_AMJ","MOCI_OND","Dist_Breed", 0,            "Coy_Breed","eSeal",    #Breed_A
  "MOCI_JFM","MOCI_AMJ","MOCI_OND","Dist_Breed", 0,            "Coy_Breed","eSeal",    #Breed_M
  "MOCI_JFM","MOCI_AMJ","MOCI_OND","Dist_Breed", 0,            "Coy_Breed","eSeal",    #Breed_P
  "MOCI_JFM","MOCI_AMJ","MOCI_OND", 0,           "Dist_Haulout",0,          0,         #HaulOut_A
  "MOCI_JFM","MOCI_AMJ","MOCI_OND", 0,           "Dist_Haulout",0,          0,         #HaulOut_M
  "MOCI_JFM","MOCI_AMJ","MOCI_OND", 0,           0,             0,          0          #HaulOut_P (unidentified)
), nrow = 6, ncol = 7, byrow = TRUE)

m.B_8site <- run_marss8("B: Breeding (18) vs Haul-out (6 grouped) [5 identifiable]",
                        Z.breed_haulout, Ut.bh, C.bh, cov_bh)
save(m.B_8site, file = "Output/m.B_8site.RData")

# ── MODEL C: Estuary/Bay vs Coastal/Exposed (6 states) — NEW ─────────────────
# Ecological grouping based on tidal shelter and oceanographic exposure:
#   Estuary/Bay (Est):   BL (Bolinas Lagoon), DE (Drakes Estero), TB (Tomales Bay)
#   Coastal/Exposed (Coast): DP (Double Point), DR (Duxbury Reef), PB (Point Bonita),
#                            PRH (Point Reyes Headlands), TP (Tomales Point)
#
# DR and PB are coastal rocky haul-outs; DR_PUP and PB_PUP map to Coast_P
# but contribute no data (all-NA), so Coast_P is identified by DP/PRH/TP pups.
#
# Covariate pooling decisions (logged for transparency):
#   Dist_Est  = mean(Dist_BL, Dist_DE, Dist_TB)
#   Dist_Coast = mean(Dist_DP, Dist_DR, Dist_PB, Dist_PRH, Dist_TP)
#   Coy_Est   = mean(Coyote_BL, Coyote_DE)  [TB has zero coyote]
#   Coy_Coast = Coyote_DP only              [PRH/TP/DR/PB all zero]
#   eSeal_Est   = eSeal at DE (estuary co-occurrence)
#   eSeal_Coast = eSeal at PRH (coastal co-occurrence)
#   Both use the same regional eSeal index; separate coefficients distinguish
#   the response at the two site types.
#
# The 7-covariate design allows a direct AICc comparison with the 7-covariate
# Model B and is structurally equivalent in parameter count for that comparison.

Z.est_coast <- factor(c(
  "Est_A","Est_M","Est_P",         # BL  (estuary/lagoon)
  "Est_A","Est_M","Est_P",         # DE  (estuary)
  "Coast_A","Coast_M","Coast_P",   # DP  (coastal breeding)
  "Coast_A","Coast_M","Coast_P",   # DR  (coastal haul-out; DR_PUP all-NA)
  "Coast_A","Coast_M","Coast_P",   # PB  (coastal haul-out; PB_PUP all-NA)
  "Coast_A","Coast_M","Coast_P",   # PRH (coastal breeding)
  "Est_A","Est_M","Est_P",         # TB  (bay/estuary)
  "Coast_A","Coast_M","Coast_P"    # TP  (coastal breeding)
))

U1.ec <- matrix(c("t1_EstA","t1_EstM","t1_EstP",
                  "t1_CoastA","t1_CoastM","t1_CoastP"), 6, 1)
U2.ec <- matrix(c("t2_EstA","t2_EstM","t1_EstP",
                  "t2_CoastA","t2_CoastM","t1_CoastP"), 6, 1)
Ut.ec <- array(U2.ec, dim = c(dim(U1.ec), TT))
Ut.ec[, , 1:breakpoint] <- U1.ec

# 9-covariate grouped matrix (class-specific MOCI, grouped Dist, site Coy, split eSeal)
cov_ec <- rbind(
  MOCI_JFM   = cov["MOCI_JFM", ],
  MOCI_AMJ   = cov["MOCI_AMJ", ],
  MOCI_OND   = cov["MOCI_OND", ],
  Dist_Est   = colMeans(cov[c("Dist_BL","Dist_DE","Dist_TB"), ]),
  Dist_Coast = colMeans(cov[c("Dist_DP","Dist_DR","Dist_PB","Dist_PRH","Dist_TP"), ]),
  Coy_Est    = colMeans(cov[c("Coyote_BL","Coyote_DE"), ]),
  Coy_Coast  = cov["Coyote_DP", ],
  eSeal_Est  = cov["eSeal_Sum_Imm_MaxCount", ],
  eSeal_Coast = cov["eSeal_Sum_Imm_MaxCount", ]
)

C.ec <- matrix(list(
  # Each row = one latent state; columns = 9 grouped covariates in cov_ec row order
  # MOCI_JFM  MOCI_AMJ  MOCI_OND  Dist_Est   Dist_Coast  Coy_Est    Coy_Coast  eSeal_Est  eSeal_Coast
  "moci_A", "moci_A", "moci_A", "dst_E",   0,          "coy_E",   0,         "es_E",    0,          # Est_A
  "moci_M", "moci_M", "moci_M", "dst_E",   0,          "coy_E",   0,         "es_E",    0,          # Est_M
  "moci_P", "moci_P", "moci_P", "dst_E",   0,          "coy_E",   0,         "es_E",    0,          # Est_P
  "moci_A", "moci_A", "moci_A", 0,        "dst_C",     0,        "coy_C",    0,        "es_C",      # Coast_A
  "moci_M", "moci_M", "moci_M", 0,        "dst_C",     0,        "coy_C",    0,        "es_C",      # Coast_M
  "moci_P", "moci_P", "moci_P", 0,        "dst_C",     0,        "coy_C",    0,        "es_C"       # Coast_P
), nrow = 6, ncol = 9, byrow = TRUE)

m.C_8site <- run_marss8("C: Estuary/Bay vs Coastal/Exposed (6 grouped states)",
                        Z.est_coast, Ut.ec, C.ec, cov_ec)
save(m.C_8site, file = "Output/m.C_8site.RData")

# ── Model comparison ─────────────────────────────────────────────────────────
# NOTE on AICc comparability: Models A and C share comparable covariate designs
# (full site-specific in A; ecologically-grouped in C). Model B collapses both
# state structure AND covariates; its larger AICc gap partly reflects the
# coarser covariate specification, not purely the state grouping.
df_aic_8site <- data.frame(
  model = c(
    "A: 24 independent states (site×class)",
    "B: Breeding (18) vs Haul-out (6) [5 identifiable]",
    "C: Estuary/Bay (9) vs Coastal (15) [6 states]"),
  AICc       = c(m.A_8site$AICc,  m.B_8site$AICc,  m.C_8site$AICc),
  num_params = c(m.A_8site$num.params, m.B_8site$num.params, m.C_8site$num.params),
  logLik     = c(m.A_8site$logLik, m.B_8site$logLik, m.C_8site$logLik)
) |>
  transform(deltaAICc = AICc - min(AICc)) |>
  (\(d) d[order(d$AICc), ])()

cat("\n── Model Comparison ─────────────────────────────────────────────────────\n")
print(df_aic_8site)
write.csv(df_aic_8site, "Output/model_comparison_8site_AIC.csv", row.names = FALSE)

# ── 89% CIs for the best (independent) model ─────────────────────────────────
cat("\n── Computing 89% CIs for m.A_8site (alpha = 0.11) ──────────────────────\n")
CIs_8site <- MARSSparamCIs(m.A_8site, alpha = 0.11)
save(CIs_8site, file = "Output/CIs_8site.RData")
save(df_aic_8site, file = "Output/df_aic_8site.RData")

cat("\nObjects: m.A_8site, m.B_8site, m.C_8site, CIs_8site, df_aic_8site\n")
cat("Next: source('Code/11_marss8_plots.R')\n")
