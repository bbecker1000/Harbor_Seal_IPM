# ============================================================================
# 03_marss6_models.R
# ----------------------------------------------------------------------------
# 6-SITE MARSS — FRAMEWORK / MODEL-SELECTION ANALYSIS (1997–2025)
# Purpose: establish the latent-state structure that justifies the IPM design.
# Compares four spatial-population hypotheses via AICc, then computes 89% CIs
# for the best model. These results set up the framework; the IPM does the
# demographic heavy lifting.
#
# Hypotheses (Z structures), all sharing the class-specific C matrix and a
# 2004 breakpoint in adult/molt growth (pups constant):
#   A  Independent     — 18 states, one per site×class   (expected best)
#   B  Estuary/Outer   — 6 states  (BL/DE/TB vs DP/PRH/TP × class)
#   C  One population   — 3 states  (class only, pooled across sites)
#   D  Molt separate    — 3 states  (Breed_A, Molt, Breed_P)
#
# Prereqs (run in order):
#   source("Code/01_data_prep_6site.R")   # -> dat (18 x 29), years
#   source("Code/02_covariates_6site.R")  # -> cov_t_scaled (16 x 29)
#
# Produces / saves:
#   m.A.indep, m.B.estuary, m.C.one, m.D.molt
#   BESTMODEL (= m.A.indep), CIs  -> Output/marss6_best.RData
#   df_aic_6site -> Output/model_comparison_6site_AIC.csv
#
# Next: source("Code/04_marss6_plots.R")
# ============================================================================

library(MARSS)
library(tidyverse)
library(broom)   # tidy() method for marssMLE/marssParamCIs (used by plots too)

dir.create("Output", showWarnings = FALSE)

# ── Guards: required inputs present and correctly shaped ─────────────────────
if (!exists("dat"))          stop("Run 01_data_prep_6site.R first (need `dat`).")
if (!exists("cov_t_scaled")) stop("Run 02_covariates_6site.R first (need `cov_t_scaled`).")

SITES_6     <- c("BL", "DE", "DP", "PRH", "TB", "TP")
CANON_STATE <- as.vector(t(outer(SITES_6, c("ADULT","MOLTING","PUP"),
                                 paste, sep = "_")))
CANON_COV   <- c("MOCI_JFM","MOCI_AMJ","MOCI_OND",
                 paste0("Dist_",   SITES_6),
                 paste0("Coyote_", SITES_6),
                 "eSeal_Sum_Imm_MaxCount")

stopifnot(
  nrow(dat) == 18L, nrow(cov_t_scaled) == 16L,
  ncol(dat) == ncol(cov_t_scaled),
  identical(rownames(dat), CANON_STATE),         # C-matrix rows depend on this
  identical(rownames(cov_t_scaled), CANON_COV)   # C-matrix cols depend on this
)

cov_t_marss <- cov_t_scaled
TT          <- ncol(dat)
breakpoint  <- ceiling(TT / 4)
cat(sprintf("Data %d x %d | breakpoint col %d = year %d\n",
            nrow(dat), TT, breakpoint, years[breakpoint]))

# ── Time-varying U: adults & molt shift at 2004; pups constant ──────────────
# Class order within each site is ADULT, MOLTING, PUP (canonical).
make_Ut <- function(u1_vec, u2_vec, n_blocks) {
  U1 <- matrix(rep(u1_vec, n_blocks), ncol = 1)
  U2 <- matrix(rep(u2_vec, n_blocks), ncol = 1)
  Ut <- array(U2, dim = c(nrow(U1), 1, TT))
  Ut[, , 1:breakpoint] <- U1
  Ut
}
Ut.Class <- make_Ut(c("t1_A","t1_M","t1_P"), c("t2_A","t2_M","t1_P"), 6)

# ── Shared C matrix: 18 states × 16 covariates ──────────────────────────────
# Class-specific MOCI (A/M/P); site-specific Dist & Coy; eSeal at DE & PRH only.
# Column order matches CANON_COV exactly.
C.model.new <- matrix(list(
  #             JFM           AMJ           OND          BL        DE        DP        PRH        TB        TP        cBL       cDE       cDP       cPRH cTB cTP eSeal
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A","Dist_BL", 0,        0,        0,         0,        0,        "Coy_BL", 0,        0,        0,   0,  0,  0,        #BL_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M","Dist_BL", 0,        0,        0,         0,        0,        "Coy_BL", 0,        0,        0,   0,  0,  0,        #BL_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P","Dist_BL", 0,        0,        0,         0,        0,        "Coy_BL", 0,        0,        0,   0,  0,  0,        #BL_P
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A",0,        "Dist_DE", 0,        0,         0,        0,        0,        "Coy_DE", 0,        0,   0,  0,  "ES-DE",  #DE_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M",0,        "Dist_DE", 0,        0,         0,        0,        0,        "Coy_DE", 0,        0,   0,  0,  "ES-DE",  #DE_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P",0,        "Dist_DE", 0,        0,         0,        0,        0,        "Coy_DE", 0,        0,   0,  0,  "ES-DE",  #DE_P
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A",0,        0,        "Dist_DP", 0,         0,        0,        0,        0,        "Coy_DP", 0,   0,  0,  0,        #DP_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M",0,        0,        "Dist_DP", 0,         0,        0,        0,        0,        "Coy_DP", 0,   0,  0,  0,        #DP_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P",0,        0,        "Dist_DP", 0,         0,        0,        0,        0,        "Coy_DP", 0,   0,  0,  0,        #DP_P
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A",0,        0,        0,        "Dist_PRH", 0,        0,        0,        0,        0,        0,   0,  0,  "ES-PRH", #PRH_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M",0,        0,        0,        "Dist_PRH", 0,        0,        0,        0,        0,        0,   0,  0,  "ES-PRH", #PRH_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P",0,        0,        0,        "Dist_PRH", 0,        0,        0,        0,        0,        0,   0,  0,  "ES-PRH", #PRH_P
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A",0,        0,        0,        0,         "Dist_TB", 0,        0,        0,        0,        0,   0,  0,  0,        #TB_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M",0,        0,        0,        0,         "Dist_TB", 0,        0,        0,        0,        0,   0,  0,  0,        #TB_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P",0,        0,        0,        0,         "Dist_TB", 0,        0,        0,        0,        0,   0,  0,  0,        #TB_P
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A",0,        0,        0,        0,         0,        "Dist_TP", 0,        0,        0,        0,   0,  0,  0,        #TP_A
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M",0,        0,        0,        0,         0,        "Dist_TP", 0,        0,        0,        0,   0,  0,  0,        #TP_M
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P",0,        0,        0,        0,         0,        "Dist_TP", 0,        0,        0,        0,   0,  0,  0         #TP_P
), nrow = 18, ncol = 16, byrow = TRUE)

stopifnot(dim(C.model.new) == c(18L, 16L))

# ── Common MARSS settings ───────────────────────────────────────────────────
R_fixed <- diag(0.025, 18)
ctrl    <- list(maxit = 5000, safe = TRUE, trace = 0, allow.degen = TRUE)

run_marss <- function(label, Z, U, C, cmat) {
  cat(sprintf("\n── %s ─────────────────────────────────\n", label))
  t0 <- Sys.time()
  m <- MARSS(dat, model = list(Z = Z, U = U, R = R_fixed,
                               Q = "diagonal and equal", B = "identity",
                               C = C, c = cmat, tinitx = 1),
             control = ctrl)
  cat(sprintf("  %s | %.1fs | AICc=%.2f | params=%d\n",
              label, as.numeric(difftime(Sys.time(), t0, units="secs")),
              m$AICc, m$num.params))
  m
}

# ── Model A: 18 independent states (expected best) ──────────────────────────
m.A.indep <- run_marss("A: Independent (18 states)",
                       factor(1:18), Ut.Class, C.model.new, cov_t_marss)
save(m.A.indep, file = "Output/m.A.indep_6site.RData")

# ── Model B: Estuary (BL/DE/TB) vs Outer (DP/PRH/TP), 6 states ──────────────
Z.estuary <- factor(c(
  "Est_A","Est_M","Est_P",   # BL
  "Est_A","Est_M","Est_P",   # DE
  "Out_A","Out_M","Out_P",   # DP
  "Out_A","Out_M","Out_P",   # PRH
  "Est_A","Est_M","Est_P",   # TB
  "Out_A","Out_M","Out_P"))  # TP
Ut.estuary <- make_Ut(c("t1_Est_A","t1_Est_M","t1_Est_P",
                        "t1_Out_A","t1_Out_M","t1_Out_P"),
                      c("t2_Est_A","t2_Est_M","t1_Est_P",
                        "t2_Out_A","t2_Out_M","t1_Out_P"), 1)
# Grouped C: class-only effects (MOCI/Dist/Coy/ES pooled within group×class).
C.model.estuary <- matrix(list(
  "MOCI_A","MOCI_A","MOCI_A","Dist_A","Dist_A",0,       0,       "Dist_A",0,       "Coy_A","Coy_A",0,      0,      "Coy_A",0,      "ES_A",  #Est_A
  "MOCI_M","MOCI_M","MOCI_M","Dist_M","Dist_M",0,       0,       "Dist_M",0,       "Coy_M","Coy_M",0,      0,      "Coy_M",0,      "ES_M",  #Est_M
  "MOCI_P","MOCI_P","MOCI_P","Dist_P","Dist_P",0,       0,       "Dist_P",0,       "Coy_P","Coy_P",0,      0,      "Coy_P",0,      "ES_P",  #Est_P
  "MOCI_A","MOCI_A","MOCI_A",0,       0,       "Dist_A","Dist_A",0,       "Dist_A",0,      0,      "Coy_A","Coy_A",0,      "Coy_A","ES_A",  #Out_A
  "MOCI_M","MOCI_M","MOCI_M",0,       0,       "Dist_M","Dist_M",0,       "Dist_M",0,      0,      "Coy_M","Coy_M",0,      "Coy_M","ES_M",  #Out_M
  "MOCI_P","MOCI_P","MOCI_P",0,       0,       "Dist_P","Dist_P",0,       "Dist_P",0,      0,      "Coy_P","Coy_P",0,      "Coy_P","ES_P"   #Out_P
), nrow = 6, ncol = 16, byrow = TRUE)
m.B.estuary <- run_marss("B: Estuary vs Outer (6 states)",
                         Z.estuary, Ut.estuary, C.model.estuary, cov_t_marss)
save(m.B.estuary, file = "Output/m.B.estuary_6site.RData")

# ── Model C: one population, 3 class states ─────────────────────────────────
Z.one  <- factor(rep(c("A","M","P"), 6))
Ut.one <- make_Ut(c("t1_A","t1_M","t1_P"), c("t2_A","t2_M","t1_P"), 1)
C.model.one <- matrix(list(
  "MOCI_JFM_A","MOCI_AMJ_A","MOCI_OND_A","Dist_A","Dist_A","Dist_A","Dist_A","Dist_A","Dist_A","Coy_A","Coy_A","Coy_A","Coy_A","Coy_A","Coy_A","ES_A",
  "MOCI_JFM_M","MOCI_AMJ_M","MOCI_OND_M","Dist_M","Dist_M","Dist_M","Dist_M","Dist_M","Dist_M","Coy_M","Coy_M","Coy_M","Coy_M","Coy_M","Coy_M","ES_M",
  "MOCI_JFM_P","MOCI_AMJ_P","MOCI_OND_P","Dist_P","Dist_P","Dist_P","Dist_P","Dist_P","Dist_P","Coy_P","Coy_P","Coy_P","Coy_P","Coy_P","Coy_P","ES_P"
), nrow = 3, ncol = 16, byrow = TRUE)
m.C.one <- run_marss("C: One population (3 states)",
                     Z.one, Ut.one, C.model.one, cov_t_marss)
save(m.C.one, file = "Output/m.C.one_6site.RData")

# ── Model D: molt separate from adults/pups, 3 states ───────────────────────
Z.molt  <- factor(rep(c("Breed_A","Molt","Breed_P"), 6))
Ut.molt <- make_Ut(c("t1_Breed_A","t1_Molt","t1_Breed_P"),
                   c("t2_Breed_A","t2_Molt","t1_Breed_P"), 1)
C.model.molt <- matrix(list(
  "MOCI_JFM_Breed","MOCI_AMJ_Breed","MOCI_OND_Breed","Dist_Breed","Dist_Breed","Dist_Breed","Dist_Breed","Dist_Breed","Dist_Breed","Coy_Breed","Coy_Breed","Coy_Breed","Coy_Breed","Coy_Breed","Coy_Breed","ES_Breed",
  "MOCI_JFM_Molt", "MOCI_AMJ_Molt", "MOCI_OND_Molt", "Dist_Molt", "Dist_Molt", "Dist_Molt", "Dist_Molt", "Dist_Molt", "Dist_Molt", "Coy_Molt", "Coy_Molt", "Coy_Molt", "Coy_Molt", "Coy_Molt", "Coy_Molt", "ES_Molt",
  "MOCI_JFM_Breed","MOCI_AMJ_Breed","MOCI_OND_Breed","Dist_Breed","Dist_Breed","Dist_Breed","Dist_Breed","Dist_Breed","Dist_Breed","Coy_Breed","Coy_Breed","Coy_Breed","Coy_Breed","Coy_Breed","Coy_Breed","ES_Breed"
), nrow = 3, ncol = 16, byrow = TRUE)
m.D.molt <- run_marss("D: Molt separate (3 states)",
                      Z.molt, Ut.molt, C.model.molt, cov_t_marss)
save(m.D.molt, file = "Output/m.D.molt_6site.RData")

# ── Model comparison (single table) ─────────────────────────────────────────
df_aic_6site <- tibble(
  model      = c("A: Independent (18 states)",
                 "B: Estuary vs Outer (6 states)",
                 "C: One population (3 states)",
                 "D: Molt separate (3 states)"),
  AICc       = c(m.A.indep$AICc, m.B.estuary$AICc, m.C.one$AICc, m.D.molt$AICc),
  num_params = c(m.A.indep$num.params, m.B.estuary$num.params,
                 m.C.one$num.params, m.D.molt$num.params),
  logLik     = c(m.A.indep$logLik, m.B.estuary$logLik,
                 m.C.one$logLik, m.D.molt$logLik)
) %>%
  mutate(deltaAICc = AICc - min(AICc)) %>%
  arrange(AICc)

cat("\n── Model comparison ────────────────────────────────────────────────\n")
print(df_aic_6site)
write_csv(df_aic_6site, "Output/model_comparison_6site_AIC.csv")

# ── Best model + 89% CIs (consumed by 04_marss6_plots.R) ────────────────────
# Model A (independent states) is the framework model carried forward; it
# justifies treating sites independently in the IPM.
BESTMODEL <- m.A.indep
cat("\n── Computing 89% CIs for BESTMODEL (alpha = 0.11) ──────────────────\n")
CIs <- MARSSparamCIs(BESTMODEL, alpha = 0.11)

save(BESTMODEL, CIs, df_aic_6site, file = "Output/marss6_best.RData")
cat("Saved: Output/marss6_best.RData (BESTMODEL, CIs, df_aic_6site)\n")

cat("\nObjects created: m.A.indep, m.B.estuary, m.C.one, m.D.molt, BESTMODEL, CIs\n")
cat("Next: source(\"Code/04_marss6_plots.R\")\n")