# ============================================================================
# 7b_models_legacy_MARSS.R
# ============================================================================
# MARSS TREND ANALYSIS — LONG-TERM TIME SERIES (1975–2025)
# No covariates. Purpose: establish pre-1997 baseline and long-term trajectory.
#
# Source first:
#   source("7a_data_prep_legacy_MARSS.R")
#
# Model A: All independent states, time-varying U (2004 breakpoint)
# Model B: All independent states, constant U (no breakpoint)
# Model comparison via AICc; CIs computed for winning model.
#
# Produces: m.LA, m.LB, CIs_legacy, df_aic_legacy
# ============================================================================

library(MARSS)

if (!exists("dat_legacy"))  stop("Run 7a_data_prep_legacy_MARSS.R first")

n_states  <- nrow(dat_legacy)
TT        <- ncol(dat_legacy)
cat("States:", n_states, "| Years:", TT,
    "(", min(years_legacy), "–", max(years_legacy), ")\n")

# ── Helper: parse class from state name ──────────────────────────────────────
get_class <- function(nm) {
  if (grepl("MOLTING", nm, ignore.case = TRUE)) "M"
  else if (grepl("ADULT", nm, ignore.case = TRUE)) "A"
  else "P"
}
classes <- sapply(state_names_legacy, get_class)

# ── Time-varying U: 2-phase model with 2004 breakpoint ───────────────────────
# Phase 1 (pre-2004): shared growth rate within class (t1_A, t1_M, t1_P)
# Phase 2 (2004+):    shared growth rate within class (t2_A, t2_M, t2_P)
# All three classes now allowed to shift at the 2004 breakpoint.
# Pups previously held constant (t1_P throughout) based on 8-site analysis
# (1997-2025 only); the longer legacy time series justifies estimating a
# separate Phase 2 pup trend.

bp_year <- 2004
bp_col  <- which(years_legacy == bp_year)

if (length(bp_col) == 0) {
  warning("Breakpoint year ", bp_year, " not in years_legacy — using constant U")
  bp_col <- TT + 1
}
cat("2004 breakpoint at column", bp_col,
    "(year", years_legacy[min(bp_col, TT)], ")\n")

U1 <- matrix(paste0("t1_", classes), n_states, 1,
             dimnames = list(state_names_legacy, NULL))

U2 <- matrix(paste0("t2_", classes), n_states, 1,  # all classes shift
             dimnames = list(state_names_legacy, NULL))

Ut <- array(U2, dim = c(n_states, 1, TT))
if (bp_col > 1) Ut[, , 1:(bp_col - 1)] <- U1

cat("Phase 1 U parameters:\n"); print(unique(U1))
cat("Phase 2 U parameters:\n"); print(unique(U2))

# ── Shared model components ───────────────────────────────────────────────────
common_model <- list(
  Z       = factor(seq_len(n_states)),  # all independent states
  R       = diag(0.025, n_states),      # fixed obs error (same as 8-site)
  Q       = "diagonal and equal",       # shared process variance
  B       = "identity",                 # random walk with drift
  C       = "zero",                     # no covariates
  c       = "zero",
  tinitx  = 1
)
ctrl <- list(maxit = 5000, safe = TRUE, trace = 0, allow.degen = TRUE)

# ── Model A: time-varying U (2-phase, 2004 breakpoint) ───────────────────────
cat("\n── Model A: independent states + 2-phase U (breakpoint", bp_year, ") ────\n")
t0 <- Sys.time()
m.LA <- MARSS(dat_legacy,
              model   = c(common_model, list(U = Ut)),
              control = ctrl)
cat("Run time:", round(difftime(Sys.time(), t0, units = "secs"), 1), "sec\n")
cat("Convergence:", m.LA$convergence, "| Iterations:", m.LA$numIter, "\n")
cat("AICc:", round(m.LA$AICc, 2), "| logLik:", round(m.LA$logLik, 3), "\n")
save(m.LA, file = "Output/m.LA_legacy.RData")

# ── Model B: constant U (no breakpoint) ──────────────────────────────────────
cat("\n── Model B: independent states + constant U ─────────────────────────────\n")
U_const <- matrix(paste0("u_", classes), n_states, 1)

t0 <- Sys.time()
m.LB <- MARSS(dat_legacy,
              model   = c(common_model, list(U = U_const)),
              control = ctrl)
cat("Run time:", round(difftime(Sys.time(), t0, units = "secs"), 1), "sec\n")
cat("Convergence:", m.LB$convergence, "| Iterations:", m.LB$numIter, "\n")
cat("AICc:", round(m.LB$AICc, 2), "| logLik:", round(m.LB$logLik, 3), "\n")
save(m.LB, file = "Output/m.LB_legacy.RData")

# ── Model comparison ─────────────────────────────────────────────────────────
df_aic_legacy <- data.frame(
  Model      = c("A: 2-phase U (2004 breakpoint)",
                 "B: Constant U (no breakpoint)"),
  AICc       = c(m.LA$AICc,       m.LB$AICc),
  num_params = c(m.LA$num.params, m.LB$num.params),
  logLik     = c(m.LA$logLik,     m.LB$logLik)
) %>%
  mutate(deltaAICc = AICc - min(AICc)) %>%
  arrange(AICc)

cat("\n── Model Comparison ──────────────────────────────────────────────────────\n")
print(df_aic_legacy)
write.csv(df_aic_legacy, "Output/model_comparison_legacy_AIC.csv", row.names = FALSE)

BEST_LEGACY <- if (m.LA$AICc <= m.LB$AICc) m.LA else m.LB
cat("\nBest model:", df_aic_legacy$Model[1], "\n")

# ── Growth rate summary ───────────────────────────────────────────────────────
cat("\n── Estimated trend parameters (best model) ──────────────────────────────\n")
u_ests <- BEST_LEGACY$par$U
cat("Phase 1 adult growth rate (t1_A):",
     round(u_ests["t1_A", 1], 4), "\n")
cat("Phase 1 molt  growth rate (t1_M):",
    round(u_ests["t1_M", 1], 4), "\n")
cat("Phase 1 pup   growth rate (t1_P):",
    round(u_ests["t1_P", 1], 4), "\n")
cat("Phase 2 adult growth rate (t2_A):",
    round(u_ests["t2_A", 1], 4), "\n")
cat("Phase 2 molt  growth rate (t2_M):",
    round(u_ests["t2_M", 1], 4), "\n")
cat("Phase 2 pup   growth rate (t2_P):",
    round(u_ests["t2_P", 1], 4), "\n")
cat("Process variance (Q):",
    round(BEST_LEGACY$par$Q[1], 5), "\n")

# ── 89% CIs for best model ────────────────────────────────────────────────────
cat("\n── Computing 89% CIs ────────────────────────────────────────────────────\n")
CIs_legacy <- MARSSparamCIs(BEST_LEGACY, alpha = 0.11)  # alpha=0.11 → 89% CI
save(CIs_legacy, file = "Output/CIs_legacy.RData")

cat("\nU parameter CIs:\n")
u_ci <- tidy(CIs_legacy) %>% filter(grepl("^U\\.", term))
print(u_ci)

cat("\nObjects created: m.LA, m.LB, CIs_legacy, df_aic_legacy\n")
cat("Source 14_marss_legacy_plots.R to generate plots.\n")

