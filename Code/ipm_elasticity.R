#####

library(posterior)
library(dplyr)

drws <- out$fit$draws(format = "df")

# ── 1. AGGREGATE LAMBDA FROM N STATE VARIABLES ──────────────────────────────
# N[site, year], years 1–29 = 1997–2025
# Total N per year = sum across sites of (N_adult_F + N_adult_M + N_juv_F + N_juv_M + N_pup)
# We'll use adult F + adult M + juv (combined) + pup as proxies since juv_F/M aren't separated in draws

S <- 6; T <- 29

get_N_total <- function(drws, t) {
  aF  <- rowMeans(drws[, paste0("N_adult_F[", 1:S, ",", t, "]")])
  aM  <- rowMeans(drws[, paste0("N_adult_M[", 1:S, ",", t, "]")])
  # Check if juv_F and juv_M exist separately
  jF_vars <- paste0("N_juv_F[", 1:S, ",", t, "]")
  jM_vars <- paste0("N_juv_M[", 1:S, ",", t, "]")
  pup_vars <- paste0("N_pup[",   1:S, ",", t, "]")
  jF  <- if(all(jF_vars %in% names(drws))) rowMeans(drws[, jF_vars]) else 0
  jM  <- if(all(jM_vars %in% names(drws))) rowMeans(drws[, jM_vars]) else 0
  pup <- if(all(pup_vars %in% names(drws))) rowMeans(drws[, pup_vars]) else 0
  aF + aM + jF + jM + pup
}

N_by_year <- sapply(1:T, get_N_total, drws = drws)  # draws × years matrix

# Lambda = N[t+1] / N[t], years 1998–2025 = 28 transitions
lambda_draws <- N_by_year[, 2:T] / N_by_year[, 1:(T-1)]  # draws × 28 years

# Summary
lambda_annual_median <- apply(lambda_draws, 2, median)
lambda_all           <- as.vector(lambda_draws)

cat("── Aggregate λ summary ──\n")
cat("Mean λ:           ", round(mean(lambda_annual_median), 3), "\n")
cat("89% CI:            ", round(quantile(lambda_annual_median, 0.055), 3), "–",
    round(quantile(lambda_annual_median, 0.945), 3), "\n")
cat("Range (min–max):  ", round(min(lambda_annual_median), 3), "–",
    round(max(lambda_annual_median), 3), "\n")
cat("Years above 1.0:  ", sum(lambda_annual_median > 1), "of 28\n")
cat("Years below 1.0:  ", sum(lambda_annual_median < 1), "of 28\n")

# ── 2. ELASTICITY FROM POSTERIOR MEDIAN VITAL RATES ─────────────────────────
# Lefkovitch matrix (females only, pre-breeding census)
# States: pup, juv_F, adult_F  (male states contribute symmetrically but
# elasticity is typically reported for the female sub-matrix)

phi_pup   <- median(drws$phi_pup_base)
phi_juv   <- median(drws$phi_juv_base)
phi_aF    <- median(drws$phi_adult_F_base)
f_avg     <- median(drws$avg_fecundity)
pf        <- 0.50   # prop female

# Lefkovitch (3×3 female sub-matrix):
#           pup      juv_F    adult_F
# pup    [  0        0        f_avg  ]
# juv_F  [ phi_pup*pf  2/3*phi_juv  0 ]
# adult_F[  0        1/3*phi_juv  phi_aF ]

A <- matrix(c(
  0,            0,              f_avg,
  phi_pup * pf, (2/3)*phi_juv,  0,
  0,            (1/3)*phi_juv,  phi_aF
), nrow=3, byrow=TRUE)

# Eigenanalysis
ev    <- eigen(A)
lambda_det <- Re(ev$values[1])
w     <- Re(ev$vectors[, 1]); w <- w / sum(w)   # stable stage (col)
vt    <- Re(eigen(t(A))$vectors[, 1])            # reproductive value (row)
vt    <- vt / vt[1]

# Elasticity = (v_i * w_j / <v,w>) * (a_ij / lambda)
vw    <- sum(vt * w)
elas  <- (outer(vt, w) / vw) * (A / lambda_det)

cat("\n── Deterministic λ from posterior medians ──\n")
cat("λ =", round(lambda_det, 4), "\n")

cat("\n── Elasticity matrix (rows=to, cols=from) ──\n")
rownames(elas) <- colnames(elas) <- c("pup","juv_F","adult_F")
print(round(elas, 4))

cat("\n── Summed elasticities by vital rate ──\n")
cat("ε_adult survival (phi_aF):  ", round(elas["adult_F","adult_F"], 3), "\n")
cat("ε_juv survival (phi_juv):   ", round(elas["adult_F","juv_F"] + elas["juv_F","juv_F"], 3), "\n")
cat("ε_pup survival (phi_pup):   ", round(elas["juv_F","pup"], 3), "\n")
cat("ε_fecundity (f_avg):        ", round(elas["pup","adult_F"], 3), "\n")