# ============================================================================
# HARBOR SEAL IPM v3.2 — 00_load.R
# Loads all project packages and prints a startup banner.
# Called automatically from .Rprofile, OR run manually:
#   source("00_load.R")
# ============================================================================

.pkgs <- c("plyr", "cmdstanr", "posterior", "bayesplot", "tidyverse", "lubridate", "patchwork", "scales", "viridis", "RColorBrewer", "here", "readxl", "loo")

.failed <- character(0)
for (.p in .pkgs) {
  .ok <- suppressPackageStartupMessages(
    tryCatch({ library(.p, character.only=TRUE, quietly=TRUE); TRUE },
             error = function(e) FALSE))
  if (isFALSE(.ok)) .failed <- c(.failed, .p)
}

# Create standard output directories
for (.d in c("Output", "Output/Plots", "Code", "Data"))
  if (!dir.exists(.d)) dir.create(.d, recursive=TRUE, showWarnings=FALSE)

# Banner
cat("\n══════════════════════════════════════════════════════\n")
cat("  Harbor Seal IPM v3.2\n")
if (length(.failed) > 0) {
  cat("  Missing packages — run source(\"harbor_seal_setup.R\"):\n")
  cat("   ", paste(.failed, collapse=", "), "\n")
} else {
  .cs <- tryCatch(cmdstanr::cmdstan_version(), error=function(e) "not found")
  cat(sprintf("  CmdStan  : v%s\n", .cs))
  cat(sprintf("  tidyverse: v%s  |  posterior: v%s\n",
              as.character(packageVersion("tidyverse")),
              as.character(packageVersion("posterior"))))
  cat(sprintf("  mc.cores : %d\n", getOption("mc.cores", 1L)))
  cat(sprintf("  Loaded   : %s\n", paste(.pkgs[!.pkgs %in% .failed], collapse=", ")))
}
cat("\n  Quick start:\n")
cat("    source(\"Code/harbor_seal_ipm_v3.2.R\")\n")
cat("    out <- load_seal_results(\"IPM_v3.2_real\")\n")
cat("══════════════════════════════════════════════════════\n\n")
rm(.pkgs, .failed, .ok, .p, .cs, .d)

