# Harbor Seal IPM v3.2 — restore_packages.R
# Run this after any container restart if packages are missing
# Usage: source("restore_packages.R")

cat("Restoring packages from renv.lock...\n")
renv::restore(prompt=FALSE)
cat("Done — restarting R...\n")
.rs.restartR()
