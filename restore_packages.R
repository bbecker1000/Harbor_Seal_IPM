# restore_packages.R — install project packages directly
# Bypasses renv restore to avoid V8/system library issues

options(repos = c(
  RSPM = "https://packagemanager.posit.co/cran/__linux__/noble/latest",
  CRAN = "https://cloud.r-project.org"))

pkgs <- c("posterior","bayesplot","tidyverse","lubridate","plyr",
          "patchwork","scales","viridis","RColorBrewer","here",
          "readxl","loo","mgcv","knitr","kableExtra","abind","progress")

cat("Installing packages...\n")
install.packages(pkgs)
install.packages("cmdstanr",
  repos=c("https://mc-stan.org/r-packages/","https://cloud.r-project.org"))

cat("Done — restarting R...\n")
.rs.restartR()
