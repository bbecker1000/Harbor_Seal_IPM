# Harbor Seal IPM v3.2 — .Rprofile
# Lightweight startup — loads minimum to avoid crash
# Run source("00_load.R") manually to load all packages

# renv
try(source("renv/activate.R"), silent=TRUE)

# Options
options(warn=1, scipen=999)
options(mc.cores=4L)

# CmdStan path
if (requireNamespace("cmdstanr", quietly=TRUE))
  try(cmdstanr::set_cmdstan_path("~/.cmdstan/cmdstan-2.38.0"), silent=TRUE)

# Notify — do NOT auto-source 00_load.R here
if (interactive())
  cat("\nHarbor Seal IPM v3.2 — run source(\"00_load.R\") to load packages\n\n")
