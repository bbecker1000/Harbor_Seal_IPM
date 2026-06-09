# ============================================================================
# HARBOR SEAL IPM v3.2 — PROJECT SETUP
# Run once per machine. Works with or without renv already active.
#
# Usage:  source("harbor_seal_setup.R")
# ============================================================================

cat("\n══════════════════════════════════════════════════════════════\n")
cat("  HARBOR SEAL IPM v3.2 — PROJECT SETUP\n")
cat("══════════════════════════════════════════════════════════════\n\n")
cat("R version:    ", R.version.string, "\n")
cat("Working dir:  ", getwd(), "\n")
cat("Library path: ", .libPaths()[1], "\n\n")

# ── Detect environment ────────────────────────────────────────────────────────
renv_active <- any(grepl("renv", .libPaths()))
is_linux    <- .Platform$OS.type == "unix" && !grepl("darwin", R.version$os)
cat(sprintf("renv active:  %s\n", renv_active))
cat(sprintf("Platform:     %s\n\n", if (is_linux) "Linux" else R.version$os))

#--
# After Step 1 confirms renv is installed, add this:
if (file.exists("renv.lock")) {
  cat("  renv.lock found — restoring all packages from lockfile...\n")
  cat("  (This reinstalls exact versions used previously)\n")
  renv::restore(prompt = FALSE)
  cat("  Restore complete\n")
} else {
  cat("  No renv.lock found — will install fresh and snapshot at end\n")
}

# ── Repository ────────────────────────────────────────────────────────────────
if (is_linux) {
  repos <- c(
    RSPM = "https://packagemanager.posit.co/cran/__linux__/noble/latest",
    CRAN = "https://cloud.r-project.org"
  )
  cat("Using Posit Package Manager (pre-compiled Ubuntu Noble binaries)\n\n")
} else {
  repos <- c(CRAN = "https://cloud.r-project.org")
}
options(repos = repos)

# ── Detect allocated cores ────────────────────────────────────────────────────
get_safe_cores <- function() {
  tryCatch({
    if (file.exists("/sys/fs/cgroup/cpu.max")) {
      p <- strsplit(readLines("/sys/fs/cgroup/cpu.max", warn=FALSE)[1], " ")[[1]]
      if (p[1] != "max") return(max(1L, as.integer(as.numeric(p[1])/as.numeric(p[2]))))
    } else if (file.exists("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")) {
      q  <- as.numeric(readLines("/sys/fs/cgroup/cpu/cpu.cfs_quota_us",  warn=FALSE)[1])
      pe <- as.numeric(readLines("/sys/fs/cgroup/cpu/cpu.cfs_period_us", warn=FALSE)[1])
      if (q > 0) return(max(1L, as.integer(q/pe)))
    }
  }, error = function(e) NULL)
  env <- Sys.getenv("NSLOTS")
  if (nzchar(env)) return(max(1L, as.integer(env)))
  max(1L, parallel::detectCores() - 1L)
}
n_cores <- get_safe_cores()
cat(sprintf("Allocated cores: %d\n\n", n_cores))

# ── Install helper ────────────────────────────────────────────────────────────
pkg_ok <- function(p) requireNamespace(p, quietly = TRUE)

install_pkgs <- function(pkgs, extra_repos = NULL) {
  missing <- pkgs[!sapply(pkgs, pkg_ok)]
  if (length(missing) == 0) return(invisible(NULL))
  all_repos <- if (is.null(extra_repos)) repos else c(extra_repos, repos)
  cat(sprintf("  Installing: %s\n", paste(missing, collapse=", ")))
  if (renv_active) {
    tryCatch(
      renv::install(missing, repos=all_repos, prompt=FALSE),
      error = function(e) install.packages(missing, repos=all_repos,
                                           type="binary", dependencies=TRUE, Ncpus=n_cores)
    )
  } else {
    install.packages(missing, repos=all_repos, dependencies=TRUE, Ncpus=n_cores)
  }
  still <- missing[!sapply(missing, pkg_ok)]
  if (length(still) > 0) cat("  FAILED:", paste(still, collapse=", "), "\n")
  else                    cat("  OK\n")
}

# ── Step 0: Git configuration ─────────────────────────────────────────────────
cat("── Step 0/5: Git configuration ──────────────────────────────\n")

system('git config --global credential.helper store')
system('git config --global user.email "your@email.com"')   # ← edit
system('git config --global user.name  "Your Name"')        # ← edit

# Confirm
cat(system('git config --global --list', intern=TRUE), sep="\n")
cat("  Git configured\n")
cat("  NOTE: run 'git push' once manually to cache your token\n\n")


# ── 1. renv ───────────────────────────────────────────────────────────────────
cat("── Step 1/5: renv ───────────────────────────────────────────\n")
if (!pkg_ok("renv")) install.packages("renv", repos=repos, quiet=TRUE)
cat(sprintf("  renv v%s\n", as.character(packageVersion("renv"))))

# ── 2. CRAN packages ──────────────────────────────────────────────────────────
cat("\n── Step 2/5: CRAN packages ──────────────────────────────────\n")

all_pkgs <- c(
  # ── Stan modelling ───────────────────────────────────────────────────────────
  "posterior",       # draws_df / extract_variable()
  "bayesplot",       # MCMC diagnostics, trace plots, PPCs
  # ── Tidyverse suite (ggplot2, dplyr, tidyr, purrr, readr, tibble, stringr…) ─
  "tidyverse",
  # ── Date handling (also in tidyverse but loaded explicitly in data prep) ─────
  "lubridate",       # ymd(), year(), yday() — Phoca data prep
  # ── Data manipulation (load BEFORE dplyr to avoid masking conflicts) ─────────
  "plyr",            # ddply() — Phoca data prep; must precede dplyr
  # ── Plot layout + colour ─────────────────────────────────────────────────────
  "patchwork",       # wrap_plots(), plot_annotation()
  "scales",          # scales::squish, label helpers
  "viridis",         # scale_fill_viridis_c/d
  "RColorBrewer",    # scale_color_brewer Dark2 / Set1
  # ── Tables ───────────────────────────────────────────────────────────────────
  "knitr",           # kable()
  "kableExtra",      # formatted HTML/LaTeX tables
  # ── File I/O ─────────────────────────────────────────────────────────────────
  "here",            # here::here() project-relative paths
  "readxl",          # read_excel() for count data spreadsheets
  # ── Model utilities ──────────────────────────────────────────────────────────
  "loo",             # LOO-CV model comparison
  # ── GAM (phenology paper) ────────────────────────────────────────────────────
  "mgcv",            # gam(), s() — hierarchical GAM for phenology timing
  # ── General utilities ────────────────────────────────────────────────────────
  "abind",           # array binding for posterior arrays
  "progress"         # progress bars in projection loops
)

already    <- all_pkgs[sapply(all_pkgs, pkg_ok)]
to_install <- all_pkgs[!sapply(all_pkgs, pkg_ok)]

if (length(already)    > 0) cat("  Already installed:", paste(already, collapse=", "), "\n")
if (length(to_install) > 0) install_pkgs(to_install)

# ── 3. cmdstanr ───────────────────────────────────────────────────────────────
cat("\n── Step 3/5: cmdstanr ───────────────────────────────────────\n")
cmdstanr_repos <- c("https://mc-stan.org/r-packages/", repos)
if (!pkg_ok("cmdstanr")) {
  cat("  Installing cmdstanr...\n")
  if (renv_active) {
    tryCatch(
      renv::install("cmdstanr", repos=cmdstanr_repos, prompt=FALSE),
      error = function(e) install.packages("cmdstanr", repos=cmdstanr_repos, dependencies=TRUE)
    )
  } else {
    install.packages("cmdstanr", repos=cmdstanr_repos, dependencies=TRUE)
  }
}
if (!pkg_ok("cmdstanr")) stop("cmdstanr install failed — try renv::install('cmdstanr')")
suppressPackageStartupMessages(library(cmdstanr))

# ── 4. CmdStan ────────────────────────────────────────────────────────────────
cat("\n── Step 4/5: CmdStan ────────────────────────────────────────\n")
cs_ok <- tryCatch({
  cat(sprintf("  Already installed: v%s\n  Path: %s\n",
              cmdstan_version(), cmdstan_path()))
  TRUE
}, error = function(e) FALSE)

if (!cs_ok) {
  cat(sprintf("  Compiling on %d core(s) — 5-15 min...\n", n_cores))
  install_cmdstan(cores=n_cores)
  cat(sprintf("  Installed: v%s\n", cmdstan_version()))
}
tryCatch(check_cmdstan_toolchain(fix=TRUE, quiet=TRUE),
         error = function(e) message("  WARNING: toolchain — ", conditionMessage(e)))
cmdstan_installed_path <- tryCatch(cmdstan_path(), error=function(e) NULL)

# ── 5. Write .Rprofile + 00_load.R ───────────────────────────────────────────
cat("\n── Step 5/5: Writing .Rprofile and 00_load.R ────────────────\n")

# Packages loaded at every session start
# plyr MUST come before dplyr (both loaded by tidyverse but plyr first avoids
# masking warnings on summarise/mutate)
AUTOLOAD_PKGS <- c(
  "plyr",          # before tidyverse/dplyr to avoid masking
  "cmdstanr",      # Stan interface; sets CmdStan PATH
  "posterior",     # draws_df, extract_variable()
  "bayesplot",     # trace plots, PPCs
  "tidyverse",     # ggplot2, dplyr, tidyr, purrr, readr, tibble, stringr
  "lubridate",     # date helpers
  "patchwork",     # multi-panel figures
  "scales",        # scales::squish, axis helpers
  "viridis",       # colour palettes
  "RColorBrewer",  # brewer palettes
  "here",          # project-relative paths
  "readxl",        # read_excel()
  "loo"            # model comparison
)
pkg_str <- paste0('"', AUTOLOAD_PKGS, '"', collapse=", ")

# ── .Rprofile ─────────────────────────────────────────────────────────────────
rprofile_path <- file.path(getwd(), ".Rprofile")
existing      <- if (file.exists(rprofile_path)) readLines(rprofile_path, warn=FALSE) else character(0)
renv_line <- if (any(grepl("renv/activate", existing))) {
  existing[grep("renv/activate", existing)[1]]
} else {
  'if (file.exists("renv/activate.R")) source("renv/activate.R")'
}
cmdstan_line <- if (!is.null(cmdstan_installed_path)) {
  paste0('if (requireNamespace("cmdstanr",quietly=TRUE))\n',
         '  tryCatch(cmdstanr::set_cmdstan_path("',
         cmdstan_installed_path, '"), error=function(e) NULL)')
} else {
  '# CmdStan path unknown at setup time'
}

rprofile_text <- paste0(
  '# ============================================================================
# HARBOR SEAL IPM v3.2 — .Rprofile  (generated by harbor_seal_setup.R)
# Auto-sourced by RStudio at session start for this project.
# ============================================================================

# ── renv: activate locked package library ────────────────────────────────────
', renv_line, '

# ── Repository ────────────────────────────────────────────────────────────────
local({
  linux <- .Platform$OS.type == "unix" && !grepl("darwin", R.version$os)
  options(repos = if (linux)
    c(RSPM = "https://packagemanager.posit.co/cran/__linux__/noble/latest",
      CRAN = "https://cloud.r-project.org")
  else
    c(CRAN = "https://cloud.r-project.org")
  )
})

# ── Session options ───────────────────────────────────────────────────────────
options(warn=1, scipen=999)
options(mc.cores = ', n_cores, 'L)   # parallel Stan chains

# ── CmdStan path ──────────────────────────────────────────────────────────────
', cmdstan_line, '

# ── Autoload packages ─────────────────────────────────────────────────────────
# Sourcing 00_load.R handles both .Rprofile-triggered and manual starts.
if (file.exists("00_load.R")) {
  source("00_load.R")
} else {
  # Inline fallback if 00_load.R is missing
  .hs_pkgs <- c(', pkg_str, ')
  for (.p in .hs_pkgs)
    suppressPackageStartupMessages(
      tryCatch(library(.p, character.only=TRUE, quietly=TRUE),
               error=function(e) NULL))
  rm(.hs_pkgs, .p)
}
')

writeLines(rprofile_text, rprofile_path)
cat("  .Rprofile written:", rprofile_path, "\n")

# ── 00_load.R — explicit session starter (more reliable than .Rprofile) ───────
# Can be sourced manually: source("00_load.R")
# Also sourced automatically from .Rprofile above.
load_script <- paste0(
  '# ============================================================================
# HARBOR SEAL IPM v3.2 — 00_load.R
# Loads all project packages and prints a startup banner.
# Called automatically from .Rprofile, OR run manually:
#   source("00_load.R")
# ============================================================================

.pkgs <- c(', pkg_str, ')

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
cat("\\n══════════════════════════════════════════════════════\\n")
cat("  Harbor Seal IPM v3.2\\n")
if (length(.failed) > 0) {
  cat("  Missing packages — run source(\\"harbor_seal_setup.R\\"):\\n")
  cat("   ", paste(.failed, collapse=", "), "\\n")
} else {
  .cs <- tryCatch(cmdstanr::cmdstan_version(), error=function(e) "not found")
  cat(sprintf("  CmdStan  : v%s\\n", .cs))
  cat(sprintf("  tidyverse: v%s  |  posterior: v%s\\n",
              as.character(packageVersion("tidyverse")),
              as.character(packageVersion("posterior"))))
  cat(sprintf("  mc.cores : %d\\n", getOption("mc.cores", 1L)))
  cat(sprintf("  Loaded   : %s\\n", paste(.pkgs[!.pkgs %%in%% .failed], collapse=", ")))
}
cat("\\n  Quick start:\\n")
cat("    source(\\"Code/harbor_seal_ipm_v3.2.R\\")\\n")
cat("    out <- load_seal_results(\\"IPM_v3.2_real\\")\\n")
cat("══════════════════════════════════════════════════════\\n\\n")
rm(.pkgs, .failed, .ok, .p, .cs, .d)
')

# Fix the %%in%% escaping (Python heredoc artifact)
load_script <- gsub("%%in%%", "%in%", load_script, fixed=TRUE)
writeLines(load_script, "00_load.R")
cat("  00_load.R written:", file.path(getwd(), "00_load.R"), "\n")

# ── Lock versions ─────────────────────────────────────────────────────────────
cat("\n  Locking package versions...\n")
tryCatch({
  if (file.exists("renv.lock")) renv::snapshot(type="implicit", prompt=FALSE)
  else                          renv::snapshot(prompt=FALSE)
  cat("  renv.lock updated\n")
}, error = function(e) cat("  renv snapshot skipped:", conditionMessage(e), "\n"))


# Ensure dplyr verbs are not masked by stats or plyr
filter  <- dplyr::filter
select  <- dplyr::select
mutate  <- dplyr::mutate
arrange <- dplyr::arrange



# ── Done ──────────────────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════\n")
cat("  SETUP COMPLETE\n\n")
cat("  → Restart R (Session ▸ Restart R)\n")
cat("  → Packages autoload via .Rprofile → 00_load.R\n")
cat("  → If autoload doesn't fire: source('00_load.R')\n\n")
cat("  Packages installed:\n")
cat("  ", paste(c("cmdstanr", all_pkgs), collapse=", "), "\n\n")
cat("  On a new machine: renv::restore() then source('harbor_seal_setup.R')\n")
cat("══════════════════════════════════════════════════════════════\n\n")