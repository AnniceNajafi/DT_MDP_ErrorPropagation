#' dependencies.R
#' ----------------------------------------------------------------------
#' Installs the CRAN packages needed to reproduce the paper.
#'
#'   Rscript dependencies.R            # install missing packages (latest CRAN)
#'   Rscript dependencies.R --exact    # install the exact versions used for
#'                                     # the paper (needs the 'remotes' package)
#'   Rscript dependencies.R --check    # only report installed vs. reference versions
#'
#' Reference environment: R 4.5.1 (2025-06-13 ucrt), Windows 11 x64.
#' renv is optional: `renv::init(bare = TRUE); source("dependencies.R");
#' renv::snapshot()` creates a lockfile if you want one.
#' ----------------------------------------------------------------------

# Versions used to produce the results in the paper (packageVersion())
required <- c(
  depmixS4       = "1.5.1",   # HMM fitting (EM, Viterbi, forward-backward)
  clue           = "0.3.66",  # Hungarian assignment (solve_LSAP)
  dplyr          = "1.1.4",
  tibble         = "3.3.0",
  expm           = "1.0.0",   # matrix exponential / logarithm
  Matrix         = "1.7.3",   # (dependency of expm)
  RMCDA          = "0.3",     # TOPSIS baseline (apply.TOPSIS)
  ggplot2        = "4.0.3",
  reshape2       = "1.4.4",
  gridExtra      = "2.3",
  patchwork      = "1.3.2",
  scales         = "1.4.0",
  microbenchmark = "1.5.0"    # computational-cost table
)
# Only needed to knit the notebooks in notebooks/ and tep/ (optional)
optional <- c(
  knitr      = "1.50",
  rmarkdown  = "2.30",        # plus pandoc (tested with 3.8.3)
  tidyr      = "1.3.1",
  kableExtra = "1.4.0"
)

args  <- commandArgs(trailingOnly = TRUE)
exact <- "--exact" %in% args
check <- "--check" %in% args
pkgs  <- c(required, optional)

installed_version <- function(p)
  if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_

if (check) {
  cur <- vapply(names(pkgs), installed_version, character(1))
  print(data.frame(package = names(pkgs), reference = unname(pkgs),
                   installed = unname(cur),
                   required = names(pkgs) %in% names(required)), row.names = FALSE)
  quit(save = "no", status = if (any(is.na(cur[names(required)]))) 1 else 0)
}

repos <- getOption("repos")
if (is.null(repos) || identical(unname(repos["CRAN"]), "@CRAN@"))
  options(repos = c(CRAN = "https://cloud.r-project.org"))

if (exact) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  for (p in names(pkgs)) {
    if (!identical(installed_version(p), pkgs[[p]]))
      remotes::install_version(p, version = pkgs[[p]], upgrade = "never")
  }
} else {
  missing <- names(pkgs)[is.na(vapply(names(pkgs), installed_version, character(1)))]
  if (length(missing)) install.packages(missing) else cat("All packages already installed.\n")
}
