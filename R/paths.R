#' paths.R
#' ----------------------------------------------------------------------
#' Path and run-mode helpers shared by every script in the repository.
#'
#' All paths are resolved relative to the repository root. The root is
#' taken from the environment variable DTMDP_ROOT (set by run_all.R);
#' if it is not set, the current working directory and its parents are
#' searched for the directory that contains run_all.R and R/paths.R.
#'
#' Every entry-point script starts with the same two-line bootstrap:
#'
#'   local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
#'     while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
#'     source(file.path(r, "R", "paths.R")) })
#'
#' Environment variables understood by the scripts:
#'   DTMDP_ROOT             repository root (absolute path)
#'   DTMDP_QUICK            "1" = quick mode (reduced Monte Carlo sizes)
#'   DTMDP_OUT              output directory (default <root>/output)
#'   DTMDP_REUSE_MATRICES   "1" = reuse output/data_driven/data_driven_matrices.rds
#'                          instead of re-fitting the HMMs (run_all.R sets this
#'                          after the pipeline step has regenerated the file)
#'   DTMDP_TEP_DATA         directory with the Tennessee Eastman .dat files
#'                          (default <root>/tep/TEP_data)
#' ----------------------------------------------------------------------

local({
  root <- Sys.getenv("DTMDP_ROOT")
  if (!nzchar(root)) {
    root <- getwd()
    while (!(file.exists(file.path(root, "R", "paths.R")) &&
             file.exists(file.path(root, "run_all.R"))) &&
           dirname(root) != root) {
      root <- dirname(root)
    }
    if (!file.exists(file.path(root, "R", "paths.R")))
      stop("Cannot locate the repository root. Run from inside the repository ",
           "or set the DTMDP_ROOT environment variable.")
  }
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  Sys.setenv(DTMDP_ROOT = root)
  assign("DTMDP_ROOT", root, envir = globalenv())
})

#' Path inside the repository.
repo_path <- function(...) file.path(DTMDP_ROOT, ...)

#' Path inside the output directory; parent directories are created.
out_path <- function(...) {
  base <- Sys.getenv("DTMDP_OUT")
  if (!nzchar(base)) base <- file.path(DTMDP_ROOT, "output")
  p <- file.path(base, ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}

#' TRUE in quick mode (DTMDP_QUICK=1).
QUICK <- tolower(Sys.getenv("DTMDP_QUICK", "0")) %in% c("1", "true", "yes")

#' Monte Carlo size: `full` in full mode, `quick` in quick mode.
nq <- function(full, quick = max(20L, as.integer(round(full / 10)))) {
  if (QUICK) quick else full
}

#' Source the full model (pipeline -> MDP -> POMDP) once per R session.
load_core <- function() {
  if (!isTRUE(get0(".dtmdp_core_loaded", envir = globalenv(), inherits = FALSE)))
    source(repo_path("R", "POMDPFormulation.R"), local = globalenv())
  invisible(TRUE)
}

if (QUICK && !isTRUE(getOption("dtmdp.quick_msg"))) {
  cat("[quick mode: reduced Monte Carlo sizes; numbers will differ from the paper]\n")
  options(dtmdp.quick_msg = TRUE)
}
