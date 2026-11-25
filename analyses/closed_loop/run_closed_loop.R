#' run_closed_loop.R -- closed-loop validation on the six-module twin
#' (Reviewer 2, comments 6 and 8; see README.md in this folder).
#'
#' Runs, each in a fresh R process, with logs in output/logs/closed_loop_*.log:
#'   cl_setup.R        decision model + causal decoder HMM -> output/closed_loop/setup.rds
#'   verify.R          stepwise twin == SimulationPipeline; calibration diagnostics
#'   cl_repair.R       repair-time experiments (default, W_bc = 2 s, W_rp = 4 s)
#'   cl_closed_loop.R  closed-loop controller comparison (default + 4 variants)
#'
#' Usage (repository root):
#'   Rscript analyses/closed_loop/run_closed_loop.R           # 500 / 1000 seeds (paper)
#'   Rscript analyses/closed_loop/run_closed_loop.R --quick   # 20 seeds, smoke test
#' DTMDP_CORES sets the number of parallel workers.

local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
if ("--quick" %in% commandArgs(trailingOnly = TRUE)) QUICK <- TRUE

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
code    <- repo_path("analyses", "closed_loop")
n_rep   <- nq(500L, 20L)
n_cl    <- nq(1000L, 20L)

jobs <- list(
  list("setup",        "cl_setup.R",       character()),
  list("verify",       "verify.R",         character()),
  list("repair",       "cl_repair.R",      c(CL_N_REPAIR = n_rep)),
  list("repair_Wbc2s", "cl_repair.R",      c(CL_N_REPAIR = n_rep, CL_CFG = "W_bc=100", CL_TAG = "_Wbc2s")),
  list("repair_Wrp4s", "cl_repair.R",      c(CL_N_REPAIR = n_rep, CL_CFG = "W_rp=200", CL_TAG = "_Wrp4s")),
  list("cl",           "cl_closed_loop.R", c(CL_N = n_cl)),
  list("cl_aonly",     "cl_closed_loop.R", c(CL_N = n_cl, CL_VARIANT = "aonly")),
  list("cl_rawdecoder","cl_closed_loop.R", c(CL_N = n_cl, CL_DEC_EPS = "0", CL_TAG = "_rawdecoder")),
  list("cl_Wbc2s",     "cl_closed_loop.R", c(CL_N = n_cl, CL_CFG = "W_bc=100", CL_TAG = "_Wbc2s")),
  list("cl_Wrp4s",     "cl_closed_loop.R", c(CL_N = n_cl, CL_CFG = "W_rp=200", CL_TAG = "_Wrp4s"))
)

reset <- c("CL_N", "CL_N_REPAIR", "CL_CFG", "CL_TAG", "CL_VARIANT", "CL_DEC_EPS")
for (j in jobs) {
  cat(sprintf("[%s] closed_loop/%-14s ... ", format(Sys.time(), "%H:%M:%S"), j[[1]]))
  Sys.unsetenv(reset)
  if (length(j[[3]])) do.call(Sys.setenv, as.list(j[[3]]))
  log <- out_path("logs", paste0("closed_loop_", j[[1]], ".log"))
  t0 <- Sys.time()
  st <- system2(rscript, c("--no-save", "--no-restore", shQuote(file.path(code, j[[2]]))),
                stdout = log, stderr = log)
  cat(if (st == 0) "ok" else "FAILED", sprintf("(%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (st != 0) stop("closed-loop step ", j[[1]], " failed; see ", log)
}
Sys.unsetenv(reset)
cat("Closed-loop results in", out_path("closed_loop"), "\n")
