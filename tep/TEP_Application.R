#' TEP_Application.R
#' ──────────────────────────────────────────────────────────────────────
#' End-to-end driver: applies the ARX → HMM → MDP/POMDP pipeline to the
#' Tennessee Eastman benchmark (paper Section "Transfer to the Tennessee
#' Eastman Process").
#'
#'   1. Source tep/TEP_Pipeline.R → A_NoAction, mcda_confusion, mcda_accuracy
#'                                  derived from real TEP fault data
#'   2. Source R/POMDPFormulation.R with options(dtmdp.use_existing_matrices
#'      = TRUE, dtmdp.dt = 180), so R/BaseMDPFormulation.R uses the TEP
#'      parameters instead of re-running the synthetic pipeline
#'   3. Solve MDP via value iteration and POMDP via PBVI
#'   4. Compare TEP-derived versus synthetic policy + value
#'
#' Output: output/tep/TEP_results.rds (used by analyses/revision_analyses.R,
#' section 7, and tep/TEP_Application_Notebook.Rmd).
#'
#' Requires the TEP data files in tep/TEP_data/ (see tep/download_tep_data.R);
#' if they are absent the script prints a message and exits without error.
#' ──────────────────────────────────────────────────────────────────────

local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })

tep_dir <- Sys.getenv("DTMDP_TEP_DATA", repo_path("tep", "TEP_data"))
tep_files_needed <- c("d00.dat", "d00_te.dat", "d01_te.dat", "d11_te.dat", "d13_te.dat")
if (!all(file.exists(file.path(tep_dir, tep_files_needed)))) {
  message("Tennessee Eastman data not found in ", tep_dir, ".\n",
          "Run `Rscript tep/download_tep_data.R` (or set DTMDP_TEP_DATA) and re-run. ",
          "Skipping the TEP application.")
  if (!interactive()) quit(save = "no", status = 0)
  stop("TEP data missing", call. = FALSE)
}

cat("\n══════════════════════════════════════════════════════════\n")
cat(" TEP APPLICATION OF ARX → HMM → MDP/POMDP PIPELINE\n")
cat("══════════════════════════════════════════════════════════\n\n")

# Synthetic-testbed values for the comparison table (diag of O and A_NoAction
# from R/SimulationPipeline.R; read from output/ if the pipeline has been run)
dd_file <- out_path("data_driven", "data_driven_matrices.rds")
if (file.exists(dd_file)) {
  dd <- readRDS(dd_file)
  acc_synth  <- dd$mcda_accuracy
  diag_synth <- diag(dd$A_NoAction)
} else {
  acc_synth  <- c(Nominal = 0.942, SensorNoisy = 0.210, DynamicsOff = 0.799, Drift = 0.996)
  diag_synth <- c(Nominal = 0.997, SensorNoisy = 0.985, DynamicsOff = 0.907, Drift = 0.856)
}

# -----------------------------------------------------------------------
# Step 1: TEP-derived parameters (A_NoAction, mcda_confusion, mcda_accuracy)
# -----------------------------------------------------------------------
source(repo_path("tep", "TEP_Pipeline.R"))

A_NoAction_TEP     <- A_NoAction
mcda_confusion_TEP <- mcda_confusion
mcda_accuracy_TEP  <- mcda_accuracy

# -----------------------------------------------------------------------
# Step 2-3: MDP + POMDP on the TEP parameters.
#   dtmdp.use_existing_matrices: R/BaseMDPFormulation.R does not source the
#     synthetic R/SimulationPipeline.R and uses A_NoAction / mcda_confusion
#     from the global environment.
#   dtmdp.dt = 180: TEP samples every 3 min (180 s).
#   V_ct falls back to the discrete V* if (r I - Q_pi) is singular
#     (repair-modified rows have near-zero off-diagonals at this dt).
# -----------------------------------------------------------------------
options(dtmdp.use_existing_matrices = TRUE, dtmdp.dt = 180)
cat("\n── Sourcing R/POMDPFormulation.R with TEP parameters ──\n")
source(repo_path("R", "POMDPFormulation.R"))
options(dtmdp.use_existing_matrices = NULL, dtmdp.dt = NULL)

cat("\n── MDP results on TEP-derived parameters ──\n")
cat("Optimal policy (state → action):\n")
print(vi$policy)
cat("\nOptimal MDP values V*(s):\n")
print(round(vi$V, 4))
cat("\nContinuous-time Bellman values V_ct(s):\n")
print(round(V_ct, 4))

# -----------------------------------------------------------------------
# Step 4: Final summary table comparing TEP vs synthetic
# -----------------------------------------------------------------------
cat("\n══════════════════════════════════════════════════════════\n")
cat(" FINAL SUMMARY: TEP vs SYNTHETIC PIPELINE\n")
cat("══════════════════════════════════════════════════════════\n")

cat("\nState-by-state classification accuracy (HMM Viterbi vs ground truth):\n")
acc_tbl <- rbind(synthetic = acc_synth, TEP = mcda_accuracy_TEP)
print(round(acc_tbl, 4))

cat("\nA_NoAction diagonal (regime persistence per dt):\n")
diag_tbl <- rbind(synthetic = diag_synth, TEP = diag(A_NoAction_TEP))
print(round(diag_tbl, 4))

cat("\nOptimal MDP policy on TEP-derived dynamics:\n")
print(vi$policy)

cat("\nMDP value V*(s):\n")
print(round(vi$V, 3))

# Persist key results to disk for reuse
saveRDS(list(
  A_NoAction     = A_NoAction_TEP,
  mcda_confusion = mcda_confusion_TEP,
  mcda_accuracy  = mcda_accuracy_TEP,
  policy         = vi$policy,
  V_star         = vi$V,
  V_ct           = V_ct
), file = out_path("tep", "TEP_results.rds"))
cat("\nSaved output/tep/TEP_results.rds\n")

cat("\n══════════════════════════════════════════════════════════\n")
cat(" TEP application complete.\n")
cat("══════════════════════════════════════════════════════════\n")
