#!/usr/bin/env Rscript
#' run_all.R
#' ----------------------------------------------------------------------
#' Regenerates, from a fresh clone, the data-driven matrices, the MDP and
#' POMDP solutions, all tables, the revision analyses and all figures of
#'
#'   Najafi et al., "Optimal Sequential Decision-Making for Error
#'   Propagation Mitigation in Modular Digital Twins".
#'
#' Usage (from the repository root, or from anywhere with the path):
#'   Rscript run_all.R                 # full mode (paper settings)
#'   Rscript run_all.R --quick         # quick mode (reduced Monte Carlo sizes)
#'   Rscript run_all.R --only=pipeline,model,table5
#'   Rscript run_all.R --skip=rl_sweep,tep
#'   Rscript run_all.R --notebooks     # additionally knit notebooks/*.Rmd
#'   Rscript run_all.R --download-tep  # fetch the Tennessee Eastman data first
#'   Rscript run_all.R --list          # list the steps and exit
#' Quick mode can also be selected with the environment variable DTMDP_QUICK=1.
#'
#' Every step runs in a fresh R process (Rscript) with DTMDP_ROOT set
#' to the repository root (no setwd needed); its console output goes to
#' output/logs/<step>.log. All results are written below output/.
#' The first step regenerates output/data_driven/data_driven_matrices.rds;
#' later steps reuse it (DTMDP_REUSE_MATRICES=1) instead of re-fitting the
#' HMMs, which gives identical matrices.
#' ----------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
root <- if (length(file_arg)) dirname(normalizePath(file_arg[1], winslash = "/")) else getwd()
root <- normalizePath(root, winslash = "/")
stopifnot(file.exists(file.path(root, "R", "paths.R")))

targs <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% targs ||
  tolower(Sys.getenv("DTMDP_QUICK", "0")) %in% c("1", "true", "yes")
if ("--full" %in% targs) quick <- FALSE
get_list <- function(flag) {
  a <- targs[startsWith(targs, paste0(flag, "="))]
  if (!length(a)) return(character(0))
  strsplit(sub(paste0("^", flag, "="), "", a[1]), ",")[[1]]
}
only <- get_list("--only"); skip <- get_list("--skip")

Sys.setenv(DTMDP_ROOT = root, DTMDP_QUICK = if (quick) "1" else "0")
source(file.path(root, "R", "paths.R"))
log_dir <- dirname(out_path("logs", "x"))

# ---------------------------------------------------------------------
# Steps: id, script, description
# ---------------------------------------------------------------------
steps <- list(
  list("pipeline",    "R/SimulationPipeline.R",
       "physics -> ARX -> HMM: data-driven A_NoAction and O (mcda_confusion)"),
  list("model",       "analyses/model_summary.R",
       "MDP (value iteration) and POMDP (PBVI) solutions"),
  list("arx",         "analyses/arx_validation.R",
       "ARX surrogate held-out fit (Table: ARX surrogates)"),
  list("mdp_valid",   "analyses/mdp_validation.R",
       "Bellman vs Gillespie Monte Carlo validation"),
  list("pbvi_conv",   "analyses/Revision_PBVI_Convergence.R",
       "PBVI convergence trace (Supplementary)"),
  list("kstep",       "analyses/Revision_KStepHeuristicBaseline.R",
       "k-step heuristic sweep, k in {1,2,3,5,10,25}"),
  list("table5",      "analyses/table5_final.R",
       "policy summary table + Wilcoxon / Cliff's delta tests"),
  list("state_diag",  "analyses/state_diagnostics.R",
       "occupancy over time and POMDP mismatch decomposition (Section 5.4, Supplementary S5.2)"),
  list("comp_cost",   "analyses/Revision_ComputationalCost.R",
       "computational cost table (Supplementary)"),
  list("reward_sens", "analyses/Revision_RewardSensitivity.R",
       "reward-perturbation sensitivity (Supplementary)"),
  list("tep",         "tep/TEP_Application.R",
       "Tennessee Eastman transfer (skipped if tep/TEP_data/ is absent)"),
  list("revision",    "analyses/revision_analyses.R",
       "CT approximation, PBVI init, myopic/QMDP, causal O, misspecification, slow repair"),
  list("rl_sweep",    c("analyses/rl_sweep/00_setup_model.R",
                        "analyses/rl_sweep/02_run_sweep.R",
                        "analyses/rl_sweep/04_summarize.R"),
       "RL tuning / budget / seed sweep (Supplementary; parallel)"),
  list("fig01",       "figures/fig01_trajectory.R",        "fig01_trajectory.pdf"),
  list("fig03",       "figures/fig03_reward_comparison.R",
       "fig03_reward_comparison.pdf (from the table5 trajectories)"),
  list("fig04",       "figures/fig04_sensitivity_analysis.R",
       "fig04_sensitivity_analysis.pdf (sweeps + table5 means)"),
  list("fig06",       "figures/fig06_action_mismatch.R",
       "fig06_action_mismatch.pdf (from the table5 trajectories)"),
  list("figs_mcda",   "figures/regen_figures.R",
       "fig02_state_occupancy.pdf, fig07_discount_factor.pdf"),
  list("fig05",       "figures/fig05_belief_trajectory.R", "fig05_belief_trajectory.pdf"),
  list("replot",      "figures/Revision_Replot.R",
       "PBVI_Convergence.pdf, RewardSensitivity_Combined.pdf")
)
if (dir.exists(file.path(root, "analyses", "closed_loop")) &&
    file.exists(file.path(root, "analyses", "closed_loop", "run_closed_loop.R"))) {
  steps <- append(steps, list(list("closed_loop", "analyses/closed_loop/run_closed_loop.R",
                                   "closed-loop digital-twin check")), after = 11)
}
ids <- vapply(steps, `[[`, "", 1)

if ("--list" %in% targs) {
  for (s in steps) cat(sprintf("  %-12s %s\n               %s\n", s[[1]], paste(s[[2]], collapse = " + "), s[[3]]))
  quit(save = "no")
}
bad <- setdiff(c(only, skip), ids)
if (length(bad)) stop("Unknown step(s): ", paste(bad, collapse = ", "),
                      ". Use --list to see the step ids.")
run_ids <- if (length(only)) intersect(ids, only) else ids
run_ids <- setdiff(run_ids, skip)

# ---------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
cat(sprintf("Repository root: %s\nMode: %s\nSteps: %s\n\n", root,
            if (quick) "QUICK (reduced Monte Carlo sizes)" else "FULL (paper settings)",
            paste(run_ids, collapse = ", ")))

# Reuse cached matrices only if the pipeline step is not part of this run
# (or after it has run successfully).
if (!"pipeline" %in% run_ids &&
    file.exists(out_path("data_driven", "data_driven_matrices.rds")))
  Sys.setenv(DTMDP_REUSE_MATRICES = "1")

summary_rows <- list()
for (s in steps[match(run_ids, ids)]) {
  id <- s[[1]]; scripts <- s[[2]]
  logf <- file.path(log_dir, paste0(id, ".log"))
  if (file.exists(logf)) file.remove(logf)
  cat(sprintf("[%s] %-12s %s ... ", format(Sys.time(), "%H:%M:%S"), id, s[[3]]))
  t0 <- Sys.time(); status <- 0L
  if (id == "tep" && "--download-tep" %in% targs)
    scripts <- c("tep/download_tep_data.R", scripts)
  for (sc in scripts) {
    status <- system2(rscript, c("--no-save", "--no-restore", shQuote(file.path(root, sc))),
                      stdout = logf, stderr = logf, wait = TRUE)
    # system2 overwrites the log on each call; keep the logs of multi-script steps
    if (length(scripts) > 1) file.copy(logf, file.path(log_dir, paste0(id, "_", basename(sc), ".log")),
                                       overwrite = TRUE)
    if (!identical(as.integer(status), 0L)) break
  }
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  ok <- identical(as.integer(status), 0L)
  cat(sprintf("%s (%.1f min)\n", if (ok) "ok" else sprintf("FAILED (status %s, see %s)", status, logf), el))
  summary_rows[[id]] <- data.frame(step = id, script = paste(scripts, collapse = " + "),
                                   status = if (ok) "ok" else "failed", minutes = round(el, 2))
  if (id == "pipeline") {
    if (!ok) stop("The data-driven pipeline failed; later steps depend on it. See ", logf)
    Sys.setenv(DTMDP_REUSE_MATRICES = "1")
  }
}

# Optional: knit the notebooks (HTML into output/notebooks/)
if ("--notebooks" %in% targs) {
  nb <- c(list.files(file.path(root, "notebooks"), "\\.Rmd$", full.names = TRUE),
          if (file.exists(out_path("tep", "TEP_results.rds")))
            file.path(root, "tep", "TEP_Application_Notebook.Rmd"))
  for (f in nb) {
    cat(sprintf("[%s] knit %s ... ", format(Sys.time(), "%H:%M:%S"), basename(f)))
    logf <- file.path(log_dir, paste0("notebook_", sub("\\.Rmd$", "", basename(f)), ".log"))
    expr <- sprintf("rmarkdown::render('%s', output_dir = '%s', intermediates_dir = tempdir(), quiet = TRUE)",
                    f, dirname(out_path("notebooks", "x")))
    st <- system2(rscript, c("-e", shQuote(expr)), stdout = logf, stderr = logf)
    cat(if (identical(as.integer(st), 0L)) "ok\n" else sprintf("FAILED (see %s)\n", logf))
    summary_rows[[basename(f)]] <- data.frame(step = basename(f), script = f,
      status = if (identical(as.integer(st), 0L)) "ok" else "failed", minutes = NA)
  }
}

summ <- do.call(rbind, summary_rows)
write.csv(summ, file.path(log_dir, "run_all_summary.csv"), row.names = FALSE)
writeLines(c(sprintf("mode: %s", if (quick) "quick" else "full"),
             capture.output(sessionInfo())), file.path(log_dir, "sessionInfo.txt"))

dd <- out_path("data_driven", "data_driven_matrices.rds")
if (file.exists(dd)) {
  O <- readRDS(dd)$mcda_confusion
  cat(sprintf("\nData-driven O: diag = %s (Nominal accuracy %.3f)\n",
              paste(sprintf("%.3f", diag(O)), collapse = ", "), O["Nominal", "Nominal"]))
}
cat(sprintf("\n%d/%d steps succeeded. Results in %s\n", sum(summ$status == "ok"),
            nrow(summ), dirname(log_dir)))
if (any(summ$status != "ok")) quit(save = "no", status = 1)
