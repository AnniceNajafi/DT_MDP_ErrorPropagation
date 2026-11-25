#' TEP_Pipeline.R
#' ──────────────────────────────────────────────────────────────────────
#' Tennessee Eastman Process application of the ARX → HMM → MDP/POMDP
#' pipeline.  Drop-in replacement for SimulationPipeline.R: produces the
#' same downstream objects (A_NoAction, mcda_confusion, mcda_accuracy,
#' states_pipeline) so that BaseMDPFormulation.R / POMDPFormulation.R
#' can be sourced without modification.
#'
#' Data: Russell-Chiang-Braatz TEP benchmark, placed in tep/TEP_data/ by
#' tep/download_tep_data.R (or in the directory given by DTMDP_TEP_DATA).
#'   d00.dat      52 × 500  normal-operation training (transposed)
#'   d00_te.dat   960 × 52  normal-operation testing
#'   dXX_te.dat   960 × 52  fault-XX testing; fault begins at sample 161
#'
#' Regime mapping (see the TEP fault-detection literature):
#'   Nominal     ← d00_te                   pre-fault baseline
#'   SensorNoisy ← d11_te samples 161-960   IDV(11) random reactor cool-water T
#'   DynamicsOff ← d01_te samples 161-960   IDV(1)  step in A/C feed ratio
#'   Drift       ← d13_te samples 161-960   IDV(13) slow drift in reaction kinetics
#'
#' IDV(1), IDV(11), IDV(13) chosen because they (i) produce clearly distinct
#' residual signatures for HMM separability, (ii) physically correspond to
#' the SensorNoisy/DynamicsOff/Drift archetypes from the synthetic testbed.
#'
#' Surrogate model: per-output ARX driven by the 11 manipulated variables
#' XMV(1)..XMV(11), trained on d00.dat (normal operation).  Residuals
#' = actual XMEAS - ARX prediction.  Eight residual channels picked to
#' span the TEP unit-op cascade (reactor, separator, stripper, compressor).
#' ──────────────────────────────────────────────────────────────────────

local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })

suppressPackageStartupMessages({
  library(depmixS4)
  library(clue)
  library(dplyr)
  library(tibble)
})

cat("── TEP_Pipeline: extracting parameters from Tennessee Eastman data ──\n")

# Per-line file logger so we can localize a segfault that wipes stdout buffers
log_file <- out_path("logs", "TEP_progress.log")
unlink(log_file)
tep_log <- function(msg) {
  con <- file(log_file, open = "a")
  writeLines(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg), con)
  close(con)
}
tep_log("pipeline start")

# =============================================================================
# 1. PATHS & CONSTANTS
# =============================================================================

tep_dir <- Sys.getenv("DTMDP_TEP_DATA", repo_path("tep", "TEP_data"))
tep_files_needed <- c("d00.dat", "d00_te.dat", "d01_te.dat", "d11_te.dat", "d13_te.dat")
if (!all(file.exists(file.path(tep_dir, tep_files_needed))))
  stop("Tennessee Eastman data not found in ", tep_dir,
       " -- run `Rscript tep/download_tep_data.R` first.")

# 41 measured + 11 manipulated = 52 columns
n_xmeas <- 41L
n_xmv   <- 11L
n_vars  <- n_xmeas + n_xmv

# Eight measurement channels chosen to span the TEP unit-op cascade
# (reactor → separator → stripper → compressor → product)
sensor_idx <- c(
  reactor_pressure = 7L,
  reactor_level    = 8L,
  reactor_temp     = 9L,
  sep_temp         = 11L,
  sep_level        = 12L,
  stripper_temp    = 18L,
  comp_work        = 20L,
  reactor_cool_T   = 21L
)

# Canonical state ordering (must match MDP — same as synthetic pipeline)
states_pipeline <- c("Nominal", "SensorNoisy", "DynamicsOff", "Drift")
K_pipe <- length(states_pipeline)

# =============================================================================
# 2. FILE LOADERS
# =============================================================================

read_tep_standard <- function(path) {
  m <- as.matrix(read.table(path))
  stopifnot(ncol(m) == n_vars)
  unname(m)
}

read_tep_transposed <- function(path) {
  m <- as.matrix(read.table(path))
  stopifnot(nrow(m) == n_vars)
  unname(t(m))
}

cat("  Loading TEP files...\n")
D_norm_train <- read_tep_transposed(file.path(tep_dir, "d00.dat"))     # 500 × 52
D_norm_test  <- read_tep_standard  (file.path(tep_dir, "d00_te.dat"))  # 960 × 52
D_f01        <- read_tep_standard  (file.path(tep_dir, "d01_te.dat"))  # IDV(1)
D_f11        <- read_tep_standard  (file.path(tep_dir, "d11_te.dat"))  # IDV(11)
D_f13        <- read_tep_standard  (file.path(tep_dir, "d13_te.dat"))  # IDV(13)

cat(sprintf("    d00 train:   %d × %d\n", nrow(D_norm_train), ncol(D_norm_train)))
cat(sprintf("    d00 test:    %d × %d\n", nrow(D_norm_test),  ncol(D_norm_test)))
cat(sprintf("    d01/d11/d13: %d × %d each\n", nrow(D_f01), ncol(D_f01)))

# =============================================================================
# 3. ARX SURROGATE FITTING (driven by XMV inputs, trained on normal data)
# =============================================================================

# Variable layout in TEP files: XMEAS(1..41) then XMV(1..11)
xmeas_cols <- 1:n_xmeas
xmv_cols   <- (n_xmeas + 1):n_vars

fit_arx_miso <- function(y, U, na, nb) {
  # y : vector length n  (output)
  # U : matrix n × m     (m manipulated inputs)
  # na: AR order, nb: input order
  n <- length(y)
  m <- ncol(U)
  p <- max(na, nb)
  stopifnot(n > p)

  rows  <- vector("list", n - p)
  Tvec  <- numeric(n - p)
  for (k in (p + 1):n) {
    phi_y <- sapply(1:na, function(i) y[k - i])
    phi_u <- as.numeric(sapply(1:nb, function(j) U[k - j, ]))  # length m*nb
    rows[[k - p]] <- c(phi_y, phi_u)
    Tvec[k - p]   <- y[k]
  }
  Phi <- do.call(rbind, rows)
  # Ridge regularization for numerical stability across all 41 outputs
  lambda <- 1e-4
  XtX    <- t(Phi) %*% Phi + lambda * diag(ncol(Phi))
  theta  <- as.numeric(solve(XtX, t(Phi) %*% Tvec))
  list(theta = theta, na = na, nb = nb, m = m)
}

predict_arx_miso_onestep <- function(y, U, model) {
  # One-step-ahead prediction:  yhat[k] uses actual y[k-1..k-na] (not yhat).
  # This is the standard residual-generator form for fault detection in
  # closed-loop systems like TEP.  Free-running simulation diverges because
  # XMV inputs are themselves controller responses.
  na <- model$na; nb <- model$nb; m <- model$m
  theta <- model$theta
  n <- length(y)
  p <- max(na, nb)
  stopifnot(nrow(U) == n, n > p)

  yhat <- numeric(n)
  yhat[1:p] <- y[1:p]
  for (k in (p + 1):n) {
    phi_y <- sapply(1:na, function(i) y[k - i])           # ACTUAL past y
    phi_u <- as.numeric(sapply(1:nb, function(j) U[k - j, ]))
    yhat[k] <- sum(theta * c(phi_y, phi_u))
  }
  yhat
}

# Train one ARX surrogate per selected sensor on the normal-operation data
arx_na <- 3L
arx_nb <- 2L
U_train <- D_norm_train[, xmv_cols, drop = FALSE]

cat("  Fitting ARX surrogates on d00 (normal training data)...\n")
arx_models <- lapply(sensor_idx, function(idx) {
  fit_arx_miso(D_norm_train[, idx], U_train, arx_na, arx_nb)
})

# =============================================================================
# 4. RESIDUAL EXTRACTION (apply trained ARX to a test trajectory)
# =============================================================================

residuals_for_block <- function(D_block) {
  # D_block : n × 52 matrix
  # returns : n × length(sensor_idx) residual matrix
  U_b <- D_block[, xmv_cols, drop = FALSE]
  res <- matrix(0, nrow = nrow(D_block), ncol = length(sensor_idx))
  colnames(res) <- names(sensor_idx)
  for (j in seq_along(sensor_idx)) {
    idx     <- sensor_idx[[j]]
    y_actual <- D_block[, idx]
    y_hat    <- predict_arx_miso_onestep(y_actual, U_b, arx_models[[j]])
    res[, j] <- y_actual - y_hat
  }
  res
}

# =============================================================================
# 5. MULTI-REGIME TRAJECTORY CONSTRUCTION
# =============================================================================
#
# Each test file contains a fault that begins at sample 161 (after 8 hours of
# normal operation at 3-min sampling).  We exploit this: the first 160 samples
# of any test file are Nominal, samples 161+ are the labeled fault.
#
# To match the synthetic pipeline (Nominal → SensorNoisy → DynamicsOff → Drift,
# one trajectory per "seed"), we build a synthetic trajectory by stitching
# windows from different test files.

build_trajectory <- function(seed,
                             nom_len = 500L,
                             fault_len = 500L) {
  set.seed(seed)
  start_nom <- sample.int(nrow(D_norm_test) - nom_len, 1)
  blk_nom   <- D_norm_test[start_nom:(start_nom + nom_len - 1), , drop = FALSE]

  pick_fault <- function(D_fault) {
    onset <- 161L
    max_start <- nrow(D_fault) - fault_len
    start_f <- sample(onset:max_start, 1)
    D_fault[start_f:(start_f + fault_len - 1), , drop = FALSE]
  }
  blk_sn <- pick_fault(D_f11)
  blk_do <- pick_fault(D_f01)
  blk_dr <- pick_fault(D_f13)

  D_traj <- rbind(blk_nom, blk_sn, blk_do, blk_dr)
  gt <- c(rep("Nominal",     nom_len),
          rep("SensorNoisy", fault_len),
          rep("DynamicsOff", fault_len),
          rep("Drift",       fault_len))
  list(D = D_traj, gt = gt)
}

# =============================================================================
# 6. HMM TRANSITION EXTRACTOR (re-used verbatim from SimulationPipeline.R)
# =============================================================================

extract_hmm_transition <- function(fm) {
  nst <- fm@nstates
  trans <- matrix(0, nst, nst)
  for (i in 1:nst) trans[i, ] <- getpars(fm@transition[[i]])
  trans
}

# =============================================================================
# 7. MONTE CARLO LOOP
# =============================================================================

n_mc_seeds      <- 10L
base_seed       <- 11L
min_map_quality <- 0.55

trans_accum  <- matrix(0, K_pipe, K_pipe,
                       dimnames = list(states_pipeline, states_pipeline))
n_good_seeds <- 0L

conf_counts  <- matrix(0, K_pipe, K_pipe,
                       dimnames = list(states_pipeline, states_pipeline))

for (seed_offset in 0:(n_mc_seeds - 1L)) {
  seed_i <- base_seed + seed_offset

  tep_log(sprintf("seed %d: start", seed_i))
  traj   <- build_trajectory(seed_i)
  res_i  <- residuals_for_block(traj$D)        # n × 8
  gt_i   <- traj$gt
  tep_log(sprintf("seed %d: residuals done (n=%d)", seed_i, nrow(res_i)))

  # Feature engineering: raw residuals + rolling SD features.
  # Rolling SD over a 25-sample window discriminates "noisy" (SensorNoisy)
  # from "calm" (Nominal) regimes — the instantaneous residual values are
  # similar in mean but the *variability* is the distinguishing signature.
  feat_i <- as.data.frame(res_i)

  rolling_sd <- function(x, w = 25L) {
    n <- length(x)
    out <- numeric(n)
    for (k in seq_len(n)) {
      lo <- max(1L, k - w + 1L)
      out[k] <- sd(x[lo:k])
    }
    out[is.na(out)] <- 0
    out
  }
  for (nm in colnames(res_i)) {
    feat_i[[paste0(nm, "_rsd")]] <- rolling_sd(res_i[, nm])
  }
  feat_i$res_abs_sum <- rowSums(abs(res_i))
  feat_i$res_rms     <- sqrt(rowMeans(res_i^2))

  # Standardize each feature (z-score)
  feat_i <- as.data.frame(scale(feat_i))

  fams <- rep(list(gaussian()), ncol(feat_i))
  resp_formulas <- lapply(names(feat_i),
                          function(nm) as.formula(paste0(nm, " ~ 1")))

  # Multi-start: rebuild the depmix object each time (avoids depmixS4 segfault
  # when the same S4 object is refit) and keep the highest logLik fit.
  fm_candidates <- vector("list", 5L)
  for (rs in 1:5L) {
    set.seed(seed_i * 100L + rs)
    mod_rs <- depmixS4::depmix(
      response = resp_formulas, data = feat_i,
      nstates = K_pipe, family = fams
    )
    fm_candidates[[rs]] <- tryCatch(
      depmixS4::fit(mod_rs, verbose = FALSE,
                    emcontrol = depmixS4::em.control(
                      maxit = 200, tol = 1e-6, classification = "soft")),
      error = function(e) NULL
    )
  }
  fm_candidates <- Filter(Negate(is.null), fm_candidates)
  if (length(fm_candidates) == 0L) {
    cat(sprintf("  [seed %d] HMM fit failed (all 5 restarts) — skipping\n",
                seed_i))
    next
  }
  ll   <- vapply(fm_candidates,
                 function(f) tryCatch(as.numeric(logLik(f)),
                                      error = function(e) -Inf),
                 numeric(1))
  fm_i <- fm_candidates[[which.max(ll)]]

  vpath_i <- depmixS4::viterbi(fm_i)$state
  tab_i   <- table(State = vpath_i, Regime = gt_i)
  C_i     <- as.matrix(tab_i)
  nr <- nrow(C_i); nc <- ncol(C_i)
  if (nr > nc) C_i <- cbind(C_i, matrix(0, nrow = nr, ncol = nr - nc))
  if (nc > nr) C_i <- rbind(C_i, matrix(0, nrow = nc - nr, ncol = nc))
  perm_i <- clue::solve_LSAP(max(C_i) - C_i)
  regime_names_i <- c(colnames(tab_i),
                      rep("<pad>", ncol(C_i) - ncol(tab_i)))
  map_i <- regime_names_i[perm_i][seq_len(nrow(tab_i))]
  names(map_i) <- rownames(tab_i)

  mapped_path <- map_i[as.character(vpath_i)]
  map_quality <- mean(mapped_path == gt_i, na.rm = TRUE)

  # Accumulate confusion counts (every seed)
  for (tt in seq_along(mapped_path)) {
    true_s <- gt_i[tt]
    pred_s <- mapped_path[tt]
    if (true_s %in% states_pipeline && pred_s %in% states_pipeline) {
      conf_counts[true_s, pred_s] <- conf_counts[true_s, pred_s] + 1
    }
  }

  # Reorder HMM transition matrix to canonical regime order
  hmm_trans_raw <- extract_hmm_transition(fm_i)
  perm_vec <- as.integer(perm_i)[1:K_pipe]
  hmm_trans_reordered <- matrix(0, K_pipe, K_pipe)
  for (hi in 1:K_pipe) {
    for (hj in 1:K_pipe) {
      ri <- perm_vec[hi]
      rj <- perm_vec[hj]
      hmm_trans_reordered[ri, rj] <- hmm_trans_raw[hi, hj]
    }
  }
  hmm_trans_reordered <- hmm_trans_reordered / rowSums(hmm_trans_reordered)

  diag_ok <- all(diag(hmm_trans_reordered) > 0.1)
  if (map_quality >= min_map_quality && diag_ok) {
    trans_accum  <- trans_accum + hmm_trans_reordered
    n_good_seeds <- n_good_seeds + 1L
    cat(sprintf("  [seed %d] GOOD (%.1f%%) — diag(A_hmm): %s\n",
                seed_i, 100 * map_quality,
                paste(round(diag(hmm_trans_reordered), 4), collapse = ", ")))
  } else {
    reason <- if (!diag_ok) "zero-diagonal state" else "low accuracy"
    cat(sprintf("  [seed %d] SKIP (%s, %.1f%%) — diag: %s\n",
                seed_i, reason, 100 * map_quality,
                paste(round(diag(hmm_trans_reordered), 4), collapse = ", ")))
  }
}

cat(sprintf("  Used %d / %d seeds for transition matrix\n",
            n_good_seeds, n_mc_seeds))

# =============================================================================
# 8. NORMALIZE → A_NoAction, mcda_confusion, mcda_accuracy
# =============================================================================

if (n_good_seeds == 0L) {
  stop("No usable seeds — every HMM fit failed the quality filter. ",
       "Inspect residual features or relax min_map_quality.")
}

A_NoAction <- trans_accum / n_good_seeds
eps <- 1e-3
A_NoAction <- A_NoAction + eps
A_NoAction <- A_NoAction / rowSums(A_NoAction)
dimnames(A_NoAction) <- list(states_pipeline, states_pipeline)

mcda_confusion <- conf_counts / rowSums(conf_counts)
dimnames(mcda_confusion) <- list(states_pipeline, states_pipeline)

mcda_accuracy <- diag(mcda_confusion)
names(mcda_accuracy) <- states_pipeline

# =============================================================================
# 9. SANITY CHECKS & REPORT
# =============================================================================

cat("\n── TEP-derived A_NoAction ──\n")
print(round(A_NoAction, 6))
cat("Row sums:", round(rowSums(A_NoAction), 6), "\n")
cat("Diagonal:", round(diag(A_NoAction), 4), "\n")

cat("\n── TEP-derived mcda_confusion ──\n")
print(round(mcda_confusion, 4))
cat("Row sums:", round(rowSums(mcda_confusion), 6), "\n")

cat("\n── mcda_accuracy (diagonal of confusion) ──\n")
print(round(mcda_accuracy, 4))

stopifnot(all(abs(rowSums(A_NoAction)     - 1) < 1e-6))
stopifnot(all(abs(rowSums(mcda_confusion) - 1) < 1e-6))

cat("\n── Comparison to synthetic-testbed values (this paper) ──\n")
# Synthetic-testbed accuracies = diag of O from R/SimulationPipeline.R
# (read from output/data_driven/ if available; values of the paper otherwise)
dd_file <- out_path("data_driven", "data_driven_matrices.rds")
acc_synth <- if (file.exists(dd_file)) readRDS(dd_file)$mcda_accuracy else
  c(Nominal = 0.942, SensorNoisy = 0.210, DynamicsOff = 0.799, Drift = 0.996)
cat("mcda_accuracy — synthetic:", round(acc_synth, 4),
    "\n                       TEP:", round(mcda_accuracy, 4), "\n")

cat("\n── TEP_Pipeline complete. ──\n")
