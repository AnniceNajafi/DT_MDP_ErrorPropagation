#' SimulationPipeline.R
#' ──────────────────────────────────────────────────────────────────────
#' Data-driven extraction of MDP/POMDP parameters from the
#' physics → ARX → HMM pipeline.
#'
#' Runs Monte Carlo simulation (10 seeds) to derive:
#'   - A_NoAction  : 4×4 baseline transition matrix
#'   - mcda_confusion : 4×4 observation/confusion matrix
#'   - mcda_accuracy  : per-state classification accuracy (diagonal)
#'
#' Adapted from ConfusionMatrix.R.  Designed to be source()'d by
#' BaseMDPFormulation.R and other downstream scripts; can also be run on
#' its own (Rscript R/SimulationPipeline.R from the repository root).
#'
#' Writes output/data_driven/{A_NoAction.csv, O_confusion.csv,
#' data_driven_matrices.rds}. If DTMDP_REUSE_MATRICES=1 and the .rds
#' exists, the Monte Carlo/HMM loop (sections 5-6) is skipped and the
#' cached matrices are loaded instead; all helper functions are still
#' defined.
#'
#' Authors@R: person("Annice", "Najafi")
#' ──────────────────────────────────────────────────────────────────────

# Repository-relative paths (see R/paths.R)
if (!exists("repo_path")) local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })

suppressPackageStartupMessages({
  library(depmixS4)
  library(clue)
  library(dplyr)
  library(tibble)
})

cat("── SimulationPipeline: extracting data-driven parameters ──\n")

# =============================================================================
# 1. SIMULATION PARAMETERS
# =============================================================================

pipe_params <- list(
  dt   = 0.02,
  T    = 20,
  tau_cmd      = 0.15,
  cmd_proc_std = 0.005,
  m = 1.0, c = 0.4, k = 4.0, plant_proc_std = 0.04,
  tau_sa = 0.08, sa_meas_std = 0.03, sa_bias = 0.0,
  tau_sb = 0.25, sb_meas_std = 0.02, sb_bias = 0.01,
  alpha = 0.6, w1 = 1.0, w2 = 0.2,
  train_frac = 0.5, na_cmd = 1L, nb_cmd = 1L,
  na_plant = 2L, nb_plant = 2L,
  na_sens  = 1L, nb_sens  = 1L,
  K = 4
)

pipe_schedule <- tibble::tribble(
  ~regime,        ~t_start, ~t_end,
  "Nominal",         0.0,      5.0,
  "SensorNoisy",     5.0,     10.0,
  "DynamicsOff",    10.0,     14.5,
  "Drift",          14.5,     pipe_params$T
)

# Canonical state ordering (must match MDP)
states_pipeline <- c("Nominal", "SensorNoisy", "DynamicsOff", "Drift")
K_pipe <- length(states_pipeline)

# =============================================================================
# 2. PHYSICS + ARX HELPER FUNCTIONS
# =============================================================================

make_u_cmd_pipe <- function(n, amp = 1.0, dt = 0.02) {
  t <- (0:(n-1)) * dt
  u <- 0.6 * sin(2*pi*0.25*t) + 0.4 * sign(sin(2*pi*0.05*t))
  step_every <- as.integer(2.5/dt)
  for (i in seq(1, n, by = step_every)) u[i:n] <- u[i:n] + 0.2 * rnorm(1)
  amp * u
}

lpf_first_order_pipe <- function(u, tau, dt, proc_std = 0.0) {
  n <- length(u); y <- numeric(n)
  a <- exp(-dt / max(tau, 1e-6))
  for (k in 2:n) y[k] <- a*y[k-1] + (1 - a)*u[k-1] + proc_std*rnorm(1)
  y
}

plant_msd_pipe <- function(u, dt, m, c, k, proc_std = 0.0) {
  n <- length(u); x <- numeric(n); v <- numeric(n)
  for (kk in 2:n) {
    acc <- (u[kk-1] - c*v[kk-1] - k*x[kk-1]) / m + proc_std*rnorm(1)
    v[kk] <- v[kk-1] + dt*acc
    x[kk] <- x[kk-1] + dt*v[kk-1]
  }
  list(pos = x, vel = v)
}

sensor_first_order_pipe <- function(xpos, dt, tau, meas_std = 0.0, bias = 0.0) {
  y <- lpf_first_order_pipe(xpos, tau, dt, 0)
  y + bias + rnorm(length(y), sd = meas_std)
}

fit_arx_pipe <- function(y, u, na, nb) {
  stopifnot(length(y) == length(u))
  kmax <- length(y); p <- max(na, nb)
  if (kmax <= p) stop("ARX series too short.")
  rows <- vector("list", kmax - p); Tvec <- numeric(kmax - p)
  for (k in (p+1):kmax) {
    phi <- c(sapply(1:na, function(i) y[k-i]),
             sapply(1:nb, function(j) u[k-j]))
    rows[[k - p]] <- phi
    Tvec[k - p]  <- y[k]
  }
  Phi   <- do.call(rbind, rows)
  theta <- as.numeric(solve(t(Phi) %*% Phi, t(Phi) %*% Tvec))
  list(theta = theta, na = na, nb = nb)
}

sim_arx_pipe <- function(u, y_seed, model) {
  na <- model$na; nb <- model$nb; theta <- model$theta
  n <- length(u); yhat <- numeric(n)
  stopifnot(length(y_seed) >= na)
  yhat[1:na] <- tail(y_seed, na)
  p <- max(na, nb)
  for (k in (p+1):n) {
    phi <- c(sapply(1:na, function(i) yhat[k-i]),
             sapply(1:nb, function(j) u[k-j]))
    yhat[k] <- sum(theta * phi)
  }
  yhat
}

# =============================================================================
# 3. FULL SIMULATION WITH FAULT INJECTION
# =============================================================================

simulate_with_faults_pipe <- function(params, schedule, seed) {
  set.seed(seed)
  n <- as.integer(params$T / params$dt)
  t <- (0:(n-1)) * params$dt
  u_cmd <- make_u_cmd_pipe(n, dt = params$dt)

  u_act_P <- lpf_first_order_pipe(u_cmd, params$tau_cmd, params$dt, params$cmd_proc_std)
  pl      <- plant_msd_pipe(u_act_P, params$dt, params$m, params$c, params$k, params$plant_proc_std)
  x_pos_P <- pl$pos; x_vel_P <- pl$vel
  yA_P    <- sensor_first_order_pipe(x_pos_P, params$dt, params$tau_sa, params$sa_meas_std, params$sa_bias)
  yB_P    <- sensor_first_order_pipe(x_pos_P, params$dt, params$tau_sb, params$sb_meas_std, params$sb_bias)
  y_fused_P <- params$alpha * yA_P + (1 - params$alpha) * yB_P
  z_P <- params$w1 * x_pos_P + params$w2 * x_vel_P

  split <- as.integer(params$train_frac * n)
  arx_cmd   <- fit_arx_pipe(u_act_P[1:split], u_cmd[1:split],     params$na_cmd,   params$nb_cmd)
  arx_plant <- fit_arx_pipe(x_pos_P[1:split], u_act_P[1:split],   params$na_plant, params$nb_plant)
  arx_sa    <- fit_arx_pipe(yA_P[1:split],    x_pos_P[1:split],   params$na_sens,  params$nb_sens)
  arx_sb    <- fit_arx_pipe(yB_P[1:split],    x_pos_P[1:split],   params$na_sens,  params$nb_sens)

  u_act_A <- sim_arx_pipe(u_cmd,   y_seed = u_act_P[1:params$na_cmd],   arx_cmd)
  xpos_A  <- sim_arx_pipe(u_act_A, y_seed = x_pos_P[1:params$na_plant], arx_plant)
  xvel_A  <- c(0, diff(xpos_A) / params$dt)
  yA_A    <- sim_arx_pipe(xpos_A,  y_seed = yA_P[1:params$na_sens],     arx_sa)
  yB_A    <- sim_arx_pipe(xpos_A,  y_seed = yB_P[1:params$na_sens],     arx_sb)
  y_fused_A <- params$alpha * yA_A + (1 - params$alpha) * yB_A
  z_A <- params$w1 * xpos_A + params$w2 * xvel_A

  regime_by_t <- rep("Nominal", n)
  for (r in seq_len(nrow(schedule))) {
    idx <- which(t >= schedule$t_start[r] & t < schedule$t_end[r])
    regime_by_t[idx] <- schedule$regime[r]
  }

  x_pos_A_fault <- xpos_A
  yA_A_fault    <- yA_A
  yB_A_fault    <- yB_A

  for (k in seq_len(n)) {
    reg <- regime_by_t[k]
    if (reg == "SensorNoisy") {
      yA_A_fault[k] <- yA_A_fault[k] + rnorm(1, 0, 0.06)
      yB_A_fault[k] <- yB_A_fault[k] + rnorm(1, 0, 0.05)
    } else if (reg == "DynamicsOff") {
      x_pos_A_fault[k] <- 0.8 * x_pos_A_fault[k]
    } else if (reg == "Drift") {
      yA_A_fault[k] <- yA_A_fault[k] + 0.0007 * (k - 1)
    }
  }
  y_fused_A_fault <- params$alpha * yA_A_fault + (1 - params$alpha) * yB_A_fault
  z_A_fault <- params$w1 * x_pos_A_fault + params$w2 * c(0, diff(x_pos_A_fault) / params$dt)

  E <- tibble::tibble(
    t   = t,
    e0  = u_act_A          - u_act_P,
    e1  = x_pos_A_fault    - x_pos_P,
    e2a = yA_A_fault       - yA_P,
    e2b = yB_A_fault       - yB_P,
    e3  = y_fused_A_fault  - y_fused_P,
    e4  = z_A_fault        - z_P
  )

  list(E = E, t = t, regime_t = regime_by_t)
}

# =============================================================================
# 4. HELPER: Extract HMM transition matrix from depmixS4 fitted model
# =============================================================================
#
# The HMM's EM-learned transition matrix gives P(s_{t+1}|s_t) — the proper
# probabilistic estimate of state dynamics.  This is more principled than
# counting Viterbi-path transitions, which are biased by the fixed schedule
# (e.g. the last regime appears absorbing because the simulation ends there).

extract_hmm_transition <- function(fm) {
  nst <- fm@nstates
  trans <- matrix(0, nst, nst)
  for (i in 1:nst) {
    # getpars on a transInit object returns the actual probabilities directly
    trans[i, ] <- getpars(fm@transition[[i]])
  }
  trans
}

dd_cache_file <- out_path("data_driven", "data_driven_matrices.rds")
reuse_matrices <- identical(Sys.getenv("DTMDP_REUSE_MATRICES"), "1") &&
                  file.exists(dd_cache_file)

if (reuse_matrices) {
  cat("  Reusing cached data-driven matrices:", dd_cache_file, "\n")
  dd_cached      <- readRDS(dd_cache_file)
  A_NoAction     <- dd_cached$A_NoAction
  mcda_confusion <- dd_cached$mcda_confusion
  mcda_accuracy  <- dd_cached$mcda_accuracy
  n_good_seeds   <- dd_cached$n_good_seeds
} else {

# =============================================================================
# 5. MONTE CARLO LOOP: 10 seeds → transition matrices & confusion counts
# =============================================================================

n_mc_seeds  <- 10
base_seed   <- 11L
min_map_quality <- 0.55  # minimum mapping accuracy to include a seed

# Accumulate HMM-learned transition matrices (from quality-filtered seeds)
trans_accum <- matrix(0, K_pipe, K_pipe,
                      dimnames = list(states_pipeline, states_pipeline))
n_good_seeds <- 0

# Accumulate confusion counts (Viterbi vs ground truth — all seeds)
conf_counts <- matrix(0, K_pipe, K_pipe,
                      dimnames = list(states_pipeline, states_pipeline))

for (seed_offset in 0:(n_mc_seeds - 1)) {
  seed_i <- base_seed + seed_offset

  sim_i <- simulate_with_faults_pipe(pipe_params, pipe_schedule, seed_i)
  E_i   <- sim_i$E
  gt_i  <- sim_i$regime_t

  # 10 z-scored residual features
  feat_i <- E_i %>%
    dplyr::transmute(
      t,
      e0 = e0, e1 = e1, e2a = e2a, e2b = e2b, e3 = e3, e4 = e4,
      e_abs_sum = abs(e0) + abs(e1) + abs(e2a) + abs(e2b) + abs(e3) + abs(e4),
      e_rms     = sqrt((e0^2 + e1^2 + e2a^2 + e2b^2 + e3^2 + e4^2) / 6),
      e1_diff   = c(0, diff(e1)),
      e3_diff   = c(0, diff(e3))
    ) %>%
    dplyr::mutate(dplyr::across(-t, ~ as.numeric(scale(.x))))

  Ez_i <- feat_i %>% dplyr::select(-t) %>% as.data.frame()
  fams <- rep(list(gaussian()), ncol(Ez_i))
  mod_i <- depmixS4::depmix(
    response = lapply(names(Ez_i), function(nm) as.formula(paste0(nm, "~ 1"))),
    data = Ez_i, nstates = pipe_params$K, family = fams
  )
  fm_i <- tryCatch(depmixS4::fit(mod_i, verbose = FALSE), error = function(e) NULL)
  if (is.null(fm_i)) {
    cat(sprintf("  [seed %d] HMM fit failed — skipping\n", seed_i))
    next
  }

  # Viterbi + Hungarian mapping
  vpath_i <- depmixS4::viterbi(fm_i)$state
  tab_i   <- table(State = vpath_i, Regime = gt_i)
  C_i     <- as.matrix(tab_i)
  nr <- nrow(C_i); nc <- ncol(C_i)
  if (nr > nc) C_i <- cbind(C_i, matrix(0, nrow = nr, ncol = nr - nc))
  if (nc > nr) C_i <- rbind(C_i, matrix(0, nrow = nc - nr, ncol = nc))
  perm_i <- clue::solve_LSAP(max(C_i) - C_i)
  regime_names_i <- c(colnames(tab_i), rep("<pad>", ncol(C_i) - ncol(tab_i)))
  map_i <- regime_names_i[perm_i][seq_len(nrow(tab_i))]
  names(map_i) <- rownames(tab_i)

  # Mapping quality = fraction of timesteps correctly mapped
  mapped_path <- map_i[as.character(vpath_i)]
  map_quality <- mean(mapped_path == gt_i, na.rm = TRUE)

  # --- Confusion matrix: always accumulate (Viterbi path vs ground truth) ---
  for (tt in seq_along(mapped_path)) {
    true_s <- gt_i[tt]
    pred_s <- mapped_path[tt]
    if (true_s %in% states_pipeline && pred_s %in% states_pipeline) {
      conf_counts[true_s, pred_s] <- conf_counts[true_s, pred_s] + 1
    }
  }

  # --- Transition matrix: extract from HMM model if mapping is good ---
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

  # Quality gate: overall mapping accuracy AND all diagonal entries > 0.1
  # (diagonal = 0 means the HMM state mapped to that regime is spurious)
  diag_ok <- all(diag(hmm_trans_reordered) > 0.1)
  if (map_quality >= min_map_quality && diag_ok) {
    trans_accum <- trans_accum + hmm_trans_reordered
    n_good_seeds <- n_good_seeds + 1
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

cat(sprintf("  Used %d / %d seeds for transition matrix\n", n_good_seeds, n_mc_seeds))

# =============================================================================
# 6. NORMALIZE → A_NoAction, mcda_confusion, mcda_accuracy
# =============================================================================

# Average HMM-learned transition matrices across quality-filtered seeds
A_NoAction <- trans_accum / n_good_seeds

# Laplace smoothing: ensure no entry is 0 or 1 (prevents absorbing states)
# This is critical because the MDP requires all transitions to be possible
eps <- 1e-3
A_NoAction <- A_NoAction + eps
A_NoAction <- A_NoAction / rowSums(A_NoAction)
dimnames(A_NoAction) <- list(states_pipeline, states_pipeline)

# Confusion matrix (from all seeds — represents actual classification quality)
mcda_confusion <- conf_counts / rowSums(conf_counts)
dimnames(mcda_confusion) <- list(states_pipeline, states_pipeline)

mcda_accuracy <- diag(mcda_confusion)
names(mcda_accuracy) <- states_pipeline


# Save the regenerated matrices (repository-relative output directory)
saveRDS(list(A_NoAction = A_NoAction, mcda_confusion = mcda_confusion,
             mcda_accuracy = mcda_accuracy, n_good_seeds = n_good_seeds,
             conf_counts = conf_counts, n_mc_seeds = n_mc_seeds,
             base_seed = base_seed),
        dd_cache_file)
write.csv(round(A_NoAction, 6), out_path("data_driven", "A_NoAction.csv"))
write.csv(round(mcda_confusion, 6), out_path("data_driven", "O_confusion.csv"))
cat("  -> Saved output/data_driven/{A_NoAction.csv, O_confusion.csv, data_driven_matrices.rds}\n")

} # end if (reuse_matrices) ... else

# =============================================================================
# 7. SANITY CHECKS
# =============================================================================

cat("\n── Data-driven A_NoAction ──\n")
print(round(A_NoAction, 6))
cat("Row sums:", round(rowSums(A_NoAction), 6), "\n")
cat("Diagonal:", round(diag(A_NoAction), 4), "\n")

cat("\n── Data-driven mcda_confusion ──\n")
print(round(mcda_confusion, 4))
cat("Row sums:", round(rowSums(mcda_confusion), 6), "\n")

cat("\n── mcda_accuracy (diagonal of confusion) ──\n")
print(round(mcda_accuracy, 4))

# Validate thresholds
stopifnot(all(abs(rowSums(A_NoAction) - 1) < 1e-6))
stopifnot(all(abs(rowSums(mcda_confusion) - 1) < 1e-6))

# Transition matrix: diagonals should dominate (regimes are persistent)
if (!all(diag(A_NoAction) > 0.8)) {
  warning("Some A_NoAction diagonal entries < 0.8: ",
          paste(round(diag(A_NoAction), 4), collapse = ", "))
}
# Confusion matrix: each regime should be classified better than chance (>0.25 for K=4)
if (!all(diag(mcda_confusion) > 0.25)) {
  warning("Some mcda_confusion diagonal entries < 0.25: ",
          paste(round(diag(mcda_confusion), 4), collapse = ", "))
}

# Comparison to old hardcoded values
cat("\n── Comparison to old hardcoded values ──\n")
A_old <- matrix(c(
  0.985, 0.007, 0.005, 0.003,
  0.005, 0.985, 0.005, 0.005,
  0.003, 0.005, 0.987, 0.005,
  0.003, 0.007, 0.005, 0.985
), nrow = 4, byrow = TRUE, dimnames = list(states_pipeline, states_pipeline))

cat("A_NoAction diagonal — old:", round(diag(A_old), 4),
    " new:", round(diag(A_NoAction), 4), "\n")

acc_old <- c(Nominal = 0.89, SensorNoisy = 0.70, DynamicsOff = 0.99, Drift = 0.99)
cat("mcda_accuracy — old:", round(acc_old, 4),
    " new:", round(mcda_accuracy, 4), "\n")

cat("\n── SimulationPipeline complete. ──\n")
