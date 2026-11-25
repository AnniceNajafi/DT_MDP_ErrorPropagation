#' BaseMDPFormulation.R
#' ──────────────────────────────────────────────────────────────────────
#' MDP definition, Value Iteration, continuous-time Bellman solution,
#' Gillespie simulation functions and policy definitions.
#'
#' The Monte Carlo Bellman-vs-Gillespie validation and its three figures
#' (formerly sections 6-7 of this file) are in analyses/mdp_validation.R,
#' so that sourcing this file only builds the model (no PDFs are written).
#'
#' Sourced by R/POMDPFormulation.R. Options:
#'   options(dtmdp.use_existing_matrices = TRUE)  use A_NoAction /
#'       mcda_confusion already in the global environment instead of
#'       sourcing R/SimulationPipeline.R (used by tep/TEP_Application.R)
#'   options(dtmdp.dt = 180)  override the 0.02 s time step (TEP)
#'
#' Authors@R: person("Annice", "Najafi")
#' ──────────────────────────────────────────────────────────────────────

# Repository-relative paths (see R/paths.R)
if (!exists("repo_path")) local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(expm)         # matrix exponential
  library(gridExtra)    # multi-panel figures
  library(RMCDA)        # real TOPSIS baseline (companion paper [4])
})

# --- Data-driven parameters from physics → ARX → HMM pipeline ---
if (!isTRUE(getOption("dtmdp.use_existing_matrices"))) {
  source(repo_path("R", "SimulationPipeline.R"), local = globalenv())
}
# Exports: A_NoAction, mcda_confusion, mcda_accuracy

set.seed(42)

# =====================================================================
# 1. MDP DEFINITION
# =====================================================================

states  <- c("Nominal", "SensorNoisy", "DynamicsOff", "Drift")
actions <- c("NoAction", "DwSensorA", "DwSensorB",
             "IncFilter", "ReidentPlant", "BiasCorrect")
K  <- length(states)
M  <- length(actions)
dt <- getOption("dtmdp.dt", 0.02)   # discrete time-step (s), 50 Hz sampling

# A_NoAction is now provided by SimulationPipeline.R (data-driven)

# --- Repair probabilities ---
rho <- matrix(0, nrow = M, ncol = K,
              dimnames = list(actions, states))
rho["DwSensorA",    "SensorNoisy"]  <- 0.60
rho["DwSensorB",    "SensorNoisy"]  <- 0.60
rho["IncFilter",    "SensorNoisy"]  <- 0.40
rho["ReidentPlant", "DynamicsOff"]  <- 0.70
rho["BiasCorrect",  "Drift"]        <- 0.80

# --- Action-dependent transition matrices ---
A <- list()
for (a in actions) {
  A[[a]] <- A_NoAction
  for (s in 1:K) {
    if (rho[a, s] > 0) {
      r <- rho[a, s]
      A[[a]][s, ] <- (1 - r) * A_NoAction[s, ]
      A[[a]][s, 1] <- A[[a]][s, 1] + r
    }
  }
}

# --- Reward matrix R(action, state) ---
R_mat <- matrix(c(
  # Nominal  SensorNoisy  DynamicsOff  Drift
  +1.0,      -0.5,        -1.0,       -0.8,   # NoAction
  -0.1,      +0.7,        -0.2,       -0.1,   # DwSensorA
  -0.1,      +0.6,        -0.2,       -0.1,   # DwSensorB
  -0.1,      +0.4,         0.0,       +0.1,   # IncFilter
  -0.5,      -0.3,        +0.9,       -0.2,   # ReidentPlant
  -0.2,       0.0,        +0.1,       +0.8    # BiasCorrect
), nrow = M, byrow = TRUE,
dimnames = list(actions, states))

# --- Discount factor ---
gamma_disc <- 0.99

# =====================================================================
# 2. VALUE ITERATION  (Discrete-time Bellman)
# =====================================================================

value_iteration <- function(A_list, R, gamma, tol = 1e-10, max_iter = 10000) {
  K <- ncol(R);  M_ <- nrow(R)
  act <- rownames(R);  st <- colnames(R)
  V <- rep(0, K);  names(V) <- st
  policy <- character(K);  names(policy) <- st
  history <- matrix(NA, max_iter, K, dimnames = list(NULL, st))

  for (it in 1:max_iter) {
    V_old <- V
    Q <- matrix(NA, M_, K, dimnames = list(act, st))
    for (ai in 1:M_) {
      a <- act[ai]
      Q[ai, ] <- R[ai, ] + gamma * (A_list[[a]] %*% V)
    }
    for (s in 1:K) {
      best   <- which.max(Q[, s])
      V[s]      <- Q[best, s]
      policy[s] <- act[best]
    }
    history[it, ] <- V
    if (max(abs(V - V_old)) < tol) break
  }
  list(V = V, policy = policy, Q = Q,
       history = history[1:it, , drop = FALSE], iterations = it)
}

vi <- value_iteration(A, R_mat, gamma_disc)
cat("Value Iteration converged in", vi$iterations, "iterations.\n")
cat("Optimal policy:\n"); print(vi$policy)
cat("Optimal values:\n"); print(round(vi$V, 4))

# =====================================================================
# 3. CONTINUOUS-TIME FORMULATION
# =====================================================================

# Rate matrices  Q = (A - I) / dt
Q_rate <- lapply(A, function(Aa) {
  Qa <- (Aa - diag(K)) / dt
  dimnames(Qa) <- list(states, states)
  Qa
})

r_disc <- -log(gamma_disc) / dt     # continuous discount rate
c_rate <- R_mat / dt                 # reward rate (per second)

# Continuous-time Bellman solution  V = (rI - Q_pi)^{-1} c_pi
Q_pi <- matrix(0, K, K, dimnames = list(states, states))
c_pi <- numeric(K); names(c_pi) <- states

for (s in 1:K) {
  a <- vi$policy[states[s]]
  Q_pi[s, ] <- Q_rate[[a]][s, ]
  c_pi[s]   <- c_rate[a, s]
}

# (tryCatch: for the TEP time step the system can be numerically singular;
#  V_ct is not used downstream, so fall back to the discrete V*.)
V_ct <- tryCatch(solve(r_disc * diag(K) - Q_pi, c_pi), error = function(e) {
  cat("  (V_ct solve singular -- using discrete V* as fallback)\n"); vi$V })
names(V_ct) <- states
cat("\nContinuous-time Bellman values:\n"); print(round(V_ct, 4))

# =====================================================================
# 4. GILLESPIE SIMULATION FUNCTIONS
# =====================================================================

# --- Gillespie under a fixed deterministic policy ---
gillespie_mdp <- function(Q_list, policy, R_rate, s0, T_end, r_disc) {
  K  <- nrow(Q_list[[1]])
  sn <- rownames(Q_list[[1]])

  cap <- 2000
  times  <- numeric(cap)
  st_seq <- integer(cap)
  idx    <- 1
  t_now  <- 0
  s      <- s0

  times[1]  <- 0
  st_seq[1] <- s

  while (t_now < T_end) {
    a  <- policy[sn[s]]
    Qa <- Q_list[[a]]
    ro <- pmax(Qa[s, ], 0);  ro[s] <- 0
    lam <- sum(ro)
    if (lam <= 0) break

    tau <- rexp(1, lam)
    t_now <- t_now + tau
    if (t_now > T_end) break

    s <- sample.int(K, 1, prob = ro)

    idx <- idx + 1
    if (idx > length(times)) {
      times  <- c(times,  numeric(cap))
      st_seq <- c(st_seq, integer(cap))
    }
    times[idx]  <- t_now
    st_seq[idx] <- s
  }

  times  <- times[1:idx]
  st_seq <- st_seq[1:idx]

  # Discounted return
  disc_ret <- 0
  for (i in seq_along(times)) {
    a_i     <- policy[sn[st_seq[i]]]
    c_i     <- R_rate[a_i, st_seq[i]]
    t_start <- times[i]
    t_end_i <- if (i < idx) times[i + 1] else T_end
    dt_seg  <- t_end_i - t_start
    if (dt_seg > 0 && r_disc > 0)
      disc_ret <- disc_ret +
        (c_i / r_disc) * exp(-r_disc * t_start) * (1 - exp(-r_disc * dt_seg))
  }

  list(times = times, states = st_seq, disc_return = disc_ret)
}

# --- Gillespie with random policy ---
gillespie_random <- function(Q_list, R_rate, s0, T_end, r_disc) {
  K  <- nrow(Q_list[[1]])
  sn <- rownames(Q_list[[1]])
  act_names <- names(Q_list)
  t_now <- 0;  s <- s0;  disc_ret <- 0

  while (t_now < T_end) {
    a  <- sample(act_names, 1)
    Qa <- Q_list[[a]]
    ro <- pmax(Qa[s, ], 0);  ro[s] <- 0
    lam <- sum(ro);  if (lam <= 0) break

    tau <- rexp(1, lam)
    t_end_seg <- min(t_now + tau, T_end)
    dt_seg    <- t_end_seg - t_now

    ci <- R_rate[a, s]
    disc_ret <- disc_ret +
      (ci / r_disc) * exp(-r_disc * t_now) * (1 - exp(-r_disc * dt_seg))

    if (t_now + tau >= T_end) break
    t_now <- t_now + tau
    s <- sample.int(K, 1, prob = ro)
  }
  disc_ret
}

# --- Gillespie with real MCDA/TOPSIS policy (companion paper [4]) ---
# At each transition event, an EMA posterior over regimes is updated
# from the noisy observation z. Five criteria (KPI_gain, Module_gain,
# Drift_gain, Stability, Cost) are computed per repair action from the
# posterior, and RMCDA::apply.TOPSIS ranks the five repair actions.
# A gate forces NoAction when the posterior is confidently Nominal and
# KPI proxy is low. Time-based dwell hysteresis (dwell_sec) prevents
# rapid action switching. Criteria weights, action costs, and stability
# penalties match MCDA_intervention.R from the companion HMM paper.
#
# Signature is kept identical to the legacy shortcut so existing callers
# do not need to change. The mcda_action argument is ignored (TOPSIS
# computes actions from criteria; it is retained only for backward
# compatibility with notebooks passing vi$policy at call sites).
gillespie_mcda <- function(Q_list, mcda_action, confusion, dwell_sec,
                           R_rate, s0, T_end, r_disc,
                           alpha_ema = 0.3) {
  K  <- nrow(Q_list[[1]])
  sn <- rownames(Q_list[[1]])

  topsis_actions <- c("DwSensorA", "DwSensorB", "IncFilter",
                      "ReidentPlant", "BiasCorrect")
  stability_penalty <- c(0.08, 0.08, 0.05, 0.15, 0.12)
  action_cost       <- c(0.10, 0.10, 0.15, 0.40, 0.25)
  topsis_weights    <- c(KPI_gain = 0.35, Module_gain = 0.25,
                         Drift_gain = 0.15, Stability = 0.10, Cost = 0.15)

  cap     <- 2000
  times   <- numeric(cap)
  st_seq  <- integer(cap)
  act_seq <- character(cap)
  idx     <- 1
  t_now   <- 0
  s       <- s0

  posterior <- rep(0, K)
  z0 <- sample.int(K, 1, prob = confusion[s0, ])
  posterior[z0] <- 1.0

  current_action   <- "NoAction"
  last_switch_time <- 0

  times[1]   <- 0
  st_seq[1]  <- s
  act_seq[1] <- current_action

  while (t_now < T_end) {
    a  <- current_action
    Qa <- Q_list[[a]]
    ro <- pmax(Qa[s, ], 0); ro[s] <- 0
    lam <- sum(ro)
    if (lam <= 0) break

    tau   <- rexp(1, lam)
    t_now <- t_now + tau
    if (t_now > T_end) break

    s <- sample.int(K, 1, prob = ro)
    z <- sample.int(K, 1, prob = confusion[s, ])

    indicator <- rep(0, K); indicator[z] <- 1.0
    posterior <- alpha_ema * indicator + (1 - alpha_ema) * posterior
    posterior <- posterior / sum(posterior)

    pN <- posterior[1]; pS <- posterior[2]
    pD <- posterior[3]; pR <- posterior[4]

    sensor_n <- pS
    plant_n  <- pD
    drift_n  <- pR
    kpi_n    <- 1 - pN
    noise_n  <- sensor_n

    kpi_gain <- c(
      pS * kpi_n * sensor_n * 1.2,
      pS * kpi_n * sensor_n * 0.6,
      pS * kpi_n * noise_n,
      pD * kpi_n * plant_n,
      (pR + 0.5 * pD) * kpi_n * (drift_n + 0.3 * plant_n)
    )
    module_gain <- c(
      pS * sensor_n * 1.2,
      pS * sensor_n * 0.6,
      pS * sensor_n * noise_n,
      pD * plant_n,
      (pR + 0.5 * pD) * (drift_n + 0.3 * plant_n)
    )
    drift_gain <- c(
      pS * drift_n * 1.1,
      pS * drift_n * 0.7,
      pS * drift_n * 0.6,
      pD * drift_n * 0.4,
      (pR + 0.2 * pS) * drift_n
    )

    dm <- cbind(
      KPI_gain    = kpi_gain,
      Module_gain = module_gain,
      Drift_gain  = drift_gain,
      Stability   = -stability_penalty,
      Cost        = -action_cost
    )
    rownames(dm) <- topsis_actions

    topsis_scores <- suppressWarnings(
      as.numeric(RMCDA::apply.TOPSIS(dm, w = as.numeric(topsis_weights)))
    )
    if (all(is.na(topsis_scores)) || all(topsis_scores == 0)) {
      desired_action <- "NoAction"
    } else {
      desired_action <- topsis_actions[which.max(topsis_scores)]
    }

    # Gate: confidently Nominal + low KPI proxy => NoAction
    if (kpi_n < 0.25 && pN > 0.75) {
      desired_action <- "NoAction"
    }

    if (desired_action != current_action &&
        (t_now - last_switch_time) >= dwell_sec) {
      current_action   <- desired_action
      last_switch_time <- t_now
    }

    idx <- idx + 1
    if (idx > length(times)) {
      times   <- c(times,   numeric(cap))
      st_seq  <- c(st_seq,  integer(cap))
      act_seq <- c(act_seq, character(cap))
    }
    times[idx]   <- t_now
    st_seq[idx]  <- s
    act_seq[idx] <- current_action
  }

  times   <- times[1:idx]
  st_seq  <- st_seq[1:idx]
  act_seq <- act_seq[1:idx]

  disc_ret <- 0
  for (i in seq_along(times)) {
    a_i     <- act_seq[i]
    c_i     <- R_rate[a_i, st_seq[i]]
    t_start <- times[i]
    t_end_i <- if (i < idx) times[i + 1] else T_end
    dt_seg  <- t_end_i - t_start
    if (dt_seg > 0 && r_disc > 0)
      disc_ret <- disc_ret +
        (c_i / r_disc) * exp(-r_disc * t_start) * (1 - exp(-r_disc * dt_seg))
  }

  list(times = times, states = st_seq, actions = act_seq, disc_return = disc_ret)
}

# --- Legacy shortcut: MDP-optimal policy applied to hard classifications ---
# Retained for ablation / comparison purposes. Not used by the figure
# notebooks unless explicitly called. Uses the same signature as
# gillespie_mcda so it is a drop-in ablation swap.
gillespie_hard_mcda <- function(Q_list, mcda_action, confusion, dwell_sec,
                                R_rate, s0, T_end, r_disc) {
  K  <- nrow(Q_list[[1]])
  sn <- rownames(Q_list[[1]])

  cap <- 2000
  times   <- numeric(cap)
  st_seq  <- integer(cap)
  act_seq <- character(cap)
  idx     <- 1
  t_now   <- 0
  s       <- s0

  times[1]  <- 0
  st_seq[1] <- s

  s_perceived <- sample.int(K, 1, prob = confusion[s, ])
  current_action <- mcda_action[sn[s_perceived]]
  act_seq[1] <- current_action
  last_switch_time <- 0

  while (t_now < T_end) {
    a  <- current_action
    Qa <- Q_list[[a]]
    ro <- pmax(Qa[s, ], 0);  ro[s] <- 0
    lam <- sum(ro)
    if (lam <= 0) break

    tau <- rexp(1, lam)
    t_now <- t_now + tau
    if (t_now > T_end) break

    s <- sample.int(K, 1, prob = ro)

    s_perceived <- sample.int(K, 1, prob = confusion[s, ])
    desired_action <- mcda_action[sn[s_perceived]]

    if (desired_action != current_action &&
        (t_now - last_switch_time) >= dwell_sec) {
      current_action <- desired_action
      last_switch_time <- t_now
    }

    idx <- idx + 1
    if (idx > length(times)) {
      times   <- c(times,   numeric(cap))
      st_seq  <- c(st_seq,  integer(cap))
      act_seq <- c(act_seq, character(cap))
    }
    times[idx]   <- t_now
    st_seq[idx]  <- s
    act_seq[idx] <- current_action
  }

  times   <- times[1:idx]
  st_seq  <- st_seq[1:idx]
  act_seq <- act_seq[1:idx]

  disc_ret <- 0
  for (i in seq_along(times)) {
    a_i     <- act_seq[i]
    c_i     <- R_rate[a_i, st_seq[i]]
    t_start <- times[i]
    t_end_i <- if (i < idx) times[i + 1] else T_end
    dt_seg  <- t_end_i - t_start
    if (dt_seg > 0 && r_disc > 0)
      disc_ret <- disc_ret +
        (c_i / r_disc) * exp(-r_disc * t_start) * (1 - exp(-r_disc * dt_seg))
  }

  list(times = times, states = st_seq, actions = act_seq, disc_return = disc_ret)
}

# --- Cumulative reward at specified time points (for a deterministic-policy trajectory) ---
cum_reward_at_pts <- function(tr, t_pts, policy, c_rate, r_disc) {
  cr <- 0;  prev_t <- 0;  prev_s <- tr$states[1];  ti <- 1
  out <- numeric(length(t_pts))
  for (b in seq_along(t_pts)) {
    while (ti < length(tr$times) && tr$times[ti + 1] <= t_pts[b]) {
      seg <- tr$times[ti + 1] - prev_t
      a_i <- policy[states[prev_s]]
      ci  <- c_rate[a_i, prev_s]
      cr  <- cr + (ci / r_disc) *
             exp(-r_disc * prev_t) * (1 - exp(-r_disc * seg))
      prev_t <- tr$times[ti + 1];  ti <- ti + 1;  prev_s <- tr$states[ti]
    }
    seg <- t_pts[b] - prev_t
    if (seg > 0) {
      a_i <- policy[states[prev_s]]
      ci  <- c_rate[a_i, prev_s]
      out[b] <- cr + (ci / r_disc) *
        exp(-r_disc * prev_t) * (1 - exp(-r_disc * seg))
    } else {
      out[b] <- cr
    }
  }
  out
}

# --- Cumulative reward at time points for MCDA trajectory ---
# Uses the stored action sequence (tr$actions) rather than re-deriving from true state
cum_reward_at_pts_mcda <- function(tr, t_pts, mcda_act, c_rate, r_disc) {
  out <- numeric(length(t_pts))
  idx <- 1
  cum <- 0
  for (p in seq_along(t_pts)) {
    while (idx < length(tr$times) && tr$times[idx + 1] <= t_pts[p]) {
      a_i     <- tr$actions[idx]
      t_start <- tr$times[idx]
      t_end_i <- tr$times[idx + 1]
      dt_seg  <- t_end_i - t_start
      c_i     <- c_rate[a_i, tr$states[idx]]
      if (dt_seg > 0 && r_disc > 0)
        cum <- cum + (c_i / r_disc) * exp(-r_disc * t_start) *
                     (1 - exp(-r_disc * dt_seg))
      idx <- idx + 1
    }
    a_i     <- tr$actions[idx]
    c_i     <- c_rate[a_i, tr$states[idx]]
    t_start <- tr$times[idx]
    dt_seg  <- t_pts[p] - t_start
    if (dt_seg > 0 && r_disc > 0)
      out[p] <- cum + (c_i / r_disc) * exp(-r_disc * t_start) *
                      (1 - exp(-r_disc * dt_seg))
    else
      out[p] <- cum
  }
  out
}

# --- Bin a trajectory into state indices at given time points ---
bin_trajectory <- function(tr, t_bins, K) {
  out <- integer(length(t_bins))
  for (b in seq_along(t_bins)) {
    j <- max(which(tr$times <= t_bins[b]))
    out[b] <- tr$states[j]
  }
  out
}

# --- Fraction of time in Nominal for one trajectory ---
# (n_ev > 1 guard: a trajectory with no transition before T_end has a
#  single event; the original 1:(n_ev - 1) loop indexed states[0].)
frac_in_nominal <- function(tr, T_end) {
  occ <- 0
  n_ev <- length(tr$times)
  if (n_ev > 1) {
    for (i in seq_len(n_ev - 1)) {
      if (tr$states[i] == 1)
        occ <- occ + (tr$times[i + 1] - tr$times[i])
    }
  }
  dt_last <- T_end - tr$times[n_ev]
  if (tr$states[n_ev] == 1) occ <- occ + dt_last
  occ / T_end
}

# =====================================================================
# 5. POLICY DEFINITIONS
# =====================================================================

# No-intervention policy
policy_none <- setNames(rep("NoAction", K), states)

# MCDA policy parameters
# mcda_accuracy and mcda_confusion are now provided by SimulationPipeline.R

mcda_ideal_action <- vi$policy
mcda_dwell_time   <- 0.50

# Reward-rate vector under optimal policy (for theoretical curves)
R_pi_vec <- sapply(1:K, function(s) c_rate[vi$policy[states[s]], s])

cat("\n══ BaseMDPFormulation (model definition) complete. ══\n")
