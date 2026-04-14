#MDP_POMDP_Formulation.R
#Main formulation script: MDP definition, Value Iteration,
#continuous-time Bellman solution, Gillespie simulation functions,
#POMDP observation model, PBVI solver, and belief tracking.
#
#This script sources SimulationPipeline.R for data-driven parameters
#and exports all model objects and functions needed by figure notebooks.
#
#Authors: Annice Najafi

library(ggplot2)
library(reshape2)
library(expm)
library(gridExtra)
library(RMCDA)     #real TOPSIS baseline (companion paper [4])

#Data-driven parameters from physics -> ARX -> HMM pipeline
source("SimulationPipeline.R")
#Exports: A_NoAction, mcda_confusion, mcda_accuracy

set.seed(42)

#=====================================================================
#1. MDP DEFINITION
#=====================================================================

states  <- c("Nominal", "SensorNoisy", "DynamicsOff", "Drift")
actions <- c("NoAction", "DwSensorA", "DwSensorB",
             "IncFilter", "ReidentPlant", "BiasCorrect")
K  <- length(states)
M  <- length(actions)
dt <- 0.02  #discrete time-step (s), 50 Hz sampling

#A_NoAction is provided by SimulationPipeline.R (data-driven)

#Repair probabilities rho(action, state)
#Each entry is the probability that action a repairs fault state s
rho <- matrix(0, nrow = M, ncol = K,
              dimnames = list(actions, states))
rho["DwSensorA",    "SensorNoisy"]  <- 0.60
rho["DwSensorB",    "SensorNoisy"]  <- 0.60
rho["IncFilter",    "SensorNoisy"]  <- 0.40
rho["ReidentPlant", "DynamicsOff"]  <- 0.70
rho["BiasCorrect",  "Drift"]        <- 0.80

#Action-dependent transition matrices A(a)
#When action a has repair probability rho[a,s] > 0 for state s,
#the row is modified: with probability rho the system returns to
#Nominal, with probability (1-rho) it follows the natural dynamics.
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

#Reward matrix R(action, state)
#Rows = actions, Columns = states
#Positive rewards for correct interventions, penalties for
#wrong actions or unaddressed faults
R_mat <- matrix(c(
  #Nominal  SensorNoisy  DynamicsOff  Drift
  +1.0,      -0.5,        -1.0,       -0.8,   #NoAction
  -0.1,      +0.7,        -0.2,       -0.1,   #DwSensorA
  -0.1,      +0.6,        -0.2,       -0.1,   #DwSensorB
  -0.1,      +0.4,         0.0,       +0.1,   #IncFilter
  -0.5,      -0.3,        +0.9,       -0.2,   #ReidentPlant
  -0.2,       0.0,        +0.1,       +0.8    #BiasCorrect
), nrow = M, byrow = TRUE,
dimnames = list(actions, states))

#Discount factor
gamma_disc <- 0.99

#=====================================================================
#2. VALUE ITERATION (Discrete-time Bellman)
#=====================================================================
#Solves V*(s) = max_a { R(a,s) + gamma * sum_s' A(a)[s,s'] V*(s') }
#Returns optimal value function V*, policy pi*, Q-values, and
#convergence history.

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
      best      <- which.max(Q[, s])
      V[s]      <- Q[best, s]
      policy[s] <- act[best]
    }
    history[it, ] <- V
    if (max(abs(V - V_old)) < tol) break
  }
  list(V = V, policy = policy, Q = Q,
       history = history[1:it, , drop = FALSE], iterations = it)
}

#Solve the MDP
vi <- value_iteration(A, R_mat, gamma_disc)
print(paste0("Value Iteration converged in ", vi$iterations, " iterations."))
print(paste0("Optimal policy:"))
print(vi$policy)
print(paste0("Optimal values:"))
print(round(vi$V, 4))

#=====================================================================
#3. CONTINUOUS-TIME FORMULATION
#=====================================================================
#Convert discrete-time MDP to continuous-time representation.
#Rate matrix: Q = (A - I) / dt
#Discount rate: r = -ln(gamma) / dt
#Reward rate: c = R / dt
#Continuous-time Bellman: V = (r*I - Q_pi)^{-1} c_pi

#Rate matrices Q(a) = (A(a) - I) / dt
Q_rate <- lapply(A, function(Aa) {
  Qa <- (Aa - diag(K)) / dt
  dimnames(Qa) <- list(states, states)
  Qa
})

r_disc <- -log(gamma_disc) / dt     #continuous discount rate
c_rate <- R_mat / dt                 #reward rate (per second)

#Build Q_pi and c_pi under the optimal policy
Q_pi <- matrix(0, K, K, dimnames = list(states, states))
c_pi <- numeric(K); names(c_pi) <- states

for (s in 1:K) {
  a <- vi$policy[states[s]]
  Q_pi[s, ] <- Q_rate[[a]][s, ]
  c_pi[s]   <- c_rate[a, s]
}

#Solve continuous-time Bellman equation
V_ct <- solve(r_disc * diag(K) - Q_pi, c_pi)
names(V_ct) <- states
print(paste0("Continuous-time Bellman values:"))
print(round(V_ct, 4))

#Reward-rate vector under optimal policy (for theoretical curves)
R_pi_vec <- sapply(1:K, function(s) c_rate[vi$policy[states[s]], s])

#=====================================================================
#4. GILLESPIE SIMULATION FUNCTIONS
#=====================================================================

#Gillespie SSA under a fixed deterministic policy
#Simulates continuous-time Markov chain using the Gillespie
#Stochastic Simulation Algorithm. At each event:
#  1. Sample sojourn time tau ~ Exp(total exit rate)
#  2. Sample next state proportional to exit rates
#  3. Accumulate discounted return over the sojourn
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

  #Discounted return
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

#Gillespie with random policy
#At each sojourn, a uniformly random action is selected
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

#Gillespie with real MCDA/TOPSIS policy (companion paper [4]).
#At each transition event, an EMA posterior over regimes is updated
#from the noisy observation z. Five criteria (KPI_gain, Module_gain,
#Drift_gain, Stability, Cost) are computed per repair action from the
#posterior, and RMCDA::apply.TOPSIS ranks the five repair actions.
#A gate forces NoAction when the posterior is confidently Nominal and
#KPI proxy is low. Time-based dwell hysteresis (dwell_sec) prevents
#rapid switching. Weights, action costs, and stability penalties match
#MCDA_intervention.R from the companion HMM paper.
#
#The mcda_action argument is retained for signature compatibility with
#legacy callers but is unused; TOPSIS computes actions from criteria.
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

    #Gate: confidently Nominal + low KPI proxy => NoAction
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

#Legacy shortcut: MDP-optimal policy applied to hard classifications.
#Retained for ablation purposes only; NOT used by the figure notebooks.
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

#=====================================================================
#5. CUMULATIVE REWARD HELPER FUNCTIONS
#=====================================================================

#Cumulative reward at specified time points for a deterministic-policy trajectory
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

#Cumulative reward at time points for MCDA trajectory
#Uses the stored action sequence (tr$actions)
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

#Cumulative reward at time points for POMDP/action-stored trajectory
cum_reward_at_pts_pomdp <- function(tr, t_pts, c_rate, r_disc) {
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

#=====================================================================
#6. STATE-OCCUPANCY AND SOJOURN HELPERS
#=====================================================================

#Bin a trajectory into state indices at given time points
bin_trajectory <- function(tr, t_bins, K) {
  out <- integer(length(t_bins))
  for (b in seq_along(t_bins)) {
    j <- max(which(tr$times <= t_bins[b]))
    out[b] <- tr$states[j]
  }
  out
}

#Fraction of time in Nominal for one trajectory
frac_in_nominal <- function(tr, T_end) {
  occ <- 0
  n_ev <- length(tr$times)
  if (n_ev > 1) {
    for (i in 1:(n_ev - 1)) {
      if (tr$states[i] == 1)
        occ <- occ + (tr$times[i + 1] - tr$times[i])
    }
  }
  dt_last <- T_end - tr$times[n_ev]
  if (tr$states[n_ev] == 1) occ <- occ + dt_last
  occ / T_end
}

#=====================================================================
#7. POLICY DEFINITIONS
#=====================================================================

#No-intervention policy: always NoAction
policy_none <- setNames(rep("NoAction", K), states)

#MCDA policy parameters
#mcda_accuracy and mcda_confusion are provided by SimulationPipeline.R
mcda_ideal_action <- vi$policy
mcda_dwell_time   <- 0.50

#=====================================================================
#8. OBSERVATION MODEL (POMDP)
#=====================================================================
#The observation matrix O_mat[s, z] = P(observe z | true state s)
#Reuses the MCDA confusion matrix as the noisy sensor model.

O_mat <- mcda_confusion
print(paste0("Observation matrix (rows = true state, cols = perceived state):"))
print(round(O_mat, 3))

#=====================================================================
#9. BELIEF POINT GENERATION
#=====================================================================
#Generates a set of belief points on the K-simplex for PBVI.
#Includes: corner beliefs (pure states), edge midpoints,
#face centres (triplets), centroid, and random Dirichlet samples.

generate_belief_points <- function(K, n_random = 185) {
  pts <- list()

  #Corner beliefs (pure states)
  for (i in 1:K) {
    b <- rep(0, K); b[i] <- 1
    pts[[length(pts) + 1]] <- b
  }

  #Edge midpoints
  for (i in 1:(K - 1)) {
    for (j in (i + 1):K) {
      b <- rep(0, K); b[i] <- 0.5; b[j] <- 0.5
      pts[[length(pts) + 1]] <- b
    }
  }

  #Face centres (triplets)
  if (K >= 3) {
    combs <- combn(K, 3)
    for (c_idx in 1:ncol(combs)) {
      b <- rep(0, K); b[combs[, c_idx]] <- 1 / 3
      pts[[length(pts) + 1]] <- b
    }
  }

  #Centroid
  pts[[length(pts) + 1]] <- rep(1 / K, K)

  #Random Dirichlet(1,...,1) samples
  for (r in 1:n_random) {
    g <- rgamma(K, 1, 1)
    pts[[length(pts) + 1]] <- g / sum(g)
  }

  pts
}

#=====================================================================
#10. PBVI SOLVER (Point-Based Value Iteration)
#=====================================================================
#Iteratively computes alpha-vectors that parameterize the value
#function over the belief simplex. For each belief point, finds the
#best action and backs up through the observation model.
#Warm-starts with MDP value function for faster convergence.

pbvi_solve <- function(A_list, O_mat, R_mat, gamma, belief_pts,
                       max_horizon = 2000, tol = 1e-6, V_init = NULL) {
  K <- ncol(O_mat)
  act <- rownames(R_mat)
  M <- length(act)
  n_pts <- length(belief_pts)

  #Stack beliefs into a matrix (n_pts x K)
  B <- do.call(rbind, belief_pts)

  #Initialise alpha vectors: use MDP values if available
  if (!is.null(V_init)) {
    Gamma_mat <- matrix(V_init, nrow = 1, ncol = K)
    Gamma_labels <- act[1]
  } else {
    Gamma_mat <- matrix(0, nrow = 1, ncol = K)
    Gamma_labels <- act[1]
  }

  for (iter in 1:max_horizon) {
    n_alpha <- nrow(Gamma_mat)

    #Precompute backed-up vectors for each (action, observation)
    backed <- vector("list", M)
    for (ai in 1:M) {
      a <- act[ai]
      backed[[ai]] <- vector("list", K)
      for (z in 1:K) {
        #weighted[k, j] = Gamma_mat[j, k] * O_mat[k, z]
        weighted <- t(Gamma_mat) * O_mat[, z]
        #backed[[ai]][[z]] is K x n_alpha
        backed[[ai]][[z]] <- gamma * (A_list[[a]] %*% weighted)
      }
    }

    Gamma_new   <- matrix(NA, n_pts, K)
    Labels_new  <- character(n_pts)

    for (bi in 1:n_pts) {
      b <- B[bi, ]
      best_val <- -Inf; best_alpha <- NULL; best_action <- NULL

      for (ai in 1:M) {
        alpha_a <- R_mat[ai, ]
        for (z in 1:K) {
          scores <- as.numeric(crossprod(backed[[ai]][[z]], b))
          best_gi <- which.max(scores)
          alpha_a <- alpha_a + backed[[ai]][[z]][, best_gi]
        }
        val <- sum(alpha_a * b)
        if (val > best_val) {
          best_val <- val
          best_alpha <- alpha_a
          best_action <- act[ai]
        }
      }
      Gamma_new[bi, ]  <- best_alpha
      Labels_new[bi]   <- best_action
    }

    #Convergence check
    V_new <- rowSums(Gamma_new * B)
    V_old <- apply(B %*% t(Gamma_mat), 1, max)
    max_diff <- max(abs(V_new - V_old))

    #Prune dominated alpha vectors
    keep <- !duplicated(round(Gamma_new, 10))
    Gamma_mat    <- Gamma_new[keep, , drop = FALSE]
    Gamma_labels <- Labels_new[keep]

    if (max_diff < tol) {
      print(paste0("PBVI converged at iteration ", iter,
                   " with max Bellman residual ", signif(max_diff, 3)))
      break
    }
    if (iter %% 100 == 0)
      print(paste0("PBVI iteration ", iter,
                   ": max diff = ", signif(max_diff, 4),
                   ", n_alpha = ", nrow(Gamma_mat)))
  }

  if (iter == max_horizon)
    print(paste0("PBVI reached max horizon (", max_horizon,
                 ") with max diff = ", signif(max_diff, 4)))

  Gamma_list <- lapply(1:nrow(Gamma_mat), function(i) Gamma_mat[i, ])
  list(Gamma = Gamma_list, actions = Gamma_labels, iterations = iter)
}

#=====================================================================
#11. POMDP ACTION SELECTION AND BELIEF UPDATE
#=====================================================================

#Select action that maximizes alpha-vector dot product with belief
pomdp_action <- function(b, Gamma, action_labels) {
  vals <- sapply(Gamma, function(alpha) sum(alpha * b))
  action_labels[which.max(vals)]
}

#Bayesian belief update: predict with transition matrix,
#then condition on observation likelihood
belief_update <- function(b, a, z, A_trans, O_mat) {
  b_pred <- as.numeric(b %*% A_trans)
  b_new  <- b_pred * O_mat[, z]
  s <- sum(b_new)
  if (s > 0) b_new / s else rep(1 / length(b), length(b))
}

#=====================================================================
#12. GILLESPIE SSA UNDER POMDP POLICY
#=====================================================================
#The true state evolves as a hidden Markov chain. At each transition
#the agent receives a noisy observation, updates its belief via
#Bayesian filtering, and selects an action using the PBVI solution.

gillespie_pomdp <- function(Q_list, A_list, pomdp_Gamma, pomdp_actions,
                            O_mat, R_rate, s0, T_end, r_disc) {
  K  <- nrow(Q_list[[1]])
  sn <- rownames(Q_list[[1]])

  cap     <- 2000
  times   <- numeric(cap)
  st_seq  <- integer(cap)
  act_seq <- character(cap)
  idx     <- 1
  t_now   <- 0
  s       <- s0

  #Initial belief: certain about starting state
  b <- rep(0, K); b[s0] <- 1
  current_action <- pomdp_action(b, pomdp_Gamma, pomdp_actions)

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

    #True state transition (hidden from agent)
    s <- sample.int(K, 1, prob = ro)

    #Generate noisy observation
    z <- sample.int(K, 1, prob = O_mat[s, ])

    #Belief update using continuous-time transition matrix
    A_tau <- expm(Q_list[[a]] * tau)
    b <- belief_update(b, a, z, A_tau, O_mat)

    #Select action based on updated belief
    current_action <- pomdp_action(b, pomdp_Gamma, pomdp_actions)

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

  #Compute discounted return using actual actions taken
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

  list(times = times, states = st_seq, actions = act_seq,
       disc_return = disc_ret)
}

#=====================================================================
#13. SOLVE POMDP (PBVI with MDP warm start)
#=====================================================================

print(paste0("Solving POMDP (PBVI)..."))
belief_pts <- generate_belief_points(K, n_random = 185)
print(paste0("Generated ", length(belief_pts), " belief points on the K=4 simplex."))

pomdp <- pbvi_solve(A, O_mat, R_mat, gamma_disc, belief_pts,
                    max_horizon = 2000, tol = 0.001, V_init = vi$V)

print(paste0("PBVI completed in ", pomdp$iterations, " iterations with ",
             length(pomdp$Gamma), " alpha vectors."))

#Print policy at corner beliefs (should match MDP)
print(paste0("POMDP policy at pure-state beliefs:"))
for (s in 1:K) {
  b_corner <- rep(0, K); b_corner[s] <- 1
  a_pomdp <- pomdp_action(b_corner, pomdp$Gamma, pomdp$actions)
  print(paste0("  b = ", states[s], " -> action: ", a_pomdp,
               " (MDP: ", vi$policy[states[s]], ")"))
}

print(paste0("MDP_POMDP_Formulation complete."))
