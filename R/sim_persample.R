#' sim_persample.R
#' ----------------------------------------------------------------------
#' Per-sample policy evaluation for the decision layer.
#'
#' The HMM decoder labels every residual sample (dt = 20 ms), and PBVI
#' solves the discrete-time POMDP in which an observation arrives at every
#' step. Policies that act on decoder output (POMDP, myopic belief, QMDP,
#' RL, k-Step, MCDA) are therefore evaluated here step by step:
#'
#'   at step k = 0, ..., H-1  (t = k dt, H = T_end / dt):
#'     a_k   = controller action (from the true state for the MDP oracle,
#'             otherwise from the observation history)
#'     G    += gamma^k R(a_k, s_k)
#'     s_k+1 ~ A^(a_k)[s_k, ]
#'     z_k+1 ~ O[s_k+1, ]          (decoder label of the next sample)
#'     controller update with (a_k, z_k+1)
#'
#' The return G is the discounted sum the MDP and the POMDP optimize, so
#' simulated means are directly comparable with V*(s0) and with the PBVI
#' value at the initial belief. The event-driven Gillespie simulators of
#' R/BaseMDPFormulation.R remain in use for the continuous-time checks of
#' the fully observed MDP (analyses/mdp_validation.R, figures/fig01).
#'
#' Requires K, states, actions, R_mat, gamma_disc (R/BaseMDPFormulation.R).
#' ----------------------------------------------------------------------

#' Simulate one trajectory.
#' @param ctl  controller: list(init(s0, O_true) -> state env,
#'             act(st, s) -> action name, update(st, a, z, t) -> state env)
#' @param A_list named list of K x K transition matrices (one per action)
#' @param O_true observation matrix that generates the decoder labels
#' @param t_pts optional time points for the cumulative discounted reward
#' @return list(times, states, actions) at the steps where the state or the
#'   action changes (the format of the Gillespie trajectories, so the metric
#'   helpers apply unchanged), disc_return, and cum (at t_pts) if requested.
sim_ps <- function(ctl, A_list, O_true, s0 = 1L, T_end = 20, dt = 0.02,
                   gamma = gamma_disc, R = R_mat, t_pts = NULL,
                   record_belief = FALSE) {
  H <- as.integer(round(T_end / dt))
  st <- ctl$init(s0, O_true)
  s <- s0; G <- 0; disc <- 1
  s_seq <- integer(H); a_seq <- character(H); cumr <- numeric(H)
  bel <- if (record_belief) matrix(NA_real_, H, K) else NULL
  for (k in seq_len(H)) {
    if (record_belief && !is.null(st$b)) bel[k, ] <- st$b
    a <- ctl$act(st, s)
    s_seq[k] <- s; a_seq[k] <- a
    G <- G + disc * R[a, s]; cumr[k] <- G; disc <- disc * gamma
    s <- sample.int(K, 1L, prob = A_list[[a]][s, ])
    z <- sample.int(K, 1L, prob = O_true[s, ])
    st <- ctl$update(st, a, z, k * dt)
  }
  chg <- c(1L, which(s_seq[-1L] != s_seq[-H] | a_seq[-1L] != a_seq[-H]) + 1L)
  out <- list(times = (chg - 1L) * dt, states = s_seq[chg], actions = a_seq[chg],
              disc_return = G)
  if (!is.null(t_pts)) {
    # reward accrued by time t: steps k with k*dt < t
    n_k <- pmin(H, ceiling(t_pts / dt - 1e-9))
    out$cum <- ifelse(n_k > 0, cumr[pmax(n_k, 1L)], 0)
  }
  if (record_belief) { out$belief <- bel; out$step_states <- s_seq; out$step_actions <- a_seq }
  out
}

# ---------------------------------------------------------------------
# Controllers
# ---------------------------------------------------------------------

#' Fully observed state-feedback policy (MDP oracle, No Intervention).
ctl_state <- function(policy) {
  list(init = function(s0, O_true) new.env(),
       act = function(st, s) unname(policy[states[s]]),
       update = function(st, a, z, t) st)
}

#' Belief-feedback policy: Bayesian filter b' ~ O_bel[, z] * (b A_bel^(a)),
#' started at the known initial state; `choose(b)` returns an action name.
ctl_belief <- function(choose, A_bel = A, O_bel = O_mat) {
  list(init = function(s0, O_true) { e <- new.env(); e$b <- replace(numeric(K), s0, 1); e },
       act = function(st, s) choose(st$b),
       update = function(st, a, z, t) {
         bp <- as.numeric(st$b %*% A_bel[[a]]) * O_bel[, z]
         sb <- sum(bp)
         st$b <- if (sb > 0) bp / sb else rep(1 / K, K)
         st
       })
}

#' k-Step consecutive-label heuristic (R/kstep_policy.R, per decoder label):
#' NoAction until the label has been the same non-Nominal regime for k
#' consecutive samples, then that regime's canonical repair; a Nominal
#' label resets; switches are separated by at least dwell_sec.
ctl_kstep <- function(k_threshold = 1L, dwell_sec = 0.5,
                      repair = c(Nominal = "NoAction", SensorNoisy = "DwSensorA",
                                 DynamicsOff = "ReidentPlant", Drift = "BiasCorrect")) {
  list(init = function(s0, O_true) {
         e <- new.env(); e$a <- "NoAction"; e$lab <- 1L; e$cnt <- 0L; e$t_sw <- 0; e
       },
       act = function(st, s) st$a,
       update = function(st, a, z, t) {
         if (z == 1L) { st$cnt <- 0L; want <- "NoAction" } else {
           if (z == st$lab) st$cnt <- st$cnt + 1L else { st$lab <- z; st$cnt <- 1L }
           want <- if (st$cnt >= k_threshold) unname(repair[states[z]]) else "NoAction"
         }
         if (want != st$a && (t - st$t_sw) >= dwell_sec - 1e-9) { st$a <- want; st$t_sw <- t }
         st
       })
}

#' TOPSIS decision of the MCDA baseline for a smoothed label distribution
#' (criteria of Supplementary Section S5.6; identical to gillespie_mcda()).
mcda_topsis_action <- function(posterior) {
  topsis_actions <- c("DwSensorA", "DwSensorB", "IncFilter", "ReidentPlant", "BiasCorrect")
  stability_penalty <- c(0.08, 0.08, 0.05, 0.15, 0.12)
  action_cost       <- c(0.10, 0.10, 0.15, 0.40, 0.25)
  topsis_weights    <- c(0.35, 0.25, 0.15, 0.10, 0.15)
  pN <- posterior[1]; pS <- posterior[2]; pD <- posterior[3]; pR <- posterior[4]
  sensor_n <- pS; plant_n <- pD; drift_n <- pR; kpi_n <- 1 - pN; noise_n <- sensor_n
  kpi_gain <- c(pS * kpi_n * sensor_n * 1.2, pS * kpi_n * sensor_n * 0.6,
                pS * kpi_n * noise_n, pD * kpi_n * plant_n,
                (pR + 0.5 * pD) * kpi_n * (drift_n + 0.3 * plant_n))
  module_gain <- c(pS * sensor_n * 1.2, pS * sensor_n * 0.6, pS * sensor_n * noise_n,
                   pD * plant_n, (pR + 0.5 * pD) * (drift_n + 0.3 * plant_n))
  drift_gain <- c(pS * drift_n * 1.1, pS * drift_n * 0.7, pS * drift_n * 0.6,
                  pD * drift_n * 0.4, (pR + 0.2 * pS) * drift_n)
  dm <- cbind(KPI_gain = kpi_gain, Module_gain = module_gain, Drift_gain = drift_gain,
              Stability = -stability_penalty, Cost = -action_cost)
  rownames(dm) <- topsis_actions
  sc <- suppressWarnings(as.numeric(RMCDA::apply.TOPSIS(dm, w = topsis_weights)))
  desired <- if (all(is.na(sc)) || all(sc == 0)) "NoAction" else topsis_actions[which.max(sc)]
  if (kpi_n < 0.25 && pN > 0.75) desired <- "NoAction"
  desired
}

#' MCDA baseline: exponential smoothing (alpha_ema) of the one-hot decoder
#' labels, TOPSIS ranking, dwell-time hysteresis (as gillespie_mcda()).
ctl_mcda <- function(dwell_sec = 0.5, alpha_ema = 0.3) {
  list(init = function(s0, O_true) {
         e <- new.env(); e$p <- numeric(K)
         e$p[sample.int(K, 1L, prob = O_true[s0, ])] <- 1
         e$a <- "NoAction"; e$t_sw <- 0; e
       },
       act = function(st, s) st$a,
       update = function(st, a, z, t) {
         ind <- numeric(K); ind[z] <- 1
         st$p <- alpha_ema * ind + (1 - alpha_ema) * st$p; st$p <- st$p / sum(st$p)
         want <- mcda_topsis_action(st$p)
         if (want != st$a && (t - st$t_sw) >= dwell_sec - 1e-9) { st$a <- want; st$t_sw <- t }
         st
       })
}

#' Greedy action of an alpha-vector set (PBVI) at belief b.
alpha_choose <- function(Gamma, labels) {
  G <- do.call(rbind, Gamma)
  function(b) labels[which.max(G %*% b)]
}
