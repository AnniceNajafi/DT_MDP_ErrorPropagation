#' kstep_policy.R
#' ──────────────────────────────────────────────────────────────────────
#' k-step consecutive-label heuristic baseline (Reviewer comment 5 of the
#' first review round): apply NoAction until the observed HMM label has
#' been non-Nominal for k consecutive observations, then apply the
#' canonical repair for that label. Function definitions only; the k
#' sweep is in analyses/Revision_KStepHeuristicBaseline.R.
#'
#' Requires the model objects of R/BaseMDPFormulation.R in the global
#' environment.
#'
#' Authors@R: person("Annice", "Najafi")
#' ──────────────────────────────────────────────────────────────────────

# =====================================================================
# 1. CANONICAL REPAIR MAP (from MDP-optimal policy at corner beliefs)
# =====================================================================
# For each observed non-Nominal label, the physically targeted repair.

canonical_repair <- c(
  "Nominal"     = "NoAction",
  "SensorNoisy" = "DwSensorA",
  "DynamicsOff" = "ReidentPlant",
  "Drift"       = "BiasCorrect"
)

# =====================================================================
# 2. GILLESPIE SSA UNDER k-STEP HEURISTIC POLICY
# =====================================================================
# Mirrors the structure of gillespie_mcda() from BaseMDPFormulation.R
# so that the act_seq is stored for cum_reward_at_pts_mcda() reuse.

gillespie_kstep <- function(Q_list, O_mat, R_rate, s0, T_end, r_disc,
                            k_threshold = 2, dwell_sec = 0.5) {
  K  <- nrow(Q_list[[1]])
  sn <- rownames(Q_list[[1]])

  cap     <- 2000
  times   <- numeric(cap)
  st_seq  <- integer(cap)
  act_seq <- character(cap)
  idx     <- 1
  t_now   <- 0
  s       <- s0

  current_action   <- "NoAction"
  consecutive_label <- 1L           # index of last non-Nominal label
  consecutive_count <- 0L
  last_switch_time  <- 0

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

    # True state transition (hidden from the operator)
    s <- sample.int(K, 1, prob = ro)

    # Noisy HMM/Viterbi observation at the event
    z <- sample.int(K, 1, prob = O_mat[s, ])

    # k-step consecutive non-Nominal counter
    if (z == 1) {
      # Nominal observation resets the counter and returns to NoAction
      consecutive_count <- 0L
      desired_action    <- "NoAction"
    } else {
      if (z == consecutive_label) {
        consecutive_count <- consecutive_count + 1L
      } else {
        consecutive_label <- z
        consecutive_count <- 1L
      }
      if (consecutive_count >= k_threshold) {
        desired_action <- canonical_repair[sn[z]]
      } else {
        desired_action <- "NoAction"
      }
    }

    # Dwell hysteresis to prevent rapid flip-flopping
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

  list(times = times, states = st_seq, actions = act_seq,
       disc_return = disc_ret)
}
