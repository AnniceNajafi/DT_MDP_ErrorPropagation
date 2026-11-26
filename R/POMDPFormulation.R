#' POMDPFormulation.R
#' ──────────────────────────────────────────────────────────────────────
#' POMDP (Partially Observable MDP) solver and simulation functions.
#' The POMDP operates on observations rather than direct state knowledge,
#' maintaining a belief vector updated via Bayesian filtering.
#'
#' Sources BaseMDPFormulation.R for model definitions (states, actions,
#' transition matrices, reward matrices, discount factors, etc.).
#'
#' Implements:
#'   - Observation model (from MCDA confusion matrix)
#'   - Belief point generation on the K-simplex
#'   - PBVI (Point-Based Value Iteration) solver
#'   - Belief update (Bayesian filtering with continuous-time prediction)
#'   - Gillespie SSA simulator under POMDP policy
#'   - Cumulative reward computation for POMDP trajectories
#'
#' PBVI settings used for the paper (Section 3 / Supplementary S3):
#' 185 random Dirichlet belief points drawn after set.seed(2024) plus the
#' 15 structured points, MDP warm start, fixed budget of 500 iterations
#' (tol = 0.001 is not reached within 500 iterations with the data-driven O; the
#' successive difference is 2.2e-3 at iteration 500, see
#' analyses/Revision_PBVI_Convergence.R and analyses/pbvi_bounds_and_gap.R).
#'
#' Authors@R: person("Annice", "Najafi")
#' ──────────────────────────────────────────────────────────────────────

# Repository-relative paths (see R/paths.R)
if (!exists("repo_path")) local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })

# Source base MDP model (which sources R/SimulationPipeline.R)
source(repo_path("R", "BaseMDPFormulation.R"), local = globalenv())

# =====================================================================
# 1. OBSERVATION MODEL
# =====================================================================
# The observation matrix O_mat[s, z] = P(observe z | true state s)
# Reuses the MCDA confusion matrix as the noisy sensor model.

O_mat <- mcda_confusion
cat("\nObservation matrix (rows = true state, cols = perceived state):\n")
print(round(O_mat, 3))

# =====================================================================
# 2. BELIEF POINT GENERATION
# =====================================================================

generate_belief_points <- function(K, n_random = 185) {
  pts <- list()

  # Corner beliefs (pure states)
  for (i in 1:K) {
    b <- rep(0, K); b[i] <- 1
    pts[[length(pts) + 1]] <- b
  }

  # Edge midpoints
  for (i in 1:(K - 1)) {
    for (j in (i + 1):K) {
      b <- rep(0, K); b[i] <- 0.5; b[j] <- 0.5
      pts[[length(pts) + 1]] <- b
    }
  }

  # Face centres (triplets)
  if (K >= 3) {
    combs <- combn(K, 3)
    for (c_idx in 1:ncol(combs)) {
      b <- rep(0, K); b[combs[, c_idx]] <- 1 / 3
      pts[[length(pts) + 1]] <- b
    }
  }

  # Centroid
  pts[[length(pts) + 1]] <- rep(1 / K, K)

  # Random Dirichlet(1,...,1) samples
  for (r in 1:n_random) {
    g <- rgamma(K, 1, 1)
    pts[[length(pts) + 1]] <- g / sum(g)
  }

  pts
}

# =====================================================================
# 3. PBVI SOLVER (Point-Based Value Iteration)
# =====================================================================

pbvi_solve <- function(A_list, O_mat, R_mat, gamma, belief_pts,
                       max_horizon = 2000, tol = 1e-6, V_init = NULL) {
  K <- ncol(O_mat)
  act <- rownames(R_mat)
  M <- length(act)
  n_pts <- length(belief_pts)

  # Stack beliefs into a matrix (n_pts x K)
  B <- do.call(rbind, belief_pts)

  # Initialise alpha vectors: use MDP values if available
  if (!is.null(V_init)) {
    Gamma_mat <- matrix(V_init, nrow = 1, ncol = K)
    Gamma_labels <- act[1]
  } else {
    Gamma_mat <- matrix(0, nrow = 1, ncol = K)
    Gamma_labels <- act[1]
  }

  for (iter in 1:max_horizon) {
    n_alpha <- nrow(Gamma_mat)

    # Precompute backed-up vectors for each (action, observation)
    backed <- vector("list", M)
    for (ai in 1:M) {
      a <- act[ai]
      backed[[ai]] <- vector("list", K)
      for (z in 1:K) {
        # weighted[k, j] = Gamma_mat[j, k] * O_mat[k, z]
        weighted <- t(Gamma_mat) * O_mat[, z]
        # backed[[ai]][[z]] is K x n_alpha
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

    # Convergence check
    V_new <- rowSums(Gamma_new * B)
    V_old <- apply(B %*% t(Gamma_mat), 1, max)
    max_diff <- max(abs(V_new - V_old))

    # Prune dominated alpha vectors
    keep <- !duplicated(round(Gamma_new, 10))
    Gamma_mat    <- Gamma_new[keep, , drop = FALSE]
    Gamma_labels <- Labels_new[keep]

    if (max_diff < tol) {
      cat("  PBVI converged at iteration", iter,
          "with max Bellman residual", signif(max_diff, 3), "\n")
      break
    }
    if (iter %% 100 == 0)
      cat("  PBVI iteration", iter,
          ": max diff =", signif(max_diff, 4),
          ", n_alpha =", nrow(Gamma_mat), "\n")
  }

  if (iter == max_horizon)
    cat("  PBVI reached max horizon (", max_horizon,
        ") with max diff =", signif(max_diff, 4), "\n")

  Gamma_list <- lapply(1:nrow(Gamma_mat), function(i) Gamma_mat[i, ])
  list(Gamma = Gamma_list, actions = Gamma_labels, iterations = iter)
}

# =====================================================================
# 4. POMDP ACTION SELECTION & BELIEF UPDATE
# =====================================================================

pomdp_action <- function(b, Gamma, action_labels) {
  vals <- sapply(Gamma, function(alpha) sum(alpha * b))
  action_labels[which.max(vals)]
}

belief_update <- function(b, a, z, A_trans, O_mat) {
  b_pred <- as.numeric(b %*% A_trans)
  b_new  <- b_pred * O_mat[, z]
  s <- sum(b_new)
  if (s > 0) b_new / s else rep(1 / length(b), length(b))
}

# =====================================================================
# 5. GILLESPIE SSA UNDER POMDP POLICY
# =====================================================================

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

  # Initial belief: certain about starting state
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

    # True state transition (hidden from agent)
    s <- sample.int(K, 1, prob = ro)

    # Generate noisy observation
    z <- sample.int(K, 1, prob = O_mat[s, ])

    # Belief update using continuous-time transition matrix
    A_tau <- expm(Q_list[[a]] * tau)
    b <- belief_update(b, a, z, A_tau, O_mat)

    # Select action based on updated belief
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

  # Compute discounted return using actual actions taken
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

# =====================================================================
# 6. CUMULATIVE REWARD AT TIME POINTS (POMDP)
# =====================================================================

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

# =====================================================================
# 7. SOLVE POMDP (PBVI with MDP warm start)
# =====================================================================

cat("\n── Solving POMDP (PBVI) ──\n")
set.seed(2024)   # belief points (Dirichlet samples) used for the paper
belief_pts <- generate_belief_points(K, n_random = 185)
cat("  Generated", length(belief_pts), "belief points on the K=4 simplex.\n")

pbvi_horizon <- getOption("dtmdp.pbvi_horizon", 500)
pomdp <- pbvi_solve(A, O_mat, R_mat, gamma_disc, belief_pts,
                    max_horizon = pbvi_horizon, tol = 0.001, V_init = vi$V)

cat("  PBVI completed in", pomdp$iterations, "iterations with",
    length(pomdp$Gamma), "alpha vectors.\n")

# Print policy at corner beliefs (should match MDP)
cat("\nPOMDP policy at pure-state beliefs:\n")
for (s in 1:K) {
  b_corner <- rep(0, K); b_corner[s] <- 1
  a_pomdp <- pomdp_action(b_corner, pomdp$Gamma, pomdp$actions)
  cat("  b =", states[s], "-> action:", a_pomdp,
      "(MDP:", vi$policy[states[s]], ")\n")
}

source(repo_path("R", "sim_persample.R"), local = globalenv())  # per-sample policy evaluation
.dtmdp_core_loaded <- TRUE
cat("\n══ POMDPFormulation complete. ══\n")
