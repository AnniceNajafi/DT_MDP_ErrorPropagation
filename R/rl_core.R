# rl_core.R
# Standalone extraction of the RL environment, Q-learning, REINFORCE and
# evaluation code from the model-free RL notebook
# (notebooks/ModelFreeRL_Significance.Rmd).
# Algorithms are identical to the Rmd; only the swept settings
# (grid delta, learning rate, episodes, training seed) are made arguments.
# Requires model objects (K, M, actions, A, R_mat, gamma_disc, O_mat, ...)
# in the global environment (from R/POMDPFormulation.R, or loaded from
# output/rl_sweep/model_objects.rds by the sweep workers).

suppressPackageStartupMessages(library(expm))

# ---------------- RL environment: per-sample POMDP (dt = 20 ms) ----------------
# One step = one decoder sample: reward R(a, s), transition s' ~ A^(a)[s, ],
# label z ~ O[s', ], Bayesian belief update with A^(a) and O (the model PBVI
# solves; see R/sim_persample.R). An episode of T_end seconds has
# T_end / dt steps. `reward` is R(a, s); `reward_disc` is gamma^k R(a, s),
# discounted to the start of the episode.
rl_env_new <- function(s0, T_end, dt = 0.02) {
  b <- rep(0, K); b[s0] <- 1
  list(s = s0, b = b, k = 0L, H = as.integer(round(T_end / dt)), done = FALSE)
}

rl_env_step <- function(env, action_idx) {
  a <- actions[action_idx]
  s <- env$s
  rew <- R_mat[a, s]
  rew_disc <- gamma_disc^env$k * rew
  s_new <- sample.int(K, 1, prob = A[[a]][s, ])
  z <- sample.int(K, 1, prob = O_mat[s_new, ])
  b_new <- as.numeric(env$b %*% A[[a]]) * O_mat[, z]
  b_sum <- sum(b_new)
  env$b <- if (b_sum > 0) b_new / b_sum else rep(1 / K, K)
  env$s <- s_new; env$k <- env$k + 1L
  env$done <- env$k >= env$H
  list(env = env, obs = z, reward = rew, reward_disc = rew_disc, done = env$done)
}

# ---------------- Belief grid (chunk belief-grid) ----------------
generate_simplex_grid <- function(K, delta = 0.1) {
  n <- round(1 / delta)
  grid_pts <- list()
  gen <- function(remaining, depth, current) {
    if (depth == K) {
      current[K] <- remaining
      grid_pts[[length(grid_pts) + 1]] <<- current / n
      return()
    }
    for (v in 0:remaining) {
      current[depth] <- v
      gen(remaining - v, depth + 1, current)
    }
  }
  gen(n, 1, numeric(K))
  do.call(rbind, grid_pts)
}

nearest_grid <- function(b, grid) {
  dists <- rowSums(abs(sweep(grid, 2, b)))
  which.min(dists)
}

# ---------------- Q-learning training (chunk qlearning-train) ----------------
train_qlearning <- function(delta = 0.1, alpha_ql = 0.01, n_episodes_ql = 5000,
                            seed = 123, gamma_ql = 0.99, eps_start = 1.0,
                            eps_end = 0.05, T_rl = 15) {
  belief_grid <- generate_simplex_grid(K, delta = delta)
  n_grid <- nrow(belief_grid)
  eps_schedule <- exp(seq(log(eps_start), log(eps_end), length.out = n_episodes_ql))
  Q_table <- matrix(0, nrow = n_grid, ncol = M, dimnames = list(NULL, actions))
  episode_returns_ql <- numeric(n_episodes_ql)
  set.seed(seed)
  for (ep in 1:n_episodes_ql) {
    s0 <- sample.int(K, 1, prob = c(0.4, 0.2, 0.2, 0.2))
    env <- rl_env_new(s0, T_rl)
    eps <- eps_schedule[ep]
    total_return <- 0
    g <- nearest_grid(env$b, belief_grid)
    while (!env$done) {
      if (runif(1) < eps) a_idx <- sample.int(M, 1) else a_idx <- which.max(Q_table[g, ])
      step <- rl_env_step(env, a_idx)
      env <- step$env
      r <- step$reward
      total_return <- total_return + step$reward_disc
      if (!step$done) {
        g_next <- nearest_grid(env$b, belief_grid)
        td_target <- r + gamma_ql * max(Q_table[g_next, ])
        Q_table[g, a_idx] <- Q_table[g, a_idx] + alpha_ql * (td_target - Q_table[g, a_idx])
        g <- g_next
      } else {
        Q_table[g, a_idx] <- Q_table[g, a_idx] + alpha_ql * (r - Q_table[g, a_idx])  # terminal step
      }
    }
    episode_returns_ql[ep] <- total_return
  }
  list(Q_table = Q_table, belief_grid = belief_grid, returns = episode_returns_ql,
       visited = sum(rowSums(Q_table != 0) > 0), n_grid = n_grid)
}

# ---------------- REINFORCE (chunks reinforce-features / reinforce-train) ----------------
belief_features <- function(b) {
  phi <- b
  for (i in 1:3) for (j in (i + 1):4) phi <- c(phi, b[i] * b[j])
  ent <- -sum(b * log(b + 1e-10))
  c(phi, ent)
}
n_features <- 11

softmax_policy <- function(theta, phi) {
  logits <- theta %*% phi
  logits <- logits - max(logits)
  probs <- exp(logits) / sum(exp(logits))
  as.numeric(probs)
}

softmax_log_grad <- function(theta, phi, a_idx) {
  probs <- softmax_policy(theta, phi)
  grad <- matrix(0, nrow = M, ncol = n_features)
  for (a in 1:M) {
    indicator <- if (a == a_idx) 1 else 0
    grad[a, ] <- (indicator - probs[a]) * phi
  }
  grad
}

train_reinforce <- function(lr_re = 1e-4, n_episodes_re = 10000, seed = 456,
                            T_rl_re = 15, baseline_decay = 0.99) {
  set.seed(seed)
  theta_re <- matrix(rnorm(M * n_features, sd = 0.01), nrow = M, ncol = n_features)
  baseline <- 0
  episode_returns_re <- numeric(n_episodes_re)
  for (ep in 1:n_episodes_re) {
    s0 <- sample.int(K, 1, prob = c(0.4, 0.2, 0.2, 0.2))
    env <- rl_env_new(s0, T_rl_re)
    H <- env$H
    Phi <- matrix(0, H, n_features); P <- matrix(0, H, M); a_taken <- integer(H)
    rewards <- numeric(H); step_count <- 0
    while (!env$done) {
      phi <- belief_features(env$b)
      probs <- softmax_policy(theta_re, phi)
      a_idx <- sample.int(M, 1, prob = probs)
      step_count <- step_count + 1
      Phi[step_count, ] <- phi; P[step_count, ] <- probs; a_taken[step_count] <- a_idx
      step <- rl_env_step(env, a_idx)
      env <- step$env
      rewards[step_count] <- step$reward_disc
    }
    n <- step_count
    G <- sum(rewards[1:n])
    episode_returns_re[ep] <- G
    baseline <- baseline_decay * baseline + (1 - baseline_decay) * G
    # sum_t (G_t - baseline) * grad log pi(a_t | b_t), G_t = sum_{j >= t} rewards_j
    G_t <- rev(cumsum(rev(rewards[1:n])))
    onehot <- matrix(0, n, M); onehot[cbind(1:n, a_taken[1:n])] <- 1
    grad_sum <- t((onehot - P[1:n, , drop = FALSE]) * (G_t - baseline)) %*% Phi[1:n, , drop = FALSE]
    theta_re <- theta_re + lr_re * grad_sum
  }
  list(theta = theta_re, returns = episode_returns_re)
}

# ---------------- Gillespie evaluation (chunk gillespie-learned) ----------------
# Generic greedy belief-policy simulator; `select_action(b)` returns an action
# name. Body is identical to gillespie_qlearning / gillespie_reinforce /
# gillespie_pomdp except for the action-selection call.
gillespie_belief_policy <- function(Q_list, select_action, O_mat, R_rate,
                                    s0, T_end, r_disc) {
  K <- nrow(Q_list[[1]]); cap <- 2000
  times <- numeric(cap); st_seq <- integer(cap); act_seq <- character(cap)
  idx <- 1; t_now <- 0; s <- s0
  b <- rep(0, K); b[s0] <- 1
  current_action <- select_action(b)
  times[1] <- 0; st_seq[1] <- s; act_seq[1] <- current_action
  while (t_now < T_end) {
    a <- current_action; Qa <- Q_list[[a]]
    ro <- pmax(Qa[s, ], 0); ro[s] <- 0; lam <- sum(ro)
    if (lam <= 0) break
    tau <- rexp(1, lam); t_now <- t_now + tau
    if (t_now > T_end) break
    s <- sample.int(K, 1, prob = ro)
    z <- sample.int(K, 1, prob = O_mat[s, ])
    A_tau <- expm(Qa * tau)
    b_pred <- as.numeric(b %*% A_tau)
    b_new <- b_pred * O_mat[, z]
    b_sum <- sum(b_new)
    if (b_sum > 0) b <- b_new / b_sum else b <- rep(1 / K, K)
    current_action <- select_action(b)
    idx <- idx + 1
    if (idx > length(times)) {
      times <- c(times, numeric(cap)); st_seq <- c(st_seq, integer(cap))
      act_seq <- c(act_seq, character(cap))
    }
    times[idx] <- t_now; st_seq[idx] <- s; act_seq[idx] <- current_action
  }
  times <- times[1:idx]; st_seq <- st_seq[1:idx]; act_seq <- act_seq[1:idx]
  disc_ret <- 0
  for (i in seq_along(times)) {
    a_i <- act_seq[i]; c_i <- R_rate[a_i, st_seq[i]]
    t_start <- times[i]; t_end_i <- if (i < idx) times[i + 1] else T_end
    dt_seg <- t_end_i - t_start
    if (dt_seg > 0 && r_disc > 0)
      disc_ret <- disc_ret + (c_i / r_disc) * exp(-r_disc * t_start) *
                  (1 - exp(-r_disc * dt_seg))
  }
  list(times = times, states = st_seq, actions = act_seq, disc_return = disc_ret)
}

ql_selector <- function(Q_table, belief_grid)
  function(b) actions[which.max(Q_table[nearest_grid(b, belief_grid), ])]
re_selector <- function(theta)
  function(b) actions[which.max(softmax_policy(theta, belief_features(b)))]
pbvi_selector <- function(Gamma, labels)
  function(b) {
    vals <- sapply(Gamma, function(alpha) sum(alpha * b))
    labels[which.max(vals)]
  }

# Evaluate with common random numbers: trajectory n uses set.seed(eval_seed0 + n).
# Per-sample evaluation (R/sim_persample.R), identical to the evaluation of
# PBVI and the other belief policies.
evaluate_policy <- function(select_action, N_eval = 500, s0 = 1, T_end = 15,
                            eval_seed0 = 900000) {
  ctl <- ctl_belief(select_action, A, O_mat)
  vapply(1:N_eval, function(n) {
    set.seed(eval_seed0 + n)
    sim_ps(ctl, A, O_mat, s0 = s0, T_end = T_end, gamma = gamma_disc, R = R_mat)$disc_return
  }, numeric(1))
}
