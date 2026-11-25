# rl_core.R
# Standalone extraction of the RL environment, Q-learning, REINFORCE and
# Gillespie evaluation code from the model-free RL notebook
# (notebooks/ModelFreeRL_Significance.Rmd).
# Algorithms are identical to the Rmd; only the swept settings
# (grid delta, learning rate, episodes, training seed) are made arguments.
# Requires model objects (K, M, actions, Q_rate, c_rate, r_disc, O_mat, ...)
# in the global environment (from R/POMDPFormulation.R, or loaded from
# output/rl_sweep/model_objects.rds by the sweep workers).

suppressPackageStartupMessages(library(expm))

# ---------------- RL environment (chunk rl-environment) ----------------
rl_env_new <- function(s0, T_end) {
  b <- rep(0, K); b[s0] <- 1
  list(s = s0, b = b, t_now = 0, T_end = T_end, done = FALSE)
}

rl_env_step <- function(env, action_idx) {
  a <- actions[action_idx]
  Qa <- Q_rate[[a]]
  s <- env$s
  ro <- pmax(Qa[s, ], 0); ro[s] <- 0
  lam <- sum(ro)
  if (lam <= 0 || env$done) {
    dt_rem <- env$T_end - env$t_now
    rew <- 0
    if (dt_rem > 0 && r_disc > 0)
      rew <- (c_rate[a, s] / r_disc) * exp(-r_disc * env$t_now) *
             (1 - exp(-r_disc * dt_rem))
    env$done <- TRUE
    return(list(env = env, obs = s, reward = rew, done = TRUE))
  }
  tau <- rexp(1, lam)
  if (env$t_now + tau > env$T_end) {
    dt_rem <- env$T_end - env$t_now
    rew <- 0
    if (dt_rem > 0 && r_disc > 0)
      rew <- (c_rate[a, s] / r_disc) * exp(-r_disc * env$t_now) *
             (1 - exp(-r_disc * dt_rem))
    env$t_now <- env$T_end
    env$done <- TRUE
    return(list(env = env, obs = s, reward = rew, done = TRUE))
  }
  rew <- 0
  if (tau > 0 && r_disc > 0)
    rew <- (c_rate[a, s] / r_disc) * exp(-r_disc * env$t_now) *
           (1 - exp(-r_disc * tau))
  env$t_now <- env$t_now + tau
  s_new <- sample.int(K, 1, prob = ro)
  env$s <- s_new
  z <- sample.int(K, 1, prob = O_mat[s_new, ])
  A_tau <- expm(Qa * tau)
  b_pred <- as.numeric(env$b %*% A_tau)
  b_new <- b_pred * O_mat[, z]
  b_sum <- sum(b_new)
  if (b_sum > 0) env$b <- b_new / b_sum else env$b <- rep(1 / K, K)
  list(env = env, obs = z, reward = rew, done = FALSE)
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
train_qlearning <- function(delta = 0.1, alpha_ql = 0.1, n_episodes_ql = 5000,
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
      total_return <- total_return + r
      if (!step$done) {
        g_next <- nearest_grid(env$b, belief_grid)
        td_target <- r + gamma_ql * max(Q_table[g_next, ])
        Q_table[g, a_idx] <- Q_table[g, a_idx] + alpha_ql * (td_target - Q_table[g, a_idx])
        g <- g_next
      } else {
        Q_table[g, a_idx] <- Q_table[g, a_idx] + alpha_ql * (r - Q_table[g, a_idx])
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

train_reinforce <- function(lr_re = 0.001, n_episodes_re = 10000, seed = 456,
                            T_rl_re = 15, baseline_decay = 0.99) {
  set.seed(seed)
  theta_re <- matrix(rnorm(M * n_features, sd = 0.01), nrow = M, ncol = n_features)
  baseline <- 0
  episode_returns_re <- numeric(n_episodes_re)
  for (ep in 1:n_episodes_re) {
    s0 <- sample.int(K, 1, prob = c(0.4, 0.2, 0.2, 0.2))
    env <- rl_env_new(s0, T_rl_re)
    log_grads <- list(); rewards <- numeric(0); step_count <- 0
    while (!env$done) {
      phi <- belief_features(env$b)
      probs <- softmax_policy(theta_re, phi)
      a_idx <- sample.int(M, 1, prob = probs)
      step_count <- step_count + 1
      log_grads[[step_count]] <- softmax_log_grad(theta_re, phi, a_idx)
      step <- rl_env_step(env, a_idx)
      env <- step$env
      rewards <- c(rewards, step$reward)
    }
    G <- sum(rewards)
    episode_returns_re[ep] <- G
    baseline <- baseline_decay * baseline + (1 - baseline_decay) * G
    grad_sum <- matrix(0, nrow = M, ncol = n_features)
    for (t in seq_along(log_grads)) {
      G_t <- sum(rewards[t:length(rewards)])
      grad_sum <- grad_sum + (G_t - baseline) * log_grads[[t]]
    }
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

# Evaluate with common random numbers: trajectory n uses set.seed(eval_seed0 + n)
evaluate_policy <- function(select_action, N_eval = 500, s0 = 1, T_end = 15,
                            eval_seed0 = 900000) {
  vapply(1:N_eval, function(n) {
    set.seed(eval_seed0 + n)
    gillespie_belief_policy(Q_rate, select_action, O_mat, c_rate,
                            s0, T_end, r_disc)$disc_return
  }, numeric(1))
}
