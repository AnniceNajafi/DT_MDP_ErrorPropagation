# 02_run_sweep.R -- RL sensitivity sweep (Q-learning, REINFORCE) vs PBVI
local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
out_dir <- dirname(out_path("rl_sweep", "x"))     # output/rl_sweep
rl_core <- repo_path("R", "rl_core.R")
library(parallel)

N_EVAL <- nq(500, 50); T_EVAL <- 15; S0_EVAL <- 1; EVAL_SEED0 <- 900000      # primary (Rmd setting)
N_EVAL2 <- nq(1000, 50); T_EVAL2 <- 20; EVAL_SEED0_2 <- 800000          # secondary (paper setting)
N_SEEDS <- nq(5L, 1L)

# Full factorial design of the paper (130 runs); quick mode trains only the
# default configurations with reduced budgets.
ql_grid <- if (QUICK) expand.grid(delta = 0.1, lr = 0.01, episodes = 500, rep = 1:N_SEEDS) else
  expand.grid(delta = c(0.2, 0.1, 0.05), lr = c(0.005, 0.01, 0.02),
              episodes = c(5000, 20000), rep = 1:N_SEEDS)
ql_tasks <- lapply(seq_len(nrow(ql_grid)), function(i)
  list(algo = "QL", delta = ql_grid$delta[i], lr = ql_grid$lr[i],
       episodes = ql_grid$episodes[i], seed = 122 + ql_grid$rep[i], rep = ql_grid$rep[i]))
re_grid <- if (QUICK) expand.grid(lr = 1e-4, episodes = 1000, rep = 1:N_SEEDS) else
  expand.grid(lr = c(1e-5, 3e-5, 1e-4, 3e-4),
              episodes = c(10000, 30000), rep = 1:N_SEEDS)
re_tasks <- lapply(seq_len(nrow(re_grid)), function(i)
  list(algo = "REINFORCE", delta = NA, lr = re_grid$lr[i],
       episodes = re_grid$episodes[i], seed = 455 + re_grid$rep[i], rep = re_grid$rep[i]))
tasks <- c(re_tasks[order(-sapply(re_tasks, `[[`, "episodes"))],
           ql_tasks[order(-sapply(ql_tasks, `[[`, "episodes"))])
cat("Total training runs:", length(tasks), "\n")

worker_init <- function(out_dir, rl_core) {
  mo <- readRDS(file.path(out_dir, "model_objects.rds")); list2env(mo, globalenv())
  source(file.path(dirname(rl_core), "sim_persample.R"))  # sim_ps, ctl_belief
  source(rl_core)
  NULL
}

run_task <- function(tk, N_EVAL, T_EVAL, S0_EVAL, EVAL_SEED0, N_EVAL2, T_EVAL2, EVAL_SEED0_2) {
  t0 <- proc.time()[["elapsed"]]
  if (tk$algo == "QL") {
    fit <- train_qlearning(delta = tk$delta, alpha_ql = tk$lr,
                           n_episodes_ql = tk$episodes, seed = tk$seed)
    sel <- ql_selector(fit$Q_table, fit$belief_grid)
    extra <- list(n_grid = fit$n_grid, visited = fit$visited)
  } else {
    fit <- train_reinforce(lr_re = tk$lr, n_episodes_re = tk$episodes, seed = tk$seed)
    sel <- re_selector(fit$theta)
    extra <- list(n_grid = NA, visited = NA)
  }
  train_time <- proc.time()[["elapsed"]] - t0
  ev  <- evaluate_policy(sel, N_eval = N_EVAL, s0 = S0_EVAL, T_end = T_EVAL, eval_seed0 = EVAL_SEED0)
  ev2 <- evaluate_policy(sel, N_eval = N_EVAL2, s0 = S0_EVAL, T_end = T_EVAL2, eval_seed0 = EVAL_SEED0_2)
  # greedy action at the four corner beliefs
  corner <- sapply(1:K, function(s) { b <- rep(0, K); b[s] <- 1; sel(b) })
  R <- fit$returns; n <- length(R); k <- floor(0.1 * n)
  c(tk, extra, list(train_time = train_time, eval = ev, eval2 = ev2,
                    corner = paste(corner, collapse = "/"),
                    last10 = mean(R[(n - k + 1):n]), prev10 = mean(R[(n - 2 * k + 1):(n - k)]),
                    returns = R))
}

t_start <- Sys.time()
n_workers <- as.integer(Sys.getenv("DTMDP_CORES", min(28, max(1, detectCores() - 2))))
cl <- makeCluster(n_workers)
clusterCall(cl, worker_init, out_dir, rl_core)
res <- clusterApplyLB(cl, tasks, run_task, N_EVAL, T_EVAL, S0_EVAL, EVAL_SEED0,
                      N_EVAL2, T_EVAL2, EVAL_SEED0_2)
# PBVI reference on identical evaluation seeds
pbvi_ev <- clusterCall(cl[1], function(...) {
  sel <- pbvi_selector(pomdp$Gamma, pomdp$actions)
  list(ev = evaluate_policy(sel, N_eval = ..1, s0 = ..2, T_end = ..3, eval_seed0 = ..4),
       ev2 = evaluate_policy(sel, N_eval = ..5, s0 = ..2, T_end = ..6, eval_seed0 = ..7),
       corner = paste(sapply(1:K, function(s) { b <- rep(0, K); b[s] <- 1; sel(b) }), collapse = "/"))
}, N_EVAL, S0_EVAL, T_EVAL, EVAL_SEED0, N_EVAL2, T_EVAL2, EVAL_SEED0_2)[[1]]
stopCluster(cl)
cat("Sweep wall-clock:", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 2), "min\n")
saveRDS(list(res = res, pbvi = pbvi_ev), file.path(out_dir, "sweep_results_raw.rds"))

# ---------------- Tables ----------------
mo <- readRDS(file.path(out_dir, "model_objects.rds"))
pb <- pbvi_ev$ev; pb2 <- pbvi_ev$ev2
per_run <- do.call(rbind, lapply(res, function(r) {
  d <- pb - r$eval; d2 <- pb2 - r$eval2
  data.frame(algo = r$algo, delta = r$delta, n_grid = r$n_grid, grid_cells_visited = r$visited,
             lr = r$lr, episodes = r$episodes, rep = r$rep, train_seed = r$seed,
             train_time_s = round(r$train_time, 2),
             eval_mean_T15 = mean(r$eval), eval_sd_T15 = sd(r$eval),
             pbvi_minus_agent_T15 = mean(d), paired_se_T15 = sd(d) / sqrt(length(d)),
             eval_mean_T20 = mean(r$eval2), eval_sd_T20 = sd(r$eval2),
             pbvi_minus_agent_T20 = mean(d2), paired_se_T20 = sd(d2) / sqrt(length(d2)),
             frac_traj_identical_to_pbvi_T15 = mean(abs(d) < 1e-9),
             train_last10pct = r$last10, train_prev10pct = r$prev10,
             train_last_minus_prev = r$last10 - r$prev10,
             corner_actions = r$corner)
}))
pbvi_row <- data.frame(policy = "PBVI", eval_mean_T15 = mean(pb), eval_sd_T15 = sd(pb),
                       eval_mean_T20 = mean(pb2), eval_sd_T20 = sd(pb2),
                       n_alpha = length(mo$pomdp$Gamma), iterations = mo$pomdp$iterations,
                       solve_time_s = round(mo$pbvi_time, 1), corner_actions = pbvi_ev$corner)
write.csv(per_run, file.path(out_dir, "per_run_results.csv"), row.names = FALSE)
write.csv(pbvi_row, file.path(out_dir, "pbvi_reference.csv"), row.names = FALSE)

key <- interaction(per_run$algo, per_run$delta, per_run$lr, per_run$episodes, drop = TRUE)
summ <- do.call(rbind, lapply(split(per_run, key), function(d) {
  b <- which.max(d$eval_mean_T15)
  data.frame(algo = d$algo[1], delta = d$delta[1], n_grid = d$n_grid[1], lr = d$lr[1],
             episodes = d$episodes[1], n_seeds = nrow(d),
             mean_eval_T15 = mean(d$eval_mean_T15), sd_across_seeds_T15 = sd(d$eval_mean_T15),
             min_seed_T15 = min(d$eval_mean_T15), best_seed_T15 = max(d$eval_mean_T15),
             best_seed_id = d$train_seed[b],
             mean_gap_to_pbvi_T15 = mean(d$pbvi_minus_agent_T15),
             mean_eval_T20 = mean(d$eval_mean_T20), sd_across_seeds_T20 = sd(d$eval_mean_T20),
             best_seed_T20 = max(d$eval_mean_T20),
             mean_train_time_s = mean(d$train_time_s),
             mean_last_minus_prev10 = mean(d$train_last_minus_prev))
}))
summ <- summ[order(summ$algo, summ$delta, summ$lr, summ$episodes), ]
write.csv(summ, file.path(out_dir, "config_summary.csv"), row.names = FALSE)

# eval returns per trajectory (long), for paired re-analysis
long <- do.call(rbind, lapply(seq_along(res), function(i) {
  r <- res[[i]]
  data.frame(algo = r$algo, delta = r$delta, lr = r$lr, episodes = r$episodes,
             train_seed = r$seed, traj = seq_along(r$eval), return_T15 = r$eval)
}))
long <- rbind(long, data.frame(algo = "PBVI", delta = NA, lr = NA, episodes = NA,
                               train_seed = NA, traj = seq_along(pb), return_T15 = pb))
write.csv(long, file.path(out_dir, "eval_returns_per_trajectory_T15.csv"), row.names = FALSE)

# Learning curves (moving average, window 200) for default configurations
ma <- function(x, w = 200) as.numeric(stats::filter(x, rep(1 / w, w), sides = 1))
is_default <- function(r) (r$algo == "QL" && r$delta == 0.1 && r$lr == 0.01) ||
                          (r$algo == "REINFORCE" && r$lr == 1e-4)
lc <- do.call(rbind, lapply(Filter(is_default, res), function(r)
  data.frame(algo = r$algo, delta = r$delta, lr = r$lr, episodes = r$episodes,
             train_seed = r$seed, episode = seq_along(r$returns),
             return = r$returns, ma200 = ma(r$returns))))
write.csv(lc, file.path(out_dir, "learning_curves_default.csv"), row.names = FALSE)
print(summ, digits = 4)
print(pbvi_row)
