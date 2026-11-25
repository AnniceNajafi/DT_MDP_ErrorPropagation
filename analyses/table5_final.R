# table5_final.R  (from table4_rerun.R)
# Changes: k-Step uses the k selected by the held-out sweep (KStepHeuristic_Summary.csv);
# all policies are evaluated per sample (R/sim_persample.R);
# adds Myopic-belief and QMDP rows; stores cumulative-reward curves and per-trajectory
# fraction in Nominal for Fig. 3 so the figure and Table 5 come from the same trajectories.
# Regenerates Table 4 (policy summary) and the Wilcoxon / Cliff's delta table
# with the current code. Uses the original model code unchanged:
#   R/BaseMDPFormulation.R and R/POMDPFormulation.R (model + PBVI solve),
#   R/kstep_policy.R (gillespie_kstep), R/rl_core.R (Q-learning, REINFORCE),
#   and the metric definitions of the earlier Revision_PolicySummaryTable.R.
# Differences from Revision_PolicySummaryTable.R:
#   * PBVI is solved with 500 iterations (as stated in the paper).
#   * RL rows: 5 agents per method trained at the paper's default settings
#     (training seeds 123-127 / 456-460); trajectory n is run by agent ((n-1) %% 5) + 1.
#   * Common random numbers: trajectory n of every policy uses set.seed(700000 + n).
#   * The same trajectories feed the Wilcoxon / Cliff's delta table.
# Run from the repository root:  Rscript analyses/table5_final.R

local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
out_dir <- dirname(out_path("table5_final", "x"))
rl_core <- repo_path("R", "rl_core.R")
suppressPackageStartupMessages({ library(depmixS4); library(clue); library(expm) })

# Model: R/POMDPFormulation.R builds the data-driven model (sections 1-5 of
# R/BaseMDPFormulation.R, frac_in_nominal with the single-event guard) and
# solves PBVI exactly as this script originally did (set.seed(2024) belief
# points, 500 iterations, tol = 0.001, MDP warm start).
load_core()
source(repo_path("R", "kstep_policy.R"))   # gillespie_kstep

# Quick mode: fewer RL agents / episodes and fewer trajectories.
N_AGENTS <- nq(5L, 2L)
QL_EPISODES <- nq(5000L, 500L)
RE_EPISODES <- nq(10000L, 1000L)

source(rl_core)
cat("Training RL agents at default settings...\n")
ql_agents <- lapply(123L + seq_len(N_AGENTS) - 1L, function(sd) {
  f <- train_qlearning(delta = 0.1, alpha_ql = 0.01, n_episodes_ql = QL_EPISODES, seed = sd)
  ql_selector(f$Q_table, f$belief_grid)
})
re_agents <- lapply(456L + seq_len(N_AGENTS) - 1L, function(sd) {
  f <- train_reinforce(lr_re = 1e-4, n_episodes_re = RE_EPISODES, seed = sd)
  re_selector(f$theta)
})

# ---- metric helpers (copied from Revision_PolicySummaryTable.R) ----
first_detection_delay <- function(tr, target_state, canonical_action) {
  hit <- which(tr$states == target_state)
  if (length(hit) == 0) return(NA_real_)
  t_enter <- tr$times[hit[1]]
  applied <- which(tr$actions == canonical_action & tr$times >= t_enter)
  if (length(applied) == 0) return(NA_real_)
  tr$times[applied[1]] - t_enter
}
mismatch_rate <- function(tr, mdp_policy, T_end) {
  n <- length(tr$times); wrong_time <- 0
  for (i in seq_len(n - 1)) {
    if (tr$actions[i] != mdp_policy[states[tr$states[i]]])
      wrong_time <- wrong_time + (tr$times[i + 1] - tr$times[i])
  }
  if (tr$actions[n] != mdp_policy[states[tr$states[n]]])
    wrong_time <- wrong_time + (T_end - tr$times[n])
  wrong_time / T_end
}
gillespie_det_wrap <- function(Q_list, policy, R_rate, s0, T_end, r_disc) {
  tr <- gillespie_mdp(Q_list, policy, R_rate, s0, T_end, r_disc)
  list(times = tr$times, states = tr$states, actions = policy[states[tr$states]],
       disc_return = tr$disc_return)
}
canonical_repair_tbl <- c(SensorNoisy = "DwSensorA", DynamicsOff = "ReidentPlant",
                          Drift = "BiasCorrect")

N_mc <- nq(2000, 200); T_end <- 20; seed0 <- 700000
t_pts <- seq(0, T_end, by = 0.2)
ch_myo  <- function(b) actions[which.max(R_mat %*% b)]
ch_qmdp <- function(b) actions[which.max(vi$Q %*% b)]
run_policy <- function(label, sim_fn) {
  cat("--", label, "\n")
  rets <- fracs <- mism <- numeric(N_mc)
  curves <- matrix(NA_real_, N_mc, length(t_pts))
  delays <- matrix(NA, N_mc, 3, dimnames = list(NULL, names(canonical_repair_tbl)))
  for (n in seq_len(N_mc)) {
    set.seed(seed0 + n)
    tr <- sim_fn(n)
    rets[n] <- tr$disc_return
    curves[n, ] <- tr$cum
    fracs[n] <- frac_in_nominal(tr, T_end)
    for (g in names(canonical_repair_tbl))
      delays[n, g] <- first_detection_delay(tr, match(g, states), canonical_repair_tbl[g])
    mism[n] <- mismatch_rate(tr, vi$policy, T_end)
  }
  list(label = label, rets = rets, curves = curves, fracs = fracs, mism = mism,
       row = data.frame(Policy = label, MeanReward = mean(rets), SDReward = sd(rets),
                        SEReward = sd(rets) / sqrt(N_mc), FracNominal = mean(fracs),
                        Delay_SensorNoisy = mean(delays[, 1], na.rm = TRUE),
                        Delay_DynamicsOff = mean(delays[, 2], na.rm = TRUE),
                        Delay_Drift = mean(delays[, 3], na.rm = TRUE),
                        MismatchPct = 100 * mean(mism)))
}

# k-Step threshold: best k of the held-out sweep (analyses/Revision_KStepHeuristicBaseline.R)
k_file <- out_path("tables", "KStepHeuristic_Summary.csv")
K_STEP <- if (file.exists(k_file)) { ks <- read.csv(k_file); ks$k[which.max(ks$mean_return)] } else 1L
cat("k-Step threshold k =", K_STEP, "
")

# Per-sample evaluation (R/sim_persample.R): every policy sees one decoder
# label per 20 ms sample; trajectory n of every policy uses set.seed(seed0 + n).
ps <- function(ctl) function(n) sim_ps(ctl, A, O_mat, s0 = 1L, T_end = T_end, t_pts = t_pts)
res <- list(
  run_policy("Optimal MDP", ps(ctl_state(vi$policy))),
  run_policy("POMDP (PBVI)", ps(ctl_belief(alpha_choose(pomdp$Gamma, pomdp$actions)))),
  run_policy("Q-Learning", function(n) ps(ctl_belief(ql_agents[[(n - 1) %% N_AGENTS + 1]]))(n)),
  run_policy("REINFORCE", function(n) ps(ctl_belief(re_agents[[(n - 1) %% N_AGENTS + 1]]))(n)),
  run_policy("Myopic belief", ps(ctl_belief(ch_myo))),
  run_policy("QMDP", ps(ctl_belief(ch_qmdp))),
  run_policy(sprintf("k-Step (k=%d)", K_STEP), ps(ctl_kstep(K_STEP, dwell_sec = 0.5))),
  run_policy("MCDA (TOPSIS)", ps(ctl_mcda(mcda_dwell_time))),
  run_policy("No Intervention", ps(ctl_state(policy_none)))
)
tbl <- do.call(rbind, lapply(res, `[[`, "row"))
saveRDS(list(t_pts = t_pts, res = res), file.path(out_dir, "table5_runs.rds"))
write.csv(tbl, file.path(out_dir, "PolicySummaryTable.csv"), row.names = FALSE)
R <- sapply(res, `[[`, "rets"); colnames(R) <- tbl$Policy
write.csv(R, file.path(out_dir, "returns_per_trajectory.csv"), row.names = FALSE)

# ---- Wilcoxon (one-sided, adjacent pairs in order of mean) + Cliff's delta ----
cliff <- function(x, y) mean(sign(outer(x, y, "-")))
ord <- order(-tbl$MeanReward)
st <- do.call(rbind, lapply(seq_len(length(ord) - 1), function(i) {
  a <- ord[i]; b <- ord[i + 1]
  data.frame(Comparison = paste(tbl$Policy[a], "vs", tbl$Policy[b]),
             p = wilcox.test(R[, a], R[, b], alternative = "greater")$p.value,
             paired_p = wilcox.test(R[, a], R[, b], paired = TRUE, alternative = "greater")$p.value,
             cliff_d = cliff(R[, a], R[, b]))
}))
write.csv(st, file.path(out_dir, "stat_tests.csv"), row.names = FALSE)
options(width = 200)
print(tbl, digits = 4); print(st, digits = 3)
