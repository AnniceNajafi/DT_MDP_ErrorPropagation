#' Revision_RewardSensitivity.R
#' ──────────────────────────────────────────────────────────────────────
#' Reviewer Comment 2 (first review round):
#'   "The reward values are presented without sensitivity analysis on
#'    reward values. The repair-probability sweep is useful but does not
#'    cover reward structure variation, please add."
#'
#' Three perturbations of the reward matrix R(a, s)
#' (Supplementary Section "Reward perturbations", Fig. RewardSensitivity):
#'   (A) Uniform scaling        R' = alpha * R               alpha in {0.5, 1, 2}
#'   (B) Penalty-only scaling   R' = beta * R_neg + R_pos    beta  in {0.5, 1, 2}
#'   (C) Gaussian noise         R' = R + N(0, sigma^2)       sigma in {0.05, 0.10, 0.20},
#'                                                           20 replicates per sigma
#'
#' For each perturbation we
#'   1. rebuild the reward rate c_rate' = R' / dt,
#'   2. re-solve the MDP by value iteration,
#'   3. re-solve the POMDP by PBVI (MDP warm start, up to 2000 iterations),
#'   4. simulate N_mc trajectories (T = 20 s, start in Nominal) for
#'      MDP, POMDP, k-Step (k = 2), MCDA (TOPSIS) and No Intervention.
#'
#' This file merges the former Revision_RewardSensitivity.R (paper
#' settings: 200 trajectories for A/B, 100 for C, 20 replicates) and the
#' data part of Revision_Runner.R (policy labels / palette). The figure
#' is drawn by figures/Revision_Replot.R.
#'
#' Outputs:
#'   - output/analyses/Revision_Sensitivity_df{A,B,C}.rds
#'   - output/tables/RewardSensitivity_{Alpha,Beta,Sigma}.csv
#'
#' Authors@R: person("Annice", "Najafi")
#' ──────────────────────────────────────────────────────────────────────

local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
load_core()                                  # pipeline -> MDP -> POMDP
source(repo_path("R", "kstep_policy.R"))     # gillespie_kstep()

set.seed(42)

N_mc_AB      <- nq(200, 20)     # trajectories per policy, sweeps A and B
N_mc_C       <- nq(100, 20)     # trajectories per policy, sweep C
n_rep        <- nq(20, 2)       # Gaussian replicates per sigma
pbvi_horizon <- nq(2000, 100)   # PBVI iterations per re-solve (as in the original script)
k_file <- out_path("tables", "KStepHeuristic_Summary.csv")   # k of the policy summary table
k_step <- if (file.exists(k_file)) { ks <- read.csv(k_file); as.integer(ks$k[which.max(ks$mean_return)]) } else 1L
T_end        <- 20

policy_labels <- c("Optimal MDP", "POMDP (PBVI)", sprintf("k-Step (k=%d)", k_step),
                   "MCDA (TOPSIS)", "No Intervention")

# =====================================================================
# 1. HELPER: re-solve MDP and POMDP under a perturbed reward matrix
# =====================================================================

solve_under_reward <- function(R_pert) {
  vi_new <- value_iteration(A, R_pert, gamma_disc,
                            tol = 1e-10, max_iter = 10000)
  pbvi_new <- pbvi_solve(A, O_mat, R_pert, gamma_disc, belief_pts,
                         max_horizon = pbvi_horizon, tol = 0.001,
                         V_init = vi_new$V)
  list(vi = vi_new, pbvi = pbvi_new, c_rate = R_pert / dt)
}

# =====================================================================
# 2. HELPER: Monte Carlo evaluation of the five policies
# =====================================================================

evaluate_policies <- function(sol, N_mc) {
  rets <- matrix(NA_real_, N_mc, 5)
  for (n in seq_len(N_mc)) {
    s0 <- 1L
    # per-sample evaluation (R/sim_persample.R) under the perturbed reward
    Rp <- sol$c_rate * dt
    rets[n, 1] <- sim_ps(ctl_state(sol$vi$policy), A, O_mat, s0, T_end, R = Rp)$disc_return
    rets[n, 2] <- sim_ps(ctl_belief(alpha_choose(sol$pbvi$Gamma, sol$pbvi$actions)), A, O_mat,
                         s0, T_end, R = Rp)$disc_return
    rets[n, 4] <- sim_ps(ctl_mcda(mcda_dwell_time), A, O_mat, s0, T_end, R = Rp)$disc_return
    rets[n, 3] <- sim_ps(ctl_kstep(k_step, dwell_sec = 0.5), A, O_mat, s0, T_end, R = Rp)$disc_return
    rets[n, 5] <- sim_ps(ctl_state(policy_none), A, O_mat, s0, T_end, R = Rp)$disc_return
  }
  data.frame(policy = policy_labels,
             mean   = colMeans(rets),
             sd     = apply(rets, 2, sd))
}

# =====================================================================
# 3. SWEEPS
# =====================================================================

cat("\n── Reward sensitivity A: uniform scaling alpha ──\n")
dfA <- data.frame()
for (al in c(0.5, 1.0, 2.0)) {
  cat(sprintf("   alpha = %.2f\n", al))
  d <- evaluate_policies(solve_under_reward(al * R_mat), N_mc_AB)
  d$level <- as.character(al); dfA <- rbind(dfA, d)
}

cat("\n── Reward sensitivity B: penalty-only scaling beta ──\n")
dfB <- data.frame()
for (be in c(0.5, 1.0, 2.0)) {
  cat(sprintf("   beta = %.2f\n", be))
  R_p <- R_mat; neg <- R_p < 0; R_p[neg] <- be * R_p[neg]
  d <- evaluate_policies(solve_under_reward(R_p), N_mc_AB)
  d$level <- as.character(be); dfB <- rbind(dfB, d)
}

cat(sprintf("\n── Reward sensitivity C: Gaussian sigma (%d replicates) ──\n", n_rep))
dfC <- data.frame()
for (sg in c(0.05, 0.10, 0.20)) {
  cat(sprintf("   sigma = %.2f\n", sg))
  for (rep in seq_len(n_rep)) {
    R_p <- R_mat + matrix(rnorm(length(R_mat), 0, sg),
                          nrow(R_mat), ncol(R_mat), dimnames = dimnames(R_mat))
    d <- evaluate_policies(solve_under_reward(R_p), N_mc_C)
    d$level <- as.character(sg); d$rep <- rep
    dfC <- rbind(dfC, d)
  }
}

# =====================================================================
# 4. SAVE
# =====================================================================

for (nm in c("dfA", "dfB", "dfC")) {
  d <- get(nm); d$policy <- factor(d$policy, levels = policy_labels)
  assign(nm, d)
}
saveRDS(dfA, out_path("analyses", "Revision_Sensitivity_dfA.rds"))
saveRDS(dfB, out_path("analyses", "Revision_Sensitivity_dfB.rds"))
saveRDS(dfC, out_path("analyses", "Revision_Sensitivity_dfC.rds"))
write.csv(dfA, out_path("tables", "RewardSensitivity_Alpha.csv"), row.names = FALSE)
write.csv(dfB, out_path("tables", "RewardSensitivity_Beta.csv"),  row.names = FALSE)
write.csv(dfC, out_path("tables", "RewardSensitivity_Sigma.csv"), row.names = FALSE)
print(dfA); print(dfB)
print(aggregate(mean ~ policy + level, data = dfC, FUN = mean))

cat("\n══ Revision_RewardSensitivity complete. ══\n")
