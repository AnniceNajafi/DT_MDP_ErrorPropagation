#' Revision_KStepHeuristicBaseline.R
#' ──────────────────────────────────────────────────────────────────────
#' Reviewer Comment 5:
#'   "A better baseline would be intervene if the HMM label is
#'    non-nominal for k consecutive steps."
#'
#' This script implements a reactive-operator baseline that:
#'   (1) receives the noisy HMM/Viterbi observation at each transition,
#'   (2) applies NoAction until the observed label has been non-Nominal
#'       for k consecutive observations,
#'   (3) upon threshold crossing, applies the canonical repair for the
#'       dominant non-Nominal label (DwSensorA for SensorNoisy,
#'       ReidentPlant for DynamicsOff, BiasCorrect for Drift).
#'
#' We then sweep k in {1, 2, 3, 5}, keeping all other environment
#' parameters (Q_rate, O_mat, c_rate, r_disc) identical to those used
#' by BaseMDPFormulation.R and POMDPFormulation.R.
#'
#' Outputs:
#'   - output/figures/KStepHeuristic_RewardSweep.pdf  mean discounted return vs k
#'   - output/tables/KStepHeuristic_Summary.csv      per-k mean, SD, frac-Nominal
#'
#' Authors@R: person("Annice", "Najafi")
#' ──────────────────────────────────────────────────────────────────────

local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
load_core()                                  # pipeline -> MDP -> POMDP
source(repo_path("R", "kstep_policy.R"))     # canonical_repair, gillespie_kstep

library(ggplot2)

set.seed(42)

# (frac_in_nominal with the single-event guard is defined in
#  R/BaseMDPFormulation.R.)

# =====================================================================
# 3. SWEEP k IN {1, 2, 3, 5, 10, 25} (k counts decoder samples of 20 ms)
# =====================================================================

k_values <- c(1L, 2L, 3L, 5L, 10L, 25L)
N_traj   <- nq(1000, 100)
T_end    <- 20

results <- data.frame(
  k            = integer(0),
  mean_return  = numeric(0),
  sd_return    = numeric(0),
  frac_nominal = numeric(0)
)

for (k_t in k_values) {
  cat(sprintf("── k = %d  (N = %d trajectories) ──\n", k_t, N_traj))
  rets  <- numeric(N_traj)
  fracs <- numeric(N_traj)
  for (n in seq_len(N_traj)) {
    # per-sample evaluation (R/sim_persample.R); seeds separate from Table 5
    tr <- sim_ps(ctl_kstep(k_t, dwell_sec = 0.5), A, O_mat, s0 = 1L, T_end = T_end)
    rets[n]  <- tr$disc_return
    fracs[n] <- frac_in_nominal(tr, T_end)
  }
  results <- rbind(results, data.frame(
    k            = k_t,
    mean_return  = mean(rets),
    sd_return    = sd(rets),
    frac_nominal = mean(fracs)
  ))
  cat(sprintf("   mean reward = %.2f  (SD %.2f)   frac-Nominal = %.3f\n",
              mean(rets), sd(rets), mean(fracs)))
}

best_k <- results$k[which.max(results$mean_return)]
cat(sprintf("\nBest k = %d   (mean reward %.2f)\n",
            best_k, max(results$mean_return)))

write.csv(results, out_path("tables", "KStepHeuristic_Summary.csv"), row.names = FALSE)

# =====================================================================
# 4. FIGURE: mean return vs k
# =====================================================================

fig <- ggplot(results, aes(factor(k), mean_return)) +
  geom_col(fill = "#5B8FF9", width = 0.5) +
  geom_errorbar(aes(ymin = mean_return - sd_return,
                    ymax = mean_return + sd_return),
                width = 0.15) +
  labs(x = "k (consecutive non-Nominal decoder samples)",
       y = "Mean discounted return",
       title = "k-step consecutive heuristic baseline") +
  theme_minimal(base_size = 13)

pdf(out_path("figures", "KStepHeuristic_RewardSweep.pdf"), width = 7, height = 4.5)
print(fig)
dev.off()
cat("  -> Saved KStepHeuristic_RewardSweep.pdf\n")

# Exported summary of the best k
kstep_best_mean <- max(results$mean_return)
kstep_best_sd   <- results$sd_return[which.max(results$mean_return)]
kstep_best_k    <- best_k
kstep_frac_nom  <- results$frac_nominal[which.max(results$mean_return)]

cat("\n══ Revision_KStepHeuristicBaseline complete. ══\n")
