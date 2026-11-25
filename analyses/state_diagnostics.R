# State diagnostics quoted in the text:
#  (1) regime occupancy over time under four controllers (Supplementary Section S5.2, Fig. S3);
#  (2) where the POMDP's action mismatches occur (main Section 5.4).
# Run from the repository root:  Rscript analyses/state_diagnostics.R
# Writes output/tables/State_Diagnostics.txt
Sys.setenv(DTMDP_REUSE_MATRICES = "1"); source(file.path(getwd(), "R", "paths.R"))
suppressPackageStartupMessages({ library(depmixS4); library(clue); library(expm) })
invisible(capture.output(load_core()))
out_file <- file.path("output", "tables", "State_Diagnostics.txt")
sink(out_file, split = TRUE)

ev <- eigen(t(A[["NoAction"]])); st <- Re(ev$vectors[, 1]); st <- st / sum(st)
cat("Stationary distribution of A(NoAction):", round(st, 3), "\n\n")

# (1) occupancy over time: 600 runs of 15 s from Nominal
ctls <- list(MDP = ctl_state(vi$policy), POMDP = ctl_belief(alpha_choose(pomdp$Gamma, pomdp$actions)),
             MCDA = ctl_mcda(mcda_dwell_time), None = ctl_state(setNames(rep("NoAction", K), states)))
cat("Occupancy over time (600 runs, 15 s, start in Nominal)\n")
for (nm in names(ctls)) {
  Mo <- matrix(0, 750, 4)
  for (i in 1:600) {
    set.seed(900000 + i)
    tr <- sim_ps(ctls[[nm]], A, O_mat, 1L, 15, record_belief = TRUE)
    ss <- tr$step_states[1:750]; Mo[cbind(1:750, ss)] <- Mo[cbind(1:750, ss)] + 1
  }
  Mo <- Mo / 600; late <- Mo[251:750, ]
  cat(sprintf("%-6s Nominal after 5 s: mean %.3f, range %.3f-%.3f | at t = 15 s: %s\n", nm,
              mean(late[, 1]), min(late[, 1]), max(late[, 1]), paste(round(Mo[750, ], 3), collapse = " ")))
}

# (2) POMDP mismatch decomposition on the seeds of Table 5
ctl <- ctl_belief(alpha_choose(pomdp$Gamma, pomdp$actions))
tot <- matrix(0, 4, 2); rets <- numeric(2000)
for (n in 1:2000) {
  set.seed(700000 + n)
  tr <- sim_ps(ctl, A, O_mat, s0 = 1L, T_end = 20, t_pts = seq(0, 20, by = 0.2), record_belief = TRUE)
  rets[n] <- tr$disc_return
  s <- tr$step_states; a <- tr$step_actions; wrong <- a != vi$policy[states[s]]
  for (j in 1:4) { tot[j, 1] <- tot[j, 1] + sum(s == j); tot[j, 2] <- tot[j, 2] + sum(wrong & s == j) }
}
cat("\nPOMDP on the Table 5 seeds: mean return", round(mean(rets), 3), "\n")
cat("Overall mismatch rate:", round(sum(tot[, 2]) / sum(tot[, 1]), 3), "\n")
cat("Share of mismatched samples by regime (", paste(states, collapse = ", "), "):", round(tot[, 2] / sum(tot[, 2]), 3), "\n")
cat("Mismatch rate within each regime:", round(tot[, 2] / tot[, 1], 3), "\n")
sink()
