# Policies under a reward that depends on the restoration of fidelity
#  The reward of the paper, R(a, s), scores whether the action matches the regime. Here it is
#  replaced by R_fid(a, s) = r(s) - c(a), with r(Nominal) = +1, r(s) = -1 otherwise, c(NoAction) = 0
#  and c(a) = c for every intervention. An intervention then pays off only through the repair it
#  causes in later steps. MDP, POMDP (PBVI solved to convergence), myopic belief, QMDP, k-Step,
#  MCDA and No Intervention are evaluated with the per-sample simulator of Table 5, for the repair
#  probabilities of Table 3 and for repair probabilities scaled by 0.1, and for c in {0.05, 0.1, 0.2}.
#  (Discussion; Supplementary Section S6.10)
# Run from the repository root:  Rscript analyses/fidelity_reward.R
# Writes output/tables/Fidelity_Reward.txt and output/tables/Fidelity_Reward.csv
Sys.setenv(DTMDP_REUSE_MATRICES = "1"); source(file.path(getwd(), "R", "paths.R"))
suppressPackageStartupMessages({ library(depmixS4); library(clue); library(expm) })
invisible(capture.output(load_core()))
out_file <- out_path("tables", "Fidelity_Reward.txt")
sink(out_file, split = TRUE)

N <- nq(1000, 100); seed0 <- 700000; T_end <- 20
k_file <- out_path("tables", "KStepHeuristic_Summary.csv")
K_STEP <- if (file.exists(k_file)) { ks <- read.csv(k_file); ks$k[which.max(ks$mean_return)] } else 1L

build_A <- function(scale) { Al <- list(); for (a in actions) { Al[[a]] <- A_NoAction
  for (s in 1:K) if (rho[a, s] > 0) { r <- scale * rho[a, s]
    Al[[a]][s, ] <- (1 - r) * A_NoAction[s, ]; Al[[a]][s, 1] <- Al[[a]][s, 1] + r } }; Al }
fid_reward <- function(c) { Rf <- t(sapply(actions, function(a) c(1, -1, -1, -1) - if (a == "NoAction") 0 else c))
  dimnames(Rf) <- dimnames(R_mat); Rf }

evaluate <- function(scale, c, policies) {
  Al <- build_A(scale); Rf <- fid_reward(c)
  vif <- value_iteration(Al, Rf, gamma_disc)
  msg <- capture.output(P <- pbvi_solve(Al, O_mat, Rf, gamma_disc, belief_pts,
                                        max_horizon = 5000, tol = 1e-6, V_init = vif$V))
  ctl <- list(
    "Optimal MDP"     = function() ctl_state(vif$policy),
    "POMDP (PBVI)"    = function() ctl_belief(alpha_choose(P$Gamma, P$actions), Al, O_mat),
    "Myopic belief"   = function() ctl_belief(function(b) actions[which.max(Rf %*% b)], Al, O_mat),
    "QMDP"            = function() ctl_belief(function(b) actions[which.max(vif$Q %*% b)], Al, O_mat),
    "k-Step"          = function() ctl_kstep(K_STEP, dwell_sec = 0.5),
    "MCDA (TOPSIS)"   = function() ctl_mcda(mcda_dwell_time),
    "No Intervention" = function() ctl_state(setNames(rep("NoAction", K), states)))[policies]
  cat(sprintf("\n== repair probabilities x %.2f, intervention cost %.2f ; MDP policy: %s ; %s\n",
              scale, c, paste(vif$policy, collapse = ", "), tail(msg, 1)))
  rows <- lapply(names(ctl), function(nm) {
    r <- t(sapply(1:N, function(n) { set.seed(seed0 + n)
      tr <- sim_ps(ctl[[nm]](), Al, O_mat, s0 = 1L, T_end = T_end, R = Rf)
      c(tr$disc_return, frac_in_nominal(tr, T_end)) }))
    data.frame(RepairScale = scale, Cost = c, Policy = nm, MeanReturn = mean(r[, 1]),
               SE = sd(r[, 1]) / sqrt(N), FracNominal = mean(r[, 2])) })
  tab <- do.call(rbind, rows)
  tab$PctOfMDP <- 100 * tab$MeanReturn / tab$MeanReturn[tab$Policy == "Optimal MDP"]
  print(format(tab[, -(1:2)], digits = 4), row.names = FALSE)
  tab
}

all_pol <- c("Optimal MDP", "POMDP (PBVI)", "Myopic belief", "QMDP", "k-Step", "MCDA (TOPSIS)", "No Intervention")
res <- rbind(evaluate(1,   0.10, all_pol),
             evaluate(0.1, 0.10, all_pol),
             evaluate(1,   0.05, all_pol[1:4]),
             evaluate(1,   0.20, all_pol[1:4]))
write.csv(res, out_path("tables", "Fidelity_Reward.csv"), row.names = FALSE)
sink()
