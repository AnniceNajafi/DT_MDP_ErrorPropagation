# PBVI solved to convergence and the remaining MDP-POMDP gap at perfect SensorNoisy accuracy
#  (1) warm-start and lower-bound-start PBVI on the belief set of the main solver, iterated until
#      successive value functions differ by less than 1e-6; value at b0 = Nominal, the
#      value-of-information bracket with the fast informed bound, and the policy's returns
#      compared with the 500-iteration solution (Sections 4.8, 5.3, 5.6.1; Supplementary S4, S6.2);
#  (2) POMDP returns when the SensorNoisy row of O is made perfect, and when the Nominal or the
#      DynamicsOff row is made perfect as well (Section 5.3; Supplementary S5.2, S5.4).
# Run from the repository root:  Rscript analyses/pbvi_bounds_and_gap.R
# Writes output/tables/PBVI_Bounds_and_Gap.txt
Sys.setenv(DTMDP_REUSE_MATRICES = "1"); source(file.path(getwd(), "R", "paths.R"))
suppressPackageStartupMessages({ library(depmixS4); library(clue); library(expm) })
invisible(capture.output(load_core()))
out_file <- out_path("tables", "PBVI_Bounds_and_Gap.txt")
sink(out_file, split = TRUE)

b0 <- c(1, 0, 0, 0)
val_at <- function(P, b) max(do.call(rbind, P$Gamma) %*% b)
N <- nq(1000, 100); seed0 <- 700000; T_end <- 20
ret <- function(ctl, O_true) sapply(1:N, function(n) { set.seed(seed0 + n); sim_ps(ctl, A, O_true, 1L, T_end)$disc_return })

# (1) converged PBVI and bounds
cat("==== (1) PBVI to convergence (tol 1e-6) on the main belief set (", length(belief_pts), " points) ====\n", sep = "")
cvg <- function(V0) { msg <- capture.output(P <- pbvi_solve(A, O_mat, R_mat, gamma_disc, belief_pts,
                       max_horizon = 5000, tol = 1e-6, V_init = V0)); list(P = P, last = tail(msg, 1)) }
w <- cvg(vi$V); l <- cvg(rep(min(R_mat) / (1 - gamma_disc), K))
cat("warm start:        ", w$last, "\n"); cat("lower-bound start: ", l$last, "\n")
fib <- function(A, O, R, g, tol = 1e-9) { Qf <- R; repeat { Qn <- Qf
  for (ai in 1:M) { a <- actions[ai]; for (s in 1:K) { acc <- 0
    for (z in 1:K) acc <- acc + max(sapply(1:M, function(a2) sum(A[[a]][s, ] * O[, z] * Qf[a2, ])))
    Qn[ai, s] <- R[ai, s] + g * acc } }
  if (max(abs(Qn - Qf)) < tol) break; Qf <- Qn }; Qf }
UB <- max(fib(A, O_mat, R_mat, gamma_disc) %*% b0); Vm <- sum(vi$V * b0)
cat(sprintf("b0 = Nominal: V_MDP = %.3f  fast informed bound = %.3f  converged PBVI: warm = %.3f, lower-bound start = %.3f (500 iterations: warm = %.3f)\n",
            Vm, UB, val_at(w$P, b0), val_at(l$P, b0), val_at(pomdp, b0)))
LBv <- max(val_at(w$P, b0), val_at(l$P, b0))
cat(sprintf("VoI bracket: [%.2f, %.2f]\n", Vm - UB, Vm - val_at(l$P, b0)))
set.seed(99); Bt <- matrix(rgamma(5000 * K, 1, 1), ncol = K); Bt <- Bt / rowSums(Bt)
ch <- function(P) function(b) alpha_choose(P$Gamma, P$actions)(b)
agree <- function(P1, P2) mean(apply(Bt, 1, function(b) ch(P1)(b) == ch(P2)(b)))
cat(sprintf("greedy-action agreement on 5000 random beliefs: converged warm vs 500-iteration %.4f; converged warm vs converged lower-bound start %.4f\n",
            agree(w$P, pomdp), agree(w$P, l$P)))
r500 <- ret(ctl_belief(alpha_choose(pomdp$Gamma, pomdp$actions)), O_mat)
rcv  <- ret(ctl_belief(alpha_choose(w$P$Gamma, w$P$actions)), O_mat)
d <- rcv - r500
cat(sprintf("returns over %d trajectories (Table 5 seeds): 500-iteration %.3f, converged %.3f, paired difference %.3f (SE %.3f)\n",
            N, mean(r500), mean(rcv), mean(d), sd(d) / sqrt(N)))

# (2) where the gap at p_SN = 1 comes from
cat("\n==== (2) Gap at perfect SensorNoisy accuracy (PBVI re-solved, 500 iterations, warm start) ====\n")
O_sn <- O_mat; O_sn[2, ] <- c(0, 1, 0, 0)
O_sn_do <- O_sn; O_sn_do[3, ] <- c(0, 0, 1, 0)
O_sn_nom <- O_sn; O_sn_nom[1, ] <- c(1, 0, 0, 0)
cases <- list("SensorNoisy row perfect" = O_sn, "+ DynamicsOff row perfect" = O_sn_do, "+ Nominal row perfect" = O_sn_nom)
r_mdp <- ret(ctl_state(vi$policy), O_mat)
cat(sprintf("%-28s %.3f (SE %.3f)\n", "Optimal MDP", mean(r_mdp), sd(r_mdp) / sqrt(N)))
for (nm in names(cases)) { OO <- cases[[nm]]
  invisible(capture.output(P <- pbvi_solve(A, OO, R_mat, gamma_disc, belief_pts, max_horizon = pbvi_horizon, tol = 0.001, V_init = vi$V)))
  r <- ret(ctl_belief(alpha_choose(P$Gamma, P$actions), A, OO), OO)
  cat(sprintf("%-28s %.3f (SE %.3f)\n", nm, mean(r), sd(r) / sqrt(N))) }
cat("Nominal labelled SensorNoisy:", O_mat[1, 2], "; stationary occupancy of A(NoAction):",
    round(Re(eigen(t(A[["NoAction"]]))$vectors[, 1] / sum(Re(eigen(t(A[["NoAction"]]))$vectors[, 1]))), 3), "\n")
sink()
