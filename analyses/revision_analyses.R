# Revision analyses for JRESS-D-26-03120 (second review round)
# Produces the numbers of: Supplementary Table "Discrepancy between A and
# exp(Q dt)" (section 1), PBVI-initialisation check (2), myopic-belief /
# QMDP / LB-initialised PBVI rows (3), exact-log generator ranking (4),
# held-out causal confusion matrix (5), misspecified / correlated
# observation channels (6), TEP policy evaluation (7, needs
# output/tep/TEP_results.rds from tep/TEP_Application.R; skipped otherwise),
# and the slow-repair sweep (8).
# Output: output/analyses/revision_analyses/{revision_results.txt, causal_confusion.rds}
# Run from the repository root:  Rscript analyses/revision_analyses.R

local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
suppressPackageStartupMessages({ library(expm); library(depmixS4); library(clue); library(dplyr) })

# Monte Carlo sizes (paper values; reduced in quick mode)
N1000 <- nq(1000, 100); N500 <- nq(500, 50); N300 <- nq(300, 30)
IT500 <- nq(500, 100);  IT300 <- nq(300, 60); IT200 <- nq(200, 50)
N_CAUSAL_SEEDS <- nq(10L, 3L)

out_dir <- dirname(out_path("analyses", "revision_analyses", "revision_results.txt"))
logf <- file.path(out_dir, "revision_results.txt")
sink(logf, split = TRUE)

source(repo_path("R", "SimulationPipeline.R"))
cat("\n==== Viterbi (in-sample) confusion used as O ====\n"); print(round(mcda_confusion, 4))

states  <- c("Nominal", "SensorNoisy", "DynamicsOff", "Drift")
actions <- c("NoAction", "DwSensorA", "DwSensorB", "IncFilter", "ReidentPlant", "BiasCorrect")
K <- 4; M <- 6; dt <- 0.02; gamma_disc <- 0.99
rho <- matrix(0, M, K, dimnames = list(actions, states))
rho["DwSensorA","SensorNoisy"] <- .6; rho["DwSensorB","SensorNoisy"] <- .6
rho["IncFilter","SensorNoisy"] <- .4; rho["ReidentPlant","DynamicsOff"] <- .7
rho["BiasCorrect","Drift"] <- .8
build_A <- function(A0, rho) { A <- list(); for (a in actions) { A[[a]] <- A0
  for (s in 1:K) if (rho[a,s] > 0) { A[[a]][s,] <- (1-rho[a,s])*A0[s,]; A[[a]][s,1] <- A[[a]][s,1] + rho[a,s] } }; A }
A <- build_A(A_NoAction, rho)
R_mat <- matrix(c(1,-.5,-1,-.8, -.1,.7,-.2,-.1, -.1,.6,-.2,-.1, -.1,.4,0,.1, -.5,-.3,.9,-.2, -.2,0,.1,.8),
                M, byrow = TRUE, dimnames = list(actions, states))
vi_solve <- function(A, R, g, tol = 1e-10) { V <- rep(0,K)
  repeat { Q <- t(sapply(actions, function(a) R[a,] + g * as.numeric(A[[a]] %*% V)))
    Vn <- apply(Q, 2, max); if (max(abs(Vn-V)) < tol) break; V <- Vn }
  list(V = Vn, Q = Q, policy = setNames(actions[apply(Q, 2, which.max)], states)) }
vi <- vi_solve(A, R_mat, gamma_disc)
cat("\nMDP policy:", vi$policy, "\nV*:", round(vi$V, 3), "\n")

# ---------------------------------------------------------------- 1. CTMC approximation error
cat("\n==== 1. First-order generator vs exact embedding ====\n")
Q_rate <- lapply(A, function(Aa) { Q <- (Aa - diag(K))/dt; dimnames(Q) <- list(states, states); Q })
r_disc <- -log(gamma_disc)/dt; c_rate <- R_mat/dt
ct_err <- t(sapply(actions, function(a) {
  E1 <- expm(Q_rate[[a]] * dt)
  L <- tryCatch(Re(logm(A[[a]])) / dt, error = function(e) matrix(NA, K, K))
  offneg <- if (any(is.na(L))) NA else min(L[row(L) != col(L)])
  c(max_abs_A_minus_expQdt = max(abs(A[[a]] - E1)),
    max_rel_row_TV = max(rowSums(abs(A[[a]] - E1))/2),
    logm_min_offdiag = offneg)
}))
print(signif(ct_err, 3))
# Effective discrete kernel implied by the approximate generator, and policy/value under it
A_eff <- lapply(Q_rate, function(Q) { E <- expm(Q*dt); E[E < 0] <- 0; E/rowSums(E) })
vi_eff <- vi_solve(A_eff, R_mat, gamma_disc)
cat("Policy under expm(Q dt):", vi_eff$policy, " identical:", all(vi_eff$policy == vi$policy), "\n")
cat("max |V*_A - V*_expm| =", signif(max(abs(vi$V - vi_eff$V)), 3),
    " (relative", signif(max(abs(vi$V - vi_eff$V)/abs(vi$V)), 3), ")\n")
# Exact generator (logm), if it is a valid generator after clipping
Q_exact <- lapply(A, function(Aa) { L <- Re(logm(Aa))/dt; off <- L; diag(off) <- 0; off[off < 0] <- 0
  diag(off) <- -rowSums(off); dimnames(off) <- list(states, states); off })
ct_value <- function(Ql, pol) { Qp <- t(sapply(1:K, function(s) Ql[[pol[s]]][s,]))
  cp <- sapply(1:K, function(s) c_rate[pol[s], s]); solve(r_disc*diag(K) - Qp, cp) }
pol_none <- setNames(rep("NoAction", K), states)
cat("CT value, optimal policy: first-order", round(ct_value(Q_rate, vi$policy),3), " exact-log", round(ct_value(Q_exact, vi$policy),3), "\n")
cat("CT value, NoAction:       first-order", round(ct_value(Q_rate, pol_none),3), " exact-log", round(ct_value(Q_exact, pol_none),3), "\n")
# CT policy improvement under each generator
ct_opt <- function(Ql) { pol <- vi$policy; repeat { V <- ct_value(Ql, pol)
  newp <- sapply(1:K, function(s) actions[which.max(sapply(actions, function(a) c_rate[a,s] + sum(Ql[[a]][s,]*V)))])
  if (all(newp == pol)) break; pol <- setNames(newp, states) }; pol }
cat("CT-optimal policy (first-order):", ct_opt(Q_rate), "\nCT-optimal policy (exact log):  ", ct_opt(Q_exact), "\n")

# ---------------------------------------------------------------- POMDP machinery
O_mat <- mcda_confusion
gen_bp <- function(n_random = 185) { pts <- list(); for (i in 1:K) { b <- rep(0,K); b[i] <- 1; pts[[length(pts)+1]] <- b }
  for (i in 1:(K-1)) for (j in (i+1):K) { b <- rep(0,K); b[c(i,j)] <- .5; pts[[length(pts)+1]] <- b }
  cb <- combn(K,3); for (c in 1:ncol(cb)) { b <- rep(0,K); b[cb[,c]] <- 1/3; pts[[length(pts)+1]] <- b }
  pts[[length(pts)+1]] <- rep(1/K,K); for (r in 1:n_random) { g <- rgamma(K,1,1); pts[[length(pts)+1]] <- g/sum(g) }; do.call(rbind, pts) }
pbvi <- function(A, O, R, g, B, iters = 500, init) {
  G <- matrix(init, 1, K); lab <- "NoAction"; res <- numeric(iters)
  for (it in 1:iters) {
    backed <- lapply(actions, function(a) lapply(1:K, function(z) g * (A[[a]] %*% (t(G) * O[, z]))))
    Gn <- matrix(NA, nrow(B), K); Ln <- character(nrow(B))
    for (bi in 1:nrow(B)) { b <- B[bi,]; bv <- -Inf
      for (ai in 1:M) { al <- R[ai,]; for (z in 1:K) { bz <- backed[[ai]][[z]]; al <- al + bz[, which.max(crossprod(bz, b))] }
        v <- sum(al*b); if (v > bv) { bv <- v; ba <- al; bl <- actions[ai] } }
      Gn[bi,] <- ba; Ln[bi] <- bl }
    res[it] <- max(abs(rowSums(Gn*B) - apply(B %*% t(G), 1, max)))
    keep <- !duplicated(round(Gn, 10)); G <- Gn[keep,,drop=FALSE]; lab <- Ln[keep] }
  list(G = G, lab = lab, res = res) }
pol_b <- function(b, P) P$lab[which.max(P$G %*% b)]
val_b <- function(b, P) max(P$G %*% b)

set.seed(7); B <- gen_bp()
cat("\n==== 2. PBVI initialisation ====\n")
t0 <- Sys.time(); P_warm <- pbvi(A, O_mat, R_mat, gamma_disc, B, IT500, vi$V); t_warm <- Sys.time()-t0
Rmin_lb <- min(R_mat)/(1-gamma_disc)
P_lb <- pbvi(A, O_mat, R_mat, gamma_disc, B, IT500, rep(Rmin_lb, K))
P_zero <- pbvi(A, O_mat, R_mat, gamma_disc, B, IT500, rep(0, K))
set.seed(99); Btest <- matrix(rgamma(5000*K,1,1), ncol = K); Btest <- Btest/rowSums(Btest)
agree <- function(P1, P2, BB) mean(apply(BB, 1, pol_b, P1) == apply(BB, 1, pol_b, P2))
vdiff <- function(P1, P2, BB) max(abs(apply(BB,1,val_b,P1) - apply(BB,1,val_b,P2)))
cat("alpha-vectors: warm", nrow(P_warm$G), " LB", nrow(P_lb$G), " zero", nrow(P_zero$G), "\n")
cat("Policy agreement on 5000 random beliefs: warm vs LB", agree(P_warm,P_lb,Btest), " warm vs zero", agree(P_warm,P_zero,Btest), "\n")
cat("Policy agreement on grid beliefs: warm vs LB", agree(P_warm,P_lb,B), "\n")
cat("max |V_warm - V_LB| on random beliefs:", signif(vdiff(P_warm,P_lb,Btest),3), " on grid:", signif(vdiff(P_warm,P_lb,B),3), "\n")
cat("Final residual: warm", signif(tail(P_warm$res,1),3), " LB", signif(tail(P_lb$res,1),3), "\n")
cat("gamma^500 =", signif(gamma_disc^500,3), " ; initial upper-bound gap bound (Vmdp - Rmin/(1-g)) =", round(max(vi$V)-Rmin_lb,2), "\n")
cat("Corner actions warm:", sapply(1:K, function(s) pol_b(diag(K)[s,], P_warm)), "\n")
# Bounds at b0 = Nominal
fib <- function(A, O, R, g, tol = 1e-9) { Qf <- R; repeat { Qn <- Qf
  for (ai in 1:M) { a <- actions[ai]; for (s in 1:K) { acc <- 0
    for (z in 1:K) acc <- acc + max(sapply(1:M, function(a2) sum(A[[a]][s,] * O[,z] * Qf[a2,])))
    Qn[ai,s] <- R[ai,s] + g*acc } }
  if (max(abs(Qn-Qf)) < tol) break; Qf <- Qn }; Qf }
Qfib <- fib(A, O_mat, R_mat, gamma_disc)
b0 <- c(1,0,0,0)
UB <- max(Qfib %*% b0); LB <- val_b(b0, P_lb); Vm <- vi$V[1]
cat(sprintf("b0=Nominal: V_MDP=%.3f  FIB upper bound=%.3f  PBVI(LB-init) lower bound=%.3f  PBVI(warm)=%.3f\n", Vm, UB, LB, val_b(b0,P_warm)))
cat(sprintf("VoI bracket (discrete-time values): [%.3f, %.3f]\n", Vm-UB, Vm-LB))
bu <- rep(.25,4); cat(sprintf("b=uniform: V_MDP(b)=%.3f FIB=%.3f LB=%.3f\n", sum(vi$V*bu), max(Qfib%*%bu), val_b(bu,P_lb)))

# ---------------------------------------------------------------- Gillespie POMDP with pluggable obs model
gill <- function(Ql, chooser, O_true, O_belief, s0, T_end, kappa = 0, seed) {
  set.seed(seed); s <- s0; b <- rep(0,K); b[s0] <- 1; a <- chooser(b); t <- 0
  tt <- 0; ss <- s; aa <- a; prev_err <- NA
  repeat { ro <- pmax(Ql[[a]][s,],0); ro[s] <- 0; lam <- sum(ro); if (lam <= 0) break
    tau <- rexp(1, lam); t <- t + tau; if (t > T_end) break
    s <- sample.int(K, 1, prob = ro)
    if (!is.na(prev_err) && runif(1) < kappa && prev_err != s) z <- prev_err else z <- sample.int(K, 1, prob = O_true[s,])
    prev_err <- if (z != s) z else NA
    bp <- as.numeric(b %*% expm(Ql[[a]]*tau)) * O_belief[, z]; b <- if (sum(bp) > 0) bp/sum(bp) else rep(1/K,K)
    a <- chooser(b); tt <- c(tt,t); ss <- c(ss,s); aa <- c(aa,a) }
  te <- c(tt[-1], T_end); d <- te - tt
  ret <- sum(c_rate[cbind(match(aa, actions), ss)] / r_disc * exp(-r_disc*tt) * (1-exp(-r_disc*d)))
  fn <- sum(d[ss == 1])/T_end; mm <- sum(d[aa != vi$policy[ss]])/T_end
  c(ret = ret, fnom = fn, mism = mm) }
run <- function(chooser, N = N1000, Ql = Q_rate, O_true = O_mat, O_bel = O_mat, kappa = 0) {
  r <- t(sapply(1:N, function(n) gill(Ql, chooser, O_true, O_bel, 1, 20, kappa, 5000+n)))
  c(mean = mean(r[,1]), sd = sd(r[,1]), fnom = mean(r[,2]), mism = 100*mean(r[,3])) }
ch_warm <- function(b) pol_b(b, P_warm)
ch_lb   <- function(b) pol_b(b, P_lb)
ch_myo  <- function(b) actions[which.max(R_mat %*% b)]
ch_qmdp <- function(b) actions[which.max(vi$Q %*% b)]
ch_mdp_oracle <- NULL

cat("\n==== 3. Policy evaluation (1000 traj, T=20 s, common seeds) ====\n")
res <- rbind(
  POMDP_warm = run(ch_warm), POMDP_LBinit = run(ch_lb),
  QMDP = run(ch_qmdp), Myopic_belief = run(ch_myo),
  NoAction = run(function(b) "NoAction"))
print(round(res, 3))
# MDP oracle via full-state info: emulate with identity O
res_mdp <- run(function(b) vi$policy[which.max(b)], O_true = diag(K), O_bel = diag(K))
cat("MDP (perfect obs):", round(res_mdp, 3), "\n")

cat("\n==== 4. Exact-log generator: ranking check ====\n")
res_exact <- rbind(MDP = run(function(b) vi$policy[which.max(b)], Ql = Q_exact, O_true = diag(K), O_bel = diag(K)),
  POMDP = run(ch_warm, Ql = Q_exact), Myopic = run(ch_myo, Ql = Q_exact), NoAction = run(function(b) "NoAction", Ql = Q_exact))
print(round(res_exact, 3))

# ---------------------------------------------------------------- 5. causal HMM on held-out seeds
cat("\n==== 5. Causal (forward-filter) HMM on held-out residual sequences ====\n")
feat_fun <- function(E) { f <- E %>% transmute(e0, e1, e2a, e2b, e3, e4,
  e_abs_sum = abs(e0)+abs(e1)+abs(e2a)+abs(e2b)+abs(e3)+abs(e4),
  e_rms = sqrt((e0^2+e1^2+e2a^2+e2b^2+e3^2+e4^2)/6), e1_diff = c(0,diff(e1)), e3_diff = c(0,diff(e3))); as.data.frame(f) }
mk_mod <- function(df) depmix(response = lapply(names(df), function(nm) as.formula(paste0(nm,"~1"))),
                              data = df, nstates = 4, family = rep(list(gaussian()), ncol(df)))
conf_c <- matrix(0,K,K,dimnames=list(states,states)); conf_v <- conf_c; conf_in <- conf_c
runs_err <- c(); ac1 <- c(); n_ok <- 0
for (i in 0:(N_CAUSAL_SEEDS - 1)) {
  tr <- simulate_with_faults_pipe(pipe_params, pipe_schedule, 11+i)
  te <- simulate_with_faults_pipe(pipe_params, pipe_schedule, 111+i)
  Ftr <- feat_fun(tr$E); mu <- colMeans(Ftr); sdv <- apply(Ftr,2,sd)
  Ztr <- as.data.frame(scale(Ftr, mu, sdv)); Zte <- as.data.frame(scale(feat_fun(te$E), mu, sdv))
  fm <- tryCatch(fit(mk_mod(Ztr), verbose = FALSE), error = function(e) NULL); if (is.null(fm)) next
  vp <- viterbi(fm)$state; tab <- as.matrix(table(factor(vp,1:4), factor(tr$regime_t, states)))
  perm <- as.integer(solve_LSAP(max(tab)-tab)); map <- states[perm]
  for (t in seq_along(vp)) conf_in[tr$regime_t[t], map[vp[t]]] <- conf_in[tr$regime_t[t], map[vp[t]]] + 1
  mte <- setpars(mk_mod(Zte), getpars(fm))
  fb <- forwardbackward(mte); al <- fb$alpha; al <- al/rowSums(al)
  zc <- map[apply(al,1,which.max)]; zv <- map[viterbi(mte)$state]
  for (t in seq_along(zc)) { conf_c[te$regime_t[t], zc[t]] <- conf_c[te$regime_t[t], zc[t]] + 1
                             conf_v[te$regime_t[t], zv[t]] <- conf_v[te$regime_t[t], zv[t]] + 1 }
  err <- as.integer(zc != te$regime_t); rl <- rle(err); runs_err <- c(runs_err, rl$lengths[rl$values == 1])
  ac1 <- c(ac1, cor(err[-1], err[-length(err)])); n_ok <- n_ok + 1 }
O_causal <- conf_c/rowSums(conf_c)
cat("seeds used:", n_ok, "\nIn-sample Viterbi confusion (this run):\n"); print(round(conf_in/rowSums(conf_in),3))
cat("Held-out Viterbi confusion:\n"); print(round(conf_v/rowSums(conf_v),3))
cat("Held-out CAUSAL filtered confusion:\n"); print(round(O_causal,3))
cat(sprintf("Error run lengths (steps): mean %.1f, median %.0f, 90th pct %.0f; mean lag-1 autocorr of error indicator %.3f\n",
  mean(runs_err), median(runs_err), quantile(runs_err,.9), mean(ac1, na.rm=TRUE)))
saveRDS(list(O_causal = O_causal, conf_v = conf_v, conf_in = conf_in), file.path(out_dir, "causal_confusion.rds"))

cat("\n==== 6. Misspecified / correlated observation models (belief still uses O_viterbi) ====\n")
O_sn <- function(p) { O <- O_mat; O[2,] <- O[2,]*(1-p)/(1-O[2,2]); O[2,2] <- p; O }
mis <- rbind(
  matched = run(ch_warm, N = N500),
  causal_true = run(ch_warm, N = N500, O_true = O_causal),
  causal_true_myopic = run(ch_myo, N = N500, O_true = O_causal),
  SN_acc_0.10 = run(ch_warm, N = N500, O_true = O_sn(.10)),
  SN_acc_0.40 = run(ch_warm, N = N500, O_true = O_sn(.40)),
  sticky_k0.5 = run(ch_warm, N = N500, kappa = .5),
  sticky_k0.8 = run(ch_warm, N = N500, kappa = .8),
  sticky_k0.8_myopic = run(ch_myo, N = N500, kappa = .8))
print(round(mis, 3))
# PBVI re-solved with the causal O, evaluated under causal O
P_c <- pbvi(A, O_causal, R_mat, gamma_disc, B, IT300, vi$V)
cat("PBVI re-solved with causal O, evaluated with causal O:", round(run(function(b) pol_b(b,P_c), N = N500, O_true = O_causal, O_bel = O_causal),3), "\n")

cat("\n==== 7. TEP matrices ====\n")
tep_file <- out_path("tep", "TEP_results.rds")
if (!file.exists(tep_file)) {
  cat("TEP results not found (", tep_file, ") -- run tep/TEP_Application.R first; section skipped.\n")
} else {
tep <- readRDS(tep_file); At <- build_A(tep$A_NoAction, rho); vit <- vi_solve(At, R_mat, gamma_disc)
cat("TEP MDP policy:", vit$policy, "\n"); print(round(tep$A_NoAction,4)); print(round(tep$mcda_confusion,3))
Qt <- lapply(At, function(Aa) (Aa-diag(K))/dt); Pt <- pbvi(At, tep$mcda_confusion, R_mat, gamma_disc, B, IT300, vit$V)
vi_bak <- vi; vi <- vit
rt <- rbind(MDP = run(function(b) vit$policy[which.max(b)], N = N500, Ql = Qt, O_true = diag(K), O_bel = diag(K)),
  POMDP = run(function(b) pol_b(b,Pt), N = N500, Ql = Qt, O_true = tep$mcda_confusion, O_bel = tep$mcda_confusion),
  Myopic = run(ch_myo, N = N500, Ql = Qt, O_true = tep$mcda_confusion, O_bel = tep$mcda_confusion),
  NoAction = run(function(b) "NoAction", N = N500, Ql = Qt, O_true = tep$mcda_confusion, O_bel = tep$mcda_confusion))
vi <- vi_bak; print(round(rt,3))
}



cat("\n==== 8. Slow-repair sweep (repair prob scaled by lambda) ====\n")
for (lam in c(0.01, 0.02, 0.05, 0.1)) {
  Al <- build_A(A_NoAction, rho*lam); vil <- vi_solve(Al, R_mat, gamma_disc)
  Ql <- lapply(Al, function(Aa) (Aa-diag(K))/dt); Pl <- pbvi(Al, O_mat, R_mat, gamma_disc, B, IT200, vil$V)
  vi_bak <- vi; vi <- vil
  r <- rbind(MDP = run(function(b) vil$policy[which.max(b)], N = N300, Ql = Ql, O_true = diag(K), O_bel = diag(K)),
             POMDP = run(function(b) pol_b(b,Pl), N = N300, Ql = Ql),
             Myopic = run(ch_myo, N = N300, Ql = Ql),
             NoAction = run(function(b) "NoAction", N = N300, Ql = Ql))
  vi <- vi_bak
  cat(sprintf("lambda=%.2f  mean recovery (s) for rho in [0.4,0.8]: %.2f-%.2f ; policy: %s\n", lam,
      dt/(0.8*lam), dt/(0.4*lam), paste(vil$policy, collapse=",")))
  print(round(r,3))
}
sink()
