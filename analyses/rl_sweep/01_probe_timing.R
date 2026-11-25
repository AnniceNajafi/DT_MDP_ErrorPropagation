# 01_probe_timing.R -- optional timing / equivalence probe (not run by run_all.R)
local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
out_dir <- dirname(out_path("rl_sweep", "x"))     # output/rl_sweep
rl_core <- repo_path("R", "rl_core.R")
mo <- readRDS(file.path(out_dir, "model_objects.rds")); list2env(mo, globalenv())
source(rl_core)

# Sanity: generic simulator == POMDPFormulation.R gillespie_pomdp on same seed
L <- readLines(repo_path("R", "POMDPFormulation.R"), encoding = "UTF-8")
s <- grep("^# 4[.] POMDP ACTION", L); e <- grep("^# 6[.] CUMULATIVE", L) - 2
eval(parse(text = L[s:e]))
d <- numeric(20)
for (n in 1:20) {
  set.seed(900000 + n)
  r1 <- gillespie_pomdp(Q_rate, A, pomdp$Gamma, pomdp$actions, O_mat, c_rate, 1, 15, r_disc)$disc_return
  set.seed(900000 + n)
  r2 <- gillespie_belief_policy(Q_rate, pbvi_selector(pomdp$Gamma, pomdp$actions), O_mat, c_rate, 1, 15, r_disc)$disc_return
  d[n] <- r1 - r2
}
cat("max |gillespie_pomdp - generic|:", max(abs(d)), "\n")

t0 <- Sys.time(); q <- train_qlearning(0.1, 0.1, 200, 123); t1 <- Sys.time()
cat("QL 200 ep:", round(as.numeric(difftime(t1, t0, units = "secs")), 2), "s; mean ret", mean(q$returns), "\n")
t0 <- Sys.time(); q <- train_qlearning(0.05, 0.1, 200, 123); t1 <- Sys.time()
cat("QL d=0.05 200 ep:", round(as.numeric(difftime(t1, t0, units = "secs")), 2), "s\n")
t0 <- Sys.time(); r <- train_reinforce(0.001, 200, 456); t1 <- Sys.time()
cat("RE 200 ep:", round(as.numeric(difftime(t1, t0, units = "secs")), 2), "s\n")
t0 <- Sys.time(); ev <- evaluate_policy(pbvi_selector(pomdp$Gamma, pomdp$actions), N_eval = 100); t1 <- Sys.time()
cat("PBVI eval 100:", round(as.numeric(difftime(t1, t0, units = "secs")), 2), "s; mean", mean(ev), "\n")
t0 <- Sys.time(); ev <- evaluate_policy(ql_selector(q$Q_table, q$belief_grid), N_eval = 100); t1 <- Sys.time()
cat("QL eval 100:", round(as.numeric(difftime(t1, t0, units = "secs")), 2), "s\n")
