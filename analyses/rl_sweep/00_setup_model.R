# 00_setup_model.R -- builds the model objects for the RL sweep
# (Supplementary "Tuning, budget, and variance") and saves them to
# output/rl_sweep/model_objects.rds. PBVI: V_init = vi$V, max_horizon = 500.
local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
out_dir <- dirname(out_path("rl_sweep", "x"))     # output/rl_sweep
rl_core <- repo_path("R", "rl_core.R")
suppressPackageStartupMessages({ library(depmixS4); library(clue) })

# Model objects from R/POMDPFormulation.R (pipeline -> MDP -> POMDP); PBVI is
# re-solved below (identical settings) only to record its wall-clock time.
load_core()

set.seed(2024)  # belief points for PBVI (Dirichlet samples)
belief_pts <- generate_belief_points(K, n_random = 185)
t0 <- Sys.time()
pomdp <- pbvi_solve(A, O_mat, R_mat, gamma_disc, belief_pts,
                    max_horizon = 500, tol = 0.001, V_init = vi$V)
pbvi_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat("PBVI:", pomdp$iterations, "iterations,", length(pomdp$Gamma),
    "alpha vectors,", round(pbvi_time, 1), "s\n")

saveRDS(list(A = A, R_mat = R_mat, gamma_disc = gamma_disc, Q_rate = Q_rate,
             c_rate = c_rate, r_disc = r_disc, O_mat = O_mat, vi = vi,
             K = K, M = M, states = states, actions = actions,
             pomdp = pomdp, pbvi_time = pbvi_time),
        file.path(out_dir, "model_objects.rds"))
cat("Saved model objects.\n")
print(round(O_mat, 3)); print(vi$policy)
