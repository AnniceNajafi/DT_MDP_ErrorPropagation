local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
out_dir <- dirname(out_path("rl_sweep", "x"))     # output/rl_sweep
rl_core <- repo_path("R", "rl_core.R")
options(width = 250)
s <- read.csv(file.path(out_dir, "config_summary.csv"))
print(s[, c("algo","delta","lr","episodes","mean_eval_T15","sd_across_seeds_T15","min_seed_T15","best_seed_T15","mean_gap_to_pbvi_T15","mean_eval_T20","sd_across_seeds_T20","mean_train_time_s","mean_last_minus_prev10")], digits = 4)
p <- read.csv(file.path(out_dir, "per_run_results.csv"))
sel <- p$eval_mean_T15 < 85 | (p$algo == "QL" & p$delta == 0.1 & p$lr == 0.01) | (p$algo == "REINFORCE" & p$lr == 1e-4)
print(p[sel, c("algo","delta","lr","episodes","train_seed","eval_mean_T15","pbvi_minus_agent_T15","paired_se_T15","frac_traj_identical_to_pbvi_T15","grid_cells_visited","n_grid","train_last10pct","train_prev10pct","corner_actions")], digits = 4)
cat("best runs:\n")
print(head(p[order(-p$eval_mean_T15), c("algo","delta","lr","episodes","train_seed","eval_mean_T15","pbvi_minus_agent_T15","paired_se_T15","frac_traj_identical_to_pbvi_T15")], 10), digits = 4)
cat("runs:", nrow(p), " frac within 1 of PBVI:", mean(p$pbvi_minus_agent_T15 < 1),
    " frac significantly worse (gap>2SE):", mean(p$pbvi_minus_agent_T15 > 2 * p$paired_se_T15),
    " frac significantly better:", mean(p$pbvi_minus_agent_T15 < -2 * p$paired_se_T15), "\n")
print(table(p$algo, p$corner_actions))
print(tapply(p$pbvi_minus_agent_T15 > 2 * p$paired_se_T15, p$algo, mean))
print(read.csv(file.path(out_dir, "pbvi_reference.csv")))
cat("PBVI SE of mean T15:", 14.72 / sqrt(500), "\n")
