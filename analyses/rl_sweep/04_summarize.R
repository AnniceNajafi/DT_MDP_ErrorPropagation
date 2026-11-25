# Rebuild config_summary.csv from per_run_results.csv (the version written by
# 02_run_sweep.R dropped REINFORCE rows because delta = NA broke interaction()).
local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
out_dir <- dirname(out_path("rl_sweep", "x"))     # output/rl_sweep
rl_core <- repo_path("R", "rl_core.R")
options(width = 250)
p <- read.csv(file.path(out_dir, "per_run_results.csv"))
pb <- read.csv(file.path(out_dir, "pbvi_reference.csv"))
p$delta_key <- ifelse(is.na(p$delta), "-", format(p$delta))
key <- paste(p$algo, p$delta_key, p$lr, p$episodes, sep = "|")
summ <- do.call(rbind, lapply(split(p, key), function(d) {
  b <- which.max(d$eval_mean_T15)
  data.frame(algo = d$algo[1], delta = d$delta[1], n_grid = d$n_grid[1], lr = d$lr[1],
             episodes = d$episodes[1], n_seeds = nrow(d),
             mean_eval_T15 = mean(d$eval_mean_T15), sd_across_seeds_T15 = sd(d$eval_mean_T15),
             median_seed_T15 = median(d$eval_mean_T15),
             min_seed_T15 = min(d$eval_mean_T15), best_seed_T15 = max(d$eval_mean_T15),
             best_seed_id = d$train_seed[b],
             mean_gap_to_pbvi_T15 = mean(d$pbvi_minus_agent_T15),
             n_seeds_sig_worse_than_pbvi = sum(d$pbvi_minus_agent_T15 > 2 * d$paired_se_T15),
             mean_eval_T20 = mean(d$eval_mean_T20), sd_across_seeds_T20 = sd(d$eval_mean_T20),
             best_seed_T20 = max(d$eval_mean_T20),
             mean_train_time_s = mean(d$train_time_s),
             mean_train_last10pct = mean(d$train_last10pct),
             mean_train_prev10pct = mean(d$train_prev10pct),
             mean_last_minus_prev10 = mean(d$train_last_minus_prev))
}))
summ <- summ[order(summ$algo, summ$delta, summ$lr, summ$episodes), ]
write.csv(summ, file.path(out_dir, "config_summary.csv"), row.names = FALSE)
print(summ, digits = 4)

cat("\nPooled over all runs per algorithm:\n")
for (a in unique(p$algo)) {
  d <- p[p$algo == a, ]
  cat(a, ": runs", nrow(d), " median", round(median(d$eval_mean_T15), 2),
      " sig.worse", sum(d$pbvi_minus_agent_T15 > 2 * d$paired_se_T15),
      " sig.better", sum(d$pbvi_minus_agent_T15 < -2 * d$paired_se_T15),
      " best", round(max(d$eval_mean_T15), 2), "\n")
}
# Default-config convergence diagnostic
cat("\nDefault configs, per seed (last 10% vs preceding 10% training return):\n")
dd <- p[(p$algo == "QL" & p$delta_key == "0.1" & p$lr == 0.01) |
        (p$algo == "REINFORCE" & p$lr == 1e-4),
        c("algo","lr","episodes","train_seed","train_prev10pct","train_last10pct","train_last_minus_prev","eval_mean_T15")]
print(dd, digits = 4)
