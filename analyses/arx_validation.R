# arx_validation.R
# ARX surrogate validation on the held-out half of clean (fault-free) runs
# (paper Table "ARX surrogates of the six-module digital twin", held-out fit %).
# Output: output/tables/arx_validation.csv (per seed) and
#         output/tables/arx_validation_summary.csv (mean / min over seeds 11-20).

local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
suppressPackageStartupMessages({ library(dplyr) })

# reuse helper functions only (avoid running the HMM loop)
src <- readLines(repo_path("R", "SimulationPipeline.R"), encoding = "UTF-8")
cut <- grep("^# 3. FULL SIMULATION WITH FAULT INJECTION", src)
eval(parse(text = src[1:(cut - 1)], encoding = "UTF-8"))
p <- pipe_params
fitpct <- function(y, yhat) 100 * (1 - sqrt(sum((y - yhat)^2)) / sqrt(sum((y - mean(y))^2)))
res <- list()
for (seed in 11:20) {
  set.seed(seed)
  n <- as.integer(p$T / p$dt); u_cmd <- make_u_cmd_pipe(n, dt = p$dt)
  ua <- lpf_first_order_pipe(u_cmd, p$tau_cmd, p$dt, p$cmd_proc_std)
  pl <- plant_msd_pipe(ua, p$dt, p$m, p$c, p$k, p$plant_proc_std); x <- pl$pos
  yA <- sensor_first_order_pipe(x, p$dt, p$tau_sa, p$sa_meas_std, p$sa_bias)
  yB <- sensor_first_order_pipe(x, p$dt, p$tau_sb, p$sb_meas_std, p$sb_bias)
  sp <- as.integer(p$train_frac * n); va <- (sp + 1):n
  m0 <- fit_arx_pipe(ua[1:sp], u_cmd[1:sp], p$na_cmd, p$nb_cmd)
  m1 <- fit_arx_pipe(x[1:sp], ua[1:sp], p$na_plant, p$nb_plant)
  m2a <- fit_arx_pipe(yA[1:sp], x[1:sp], p$na_sens, p$nb_sens)
  m2b <- fit_arx_pipe(yB[1:sp], x[1:sp], p$na_sens, p$nb_sens)
  # free-run cascade (surrogates fed by upstream surrogate predictions)
  uaA <- sim_arx_pipe(u_cmd, ua[1:p$na_cmd], m0)
  xA <- sim_arx_pipe(uaA, x[1:p$na_plant], m1)
  yAA <- sim_arx_pipe(xA, yA[1:p$na_sens], m2a); yBA <- sim_arx_pipe(xA, yB[1:p$na_sens], m2b)
  # module-isolated free run (fed by true upstream signal)
  xI <- sim_arx_pipe(ua, x[1:p$na_plant], m1)
  res[[length(res) + 1]] <- data.frame(seed = seed,
    M0 = fitpct(ua[va], uaA[va]), M1_cascade = fitpct(x[va], xA[va]), M1_isolated = fitpct(x[va], xI[va]),
    M2a_cascade = fitpct(yA[va], yAA[va]), M2b_cascade = fitpct(yB[va], yBA[va]))
}
r <- do.call(rbind, res)
print(round(r, 1)); cat("\nmean:\n"); print(round(colMeans(r[,-1]), 1)); cat("min:\n"); print(round(apply(r[,-1], 2, min), 1))
write.csv(r, out_path("tables", "arx_validation.csv"), row.names = FALSE)
write.csv(data.frame(module = names(r)[-1], mean_fit_pct = colMeans(r[, -1]),
                     min_fit_pct = apply(r[, -1], 2, min)),
          out_path("tables", "arx_validation_summary.csv"), row.names = FALSE)
