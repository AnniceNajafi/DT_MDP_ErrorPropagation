# verify.R: the stepwise twin (no actions, canonical schedule, pipeline RNG order) reproduces
# SimulationPipeline::simulate_with_faults_pipe residuals exactly; plus calibration diagnostics.
suppressPackageStartupMessages({ library(depmixS4); library(dplyr) })
local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
cl_code <- repo_path("analyses", "closed_loop")   # scripts
cl_dir  <- dirname(out_path("closed_loop", "x"))  # results (setup.rds, csv, rds)
L <- readLines(repo_path("R", "SimulationPipeline.R"))
eval(parse(text = L[grep("^suppressPackageStartupMessages", L)[1]:(grep("^# 4[.] HELPER", L) - 2)]))
source(file.path(cl_code, "cl_sim.R"))
S <- load_setup(file.path(cl_dir, "setup.rds"), as.numeric(Sys.getenv("CL_DEC_EPS", "1e-3"))); cfg <- cl_cfg_default
epi <- data.frame(type = c(2, 3, 4), k_on = c(251, 501, 726), k_off = c(500, 725, NA))
for (sd in c(11, 111)) {
  ref <- simulate_with_faults_pipe(pipe_params, pipe_schedule, sd)
  W <- make_world(sd, pipe_params, cfg, canonical_noise = TRUE)
  o <- run_twin(W, S, make_ctrl("NoAction", S), cfg, mode = "fixed", episodes = epi, record = TRUE)
  Er <- as.matrix(ref$E[, -1])
  cat("seed", sd, " max|E_step - E_pipeline| per channel:", signif(apply(abs(o$E - Er), 2, max), 3), "\n")
  cat("  ground truth agrees with pipeline schedule:", mean(ST[o$s] == ref$regime_t), "\n")
}
# calibration diagnostics over 200 held-out worlds
gains <- c(); e2a_nom <- c(); refit <- c()
for (sd in 300001:300200) {
  W <- make_world(sd, pipe_params, cfg)
  gains <- c(gains, arx_dc_gain(W$arx_plant$theta) / W$G_true)
  o <- run_twin(W, S, make_ctrl("NoAction", S), cfg, mode = "fixed", episodes = data.frame(type = 2, k_on = 2000, k_off = NA), record = TRUE)
  e2a_nom <- rbind(e2a_nom, c(mean = mean(o$E[501:1000, 3]), sd = sd(o$E[501:1000, 3])))
  for (k0 in c(200, 500, 800)) { idx <- (k0 - cfg$W_rp + 1):k0
    th <- fit_arx_pipe(W$x_pos_P[idx], W$u_act_P[idx], 2L, 2L)$theta
    refit <- c(refit, if (arx_stable(th)) arx_dc_gain(th) / W$G_true else NA) }
}
cat("Calibrated plant ARX DC gain / true gain: quantiles", round(quantile(gains, c(0, .05, .5, .95, 1)), 4), "\n")
cat("  fraction within +-5%:", mean(abs(gains - 1) <= 0.05), "\n")
cat("2-s window refit gain ratio: quantiles", round(quantile(refit, c(0, .05, .5, .95, 1), na.rm = TRUE), 4),
    " unstable:", mean(is.na(refit)), " within 5%:", mean(abs(refit - 1) <= 0.05, na.rm = TRUE), "\n")
cat("Nominal e2a (t in [10,20) s, no fault): mean of per-run means", round(mean(e2a_nom[, 1]), 4),
    " SD of per-run means", round(sd(e2a_nom[, 1]), 4), " mean within-run SD", round(mean(e2a_nom[, 2]), 4), "\n")
