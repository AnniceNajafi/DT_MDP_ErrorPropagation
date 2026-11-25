# cl_repair.R -- empirical repair times / per-step repair probabilities (Reviewer 2 #6)
# For each (action, fault) pair: fault injected at t = 6 s (k_on = 301) and persisting until repaired;
# the action is held from k_start until the episode is repaired (then NoAction).
# k_start = k_on ("onset") or k_on + U(0,2) s ("random").  N = 200 held-out seeds (400001-400200).
suppressPackageStartupMessages({ library(parallel) })
local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
cl_code <- repo_path("analyses", "closed_loop")   # scripts
cl_dir  <- dirname(out_path("closed_loop", "x"))  # results (setup.rds, csv, rds)
source(file.path(cl_code, "cl_sim.R"))   # ST, ACT, cfg in the master
tag <- Sys.getenv("CL_TAG", "")   # output-file suffix for sensitivity runs
N <- as.integer(Sys.getenv("CL_N_REPAIR", "200")); seeds <- 400000L + seq_len(N)
pairs <- data.frame(
  label  = c("DwSensorA|SensorNoisy", "DwSensorB|SensorNoisy", "IncFilter|SensorNoisy", "ReidentPlant|DynamicsOff",
             "BiasCorrect|Drift", "DwSensorA|SensorNoisy(A-only noise)",
             "NoAction|SensorNoisy", "NoAction|DynamicsOff", "NoAction|Drift", "NoAction|SensorNoisy(A-only noise)"),
  action = c(2, 3, 4, 5, 6, 2, 1, 1, 1, 1), fault = c(2, 2, 2, 3, 4, 2, 2, 3, 4, 2),
  rho_assumed = c(0.6, 0.6, 0.4, 0.7, 0.8, 0.6, NA, NA, NA, NA),
  a_only = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE))
k_on <- 301L

one <- function(seed, pairs, k_on) {
  P <- S$pipe_params
  out <- list()
  for (i in seq_len(nrow(pairs))) for (mode in c("onset", "random")) {
    cfg <- cfg_override(cl_cfg_default); if (pairs$a_only[i]) cfg$sn_sdB <- 0
    W <- make_world(seed, P, cfg)
    k_start <- if (mode == "onset") k_on else k_on + round(W$delay / cfg$dt)
    ctrl <- make_ctrl("Scripted", S, fixed = list(a = pairs$action[i], k_start = k_start))
    o <- run_twin(W, S, ctrl, cfg, mode = "fixed", episodes = data.frame(type = pairs$fault[i], k_on = k_on, k_off = NA), record = TRUE)
    e <- o$eps[1, ]
    rep_ok <- e$how == "repaired"
    T_rep <- if (rep_ok) (e$k_end - k_start) * cfg$dt else NA
    steps_risk <- if (rep_ok) e$k_end - k_start else W$n - k_start + 1
    # decoder: first k >= k_start starting a run of `sustain` consecutive Nominal labels
    zN <- o$z == 1L; kk <- NA
    for (k in k_start:(W$n - cfg$sustain + 1)) if (all(zN[k:(k + cfg$sustain - 1)])) { kk <- k; break }
    t_dec <- if (is.na(kk)) NA else (kk - k_start) * cfg$dt
    dec_rec_1s <- NA
    if (rep_ok && e$k_end <= W$n - cfg$sustain + 1) { dec_rec_1s <- FALSE
      for (k in e$k_end:min(e$k_end + 50L, W$n - cfg$sustain + 1)) if (all(zN[k:(k + cfg$sustain - 1)])) { dec_rec_1s <- TRUE; break } }
    dec_nom_frac_after_rep <- if (rep_ok && e$k_end + 49 <= W$n) mean(zN[e$k_end:(e$k_end + 49)]) else NA
    relapse <- if (rep_ok) any(o$s[e$k_end:W$n] != 1L) else NA
    out[[length(out) + 1]] <- data.frame(seed = seed, pair = pairs$label[i], mode = mode, k_start = k_start,
      repaired = rep_ok, T_rep = T_rep, steps_risk = steps_risk, t_dec = t_dec, dec_rec_1s = dec_rec_1s, dec_nom_frac_after_rep = dec_nom_frac_after_rep,
      dec_nom_frac_1s = mean(zN[k_start:(k_start + 49)]), relapse = relapse,
      rp_ok = o$rp_fits[["ok"]], rp_fail = o$rp_fits[["fail"]])
  }
  do.call(rbind, out)
}

ncore <- as.integer(Sys.getenv("DTMDP_CORES", min(8L, max(1L, detectCores() - 2L))))
cl <- makeCluster(ncore)
clusterExport(cl, c("cl_dir", "cl_code"))
invisible(clusterCall(cl, function(e) Sys.setenv(CL_CFG = e), Sys.getenv("CL_CFG", "")))
invisible(clusterEvalQ(cl, { source(file.path(cl_code, "cl_sim.R")); S <- load_setup(file.path(cl_dir, "setup.rds"), as.numeric(Sys.getenv("CL_DEC_EPS", "1e-3"))); NULL }))
t0 <- Sys.time()
res <- do.call(rbind, parLapply(cl, seeds, one, pairs = pairs, k_on = k_on))
stopCluster(cl)
cat("elapsed:", format(Sys.time() - t0), "\n")
write.csv(res, file.path(cl_dir, paste0("repair_runs", tag, ".csv")), row.names = FALSE)

dt <- 0.02
summ <- do.call(rbind, lapply(split(res, list(res$pair, res$mode), drop = TRUE), function(d) {
  p <- pairs[pairs$label == d$pair[1], ]
  rep <- d$repaired
  mT <- if (any(rep)) mean(d$T_rep[rep]) else NA
  rel <- d$t_dec - d$T_rep
  data.frame(pair = d$pair[1], mode = d$mode[1], n = nrow(d), frac_repaired = mean(rep),
    mean_T_rep_s = mT, sd_T_rep_s = if (sum(rep) > 1) sd(d$T_rep[rep]) else NA,
    median_T_rep_s = if (any(rep)) median(d$T_rep[rep]) else NA,
    rho_hat_dt_over_mean = if (any(rep)) dt / mT else 0,
    rho_hat_censored_MLE = sum(rep) / sum(d$steps_risk),
    rho_assumed = p$rho_assumed,
    mean_recovery_assumed_s = dt / p$rho_assumed,
    relapse_frac = if (any(rep)) mean(d$relapse[rep]) else NA,
    dec_nominal_within_1s_of_start = mean(!is.na(d$t_dec) & d$t_dec <= 1),
    dec_nominal_within_1s_of_repair = if (any(rep)) mean(d$dec_rec_1s[rep]) else NA,
    dec_nominal_before_repair = if (any(rep)) mean(!is.na(rel[rep]) & rel[rep] < 0) else NA,
    dec_nominal_ever = mean(!is.na(d$t_dec)),
    median_t_dec_s = median(d$t_dec, na.rm = TRUE),
    dec_nom_frac_first_1s = mean(d$dec_nom_frac_1s),
    dec_nom_frac_1s_after_repair = if (any(rep)) mean(d$dec_nom_frac_after_rep[rep], na.rm = TRUE) else NA,
    rp_fail_per_run = mean(d$rp_fail))
}))
summ <- summ[order(summ$mode, summ$pair), ]
write.csv(summ, file.path(cl_dir, paste0("repair_summary", tag, ".csv")), row.names = FALSE)
options(width = 250); print(summ, digits = 4, row.names = FALSE)
