# cl_closed_loop.R -- closed-loop controller comparison on the six-module twin (Reviewer 2 #8)
# Held-out seeds 500001..500000+N (never used for ARX/HMM decoder fitting: decoder trained on 11-20).
# Common random numbers: run n of every controller uses the same physical asset, ARX twin calibration,
# noise sequences and fault-type / inter-fault-gap draws (make_world(seed)).
suppressPackageStartupMessages({ library(parallel) })
local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
cl_code <- repo_path("analyses", "closed_loop")   # scripts
cl_dir  <- dirname(out_path("closed_loop", "x"))  # results (setup.rds, csv, rds)
source(file.path(cl_code, "cl_sim.R"))   # ST, ACT, cfg in the master
tag <- Sys.getenv("CL_TAG", "")   # output-file suffix for sensitivity runs
N <- as.integer(Sys.getenv("CL_N", "500")); seeds <- 500000L + seq_len(N)
ctrls <- c("Oracle", "POMDP", "Myopic", "MCDA", "NoAction")
variant <- Sys.getenv("CL_VARIANT", "main")   # "main" or "aonly" (SensorNoisy noise on sensor A only)

one <- function(seed, ctrls, variant) {
  cfg <- cfg_override(cl_cfg_default); if (variant == "aonly") cfg$sn_sdB <- 0
  W <- make_world(seed, S$pipe_params, cfg)
  rows <- list(); confs <- list()
  for (cn in ctrls) {
    o <- run_twin(W, S, make_ctrl(cn, S), cfg, mode = "random")
    e <- o$eps
    rt <- function(ty) { d <- e[e$type == ty & e$how == "repaired", ]; if (nrow(d)) mean(d$k_end - d$k_on) * cfg$dt else NA }
    rows[[cn]] <- data.frame(seed = seed, controller = cn, ret = o$ret, fnom = o$fnom, mism = o$mism, acc_dec = o$acc_dec,
      n_ep = nrow(e), n_rep = sum(e$how == "repaired"), n_self = sum(e$how == "selfclear"), n_cens = sum(e$how == "censored"),
      n_SN = sum(e$type == "SensorNoisy"), n_DO = sum(e$type == "DynamicsOff"), n_DR = sum(e$type == "Drift"),
      rep_SN = sum(e$type == "SensorNoisy" & e$how == "repaired"), rep_DO = sum(e$type == "DynamicsOff" & e$how == "repaired"),
      rep_DR = sum(e$type == "Drift" & e$how == "repaired"),
      Trep_SN = rt("SensorNoisy"), Trep_DO = rt("DynamicsOff"), Trep_DR = rt("Drift"),
      rp_ok = o$rp_fits[["ok"]], rp_fail = o$rp_fits[["fail"]], bc_installs = o$bc_installs,
      t(setNames(o$act_frac, paste0("frac_", ACT))))
    confs[[cn]] <- unclass(o$conf)
  }
  list(rows = do.call(rbind, rows), confs = confs)
}

ncore <- as.integer(Sys.getenv("DTMDP_CORES", min(10L, max(1L, detectCores() - 2L))))
cl <- makeCluster(ncore)
clusterExport(cl, c("cl_dir", "cl_code"))
invisible(clusterCall(cl, function(e) Sys.setenv(CL_CFG = e), Sys.getenv("CL_CFG", "")))
invisible(clusterEvalQ(cl, { source(file.path(cl_code, "cl_sim.R")); S <- load_setup(file.path(cl_dir, "setup.rds"), as.numeric(Sys.getenv("CL_DEC_EPS", "1e-3"))); NULL }))
t0 <- Sys.time()
out <- parLapply(cl, seeds, one, ctrls = ctrls, variant = variant)
stopCluster(cl)
cat("elapsed:", format(Sys.time() - t0), "\n")
runs <- do.call(rbind, lapply(out, `[[`, "rows"))
sfx <- paste0(if (variant == "main") "" else paste0("_", variant), tag)
write.csv(runs, file.path(cl_dir, paste0("closed_loop_runs", sfx, ".csv")), row.names = FALSE)

lab <- c(Oracle = "Oracle MDP (true regime)", POMDP = "POMDP (PBVI)", Myopic = "Myopic belief",
         MCDA = "MCDA (TOPSIS)", NoAction = "No Intervention")
summ <- do.call(rbind, lapply(ctrls, function(cn) { d <- runs[runs$controller == cn, ]
  data.frame(Controller = lab[cn], N = nrow(d), MeanReturn = mean(d$ret), SDReturn = sd(d$ret), SEReturn = sd(d$ret) / sqrt(nrow(d)),
    FracNominal = mean(d$fnom), MismatchPct = 100 * mean(d$mism), DecoderAcc = mean(d$acc_dec),
    FaultsPerRun = mean(d$n_ep), FracRepaired = sum(d$n_rep) / sum(d$n_ep), FracSelfCleared = sum(d$n_self) / sum(d$n_ep),
    RepFrac_SN = sum(d$rep_SN) / sum(d$n_SN), RepFrac_DO = sum(d$rep_DO) / sum(d$n_DO), RepFrac_DR = sum(d$rep_DR) / sum(d$n_DR),
    MeanTrep_SN = mean(d$Trep_SN, na.rm = TRUE), MeanTrep_DO = mean(d$Trep_DO, na.rm = TRUE), MeanTrep_DR = mean(d$Trep_DR, na.rm = TRUE),
    RefitFailPerRun = mean(d$rp_fail),
    t(colMeans(d[, grep("^frac_", names(d))]))) }))
summ <- summ[order(-summ$MeanReturn), ]
write.csv(summ, file.path(cl_dir, paste0("closed_loop_summary", sfx, ".csv")), row.names = FALSE)

# paired differences (common random numbers) between adjacent controllers in the ranking + vs POMDP
R <- sapply(ctrls, function(cn) runs$ret[runs$controller == cn])
cliff <- function(x, y) mean(sign(outer(x, y, "-")))
pr <- list(); ord <- match(names(lab)[match(summ$Controller, lab)], ctrls)
cmp <- unique(rbind(cbind(ord[-length(ord)], ord[-1]), cbind(match("POMDP", ctrls), match(c("Myopic", "MCDA", "NoAction"), ctrls))))
for (i in seq_len(nrow(cmp))) { a <- cmp[i, 1]; b <- cmp[i, 2]; dlt <- R[, a] - R[, b]
  pr[[i]] <- data.frame(Comparison = paste(ctrls[a], "-", ctrls[b]), MeanDiff = mean(dlt), SEDiff = sd(dlt) / sqrt(length(dlt)),
    WilcoxPairedP = suppressWarnings(wilcox.test(R[, a], R[, b], paired = TRUE)$p.value), CliffDelta = cliff(R[, a], R[, b])) }
pr <- do.call(rbind, pr)
write.csv(pr, file.path(cl_dir, paste0("closed_loop_paired", sfx, ".csv")), row.names = FALSE)

conf <- lapply(ctrls, function(cn) { M <- Reduce(`+`, lapply(out, function(x) x$confs[[cn]])); dimnames(M) <- list(ST, ST); M })
names(conf) <- ctrls
saveRDS(list(runs = runs, summary = summ, paired = pr, decoder_confusion_counts = conf,
             cfg = cfg_override(cl_cfg_default), cfg_spec = Sys.getenv("CL_CFG", ""), dec_eps = Sys.getenv("CL_DEC_EPS", "1e-3"), variant = variant, seeds = seeds),
        file.path(cl_dir, paste0("closed_loop_results", sfx, ".rds")))
options(width = 250)
print(summ, digits = 5, row.names = FALSE); print(pr, digits = 4, row.names = FALSE)
for (cn in ctrls) { cat("\nDecoder confusion (rows true, cols label), controller", cn, "\n"); print(round(conf[[cn]] / rowSums(conf[[cn]]), 3)) }
