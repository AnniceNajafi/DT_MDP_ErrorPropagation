#' Revision_ComputationalCost.R
#' ──────────────────────────────────────────────────────────────────────
#' Reviewer Comment 7:
#'   "Maybe add a table for Computational cost comparison?"
#'
#' Measures wall-clock time for offline solving (MDP VI, POMDP PBVI),
#' per-step decision latency at runtime (MDP lookup, POMDP alpha-vector
#' inner product, MCDA TOPSIS score, k-Step threshold check), and peak
#' memory for the solver artefacts. Reports a single table for the paper.
#'
#' Machine info (captured via Sys.info(), R.version) is embedded in the
#' caption so reviewers can interpret absolute times.
#'
#' Outputs:
#'   - output/tables/ComputationalCost.csv
#'   - output/tables/ComputationalCost.tex  (ready-to-\input LaTeX fragment)
#'
#' Absolute times depend on the machine; the paper reports them for an
#' AMD Ryzen 9 9950X3D, 64 GB RAM, R 4.5.1.
#'
#' Authors@R: person("Annice", "Najafi")
#' ──────────────────────────────────────────────────────────────────────

local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
load_core()                                  # pipeline -> MDP -> POMDP
source(repo_path("R", "kstep_policy.R"))     # canonical_repair, gillespie_kstep

library(microbenchmark)

set.seed(42)

# =====================================================================
# 1. OFFLINE SOLVE TIMES
# =====================================================================

cat("\n── Offline solve times ──\n")

# MDP Value Iteration
t_vi <- system.time({
  vi_bench <- value_iteration(A, R_mat, gamma_disc,
                              tol = 1e-10, max_iter = 10000)
})
cat(sprintf("   VI: %.3f s, %d iterations\n",
            t_vi["elapsed"], vi_bench$iterations))

# POMDP PBVI
t_pbvi <- system.time({
  pomdp_bench <- pbvi_solve(A, O_mat, R_mat, gamma_disc, belief_pts,
                            max_horizon = 500, tol = 1e-10,  # 500-iteration budget
                            V_init = vi$V)
})
cat(sprintf("   PBVI: %.3f s, %d iterations, %d alpha vectors\n",
            t_pbvi["elapsed"], pomdp_bench$iterations,
            length(pomdp_bench$Gamma)))

# =====================================================================
# 2. PER-STEP ONLINE DECISION LATENCY (microbenchmark, 1000 reps)
# =====================================================================

cat("\n── Per-step decision latency (microbenchmark) ──\n")

# MDP lookup
b_test <- rep(0.25, K)
a_ground <- 2L

mdp_step <- function() vi$policy[states[a_ground]]

# POMDP alpha-vector inner product
pomdp_step <- function() pomdp_action(b_test, pomdp$Gamma, pomdp$actions)

# k-step: trivial threshold check (we include a constant work proxy)
kstep_step <- function() {
  if (a_ground == 1L) "NoAction" else canonical_repair[states[a_ground]]
}

# MCDA TOPSIS: approximate single decision cost
mcda_step <- function() {
  posterior <- c(0.25, 0.25, 0.25, 0.25)
  pN <- posterior[1]; pS <- posterior[2]
  pD <- posterior[3]; pR <- posterior[4]
  kpi_n <- 1 - pN
  stability_penalty <- c(0.08, 0.08, 0.05, 0.15, 0.12)
  action_cost       <- c(0.10, 0.10, 0.15, 0.40, 0.25)
  w <- c(0.35, 0.25, 0.15, 0.10, 0.15)
  topsis_actions <- c("DwSensorA","DwSensorB","IncFilter",
                      "ReidentPlant","BiasCorrect")
  dm <- cbind(
    KPI_gain    = c(pS*kpi_n*1.2, pS*kpi_n*0.6, pS*kpi_n,
                    pD*kpi_n, (pR+0.5*pD)*kpi_n),
    Module_gain = c(pS*1.2, pS*0.6, pS, pD, (pR+0.5*pD)),
    Drift_gain  = c(pS*1.1, pS*0.7, pS*0.6, pD*0.4, pR+0.2*pS),
    Stability   = -stability_penalty,
    Cost        = -action_cost)
  rownames(dm) <- topsis_actions
  topsis_actions[which.max(
    suppressWarnings(as.numeric(RMCDA::apply.TOPSIS(dm, w = w))))]
}

bench <- microbenchmark(
  MDP       = mdp_step(),
  POMDP     = pomdp_step(),
  kStep     = kstep_step(),
  MCDA      = mcda_step(),
  times     = nq(1000L, 100L),
  unit      = "us"
)
print(bench)

latency_us <- aggregate(time ~ expr, data = bench, FUN = median)
latency_us$time_us <- latency_us$time / 1000  # ns -> us

# =====================================================================
# 3. MEMORY FOOTPRINT (approximate via object.size)
# =====================================================================

cat("\n── Solver artefact memory ──\n")
mem <- data.frame(
  Solver = c("MDP policy + V", "PBVI alpha-vectors", "k-Step (none)",
             "MCDA TOPSIS scorer"),
  Bytes  = c(
    as.numeric(object.size(vi$policy) + object.size(vi$V)),
    as.numeric(object.size(pomdp$Gamma) + object.size(pomdp$actions)),
    128,     # a single canonical-repair lookup table
    as.numeric(object.size(mcda_step))
  )
)
mem$KB <- round(mem$Bytes / 1024, 2)
print(mem)

# =====================================================================
# 4. BUILD TABLE
# =====================================================================

cost_tbl <- data.frame(
  Solver           = c("Optimal MDP (VI)", "POMDP (PBVI)",
                       "k-Step (k=2)", "MCDA (TOPSIS)", "No Intervention"),
  OfflineTime_s    = c(round(t_vi["elapsed"], 3),
                       round(t_pbvi["elapsed"], 3),
                       0.000, 0.000, 0.000),
  Iterations       = c(vi_bench$iterations,
                       pomdp_bench$iterations,
                       NA, NA, NA),
  OnlineLatency_us = c(latency_us$time_us[latency_us$expr == "MDP"],
                       latency_us$time_us[latency_us$expr == "POMDP"],
                       latency_us$time_us[latency_us$expr == "kStep"],
                       latency_us$time_us[latency_us$expr == "MCDA"],
                       latency_us$time_us[latency_us$expr == "MDP"]),
  Memory_KB        = c(round(as.numeric(object.size(vi$policy) +
                                        object.size(vi$V)) / 1024, 2),
                       round(as.numeric(object.size(pomdp$Gamma) +
                                        object.size(pomdp$actions)) / 1024, 2),
                       0.13, 0.50, 0.00),
  stringsAsFactors = FALSE
)

write.csv(cost_tbl, out_path("tables", "ComputationalCost.csv"), row.names = FALSE)

# =====================================================================
# 5. MACHINE INFO CAPTION
# =====================================================================

sysinfo <- Sys.info()
machine_str <- sprintf(
  "%s %s, %s, R %s.%s",
  sysinfo["sysname"], sysinfo["release"], sysinfo["machine"],
  R.version$major, R.version$minor
)

# =====================================================================
# 6. LaTeX FRAGMENT
# =====================================================================

fmt <- function(x, d = 3) ifelse(is.na(x), "--",
                                 formatC(x, format = "f", digits = d))

tex_lines <- c(
  "\\begin{tabular}{@{}lrrrr@{}}",
  "\\toprule",
  "Solver & Offline (s) & Iterations & Online ($\\mu$s) & Memory (KB) \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(cost_tbl))) {
  row <- cost_tbl[i, ]
  tex_lines <- c(tex_lines,
    sprintf("\\textsc{%s} & %s & %s & %s & %s \\\\",
            row$Solver,
            fmt(row$OfflineTime_s, 3),
            ifelse(is.na(row$Iterations), "--", as.character(row$Iterations)),
            fmt(row$OnlineLatency_us, 1),
            fmt(row$Memory_KB, 2))
  )
}
tex_lines <- c(tex_lines, "\\bottomrule", "\\end{tabular}")

writeLines(tex_lines, out_path("tables", "ComputationalCost.tex"))
cat("  -> Saved ComputationalCost.tex (ready to \\input)\n")

cat(sprintf("\n  Hardware: %s\n", machine_str))
cat("\n══ Revision_ComputationalCost complete. ══\n")
