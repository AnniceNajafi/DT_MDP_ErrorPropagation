#' hmm_bic.R
#' BIC of the residual HMM for K = 2,...,7 states on the ten seeds of
#' R/SimulationPipeline.R (same simulation, features and depmixS4 fit).
#' Run from the repository root:  Rscript analyses/hmm_bic.R
#' Writes output/analyses/hmm_bic.csv.
suppressPackageStartupMessages({ library(dplyr); library(depmixS4) })
source(file.path(getwd(), "R", "paths.R"))

# Function definitions, parameters and fault schedule from the pipeline
# (top-level assignments only; the Monte Carlo loop is not run).
exprs <- parse(repo_path("R", "SimulationPipeline.R"))
keep <- vapply(exprs, function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
                 is.name(e[[2]]) && (grepl("_pipe$", as.character(e[[2]])) ||
                 as.character(e[[2]]) %in% c("pipe_params", "pipe_schedule", "states_pipeline")),
               logical(1))
for (e in exprs[keep]) eval(e, globalenv())

features <- function(E) {
  E %>% transmute(
    e0, e1, e2a, e2b, e3, e4,
    e_abs_sum = abs(e0) + abs(e1) + abs(e2a) + abs(e2b) + abs(e3) + abs(e4),
    e_rms     = sqrt((e0^2 + e1^2 + e2a^2 + e2b^2 + e3^2 + e4^2) / 6),
    e1_diff   = c(0, diff(e1)),
    e3_diff   = c(0, diff(e3))) %>%
    mutate(across(everything(), ~ as.numeric(scale(.x)))) %>% as.data.frame()
}

seeds <- 11L + 0:9
Ks <- 2:7
res <- expand.grid(K = Ks, seed = seeds)
res$BIC <- NA_real_
for (i in seq_len(nrow(res))) {
  sim <- simulate_with_faults_pipe(pipe_params, pipe_schedule, res$seed[i])  # sets the seed
  Ez <- features(sim$E)
  mod <- depmix(lapply(names(Ez), function(nm) as.formula(paste0(nm, " ~ 1"))),
                data = Ez, nstates = res$K[i], family = rep(list(gaussian()), ncol(Ez)))
  fm <- tryCatch(fit(mod, verbose = FALSE), error = function(e) NULL)
  if (!is.null(fm)) res$BIC[i] <- BIC(fm)
}

out <- out_path("analyses", "hmm_bic.csv")
write.csv(res, out, row.names = FALSE)

best <- res %>% group_by(seed) %>% slice_min(BIC, n = 1, with_ties = FALSE)
print(res %>% group_by(K) %>% summarise(mean_BIC = mean(BIC, na.rm = TRUE),
                                         n_fit = sum(!is.na(BIC))))
print(table(best_K = best$K))
