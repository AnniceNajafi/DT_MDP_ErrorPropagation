# cl_setup.R -- closed-loop experiment, step 0
# (1) loads the paper's decision model unchanged (R/POMDPFormulation.R via load_core():
#     PBVI 500 iterations, seed 2024, 185 random belief points)
# (2) fits ONE causal decoder HMM (4-state Gaussian, 10 features) on the pooled
#     training seeds 11-20 (canonical schedule), train-pooled standardisation,
#     Hungarian mapping, best of several EM restarts
# (3) validates an online forward filter against depmixS4::forwardbackward and
#     reports the held-out (seeds 111-120) causal confusion
# Output: output/closed_loop/setup.rds
local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
cl_code <- repo_path("analyses", "closed_loop")   # scripts
cl_dir  <- dirname(out_path("closed_loop", "x"))  # results (setup.rds, csv, rds)
# decision model exactly as in the paper: pipeline -> MDP -> POMDP
# (R/POMDPFormulation.R: seed 2024, 185 random belief points, PBVI 500 iterations, MDP warm start)
suppressPackageStartupMessages({ library(depmixS4); library(clue); library(expm); library(dplyr) })
load_core()
G_pomdp <- do.call(rbind, pomdp$Gamma)
cat("PBVI alpha vectors:", nrow(G_pomdp), "\n")
cat("POMDP corner actions:", sapply(1:K, function(s) pomdp$actions[which.max(G_pomdp[, s])]), "\n")

# ---------------- decoder HMM ----------------
feat_fun <- function(E) { f <- E %>% transmute(e0, e1, e2a, e2b, e3, e4,
  e_abs_sum = abs(e0)+abs(e1)+abs(e2a)+abs(e2b)+abs(e3)+abs(e4),
  e_rms = sqrt((e0^2+e1^2+e2a^2+e2b^2+e3^2+e4^2)/6), e1_diff = c(0,diff(e1)), e3_diff = c(0,diff(e3))); as.data.frame(f) }
train_seeds <- 11:20; test_seeds <- 111:120
tr_sims <- lapply(train_seeds, function(s) simulate_with_faults_pipe(pipe_params, pipe_schedule, s))
Ftr <- do.call(rbind, lapply(tr_sims, function(x) feat_fun(x$E)))
gt_tr <- unlist(lapply(tr_sims, `[[`, "regime_t"))
mu <- colMeans(Ftr); sdv <- apply(Ftr, 2, sd)
clipz <- 6   # standardised features clipped to +-6 (pooled EM underflows on the 20-40 SD spikes at DynamicsOff edges)
Ztr <- as.data.frame(pmin(pmax(scale(Ftr, mu, sdv), -clipz), clipz))
mk_mod <- function(df, nt) depmix(response = lapply(names(df), function(nm) as.formula(paste0(nm, "~1"))),
                                  data = df, nstates = 4, family = rep(list(gaussian()), ncol(df)), ntimes = nt)
best <- NULL
for (r in 1:6) {
  set.seed(1000 + r)
  # depmixS4's default starting values are infeasible on the pooled data, so EM is started from
  # k-means clusters (unsupervised; restarts differ in the k-means seed)
  km <- kmeans(Ztr, 4, nstart = 5, iter.max = 50)
  m0 <- mk_mod(Ztr, rep(1000, length(train_seeds)))
  st0 <- c(rep(0.25, 4), as.numeric(t(matrix(0.01, 4, 4) + diag(0.96, 4))))
  for (s in 1:4) { Zs <- Ztr[km$cluster == s, , drop = FALSE]
    st0 <- c(st0, as.numeric(rbind(colMeans(Zs), pmax(apply(Zs, 2, sd), 0.3)))) }
  m0 <- setpars(m0, st0)
  fm <- tryCatch(fit(m0, verbose = FALSE, emcontrol = em.control(maxit = 300)), error = function(e) { print(e); NULL })
  if (is.null(fm)) next
  vp <- posterior(fm, type = "viterbi")$state
  tab <- as.matrix(table(factor(vp, 1:4), factor(gt_tr, states)))
  perm <- as.integer(solve_LSAP(max(tab) - tab)); map <- states[perm]
  acc <- mean(map[vp] == gt_tr)
  cat(sprintf("restart %d: logLik %.1f  train Viterbi acc %.3f  map %s\n", r, logLik(fm), acc, paste(map, collapse = ",")))
  if (is.null(best) || logLik(fm) > best$ll) best <- list(fm = fm, ll = as.numeric(logLik(fm)), map = map, acc = acc)
}
fm <- best$fm; map <- best$map
# extract parameters
Atr <- t(sapply(1:4, function(i) getpars(fm@transition[[i]])))
init <- getpars(fm@prior)
nf <- ncol(Ztr)
MU <- matrix(NA, 4, nf); SD <- matrix(NA, 4, nf)
for (s in 1:4) for (j in 1:nf) { p <- getpars(fm@response[[s]][[j]]); MU[s, j] <- p[1]; SD[s, j] <- p[2] }
dec <- list(mu = mu, sdv = sdv, Atr = Atr, init = init, MU = MU, SD = SD, map = map,
            map_idx = match(map, states), feat_names = names(Ztr), clipz = clipz)
cat("Decoder transition matrix (HMM-state order):\n"); print(round(Atr, 4)); cat("map:", map, "\n")

# online forward filter (log-domain emissions)
dec_emis_log <- function(dec, x) {  # x standardised feature vector
  rowSums(dnorm(matrix(x, 4, length(x), byrow = TRUE), dec$MU, dec$SD, log = TRUE))
}
fwd_offline <- function(dec, Z) {
  n <- nrow(Z); al <- matrix(0, n, 4); prev <- NULL
  for (k in 1:n) { le <- dec_emis_log(dec, Z[k, ]); pr <- if (k == 1) dec$init else as.numeric(prev %*% dec$Atr)
    w <- log(pmax(pr, 1e-300)) + le; w <- exp(w - max(w)); al[k, ] <- w / sum(w); prev <- al[k, ] }
  al }
conf_c <- matrix(0, 4, 4, dimnames = list(states, states)); agree <- c()
for (s in test_seeds) {
  te <- simulate_with_faults_pipe(pipe_params, pipe_schedule, s)
  Zte <- pmin(pmax(as.matrix(scale(feat_fun(te$E), mu, sdv)), -clipz), clipz)
  al <- fwd_offline(dec, Zte)
  mte <- setpars(mk_mod(as.data.frame(Zte), 1000), getpars(fm))
  fb <- forwardbackward(mte); al2 <- fb$alpha / rowSums(fb$alpha)
  agree <- c(agree, mean(apply(al, 1, which.max) == apply(al2, 1, which.max)), max(abs(al - al2)))
  zc <- map[apply(al, 1, which.max)]
  for (t in seq_along(zc)) conf_c[te$regime_t[t], zc[t]] <- conf_c[te$regime_t[t], zc[t]] + 1
}
cat("Online filter vs depmixS4 forwardbackward: label agreement / max|diff| per seed:\n"); print(matrix(agree, ncol = 2, byrow = TRUE))
O_causal_pooled <- conf_c / rowSums(conf_c)
cat("Held-out causal confusion of the pooled decoder (seeds 111-120):\n"); print(round(O_causal_pooled, 3))
cat("O_mat used by the POMDP (paper):\n"); print(round(O_mat, 3))

saveRDS(list(A = A, A_NoAction = A_NoAction, O_mat = O_mat, R_mat = R_mat, vi = vi, G_pomdp = G_pomdp,
             pomdp_actions = pomdp$actions, states = states, actions = actions, gamma_disc = gamma_disc,
             pipe_params = pipe_params, pipe_schedule = pipe_schedule, dec = dec,
             O_causal_pooled = O_causal_pooled, decoder_train_acc = best$acc, rho = rho,
             mcda_dwell_time = mcda_dwell_time),
        file.path(cl_dir, "setup.rds"))
cat("saved setup.rds\n")
