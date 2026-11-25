# cl_sim.R -- time-stepped six-module simulation with interventions acting on the twin,
# online causal HMM decoder, and closed-loop controllers.
# Physics/ARX helper functions are copied verbatim from ErrorPropagation_MDP/SimulationPipeline.R.

# ------------------------------------------------------------------ copied helpers
make_u_cmd_pipe <- function(n, amp = 1.0, dt = 0.02) {
  t <- (0:(n-1)) * dt
  u <- 0.6 * sin(2*pi*0.25*t) + 0.4 * sign(sin(2*pi*0.05*t))
  step_every <- as.integer(2.5/dt)
  for (i in seq(1, n, by = step_every)) u[i:n] <- u[i:n] + 0.2 * rnorm(1)
  amp * u
}
lpf_first_order_pipe <- function(u, tau, dt, proc_std = 0.0) {
  n <- length(u); y <- numeric(n)
  a <- exp(-dt / max(tau, 1e-6))
  for (k in 2:n) y[k] <- a*y[k-1] + (1 - a)*u[k-1] + proc_std*rnorm(1)
  y
}
plant_msd_pipe <- function(u, dt, m, c, k, proc_std = 0.0) {
  n <- length(u); x <- numeric(n); v <- numeric(n)
  for (kk in 2:n) {
    acc <- (u[kk-1] - c*v[kk-1] - k*x[kk-1]) / m + proc_std*rnorm(1)
    v[kk] <- v[kk-1] + dt*acc
    x[kk] <- x[kk-1] + dt*v[kk-1]
  }
  list(pos = x, vel = v)
}
sensor_first_order_pipe <- function(xpos, dt, tau, meas_std = 0.0, bias = 0.0) {
  y <- lpf_first_order_pipe(xpos, tau, dt, 0)
  y + bias + rnorm(length(y), sd = meas_std)
}
fit_arx_pipe <- function(y, u, na, nb) {
  stopifnot(length(y) == length(u))
  kmax <- length(y); p <- max(na, nb)
  if (kmax <= p) stop("ARX series too short.")
  rows <- vector("list", kmax - p); Tvec <- numeric(kmax - p)
  for (k in (p+1):kmax) {
    phi <- c(sapply(1:na, function(i) y[k-i]), sapply(1:nb, function(j) u[k-j]))
    rows[[k - p]] <- phi; Tvec[k - p] <- y[k]
  }
  Phi <- do.call(rbind, rows)
  theta <- as.numeric(solve(t(Phi) %*% Phi, t(Phi) %*% Tvec))
  list(theta = theta, na = na, nb = nb)
}
sim_arx_pipe <- function(u, y_seed, model) {
  na <- model$na; nb <- model$nb; theta <- model$theta
  n <- length(u); yhat <- numeric(n)
  yhat[1:na] <- tail(y_seed, na); p <- max(na, nb)
  for (k in (p+1):n) {
    phi <- c(sapply(1:na, function(i) yhat[k-i]), sapply(1:nb, function(j) u[k-j]))
    yhat[k] <- sum(theta * phi)
  }
  yhat
}

# ------------------------------------------------------------------ configuration
cl_cfg_default <- list(
  dt = 0.02, n = 1000L,
  alpha_nom = 0.6, alpha_dwA = 0.1, alpha_dwB = 0.9,   # fusion weight on sensor A (DwSensorB -> weight on B = 0.1)
  tau_f = 0.10,                                        # IncFilter: extra 1st-order LPF on the twin sensor channels (s)
  sn_sdA = 0.06, sn_sdB = 0.05,                        # SensorNoisy fault (as pipeline)
  do_scale = 0.8,                                      # DynamicsOff fault (as pipeline)
  drift_rate = 0.0007, drift_k0 = 725,                 # Drift: 0.0007*(k - k_on + 725); identical to pipeline for k_on = 726
  W_rp = 100L,                                         # ReidentPlant data window (samples) = 2 s
  rp_gain_tol = 0.05,                                  # accept refit iff |G_hat/G_true - 1| <= 5 % and AR part stable
  W_bc = 50L,                                          # BiasCorrect data window (samples) = 1 s
  tol_b = 0.05,                                        # Drift present iff |net sensor-A offset error| > tol_b
  W_sn = 25L, thr_sn = 0.02,                           # SensorNoisy repaired iff RMS of fault-induced fused noise over the
                                                       #   last W_sn samples (one constant non-nominal config) <= thr_sn
  Dmax = 400L,                                         # max fault duration (samples) -> self-clears (8 s)
  gap_lo = 1, gap_hi = 4,                              # Nominal sojourn before/between faults ~ U(gap_lo, gap_hi) s
  sustain = 25L                                        # decoder "back to Nominal" = 25 consecutive Nominal labels (0.5 s)
)
ST <- c("Nominal", "SensorNoisy", "DynamicsOff", "Drift")
ACT <- c("NoAction", "DwSensorA", "DwSensorB", "IncFilter", "ReidentPlant", "BiasCorrect")

# ------------------------------------------------------------------ per-seed physical asset + ARX twin calibration
make_world <- function(seed, P, cfg, canonical_noise = FALSE) {
  set.seed(seed)
  n <- cfg$n; dt <- cfg$dt
  u_cmd <- make_u_cmd_pipe(n, dt = dt)
  u_act_P <- lpf_first_order_pipe(u_cmd, P$tau_cmd, dt, P$cmd_proc_std)
  pl <- plant_msd_pipe(u_act_P, dt, P$m, P$c, P$k, P$plant_proc_std)
  x_pos_P <- pl$pos; x_vel_P <- pl$vel
  yA_P <- sensor_first_order_pipe(x_pos_P, dt, P$tau_sa, P$sa_meas_std, P$sa_bias)
  yB_P <- sensor_first_order_pipe(x_pos_P, dt, P$tau_sb, P$sb_meas_std, P$sb_bias)
  y_fused_P <- P$alpha * yA_P + (1 - P$alpha) * yB_P
  z_P <- P$w1 * x_pos_P + P$w2 * x_vel_P
  split <- as.integer(P$train_frac * n)
  arx_cmd   <- fit_arx_pipe(u_act_P[1:split], u_cmd[1:split], P$na_cmd, P$nb_cmd)
  arx_plant <- fit_arx_pipe(x_pos_P[1:split], u_act_P[1:split], P$na_plant, P$nb_plant)
  arx_sa    <- fit_arx_pipe(yA_P[1:split], x_pos_P[1:split], P$na_sens, P$nb_sens)
  arx_sb    <- fit_arx_pipe(yB_P[1:split], x_pos_P[1:split], P$na_sens, P$nb_sens)
  u_act_A <- sim_arx_pipe(u_cmd, y_seed = u_act_P[1:P$na_cmd], arx_cmd)
  if (canonical_noise) {   # reproduce the pipeline's RNG order (noise drawn only on SN samples 251..500)
    nA <- nB <- numeric(n); for (k in 251:500) { nA[k] <- rnorm(1, 0, 0.06); nB[k] <- rnorm(1, 0, 0.05) }
    nA <- nA / 0.06; nB <- nB / 0.05                         # stored as unit-variance noise
    gaps <- types <- NULL; delay <- 0
  } else {
    nA <- rnorm(n); nB <- rnorm(n)
    gaps <- runif(30, cfg$gap_lo, cfg$gap_hi); types <- sample(2:4, 30, replace = TRUE)
    delay <- runif(1, 0, 2)
  }
  list(n = n, u_act_P = u_act_P, x_pos_P = x_pos_P, yA_P = yA_P, yB_P = yB_P, y_fused_P = y_fused_P, z_P = z_P,
       u_act_A = u_act_A, arx_plant = arx_plant, arx_sa = arx_sa, arx_sb = arx_sb,
       G_true = 1 / P$k, nA = nA, nB = nB, gaps = gaps, types = types, delay = delay)
}

arx_dc_gain <- function(th) sum(th[3:4]) / (1 - sum(th[1:2]))
arx_stable <- function(th) all(Mod(polyroot(c(1, -th[1], -th[2]))) > 1)   # roots of 1 - a1 z - a2 z^2 outside unit circle

# ------------------------------------------------------------------ decoder (online forward filter)
dec_step <- function(dec, x, prev) {
  xs <- pmin(pmax((x - dec$mu) / dec$sdv, -dec$clipz), dec$clipz)
  le <- rowSums(dnorm(matrix(xs, 4, length(xs), byrow = TRUE), dec$MU, dec$SD, log = TRUE))
  pr <- if (is.null(prev)) dec$init else as.numeric(prev %*% dec$Atr)
  w <- log(pmax(pr, 1e-300)) + le; w <- exp(w - max(w)); w / sum(w)
}

# ------------------------------------------------------------------ controllers
# controller = list(init = function() state, act = function(st, z, s_true, k, info) list(a = action index, st = state))
make_ctrl <- function(type, S, fixed = NULL) {
  A <- S$A; O <- S$O_mat; R <- S$R_mat; G <- S$G_pomdp; lab <- S$pomdp_actions; pol <- S$vi$policy
  bel_upd <- function(b, a, z) { b <- as.numeric(b %*% A[[a]]) * O[, z]; if (sum(b) > 0) b / sum(b) else rep(.25, 4) }
  switch(type,
    NoAction = list(init = function() NULL, act = function(st, z, s, k, info) list(a = 1L, st = st)),
    Oracle   = list(init = function() NULL, act = function(st, z, s, k, info) list(a = match(pol[s], ACT), st = st)),
    POMDP = list(init = function() list(b = c(1, 0, 0, 0), a = 1L),
      act = function(st, z, s, k, info) {
        b <- bel_upd(st$b, st$a, z); a <- match(lab[which.max(G %*% b)], ACT); list(a = a, st = list(b = b, a = a)) }),
    Myopic = list(init = function() list(b = c(1, 0, 0, 0), a = 1L),
      act = function(st, z, s, k, info) {
        b <- bel_upd(st$b, st$a, z); a <- as.integer(which.max(R %*% b)); list(a = a, st = list(b = b, a = a)) }),
    MCDA = {
      ta <- c(2L, 3L, 4L, 5L, 6L)
      stab <- c(0.08, 0.08, 0.05, 0.15, 0.12); cost <- c(0.10, 0.10, 0.15, 0.40, 0.25)
      wts <- c(0.35, 0.25, 0.15, 0.10, 0.15); dwell <- S$mcda_dwell_time; aema <- 0.3
      topsis <- function(Am, w) {   # identical to RMCDA::apply.TOPSIS
        nA <- t(t(Am) / sqrt(colSums(Am^2))); nA <- nA * w
        Sp <- sqrt(rowSums((t(t(nA) - apply(nA, 2, max)))^2)); Sn <- sqrt(rowSums((t(t(nA) - apply(nA, 2, min)))^2))
        Sn / (Sp + Sn) }
      list(init = function() list(post = NULL, cur = 1L, last = 0),
        act = function(st, z, s, k, info) {
          t_now <- (k - 1) * 0.02
          if (is.null(st$post)) { p <- rep(0, 4); p[z] <- 1; st$post <- p; return(list(a = st$cur, st = st)) }
          ind <- rep(0, 4); ind[z] <- 1; p <- aema * ind + (1 - aema) * st$post; p <- p / sum(p); st$post <- p
          if ((t_now - st$last) >= dwell) {   # desired action only matters when the dwell allows a switch
            pN <- p[1]; pS <- p[2]; pD <- p[3]; pR <- p[4]; sn <- pS; pn <- pD; dn <- pR; kn <- 1 - pN; nn <- sn
            kg <- c(pS*kn*sn*1.2, pS*kn*sn*0.6, pS*kn*nn, pD*kn*pn, (pR+0.5*pD)*kn*(dn+0.3*pn))
            mg <- c(pS*sn*1.2, pS*sn*0.6, pS*sn*nn, pD*pn, (pR+0.5*pD)*(dn+0.3*pn))
            dg <- c(pS*dn*1.1, pS*dn*0.7, pS*dn*0.6, pD*dn*0.4, (pR+0.2*pS)*dn)
            dm <- cbind(kg, mg, dg, -stab, -cost)
            sc <- suppressWarnings(topsis(dm, wts))
            des <- if (all(is.na(sc)) || all(sc == 0, na.rm = TRUE)) 1L else ta[which.max(sc)]
            if (kn < 0.25 && pN > 0.75) des <- 1L
            if (des != st$cur) { st$cur <- des; st$last <- t_now }
          }
          list(a = st$cur, st = st) })
    },
    Scripted = list(init = function() NULL,     # fixed$a held from fixed$k_start until the fault episode is repaired
      act = function(st, z, s, k, info) {
        a <- if (k >= fixed$k_start && !info$episode_done) fixed$a else 1L
        list(a = a, st = st) })
  )
}

# ------------------------------------------------------------------ main simulator
# mode "random": randomised schedule (W$gaps, W$types); faults persist until repaired or Dmax (self-clear)
# mode "fixed":  episodes = data.frame(type, k_on, k_off); k_off = NA -> persists until repaired / end of run
run_twin <- function(W, S, ctrl, cfg, mode = "random", episodes = NULL, record = FALSE) {
  n <- W$n; dt <- cfg$dt; g <- S$gamma_disc; R <- S$R_mat; pol_idx <- match(S$vi$policy, ACT); dec <- S$dec
  th_p <- W$arx_plant$theta; th_a <- W$arx_sa$theta; th_b <- W$arx_sb$theta
  xin <- numeric(n); yAi <- numeric(n); yBi <- numeric(n)
  xin[1:2] <- W$x_pos_P[1:2]; yAi[1] <- W$yA_P[1]; yBi[1] <- W$yB_P[1]
  a_f <- exp(-dt / cfg$tau_f)
  sn_on <- FALSE; do_on <- FALSE; drift_src <- integer(0); corr_int <- 0; corr_slope <- 0
  ep <- NULL; ep_log <- list(); gi <- 1L
  next_onset <- if (mode == "random") 1L + round(W$gaps[1] / dt) else NA
  a_prev <- 1L; rp_cnt <- 0L; bc_cnt <- 0L; cfg_since <- 1L; cfg_key_prev <- "0.6 FALSE"
  fA <- 0; fB <- 0; nfA <- 0; nfB <- 0; nf_hist <- numeric(n)
  e2a_hist <- numeric(n); x_rep_prev <- 0; e1p <- 0; e3p <- 0
  post <- NULL; cst <- ctrl$init()
  ret <- 0; s_seq <- integer(n); a_seq <- integer(n); z_seq <- integer(n); off_seq <- numeric(n)
  rp_fits <- c(ok = 0L, fail = 0L); bc_installs <- 0L
  if (record) E_rec <- matrix(0, n, 6)
  drift_offset <- function(k) if (length(drift_src)) sum(cfg$drift_rate * (k - drift_src + cfg$drift_k0)) else 0
  start_ep <- function(type, k) { ep <<- list(type = type, k_on = k)
    if (type == 2) sn_on <<- TRUE else if (type == 3) do_on <<- TRUE else drift_src <<- c(drift_src, k) }
  end_ep <- function(k, how) {
    if (ep$type == 2) sn_on <<- FALSE
    if (ep$type == 3) do_on <<- FALSE
    if (ep$type == 4 && how == "selfclear") drift_src <<- drift_src[drift_src != ep$k_on]
    ep_log[[length(ep_log) + 1]] <<- data.frame(type = ST[ep$type], k_on = ep$k_on, k_end = k, how = how)
    ep <<- NULL
    if (mode == "random") { gi <<- gi + 1L; next_onset <<- k + round(W$gaps[gi] / dt) }
  }
  for (k in 1:n) {
    # ---- 1. fault schedule
    if (mode == "random") {
      if (is.null(ep) && k == next_onset) start_ep(W$types[gi], k)
      if (!is.null(ep) && (k - ep$k_on) >= cfg$Dmax) end_ep(k, "selfclear")
    } else {
      for (r in seq_len(nrow(episodes))) {
        if (!is.na(episodes$k_off[r]) && k == episodes$k_off[r] + 1 && !is.null(ep) && ep$k_on == episodes$k_on[r]) end_ep(k, "selfclear")
        if (k == episodes$k_on[r]) start_ep(episodes$type[r], k)
      }
    }
    # ---- 2. twin modules (configuration set by a_prev, the action chosen at k-1)
    if (k >= 3) xin[k] <- th_p[1]*xin[k-1] + th_p[2]*xin[k-2] + th_p[3]*W$u_act_A[k-1] + th_p[4]*W$u_act_A[k-2]
    if (k >= 2) { yAi[k] <- th_a[1]*yAi[k-1] + th_a[2]*xin[k-1]; yBi[k] <- th_b[1]*yBi[k-1] + th_b[2]*xin[k-1] }
    x_rep <- if (do_on) cfg$do_scale * xin[k] else xin[k]
    nAk <- if (sn_on) cfg$sn_sdA * W$nA[k] else 0; nBk <- if (sn_on) cfg$sn_sdB * W$nB[k] else 0
    D_k <- drift_offset(k); off_err <- D_k - corr_int; off_seq[k] <- off_err
    yA_raw <- yAi[k] + nAk + off_err
    yB_raw <- yBi[k] + nBk
    filt_on <- (a_prev == 4L)
    if (filt_on) { fA <- a_f*fA + (1-a_f)*yA_raw; fB <- a_f*fB + (1-a_f)*yB_raw
                   nfA <- a_f*nfA + (1-a_f)*nAk; nfB <- a_f*nfB + (1-a_f)*nBk
    } else { fA <- yA_raw; fB <- yB_raw; nfA <- nAk; nfB <- nBk }   # filter state tracks the raw channel when off
    alpha <- if (a_prev == 2L) cfg$alpha_dwA else if (a_prev == 3L) cfg$alpha_dwB else cfg$alpha_nom
    fused <- alpha * fA + (1 - alpha) * fB
    nf_hist[k] <- alpha * nfA + (1 - alpha) * nfB          # fault-induced component of the fused output
    vel <- if (k == 1) 0 else (x_rep - x_rep_prev) / dt
    zA <- x_rep + 0.2 * vel
    e <- c(W$u_act_A[k] - W$u_act_P[k], x_rep - W$x_pos_P[k], fA - W$yA_P[k], fB - W$yB_P[k],
           fused - W$y_fused_P[k], zA - W$z_P[k])
    e2a_hist[k] <- e[3]
    if (record) E_rec[k, ] <- e
    feat <- c(e, sum(abs(e)), sqrt(sum(e^2) / 6), if (k == 1) 0 else e[2] - e1p, if (k == 1) 0 else e[5] - e3p)
    e1p <- e[2]; e3p <- e[5]; x_rep_prev <- x_rep
    # ---- 3. causal decoder
    post <- dec_step(dec, feat, post); z <- dec$map_idx[which.max(post)]
    # ---- 4. ground-truth regime at k
    s <- if (do_on) 3L else if (sn_on) 2L else if (abs(off_err) > cfg$tol_b) 4L else 1L
    if (!is.null(ep) && ep$type == 4 && abs(off_err) <= cfg$tol_b) end_ep(k, "repaired")
    cfg_key <- paste(alpha, filt_on)
    if (cfg_key != cfg_key_prev) cfg_since <- k
    cfg_key_prev <- cfg_key
    sn_rep_now <- FALSE
    if (!is.null(ep) && ep$type == 2 && cfg_key != "0.6 FALSE" && (k - max(cfg_since, ep$k_on) + 1) >= cfg$W_sn)   # window inside the episode, one config
      sn_rep_now <- sqrt(mean(nf_hist[(k - cfg$W_sn + 1):k]^2)) <= cfg$thr_sn
    # ---- 5. controller acts on z_k (oracle on s_k); reward R(a_k, s_k)
    info <- list(episode_done = is.null(ep) && length(ep_log) > 0)
    o <- ctrl$act(cst, z, s, k, info); a <- o$a; cst <- o$st
    ret <- ret + g^(k - 1) * R[a, s]
    s_seq[k] <- s; a_seq[k] <- a; z_seq[k] <- z
    if (sn_rep_now) end_ep(k + 1, "repaired")             # noise source cleared from k+1
    # ---- 6. action progress / installs (effective from k+1)
    if (a == 5L) { rp_cnt <- rp_cnt + 1L
      if (rp_cnt >= cfg$W_rp) {
        idx <- (k - cfg$W_rp + 1):k
        th <- tryCatch(fit_arx_pipe(W$x_pos_P[idx], W$u_act_P[idx], 2L, 2L)$theta, error = function(e) NULL)
        ok <- !is.null(th) && arx_stable(th) && abs(arx_dc_gain(th) / W$G_true - 1) <= cfg$rp_gain_tol
        if (ok) { th_p <- th; rp_fits["ok"] <- rp_fits["ok"] + 1L
                  if (do_on) { if (!is.null(ep) && ep$type == 3) end_ep(k + 1, "repaired") else do_on <- FALSE }
        } else rp_fits["fail"] <- rp_fits["fail"] + 1L
        rp_cnt <- 0L } } else rp_cnt <- 0L
    if (a == 6L) { bc_cnt <- bc_cnt + 1L
      if (bc_cnt >= cfg$W_bc) {
        idx <- (k - cfg$W_bc + 1):k
        cf <- lm.fit(cbind(1, idx - k), e2a_hist[idx])$coefficients   # offset at k, slope per sample
        corr_int <- corr_int + cf[[1]]; corr_slope <- corr_slope + cf[[2]]; bc_installs <- bc_installs + 1L
        bc_cnt <- 0L } } else bc_cnt <- 0L
    corr_int <- corr_int + corr_slope                     # correction extrapolated to k+1
    a_prev <- a
  }
  if (!is.null(ep)) ep_log[[length(ep_log) + 1]] <- data.frame(type = ST[ep$type], k_on = ep$k_on, k_end = NA, how = "censored")
  eps <- if (length(ep_log)) do.call(rbind, ep_log) else
    data.frame(type = character(), k_on = integer(), k_end = integer(), how = character())
  out <- list(ret = ret, fnom = mean(s_seq == 1L), mism = mean(a_seq != pol_idx[s_seq]),
              acc_dec = mean(z_seq == s_seq), eps = eps, rp_fits = rp_fits, bc_installs = bc_installs,
              act_frac = tabulate(a_seq, 6) / n, conf = table(factor(s_seq, 1:4), factor(z_seq, 1:4)))
  if (record) { out$s <- s_seq; out$a <- a_seq; out$z <- z_seq; out$off <- off_seq; out$E <- E_rec; out$nf <- nf_hist }
  out
}

# ------------------------------------------------------------------ decoder transition smoothing
# The decoder HMM is trained on the canonical schedule Nominal->SensorNoisy->DynamicsOff->Drift, so EM learns
# exact zeros (Drift absorbing; nothing returns to Nominal from DynamicsOff/Drift).  A causal filter with zero
# transition probabilities can never report a repair.  Default: the same Laplace smoothing the paper applies to
# A_NoAction (eps = 1e-3, rows renormalised).  eps = 0 reproduces the raw EM decoder.
load_setup <- function(path, dec_eps = 1e-3) {
  S <- readRDS(path)
  S$dec$Atr_raw <- S$dec$Atr
  if (dec_eps > 0) { A <- S$dec$Atr + dec_eps; S$dec$Atr <- A / rowSums(A) }
  S
}

# ------------------------------------------------------------------ sensitivity overrides, e.g. CL_CFG="W_bc=100,W_rp=200"
cfg_override <- function(cfg, spec = Sys.getenv("CL_CFG", "")) {
  if (nzchar(spec)) for (kv in strsplit(strsplit(spec, ",")[[1]], "=")) cfg[[kv[1]]] <- as.numeric(kv[2])
  cfg
}
