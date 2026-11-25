#' mdp_validation.R
#' ──────────────────────────────────────────────────────────────────────
#' Monte Carlo validation of the continuous-time Bellman solution by the
#' Gillespie simulator (formerly sections 6-7 of BaseMDPFormulation.R).
#'
#' Outputs (output/):
#'   tables/MDP_Bellman_vs_Gillespie.csv
#'   figures/validation/SingleTrajectory.pdf
#'   figures/validation/ManyTrajectories_vs_Bellman.pdf
#'   figures/validation/CumulativeReward_AllStates.pdf
#'
#' Authors@R: person("Annice", "Najafi")
#' ──────────────────────────────────────────────────────────────────────

local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
source(repo_path("R", "BaseMDPFormulation.R"))
set.seed(42)

# =====================================================================
# 6. MONTE CARLO VALIDATION  (Gillespie vs Bellman)
# =====================================================================

N_sim <- nq(3000, 300)
cat(sprintf("\n── Monte Carlo validation (%d trajectories per state) ──\n", N_sim))
T_end <- 20

V_sim <- matrix(NA, N_sim, K, dimnames = list(NULL, states))
for (s0 in 1:K) {
  for (n in 1:N_sim) {
    tr <- gillespie_mdp(Q_rate, vi$policy, c_rate, s0, T_end, r_disc)
    V_sim[n, s0] <- tr$disc_return
  }
}
V_hat <- colMeans(V_sim)
V_se  <- apply(V_sim, 2, sd) / sqrt(N_sim)

comparison <- data.frame(
  State    = states,
  Bellman  = round(V_ct, 4),
  MC_Mean  = round(V_hat, 4),
  StdErr   = round(V_se, 4),
  RelError = sprintf("%.2f%%", 100 * abs(V_hat - V_ct) / abs(V_ct))
)
print(comparison)
write.csv(comparison, out_path("tables", "MDP_Bellman_vs_Gillespie.csv"), row.names = FALSE)

# =====================================================================
# 7. FIGURES
# =====================================================================

# ------------------------------------------------------------------
# Figure 1: Single trajectory from Nominal  (3-panel)
# ------------------------------------------------------------------
cat("\n── Figure 1: Single trajectory ──\n")

set.seed(314)
tr_ex <- gillespie_mdp(Q_rate, vi$policy, c_rate, 1, 10.0, r_disc)

t_plot <- seq(0, 10, by = 0.005)
st_at_t    <- integer(length(t_plot))
rew_at_t   <- numeric(length(t_plot))
cum_at_t   <- numeric(length(t_plot))
cr <- 0;  prev_t <- 0;  prev_s <- tr_ex$states[1];  ti <- 1

for (i in seq_along(t_plot)) {
  while (ti < length(tr_ex$times) && tr_ex$times[ti + 1] <= t_plot[i]) {
    seg <- tr_ex$times[ti + 1] - prev_t
    a_i <- vi$policy[states[prev_s]]
    ci  <- c_rate[a_i, prev_s]
    cr  <- cr + (ci / r_disc) *
           exp(-r_disc * prev_t) * (1 - exp(-r_disc * seg))
    prev_t <- tr_ex$times[ti + 1];  ti <- ti + 1;  prev_s <- tr_ex$states[ti]
  }
  seg <- t_plot[i] - prev_t
  a_now <- vi$policy[states[prev_s]]
  rew_at_t[i] <- c_rate[a_now, prev_s]
  st_at_t[i]  <- prev_s
  if (seg > 0 && r_disc > 0) {
    cum_at_t[i] <- cr + (rew_at_t[i] / r_disc) *
      exp(-r_disc * prev_t) * (1 - exp(-r_disc * seg))
  } else {
    cum_at_t[i] <- cr
  }
}

df_ex <- data.frame(
  time   = t_plot,
  State  = factor(states[st_at_t], levels = states),
  Reward_Rate = rew_at_t,
  Cumulative  = cum_at_t
)

p1 <- ggplot(df_ex, aes(time, State, group = 1)) +
  geom_step(color = "#2c3e50", linewidth = 1) +
  labs(title = "Single Trajectory Starting from Nominal Under Optimal Policy",
       subtitle = "Top: which state the system is in over time",
       x = NULL, y = "State") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_blank())

p2 <- ggplot(df_ex, aes(time, Reward_Rate, group = 1)) +
  geom_step(color = "#e74c3c", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  labs(subtitle = "Middle: instantaneous reward rate c(s, pi*(s))",
       x = NULL, y = "Reward rate\n(per second)") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_blank())

p3 <- ggplot(df_ex, aes(time, Cumulative)) +
  geom_line(color = "coral", linewidth = 1) +
  geom_hline(yintercept = V_ct["Nominal"], linetype = "dashed",
             color = "steelblue", linewidth = 1) +
  annotate("text", x = 9, y = V_ct["Nominal"],
           label = paste0("V*(Nominal) = ", round(V_ct["Nominal"], 2)),
           vjust = -0.7, color = "steelblue", size = 4.5) +
  labs(subtitle = "Bottom: cumulative discounted reward",
       x = "Time (s)", y = "Cumulative\ndiscounted reward") +
  theme_minimal(base_size = 13)

pdf(out_path("figures", "validation", "SingleTrajectory.pdf"), width = 10, height = 9)
gridExtra::grid.arrange(p1, p2, p3, ncol = 1, heights = c(1, 0.8, 1))
dev.off()
cat("  -> Saved SingleTrajectory.pdf\n")

# ------------------------------------------------------------------
# Figure 2: Many trajectories vs Bellman (MC mean + theory)
# ------------------------------------------------------------------
cat("── Figure 2: Many trajectories vs Bellman ──\n")

p0_nom  <- c(1, 0, 0, 0)
t_fine  <- seq(0, 10, by = 0.01)
cum_th_nom <- numeric(length(t_fine))
for (i in 2:length(t_fine)) {
  dt_i <- t_fine[i] - t_fine[i - 1]
  p_t  <- as.numeric(p0_nom %*% expm(Q_pi * t_fine[i]))
  cum_th_nom[i] <- cum_th_nom[i - 1] +
    exp(-r_disc * t_fine[i]) * sum(p_t * R_pi_vec) * dt_i
}

N_thin <- 30
N_mean <- nq(2000, 200)
t_pts  <- seq(0, 10, by = 0.1)

set.seed(2025)
df_thin <- do.call(rbind, lapply(1:N_thin, function(n) {
  tr <- gillespie_mdp(Q_rate, vi$policy, c_rate, 1, 10, r_disc)
  data.frame(time = t_pts,
             Cumulative = cum_reward_at_pts(tr, t_pts, vi$policy, c_rate, r_disc),
             id = factor(n))
}))

cum_all <- matrix(NA, N_mean, length(t_pts))
for (n in 1:N_mean) {
  tr <- gillespie_mdp(Q_rate, vi$policy, c_rate, 1, 10, r_disc)
  cum_all[n, ] <- cum_reward_at_pts(tr, t_pts, vi$policy, c_rate, r_disc)
}
mu_nom  <- colMeans(cum_all)
q25_nom <- apply(cum_all, 2, quantile, 0.25)
q75_nom <- apply(cum_all, 2, quantile, 0.75)

th_nom_interp <- approx(t_fine, cum_th_nom, xout = t_pts)$y

fig2 <- ggplot() +
  #geom_line(data = df_thin,
  #          aes(time, Cumulative, group = id),
  #          color = "grey70", linewidth = 0.3, alpha = 0.5) +
  geom_point(aes(t_pts, mu_nom, color = "Gillespie MC mean"),
            stroke=1.5, shape=8) +
  geom_line(aes(t_pts, th_nom_interp, color = "Bellman (theoretical)"),
            linewidth = 1.5, linetype = "dashed") +
  geom_hline(yintercept = V_ct["Nominal"], linetype = "dotted",
             color = "steelblue", linewidth = 0.7) +
  labs(title = "Cumulative Discounted Reward - System Starting from Nominal",
       x = "Time (s)", y = "Cumulative Discounted Reward", color = NULL) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("Bellman (theoretical)" = "steelblue",
                                "Gillespie MC mean" = "#DDAED3"))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text = element_text(size = 15),
    text = element_text(size = 18)
  )

pdf(out_path("figures", "validation", "ManyTrajectories_vs_Bellman.pdf"), width = 10, height = 6)
print(fig2)
dev.off()
cat("  -> Saved ManyTrajectories_vs_Bellman.pdf\n")

# ------------------------------------------------------------------
# Figure 3: Cumulative reward by starting state (theory vs MC)
# ------------------------------------------------------------------
cat("── Figure 3: All starting states ──\n")

cum_th_all <- list()
for (s0 in 1:K) {
  p0_s <- numeric(K);  p0_s[s0] <- 1
  cum_s <- numeric(length(t_fine))
  for (i in 2:length(t_fine)) {
    dt_i <- t_fine[i] - t_fine[i - 1]
    p_t  <- as.numeric(p0_s %*% expm(Q_pi * t_fine[i]))
    cum_s[i] <- cum_s[i - 1] +
      exp(-r_disc * t_fine[i]) * sum(p_t * R_pi_vec) * dt_i
  }
  cum_th_all[[s0]] <- cum_s
}

N_per <- nq(1500, 150)
cum_mc_all <- list()
for (s0 in 1:K) {
  mat <- matrix(NA, N_per, length(t_pts))
  for (n in 1:N_per) {
    tr <- gillespie_mdp(Q_rate, vi$policy, c_rate, s0, 10, r_disc)
    mat[n, ] <- cum_reward_at_pts(tr, t_pts, vi$policy, c_rate, r_disc)
  }
  cum_mc_all[[s0]] <- colMeans(mat)
}

df_theory_all <- do.call(rbind, lapply(1:K, function(s0) {
  data.frame(time = t_fine, Cumulative = cum_th_all[[s0]], State = states[s0])
}))
df_mc_all <- do.call(rbind, lapply(1:K, function(s0) {
  data.frame(time = t_pts, Cumulative = cum_mc_all[[s0]], State = states[s0])
}))

fig3 <- ggplot() +
  geom_line(data = df_theory_all,
            aes(time, Cumulative, color = State),
            linewidth = 1.2, linetype = "dashed") +
  geom_line(data = df_mc_all,
            aes(time, Cumulative, color = State),
            linewidth = 0.9, alpha = 0.8) +
  labs(title = "Cumulative Discounted Reward by Starting State",
       subtitle = "Dashed = Bellman theory  |  Solid = Gillespie MC mean",
       x = "Time (s)", y = "Cumulative Discounted Reward") +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Set1")

pdf(out_path("figures", "validation", "CumulativeReward_AllStates.pdf"), width = 10, height = 6)
print(fig3)
dev.off()
cat("  -> Saved CumulativeReward_AllStates.pdf\n")

cat("\n══ mdp_validation complete. ══\n")
