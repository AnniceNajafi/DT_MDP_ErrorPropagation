#fig04_sensitivity_analysis.R: Fig. fig:sensitivity. Panels A and B rerun the sweeps of the
#earlier figures/regen_figures.R (same seeds and settings: 300 trajectories, T = 15 s, PBVI with
#80 random belief points and 300 iterations); panel C uses the Table 5 means
#(analyses/table5_final.R, 2,000 trajectories, T = 20 s; run it first). Larger fonts and a shared legend.
#Usage (repository root): Rscript figures/fig04_sensitivity_analysis.R
local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
fig_out  <- out_path("figures", "fig04_sensitivity_analysis.pdf")
cache    <- out_path("analyses", "fig04_sweeps.rds")
# The sweeps are re-run unless DTMDP_REUSE_SIM=1 and the cache exists.
if (!identical(Sys.getenv("DTMDP_REUSE_SIM"), "1") && file.exists(cache)) invisible(file.remove(cache))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(scales) })
if (!file.exists(cache)) {
  suppressPackageStartupMessages(load_core())
#Panel A: SensorNoisy classification sweep
set.seed(42)
p_sn_levels <- if (QUICK) seq(0.20, 1.00, by = 0.20) else seq(0.20, 1.00, by = 0.05)
N_sweep <- nq(300, 30); T_sweep <- 15
obs_results <- data.frame()

for (p_sn in p_sn_levels) {
  cat("  p_SN =", p_sn, "\n")
  O_mod <- O_mat
  O_mod[2, ] <- (1 - p_sn) / (K - 1); O_mod[2, 2] <- p_sn
  bp <- generate_belief_points(K, n_random = 80)
  pomdp_mod <- pbvi_solve(A, O_mod, R_mat, gamma_disc, bp,
                          max_horizon = nq(300, 60), tol = 0.001, V_init = vi$V)
  for (pol_name in c("Optimal MDP", "POMDP (PBVI)", "MCDA (TOPSIS)", "No Intervention")) {
    rews <- numeric(N_sweep)
    for (i in 1:N_sweep) {
      if (pol_name == "Optimal MDP") {
        tr <- gillespie_mdp(Q_rate, vi$policy, c_rate, 1, T_sweep, r_disc)
      } else if (pol_name == "POMDP (PBVI)") {
        tr <- gillespie_pomdp(Q_rate, A, pomdp_mod$Gamma, pomdp_mod$actions,
                              O_mod, c_rate, 1, T_sweep, r_disc)
      } else if (pol_name == "MCDA (TOPSIS)") {
        tr <- gillespie_mcda(Q_rate, mcda_ideal_action, O_mod, mcda_dwell_time,
                             c_rate, 1, T_sweep, r_disc)
      } else {
        tr <- gillespie_mdp(Q_rate, policy_none, c_rate, 1, T_sweep, r_disc)
      }
      rews[i] <- tr$disc_return
    }
    obs_results <- rbind(obs_results, data.frame(
      p_sn = p_sn, policy = pol_name,
      mean_reward = mean(rews), se = sd(rews) / sqrt(N_sweep)))
  }
}
obs_results$policy <- factor(obs_results$policy,
  levels = c("Optimal MDP", "POMDP (PBVI)", "MCDA (TOPSIS)", "No Intervention"))

#Panel B: Repair-probability sweep
set.seed(42)
rho_base <- matrix(0, nrow = M, ncol = K, dimnames = list(actions, states))
rho_base["DwSensorA", "SensorNoisy"] <- 0.60
rho_base["DwSensorB", "SensorNoisy"] <- 0.60
rho_base["IncFilter", "SensorNoisy"] <- 0.40
rho_base["ReidentPlant", "DynamicsOff"] <- 0.70
rho_base["BiasCorrect", "Drift"] <- 0.80

lambda_levels <- if (QUICK) c(0.1, 0.5, 1.0) else seq(0.1, 1.0, by = 0.1)
N_rep <- nq(300, 30); T_rep <- 15
repair_results <- data.frame()

for (lam in lambda_levels) {
  cat("  lambda =", lam, "\n")
  rho_s <- rho_base * lam
  A_new <- list()
  for (a in actions) {
    A_new[[a]] <- A_NoAction
    for (s in 1:K) {
      if (rho_s[a, s] > 0) {
        r <- rho_s[a, s]
        A_new[[a]][s, ] <- (1 - r) * A_NoAction[s, ]
        A_new[[a]][s, 1] <- A_new[[a]][s, 1] + r
      }
    }
  }
  Q_new <- lapply(A_new, function(Aa) {
    Qa <- (Aa - diag(K)) / dt; dimnames(Qa) <- list(states, states); Qa
  })
  vi_new <- value_iteration(A_new, R_mat, gamma_disc)
  c_rate_new <- R_mat / dt
  bp <- generate_belief_points(K, n_random = 80)
  pomdp_new <- pbvi_solve(A_new, O_mat, R_mat, gamma_disc, bp,
                          max_horizon = nq(300, 60), tol = 0.001, V_init = vi_new$V)
  for (pol_name in c("Optimal MDP", "POMDP (PBVI)", "MCDA (TOPSIS)", "No Intervention")) {
    rews <- numeric(N_rep)
    for (i in 1:N_rep) {
      if (pol_name == "Optimal MDP") {
        tr <- gillespie_mdp(Q_new, vi_new$policy, c_rate_new, 1, T_rep, r_disc)
      } else if (pol_name == "POMDP (PBVI)") {
        tr <- gillespie_pomdp(Q_new, A_new, pomdp_new$Gamma, pomdp_new$actions,
                              O_mat, c_rate_new, 1, T_rep, r_disc)
      } else if (pol_name == "MCDA (TOPSIS)") {
        tr <- gillespie_mcda(Q_new, vi_new$policy, O_mat, mcda_dwell_time,
                             c_rate_new, 1, T_rep, r_disc)
      } else {
        policy_n <- setNames(rep("NoAction", K), states)
        tr <- gillespie_mdp(Q_new, policy_n, c_rate_new, 1, T_rep, r_disc)
      }
      rews[i] <- tr$disc_return
    }
    repair_results <- rbind(repair_results, data.frame(
      lambda = lam, policy = pol_name,
      mean_reward = mean(rews), se = sd(rews) / sqrt(N_rep)))
  }
}
repair_results$policy <- factor(repair_results$policy,
  levels = c("Optimal MDP", "POMDP (PBVI)", "MCDA (TOPSIS)", "No Intervention"))
  saveRDS(list(obs_results = obs_results, repair_results = repair_results), cache)
}
sw <- readRDS(cache); obs_results <- sw$obs_results; repair_results <- sw$repair_results
pol <- c("Optimal MDP", "POMDP (PBVI)", "MCDA (TOPSIS)", "No Intervention")
palette_policy <- c("Optimal MDP" = "#213C51", "POMDP (PBVI)" = "#574964",
                    "MCDA (TOPSIS)" = "#5A7863", "No Intervention" = "darkred")
tb <- read.csv(out_path("table5_final", "PolicySummaryTable.csv"))
voi_data <- data.frame(policy = factor(pol, levels = pol),
                       mean_ret = tb$MeanReward[match(pol, tb$Policy)],
                       se = tb$SEReward[match(pol, tb$Policy)])
theme_paper <- theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"),
  plot.title = element_text(hjust = 0.5, size = 14),
  plot.tag = element_text(size = 15, face = "bold"),
  axis.text = element_text(size = 11.5, colour = "black"), axis.title = element_text(size = 12.5),
  text = element_text(size = 12.5), legend.position = "bottom",
  legend.title = element_text(size = 12.5), legend.text = element_text(size = 12),
  legend.key.width = unit(1.1, "cm"))
p4a <- ggplot(obs_results, aes(p_sn, mean_reward, colour = policy, fill = policy)) +
  geom_ribbon(aes(ymin = mean_reward - 1.96 * se, ymax = mean_reward + 1.96 * se), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) + geom_point(size = 1.6) +
  scale_colour_manual(values = palette_policy) + scale_fill_manual(values = palette_policy) +
  scale_x_continuous(breaks = seq(0.2, 1.0, 0.2), labels = percent) +
  labs(tag = "A", title = "Observation quality", x = expression(p[SN] ~ "(SensorNoisy accuracy)"),
       y = "Mean discounted return", colour = "Policy", fill = "Policy") + theme_paper
p4b <- ggplot(repair_results, aes(lambda, mean_reward, colour = policy, fill = policy)) +
  geom_ribbon(aes(ymin = mean_reward - 1.96 * se, ymax = mean_reward + 1.96 * se), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) + geom_point(size = 1.6) +
  scale_colour_manual(values = palette_policy) + scale_fill_manual(values = palette_policy) +
  scale_x_continuous(breaks = seq(0.2, 1.0, 0.2)) +
  labs(tag = "B", title = "Repair probability scale", x = expression(lambda ~ "(scale on" ~ rho[a*","*s] * ")"),
       y = "Mean discounted return", colour = "Policy", fill = "Policy") + theme_paper
p4c <- ggplot(voi_data, aes(policy, mean_ret, fill = policy)) +
  geom_col(width = 0.6, alpha = 0.9) +
  geom_errorbar(aes(ymin = mean_ret - 1.96 * se, ymax = mean_ret + 1.96 * se), width = 0.15) +
  geom_text(aes(label = sprintf("%.1f", mean_ret)), vjust = -0.9, size = 4.3, fontface = "bold") +
  scale_fill_manual(values = palette_policy, guide = "none") +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(tag = "C", title = "Mean discounted return (Table 5 runs)", x = NULL, y = "Mean discounted return") +
  theme_paper
p <- ((p4a | p4b) / p4c) + plot_layout(heights = c(1, 0.85), guides = "collect") &
  theme(legend.position = "bottom")
ggsave(fig_out, p, width = 10, height = 8, device = if (capabilities("cairo")) cairo_pdf else "pdf")
cat("saved", fig_out, "
")
print(obs_results); print(repair_results)
