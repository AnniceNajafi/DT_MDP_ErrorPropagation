#Regenerate the paper figures that depend on MCDA, using real TOPSIS.
#Outputs (output/figures/):
#  fig02_state_occupancy.pdf       Supplementary fig:state_occupancy
#  fig07_discount_factor.pdf       Supplementary fig:discount_factor
#The Fig03, Fig04 and Fig06 blocks that used to be here are superseded by
#figures/fig03_reward_comparison.R, figures/fig04_sensitivity_analysis.R and
#figures/fig06_action_mismatch.R (which use the policy-summary trajectories).
#Usage (repository root): Rscript figures/regen_figures.R

local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
suppressPackageStartupMessages(load_core())
library(scales)

fig_dir <- dirname(out_path("figures", "x"))

palette_states <- c("Nominal" = "darkred", "SensorNoisy" = "#213C51",
                    "DynamicsOff" = "#574964", "Drift" = "#5A7863")
palette_policy <- c("Optimal MDP" = "#213C51", "POMDP (PBVI)" = "#574964",
                    "MCDA (TOPSIS)" = "#5A7863", "No Intervention" = "darkred")

theme_paper <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  plot.title = element_text(hjust = 0.5, size = 20),
  axis.text = element_text(size = 15),
  text = element_text(size = 18),
  legend.position = "bottom",
  legend.text = element_text(size = 14)
)

safe_frac_in_nominal <- function(tr, T_end) {
  n_ev <- length(tr$times)
  if (n_ev == 0) return(NA_real_)
  if (n_ev == 1) return(if (tr$states[1] == 1) 1.0 else 0.0)
  occ <- 0
  for (i in 1:(n_ev - 1))
    if (tr$states[i] == 1) occ <- occ + (tr$times[i + 1] - tr$times[i])
  dt_last <- T_end - tr$times[n_ev]
  if (tr$states[n_ev] == 1) occ <- occ + dt_last
  occ / T_end
}

#=====================================================================
#Fig02: State-occupancy evolution (fig02_state_occupancy.pdf, 10 x 8)
#=====================================================================
cat("\n--- Fig02 state occupancy ---\n")
set.seed(42)
N_occ <- nq(2000, 200); t_bins <- seq(0, 15, by = 0.1)
occ_data <- data.frame()

for (pol_name in c("Optimal MDP", "POMDP (PBVI)", "MCDA (TOPSIS)", "No Intervention")) {
  cat("  policy:", pol_name, "\n")
  state_counts <- matrix(0, length(t_bins), K)
  for (i in 1:N_occ) {
    if (pol_name == "Optimal MDP") {
      tr <- gillespie_mdp(Q_rate, vi$policy, c_rate, 1, 15, r_disc)
    } else if (pol_name == "POMDP (PBVI)") {
      tr <- gillespie_pomdp(Q_rate, A, pomdp$Gamma, pomdp$actions, O_mat, c_rate, 1, 15, r_disc)
    } else if (pol_name == "MCDA (TOPSIS)") {
      tr <- gillespie_mcda(Q_rate, mcda_ideal_action, O_mat, mcda_dwell_time, c_rate, 1, 15, r_disc)
    } else {
      tr <- gillespie_mdp(Q_rate, policy_none, c_rate, 1, 15, r_disc)
    }
    binned <- bin_trajectory(tr, t_bins, K)
    for (j in seq_along(t_bins)) state_counts[j, binned[j]] <- state_counts[j, binned[j]] + 1
  }
  state_probs <- state_counts / N_occ
  for (s in 1:K) {
    occ_data <- rbind(occ_data, data.frame(
      time = t_bins, probability = state_probs[, s],
      state = states[s], policy = pol_name))
  }
}

occ_data$policy <- factor(occ_data$policy,
  levels = c("Optimal MDP", "POMDP (PBVI)", "MCDA (TOPSIS)", "No Intervention"))
occ_data$state <- factor(occ_data$state, levels = states)

p_fig02 <- ggplot(occ_data, aes(x = time, y = probability, colour = state)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ policy, ncol = 2) +
  scale_colour_manual(values = palette_states) +
  labs(title = "State-Occupancy Evolution Under Four Intervention Strategies",
       x = "Time (s)", y = "P(state)", colour = "State") +
  theme_paper + theme(legend.position = "right")

ggsave(file.path(fig_dir, "fig02_state_occupancy.pdf"), p_fig02,
       width = 10, height = 8, device = "pdf")
cat("  saved fig02_state_occupancy.pdf\n")

#=====================================================================
#Fig07: Discount factor (fig07_discount_factor.pdf, 9 x 5.5)
#=====================================================================
cat("\n--- Fig07 discount factor ---\n")
set.seed(42)
gamma_levels <- c(0.90, 0.95, 0.99, 0.999)
N_gam <- nq(300, 30); T_gam <- 15
gamma_results <- data.frame()
for (gam in gamma_levels) {
  cat("  gamma =", gam, "\n")
  vi_g <- value_iteration(A, R_mat, gam)
  r_disc_g <- -log(gam) / dt
  c_rate_g <- R_mat / dt
  bp <- generate_belief_points(K, n_random = 80)
  pomdp_g <- pbvi_solve(A, O_mat, R_mat, gam, bp,
                        max_horizon = nq(300, 60), tol = 0.001, V_init = vi_g$V)
  policy_n <- setNames(rep("NoAction", K), states)
  for (pol_name in c("Optimal MDP", "POMDP (PBVI)", "MCDA (TOPSIS)", "No Intervention")) {
    fracs <- numeric(N_gam)
    for (i in 1:N_gam) {
      if (pol_name == "Optimal MDP") {
        tr <- gillespie_mdp(Q_rate, vi_g$policy, c_rate_g, 1, T_gam, r_disc_g)
      } else if (pol_name == "POMDP (PBVI)") {
        tr <- gillespie_pomdp(Q_rate, A, pomdp_g$Gamma, pomdp_g$actions,
                              O_mat, c_rate_g, 1, T_gam, r_disc_g)
      } else if (pol_name == "MCDA (TOPSIS)") {
        tr <- gillespie_mcda(Q_rate, vi_g$policy, O_mat, mcda_dwell_time,
                             c_rate_g, 1, T_gam, r_disc_g)
      } else {
        tr <- gillespie_mdp(Q_rate, policy_n, c_rate_g, 1, T_gam, r_disc_g)
      }
      fracs[i] <- safe_frac_in_nominal(tr, T_gam)
    }
    gamma_results <- rbind(gamma_results, data.frame(
      gamma = gam, policy = pol_name,
      mean_frac = mean(fracs, na.rm = TRUE),
      se = sd(fracs, na.rm = TRUE) / sqrt(sum(!is.na(fracs)))))
  }
}
gamma_results$policy <- factor(gamma_results$policy,
  levels = c("Optimal MDP", "POMDP (PBVI)", "MCDA (TOPSIS)", "No Intervention"))

p_fig07 <- ggplot(gamma_results, aes(x = factor(gamma), y = mean_frac, colour = policy, group = policy)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_frac - 1.96 * se, ymax = mean_frac + 1.96 * se), width = 0.1) +
  scale_colour_manual(values = palette_policy) +
  scale_y_continuous(labels = percent) +
  labs(title = "Fraction in Nominal Across Discount Factors",
       x = expression(gamma), y = "Mean Fraction in Nominal", colour = "Policy") +
  theme_paper

ggsave(file.path(fig_dir, "fig07_discount_factor.pdf"), p_fig07,
       width = 9, height = 5.5, device = "pdf")
cat("  saved fig07_discount_factor.pdf\n")

cat("\n=== All affected figures regenerated ===\n")
cat("Output directory:", fig_dir, "\n")
