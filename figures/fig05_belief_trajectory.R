#fig05_belief_trajectory.R
#Supplementary fig:belief_trajectory (POMDP belief trajectory).
#Extracted from notebooks/Fig05_POMDP_Belief.Rmd (same code and seeds); writes output/figures/fig05_belief_trajectory.pdf.
#Usage (repository root): Rscript figures/fig05_belief_trajectory.R

local({ r <- Sys.getenv("DTMDP_ROOT"); if (!nzchar(r)) { r <- getwd()
  while (!file.exists(file.path(r, "R", "paths.R")) && dirname(r) != r) r <- dirname(r) }
  source(file.path(r, "R", "paths.R")) })
suppressPackageStartupMessages(load_core())

# ---- chunk: source-formulation ----
#Source the pipeline (BaseMDP -> POMDP)
library(scales)

palette_states <- c("Nominal" = "darkred", "SensorNoisy" = "#213C51",
                    "DynamicsOff" = "#574964", "Drift" = "#5A7863")
palette_actions <- c("NoAction" = "#7EACB5", "DwSensorA" = "#355872",
                     "DwSensorB" = "#213C51", "IncFilter" = "#574964",
                     "ReidentPlant" = "#5A7863", "BiasCorrect" = "darkred")

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

# ---- chunk: belief-tracker ----
#Per-sample POMDP trajectory with the belief at every 20 ms sample (R/sim_persample.R)
pomdp_beliefs_ps <- function(T_end) {
  tr <- sim_ps(ctl_belief(alpha_choose(pomdp$Gamma, pomdp$actions)), A, O_mat, 1L, T_end,
               record_belief = TRUE)
  H <- length(tr$step_states)
  list(times = (seq_len(H) - 1) * dt, states = tr$step_states,
       actions = tr$step_actions, beliefs = tr$belief)
}

# ---- chunk: search-trajectory ----
#Search across 500 seeds for a trajectory that visits all 4 states within 10s
best_tr <- NULL; best_n <- 0
for (seed in 1:500) {
  set.seed(seed)
  tr <- pomdp_beliefs_ps(10)
  n_u <- length(unique(tr$states))
  if (n_u > best_n) { best_n <- n_u; best_tr <- tr }
  if (n_u >= 4) break
}
print(paste0("Selected seed: ", seed, " with ", best_n, " unique states, ",
             length(best_tr$times), " samples"))

# ---- chunk: plot-figure5 ----
#Panel A: Belief vector evolution
bel_df <- data.frame(time = best_tr$times,
  Nominal = best_tr$beliefs[,1], SensorNoisy = best_tr$beliefs[,2],
  DynamicsOff = best_tr$beliefs[,3], Drift = best_tr$beliefs[,4])
bel_long <- melt(bel_df, id.vars = "time", variable.name = "State", value.name = "Belief")

#Panel B: True hidden state
true_df <- data.frame(time = best_tr$times,
                      state = factor(states[best_tr$states], levels = states))

#Panel C: Actions selected
act_df <- data.frame(time = best_tr$times,
                     action = factor(best_tr$actions, levels = actions))

p5a <- ggplot(bel_long, aes(x = time, y = Belief, colour = State)) +
  geom_step(linewidth = 0.9) +
  scale_colour_manual(values = palette_states) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(tag = "A", title = "POMDP Belief Vector Evolution", x = NULL, y = "Belief b(s)") +
  theme_paper + theme(legend.position = "right")

p5b <- ggplot(true_df, aes(x = time, y = state, colour = state)) +
  geom_point(size = 0.8) +
  scale_colour_manual(values = palette_states, guide = "none") +
  labs(tag = "B", title = "True Hidden State", x = NULL, y = NULL) +
  theme_paper + theme(legend.position = "none")

p5c <- ggplot(act_df, aes(x = time, y = action, colour = action)) +
  geom_point(size = 0.8) +
  scale_colour_manual(values = palette_actions, guide = "none") +
  labs(tag = "C", title = "Actions Selected", x = "Time (s)", y = NULL) +
  theme_paper + theme(legend.position = "none")

pdf(out_path("figures", "fig05_belief_trajectory.pdf"), width = 10, height = 9)
grid.arrange(p5a, p5b, p5c, ncol = 1, heights = c(2.5, 1.2, 1.2))
invisible(dev.off())
cat("  saved output/figures/fig05_belief_trajectory.pdf\n")
